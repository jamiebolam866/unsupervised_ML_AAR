# unsupervivsed ML AAR utils
# restarting the functions for our acc workflows, based on tidy data formats

# functions to segment continuous ACC data into bouts
#  (and identify distinct bouts based on behaviours if labelled data provided)

# segment_bouts.R --------------------------------------------------------------
# Turns a raw, unsegmented sensor stream (one row per sample, a real `time`
# column, id, optionally a behaviour `label`) into the same bout-tsibble shape
# as vultures_tsbl: key = c(id, bout), index = sample_in_bout. Two functions:
#
#   add_segments()  finds CONTIGUOUS RUNS per individual: a new run starts
#                   whenever there's a time gap of at least `min_gap_s`, or
#                   (if a label column is used) the label changes -- including
#                   to/from NA, which counts as its own run. Runs are whatever
#                   length they naturally are; nothing is fixed-length yet.
#                   Adds a `segment` column (restarts at 1 per individual, for
#                   readability when used standalone) and returns the data
#                   GROUPED by (id, segment).
#
#   make_bouts()    chops each contiguous run into FIXED-length bouts (`nrow`
#                   samples, or `duration_s` seconds via `fs`), drops any
#                   leftover partial bout, and returns a tsibble with `bout`
#                   (ONE running count across the whole output), `sample_in_bout`
#                   (position within the bout), `t_in_bout` (actual elapsed
#                   time since the bout's first sample -- not an idealised
#                   1/fs grid, since real timestamps are available), and a
#                   `segment` column (which contiguous run the bout came from
#                   -- the same leakage-avoidance signal we found in the
#                   vultures data: bouts from the same segment are
#                   near-duplicates, so split train/test by segment or
#                   individual, not by bout).
#
# Respects whatever grouping the data already carries: make_bouts() always
# adds `id_col` to the grouping if it's not already there, but otherwise uses
# exactly what it's given. This is how label-awareness is opted into, with no
# `label` argument on make_bouts() itself and no cost to using it on
# unlabelled data:
#   raw |> add_segments(min_gap_s = 5, label_col = "label") |> make_bouts(...)
#     -- bouts never straddle a label change (or a gap)
#   raw |> make_bouts(min_gap_s = 5, ...)
#     -- label ignored entirely; bouts only respect id and time gaps
# make_bouts() ALWAYS does its own gap check too (using its own `min_gap_s`,
# which need not match add_segments()'s), so it's safe even called directly,
# without add_segments(), on a stream that has real gaps in it.
#
# `time` can be POSIXct (full date-time), hms/difftime (time-of-day only, as
# in the cats dataset), or plain numeric seconds. A POSIXct column handles a
# recording that spans midnight correctly with no special handling needed --
# it's just a continuous timeline, midnight isn't special to it. A bare
# time-of-day (hms) column cannot resolve that (there is no date to say which
# day), so if time appears to run backwards, that's an error asking for a full
# date-time column instead, rather than silently guessing.
#
# Requires: dplyr (>= 1.1.0, for .by/pick), vctrs, rlang, tsibble.

# Per-individual (or per-existing-group) run id: starts a new run wherever the
# gap to the previous sample is >= min_gap_s, or (if label given) label changes.
# Assumes time_s is already known to be non-decreasing (see .enforce_time_order()
# below, which checks this BEFORE any sorting happens) -- the check here is
# just defense in depth for direct callers.
.mark_runs <- function(time_s, min_gap_s, label = NULL) {
    n <- length(time_s)
    if (n == 0) {
        return(integer(0))
    }
    dt <- diff(time_s)
    if (any(dt < 0)) {
        rlang::abort(
            "Time goes backwards -- .mark_runs() requires already-sorted, monotonic input."
        )
    }
    run_break <- c(TRUE, dt >= min_gap_s)
    if (!is.null(label)) {
        run_break <- run_break |
            c(FALSE, diff(dplyr::consecutive_id(label)) > 0)
    }
    cumsum(run_break)
}

# Seconds, from POSIXct, hms/difftime, or plain numeric.
.time_to_seconds <- function(x) {
    if (inherits(x, "POSIXct")) {
        return(as.numeric(x))
    }
    if (inherits(x, "difftime")) {
        return(as.numeric(x, units = "secs"))
    }
    if (is.numeric(x)) {
        return(x)
    }
    rlang::abort(
        "`time` must be POSIXct, hms/difftime, or plain numeric seconds."
    )
}

# Checks monotonicity of time WITHIN each id_col group, in the ROW ORDER
# `data` WAS GIVEN IN -- crucially, BEFORE any sorting, since sorting would
# silently "fix" the exact problem this is meant to catch (and would be
# actively wrong for a genuine day-rollover in a bare time-of-day column: it
# would put the post-midnight rows before the pre-midnight rows of the same
# continuous recording). Errors clearly if time ever decreases; otherwise
# returns `data`, ungrouped, sorted by (id_col, time_col), with a `..t_s..`
# column (time in seconds) added for callers to reuse.
.enforce_time_order <- function(data, id_col, time_col) {
    data <- dplyr::ungroup(data)
    data$..t_s.. <- .time_to_seconds(data[[time_col]])
    dt <- data |>
        dplyr::mutate(
            ..dt.. = ..t_s.. - dplyr::lag(..t_s..),
            .by = dplyr::all_of(id_col)
        ) |>
        dplyr::pull(..dt..)
    if (any(dt < 0, na.rm = TRUE)) {
        rlang::abort(c(
            "Time goes backwards within an individual's data (checked in the row order given, before any sorting).",
            i = "Rows must already be in increasing time order within each individual.",
            i = "If this is a genuine day-boundary crossing, supply a full date-time (POSIXct) column, not a bare time-of-day (hms) column -- a bare time-of-day has no way to tell a real gap from a day rollover."
        ))
    }
    dplyr::arrange(data, dplyr::across(dplyr::all_of(c(id_col, time_col))))
}

.check_min_gap_s <- function(min_gap_s) {
    if (
        missing(min_gap_s) ||
            is.null(min_gap_s) ||
            !is.numeric(min_gap_s) ||
            length(min_gap_s) != 1 ||
            is.na(min_gap_s) ||
            min_gap_s <= 0
    ) {
        rlang::abort(c(
            "`min_gap_s` is required.",
            i = "The minimum time gap (seconds) between consecutive samples that marks a new run -- any gap this size or larger starts a new one."
        ))
    }
}

#' Find contiguous recording runs per individual
#'
#' @param data      A data frame/tibble with (at least) `id_col` and `time_col`.
#' @param min_gap_s Required. Minimum time gap, in seconds, between consecutive
#'                  samples that counts as a break between runs.
#' @param id_col    Name of the individual-identifier column.
#' @param time_col  Name of the time column (POSIXct, hms/difftime, or numeric
#'                  seconds).
#' @param label_col Name of a behaviour-label column. If present, a label
#'                  change (including to/from NA) also breaks a run. Set to
#'                  NULL to ignore label even if the column exists; if the
#'                  named column doesn't exist, it's skipped automatically
#'                  (with a message) -- segmentation still works, on time
#'                  gaps alone.
#' @param verbose   Print progress messages.
#'
#' @return `data`, sorted by (id_col, time_col), with a `segment` column
#'   (restarts at 1 per individual), returned GROUPED by (id_col, segment) so
#'   it flows straight into make_bouts().
#' @examples
#' cats |> add_segments(min_gap_s = 5, label_col = "label")
add_segments <- function(
    data,
    min_gap_s,
    id_col = "id",
    time_col = "time",
    label_col = "label",
    verbose = TRUE
) {
    say <- function(...) if (verbose) rlang::inform(paste0(...))

    if (!is.data.frame(data)) {
        rlang::abort("`data` must be a data frame/tibble.")
    }
    .check_min_gap_s(min_gap_s)
    missing_cols <- setdiff(c(id_col, time_col), names(data))
    if (length(missing_cols)) {
        rlang::abort(paste0(
            "Column(s) not found: ",
            paste(missing_cols, collapse = ", ")
        ))
    }
    if ("segment" %in% names(data)) {
        rlang::abort(
            "`data` already has a `segment` column; rename or remove it first."
        )
    }

    use_label <- !is.null(label_col) && label_col %in% names(data)
    if (!is.null(label_col) && !use_label) {
        say("No `", label_col, "` column found; segmenting on time gaps only.")
    }

    df <- .enforce_time_order(data, id_col, time_col)

    df <- df |>
        dplyr::mutate(
            segment = .mark_runs(
                ..t_s..,
                min_gap_s,
                if (use_label) .data[[label_col]] else NULL
            ),
            .by = dplyr::all_of(id_col)
        ) |>
        dplyr::select(-..t_s..)

    n_seg <- dplyr::n_distinct(df[c(id_col, "segment")])
    say(
        format(n_seg, big.mark = ","),
        " contiguous run(s) found across ",
        format(dplyr::n_distinct(df[id_col]), big.mark = ","),
        " individual(s)."
    )

    dplyr::group_by(df, dplyr::across(dplyr::all_of(c(id_col, "segment"))))
}

#' Chop contiguous runs into fixed-length bouts, as a tsibble
#'
#' Respects any existing grouping on `data` (e.g. from add_segments()); always
#' adds `id_col` to the grouping if it isn't already there. Also always does
#' its own time-gap check within each such group, so it's safe to call
#' directly (without add_segments()) on a stream that may contain real gaps.
#'
#' @param data       A data frame/tibble, optionally grouped (e.g. by
#'                   add_segments()). Must have `id_col` and `time_col`.
#' @param min_gap_s  Required. As in add_segments() -- this function's OWN
#'                   safety check; it need not match a `min_gap_s` used
#'                   upstream in add_segments().
#' @param nrow       Bout length in samples. Give this or `duration_s`, not
#'                   both.
#' @param duration_s Bout length in seconds. Needs `fs`; converted to samples
#'                   as `round(duration_s * fs)`.
#' @param fs         Sampling rate in Hz, required with `duration_s`.
#' @param id_col     Name of the individual-identifier column.
#' @param time_col   Name of the time column (POSIXct, hms/difftime, or
#'                   numeric seconds). Kept in the output alongside the
#'                   derived `t_in_bout`.
#' @param verbose    Print progress messages.
#'
#' @return A tsibble, key = c(id_col, "bout"), index = "sample_in_bout", with
#'   columns id, bout (one running count across the WHOLE output),
#'   sample_in_bout, t_in_bout (actual elapsed time since the bout's first
#'   sample), segment (which contiguous run the bout came from), the original
#'   `time_col`, and every other original column carried through unchanged.
#'   Any leftover partial bout at the end of a run is dropped.
#' @examples
#' cats |> add_segments(min_gap_s = 5, label_col = "label") |>
#'   make_bouts(min_gap_s = 5, duration_s = 2, fs = 25)
#' cats |> make_bouts(min_gap_s = 5, nrow = 100)   # label ignored entirely
make_bouts <- function(
    data,
    min_gap_s,
    nrow = NULL,
    duration_s = NULL,
    fs = NULL,
    id_col = "id",
    time_col = "time",
    verbose = TRUE
) {
    say <- function(...) if (verbose) rlang::inform(paste0(...))

    if (!is.data.frame(data)) {
        rlang::abort("`data` must be a data frame/tibble.")
    }
    .check_min_gap_s(min_gap_s)
    missing_cols <- setdiff(c(id_col, time_col), names(data))
    if (length(missing_cols)) {
        rlang::abort(paste0(
            "Column(s) not found: ",
            paste(missing_cols, collapse = ", ")
        ))
    }
    # `segment` is allowed to already exist (it's how add_segments() output flows
    # in) and is recomputed/overwritten below; the others never are.
    clash <- intersect(c("bout", "sample_in_bout", "t_in_bout"), names(data))
    if (length(clash)) {
        rlang::abort(paste0(
            "`data` already has column(s) this function creates: ",
            paste(clash, collapse = ", "),
            "; rename or remove them first."
        ))
    }

    if (is.null(nrow) && is.null(duration_s)) {
        rlang::abort(
            "Give either `nrow` (samples per bout) or `duration_s` (+ `fs`)."
        )
    }
    if (!is.null(nrow) && !is.null(duration_s)) {
        rlang::abort("Give either `nrow` or `duration_s`, not both.")
    }
    if (!is.null(duration_s)) {
        if (
            is.null(fs) ||
                !is.numeric(fs) ||
                length(fs) != 1 ||
                is.na(fs) ||
                fs <= 0
        ) {
            rlang::abort("`duration_s` needs `fs`, the sampling rate in Hz.")
        }
        target_n <- round(duration_s * fs)
        say(
            "Bout length: ",
            target_n,
            " samples (",
            round(target_n / fs, 2),
            "s at ",
            fs,
            " Hz; requested ",
            duration_s,
            "s)."
        )
    } else {
        if (
            !is.numeric(nrow) ||
                length(nrow) != 1 ||
                is.na(nrow) ||
                nrow != round(nrow)
        ) {
            rlang::abort("`nrow` must be a single whole number of samples.")
        }
        target_n <- nrow
        say("Bout length: ", target_n, " samples.")
    }
    if (target_n < 2) {
        rlang::abort("The bout length must resolve to at least 2 samples.")
    }
    target_n <- as.integer(target_n)

    base_groups <- dplyr::group_vars(data)
    groups <- union(base_groups, id_col)

    df <- .enforce_time_order(data, id_col, time_col)

    # own gap check, within whatever grouping already exists (id, plus segment if present)
    df <- df |>
        dplyr::mutate(
            ..gap_id.. = .mark_runs(..t_s.., min_gap_s),
            .by = dplyr::all_of(groups)
        )

    atomic <- union(groups, "..gap_id..")
    n_runs_before <- nrow(dplyr::distinct(df[atomic]))

    # chop each atomic contiguous run into fixed-length blocks; drop the remainder
    df <- df |>
        dplyr::mutate(
            ..pos.. = seq_len(dplyr::n()),
            ..local_bout.. = ((..pos.. - 1L) %/% target_n) + 1L,
            sample_in_bout = ((..pos.. - 1L) %% target_n) + 1L,
            ..n_full.. = dplyr::n() %/% target_n,
            .by = dplyr::all_of(atomic)
        ) |>
        dplyr::filter(..local_bout.. <= ..n_full..)

    n_dropped <- nrow(data) - nrow(df)
    if (n_dropped > 0) {
        say(
            format(n_dropped, big.mark = ","),
            " row(s) dropped as leftover partial bout(s)."
        )
    }
    if (nrow(df) == 0) {
        rlang::abort(
            "No complete bouts of that length were found (every run was shorter than one bout)."
        )
    }

    bout_key <- union(atomic, "..local_bout..")

    df <- df |>
        dplyr::mutate(
            t_in_bout = ..t_s.. - ..t_s..[1],
            .by = dplyr::all_of(bout_key)
        ) |>
        dplyr::mutate(
            bout = vctrs::vec_group_id(dplyr::pick(dplyr::all_of(bout_key))),
            segment = vctrs::vec_group_id(dplyr::pick(dplyr::all_of(atomic)))
        ) |>
        dplyr::select(
            -dplyr::any_of(c(
                "..t_s..",
                "..gap_id..",
                "..pos..",
                "..local_bout..",
                "..n_full.."
            ))
        )

    n_bouts <- dplyr::n_distinct(df$bout)
    n_dropped_runs <- n_runs_before - dplyr::n_distinct(df$segment)
    if (n_dropped_runs > 0) {
        say(
            format(n_dropped_runs, big.mark = ","),
            " contiguous run(s) were entirely shorter than one bout and were dropped."
        )
    }
    say(
        "Result: ",
        format(n_bouts, big.mark = ","),
        " bout(s), ",
        format(nrow(df), big.mark = ","),
        " rows."
    )

    df <- df |>
        dplyr::relocate(
            dplyr::all_of(id_col),
            bout,
            sample_in_bout,
            t_in_bout,
            segment,
            dplyr::all_of(time_col)
        )

    out <- tsibble::as_tsibble(
        df,
        key = c(id_col, "bout"),
        index = "sample_in_bout"
    )
    if (any(tsibble::has_gaps(out)$.gaps)) {
        rlang::warn(
            "Unexpected gaps in the resulting tsibble -- please report this; the result may not be safe to use as-is."
        )
    }
    out
}

# plot_signal.R ---------------------------------------------------------------
# One entry point for interactive dygraph views of a tidy sensor tsibble
# (id, bout, <index>, acc_x, acc_y, acc_z, attachment, label, ...).
#
# Requires: tsibble, tibble, tidyselect, rlang, purrr, xts, dygraphs
#
# Three views, chosen with `mode`:
#   "time"        one individual, chosen streams over time
#   "behaviours"  one individual, data arranged by behaviour, one dyEvent per
#                 behaviour (unlabelled rows go in their own NA bucket, last)
#   "individuals" one behaviour, every individual with that label, one dyEvent
#                 per individual
#
# x-axis rules
#   * "time" mode + a real clock index (POSIXct/Date): the timestamp.
#   * Otherwise (no clock, or any mode that rearranges the data by label/id):
#     a stitched axis. Rows are ordered by the non-id key columns (e.g. bout)
#     and then the index (e.g. sample_in_bout), and laid end to end as
#     1, 2, 3, ... samples, or as seconds if `fs` is given. This is recorded
#     time only; gaps between bouts are not represented.
#
# Hover tooltip (modes "behaviours" and "individuals")
#   The legend's x-line shows the stitched position, the behaviour (or
#   individual) at that point, and where the point came from in the original
#   data (e.g. "bout 107, sample_in_bout 45", or the real timestamp).
#
# Colours: dyEvent lines and labels are black. Series use a colour-blind-safe
#   palette chosen to keep >= 3:1 contrast on a white background.

.series_cols <- c(
    "#0072B2",
    "#D55E00",
    "#009E73",
    "#AA3377",
    "#332288",
    "#999933",
    "#7A3B00"
)
.na_label <- "NA (unlabelled)"

#' Interactive dygraph views of a tidy sensor tsibble
#'
#' @param data      A tsibble (see the tidy ACC format).
#' @param cols      The data streams to draw, as tidyselect: `acc_x`,
#'                  `c(acc_x, acc_y, acc_z)`, `acc_x:acc_z`,
#'                  `starts_with("acc_")`, or any other numeric column.
#' @param mode      `"time"`, `"behaviours"` or `"individuals"` (see above).
#' @param id        Individual to plot. Required for `"time"` and
#'                  `"behaviours"` unless the data hold a single individual.
#'                  In `"individuals"` mode it is an optional subset of ids.
#' @param behaviour A single label value. Required for `"individuals"` mode.
#'                  Use `NA` to view unlabelled rows.
#' @param fs        Sampling rate in Hz. Only used for stitched axes, where it
#'                  turns the sample counter into seconds.
#' @param id_col,label_col Names of the id and label columns.
#' @param main      Plot title. A sensible default is built if `NULL`.
#'
#' @return A dygraphs htmlwidget.
#' @examples
#' plot_signal(vultures_tsbl, acc_x, id = "202367", fs = 20)
#' plot_signal(vultures_tsbl, c(acc_x, acc_y, acc_z),
#'             mode = "behaviours", id = "202367", fs = 20)
#' plot_signal(vultures_tsbl, acc_x:acc_z,
#'             mode = "individuals", behaviour = "Standing", fs = 20)
plot_signal <- function(
    data,
    cols,
    mode = c("time", "behaviours", "individuals"),
    id = NULL,
    behaviour = NULL,
    fs = NULL,
    id_col = "id",
    label_col = "label",
    main = NULL
) {
    mode <- match.arg(mode)

    if (!tsibble::is_tsibble(data)) {
        rlang::abort("`data` must be a tsibble in the tidy sensor format.")
    }
    if (missing(cols)) {
        rlang::abort(c(
            "Say which data streams to plot.",
            i = "e.g. `cols = acc_x` or `cols = c(acc_x, acc_y, acc_z)`."
        ))
    }

    idx_col <- tsibble::index_var(data)
    keys <- tsibble::key_vars(data)
    df <- tibble::as_tibble(data)

    # ---- which streams ---------------------------------------------------------
    sig <- names(tidyselect::eval_select(rlang::enquo(cols), df))
    if (length(sig) == 0) {
        rlang::abort("`cols` selected no columns.")
    }
    bad <- sig[!vapply(df[sig], is.numeric, logical(1))]
    if (length(bad)) {
        rlang::abort(paste0(
            "These `cols` are not numeric: ",
            paste(bad, collapse = ", ")
        ))
    }

    # ---- column checks ---------------------------------------------------------
    if (!id_col %in% names(df)) {
        rlang::abort(paste0(
            "No id column called `",
            id_col,
            "`. Set `id_col`."
        ))
    }
    if (mode != "time" && !label_col %in% names(df)) {
        rlang::abort(paste0(
            "Mode \"",
            mode,
            "\" needs a label column; `",
            label_col,
            "` not found. Set `label_col`."
        ))
    }

    has_clock <- inherits(df[[idx_col]], c("POSIXct", "Date"))
    ord <- c(setdiff(keys, id_col), idx_col) # within-individual ordering
    sort_rows <- function(d, ...) {
        d[do.call(order, c(list(...), unname(as.list(d[ord])))), , drop = FALSE]
    }

    blocks <- NULL # per-row block name, for dyEvent modes
    block_levels <- NULL # block order along the x-axis

    # ---- build the data for each mode -----------------------------------------
    if (mode == "time") {
        id <- .pick_one_id(df, id, id_col)
        d <- sort_rows(df[as.character(df[[id_col]]) == id, , drop = FALSE])
        default_main <- paste0("ID ", id)
    } else if (mode == "behaviours") {
        id <- .pick_one_id(df, id, id_col)
        d <- sort_rows(df[as.character(df[[id_col]]) == id, , drop = FALSE])
        lab <- as.character(d[[label_col]])

        known <- if (is.factor(d[[label_col]])) {
            levels(droplevels(d[[label_col]]))
        } else {
            sort(unique(lab[!is.na(lab)]))
        }
        block_levels <- c(known, if (anyNA(lab)) .na_label)
        blocks <- ifelse(is.na(lab), .na_label, lab)

        o <- order(match(blocks, block_levels)) # stable: keeps bout/sample order
        d <- d[o, , drop = FALSE]
        blocks <- blocks[o]
        default_main <- paste0(
            "ID ",
            id,
            " \u00b7 ",
            if (length(block_levels) == 1) {
                block_levels
            } else {
                "arranged by behaviour"
            }
        )
    } else {
        # "individuals"
        if (length(behaviour) != 1) {
            labs <- unique(as.character(df[[label_col]]))
            rlang::abort(c(
                "Supply a single `behaviour` for mode \"individuals\".",
                i = paste0(
                    "Labels in the data: ",
                    paste(ifelse(is.na(labs), "NA", labs), collapse = ", ")
                )
            ))
        }
        lab <- df[[label_col]]
        keep <- if (is.na(behaviour)) {
            is.na(lab)
        } else {
            !is.na(lab) & as.character(lab) == as.character(behaviour)
        }
        d <- df[keep, , drop = FALSE]
        if (!is.null(id)) {
            d <- d[
                as.character(d[[id_col]]) %in% as.character(id),
                ,
                drop = FALSE
            ]
        }
        if (nrow(d) == 0) {
            rlang::abort("No rows match that `behaviour` (and `id`, if given).")
        }

        present <- unique(as.character(d[[id_col]]))
        block_levels <- if (is.null(id)) {
            sort(present)
        } else {
            intersect(as.character(id), present)
        }
        d <- sort_rows(d, match(as.character(d[[id_col]]), block_levels))
        blocks <- as.character(d[[id_col]])
        default_main <- paste0(
            "Behaviour: ",
            if (is.na(behaviour)) "unlabelled" else behaviour,
            " \u00b7 ",
            if (length(block_levels) == 1) {
                paste0("ID ", block_levels)
            } else {
                "all individuals"
            }
        )
    }

    if (nrow(d) == 0) {
        rlang::abort("No rows to plot for that selection.")
    }
    if (nrow(d) > 2e5) {
        rlang::inform(paste0(
            format(nrow(d), big.mark = ","),
            " points; the dygraph may be slow to draw."
        ))
    }

    # ---- x-axis ----------------------------------------------------------------
    use_clock <- mode == "time" && has_clock
    if (use_clock) {
        x <- d[[idx_col]]
        xlab <- "Time"
    } else {
        n <- nrow(d)
        x <- if (is.null(fs)) seq_len(n) else (seq_len(n) - 1) / fs
        xlab <- if (is.null(fs)) "Stitched sample" else "Stitched time (s)"
    }

    # ---- dygraph ---------------------------------------------------------------
    mat <- as.matrix(d[sig])
    main <- main %||% default_main

    g <- if (use_clock) {
        dygraphs::dygraph(xts::xts(mat, order.by = x), main = main, xlab = xlab)
    } else {
        dygraphs::dygraph(
            data.frame(x = x, mat, check.names = FALSE),
            main = main,
            xlab = xlab
        )
    }

    g <- g |>
        dygraphs::dyOptions(
            colors = rep_len(.series_cols, length(sig)),
            useDataTimezone = use_clock
        ) |>
        dygraphs::dyRangeSelector()

    # ---- dyEvent per block (behaviour or individual) ---------------------------
    if (!is.null(blocks)) {
        pos <- x[match(block_levels, blocks)]

        g <- purrr::reduce(
            seq_along(block_levels),
            function(g, i) {
                dygraphs::dyEvent(
                    g,
                    x = pos[i],
                    label = block_levels[i],
                    labelLoc = "top",
                    color = "black",
                    strokePattern = "dashed"
                )
            },
            .init = g
        )

        # Hover text: position, block, and provenance in the original data.
        # rangePad keeps an event on the very first point off the axis edge.
        tip_levels <- if (mode == "individuals") {
            paste0("ID ", block_levels)
        } else {
            block_levels
        }
        g <- dygraphs::dyAxis(
            g,
            "x",
            rangePad = 10,
            valueFormatter = .block_tooltip_js(
                d = d,
                ord = ord,
                blocks = blocks,
                block_levels = block_levels,
                tip_levels = tip_levels,
                fs = fs
            )
        )
    }

    g
}

# ---- helpers -----------------------------------------------------------------

# JS x-axis valueFormatter for the legend. The per-row lookup tables are built
# once (closure) and indexed by row, so hovering stays cheap on large plots.
.block_tooltip_js <- function(d, ord, blocks, block_levels, tip_levels, fs) {
    extra <- lapply(d[ord], function(v) {
        if (inherits(v, "POSIXct")) {
            format(v, "%Y-%m-%d %H:%M:%OS2")
        } else if (inherits(v, "Date")) {
            format(v)
        } else {
            v
        }
    })
    meta <- list(
        levels = tip_levels,
        block = match(blocks, block_levels) - 1L,
        extra = extra,
        pre = jsonlite::unbox(if (is.null(fs)) "sample " else ""),
        suf = jsonlite::unbox(if (is.null(fs)) "" else " s"),
        dp = jsonlite::unbox(if (is.null(fs)) 0L else 2L)
    )
    meta_json <- jsonlite::toJSON(meta, digits = NA, na = "string")

    htmlwidgets::JS(paste0(
        "(function() {",
        "  var m = ",
        meta_json,
        ";",
        "  var keys = Object.keys(m.extra);",
        "  return function(x, opts, seriesName, g, row, col) {",
        "    var xs = m.dp === 0 ? String(Math.round(x)) : Number(x).toFixed(m.dp);",
        "    var where = keys.map(function(k) { return k + ' ' + m.extra[k][row]; }).join(', ');",
        "    return m.pre + xs + m.suf + ' \\u00b7 ' + m.levels[m.block[row]] +",
        "           (where ? ' \\u00b7 ' + where : '');",
        "  };",
        "})()"
    ))
}

.pick_one_id <- function(df, id, id_col) {
    ids <- unique(as.character(df[[id_col]]))
    shown <- paste0(
        paste(utils::head(ids, 10), collapse = ", "),
        if (length(ids) > 10) ", ..." else ""
    )
    if (is.null(id)) {
        if (length(ids) == 1) {
            return(ids)
        }
        rlang::abort(c("Please supply `id`.", i = paste0("Available: ", shown)))
    }
    id <- as.character(id)
    if (length(id) != 1 || !id %in% ids) {
        rlang::abort(c(
            "`id` must be a single id present in the data.",
            i = paste0("Available: ", shown)
        ))
    }
    id
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# calculate_features.R --------------------------------------------------------
# Bout-level feature calculation for a tidy sensor tsibble
# (key = c(id, bout), index = sample_in_bout or datetime, one column per stream).
# Modelled on theft::calculate_features(), but multivariate: any sensor columns
# can be chosen, and features that combine axes (ODBA, VeDBA, pitch, cor_x_y)
# are supported.
#
# Requires: tsibble, tibble, dplyr, tidyr, tidyselect, rlang, vctrs, zoo,
#           entropy, Rcatch22
#
# Named calc_features(), not calculate_features(), so it can be loaded
# alongside {theft} without masking theft::calculate_features().
#
# =============================================================================
# HOW FEATURES ARE ORGANISED: a library, with sets as membership lists
# =============================================================================
# Every distinct feature is computed ONCE, under one canonical column name:
#
#     <family>__<stream>__<feature>
#
#   family   how it is computed: stat, seq, pair, dba, orient, freq, catch22, custom
#   stream   the column(s) it is computed from, joined with "_" for features
#            that use several (e.g. acc_x_acc_y_acc_z); always DEVICE names
#   feature  the feature itself (mean, sd, ODBA, pitch_mean, ...)
#
#   e.g. stat__acc_z__sd, pair__acc_x_acc_y__cor, dba__acc_x_acc_y_acc_z__ODBA,
#        orient__acc_x_acc_y_acc_z__pitch_mean, catch22__acc_x__DN_HistogramMode_5
#
# A published feature set is a LIST of library features. Requesting several
# sets returns their union with no duplicated columns: collins and
# harvey_carroll share 9 features, rabc_time2 contains all of rabc_time1, and
# so on. Which set(s) each column belongs to is recorded in a lookup table,
# returned as the "feature_info" attribute of the output (one row per column x
# set) and read with feature_info(x). project_features() and
# explore_features() use it, so `feature_set = "collins"` still works there.
# feature_sets() lists every set and what it contains.
#
# =============================================================================
# SENSORS AND AXIS ROLES (`sensors` argument)
# =============================================================================
# A small table, one row per sensor column, passed explicitly on each call:
#
#   sensor  type           column  role   sign
#   acc     accelerometer  acc_x   sway      1
#   acc     accelerometer  acc_y   surge     1
#   acc     accelerometer  acc_z   heave     1
#
#   sensor  groups columns into one physical sensor (free text)
#   type    "accelerometer" is recognised by the accelerometer sets (DBA, pitch,
#           roll ...). Anything else ("gyroscope", "temperature") is allowed and
#           is used only by the generic sets (catch22, stats) and custom features.
#   role    optional, free text: name the body axis however suits the sensor
#           (surge/sway/heave for accelerometers, roll/pitch/yaw for gyros...).
#           The published accelerometer sets that are defined on body axes look
#           for roles named exactly "surge", "sway" and "heave".
#   sign    optional, 1 (default) or -1. -1 means the device axis points the
#           opposite way to the role's positive direction (surge: towards the
#           head; heave: towards the back). Signs are applied to that column
#           BEFORE any feature is computed, so every feature of that column uses
#           the corrected values; a message lists flipped columns.
#
# `sensors = NULL` (default): columns named acc_x / acc_y / acc_z among `cols`
#   are treated as one accelerometer with no roles -- the previous behaviour.
#
# One mapping per call. Roles that change by individual, attachment or over
# time (tag slippage) are not handled here: calc_features() never re-orients
# data. That is the job of separate, explicit steps on the per-sample tsibble
# (planned: check_orientation() / reorient()), after which one mapping holds.
#
# =============================================================================
# FEATURE SETS (see feature_sets() for the exact contents)
# =============================================================================
# Generic, any sensor:
#   catch22        Rcatch22 on every selected stream; catch24 = TRUE adds
#                  DN_Mean, DN_Spread_Std.
#   stats          every stat-family feature (below) on every selected stream.
#   eikelboom      Eikelboom et al. 2026 IMU features, on every selected stream
#                  (acc, gyro, magnetometer ...): summary -- min, q25, q50, q75,
#                  max, mean, sd, skew, kurtosis; sequential -- acf1, diff_acf1,
#                  diff_sd, mean_crossings, arv; frequency -- fft_amp1..5,
#                  fft_freq1..5, fft_energy, fft_entropy. 26 per stream (78 on a
#                  tri-axial accelerometer). Needs `fs`.
#
# Accelerometer, device axes (no roles needed):
#   rabc_time1     rabc time domain: per-axis mean, variance, sd, max, min,
#                  range; ODBA. Exact parity with rabc.
#   rabc_time2     rabc_time1 + per-axis norm; per-pair cov, cor, meandiff,
#                  sddiff; per-axis varsba, vardba, maxdba. Exact parity.
#   rabc_freq      per-axis freqmain, freqamp, entropy. Exact parity.
#   yu_full        Yu et al. 2021 Table 1, 78 features on a tri-axial tag:
#                  per-axis mean, variance, sd, cv, skew, kurtosis, max, min,
#                  range, norm, absdiff_mean, diff_var, q25, q50, q75; per-pair
#                  cov, cor, meandiff, sddiff; per-axis varsba, vardba, meandba,
#                  maxdba; ODBA; pitch_mean, roll_mean; per-axis freqmain,
#                  freqamp. Pitch and roll need roles: without them they are
#                  skipped (with a message), leaving 76.
#   arrieta        Arrieta 2025 Table 3.1: per-axis kurtosis, skew, static
#                  acceleration (= mean), dynamic acceleration (= meandba);
#                  ODBA; VeDBA. 14 on a tri-axial tag.
#
# Accelerometer, body axes (roles surge/sway/heave REQUIRED; error if missing):
#   yu_reduced     Yu et al. 2021 simplified set (backpacks, collars, legs):
#                  mean and sd of surge and heave; ODBA. 5 features.
#   yu_reduced_ear Yu et al. 2021 ear-tag variant: mean and sd of heave, heave
#                  freqmain; ODBA. 4 features.
#   collins        Collins et al. 2015: mean and sd of surge, sway, heave;
#                  pitch_mean, roll_mean; ODBA; VeDBA. 10 features.
#   harvey_carroll Harvey-Carroll et al. 2024: mean, sd, min, max of surge,
#                  sway, heave; pitch_mean, pitch_sd, roll_mean, roll_sd; ODBA;
#                  VeDBA. 18 features. Strict superset of collins.
#
# With fewer than three accelerometer axes selected, nothing errors: ODBA /
# VeDBA become PDBA_<axes> / PVeDBA_<axes> (e.g. PDBA_xz) and pair features are
# limited to the axes present. Sets needing a role that is not selected error.
#
# =============================================================================
# FEATURE DEFINITIONS (harmonised: one definition per feature, per bout)
# =============================================================================
#   stat   mean, variance (n-1), sd (n-1), cv = sd / mean, skew = m3 / m2^1.5,
#          kurtosis = m4 / m2^2 - 3 (excess), max, min, range, norm =
#          sqrt(sum(a^2)), absdiff_mean = mean(|diff(a)|), diff_var =
#          var(diff(a)), q25 / q50 / q75 (quantile type 7).
#   seq    acf1 = lag-1 autocorrelation (as stats::acf / tsfeatures x_acf1),
#          diff_acf1 = acf1 of diff(a) (tsfeatures diff1_acf1), diff_sd =
#          sd(diff(a)), mean_crossings = number of times a crosses its own
#          mean (sign changes of a - mean(a)), arv = average rectified value
#          mean(|a|).
#   pair   cov, cor (Pearson), meandiff = mean(a - b), sddiff = sd(a - b),
#          in the order xy, yz, xz.
#   dba    static acceleration s = centred rolling mean over `winlen_dba`;
#          dynamic d = a - s. Per axis: varsba = var(s), vardba = var(|d|),
#          meandba = mean(|d|), maxdba = max(|d|). Vector: ODBA = mean of
#          sum(|d|) over axes; VeDBA = mean of sqrt(sum(d^2)). Edge samples
#          (NA static) are excluded.
#   orient per sample, in degrees: pitch = atan(surge / sqrt(sway^2 + heave^2)),
#          roll = atan(sway / sqrt(surge^2 + heave^2)) (Collins 2015 eqs);
#          then pitch_mean, pitch_sd, roll_mean, roll_sd over the bout.
#   freq   detrended FFT: freqmain = bin * fs / n (rabc convention: one bin
#          high, i.e. 0.2 Hz high for 100 samples at 20 Hz), freqamp = amplitude
#          of that bin, entropy = spectral entropy^2 / (n / 2).
#          fft_* (eikelboom): FFT of a - mean(a), one-sided, DC excluded.
#          fft_amp1..5 = the five largest single-sided amplitudes (2|X_k|/n,
#          in the signal's units; |X_k|/n at Nyquist), descending; fft_freq1..5
#          = their frequencies k * fs / n (true bin frequency, NOT the rabc
#          convention); fft_energy = sum over non-DC bins of |X_k|^2 / n^2
#          (Bao & Intille 2004); fft_entropy = Shannon entropy of the normalised
#          one-sided power spectrum, divided by log(number of bins), so 0-1.
#
# =============================================================================
# HOW EACH PAPER DIFFERS FROM THESE DEFINITIONS, AND HOW TO MATCH IT
# =============================================================================
# Yu et al. 2021 (yu_full, yu_reduced, yu_reduced_ear)
#   * Computed per bout on device axes, as here.
#   * "Mean difference of continuous points" is ambiguous: read literally,
#     mean(diff(a)) = (last - first) / (n - 1), which carries almost no
#     information. We use mean(|diff(a)|). For the literal version:
#     features = list(diff_mean = \(v) mean(diff(v))).
#   * Pitch / roll: Yu derives them from the bout means of the axes (Table 2
#     lists "mean" as their dependency); we average per-sample angles
#     (Collins). Nearly identical unless the bout is very dynamic. To match:
#     cross_feature(c("acc_y", "acc_x", "acc_z"),
#       \(s, w, h) atan(mean(s) / sqrt(mean(w)^2 + mean(h)^2)) * 180 / pi)
#   * Skew / kurtosis estimator unspecified (R default via e1071 is type 3).
#     For another type: features = list(kurt3 = \(v) e1071::kurtosis(v, type = 3)).
#   * DBA smoothing window is not reported; set winlen_dba(_s) to the value used.
#
# Collins et al. 2015 (collins)
#   * Metrics computed over 1 s MOVING windows, not per bout. Here each bout is
#     summarised once; for 1 s bouts the two coincide: make_bouts(duration_s = 1).
#   * Static acceleration smoothed over 1 s: winlen_dba_s = 1.
#
# Harvey-Carroll et al. 2024 (harvey_carroll)
#   * Per window (5 s at 50 Hz best for discrete behaviours, 7 s at 10 Hz for
#     activity levels): make_bouts(duration_s = 5).
#   * Static-acceleration window not stated in the main text (Additional file 2);
#     set winlen_dba(_s) accordingly.
#
# Arrieta 2025 (arrieta)
#   * Features on 5 s RIGHT-aligned rolling windows at every sample, fed in
#     alongside the raw signal after resampling to 10 or 25 Hz. Here: per bout.
#   * Static acceleration = mean over the whole 5 s window, and dynamic = raw
#     minus that mean. Our meandba uses a centred smoothing window instead; to
#     match Arrieta, use 5 s bouts and winlen_dba = the bout length in samples
#     (static then = bout mean, with edge NAs; or write it directly:
#     features = list(dyn = \(v) mean(abs(v - mean(v))))).
#   * Table 3.1 formula typos (ODBA's third term |az - ay|, VeDBA without mean
#     subtraction, malformed skew denominator): standard definitions used.
#   * Kurtosis/skew use the bias-adjusted (k2-based) estimators; we use the
#     moment estimators above (differences are negligible at n >= 50).
#
# Eikelboom et al. 2026 (eikelboom)
#   * The paper lists the features but not their formulas; definitions follow
#     the sources it cites: Hyndman et al. 2020 (tsfeatures) for acf1 /
#     diff_acf1, Bao & Intille 2004 for spectral energy and entropy.
#   * Computed per input variable: 3 per feature for tri-axial data, 1 for the
#     "dimension-reduced" (magnitude) input, 21 for the "blown-up" input. To
#     reproduce the magnitude version, add a column first, e.g.
#     mutate(acc_mag = sqrt(acc_x^2 + acc_y^2 + acc_z^2)), and select it.
#   * fft_energy excluding DC equals the population variance of the bout
#     exactly (Parseval), so it is redundant with variance / sd. If their
#     implementation kept the DC term, the equivalent is mean(a^2):
#     features = list(energy_dc = \(v) mean(v^2)).
#   * "Average rectified value" read as mean(|a|), the standard definition;
#     on raw acceleration this mostly reflects posture (gravity). For the
#     dynamic version: features = list(arv_dyn = \(v) mean(abs(v - mean(v)))).
#   * Entropy normalised to 0-1 (divided by log of the number of bins) so it
#     is comparable across bout lengths; un-normalised = fft_entropy * log(K).
#
# rabc (rabc_time1, rabc_time2, rabc_freq): exact parity with the original
#   functions, checked numerically (rolling mean, ODBA, freq conventions).
#
# =============================================================================
# OUTPUT, MISSING VALUES, CUSTOM FEATURES
# =============================================================================
# A plain tibble, one row per bout: key columns (e.g. id, bout), then context
# columns (any non-selected column constant within every bout, e.g.
# attachment, label), then feature columns, in the order of the requested sets
# (each feature at its first appearance), then custom features. Attribute
# "feature_info" holds the set membership; read it with feature_info(x).
# feature_long() pivots to a long theft-style layout.
#
# Missing values: a feature is NA for a bout if any column it uses has an NA /
# non-finite value in that bout (per stream for single-stream features; any of
# the axes for ODBA, VeDBA, pitch, roll, pair features), or if it is undefined
# (cor of a constant signal, cv with mean 0). One message per family reports
# how many bouts are affected.
#
# Custom features (`features`, as in theft::calculate_features())
#   * a plain function -> applied to every selected stream separately:
#     list(mean = mean, sd = sd) -> custom__<stream>__mean, ...
#   * cross_feature(cols, fn, window = FALSE) -> one feature from several
#     streams; arguments are POSITIONAL, matched to `cols` in the order given.
#     window = TRUE appends the resolved smoothing window (samples) as one more
#     argument. roll_static(v, winlen) is the exact static-acceleration filter
#     the dba family uses. Column: custom__<col1>_<col2>..__<name>.
#   Every function must return one number; one that errors or returns
#   something else yields NA for that bout.
#
# Dynamic-acceleration window (any dba-family feature, or cross_feature(window = TRUE)):
#   winlen_dba_s (seconds, needs fs; rounded to an odd number of samples:
#   2 * floor(winlen_dba_s * fs / 2) + 1) OR winlen_dba (samples). No default.
#
# `seed`: every feature implemented so far is deterministic. The argument
#   exists for parity with theft::calculate_features() and future sets.

# =============================================================================
# Feature registry
# =============================================================================

.stat_features <- c(
    "mean",
    "variance",
    "sd",
    "cv",
    "skew",
    "kurtosis",
    "max",
    "min",
    "range",
    "norm",
    "absdiff_mean",
    "diff_var",
    "q25",
    "q50",
    "q75"
)
.pair_features <- c("cov", "cor", "meandiff", "sddiff")
.seq_features <- c("acf1", "diff_acf1", "diff_sd", "mean_crossings", "arv")
.fft_features <- c(
    paste0("fft_amp", 1:5),
    paste0("fft_freq", 1:5),
    "fft_energy",
    "fft_entropy"
)
.body_roles <- c("surge", "sway", "heave")

.reg <- function(set, family, feature, target) {
    tibble::tibble(
        feature_set = set,
        family = family,
        feature = feature,
        target = target
    )
}

# One row per (set, feature, target), in each set's documented order.
# target: "stream"  every selected stream
#         "axis"    every selected axis of each accelerometer
#         "pair"    every pair of selected axes of each accelerometer
#         "vector"  each accelerometer as a whole (ODBA, VeDBA)
#         "role:<r>" the axis with role <r>
#         "roles"   the surge, sway and heave axes together (pitch, roll)
.feature_registry <- function() {
    rabc1 <- list(
        .reg(
            "rabc_time1",
            "stat",
            c("mean", "variance", "sd", "max", "min", "range"),
            "axis"
        ),
        .reg("rabc_time1", "dba", "ODBA", "vector")
    )
    rabc2 <- c(
        lapply(rabc1, function(r) {
            r$feature_set <- "rabc_time2"
            r
        }),
        list(
            .reg("rabc_time2", "stat", "norm", "axis"),
            .reg("rabc_time2", "pair", .pair_features, "pair"),
            .reg("rabc_time2", "dba", c("varsba", "vardba", "maxdba"), "axis")
        )
    )
    role_stats <- function(set, roles, feats) {
        lapply(feats, function(f) {
            dplyr::bind_rows(lapply(roles, function(r) {
                .reg(set, "stat", f, paste0("role:", r))
            }))
        })
    }

    dplyr::bind_rows(
        .reg("catch22", "catch22", "catch22", "stream"),
        .reg("stats", "stat", .stat_features, "stream"),
        rabc1,
        rabc2,
        .reg("rabc_freq", "freq", c("freqmain", "freqamp", "entropy"), "axis"),
        # Yu et al. 2021, Table 1 order
        .reg(
            "yu_full",
            "stat",
            c(
                "mean",
                "variance",
                "sd",
                "cv",
                "skew",
                "kurtosis",
                "max",
                "min",
                "range",
                "norm"
            ),
            "axis"
        ),
        .reg("yu_full", "pair", .pair_features, "pair"),
        .reg(
            "yu_full",
            "dba",
            c("varsba", "vardba", "meandba", "maxdba"),
            "axis"
        ),
        .reg("yu_full", "dba", "ODBA", "vector"),
        .reg("yu_full", "orient", c("pitch_mean", "roll_mean"), "roles"),
        .reg("yu_full", "stat", c("absdiff_mean", "diff_var"), "axis"),
        .reg("yu_full", "freq", c("freqmain", "freqamp"), "axis"),
        .reg("yu_full", "stat", c("q25", "q50", "q75"), "axis"),
        # Yu et al. 2021 simplified sets
        role_stats("yu_reduced", c("surge", "heave"), c("mean", "sd")),
        .reg("yu_reduced", "dba", "ODBA", "vector"),
        role_stats("yu_reduced_ear", "heave", c("mean", "sd")),
        .reg("yu_reduced_ear", "dba", "ODBA", "vector"),
        .reg("yu_reduced_ear", "freq", "freqmain", "role:heave"),
        # Collins et al. 2015
        role_stats("collins", .body_roles, c("mean", "sd")),
        .reg("collins", "orient", c("pitch_mean", "roll_mean"), "roles"),
        .reg("collins", "dba", c("ODBA", "VeDBA"), "vector"),
        # Harvey-Carroll et al. 2024
        .reg("harvey_carroll", "dba", c("ODBA", "VeDBA"), "vector"),
        role_stats(
            "harvey_carroll",
            .body_roles,
            c("mean", "sd", "min", "max")
        ),
        .reg(
            "harvey_carroll",
            "orient",
            c("pitch_mean", "pitch_sd", "roll_mean", "roll_sd"),
            "roles"
        ),
        # Arrieta 2025, Table 3.1
        .reg("arrieta", "stat", c("kurtosis", "skew", "mean"), "axis"),
        .reg("arrieta", "dba", "meandba", "axis"),
        .reg("arrieta", "dba", c("ODBA", "VeDBA"), "vector"),
        # Eikelboom et al. 2026: summary, sequential, frequency-domain, per stream
        .reg(
            "eikelboom",
            "stat",
            c(
                "min",
                "q25",
                "q50",
                "q75",
                "max",
                "mean",
                "sd",
                "skew",
                "kurtosis"
            ),
            "stream"
        ),
        .reg("eikelboom", "seq", .seq_features, "stream"),
        .reg("eikelboom", "freq", .fft_features, "stream")
    )
}

.set_info <- tibble::tribble(
    ~feature_set     , ~needs                                   , ~source                                            ,
    "catch22"        , "any sensor"                             , "Lubba et al. 2019 (Rcatch22)"                     ,
    "stats"          , "any sensor"                             , "generic per-stream summaries"                     ,
    "rabc_time1"     , "accelerometer"                          , "rabc (Yu & Klaassen); exact parity"               ,
    "rabc_time2"     , "accelerometer"                          , "rabc expanded; exact parity"                      ,
    "rabc_freq"      , "accelerometer"                          , "rabc; exact parity"                               ,
    "yu_full"        , "accelerometer (+ roles for pitch/roll)" , "Yu et al. 2021, Table 1"                          ,
    "yu_reduced"     , "accelerometer + roles"                  , "Yu et al. 2021, simplified (backpack/collar/leg)" ,
    "yu_reduced_ear" , "accelerometer + roles"                  , "Yu et al. 2021, simplified (ear tag)"             ,
    "collins"        , "accelerometer + roles"                  , "Collins et al. 2015"                              ,
    "harvey_carroll" , "accelerometer + roles"                  , "Harvey-Carroll et al. 2024"                       ,
    "arrieta"        , "accelerometer"                          , "Arrieta 2025, Table 3.1"                          ,
    "eikelboom"      , "any sensor (needs fs)"                  , "Eikelboom et al. 2026, IMU features"
)
.feature_set_names <- .set_info$feature_set
.role_sets <- c("yu_reduced", "yu_reduced_ear", "collins", "harvey_carroll")

#' List the available feature sets and what each contains
#'
#' @param set Optional character vector of set names to show.
#' @return A tibble, one row per (set, family, feature, applies_to), plus where
#'   the set comes from and what it needs.
#' @examples
#' feature_sets()
#' feature_sets("collins")
#' feature_sets() |> dplyr::count(feature_set)
feature_sets <- function(set = NULL) {
    applies <- c(
        stream = "every selected stream",
        axis = "each accelerometer axis",
        pair = "each pair of accelerometer axes",
        vector = "accelerometer as a whole",
        roles = "surge + sway + heave together"
    )
    reg <- .feature_registry()
    reg$applies_to <- ifelse(
        startsWith(reg$target, "role:"),
        paste0("the ", sub("role:", "", reg$target, fixed = TRUE), " axis"),
        applies[reg$target]
    )
    reg$feature[
        reg$family == "catch22"
    ] <- "22 catch22 features (24 with catch24 = TRUE)"
    out <- dplyr::left_join(
        reg[c("feature_set", "family", "feature", "applies_to")],
        .set_info,
        by = "feature_set"
    )
    if (!is.null(set)) {
        unknown <- setdiff(set, .feature_set_names)
        if (length(unknown)) {
            rlang::abort(paste0(
                "Unknown set(s): ",
                paste(unknown, collapse = ", ")
            ))
        }
        out <- out[out$feature_set %in% set, ]
    }
    out
}

#' Feature-set membership of a calc_features() output
#'
#' @param x Output of calc_features().
#' @return A tibble, one row per (column, feature_set): column, feature_set,
#'   family, stream, feature. NULL (with a message) if the attribute was lost,
#'   e.g. after bind_rows() or a round trip through a file.
feature_info <- function(x) {
    info <- attr(x, "feature_info", exact = TRUE)
    if (is.null(info)) {
        rlang::inform(
            "No feature_info attribute (lost by e.g. bind_rows(), or saved without it)."
        )
        return(NULL)
    }
    info[info$column %in% names(x), ]
}

# =============================================================================
# Main function
# =============================================================================

#' Calculate bout-level features from a sensor tsibble
#'
#' @param data        A tsibble. Every key combination (e.g. id + bout) is one
#'                    bout and yields one row of features.
#' @param cols        Sensor streams to use, as tidyselect: `acc_x`,
#'                    `c(acc_x, acc_y, acc_z)`, `starts_with("acc_")`, ...
#' @param feature_set One or more set names (see feature_sets()). Can be
#'                    `NULL` if `features` is given.
#' @param catch24     Also compute catch24 when catch22 is requested.
#' @param features    Optional named list of custom features: plain functions
#'                    (applied per stream) and/or cross_feature() entries.
#' @param sensors     Optional sensor table (sensor, type, column, role, sign);
#'                    see the file header. NULL treats acc_x/acc_y/acc_z as one
#'                    accelerometer without roles.
#' @param winlen_dba_s Static-acceleration window in SECONDS (needs `fs`).
#' @param winlen_dba  The same window in SAMPLES. Give one of the two whenever a
#'                    dba-family feature or cross_feature(window = TRUE) is used.
#' @param fs          Sampling rate in Hz. Needed for freq features and for
#'                    `winlen_dba_s`.
#' @param seed        Integer for set.seed(), or NULL to leave the RNG alone.
#' @param verbose     Print progress messages.
#'
#' @return A tibble with one row per bout and a "feature_info" attribute
#'   (see the file header and feature_info()).
#' @examples
#' vulture_sensors <- tibble::tribble(
#'   ~sensor, ~type,           ~column, ~role,   ~sign,
#'   "acc",   "accelerometer", "acc_x", "sway",  1,
#'   "acc",   "accelerometer", "acc_y", "surge", 1,
#'   "acc",   "accelerometer", "acc_z", "heave", 1
#' )
#' calc_features(vultures_tsbl, c(acc_x, acc_y, acc_z),
#'               feature_set = c("collins", "harvey_carroll", "yu_full"),
#'               sensors = vulture_sensors, winlen_dba_s = 1, fs = 20)
calc_features <- function(
    data,
    cols,
    feature_set = "catch22",
    catch24 = FALSE,
    features = NULL,
    sensors = NULL,
    winlen_dba_s = NULL,
    winlen_dba = NULL,
    fs = NULL,
    seed = 123,
    verbose = TRUE
) {
    say <- function(...) if (verbose) rlang::inform(paste0(...))

    # ---- input checks ----------------------------------------------------------
    if (!tsibble::is_tsibble(data)) {
        rlang::abort("`data` must be a tsibble in the tidy sensor format.")
    }
    if (missing(cols)) {
        rlang::abort(c(
            "Say which sensor columns to use.",
            i = "e.g. `cols = c(acc_x, acc_y, acc_z)`."
        ))
    }

    feature_set <- unique(tolower(feature_set))
    unknown <- setdiff(feature_set, .feature_set_names)
    if (length(unknown)) {
        rlang::abort(c(
            paste0("Unknown `feature_set`: ", paste(unknown, collapse = ", ")),
            i = paste0(
                "Available: ",
                paste(.feature_set_names, collapse = ", ")
            ),
            i = "See feature_sets() for what each contains."
        ))
    }
    features <- .validate_custom_features(features)
    uses_custom <- !is.null(features)
    uses_custom_window <- uses_custom &&
        any(vapply(
            features,
            function(e) inherits(e, "cross_feature") && isTRUE(e$window),
            logical(1)
        ))
    if (length(feature_set) == 0 && !uses_custom) {
        rlang::abort(c(
            "Nothing to calculate.",
            i = "Set `feature_set` (see feature_sets()) and/or `features`, a named list of functions."
        ))
    }

    keys <- tsibble::key_vars(data)
    idx <- tsibble::index_var(data)
    if (length(keys) == 0) {
        rlang::abort(c(
            "`data` needs key column(s) that identify each bout.",
            i = "e.g. `as_tsibble(x, key = c(id, bout), index = sample_in_bout)`."
        ))
    }
    df <- tibble::as_tibble(data)

    sig <- names(tidyselect::eval_select(rlang::enquo(cols), df))
    if (length(sig) == 0) {
        rlang::abort("`cols` selected no columns.")
    }
    if (any(sig %in% c(keys, idx))) {
        rlang::abort("`cols` must be sensor columns, not key or index columns.")
    }
    bad_type <- sig[!vapply(df[sig], is.numeric, logical(1))]
    if (length(bad_type)) {
        rlang::abort(paste0(
            "These `cols` are not numeric: ",
            paste(bad_type, collapse = ", ")
        ))
    }

    if (uses_custom) {
        is_cross <- vapply(
            features,
            inherits,
            logical(1),
            what = "cross_feature"
        )
        for (nm in names(features)[is_cross]) {
            missing_cols <- setdiff(features[[nm]]$cols, sig)
            if (length(missing_cols)) {
                rlang::abort(c(
                    paste0(
                        "features$",
                        nm,
                        " (cross_feature) needs column(s) not in `cols`: ",
                        paste(missing_cols, collapse = ", ")
                    ),
                    i = "Add them to `cols`, or remove this feature."
                ))
            }
        }
    }
    if (isTRUE(catch24) && !"catch22" %in% feature_set) {
        rlang::warn(
            "`catch24 = TRUE` is ignored because \"catch22\" is not in `feature_set`."
        )
    }

    # ---- sensors, roles and the concrete features to compute ---------------------
    sens <- .resolve_sensors(sensors, sig, names(df), feature_set)
    plan <- .plan_features(feature_set, sig, sens$acc)

    uses_window <- any(plan$specs$family == "dba") || uses_custom_window
    uses_freq <- any(plan$specs$family == "freq")

    fs_ok <- !is.null(fs) &&
        is.numeric(fs) &&
        length(fs) == 1 &&
        !is.na(fs) &&
        fs > 0
    if (!is.null(fs) && !fs_ok) {
        rlang::abort(
            "`fs` must be a single positive number (sampling rate in Hz)."
        )
    }
    if (uses_freq && !fs_ok) {
        rlang::abort(
            "Frequency features (rabc_freq, yu_full, yu_reduced_ear, eikelboom) need `fs`, the sampling rate in Hz."
        )
    }
    win_requested <- winlen_dba_s
    if (uses_window) {
        winlen_dba <- .resolve_winlen(winlen_dba, winlen_dba_s, fs, fs_ok)
    } else if (!is.null(winlen_dba) || !is.null(winlen_dba_s)) {
        say(
            "`winlen_dba` / `winlen_dba_s` are ignored: not used by any requested feature."
        )
    }

    # messages (only after every argument check has passed)
    if (sens$defaulted && plan$uses_acc) {
        say(
            "No `sensors` table: treating ",
            paste(sens$acc[[1]]$cols, collapse = ", "),
            " as one accelerometer without roles."
        )
    }
    if (length(sens$flipped)) {
        say(
            "Sign-flipped before computing: ",
            paste(sens$flipped, collapse = ", ")
        )
    }
    if (length(sens$unused)) {
        say(
            "In `sensors` but not selected in `cols` (ignored): ",
            paste(sens$unused, collapse = ", ")
        )
    }
    for (m in plan$messages) {
        say(m)
    }
    if (!is.null(seed)) {
        set.seed(seed)
    }

    # ---- bout structure ----------------------------------------------------------
    df <- dplyr::arrange(df, dplyr::across(dplyr::all_of(c(keys, idx))))
    gid <- vctrs::vec_group_id(df[keys]) # consecutive after the sort above
    n_bouts <- max(gid)
    first_idx <- match(seq_len(n_bouts), gid)

    streams <- lapply(stats::setNames(sig, sig), function(col) {
        unname(split(df[[col]] * sens$sign[[col]], gid))
    })
    bad <- lapply(streams, function(s) {
        vapply(s, function(v) any(!is.finite(v)), logical(1))
    })

    if (uses_window) {
        n_per_bout <- lengths(streams[[1]])
        say(.window_note(
            winlen_dba,
            win_requested,
            if (fs_ok) fs,
            stats::median(n_per_bout)
        ))
        shortest <- min(n_per_bout)
        if (shortest < winlen_dba) {
            rlang::warn(paste0(
                "Some bouts have fewer samples (",
                shortest,
                ") than the window (",
                winlen_dba,
                "); for those bouts the dba-family features",
                if (uses_custom_window) {
                    " and window-using custom features"
                } else {
                    ""
                },
                " will be NA."
            ))
        }
    }

    # ---- key + context columns ---------------------------------------------------
    ctx_candidates <- setdiff(names(df), c(keys, idx, sig))
    is_const <- function(col) {
        v <- df[[col]]
        first <- v[first_idx][gid]
        all(((v == first) %in% TRUE) | (is.na(v) & is.na(first)))
    }
    const <- vapply(ctx_candidates, is_const, logical(1))
    varying <- ctx_candidates[!const]
    is_num <- vapply(df[ctx_candidates], is.numeric, logical(1))
    vary_num <- varying[is_num[varying]]
    vary_oth <- setdiff(varying, vary_num)
    if (length(vary_oth)) {
        rlang::warn(paste0(
            "Dropped context column(s) that change within a bout: ",
            paste(vary_oth, collapse = ", "),
            ". Bout-level values would be ambiguous."
        ))
    }
    if (length(vary_num)) {
        say(
            "Dropped numeric column(s) that vary within bouts (treated as unselected sensor streams): ",
            paste(vary_num, collapse = ", ")
        )
    }
    bout_tbl <- df[first_idx, c(keys, ctx_candidates[const]), drop = FALSE]

    n_req <- nrow(plan$membership)
    n_uni <- length(unique(plan$specs$column[plan$specs$family != "catch22"]))
    say(
        "Calculating features for ",
        format(n_bouts, big.mark = ","),
        " bouts: ",
        paste(c(feature_set, if (uses_custom) "custom"), collapse = ", "),
        if (
            length(feature_set) > 1 &&
                n_uni < sum(plan$membership$family != "catch22")
        ) {
            paste0(
                " (",
                sum(plan$membership$family != "catch22") - n_uni,
                " shared features computed once)"
            )
        } else {
            ""
        }
    )

    # ---- compute, family by family -------------------------------------------------
    sp <- plan$specs
    res <- list()
    if (any(sp$family == "stat")) {
        res <- c(res, .run_stat(sp[sp$family == "stat", ], streams, bad, say))
    }
    if (any(sp$family == "seq")) {
        res <- c(res, .run_seq(sp[sp$family == "seq", ], streams, bad, say))
    }
    if (any(sp$family == "pair")) {
        res <- c(res, .run_pair(sp[sp$family == "pair", ], streams, bad, say))
    }
    if (any(sp$family == "dba")) {
        res <- c(
            res,
            .run_dba(
                sp[sp$family == "dba", ],
                sens$acc,
                streams,
                bad,
                winlen_dba,
                say
            )
        )
    }
    if (any(sp$family == "orient")) {
        res <- c(
            res,
            .run_orient(sp[sp$family == "orient", ], streams, bad, say)
        )
    }
    if (any(sp$family == "freq")) {
        res <- c(
            res,
            .run_freq(sp[sp$family == "freq", ], streams, bad, fs, say)
        )
    }

    info <- plan$membership[plan$membership$family != "catch22", ]
    if (any(sp$family == "catch22")) {
        c22 <- .run_catch22(streams, bad, catch24, say)
        res <- c(res, c22)
        c22_info <- tibble::tibble(column = names(c22))
        parts <- strsplit(c22_info$column, "__", fixed = TRUE)
        c22_info$feature_set <- "catch22"
        c22_info$family <- "catch22"
        c22_info$stream <- vapply(parts, `[`, character(1), 2)
        c22_info$feature <- vapply(parts, `[`, character(1), 3)
        info <- .ordered_info(info, c22_info, plan$membership)
    }
    if (uses_custom) {
        cust <- .run_custom(
            streams,
            bad,
            features,
            if (uses_window) winlen_dba,
            say
        )
        res <- c(res, as.list(cust))
        parts <- strsplit(names(cust), "__", fixed = TRUE)
        info <- dplyr::bind_rows(
            info,
            tibble::tibble(
                column = names(cust),
                feature_set = "custom",
                family = "custom",
                stream = vapply(parts, `[`, character(1), 2),
                feature = vapply(parts, `[`, character(1), 3)
            )
        )
    }

    col_order <- unique(info$column)
    out <- dplyr::bind_cols(bout_tbl, tibble::as_tibble(res[col_order]))
    attr(out, "feature_info") <- info
    out
}

#' Pivot the wide feature tibble into a long, theft-style layout
#'
#' Splits `<family>__<stream>__<feature>` column names into three columns.
#' Every other column (keys, context) is kept as an identifier. For set
#' membership, join feature_info(x) on `column`.
#' @param x Output of calc_features().
feature_long <- function(x) {
    n_delim <- lengths(regmatches(
        names(x),
        gregexpr("__", names(x), fixed = TRUE)
    ))
    feat <- names(x)[n_delim == 2]
    if (length(feat) == 0) {
        rlang::abort("No `<family>__<stream>__<feature>` columns found.")
    }
    tidyr::pivot_longer(
        x,
        cols = tidyselect::all_of(feat),
        names_to = c("family", "stream", "feature"),
        names_sep = "__",
        values_to = "value"
    )
}

# =============================================================================
# Planning: sensors -> concrete feature specs
# =============================================================================

# Validate the sensors table; return accelerometer descriptions, signs per
# selected column, and bookkeeping for messages.
.resolve_sensors <- function(sensors, sig, all_cols, feature_set) {
    sign <- stats::setNames(rep(1, length(sig)), sig)
    if (is.null(sensors)) {
        acc_cols <- intersect(c("acc_x", "acc_y", "acc_z"), sig)
        acc <- if (length(acc_cols)) {
            list(.acc_entry(
                "acc",
                c("acc_x", "acc_y", "acc_z"),
                acc_cols,
                roles = character()
            ))
        } else {
            list()
        }
        return(list(
            acc = acc,
            sign = as.list(sign),
            flipped = character(),
            unused = character(),
            defaulted = length(acc) > 0
        ))
    }

    if (!is.data.frame(sensors)) {
        rlang::abort(
            "`sensors` must be a data frame (sensor, type, column, role, sign)."
        )
    }
    allowed <- c("sensor", "type", "column", "role", "sign")
    extra <- setdiff(names(sensors), allowed)
    if (length(extra)) {
        rlang::abort(c(
            paste0(
                "`sensors` has unsupported column(s): ",
                paste(extra, collapse = ", ")
            ),
            i = "Allowed: sensor, type, column, role, sign.",
            i = "Mappings that vary by individual or over time are not supported: calc_features() never re-orients data. Give one mapping per call."
        ))
    }
    need <- setdiff(c("sensor", "type", "column"), names(sensors))
    if (length(need)) {
        rlang::abort(paste0(
            "`sensors` is missing column(s): ",
            paste(need, collapse = ", ")
        ))
    }
    s <- tibble::as_tibble(sensors)
    if (!"role" %in% names(s)) {
        s$role <- NA_character_
    }
    if (!"sign" %in% names(s)) {
        s$sign <- 1
    }
    s$sign[is.na(s$sign)] <- 1
    if (anyNA(s$sensor) || anyNA(s$type) || anyNA(s$column)) {
        rlang::abort("`sensors`: sensor, type and column must not be NA.")
    }
    if (anyDuplicated(s$column)) {
        rlang::abort(paste0(
            "`sensors`: each column may appear once; duplicated: ",
            paste(unique(s$column[duplicated(s$column)]), collapse = ", ")
        ))
    }
    if (!all(s$sign %in% c(1, -1))) {
        rlang::abort("`sensors$sign` must be 1 or -1.")
    }
    missing_cols <- setdiff(s$column, all_cols)
    if (length(missing_cols)) {
        rlang::abort(paste0(
            "`sensors` lists column(s) not in `data`: ",
            paste(missing_cols, collapse = ", ")
        ))
    }
    dup_role <- s[!is.na(s$role), ] |>
        dplyr::count(.data$sensor, .data$role) |>
        dplyr::filter(.data$n > 1)
    if (nrow(dup_role)) {
        rlang::abort(paste0(
            "`sensors`: a role may be used once per sensor; repeated: ",
            paste0(dup_role$sensor, ":", dup_role$role, collapse = ", ")
        ))
    }

    used <- s[s$column %in% sig, ]
    for (i in seq_len(nrow(used))) {
        sign[[used$column[i]]] <- used$sign[i]
    }

    acc <- list()
    for (nm in unique(s$sensor[tolower(s$type) == "accelerometer"])) {
        rows <- s[s$sensor == nm, ]
        pres <- rows$column[rows$column %in% sig]
        if (!length(pres)) {
            next
        }
        rr <- rows[rows$column %in% sig & !is.na(rows$role), ]
        roles <- stats::setNames(rr$column, tolower(rr$role))
        acc[[length(acc) + 1]] <- .acc_entry(nm, rows$column, pres, roles)
    }
    list(
        acc = acc,
        sign = as.list(sign),
        flipped = used$column[used$sign == -1],
        unused = setdiff(s$column, sig),
        defaulted = FALSE
    )
}

.acc_entry <- function(sensor, all_cols, cols, roles) {
    list(
        sensor = sensor,
        all_cols = all_cols,
        cols = all_cols[all_cols %in% cols],
        roles = roles,
        suffix = .axis_suffix(all_cols)
    )
}

# Short axis labels for PDBA names: strip the common prefix (acc_x -> x).
.axis_suffix <- function(cols) {
    if (length(cols) < 2) {
        return(stats::setNames(cols, cols))
    }
    ch <- strsplit(cols, "")
    len <- min(lengths(ch))
    k <- 0
    while (
        k < len && length(unique(vapply(ch, `[`, character(1), k + 1))) == 1
    ) {
        k <- k + 1
    }
    suf <- substring(cols, k + 1)
    if (any(suf == "") || anyDuplicated(suf)) {
        suf <- cols
    }
    stats::setNames(suf, cols)
}

# Pairs in rabc order (xy, yz, xz) for a three-axis sensor; otherwise combn order.
.acc_pairs <- function(a) {
    cand <- if (length(a$all_cols) == 3) {
        list(a$all_cols[1:2], a$all_cols[2:3], a$all_cols[c(1, 3)])
    } else if (length(a$all_cols) >= 2) {
        utils::combn(a$all_cols, 2, simplify = FALSE)
    } else {
        list()
    }
    Filter(function(p) all(p %in% a$cols), cand)
}

.vector_name <- function(feature, a) {
    if (length(a$cols) == 3 && length(a$all_cols) == 3) {
        return(feature)
    }
    prefix <- if (feature == "ODBA") "PDBA" else "PVeDBA"
    paste0(prefix, "_", paste(a$suffix[a$cols], collapse = ""))
}

.spec <- function(
    set,
    family,
    feature,
    stream,
    cols,
    out_feature = feature,
    sensor = NA_character_
) {
    base <- feature # tibble() evaluates in order, so capture before `feature` is reused
    tibble::tibble(
        feature_set = set,
        family = family,
        feature = out_feature,
        stream = stream,
        cols = list(cols),
        sensor = sensor,
        base_feature = base
    )
}

# Expand the requested sets into concrete specs; errors for unbuildable sets.
.plan_features <- function(feature_set, sig, acc) {
    if (!length(feature_set)) {
        # custom features only
        empty <- tibble::tibble(
            column = character(),
            feature_set = character(),
            family = character(),
            stream = character(),
            feature = character()
        )
        return(list(
            specs = empty,
            membership = empty,
            messages = character(),
            uses_acc = FALSE
        ))
    }
    reg <- .feature_registry()
    reg <- reg[reg$feature_set %in% feature_set, ]
    reg <- reg[order(match(reg$feature_set, feature_set)), ] # stable: keeps each set's order

    acc_sets <- setdiff(feature_set, c("catch22", "stats", "eikelboom"))
    if (length(acc_sets) && !length(acc)) {
        rlang::abort(c(
            paste0(
                "These sets need an accelerometer: ",
                paste(acc_sets, collapse = ", ")
            ),
            i = "Select accelerometer columns in `cols`, and declare them in `sensors` with type = \"accelerometer\" (or name them acc_x / acc_y / acc_z)."
        ))
    }

    # roles required by role-defined sets
    have_roles <- unique(unlist(lapply(acc, function(a) names(a$roles))))
    for (set in intersect(feature_set, .role_sets)) {
        tg <- reg$target[reg$feature_set == set]
        need <- unique(c(
            sub("role:", "", tg[startsWith(tg, "role:")], fixed = TRUE),
            if ("roles" %in% tg) .body_roles
        ))
        miss <- setdiff(need, have_roles)
        if (length(miss)) {
            rlang::abort(c(
                paste0(
                    "Set \"",
                    set,
                    "\" is defined on body axes and needs role(s): ",
                    paste(miss, collapse = ", ")
                ),
                i = "Pass `sensors` with a `role` column naming surge / sway / heave, and select those columns in `cols`.",
                i = "e.g. tibble::tribble(~sensor, ~type, ~column, ~role, \"acc\", \"accelerometer\", \"acc_z\", \"heave\", ...)"
            ))
        }
    }

    specs <- list()
    msgs <- character()
    for (i in seq_len(nrow(reg))) {
        r <- reg[i, ]
        tg <- r$target
        new <- if (tg == "stream") {
            lapply(sig, function(s) {
                .spec(r$feature_set, r$family, r$feature, s, s)
            })
        } else if (tg == "axis") {
            unlist(
                lapply(acc, function(a) {
                    lapply(a$cols, function(s) {
                        .spec(
                            r$feature_set,
                            r$family,
                            r$feature,
                            s,
                            s,
                            sensor = a$sensor
                        )
                    })
                }),
                recursive = FALSE
            )
        } else if (tg == "pair") {
            unlist(
                lapply(acc, function(a) {
                    lapply(.acc_pairs(a), function(p) {
                        .spec(
                            r$feature_set,
                            r$family,
                            r$feature,
                            paste(p, collapse = "_"),
                            p,
                            sensor = a$sensor
                        )
                    })
                }),
                recursive = FALSE
            )
        } else if (tg == "vector") {
            lapply(acc, function(a) {
                .spec(
                    r$feature_set,
                    r$family,
                    r$feature,
                    paste(a$cols, collapse = "_"),
                    a$cols,
                    .vector_name(r$feature, a),
                    a$sensor
                )
            })
        } else if (startsWith(tg, "role:")) {
            role <- sub("role:", "", tg, fixed = TRUE)
            lapply(
                Filter(function(a) role %in% names(a$roles), acc),
                function(a) {
                    .spec(
                        r$feature_set,
                        r$family,
                        r$feature,
                        a$roles[[role]],
                        a$roles[[role]],
                        sensor = a$sensor
                    )
                }
            )
        } else if (tg == "roles") {
            ok <- Filter(function(a) all(.body_roles %in% names(a$roles)), acc)
            if (!length(ok)) {
                msgs <- c(
                    msgs,
                    paste0(
                        r$feature_set,
                        ": ",
                        r$feature,
                        " skipped -- needs roles surge, sway and heave (see `sensors`)."
                    )
                )
            }
            lapply(ok, function(a) {
                rc <- unname(a$roles[.body_roles])
                .spec(
                    r$feature_set,
                    r$family,
                    r$feature,
                    paste(a$cols[a$cols %in% rc], collapse = "_"),
                    rc,
                    sensor = a$sensor
                )
            })
        }
        specs <- c(specs, new)
    }
    specs <- dplyr::bind_rows(specs)
    if (!nrow(specs)) {
        rlang::abort(
            "The requested sets produce no features with the selected columns."
        )
    }
    specs$column <- ifelse(
        specs$family == "catch22",
        NA_character_,
        paste(specs$family, specs$stream, specs$feature, sep = "__")
    )

    membership <- unique(specs[c(
        "column",
        "feature_set",
        "family",
        "stream",
        "feature"
    )])
    uni <- specs[
        !duplicated(ifelse(
            is.na(specs$column),
            paste0("c22:", specs$stream),
            specs$column
        )),
    ]
    list(
        specs = uni,
        membership = membership,
        messages = unique(msgs),
        uses_acc = any(specs$family %in% c("dba", "orient", "pair")) ||
            length(acc_sets) > 0
    )
}

# Put catch22 rows where catch22 sits in the requested set order.
.ordered_info <- function(info, c22_info, membership) {
    set_order <- unique(membership$feature_set)
    rows <- lapply(set_order, function(s) {
        if (s == "catch22") c22_info else info[info$feature_set == s, ]
    })
    dplyr::bind_rows(rows)
}

# =============================================================================
# Runners: each returns a named list of numeric columns (one value per bout)
# =============================================================================

.cols_to_list <- function(mat) {
    stats::setNames(
        lapply(seq_len(ncol(mat)), function(j) mat[, j]),
        colnames(mat)
    )
}

.report_family <- function(cols, label, say) {
    if (!length(cols)) {
        return(invisible(0))
    }
    m <- do.call(cbind, cols)
    .report_na(m, label, say)
}

.run_stat <- function(specs, streams, bad, say) {
    out <- list()
    for (s in unique(specs$stream)) {
        need <- unique(specs$feature[specs$stream == s])
        S <- streams[[s]]
        mat <- .fill_matrix(
            seq_along(S),
            which(!bad[[s]]),
            function(b) .stat_all(S[[b]])[need],
            need
        )
        colnames(mat) <- paste0("stat__", s, "__", need)
        out <- c(out, .cols_to_list(mat))
    }
    .report_family(out, "stat", say)
    out
}

.stat_all <- function(v) {
    m <- mean(v)
    d <- v - m
    m2 <- mean(d^2)
    s2 <- stats::var(v)
    dv <- diff(v)
    mx <- max(v)
    mn <- min(v)
    q <- stats::quantile(v, c(0.25, 0.5, 0.75), names = FALSE, type = 7)
    c(
        mean = m,
        variance = s2,
        sd = sqrt(s2),
        cv = sqrt(s2) / m,
        skew = mean(d^3) / m2^1.5,
        kurtosis = mean(d^4) / m2^2 - 3,
        max = mx,
        min = mn,
        range = mx - mn,
        norm = sqrt(sum(v^2)),
        absdiff_mean = mean(abs(dv)),
        diff_var = stats::var(dv),
        q25 = q[1],
        q50 = q[2],
        q75 = q[3]
    )
}

.run_seq <- function(specs, streams, bad, say) {
    out <- list()
    for (s in unique(specs$stream)) {
        need <- unique(specs$feature[specs$stream == s])
        S <- streams[[s]]
        mat <- .fill_matrix(
            seq_along(S),
            which(!bad[[s]]),
            function(b) .seq_all(S[[b]])[need],
            need
        )
        colnames(mat) <- paste0("seq__", s, "__", need)
        out <- c(out, .cols_to_list(mat))
    }
    .report_family(out, "seq", say)
    out
}

# lag-1 autocorrelation, as stats::acf(x, lag.max = 1)$acf[2]
.acf1 <- function(x) {
    d <- x - mean(x)
    sum(d[-1] * d[-length(d)]) / sum(d^2)
}

.seq_all <- function(v) {
    dv <- diff(v)
    sgn <- sign(v - mean(v))
    sgn <- sgn[sgn != 0]
    c(
        acf1 = .acf1(v),
        diff_acf1 = .acf1(dv),
        diff_sd = stats::sd(dv),
        mean_crossings = sum(diff(sgn) != 0),
        arv = mean(abs(v))
    )
}

# Eikelboom frequency features: FFT of the de-meaned bout, one-sided, no DC.
.fft_all <- function(v, fs) {
    n <- length(v)
    K <- floor(n / 2)
    if (K < 1) {
        return(stats::setNames(
            rep(NA_real_, length(.fft_features)),
            .fft_features
        ))
    }
    X <- stats::fft(v - mean(v))
    mod <- Mod(X[2:(K + 1)]) # bins k = 1..K
    amp <- 2 * mod / n
    if (n %% 2 == 0) {
        amp[K] <- mod[K] / n
    } # Nyquist bin is not doubled
    ord <- order(amp, decreasing = TRUE)[1:5]
    top <- amp[ord]
    frq <- ord * fs / n
    pw <- mod^2
    p <- pw / sum(pw)
    ent <- if (K > 1) -sum(p[p > 0] * log(p[p > 0])) / log(K) else NA_real_
    c(
        stats::setNames(top, paste0("fft_amp", 1:5)),
        stats::setNames(frq, paste0("fft_freq", 1:5)),
        fft_energy = sum(Mod(X[-1])^2) / n^2,
        fft_entropy = ent
    )
}

.run_pair <- function(specs, streams, bad, say) {
    out <- list()
    for (k in which(!duplicated(specs$stream))) {
        p <- specs$cols[[k]]
        st <- specs$stream[k]
        need <- unique(specs$feature[specs$stream == st])
        A <- streams[[p[1]]]
        B <- streams[[p[2]]]
        ok <- which(!(bad[[p[1]]] | bad[[p[2]]]))
        mat <- .fill_matrix(
            seq_along(A),
            ok,
            function(b) {
                a <- A[[b]]
                bb <- B[[b]]
                suppressWarnings(c(
                    cov = stats::cov(a, bb),
                    cor = stats::cor(a, bb),
                    meandiff = mean(a - bb),
                    sddiff = stats::sd(a - bb)
                )[need])
            },
            need
        )
        colnames(mat) <- paste0("pair__", st, "__", need)
        out <- c(out, .cols_to_list(mat))
    }
    .report_family(out, "pair", say)
    out
}

.run_dba <- function(specs, acc, streams, bad, winlen, say) {
    out <- list()
    for (a in acc) {
        sp <- specs[specs$sensor %in% a$sensor, ]
        if (!nrow(sp)) {
            next
        }
        cols <- a$cols
        n <- length(streams[[cols[1]]])
        fn <- function(b) {
            ax <- lapply(cols, function(cl) streams[[cl]][[b]])
            sba <- lapply(ax, .roll_mean, w = winlen)
            dyn <- Map(`-`, ax, sba)
            adb <- lapply(dyn, abs)
            per <- suppressWarnings(unlist(lapply(seq_along(cols), function(j) {
                c(
                    varsba = stats::var(sba[[j]], na.rm = TRUE),
                    vardba = stats::var(adb[[j]], na.rm = TRUE),
                    meandba = mean(adb[[j]], na.rm = TRUE),
                    maxdba = max(adb[[j]], na.rm = TRUE)
                )
            })))
            names(per) <- paste0(rep(cols, each = 4), "__", names(per))
            vec <- c(
                ODBA = mean(Reduce(`+`, adb), na.rm = TRUE),
                VeDBA = mean(
                    sqrt(Reduce(`+`, lapply(dyn, `^`, 2))),
                    na.rm = TRUE
                )
            )
            c(per, vec)
        }
        tmpl <- c(
            paste0(
                rep(cols, each = 4),
                "__",
                c("varsba", "vardba", "meandba", "maxdba")
            ),
            "ODBA",
            "VeDBA"
        )
        mat <- .fill_matrix(seq_len(n), seq_len(n), fn, tmpl)

        # mask: per-axis features by that axis, vector features by any axis
        any_bad <- Reduce(`|`, bad[cols])
        for (j in seq_along(cols)) {
            mat[bad[[cols[j]]], grepl(paste0("^", cols[j], "__"), tmpl)] <- NA
        }
        mat[any_bad, c("ODBA", "VeDBA")] <- NA

        stream_all <- paste(cols, collapse = "_")
        for (k in seq_len(nrow(sp))) {
            f <- sp$base_feature[k]
            src <- if (f %in% c("ODBA", "VeDBA")) {
                f
            } else {
                paste0(sp$stream[k], "__", f)
            }
            out[[paste0("dba__", sp$stream[k], "__", sp$feature[k])]] <- mat[,
                src
            ]
        }
    }
    .report_family(out, "dba", say)
    out
}

.run_orient <- function(specs, streams, bad, say) {
    out <- list()
    for (k in which(!duplicated(specs$stream))) {
        rc <- specs$cols[[k]] # surge, sway, heave columns
        st <- specs$stream[k]
        need <- unique(specs$feature[specs$stream == st])
        n <- length(streams[[rc[1]]])
        ok <- which(!Reduce(`|`, bad[rc]))
        mat <- .fill_matrix(
            seq_len(n),
            ok,
            function(b) {
                su <- streams[[rc[1]]][[b]]
                sw <- streams[[rc[2]]][[b]]
                hv <- streams[[rc[3]]][[b]]
                pitch <- atan(su / sqrt(sw^2 + hv^2)) * 180 / pi
                roll <- atan(sw / sqrt(su^2 + hv^2)) * 180 / pi
                c(
                    pitch_mean = mean(pitch),
                    pitch_sd = stats::sd(pitch),
                    roll_mean = mean(roll),
                    roll_sd = stats::sd(roll)
                )[need]
            },
            need
        )
        colnames(mat) <- paste0("orient__", st, "__", need)
        out <- c(out, .cols_to_list(mat))
    }
    .report_family(out, "orient", say)
    out
}

.run_freq <- function(specs, streams, bad, fs, say) {
    out <- list()
    for (s in unique(specs$stream)) {
        need <- unique(specs$feature[specs$stream == s])
        S <- streams[[s]]
        mat <- .fill_matrix(
            seq_along(S),
            which(!bad[[s]]),
            function(b) {
                v <- S[[b]]
                out <- numeric(0)
                if (any(need %in% c("freqmain", "freqamp", "entropy"))) {
                    r <- unname(suppressWarnings(.max_freq_amp(v)))
                    out <- c(
                        freqmain = r[1] * fs / length(v),
                        freqamp = r[2],
                        entropy = r[3]
                    )
                }
                if (any(need %in% .fft_features)) {
                    out <- c(out, .fft_all(v, fs))
                }
                out[need]
            },
            need
        )
        colnames(mat) <- paste0("freq__", s, "__", need)
        out <- c(out, .cols_to_list(mat))
    }
    .report_family(out, "freq", say)
    out
}

.run_catch22 <- function(streams, bad, catch24, say) {
    fn <- function(v) .catch22_vec(v, catch24)
    tmpl <- names(fn(sin(seq_len(100))))
    out <- list()
    for (col in names(streams)) {
        S <- streams[[col]]
        mat <- .fill_matrix(
            seq_along(S),
            which(!bad[[col]]),
            function(b) fn(S[[b]]),
            tmpl
        )
        colnames(mat) <- paste0("catch22__", col, "__", tmpl)
        out <- c(out, .cols_to_list(mat))
    }
    .report_family(out, "catch22", say)
    out
}

.run_custom <- function(streams, bad, features, winlen_dba, say) {
    is_cross <- vapply(features, inherits, logical(1), what = "cross_feature")
    plain <- features[!is_cross]
    cross <- features[is_cross]
    out <- list()

    # plain functions: applied to every selected stream separately, like catch22
    if (length(plain)) {
        nms <- names(plain)
        fn_all <- function(v) {
            vapply(
                plain,
                function(f) {
                    r <- tryCatch(f(v), error = function(e) NA_real_)
                    if (length(r) != 1 || !is.numeric(r)) {
                        NA_real_
                    } else {
                        as.numeric(r)
                    }
                },
                numeric(1)
            )
        }

        out <- c(
            out,
            lapply(names(streams), function(col) {
                S <- streams[[col]]
                mat <- .fill_matrix(
                    seq_along(S),
                    which(!bad[[col]]),
                    function(b) fn_all(S[[b]]),
                    nms
                )
                colnames(mat) <- paste0("custom__", col, "__", nms)
                .report_na(mat, paste0("custom (", col, ")"), say)
                tibble::as_tibble(mat)
            })
        )
    }

    # cross_feature(): one function call per bout, arguments positional over its cols
    if (length(cross)) {
        n <- length(streams[[1]])
        out <- c(
            out,
            lapply(names(cross), function(nm) {
                e <- cross[[nm]]
                cls <- e$cols
                fn <- function(b) {
                    args <- lapply(cls, function(col) streams[[col]][[b]])
                    if (isTRUE(e$window)) {
                        args <- c(args, list(winlen_dba))
                    }
                    r <- tryCatch(do.call(e$fn, args), error = function(err) {
                        NA_real_
                    })
                    if (length(r) != 1 || !is.numeric(r)) {
                        NA_real_
                    } else {
                        as.numeric(r)
                    }
                }
                bad_any <- Reduce(`|`, bad[cls])
                mat <- .fill_matrix(seq_len(n), which(!bad_any), fn, nm)
                colnames(mat) <- paste0(
                    "custom__",
                    paste(cls, collapse = "_"),
                    "__",
                    nm
                )
                .report_na(
                    mat,
                    paste0(
                        "custom (",
                        nm,
                        ", using ",
                        paste(cls, collapse = ", "),
                        ")"
                    ),
                    say
                )
                tibble::as_tibble(mat)
            })
        )
    }

    dplyr::bind_cols(out)
}

# =============================================================================
# Per-bout helpers
# =============================================================================

.catch22_vec <- function(v, catch24) {
    o <- Rcatch22::catch22_all(v, catch24 = catch24)
    stats::setNames(o$values, o$names)
}

# Centred rolling mean used to split static from dynamic acceleration.
# Same result as zoo::rollapply(v, w, mean, fill = NA, align = "center"),
# via cumulative sums (much faster; checked against zoo for odd and even w).
.roll_mean <- function(v, w) {
    n <- length(v)
    if (n < w) {
        return(rep(NA_real_, n))
    }
    cs <- cumsum(c(0, v))
    m <- (cs[(w + 1):(n + 1)] - cs[1:(n - w + 1)]) / w
    left <- floor((w - 1) / 2)
    out <- rep(NA_real_, n)
    out[(left + 1):(left + length(m))] <- m
    out
}

# Dominant frequency bin, its amplitude, and a spectral-entropy term for one
# detrended series. Follows rabc::max_freq_amp; ties take the first maximum.
.max_freq_amp <- function(v) {
    n <- length(v)
    half <- floor(n / 2)
    if (half < 1) {
        return(rep(NA_real_, 3))
    }
    freq <- abs(stats::fft(stats::lm(as.numeric(v) ~ seq_len(n))$residuals))
    ind <- which.max(freq[seq_len(half)])
    c(ind, freq[ind], entropy::entropy(freq[seq_len(half)])^2 / half)
}

# Resolve the dynamic-acceleration window to a number of samples.
.resolve_winlen <- function(w, w_s, fs, fs_ok) {
    if (!is.null(w) && !is.null(w_s)) {
        rlang::abort(
            "Give either `winlen_dba_s` (seconds) or `winlen_dba` (samples), not both."
        )
    }
    if (is.null(w) && is.null(w_s)) {
        rlang::abort(c(
            "A smoothing window is needed (dba-family features such as ODBA / VeDBA, and/or cross_feature(window = TRUE)).",
            i = "Set it in seconds (needs `fs`): `winlen_dba_s = 1, fs = 20`.",
            i = "Or in samples: `winlen_dba = 21`."
        ))
    }
    if (!is.null(w_s)) {
        if (!is.numeric(w_s) || length(w_s) != 1 || is.na(w_s) || w_s <= 0) {
            rlang::abort(
                "`winlen_dba_s` must be a single positive number of seconds."
            )
        }
        if (!fs_ok) {
            rlang::abort(
                "`winlen_dba_s` needs `fs`, the sampling rate in Hz (e.g. `fs = 20`)."
            )
        }
        n <- 2 * floor(w_s * fs / 2 + 1e-9) + 1 # odd, so the window centres on a sample
        if (n < 3) {
            rlang::abort(paste0(
                w_s,
                " s at ",
                fs,
                " Hz is under 3 samples; use a longer `winlen_dba_s`."
            ))
        }
        return(as.integer(n))
    }
    if (
        !is.numeric(w) || length(w) != 1 || is.na(w) || w < 2 || w != round(w)
    ) {
        rlang::abort(
            "`winlen_dba` must be a whole number of samples >= 2 (e.g. `winlen_dba = 21`)."
        )
    }
    as.integer(w)
}

# One-line description of the resolved window and its edge cost.
.window_note <- function(n, requested_s, fs, median_len) {
    win <- if (is.null(fs)) {
        paste0(n, " samples")
    } else {
        paste0(
            n,
            " samples (",
            format(round(n / fs, 2), nsmall = 2),
            " s at ",
            fs,
            " Hz",
            if (!is.null(requested_s)) {
                paste0("; requested ", requested_s, " s")
            },
            ")"
        )
    }
    lost <- min(1, (n - 1) / median_len)
    paste0(
        "Smoothing window: ",
        win,
        ". Edge NAs remove ",
        round(100 * lost),
        "% of a median-length bout (",
        median_len,
        " samples)."
    )
}

# Bouts x features matrix; rows not in `ok` stay NA; non-finite results -> NA.
.fill_matrix <- function(items, ok, fn, nms) {
    mat <- matrix(
        NA_real_,
        nrow = length(items),
        ncol = length(nms),
        dimnames = list(NULL, nms)
    )
    if (length(ok)) {
        mat[ok, ] <- do.call(rbind, lapply(items[ok], fn))
    }
    mat[!is.finite(mat)] <- NA_real_
    mat
}

.report_na <- function(mat, label, say) {
    n_na <- sum(!stats::complete.cases(mat))
    if (n_na > 0) {
        say(
            label,
            ": ",
            format(n_na, big.mark = ","),
            " of ",
            format(nrow(mat), big.mark = ","),
            " bouts have at least one NA feature (missing/non-finite input or an undefined feature)."
        )
    }
    invisible(n_na)
}

# =============================================================================
# Custom-feature helpers
# =============================================================================

# A single feature computed from SEVERAL streams together (see file header).
# `fn`'s arguments are matched POSITIONALLY to `cols`, in the order given here.
# `window = TRUE` appends the resolved smoothing window (samples) as one more
# positional argument, after the columns: fn(<col1>, ..., <colN>, winlen).
cross_feature <- function(cols, fn, window = FALSE) {
    if (!is.character(cols) || length(cols) == 0 || anyNA(cols)) {
        rlang::abort(
            "`cross_feature()`: `cols` must be a non-empty character vector of column names."
        )
    }
    if (!is.function(fn)) {
        rlang::abort("`cross_feature()`: `fn` must be a function.")
    }
    if (!is.logical(window) || length(window) != 1 || is.na(window)) {
        rlang::abort("`cross_feature()`: `window` must be TRUE or FALSE.")
    }
    structure(
        list(cols = cols, fn = fn, window = window),
        class = "cross_feature"
    )
}

# The exact centred rolling mean the dba family uses to split static from
# dynamic acceleration -- exported so a custom cross_feature(window = TRUE) can
# build the same split. NA at both edges ((winlen - 1) / 2 samples per side).
roll_static <- function(v, winlen) .roll_mean(v, winlen)

# A named list of functions and/or cross_feature() entries, or NULL.
.validate_custom_features <- function(features) {
    if (is.null(features)) {
        return(NULL)
    }
    if (!is.list(features) || length(features) == 0) {
        rlang::abort(
            "`features` must be a named list of functions and/or cross_feature() entries, e.g. `list(mean = mean, sd = sd)`."
        )
    }
    nms <- names(features)
    if (is.null(nms) || any(nms == "") || anyDuplicated(nms)) {
        rlang::abort(
            "`features` must be a named list with unique, non-empty names."
        )
    }
    is_cross <- vapply(features, inherits, logical(1), what = "cross_feature")
    is_fn <- vapply(features, is.function, logical(1))
    bad <- nms[!is_cross & !is_fn]
    if (length(bad)) {
        rlang::abort(paste0(
            "These `features` entries are neither functions nor cross_feature(): ",
            paste(bad, collapse = ", ")
        ))
    }
    for (i in which(is_cross)) {
        ent <- features[[i]]
        fm <- formals(ent$fn)
        need <- length(ent$cols) + if (isTRUE(ent$window)) 1L else 0L
        if (!"..." %in% names(fm) && length(fm) < need) {
            rlang::abort(paste0(
                "features$",
                nms[i],
                ": cross_feature() needs ",
                need,
                " argument(s) (",
                length(ent$cols),
                " column(s)",
                if (isTRUE(ent$window)) " + the smoothing window" else "",
                ") but its function takes only ",
                length(fm),
                "."
            ))
        }
    }
    features
}

# project_features.R -----------------------------------------------------------
# Dimension reduction (project_features()) and its plot() method for the wide
# tibble calc_features() returns: one row per bout, feature columns named
# <family>__<stream>__<feature> (older outputs: <feature_set>__<stream>__<feature>).
#
# Adapted from theftdlc::project() and theftdlc:::plot.feature_projection().
# theftdlc's calculate_features() returns a LONG feature_calculations object
# (id, group, feature_set, names, values) with one measured variable, so
# project() first pivots it wide. calc_features() is already wide -- same
# principle as the tsibble itself: genuinely multivariate, no fold/unfold --
# so here we just select the feature columns and normalise each in place;
# there is no long-to-wide reshape at all.
#
# Requires: normaliseR, tibble, dplyr, rlang, stats, and, depending on
# low_dim_method: Rtsne (tSNE), MASS (KruskalMDS, SammonMDS), umap (UMAP).
# ggplot2 for plot.calc_projection().
#
# Deliberate deviations from theftdlc, needed so every low_dim_method produces
# the SAME shape and one simple plot() method can serve all of them:
#   * PCA: theftdlc's project() puts the raw prcomp object in BOTH ModelFit
#     and ProjectedData, and its plot() then re-derives 2D coordinates via
#     broom::augment()/tidy(). We extract the first two PC scores into
#     ProjectedData ($id/$.fitted1/$.fitted2) directly and keep the full
#     prcomp object only in ModelFit -- no broom dependency needed.
#   * KruskalMDS / SammonMDS: theftdlc's project() attaches the row id to
#     `fits$id` (the MODEL FIT, i.e. the isoMDS/sammon() return value) rather
#     than `projected$id` -- an apparent bug: ProjectedData for these two
#     methods ends up with no id column at all, and its own plot() method
#     (which reads $ProjectedData$.fitted1) would also fail on them, since
#     those columns are unnamed there too. Both are fixed here.
#
# Column selection: `feature_set` filters by set MEMBERSHIP, read from the
# "feature_info" attribute calc_features() attaches (a column can belong to
# several sets, e.g. ODBA is in collins, harvey_carroll, yu_full ...). So
# `feature_set = "collins"` gives exactly Collins' 10 features, and
# `feature_set = c("collins", "arrieta")` gives their union, no duplicates.
# For outputs without the attribute (older versions of calc_features(), or a
# tibble rebuilt from a CSV), the first part of the column name is used as
# the set, as before. NULL uses every feature column present. `cols` selects exact feature column
# names instead (e.g. for a UI that lets someone tick individual features, not
# whole sets); give one or the other, not both. `Data` in the returned object
# is always the untouched, full `data` you passed in regardless of which
# columns were modelled -- so plot()'s `colour_by` can reference ANY original
# column, feature or not, even one that wasn't part of this particular run.

# Split calc_features() output into feature columns (matched by the
# <feature_set>__<stream>__<feature> naming, i.e. exactly two "__" delimiters)
# and everything else (key + context columns).
.split_feature_cols <- function(nms) {
    n_delim <- lengths(regmatches(nms, gregexpr("__", nms, fixed = TRUE)))
    is_feat <- n_delim == 2
    info <- NULL
    if (any(is_feat)) {
        parts <- strsplit(nms[is_feat], "__", fixed = TRUE)
        info <- tibble::tibble(
            column = nms[is_feat],
            feature_set = vapply(parts, `[`, character(1), 1),
            stream = vapply(parts, `[`, character(1), 2),
            feature = vapply(parts, `[`, character(1), 3)
        )
    }
    list(feature = nms[is_feat], meta = nms[!is_feat], info = info)
}

# Which feature set(s) each feature column belongs to: from the
# "feature_info" attribute when present, else the first part of the name.
.feature_membership <- function(data) {
    parts <- .split_feature_cols(names(data))
    info <- attr(data, "feature_info", exact = TRUE)
    if (!is.null(info) && all(c("column", "feature_set") %in% names(info))) {
        info <- info[info$column %in% parts$feature, c("column", "feature_set")]
        if (nrow(info)) return(unique(info))
    }
    if (is.null(parts$info)) {
        return(tibble::tibble(column = character(), feature_set = character()))
    }
    parts$info[c("column", "feature_set")]
}

#' Reduce a calc_features() feature tibble to two dimensions
#'
#' @param data           Output of calc_features(): one row per bout/sample,
#'                       feature columns named <feature_set>__<stream>__<feature>.
#' @param feature_set    Character vector of feature_set(s) to project (e.g.
#'                       "catch22", or c("catch22", "rabc_time2")). NULL
#'                       (default) uses every feature column present. Give
#'                       this or `cols`, not both.
#' @param cols           Character vector of EXACT feature column names to
#'                       project (e.g. for per-feature, not per-set, selection).
#'                       Give this or `feature_set`, not both.
#' @param norm_method    One of "zScore", "Sigmoid", "RobustSigmoid", "MinMax",
#'                       "MaxAbs" (normaliseR::normalise()), applied per feature,
#'                       across bouts.
#' @param unit_int       Also rescale into the unit interval [0, 1] afterwards.
#' @param low_dim_method One of "PCA", "tSNE", "ClassicalMDS", "KruskalMDS",
#'                       "SammonMDS", "UMAP".
#' @param na_removal     "feature" (default): drop any feature column with an
#'                       NA, keeping every bout. "sample": keep every feature
#'                       column, drop any bout with an NA.
#' @param seed           Integer for set.seed(), used by the stochastic methods
#'                       (tSNE, UMAP).
#' @param verbose        Print progress / NA-removal messages.
#' @param ...            Passed on to the underlying method: stats::prcomp(),
#'                       Rtsne::Rtsne(), stats::cmdscale(), MASS::isoMDS(),
#'                       MASS::sammon(), or umap::umap().
#'
#' @return An object of class "calc_projection", a list with:
#'   Data           the `data` you passed in, unmodified
#'   Meta           the non-feature (key/context) columns that survived NA
#'                  removal, plus `..row_id..` (matches ProjectedData$id)
#'   ModelData      the normalised, NA-handled feature matrix actually modelled
#'   ProjectedData  tibble(id, .fitted1, .fitted2), one row per surviving bout
#'   ModelFit       the raw fit object (prcomp / Rtsne / cmdscale / isoMDS /
#'                  sammon / umap)
#'   LowDimMethod   the `low_dim_method` used
#' @examples
#' calc_features(vultures_tsbl, c(acc_x, acc_y, acc_z), feature_set = "catch22") |>
#'   project_features(norm_method = "RobustSigmoid", unit_int = TRUE,
#'                    low_dim_method = "PCA") |>
#'   plot()
project_features <- function(
    data,
    feature_set = NULL,
    cols = NULL,
    norm_method = c("zScore", "Sigmoid", "RobustSigmoid", "MinMax", "MaxAbs"),
    unit_int = FALSE,
    low_dim_method = c(
        "PCA",
        "tSNE",
        "ClassicalMDS",
        "KruskalMDS",
        "SammonMDS",
        "UMAP"
    ),
    na_removal = c("feature", "sample"),
    seed = 123,
    verbose = TRUE,
    ...
) {
    say <- function(...) if (verbose) rlang::inform(paste0(...))

    if (!is.data.frame(data)) {
        rlang::abort(
            "`data` must be a tibble/data.frame -- the output of calc_features()."
        )
    }
    norm_method <- match.arg(norm_method)
    low_dim_method <- match.arg(low_dim_method)
    na_removal <- match.arg(na_removal)

    # ---- which columns are features, and which feature_set(s) do they belong to --
    parts <- .split_feature_cols(names(data))
    if (length(parts$feature) == 0) {
        rlang::abort(
            "No `<family>__<stream>__<feature>` columns found. Is `data` the output of calc_features()?"
        )
    }
    if (!is.null(feature_set) && !is.null(cols)) {
        rlang::abort("Give either `feature_set` or `cols`, not both.")
    }
    if (!is.null(cols)) {
        unknown <- setdiff(cols, parts$feature)
        if (length(unknown)) {
            rlang::abort(c(
                paste0(
                    "`cols` not found among feature columns: ",
                    paste(unknown, collapse = ", ")
                ),
                i = "Feature columns are named <family>__<stream>__<feature>."
            ))
        }
        parts$info <- parts$info[parts$info$column %in% cols, ]
        parts$feature <- cols
    } else if (!is.null(feature_set)) {
        member <- .feature_membership(data)
        available <- unique(member$feature_set)
        unknown <- setdiff(feature_set, available)
        if (length(unknown)) {
            rlang::abort(c(
                paste0(
                    "`feature_set` not present in `data`: ",
                    paste(unknown, collapse = ", ")
                ),
                i = paste0("Available: ", paste(available, collapse = ", "))
            ))
        }
        keep_cols <- unique(member$column[member$feature_set %in% feature_set])
        parts$feature <- parts$feature[parts$feature %in% keep_cols] # original column order
        parts$info <- parts$info[parts$info$column %in% parts$feature, ]
    }
    if (length(parts$feature) < 2) {
        rlang::abort(
            "Need at least 2 feature columns to project to 2 dimensions."
        )
    }

    # ---- normalise every selected feature column in place; no reshape needed -----
    meta <- tibble::as_tibble(data[parts$meta])
    meta$..row_id.. <- as.character(seq_len(nrow(data)))

    wide <- as.data.frame(data[parts$feature])
    wide[] <- lapply(wide, function(v) {
        normaliseR::normalise(
            as.numeric(v),
            norm_method = norm_method,
            unit_int = unit_int
        )
    })
    rownames(wide) <- meta$..row_id..

    n_features <- ncol(wide)
    n_samples <- nrow(wide)

    if (na_removal == "feature") {
        keep_col <- vapply(wide, function(v) !anyNA(v), logical(1))
        wide <- wide[, keep_col, drop = FALSE]
    } else {
        keep_row <- stats::complete.cases(wide)
        wide <- wide[keep_row, , drop = FALSE]
        meta <- meta[keep_row, , drop = FALSE]
    }

    n_features_omitted <- n_features - ncol(wide)
    n_samples_omitted <- n_samples - nrow(wide)
    if (n_features_omitted > 0) {
        say(n_features_omitted, " feature(s) omitted due to NAs.")
    }
    if (n_samples_omitted > 0) {
        say(n_samples_omitted, " bout(s) omitted due to NAs.")
    }
    if (ncol(wide) < 2) {
        rlang::abort(
            "Fewer than 2 features survive NA removal; try `na_removal = \"sample\"`, or a different `feature_set`."
        )
    }
    if (nrow(wide) < 3) {
        rlang::abort("Fewer than 3 bouts survive NA removal.")
    }

    # ---- dimension reduction -------------------------------------------------------
    set.seed(seed)
    ids <- rownames(wide)

    if (low_dim_method == "PCA") {
        fit <- stats::prcomp(wide, center = FALSE, scale. = FALSE, ...)
        projected <- tibble::tibble(
            id = ids,
            .fitted1 = fit$x[, 1],
            .fitted2 = fit$x[, 2]
        )
    } else if (low_dim_method == "tSNE") {
        fit <- Rtsne::Rtsne(
            as.matrix(wide),
            dims = 2,
            check_duplicates = FALSE,
            ...
        )
        projected <- tibble::tibble(
            id = ids,
            .fitted1 = fit$Y[, 1],
            .fitted2 = fit$Y[, 2]
        )
    } else if (low_dim_method == "ClassicalMDS") {
        fit <- stats::cmdscale(stats::dist(wide), k = 2, ...)
        projected <- tibble::tibble(
            id = ids,
            .fitted1 = fit[, 1],
            .fitted2 = fit[, 2]
        )
    } else if (low_dim_method == "KruskalMDS") {
        fit <- MASS::isoMDS(stats::dist(wide), k = 2, ...)
        projected <- tibble::tibble(
            id = ids,
            .fitted1 = fit$points[, 1],
            .fitted2 = fit$points[, 2]
        )
    } else if (low_dim_method == "SammonMDS") {
        fit <- MASS::sammon(stats::dist(wide), k = 2, ...)
        projected <- tibble::tibble(
            id = ids,
            .fitted1 = fit$points[, 1],
            .fitted2 = fit$points[, 2]
        )
    } else {
        fit <- umap::umap(as.matrix(wide), n_components = 2, ...)
        projected <- tibble::tibble(
            id = ids,
            .fitted1 = fit$layout[, 1],
            .fitted2 = fit$layout[, 2]
        )
    }

    structure(
        list(
            Data = data,
            Meta = meta,
            ModelData = wide,
            ProjectedData = projected,
            ModelFit = fit,
            LowDimMethod = low_dim_method
        ),
        class = "calc_projection"
    )
}

#' Plot a calc_projection object
#'
#' @param x            A "calc_projection" object from project_features().
#' @param colour_by    Name of a column to colour points by: a column from
#'                     `x$Meta` (context, e.g. "label", "id" -- discrete
#'                     palette), OR the name of any column in `x$Data` (e.g. a
#'                     feature, whether or not it was used in this particular
#'                     projection -- continuous gradient). Defaults to "label"
#'                     if present, else no colour.
#' @param show_covariance Draw a covariance ellipse per colour group. Ignored
#'                     if `colour_by` is NULL or continuous (a feature value).
#' @param ...          Unused; present for S3 consistency.
#' @return A ggplot object.
#' @export
plot.calc_projection <- function(
    x,
    colour_by = NULL,
    show_covariance = TRUE,
    ...
) {
    if (!inherits(x, "calc_projection")) {
        rlang::abort(
            "`x` must be a calc_projection object (see project_features())."
        )
    }

    fits <- x$ProjectedData
    meta_cols <- setdiff(names(x$Meta), "..row_id..")

    if (is.null(colour_by) && "label" %in% meta_cols) {
        colour_by <- "label"
    }

    is_continuous <- FALSE
    if (!is.null(colour_by)) {
        if (colour_by %in% meta_cols) {
            grp <- stats::setNames(
                x$Meta[c("..row_id..", colour_by)],
                c("id", "group_id")
            )
            fits <- dplyr::inner_join(fits, grp, by = "id")
            fits$group_id <- factor(fits$group_id)
        } else if (colour_by %in% names(x$Data)) {
            row_idx <- as.integer(x$Meta$..row_id..)
            val <- x$Data[[colour_by]][row_idx]
            fits$group_id <- val[match(fits$id, x$Meta$..row_id..)]
            is_continuous <- is.numeric(fits$group_id)
            if (!is_continuous) fits$group_id <- factor(fits$group_id)
        } else {
            rlang::abort(c(
                paste0("`colour_by` not found: ", colour_by),
                i = paste0(
                    "Context columns: ",
                    paste(meta_cols, collapse = ", ")
                ),
                i = "Or any column name from the original data passed to project_features() (e.g. a feature, for a continuous gradient)."
            ))
        }
    }

    if (x$LowDimMethod == "PCA") {
        pct <- round(100 * x$ModelFit$sdev^2 / sum(x$ModelFit$sdev^2))
        xlab <- paste0("PC1 (", pct[1], "%)")
        ylab <- paste0("PC2 (", pct[2], "%)")
    } else {
        xlab <- "Dimension 1"
        ylab <- "Dimension 2"
    }

    pt_size <- if (nrow(fits) > 200) 1.5 else 2.25
    p <- ggplot2::ggplot(
        fits,
        ggplot2::aes(x = .data$.fitted1, y = .data$.fitted2)
    )

    if (!is.null(colour_by) && is_continuous) {
        p <- p +
            ggplot2::geom_point(
                ggplot2::aes(colour = .data$group_id),
                size = pt_size
            ) +
            ggplot2::scale_colour_viridis_c() +
            ggplot2::labs(colour = colour_by)
    } else if (!is.null(colour_by)) {
        n_levels <- nlevels(fits$group_id)
        use_viridis <- n_levels > 8 # RColorBrewer's Dark2 supports at most 8 distinct colours;
        # beyond that it silently assigns NA, and ggplot2 then drops
        # those points entirely -- viridis scales to any n instead.
        if (isTRUE(show_covariance)) {
            p <- p +
                ggplot2::stat_ellipse(
                    ggplot2::aes(fill = .data$group_id),
                    geom = "polygon",
                    alpha = 0.2
                ) +
                ggplot2::guides(fill = "none") +
                (if (use_viridis) {
                    ggplot2::scale_fill_viridis_d()
                } else {
                    ggplot2::scale_fill_brewer(palette = "Dark2")
                })
        }
        p <- p +
            ggplot2::geom_point(
                ggplot2::aes(colour = .data$group_id),
                size = pt_size
            ) +
            (if (use_viridis) {
                ggplot2::scale_colour_viridis_d()
            } else {
                ggplot2::scale_colour_brewer(palette = "Dark2")
            }) +
            ggplot2::labs(colour = colour_by)
    } else {
        p <- p + ggplot2::geom_point(size = pt_size, colour = "black")
    }

    p +
        ggplot2::labs(
            title = "Low-dimensional projection of bouts",
            subtitle = paste0("Dimension reduction: ", x$LowDimMethod),
            x = xlab,
            y = ylab
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            panel.grid.minor = ggplot2::element_blank(),
            legend.position = "bottom"
        )
}

# explore_features.R -----------------------------------------------------------
# Interactive Shiny app wrapping project_features() / plot.calc_projection(),
# inspired by rabc::plot_UMAP() but adapted to calc_features()'s single wide
# tibble and project_features()'s full range of low_dim_method options.
#
# Requires: shiny, bslib, plus everything project_features.R requires. Source
# project_features.R (and calculate_features.R, which it needs) first.
#
# Structure (one page, not rabc's three tabs -- merged into one control panel):
#   Sidebar -- feature sets (checkboxes, populated from the feature_set
#     prefixes actually present in `data`'s column names), an optional
#     fine-tune multi-select for individual features (grouped by set), then
#     normalisation / NA-handling / low_dim_method controls, and a "Project"
#     button. Colour-by and the covariance-ellipse toggle sit below the
#     button: they restyle the CURRENT projection instantly (cheap), while
#     everything above the button requires a fresh Project click, since
#     re-running tSNE/UMAP on every checkbox tick would be slow and annoying.
#   Main panel -- the plot, reusing plot.calc_projection() as-is (no
#     duplicated plotting code), plus a small summary of the current
#     settings and any NA-removal messages from project_features().
#
# "Colour by" lists BOTH context columns (id, label, attachment, segment, ...
#   -- discrete palette) and every feature column in the FULL original
#   dataset (continuous gradient), not just whichever ones are currently
#   toggled on for the projection -- project_features() always keeps `Data`
#   as the untouched original tibble for exactly this: seeing how a feature
#   maps onto a projection built from a DIFFERENT set of features.
#
# Feature/hyperparameter scope for v1: per-method hyperparameters (UMAP's
# n_neighbors, tSNE's perplexity, ...) are deliberately not exposed yet --
# each method's defaults are used via project_features()'s `...`. Add a
# conditionalPanel per method once the feature-selection workflow is settled.

#' Launch the feature projection explorer
#'
#' @param data Output of calc_features(): one row per bout, feature columns
#'   named <family>__<stream>__<feature>. Feature-set checkboxes come from
#'   its "feature_info" attribute (a feature can be in several sets).
#' @return A shiny.appobj (call the function to launch it).
#' @examples
#' calc_features(vultures_tsbl, c(acc_x, acc_y, acc_z),
#'              feature_set = c("catch22", "rabc_time2")) |>
#'   explore_features()
explore_features <- function(data) {
    if (
        !exists("project_features", mode = "function") ||
            !exists(".split_feature_cols", mode = "function")
    ) {
        rlang::abort(
            "Source project_features.R first (this builds on project_features() and plot.calc_projection())."
        )
    }
    if (!is.data.frame(data)) {
        rlang::abort(
            "`data` must be a tibble/data.frame -- the output of calc_features()."
        )
    }

    parts <- .split_feature_cols(names(data))
    if (length(parts$feature) == 0) {
        rlang::abort(
            "No `<family>__<stream>__<feature>` columns found. Is `data` the output of calc_features()?"
        )
    }

    member <- .feature_membership(data) # column x feature_set (from feature_info)
    feature_sets <- unique(member$feature_set)
    meta_cols <- parts$meta
    # fine-tune list grouped by family (first part of the name): a feature shared
    # by several sets appears once
    families <- unique(parts$info$feature_set)
    by_family <- split(
        parts$info$column,
        factor(parts$info$feature_set, levels = families)
    )
    colour_choices <- list(Context = meta_cols, Features = sort(parts$feature))
    default_colour <- if ("label" %in% meta_cols) "label" else meta_cols[1]

    ui <- bslib::page_sidebar(
        title = "Feature projection explorer",
        sidebar = bslib::sidebar(
            width = 340,
            shiny::checkboxGroupInput(
                "feature_sets",
                "Feature sets",
                choices = feature_sets,
                selected = feature_sets
            ),
            shiny::checkboxInput(
                "fine_tune",
                "Fine-tune individual features",
                value = FALSE
            ),
            shiny::conditionalPanel(
                "input.fine_tune",
                shiny::selectizeInput(
                    "features_custom",
                    NULL,
                    choices = by_family,
                    multiple = TRUE,
                    options = list(
                        placeholder = "Defaults to all features in the checked sets above"
                    )
                )
            ),
            shiny::hr(),
            shiny::selectInput(
                "norm_method",
                "Normalisation",
                choices = c(
                    "zScore",
                    "Sigmoid",
                    "RobustSigmoid",
                    "MinMax",
                    "MaxAbs"
                )
            ),
            shiny::checkboxInput(
                "unit_int",
                "Rescale to [0, 1] afterwards",
                value = FALSE
            ),
            shiny::selectInput(
                "na_removal",
                "NA handling",
                choices = c(
                    "Drop features with any NA" = "feature",
                    "Drop bouts with any NA" = "sample"
                )
            ),
            shiny::selectInput(
                "low_dim_method",
                "Method",
                choices = c(
                    "PCA",
                    "tSNE",
                    "ClassicalMDS",
                    "KruskalMDS",
                    "SammonMDS",
                    "UMAP"
                )
            ),
            shiny::hr(),
            shiny::actionButton("go", "Project", class = "btn-primary"),
            shiny::hr(),
            shiny::selectInput(
                "colour_by",
                "Colour by",
                choices = colour_choices,
                selected = default_colour
            ),
            shiny::checkboxInput(
                "show_covariance",
                "Show covariance ellipse",
                value = TRUE
            )
        ),
        shiny::uiOutput("summary"),
        shiny::plotOutput("proj_plot", height = "600px")
    )

    server <- function(input, output, session) {
        chosen_cols <- shiny::reactive({
            if (isTRUE(input$fine_tune) && length(input$features_custom) > 0) {
                input$features_custom
            } else {
                unique(member$column[
                    member$feature_set %in% input$feature_sets
                ])
            }
        })

        projection <- shiny::eventReactive(input$go, {
            cols <- chosen_cols()
            shiny::validate(shiny::need(
                length(cols) >= 2,
                "Select at least 2 features to project."
            ))

            msgs <- character()
            proj <- tryCatch(
                withCallingHandlers(
                    project_features(
                        data,
                        cols = cols,
                        norm_method = input$norm_method,
                        unit_int = input$unit_int,
                        low_dim_method = input$low_dim_method,
                        na_removal = input$na_removal,
                        verbose = TRUE
                    ),
                    message = function(m) {
                        msgs <<- c(msgs, trimws(conditionMessage(m)))
                        invokeRestart("muffleMessage")
                    }
                ),
                error = function(e) {
                    shiny::validate(shiny::need(FALSE, conditionMessage(e)))
                }
            )
            list(proj = proj, messages = msgs)
        })

        output$proj_plot <- shiny::renderPlot(
            {
                shiny::req(projection())
                plot(
                    projection()$proj,
                    colour_by = input$colour_by,
                    show_covariance = input$show_covariance
                )
            },
            height = function() session$clientData$output_proj_plot_width * 0.75
        )

        output$summary <- shiny::renderUI({
            shiny::req(projection())
            p <- projection()$proj
            shiny::tagList(
                shiny::p(
                    shiny::strong(paste0(
                        nrow(p$ModelData),
                        " bouts \u00d7 ",
                        ncol(p$ModelData),
                        " features"
                    )),
                    " -- ",
                    p$LowDimMethod,
                    ", ",
                    input$norm_method,
                    if (isTRUE(input$unit_int)) " (rescaled to [0, 1])" else ""
                ),
                if (length(projection()$messages)) {
                    shiny::tags$ul(lapply(
                        projection()$messages,
                        shiny::tags$li
                    ))
                }
            )
        })
    }

    shiny::shinyApp(ui, server)
}

# winlen_sensitivity.R --------------------------------------------------------
# How sensitive is ODBA (or PDBA) to the smoothing window?
#
# Runs calc_features(feature_set = "rabc_time1") over a grid of windows
# given in SECONDS, then asks three questions of the resulting ODBA values:
#
#   1. Magnitude   How does ODBA change with the window, per behaviour?
#                  (median and IQR by label; median relative to a reference)
#   2. Separation  Does ODBA still tell behaviours apart?
#                  * Kruskal-Wallis epsilon^2 across all labels (0-1)
#                  * pairwise separability = max(AUC, 1 - AUC) for every pair
#                    of labels (0.5 = no separation, 1 = perfect); the median
#                    and minimum over pairs are summarised
#   3. Stability   Does the ORDER of bouts by ODBA change? Spearman rank
#                  correlation of every window with the reference window.
#                  Also reported: the share of each bout lost to edge NAs.
#
# Requires: calculate_features.R sourced first (for calc_features()), plus dplyr,
#           tidyr, purrr,
#           tidyselect, tsibble, rlang; ggplot2 for the plots (optional).

#' Sensitivity of ODBA/PDBA to the smoothing window
#'
#' @param data      A tsibble in the tidy sensor format (see calc_features).
#' @param cols      Accelerometer axes, tidyselect. With fewer than three axes
#'                  PDBA is analysed instead of ODBA, as in calc_features().
#' @param fs        Sampling rate in Hz (required).
#' @param winlen_s  Windows to test, in seconds. Windows that resolve to the
#'                  same odd number of samples are collapsed, and windows as
#'                  long as a typical bout are dropped (both with a message).
#' @param ref_s     Reference window (seconds) for the relative-magnitude and
#'                  rank-correlation comparisons. Added to the grid if missing.
#' @param label_col Name of the behaviour label column. If absent, only
#'                  magnitude and stability are analysed.
#' @param sensors   Optional sensor table passed to calc_features() (see
#'                  calculate_features.R); NULL treats acc_x/y/z as one
#'                  accelerometer.
#' @param min_n     Labels with fewer bouts than this are left out of the
#'                  separation statistics (noisy with few bouts).
#' @param plot      Build ggplot2 figures (needs ggplot2).
#' @param verbose   Print progress messages and the overview table.
#'
#' @return A list: `overview` (one row per window), `by_label`, `pairs`
#'   (pairwise separability), `rank_cor` (window x window Spearman matrix),
#'   `data` (bout-level values for every window), `plots`, `metric`.
#' @examples
#' res <- winlen_sensitivity(vultures_tsbl, fs = 20)
#' res$overview
#' res$plots$separation
winlen_sensitivity <- function(
    data,
    cols = c(acc_x, acc_y, acc_z),
    fs,
    winlen_s = c(0.25, 0.5, 1, 1.5, 2, 3),
    ref_s = 1,
    label_col = "label",
    min_n = 10,
    sensors = NULL,
    plot = TRUE,
    verbose = TRUE
) {
    say <- function(...) if (verbose) rlang::inform(paste0(...))

    # ---- checks ----------------------------------------------------------------
    if (!exists("calc_features", mode = "function")) {
        rlang::abort(
            "Source calculate_features.R first (this builds on calc_features())."
        )
    }
    if (!tsibble::is_tsibble(data)) {
        rlang::abort("`data` must be a tsibble in the tidy sensor format.")
    }
    if (
        missing(fs) ||
            !is.numeric(fs) ||
            length(fs) != 1 ||
            is.na(fs) ||
            fs <= 0
    ) {
        rlang::abort(
            "`fs`, the sampling rate in Hz, is required (e.g. `fs = 20`)."
        )
    }
    if (!is.numeric(winlen_s) || any(is.na(winlen_s)) || any(winlen_s <= 0)) {
        rlang::abort("`winlen_s` must be positive numbers of seconds.")
    }

    keys <- tsibble::key_vars(data)
    df <- tibble::as_tibble(data)
    sig <- names(tidyselect::eval_select(rlang::enquo(cols), df))
    has_label <- label_col %in% names(df)
    if (!has_label) {
        say(
            "No `",
            label_col,
            "` column: separation statistics will be skipped."
        )
    }

    n_per_bout <- dplyr::count(df, dplyr::across(dplyr::all_of(keys)))$n
    med_len <- stats::median(n_per_bout)

    # ---- window grid -----------------------------------------------------------
    ws <- sort(unique(c(winlen_s, ref_s)))
    n_s <- vapply(
        ws,
        function(s) {
            tryCatch(.resolve_winlen(NULL, s, fs, TRUE), error = function(e) {
                NA_integer_
            })
        },
        integer(1)
    )

    too_short <- is.na(n_s) # under 3 samples at this fs
    pref <- order(ws != ref_s) # reference first, so it survives de-duplication
    first_hit <- logical(length(ws))
    first_hit[pref] <- !duplicated(n_s[pref])
    dup <- !first_hit & !too_short
    too_long <- !is.na(n_s) & n_s >= med_len # window as long as the bout: all NA

    if (any(ws[too_short] == ref_s) || any(ws[too_long] == ref_s)) {
        rlang::abort(
            "The reference window is unusable at this `fs` / bout length; choose another `ref_s`."
        )
    }
    if (any(too_short)) {
        say(
            "Dropped (under 3 samples at ",
            fs,
            " Hz): ",
            paste(ws[too_short], collapse = ", "),
            " s"
        )
    }
    if (any(dup)) {
        say(
            "Dropped (same odd sample count as another window): ",
            paste(ws[dup], collapse = ", "),
            " s"
        )
    }
    if (any(too_long)) {
        say(
            "Dropped (as long as a median bout, ",
            med_len,
            " samples): ",
            paste(ws[too_long], collapse = ", "),
            " s"
        )
    }
    use <- !too_short & !dup & !too_long
    ws <- ws[use]
    n_s <- n_s[use]
    if (length(ws) < 2) {
        rlang::abort("Fewer than two usable windows; widen `winlen_s`.")
    }
    say(
        "Windows tested: ",
        paste0(ws, " s (", n_s, " samples)", collapse = ", ")
    )

    # ---- ODBA / PDBA at every window -------------------------------------------
    long <- purrr::map2(ws, n_s, function(s, k) {
        f <- calc_features(
            data,
            dplyr::all_of(sig),
            feature_set = "rabc_time1",
            sensors = sensors,
            winlen_dba_s = s,
            fs = fs,
            verbose = FALSE
        )
        dyn <- grep("^dba__.*__(ODBA|PDBA_.+)$", names(f), value = TRUE)
        out <- f[keys]
        if (has_label) {
            out$label <- as.character(f[[label_col]])
        }
        out$winlen_s <- s
        out$winlen_samples <- k
        out$metric <- sub(".*__", "", dyn)
        out$dba <- f[[dyn]]
        out
    }) |>
        purrr::list_rbind()
    metric <- unique(long$metric)

    # ---- 1. magnitude ------------------------------------------------------------
    by_label <- NULL
    if (has_label) {
        by_label <- long |>
            dplyr::filter(!is.na(dba), !is.na(label)) |>
            dplyr::summarise(
                n = dplyr::n(),
                median = stats::median(dba),
                q25 = stats::quantile(dba, 0.25),
                q75 = stats::quantile(dba, 0.75),
                .by = c(label, winlen_s, winlen_samples)
            ) |>
            dplyr::mutate(
                median_rel_ref = median / median[winlen_s == ref_s],
                .by = label
            ) |>
            dplyr::arrange(label, winlen_s)
    }

    # ---- 3. stability (rank correlation between windows) ---------------------------
    wide <- tidyr::pivot_wider(
        long[c(keys, "winlen_s", "dba")],
        id_cols = dplyr::all_of(keys),
        names_from = winlen_s,
        values_from = dba
    )
    rank_cor <- stats::cor(
        as.matrix(wide[as.character(ws)]),
        method = "spearman",
        use = "pairwise.complete.obs"
    )

    # ---- 2. separation -------------------------------------------------------------
    pairs <- NULL
    sep <- tibble::tibble(
        winlen_s = ws,
        epsilon2 = NA_real_,
        median_pair_sep = NA_real_,
        min_pair_sep = NA_real_
    )
    if (has_label) {
        res <- lapply(ws, function(s) {
            .separation(long[long$winlen_s == s, ], min_n)
        })
        sep$epsilon2 <- vapply(res, function(r) r$eps2, numeric(1))
        pairs <- purrr::map2(res, ws, function(r, s) {
            if (nrow(r$pairs)) dplyr::mutate(r$pairs, winlen_s = s)
        }) |>
            purrr::list_rbind()
        if (nrow(pairs)) {
            ps <- pairs |>
                dplyr::summarise(
                    median_pair_sep = stats::median(sep),
                    min_pair_sep = min(sep),
                    .by = winlen_s
                )
            sep$median_pair_sep <- ps$median_pair_sep[match(ws, ps$winlen_s)]
            sep$min_pair_sep <- ps$min_pair_sep[match(ws, ps$winlen_s)]
        }
    }

    overview <- tibble::tibble(
        winlen_s = ws,
        samples = n_s,
        realised_s = n_s / fs,
        edge_loss_pct = pmin(100, 100 * (n_s - 1) / med_len),
        median_dba = vapply(
            ws,
            function(s) {
                stats::median(long$dba[long$winlen_s == s], na.rm = TRUE)
            },
            numeric(1)
        ),
        rank_cor_ref = unname(rank_cor[, as.character(ref_s)])
    ) |>
        dplyr::left_join(sep, by = "winlen_s") |>
        dplyr::mutate(is_ref = winlen_s == ref_s)
    names(overview)[names(overview) == "median_dba"] <- paste0(
        "median_",
        metric
    )

    if (verbose) {
        rlang::inform(paste0(
            "Overview (",
            metric,
            "; reference window = ",
            ref_s,
            " s):"
        ))
        print(overview, n = Inf, width = Inf)
    }

    plots <- NULL
    if (plot) {
        if (requireNamespace("ggplot2", quietly = TRUE)) {
            plots <- .plot_winlen(by_label, overview, ref_s, metric)
        } else {
            say("ggplot2 not installed; skipping plots.")
        }
    }

    list(
        overview = overview,
        by_label = by_label,
        pairs = pairs,
        rank_cor = rank_cor,
        data = long,
        plots = plots,
        metric = metric
    )
}

# ---- helpers -----------------------------------------------------------------

# Direction-free separability of two samples: max(AUC, 1 - AUC).
# AUC = P(a > b) via the rank-sum (Mann-Whitney U) identity.
.separability <- function(a, b) {
    r <- rank(c(a, b))
    n1 <- length(a)
    n2 <- length(b)
    auc <- (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * n2)
    max(auc, 1 - auc)
}

# Kruskal-Wallis epsilon^2 (0-1) across labels, plus pairwise separability.
.separation <- function(d, min_n) {
    d <- d[!is.na(d$dba) & !is.na(d$label), ]
    cnt <- table(d$label)
    ok <- sort(names(cnt)[cnt >= min_n])
    d <- d[d$label %in% ok, ]
    if (length(ok) < 2) {
        return(list(eps2 = NA_real_, pairs = tibble::tibble()))
    }

    H <- unname(stats::kruskal.test(dba ~ label, data = d)$statistic)
    n <- nrow(d)
    eps2 <- H * (n + 1) / (n^2 - 1) # Tomczak & Tomczak (2014)

    pairs <- utils::combn(ok, 2, simplify = FALSE) |>
        purrr::map(function(p) {
            tibble::tibble(
                label_a = p[1],
                label_b = p[2],
                sep = .separability(
                    d$dba[d$label == p[1]],
                    d$dba[d$label == p[2]]
                )
            )
        }) |>
        purrr::list_rbind()
    list(eps2 = eps2, pairs = pairs)
}

.plot_winlen <- function(by_label, overview, ref_s, metric) {
    pal <- c(blue = "#0072B2", vermillion = "#D55E00", green = "#009E73")
    vref <- ggplot2::geom_vline(
        xintercept = ref_s,
        linetype = "dotted",
        colour = "black"
    )
    theme <- ggplot2::theme_bw()
    xs <- ggplot2::scale_x_continuous(name = "Smoothing window (s)")
    out <- list()

    if (!is.null(by_label)) {
        out$odba_by_label <- ggplot2::ggplot(
            by_label,
            ggplot2::aes(winlen_s, median)
        ) +
            ggplot2::geom_ribbon(
                ggplot2::aes(ymin = q25, ymax = q75),
                fill = pal[["blue"]],
                alpha = 0.2
            ) +
            ggplot2::geom_line(colour = pal[["blue"]]) +
            ggplot2::geom_point(colour = pal[["blue"]]) +
            vref +
            ggplot2::facet_wrap(~label, scales = "free_y") +
            xs +
            ggplot2::labs(
                y = paste0("Median ", metric, " (IQR shaded)"),
                title = paste0(
                    metric,
                    " versus smoothing window, by behaviour"
                ),
                subtitle = "Dotted line = reference window"
            ) +
            theme
    }

    sep_long <- overview |>
        dplyr::select(winlen_s, epsilon2, median_pair_sep, min_pair_sep) |>
        tidyr::pivot_longer(
            -winlen_s,
            names_to = "stat",
            values_to = "value"
        ) |>
        dplyr::mutate(
            stat = factor(
                stat,
                levels = c("epsilon2", "median_pair_sep", "min_pair_sep"),
                labels = c(
                    "Kruskal-Wallis epsilon-squared (all labels)",
                    "Median pairwise separability",
                    "Minimum pairwise separability"
                )
            )
        )
    if (any(!is.na(sep_long$value))) {
        out$separation <- ggplot2::ggplot(
            sep_long,
            ggplot2::aes(winlen_s, value)
        ) +
            ggplot2::geom_line(colour = pal[["vermillion"]]) +
            ggplot2::geom_point(colour = pal[["vermillion"]]) +
            vref +
            ggplot2::facet_wrap(~stat, ncol = 1, scales = "free_y") +
            xs +
            ggplot2::labs(
                y = NULL,
                title = paste0("Does ", metric, " still separate behaviours?"),
                subtitle = "Higher = better separation. Pairwise: 0.5 = none, 1 = perfect"
            ) +
            theme
    }

    st_long <- overview |>
        dplyr::select(winlen_s, rank_cor_ref, edge_loss_pct) |>
        tidyr::pivot_longer(
            -winlen_s,
            names_to = "stat",
            values_to = "value"
        ) |>
        dplyr::mutate(
            stat = factor(
                stat,
                levels = c("rank_cor_ref", "edge_loss_pct"),
                labels = c(
                    paste0("Spearman rho with the ", ref_s, " s window"),
                    "Share of a median bout lost to edge NAs (%)"
                )
            )
        )
    out$stability_cost <- ggplot2::ggplot(
        st_long,
        ggplot2::aes(winlen_s, value)
    ) +
        ggplot2::geom_line(colour = pal[["green"]]) +
        ggplot2::geom_point(colour = pal[["green"]]) +
        vref +
        ggplot2::facet_wrap(~stat, ncol = 1, scales = "free_y") +
        xs +
        ggplot2::labs(
            y = NULL,
            title = "Stability and cost of the window",
            subtitle = "Rank agreement of bouts with the reference, and data lost at the bout edges"
        ) +
        theme
    out
}

# screen_features.R ------------------------------------------------------------
# Per-feature screening of a calc_features() output, before feature selection.
# Scores every feature on several complementary criteria, with or without
# behaviour labels, and returns tidy tables you can sort, filter and plot.
#
# Requires: dplyr, tibble, tidyr, rlang, stats, normaliseR, diptest, mclust,
#           RANN; ggplot2 for plot(). (sparcl is not needed: sparse k-means is
#           re-implemented here, see .sparse_kmeans().)
#
# =============================================================================
# CRITERIA
# =============================================================================
# Label-free, K-free (one value per feature)
#   redundancy   Features are clustered on |Spearman correlation| (complete
#                linkage; each column ranked once, so with missing values the
#                ranks use each feature's own non-missing bouts) and cut so every pair inside a group has |rho| >=
#                `cor_threshold` (default 0.9, rabc's cutoff). One
#                representative per group: the member with the highest mean
#                |rho| to the rest. Downstream multivariate criteria
#                (sparse k-means, set-level stability, the Laplacian graph)
#                use representatives only, so redundant copies of the same
#                signal do not get extra votes.
#   dip_stat, dip_p
#                Hartigan's dip test of unimodality on the raw values
#                (diptest). Small p = evidence of more than one mode.
#   gmm_components
#                Number of Gaussian mixture components chosen by BIC (mclust,
#                1-5, unequal variances). More sensitive than the dip test to
#                small modes from rare behaviours, but a skewed unimodal
#                feature can also need 2 components, so this is a count of
#                components, not of modes.
#   within_id_multimodal
#                Share of individuals (with >= `min_bouts_id` bouts) whose own
#                distribution of the feature is multimodal (dip p < 0.05).
#                Separates "modes from behaviour" (present within animals)
#                from "modes from differences between animals or tags".
#                Only informative for individuals that perform several
#                behaviours: an animal recorded doing one thing is unimodal
#                for every feature, however useful the feature is.
#   laplacian    Laplacian score (He, Cai & Niyogi 2005), as in
#                Rdimtools::do.lscore: k-nearest-neighbour graph of bouts
#                (union-symmetrised, heat kernel exp(-d^2 / t)) built on the
#                normalised redundancy representatives; t = mean squared
#                distance to the k neighbours. LOWER = the feature varies more
#                smoothly over the graph, i.e. better preserves local
#                structure. Built with a sparse neighbour search, so it scales
#                to large datasets (Rdimtools builds dense n x n matrices).
#
# Label-free, depends on K (one value per feature per K)
#   sparse_weight
#                Feature weight from sparse k-means (Witten & Tibshirani 2010)
#                on the normalised representatives. 0 = not used. The L1 bound
#                is tuned once by permutation (gap statistic) at the middle K,
#                or set via `sparse_wbound`. Non-representative features get NA
#                (see their group). Same algorithm as
#                sparcl::KMeansSparseCluster / .permute, step for step and with
#                the same random-number use (results match sparcl exactly for
#                a given seed), but the nearest-centre step uses bouts x K
#                distances instead of sparcl's full bouts x bouts distance
#                matrix, which made it ~100x slower on 5,000+ bouts. Like
#                sparcl, the inner k-means updates stop at 10 iterations; the
#                resulting "did not converge" warnings (which sparcl also
#                emits) are silenced.
#   stability    How well a 1-D k-means clustering of this feature transfers
#                between individuals: clusters are fitted on the bouts of
#                some individuals, held-out individuals are assigned to the
#                nearest centre, and that assignment is compared (adjusted
#                Rand index) with clustering the held-out bouts directly.
#                Folds are by individual (`n_folds`). 1 = identical, 0 =
#                chance. Low values flag features whose cut-points differ
#                between animals (e.g. posture offsets from tag placement).
#                CAUTION: high stability does NOT mean the feature has
#                structure. Pure noise with the same distribution in every
#                animal is cut in the same places, so it scores high too.
#                Read it as a filter (low = individual-specific, distrust)
#                alongside multimodality or labelled scores, never alone.
#   Set-level stability (same procedure, all representatives together) is
#   returned per K in $stability, to see which numbers of clusters
#   reproduce across individuals. Read it as a PLATEAU, not a peak: merging
#   real clusters into coarser groups is just as reproducible, so K = 2 is
#   often as stable as the true K (Ben-David, von Luxburg & Pal 2006). In
#   simulations with 6 true clusters the single highest value picked K = 2 in
#   3 of 4 cases; the largest K whose stability stays within `stable_tol`
#   (default 0.05) of the maximum gave 7, 6, 6, 4. That "largest stable K" is
#   what print() reports and what the boundary warning uses.
#
# Labelled (only if `label_col` exists; labels with < `min_n` bouts dropped)
#   fisher       Fisher score: sum_c n_c (mu_c - mu)^2 / sum_c n_c var_c (sample
#                variance), as Rdimtools::do.fscore.
#   kw_epsilon2  Kruskal-Wallis epsilon^2 (rank-based, 0-1).
#   pair_sep_median, pair_sep_min
#                Median / worst separability over behaviour pairs, where
#                separability = max(AUC, 1 - AUC) (0.5 none, 1 perfect).
#                The full feature x pair table is returned in $pairs.
#
# =============================================================================
# K (number of clusters) for the K-dependent criteria
# =============================================================================
#   k = NULL  sweep 2 .. min(12, floor(n_bouts / 50)). 12 is about twice a
#             typical labelled ethogram; the size cap keeps clusters large
#             enough (~50 bouts on average) for stable centres. Labels are
#             never used to choose the range, so labelled and wild data are
#             treated alike.
#   k = 6     a single K.   k = 3:20   your own range.
#   A warning is given if the largest stable K is the top of the range (the
#   range may be too narrow).
#
# =============================================================================
# OUTPUT: a "feature_screen" object (a list of tibbles)
# =============================================================================
#   $summary    one row per feature: identity, redundancy group, every K-free
#               score, and K-dependent scores summarised over K
#               (median / share selected / min).
#   $scores     long: feature, criterion, k (NA if K-free), value.
#   $pairs      feature x behaviour pair separability (labelled data only).
#   $groups     redundancy groups: group, feature, is_representative, rho.
#   $stability  set-level stability per K (mean and sd over folds).
#   $settings   the arguments used.
#   plot(x, type = "stability" | "pairs" | "redundancy")
#
# Normalisation for the clustering/graph criteria: `norm_method` (default
# "RobustSigmoid", via normaliseR) so heavy-tailed features (cv, kurtosis)
# cannot dominate k-means through a few extreme bouts. NAs are replaced by the
# feature's median after normalisation. Dip test, mixture models and labelled
# scores use the raw values (rank-based labelled scores are unaffected).

#' Screen features on several unsupervised (and optional labelled) criteria
#'
#' @param data          Output of calc_features().
#' @param feature_set   Optional set name(s) to screen (via feature_info).
#' @param cols          Optional exact feature columns. Give this or
#'                      `feature_set`; NULL for both = every feature column.
#' @param label_col     Behaviour label column; labelled scores are skipped if
#'                      absent.
#' @param id_col        Individual column, used for folds and within-individual
#'                      multimodality.
#' @param k             NULL (default sweep), a single K, or a vector of K.
#' @param cor_threshold |Spearman rho| above which features count as redundant.
#' @param norm_method   Normalisation for the clustering/graph criteria; one of
#'                      "RobustSigmoid", "zScore", "Sigmoid", "MinMax", "MaxAbs".
#' @param n_folds       Folds (by individual) for the stability criteria.
#' @param min_bouts_id  Minimum bouts for an individual to enter
#'                      within_id_multimodal.
#' @param min_n         Minimum bouts for a label to enter the labelled scores.
#' @param laplacian_k   Neighbours in the Laplacian-score graph.
#' @param sparse_wbound L1 bound for sparse k-means; NULL tunes it by
#'                      permutation at the middle K.
#' @param sparse_nperms Permutations for that tuning.
#' @param feature_stability Compute per-feature stability (the slowest part).
#' @param stable_tol    Tolerance below the maximum set-level stability within
#'                      which a K counts as "stable" (see file header).
#' @param max_bouts     Optional cap on bouts, sampled evenly across
#'                      individuals, for very large datasets.
#' @param seed          Random seed (k-means starts, folds, permutations).
#' @param verbose       Progress messages.
#' @return A "feature_screen" object (see file header).
#' @examples
#' scr <- screen_features(vulture_features)
#' scr$summary |> dplyr::arrange(laplacian)
#' plot(scr, "stability")
screen_features <- function(
    data,
    feature_set = NULL,
    cols = NULL,
    label_col = "label",
    id_col = "id",
    k = NULL,
    cor_threshold = 0.9,
    norm_method = c("RobustSigmoid", "zScore", "Sigmoid", "MinMax", "MaxAbs"),
    n_folds = 5,
    min_bouts_id = 30,
    min_n = 10,
    laplacian_k = 10,
    sparse_wbound = NULL,
    sparse_nperms = 10,
    feature_stability = TRUE,
    stable_tol = 0.05,
    max_bouts = NULL,
    seed = 123,
    verbose = TRUE
) {
    say <- function(...) if (verbose) rlang::inform(paste0(...))
    t0 <- Sys.time()
    norm_method <- match.arg(norm_method)

    # ---- checks and feature columns ----------------------------------------------
    if (!is.data.frame(data)) {
        rlang::abort("`data` must be the tibble returned by calc_features().")
    }
    if (!id_col %in% names(data)) {
        rlang::abort(paste0(
            "No `",
            id_col,
            "` column (individuals); set `id_col`."
        ))
    }
    if (
        !is.numeric(cor_threshold) ||
            length(cor_threshold) != 1 ||
            cor_threshold <= 0 ||
            cor_threshold >= 1
    ) {
        rlang::abort("`cor_threshold` must be a single number between 0 and 1.")
    }
    feat <- .screen_cols(data, feature_set, cols)
    if (length(feat) < 2) {
        rlang::abort("Need at least 2 feature columns to screen.")
    }

    set.seed(seed)
    df <- tibble::as_tibble(data)
    if (!is.null(max_bouts) && nrow(df) > max_bouts) {
        df <- .sample_by_id(df, id_col, max_bouts)
        say(
            "Using ",
            format(nrow(df), big.mark = ","),
            " bouts sampled evenly across individuals (max_bouts)."
        )
    }
    n <- nrow(df)
    ids <- as.character(df[[id_col]])

    X <- as.matrix(df[feat])
    storage.mode(X) <- "double"
    X[!is.finite(X)] <- NA
    n_ok <- colSums(!is.na(X))
    sds <- apply(X, 2, stats::sd, na.rm = TRUE)
    unusable <- feat[n_ok < 10 | is.na(sds) | sds == 0]
    if (length(unusable)) {
        say(
            length(unusable),
            " feature(s) constant or (almost) all NA, not scored: ",
            paste(utils::head(unusable, 5), collapse = ", "),
            if (length(unusable) > 5) ", ..." else ""
        )
    }
    use <- setdiff(feat, unusable)
    if (length(use) < 2) {
        rlang::abort("Fewer than 2 usable (non-constant) features.")
    }
    X <- X[, use, drop = FALSE]

    # normalised matrix for clustering / graph criteria
    Xn <- apply(X, 2, function(v) {
        z <- normaliseR::normalise(as.numeric(v), norm_method = norm_method)
        z[!is.finite(z)] <- NA
        z[is.na(z)] <- stats::median(z, na.rm = TRUE)
        z
    })
    Xn[!is.finite(Xn)] <- 0

    # ---- K range --------------------------------------------------------------------
    if (is.null(k)) {
        k_max <- min(12, floor(n / 50))
        if (k_max < 2) {
            rlang::abort(
                "Too few bouts for the default K range (needs >= 100); give `k` explicitly."
            )
        }
        k <- 2:k_max
    }
    k <- sort(unique(as.integer(k)))
    if (any(is.na(k)) || any(k < 2) || any(k >= n)) {
        rlang::abort(
            "`k` must be whole numbers >= 2 and below the number of bouts."
        )
    }
    say(
        "Screening ",
        length(use),
        " features on ",
        format(n, big.mark = ","),
        " bouts from ",
        length(unique(ids)),
        " individuals; K = ",
        if (length(k) > 1) paste0(min(k), "-", max(k)) else k,
        "."
    )

    scores <- list()
    add <- function(criterion, values, kk = NA_integer_) {
        scores[[length(scores) + 1]] <<- tibble::tibble(
            feature = names(values),
            criterion = criterion,
            k = as.integer(kk),
            value = unname(as.numeric(values))
        )
    }

    # ---- 1. redundancy ----------------------------------------------------------------
    say("Redundancy groups (|Spearman rho| >= ", cor_threshold, ") ...")
    red <- .redundancy_groups(X, cor_threshold)
    reps <- red$feature[red$is_representative]
    say("  ", length(use), " features -> ", length(reps), " groups.")

    # ---- 2. multimodality -------------------------------------------------------------
    say("Multimodality (dip test, mixture models, within individuals) ...")
    dip <- t(vapply(
        use,
        function(f) {
            x <- X[, f]
            x <- x[!is.na(x)]
            d <- suppressWarnings(diptest::dip.test(x))
            c(stat = unname(d$statistic), p = d$p.value)
        },
        numeric(2)
    ))
    add("dip_stat", dip[, "stat"])
    add("dip_p", dip[, "p"])
    add(
        "gmm_components",
        vapply(use, function(f) .gmm_components(X[, f]), numeric(1))
    )
    id_tab <- table(ids)
    big_ids <- names(id_tab)[id_tab >= min_bouts_id]
    if (length(big_ids)) {
        add(
            "within_id_multimodal",
            vapply(
                use,
                function(f) {
                    p <- vapply(
                        big_ids,
                        function(i) {
                            x <- X[ids == i, f]
                            x <- x[!is.na(x)]
                            if (length(x) < min_bouts_id || stats::sd(x) == 0) {
                                return(NA_real_)
                            }
                            suppressWarnings(diptest::dip.test(x)$p.value)
                        },
                        numeric(1)
                    )
                    if (all(is.na(p))) {
                        NA_real_
                    } else {
                        mean(p < 0.05, na.rm = TRUE)
                    }
                },
                numeric(1)
            )
        )
    } else {
        say(
            "  No individual has >= ",
            min_bouts_id,
            " bouts; within_id_multimodal skipped."
        )
    }

    # ---- 3. Laplacian score ----------------------------------------------------------
    say("Laplacian score (", laplacian_k, "-nearest-neighbour graph) ...")
    add("laplacian", .laplacian_score(Xn, graph_cols = reps, k = laplacian_k))

    # ---- 4. sparse k-means -------------------------------------------------------------
    say("Sparse k-means over K ...")
    Xr <- scale(Xn[, reps, drop = FALSE])
    Xr[!is.finite(Xr)] <- 0
    wb <- sparse_wbound
    if (is.null(wb)) {
        k_mid <- k[ceiling(length(k) / 2)]
        perm <- suppressWarnings(.sparse_kmeans_permute(
            Xr,
            K = k_mid,
            nperms = sparse_nperms
        ))
        wb <- perm$bestw
        say(
            "  L1 bound tuned by permutation at K = ",
            k_mid,
            ": ",
            round(wb, 2),
            " (range 1 to ",
            round(sqrt(ncol(Xr)), 2),
            ")."
        )
    }
    for (kk in k) {
        fit <- suppressWarnings(.sparse_kmeans(Xr, K = kk, wbounds = wb))
        add("sparse_weight", stats::setNames(fit[[1]]$ws, reps), kk)
    }

    # ---- 5. stability across individuals ------------------------------------------------
    folds <- .id_folds(ids, n_folds)
    stab_k <- NULL
    if (is.null(folds)) {
        say("Fewer than 2 individuals: stability skipped.")
    } else {
        say(
            "Stability across individuals (",
            length(unique(folds)),
            " folds by individual) ..."
        )
        stab_k <- dplyr::bind_rows(lapply(k, function(kk) {
            a <- vapply(
                sort(unique(folds)),
                function(f) .transfer_ari(Xr, folds == f, kk, nstart = 10),
                numeric(1)
            )
            tibble::tibble(
                k = kk,
                stability = mean(a, na.rm = TRUE),
                sd = stats::sd(a, na.rm = TRUE)
            )
        }))
        stab_k$near_max <- stab_k$stability >=
            max(stab_k$stability, na.rm = TRUE) - stable_tol
        largest_stable <- max(stab_k$k[stab_k$near_max %in% TRUE])
        if (length(k) > 2 && largest_stable == max(k)) {
            rlang::warn(paste0(
                "Stability is still near its maximum at the largest K tested (",
                max(k),
                "); the range may be too narrow. Try a wider `k`."
            ))
        }
        if (isTRUE(feature_stability)) {
            for (kk in k) {
                add(
                    "stability",
                    vapply(
                        use,
                        function(f) {
                            x <- Xn[, f, drop = FALSE]
                            if (length(unique(x)) <= kk) {
                                return(NA_real_)
                            }
                            mean(
                                vapply(
                                    sort(unique(folds)),
                                    function(fo) {
                                        .transfer_ari(
                                            x,
                                            folds == fo,
                                            kk,
                                            nstart = 3
                                        )
                                    },
                                    numeric(1)
                                ),
                                na.rm = TRUE
                            )
                        },
                        numeric(1)
                    ),
                    kk
                )
            }
        }
    }

    # ---- 6. labelled scores --------------------------------------------------------------
    pairs <- NULL
    has_label <- !is.null(label_col) && label_col %in% names(df)
    n_labels <- NA_integer_
    if (has_label) {
        lab <- as.character(df[[label_col]])
        cnt <- table(lab[!is.na(lab)])
        keep_lab <- sort(names(cnt)[cnt >= min_n])
        if (length(keep_lab) >= 2) {
            n_labels <- length(keep_lab)
            say("Labelled scores on ", n_labels, " behaviours ...")
            in_lab <- !is.na(lab) & lab %in% keep_lab
            L <- lab[in_lab]
            lab_sc <- t(vapply(
                use,
                function(f) .label_scores(X[in_lab, f], L),
                numeric(2)
            ))
            add("fisher", lab_sc[, "fisher"])
            add("kw_epsilon2", lab_sc[, "kw_epsilon2"])
            prs <- utils::combn(keep_lab, 2, simplify = FALSE)
            pairs <- dplyr::bind_rows(lapply(use, function(f) {
                x <- X[in_lab, f]
                tibble::tibble(
                    feature = f,
                    label_a = vapply(prs, `[`, character(1), 1),
                    label_b = vapply(prs, `[`, character(1), 2),
                    separability = vapply(
                        prs,
                        function(p) .separability(x[L == p[1]], x[L == p[2]]),
                        numeric(1)
                    )
                )
            }))
            ps <- pairs |>
                dplyr::summarise(
                    med = stats::median(.data$separability, na.rm = TRUE),
                    min = min(.data$separability, na.rm = TRUE),
                    .by = "feature"
                )
            add("pair_sep_median", stats::setNames(ps$med, ps$feature))
            add("pair_sep_min", stats::setNames(ps$min, ps$feature))
        } else {
            say(
                "Fewer than 2 behaviours with >= ",
                min_n,
                " bouts; labelled scores skipped."
            )
        }
    }

    # ---- assemble -------------------------------------------------------------------------
    scores <- dplyr::bind_rows(scores)
    out <- structure(
        list(
            summary = .screen_summary(scores, red, feat, unusable, data),
            scores = scores,
            pairs = pairs,
            groups = red,
            stability = stab_k,
            settings = list(
                k = k,
                cor_threshold = cor_threshold,
                norm_method = norm_method,
                stable_tol = stable_tol,
                n_folds = n_folds,
                laplacian_k = laplacian_k,
                sparse_wbound = wb,
                n_bouts = n,
                n_individuals = length(unique(ids)),
                n_labels = n_labels,
                label_col = if (has_label) label_col else NA_character_,
                seed = seed
            )
        ),
        class = "feature_screen"
    )
    say(
        "Done in ",
        round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1),
        " s."
    )
    out
}

# =============================================================================
# Criteria helpers
# =============================================================================

# Complete-linkage clustering on 1 - |Spearman rho|; every pair within a group
# has |rho| >= threshold. Representative = highest mean |rho| to the others.
.redundancy_groups <- function(X, threshold) {
    R <- apply(X, 2, rank, na.last = "keep", ties.method = "average") # rank once, then Pearson
    r <- suppressWarnings(stats::cor(R, use = "pairwise.complete.obs"))
    r[is.na(r)] <- 0
    a <- abs(r)
    grp <- if (ncol(X) > 1) {
        stats::cutree(
            stats::hclust(stats::as.dist(1 - a), method = "complete"),
            h = 1 - threshold
        )
    } else {
        stats::setNames(1L, colnames(X))
    }
    # renumber groups in column order
    grp <- match(grp, unique(grp))
    res <- tibble::tibble(feature = colnames(X), group = grp)
    res$group_size <- as.integer(table(grp)[as.character(grp)])
    res$mean_abs_rho <- vapply(
        seq_len(ncol(X)),
        function(j) {
            m <- which(grp == grp[j] & seq_len(ncol(X)) != j)
            if (length(m)) mean(a[j, m]) else NA_real_
        },
        numeric(1)
    )
    res$is_representative <- FALSE
    for (g in unique(grp)) {
        m <- which(grp == g)
        best <- if (length(m) == 1) m else m[which.max(res$mean_abs_rho[m])]
        res$is_representative[best] <- TRUE
    }
    res$representative <- res$feature[res$is_representative][match(
        res$group,
        res$group[res$is_representative]
    )]
    res
}

.gmm_components <- function(x, max_n = 3000) {
    x <- x[!is.na(x)]
    if (length(x) > max_n) {
        x <- sample(x, max_n)
    }
    if (length(unique(x)) < 5) {
        return(NA_real_)
    }
    # mclustBIC directly: Mclust() only works with mclust attached
    bic <- tryCatch(
        .quiet(mclust::mclustBIC(
            x,
            G = 1:5,
            modelNames = "V",
            verbose = FALSE
        )),
        error = function(e) NULL
    )
    if (is.null(bic) || all(is.na(bic))) {
        return(NA_real_)
    }
    as.numeric(which.max(bic[, 1])) # mclust's BIC: larger is better
}

# Laplacian score (He et al. 2005), as Rdimtools::do.lscore but with a sparse
# k-NN graph: S_ij = exp(-d_ij^2 / t) for union-symmetrised neighbours.
.laplacian_score <- function(Xn, graph_cols, k = 10, t = NULL) {
    G <- Xn[, graph_cols, drop = FALSE]
    n <- nrow(G)
    kk <- min(k + 1, n)
    nn <- RANN::nn2(G, k = kk)
    i <- rep(seq_len(n), kk - 1)
    j <- as.vector(nn$nn.idx[, -1, drop = FALSE])
    d2 <- as.vector(nn$nn.dists[, -1, drop = FALSE])^2
    if (is.null(t)) {
        t <- mean(d2)
    }
    if (!is.finite(t) || t <= 0) {
        t <- 1
    }
    # union symmetrisation: keep each unordered pair once
    key <- paste(pmin(i, j), pmax(i, j))
    keep <- !duplicated(key) & i != j
    i <- i[keep]
    j <- j[keep]
    w <- exp(-d2[keep] / t)
    deg <- tabulate(i, n) * 0
    deg <- deg +
        vapply(
            split(c(w, w), factor(c(i, j), levels = seq_len(n))),
            sum,
            numeric(1)
        )
    vapply(
        colnames(Xn),
        function(f) {
            fr <- Xn[, f]
            ft <- fr - sum(fr * deg) / sum(deg)
            num <- sum(w * (ft[i] - ft[j])^2) # f' L f for the symmetric graph
            den <- sum(deg * ft^2)
            if (den <= 0) NA_real_ else num / den
        },
        numeric(1)
    )
}

# ---- sparse k-means (Witten & Tibshirani 2010), following sparcl 1.0.4 -------------
# Same steps, defaults and random-number use as sparcl::KMeansSparseCluster and
# KMeansSparseCluster.permute; only the nearest-centre computation differs
# (bouts x K instead of a full bouts x bouts distance matrix).
.skm_wcss <- function(x, Cs) {
    wcss <- numeric(ncol(x))
    for (k in unique(Cs)) {
        w <- Cs == k
        if (sum(w) > 1) {
            wcss <- wcss +
                colSums(
                    scale(x[w, , drop = FALSE], center = TRUE, scale = FALSE)^2
                )
        }
    }
    bcss <- colSums(scale(x, center = TRUE, scale = FALSE)^2) - wcss
    list(wcss.perfeature = wcss, bcss.perfeature = bcss)
}
.skm_soft <- function(x, d) sign(x) * pmax(0, abs(x) - d)
.skm_l2n <- function(v) sqrt(sum(v^2))
.skm_binary_search <- function(argu, sumabs) {
    if (.skm_l2n(argu) == 0 || sum(abs(argu / .skm_l2n(argu))) <= sumabs) {
        return(0)
    }
    lam1 <- 0
    lam2 <- max(abs(argu)) - 1e-5
    iter <- 1
    while (iter <= 15 && (lam2 - lam1) > 1e-4) {
        su <- .skm_soft(argu, (lam1 + lam2) / 2)
        if (sum(abs(su / .skm_l2n(su))) < sumabs) {
            lam2 <- (lam1 + lam2) / 2
        } else {
            lam1 <- (lam1 + lam2) / 2
        }
        iter <- iter + 1
    }
    (lam1 + lam2) / 2
}
.skm_update_ws <- function(x, Cs, l1bound) {
    wcss <- .skm_wcss(x, Cs)$wcss.perfeature
    tss <- .skm_wcss(x, rep(1, nrow(x)))$wcss.perfeature
    lam <- .skm_binary_search(-wcss + tss, l1bound)
    ws <- .skm_soft(-wcss + tss, lam)
    ws / .skm_l2n(ws)
}
.skm_update_cs <- function(x, K, ws, Cs) {
    x <- x[, ws != 0, drop = FALSE]
    z <- sweep(x, 2, sqrt(ws[ws != 0]), "*")
    mus <- NULL
    if (!is.null(Cs)) {
        for (k in unique(Cs)) {
            if (sum(Cs == k) > 1) {
                mus <- rbind(mus, colMeans(z[Cs == k, , drop = FALSE]))
            }
            if (sum(Cs == k) == 1) mus <- rbind(mus, z[Cs == k, ])
        }
    }
    if (is.null(mus)) {
        return(stats::kmeans(z, centers = K, nstart = 10)$cluster)
    }
    d2 <- vapply(
        seq_len(nrow(mus)),
        function(c) rowSums(sweep(z, 2, mus[c, ])^2),
        numeric(nrow(z))
    )
    nearest <- if (is.matrix(d2)) max.col(-d2, ties.method = "first") else 1L
    km <- if (length(unique(nearest)) == K) {
        stats::kmeans(z, centers = mus)
    } else {
        stats::kmeans(z, centers = K, nstart = 10)
    }
    km$cluster
}
.sparse_kmeans <- function(x, K, wbounds = NULL, nstart = 20, maxiter = 6) {
    if (is.null(wbounds)) {
        wbounds <- seq(1.1, sqrt(ncol(x)), len = 20)
    }
    if (min(wbounds) <= 1) {
        rlang::abort("Sparse k-means: wbounds must be > 1.")
    }
    Cs <- stats::kmeans(x, centers = K, nstart = nstart)$cluster
    out <- list()
    for (i in seq_along(wbounds)) {
        ws <- rep(1 / sqrt(ncol(x)), ncol(x))
        ws.old <- stats::rnorm(ncol(x))
        niter <- 0
        while (
            (sum(abs(ws - ws.old)) / sum(abs(ws.old))) > 1e-4 && niter < maxiter
        ) {
            niter <- niter + 1
            ws.old <- ws
            if (niter > 1) {
                Cs <- .skm_update_cs(x, K, ws, Cs)
            }
            ws <- .skm_update_ws(x, Cs, wbounds[i])
        }
        out[[i]] <- list(ws = ws, Cs = Cs, wbound = wbounds[i])
    }
    out
}
.sparse_kmeans_permute <- function(
    x,
    K,
    nperms = 25,
    wbounds = NULL,
    nvals = 10
) {
    if (is.null(wbounds)) {
        wbounds <- exp(seq(log(1.2), log(sqrt(ncol(x)) * 0.9), len = nvals))
    }
    permx <- lapply(seq_len(nperms), function(i) apply(x, 2, sample))
    out <- .sparse_kmeans(x, K, wbounds)
    tots <- vapply(
        out,
        function(o) sum(o$ws * .skm_wcss(x, o$Cs)$bcss.perfeature),
        numeric(1)
    )
    permtots <- vapply(
        permx,
        function(px) {
            po <- .sparse_kmeans(px, K, wbounds)
            vapply(
                po,
                function(o) sum(o$ws * .skm_wcss(px, o$Cs)$bcss.perfeature),
                numeric(1)
            )
        },
        numeric(length(wbounds))
    )
    permtots <- matrix(permtots, nrow = length(wbounds))
    gaps <- log(tots) - rowMeans(log(permtots))
    list(
        tots = tots,
        permtots = permtots,
        gaps = gaps,
        wbounds = wbounds,
        nnonzerows = vapply(out, function(o) sum(o$ws != 0), numeric(1)),
        bestw = wbounds[which.max(gaps)]
    )
}

# Assign individuals to folds with roughly equal numbers of bouts.
.id_folds <- function(ids, n_folds) {
    u <- unique(ids)
    if (length(u) < 2) {
        return(NULL)
    }
    nf <- min(n_folds, length(u))
    sz <- table(ids)[u]
    ord <- u[order(-as.numeric(sz), stats::runif(length(u)))]
    load <- numeric(nf)
    f_of <- stats::setNames(integer(length(u)), ord)
    for (x in ord) {
        b <- which.min(load)
        f_of[x] <- b
        load[b] <- load[b] + sz[[x]]
    }
    unname(f_of[ids])
}

# Cluster the training rows, assign the test rows to the nearest centre, and
# compare with clustering the test rows directly (adjusted Rand index).
.transfer_ari <- function(X, test, k, nstart = 10) {
    tr <- X[!test, , drop = FALSE]
    te <- X[test, , drop = FALSE]
    if (nrow(te) <= k || nrow(tr) <= k) {
        return(NA_real_)
    }
    km_tr <- tryCatch(
        suppressWarnings(stats::kmeans(tr, k, nstart = nstart, iter.max = 50)),
        error = function(e) NULL
    )
    km_te <- tryCatch(
        suppressWarnings(stats::kmeans(te, k, nstart = nstart, iter.max = 50)),
        error = function(e) NULL
    )
    if (is.null(km_tr) || is.null(km_te)) {
        return(NA_real_)
    }
    cen <- km_tr$centers
    d <- vapply(
        seq_len(nrow(cen)),
        function(c) rowSums(sweep(te, 2, cen[c, ])^2),
        numeric(nrow(te))
    )
    pred <- if (is.matrix(d)) {
        max.col(-d, ties.method = "first")
    } else {
        which.min(d)
    }
    .ari(pred, km_te$cluster)
}

.ari <- function(a, b) {
    tab <- table(a, b)
    ch2 <- function(x) sum(x * (x - 1) / 2)
    s_ij <- ch2(tab)
    s_a <- ch2(rowSums(tab))
    s_b <- ch2(colSums(tab))
    nn <- ch2(sum(tab))
    expct <- s_a * s_b / nn
    mx <- (s_a + s_b) / 2
    if (mx == expct) {
        return(NA_real_)
    }
    (s_ij - expct) / (mx - expct)
}

.label_scores <- function(x, L) {
    ok <- !is.na(x)
    x <- x[ok]
    L <- L[ok]
    if (length(unique(L)) < 2) {
        return(c(fisher = NA_real_, kw_epsilon2 = NA_real_))
    }
    mu <- mean(x)
    grp <- split(x, L)
    num <- sum(vapply(
        grp,
        function(g) length(g) * (mean(g) - mu)^2,
        numeric(1)
    ))
    den <- sum(vapply(
        grp,
        function(g) if (length(g) > 1) length(g) * stats::var(g) else 0,
        numeric(1)
    )) # as Rdimtools::do.fscore
    H <- suppressWarnings(unname(stats::kruskal.test(x, factor(L))$statistic))
    nn <- length(x)
    c(
        fisher = if (den > 0) num / den else NA_real_,
        kw_epsilon2 = H * (nn + 1) / (nn^2 - 1)
    )
}

# max(AUC, 1 - AUC) via the rank-sum identity (same as winlen_sensitivity()).
.separability <- function(a, b) {
    a <- a[!is.na(a)]
    b <- b[!is.na(b)]
    if (!length(a) || !length(b)) {
        return(NA_real_)
    }
    r <- rank(c(a, b))
    n1 <- length(a)
    n2 <- length(b)
    auc <- (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * n2)
    max(auc, 1 - auc)
}

# =============================================================================
# Plumbing
# =============================================================================

# Feature columns: exactly two "__" in the name; `feature_set` via the
# feature_info attribute (fallback: first part of the name).
.screen_cols <- function(data, feature_set, cols) {
    nms <- names(data)
    feat <- nms[
        lengths(regmatches(nms, gregexpr("__", nms, fixed = TRUE))) == 2
    ]
    if (!length(feat)) {
        rlang::abort(
            "No `<family>__<stream>__<feature>` columns: is `data` the output of calc_features()?"
        )
    }
    if (!is.null(feature_set) && !is.null(cols)) {
        rlang::abort("Give either `feature_set` or `cols`, not both.")
    }
    if (!is.null(cols)) {
        bad <- setdiff(cols, feat)
        if (length(bad)) {
            rlang::abort(paste0(
                "`cols` not among feature columns: ",
                paste(bad, collapse = ", ")
            ))
        }
        return(cols)
    }
    if (!is.null(feature_set)) {
        info <- attr(data, "feature_info", exact = TRUE)
        member <- if (!is.null(info)) {
            unique(info[info$column %in% feat, c("column", "feature_set")])
        } else {
            tibble::tibble(column = feat, feature_set = sub("__.*", "", feat))
        }
        bad <- setdiff(feature_set, member$feature_set)
        if (length(bad)) {
            rlang::abort(c(
                paste0(
                    "`feature_set` not in `data`: ",
                    paste(bad, collapse = ", ")
                ),
                i = paste0(
                    "Available: ",
                    paste(unique(member$feature_set), collapse = ", ")
                )
            ))
        }
        keep <- unique(member$column[member$feature_set %in% feature_set])
        return(feat[feat %in% keep])
    }
    feat
}

.sample_by_id <- function(df, id_col, max_bouts) {
    ids <- as.character(df[[id_col]])
    per <- split(seq_len(nrow(df)), ids)
    cap <- max(1, floor(max_bouts / length(per)))
    rows <- unlist(lapply(per, function(r) {
        if (length(r) > cap) sample(r, cap) else r
    }))
    df[sort(rows), ]
}

.quiet <- function(expr) {
    out <- NULL
    utils::capture.output(out <- suppressWarnings(suppressMessages(expr)))
    out
}

.screen_summary <- function(scores, red, feat, unusable, data) {
    parts <- strsplit(feat, "__", fixed = TRUE)
    base <- tibble::tibble(
        feature = feat,
        family = vapply(parts, `[`, character(1), 1),
        stream = vapply(parts, `[`, character(1), 2),
        name = vapply(parts, `[`, character(1), 3)
    )
    info <- attr(data, "feature_info", exact = TRUE)
    if (!is.null(info)) {
        sets <- info |>
            dplyr::summarise(
                feature_sets = paste(
                    unique(.data$feature_set),
                    collapse = ", "
                ),
                .by = "column"
            )
        base$feature_sets <- sets$feature_sets[match(base$feature, sets$column)]
    }
    base$usable <- !base$feature %in% unusable
    base <- dplyr::left_join(
        base,
        red[c(
            "feature",
            "group",
            "group_size",
            "is_representative",
            "representative"
        )],
        by = "feature"
    )

    kfree <- scores[is.na(scores$k), ]
    wide <- tidyr::pivot_wider(
        kfree[c("feature", "criterion", "value")],
        names_from = "criterion",
        values_from = "value"
    )
    base <- dplyr::left_join(base, wide, by = "feature")

    kdep <- scores[!is.na(scores$k), ]
    if (nrow(kdep)) {
        kd <- kdep |>
            dplyr::summarise(
                sparse_weight_median = if (
                    any(.data$criterion == "sparse_weight")
                ) {
                    stats::median(
                        .data$value[.data$criterion == "sparse_weight"],
                        na.rm = TRUE
                    )
                } else {
                    NA_real_
                },
                sparse_selected_prop = if (
                    any(.data$criterion == "sparse_weight")
                ) {
                    mean(
                        .data$value[.data$criterion == "sparse_weight"] > 0,
                        na.rm = TRUE
                    )
                } else {
                    NA_real_
                },
                stability_median = if (any(.data$criterion == "stability")) {
                    stats::median(
                        .data$value[.data$criterion == "stability"],
                        na.rm = TRUE
                    )
                } else {
                    NA_real_
                },
                stability_min = if (any(.data$criterion == "stability")) {
                    suppressWarnings(min(
                        .data$value[.data$criterion == "stability"],
                        na.rm = TRUE
                    ))
                } else {
                    NA_real_
                },
                .by = "feature"
            )
        kd$stability_min[!is.finite(kd$stability_min)] <- NA
        base <- dplyr::left_join(base, kd, by = "feature")
    }
    base
}

#' @export
print.feature_screen <- function(x, ...) {
    s <- x$settings
    cat(
        "<feature_screen>",
        nrow(x$summary),
        "features,",
        format(s$n_bouts, big.mark = ","),
        "bouts,",
        s$n_individuals,
        "individuals\n"
    )
    cat(
        "  Redundancy groups (|rho| >=",
        s$cor_threshold,
        "):",
        length(unique(stats::na.omit(x$groups$group))),
        "\n"
    )
    cat(
        "  K tested:",
        paste(range(s$k), collapse = "-"),
        " | labelled scores:",
        if (is.na(s$n_labels)) "no" else paste(s$n_labels, "behaviours"),
        "\n"
    )
    if (!is.null(x$stability)) {
        near <- x$stability$k[x$stability$near_max %in% TRUE]
        cat(
            "  Stable K (within",
            s$stable_tol,
            "of max stability):",
            paste(near, collapse = ", "),
            "-> largest stable K =",
            max(near),
            "\n"
        )
    }
    cat(
        "  Tables: $summary, $scores, $pairs, $groups, $stability. plot(x, 'stability' | 'pairs' | 'redundancy')\n"
    )
    invisible(x)
}

#' Plot a feature_screen
#'
#' @param x    A feature_screen.
#' @param type "stability" (set-level stability against K), "pairs" (feature x
#'             behaviour-pair separability heatmap, labelled data) or
#'             "redundancy" (group sizes).
#' @param top  For "pairs": number of features shown (highest median
#'             separability).
#' @param ...  Unused.
plot.feature_screen <- function(
    x,
    type = c("stability", "pairs", "redundancy"),
    top = 25,
    ...
) {
    type <- match.arg(type)
    th <- ggplot2::theme_bw() +
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
    if (type == "stability") {
        if (is.null(x$stability)) {
            rlang::abort("No stability results (fewer than 2 individuals).")
        }
        top <- max(x$stability$stability, na.rm = TRUE)
        p <- ggplot2::ggplot(
            x$stability,
            ggplot2::aes(.data$k, .data$stability)
        ) +
            ggplot2::annotate(
                "rect",
                xmin = -Inf,
                xmax = Inf,
                ymin = top - x$settings$stable_tol,
                ymax = top,
                fill = "#0072B2",
                alpha = 0.08
            ) +
            ggplot2::geom_errorbar(
                ggplot2::aes(
                    ymin = .data$stability - .data$sd,
                    ymax = .data$stability + .data$sd
                ),
                width = 0.15,
                colour = "grey60"
            ) +
            ggplot2::geom_line(colour = "#0072B2") +
            ggplot2::geom_point(colour = "#0072B2", size = 2) +
            ggplot2::scale_x_continuous(breaks = x$settings$k) +
            ggplot2::labs(
                x = "Number of clusters (K)",
                y = "Stability across individuals (ARI)",
                title = "Which numbers of clusters reproduce across individuals?",
                subtitle = paste0(
                    "Mean +/- sd over individual folds. Shaded band: within ",
                    x$settings$stable_tol,
                    " of the maximum.\nRead the largest K inside the band, not the peak."
                )
            ) +
            th
        if (!is.na(x$settings$n_labels)) {
            p <- p +
                ggplot2::geom_vline(
                    xintercept = x$settings$n_labels,
                    linetype = "dashed"
                ) +
                ggplot2::annotate(
                    "text",
                    x = x$settings$n_labels,
                    y = max(x$stability$stability, na.rm = TRUE),
                    label = " labelled behaviours",
                    hjust = 0,
                    vjust = 1,
                    size = 3.2
                )
        }
        return(p)
    }
    if (type == "pairs") {
        if (is.null(x$pairs)) {
            rlang::abort(
                "No labelled results: `label_col` absent or too few labelled bouts."
            )
        }
        keep <- x$summary |>
            dplyr::filter(!is.na(.data$pair_sep_median)) |>
            dplyr::arrange(dplyr::desc(.data$pair_sep_median)) |>
            utils::head(top)
        d <- x$pairs |>
            dplyr::filter(.data$feature %in% keep$feature) |>
            dplyr::mutate(
                pair = paste(.data$label_a, "vs", .data$label_b),
                feature = factor(.data$feature, levels = rev(keep$feature))
            )
        return(
            ggplot2::ggplot(
                d,
                ggplot2::aes(
                    .data$pair,
                    .data$feature,
                    fill = .data$separability
                )
            ) +
                ggplot2::geom_tile(colour = "white") +
                ggplot2::scale_fill_viridis_c(
                    limits = c(0.5, 1),
                    name = "Separability"
                ) +
                ggplot2::labs(
                    x = NULL,
                    y = NULL,
                    title = "Which features separate which behaviours?",
                    subtitle = paste0(
                        "Top ",
                        nrow(keep),
                        " features by median separability (0.5 = none, 1 = perfect)"
                    )
                ) +
                th +
                ggplot2::theme(
                    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
                )
        )
    }
    g <- x$groups |>
        dplyr::filter(.data$group_size > 1) |>
        dplyr::distinct(.data$group, .data$group_size, .data$representative)
    if (!nrow(g)) {
        rlang::abort("No redundant features at this threshold.")
    }
    ggplot2::ggplot(
        g,
        ggplot2::aes(
            stats::reorder(.data$representative, .data$group_size),
            .data$group_size
        )
    ) +
        ggplot2::geom_col(fill = "grey40") +
        ggplot2::coord_flip() +
        ggplot2::labs(
            x = "Group representative",
            y = "Features in group",
            title = paste0(
                "Redundant feature groups (|Spearman rho| >= ",
                x$settings$cor_threshold,
                ")"
            )
        ) +
        th
}
