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
    run_break <- run_break | c(FALSE, diff(dplyr::consecutive_id(label)) > 0)
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
      is.null(fs) || !is.numeric(fs) || length(fs) != 1 || is.na(fs) || fs <= 0
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

# now to make an interactive dygraph to visualise data
# inspired by rabc's plot_acc
# How to have multiple visualisation options all rolled into one function
# 1) visualise a specified individual's labels/behaviours - each marked by a dyEvent
# if the data is split into bouts and lacks a timestamp, the timestamp must be sample_in_bout
# and separate bouts can be stiched to follow one another
# arranged by behaviour? or left to the user to specify based on if they want to see it by timestamp or not
# if the data has a timestamp, can use that as the x_axis
# 2) visualise a specified behaviour of multiple individual's - dyEvent marking id instead of label

# and ability to specify which time series axes to plot
# (must be named - e.g. acc_x, acc_y, etc.
#  Our users might have other data streams like magnetometers)

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

# function
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
    rlang::abort(paste0("No id column called `", id_col, "`. Set `id_col`."))
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
      if (length(block_levels) == 1) block_levels else "arranged by behaviour"
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
      d <- d[as.character(d[[id_col]]) %in% as.character(id), , drop = FALSE]
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


# calc features 2
# Bout-level feature calculation for a tidy sensor tsibble
# (key = c(id, bout), index = sample_in_bout or datetime, one column per stream).
# Modelled on theft::calculate_features(), but multivariate: any sensor columns
# can be chosen, and feature sets that combine axes (ODBA, cor_x_y, ...) are
# supported.
#
# Requires: tsibble, tibble, dplyr, tidyr, tidyselect, rlang, vctrs, zoo,
#           entropy, Rcatch22
#
# Named calc_features(), not calculate_features(), so it can be loaded
# alongside {theft} without masking theft::calculate_features().
#
# Feature sets
#   catch22     Rcatch22 features on EVERY selected stream (any numeric
#               column). `catch24 = TRUE` adds DN_Mean and DN_Spread_Std.
#   rabc_time1  Original rabc time-domain set: per-axis mean, variance, sd,
#               max, min, range, plus ODBA.
#   rabc_time2  Expanded set = rabc_time1 + per-axis norm, cross-axis cov / cor /
#               meandiff / sddiff, and per-axis varsba / vardba / maxdba.
#               It CONTAINS rabc_time1, so if both are requested only
#               rabc_time2 is computed (with a message).
#   rabc_freq   Per-axis dominant frequency (freqmain), its amplitude
#               (freqamp) and a spectral entropy term (entropy).
#
# Accelerometer axes for the rabc sets
#   The rabc sets read the columns named in `acc_cols` (default acc_x, acc_y,
#   acc_z) from the selected `cols`. Missing axes do not cause an error:
#     * x, y and z present -> ODBA
#     * two axes present   -> PDBA_xy / PDBA_xz / PDBA_yz in place of ODBA, and
#                             only the cross-axis features for that pair
#     * one axis present   -> PDBA_x / PDBA_y / PDBA_z, no cross-axis features
#   An error is raised only if NONE of the axes is selected.
#
# Output (a plain tibble, one row per bout)
#   key columns (e.g. id, bout), then context columns (e.g. attachment, label;
#   any non-selected column that is constant within every bout), then feature
#   columns named   <feature_set>__<stream>__<feature>
#   e.g. catch22__acc_x__DN_HistogramMode_5, rabc_time2__acc_x_acc_y__cor,
#        rabc_time1__acc_x_acc_y_acc_z__ODBA, rabc_time1__acc_x_acc_z__PDBA_xz
#   custom      (only if `features` is given) EVERY function in `features` on
#               EVERY selected stream: custom__<stream>__<name>.
#   Use feature_long() to pivot into a long theft-style layout.
#
# Dynamic-acceleration window (rabc_time sets, and any cross_feature(window = TRUE))
#   ODBA/PDBA, varsba, vardba and maxdba split each axis into a static part (a
#   centred moving average) and a dynamic part (raw minus static). Set the
#   window with EITHER
#     winlen_dba_s  seconds (needs `fs`); the recommended, dataset-independent
#                   way. Converted to an ODD number of samples:
#                   n = 2 * floor(winlen_dba_s * fs / 2) + 1, so the window is
#                   centred on a sample. e.g. 1 s at 20 Hz -> 21 samples.
#     winlen_dba    samples; use to reproduce earlier rabc results exactly.
#   Giving both is an error; there is no default. The resolved window and the
#   share of a median-length bout lost to edge NAs are always reported.
#
# Missing values: a bout with any NA / non-finite value in a stream a set
#   needs gets NA for that set's features (per stream for catch22; for the whole
#   set for rabc sets, since features combine axes). Features that are
#   undefined (e.g. a correlation on a constant signal) are also NA. A message
#   reports how many bouts are affected.
#
# Custom features (`features`, as in theft::calculate_features())
#   A named list of entries, each producing ONE feature column, computed per
#   bout. Two kinds of entry:
#     * a plain function -> applied to EVERY selected stream separately
#       (like catch22): `list(mean = mean, sd = sd)` gives one `mean`/`sd`
#       column per stream. Column name: custom__<stream>__<name>.
#     * cross_feature(cols, fn, window = FALSE) -> a feature computed from
#       SEVERAL streams together, e.g. a hand-written ODBA:
#         cross_feature(c("acc_x", "acc_y", "acc_z"), function(x, y, z) {...})
#       `cols` must already be among the selected `cols`. Arguments are POSITIONAL,
#       matched to `cols` in the order given there (names x/y/z above are just
#       the parameter names chosen for readability; acc_x is always first
#       because it's listed first) -- NOT by matching parameter names to column
#       names, and not by each column's position in `cols`/the data. One output
#       column, named after every column it uses in that order:
#       custom__acc_x_acc_y_acc_z__<name>.
#       `window = TRUE` adds ONE more positional argument after the columns:
#       the resolved smoothing window, in SAMPLES -- the same value rabc_time1/2
#       use for ODBA/PDBA (see "Dynamic-acceleration window" below). This lets a
#       custom feature reuse the exact same static/dynamic split, e.g. a
#       hand-written ODBA that matches rabc_time1's:
#         cross_feature(c("acc_x", "acc_y", "acc_z"),
#           function(x, y, z, winlen) {
#             sba <- function(v) roll_static(v, winlen)   # see helper below
#             mean(abs(x - sba(x)) + abs(y - sba(y)) + abs(z - sba(z)), na.rm = TRUE)
#           }, window = TRUE)
#       Requesting this makes `winlen_dba`/`winlen_dba_s` required even if no
#       rabc_time set is requested, and the SAME window is then used for both
#       rabc_time1/2 and every window-using custom feature -- there is no way
#       to give a custom feature a different window from rabc_time in one call.
#       roll_static(v, winlen) is exported: the exact centred rolling mean
#       rabc_time1/2 use internally, for building a matching custom feature.
#   Every function (plain or inside cross_feature) must return a single number.
#   One that errors, or does not, yields NA for that stream/bout (folded into
#   the same NA count as missing input; there is no separate error report).
#   Plain-function columns come first (grouped per stream), then cross_feature
#   columns, in the order given in `features`.
#   `feature_set` can be empty (`feature_set = NULL`) if only custom features
#   are wanted; at least one of `feature_set` / `features` is required.
#
# `seed`: every set implemented so far is deterministic. The argument exists
#   for parity with theft::calculate_features() and future stochastic sets.

.feature_sets <- c("catch22", "rabc_time1", "rabc_time2", "rabc_freq")

#' Calculate bout-level features from a sensor tsibble
#'
#' @param data        A tsibble. Every key combination (e.g. id + bout) is one
#'                    bout and yields one row of features.
#' @param cols        Sensor streams to use, as tidyselect: `acc_x`,
#'                    `c(acc_x, acc_y, acc_z)`, `starts_with("acc_")`, ...
#' @param feature_set One or more of "catch22", "rabc_time1", "rabc_time2",
#'                    "rabc_freq". Can be `NULL`/empty if `features` is given.
#' @param catch24     Also compute catch24 (adds DN_Mean, DN_Spread_Std) when
#'                    catch22 is requested.
#' @param winlen_dba_s Rolling-mean window in SECONDS used to separate static from
#'                    dynamic acceleration (rabc_time1/2). Needs `fs`. Rounded to
#'                    an odd number of samples. Give this or `winlen_dba`.
#' @param winlen_dba  The same window in SAMPLES. Give this or `winlen_dba_s`.
#' @param fs          Sampling rate in Hz. Required for rabc_freq and for
#'                    `winlen_dba_s`.
#' @param features    Optional named list of custom features (as in
#'                    theft::calculate_features()). A plain function, e.g.
#'                    `list(mean = mean, sd = sd)`, is applied to every selected
#'                    stream separately. `cross_feature()` computes one feature
#'                    from several streams together; see its own help and the
#'                    file header for details, column naming, and the `window`
#'                    option (a hand-written ODBA using the same window as
#'                    rabc_time1/2).
#' @param acc_cols    Named character vector mapping axes to column names.
#'                    Change it if your axes are not called acc_x/acc_y/acc_z.
#' @param seed        Integer for set.seed(), or NULL to leave the RNG alone.
#' @param verbose     Print progress messages.
#'
#' @return A tibble with one row per bout (see file header for the layout).
#' @examples
#' calc_features(vultures_tsbl, c(acc_x, acc_y, acc_z),
#'               feature_set = c("catch22", "rabc_time2", "rabc_freq"),
#'               winlen_dba_s = 1, fs = 20)
#' calc_features(vultures_tsbl, acc_x, feature_set = "catch22",
#'               features = list(mean = mean, sd = sd))
#' calc_features(vultures_tsbl, c(acc_x, acc_y, acc_z), feature_set = NULL,
#'               features = list(my_odba = cross_feature(
#'                 c("acc_x", "acc_y", "acc_z"),
#'                 function(x, y, z) mean(abs(x) + abs(y) + abs(z))
#'               )))
#' # a hand-written ODBA that reuses rabc_time1/2's own smoothing window:
#' calc_features(vultures_tsbl, c(acc_x, acc_y, acc_z), feature_set = NULL,
#'               winlen_dba_s = 1, fs = 20,
#'               features = list(my_odba = cross_feature(
#'                 c("acc_x", "acc_y", "acc_z"),
#'                 function(x, y, z, winlen) {
#'                   sba <- function(v) roll_static(v, winlen)
#'                   mean(abs(x - sba(x)) + abs(y - sba(y)) + abs(z - sba(z)), na.rm = TRUE)
#'                 }, window = TRUE
#'               )))
calc_features <- function(
  data,
  cols,
  feature_set = "catch22",
  catch24 = FALSE,
  features = NULL,
  winlen_dba_s = NULL,
  winlen_dba = NULL,
  fs = NULL,
  acc_cols = c(x = "acc_x", y = "acc_y", z = "acc_z"),
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
  unknown <- setdiff(feature_set, .feature_sets)
  if (length(unknown)) {
    rlang::abort(c(
      paste0("Unknown `feature_set`: ", paste(unknown, collapse = ", ")),
      i = paste0("Available: ", paste(.feature_sets, collapse = ", "))
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
      i = "Set `feature_set` (e.g. \"catch22\") and/or `features`, a named list of functions."
    ))
  }
  if (all(c("rabc_time1", "rabc_time2") %in% feature_set)) {
    say(
      "rabc_time2 already contains every rabc_time1 feature; computing rabc_time2 only."
    )
    feature_set <- setdiff(feature_set, "rabc_time1")
  }
  uses_time <- any(c("rabc_time1", "rabc_time2") %in% feature_set)
  uses_freq <- "rabc_freq" %in% feature_set
  uses_rabc <- uses_time || uses_freq
  uses_window <- uses_time || uses_custom_window

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
    is_cross <- vapply(features, inherits, logical(1), what = "cross_feature")
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

  # ---- sampling rate and dynamic-acceleration window (validated before any messages)
  fs_ok <- !is.null(fs) &&
    is.numeric(fs) &&
    length(fs) == 1 &&
    !is.na(fs) &&
    fs > 0
  if (uses_freq && !fs_ok) {
    rlang::abort(
      "rabc_freq needs `fs`, the sampling rate in Hz (e.g. `fs = 20`)."
    )
  }
  if (!is.null(fs) && !fs_ok) {
    rlang::abort("`fs` must be a single positive number (sampling rate in Hz).")
  }
  win_requested <- winlen_dba_s
  if (uses_window) {
    winlen_dba <- .resolve_winlen(winlen_dba, winlen_dba_s, fs, fs_ok)
  } else if (!is.null(winlen_dba) || !is.null(winlen_dba_s)) {
    say(
      "`winlen_dba` / `winlen_dba_s` are ignored: not used by any requested feature."
    )
  }

  # ---- accelerometer axes for the rabc sets ------------------------------------
  axes <- NULL
  if (uses_rabc) {
    axes <- .resolve_axes(sig, acc_cols)
    used <- names(axes)
    if (length(axes) < 3 && uses_time) {
      pair_txt <- if (length(axes) == 2) {
        paste0("PDBA_", paste(used, collapse = ""))
      } else {
        paste0("PDBA_", used)
      }
      say(
        "Axes available: ",
        paste(unname(axes), collapse = ", "),
        ". ODBA is replaced by ",
        pair_txt,
        "; cross-axis features are limited to the axes present."
      )
    }
    ignored <- setdiff(sig, axes)
    if (length(ignored) && !"catch22" %in% feature_set) {
      say(
        "rabc sets use accelerometer axes only; ignoring: ",
        paste(ignored, collapse = ", ")
      )
    }
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # ---- bout structure ----------------------------------------------------------
  df <- dplyr::arrange(df, dplyr::across(dplyr::all_of(c(keys, idx))))
  gid <- vctrs::vec_group_id(df[keys]) # consecutive after the sort above
  n_bouts <- max(gid)
  first_idx <- match(seq_len(n_bouts), gid)

  need <- unique(c(
    if ("catch22" %in% feature_set || uses_custom) sig,
    if (uses_rabc) unname(axes)
  ))
  streams <- lapply(stats::setNames(need, need), function(col) {
    unname(split(df[[col]], gid))
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
      affected <- c(
        if (uses_time) {
          paste0(
            "the rolling-window features (ODBA/PDBA",
            if ("rabc_time2" %in% feature_set) {
              ", varsba, vardba, maxdba"
            } else {
              ""
            },
            ")"
          )
        },
        if (uses_custom_window) "window-using custom feature(s)"
      )
      rlang::warn(paste0(
        "Some bouts have fewer samples (",
        shortest,
        ") than the window (",
        winlen_dba,
        "); for those bouts ",
        paste(affected, collapse = " and "),
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

  say(
    "Calculating features for ",
    format(n_bouts, big.mark = ","),
    " bouts: ",
    paste(c(feature_set, if (uses_custom) "custom"), collapse = ", ")
  )

  # ---- feature sets ------------------------------------------------------------
  results <- list()
  if ("catch22" %in% feature_set) {
    results$catch22 <- .run_catch22(streams[sig], bad[sig], catch24, say)
  }
  if ("rabc_time1" %in% feature_set) {
    results$rabc_time1 <- .run_rabc_time(
      streams,
      bad,
      axes,
      winlen_dba,
      expanded = FALSE,
      set = "rabc_time1",
      say = say
    )
  }
  if ("rabc_time2" %in% feature_set) {
    results$rabc_time2 <- .run_rabc_time(
      streams,
      bad,
      axes,
      winlen_dba,
      expanded = TRUE,
      set = "rabc_time2",
      say = say
    )
  }
  if (uses_freq) {
    results$rabc_freq <- .run_rabc_freq(streams, bad, axes, fs, say)
  }
  if (uses_custom) {
    results$custom <- .run_custom(
      streams[sig],
      bad[sig],
      features,
      if (uses_window) winlen_dba,
      say
    )
  }
  results <- results[c(feature_set, if (uses_custom) "custom")] # keep the order requested

  dplyr::bind_cols(c(list(bout_tbl), unname(results)))
}

#' Pivot the wide feature tibble into a long, theft-style layout
#'
#' Splits `<feature_set>__<stream>__<feature>` column names back into three
#' columns. Every other column (keys, context) is kept as an identifier.
#' @param x Output of calc_features().
feature_long <- function(x) {
  n_delim <- lengths(regmatches(
    names(x),
    gregexpr("__", names(x), fixed = TRUE)
  ))
  feat <- names(x)[n_delim == 2]
  if (length(feat) == 0) {
    rlang::abort("No `<feature_set>__<stream>__<feature>` columns found.")
  }
  tidyr::pivot_longer(
    x,
    cols = tidyselect::all_of(feat),
    names_to = c("feature_set", "stream", "feature"),
    names_sep = "__",
    values_to = "value"
  )
}

# ---- runners ---------------------------------------------------------------------

.run_catch22 <- function(streams, bad, catch24, say) {
  fn <- function(v) .catch22_vec(v, catch24)
  tmpl <- names(fn(sin(seq_len(100))))

  per_stream <- lapply(names(streams), function(col) {
    S <- streams[[col]]
    mat <- .fill_matrix(
      seq_along(S),
      which(!bad[[col]]),
      function(b) fn(S[[b]]),
      tmpl
    )
    colnames(mat) <- paste0("catch22__", col, "__", tmpl)
    .report_na(mat, paste0("catch22 (", col, ")"), say)
    tibble::as_tibble(mat)
  })
  dplyr::bind_cols(per_stream)
}

.run_rabc_time <- function(streams, bad, axes, winlen, expanded, set, say) {
  n <- length(streams[[unname(axes[1])]])
  fn <- function(b) {
    ax <- lapply(axes, function(col) streams[[col]][[b]])
    .rabc_time_bout(ax, axes, winlen, expanded)
  }
  dummy <- lapply(axes, function(.) sin(seq_len(max(2 * winlen, 20)) * 0.7))
  tmpl <- names(.rabc_time_bout(dummy, axes, winlen, expanded))

  bad_any <- Reduce(`|`, bad[unname(axes)])
  mat <- .fill_matrix(seq_len(n), which(!bad_any), fn, tmpl)
  colnames(mat) <- paste0(set, "__", tmpl)
  .report_na(mat, set, say)
  tibble::as_tibble(mat)
}

.run_rabc_freq <- function(streams, bad, axes, fs, say) {
  n <- length(streams[[unname(axes[1])]])
  fn <- function(b) {
    ax <- lapply(axes, function(col) streams[[col]][[b]])
    .rabc_freq_bout(ax, axes, fs)
  }
  dummy <- lapply(axes, function(.) sin(seq_len(40) * 0.7))
  tmpl <- names(.rabc_freq_bout(dummy, axes, fs))

  bad_any <- Reduce(`|`, bad[unname(axes)])
  mat <- .fill_matrix(seq_len(n), which(!bad_any), fn, tmpl)
  colnames(mat) <- paste0("rabc_freq__", tmpl)
  .report_na(mat, "rabc_freq", say)
  tibble::as_tibble(mat)
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
          if (length(r) != 1 || !is.numeric(r)) NA_real_ else as.numeric(r)
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
          r <- tryCatch(do.call(e$fn, args), error = function(err) NA_real_)
          if (length(r) != 1 || !is.numeric(r)) NA_real_ else as.numeric(r)
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
          paste0("custom (", nm, ", using ", paste(cls, collapse = ", "), ")"),
          say
        )
        tibble::as_tibble(mat)
      })
    )
  }

  dplyr::bind_cols(out)
}

# ---- per-bout feature functions ----------------------------------------------------

.catch22_vec <- function(v, catch24) {
  o <- Rcatch22::catch22_all(v, catch24 = catch24)
  stats::setNames(o$values, o$names)
}

# Rolling mean used to split static (sba) from dynamic (dba) acceleration.
.roll_mean <- function(v, w) {
  if (length(v) < w) {
    return(rep(NA_real_, length(v)))
  }
  zoo::rollapply(v, width = w, FUN = mean, fill = NA, align = "center")
}

# ax:   named list of numeric vectors for the axes present (names in x, y, z)
# cmap: named character vector mapping those axes to column names
.rabc_time_bout <- function(ax, cmap, winlen, expanded) {
  ns <- names(ax)
  st <- unname(cmap[ns])
  per <- function(f) unname(vapply(ax, f, numeric(1)))
  nm <- function(vals, streams, feature) {
    if (length(vals) == 0) {
      return(numeric(0))
    }
    stats::setNames(unname(vals), paste0(streams, "__", feature))
  }

  sba <- lapply(ax, .roll_mean, w = winlen)
  dba <- Map(function(v, s) abs(v - s), ax, sba)

  # ODBA with all three axes; PDBA_<axes> when only some are available
  dyn <- mean(Reduce(`+`, dba), na.rm = TRUE)
  dyn_name <- if (length(ns) == 3) {
    "ODBA"
  } else {
    paste0("PDBA_", paste(ns, collapse = ""))
  }
  dyn_stream <- paste(st, collapse = "_")

  mx <- per(max)
  mn <- per(min)
  out <- c(
    nm(per(mean), st, "mean"),
    nm(per(stats::var), st, "variance"),
    nm(per(stats::sd), st, "sd"),
    nm(mx, st, "max"),
    nm(mn, st, "min"),
    nm(mx - mn, st, "range"),
    nm(dyn, dyn_stream, dyn_name)
  )
  if (!expanded) {
    return(out)
  }

  # cross-axis pairs, direction as in the original code: xy, yz, xz
  prs <- Filter(
    function(p) all(p %in% ns),
    list(c("x", "y"), c("y", "z"), c("x", "z"))
  )
  ps <- vapply(prs, function(p) paste(cmap[p], collapse = "_"), character(1))
  pair <- function(f) {
    vapply(prs, function(p) f(ax[[p[1]]], ax[[p[2]]]), numeric(1))
  }

  suppressWarnings(c(
    out,
    nm(per(function(v) sqrt(sum(v^2))), st, "norm"),
    nm(pair(stats::cov), ps, "cov"),
    nm(pair(stats::cor), ps, "cor"),
    nm(pair(function(a, b) mean(a - b)), ps, "meandiff"),
    nm(pair(function(a, b) stats::sd(a - b)), ps, "sddiff"),
    nm(vapply(sba, stats::var, numeric(1), na.rm = TRUE), st, "varsba"),
    nm(vapply(dba, stats::var, numeric(1), na.rm = TRUE), st, "vardba"),
    nm(vapply(dba, max, numeric(1), na.rm = TRUE), st, "maxdba")
  ))
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

.rabc_freq_bout <- function(ax, cmap, fs) {
  st <- unname(cmap[names(ax)])
  res <- suppressWarnings(vapply(ax, .max_freq_amp, numeric(3)))
  res <- matrix(res, nrow = 3)
  fi <- fs / length(ax[[1]]) # as in rabc: frequency = bin * fs / n
  c(
    stats::setNames(res[1, ] * fi, paste0(st, "__freqmain")),
    stats::setNames(res[2, ], paste0(st, "__freqamp")),
    stats::setNames(res[3, ], paste0(st, "__entropy"))
  )
}

# ---- helpers -------------------------------------------------------------------

# Resolve the dynamic-acceleration window to a number of samples.
.resolve_winlen <- function(w, w_s, fs, fs_ok) {
  if (!is.null(w) && !is.null(w_s)) {
    rlang::abort(
      "Give either `winlen_dba_s` (seconds) or `winlen_dba` (samples), not both."
    )
  }
  if (is.null(w) && is.null(w_s)) {
    rlang::abort(c(
      "A smoothing window is needed (used by rabc_time1/2 and/or any cross_feature(window = TRUE)).",
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
  if (!is.numeric(w) || length(w) != 1 || is.na(w) || w < 2 || w != round(w)) {
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
      if (!is.null(requested_s)) paste0("; requested ", requested_s, " s"),
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
  mat <- matrix(NA_real_, nrow = length(items), ncol = length(nms))
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

# Which of the requested axis columns are among the selected streams?
.resolve_axes <- function(sig, acc_cols) {
  if (is.null(names(acc_cols)) || !all(names(acc_cols) %in% c("x", "y", "z"))) {
    rlang::abort(
      "`acc_cols` must be a named character vector with names from x, y, z."
    )
  }
  present <- acc_cols[acc_cols %in% sig]
  if (length(present) == 0) {
    rlang::abort(c(
      "The rabc feature sets need at least one accelerometer axis among `cols`.",
      i = paste0(
        "Looked for: ",
        paste(acc_cols, collapse = ", "),
        " (change `acc_cols` if yours are named differently)."
      ),
      i = paste0("Selected: ", paste(sig, collapse = ", "))
    ))
  }
  present[intersect(c("x", "y", "z"), names(present))]
}

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

# The exact centred rolling mean rabc_time1/2 use to split static from dynamic
# acceleration -- exported so a custom cross_feature(window = TRUE) can build
# the same split. NA at both edges ((winlen - 1) / 2 samples per side).
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
    missing(fs) || !is.numeric(fs) || length(fs) != 1 || is.na(fs) || fs <= 0
  ) {
    rlang::abort("`fs`, the sampling rate in Hz, is required (e.g. `fs = 20`).")
  }
  if (!is.numeric(winlen_s) || any(is.na(winlen_s)) || any(winlen_s <= 0)) {
    rlang::abort("`winlen_s` must be positive numbers of seconds.")
  }

  keys <- tsibble::key_vars(data)
  df <- tibble::as_tibble(data)
  sig <- names(tidyselect::eval_select(rlang::enquo(cols), df))
  has_label <- label_col %in% names(df)
  if (!has_label) {
    say("No `", label_col, "` column: separation statistics will be skipped.")
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
  say("Windows tested: ", paste0(ws, " s (", n_s, " samples)", collapse = ", "))

  # ---- ODBA / PDBA at every window -------------------------------------------
  long <- purrr::map2(ws, n_s, function(s, k) {
    f <- calc_features(
      data,
      dplyr::all_of(sig),
      feature_set = "rabc_time1",
      winlen_dba_s = s,
      fs = fs,
      verbose = FALSE
    )
    dyn <- grep("__(ODBA|PDBA_[xyz]+)$", names(f), value = TRUE)
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
      function(s) stats::median(long$dba[long$winlen_s == s], na.rm = TRUE),
      numeric(1)
    ),
    rank_cor_ref = unname(rank_cor[, as.character(ref_s)])
  ) |>
    dplyr::left_join(sep, by = "winlen_s") |>
    dplyr::mutate(is_ref = winlen_s == ref_s)
  names(overview)[names(overview) == "median_dba"] <- paste0("median_", metric)

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
        sep = .separability(d$dba[d$label == p[1]], d$dba[d$label == p[2]])
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
        title = paste0(metric, " versus smoothing window, by behaviour"),
        subtitle = "Dotted line = reference window"
      ) +
      theme
  }

  sep_long <- overview |>
    dplyr::select(winlen_s, epsilon2, median_pair_sep, min_pair_sep) |>
    tidyr::pivot_longer(-winlen_s, names_to = "stat", values_to = "value") |>
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
    out$separation <- ggplot2::ggplot(sep_long, ggplot2::aes(winlen_s, value)) +
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
    tidyr::pivot_longer(-winlen_s, names_to = "stat", values_to = "value") |>
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

# project_features.R -----------------------------------------------------------
# Dimension reduction (project_features()) and its plot() method for the wide
# tibble calc_features() returns: one row per bout, feature columns named
# <feature_set>__<stream>__<feature>.
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
# Column selection: `feature_set` filters by the SAME prefix calc_features()
# writes into every feature column name (e.g. "catch22", "rabc_time2"), so
# projecting a single feature set, or a chosen few, is `feature_set = "catch22"`
# / `feature_set = c("catch22", "rabc_time2")`. NULL uses every feature column
# present, whichever set(s) they came from.

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

#' Reduce a calc_features() feature tibble to two dimensions
#'
#' @param data           Output of calc_features(): one row per bout/sample,
#'                       feature columns named <feature_set>__<stream>__<feature>.
#' @param feature_set    Character vector of feature_set(s) to project (e.g.
#'                       "catch22", or c("catch22", "rabc_time2")). NULL
#'                       (default) uses every feature column present.
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
      "No `<feature_set>__<stream>__<feature>` columns found. Is `data` the output of calc_features()?"
    )
  }
  if (!is.null(feature_set)) {
    available <- unique(parts$info$feature_set)
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
    keep <- parts$info$feature_set %in% feature_set
    parts$info <- parts$info[keep, ]
    parts$feature <- parts$info$column
  }
  if (length(parts$feature) < 2) {
    rlang::abort("Need at least 2 feature columns to project to 2 dimensions.")
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
#' @param colour_by    Name of a column in `x$Meta` to colour points by (e.g.
#'                     "label", "id"). Defaults to "label" if present, else no
#'                     colour.
#' @param show_covariance Draw a covariance ellipse per colour group. Ignored
#'                     if `colour_by` is NULL.
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

  if (is.null(colour_by) && "label" %in% names(x$Meta)) {
    colour_by <- "label"
  }
  if (!is.null(colour_by) && !colour_by %in% names(x$Meta)) {
    rlang::abort(c(
      paste0("`colour_by` not found: ", colour_by),
      i = paste0(
        "Available: ",
        paste(setdiff(names(x$Meta), "..row_id.."), collapse = ", ")
      )
    ))
  }
  if (!is.null(colour_by)) {
    grp <- stats::setNames(
      x$Meta[c("..row_id..", colour_by)],
      c("id", "group_id")
    )
    fits <- dplyr::inner_join(fits, grp, by = "id")
    fits$group_id <- factor(fits$group_id)
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

  if (!is.null(colour_by)) {
    if (isTRUE(show_covariance)) {
      p <- p +
        ggplot2::stat_ellipse(
          ggplot2::aes(fill = .data$group_id),
          geom = "polygon",
          alpha = 0.2
        ) +
        ggplot2::guides(fill = "none") +
        ggplot2::scale_fill_brewer(palette = "Dark2")
    }
    p <- p +
      ggplot2::geom_point(
        ggplot2::aes(colour = .data$group_id),
        size = pt_size
      ) +
      ggplot2::scale_colour_brewer(palette = "Dark2") +
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

# misc / old functions
# old calc_features before adding custom features and custom cross-features

# # calc_features --------------------------------------------------------
# # Bout-level feature calculation for a tidy sensor tsibble
# # (key = c(id, bout), index = sample_in_bout or datetime, one column per stream).
# # Modelled on theft::calc_features(), but multivariate: any sensor columns
# # can be chosen, and feature sets that combine axes (ODBA, cor_x_y, ...) are
# # supported.
# #
# # Requires: tsibble, tibble, dplyr, tidyr, tidyselect, rlang, vctrs, zoo,
# #           entropy, Rcatch22
# #
# # NOTE: this is deliberately named calc_features(), as requested. If you
# # also attach {theft}, whichever is attached/sourced last masks the other; call
# # theft's version as theft::calc_features().
# #
# # Feature sets
# #   catch22     Rcatch22 features on EVERY selected stream (any numeric
# #               column). `catch24 = TRUE` adds DN_Mean and DN_Spread_Std.
# #   rabc_time1  Original rabc time-domain set: per-axis mean, variance, sd,
# #               max, min, range, plus ODBA.
# #   rabc_time2  Expanded set = rabc_time1 + per-axis norm, cross-axis cov / cor /
# #               meandiff / sddiff, and per-axis varsba / vardba / maxdba.
# #               It CONTAINS rabc_time1, so if both are requested only
# #               rabc_time2 is computed (with a message).
# #   rabc_freq   Per-axis dominant frequency (freqmain), its amplitude
# #               (freqamp) and a spectral entropy term (entropy).
# #
# # Accelerometer axes for the rabc sets
# #   The rabc sets read the columns named in `acc_cols` (default acc_x, acc_y,
# #   acc_z) from the selected `cols`. Missing axes do not cause an error:
# #     * x, y and z present -> ODBA
# #     * two axes present   -> PDBA_xy / PDBA_xz / PDBA_yz in place of ODBA, and
# #                             only the cross-axis features for that pair
# #     * one axis present   -> PDBA_x / PDBA_y / PDBA_z, no cross-axis features
# #   An error is raised only if NONE of the axes is selected.
# #
# # Output (a plain tibble, one row per bout)
# #   key columns (e.g. id, bout), then context columns (e.g. attachment, label;
# #   any non-selected column that is constant within every bout), then feature
# #   columns named   <feature_set>__<stream>__<feature>
# #   e.g. catch22__acc_x__DN_HistogramMode_5, rabc_time2__acc_x_acc_y__cor,
# #        rabc_time1__acc_x_acc_y_acc_z__ODBA, rabc_time1__acc_x_acc_z__PDBA_xz
# #   Use feature_long() to pivot into a long theft-style layout.
# #
# # Dynamic-acceleration window (rabc_time sets)
# #   ODBA/PDBA, varsba, vardba and maxdba split each axis into a static part (a
# #   centred moving average) and a dynamic part (raw minus static). Set the
# #   window with EITHER
# #     winlen_dba_s  seconds (needs `fs`); the recommended, dataset-independent
# #                   way. Converted to an ODD number of samples:
# #                   n = 2 * floor(winlen_dba_s * fs / 2) + 1, so the window is
# #                   centred on a sample. e.g. 1 s at 20 Hz -> 21 samples.
# #     winlen_dba    samples; use to reproduce earlier rabc results exactly.
# #   Giving both is an error; there is no default. The resolved window and the
# #   share of a median-length bout lost to edge NAs are always reported.
# #
# # Missing values: a bout with any NA / non-finite value in a stream a set
# #   needs gets NA for that set's features (per stream for catch22; for the whole
# #   set for rabc sets, since features combine axes). Features that are
# #   undefined (e.g. a correlation on a constant signal) are also NA. A message
# #   reports how many bouts are affected.
# #
# # `seed`: every set implemented so far is deterministic. The argument exists
# #   for parity with theft::calc_features() and future stochastic sets.

# .feature_sets <- c("catch22", "rabc_time1", "rabc_time2", "rabc_freq")

# #' Calculate bout-level features from a sensor tsibble
# #'
# #' @param data        A tsibble. Every key combination (e.g. id + bout) is one
# #'                    bout and yields one row of features.
# #' @param cols        Sensor streams to use, as tidyselect: `acc_x`,
# #'                    `c(acc_x, acc_y, acc_z)`, `starts_with("acc_")`, ...
# #' @param feature_set One or more of "catch22", "rabc_time1", "rabc_time2",
# #'                    "rabc_freq".
# #' @param catch24     Also compute catch24 (adds DN_Mean, DN_Spread_Std) when
# #'                    catch22 is requested.
# #' @param winlen_dba_s Rolling-mean window in SECONDS used to separate static from
# #'                    dynamic acceleration (rabc_time1/2). Needs `fs`. Rounded to
# #'                    an odd number of samples. Give this or `winlen_dba`.
# #' @param winlen_dba  The same window in SAMPLES. Give this or `winlen_dba_s`.
# #' @param fs          Sampling rate in Hz. Required for rabc_freq and for
# #'                    `winlen_dba_s`.
# #' @param acc_cols    Named character vector mapping axes to column names.
# #'                    Change it if your axes are not called acc_x/acc_y/acc_z.
# #' @param seed        Integer for set.seed(), or NULL to leave the RNG alone.
# #' @param verbose     Print progress messages.
# #'
# #' @return A tibble with one row per bout (see file header for the layout).
# #' @examples

# calc_features <- function(
#   data,
#   cols,
#   feature_set = "catch22",
#   catch24 = FALSE,
#   winlen_dba_s = NULL,
#   winlen_dba = NULL,
#   fs = NULL,
#   acc_cols = c(x = "acc_x", y = "acc_y", z = "acc_z"),
#   seed = 123,
#   verbose = TRUE
# ) {
#   say <- function(...) if (verbose) rlang::inform(paste0(...))

#   # ---- input checks ----------------------------------------------------------
#   if (!tsibble::is_tsibble(data)) {
#     rlang::abort("`data` must be a tsibble in the tidy sensor format.")
#   }
#   if (missing(cols)) {
#     rlang::abort(c(
#       "Say which sensor columns to use.",
#       i = "e.g. `cols = c(acc_x, acc_y, acc_z)`."
#     ))
#   }

#   feature_set <- unique(tolower(feature_set))
#   unknown <- setdiff(feature_set, .feature_sets)
#   if (length(unknown)) {
#     rlang::abort(c(
#       paste0("Unknown `feature_set`: ", paste(unknown, collapse = ", ")),
#       i = paste0("Available: ", paste(.feature_sets, collapse = ", "))
#     ))
#   }
#   if (all(c("rabc_time1", "rabc_time2") %in% feature_set)) {
#     say(
#       "rabc_time2 already contains every rabc_time1 feature; computing rabc_time2 only."
#     )
#     feature_set <- setdiff(feature_set, "rabc_time1")
#   }
#   uses_time <- any(c("rabc_time1", "rabc_time2") %in% feature_set)
#   uses_freq <- "rabc_freq" %in% feature_set
#   uses_rabc <- uses_time || uses_freq

#   keys <- tsibble::key_vars(data)
#   idx <- tsibble::index_var(data)
#   if (length(keys) == 0) {
#     rlang::abort(c(
#       "`data` needs key column(s) that identify each bout.",
#       i = "e.g. `as_tsibble(x, key = c(id, bout), index = sample_in_bout)`."
#     ))
#   }
#   df <- tibble::as_tibble(data)

#   sig <- names(tidyselect::eval_select(rlang::enquo(cols), df))
#   if (length(sig) == 0) {
#     rlang::abort("`cols` selected no columns.")
#   }
#   if (any(sig %in% c(keys, idx))) {
#     rlang::abort("`cols` must be sensor columns, not key or index columns.")
#   }
#   bad_type <- sig[!vapply(df[sig], is.numeric, logical(1))]
#   if (length(bad_type)) {
#     rlang::abort(paste0(
#       "These `cols` are not numeric: ",
#       paste(bad_type, collapse = ", ")
#     ))
#   }

#   if (isTRUE(catch24) && !"catch22" %in% feature_set) {
#     rlang::warn(
#       "`catch24 = TRUE` is ignored because \"catch22\" is not in `feature_set`."
#     )
#   }

#   # ---- sampling rate and dynamic-acceleration window (validated before any messages)
#   fs_ok <- !is.null(fs) &&
#     is.numeric(fs) &&
#     length(fs) == 1 &&
#     !is.na(fs) &&
#     fs > 0
#   if (uses_freq && !fs_ok) {
#     rlang::abort(
#       "rabc_freq needs `fs`, the sampling rate in Hz (e.g. `fs = 20`)."
#     )
#   }
#   if (!is.null(fs) && !fs_ok) {
#     rlang::abort("`fs` must be a single positive number (sampling rate in Hz).")
#   }
#   win_requested <- winlen_dba_s
#   if (uses_time) {
#     winlen_dba <- .resolve_winlen(winlen_dba, winlen_dba_s, fs, fs_ok)
#   } else if (!is.null(winlen_dba) || !is.null(winlen_dba_s)) {
#     say(
#       "`winlen_dba` / `winlen_dba_s` are ignored: no rabc_time set requested."
#     )
#   }

#   # ---- accelerometer axes for the rabc sets ------------------------------------
#   axes <- NULL
#   if (uses_rabc) {
#     axes <- .resolve_axes(sig, acc_cols)
#     used <- names(axes)
#     if (length(axes) < 3 && uses_time) {
#       pair_txt <- if (length(axes) == 2) {
#         paste0("PDBA_", paste(used, collapse = ""))
#       } else {
#         paste0("PDBA_", used)
#       }
#       say(
#         "Axes available: ",
#         paste(unname(axes), collapse = ", "),
#         ". ODBA is replaced by ",
#         pair_txt,
#         "; cross-axis features are limited to the axes present."
#       )
#     }
#     ignored <- setdiff(sig, axes)
#     if (length(ignored) && !"catch22" %in% feature_set) {
#       say(
#         "rabc sets use accelerometer axes only; ignoring: ",
#         paste(ignored, collapse = ", ")
#       )
#     }
#   }
#   if (!is.null(seed)) {
#     set.seed(seed)
#   }

#   # ---- bout structure ----------------------------------------------------------
#   df <- dplyr::arrange(df, dplyr::across(dplyr::all_of(c(keys, idx))))
#   gid <- vctrs::vec_group_id(df[keys]) # consecutive after the sort above
#   n_bouts <- max(gid)
#   first_idx <- match(seq_len(n_bouts), gid)

#   need <- unique(c(
#     if ("catch22" %in% feature_set) sig,
#     if (uses_rabc) unname(axes)
#   ))
#   streams <- lapply(stats::setNames(need, need), function(col) {
#     unname(split(df[[col]], gid))
#   })
#   bad <- lapply(streams, function(s) {
#     vapply(s, function(v) any(!is.finite(v)), logical(1))
#   })

#   if (uses_time) {
#     n_per_bout <- lengths(streams[[unname(axes[1])]])
#     say(.window_note(
#       winlen_dba,
#       win_requested,
#       if (fs_ok) fs,
#       stats::median(n_per_bout)
#     ))
#     shortest <- min(n_per_bout)
#     if (shortest < winlen_dba) {
#       rlang::warn(paste0(
#         "Some bouts have fewer samples (",
#         shortest,
#         ") than the window (",
#         winlen_dba,
#         "); for those bouts the rolling-window features (ODBA/PDBA",
#         if ("rabc_time2" %in% feature_set) ", varsba, vardba, maxdba" else "",
#         ") will be NA."
#       ))
#     }
#   }

#   # ---- key + context columns ---------------------------------------------------
#   ctx_candidates <- setdiff(names(df), c(keys, idx, sig))
#   is_const <- function(col) {
#     v <- df[[col]]
#     first <- v[first_idx][gid]
#     all(((v == first) %in% TRUE) | (is.na(v) & is.na(first)))
#   }
#   const <- vapply(ctx_candidates, is_const, logical(1))
#   varying <- ctx_candidates[!const]
#   is_num <- vapply(df[ctx_candidates], is.numeric, logical(1))
#   vary_num <- varying[is_num[varying]]
#   vary_oth <- setdiff(varying, vary_num)
#   if (length(vary_oth)) {
#     rlang::warn(paste0(
#       "Dropped context column(s) that change within a bout: ",
#       paste(vary_oth, collapse = ", "),
#       ". Bout-level values would be ambiguous."
#     ))
#   }
#   if (length(vary_num)) {
#     say(
#       "Dropped numeric column(s) that vary within bouts (treated as unselected sensor streams): ",
#       paste(vary_num, collapse = ", ")
#     )
#   }
#   bout_tbl <- df[first_idx, c(keys, ctx_candidates[const]), drop = FALSE]

#   say(
#     "Calculating features for ",
#     format(n_bouts, big.mark = ","),
#     " bouts: ",
#     paste(feature_set, collapse = ", ")
#   )

#   # ---- feature sets ------------------------------------------------------------
#   results <- list()
#   if ("catch22" %in% feature_set) {
#     results$catch22 <- .run_catch22(streams[sig], bad[sig], catch24, say)
#   }
#   if ("rabc_time1" %in% feature_set) {
#     results$rabc_time1 <- .run_rabc_time(
#       streams,
#       bad,
#       axes,
#       winlen_dba,
#       expanded = FALSE,
#       set = "rabc_time1",
#       say = say
#     )
#   }
#   if ("rabc_time2" %in% feature_set) {
#     results$rabc_time2 <- .run_rabc_time(
#       streams,
#       bad,
#       axes,
#       winlen_dba,
#       expanded = TRUE,
#       set = "rabc_time2",
#       say = say
#     )
#   }
#   if (uses_freq) {
#     results$rabc_freq <- .run_rabc_freq(streams, bad, axes, fs, say)
#   }
#   results <- results[feature_set] # keep the order the user asked for

#   dplyr::bind_cols(c(list(bout_tbl), unname(results)))
# }

# #' Pivot the wide feature tibble into a long, theft-style layout
# #'
# #' Splits `<feature_set>__<stream>__<feature>` column names back into three
# #' columns. Every other column (keys, context) is kept as an identifier.
# #' @param x Output of calc_features().
# feature_long <- function(x) {
#   n_delim <- lengths(regmatches(
#     names(x),
#     gregexpr("__", names(x), fixed = TRUE)
#   ))
#   feat <- names(x)[n_delim == 2]
#   if (length(feat) == 0) {
#     rlang::abort("No `<feature_set>__<stream>__<feature>` columns found.")
#   }
#   tidyr::pivot_longer(
#     x,
#     cols = tidyselect::all_of(feat),
#     names_to = c("feature_set", "stream", "feature"),
#     names_sep = "__",
#     values_to = "value"
#   )
# }

# # ---- runners ---------------------------------------------------------------------

# .run_catch22 <- function(streams, bad, catch24, say) {
#   fn <- function(v) .catch22_vec(v, catch24)
#   tmpl <- names(fn(sin(seq_len(100))))

#   per_stream <- lapply(names(streams), function(col) {
#     S <- streams[[col]]
#     mat <- .fill_matrix(
#       seq_along(S),
#       which(!bad[[col]]),
#       function(b) fn(S[[b]]),
#       tmpl
#     )
#     colnames(mat) <- paste0("catch22__", col, "__", tmpl)
#     .report_na(mat, paste0("catch22 (", col, ")"), say)
#     tibble::as_tibble(mat)
#   })
#   dplyr::bind_cols(per_stream)
# }

# .run_rabc_time <- function(streams, bad, axes, winlen, expanded, set, say) {
#   n <- length(streams[[unname(axes[1])]])
#   fn <- function(b) {
#     ax <- lapply(axes, function(col) streams[[col]][[b]])
#     .rabc_time_bout(ax, axes, winlen, expanded)
#   }
#   dummy <- lapply(axes, function(.) sin(seq_len(max(2 * winlen, 20)) * 0.7))
#   tmpl <- names(.rabc_time_bout(dummy, axes, winlen, expanded))

#   bad_any <- Reduce(`|`, bad[unname(axes)])
#   mat <- .fill_matrix(seq_len(n), which(!bad_any), fn, tmpl)
#   colnames(mat) <- paste0(set, "__", tmpl)
#   .report_na(mat, set, say)
#   tibble::as_tibble(mat)
# }

# .run_rabc_freq <- function(streams, bad, axes, fs, say) {
#   n <- length(streams[[unname(axes[1])]])
#   fn <- function(b) {
#     ax <- lapply(axes, function(col) streams[[col]][[b]])
#     .rabc_freq_bout(ax, axes, fs)
#   }
#   dummy <- lapply(axes, function(.) sin(seq_len(40) * 0.7))
#   tmpl <- names(.rabc_freq_bout(dummy, axes, fs))

#   bad_any <- Reduce(`|`, bad[unname(axes)])
#   mat <- .fill_matrix(seq_len(n), which(!bad_any), fn, tmpl)
#   colnames(mat) <- paste0("rabc_freq__", tmpl)
#   .report_na(mat, "rabc_freq", say)
#   tibble::as_tibble(mat)
# }

# # ---- per-bout feature functions ----------------------------------------------------

# .catch22_vec <- function(v, catch24) {
#   o <- Rcatch22::catch22_all(v, catch24 = catch24)
#   stats::setNames(o$values, o$names)
# }

# # Rolling mean used to split static (sba) from dynamic (dba) acceleration.
# .roll_mean <- function(v, w) {
#   if (length(v) < w) {
#     return(rep(NA_real_, length(v)))
#   }
#   zoo::rollapply(v, width = w, FUN = mean, fill = NA, align = "center")
# }

# # ax:   named list of numeric vectors for the axes present (names in x, y, z)
# # cmap: named character vector mapping those axes to column names
# .rabc_time_bout <- function(ax, cmap, winlen, expanded) {
#   ns <- names(ax)
#   st <- unname(cmap[ns])
#   per <- function(f) unname(vapply(ax, f, numeric(1)))
#   nm <- function(vals, streams, feature) {
#     if (length(vals) == 0) {
#       return(numeric(0))
#     }
#     stats::setNames(unname(vals), paste0(streams, "__", feature))
#   }

#   sba <- lapply(ax, .roll_mean, w = winlen)
#   dba <- Map(function(v, s) abs(v - s), ax, sba)

#   # ODBA with all three axes; PDBA_<axes> when only some are available
#   dyn <- mean(Reduce(`+`, dba), na.rm = TRUE)
#   dyn_name <- if (length(ns) == 3) {
#     "ODBA"
#   } else {
#     paste0("PDBA_", paste(ns, collapse = ""))
#   }
#   dyn_stream <- paste(st, collapse = "_")

#   mx <- per(max)
#   mn <- per(min)
#   out <- c(
#     nm(per(mean), st, "mean"),
#     nm(per(stats::var), st, "variance"),
#     nm(per(stats::sd), st, "sd"),
#     nm(mx, st, "max"),
#     nm(mn, st, "min"),
#     nm(mx - mn, st, "range"),
#     nm(dyn, dyn_stream, dyn_name)
#   )
#   if (!expanded) {
#     return(out)
#   }

#   # cross-axis pairs, direction as in the original code: xy, yz, xz
#   prs <- Filter(
#     function(p) all(p %in% ns),
#     list(c("x", "y"), c("y", "z"), c("x", "z"))
#   )
#   ps <- vapply(prs, function(p) paste(cmap[p], collapse = "_"), character(1))
#   pair <- function(f) {
#     vapply(prs, function(p) f(ax[[p[1]]], ax[[p[2]]]), numeric(1))
#   }

#   suppressWarnings(c(
#     out,
#     nm(per(function(v) sqrt(sum(v^2))), st, "norm"),
#     nm(pair(stats::cov), ps, "cov"),
#     nm(pair(stats::cor), ps, "cor"),
#     nm(pair(function(a, b) mean(a - b)), ps, "meandiff"),
#     nm(pair(function(a, b) stats::sd(a - b)), ps, "sddiff"),
#     nm(vapply(sba, stats::var, numeric(1), na.rm = TRUE), st, "varsba"),
#     nm(vapply(dba, stats::var, numeric(1), na.rm = TRUE), st, "vardba"),
#     nm(vapply(dba, max, numeric(1), na.rm = TRUE), st, "maxdba")
#   ))
# }

# # Dominant frequency bin, its amplitude, and a spectral-entropy term for one
# # detrended series. Follows rabc::max_freq_amp; ties take the first maximum.
# .max_freq_amp <- function(v) {
#   n <- length(v)
#   half <- floor(n / 2)
#   if (half < 1) {
#     return(rep(NA_real_, 3))
#   }
#   freq <- abs(stats::fft(stats::lm(as.numeric(v) ~ seq_len(n))$residuals))
#   ind <- which.max(freq[seq_len(half)])
#   c(ind, freq[ind], entropy::entropy(freq[seq_len(half)])^2 / half)
# }

# .rabc_freq_bout <- function(ax, cmap, fs) {
#   st <- unname(cmap[names(ax)])
#   res <- suppressWarnings(vapply(ax, .max_freq_amp, numeric(3)))
#   res <- matrix(res, nrow = 3)
#   fi <- fs / length(ax[[1]]) # as in rabc: frequency = bin * fs / n
#   c(
#     stats::setNames(res[1, ] * fi, paste0(st, "__freqmain")),
#     stats::setNames(res[2, ], paste0(st, "__freqamp")),
#     stats::setNames(res[3, ], paste0(st, "__entropy"))
#   )
# }

# # ---- helpers -------------------------------------------------------------------

# # Resolve the dynamic-acceleration window to a number of samples.
# .resolve_winlen <- function(w, w_s, fs, fs_ok) {
#   if (!is.null(w) && !is.null(w_s)) {
#     rlang::abort(
#       "Give either `winlen_dba_s` (seconds) or `winlen_dba` (samples), not both."
#     )
#   }
#   if (is.null(w) && is.null(w_s)) {
#     rlang::abort(c(
#       "The rabc_time sets need a window for the static-acceleration estimate.",
#       i = "Set it in seconds (needs `fs`): `winlen_dba_s = 1, fs = 20`.",
#       i = "Or in samples: `winlen_dba = 21`."
#     ))
#   }
#   if (!is.null(w_s)) {
#     if (!is.numeric(w_s) || length(w_s) != 1 || is.na(w_s) || w_s <= 0) {
#       rlang::abort(
#         "`winlen_dba_s` must be a single positive number of seconds."
#       )
#     }
#     if (!fs_ok) {
#       rlang::abort(
#         "`winlen_dba_s` needs `fs`, the sampling rate in Hz (e.g. `fs = 20`)."
#       )
#     }
#     n <- 2 * floor(w_s * fs / 2 + 1e-9) + 1 # odd, so the window centres on a sample
#     if (n < 3) {
#       rlang::abort(paste0(
#         w_s,
#         " s at ",
#         fs,
#         " Hz is under 3 samples; use a longer `winlen_dba_s`."
#       ))
#     }
#     return(as.integer(n))
#   }
#   if (!is.numeric(w) || length(w) != 1 || is.na(w) || w < 2 || w != round(w)) {
#     rlang::abort(
#       "`winlen_dba` must be a whole number of samples >= 2 (e.g. `winlen_dba = 21`)."
#     )
#   }
#   as.integer(w)
# }

# # One-line description of the resolved window and its edge cost.
# .window_note <- function(n, requested_s, fs, median_len) {
#   win <- if (is.null(fs)) {
#     paste0(n, " samples")
#   } else {
#     paste0(
#       n,
#       " samples (",
#       format(round(n / fs, 2), nsmall = 2),
#       " s at ",
#       fs,
#       " Hz",
#       if (!is.null(requested_s)) paste0("; requested ", requested_s, " s"),
#       ")"
#     )
#   }
#   lost <- min(1, (n - 1) / median_len)
#   paste0(
#     "rabc_time window: ",
#     win,
#     ". Edge NAs remove ",
#     round(100 * lost),
#     "% of a median-length bout (",
#     median_len,
#     " samples)."
#   )
# }

# # Bouts x features matrix; rows not in `ok` stay NA; non-finite results -> NA.
# .fill_matrix <- function(items, ok, fn, nms) {
#   mat <- matrix(NA_real_, nrow = length(items), ncol = length(nms))
#   if (length(ok)) {
#     mat[ok, ] <- do.call(rbind, lapply(items[ok], fn))
#   }
#   mat[!is.finite(mat)] <- NA_real_
#   mat
# }

# .report_na <- function(mat, label, say) {
#   n_na <- sum(!stats::complete.cases(mat))
#   if (n_na > 0) {
#     say(
#       label,
#       ": ",
#       format(n_na, big.mark = ","),
#       " of ",
#       format(nrow(mat), big.mark = ","),
#       " bouts have at least one NA feature (missing/non-finite input or an undefined feature)."
#     )
#   }
#   invisible(n_na)
# }

# # Which of the requested axis columns are among the selected streams?
# .resolve_axes <- function(sig, acc_cols) {
#   if (is.null(names(acc_cols)) || !all(names(acc_cols) %in% c("x", "y", "z"))) {
#     rlang::abort(
#       "`acc_cols` must be a named character vector with names from x, y, z."
#     )
#   }
#   present <- acc_cols[acc_cols %in% sig]
#   if (length(present) == 0) {
#     rlang::abort(c(
#       "The rabc feature sets need at least one accelerometer axis among `cols`.",
#       i = paste0(
#         "Looked for: ",
#         paste(acc_cols, collapse = ", "),
#         " (change `acc_cols` if yours are named differently)."
#       ),
#       i = paste0("Selected: ", paste(sig, collapse = ", "))
#     ))
#   }
#   present[intersect(c("x", "y", "z"), names(present))]
# }
