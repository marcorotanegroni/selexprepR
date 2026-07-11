.qc_plot_colours <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9",
    "#000000")
.qc_plot_unavailable <- function(title, message) {
    graphics::plot.new()
    graphics::title(main = title)
    graphics::text(0.5, 0.5, message, cex = 0.8, col = "#555555")
    graphics::box(col = "#B3B3B3")
}
.qc_plot_named_values <- function(values, rounds) {
    if (is.null(values)) {
        return(stats::setNames(rep(NA_real_, length(rounds)), rounds))
    }
    value_names <- names(values)
    values <- as.numeric(values)
    if (is.null(value_names) && length(values) == length(rounds)) {
        value_names <- rounds
    }
    result <- stats::setNames(rep(NA_real_, length(rounds)), rounds)
    if (!is.null(value_names)) {
        matched <- match(rounds, value_names)
        available <- !is.na(matched)
        result[available] <- values[matched[available]]
    }
    result
}
.qc_plot_provenance <- function(experiment, rounds) {
    if (is.null(experiment)) {
        return(list(input = NULL, output = NULL, report = NULL))
    }
    metadata <- S4Vectors::metadata(experiment)
    extraction <- metadata$extraction
    if (!inherits(extraction, "selexprep_extraction")) {
        extraction <- NULL
    }
    list(input = .qc_plot_named_values(extraction$input_reads, rounds), output = .qc_plot_named_values(extraction$output_reads,
        rounds), report = metadata$library_report)
}
.qc_plot_retention <- function(x, provenance, rounds) {
    complete <- !is.null(provenance$input) && all(is.finite(provenance$input)) &&
        all(is.finite(provenance$output))
    if (complete) {
        values <- rbind(input = provenance$input, extracted = provenance$output)
        graphics::barplot(values, beside = TRUE, names.arg = rounds, col = .qc_plot_colours[1:2],
            border = NA, las = 2, cex.names = 0.75, ylab = "Reads", main = "Read retention")
        graphics::legend("topright", legend = rownames(values), fill = .qc_plot_colours[1:2],
            border = NA, bty = "n", cex = 0.75)
        return(invisible(NULL))
    }
    graphics::barplot(as.numeric(x$per_round$n_reads), names.arg = rounds, col = .qc_plot_colours[[2L]],
        border = NA, las = 2, cex.names = 0.75, ylab = "Reads", main = "Counted reads")
    graphics::mtext("Input-read provenance unavailable", side = 3, line = 0.1,
        cex = 0.65, col = "#555555")
    invisible(NULL)
}
.qc_plot_primer_match <- function(provenance, rounds) {
    complete <- !is.null(provenance$input) && all(is.finite(provenance$input)) &&
        all(is.finite(provenance$output))
    if (!complete) {
        .qc_plot_unavailable("Primer-match proxy", "Extraction provenance\nis unavailable")
        return(invisible(NULL))
    }
    retention <- ifelse(provenance$input > 0, provenance$output/provenance$input,
        NA_real_)
    positions <- seq_along(rounds)
    graphics::plot(positions, retention, type = "b", pch = 16, xaxt = "n", ylim = c(0,
        1.05), xlab = "Round", ylab = "Fraction", main = "Primer-match / retention proxy",
        col = .qc_plot_colours[[1L]])
    graphics::axis(1, at = positions, labels = rounds, las = 2, cex.axis = 0.75)
    legend_labels <- "extracted / input"
    legend_colours <- .qc_plot_colours[[1L]]
    legend_lty <- 1
    legend_pch <- 16
    report <- provenance$report
    if (inherits(report, "selexprep_library_report")) {
        references <- c(`5' report match` = report$match_rate_5p, `3' report match` = report$match_rate_3p)
        keep <- is.finite(references)
        reference_colours <- .qc_plot_colours[2:3][keep]
        references <- references[keep]
        if (length(references)) {
            for (index in seq_along(references)) {
                graphics::abline(h = references[[index]], col = reference_colours[[index]],
                  lty = 2)
            }
            legend_labels <- c(legend_labels, names(references))
            legend_colours <- c(legend_colours, reference_colours)
            legend_lty <- c(legend_lty, rep(2, length(references)))
            legend_pch <- c(legend_pch, rep(NA_integer_, length(references)))
        }
    }
    graphics::legend("bottomleft", legend = legend_labels, col = legend_colours,
        lty = legend_lty, pch = legend_pch, bty = "n", cex = 0.7)
    invisible(NULL)
}
.qc_plot_length_distribution <- function(experiment, rounds) {
    if (is.null(experiment)) {
        .qc_plot_unavailable("Random-region lengths", "Count experiment\nis unavailable")
        return(invisible(NULL))
    }
    counts <- .counts_by_round(experiment)
    if (!all(rounds %in% names(counts))) {
        stop("QC and experiment round names do not agree.", call. = FALSE)
    }
    counts <- counts[rounds]
    distributions <- lapply(counts, function(round_counts) {
        if (!length(round_counts)) {
            return(numeric())
        }
        tapply(as.numeric(round_counts), nchar(names(round_counts)), sum)
    })
    lengths <- sort(unique(as.integer(unlist(lapply(distributions, names), use.names = FALSE))))
    if (!length(lengths)) {
        .qc_plot_unavailable("Random-region lengths", "No counted sequences")
        return(invisible(NULL))
    }
    fractions <- matrix(0, nrow = length(lengths), ncol = length(rounds), dimnames = list(lengths,
        rounds))
    for (index in seq_along(distributions)) {
        distribution <- distributions[[index]]
        if (!length(distribution)) {
            next
        }
        rows <- match(as.integer(names(distribution)), lengths)
        fractions[rows, index] <- as.numeric(distribution)
        total <- sum(fractions[, index])
        if (total > 0) {
            fractions[, index] <- fractions[, index]/total
        }
    }
    x_limits <- range(lengths)
    if (!diff(x_limits)) {
        x_limits <- x_limits + c(-0.5, 0.5)
    }
    colours <- rep(.qc_plot_colours, length.out = length(rounds))
    graphics::matplot(lengths, fractions, type = "o", lty = rep(1:6, length.out = length(rounds)),
        pch = rep(15:20, length.out = length(rounds)), col = colours, xlim = x_limits,
        xlab = "Random-region length", ylab = "Fraction of reads", main = "Length distribution")
    graphics::legend("topright", legend = rounds, col = colours, lty = rep(1:6,
        length.out = length(rounds)), pch = rep(15:20, length.out = length(rounds)),
        bty = "n", cex = 0.65, ncol = max(1, ceiling(length(rounds)/8)))
    invisible(NULL)
}
.qc_plot_diversity <- function(x, rounds) {
    positions <- seq_along(rounds)
    graphics::barplot(as.numeric(x$per_round$n_unique), names.arg = rounds, col = .qc_plot_colours[[3L]],
        border = NA, las = 2, cex.names = 0.75, ylab = "Unique sequences", main = "Unique count")
    graphics::plot(positions, as.numeric(x$per_round$shannon_entropy_bits), type = "b",
        pch = 16, xaxt = "n", xlab = "Round", ylab = "Entropy (bits)", main = "Shannon diversity",
        col = .qc_plot_colours[[1L]])
    graphics::axis(1, at = positions, labels = rounds, las = 2, cex.axis = 0.75)
    graphics::plot(positions, as.numeric(x$per_round$top_100_coverage), type = "b",
        pch = 16, xaxt = "n", ylim = c(0, 1.05), xlab = "Round", ylab = "Fraction of reads",
        main = "Top-100 coverage", col = .qc_plot_colours[[2L]])
    graphics::axis(1, at = positions, labels = rounds, las = 2, cex.axis = 0.75)
    invisible(NULL)
}
.qc_plot_layout <- function(n_panels) {
    if (n_panels <= 1L) {
        return(c(1L, 1L))
    }
    if (n_panels <= 3L) {
        return(c(1L, n_panels))
    }
    if (n_panels <= 4L) {
        return(c(2L, 2L))
    }
    c(2L, ceiling(n_panels/2L))
}
#' Plot SELEX quality-control diagnostics
#'
#' Renders an R-native diagnostic surface corresponding to the four QC
#' artifacts produced by the original Python workflow. Read retention and its
#' primer-match proxy use extraction provenance from `y`. Sequence lengths are
#' shown as read-weighted, per-round fractions so rounds with different depths
#' remain comparable. The diversity view contains unique counts, Shannon
#' entropy, and top-100 coverage.
#'
#' When `y` is omitted, summaries already retained in `x` are still plotted.
#' Panels that require raw counts or extraction provenance state that the data
#' are unavailable instead of reconstructing them from incomplete summaries.
#'
#' @param x A `selexprep_qc` object returned by `selexprep_qc()`.
#' @param y Optional `SummarizedExperiment` used to compute `x`. It supplies
#'   count distributions and extraction provenance.
#' @param which One or more diagnostic groups: `'retention'`,
#'   `'primer_match'`, `'length'`, or `'diversity'`.
#' @param ... Named graphical parameters passed to [graphics::par()].
#'
#' @return `x`, invisibly.
#' @export
#' @examples
#' p5 <- 'GGTAATACGACTCACTATAGGG'
#' p3 <- 'CCATGCATGCATGCATGCAT'
#' pools <- list(round_00 = rep(paste0(p5, 'ACGTACGT', p3), 500))
#' result <- run_selexprep(pools, low_total_reads = 0)
#' plot(S4Vectors::metadata(result)$qc, result, which = 'diversity')
plot.selexprep_qc <- function(x, y = NULL, which = c("retention", "primer_match",
    "length", "diversity"), ...) {
    required <- c("round", "n_reads", "n_unique", "shannon_entropy_bits", "top_100_coverage")
    if (!inherits(x, "selexprep_qc") || is.null(x$per_round) || !all(required %in%
        colnames(x$per_round))) {
        stop("`x` must be a valid selexprep_qc object.", call. = FALSE)
    }
    if (!is.null(y) && !inherits(y, "SummarizedExperiment")) {
        stop("`y` must be NULL or a SummarizedExperiment.", call. = FALSE)
    }
    choices <- c("retention", "primer_match", "length", "diversity")
    which <- match.arg(which, choices, several.ok = TRUE)
    rounds <- as.character(x$per_round$round)
    if (!length(rounds)) {
        stop("`x` does not contain any rounds to plot.", call. = FALSE)
    }
    panel_count <- sum(which != "diversity") + 3L * sum(which == "diversity")
    dimensions <- .qc_plot_layout(panel_count)
    old_parameters <- graphics::par(no.readonly = TRUE)
    on.exit(suppressWarnings(graphics::par(old_parameters)), add = TRUE)
    graphics::par(mfrow = dimensions, mar = c(6, 4, 3, 1) + 0.1, mgp = c(2.3, 0.7,
        0))
    graphical_parameters <- list(...)
    if (length(graphical_parameters)) {
        do.call(graphics::par, graphical_parameters)
    }
    provenance <- .qc_plot_provenance(y, rounds)
    for (diagnostic in which) {
        switch(diagnostic, retention = .qc_plot_retention(x, provenance, rounds),
            primer_match = .qc_plot_primer_match(provenance, rounds), length = .qc_plot_length_distribution(y,
                rounds), diversity = .qc_plot_diversity(x, rounds))
    }
    invisible(x)
}
