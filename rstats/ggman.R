require(conflicted)
require(tidyverse)
# require(ggrepel)


ggmanh = function(
    x, chr = "CHR", bp = "POS", p = "P",
    col = c("#1f78b4", "#a6cee3"),
    chrlabs = c(1:10, 12, 14, 17, 20, 24, 36),
    logp = TRUE,
    xlab = "Chromosome",
    ylab = "P",
    title = NULL,
    plot = TRUE,
    ylim = NULL,
    base_size = 11,
    ...
  ) {
    columns <- c(chr, bp, p)

    if (!all(columns %in% colnames(x))) {
      missing_columns <- columns[!columns %in% colnames(x)]
      stop(paste("Error: The following columns are missing:", paste(missing_columns, collapse = ", ")))
    }

    tbl = x |>
      dplyr::select(all_of(columns)) |>
      dplyr::rename(chr = !!sym(chr), bp = !!sym(bp), p = !!sym(p)) |>
      dplyr::as_tibble()

    data_cum = tbl |>
      dplyr::group_by(chr) |>
      dplyr::summarise(max_bp = max(bp)) |>
      dplyr::mutate(bp_add = dplyr::lag(cumsum(max_bp), default = 0)) |>
      dplyr::select(chr, bp_add)

    tbl4plot = tbl |>
      dplyr::inner_join(data_cum, by = "chr") |>
      dplyr::mutate(bp_cum = bp + bp_add)

    axis_set = tbl4plot |>
      dplyr::group_by(chr) |>
      dplyr::summarize(center = mean(bp_cum)) |>
      dplyr::mutate(chr = dplyr::if_else(chr %in% chrlabs, as.character(chr), ""))

    if (is.null(ylim)) {
      ymin = 0
      ymax = max(tbl4plot$p)*1.1
      ylim = c(ymin, ymax)
    }

    g = ggplot(tbl4plot) +
      aes(x = bp_cum, y = p, color = forcats::as_factor(chr)) +
      geom_point(...) +
      scale_x_continuous(expand = c(.02, .02), label = axis_set$chr, breaks = axis_set$center) +
      scale_y_continuous(expand = c(.02, .02), limits = ylim) +
      scale_color_manual(values = rep(col, unique(length(axis_set$chr)))) +
      labs(x = xlab, y = ylab, title = title) +
      theme_classic(base_size = base_size) +
      theme(
        legend.position = "none",
        panel.grid = element_blank(),
        axis.line.x = element_blank(),
        plot.title = element_text(hjust = 0.5)
      )

    if (plot) {
      print(g)
    } else {
      g
    }
}
