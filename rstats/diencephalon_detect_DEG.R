library(conflicted)
library(tidyverse)
library(TCC)


## DEG detection ---------------------------------------------------------------

group = c(1,1,1,1,1,1,2,2,2,2,2,2)
DE_method = "edger"
FDR = 0.05

count_data = readr::read_csv("out/diencephalon_gene_count_matrix.csv") |>
  tidyr::drop_na()
count_matrix = count_data |>
  tibble::column_to_rownames(var = "gene_id") |>
  as.matrix()

tcc = new("TCC", count_matrix, group)
tcc = TCC::calcNormFactors(tcc, iteration = TRUE)
tcc = TCC::estimateDE(tcc, test.method = DE_method, FDR = FDR)
result = TCC::getResult(tcc, sort = FALSE) |> dplyr::arrange(rank)

## Write and re-read TCC result

readr::write_tsv(result, "out/diencephalon_stringtie_TCC.tsv")
sorted = readr::read_tsv("out/diencephalon_stringtie_TCC.tsv") |>
  dplyr::mutate(definedDEG = dplyr::case_when(
    q.value < 0.05 & m.value > 1 ~ "WRM",
    q.value < 0.05 & m.value < -1 ~ "BF",
    .default = "N"
  )) |>
  dplyr::arrange(desc(rank))

sorted |> dplyr::filter(definedDEG == "BF") |> nrow()
sorted |> dplyr::filter(definedDEG == "WRM") |> nrow()
# sorted |>
#   dplyr::filter(definedDEG != "N") |>
#   dplyr::mutate(gene_id = stringr::str_split(gene_id, "\\|", simplify = TRUE)[,2]) |>
#   dplyr::filter(stringr::str_detect(gene_id, "^LOC", negate = TRUE)) |>
#   dplyr::select(!(dplyr::starts_with(c("B0", "W0", "estimatedDEG")))) |>
#   dplyr::arrange(rank) |>
#   readr::write_excel_csv("docs/rnaseq_deg_list.csv")


## MA-plot

pma = ggplot(sorted) +
  aes(x = a.value, y = m.value) +
  geom_point(aes(color = definedDEG), size = 2) +
  scale_color_manual(
    values = c("BF" = "#d8b365", "WRM" = "#5ab4ac", "N" = "grey"),
    labels = c("BF" = "Highly expressed in BF", "WRM" = "Highly expressed in WRM", "N" = "No significant")
  ) +
  labs(
    x = expression(paste({log[2]}, "Mean Expression")),
    y = expression(paste({log[2]}, "Fold Change"))
  ) +
  theme_test(base_size = 12) +
  theme(
    legend.position = "inside",
    legend.justification = c(.99, .01),
    legend.background = element_rect(color = NA),
    legend.title = element_blank()
  )
pma

## Volcano plot

pvol = ggplot(sorted) +
  aes(x = m.value, y = -log10(q.value)) +
  geom_point(aes(color = definedDEG), size = 2) +
  scale_color_manual(
    values = c("BF" = "#d8b365", "WRM" = "#5ab4ac", "N" = "grey"),
    labels = c("BF" = "Highly expressed in BF", "WRM" = "Highly expressed in WRM", "N" = "No significant")
  ) +
  labs(
    x = expression(paste({log[2]}, "(Fold Change)")),
    y = expression(paste({-log[10]}, "(q-value)"))
  ) +
  theme_test(base_size = 12) +
  theme(
    legend.position = "inside",
    legend.justification = c(.99, .99),
    legend.title = element_blank()
  )
pvol
