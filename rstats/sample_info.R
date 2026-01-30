library(conflicted)
library(tidyverse)
library(readxl)


sample_20231227 = readxl::read_excel("サンプル情報/231227_DNA-EN_サンプル詳細追記.xlsx", sheet = "工作表1") |>
  dplyr::rename(sample_id = `* Sample Name`, sex = `性別`) |>
  dplyr::mutate(species = dplyr::if_else(`* Species` == "Lonchura striata", "white-rumped munia", "Bengalese finch")) |>
  dplyr::select(sample_id, sex, species)

sample_20240828 = readxl::read_excel("docs/20240828_DNA送付/20240828_麻布大_サンプルリスト.xlsx", sheet = "DNA_sample_list") |>
  dplyr::rename(sample_id = ID, species = `種`) |> 
  dplyr::mutate(sex = dplyr::case_when(
    `性別` %in% c("F", "M") ~ `性別`,
    is.na(`性別`) ~ stringr::str_extract(`麻布大学 PCRの性判定`, "F|M"),
    .default = NA
  )) |>
  dplyr::select(sample_id, sex, species)

sample_info = dplyr::bind_rows(sample_20231227, sample_20240828) |>
  dplyr::mutate(species = stringr::str_to_sentence(species)) |>
  dplyr::mutate(colony = dplyr::case_when(
    stringr::str_detect(sample_id, "WRM|SBM") ~ "Wild",
    stringr::str_detect(sample_id, "TB|WP") ~ "Tomakomai",
    stringr::str_detect(sample_id, "TO") ~ "Okanoya Lab",
    stringr::str_detect(sample_id, "BF66|BF67|AMBF") ~ "Others",
    .default = "Azabu University"
  )) |>
  dplyr::filter(species != "Scaly-breasted munia")
readr::write_tsv(sample_info, "out/sample_information_excel.tsv")
