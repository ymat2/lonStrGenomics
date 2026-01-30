lonStrDom2_gene_list = readr::read_tsv("out/lonStrDom2_gene_list.tsv")
annot = function(df, seqnames, start, end, ...) {
  gr = GenomicRanges::makeGRangesFromDataFrame(
    df,
    keep.extra.columns = TRUE,
    ignore.strand = TRUE,
    seqnames.field = seqnames,
    start.field = start,
    end.field = end
  )
  
  gene_list = lonStrDom2_gene_list |>
    tidyr::drop_na(Symbol) |>
    dplyr::rename(gene = Symbol) |>
    GenomicRanges::makeGRangesFromDataFrame(
      keep.extra.columns = TRUE,
      ignore.strand = TRUE,
      seqnames.field = "Accession",
      start.field = "Begin",
      end.field = "End"
    )
  
  df_overlap = IRanges::mergeByOverlaps(
    gr,
    gene_list,
    maxgap = -1L,
    minoverlap = 0L,
    type = "any",
    select = "all"
  ) |>
    dplyr::as_tibble() |>
    dplyr::select(
      "gr.seqnames",
      "gr.start",
      "gr.end",
      "gene"
    ) |>
    dplyr::rename(
      !!seqnames := "gr.seqnames",
      !!start := "gr.start",
      !!end := "gr.end"
    )
  
  df_merge = df |> dplyr::left_join(
    df_overlap,
    by = c(seqnames, start, end)
  )
  
  return(df_merge)
}