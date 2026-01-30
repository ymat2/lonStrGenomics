#colors = c("WRM" = "#5ab4ac", "BF" = "#d8b365", "SBM" = "#999999")
colWRM = "#5ab4ac"
colBF = "#d8b365"
colors = c("WRM" = colWRM, "BF" = colBF)
lighter = function(col, alpha = "66") {stringr::str_flatten(c(col, alpha))}

BASESIZE = 7
LABELSIZE = 10

ann_field = c(
  "Allele",
  "Annotation",
  "Putative_impact",
  "Gene_Name",
  "Gene_ID",
  "Feature_type",
  "Feature_ID",
  "Transcript_biotype",
  "Rank_total",
  "HGVS_c",
  "HGVS_p",
  "cDNA_position_cDNA_length",
  "CDS_position_CDS_length",
  "Protein_position_Protein_length",
  "Distance_to_feature",
  "Log"
)
