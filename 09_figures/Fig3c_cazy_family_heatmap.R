#!/usr/bin/env Rscript
# Fig. 3c - CAZy family counts across the eight Shewanella genomes.
# Rows = genomes, columns = CAZy families grouped by predicted primary substrate.
#
# Usage:  Rscript Fig3c_cazy_family_heatmap.R [all_CAZy_family_matrix.tsv]
#   Input: family x genome count matrix written by
#          06_cazy_kegg_transporters/build_cazy_consensus_matrix.py
# Output: Fig3c_cazy_family_heatmap.svg
#
# The four cellulase families GH6, GH7, GH9 and GH12 are absent from all eight
# genomes and therefore missing from the matrix; they are added as zero columns so
# that the lack of cellulases is visible. The GH1 and GH3 highlight, the "Focal
# CAZymes" bracket and the group labels of the published panel were added
# afterwards in the vector editor.

library(ComplexHeatmap)
library(circlize)
library(svglite)
library(grid)

args <- commandArgs(trailingOnly = TRUE)
DATA <- if (length(args) >= 1) args[1] else "all_CAZy_family_matrix.tsv"

# read data
raw <- read.table(DATA, header = TRUE, sep = "\t",
                  row.names = 1, check.names = FALSE)
# raw: family x genome; transpose so that genomes are rows
mat_t <- t(as.matrix(raw))
class(mat_t) <- "numeric"

# add the all-zero cellulase families
.cel_missing <- setdiff(c("GH6","GH7","GH9","GH12"), colnames(mat_t))
if (length(.cel_missing) > 0) {
  .zero <- matrix(0, nrow = nrow(mat_t), ncol = length(.cel_missing),
                  dimnames = list(rownames(mat_t), .cel_missing))
  mat_t <- cbind(mat_t, .zero)
}

# row order (genomes); L5_PB002 is E1
ROW_ORDER <- c("L5_PB002","M7_REF","S_frigidimarina","S_livingstonensis",
               "S_polaris","S_psychromarinicola","S_sp002836315","S_sp014164505")
mat_t <- mat_t[ROW_ORDER, ]

ROW_LABELS <- c(
  L5_PB002             = "S. vesiculosa E1",
  M7_REF               = "S. vesiculosa M7",
  S_frigidimarina      = "S. frigidimarina",
  S_livingstonensis    = "S. livingstonensis",
  S_polaris            = "S. polaris",
  S_psychromarinicola  = "S. psychromarinicola",
  S_sp002836315        = "S.sp002836315",
  S_sp014164505        = "S.sp014164505"
)

# Families grouped by predicted primary substrate (CAZy / dbCAN)
SUBSTRATE <- list(
  "Cellulose"     = c("GH5","GH6","GH7","GH9","GH12"),         # cellulases
  "β-glucan"      = c("GH1","GH3","GH16"),                     # GH16 is assigned here
  "Chitin"        = c("AA10","GH18","GH19","GH20"),
  "α-glucan"      = c("GH13","GH15","GH57","GH77","GH97"),
  "Peptidoglycan" = c("GH23","GH24","GH73","GH103"),
  "Hemicellulose" = c("GH10","GH43"),                          # xylan, arabinan
  "Other glycans" = c("GH32","GH36","GH37","GH42","GH109","GH186"),
  "Auxiliary/GT"  = c("AA1","AA2","AA3","AA7","GT4","GT26")
)
SUB_LEVELS <- names(SUBSTRATE)
fam2sub <- unlist(lapply(SUB_LEVELS, function(s)
  setNames(rep(s, length(SUBSTRATE[[s]])), SUBSTRATE[[s]])))

# column order = substrate groups, families present in the data only
fam_order <- unname(unlist(SUBSTRATE))
fam_order <- fam_order[fam_order %in% colnames(mat_t)]
mat_t <- mat_t[, fam_order]

# substrate of each column (annotation bar and column split)
fam_sub <- factor(fam2sub[colnames(mat_t)], levels = SUB_LEVELS)

# Substrate annotation colours
SUB_COLS <- c(
  "Cellulose"     = "#B0553E",
  "β-glucan"      = "#3C7A6A",
  "Chitin"        = "#2B6CB0",
  "α-glucan"      = "#E0913A",
  "Peptidoglycan" = "#8E6CAE",
  "Hemicellulose" = "#C7A33E",
  "Other glycans" = "#7FA8C9",
  "Auxiliary/GT"  = "#AEB4B8"
)

col_anno <- HeatmapAnnotation(
  Substrate = fam_sub,
  col       = list(Substrate = SUB_COLS),
  annotation_legend_param = list(
    Substrate = list(
      title     = "Predicted substrate",
      title_gp  = gpar(fontsize = 7, fontface = "bold"),
      labels_gp = gpar(fontsize = 6.5)
    )
  ),
  simple_anno_size = unit(3.5, "mm"),
  show_annotation_name = FALSE,
  border = FALSE
)

# Colour scale
vmax <- max(mat_t)
col_fun <- colorRamp2(
  c(0, 1, 4, 8, vmax),
  c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B")
)

# Heatmap
ht <- Heatmap(
  mat_t,
  name = "Count",
  col  = col_fun,

  cluster_rows    = FALSE,
  cluster_columns = FALSE,

  # row labels
  row_labels   = ROW_LABELS[ROW_ORDER],
  row_names_gp = gpar(
    fontsize = 7.5,
    fontface = ifelse(ROW_ORDER %in% c("L5_PB002","M7_REF"), "bold.italic", "italic"),
    col      = ifelse(ROW_ORDER == "L5_PB002", "#C0392B", "#2C2C2C")
  ),
  row_names_side      = "left",
  row_names_max_width = unit(6.5, "cm"),

  # column labels
  column_labels          = colnames(mat_t),
  column_names_gp        = gpar(fontsize = 6, col = "#2C2C2C"),
  column_names_rot        = 90,
  column_names_side       = "bottom",
  column_names_max_height = unit(3.5, "cm"),

  # substrate bar above the columns
  top_annotation = col_anno,

  # columns split by substrate
  column_split = fam_sub,
  column_gap   = unit(1.5, "mm"),
  column_title = NULL,

  # counts in cells; zeros are not printed
  cell_fun = function(j, i, x, y, width, height, fill) {
    val <- mat_t[i, j]
    if (val == 0) return()
    txt_col <- if (val >= 8) "#FFFFFF" else "#0D2233"
    grid.text(as.character(val), x, y,
              gp = gpar(fontsize = 5.5, col = txt_col))
  },

  rect_gp = gpar(col = "#E0EAF4", lwd = 0.4),

  width  = unit(ncol(mat_t) * 5.5, "mm"),
  height = unit(nrow(mat_t) * 9, "mm"),

  heatmap_legend_param = list(
    title         = "Count",
    title_gp      = gpar(fontsize = 7, fontface = "bold"),
    labels_gp     = gpar(fontsize = 6.5),
    legend_height = unit(4, "cm"),
    at            = c(0, 2, 4, 6, 8, 10, 12),
    border        = "#BBBBBB",
    grid_width    = unit(3.5, "mm")
  ),

  border = FALSE,
  show_heatmap_legend = TRUE
)

# Output
out_svg <- "Fig3c_cazy_family_heatmap.svg"

svglite(out_svg, width = 22, height = 6,
        system_fonts = list(sans = "Arial"),
        fix_text_size = FALSE)

draw(ht,
     padding                = unit(c(4, 4, 4, 4), "mm"),
     merge_legend           = FALSE,
     heatmap_legend_side    = "right",
     annotation_legend_side = "right")

dev.off()

cat(sprintf("written: %s\n", out_svg))
