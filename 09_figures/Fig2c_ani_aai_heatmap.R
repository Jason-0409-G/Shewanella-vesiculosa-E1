#!/usr/bin/env Rscript
# Fig. 2c - combined ANI / AAI heatmap of the eight Shewanella genomes.
# Upper triangle = ANI (FastANI), lower triangle = AAI (EzAAI), diagonal = 100.
#
# Usage:  Rscript Fig2c_ani_aai_heatmap.R [ANI_8x8.tsv] [AAI_8x8.tsv]
#   ANI_8x8.tsv  FastANI long format (query, reference, ANI, mapped, total fragments)
#   AAI_8x8.tsv  EzAAI long format with columns "Label 1", "Label 2", "AAI"
#   (both written by 03_comparative_genomics/ani_aai_8genomes.sh)
# Output: Fig2c_ani_aai_heatmap.svg
#
# The ANI and AAI triangles have separate colour scales (ANI 83-100, AAI 88-100)
# and legends. The published panel uses a single merged colour bar
# ("Identity (%)", 85-100), applied afterwards in the vector editor.

library(ComplexHeatmap)
library(circlize)
library(svglite)
library(grid)

args    <- commandArgs(trailingOnly = TRUE)
ANI_TSV <- if (length(args) >= 1) args[1] else "ANI_8x8.tsv"
AAI_TSV <- if (length(args) >= 2) args[2] else "AAI_8x8.tsv"

# strain order and display names (L5 is the internal name of E1)
STRAINS <- c("L5", "M7", "frig", "sp02", "sp14", "pola", "psyc", "livi")

STRAIN_LABELS <- c(
  L5   = "S. vesiculosa E1",
  M7   = "S. vesiculosa M7",
  frig = "S. frigidimarina",
  sp02 = "S. sp002836315",
  sp14 = "S. sp014164505",
  pola = "S. polaris",
  psyc = "S. psychromarinicola",
  livi = "S. livingstonensis"
)

# ANI
ani_raw <- read.table(
  ANI_TSV,
  header = FALSE, sep = "\t",
  col.names = c("query", "ref", "ani", "map", "total")
)
extract_strain <- function(x) sub(".*/(.*)\\.fna", "\\1", x)
ani_raw$q <- extract_strain(ani_raw$query)
ani_raw$r <- extract_strain(ani_raw$ref)

ani_mat <- matrix(NA_real_, 8, 8, dimnames = list(STRAINS, STRAINS))
for (k in seq_len(nrow(ani_raw))) {
  q <- ani_raw$q[k]; r <- ani_raw$r[k]
  if (q %in% STRAINS && r %in% STRAINS)
    ani_mat[q, r] <- round(ani_raw$ani[k], 1)
}

# AAI
aai_raw <- read.table(
  AAI_TSV,
  header = TRUE, sep = "\t", check.names = FALSE
)
colnames(aai_raw) <- make.names(colnames(aai_raw))  # column names contain spaces

aai_mat <- matrix(NA_real_, 8, 8, dimnames = list(STRAINS, STRAINS))
for (k in seq_len(nrow(aai_raw))) {
  q <- as.character(aai_raw$Label.1[k])
  r <- as.character(aai_raw$Label.2[k])
  if (q %in% STRAINS && r %in% STRAINS)
    aai_mat[q, r] <- round(aai_raw$AAI[k], 1)
}

# Combined matrix
n <- length(STRAINS)
comb_mat <- matrix(NA_real_, n, n, dimnames = list(STRAINS, STRAINS))
for (i in seq_len(n)) for (j in seq_len(n)) {
  if (i == j)      comb_mat[i, j] <- 100
  else if (i < j)  comb_mat[i, j] <- ani_mat[STRAINS[i], STRAINS[j]]  # ANI
  else             comb_mat[i, j] <- aai_mat[STRAINS[i], STRAINS[j]]  # AAI
}

# Colour scales; ANI ranges from 83.7 to 100
ani_col <- colorRamp2(
  c(83, 87, 92, 97, 100),
  c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B")
)
# AAI ranges from 88.6 to 100
aai_col <- colorRamp2(
  c(88, 91, 94, 97, 100),
  c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B")
)

# Groups
GROUPS <- c(L5   = "vesiculosa",    M7   = "vesiculosa",
            frig = "frigidimarina", sp02 = "frigidimarina", sp14 = "frigidimarina",
            pola = "other",         psyc = "other",         livi = "other")
grp_fac <- factor(GROUPS[STRAINS], levels = c("vesiculosa", "frigidimarina", "other"))

# Heatmap
ht <- Heatmap(
  comb_mat,
  name = "ANI_AAI",

  # placeholder; cell colours are drawn in cell_fun
  col = colorRamp2(c(83, 100), c("#F7FBFF", "#08306B")),

  cluster_rows    = FALSE,
  cluster_columns = FALSE,

  # row labels (E1 in red bold italic)
  row_labels   = STRAIN_LABELS[STRAINS],
  row_names_gp = gpar(
    fontsize = 8,
    fontface = ifelse(STRAINS == "L5", "bold.italic", "italic"),
    col      = ifelse(STRAINS == "L5", "#C0392B", "#2C2C2C")
  ),
  row_names_side      = "left",
  row_names_max_width = unit(4.5, "cm"),

  # column labels
  column_labels = STRAIN_LABELS[STRAINS],
  column_names_gp = gpar(
    fontsize = 8,
    fontface = ifelse(STRAINS == "L5", "bold.italic", "italic"),
    col      = ifelse(STRAINS == "L5", "#C0392B", "#2C2C2C")
  ),
  column_names_rot        = 45,
  column_names_side       = "bottom",
  column_names_max_height = unit(4.5, "cm"),

  row_title    = NULL,
  column_title = NULL,

  # upper triangle (incl. diagonal): ANI colours; lower triangle: AAI colours
  cell_fun = function(j, i, x, y, width, height, fill) {
    val <- comb_mat[i, j]
    if (is.na(val)) return()

    bg <- if (i <= j) ani_col(val) else aai_col(val)

    grid.rect(x, y, width, height,
              gp = gpar(fill = bg, col = "#FFFFFF", lwd = 0.6))

    # white text on dark cells
    txt_col <- if (val >= 97) "#FFFFFF" else "#0D2233"
    fw      <- if (i == j) "bold" else "plain"
    grid.text(sprintf("%.1f", val), x, y,
              gp = gpar(fontsize = 6.5, col = txt_col, fontface = fw))
  },

  rect_gp = gpar(col = NA),

  width  = unit(n * 11, "mm"),
  height = unit(n * 11, "mm"),

  show_heatmap_legend = FALSE,  # custom legends below
  border = FALSE
)

# Legends
ani_lgd <- Legend(
  col_fun       = ani_col,
  title         = "ANI (%)",
  title_gp      = gpar(fontsize = 7.5, fontface = "bold", col = "#08306B"),
  labels_gp     = gpar(fontsize = 6.5),
  at            = c(84, 88, 92, 96, 100),
  legend_height = unit(2.5, "cm"),
  border        = "#BBBBBB",
  grid_width    = unit(3.5, "mm")
)

aai_lgd <- Legend(
  col_fun       = aai_col,
  title         = "AAI (%)",
  title_gp      = gpar(fontsize = 7.5, fontface = "bold", col = "#08306B"),
  labels_gp     = gpar(fontsize = 6.5),
  at            = c(88, 91, 94, 97, 100),
  legend_height = unit(2.5, "cm"),
  border        = "#BBBBBB",
  grid_width    = unit(3.5, "mm")
)

# Output
out_svg <- "Fig2c_ani_aai_heatmap.svg"

svglite(out_svg, width = 9.5, height = 7.5,
        system_fonts = list(sans = "Arial"),
        fix_text_size = FALSE)

draw(ht,
     padding             = unit(c(8, 4, 6, 4), "mm"),
     heatmap_legend_list = list(ani_lgd, aai_lgd),
     heatmap_legend_side = "right",
     merge_legend        = FALSE)

# triangle labels
decorate_heatmap_body("ANI_AAI", {
  grid.text("ANI",
            x    = unit(1, "npc"),
            y    = unit(1, "npc"),
            just = c("right", "top"),
            gp   = gpar(fontsize = 9, fontface = "bold", col = "#08306B"))
  grid.text("AAI",
            x    = unit(0, "npc"),
            y    = unit(0, "npc"),
            just = c("left", "bottom"),
            gp   = gpar(fontsize = 9, fontface = "bold", col = "#08306B"))
})

dev.off()

cat(sprintf("written: %s\n", out_svg))
