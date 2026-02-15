library(readr)
library(dplyr)
library(SVbyEye)
library(pafr)
library(stringr)
library(GenomicRanges)

##### ------------------------------------
##### Set working directory
##### ------------------------------------

# setwd("/C:/Documents/YT_synteny_hoCxb")

##### ------------------------------------
##### Load all six PAF files
##### ------------------------------------

# Target: Human | Query: Mouse
mouse_human_hoxb_paf <- readPaf(
  paf.file = "QUERY_mouse_mm39_chr11_hoxb_cluster_TARGET_human_hg38_chr17_hoxb_cluster_lastz.paf",
  include.paf.tags = TRUE, restrict.paf.tags = "cg"
)
mouse_human_hoxb_paf$q.name <- "Mouse_chr11"
mouse_human_hoxb_paf$t.name <- "Human_chr17"

##### ------------------------------------
##### Load all seven custom GTF files to specify where to place rectangles for each KZFP
##### ------------------------------------

# Human hg38 annotation
humanHg38_anno <- read_tsv("hg38_hoxb_cluster.gtf",
                           comment = "#", col_names = FALSE)

# Dog canFam3 annotation
mousemm39_anno <- read_tsv("mm39_hoxb_cluster.gtf",
                            comment = "#", col_names = FALSE)

##### ------------------------------------
##### Format annotations for plot
##### ------------------------------------

# Human Annotation
humanHg38_chr17_anno <- humanHg38_anno %>%
  dplyr::mutate(gene_name = str_extract(X9, '(?<=gene_name ")[^"]+')) %>%
  dplyr::filter(., X1=='chr17'
                & X4>48526807
                & X5<48627161 ) %>%
  mutate(., chr='Human_chr17')%>%
  mutate(., start=X4-48526807) %>%
  mutate(., end=X5-48526807) %>%
  select(., chr, start, end, X7, gene_name) %>%
  arrange(., start) %>%
  dplyr::rename(., strand = X7)

# Dog Annotation
mousemm39_chr11_anno <- mousemm39_anno %>%
  dplyr::mutate(gene_name = str_extract(X9, '(?<=gene_name ")[^"]+')) %>%
  dplyr::filter(., X1=='chr11'
                & X4>96159502
                & X5<96263131 ) %>%
  mutate(., chr='Mouse_chr11')%>%
  mutate(., start=X4-96159502) %>%
  mutate(., end=X5-96159502) %>%
  select(., chr, start, end, X7, gene_name) %>%
  arrange(., start) %>%
  dplyr::rename(., strand = X7)

##### ------------------------------------
##### Bind all dataframes in order
##### ------------------------------------

# hoxb_cluster_list <- list(mouse_human_hoxb_paf = mouse_human_hoxb_paf)
# hoxb_clusterPafs <- bind_rows(hoxb_cluster_list)

hoxb_cluster_anno <- dplyr::bind_rows(humanHg38_chr17_anno, mousemm39_chr11_anno)

# Retain gene_name metadata
hoxb_cluster_anno_gr <- GenomicRanges::makeGRangesFromDataFrame(
  hoxb_cluster_anno, 
  seqnames.field = "chr",
  start.field    = "start",
  end.field      = "end",
  strand.field   = "strand",
  keep.extra.columns = TRUE   # <- key
)

# Assign the assembly order in the plot
seqnames.order <- c("Human_chr17", "Mouse_chr11")

# Check that worked
head(mouse_human_hoxb_paf)

hoxb_clusterPafs <- mouse_human_hoxb_paf

# Assemble base plot
hoxb_clusterdir <- plotAVA(
  paf.table = hoxb_clusterPafs,
  color.by = "direction",
  # color.by = "identity",
  color.palette = c("+" = "#F5DEB3", "-" = "#0077b6"),
  seqnames.order = seqnames.order
)

hoxb_clusterdir


# Add rectangle KZFP gene annotations, and optionally gene_name labels
hoxb_clusterdir_anno <- addAnnotation(
  ggplot.obj = hoxb_clusterdir, 
  annot.gr = hoxb_cluster_anno_gr, 
  shape = "rectangle",
  coordinate.space = 'self', 
  y.label.id = 'seqnames', 
  label.by = 'gene_name',
  annotation.level = 0
)

# Check plot
hoxb_clusterdir_anno

# Inspect layers to adjust gene_name label size as needed
hoxb_clusterdir_anno$layers
# 
hoxb_clusterdir_anno$layers[[4]]$aes_params$size <- 2.6
hoxb_clusterdir_anno
# 
# hoxb_clusterdir_anno$layers[[4]]$mapping$y <-
#   rlang::expr(y.offset + 0)
# 
# hoxb_clusterdir_anno


##### ------------------------------------
##### Save plot
##### ------------------------------------

# Save pdf of plot
pdf("hoxb_cluster_hg38_mm39.pdf", width = 10, height = 6)
print(hoxb_clusterdir_anno)
dev.off()

# Add title and other adjustments here
plotWithTitle <- hoxb_clusterdir_anno + labs(
  title = "HoxB Cluster Synteny in Human (hg38) and Mouse (mm39)",
  x = "HoxB Cluster Locus Size (bp)",
  # y = "Sequence Blocks",
  # color = "Strand Direction"
)

plotWithTitle

# Save with adjustments
pdf("hoxb_cluster_hg38_mm39_withTitle.pdf", width = 10, height = 7)
print(plotWithTitle)
dev.off()

##### ------------------------------------
##### View final plot
##### ------------------------------------

plotWithTitle