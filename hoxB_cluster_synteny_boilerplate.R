library(readr)
library(dplyr)
library(SVbyEye)
library(pafr)
library(stringr)
library(GenomicRanges)

##### ------------------------------------
##### Set working directory
##### ------------------------------------

setwd("/C:/Documents/YT_synteny_hoxb")

##### ------------------------------------
##### Load all six PAF files
##### ------------------------------------

# Target: Human | Query: Mouse
dog_human_paf <- readPaf(
  paf.file = "./paf/<>.paf",
  include.paf.tags = TRUE, restrict.paf.tags = "cg"
)
dog_human_paf$q.name <- "Dog_chr16"
dog_human_paf$t.name <- "Human_chr7"

##### ------------------------------------
##### Load all seven custom GTF files to specify where to place rectangles for each KZFP
##### ------------------------------------

# Human hg38 annotation
humanHg38_anno <- read_tsv("./gtf/<>.gtf",
                           comment = "#", col_names = FALSE)

# Dog canFam3 annotation
dogCanFam3_anno <- read_tsv("./gtf/<>.gtf",
                            comment = "#", col_names = FALSE)

##### ------------------------------------
##### Format annotations for plot
##### ------------------------------------
 
# Human Annotation
humanHg38_chr7_anno <- humanHg38_anno %>%
  dplyr::mutate(gene_name = str_extract(X9, '(?<=gene_name ")[^"]+')) %>%
  dplyr::filter(., X1=='chr7'
                & X4>149059641
                & X5<149507809 ) %>%
  mutate(., chr='Human_chr7')%>%
  mutate(., start=X4-149059641) %>%
  mutate(., end=X5-149059641) %>%
  select(., chr, start, end, X7, gene_name) %>%
  arrange(., start) %>%
  dplyr::rename(., strand = X7)

# Dog Annotation
dogCanFam3_chr16_anno <- dogCanFam3_anno %>%
  dplyr::mutate(gene_name = str_extract(X9, '(?<=gene_name ")[^"]+')) %>%
  dplyr::filter(., X1=='chr16'
                & X4>14134497
                & X5<14436586 ) %>%
  mutate(., chr='Dog_chr16')%>%
  mutate(., start=X4-14134497) %>%
  mutate(., end=X5-14134497) %>%
  select(., chr, start, end, X7, gene_name) %>%
  arrange(., start) %>%
  dplyr::rename(., strand = X7)

##### ------------------------------------
##### Bind all dataframes in order
##### ------------------------------------

<>_list <- list(dog_human_paf = dog_human_paf)
<>Pafs <- bind_rows(<>_list)

<>_anno <- dplyr::bind_rows(humanHg38_chr7_anno, dogCanFam3_chr16_anno)

# Retain gene_name metadata
<>_anno_gr <- GenomicRanges::makeGRangesFromDataFrame(
  <>_anno, 
  seqnames.field = "chr",
  start.field    = "start",
  end.field      = "end",
  strand.field   = "strand",
  keep.extra.columns = TRUE   # <- key
)

# Assign the assembly order in the plot
seqnames.order <- c("Human_chr7", "Dog_chr16")

# Check that worked
head(<>Pafs)

# Assemble base plot
<>dir <- plotAVA(
  paf.table = <>Pafs,
  color.by = "direction",
  # color.by = "identity",
  color.palette = c("+" = "#F5DEB3", "-" = "#C4C4C2"),
  seqnames.order = seqnames.order
)

# Add rectangle KZFP gene annotations, and optionally gene_name labels
<>dir_anno <- addAnnotation(
  ggplot.obj = <>dir, 
  annot.gr = <>_anno_gr, 
  shape = "rectangle",
  coordinate.space = 'self', 
  y.label.id = 'seqnames', 
  # label.by = 'gene_name',
  annotation.level = 0
)

# Check plot
<>dir_anno

# Inspect layers to adjust gene_name label size as needed
<>dir_anno$layers
# 
# <>dir_anno$layers[[4]]$aes_params$size <- 1.5
# <>dir_anno
# 
# <>dir_anno$layers[[4]]$mapping$y <-
#   rlang::expr(y.offset + 0)
# 
# <>dir_anno


##### ------------------------------------
##### Save plot
##### ------------------------------------

# Save pdf of plot
pdf("./plots/<>.pdf", width = 10, height = 5)
print(<>dir_anno)
dev.off()

# Add title and other adjustments here
plotWithTitle <- <>dir_anno + labs(
  title = "<>",
  x = "<> Locus Size (bp)",
  # y = "Sequence Blocks",
  color = "Strand Direction"
)

# Save with adjustments
pdf("./plots/<>.pdf", width = 10, height = 5)
print(plotWithTitle)
dev.off()

##### ------------------------------------
##### View final plot
##### ------------------------------------

plotWithTitle