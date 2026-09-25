# Data import and differential expression ----
# Packages and setup ----
library(rtracklayer)
library(tximport)
library(tidyverse)
library(DESeq2)
library(pheatmap)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggrepel)

set.seed(20250924)

# Create output directories ----
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# Import reference annotation ----
gtf <- rtracklayer::import("reference/Homo_sapiens.GRCh38.113.gtf")
txdb <- gtf[gtf$type == "transcript"]

# Create transcript-to-gene mapping ----
tx2gene <- data.frame(
  transcript = txdb$transcript_id,
  gene = txdb$gene_id)

tx2gene <- tx2gene[!duplicated(tx2gene$transcript), ]

# Define sample IDs ----
sample_ids <- c(
  "SRR1039508", "SRR1039509",
  "SRR1039512", "SRR1039513",
  "SRR1039516", "SRR1039517",
  "SRR1039520", "SRR1039521")

# Define Salmon quantification files ----
files <- file.path(
  "quant",
  paste0(sample_ids, "_quant"),
  "quant.sf")

names(files) <- sample_ids

# Import Salmon quantification ----
txi <- tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene,
  txOut = FALSE,
  countsFromAbundance = "no",
  ignoreTxVersion = TRUE)

# Define sample metadata ----
metadata <- data.frame(
  sample = sample_ids,
  donor = c(
    "N61311", "N61311",
    "N052611", "N052611",
    "N080611", "N080611",
    "N061011", "N061011"),
  treatment = c(
    "Untreated", "Dexamethasone",
    "Untreated", "Dexamethasone",
    "Untreated", "Dexamethasone",
    "Untreated", "Dexamethasone"),
  row.names = "sample")

metadata$donor <- factor(metadata$donor)
metadata$treatment <- factor(metadata$treatment, levels = c("Untreated", "Dexamethasone"))

# Verify metadata and quantification match ----
stopifnot(identical(colnames(txi$counts), rownames(metadata)))

# Create DESeq2 dataset ----
dds <- DESeqDataSetFromTximport(
  txi = txi,
  colData = metadata,
  design = ~ donor + treatment)

# Filter low-count genes ----
keep <- rowSums(counts(dds) >= 10) >= 4
dds <- dds[keep, ]

# Run DESeq2 ----
dds <- DESeq(dds)

# Extract differential expression results ----
res <- results(dds, contrast = c("treatment", "Dexamethasone", "Untreated"))

# Shrink log2 fold changes ----
res_shrunk <- lfcShrink(
  dds,
  coef = "treatment_Dexamethasone_vs_Untreated",
  type = "apeglm")

# Convert results to data frame ----
res_df <- as.data.frame(res_shrunk)
res_df$gene_id <- rownames(res_df)

# Remove genes without statistical results ----
res_df <- res_df[
  !is.na(res_df$padj) &
    !is.na(res_df$log2FoldChange),]

# Remove Ensembl version numbers ----
res_df$ensembl_id <- sub("\\..*$", "", res_df$gene_id)

# Annotate Ensembl IDs ----
gene_annotation <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = res_df$ensembl_id,
  keytype = "ENSEMBL",
  columns = c("ENSEMBL", "SYMBOL", "GENENAME", "ENTREZID"))

gene_annotation <- gene_annotation[!duplicated(gene_annotation$ENSEMBL), ]

# Add gene annotation to DESeq2 results ----
res_df <- dplyr::left_join(
  res_df,
  gene_annotation,
  by = c("ensembl_id" = "ENSEMBL"))

# Define differential expression categories ----
res_df$significance <- "Not significant"

res_df$significance[
  res_df$padj < 0.05 &
    res_df$log2FoldChange > 0] <- "Upregulated"

res_df$significance[
  res_df$padj < 0.05 &
    res_df$log2FoldChange < 0] <- "Downregulated"

# Calculate -log10 adjusted p-value ----
res_df$neg_log10_padj <- -log10(res_df$padj)

# Define gene labels ----
res_df$label <- ifelse(
  !is.na(res_df$SYMBOL) &
    res_df$SYMBOL != "",
  res_df$SYMBOL,
  res_df$ensembl_id)

# Save complete differential expression results ----
write.csv(res_df, "results/differential_expression_results.csv", 
          row.names = FALSE)

# Define stringent DEGs ----
deg_strict <- res_df %>%
  filter(padj < 0.05, abs(log2FoldChange) >= 2) %>%
  arrange(padj)

# Save stringent DEG results ----
write.csv(deg_strict, "results/stringent_DEGs.csv", row.names = FALSE)

# Variance-stabilizing transformation ----
vsd <- vst(dds, blind = FALSE)

# PCA ----
pca_data <- plotPCA(
  vsd,
  intgroup = c("donor", "treatment"),
  returnData = TRUE)

percent_var <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(
  pca_data,
  aes(
    x = PC1,
    y = PC2,
    color = treatment,
    shape = donor,
    label = donor)) +
  geom_point(size = 4) +
  ggrepel::geom_text_repel(
    size = 3.5,
    box.padding = 0.5,
    point.padding = 0.3,
    max.overlaps = Inf) +
  labs(
    title = "PCA of primary human airway smooth muscle cells",
    subtitle = "Paired design: untreated versus dexamethasone treatment",
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance"),
    color = "Treatment",
    shape = "Donor") +
  theme_classic()

ggsave(
  "figures/Figure1_PCA.png",
  plot = pca_plot,
  width = 10,
  height = 8,
  dpi = 600)

# Sample correlation ----
sample_cor <- cor(assay(vsd), method = "pearson")

sample_labels <- paste(
  metadata$donor,
  metadata$treatment,
  sep = "_")

rownames(sample_cor) <- sample_labels
colnames(sample_cor) <- sample_labels

pheatmap(
  sample_cor,
  main = "Sample correlation",
  border_color = NA,
  color = colorRampPalette(c("#2878B5", "#FFFFFF", "#D94841"))(100),
  filename = "figures/Figure2_Sample_correlation.png",
  width = 10,
  height = 8)

# Sample distance ----
sample_dist <- dist(t(assay(vsd)))
sample_dist_matrix <- as.matrix(sample_dist)

rownames(sample_dist_matrix) <- sample_labels
colnames(sample_dist_matrix) <- sample_labels

pheatmap(
  sample_dist_matrix,
  main = "Sample distance",
  border_color = NA,
  color = colorRampPalette(c("#2878B5", "#FFFFFF", "#D94841"))(100),
  filename = "figures/Figure3_Sample_distance.png",
  width = 10,
  height = 8)

# MA plot ----
png(
  "figures/Figure4_MA_plot.png",
  width = 10,
  height = 8,
  units = "in",
  res = 600)

plotMA(
  res,
  alpha = 0.05,
  main = "Dexamethasone vs Untreated")

# Volcano plot ----
volcano_labels <- res_df %>%
  filter(
    padj < 0.05,
    abs(log2FoldChange) >= 2,
    !is.na(SYMBOL),
    SYMBOL != "") %>%
  arrange(desc(log2FoldChange))

top_10_up <- volcano_labels %>%
  slice_head(n = 10)

top_10_down <- volcano_labels %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10)

volcano_labels <- dplyr::bind_rows(
  top_10_up,
  top_10_down)

x_limit <- max(abs(res_df$log2FoldChange), na.rm = TRUE)

volcano_plot <- ggplot(
  res_df,
  aes(
    x = log2FoldChange,
    y = neg_log10_padj)) +
  geom_point(
    aes(
      color = padj < 0.05 &
        abs(log2FoldChange) >= 2),
    alpha = 0.7,
    size = 1.8) +
  scale_color_manual(
    values = c(
      "FALSE" = "gray70",
      "TRUE" = "red"),
    labels = c(
      "FALSE" = "Not significant",
      "TRUE" = "Significant"),
    name = "Category") +
  geom_vline(
    xintercept = c(-2, 0, 2),
    linetype = "dashed") +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed") +
  ggrepel::geom_text_repel(
    data = volcano_labels,
    aes(label = SYMBOL),
    size = 3.5,
    box.padding = 0.6,
    point.padding = 0.4,
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.color = "black",
    segment.size = 0.4) +
  scale_x_continuous(limits = c(-x_limit, x_limit)) +
  labs(
    title = "Dexamethasone differential expression",
    x = "Shrunken log2 fold change",
    y = "-log10 adjusted p-value") +
  theme_classic()

ggsave(
  "figures/Figure5_Volcano.png",
  plot = volcano_plot,
  width = 10,
  height = 8,
  dpi = 600)

# Heatmaps and DEG tables ----

# Define sample order

sample_order <- c(
  "SRR1039508", "SRR1039512", "SRR1039516", "SRR1039520",
  "SRR1039509", "SRR1039513", "SRR1039517", "SRR1039521")

# Define treatment annotation
annotation_col <- data.frame(
  Treatment = metadata[sample_order, "treatment"])

rownames(annotation_col) <- sample_order

annotation_colors <- list(
  Treatment = c(
    "Untreated" = "grey80",
    "Dexamethasone" = "grey40"))

# Define heatmap colors

heatmap_colors <- colorRampPalette(
  c("#2878B5", "#FFFFFF", "#D94841"))(100)

# Define stringent upregulated DEGs
upregulated_degs <- res_df %>%
  filter(
    padj < 0.05,
    log2FoldChange >= 2) %>%
  arrange(desc(log2FoldChange)) %>%
  filter(!duplicated(SYMBOL)) %>%
  filter(!is.na(SYMBOL), SYMBOL != "")

# Define stringent downregulated DEGs ----
downregulated_degs <- res_df %>%
  filter(
    padj < 0.05,
    log2FoldChange <= -2) %>%
  arrange(log2FoldChange) %>%
  filter(!duplicated(SYMBOL)) %>%
  filter(!is.na(SYMBOL), SYMBOL != "")

# Save upregulated DEG table
write.csv(
  upregulated_degs,
  "results/upregulated_DEGs.csv",
  row.names = FALSE)

# Save downregulated DEG table
write.csv(
  downregulated_degs,
  "results/downregulated_DEGs.csv",
  row.names = FALSE)

# Select top 60 upregulated genes
top_upregulated <- upregulated_degs %>%
  slice_head(n = 60)

# Select top 60 downregulated genes
top_downregulated <- downregulated_degs %>%
  slice_head(n = 60)

# Create expression matrix for top upregulated genes
up_matrix <- assay(vsd)[
  top_upregulated$gene_id,
  sample_order,
  drop = FALSE]

rownames(up_matrix) <- top_upregulated$SYMBOL

# Create expression matrix for top downregulated genes
down_matrix <- assay(vsd)[
  top_downregulated$gene_id,
  sample_order,
  drop = FALSE]

rownames(down_matrix) <- top_downregulated$SYMBOL

# Plot top 60 upregulated genes
pheatmap(
  up_matrix,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  color = heatmap_colors,
  border_color = NA,
  fontsize_row = 7,
  fontsize_col = 9,
  main = "Top 60 upregulated genes",
  filename = "figures/Figure6_Heatmap_top60_upregulated.png",
  width = 10,
  height = 14)

# Plot top 60 downregulated genes
pheatmap(
  down_matrix,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  color = heatmap_colors,
  border_color = NA,
  fontsize_row = 7,
  fontsize_col = 9,
  main = "Top 60 downregulated genes",
  filename = "figures/Figure7_Heatmap_top60_downregulated.png",
  width = 10,
  height = 14)

# DAVID pathway analysis figures ----

# Import DAVID results
david_up <- read.csv(
  "DAVID_results/DAVIDFunctAnnotClusterReport_upregulated.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE)

david_down <- read.csv(
  "DAVID_results/DAVIDFunctAnnotClusterReport_downregulated.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE)

# Inspect DAVID result structure
david_up <- david_up %>%
  filter(!is.na(Category), Category != "")

david_down <- david_down %>%
  filter(!is.na(Category), Category != "")
# Prepare upregulated pathway results ----
david_up_plot <- david_up %>%
  filter(!is.na(Category), Category != "") %>%
  mutate(
    Term = sub("^.*~", "", Term),
    fold_enrichment = as.numeric(`Fold Enrichment`),
    p_value = as.numeric(`P-Value`),
    neg_log10_p = -log10(p_value)
  ) %>%
  arrange(p_value) %>%
  filter(!duplicated(Term)) %>%
  slice_head(n = 15) %>%
  mutate(Term = factor(Term, levels = rev(unique(Term))))

# Prepare downregulated pathway results ----
david_down_plot <- david_down %>%
  filter(!is.na(Category), Category != "") %>%
  mutate(
    Term = sub("^.*~", "", Term),
    fold_enrichment = as.numeric(`Fold Enrichment`),
    p_value = as.numeric(`P-Value`),
    neg_log10_p = -log10(p_value)
  ) %>%
  arrange(p_value) %>%
  filter(!duplicated(Term)) %>%
  slice_head(n = 15) %>%
  mutate(Term = factor(Term, levels = rev(unique(Term))))

# Plot upregulated pathways ----
david_up_figure <- ggplot(
  david_up_plot,
  aes(
    x = fold_enrichment,
    y = Term,
    size = Count,
    color = neg_log10_p)) +
  geom_point(alpha = 0.85) +
  scale_color_gradient(
    low = "#2878B5",
    high = "#D94841") +
  labs(
    title = "DAVID functional enrichment of upregulated genes",
    x = "Fold enrichment",
    y = NULL,
    size = "Gene count",
    color = "-log10 P value") +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0),
    axis.text.y = element_text(size = 10))

ggsave(
  "figures/Figure8_DAVID_upregulated.png",
  plot = david_up_figure,
  width = 11,
  height = 8,
  dpi = 600)

# Plot downregulated pathways ----
david_down_figure <- ggplot(
  david_down_plot,
  aes(
    x = fold_enrichment,
    y = Term,
    size = Count,
    color = neg_log10_p)) +
  geom_point(alpha = 0.85) +
  scale_color_gradient(
    low = "#2878B5",
    high = "#D94841") +
  labs(
    title = "DAVID functional enrichment of downregulated genes",
    x = "Fold enrichment",
    y = NULL,
    size = "Gene count",
    color = "-log10 P value") +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0),
    axis.text.y = element_text(size = 10))

ggsave(
  "figures/Figure9_DAVID_downregulated.png",
  plot = david_down_figure,
  width = 11,
  height = 8,
  dpi = 600)

# Literature validation and final outputs ----

# Define literature-supported genes
induced_genes <- c(
  "DUSP1", "FKBP5", "KLF15", "TSC22D3",
  "PER1", "CRISPLD2", "SERPINA3", "PTX3")

repressed_genes <- c(
  "KCTD12", "VCAM1", "TSLP", "TNFAIP6",
  "WNT2", "SOX4", "BDKRB2", "AQP3")

# Prepare normalized expression data
expression_df <- as.data.frame(assay(vsd))
expression_df$gene_id <- rownames(expression_df)

expression_df <- expression_df %>%
  dplyr::left_join(
    res_df %>% dplyr::select(gene_id, SYMBOL),
    by = "gene_id") %>%
  filter(!is.na(SYMBOL), SYMBOL != "")

expression_df <- expression_df %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

# Convert expression matrix to long format
expression_long <- expression_df %>%
  dplyr::select(SYMBOL, dplyr::all_of(sample_order)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(sample_order),
    names_to = "sample",
    values_to = "expression") %>%
  dplyr::left_join(
    metadata %>%
      tibble::rownames_to_column("sample"),
    by = "sample")

expression_long$treatment <- factor(
  expression_long$treatment,
  levels = c("Untreated", "Dexamethasone"))

# Define induced gene expression
induced_expression <- expression_long %>%
  filter(SYMBOL %in% induced_genes) %>%
  mutate(SYMBOL = factor(SYMBOL, levels = induced_genes))

# Define repressed gene expression
repressed_expression <- expression_long %>%
  filter(SYMBOL %in% repressed_genes) %>%
  mutate(SYMBOL = factor(SYMBOL, levels = repressed_genes))

# Calculate paired significance for induced genes
induced_stats <- induced_expression %>%
  group_by(SYMBOL) %>%
  summarise(
    p_value = t.test(
      expression[treatment == "Untreated"],
      expression[treatment == "Dexamethasone"],
      paired = TRUE
    )$p.value,
    .groups = "drop") %>%
  mutate(
    label = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"))

# Calculate paired significance for repressed genes
repressed_stats <- repressed_expression %>%
  group_by(SYMBOL) %>%
  summarise(
    p_value = t.test(
      expression[treatment == "Untreated"],
      expression[treatment == "Dexamethasone"],
      paired = TRUE
    )$p.value,
    .groups = "drop") %>%
  mutate(
    label = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"))

# Define annotation positions for induced genes
induced_positions <- induced_expression %>%
  group_by(SYMBOL) %>%
  summarise(
    y_position = max(expression, na.rm = TRUE),
    .groups = "drop") %>%
  left_join(induced_stats, by = "SYMBOL") %>%
  mutate(
    y_position = y_position + abs(y_position) * 0.12 + 0.2,
    xmin = 1,
    xmax = 2)

# Define annotation positions for repressed genes
repressed_positions <- repressed_expression %>%
  group_by(SYMBOL) %>%
  summarise(
    y_position = max(expression, na.rm = TRUE),
    .groups = "drop") %>%
  left_join(repressed_stats, by = "SYMBOL") %>%
  mutate(
    y_position = y_position + abs(y_position) * 0.12 + 0.2,
    xmin = 1,
    xmax = 2)

# Plot literature-supported induced genes
induced_plot <- ggplot(
  induced_expression,
  aes(x = treatment, y = expression, fill = treatment)) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA,
    color = "black") +
  geom_jitter(
    aes(shape = donor),
    width = 0.025,
    size = 1.3,
    color = "black") +
  geom_segment(
    data = induced_positions,
    aes(
      x = xmin,
      xend = xmax,
      y = y_position,
      yend = y_position),
    inherit.aes = FALSE,
    color = "black") +
  geom_text(
    data = induced_positions,
    aes(
      x = 1.5,
      y = y_position + 0.05,
      label = label),
    inherit.aes = FALSE,
    size = 3) +
  facet_wrap(
    ~SYMBOL,
    ncol = 4,
    scales = "free_y",
    axes = "all",
    axis.labels = "all") +
  scale_fill_manual(
    values = c(
      "Untreated" = "#4DBBD5",
      "Dexamethasone" = "#E64B35")) +
  labs(
    title = "Literature-supported glucocorticoid-responsive genes",
    x = NULL,
    y = "VST expression") +
  theme_classic(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0),
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.7, "lines")) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.20)))

ggsave(
  "figures/Figure10_Literature_validation_induced.png",
  plot = induced_plot,
  width = 12,
  height = 10,
  dpi = 600)

# Plot literature-supported repressed genes
repressed_plot <- ggplot(
  repressed_expression,
  aes(x = treatment, y = expression, fill = treatment)) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA,
    color = "black") +
  geom_jitter(
    aes(shape = donor),
    width = 0.025,
    size = 1.3,
    color = "black") +
  geom_segment(
    data = repressed_positions,
    aes(
      x = xmin,
      xend = xmax,
      y = y_position,
      yend = y_position),
    inherit.aes = FALSE,
    color = "black") +
  geom_text(
    data = repressed_positions,
    aes(
      x = 1.5,
      y = y_position + 0.05,
      label = label),
    inherit.aes = FALSE,
    size = 3) +
  facet_wrap(
    ~SYMBOL,
    ncol = 4,
    scales = "free_y",
    axes = "all",
    axis.labels = "all") +
  scale_fill_manual(
    values = c(
      "Untreated" = "#4DBBD5",
      "Dexamethasone" = "#E64B35")) +
  labs(
    title = "Literature-supported repressed genes",
    x = NULL,
    y = "VST expression") +
  theme_classic(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0),
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.7, "lines")) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.20)))

ggsave(
  "figures/Figure11_Literature_validation_repressed.png",
  plot = repressed_plot,
  width = 12,
  height = 10,
  dpi = 600)

# Save final significant DEG tables
write.csv(
  res_df %>% filter(significance == "Upregulated") %>% arrange(padj),
  "results/all_upregulated_genes.csv",
  row.names = FALSE)

write.csv(
  res_df %>% filter(significance == "Downregulated") %>% arrange(padj),
  "results/all_downregulated_genes.csv",
  row.names = FALSE)

# Save session information
session_info <- capture.output(sessionInfo())

writeLines(
  session_info,
  "results/sessionInfo.txt")
