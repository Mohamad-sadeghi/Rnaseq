# =============================
# 🔧 KEGG Analysis ( Counts for DE + Optional TPM for plots)
# =============================
 

# اگر نیاز بود:
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("DESeq2","clusterProfiler","enrichplot"))
# install.packages(c("readxl","openxlsx","dplyr","tibble","ggplot2","ggrepel","pheatmap","ggpubr","reshape2"))

suppressPackageStartupMessages({
  library(readxl)
  library(DESeq2)
  library(clusterProfiler)
  library(enrichplot)
  library(openxlsx)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(ggpubr)
  library(reshape2)
  library(tidyr)
})
 

# =============================
# ⚙️ CONFIGURATION BLOCK (User Editable)
# =============================

# -----------------------------
# 📁 Paths and File Settings
# -----------------------------
work_dir        <- "D:\\ArticleScripts"    # مسیر کار
output_dir      <- "D:\\ArticleScripts\\Results"

file_name       <- "Expression_matrix.xlsx"              # نام فایل اکسل
reads_sheet     <- "Reads"                       # شیت شمارش خام
tpm_sheet       <- "TPM"                         # شیت TPM (اختیاری)

# -----------------------------
# 🔬 Analysis Parameters
# -----------------------------
USE_TPM_FOR_PLOTS <- TRUE                        # استفاده از TPM برای نمودارها
VOLCANO_PVAL      <- 0.05                        # آستانه padj برای ولکانو
L2FC_THRESHOLD    <- 1                           # آستانه |log2FC| برای معنی‌داری


# ایجاد پوشه های خروجی
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "Plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "DEG_Module_Results"), showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 🧬 MODULE DEFINITIONS (User Editable)
# -----------------------------

# پاکسازی فاصله‌ها
trimK <- function(x) gsub("\\s+", "", x)

# ==========================ِDefine your module here==================================Add More Modules




Mevalonate_Module <- c("K00626","K01641","K00021","K00054","K00869","K18689",
                       "K00938","K13273","K25517","K18690","K01597","K17942","K25518","K22813","K06981","K01823","K09128","K03186")

MEP_Module <- c("K01662","K00099","K00991","K00919","K01770","K03526","K12506",
                "K03527","K01823","K14066","K00787","K00795","K00804","K13789","K13787")

 
# KEGG Pathway Modules
map00906 <- c("K00514","K02291","K02292","K02293","K02294","K06443","K06444","K08977","K09835","K09836","K09837","K09838","K09839",
              "K09840","K09841","K09842","K09843","K09844","K09845","K09846","K09847","K09879","K10027","K10208","K10209","K10210",
              "K10211","K10212","K14593","K14594","K14595","K14596","K14597","K14598","K14605","K14606","K15744","K15745","K15746",
              "K15747","K15748","K17819","K17841","K17842","K17911","K17912","K17913","K20611","K20616","K22445","K22492","K22502",
              "K23037","K25072","K25073","K25074")

map00908 <- c("K00279","K00791","K10717","K10760","K13492","K13493","K13494","K13495","K13496","K23452")

map00909 <- c("K00511","K00801","K06045","K10156","K10187","K12249","K12250","K12645","K14173","K14174","K14175","K14176","K14177",
              "K14178","K14179","K14180","K14181","K14182","K14183","K14184","K14185","K14186","K15472","K15793","K15794","K15795",
              "K15796","K15797","K15798","K15799","K15800","K15801","K15802","K15803","K15804","K15805","K15806","K15807","K15808",
              "K15809","K15810","K15811","K15812","K15813","K15814","K15815","K15816","K15817","K15818","K15819","K15820","K15821",
              "K15822","K15823","K15891","K15907","K16204","K16205","K16206","K16207","K18109","K18110","K18111","K18112","K19010",
              "K19011","K20561","K20658","K20659","K21927","K21928","K21984","K22064","K22065","K26054")

map00100 <- c("K00213","K00222","K00223","K00227","K00511","K00559","K00637","K00801","K01052","K01824","K01852","K01853",
              "K05917","K07419","K07436","K07438","K07748","K07750","K08242","K08246","K09827","K09828","K09829","K09831",
              "K09832","K12298","K13373","K14423","K14424","K14674","K19532","K21146","K23485","K23558")

map00130 <- c("K00355","K00457","K00487","K00568","K00591","K00815","K00838","K01075","K01661","K01851","K01904","K01911",
              "K02361","K02548","K02549","K02551","K02552","K03179","K03181","K03182","K03183","K03184","K03185","K03186",
              "K03809","K05357","K05928","K06125","K06126","K06127","K06134","K08680","K09833","K09834","K10106","K11782",
              "K11783","K11784","K11785","K12073","K12501","K12502","K13565","K14759","K14760","K17872","K18240","K18284",
              "K18285","K18286","K18534","K18586","K18606","K18800","K19222","K19267","K20810","K23094","K23095","K24441",
              "K24843","K24844","K24845","K28034")

map00403 <- c("K18385","K18386","K18387","K18388","K18389","K18390","K18391","K18392","K18393","K18394","K18395","K18396","K18397")

map00510 <- c("K00717","K00721","K00726","K00729","K00730","K00736","K00737","K00738","K00744","K00778","K00779","K00902","K01001",
              "K01228","K01230","K01231","K03842","K03843","K03844","K03845","K03846","K03847","K03848","K03849","K03850","K05546",
              "K07151","K07252","K07432","K07441","K07966","K07967","K07968","K09658","K09659","K09661","K11170","K12345","K12666",
              "K12667","K12668","K12669","K12670","K12691","K13748","K19478","K23741","K26501","K26502","K28202")

map00270 <- c("K00003","K00016","K00024","K00025","K00026","K00058","K00133","K00456","K00544","K00547","K00548","K00549","K00552",
              "K00558","K00640","K00641","K00651","K00772","K00789","K00797","K00802","K00811","K00812","K00813","K00815","K00816",
              "K00826","K00827","K00831","K00832","K00837","K00838","K00872","K00899","K00928","K01011","K01243","K01244","K01251",
              "K01505","K01566","K01611","K01697","K01738","K01739","K01740","K01752","K01758","K01760","K01761","K01762","K01767",
              "K01919","K01920","K03334","K05396","K05810","K05933","K05953","K07173","K08963","K08964","K08965","K08966","K08967",
              "K08968","K08969","K09758","K09880","K10150","K10764","K11204","K11205","K11358","K12339","K12524","K12525","K12526",
              "K12960","K13034","K13060","K13061","K13062","K14155","K14454","K14455","K16054","K16843","K16844","K16845","K16846",
              "K17069","K17216","K17217","K17398","K17399","K17462","K17950","K17989","K18284","K19696","K20021","K20248","K20249",
              "K20250","K20772","K21456","K21623","K22207","K22846","K22847","K22954","K22955","K22956","K22957","K22968","K23304",
              "K23370","K23975","K23976","K23977","K24034","K24042","K25035","K25316","K25317","K27857","K28205")

Taxadine_module <- c("K12921","K14039","K12923","K27503","K12924","K20512","K12926","K27849","K27850","K28367", "K20709")

Gibberellin_biosynthesis <- c("K20657","K21292","K12917","K12918","K12919","K05282","K04124","K04120","K04121","K04122","K04123")

PhenylAlanine_Module <- c("K01850","K15849","K05359","K04092","K14187","K04093","K04516","K06208","K06209",
                          "K01713","K04518","K14170","K00832","K00838")

competitor_genes <- c("K12742", "K05356", "K10960", "K15888", "K05355", "K21268", "K12504", "K02523", "K12505", "K00805", "K24873",
                      "K05355", "K21274", "K21275", "K00806", "K06447", "K11778", "K19177", "K05954", "K05955", "K15793", "K12503",
                      "K15887", "K10208", "K02291", "K00801", "K20986", "K20979", "K15086", "K22049", "K12467", "K21925", "K21926",
                      "K15087", "K15088", "K15096", "K07384", "K22050", "K15097", "K07385", "K18108", "K21925", "K22208", "K21938",
                      "K15098", "K00791", "K10760", "K18385")
# -----------------------------
# 🧪 PATHWAY STEP DEFINITIONS (User Editable)
# -----------------------------

# Pathway step rules for precursor analysis
pathway_steps <- tibble::tribble(
  ~pathway, ~step,     ~rule,
  "M00364", "Step_1",  "K01823",
  "M00364", "Step_2",  "K00795|K13789|K13787",
  
  "M00365", "Step_1",  "K01823",
  "M00365", "Step_2",  "K13787",
  
  "M00366", "Step_1",  "K01823",
  "M00366", "Step_2",  "K14066",
  "M00366", "Step_3",  "K00787",
  "M00366", "Step_4",  "K13789",
  
  "M00367", "Step_1",  "K01823",
  "M00367", "Step_2",  "K00787",
  "M00367", "Step_3",  "K00804",
  
  "M00095", "Step_1",  "K00626",
  "M00095", "Step_2",  "K01641",
  "M00095", "Step_3",  "K00021",
  "M00095", "Step_4",  "K00869",
  "M00095", "Step_5",  "K00938|K13273",
  "M00095", "Step_6",  "K01597",
  "M00095", "Step_7",  "K01823",
  
  # Archaeal mevalonate
  "M00849", "Step_1",  "K00626",
  "M00849", "Step_2",  "K01641",
  "M00849", "Step_3",  "K00021|K00054",
  "M00849", "Step_4",  "K00869|K18689",
  "M00849", "Step_5",  "K06981",
  "M00849", "Step_6",  "K01823",
  
  # Non-mevalonate
  "M00096", "Step_1",  "K01662",
  "M00096", "Step_2",  "K00099",
  "M00096", "Step_3",  "K00991|K12506",
  "M00096", "Step_4",  "K00919",
  "M00096", "Step_5",  "K01770|K12506",
  "M00096", "Step_6",  "K03526",
  "M00096", "Step_7",  "K03527",
  "M00096", "Step_8",  "K01823",
  
  # Taxadiene 
  "Taxadiene", "Step_1",  "K12921",
  "Taxadiene", "Step_2",  "K14039",
  "Taxadiene", "Step_3",  "K12923",
  "Taxadiene", "Step_4",  "K12924",
  "Taxadiene", "Step_5",  "K12926",
  "Taxadiene", "Step_6",  "K27849",
  "Taxadiene", "Step_7",  "K28368",
  "Taxadiene", "Step_8",  "K27850",
  
  # Gibberellin
  "M00927", "Step_1",  "K04120",
  "M00927", "Step_2",  "K04121",
  "M00927", "Step_3",  "K04122",
  "M00927", "Step_4",  "K04123"
)

# Module blocks for visual grouping
module_block <- c(
  M00095 = "Mevalonate  ",
  M00849 = "Mevalonate(archaea) ",
  M00096 = "Non-mevalonate",
  M00364 = "M00364 (alt)",
  M00365 = "M00365 (alt)",
  M00366 = "M00366 (alt)",
  M00367 = "M00367 (alt)",
  Taxadiene = "Taxadiene",
  M00927 = "Gibberellin" 
)

# Custom steps for bubble chart
custom_steps <- tibble::tribble(
  ~pathway, ~step,     ~rule,
  
  # Mevalonate pathway (Archaeal)
  "M00849", "Step_1",  "K00626",
  "M00849", "Step_2",  "K01641",
  "M00849", "Step_3",  "K00021|K00054",
  "M00849", "Step_4",  "K00869&K17942 | K00869&K25517&K09128&K25518&K03186 | K18689&K18690&K22813",
  "M00849", "Step_5",  "K06981",
  "M00849", "Step_6",  "K01823",
  
  # Mevalonate pathway (Plants)
  "M00095", "Step_1",  "K00626",
  "M00095", "Step_2",  "K01641",
  "M00095", "Step_3",  "K00021",
  "M00095", "Step_4",  "K00869",
  "M00095", "Step_5",  "K00938|K13273",
  "M00095", "Step_6",  "K01597",
  "M00095", "Step_7",  "K01823",
  
  # Non-mevalonate
  "M00096", "Step_1",  "K01662",
  "M00096", "Step_2",  "K00099",
  "M00096", "Step_3",  "K00991|K12506",
  "M00096", "Step_4",  "K00919",
  "M00096", "Step_5",  "K01770|K12506",
  "M00096", "Step_6",  "K03526",
  "M00096", "Step_7",  "K03527",
  "M00096", "Step_8",  "K01823",
  
  # C10-C20 isoprenoid biosynthesis pathways
  "M00364", "Step_1",  "K01823",
  "M00364", "Step_2",  "K00795|K13789|K13787",
  
  "M00365", "Step_1",  "K01823",
  "M00365", "Step_2",  "K13787",
  
  "M00366", "Step_1",  "K01823",
  "M00366", "Step_2",  "K14066",
  "M00366", "Step_3",  "K00787",
  "M00366", "Step_4",  "K13789",
  
  "M00367", "Step_1",  "K01823",
  "M00367", "Step_2",  "K00787",
  "M00367", "Step_3",  "K00804",
  
  # Taxadiene
  "Taxadiene", "Step_1",  "K12921",
  "Taxadiene", "Step_2",  "K14039",
  "Taxadiene", "Step_3",  "K12923",
  "Taxadiene", "Step_4",  "K12924",
  "Taxadiene", "Step_5",  "K12926",
  "Taxadiene", "Step_6",  "K27849",
  
  # Gibberellin
  "M00927", "Step_1",  "K04120",
  "M00927", "Step_2",  "K04121",
  "M00927", "Step_3",  "K04122",
  "M00927", "Step_4",  "K04123"
)

# Module labels for bubble chart
labels <- c(
  M00849 = "Mevalonate • M00849(archaea)",
  M00095 = "Mevalonate • M00095",
  M00096 = "Non-Mevalonate • M00096",
  M00364 = "C10-C20 isoprenoid biosynthesis bacteria",
  M00365 = "C10-C20 isoprenoid biosynthesis archaea",
  M00366 = "C10-C20 isoprenoid biosynthesis plants",
  M00367 = "C10-C20 isoprenoid biosynthesis eukaryotes",
  M00927 = "Gibberellin",
  Taxadiene = "Taxadiene"
)

# =============================
# 🚀 EXECUTION CODE (Do not modify below unless necessary)
# =============================

# Clean KO lists
Taxadine_module          <- trimK(Taxadine_module)
Mevalonate_Module        <- trimK(Mevalonate_Module)
MEP_Module               <- trimK(MEP_Module)
PhenylAlanine_Module     <- trimK(PhenylAlanine_Module)
Gibberellin_biosynthesis <- trimK(Gibberellin_biosynthesis)
competitor_genes         <- trimK(competitor_genes)
map00906                 <- trimK(map00906)
map00908                 <- trimK(map00908)
map00909                 <- trimK(map00909)
map00100                 <- trimK(map00100)
map00130                 <- trimK(map00130)
map00403                 <- trimK(map00403)
map00510                 <- trimK(map00510)
map00270                 <- trimK(map00270)

# Combine all selected KOs
selected_KOs <- unique(c(Mevalonate_Module, MEP_Module, PhenylAlanine_Module, 
                         Gibberellin_biosynthesis, Taxadine_module,
                         map00906, map00908, map00909, map00270))


# -----------------------------
# 📥 2) خواندن داده‌ها (Counts + TPM)
# -----------------------------
message("📥 Reading counts from sheet:", reads_sheet)
reads_data <- read_excel(file_name, sheet = reads_sheet)
stopifnot(ncol(reads_data) >= 2)

gene_names <- reads_data[[1]]
expr_counts <- reads_data[, -1]
colnames(expr_counts) <- colnames(reads_data)[-1]
expr_counts <- as.data.frame(lapply(expr_counts, function(x) as.integer(as.numeric(x))))
rownames(expr_counts) <- gene_names
expr_counts[is.na(expr_counts)] <- 0

# سعی برای خواندن TPM
expr_tpm <- NULL
if (USE_TPM_FOR_PLOTS) {
  tp <- tryCatch(read_excel(file_name, sheet = tpm_sheet), error = function(e) NULL)
  if (!is.null(tp)) {
    message("📥 Reading TPM from sheet:", tpm_sheet)
    gene_names_tpm <- tp[[1]]
    expr_tpm0 <- tp[, -1]
    colnames(expr_tpm0) <- colnames(tp)[-1]
    expr_tpm <- as.data.frame(lapply(expr_tpm0, function(x) as.numeric(x)))
    rownames(expr_tpm) <- gene_names_tpm
  } else {
    warning("⚠️ TPM sheet not found; plots will use VSD instead.")
    USE_TPM_FOR_PLOTS <- FALSE
  }
}

# فقط ژن‌های ماژولی که در ماتریس وجود دارند
genes_available <- intersect(selected_KOs, rownames(expr_counts))
stopifnot(length(genes_available) > 0)
expr_counts <- expr_counts[genes_available, , drop = FALSE]
if (!is.null(expr_tpm)) expr_tpm <- expr_tpm[genes_available, intersect(colnames(expr_tpm), colnames(expr_counts)), drop = FALSE]

# -----------------------------
# 👥 3) گروه‌ها و نام‌گذاری نمونه‌ها (مطابق نسخه‌ی شما)
# -----------------------------
# توجه: این بخش را مطابق داده‌های خود تنظیم کنید.
# همان ساختار قبلی شما حفظ شده است.

group_list <- list(
  "T.Mairei"    = c("X8080080",   "X8080082"),
  "T.Cuspidata" = c("X8080084", "X8080085"),
  "T.xMedia"    = c("X8080086",   "X8080087"),
 # "Corylus alevena"      = c("X3303631","X3303632"),
 # "Nicotiana benthamiana"= c("X6076955","X6076957"), 
 
  "T.Chinensis KL27" = c("X15317970","X15317979") ,
  "T.Chinensis " = c("X15317980", "X15317981") , 
  "T.wallichiana(Ly)" = c("X25203150") ,
  "T.wallichiana(Hy)" = c("X25203151") , 
  "T.Chinensis LL" = c("SRR19646424", "SRR19646422"),
  "T.Chinensis HL" = c("SRR19646426","SRR19646427")  
  
)

sample_names <- c(
  "X8080080" = "T.Mairei 1",
  "X8080082" = "T.Mairei 2",
  "X8080084" = "T.Cuspidata 1",
  "X8080085" = "T.Cuspidata 2",
  "X8080086" = "T.xMedia 1",
  "X8080087" = "T.xMedia 2",
  "X15317980" = "T.Chinensis 1",
  "X15317981" = "T.Chinensis 2",
  "X15317970" = "T.Chinensis Kl27 1",
  "X15317979" = "T.Chinensis Kl27 2",
  "SRR19646422" = "T.Chinensis  LL 1",
  "SRR19646424" = "T.Chinensis  LL 2",
  "SRR19646426" = "T.Chinensis  HL 1",
  "SRR19646427" = "T.Chinensis  HL 2",
  "X25203150" = "T.Wallichiana Ly",
  "X25203151" = "T.Wallichiana Hy",
  "X3303631"  = "C.alevena 1",
  "X3303632"  = "C.alevena 2",
  "X6076955"  = "N.benthamiana 1",
  "X6076956"  = "N.benthamiana 2",

)

selected_samples <- unlist(group_list)
conditions <- unlist(mapply(function(group, samples) rep(group, length(samples)), names(group_list), group_list))

expr_counts_sub <- expr_counts[, selected_samples, drop = FALSE]
col_data <- data.frame(row.names = selected_samples, condition = factor(conditions))

# نام‌نمایشی
rownames(col_data) <- sample_names[rownames(col_data)]
colnames(expr_counts_sub) <- sample_names[colnames(expr_counts_sub)]

# TPM زیرمجموعه برای رسم (در صورت وجود)
expr_tpm_sub <- NULL
if (!is.null(expr_tpm)) {
  # هم‌تراز با samples انتخابی و نام‌های نمایشی
  keep <- intersect(colnames(expr_tpm), selected_samples)
  expr_tpm_sub <- expr_tpm[, keep, drop = FALSE]
  colnames(expr_tpm_sub) <- sample_names[colnames(expr_tpm_sub)]
  # برای هماهنگی کامل ترتیب ستون‌ها را بر اساس col_data تنظیم کن
  expr_tpm_sub <- expr_tpm_sub[, rownames(col_data), drop = FALSE]
}

# -----------------------------
# 🔬 4) DESeq2 روی counts (استاندارد)
# -----------------------------
message("🔬 Running DESeq2 on raw counts…")
dds <- DESeqDataSetFromMatrix(countData = as.matrix(expr_counts_sub), colData = col_data, design = ~ condition)
dds <- DESeq(dds)
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)

# ماتریس برای نمودارها: TPM (ترجیحی) یا VSD
get_matrix <- function(obj) if (is.matrix(obj)) obj else assay(obj)
plot_matrix <- if (USE_TPM_FOR_PLOTS && !is.null(expr_tpm_sub)) as.matrix(expr_tpm_sub) else get_matrix(vsd)

# 🔒 Sanitize matrix for downstream plots (remove NA/Inf, drop zero-variance rows/cols)
sanitize_matrix <- function(m) {
  m <- as.matrix(m)
  storage.mode(m) <- "numeric"
  m[!is.finite(m)] <- NA
  m[is.na(m)] <- 0
  # drop rows/cols with zero variance
  if (nrow(m) > 1) {
    keep_rows <- apply(m, 1, function(r) stats::sd(r) > 0)
  } else { keep_rows <- TRUE }
  if (ncol(m) > 1) {
    keep_cols <- apply(m, 2, function(c) stats::sd(c) > 0)
  } else { keep_cols <- TRUE }
  if (any(!keep_rows)) m <- m[keep_rows, , drop = FALSE]
  if (any(!keep_cols)) m <- m[, keep_cols, drop = FALSE]
  m
}
plot_matrix <- sanitize_matrix(plot_matrix)
# Align col_data with plot_matrix columns
if (!all(colnames(plot_matrix) %in% rownames(col_data))) {
  common <- intersect(colnames(plot_matrix), rownames(col_data))
  plot_matrix <- plot_matrix[, common, drop = FALSE]
}
col_data <- col_data[colnames(plot_matrix), , drop = FALSE]

# -----------------------------
# 🔁 5) Pairwise DEG + Volcano/MA
# -----------------------------
dir.create("Plots", showWarnings = FALSE)

group_names <- levels(col_data$condition)
pairwise_contrasts <- combn(group_names, 2, simplify = FALSE)
all_sig_genes <- character(0)

for (pair in pairwise_contrasts) {
  g1 <- pair[1]; g2 <- pair[2]
  res_df <- results(dds, contrast = c("condition", g1, g2)) %>%
    as.data.frame() %>% tibble::rownames_to_column("Gene") %>% dplyr::filter(complete.cases(.))
  
  res_df$padj_plot <- ifelse(res_df$padj == 0, 1e-300, res_df$padj)
  res_df$significant <- res_df$padj < VOLCANO_PVAL & abs(res_df$log2FoldChange) > L2FC_THRESHOLD
  
  # تجمیع ژن‌های معنی‌دار برای مراحل بعدی
  all_sig_genes <- union(all_sig_genes, res_df$Gene[res_df$significant])
  
  topN <- res_df %>% dplyr::filter(significant) %>% dplyr::arrange(dplyr::desc(abs(log2FoldChange))) %>% dplyr::slice_head(n = 20)
  n_up <- sum(res_df$significant & res_df$log2FoldChange > 0, na.rm = TRUE)
  n_down <- sum(res_df$significant & res_df$log2FoldChange < 0, na.rm = TRUE)
  
  volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj_plot))) +
    geom_point(aes(color = significant), alpha = 0.7, size = 1.7) +
    geom_vline(xintercept = c(-L2FC_THRESHOLD, L2FC_THRESHOLD), linetype = "dashed", color = "grey60") +
    geom_hline(yintercept = -log10(VOLCANO_PVAL), linetype = "dashed", color = "grey60") +
    ggrepel::geom_label_repel(
      data = topN,
      aes(label = sprintf("%s (%.2f)", Gene, log2FoldChange)),
      size = 3, label.size = 0.2, box.padding = 0.05, label.padding = 0.15,
      max.overlaps = Inf, min.segment.length = 0
    ) +
    scale_color_manual(values = c("grey70", "red")) +
    labs(title = paste("Volcano:", g1, "vs", g2),
         subtitle = sprintf("Significant: Up=%d, Down=%d  |  padj<%g & |L2FC|>%g", n_up, n_down, VOLCANO_PVAL, L2FC_THRESHOLD),
         x = "log2 Fold Change", y = "-log10 Adjusted P-value", color = "significant") +
    theme_classic(base_size = 12)
  
  
  ggsave(paste0(output_dir,"/Volcano_", g1, "_vs_", g2, ".png"), plot = volcano, width = 7.5, height = 5.2, dpi = 300)
  
  # MA Plot
  res <- results(dds, contrast = c("condition", g1, g2))
  png(paste0(output_dir,"/MA_", g1, "_vs_", g2, ".png"), width = 800, height = 600, res = 120)
  DESeq2::plotMA(res, ylim = c(-5,5), alpha = VOLCANO_PVAL, main = paste("MA Plot:", g1, "vs", g2))
  dev.off()
}

# -----------------------------
# 📊 6) PCA (robust to NA/zero-variance)
# -----------------------------
message(if (USE_TPM_FOR_PLOTS && !is.null(expr_tpm_sub)) "📊 PCA on TPM matrix" else "📊 PCA on VSD (DESeq2)")
# Use sanitized matrix and align with col_data
pca_mat <- sanitize_matrix(plot_matrix)
col_data_pca <- col_data[colnames(pca_mat), , drop = FALSE]
if (ncol(pca_mat) < 2 || nrow(pca_mat) < 2) {
  warning("Not enough non-constant genes/samples for PCA; skipping.")
} else {
  pca <- prcomp(t(pca_mat), center = TRUE, scale. = TRUE)
  percentVar <- round(100 * (pca$sdev^2/sum(pca$sdev^2)))[1:2]
  pca_df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], condition = col_data_pca$condition, row.names = rownames(col_data_pca))
  p <- ggplot(pca_df, aes(PC1, PC2, color = condition, shape = condition)) +
    geom_point(size = 3) +
    xlab(paste0("PC1: ", percentVar[1], "%")) + ylab(paste0("PC2: ", percentVar[2], "%")) +
    ggtitle(if (USE_TPM_FOR_PLOTS && !is.null(expr_tpm_sub)) "PCA of TPM (module genes)" else "PCA of Module Genes (VSD)") +
    theme_classic()
  ggsave(paste(output_dir,"/PCA_Module_Genes.png"), plot = p, width = 8, height = 6)
}

# -----------------------------
# 🔥 7) Heatmaps by module + All Modules
# -----------------------------
plot_module_heatmap <- function(genes, label, mat_for_plot, col_data,
                                sig_only = TRUE, all_sig_genes = NULL,
                                topN = 40, scale_row = TRUE) {
  genes <- unique(trimK(genes))
  # Align & sanitize matrix
  vsd_mat <- sanitize_matrix(mat_for_plot)
  common <- intersect(colnames(vsd_mat), rownames(col_data))
  vsd_mat <- vsd_mat[, common, drop = FALSE]
  ann <- col_data[common, , drop = FALSE]
  
  # Choose genes (prefer significant if available)
  sel <- character(0)
  if (sig_only && !is.null(all_sig_genes)) {
    sel_sig <- intersect(trimK(all_sig_genes), genes)
    sel_sig <- sel_sig[!is.na(sel_sig)]
    sel_sig <- intersect(sel_sig, rownames(vsd_mat))
  } else sel_sig <- character(0)
  
  if (length(sel_sig) >= 2) {
    sel <- unique(sel_sig)
    title_extra <- sprintf(" (signif, n=%d)", length(sel))
  } else {
    sel_all <- intersect(genes, rownames(vsd_mat))
    if (length(sel_all) == 0) {
      message(sprintf("⛔ No genes from module '%s' found in matrix. Skipped.", label))
      return(invisible(NULL))
    }
    mat0 <- vsd_mat[sel_all, , drop = FALSE]
    vars <- apply(mat0, 1, stats::var)
    vars[!is.finite(vars)] <- 0
    keep <- names(sort(vars[vars > 0], decreasing = TRUE))
    sel <- head(keep, topN)
    if (length(sel) == 0) {
      message(sprintf("⛔ All candidate genes had zero variance for '%s'; skipped.", label))
      return(invisible(NULL))
    }
    title_extra <- sprintf(" (top-variance, n=%d)", length(sel))
  }
  
  sel <- sel[!is.na(sel) & sel %in% rownames(vsd_mat)]
  if (length(sel) == 0) {
    message(sprintf("⛔ Selection for '%s' resulted in empty set; skipped.", label))
    return(invisible(NULL))
  }
  
  mat <- vsd_mat[sel, , drop = FALSE]
  
  if (scale_row) {
    rn <- rownames(mat)
    mat <- t(scale(t(mat)))
    mat[!is.finite(mat)] <- 0
    rownames(mat) <- rn
  }
  
  do_rowclust <- nrow(mat) > 1
  do_colclust <- ncol(mat) > 1
  fn <- file.path(output_dir,   sprintf("Heatmap_%s.png", label))
  png(fn, width = 1200, height = 900, res = 150)
  pheatmap(mat,
           cluster_rows = do_rowclust,
           cluster_cols = do_colclust,
           show_rownames = TRUE,
           show_colnames = TRUE,
           annotation_col = ann,
           main =  file.path(output_dir,  sprintf("Heatmap: %s%s", label, title_extra)))
  dev.off()
  message("✅ Saved: ", fn)
}

sigK <- if (length(all_sig_genes)) trimK(all_sig_genes) else NULL

plot_module_heatmap(Mevalonate_Module,        "MVA_Module", plot_matrix, col_data, FALSE, sigK, topN = 30)
plot_module_heatmap(MEP_Module,               "MEP_Module", plot_matrix, col_data, FALSE, sigK, topN = 30)
plot_module_heatmap(PhenylAlanine_Module,     "PHE_Module", plot_matrix, col_data, FALSE, sigK, topN = 30)
plot_module_heatmap(Gibberellin_biosynthesis, "GA_Module",  plot_matrix, col_data, FALSE, sigK, topN = 30)
plot_module_heatmap(Taxadine_module,          "Taxadiene",   plot_matrix, col_data, FALSE, sigK, topN = 30)

AllModules <- unique(c(Mevalonate_Module, MEP_Module, PhenylAlanine_Module, Gibberellin_biosynthesis))
plot_module_heatmap(AllModules, "All_Modules", plot_matrix, col_data, FALSE, sigK, topN = 60)

# -----------------------------
# 🧪 8) KEGG Enrichment بر اساس KO (در صورت وجود ژن معنی‌دار)
# -----------------------------
if (length(all_sig_genes)) {
  ks <- unique(trimK(all_sig_genes))
  ko_enrich <- tryCatch(
    enrichMKEGG(
      gene          = ks,
      organism      = "ko",
      keyType       = "kegg",
      pAdjustMethod = "BH",
      pvalueCutoff  = 1
    ), error = function(e) { message("enrichMKEGG error: ", e$message); NULL }
  )
  df_enr <- tryCatch(as.data.frame(ko_enrich), error = function(e) NULL)
  
  # robust barplot wrapper (handles different package versions)
  plot_enrich_bar <- function(x, outfile, showCategory = 10, title = "KEGG Enrichment - Barplot") {
    p <- NULL
    used <- NULL
    if (requireNamespace("enrichplot", quietly = TRUE) && "barplot" %in% getNamespaceExports("enrichplot")) {
      p <- try(enrichplot::barplot(x, showCategory = showCategory, title = title) + theme_classic(), silent = TRUE)
      if (!inherits(p, "try-error")) used <- "enrichplot::barplot"
    }
    if (is.null(used) && requireNamespace("clusterProfiler", quietly = TRUE) && "barplot" %in% getNamespaceExports("clusterProfiler")) {
      p <- try(clusterProfiler::barplot(x, showCategory = showCategory, title = title) + theme_classic(), silent = TRUE)
      if (!inherits(p, "try-error")) used <- "clusterProfiler::barplot"
    }
    if (is.null(used)) {
      p <- enrichplot::dotplot(x, showCategory = showCategory, title = "KEGG Enrichment") + theme_classic()
      used <- "enrichplot::dotplot (fallback)"
    }
    message("Enrichment barplot method: ", used)
    ggsave(outfile, plot = p, width = 8, height = 6, dpi = 600)
  }
  if (!is.null(df_enr) && nrow(df_enr) > 0) {
    plot_enrich_bar(ko_enrich, outfile = "Plots/KEGG_Barplot.png", showCategory = 10)
    p_dot <- enrichplot::dotplot(ko_enrich, showCategory = 10, title = "KEGG Enrichment - Dotplot") + theme_classic()
    ggsave( paste(output_dir,"/KEGG_Dotplot.png"), plot = p_dot, width = 8, height = 6, dpi = 300)
    if (nrow(df_enr) >= 2) {
      if (requireNamespace("enrichplot", quietly = TRUE) && "pairwise_termsem" %in% getNamespaceExports("enrichplot")) {
        em <- enrichplot::pairwise_termsem(ko_enrich)
        p_emap <- enrichplot::emapplot(em, showCategory = min(20, nrow(df_enr)))
        ggsave(paste(output_dir,"/KEGG_EnrichmentMap.png"), plot = p_emap, width = 10, height = 8, dpi = 300)
      } else {
        message("pairwise_termsem not available; skipping emapplot.")
      }
    }
    openxlsx::write.xlsx(df_enr, file = file.path(output_dir, "KEGG_Enrichment_KO_Results.xlsx"), overwrite = TRUE)
  } else {
    message("No enriched KEGG terms found (or object is NULL). Skipping KEGG plots.")
  }
} else {
  message("No significant genes to enrich.")
}

# -----------------------------
# 📈 9) Module Activity per Sample (بر اساس منبع رسم)
# -----------------------------
USE_ZSCORE_FOR_ACTIVITY <- TRUE

safe_scale_rows <- function(m) {
  rn <- rownames(m)
  m  <- t(scale(t(m)))
  m[!is.finite(m)] <- NA
  rownames(m) <- rn
  m
}

compute_module_activity <- function(module_genes, vsd_matrix) {
  module_genes <- trimK(module_genes)
  found_genes  <- intersect(module_genes, rownames(vsd_matrix))
  if (length(found_genes) < 2) return(rep(NA, ncol(vsd_matrix)))
  mat <- vsd_matrix[found_genes, , drop = FALSE]
  if (isTRUE(USE_ZSCORE_FOR_ACTIVITY)) mat <- safe_scale_rows(mat)
  colMeans(mat, na.rm = TRUE)
}

vsd_mat_for_activity <- as.matrix(plot_matrix)  # TPM یا VSD
modules_list <- list(
  Mevalonate      = Mevalonate_Module,
  MEP             = MEP_Module,
  Phenylalanine   = PhenylAlanine_Module,
  GA_Biosynthesis = Gibberellin_biosynthesis,
  map_00906          = map00906,
  map_00908=        map00908,
  map_00909=         map00909,
  map_00270=         map00270
)

activity_matrix <- sapply(modules_list, compute_module_activity, vsd_matrix = vsd_mat_for_activity)
activity_df <- as.data.frame(activity_matrix)
if ("Taxadine" %in% names(activity_df)) {
  names(activity_df)[names(activity_df) == "Taxadine"] <- "Taxadiene"
}
activity_df$Sample <- rownames(col_data)
activity_df$Group  <- col_data$condition
activity_df <- activity_df[, c("Sample","Group", setdiff(colnames(activity_df), c("Sample","Group")))]
openxlsx::write.xlsx(activity_df, file = file.path(output_dir,"/Module_Activity_Scores.xlsx"), overwrite = TRUE)

# Boxplots
plots_dir <- output_dir
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)
for (module in colnames(activity_df)[!(colnames(activity_df) %in% c("Sample","Group"))]) {
  p <- ggboxplot(activity_df, x = "Group", y = module, color = "Group",
                 add = "jitter", palette = "jco", outlier.shape = NA) +
    stat_compare_means(method = "anova",
                       label.y = max(activity_df[[module]], na.rm = TRUE) * 1.05) +
    theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste("Module Activity:", module), x = "Group", y = "Activity Score")
  ggsave(filename = file.path(plots_dir, paste0("Boxplot_", module, ".png")), plot = p, width = 8, height = 5)
}

# Group-wise heatmap
group_activity_means <- activity_df %>% group_by(Group) %>% summarise(across(-Sample, ~mean(.x, na.rm = TRUE)))
activity_mat <- as.matrix(group_activity_means[, -1]); rownames(activity_mat) <- group_activity_means$Group
png(file.path(output_dir,"/Module_Activity_Heatmap.png"), width = 1000, height = 700)
pheatmap(activity_mat, cluster_rows = TRUE, cluster_cols = TRUE, display_numbers = TRUE,
         main = "Group-wise Module Activity Heatmap", angle_col = 45)
dev.off()




print("______________________________________compute Module Activity ")
# -----------------------------
# 📊 Bar chart of module activities by group (mean ± SE) - SAFE VERSION
# -----------------------------

# ⛑️ یک لابلر امن برای محور Y بساز (سازگار با همه نسخه‌های `scales`)
lab_si_safe <- function(accuracy = NULL, unit = "") {
  if (!requireNamespace("scales", quietly = TRUE)) {
    return(function(x) x)  # اگر پکیج نبود
  }
  # تلاش برای استفاده از cut_si(unit=...)
  ff <- tryCatch(
    scales::label_number(scale_cut = do.call(scales::cut_si, list(unit = unit)),
                         accuracy = accuracy),
    error = function(e) {
      # اگر شکست خورد (یا نسخه قدیمی/ناسازگار بود)، از label_number ساده استفاده کن
      scales::label_number(accuracy = accuracy, big.mark = ",")
    }
  )
  ff
}

bars_dir <- output_dir
if (!dir.exists(bars_dir)) dir.create(bars_dir, recursive = TRUE)

# بررسی وجود activity_df و داده‌های کافی
if (!exists("activity_df") || nrow(activity_df) == 0) {
  message("⏩ Skipping bar charts: activity_df not found or empty")
} else {
  message("📊 Preparing bar charts from activity_df with ", nrow(activity_df), " rows")
  
  # activity_df: data.frame با ستون‌های Sample, Group و ستون‌های ماژول‌ها
  long_act <- melt(activity_df, id.vars = c("Sample","Group"),
                   variable.name = "Module", value.name = "Activity") |>
    dplyr::filter(is.finite(Activity))
  
  # بررسی وجود داده پس از فیلتر
  if (nrow(long_act) == 0) {
    message("⏩ Skipping bar charts: No finite activity values after filtering")
  } else {
    message("📈 After filtering: ", nrow(long_act), " activity measurements")
    
    summary_act <- long_act |>
      dplyr::group_by(Group, Module) |>
      dplyr::summarise(
        n    = sum(!is.na(Activity)),
        mean = mean(Activity, na.rm = TRUE),
        sd   = sd(Activity,   na.rm = TRUE),
        se   = sd / sqrt(pmax(n, 1)),
        .groups = "drop"
      )
    
    # بررسی وجود داده در summary_act
    if (nrow(summary_act) == 0) {
      message("⏩ Skipping bar charts: summary_act is empty after summarization")
    } else {
      message("📊 Summary statistics: ", nrow(summary_act), " group-module combinations")
      
      # ترتیب گروه‌ها مطابق col_data (اگر موجود است)
      if (exists("col_data") && "condition" %in% names(col_data)) {
        summary_act$Group <- factor(summary_act$Group, levels = levels(col_data$condition))
      } else {
        summary_act$Group <- factor(summary_act$Group)
      }
      
      # ========== گزینه A: Facet با محور مستقل ==========
      create_facet_plot <- function() {
        # بررسی وجود حداقل یک ماژول
        if (length(unique(summary_act$Module)) == 0) {
          message("⏩ Skipping facet plot: No modules available")
          return(FALSE)
        }
        
        tryCatch({
          p_facet <- ggplot(summary_act, aes(x = Group, y = mean)) +
            geom_col(fill = "grey35", width = .72) +
            geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = .22) +
            facet_wrap(~ Module, scales = "free_y") +
            scale_y_continuous(labels = lab_si_safe()) +
            labs(title = "Module activity by group (mean ± SE)",
                 x = "Group", y = "Activity score") +
            theme_bw(base_size = 12) +
            theme(axis.text.x = element_text(angle = 35, hjust = 1),
                  panel.grid.minor = element_blank())
          
          ggsave(file.path(output_dir, "Module_Activity_Barplot_Facet.png"),
                 plot = p_facet, width = 12, height = 8, dpi = 600)
          message("✅ Saved: Module_Activity_Barplot_Facet.png")
          return(TRUE)
        }, error = function(e) {
          message("❌ Error in facet plot: ", e$message)
          return(FALSE)
        })
      }
      
      # ========== گزینه B: گروهی با محور لگاریتمی ==========
      create_grouped_log_plot <- function() {
        if (nrow(summary_act) == 0 || !all(summary_act$mean > 0, na.rm = TRUE)) {
          message("⏩ Skipping log plot: Not all mean values are positive")
          return(FALSE)
        }
        
        tryCatch({
          p_grouped_log <- ggplot(summary_act, aes(x = Group, y = mean, fill = Module)) +
            geom_col(position = position_dodge(.85), width = .72) +
            geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                          position = position_dodge(.85), width = .22) +
            scale_y_continuous(
              trans   = "log10",
              labels  = lab_si_safe(accuracy = 1),
              minor_breaks = NULL
            ) +
            labs(title = "Module activity by group (mean ± SE, log10)",
                 x = "Group", y = "Activity (log10)") +
            theme_bw(base_size = 12) +
            theme(axis.text.x = element_text(angle = 35, hjust = 1),
                  panel.grid.minor = element_blank())
          
          ggsave(file.path(output_dir, "Module_Activity_Barplot_Grouped_log.png"),
                 plot = p_grouped_log, width = 11, height = 6.5, dpi = 600)
          message("✅ Saved: Module_Activity_Barplot_Grouped_log.png")
          return(TRUE)
        }, error = function(e) {
          message("❌ Error in log plot: ", e$message)
          return(FALSE)
        })
      }
      
      # ========== گزینه C: نسبی نسبت به مرجع ==========
      create_relative_plot <- function() {
        baseline_group <- "T.Chinensis "
        
        # بررسی وجود گروه baseline در داده‌ها
        if (!baseline_group %in% summary_act$Group) {
          message("⏩ Skipping relative plot: Baseline group '", baseline_group, "' not found in data")
          message("   Available groups: ", paste(unique(summary_act$Group), collapse = ", "))
          return(FALSE)
        }
        
        ref_means <- dplyr::filter(summary_act, Group == baseline_group) |>
          dplyr::select(Module, ref_mean = mean)
        
        # بررسی وجود reference means
        if (nrow(ref_means) == 0) {
          message("⏩ Skipping relative plot: No reference means found for baseline group")
          return(FALSE)
        }
        
        summary_rel <- dplyr::inner_join(summary_act, ref_means, by = "Module") |>
          dplyr::mutate(mean_rel = 100 * mean / ref_mean,
                        se_rel   = 100 * se   / ref_mean)
        
        # بررسی داده‌های نسبی
        if (nrow(summary_rel) == 0) {
          message("⏩ Skipping relative plot: No data after joining with reference")
          return(FALSE)
        }
        
        tryCatch({
          p_facet_rel <- ggplot(summary_rel, aes(x = Group, y = mean_rel)) +
            geom_col(fill = "steelblue", width = .72) +
            geom_errorbar(aes(ymin = mean_rel - se_rel, ymax = mean_rel + se_rel), width = .22) +
            facet_wrap(~ Module, scales = "free_y") +
            labs(title = paste0("Relative activity to ", baseline_group, " (mean ± SE)"),
                 x = "Group", y = "Relative activity (%)") +
            theme_bw(base_size = 12) +
            theme(axis.text.x = element_text(angle = 35, hjust = 1),
                  panel.grid.minor = element_blank())
          
          ggsave(file.path(output_dir, "Module_Activity_Barplot_Facet_Relative.png"),
                 plot = p_facet_rel, width = 12, height = 8, dpi = 300)
          message("✅ Saved: Module_Activity_Barplot_Facet_Relative.png")
          return(TRUE)
        }, error = function(e) {
          message("❌ Error in relative plot: ", e$message)
          return(FALSE)
        })
      }
      
      # اجرای تمام گزینه‌ها
      create_facet_plot()
      create_grouped_log_plot()
      create_relative_plot()
    }
  }
}

# Composite Index (CI)
if (exists("activity_df") && all(c("Taxadiene","GA_Biosynthesis") %in% colnames(activity_df))) {
  message("📊 Calculating Composite Index...")
  activity_df$CI <- activity_df$Taxadiene - activity_df$GA_Biosynthesis
  
  # ذخیرهٔ کامل
  openxlsx::write.xlsx(activity_df,
                       file = file.path(output_dir, "Module_Activity_Scores_with_CI.xlsx"),
                       overwrite = TRUE)
  message("✅ Saved: Module_Activity_Scores_with_CI.xlsx")
  
  # میانگین گروهی
  group_ci <- activity_df |>
    dplyr::group_by(Group) |>
    dplyr::summarise(
      Taxadiene       = mean(Taxadiene,       na.rm=TRUE),
      GA_Biosynthesis = mean(GA_Biosynthesis, na.rm=TRUE),
      CI              = mean(CI,              na.rm=TRUE),
      .groups = "drop"
    )
  openxlsx::write.xlsx(group_ci,
                       file = file.path(output_dir, "Composite_Index_CI_GroupMeans.xlsx"),
                       overwrite = TRUE)
  message("✅ Saved: Composite_Index_CI_GroupMeans.xlsx")
} else {
  message("⏩ Skipping Composite Index: Required columns not found in activity_df")
  if (exists("activity_df")) {
    message("   Available columns: ", paste(colnames(activity_df), collapse = ", "))
  }
}
# =====================================================================
# 🌿 Precursor Pathways: Presence (step-level) + Activity (module/macro)
# Refactored, tidy, robust — drop-in replacement for the new section
# =====================================================================

print(">>>>>Check Presence (step-level) + Activity (module/macro)__________________ ")
# -----------------------------
# ⚙️ Config
# -----------------------------
# Threshold for "present" (on VSD scale). Keep 0 if presence is binary (>0).
PRESENCE_THRESH <- 0

# -----------------------------
# 🧱 Matrix prep
# -----------------------------
vsd_mat <- if (is.matrix(vsd)) vsd else tryCatch(assay(vsd), error = function(e) as.matrix(vsd))
storage.mode(vsd_mat) <- "numeric"
rownames(vsd_mat) <- trimws(rownames(vsd_mat))
samples <- colnames(vsd_mat)

# Attach groups (from existing col_data) if available
sample_to_group <- NULL
if (exists("col_data") && all(samples %in% rownames(col_data))) {
  sample_to_group <- setNames(as.character(col_data[samples, "condition"]), samples)
}

# -----------------------------
# 🔎 Rule evaluation helpers
# -----------------------------
split_or  <- function(x) strsplit(x, "\\|", perl = TRUE)[[1]] |> trimws()
split_and <- function(x) strsplit(x, "&",  perl = TRUE)[[1]] |> trimws()

eval_step_presence <- function(rule, sample_vec, thresh = PRESENCE_THRESH) {
  if (is.na(rule) || !nzchar(rule)) return(FALSE)
  combos <- split_or(rule)
  any(vapply(combos, function(cb) {
    genes <- split_and(cb)
    all(genes %in% names(sample_vec)) &&
      all(sample_vec[genes] > thresh, na.rm = TRUE)
  }, logical(1)))
}

eval_step_score <- function(rule, sample_vec) {
  if (is.na(rule) || !nzchar(rule)) return(0)
  combos <- split_or(rule)
  scores <- vapply(combos, function(cb) {
    genes <- split_and(cb)
    if (!all(genes %in% names(sample_vec))) return(0)
    min(sample_vec[genes], na.rm = TRUE)  # AND bottleneck
  }, numeric(1))
  max(scores, na.rm = TRUE)               # OR best branch
}
print(">>>> Compute Coverage ____________ ")
# -----------------------------
# 🧮 Compute coverage + activity (step & module)
# -----------------------------
compute_coverage_for_sample <- function(sample_name) {
  sample_vec <- vsd_mat[, sample_name, drop = TRUE]
  names(sample_vec) <- rownames(vsd_mat)
  pathway_steps |>
    mutate(
      covered    = vapply(rule, eval_step_presence, logical(1), sample_vec = sample_vec),
      step_score = vapply(rule, eval_step_score,  numeric(1),  sample_vec = sample_vec),
      sample     = sample_name
    )
}

coverage_df <- bind_rows(lapply(samples, compute_coverage_for_sample)) |>
  mutate(
    step_id = paste(pathway, step, sep = "_"),
    Group   = if (!is.null(sample_to_group)) sample_to_group[sample] else NA_character_
  )

pathway_activity <- coverage_df |>
  group_by(sample, pathway, Group) |>
  summarise(
    activity_raw      = mean(step_score,   na.rm = TRUE),
    coverage_fraction = mean(covered,      na.rm = TRUE),
    covered_all_steps = all(covered),
    .groups = "drop"
  ) |>
  group_by(pathway) |>
  mutate(activity_z = as.numeric(scale(activity_raw))) |>
  ungroup()




print(">>>>__ (module/macro)logic______________ ")
# -----------------------------
# 🧩 Macro logic (Upstream/Downstream/Overall)
# -----------------------------
module_complete <- coverage_df |>
  group_by(sample, pathway) |>
  summarise(complete = all(covered), .groups = "drop")

module_activity <- pathway_activity |>
  select(sample, pathway, activity_raw)

presence_summary <- module_complete |>
  tidyr::pivot_wider(names_from = pathway, values_from = complete, values_fill = FALSE) |>
  mutate(
    upstream_mev   = (`M00095` | `M00849`),
    upstream_non   = (`M00096`),
    upstream_any   = (upstream_mev | upstream_non),
    downstream_any = (`M00364` | `M00365` | `M00366` | `M00367`),
    overall_present = (upstream_any & downstream_any)
  )

activity_summary <- module_activity |>
  tidyr::pivot_wider(names_from = pathway, values_from = activity_raw, values_fill = 0) |>
  dplyr::mutate(
    mev_activity_sum     = coalesce(`M00095`, 0) + coalesce(`M00849`, 0),
    non_mev_activity     = coalesce(`M00096`, 0),
    upstream_sum         = mev_activity_sum + non_mev_activity,
    downstream_sum       = rowSums(cbind(
      coalesce(`M00364`, 0),
      coalesce(`M00365`, 0),
      coalesce(`M00366`, 0),
      coalesce(`M00367`, 0)
    ), na.rm = TRUE),
    overall_activity     = pmin(upstream_sum, downstream_sum, na.rm = TRUE)
  )

macro_summary <- presence_summary |>
  left_join(activity_summary, by = "sample") |>
  mutate(Group = if (!is.null(sample_to_group)) sample_to_group[sample] else NA_character_) |>
  select(sample, Group,
         upstream_mev, upstream_non, upstream_any,
         downstream_any, overall_present,
         mev_activity_sum, non_mev_activity, upstream_sum,
         downstream_sum, overall_activity)

# -----------------------------
# 🎨 Plots (white background, consistent themes)
# -----------------------------
print(">>>>__ Plot Coverages______________ ")
theme_heat <- theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "grey40", fill = NA, linewidth = .5))

# 1) Step presence (samples × steps)
p_presence <- ggplot(coverage_df, aes(x = sample, y = step_id, fill = covered)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c(`TRUE` = "#2E7D32", `FALSE` = "#C62828"),
                    labels = c(`TRUE` = "Covered", `FALSE` = "Missing"),
                    name   = "Presence") +
  labs(title = "Precursor Pathways • Step Coverage per Sample",
       x = "Sample", y = "Pathway • Step") +
  theme_heat
ggsave(paste(output_dir,"/Pathway_Step_Coverage_per_Sample.png"), p_presence, width = 13, height = 8, dpi = 300)

# 2) Pathway activity (z) per sample
p_activity <- ggplot(pathway_activity, aes(x = sample, y = pathway, fill = activity_z)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#D32F2F", midpoint = 0,
                       name = "Activity (z)") +
  labs(title = "Precursor Pathways • Activity per Sample (z-scored)",
       x = "Sample", y = "Pathway") +
  theme_heat
ggsave(paste(output_dir,"/Pathway_Activity_per_Sample.png"), p_activity, width = 13, height = 6.5, dpi = 300)

# 3) Step presence with clear module borders (facet by pathway)
coverage_df2 <- coverage_df |>
  mutate(
    block = module_block[pathway],
    facet = factor(paste0(pathway, " • ", block),
                   levels = c("M00095 • Mevalonate",
                              "M00849 • Mevalonate",
                              "M00096 • Non-mevalonate",
                              "M00364 •(alt)",
                              "M00365 •(alt)",
                              "M00366 •(alt)",
                              "M00367 •(alt)",
                              "Giberlline" ,              
                              "Taxadiene"              ))
  )
p_presence_blocked <- ggplot(coverage_df2, aes(x = sample, y = step, fill = covered)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c(`TRUE` = "#2E7D32", `FALSE` = "#C62828"),
                    labels = c(`TRUE` = "Covered", `FALSE` = "Missing"),
                    name   = "Presence") +
  labs(title = "Precursor Pathways • Step Coverage per Sample (modules separated)",
       x = "Sample", y = "Step") +
  facet_grid(rows = vars(facet), switch = "y", scales = "free_y", space = "free_y") +
  theme_heat +
  theme(strip.background = element_rect(fill = "grey95", color = "grey60"),
        strip.placement  = "outside",
        panel.spacing.y  = unit(6, "pt"))
ggsave(paste(output_dir,"/Pathway_Step_Coverage_per_Sample_byModule.png"), p_presence_blocked, width = 14, height = 11, dpi = 600)

# 4) Macro presence (Upstream/Downstream/Overall)
presence_long <- macro_summary |>
  transmute(sample,
            `Upstream: Mevalonate` = upstream_mev,
            `Upstream: Non-meval.` = upstream_non,
            `Downstream: Alt`      = downstream_any,
            `Overall`              = overall_present) |>
  pivot_longer(-sample, names_to = "track", values_to = "present")
p_macro_presence <- ggplot(presence_long, aes(x = sample, y = track, fill = present)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c(`TRUE` = "#2E7D32", `FALSE` = "#C62828"),
                    name = "Presence", labels = c(`TRUE` = "Yes", `FALSE` = "No")) +
  labs(title = "Isoprenoid Precursor Pathway • Macro Presence per Sample",
       x = "Sample", y = "") +
  theme_heat
ggsave(paste(output_dir,"/Precursor_Macro_Presence_per_Sample.png"), p_macro_presence, width = 12, height = 4.8, dpi = 300)

# 5) Macro activity (Upstream-best / Downstream-best / Overall) — z on samples
activity_long <- macro_summary |>
  transmute(sample,
            `Upstream (sum of Mev+Non)` = upstream_sum,
            `Downstream (best alt)`      = downstream_sum,
            `Overall (bottleneck)`       = overall_activity) |>
  pivot_longer(-sample, names_to = "component", values_to = "activity") |>
  group_by(component) |>
  mutate(activity_z = as.numeric(scale(activity))) |>
  ungroup()

p_macro_activity <- ggplot(activity_long, aes(x = sample, y = component, fill = activity_z)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#D32F2F", midpoint = 0,
                       name = "Activity (z)") +
  labs(title = "Isoprenoid Precursor Pathway • Macro Activity per Sample",
       x = "Sample", y = "") +
  theme_heat
ggsave(paste(output_dir,"/Precursor_Macro_Activity_per_Sample.png"),
       p_macro_activity, width = 12, height = 4.8, dpi = 300)

# -----------------------------
# 👥 Group-level summaries (mean over replicates)
# -----------------------------
print(">>>>__ Group Level Summaries______________ ")
if (!is.null(sample_to_group)) {
  # Step-level presence fraction & mean step score
  group_step_cov <- coverage_df |>
    group_by(Group, pathway, step) |>
    summarise(
      present_frac    = mean(covered,    na.rm = TRUE),
      step_score_mean = mean(step_score, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(step_id = paste(pathway, step, sep = "_"))
  
  p_group_step_presence <- ggplot(group_step_cov, aes(x = Group, y = step_id, fill = present_frac)) +
    geom_tile(color = "white") +
    scale_fill_gradient(limits = c(0, 1), low = "#FBE9E7", high = "#2E7D32",
                        name = "Presence\nfraction") +
    labs(title = "Step Coverage by Group (fraction of replicates present)",
         x = "Group", y = "Pathway • Step") +
    theme_heat +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  ggsave(paste(output_dir,"/Group_Step_Coverage_Fraction.png"), p_group_step_presence, width = 12, height = 8, dpi = 600)
  
  # Module activity (z on group means)
  group_pathway_activity <- pathway_activity |>
    group_by(Group, pathway) |>
    summarise(activity_raw_mean = mean(activity_raw, na.rm = TRUE), .groups = "drop") |>
    group_by(pathway) |>
    mutate(activity_z = as.numeric(scale(activity_raw_mean))) |>
    ungroup()
  
  p_group_activity <- ggplot(group_pathway_activity, aes(x = Group, y = pathway, fill = activity_z)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#D32F2F", midpoint = 0,
                         name = "Activity (z)") +
    labs(title = "Pathway Activity by Group (z-scored on group means)",
         x = "Group", y = "Pathway") +
    theme_heat +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  ggsave(paste(output_dir,"/Group_Pathway_Activity_z.png"), p_group_activity, width = 12, height = 6.5, dpi = 600)
  
  # Macro presence/activity (group level)
  presence_group <- macro_summary |>
    group_by(Group) |>
    summarise(
      `Upstream: Mevalonate` = mean(upstream_mev, na.rm = TRUE),
      `Upstream: Non-meval.` = mean(upstream_non, na.rm = TRUE),
      `Downstream: Alt`      = mean(downstream_any, na.rm = TRUE),
      `Overall`              = mean(overall_present, na.rm = TRUE),
      .groups = "drop"
    ) |>
    pivot_longer(-Group, names_to = "track", values_to = "present_frac")
  
  p_macro_presence_group <- ggplot(presence_group, aes(x = Group, y = track, fill = present_frac)) +
    geom_tile(color = "white") +
    scale_fill_gradient(limits = c(0, 1),
                        labels = scales::label_percent(accuracy = 1),
                        low = "#FFEBEE", high = "#2E7D32",
                        name = "Presence") +
    labs(title = "Precursor Macro Presence • Group-level (share of present replicates)",
         x = "Group", y = "") +
    theme_heat +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  ggsave(paste(output_dir,"/Group_Macro_Presence.png"), p_macro_presence_group, width = 12, height = 4.8, dpi = 300)
  
  activity_group <- macro_summary |>
    group_by(Group) |>
    summarise(
      `Upstream (sum of Mev+Non)` = mean(upstream_sum,       na.rm = TRUE),
      `Downstream (best alt)`     = mean(downstream_sum, na.rm = TRUE),
      `Overall (bottleneck)`      = mean(overall_activity,    na.rm = TRUE),
      .groups = "drop"
    ) |>
    pivot_longer(-Group, names_to = "component", values_to = "activity_raw_mean") |>
    group_by(component) |>
    mutate(activity_z = as.numeric(scale(activity_raw_mean))) |>
    ungroup()
  
  p_macro_activity_group <- ggplot(activity_group, aes(x = Group, y = component, fill = activity_z)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#D32F2F", midpoint = 0,
                         name = "Activity (z)") +
    labs(title = "Precursor Macro Activity • Group-level (z on group means)",
         x = "Group", y = "") +
    theme_heat +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  ggsave(paste(output_dir,"/Group_Macro_Activity.png"), p_macro_activity_group, width = 12, height = 4.8, dpi = 600)
}

# -----------------------------
# 💾 Optional: export tables
# -----------------------------
if (requireNamespace("openxlsx", quietly = TRUE)) {
  openxlsx::write.xlsx(
    list(
      StepCoverage         = coverage_df %>% select(sample, Group, pathway, step, covered, step_score),
      PathwayActivity      = pathway_activity,
      MacroSummary_Sample  = macro_summary,
      Group_Step_Coverage  = if (exists("group_step_cov")) group_step_cov else NULL,
      Group_Pathway_Act    = if (exists("group_pathway_activity")) group_pathway_activity else NULL,
      Group_Macro_Presence = if (exists("presence_group")) presence_group else NULL,
      Group_Macro_Activity = if (exists("activity_group")) activity_group else NULL
    ),
    file =   file.path(output_dir,"Precursor_Pathways_Summaries.xlsx"),
    overwrite = TRUE
  )
}

# ============================================
# 🔁 Make pathway activity consistent with "Groupwise module activity"
#    -> same compute_module_activity(), same source (vsd_mat_for_activity)
# ============================================
print(">>>>__ Makw oathway activity consistent... ... ... ______________ ")
# 1) Build gene sets per pathway from rules (union of all KOs in that pathway)
extract_genes_from_rule <- function(rule) {
  if (is.na(rule) || !nzchar(rule)) return(character())
  or_parts  <- strsplit(rule, "\\|", perl = TRUE)[[1]]
  genes_vec <- unlist(strsplit(or_parts, "&", perl = TRUE), use.names = FALSE)
  trimws(unique(genes_vec))
}

pathway_gene_sets <- pathway_steps |>
  dplyr::group_by(pathway) |>
  dplyr::summarise(
    genes = list(unique(unlist(lapply(rule, extract_genes_from_rule)))),
    .groups = "drop"
  )

# 2) Choose the same matrix used in your Groupwise block
#    (there you use: vsd_mat_for_activity <- as.matrix(plot_matrix))
activity_src <- if (exists("vsd_mat_for_activity")) {
  as.matrix(vsd_mat_for_activity)
} else if (exists("plot_matrix")) {
  as.matrix(plot_matrix)
} else {
  as.matrix(vsd_mat)  # fallback
}

# 3) Compute activity per pathway using your compute_module_activity()
#    (this respects USE_ZSCORE_FOR_ACTIVITY and safe_scale_rows you defined earlier)
pathway_activity_matrix <- sapply(
  setNames(pathway_gene_sets$genes, pathway_gene_sets$pathway),
  compute_module_activity,
  vsd_matrix = activity_src
)

# pathway_activity_df: long form, with z-score per pathway (like your plots)
pathway_activity_df <- as.data.frame(pathway_activity_matrix)
pathway_activity_df$Sample <- colnames(activity_src)
if (exists("col_data") && all(pathway_activity_df$Sample %in% rownames(col_data))) {
  pathway_activity_df$Group <- as.character(col_data[pathway_activity_df$Sample, "condition"])
} else {
  pathway_activity_df$Group <- NA_character_
}

pathway_activity_long <- tidyr::pivot_longer(
  pathway_activity_df,
  cols = -c(Sample, Group),
  names_to = "pathway",
  values_to = "activity_raw"
) |>
  dplyr::group_by(pathway) |>
  dplyr::mutate(activity_z = as.numeric(scale(activity_raw))) |>
  dplyr::ungroup()

# 4) Re-plot per-sample pathway activity (consistent now)
p_activity_cons <- ggplot(pathway_activity_long,
                          aes(x = Sample, y = pathway, fill = activity_z)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#D32F2F",
                       midpoint = 0, name = "Activity (z)") +
  labs(title = "Precursor Pathways • Activity per Sample (z-scored, consistent)",
       x = "Sample", y = "Pathway") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste(output_dir,"/Pathway_Activity_per_Sample_consistent.png"),
       p_activity_cons, width = 13, height = 6.5, dpi = 600)

# 5) Macro activity rebuilt from these consistent pathway activities
macro_from_consistent <- tidyr::pivot_wider(
  pathway_activity_long,
  id_cols = c(Sample, Group),
  names_from = pathway,
  values_from = activity_raw,
  values_fill = 0
) |>
  dplyr::mutate(
    # Upstream: sum of both mevalonate branches + non-mev
    mev_activity_sum     = coalesce(`M00095`, 0) + coalesce(`M00849`, 0),
    non_mev_activity     = coalesce(`M00096`, 0),
    upstream_sum         = mev_activity_sum + non_mev_activity,
    
    # Downstream: sum of all alternative downstream modules
    downstream_sum       = rowSums(cbind(
      coalesce(`M00364`, 0),
      coalesce(`M00365`, 0),
      coalesce(`M00366`, 0),
      coalesce(`M00367`, 0)
    ), na.rm = TRUE),
    
    # Overall bottleneck still limited by the smaller of upstream vs downstream
    overall_activity     = pmin(upstream_sum, downstream_sum, na.rm = TRUE)
  )

# sample-level macro heatmap (z per component)
activity_long_cons <- macro_from_consistent |>
  tidyr::pivot_longer(-c(Sample, Group),
                      names_to = "component", values_to = "activity") |>
  dplyr::filter(component %in% c("upstream_sum","downstream_sum","overall_activity")) |>
  dplyr::mutate(component = dplyr::recode(component,
                                          upstream_sum     = "Upstream (sum of Mev + Non-mev)",
                                          downstream_sum   = "Downstream (sum of alts)",
                                          overall_activity = "Overall (bottleneck)"
  )) |>
  dplyr::group_by(component) |>
  dplyr::mutate(activity_z = as.numeric(scale(activity))) |>
  dplyr::ungroup()

p_macro_activity_cons <- ggplot(activity_long_cons,
                                aes(x = Sample, y = component, fill = activity_z)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#D32F2F",
                       midpoint = 0, name = "Activity (z)") +
  labs(title = "Isoprenoid Precursor Pathway • Macro Activity per Sample (consistent)",
       x = "Sample", y = "") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste(output_dir,"/Precursor_Macro_Activity_per_Sample_consistent.png"),
       p_macro_activity_cons, width = 12, height = 4.8, dpi = 300)

# 6) Group-level (mean of replicates) using consistent activities
if (!all(is.na(pathway_activity_df$Group))) {
  # module level
  group_pathway_activity_cons <- pathway_activity_long |>
    dplyr::group_by(Group, pathway) |>
    dplyr::summarise(activity_raw_mean = mean(activity_raw, na.rm = TRUE), .groups = "drop") |>
    dplyr::group_by(pathway) |>
    dplyr::mutate(activity_z = as.numeric(scale(activity_raw_mean))) |>
    dplyr::ungroup()
  
  p_group_activity_cons <- ggplot(group_pathway_activity_cons,
                                  aes(x = Group, y = pathway, fill = activity_z)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#D32F2F",
                         midpoint = 0, name = "Activity (z)") +
    labs(title = "Pathway Activity by Group (z-scored on group means, consistent)",
         x = "Group", y = "Pathway") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  ggsave(paste(output_dir,"/Group_Pathway_Activity_z_consistent.png"),
         p_group_activity_cons, width = 12, height = 6.5, dpi = 300)
  
  # macro level (group mean)
  macro_group_cons <- macro_from_consistent |>
    dplyr::group_by(Group) |>
    dplyr::summarise(
      `Upstream (sum of Mev + Non-mev)` = mean(upstream_sum,   na.rm = TRUE),
      `Downstream (sum of alts)`        = mean(downstream_sum, na.rm = TRUE),
      `Overall (bottleneck)`            = mean(overall_activity, na.rm = TRUE),
      .groups = "drop"
    ) |>
    tidyr::pivot_longer(-Group, names_to = "component", values_to = "activity_raw_mean") |>
    dplyr::group_by(component) |>
    dplyr::mutate(activity_z = as.numeric(scale(activity_raw_mean))) |>
    dplyr::ungroup()
  
  p_macro_activity_group_cons <- ggplot(macro_group_cons,
                                        aes(x = Group, y = component, fill = activity_z)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#D32F2F",
                         midpoint = 0, name = "Activity (z)") +
    labs(title = "Precursor Macro Activity • Group-level (z on group means, consistent)",
         x = "Group", y = "") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  ggsave(paste(output_dir,"/Group_Macro_Activity_consistent.png"),
         p_macro_activity_group_cons, width = 12, height = 4.8, dpi = 300)
}

# =========================================================
# 🟢 Bubble chart with expression-driven sizes (robust)
# =========================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(scales)
})

# ---- knobs ----
SIZE_MODE        <- "best_combo_min"   # "fraction" | "sum_expr_union" | "mean_expr_union" | "best_combo_min" | "best_combo_sum"
PRESENCE_THRESH  <- get0("PRESENCE_THRESH", ifnotfound = 0)
MIN_PRESENT_SIZE <- 0.20               # حداقل اندازه برای present (0..1)
ABSENT_SIZE      <- 0.06               # اندازه ثابت برای absent (0..1)
USE_LOG_SIZE     <- TRUE               # برای اندازه از log1p استفاده شود؟

# ---- choose expression source (same as your activity source) ----
expr_src <- if (exists("vsd_mat_for_activity")) {
  as.matrix(vsd_mat_for_activity)
} else if (exists("plot_matrix")) {
  as.matrix(plot_matrix)
} else {
  as.matrix(vsd_mat)
}
storage.mode(expr_src) <- "numeric"

# ---- helpers ----
split_or  <- function(x) strsplit(x, "\\|", perl = TRUE)[[1]] |> trimws()
split_and <- function(x) strsplit(x, "&",  perl = TRUE)[[1]] |> trimws()
extract_union_genes <- function(rule) {
  if (is.na(rule) || !nzchar(rule)) return(character(0))
  parts <- unlist(strsplit(gsub("[&|]", ",", rule), ",", fixed = FALSE), use.names = FALSE)
  trimws(unique(parts[nchar(parts) > 0]))
}
compute_step_expr_score <- function(rule, sample_name, mode = "best_combo_min") {
  if (is.na(rule) || !nzchar(rule) || !sample_name %in% colnames(expr_src)) return(0)
  sv <- expr_src[, sample_name, drop = TRUE]; names(sv) <- rownames(expr_src)
  if (mode %in% c("best_combo_min", "best_combo_sum")) {
    combos <- split_or(rule); if (!length(combos)) return(0)
    cs <- vapply(combos, function(cb) {
      g <- intersect(split_and(cb), names(sv)); if (!length(g)) return(0)
      v <- sv[g]
      if (mode == "best_combo_min") min(v, na.rm = TRUE) else sum(v, na.rm = TRUE)
    }, numeric(1))
    return(max(cs, na.rm = TRUE))
  } else {
    g <- intersect(extract_union_genes(rule), names(sv)); if (!length(g)) return(0)
    v <- sv[g]; if (mode == "sum_expr_union") sum(v, na.rm = TRUE) else mean(v, na.rm = TRUE)
  }
}

# ---- USE group_list samples (robust) ----
stopifnot(exists("col_data"))
if (exists("group_list") && any(col_data$condition %in% group_list)) {
  selected_samples <- rownames(col_data)[col_data$condition %in% group_list]
} else if (exists("group_list") && all(group_list %in% rownames(col_data))) {
  selected_samples <- intersect(group_list, rownames(col_data))
} else {
  selected_samples <- colnames(vsd_mat)
}
# تقاطع با هر دو ماتریس تا خطای subscript نگیریم
selected_samples <- intersect(selected_samples, intersect(colnames(vsd_mat), colnames(expr_src)))
if (!length(selected_samples)) {
  stop("No overlap between selected samples and matrices. Check group_list / sample names.")
}
# ترتیب: بر اساس ترتیب گروه‌ها و سپس نام نمونه
if (exists("group_list")) {
  cf <- factor(col_data[selected_samples, "condition"], levels = group_list)
  ord <- order(cf, selected_samples)
  selected_samples <- selected_samples[ord]
}

# زیرمجموعه‌سازی همهٔ منابع
vsd_mat    <- vsd_mat[, selected_samples, drop = FALSE]
expr_src   <- expr_src[, selected_samples, drop = FALSE]
present_mat <- vsd_mat > PRESENCE_THRESH
if (exists("col_data")) col_data <- col_data[selected_samples, , drop = FALSE]

# ---- step gene unions ----
step_union_genes <- pathway_steps %>% mutate(genes_union = lapply(rule, extract_union_genes))

# ---- build bubble_expr_df from scratch (includes expr_score) ----
bubble_expr_df <- tidyr::expand_grid(
  sample = selected_samples,
  row_id = seq_len(nrow(step_union_genes))
) %>%
  mutate(
    pathway = step_union_genes$pathway[row_id],
    step    = step_union_genes$step[row_id],
    rule    = pathway_steps$rule[row_id],
    genes_u = step_union_genes$genes_union[row_id]
  ) %>%
  rowwise() %>%
  mutate(
    # fraction present (for reference)
    present_n = {
      g <- intersect(genes_u, rownames(present_mat))
      if (!length(g)) 0L else sum(present_mat[g, sample, drop = TRUE], na.rm = TRUE)
    },
    total_n   = length(genes_u),
    frac_pres = ifelse(total_n > 0, present_n / total_n, 0),
    # expression-based score for bubble sizing
    expr_score = compute_step_expr_score(rule, sample, mode = SIZE_MODE)
  ) %>%
  ungroup() %>%
  select(sample, pathway, step, frac_pres, expr_score) %>%
  left_join(coverage_df %>% select(sample, pathway, step, covered),
            by = c("sample", "pathway", "step")) %>%
  mutate(
    # order X as selected_samples
    sample = factor(sample, levels = selected_samples),
    # facet text
    block = module_block[pathway],
    facet = factor(paste0(pathway, " • ", block),
                   levels = c("M00095",
                              "M00849",
                              "M00096",
                              "M00364",
                              "M00365",
                              "M00366",
                              "M00367","M00927")),
    step = factor(step, levels = sort(unique(step)))
  )

# ---- scale to sizes & colors ----
expr_pos <- pmax(bubble_expr_df$expr_score, 0)
expr_scaled <- if (USE_LOG_SIZE) log1p(expr_pos) else expr_pos
expr_scaled <- if (max(expr_scaled, na.rm = TRUE) > 0) expr_scaled / max(expr_scaled, na.rm = TRUE) else expr_scaled
bubble_expr_df$expr_scaled <- expr_scaled

pretty_names <- c(
  M00095 = "M00095",
  M00849 = "M00849",
  M00096 = "M00096",
  M00364 = "M00364",
  M00365 = "M00365",
  M00366 = "M00366",
  M00367 = "M00367",
  M00927 = "M00927",
  Taxadiene="Taxadiene"
)

bubble_expr_df <- bubble_expr_df %>%
  mutate(
    # برچسب پیشنهادی بر اساس نگاشت
    facet_try = unname(pretty_names[as.character(pathway)]),
    # اگر نگاشت چیزی نداد، برچسب پیش‌فرض بساز
    facet_lbl = ifelse(
      is.na(facet_try) | facet_try == "",
      paste0(pathway, " • ", module_block[pathway]),
      facet_try
    ),
    facet = factor(
      facet_lbl,
      levels = c(
        pretty_names["M00095"],
        pretty_names["M00849"],
        pretty_names["M00096"],
        pretty_names["M00364"],
        pretty_names["M00365"],
        pretty_names["M00366"],
        pretty_names["M00367"],
        pretty_names["M00927"] 
        
        
      )
    )
  )

# 1) join 'covered' if missing
if (!"covered" %in% names(bubble_expr_df)) {
  bubble_expr_df <- bubble_expr_df %>%
    dplyr::left_join(coverage_df %>% dplyr::select(sample, pathway, step, covered),
                     by = c("sample","pathway","step"))
}

# 2) safe logical flag for presence
bubble_expr_df <- bubble_expr_df %>%
  dplyr::mutate(covered2 = !is.na(covered) & as.logical(covered))

# 3) (re)compute scaled expression 0..1 if missing
if (!"expr_scaled" %in% names(bubble_expr_df)) {
  expr_pos <- pmax(bubble_expr_df$expr_score, 0)
  expr_scaled <- if (exists("USE_LOG_SIZE") && USE_LOG_SIZE) log1p(expr_pos) else expr_pos
  expr_scaled <- if (max(expr_scaled, na.rm = TRUE) > 0) expr_scaled / max(expr_scaled, na.rm = TRUE) else expr_scaled
  bubble_expr_df$expr_scaled <- expr_scaled
}

# 4) compute size & fill aesthetics
if (!exists("MIN_PRESENT_SIZE")) MIN_PRESENT_SIZE <- 0.20
if (!exists("ABSENT_SIZE"))      ABSENT_SIZE      <- 0.06

bubble_expr_df <- bubble_expr_df %>%
  dplyr::mutate(
    size_present = MIN_PRESENT_SIZE + (1 - MIN_PRESENT_SIZE) * expr_scaled,
    size_var     = ifelse(covered2, size_present, ABSENT_SIZE),
    fill_val     = ifelse(covered2, expr_scaled, NA_real_)
  )

# ---- plot ----
p_bubble_expr <- ggplot(bubble_expr_df, aes(x = sample, y = step)) +
  geom_point(aes(size = size_var, fill = fill_val),
             shape = 21, color = "grey25", alpha = 0.9) +
  scale_fill_gradient(low = "#C8E6C9", high = "#1B5E20",
                      na.value = "#C62828", name = "Expression (scaled)") +
  scale_size_area(limits = c(ABSENT_SIZE, 1), max_size = 9,
                  breaks = c(ABSENT_SIZE, MIN_PRESENT_SIZE, (1+MIN_PRESENT_SIZE)/2, 1),
                  labels = c("Absent", "Present (min)", "mid", "max"),
                  name = "Bubble size") +
  scale_y_discrete(limits = function(x) rev(x)) + 
  labs(title = paste0("Precursor Pathways • Step Coverage per Sample (bubble, ", SIZE_MODE, ")"),
       x = "Sample", y = "Step") +
  facet_grid(rows = vars(facet), switch = "y", scales = "free_y", space = "free_y") +
  theme_bw(base_size = 12) +
  theme(
    strip.background   = element_rect(fill = "grey95", color = "grey60"),
    strip.placement    = "outside",
    strip.text.y.left  = element_text(angle = 0, hjust = 1,
                                      margin = margin(r = 8), size = 11),
    strip.clip         = "off",
    plot.margin        = margin(t = 10, r = 20, b = 10, l = 140),
    axis.text.x        = element_text(angle = 45, hjust = 1),
    panel.spacing.y    = unit(6, "pt"),
    legend.box         = "vertical",
    legend.spacing.y   = unit(3, "pt")
  ) +
  coord_cartesian(clip = "off")

 
ggsave(sprintf(paste(output_dir,"/Pathway_Step_Coverage_Bubble_byModule_%s.png"), SIZE_MODE),
       p_bubble_expr, width = 15, height = 11, dpi = 600)

message("✅ Precursor pathway analysis finished.")

#========================================================= 
# 📈 Multi-module bubble chart using a custom steps list 
# steps_df: data.frame(pathway, step, rule) # modules : character vector of pathway IDs to include 
# =========================================================

plot_modules_bubble <- function(modules,
                                steps_df,
                                size_mode        = get0("SIZE_MODE", ifnotfound = "best_combo_min"),
                                # size controls
                                size_strategy    = c("raw", "gamma", "zstep"),
                                size_gamma       = 0.6,
                                trim_quant       = c(0.05, 0.95),
                                rescale_within   = TRUE,   # scale within each module (facet)
                                use_log_size     = TRUE,   # log1p قبل از rescale
                                min_present_size = 0.20,
                                absent_size      = 0.06,
                                max_size_pt      = 9,
                                # matrices / samples
                                vsd_mat          = NULL,
                                expr_mat         = NULL,
                                samples          = NULL,
                                # --- 🆕 grouping / averaging ---
                                aggregate_groups = FALSE,
                                group_var        = "condition",
                                groups           = NULL,         # زیرمجموعه گروه‌ها (اختیاری)
                                group_fun        = c("mean","median"),
                                presence_rule    = c("any","majority","all"),
                                # labels & output
                                module_labels    = NULL,   # named chr: c(M00849="…", Taxadiene="…")
                                outdir           = "Plots",
                                file_tag         = NULL) {
  size_strategy <- match.arg(size_strategy)
  group_fun     <- match.arg(group_fun)
  presence_rule <- match.arg(presence_rule)
  
  stopifnot(is.data.frame(steps_df))
  stopifnot(all(c("pathway","step","rule") %in% names(steps_df)))
  
  PRESENCE_THRESH <- get0("PRESENCE_THRESH", ifnotfound = 0)
  
  # ---------- resolve matrices ----------
  resolve_matrix <- function(...) { for (m in list(...)) if (!is.null(m)) return(as.matrix(m)); NULL }
  if (is.null(vsd_mat)) {
    vsd_obj <- get0("vsd", ifnotfound = NULL)
    vsd_from_vsd <- if (!is.null(vsd_obj)) {
      tryCatch(SummarizedExperiment::assay(vsd_obj),
               error = function(e) tryCatch(as.matrix(vsd_obj), error = function(e) NULL))
    } else NULL
    vsd_mat <- resolve_matrix(get0("vsd_mat", ifnotfound = NULL),
                              vsd_from_vsd,
                              get0("plot_matrix", ifnotfound = NULL),
                              get0("vsd_mat_for_activity", ifnotfound = NULL))
  } else vsd_mat <- as.matrix(vsd_mat)
  if (is.null(vsd_mat)) stop("Cannot resolve 'vsd_mat'."); storage.mode(vsd_mat) <- "numeric"
  
  if (is.null(expr_mat)) {
    expr_mat <- resolve_matrix(get0("vsd_mat_for_activity", ifnotfound = NULL),
                               get0("plot_matrix", ifnotfound = NULL),
                               vsd_mat)
  } else expr_mat <- as.matrix(expr_mat)
  storage.mode(expr_mat) <- "numeric"
  
  # ---------- choose modules & steps ----------
  modules <- unique(as.character(modules))
  steps_df <- steps_df |>
    dplyr::filter(.data$pathway %in% modules) |>
    dplyr::mutate(pathway = as.character(.data$pathway),
                  step    = as.character(.data$step),
                  rule    = as.character(.data$rule))
  if (nrow(steps_df) == 0) stop("No steps found for the requested modules.")
  
  # ---------- samples OR groups ----------
  col_data <- get0("col_data", ifnotfound = NULL)
  
  if (!aggregate_groups) {
    # مثل قبل: محور x = نمونه‌ها
    if (is.null(samples)) {
      group_list <- get0("group_list", ifnotfound = NULL)
      if (!is.null(group_list) && !is.null(col_data) && group_var %in% names(col_data)) {
        ss <- rownames(col_data)[col_data[[group_var]] %in% group_list]
        samples <- intersect(ss, intersect(colnames(vsd_mat), colnames(expr_mat)))
      }
      if (is.null(samples) || length(samples) == 0L) {
        samples <- intersect(colnames(vsd_mat), colnames(expr_mat))
      }
    } else {
      samples <- intersect(samples, intersect(colnames(vsd_mat), colnames(expr_mat)))
    }
    if (!length(samples)) stop("No overlapping samples between matrices.")
    vsd_mat  <- vsd_mat[, samples, drop = FALSE]
    expr_mat <- expr_mat[, samples, drop = FALSE]
  } else {
    # محور x = گروه‌ها (میانگین/مدین تکرارها)
    if (is.null(col_data) || !group_var %in% names(col_data))
      stop("For aggregate_groups=TRUE you must have 'col_data' with column '", group_var, "'.")
    if (is.null(groups)) {
      groups <- get0("group_list", ifnotfound = NULL)
      if (is.null(groups)) groups <- unique(col_data[[group_var]])
    }
    groups <- intersect(as.character(groups), as.character(col_data[[group_var]]))
    if (!length(groups)) stop("No groups found in col_data$", group_var, " matching requested 'groups'.")
    
    # نگاشت sample -> group
    sample2group <- setNames(as.character(col_data[colnames(vsd_mat), group_var]),
                             colnames(vsd_mat))
    keep_smpl <- names(sample2group)[sample2group %in% groups]
    if (!length(keep_smpl)) stop("No samples in requested groups.")
    vsd_mat  <- vsd_mat[, keep_smpl, drop = FALSE]
    expr_mat <- expr_mat[, keep_smpl, drop = FALSE]
  }
  
  present_mat <- vsd_mat > PRESENCE_THRESH
  
  # ---------- helpers ----------
  split_or  <- function(x) strsplit(x, "\\|", perl = TRUE)[[1]] |> trimws()
  split_and <- function(x) strsplit(x, "&",  perl = TRUE)[[1]] |> trimws()
  extract_union_genes <- function(rule) {
    if (is.na(rule) || !nzchar(rule)) return(character(0))
    parts <- unlist(strsplit(gsub("[&|]", ",", rule), ","))
    trimws(unique(parts[nchar(parts) > 0]))
  }
  eval_step_presence <- function(rule, smp) {
    if (is.na(rule) || !nzchar(rule)) return(FALSE)
    combos <- split_or(rule)
    any(vapply(combos, function(cb) {
      g <- intersect(split_and(cb), rownames(vsd_mat))
      length(g) > 0 && all(vsd_mat[g, smp, drop = TRUE] > PRESENCE_THRESH, na.rm = TRUE)
    }, logical(1)))
  }
  compute_step_expr_score <- function(rule, smp, mode = "best_combo_min") {
    if (is.na(rule) || !nzchar(rule)) return(0)
    sv <- expr_mat[, smp, drop = TRUE]; names(sv) <- rownames(expr_mat)
    if (mode %in% c("best_combo_min", "best_combo_sum")) {
      combos <- split_or(rule); if (!length(combos)) return(0)
      cs <- vapply(combos, function(cb) {
        g <- intersect(split_and(cb), names(sv)); if (!length(g)) return(0)
        v <- sv[g]; if (mode == "best_combo_min") min(v, na.rm = TRUE) else sum(v, na.rm = TRUE)
      }, numeric(1))
      max(cs, na.rm = TRUE)
    } else {
      g <- intersect(extract_union_genes(rule), names(sv)); if (!length(g)) return(0)
      v <- sv[g]; if (mode == "sum_expr_union") sum(v, na.rm = TRUE) else mean(v, na.rm = TRUE)
    }
  }
  
  # ---------- expand all rows ----------
  steps_df$genes_union <- lapply(steps_df$rule, extract_union_genes)
  
  base_df <- tidyr::expand_grid(sample = colnames(vsd_mat), row_id = seq_len(nrow(steps_df))) |>
    dplyr::mutate(
      pathway = steps_df$pathway[row_id],
      step    = steps_df$step[row_id],
      rule    = steps_df$rule[row_id],
      genes_u = steps_df$genes_union[row_id]
    ) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      present_n = { g <- intersect(genes_u, rownames(present_mat));
      if (!length(g)) 0L else sum(present_mat[g, sample, drop = TRUE], na.rm = TRUE) },
      total_n   = length(genes_u),
      frac_pres = ifelse(total_n > 0, present_n / total_n, 0),
      expr_score = compute_step_expr_score(rule, sample, mode = size_mode),
      covered    = eval_step_presence(rule, sample)
    ) |>
    dplyr::ungroup()
  
  # ---------- (اختیاری) تجمیع گروهی ----------
  if (aggregate_groups) {
    # اضافه کردن نام گروه به ردیف‌ها
    base_df$group <- sample2group[base_df$sample]
    
    agg_fun <- if (group_fun == "mean") mean else stats::median
    
    df <- base_df |>
      dplyr::group_by(pathway, step, group) |>
      dplyr::summarise(
        expr_score = agg_fun(expr_score, na.rm = TRUE),
        # proportion of covered replicates:
        prop_cov   = mean(covered, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::mutate(
        covered = dplyr::case_when(
          presence_rule == "all"      ~ prop_cov >= 1,
          presence_rule == "majority" ~ prop_cov >= 0.5,
          TRUE                        ~ prop_cov > 0        # "any"
        ),
        sample = group                # محور x = نام گروه
      )
  } else {
    df <- base_df
  }
  
  # ---------- scale expression to [0,1] ----------
  expr_pos <- pmax(df$expr_score, 0)
  if (use_log_size) expr_pos <- log1p(expr_pos)
  
  if (rescale_within) {
    df$expr_scaled <- dplyr::group_by(df, pathway) |>
      dplyr::mutate(expr_scaled = {
        r <- range(expr_pos[dplyr::cur_group_rows()], na.rm = TRUE)
        if (is.finite(r[1]) && is.finite(r[2]) && diff(r) > 0)
          (expr_pos[dplyr::cur_group_rows()] - r[1]) / diff(r)
        else
          expr_pos[dplyr::cur_group_rows()] * 0
      }) |>
      dplyr::ungroup() |>
      dplyr::pull(expr_scaled)
  } else {
    mx <- max(expr_pos, na.rm = TRUE)
    df$expr_scaled <- if (mx > 0) expr_pos / mx else expr_pos
  }
  
  # ---------- size strategies ----------
  if (size_strategy == "raw") {
    size_base <- df$expr_scaled
    
  } else if (size_strategy == "gamma") {
    size_base <- dplyr::group_by(df, pathway) |>
      dplyr::mutate(tmp_es = {
        es <- df$expr_scaled[dplyr::cur_group_rows()]
        es_pres <- es[df$covered[dplyr::cur_group_rows()]]
        if (length(na.omit(es_pres)) > 2) {
          q <- stats::quantile(es_pres, trim_quant, na.rm = TRUE)
          es_cs <- scales::rescale(pmin(pmax(es, q[1]), q[2]), to = c(0, 1), from = q)
        } else {
          es_cs <- es
        }
        es_cs ^ size_gamma
      }) |>
      dplyr::ungroup() |>
      dplyr::pull(tmp_es)
    
  } else if (size_strategy == "zstep") {
    z_by_step <- ave(df$expr_score, interaction(df$pathway, df$step), FUN = function(v) as.numeric(scale(v)))
    z_by_step[is.na(z_by_step)] <- 0
    z_pos <- pmax(z_by_step, 0)
    size_base <- if (max(z_pos, na.rm = TRUE) > 0) z_pos / max(z_pos, na.rm = TRUE) else z_pos
  }
  
  df <- df |>
    dplyr::mutate(
      size_present = min_present_size + (1 - min_present_size) * size_base,
      size_var     = ifelse(covered, size_present, absent_size),
      fill_val     = ifelse(covered, expr_scaled, NA_real_)
    )
  
  # ---------- facet labels ----------
  if (is.null(module_labels)) module_labels <- setNames(modules, modules)
  facet_levels <- unname(module_labels[modules])
  
  # ترتیب محور x: نمونه‌ها یا گروه‌ها
  x_levels <- if (aggregate_groups) {
    # به ترتیب 'groups' بده
    groups
  } else {
    colnames(vsd_mat)
  }
  
  df <- df |>
    dplyr::mutate(
      facet  = factor(module_labels[pathway], levels = facet_levels),
      step   = factor(step, levels = stats::na.omit(unique(step[order(pathway, step)]))),
      sample = factor(sample, levels = x_levels)
    )
  
  # ---------- plot ----------
  p <- ggplot2::ggplot(df, ggplot2::aes(x = sample, y = step)) +
    ggplot2::geom_point(ggplot2::aes(size = size_var, fill = fill_val),
                        shape = 21, color = "grey25", alpha = 0.9) +
    ggplot2::scale_fill_gradient(low = "#C8E6C9", high = "#1B5E20",
                                 na.value = "#C62828", name = "Expression (scaled)") +
    ggplot2::scale_size_area(limits = c(absent_size, 1), max_size = max_size_pt,
                             breaks = c(absent_size, min_present_size, (1 + min_present_size) / 2, 1),
                             labels = c("Absent", "Present (min)", "mid", "max"),
                             name = "Bubble size") +
    scale_y_discrete(limits = function(x) rev(x)) +
    ggplot2::labs(
      title = sprintf("Step Coverage per Sample (bubble) • %s • %s%s",
                      size_mode, size_strategy,
                      if (aggregate_groups) " • grouped means" else ""),
      x = if (aggregate_groups) "Group" else "Sample",
      y = "Step"
    ) +
    ggplot2::facet_grid(rows = ggplot2::vars(facet), switch = "y",
                        scales = "free_y", space = "free_y", drop = FALSE) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      strip.background   = ggplot2::element_rect(fill = "grey95", color = "grey60"),
      strip.placement    = "outside",
      strip.text.y.left  = ggplot2::element_text(angle = 0, hjust = 1,
                                                 margin = ggplot2::margin(r = 8), size = 11),
      strip.clip         = "off",
      plot.margin        = ggplot2::margin(t = 10, r = 20, b = 10, l = 140),
      axis.text.x        = ggplot2::element_text(angle = 45, hjust = 1),
      panel.spacing.y    = grid::unit(6, "pt"),
      legend.box         = "vertical",
      legend.spacing.y   = grid::unit(3, "pt")
    ) +
    ggplot2::coord_cartesian(clip = "off")
  
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  if (is.null(file_tag)) {
    tag <- paste0(modules[1], if (length(modules) > 1) paste0("_plus", length(modules) - 1) else "")
  } else {
    tag <- file_tag
  }
  file <- file.path(outdir, sprintf("Modules_Bubble_%s_%s_%s%s.png",
                                    tag, size_mode, size_strategy,
                                    if (aggregate_groups) "_grouped" else ""))
  ggplot2::ggsave(file, p, width = 12, height = 1.8 * length(modules) + 5, dpi = 300)
  message("✅ Saved: ", file)
  
  invisible(p)
}

# -----------------------------
# 🎯 Example usage of plot_modules_bubble function
# -----------------------------

# دو ماژول کنار هم، با تقویت اختلاف‌ها
exists("vsd_mat"); exists("vsd"); exists("plot_matrix"); exists("vsd_mat_for_activity")
if (exists("vsd")) print(dim(SummarizedExperiment::assay(vsd)))

plot_modules_bubble(modules = c("M00095","M00849","M00096","M00366"), steps_df = custom_steps,
                    module_labels  = labels,
                    size_strategy  = "gamma", size_gamma = 0.66, trim_quant = c(0.1, 0.9),
                    rescale_within = TRUE,
                  #  aggregate_groups = TRUE,
                  #  group_var = "condition",
                  #  groups    = names(group_list),
                  #  presence_rule = "any",
                  #  group_fun = "mean"
)
plot_modules_bubble(modules = c("M00927","Taxadiene"), steps_df = custom_steps,
                    module_labels  = labels,
                    size_strategy  = "gamma", size_gamma = 0.66, trim_quant = c(0.1, 0.9),
                    rescale_within = TRUE,
                    #  aggregate_groups = TRUE,
                    #  group_var = "condition",
                    #  groups    = names(group_list),
                    #  presence_rule = "any",
                    #  group_fun = "mean"
)
message("🎉 Done.")
