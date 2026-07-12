# Mission 01 : Data Preparation
# --------------------------------------------------------
library(DESeq2)
library(ggplot2)
library(pheatmap)

# ===============================
#   Langkah 1: Download Data
# ===============================
# 1. Load package
library(GEOquery)

library(dplyr)
library(tidyr)
library(readr)


# 2. Download series_matrix GSE 
gse <- getGEO("GSE285595", GSEMatrix=TRUE)

# 3. Download semua file supplementary dari GSE
getGEOSuppFiles("GSE285595")

# ===============================
#   Langkah 2: Data Cleaning
# ===============================
# 1. Ambil data ekspresi dan metadata dari series matrix
gse_data <- gse[[1]]
expression_data <-exprs(gse_data)
metadata <- pData(gse_data)

# 2. Cek struktur datanya
metadata[1]
metadata[2]
metadata[3]
metadata[4]
length(metadata)
dim(expression_data)

# 3. Lihat list supplementary file yang sudah didownload
list.files("GSE285595/")

# 4. Baca isi file, dan simpan sebagai objek “counts”. 
# Sesuaikan nama file dengan output list.files tadi
counts <- read_tsv("GSE285595/GSE285595_Expression_Values_Raw_Data.txt.gz")

# 5. Lihat struktur data. Seharusnya jumlah baris bukan 0
dim(counts)
head(counts)
colnames(counts)

# 6. Ubah counts menjadi format data frame
counts <- as.data.frame(counts)

# 7. Ambil isi kolom pertama (yang berisi ID gen) dan 
# jadikan sebagai nama baris (rownames). Lalu cek strukturnya.
head(counts)
counts[,1]
rownames(counts) <- counts[,1]
counts <- counts[,-1]
dim(counts)

# ===============================
#   Langkah 3: Membuat Metadata
# ===============================
colnames(counts)

sample_names <- colnames(counts)
condition <- c("DMSO", "DMSO", "DTG", "DTG", "DTG", "DMSO")

metadata <- data.frame(
  sample    = sample_names,
  condition = condition 
)

rownames(metadata) <- sample_names

print(metadata)

#  Verifikasi dengan
#  colnames(counts) == rownames(metadata)
colnames(counts) == rownames(metadata)

# ===============================
#   Langkah 4: Distribusi Data
# ===============================
# 1. Load library untuk visualisasi
library(ggplot2)
library(tidyr)
library(dplyr)

# 2. Ubah data ke format panjang (long format) untuk ggplot
counts_long <- counts %>% 
  as.data.frame() %>% mutate(gene = rownames(.)) %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "count")

# 3. Boxplot (dalam skala log2 agar lebih jelas)
ggplot(counts_long, aes(x=sample, y=log2(count+1), fill = sample)) +
  geom_boxplot() +
  theme_minimal() +
  labs( title="Distribusi Ekspresi Gen per Sampel (Data Mentah)", 
        x = "Sampel",
        y = "Log2(count+1)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ===================================================
#   Filter Gen dengan Reads Rendah (Optional)
# ===================================================
# 1. Filter: gen dengan total reads >= 10 di semua sampel
keep <- rowSums(counts) >= 10
counts_filtered <- counts[keep,]

# 2. Bandingkan jumlah gen sebelum dan sesudah filter
dim(counts)
dim(counts_filtered)

# ======================================
#   Langkah 5: Korelasi Antar Sampel
# ======================================
# 1. Load pheatmap
library(pheatmap)

# 2. Hitung korelasi antar sampel (pakai data yang sudah difilter)
cor_matrix <- cor(counts_filtered)

# 3. Tampilkan sebagai heatmap
pheatmap(cor_matrix, 
         main = "Korelasi Antar Sampel (Data Mentah)",
         display_numbers = TRUE,
         number_format = "%.2f",
         color = colorRampPalette(c("blue", "white", "red"))(50))

# ======================================
#   Langkah 6: PCA 
# ======================================
# 1. PCA dengan data yang sudah di-log (karena data counts sangat skewed)
log_data <- log2(counts_filtered + 1)
pca_result <- prcomp(t(log_data), scale=TRUE)
pca_result
pca_result[2]

# 2. Buat data frame untuk plotting
pca_df <- data.frame(
  PC1 = pca_result$x[,1],
  PC2 = pca_result$x[,2],
  condition = metadata$condition
  
)
head(pca_df)

# 3. Hitung persentase varians
var_explained <- summary(pca_result)$importance[2,1:2]*100

# 4. Plot PCA
ggplot(pca_df, aes(x=PC1, y=PC2, colour = condition)) +
  geom_point(size=5) +
  theme_minimal() +
  labs( title = "PCA Plot(Data Mentah - Log2)", 
        x = paste0("PC1:", round(var_explained[1], 1), "% variance"),
        y = paste0("PC2:", round(var_explained[2], 1), "% variance"))+
  theme(legend.position = "bottom")      


# Mission 02 : Normalization
# --------------------------------------------------------
# =============================
#   Langkah 1: Normalisasi
# =============================
# 1. Load library
library(DESeq2)

# 2. Buat DESeq2 object (dari count matrix dan metadata)
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~condition
)

# 3. Mulai normalisasi
dds <- DESeq(dds)

# 4. Ambil size factors (faktor normalisasi)
sizeFactors(dds)

# 5. Ambil data yang sudah dinormalisasi
# Data ini bisa kita gunakan untuk visualisasi (heatmap, PCA)
rld <- rlog(dds, blind = FALSE) # -> Jika dataset kecil, karena lebih lambat
vsd <- vst(dds, blind = FALSE)  # -> Jika dataset besar, karena lebih cepat

# =====================================
#   Langkah 2: Persiapan Visualisasi
# =====================================
# 1. Ubah data normalisasi (rlog) ke format panjang
rld_data <- assay(rld)
rld_long <- rld_data %>%
  as.data.frame() %>%
  mutate(gene=rownames(.)) %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "expression") %>%
  mutate(type="Data Normalisasi (rlog)")

# 2. Gabungkan kedua data
combined_data <- bind_rows(
  counts_long %>% mutate(value = log2(count +1),
                         type = "Data Mentah"),
  rld_long %>% mutate(value = expression,
                      type = "Data Normalisasi (rlog)")
)

# 3. Visualisasi dengan Boxplot
library(ggplot2)
ggplot(combined_data, aes(x=sample, y=value, fill = type)) +
  geom_boxplot() +
  facet_wrap(~type, scales="free_y", ncol=2 ) + 
  theme_minimal() +
  labs( title = "Perbandingan Distribusi Ekspresi Gen",
        subtitle = "Data Mentah vs Data Normalisasi (rlog)",
        x = "Sampel",
        y = "Log2 Ekspresi") +
  theme(axis.text.x = element_text(angle=45, hjust=1, size=8),
        strip.text = element_text(face = "bold", size=10),
        legend.position = "none")


# Mission 03 : Differential Gene Expression (DEG) Analysis
# --------------------------------------------------------
# ======================
# Langkah 1: DESeq2
# ======================
# 1. Membuat object dds
dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = metadata,
  design = ~condition
)

# 2. Jalankan DESeq2
dds <- DESeq(dds)

# 3. Ambil hasil perbandingan
res <- results(dds, contrast = c("condition", "DTG", "DMSO"))
res <- res[order(res$padj),]
head(res)

# ======================
# Langkah 2: Filter Gen
# ======================
# 1. Ubah data jadi format dataframe
res_df <- as.data.frame(res)

# 2. Tambahkan kolom status signifikansi
res_df$significant <- ifelse(
  !is.na(res_df$padj) & 
    res_df$padj < 0.05 &
    abs(res_df$log2FoldChange) > 1,
  "Signifikan", 
  "Tidak Signfikan"
)

# 3. Tambahkan kolom arah perubahan
res_df$direction <- ifelse(
  !is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange >1, "Naik", 
  ifelse(
    !is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange < -1, "Turun",
    "Tidak Signifikan"))


# 4. Filter gen signifikan
sig_genes <- res_df[res_df$significant == "Signifikan",]
nrow(sig_genes)

# 5. Lihat gen naik dan turun
table(sig_genes$direction)

# =========================
# Langkah 3: Export Hasil
# =========================
# 1. Export hasil
write.csv(res_df, "DEG_results_GSM285595.csv")

# 2. Lihat top 10 gen naik
top_up <- head(res_df[res_df$direction == "Naik",c("log2FoldChange", "pvalue", "padj")], 10)
print(top_up)

  
# 3. Lihat top 10 gen turun
top_down <- head(res_df[res_df$direction == "Turun", 
                        c("log2FoldChange", "pvalue", "padj")], 10)
print(top_down)

# Mission 04 : Visualization
# --------------------------------------------------------
# =========================
# Langkah 1: Volcano Plot
# =========================
library(dplyr)
library(ggplot2)

# 1. Ambil 10 gen dengan padj terkecil (paling signifikan)
top_genes <- res_df %>% 
  filter(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1) %>% 
  arrange(padj) %>%
  head(10)

# 2. Tambahkan kolom label (hanya untuk gen top)
res_df$label <- ifelse(rownames(res_df) %in% 
                         rownames(top_genes), rownames(res_df), "")

# 3. Buat volcano plot dengan label
#volcano_plot_labeled <- ggplot(res_df, aes(x=log2FoldChange, y = -log10(padj), color=direction))+
  
