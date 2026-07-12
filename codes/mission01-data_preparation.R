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
