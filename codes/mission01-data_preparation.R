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
