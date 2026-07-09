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
counts[,-1]
