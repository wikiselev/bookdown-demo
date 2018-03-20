---
output: html_document
---

# Tabula Muris

## Introduction

To give you hands-on experience analyzing from start to finish a single-cell RNASeq dataset we will be using as an example, data from the [Tabula Muris](https://www.biorxiv.org/content/early/2017/12/20/237446) initial release. The Tabula Muris 
is an international collaboration with the aim to profile every cell-type in the mouse using a standardized method. They combine highthroughput but low-coverage 10X data with lower throughput
but high-coverage FACS-sorted cells + Smartseq2. 

The initial release of the data (20 Dec 2017), contain almost 100,000 cells across 20 different tissues/organs. You have been assigned one of these tissues as an example to work on over this course, and on Friday each person will have 3 minutes to present the result for their tissue. 

## Downloading the data
Unlike most single-cell RNASeq data Tabula Muris has release their data through the [figshare](https://figshare.com/) platform rather than uploading it to [GEO](https://www.ncbi.nlm.nih.gov/geo/) or [ArrayExpress](https://www.ebi.ac.uk/arrayexpress/). You can find the data by using the doi's in their paper : [10.6084/m9.figshare.5715040](https://figshare.com/articles/Single-cell_RNA-seq_data_from_Smart-seq2_sequencing_of_FACS_sorted_cells/5715040) for FACS/Smartseq2 and [10.6084/m9.figshare.5715025](https://figshare.com/articles/Single-cell_RNA-seq_data_from_microfluidic_emulsion/5715025) for 10X data. The data can be downloaded manually by clinking the doi links or by using the command-line commands below:

Terminal-based download of FACS data: 


```bash
wget https://ndownloader.figshare.com/files/10038307
unzip 10038307
wget https://ndownloader.figshare.com/files/10038310
mv 10038310 FACS_metadata.csv
wget https://ndownloader.figshare.com/files/10039267
mv 10039267 FACS_annotations.csv
```

Terminal-based download of 10X data:

```bash
wget https://ndownloader.figshare.com/files/10038325
unzip 10038325
wget https://ndownloader.figshare.com/files/10038328
mv 10038328 droplet_metadata.csv
wget https://ndownloader.figshare.com/files/10039264
mv 10039264 droplet_annotation.csv
```

Note if you download the data by hand you should unzip & rename the files as above before continuing.

You should now have two folders : "FACS" and "droplet" and one annotation and metadata file for each. To inspect these files you can use the `head` to see the top few lines of the text files (Press "q" to exit):

```bash
head -n 10 droplet_metadata.csv
```

```
## channel,mouse.id,tissue,subtissue,mouse.sex
## 10X_P4_0,3-M-8,Tongue,NA,M
## 10X_P4_1,3-M-9,Tongue,NA,M
## 10X_P4_2,3-M-8/9,Liver,hepatocytes,M
## 10X_P4_3,3-M-8,Bladder,NA,M
## 10X_P4_4,3-M-9,Bladder,NA,M
## 10X_P4_5,3-M-8,Kidney,NA,M
## 10X_P4_6,3-M-9,Kidney,NA,M
## 10X_P4_7,3-M-8,Spleen,NA,M
## 10X_P7_0,3-F-56,Liver,NA,F
```
You can also check the number of rows in each file using:

```bash
wc -l droplet_annotation.csv
```

```
## 54838 droplet_annotation.csv
```

__Exercise__
How many cells do we have annotations for from FACS? from 10X?

__Answer__
FACS : 54,838 cells
Droplet : 42,193 cells

## Reading the data (Smartseq2)

We can now read in the relevant count matrix from the comma-separated file. Then inspect the resulting dataframe:


```r
dat = read.delim("FACS/Kidney-counts.csv", sep=",", header=TRUE)
dat[1:5,1:5]
```

```
##               X A14.MAA000545.3_8_M.1.1 E1.MAA000545.3_8_M.1.1
## 1 0610005C13Rik                       0                      0
## 2 0610007C21Rik                       1                      0
## 3 0610007L01Rik                       0                      0
## 4 0610007N19Rik                       0                      0
## 5 0610007P08Rik                       0                      0
##   M4.MAA000545.3_8_M.1.1 O21.MAA000545.3_8_M.1.1
## 1                      0                       0
## 2                      0                       0
## 3                      0                       0
## 4                      0                       0
## 5                      0                       0
```
We can see that the first column in the dataframe is the gene names, so first we move these to the rownames so we have a numeric matrix:


```r
dim(dat)
```

```
## [1] 23433   866
```

```r
rownames(dat) <- dat[,1]
dat <- dat[,-1]
```

Since this is a Smartseq2 dataset it may contain spike-ins so lets check:


```r
rownames(dat)[grep("^ERCC-", rownames(dat))]
```

```
##  [1] "ERCC-00002" "ERCC-00003" "ERCC-00004" "ERCC-00009" "ERCC-00012"
##  [6] "ERCC-00013" "ERCC-00014" "ERCC-00016" "ERCC-00017" "ERCC-00019"
## [11] "ERCC-00022" "ERCC-00024" "ERCC-00025" "ERCC-00028" "ERCC-00031"
## [16] "ERCC-00033" "ERCC-00034" "ERCC-00035" "ERCC-00039" "ERCC-00040"
## [21] "ERCC-00041" "ERCC-00042" "ERCC-00043" "ERCC-00044" "ERCC-00046"
## [26] "ERCC-00048" "ERCC-00051" "ERCC-00053" "ERCC-00054" "ERCC-00057"
## [31] "ERCC-00058" "ERCC-00059" "ERCC-00060" "ERCC-00061" "ERCC-00062"
## [36] "ERCC-00067" "ERCC-00069" "ERCC-00071" "ERCC-00073" "ERCC-00074"
## [41] "ERCC-00075" "ERCC-00076" "ERCC-00077" "ERCC-00078" "ERCC-00079"
## [46] "ERCC-00081" "ERCC-00083" "ERCC-00084" "ERCC-00085" "ERCC-00086"
## [51] "ERCC-00092" "ERCC-00095" "ERCC-00096" "ERCC-00097" "ERCC-00098"
## [56] "ERCC-00099" "ERCC-00104" "ERCC-00108" "ERCC-00109" "ERCC-00111"
## [61] "ERCC-00112" "ERCC-00113" "ERCC-00116" "ERCC-00117" "ERCC-00120"
## [66] "ERCC-00123" "ERCC-00126" "ERCC-00130" "ERCC-00131" "ERCC-00134"
## [71] "ERCC-00136" "ERCC-00137" "ERCC-00138" "ERCC-00142" "ERCC-00143"
## [76] "ERCC-00144" "ERCC-00145" "ERCC-00147" "ERCC-00148" "ERCC-00150"
## [81] "ERCC-00154" "ERCC-00156" "ERCC-00157" "ERCC-00158" "ERCC-00160"
## [86] "ERCC-00162" "ERCC-00163" "ERCC-00164" "ERCC-00165" "ERCC-00168"
## [91] "ERCC-00170" "ERCC-00171"
```

Now we can extract much of the metadata for this data from the column names:


```r
cellIDs <- colnames(dat)
cell_info <- strsplit(cellIDs, "\\.")
Well <- lapply(cell_info, function(x){x[1]})
Well <- unlist(Well)
Plate <- unlist(lapply(cell_info, function(x){x[2]}))
Mouse <- unlist(lapply(cell_info, function(x){x[3]}))
```
We can check the distributions of each of these metadata classifications:


```r
summary(factor(Mouse))
```

```
## 3_10_M 3_11_M 3_38_F 3_39_F  3_8_M  3_9_M 
##    104    196    237    169     82     77
```

We can also check if any technical factors are confounded:


```r
table(Mouse, Plate)
```

```
##         Plate
## Mouse    B001717 B002775 MAA000545 MAA000752 MAA000801 MAA000922
##   3_10_M       0       0         0       104         0         0
##   3_11_M       0       0         0         0       196         0
##   3_38_F     237       0         0         0         0         0
##   3_39_F       0     169         0         0         0         0
##   3_8_M        0       0        82         0         0         0
##   3_9_M        0       0         0         0         0        77
```

Lastly we will read the computationally inferred cell-type annotation and match them to the cell in our expression matrix:


```r
ann <- read.table("FACS_annotations.csv", sep=",", header=TRUE)
ann <- ann[match(cellIDs, ann[,1]),]
celltype <- ann[,3]
```

## Building a scater object
To create a SingleCellExperiment object we must put together all the cell annotations into a single dataframe, since the experimental batch (pcr plate) is completely confounded with donor mouse we will only keep one of them.


```r
require("SingleCellExperiment")
```

```
## Loading required package: SingleCellExperiment
```

```
## Loading required package: SummarizedExperiment
```

```
## Loading required package: methods
```

```
## Loading required package: GenomicRanges
```

```
## Loading required package: stats4
```

```
## Loading required package: BiocGenerics
```

```
## Loading required package: parallel
```

```
## 
## Attaching package: 'BiocGenerics'
```

```
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
```

```
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, sd, var, xtabs
```

```
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, cbind, colMeans,
##     colnames, colSums, do.call, duplicated, eval, evalq, Filter,
##     Find, get, grep, grepl, intersect, is.unsorted, lapply,
##     lengths, Map, mapply, match, mget, order, paste, pmax,
##     pmax.int, pmin, pmin.int, Position, rank, rbind, Reduce,
##     rowMeans, rownames, rowSums, sapply, setdiff, sort, table,
##     tapply, union, unique, unsplit, which, which.max, which.min
```

```
## Loading required package: S4Vectors
```

```
## 
## Attaching package: 'S4Vectors'
```

```
## The following object is masked from 'package:base':
## 
##     expand.grid
```

```
## Loading required package: IRanges
```

```
## Loading required package: GenomeInfoDb
```

```
## Loading required package: Biobase
```

```
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
```

```
## Loading required package: DelayedArray
```

```
## Loading required package: matrixStats
```

```
## 
## Attaching package: 'matrixStats'
```

```
## The following objects are masked from 'package:Biobase':
## 
##     anyMissing, rowMedians
```

```
## 
## Attaching package: 'DelayedArray'
```

```
## The following objects are masked from 'package:matrixStats':
## 
##     colMaxs, colMins, colRanges, rowMaxs, rowMins, rowRanges
```

```
## The following object is masked from 'package:base':
## 
##     apply
```

```r
require("scater")
```

```
## Loading required package: scater
```

```
## Loading required package: ggplot2
```

```
## 
## Attaching package: 'scater'
```

```
## The following object is masked from 'package:S4Vectors':
## 
##     rename
```

```
## The following object is masked from 'package:stats':
## 
##     filter
```

```r
cell_anns <- data.frame(mouse = Mouse, well=Well, type=celltype)
rownames(cell_anns) <- colnames(dat)
sceset <- SingleCellExperiment(assays = list(counts = as.matrix(dat)), colData=cell_anns)
```

Finally if the dataset contains spike-ins we a hidden variable in the SingleCellExperiment object to track them:

```r
isSpike(sceset, "ERCC") <- grepl("ERCC-", rownames(sceset))
```


## Reading the data (10X)
Due to the large size and sparsity of 10X data (upto 90% of the expression matrix may be 0s) it is typically 
stored as a sparse matrix. The default output format for CellRanger is an .mtx file which stores this sparse 
matrix as a column of row coordinates, a column of column corodinates, and a column of expression values > 0. 
Note if you look at the .mtx file you will see two header lines followed by a line detailing the 
total number of rows, columns and counts for the full matrix. Since only the coordinates are stored in the .mtx 
file, the names of each row & column must be stored separately in the "genes.tsv" and "barcodes.tsv" files 
respectively.

We will be using the "Matrix" package to store matrices in sparse-matrix format in R.


```r
require("Matrix")
```

```
## Loading required package: Matrix
```

```
## 
## Attaching package: 'Matrix'
```

```
## The following object is masked from 'package:S4Vectors':
## 
##     expand
```

```r
cellbarcodes <- read.table("droplet/Kidney-10X_P4_5/barcodes.tsv")
genenames <- read.table("droplet/Kidney-10X_P4_5/genes.tsv")
molecules <- Matrix::readMM("droplet/Kidney-10X_P4_5/matrix.mtx")
```
Now we will add the appropriate row and column names. However, if you inspect the read cellbarcodes you will see that they are just the barcode sequence associated with each cell. This is a problem since each batch of 10X data uses the same pool of barcodes so if we need to combine data from multiple 10X batches the cellbarcodes will not be unique. Hence we will attach the batch ID to each cell barcode:

```r
head(cellbarcodes)
```

```
##                   V1
## 1 AAACCTGAGATGCCAG-1
## 2 AAACCTGAGTGTCCAT-1
## 3 AAACCTGCAAGGCTCC-1
## 4 AAACCTGTCCTTGCCA-1
## 5 AAACGGGAGCTGAACG-1
## 6 AAACGGGCAGGACCCT-1
```


```r
rownames(molecules) <- genenames[,1]
colnames(molecules) <- paste("10X_P4_5", cellbarcodes[,1], sep="_")
```
Now lets get the metadata and computational annotations for this data:


```r
meta <- read.delim("droplet_metadata.csv", sep=",", header=TRUE)
head(meta)
```

```
##    channel mouse.id  tissue   subtissue mouse.sex
## 1 10X_P4_0    3-M-8  Tongue        <NA>         M
## 2 10X_P4_1    3-M-9  Tongue        <NA>         M
## 3 10X_P4_2  3-M-8/9   Liver hepatocytes         M
## 4 10X_P4_3    3-M-8 Bladder        <NA>         M
## 5 10X_P4_4    3-M-9 Bladder        <NA>         M
## 6 10X_P4_5    3-M-8  Kidney        <NA>         M
```
Here we can see that we need to use "10X_P4_5" to find the metadata for this batch, also note that the format of the mouse ID is different in this metadata table with hyphens instead of underscores and with the gender in the middle of the ID. From checking the methods section of the accompanying paper we know that the same 8 mice were used for both droplet and plate-based techniques. So we need to fix the mouse IDs to be consistent with those used in the FACS experiments. 


```r
meta[meta$channel == "10X_P4_5",]
```

```
##    channel mouse.id tissue subtissue mouse.sex
## 6 10X_P4_5    3-M-8 Kidney      <NA>         M
```

```r
mouseID <- "3_8_M"
```
Note: depending on the tissue you have been assigned you may have 10X data from mixed samples : e.g. mouse id = 3-M-5/6. You should still reformat these to be consistent but they will not match mouse ids from the FACS data which may affect your downstream analysis. If the mice weren't from an inbred strain it would be possible to assign individual cells to a specific mouse using exonic-SNPs but that is beyond the scope of this course.


```r
ann <- read.delim("droplet_annotation.csv", sep=",", header=TRUE)
head(ann)
```

```
##                        cell  tissue cell_ontology_class
## 1 10X_P4_3_AAAGTAGAGATGCCAG Bladder    mesenchymal cell
## 2 10X_P4_3_AACCGCGTCCAACCAA Bladder    mesenchymal cell
## 3 10X_P4_3_AACTCCCGTCGGGTCT Bladder    mesenchymal cell
## 4 10X_P4_3_AACTCTTAGTTGCAGG Bladder        bladder cell
## 5 10X_P4_3_AACTCTTTCATAACCG Bladder    mesenchymal cell
## 6 10X_P4_3_AAGACCTAGATCCGAG Bladder    endothelial cell
##                      cell_ontology_term_iri cell_ontology_id
## 1 http://purl.obolibrary.org/obo/CL_0008019       CL:0008019
## 2 http://purl.obolibrary.org/obo/CL_0008019       CL:0008019
## 3 http://purl.obolibrary.org/obo/CL_0008019       CL:0008019
## 4 http://purl.obolibrary.org/obo/CL_1001319       CL:1001319
## 5 http://purl.obolibrary.org/obo/CL_0008019       CL:0008019
## 6 http://purl.obolibrary.org/obo/CL_0000115       CL:0000115
```
Again you will find a slight formating difference between the cellID in the annotation and the cellbarcodes which we will have to correct before matching them.


```r
ann[,1] <- paste(ann[,1], "-1", sep="")
ann_subset <- ann[match(colnames(molecules), ann[,1]),]
celltype <- ann_subset[,3]
```

Now lets build the cell-metadata dataframe:

```r
cell_anns <- data.frame(mouse = rep(mouseID, times=ncol(molecules)), type=celltype)
rownames(cell_anns) <- colnames(molecules);
```

__Exercise__ Repeat the above for the other 10X batches for your tissue.

__Answer__



## Building a scater object

Now that we have read the 10X data in multiple batches we need to combine them into a single SingleCellExperiment object. First we will check that the gene names are the same and in the same order across all batches:


```r
identical(rownames(molecules1), rownames(molecules2))
```

```
## [1] TRUE
```

```r
identical(rownames(molecules1), rownames(molecules3))
```

```
## [1] TRUE
```

Now we'll check that there aren't any repeated cellIDs:

```r
sum(colnames(molecules1) %in% colnames(molecules2))
```

```
## [1] 0
```

```r
sum(colnames(molecules1) %in% colnames(molecules3))
```

```
## [1] 0
```

```r
sum(colnames(molecules2) %in% colnames(molecules3))
```

```
## [1] 0
```

Everything is ok, so we can go ahead and combine them:


```r
all_molecules <- cbind(molecules1, molecules2, molecules3)
all_cell_anns <- as.data.frame(rbind(cell_anns1, cell_anns2, cell_anns3))
all_cell_anns$batch <- rep(c("10X_P4_5", "10X_P4_6","10X_P7_5"), times = c(nrow(cell_anns1), nrow(cell_anns2), nrow(cell_anns3)))
```

__Exercise__
How many cells are in the whole dataset?

__Answer__


Now build the SingleCellExperiment object. One of the advantages of the SingleCellExperiment class is that it is capable of storing data in normal matrix or sparse matrix format, as well as HDF5 format which allows large non-sparse matrices to be stored & accessed on disk in an efficient manner rather than loading the whole thing into RAM.


```r
require("SingleCellExperiment")
require("scater")
all_molecules <- as.matrix(all_molecules)
sceset <- SingleCellExperiment(assays = list(counts = as.matrix(all_molecules)), colData=all_cell_anns)
```

Since this is 10X data it will not contain spike-ins, so we just save the data:

```r
saveRDS(sceset, "kidney_droplet.rds")
```

## Advanced Exercise

Write an R function/script which will fully automate this procedure for each data-type for any tissue.
