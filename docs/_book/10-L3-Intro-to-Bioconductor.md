---
output: html_document
---

## Data Types

### What is Tidy Data?

Tidy data is a concept largely defined by Hadley Wickham [@wickham_2014]. Tidy data has the following three characteristics:

1. Each variable has its own column.
2. Each observation has its own row.
3. Each value has its own cell.

Here is an example of some tidy data:

```
##   Students   Subject Years Score
## 1     Mark     Maths     1     5
## 2     Jane   Biology     2     6
## 3 Mohammed   Physics     3     4
## 4      Tom     Maths     2     7
## 5    Celia Computing     3     9
```

Here is an example of some untidy data:

```
##    Students    Sport Category Counts
## 1      Matt   Tennis     Wins      0
## 2      Matt   Tennis   Losses      1
## 3     Ellie    Rugby     Wins      3
## 4     Ellie    Rugby   Losses      2
## 5       Tim Football     Wins      1
## 6       Tim Football   Losses      4
## 7    Louise Swimming     Wins      2
## 8    Louise Swimming   Losses      2
## 9     Kelly  Running     Wins      5
## 10    Kelly  Running   Losses      1
```

Task 1: In what ways is the untidy data not tidy? How could we make the untidy data tidy?

Tidy data is generally easier to work with than untidy data, especially if you are working with packages such as ggplot. Fortunately, packages are available to make untidy data tidy. Today we will explore a few of the functions available in the tidyr package which can be used to make untidy data tidy. If you are interested in finding out more about tidying data, we recommend reading "R for Data Science", by Garrett Grolemund and Hadley Wickham. An electronic copy is available here: http://r4ds.had.co.nz/

The untidy data above is untidy because two variables (`Wins` and `Losses`) are stored in one column (`Category`). This is a common way in which data can be untidy. To tidy this data, we need to make `Wins` and `Losses` into columns, and store the values in `Counts` in these columns. Fortunately, there is a function from the tidyverse packages to perform this operation. The function is called `spread`, and it takes two arguments, `key` and `value`. You should pass the name of the column which contains multiple variables to `key`, and pass the name of the column which contains values from multiple variables to `value`. For example:


```r
library(tidyverse)
sports<-data.frame(Students=c("Matt", "Matt", "Ellie", "Ellie", "Tim", "Tim", "Louise", "Louise", "Kelly", "Kelly"), Sport=c("Tennis","Tennis", "Rugby", "Rugby","Football", "Football","Swimming","Swimming", "Running", "Running"), Category=c("Wins", "Losses", "Wins", "Losses", "Wins", "Losses", "Wins", "Losses", "Wins", "Losses"), Counts=c(0,1,3,2,1,4,2,2,5,1))
sports
```

```
##    Students    Sport Category Counts
## 1      Matt   Tennis     Wins      0
## 2      Matt   Tennis   Losses      1
## 3     Ellie    Rugby     Wins      3
## 4     Ellie    Rugby   Losses      2
## 5       Tim Football     Wins      1
## 6       Tim Football   Losses      4
## 7    Louise Swimming     Wins      2
## 8    Louise Swimming   Losses      2
## 9     Kelly  Running     Wins      5
## 10    Kelly  Running   Losses      1
```

```r
spread(sports, key=Category, value=Counts)
```

```
##   Students    Sport Losses Wins
## 1    Ellie    Rugby      2    3
## 2    Kelly  Running      1    5
## 3   Louise Swimming      2    2
## 4     Matt   Tennis      1    0
## 5      Tim Football      4    1
```

Task 2: The dataframe `foods` defined below is untidy. Work out why and use `spread()` to tidy it


```r
foods<-data.frame(student=c("Antoinette","Antoinette","Taylor", "Taylor", "Alexa", "Alexa"), Category=c("Dinner", "Dessert", "Dinner", "Dessert", "Dinner","Dessert"), Frequency=c(3,1,4,5,2,1))
```

The other common way in which data can be untidy is if the columns are values instead of variables. For example, the dataframe below shows the percentages some students got in tests they did in May and June. The data is untidy because the columns `May` and `June` are values, not variables.


```r
percentages<-data.frame(student=c("Alejandro", "Pietro", "Jane"), "May"=c(90,12,45), "June"=c(80,30,100))
```

Fortunately, there is a function in the tidyverse packages to deal with this problem too. `gather()` takes the names of the columns which are values, the `key` and the `value` as arguments. This time, the `key` is the name of the variable with values as column names, and the `value` is the name of the variable with values spread over multiple columns. Ie:


```r
gather(percentages, "May", "June", key="Month", value = "Percentage")
```

```
##     student Month Percentage
## 1 Alejandro   May         90
## 2    Pietro   May         12
## 3      Jane   May         45
## 4 Alejandro  June         80
## 5    Pietro  June         30
## 6      Jane  June        100
```

These examples don't have much to do with single-cell RNA-seq analysis, but are designed to help illustrate the features of tidy and untidy data. You will find it much easier to analyse your single-cell RNA-seq data if your data is stored in a tidy format. Fortunately, the data structures we commonly use to facilitate single-cell RNA-seq analysis usually encourage store your data in a tidy manner.

### What is Rich Data?

If you google 'rich data', you will find lots of different definitions for this term. In this course, we will use 'rich data' to mean data which is generated by combining information from multiple sources. For example, you could make rich data by creating an object in R which contains a matrix of gene expression values across the cells in your single-cell RNA-seq experiment, but also information about how the experiment was performed. Objects of the `SingleCellExperiment` class, which we will discuss below, are an example of rich data.



### What is Bioconductor?

From [Wikipedia](https://en.wikipedia.org/wiki/Bioconductor):
[Bioconductor](https://www.bioconductor.org/) is a free, open source and open development software project for the analysis and comprehension of genomic data generated by wet lab experiments in molecular biology. Bioconductor is based primarily on the statistical R programming language, but does contain contributions in other programming languages. It has two releases each year that follow the semiannual releases of R. At any one time there is a release version, which corresponds to the released version of R, and a development version, which corresponds to the development version of R. Most users will find the release version appropriate for their needs.

We strongly recommend all new comers and even experienced high-throughput data analysts to use well developed and maintained [Bioconductor methods and classes](https://www.bioconductor.org/developers/how-to/commonMethodsAndClasses/).

### `SingleCellExperiment` class

[`SingleCellExperiment`](http://bioconductor.org/packages/SingleCellExperiment) (SCE) is a S4 class for storing data from single-cell experiments. This includes specialized methods to store and retrieve spike-in information, dimensionality reduction coordinates and size factors for each cell, along with the usual metadata for genes and libraries.

In practice, an object of this class can be created using its constructor:

```r
library(SingleCellExperiment)
counts <- matrix(rpois(100, lambda = 10), ncol=10, nrow=10)
rownames(counts) <- paste("gene", 1:10, sep = "")
colnames(counts) <- paste("cell", 1:10, sep = "")
sce <- SingleCellExperiment(
    assays = list(counts = counts),
    rowData = data.frame(gene_names = paste("gene_name", 1:10, sep = "")),
    colData = data.frame(cell_names = paste("cell_name", 1:10, sep = ""))
)
sce
```

```
## class: SingleCellExperiment 
## dim: 10 10 
## metadata(0):
## assays(1): counts
## rownames(10): gene1 gene2 ... gene9 gene10
## rowData names(1): gene_names
## colnames(10): cell1 cell2 ... cell9 cell10
## colData names(1): cell_names
## reducedDimNames(0):
## spikeNames(0):
```

In the `SingleCellExperiment`, users can assign arbitrary names to entries of assays. To assist interoperability between packages, some suggestions for what the names should be for particular types of data are provided by the authors:

* __counts__: Raw count data, e.g., number of reads or transcripts for a particular gene.
* __normcounts__: Normalized values on the same scale as the original counts. For example, counts divided by cell-specific size factors that are centred at unity.
* __logcounts__: Log-transformed counts or count-like values. In most cases, this will be defined as log-transformed normcounts, e.g., using log base 2 and a pseudo-count of 1.
* __cpm__: Counts-per-million. This is the read count for each gene in each cell, divided by the library size of each cell in millions.
* __tpm__: Transcripts-per-million. This is the number of transcripts for each gene in each cell, divided by the total number of transcripts in that cell (in millions).

Each of these suggested names has an appropriate getter/setter method for convenient manipulation of the `SingleCellExperiment`. For example, we can take the (very specifically named) `counts` slot, normalise it and assign it to `normcounts` instead:


```r
normcounts(sce) <- log2(counts(sce) + 1)
sce
```

```
## class: SingleCellExperiment 
## dim: 10 10 
## metadata(0):
## assays(2): counts normcounts
## rownames(10): gene1 gene2 ... gene9 gene10
## rowData names(1): gene_names
## colnames(10): cell1 cell2 ... cell9 cell10
## colData names(1): cell_names
## reducedDimNames(0):
## spikeNames(0):
```

```r
dim(normcounts(sce))
```

```
## [1] 10 10
```

```r
head(normcounts(sce))
```

```
##          cell1    cell2    cell3    cell4    cell5    cell6    cell7
## gene1 3.169925 3.169925 2.000000 2.584963 2.584963 3.321928 3.584963
## gene2 3.459432 1.584963 3.584963 3.807355 3.700440 3.700440 3.000000
## gene3 3.000000 3.169925 3.807355 3.169925 3.321928 3.321928 3.321928
## gene4 3.584963 3.459432 3.000000 3.807355 3.700440 3.700440 3.700440
## gene5 3.906891 3.000000 3.169925 3.321928 3.584963 3.459432 3.807355
## gene6 3.700440 3.700440 3.584963 4.000000 3.169925 3.000000 3.459432
##          cell8    cell9   cell10
## gene1 3.321928 3.807355 2.807355
## gene2 3.807355 3.700440 4.000000
## gene3 2.584963 4.000000 3.700440
## gene4 3.169925 3.584963 3.700440
## gene5 3.807355 2.584963 3.584963
## gene6 3.321928 3.459432 4.000000
```

### `scater` package

[`scater`](http://bioconductor.org/packages/scater/) is a R package for single-cell RNA-seq analysis [@McCarthy2017-kb]. The package contains several useful methods for quality control, visualisation and pre-processing of data prior to further downstream analysis.

`scater` features the following functionality:

* Automated computation of QC metrics
* Transcript quantification from read data with pseudo-alignment
* Data format standardisation
* Rich visualizations for exploratory analysis
* Seamless integration into the Bioconductor universe
* Simple normalisation methods

We highly recommend to use `scater` for all single-cell RNA-seq analyses and `scater` is the basis of the first part of the course.

As illustrated in the figure below, `scater` will help you with quality control, filtering and normalization of your expression matrix following mapping and alignment. <span style="color:red">Keep in mind that this figure represents the original version of `scater` where an `SCESet` class was used. In the newest version this figure is still correct, except that `SCESet` can be substituted with the `SingleCellExperiment` class.</span>


![](figures/scater_qc_workflow.png)


