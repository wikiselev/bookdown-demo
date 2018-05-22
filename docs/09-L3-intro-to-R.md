---
output: html_document
---



# Introduction to R/Bioconductor

## Installing packages

### CRAN

The Comprehensive R Archive Network [CRAN](https://cran.r-project.org/) is the biggest archive of R packages. There are few requirements for uploading packages besides building and installing succesfully, hence documentation and support is often minimal and figuring how to use these packages can be a challenge it itself. CRAN is the default repository R will search to find packages to install:


```r
install.packages("devtools")
require("devtools")
```

### Github

[Github](https://github.com/) isn't specific to R, any code of any type in any state can be uploaded. There is no guarantee a package uploaded to github will even install, nevermind do what it claims to do. R packages can be downloaded and installed directly from github using the "devtools" package installed above.


```r
devtools::install_github("tallulandrews/M3Drop")
```

Github is also a version control system which stores multiple versions of any package. By default the most recent "master" version of the package is installed. If you want an older version or the development branch this can be specified using the "ref" parameter:


```r
# different branch
devtools::install_github("tallulandrews/M3D", ref="nbumi")
# previous commit
devtools::install_github("tallulandrews/M3Drop", ref="434d2da28254acc8de4940c1dc3907ac72973135")
```
Note: make sure you re-install the M3Drop master branch for later in the course.

### Bioconductor
Bioconductor is a repository of R-packages specifically for biological analyses. It has the strictest requirements for submission, including installation on every platform and full documentation with a tutorial (called a vignette) explaining how the package should be used. Bioconductor also encourages utilization of standard data structures/classes and coding style/naming conventions, so that, in theory, packages and analyses can be combined into large pipelines or workflows. 



```r
source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")
```

Note: in some situations it is necessary to substitute "http://" for "https://" in the above depending on the security features of your internet connection/network. 

Bioconductor also requires creators to support their packages and has a regular 6-month release schedule. Make sure you are using the most recent release of bioconductor before trying to install packages for the course.


```r
source("https://bioconductor.org/biocLite.R")
biocLite("BiocUpgrade")
```

### Source
The final way to install packages is directly from source. In this case you have to download a fully built source code file, usually packagename.tar.gz, or clone the github repository and rebuild the package yourself. Generally this will only be done if you want to edit a package yourself, or if for some reason the former methods have failed.


```r
install.packages("M3Drop_3.05.00.tar.gz", type="source")
```

## Installation instructions:
All the packages necessary for this course are available [here](https://github.com/hemberg-lab/scRNA.seq.course/blob/master/Dockerfile). Starting from "RUN Rscript -e "install.packages('devtools')" ", run each of the commands (minus "RUN") on the command line or start an R session and run each of the commands within the quotation marks. Note the ordering of the installation is important in some cases, so make sure you run them in order from top to bottom. 

## Data-types/classes

R is a high level language so the underlying data-type is generally not important. The exception if you are accessing R data directly using another language such as C, but that is beyond the scope of this course. Instead we will consider the basic data classes: numeric, integer, logical, and character, and the higher level data class called "factor". You can check what class your data is using the "class()" function.

Aside: R can also store data as "complex" for complex numbers but generally this isn't relevant for biological analyses.

### Numeric

The "numeric" class is the default class for storing any numeric data - integers, decimal numbers, numbers in scientific notation, etc... 


```r
x = 1.141
class(x)
```

```
## [1] "numeric"
```

```r
y = 42
class(y)
```

```
## [1] "numeric"
```

```r
z = 6.02e23
class(z)
```

```
## [1] "numeric"
```

Here we see that even though R has an "integer" class and 42 could be stored more efficiently as an integer the default is to store it as "numeric". If we want 42 to be stored as an integer we must "coerce" it to that class:


```r
y = as.integer(42)
class(y)
```

```
## [1] "integer"
```


Coercion will force R to store data as a particular class, if our data is incompatible with that class it will still do it but the data will be converted to NAs:


```r
as.numeric("H")
```

```
## Warning: NAs introduced by coercion
```

```
## [1] NA
```

Above we tried to coerce "character" data, identified by the double quotation marks, into numeric data which doesn't make sense, so we triggered ("threw") an warning message. Since this is only a warning R would continue with any subsequent commands in a script/function, whereas an "error" would cause R to halt. 

### Character/String

The "character" class stores all kinds of text data. Programing convention calls data containing multiple letters a "string", thus most R functions which act on character data will refer to the data as "strings" and will often have "str" or "string" in it's name. Strings are identified by being flanked by double quotation marks, whereas variable/function names are not:

```r
x = 5

a = "x" # character "x"
a
```

```
## [1] "x"
```

```r
b = x # variable x
b
```

```
## [1] 5
```

In addition to standard alphanumeric characters, strings can also store various special characters. Special characters are identified using a backlash followed by a single character, the most relevant are the special character for tab : `\t` and new line : `\n`. To demonstrate the these special characters lets concatenate (cat) together two strings with these characters separating (sep) them:

```r
cat("Hello", "World", sep= " ")
```

```
## Hello World
```

```r
cat("Hello", "World", sep= "\t")
```

```
## Hello	World
```

```r
cat("Hello", "World", sep= "\n")
```

```
## Hello
## World
```
Note that special characters work differently in different functions. For instance the `paste` function does the same thing as `cat` but does not recognize special characters.


```r
paste("Hello", "World", sep= " ")
```

```
## [1] "Hello World"
```

```r
paste("Hello", "World", sep= "\t")
```

```
## [1] "Hello\tWorld"
```

```r
paste("Hello", "World", sep= "\n")
```

```
## [1] "Hello\nWorld"
```

Single or double backslash is also used as an `escape` character to turn off special characters or allow quotation marks to be included in strings:


```r
cat("This \"string\" contains quotation marks.")
```

```
## This "string" contains quotation marks.
```

Special characters are generally only used in pattern matching, and reading/writing data to files. For instance this is how you would read a tab-separated file into R.

```r
dat = read.delim("file.tsv", sep="\t")
```

Another special type of character data are colours. Colours can be specified in three main ways: by name from those [available](http://bxhorn.com/r-color-tables/), by red, green, blue values using the `rgb` function, and by hue (colour), saturation (colour vs white) and value (colour/white vs black) using the `hsv` function. By default rgb and hsv expect three values in 0-1 with an optional fourth value for transparency. Alternatively, sets of predetermined colours with useful properties can be loaded from many different packages with [RColorBrewer](http://colorbrewer2.org/) being one of the most popular.


```r
reds = c("red", rgb(1,0,0), hsv(0, 1, 1))
reds
```

```
## [1] "red"     "#FF0000" "#FF0000"
```

```r
barplot(c(1,1,1), col=reds, names=c("by_name", "by_rgb", "by_hsv"))
```

<img src="09-L3-intro-to-R_files/figure-html/unnamed-chunk-16-1.png" width="90%" style="display: block; margin: auto;" />

### Logical

The `logical` class stores boolean truth values, i.e. TRUE and FALSE. It is used for storing the results of logical operations and conditional statements will be coerced to this class. Most other data-types can be coerced to boolean without triggering (or "throwing") error messages, which may cause unexpected behaviour.


```r
x = TRUE
class(x)
```

```
## [1] "logical"
```

```r
y = "T"
as.logical(y)
```

```
## [1] TRUE
```

```r
z = 5
as.logical(z)
```

```
## [1] TRUE
```

```r
x = FALSE
class(x)
```

```
## [1] "logical"
```

```r
y = "F"
as.logical(y)
```

```
## [1] FALSE
```

```r
z = 0
as.logical(z)
```

```
## [1] FALSE
```
__Exercise 1__
Experiment with other character and numeric values, which are coerced to TRUE or FALSE? which are coerced to neither? 
Do you ever throw a warning/error message?

### Factors

String/Character data is very memory inefficient to store, each letter generally requires the same amount of memory as any integer. Thus when storing a vector of strings with repeated elements it is more efficient assign each element to an integer and store the vector as integers and an additional string-to-integer association table. Thus, by default R will read in text columns of a data table as factors. 


```r
str_vector = c("Apple", "Apple", "Banana", "Banana", "Banana", "Carrot", "Carrot", "Apple", "Banana") 
factored_vector = factor(str_vector)
factored_vector
```

```
## [1] Apple  Apple  Banana Banana Banana Carrot Carrot Apple  Banana
## Levels: Apple Banana Carrot
```

```r
as.numeric(factored_vector)
```

```
## [1] 1 1 2 2 2 3 3 1 2
```

The double nature of factors can cause some unintuitive behaviour. E.g. joining two factors together will convert them to the numeric form and the original strings will be lost.


```r
c(factored_vector, factored_vector)
```

```
##  [1] 1 1 2 2 2 3 3 1 2 1 1 2 2 2 3 3 1 2
```

Likewise if due to formatting issues numeric data is mistakenly interpretted as strings, then you must convert the factor back to strings before coercing to numeric values:

```r
x = c("20", "25", "23", "38", "20", "40", "25", "30")
x = factor(x)
as.numeric(x)
```

```
## [1] 1 3 2 5 1 6 3 4
```

```r
as.numeric(as.character(x))
```

```
## [1] 20 25 23 38 20 40 25 30
```

To make R read text as character data instead of factors set the environment option `stringsAsFactors=FALSE`. This must be done at the start of each R session.


```r
options(stringsAsFactors=FALSE)
```
__Exercise__
How would you use factors to create a vector of colours for an arbitrarily long vector of fruits like `str_vector` above?
__Answer__


### Checking class/type
We recommend checking your data is of the correct class after reading from files:


```r
x = 1.4
is.numeric(x)
```

```
## [1] TRUE
```

```r
is.character(x)
```

```
## [1] FALSE
```

```r
is.logical(x)
```

```
## [1] FALSE
```

```r
is.factor(x)
```

```
## [1] FALSE
```


## Basic data structures
So far we have only looked at single values and vectors. Vectors are the simplest data structure in R. They are a 1-dimensional array of data all of the same type. If the input when creating a vector is of different types it will be coerced to the data-type that is most consistent with the data.


```r
x = c("Hello", 5, TRUE)
x
```

```
## [1] "Hello" "5"     "TRUE"
```

```r
class(x)
```

```
## [1] "character"
```
Here we tried to put character, numeric and logical data into a single vector so all the values were coerced to `character` data.

A `matrix` is the two dimensional version of a vector, it also requires all data to be of the same type. 
If we combine a character vector and a numeric vector into a matrix, all the data will be coerced to characters:


```r
x = c("A", "B", "C")
y = c(1, 2, 3)
class(x)
```

```
## [1] "character"
```

```r
class(y)
```

```
## [1] "numeric"
```

```r
m = cbind(x, y)
m
```

```
##      x   y  
## [1,] "A" "1"
## [2,] "B" "2"
## [3,] "C" "3"
```

The quotation marks indicate that the numeric vector has been coerced to characters. Alternatively, to store data with columns of different data-types we can use a dataframe.


```r
z = data.frame(x, y)
z
```

```
##   x y
## 1 A 1
## 2 B 2
## 3 C 3
```

```r
class(z[,1])
```

```
## [1] "character"
```

```r
class(z[,2])
```

```
## [1] "numeric"
```

If you have set stringsAsFactors=FALSE as above you will find the first column remains characters, otherwise it will be automatically converted to a factor.


```r
options(stringsAsFactors=TRUE)
z = data.frame(x, y)
class(z[,1])
```

```
## [1] "factor"
```

Another difference between matrices and dataframes is the ability to select columns using the `$` operator:


```r
m$x # throws an error
z$x # ok
```

The final basic data structure is the `list`. Lists allow data of different types and different lengths to be stored in a single object. Each element of a list can be any other R object : data of any type, any data structure, even other lists or functions. 


```r
l = list(m, z)
ll = list(sublist=l, a_matrix=m, numeric_value=42, this_string="Hello World", even_a_function=cbind)
ll
```

```
## $sublist
## $sublist[[1]]
##      x   y  
## [1,] "A" "1"
## [2,] "B" "2"
## [3,] "C" "3"
## 
## $sublist[[2]]
##   x y
## 1 A 1
## 2 B 2
## 3 C 3
## 
## 
## $a_matrix
##      x   y  
## [1,] "A" "1"
## [2,] "B" "2"
## [3,] "C" "3"
## 
## $numeric_value
## [1] 42
## 
## $this_string
## [1] "Hello World"
## 
## $even_a_function
## function (..., deparse.level = 1) 
## .Internal(cbind(deparse.level, ...))
## <bytecode: 0x55fca094d378>
## <environment: namespace:base>
```

Lists are most commonly used when returning a large number of results from a function that do not fit into any of the previous data structures. 

## More information

You can get more information about any R commands relevant to these datatypes using by typing `?function` in an interactive session.
