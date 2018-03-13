---
output: html_document
---



## Unique Molecular Identifiers (UMIs) {#umichapter}

Thanks to Andreas Buness from [EMBL Monterotondo](https://www.embl.it/services/bioinformatics/) for collaboration on this section.

### Introduction

Unique Molecular Identifiers are short (4-10bp) random barcodes added to transcripts during reverse-transcription. They enable sequencing reads to be assigned to individual transcript molecules and thus the removal of amplification noise and biases from scRNASeq data. 

<div class="figure" style="text-align: center">
<img src="figures/UMI-Seq-protocol.png" alt="UMI sequencing protocol" width="90%" />
<p class="caption">(\#fig:intro-umi-protocol)UMI sequencing protocol</p>
</div>

When sequencing UMI containing data, techniques are used to specifically sequence only the end of the transcript containing the UMI (usually the 3' end).

### Mapping Barcodes

Since the number of unique barcodes ($4^N$, where $N$ is the length of UMI) is much smaller than the total number of molecules per cell (~$10^6$), each barcode will typically be assigned to multiple transcripts. Hence, to identify unique molecules both barcode and mapping location (transcript) must be used. The first step is to map UMI reads, for which we recommend using STAR since it is fast and outputs good quality BAM-alignments. Moreover, mapping locations can be useful for eg. identifying poorly-annotated 3' UTRs of transcripts.

UMI-sequencing typically consists of paired-end reads where one read from each pair captures the cell and UMI barcodes while the other read consists of exonic sequence from the transcript (Figure \@ref(fig:intro-umi-reads)). Note that trimming and/or filtering to remove reads containing poly-A sequence is recommended to avoid erors due to these read mapping to genes/transcripts with internal poly-A/poly-T sequences.

After processing the reads from a UMI experiment, the following conventions are often used:

1. The UMI is added to the read name of the other paired read. 

2. Reads are sorted into separate files by cell barcode
	+ For extremely large, shallow datasets, the cell barcode may be added to the read name as well to reduce the number of files.

<div class="figure" style="text-align: center">
<img src="figures/UMI-Seq-reads.png" alt="UMI sequencing reads, red lightning bolts represent different fragmentation locations" width="90%" />
<p class="caption">(\#fig:intro-umi-reads)UMI sequencing reads, red lightning bolts represent different fragmentation locations</p>
</div>

### Counting Barcodes

In theory, every unique UMI-transcript pair should represent all reads originating from a single RNA molecule. However, in practice this is frequently not the case and the most common reasons are:

1. __Different UMI does not necessarily mean different molecule__
	+ Due to PCR or sequencing errors, base-pair substitution events can result in new UMI sequences. Longer UMIs give more opportunity for errors to arise and based on estimates from cell barcodes we expect 7-10% of 10bp UMIs to contain at least one error. If not corrected for, this type of error will result in an overestimate of the number of transcripts.

2. __Different transcript does not necessarily mean different molecule__
	+ Mapping errors and/or multimapping reads may result in some UMIs being assigned to the wrong gene/transcript. This type of error will also result in an overestimate of the number of transcripts.

3. __Same UMI does not necessarily mean same molecule__
	+ Biases in UMI frequency and short UMIs can result in the same UMI being attached to different mRNA molecules from the same gene. Thus, the number of transcripts may be underestimated.

<div class="figure" style="text-align: center">
<img src="figures/UMI-Seq-errors.png" alt="Potential Errors in UMIs" width="90%" />
<p class="caption">(\#fig:intro-umi-errors)Potential Errors in UMIs</p>
</div>

### Correcting for Errors

How to best account for errors in UMIs remains an active area of research. The best approaches that we are aware of for resolving the issues mentioned above are:

1. [UMI-tools'](https://github.com/CGATOxford/UMI-tools) directional-adjacency method implements a procedure which considers both the number of mismatches and the relative frequency of similar UMIs to identify likely PCR/sequencing errors.

2. Currently an open question. The problem may be mitigated by removing UMIs with few reads to support their association with a particular transcript, or by removing all multi-mapping reads.

3. Simple saturation (aka "collision probability") correction proposed by [Grun, Kester and van Oudenaarden (2014)](http://www.nature.com/nmeth/journal/v11/n6/full/nmeth.2930.html#methods) to estimate the true number of molecules $M$:

$$M \approx -N*log(1 - \frac{n}{N})$$ 
where N = total number of unique UMI barcodes and n = number of observed barcodes.
 
An important caveat of this method is that it assumes that all UMIs are equally frequent. In most cases this is incorrect, since there is often a bias related to the GC content. 

<div class="figure" style="text-align: center">
<img src="figures/UMI-Seq-amp.png" alt="Per gene amplification rate" width="60%" />
<p class="caption">(\#fig:intro-umi-amp)Per gene amplification rate</p>
</div>

Determining how to best process and use UMIs is currently an active area of research in the bioinformatics community. We are aware of several methods that have recently been developed, including:

* [UMI-tools](https://github.com/CGATOxford/UMI-tools)
* [PoissonUMIs](https://github.com/tallulandrews/PoissonUMIs)
* [zUMIs](https://github.com/sdparekh/zUMIs)
* [dropEst](https://github.com/hms-dbmi/dropEst)

### Downstream Analysis

Current UMI platforms (DropSeq, InDrop, ICell8) exhibit low and highly variable capture efficiency as shown in the figure below. 

<div class="figure" style="text-align: center">
<img src="figures/UMI-Seq-capture.png" alt="Variability in Capture Efficiency" width="70%" />
<p class="caption">(\#fig:intro-umi-capture)Variability in Capture Efficiency</p>
</div>

This variability can introduce strong biases and it needs to be considered in downstream analysis. Recent analyses often pool cells/genes together based on cell-type or biological pathway to increase the power. Robust statistical analyses of this data is still an open research question and it remains to be determined how to best adjust for biases.

__Exercise 1__ We have provided you with UMI counts and read counts from induced pluripotent stem cells generated from three different individuals [@Tung2017-ba] (see: Chapter \@ref(exprs-qc) for details of this dataset).


```r
umi_counts <- read.table("tung/molecules.txt", sep = "\t")
read_counts <- read.table("tung/reads.txt", sep = "\t")
```
Using this data:

1. Plot the variability in capture efficiency

2. Determine the amplification rate: average number of reads per UMI.


