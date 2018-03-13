---
output: html_document
---

## Using STAR to Align Reads

Now we have trimmed our reads and established that they are of good quality, we would like to map them to a reference genome. This process is known as alignment. Some form of alignment is generally required if we want to quantify gene expression or find genes which are differentially expressed between samples.

Many tools have been developed for read alignment, but today we will focus on two. The first tool we will consider is STAR [@dobin]. For each read in our reads data, STAR tries to find the longest possible sequence which matches one or more sequences in the reference genome. For example, in the figure below, we have a read (blue) which spans two exons and an alternative splicing junction (purple). STAR finds that the first part of the read is the same as the sequence of the first exon, whilst the second part of the read matches the sequence in the second exon. Because STAR is able to recognise splicing events in this way, it is described as a 'splice aware' aligner.

![Figure 1: Diagram of how STAR performs alignments, taken from Dobin et al.](L1-images/STAR_explanation.png)

Usually STAR aligns reads to a reference genome, potentially allowing it to detect novel splicing events or chromosomal rearrangements. However, one issue with STAR is that it needs a lot of RAM, especially if your reference genome is large (eg. mouse and human). To speed up our analysis today, we will use STAR to align reads from to a reference transcriptome of 2000 transcripts. Note that this is NOT normal or recommended practice, we only do it here for reasons of time. We recommend that normally you should align to a reference genome.

Two steps are required to perform STAR alignment. In the first step, the user provides STAR with reference genome sequences (FASTA) and annotations (GTF), which STAR uses to create a genome index. In the second step, STAR maps the user's reads data to the genome index.

Let's create the index now. Remember, for reasons of time we are aligning to a transcriptome rather than a genome today, meaning we only need to provide STAR with the sequences of the transcripts we will be aligning reads to. You can obtain transcriptomes for many model organisms from Ensembl (https://www.ensembl.org/info/data/ftp/index.html).

Task 1: Execute the commands below to create the index:

```bash
mkdir indices
mkdir indices/STAR
STAR --runThreadN 4 --runMode genomeGenerate --genomeDir indices/STAR --genomeFastaFiles Share/2000_reference.transcripts.fa
```

Task 2: What does each of the options we used do? Hint: Use the STAR manual to help you (https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)

Task 3: How would the command we used in Task 1 be different if we were aligning to the genome rather than the transcriptome?

Now that we have created the index, we can perform the mapping step.

Task 4: Try to work out what command you should use to map our trimmed reads (from ERR522959) to the index you created. Use the STAR manual to help you. One you think you know the answer, check whether it matches the solution in the next section and execute the alignment.

Task 5: Try to understand the output of your alignment. Talk to one of the instructors if you need help!

### Solution for STAR Alignment

You can use the folowing commands to perform the mapping step:

```bash
mkdir results
mkdir results/STAR

STAR --runThreadN 4 --genomeDir indices/STAR --readFilesIn Share/ERR522959_1.fastq Share/ERR522959_2.fastq --outFileNamePrefix results/STAR/
```

## Kallisto and Pseudo-Alignment

STAR is a reads aligner, whereas Kallisto is a pseudo-aligner [@bray_2016]. The main difference between aligners and pseudo-aligners is that whereas aligners map reads to a reference, pseudo-aligners map k-mers to a reference.

### What is a k-mer?

A k-mer is a sequence of length k derived from a read. For example, imagine we have a read with the sequence ATCCCGGGTTAT and we want to make 7-mers from it. To do this, we would find the first 7-mer by counting the first seven bases of the read. We would find the second 7-mer by moving one base along, then counting the next seven bases. Figure 2 shows all the 7-mers that could be derived from our read:

![Figure 2: The 7-mers derived from an example read](L1-images/Kmers.png)

### Why map k-mers rather than reads?
There are two main reasons:

1. Pseudo-aligners use k-mers and a computational trick to make pseudo-alignment much faster than traditional aligners. If you are interested in how this is acheived, see (Bray et al., 2017) for details.

2. Under some circumstances, pseudo-aligners may be able to cope better with sequencing errors than traditional aligners. For example, imagine there was a sequencing error in the first base of the read above and the A was actually a T. This would impact on the pseudo-aligners ability to map the first 7-mer but none of the following 7-mers.

### Kallisto's pseudo mode

Kallisto has a specially designed mode for pseudo-aligning reads from single-cell RNA-seq experiments. Unlike STAR, Kallisto psuedo-aligns to a reference transcriptome rather than a reference genome. This means Kallisto maps reads to splice isoforms rather than genes. Mapping reads to isoforms rather than genes is especially challenging for single-cell RNA-seq for the following reasons:

 * Single-cell RNA-seq is lower coverage than bulk RNA-seq, meaning the total amount of information available from reads is reduced.
 * Many single-cell RNA-seq protocols have 3' coverage bias, meaning if two isoforms differ only at their 5' end, it might not be possible to work out which isoform the read came from.
 * Some single-cell RNA-seq protocols have short read lengths, which can also mean it is not possible to work out which isoform the read came from.

Kallisto's pseudo mode takes a slightly different approach to pseudo-alignment. Instead of aligning to isoforms, Kallisto aligns to equivalence classes. Essentially, this means if a read maps to multiple isoforms, Kallisto records the read as mapping to an equivalence class containing all the isoforms it maps to. Instead of using gene or isoform expression estimates in downstream analysis such as clustering, equivalence class counts can be used instead. Figure 3 shows a diagram which helps explain this.

![Figure 3: A diagram explaining Kallisto's Equivalence Classes, taken from https://pachterlab.github.io/kallisto/singlecell.html.](L1-images/TCC.jpg)

Today we will just perform pseudo-alignment with one cell, but Kallisto is also capable of pseudo-aligning multiple cells simultaneously and using information from UMIs. See https://pachterlab.github.io/kallisto/manual for details.

As for STAR, you will need to produce an index for Kallisto before the pseudo-alignment step.

Task 6: Use the below command to produce the Kallisto index. Use the Kallisto manual (https://pachterlab.github.io/kallisto/manual) to work out what the options do in this command.


```bash
mkdir indices/Kallisto
kallisto index -i indices/Kallisto/transcripts.idx Share/2000_reference.transcripts.fa
```

Task 7: Use the Kallisto manual to work out what command to use to perform pseudo-alignment. One you think you know the answer, check whether it matches the solution in the next section and execute the pseudo-alignment.

### Solution to Kallisto Pseudo-Alignment

Use the below command to perform pseudo-alignment


```bash
mkdir results/Kallisto
kallisto pseudo -i indices/Kallisto/transcripts.idx -o results/Kallisto -b batch.txt 
```

See https://pachterlab.github.io/kallisto/manual for instructions on creating batch.txt, or ask an instructor if you get stuck.

### Understanding the Output of Kallisto Pseudo-Alignment

The command above should produce 4 files - matrix.cells, matrix.ec, matrix.tsv and run_info.json.

* matrix.cells contains a list of cell IDs. As we only used one cell, this file should just contain "ERR522959"
* matrix.ec contains information about the equivalence classes used. The first number in each row is the equivalence class ID. The second number(s) correspond to the transcript ID(s) in that equivalence class. For example "10 1,2,3" would mean that equivalence class 10 contains transcript IDs 1,2 and 3. The ID numbers correspond to the order that the transcripts appear in reference.transcripts.fa. Zero indexing is used, meaning transcript IDs 1,2 and 3 correspond to the second, third and fourth transcripts in 2000_reference.transcripts.fa.
* matrix.tsv contains information about how many reads in each cell map to each equivalence class. The first number is the equivalence class ID, as defined in matrix.ec. The second number is the cell ID, where the cell ID corresponds to the order that the cell came in the matrix.cells file. The third number is the number of reads which fall into that equivalence class. For example, "5 1 3" means that 3 reads from cell 1 map to equivalence class 5. Note that zero indexing is used, so cell 1 corresponds to the second line of matrix.cells.
* run_info.json contains information about how Kallisto was executed and can be ignored.



