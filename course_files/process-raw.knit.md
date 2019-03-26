---
output: html_document
---


## File formats

### FastQ
FastQ is the most raw form of scRNASeq data you will encounter. All scRNASeq 
protocols are sequenced with paired-end sequencing. Barcode sequences may occur 
in one or both reads depending on the protocol employed. However, protocols 
using unique molecular identifiers (UMIs) will generally contain one read with 
the cell and UMI barcodes plus adapters but without any transcript sequence. 
Thus reads will be mapped as if they are single-end sequenced despite actually
being paired end. 

FastQ files have the format:

```eval
>ReadID
READ SEQUENCE
+
SEQUENCING QUALITY SCORES
```

### BAM

BAM file format stores mapped reads in a standard and efficient manner. The 
human-readable version is called a SAM file, while the BAM file is the highly
compressed version. BAM/SAM files contain a header which typically includes  
information on the sample preparation, sequencing and mapping; and a tab-separated row for each individual alignment of each read. 

Alignment rows employ a standard format with the following columns:

(1) QNAME : read name (generally will include UMI barcode if applicable)

(2) FLAG : number tag indicating the "type" of alignment, [link](https://broadinstitute.github.io/picard/explain-flags.html) to explanation of all possible "types"

(3) RNAME : reference sequence name (i.e. chromosome read is mapped to).

(4) POS : leftmost mapping position

(5) MAPQ : Mapping quality

(6) CIGAR : string indicating the matching/mismatching parts of the read (may include soft-clipping).

(7) RNEXT : reference name of the mate/next read

(8) PNEXT : POS for mate/next read

(9) TLEN : Template length (length of reference region the read is mapped to)

(10) SEQ : read sequence

(11) QUAL : read quality

BAM/SAM files can be converted to the other format using 'samtools':


```bash
samtools view -S -b file.sam > file.bam
samtools view -h file.bam > file.sam
```

Some sequencing facilities will automatically map your reads to the a standard 
genome and deliver either BAM or CRAM formatted files. Generally they will not 
have included ERCC sequences in the genome thus no ERCC reads will be mapped in
the BAM/CRAM file. To quantify ERCCs (or any other genetic alterations) or if 
you just want to use a different alignment algorithm than whatever is in the 
generic pipeline (often outdated), then you will need to convert the BAM/CRAM 
files back to FastQs:

BAM files can be converted to FastQ using bedtools. To ensure a single copy for
multi-mapping reads first sort by read name and remove secondary alignments
using samtools. [Picard](https://broadinstitute.github.io/picard/index.html) 
also contains a method for converting BAM to FastQ files. 


```bash
# sort reads by name
samtools sort -n original.bam -o sorted_by_name.bam
# remove secondary alignments
samtools view -b -F 256 sorted_by_name.bam -o primary_alignment_only.bam
# convert to fastq
bedtools bamtofastq -i primary_alignment_only.bam -fq read1.fq -fq2 read2.fq
```
### CRAM

[CRAM](https://www.ebi.ac.uk/ena/software/cram-usage) files are similar to BAM files only they contain information in the header 
to the reference genome used in the mapping in the header. This allow the bases
in each read that are identical to the reference to be further compressed. CRAM
also supports some lossy data compression approaches to further optimize storage
compared to BAMs. CRAMs are mainly used by the Sanger/EBI sequencing facility.

CRAM and BAM files can be interchanged using the lastest version of samtools (>=v1.0). 
However, this conversion may require downloading the reference genome into cache.
Alternatively, you may pre-download the correct reference either from metadata in the header
of the CRAM file, or from talking to whomever generated the CRAM and specify that file using '-T'
Thus we recommend setting a specific cache location prior to doing this:


```bash
export REF_CACHE=/path_to/cache_directory_for_reference_genome
samtools view -b -h -T reference_genome.fasta file.cram -o file.bam
samtools view -C -h -T reference_genome.fasta file.bam -o file.cram
```

### Mannually Inspecting files

At times it may be useful to mannual inspect files for example to check the metadata in headers that the files 
are from the correct sample. 'less' and 'more' can be used to inspect any text files from the command line.
By "pipe-ing" the output of samtools view into these commands using '|' we check each of these file types without having to save
multiple copies of each file.


```bash
less file.txt
more file.txt
# counts the number of lines in file.txt
wc -l file.txt
samtools view -h file.[cram/bam] | more
# counts the number of lines in the samtools output
samtools view -h file.[cram/bam] | wc -l
```

__Exercises__

You have been provided with a small cram file: EXAMPLE.cram 

Task 1: How was this file aligned? What software was used? What was used as the genome? (Hint: check the header)

Task 2: How many reads are unmapped/mapped? How total reads are there? How many secondary alignments are present? (Hint: use the FLAG)

Task 3: Convert the CRAM into two Fastq files. Did you get exactly one copy of each read? (name these files "10cells_read1.fastq" "10cells_read2.fastq")

If you get stuck help information for each piece of software can be displayed 
by entering running the command "naked" - e.g. 'samtools view', 'bedtools'


__Answer__



























