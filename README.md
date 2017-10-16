## De-duper

### Part 1
Write up a strategy for writing a Reference based PCR duplicate removal tool. Be sure to include an example input and output SAM file.

### Background
To prepare for downstream analysis of next-generation sequencing data it is necessary to remove PCR duplicates. PCR duplicates arise from PCR amplification during library prep. Given a sample of adapter-ligated library fragments (molecules of DNA or cDNA), if there is an over-represented library fragment, there is a greater likelihood that this fragment will end up being sequenced twice because it generated two separate clusters on the flow cell (for Illumina reads). Over representation of a fragment can occur for a number of reasons. Having too few starting fragments can lead to over-representation because more PCR cycles need to be run in order to have enough DNA to sequence. Alternatively, a large variation in fragment size can lead to PCR bias because shorter fragments amplify better than longer ones.


### What does it look like?
The PCR duplicates appear as the same sequences with the same unique molecular indices (UMIs) or inline 'randomer'. The sequences can be characterized by alignment to the same chromosome, position, and strand, of the genome. Complications that arise from sequencing can make it hard to identify PCR duplicates and remove them programatically such as soft-clipping and low quality base calls.


### What tools are available?
Various tools have been built to de-duplicate sequencing reads. [Samtools rmdup](http://www.htslib.org/doc/samtools.html) identifies duplicates by having the same alignment position from a SAM file. [Picard MarkDuplicates](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates) marks duplicates and reports the number found for both single- and paired-end data; duplicate sequences can be removed if specified. MarkDuplicates also accounts for any soft clipping in the reads. [UMI-tools](https://github.com/CGATOxford/UMI-tools) is another tool used to de-duplicate reads with the same alignment position, accounting for UMIs and any simple soft clipping of the reads. Other slower tools for duplicate removal include [Stacks clone_filter](http://catchenlab.life.illinois.edu/stacks/comp/clone_filter.php), and [FastUniq](https://bioconda.github.io/recipes/fastuniq/README.html).


### Strategize and Algorithm
Here I will begin developing a strategy for writing my own PCR duplicate removal tool in Python that de-duplicates a sorted SAM file of aligned reads.


### Starting materials
We start with a SAM file of uniquely mapped reads. First `samtools sort` will be used to sort reads in the SAM file by left-most base position, and return the sorted file in SAM format. A samtools command of this format could be used:

```
samtools sort -T sample.sort -o sample.sort.sam sample.sam
```

Make sure to douible check that this sorts reads by left-most base position

### High-level functions

A function is needed to parse the forth field of the alignment section, returning the 1-based left-most mapping position of each read.


A function will be needed to parse the CIGAR string, and determine the level of soft clipping for the reads. If (regex to match number followed by S), returning the number of bases (if any) soft clipped.


A function is needed to apply the soft-clip matching to reads that have the same left-most base position, but have soft clipped nucleotides.


### Algorithm

Take in file <example>.sam

Remove SAM header

Clean reads by removing reads with RNAME (field 3) set to "\*"

Don't consider reads with flag && 4 == True (only consider mapped reads)

Step through reads one at a time, storing each read to a dictionary with the left-most base position as they key, and a list (sequence, qscore, sclipped) as the value.

If that left-most base position is already in the dictionary, determine if the sequences are identical.

If they are, retain the sequence of highest quality. If not identical, If soft clipped, determine if there is exact match past soft clipped bases, if match retain highest quality read.

Return a file, <example>.dedup.sam
