## De-duper

### Part 1
Write up a strategy for writing a Reference based PCR duplicate removal tool. Be sure to include an example input and output SAM file.

### Background
To prepare for downstream analysis of next-generation sequencing data it is necessary to remove PCR duplicates. PCR duplicates arise from PCR amplification during library prep. Given a sample of adapter-ligated library fragments (molecules of DNA or cDNA), if there is an over-represented library fragment, there is a greater likelihood that this fragment will end up being sequenced twice because it generated two separate clusters on the flow cell (for Illumina reads). Over representation of a fragment can occur for a number of reasons. Having too few starting fragments can lead to over-representation because more PCR cycles need to be run in order to have enough DNA to sequence. Alternatively, a large variation in fragment size can lead to PCR bias because shorter fragments amplify better than longer ones.


### What does it look like?
The PCR duplicates appear as the same sequences with the same unique molecular indices (UMIs) or inline 'randomer'. The sequences can be characterized by alignment to the same chromosome, position, and strand, of the genome. Complications that arise from sequencing can make it hard to identify PCR duplicates and remove them programatically.


### What tools are available?
Various tools have been built to de-duplicate sequencing reads. [Samtools rmdup](http://www.htslib.org/doc/samtools.html) identifies duplicates by having the same alignment position from a SAM file. [Picard MarkDuplicates](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates) marks duplicates and reports the number found for both single- and paired-end data; duplicate sequences can be removed if specified. MarkDuplicates also accounts for any soft clipping in the reads. [UMI-tools](https://github.com/CGATOxford/UMI-tools) is another tool used to de-duplicate reads with the same alignment position, accounting for UMIs and any simple soft clipping of the reads. Other slower tools for duplicate removal include [Stacks clone_filter](http://catchenlab.life.illinois.edu/stacks/comp/clone_filter.php), and [FastUniq](https://bioconda.github.io/recipes/fastuniq/README.html).


### Strategize and Algorithm
Here I will begin developing a strategy for writing my own PCR duplicate removal tool in Python that de-duplicates a sorted SAM file of aligned reads.


### Begin
We start with a SAM file of uniquely mapped reads. First `samtools sort` will be used to sort the SAM file by base position, and return the sorted file in SAM format. A samtools command of this format could be used:

```
samtools sort -O sam -T sample.sort -o sample.sort.sam sample.sam
```

Now that the same file is sorted. We need to account for simple soft-clipping of the aligned reads.
