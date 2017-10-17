## De-duper

### Part 1
Write up a strategy for writing a Reference based PCR duplicate removal tool. Be sure to include an example input and output SAM file.

### Background
To prepare for downstream analysis of next-generation sequencing data it is necessary to remove PCR duplicates. PCR duplicates arise from PCR amplification during library prep. Given a sample of adapter-ligated library fragments (molecules of DNA or cDNA), if there is an over-represented library fragment, there is a greater likelihood that this fragment will end up being sequenced twice because it generated two separate clusters on the flow cell (for Illumina reads). Over representation of a fragment can occur for a number of reasons. Having too few starting fragments can lead to over-representation because more PCR cycles need to be run in order to have enough DNA to sequence. Alternatively, a large variation in fragment size can lead to PCR bias because shorter fragments amplify better than longer ones.


### What does it look like?
The PCR duplicates appear as the same sequences with the same unique molecular indices (UMIs) or inline 'randomer'. The sequences can be characterized by alignment to the same chromosome, position, and strand, of the genome. Complications that arise from sequencing can make it hard to identify PCR duplicates and remove them programatically such as soft-clipping and low quality base calls.


### What tools are available?
Various tools have been built to de-duplicate sequencing reads. [Samtools rmdup](http://www.htslib.org/doc/samtools.html) identifies duplicates by having the same alignment position from a SAM file. [Picard MarkDuplicates](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates) marks duplicates and reports the number found for both single- and paired-end data; duplicate sequences can be removed if specified. MarkDuplicates also accounts for any soft clipping in the reads. [UMI-tools](https://github.com/CGATOxford/UMI-tools) is another tool used to de-duplicate reads with the same alignment position, accounting for UMIs and any simple soft clipping of the reads. Other reference-free tools for duplicate removal include [Stacks clone_filter](http://catchenlab.life.illinois.edu/stacks/comp/clone_filter.php), and [FastUniq](https://bioconda.github.io/recipes/fastuniq/README.html), they tend to be more computationally expensive (slow).


### Strategize and Algorithm
Here I will begin developing a strategy for writing my own reference-based PCR duplicate removal tool in Python that de-duplicates a sorted SAM file of aligned reads.


### Starting materials
We start with a SAM file of uniquely mapped reads. First `samtools sort` will be used to sort reads in the SAM file by chromosome and left-most base position, and return the sorted file in SAM format. A samtools command of this format could be used:

```
samtools sort -T sample.sort -o sample.sort.sam sample.sam
```

Make sure to checkout how this sorts the reads with some UNIX commands

### High-level functions

```
get_umi: parses the header information and returns the UMI

get_pos: parse the forth field of the alignment section, returning the 1-based left-most mapping position of each read.

get_chrm: returns chromosome from alignment section (RNAME: field 3)

num_clipped: parses CIGAR string, and determine the level of soft clipping for the reads. If (regex to match number followed by S), returns the number of bases (if any) soft clipped.


```

### Algorithm

```
Take in file <example>.sam as an argument.

Step through alignments one at a time:

Store the first alignment in dict. key: (UMI, POS, CHRM), value: (alignment)
store POS of this read in a variable

if UMI in 96_UMI:

  if not soft clipped:
    add alignment if key (UMI, POS, CHRM) not in dict, and POS is not within +_ 3/4 read-length of previous POS.
    update previous POS to current POS.

  else:
    adjust POS=POS-num_clipped
    if not (UMI, POS, CHRM) in dict and POS is not within +_ 3/4 read-length of previous POS, update dict.
    update previous POS to current POS.

else:
  increment bad UMI counter

write out alignment dict to a file:  <example>.dedup.sam


```

### Additional functionality

* Account for paired-end reads (this could be done by throwing out the mate-pair of the duplicate read, a function would be needed to find the mate-pair and skip that read)

* Maybe include complex clipping?
