
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
