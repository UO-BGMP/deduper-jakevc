#!/usr/bin/env python3.6

"""
Jake VanCampen - Jan. 2018

This script sorts a sam file of uniquely mapped RNA-seq reads by base
position, then removes all PCR duplicates determined by chromosome,
umi, and position.

**** A current version of samtools is required to run this script ****
"""

import argparse
import subprocess
import re
from collections import defaultdict


def get_umi(alignment):
    '''Returns the unique molecular index (UMI) for the alignment.'''
    umi = alignment[0].split(':')[-1]
    return umi


def get_pos(alignment):
    '''Returns the 1-based, leftmost mapping position from the 4th
    field of the alignment section.
    '''
    pos = alignment[3]
    return int(pos)


def get_tlen(alignment):
    tlen = alignment[8]
    return int(tlen)


def get_chrm(alignment):
    '''Returns chromosome number from alignment section'''
    chrm = alignment[2]
    return chrm


def num_clipped(alignment):
    '''Parses the cigar string for each alignment, returning the number
    of soft-clipped bases.
    '''
    cigar = alignment[5]

    # determine number soft clipped bases
    soft_match = re.match('^\d+S', cigar)
    if soft_match:
        soft_bases = soft_match.group()

        # extract the number of soft clipped bases
        num_soft = int(re.match('\d+', soft_bases).group())

    else:
        num_soft = 0

    return num_soft


def correct_pos(alignment):
    '''Corrects mapping position by num_clipped.'''
    if num_clipped(alignment) != 0:
        # left-most mapping postition
        pos = get_pos(alignment) - num_clipped(alignment)
    else:
        pos = get_pos(alignment)
    return pos


# function for future paired-end functionality
def check_strand(alignment):
    '''Takes bitwise flag and returns + or - indicating strand'''
    FLAG = alignment[1]
    if ((FLAG & 16) != 16):
        return "+"
    else:
        return "-"


def get_args():
    '''Define and return command line options.'''
    parser = argparse.ArgumentParser(
        prog='vancampen_deduper.py',
        description='This program removes PCR duplicates from a SAM \
        file of uniquely mapped, single-end reads. The program first \
        sorts the SAM file based on 1-based, leftmost mapping position\
        using SAM tools. Please ensure that samtools is installed on\
        your machine before running this script')

    parser.add_argument('-i', '--infile',
                        help='specify input file (abs path)',
                        required=True,
                        type=argparse.FileType('rt',
                                               encoding='UTF-8 '))

    parser.add_argument('-u', '--umifile',
                        help='file containing list of known UMIs',
                        required=True,
                        type=argparse.FileType('rt',
                                               encoding='UTF-8 '))

    parser.add_argument('-p', '--paired', action='store_true',
                        help='flag indicating paired-end alignments')

    return parser.parse_args()


# return command line arguments
args = get_args()

# cache infile
infile = args.infile.name

outfile = infile.replace('.sam', '_dedup.sam')

# name the temp file
samtemp = infile + '.sort'

# name the output file
samout = infile.replace('.sam', '_sorted.sam')

# sort the samfile using samtools sort
subprocess.run(['samtools', 'sort', '-T',
               samtemp, '-o', samout, infile])

# catch paired end flag
if args.paired:
    raise ValueError('No paired-end functionality implimented')

# generate list of umi's from umi file
umifile = args.umifile.name
umilist = []
with open(umifile, 'r') as umi:
    for line in umi:
        line = line.strip().split()
        umilist += line

# match dict
match_dict = defaultdict()

# firstentry to populate dict
firstentry = True

# read SAM file alignments
with open(samout, 'r') as fh, open(outfile, 'w') as of:

    for line in fh:

        # strip whitespace and newlines
        line = line.strip()

        # write header lines straight to output file
        if line[0] == '@':
            of.write(line+'\n')

        elif firstentry:

            # split by field
            line = line.split()

            # retrieve umi
            umi = get_umi(line)

            # store current line position
            current_pos = correct_pos(line)

            # chromosome number
            chrm = get_chrm(line)

            # generate first entry
            if umi in umilist:
                match_dict[(umi, current_pos, chrm)] = line

            # endcase
            firstentry = False

            # update previous position
            prev_position = current_pos

        else:
            # split the lines by field
            line = line.split()

            # retrieve umi
            umi = get_umi(line)

            # correct alignmnet postition
            current_pos = correct_pos(line)

            # chromosome number
            chrm = get_chrm(line)

            # length of read
            tlen = get_tlen(line)

            # attributes to compare
            attributes = (umi, current_pos, chrm)

            # only add a new entry if it is valid
            # and not already present.
            if current_pos >= (tlen + prev_position) and\
               umi in umilist and\
               attributes not in match_dict:
                match_dict[attributes] = line

            # update current position
            prev_position = current_pos

    # write retained entries to output file
    for key in match_dict:
        of.write('\t'.join(match_dict[key])+'\n')
