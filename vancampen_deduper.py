#!/usr/bin/env python3.6

"""
This script removes all PCR duplicates from a SAM file
containing uniquely mapped RNA-seq reads.
"""

import argparse
import re
from collections import defaultdict


def get_umi(alignment):
    '''Returns the unique molecular index (UMI) for the alignment.'''
    umi = alignment[0].split(':')[1]
    return umi


def get_pos(alignment):
    '''Returns the 1-based, leftmost mapping position from the 4th
    field of the alignment section.
    '''
    pos = alignment[3]
    return int(pos)


def get_tlen(alignment):
    tlen = alignment[8]
    return tlen


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


def decode_flag():
    '''Verifies mapping conditions.'''


def umi_check():
    '''
    Returns TRUE if UMI in list of known UMIs.
    '''


def samtools_sort():
    '''Sorts the SAM file using samtools.'''

def get_args():
    '''Define and return command line options.'''
    parser = argparse.ArgumentParser(prog='vancampen_deduper.py',
                                     description='Reference-based PCR\
                                     deduplicator of uniquely mapped\
                                     alignments in a SAM file.')

    parser.add_argument('-i', '--infile',
                        help='specify input file (abs path)',
                        required=True,
                        type=argparse.FileType('rt',
                                               encoding='UTF-8 '))

    parser.add_argument('-u', '--umifile',
                        help='file containing list of known UMIs',
                        required=False,
                        type=argparse.FileType('rt',
                                               encoding='UTF-8 '))

    parser.add_argument('-o', '--outfile',
                        help='specify output file (path)',
                        required=False,
                        type=argparse.FileType('wt',
                                               encoding='UTF-8 '))

    parser.add_argument('-p', '--paired',
                        help='flag indicating paired-end alignments',
                        required=False,
                        type=bool
                        )
    return parser.parse_args()


# return command line arguments
args = get_args()

# cache infile
infile = args.infile.name


# ... umifile
umifile = args.umifile


# match dict
match_dict = defaultdict()

firstline = True

# read SAM file alignments
with open(infile, 'r') as fh:
    fh.readline()

    if firstline:
        for line in fh:
            line = line.strip().split('  ')

            # retrieve umi
            umi = get_umi(line)

            current_pos = correct_pos(line)

            # chromosome number
            chrm = get_chrm(line)

            match_dict[(umi, current_pos, chrm)] = line

        firstline = False
    else:
        for line in fh:
            line = line.strip().split('  ')

            # retrieve umi
            umi = get_umi(line)

            # correct alignmnet postition
            pos = correct_pos(line)

            # chromosome number
            chrm = get_chrm(line)

            # length of template
            tlen = get_tlen(line)

            if pos < (0.75*tlen + current_pos):
                match_dict[(umi, pos, chrm)] = line
            else:
                pass


for key in match_dict:
    print(f'{key}: {match_dict[key]}')
