import sys
import argparse
import itertools
import numpy as np
import time
import zipfile
import pickle
np.set_printoptions(threshold=sys.maxsize)

from numpy.lib import insert

def parse_reads_file(reads_fn):
    """
    :param reads_fn: the file containing all of the reads
    :return: outputs a list of all paired-end reads

    HINT: This might not work well if the number of reads is too large to handle in memory
    """
    try:
        with open(reads_fn, 'r') as rFile:
            print("Parsing Reads")
            first_line = True
            count = 0
            all_reads = []
            for line in rFile:
                count += 1
                if count % 1000 == 0:
                    print(count, " reads done")
                if first_line:
                    first_line = False
                    continue
                ends = line.strip().split(',')
                all_reads.append(ends)
        return all_reads
    except IOError:
        print("Could not read file: ", reads_fn)
        return None


def parse_ref_file(ref_fn):
    """
    :param ref_fn: the file containing the reference genome
    :return: a string containing the reference genome
    """
    try:
        with open(ref_fn, 'r') as gFile:
            print("Parsing Ref")
            first_line = True
            ref_genome = ''
            for line in gFile:
                if first_line:
                    first_line = False
                    continue
                ref_genome += line.strip()
        return ref_genome
    except IOError:
        print("Could not read file: ", ref_fn)
        return None


"""
    TODO: Use this space to implement any additional functions you might need

"""

# INDEL FUNCTIONS

# Adapted from: https://stackoverflow.com/questions/45706896/attributeerror-module-numpy-has-no-attribute-flip
# .flip() is not being found in my numpy module so I found an implementation for it online
def flip(m, axis):
    if not hasattr(m, 'ndim'):
        m = np.asarray(m)
    indexer = [slice(None)] * m.ndim
    try:
        indexer[axis] = slice(None, None, -1)
    except IndexError:
        raise ValueError("axis=%i is invalid for the %i-dimensional input array"
                         % (axis, m.ndim))
    return m[tuple(indexer)]

match_cost = 3
mismatch_cost = -1
gap_cost = -2

# Adapted from: https://tiefenauer.github.io/blog/smith-waterman/
# Used this resource to find a numpy implementation of the algorithm, but a lot of the code is my own
# Affine gap not implemented
# Build the scoring matrix
def build_matrix(ref, read):
    global match_cost
    global mismatch_cost
    global gap_cost

    # Zero out a 2D array
    H = np.zeros((len(ref) + 1, len(read) + 1), np.int)
    # Initialize the first row/col with values
    count = -gap_cost
    for i in range(1, len(read) + 1):
        H[0, i] = count
        count -= gap_cost
    count = -gap_cost
    for i in range(1, len(ref) + 1):
        H[i, 0] = count
        count -= gap_cost

    for row, col in itertools.product(range(1, H.shape[0]), range(1, H.shape[1])):
        match = 0
        # If the letters match, diagonal is match; else, diagonal is mismatch (SNP)
        if ref[row - 1] == read[col - 1]:
            match = H[row - 1, col - 1] + match_cost
        else:
            match = H[row - 1, col - 1] + mismatch_cost
        # Down is deletion
        deletion = H[row - 1, col] + gap_cost
        # Right is insertion
        insertion = H[row, col - 1] + gap_cost
        H[row, col] = max(match, deletion, insertion)
    return H, H[len(ref), len(read)]

# Backtrack to find the 
def backtrack(H_flip, ref_rev, read_rev, pos, index, SNP, inserts, deletes):
    global match_cost
    global mismatch_cost
    global gap_cost

    # Flip the matrix so the bottom right corner is now the top left corner
    if len(ref_rev) == 0 or len(read_rev) == 0:
        return

    # In each step we prune away the used parts of the H_flip matrix
    # In the insertion case, we only prune the leftmost column of H_flip and do not prune ref_rev, since we move right and did not use a character of ref_rev
    # In the deletion case, we only prune the top row of H_flip and do not prune read_rev, since we move down and did not use a character of read_rev

    # Test if diagonal (match or SNP) was a possible path

    if ref_rev[0] == read_rev[0]:
        if H_flip[0, 0] - match_cost == H_flip[1, 1]:
            backtrack(H_flip[1:, 1:], ref_rev[1:], read_rev[1:], pos-1, index, SNP, inserts, deletes)
    else:
        if H_flip[0, 0] - mismatch_cost == H_flip[1, 1]:
            SNP.append([ref_rev[0], read_rev[0], index + pos])
            backtrack(H_flip[1:, 1:], ref_rev[1:], read_rev[1:], pos-1, index, SNP, inserts, deletes)

    # Test if right (insertion) was a possible path
    if H_flip[0, 0] - gap_cost == H_flip[0, 1]:
        inserts.append([read_rev[0], index + pos])
        backtrack(H_flip[0:, 1:], ref_rev, read_rev[1:], pos, index, SNP, inserts, deletes)

    # Test if down (delete) was a possible path
    if H_flip[0, 0] - gap_cost == H_flip[1, 0]:
        deletes.append([ref_rev[0], index + pos])
        backtrack(H_flip[1:, 0:], ref_rev[1:], read_rev, pos-1, index, SNP, inserts, deletes)

    return

def needleman_wunsch(ref, read, index):
    H, score = build_matrix(ref, read)
    

    # Flip the matrix so the bottom right corner is now the top left corner for backtracking
    H_flip = flip(flip(H, 0), 1)
    read_rev = read[::-1]
    ref_rev = ref[::-1]

    print(H_flip[22:, 0:])

    SNP = []
    inserts = []
    deletes = []

    backtrack(H_flip, ref_rev, read_rev, len(ref_rev)-1, index, SNP, inserts, deletes)
    if len(SNP) > 20:
        SNP = []
    if len(inserts) > 20:
        inserts = []
    if len(deletes) > 20:
        deletes = []
        
    return score, SNP, inserts, deletes

# SNP FUNCTIONS

# Return position of matched 50 bp end in reference genome, otherwise return -1
def test_read(seq_read, ref_index):
    if ref_index.get(seq_read, "DNE") != "DNE":
        return ref_index[seq_read]
    else:
        return -1

# Find if this read has a match somewhere in the reference
# Match is defined by having 1 or fewer SNPs
def find_match(seq_read, ref_dict, mismatches):
    match_found = -1
    for k, v in ref_dict.items():
        for index in range(0, len(k)-1):
            if seq_read[:index] == k[:index] and seq_read[index + 1:] == k[index + 1:]:
                if mismatches.get((k[index], seq_read[index], v + index), "DNE") == "DNE":
                    mismatches[(k[index], seq_read[index], v + index)] = 1
                else: 
                    mismatches[(k[index], seq_read[index], v + index)] += 1
                match_found = v
                break
    return match_found



if __name__ == "__main__":
    ref = 'TGCACGTTAGCCGGTCCAACACGCTAAACTTGTGATACGCCAGAACAGTGGAA'
    read = 'TGCACGTTAGCCGGTCCAACACGCTTGTGATACGCCAGAACAGTGGAA'

    score, SNP, insertions, deletions = needleman_wunsch(ref, read, 67)
    print(SNP)
    print(insertions)
    print(deletions)
