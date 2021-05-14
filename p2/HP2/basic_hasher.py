import sys
import argparse
import itertools
import numpy as np
import time
import zipfile

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
    return H

# Backtrack to find the 
def backtrack(H_flip, ref_rev, read_rev, pos, SNP, insert, delete):
    global match_cost
    global mismatch_cost
    global gap_cost

    # Flip the matrix so the bottom right corner is now the top left corner
    if len(H_flip) == 1 and len(H_flip[0]) == 1:
        return

    # In each step we prune away the used parts of the H_flip matrix
    # In the insertion case, we only prune the leftmost column of H_flip and do not prune ref_rev, since we move right and did not use a character of ref_rev
    # In the deletion case, we only prune the top row of H_flip and do not prune read_rev, since we move down and did not use a character of read_rev

    # Test if diagonal (match or SNP) was a possible path
    if ref_rev[0] == read_rev[0]:
        if H_flip[0, 0] - match_cost == H_flip[1, 1]:
            backtrack(H_flip[1:, 1:], ref_rev[1:], read_rev[1:], pos-1, SNP, insert, delete)
    else:
        if H_flip[0, 0] - mismatch_cost == H_flip[1, 1]:
            SNP.append([ref_rev[0], read_rev[0], pos])
            backtrack(H_flip[1:, 1:], ref_rev[1:], read_rev[1:], pos-1, SNP, insert, delete)

    # Test if right (insertion) was a possible path
    if H_flip[0, 0] - gap_cost == H_flip[0, 1]:
        insert.append([ref_rev[0], pos])
        backtrack(H_flip[0:, 1:], ref_rev, read_rev[1:], pos, SNP, insert, delete)
    
    # Test if down (delete) was a possible path
    if H_flip[0, 0] - gap_cost == H_flip[1, 0]:
        delete.append([ref_rev[0], pos])
        backtrack(H_flip[1:, 0:], ref_rev[1:], read_rev, pos-1, SNP, insert, delete)

    return

def needleman_wunsch(ref, read, SNP_dict, insertion_dict, deletion_dict):
    H = build_matrix(ref, read)
    
    # Flip the matrix so the bottom right corner is now the top left corner for backtracking
    H_flip = flip(flip(H, 0), 1)
    read_rev = read[::-1]
    ref_rev = ref[::-1]

    backtrack(H_flip, ref_rev, read_rev, len(ref_rev)-1, SNP_dict, insertion_dict, deletion_dict)
    return

# SNP FUNCTIONS

# Return position of matched 50 bp end in reference genome, otherwise return -1
def test_read(seq_read, ref_index):
    if ref_index.get(seq_read, "DNE") != "DNE":
        return ref_index[seq_read]
    else:
        return -1

def find_match(seq_read, ref_dict, mismatches):
    match_found = -1
    for k, v in ref_dict.items():
        for index in range(0, len(k)-1):
            if seq_read[:index] == k[:index] and seq_read[index + 1:] == k[index + 1:]:
                if mismatches.get((k[index], seq_read[index], v + index), "DNE") == "DNE":
                    mismatches[(k[index], seq_read[index], v + index)] = 1
                else: 
                    mismatches[(k[index], seq_read[index], v + index)] += 1
                match_found = 1
                continue
    return match_found



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='basic_hasher.py takes in data for homework assignment 2 consisting '
                                     'of a genome and a set of reads and aligns the reads to the reference genome, '
                                     'then calls SNPS and indels based on this alignment.')
    parser.add_argument('-g', '--referenceGenome', required=True, dest='reference_file',
                        help='File containing a reference genome.')
    parser.add_argument('-r', '--reads', required=True, dest='reads_file',
                        help='File containg sequencing reads.')
    parser.add_argument('-o', '--outputFile', required=True, dest='output_file',
                        help='Output file name.')
    parser.add_argument('-t', '--outputHeader', required=True, dest='output_header',
                        help='String that needs to be outputted on the first line of the output file so that the\n'
                             'online submission system recognizes which leaderboard this file should be submitted to.\n'
                             'This HAS to be one of the following:\n'
                             '1) practice_W_3_chr_1 for 10K length genome practice data\n'
                             '2) practice_E_1_chr_1 for 1 million length genome practice data\n'
                             '3) hw2undergrad_E_2_chr_1 for project 2 undergrad for-credit data\n'
                             '4) hw2grad_M_1_chr_1 for project 2 grad for-credit data\n')
    args = parser.parse_args()
    reference_fn = args.reference_file
    reads_fn = args.reads_file

    input_reads = parse_reads_file(reads_fn)
    if input_reads is None:
        sys.exit(1)
    reference = parse_ref_file(reference_fn)
    if reference is None:
        sys.exit(1)

    """
        TODO: Call functions to do the actual read alignment here

    """
    # Remember that for each paired read, the left (first) end always precedes the right (second) end in the reference genome.
    # So, we handle each left end, but after we find (if we find) a match for the left end, we will match the right end 
    # using read_order_dict to get the right end immediately.
    
    # REFERENCE INDEX
    # Create dictionary of all substrings of length 50 (key) and their positions (value)
    ref_index = dict()
    final_index = len(reference) - len(input_reads[0][0])
    for index in range(0, final_index + 1):
        ref_index[reference[index:index+50]] = index

    # ref_index = dict()
    # for index in range(0, len(reference)):
    #     first_part = reference[index:index+10]
    #     if ref_index.get(first_part, "DNE") == "DNE":
    #         ref_index[first_part] = [index]
    #     else:
    #         ref_index[first_part].append(index)

    # TARGET DICTIONARIES
    # Preliminary SNPs
    mismatches = dict()
    # Preliminary insertions
    inserts = dict()
    # Preliminary deletes
    deletes = dict()

    # Create dictionary of all mutations (key) and the number of times they are observed (value)
    for pair in input_reads:
        # Case 1: check if left end is in reference dictionary
        first_pos = test_read(pair[0], ref_index)
        # Case 2: check if reverse of left end is in reference dictionary
        first_rev_pos = test_read(pair[0][::-1], ref_index)

        # If found a match for left end in Case 1
        if first_pos != -1:
            # If the right end of the read also matches, this read is a perfect match, continue
            if test_read(pair[1][::-1], ref_index) != -1:
                continue
            # Find the mutations in the right end of the read (90 - 110 bp ahead)
            else:
                find_match(pair[1][::-1], ref_index, mismatches)
        else:
            match_found = find_match(pair[0], ref_index, mismatches)
            
            if match_found == 1:
                # If the right end of the read also matches, this read is a perfect match, continue
                if test_read(pair[1][::-1], ref_index) != -1:
                    continue
                # Find the mutations in the right end of the read (90 - 110 bp ahead)
                else:
                    find_match(pair[1][::-1], ref_index, mismatches)
            else:
                continue

        # If found a match for left end in Case 2
        if first_rev_pos != -1:
            # If the right end of the read also matches, this read is a perfect match, continue
            if test_read(pair[1], ref_index) != -1:
                continue
            # Find the mutations in the right end of the read (90 - 110 bp ahead)
            else:
                find_match(pair[1], ref_index, mismatches)
        else:
            match_found = find_match(pair[0][::-1], ref_index, mismatches)
            
            if match_found == 1:
                # If the right end of the read also matches, this read is a perfect match, continue
                if test_read(pair[1], ref_index) != -1:
                    continue
                # Find the mutations in the right end of the read (90 - 110 bp ahead)
                else:
                    find_match(pair[1], ref_index, mismatches)
            else:
                continue


    # Iterate through the dictionary of mutations and pick only the ones that are observed above a certain threshold
    true_snps = []
    for k, v in mismatches.items():
        # Adjust what v > ? for changing confidence level of SNP (v is the number of occurrences we observed this SNP)
        if v > 1:
            true_snps.append([k[0], k[1], k[2]])

    true_snps = sorted(true_snps, key=lambda x: x[2])
    print(true_snps)




    snps = [['A', 'G', 3425]]
    insertions = [['ACGTA', 12434]]
    deletions = [['CACGG', 12]]

    output_fn = args.output_file
    zip_fn = output_fn + '.zip'
    with open(output_fn, 'w') as output_file:
        output_file.write('>' + args.output_header + '\n>SNP\n')
        for x in snps:
            output_file.write(','.join([str(u) for u in x]) + '\n')
        output_file.write('>INS\n')
        for x in insertions:
            output_file.write(','.join([str(u) for u in x]) + '\n')
        output_file.write('>DEL\n')
        for x in deletions:
            output_file.write(','.join([str(u) for u in x]) + '\n')
    with zipfile.ZipFile(zip_fn, 'w') as myzip:
        myzip.write(output_fn)
