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
# https://stackoverflow.com/questions/45706896/attributeerror-module-numpy-has-no-attribute-flip
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

def build_matrix(ref, read):
    global match_cost
    global mismatch_cost
    global gap_cost

    H = np.zeros((len(ref) + 1, len(read) + 1), np.int)
    H[0] = [0, -2, -4, -6, -8, -10, -12]
    H[:, 0] = [0, -2, -4, -6, -8, -10, -12, -14]

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

def backtrack(H_flip, ref_rev, read_rev, pos, SNP, insert, delete):
    global match_cost
    global mismatch_cost
    global gap_cost

    if len(H_flip) == 1 and len(H_flip[0]) == 1:
        return

    if ref_rev[0] == read_rev[0]:
        if H_flip[0, 0] - match_cost == H_flip[1, 1]:
            backtrack(H_flip[1:, 1:], ref_rev[1:], read_rev[1:], pos-1, SNP, insert, delete)
    else:
        if H_flip[0, 0] - mismatch_cost == H_flip[1, 1]:
            SNP.append([ref_rev[0], read_rev[0], pos])
            backtrack(H_flip[1:, 1:], ref_rev[1:], read_rev[1:], pos-1, SNP, insert, delete)

    if H_flip[0, 0] - gap_cost == H_flip[0, 1]:
        insert.append([ref_rev[0], pos])
        backtrack(H_flip[0:, 1:], ref_rev, read_rev[1:], pos, SNP, insert, delete)
    
    if H_flip[0, 0] - gap_cost == H_flip[1, 0]:
        delete.append([ref_rev[0], pos])
        backtrack(H_flip[1:, 0:], ref_rev[1:], read_rev, pos-1, SNP, insert, delete)


def traceback(H, b, b_='', old_i=0):
    # flip H to get index of **last** occurrence of H.max() with np.argmax()
    H_flip = flip(flip(H, 0), 1)
    i_, j_ = np.unravel_index(H_flip.argmax(), H_flip.shape)
    i, j = np.subtract(H.shape, (i_ + 1, j_ + 1))  # (i, j) are **last** indexes of H.max()
    if H[i, j] == 0:
        return b_, j
    b_ = b[j - 1] + '-' + b_ if old_i - i > 1 else b[j - 1] + b_
    return traceback(H[0:i, 0:j], b, b_, i)

# def smith_waterman(a, b, match_score=3, gap_cost=2):
#     a, b = a.upper(), b.upper()
#     H = matrix(a, b, match_score, gap_cost)
#     b_, pos = traceback(H, b)
#     return pos, pos + len(b_)

# def rotate(text, n):
#     first_half = text[0:len(text) - n] 
#     second_half = text[len(text) - n:]
#     return second_half + first_half

# def construct_BWT_matrix(text):
#     rotations_array = []
#     for index, character in enumerate(text):
#         rotations_array.append(rotate(text, index))
#     rotations_array = sorted(rotations_array)
#     return rotations_array

# def get_BWT_string(rotations_array):
#     BWT = ""
#     for rotation in rotations_array:
#         BWT += rotation[len(rotation)-1]
#     return BWT

# def get_BWT(reference_w_dollar):
#     rotations_array = construct_BWT_matrix(reference_w_dollar)
#     return get_BWT_string(rotations_array)



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
    ref_index = dict()
    for index in range(0, len(reference)):
        first_part = reference[index:index+10]
        if ref_index.get(first_part, "DNE") == "DNE":
            ref_index[first_part] = [index]
        else:
            ref_index[first_part].append(index)

    read_order_dict = dict()
    left_ends = []
    for input_read in input_reads:
        read_order_dict[input_read[0]] = input_read[1]
        left_ends.append(input_read[0])
    
    # for left_end in left_ends:
    #     left_end_rev = left_end[::-1]
    #     if ref_index.get(left_end[0:10], "DNE") != "DNE":
    #         for start_index in ref_index[left_end]:
    #             if ref_index == 0
    #     if ref_index.get(left_end_rev[0:10], "DNE") != "DNE":

    read = "ATGAGT"
    ref = "ATGGCGT"

    H = np.zeros((len(ref) + 1, len(read) + 1), np.int)
    H[0] = [0, -2, -4, -6, -8, -10, -12]
    H[:, 0] = [0, -2, -4, -6, -8, -10, -12, -14]

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

    SNP = []
    inserts = []
    deletes = []

    read_rev = read[::-1]
    ref_rev = ref[::-1]

    H_flip = flip(flip(H, 0), 1)

    print(H_flip[0:, 0:])

    backtrack(H_flip, ref_rev, read_rev, len(ref_rev)-1, SNP, inserts, deletes)

    print(SNP)
    print(inserts)
    print(deletes)


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
