import sys
import argparse
import itertools
import numpy as np
import time
import zipfile
import pickle

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

# Backtrack to find the different paths that yield the highest score
def backtrack(H_flip, ref_rev, read_rev, pos, index, SNP, inserts, deletes):
    global match_cost
    global mismatch_cost
    global gap_cost

    # If one of the lists is getting too large, toss this matching out
    if len(SNP) > 20 or len(inserts) > 20 or len(deletes) > 20:
        return 1

    # Flip the matrix so the bottom right corner is now the top left corner
    if len(ref_rev) == 0 or len(read_rev) == 0:
        return -1

    # In each step we prune away the used parts of the H_flip matrix
    # In the insertion case, we only prune the leftmost column of H_flip and do not prune ref_rev, since we move right and did not use a character of ref_rev
    # In the deletion case, we only prune the top row of H_flip and do not prune read_rev, since we move down and did not use a character of read_rev

    # Test if diagonal (match or SNP) was a possible path
    if ref_rev[0] == read_rev[0]:
        if H_flip[0, 0] - match_cost == H_flip[1, 1]:
            abort_status = backtrack(H_flip[1:, 1:], ref_rev[1:], read_rev[1:], pos-1, index, SNP, inserts, deletes)
            if abort_status == 1:
                return abort_status
    else:
        if H_flip[0, 0] - mismatch_cost == H_flip[1, 1]:
            SNP.append([ref_rev[0], read_rev[0], index + pos])
            abort_status = backtrack(H_flip[1:, 1:], ref_rev[1:], read_rev[1:], pos-1, index, SNP, inserts, deletes)
            if abort_status == 1:
                return abort_status

    # Test if right (insertion) was a possible path
    if H_flip[0, 0] - gap_cost == H_flip[0, 1]:
        inserts.append([read_rev[0], index + pos])
        abort_status = backtrack(H_flip[0:, 1:], ref_rev, read_rev[1:], pos, index, SNP, inserts, deletes)
        if abort_status == 1:
            return abort_status

    # Test if down (delete) was a possible path
    if H_flip[0, 0] - gap_cost == H_flip[1, 0]:
        deletes.append([ref_rev[0], index + pos])
        abort_status = backtrack(H_flip[1:, 0:], ref_rev[1:], read_rev, pos-1, index, SNP, inserts, deletes)
        if abort_status == 1:
            return abort_status

    return -1

# Executes the needleman-wunsch algorithm and returns lists of SNPs, insertions, and deletions
# Does not handle multiple base pair indels, but from my observations it does return a good chunk of them as single bp insertions or deletions
# (especially deletions) correctly, but gives false positives at positions beyond the ends of the indel
# If I had time, implement smith-waterman instead as TA Brandon Jew said it would help remove the false positives at the ends and then I could simply
# loop through the single base pair indels and build multi-base-pair indels
def needleman_wunsch(ref, read, index):
    # Get the scoring matrix and the end highest score
    H, score = build_matrix(ref, read)
    
    # Flip the matrix so the bottom right corner is now the top left corner for backtracking
    H_flip = flip(flip(H, 0), 1)
    read_rev = read[::-1]
    ref_rev = ref[::-1]

    # Store the returned targets for each category
    SNP = []
    inserts = []
    deletes = []

    # Backtrack through scoring matrix and find any target indels/SNPs and store them
    abort_status = backtrack(H_flip, ref_rev, read_rev, len(ref_rev)-1, index, SNP, inserts, deletes)
    # If we encounter an especially large list, we probably hit an alternating substring e.g. GCGCGCGCGCGCGCGCGCGC, so toss this matching
    if abort_status == 1:
        SNP = []
        inserts = []
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
# Match is defined by having 2 or fewer SNPs
def find_match(ref, seq_read, substr_ref_dict, mismatches):
    match_found = -1
    match_count = 0

    # Get the indices in the reference string where the thirds of the read string perfectly match
    indices = [-1, -1, -1]
    if substr_ref_dict.get(seq_read[0:16], "DNE") != "DNE":
        indices[0] = substr_ref_dict[seq_read[0:16]]
        match_count += 1
    if substr_ref_dict.get(seq_read[16:32], "DNE") != "DNE":
        indices[1] = substr_ref_dict[seq_read[16:32]]
        match_count += 1
    if substr_ref_dict.get(seq_read[32:48], "DNE") != "DNE":
        indices[2] = substr_ref_dict[seq_read[32:48]]
        match_count += 1

    # If we had at least one matching third, we have a potential match (SNP could lie in the other 2 thirds at most)
    if match_count > 0:
        # Need to get the full reference string to compare the read to, so we get the starting and ending indices for the full 50 bp length string
        ref_start_index = 0
        ref_end_index = 0

        # If we only have 1 matching third, this means this string has at least 2 SNPs, so we produce the starting and ending indices from one index
        if match_count == 1:
            if indices[0] != -1:
                ref_start_index = indices[0]
                ref_end_index = indices[0] + 50
            elif indices[1] != -1:
                ref_start_index = indices[1] - 16
                ref_end_index = indices[1] + 34
            elif indices[2] != -1:
                ref_start_index = indices[2] - 32
                ref_end_index = indices[2] + 18
        # If we have 2 or 3 matching thirds, we get the starting and ending indices from 2 of the matching thirds
        # Note, starting index will be either the first or second third (one of these two must have been a match)
        # Likewise with ending index and second and third thirds
        # Even though we only handle matching 48 base pairs, an SNP in the last two base pairs is still handled in the case that all thirds match
        else:
            if indices[0] != -1:
                ref_start_index = indices[0]
            elif indices[1] != -1:
                ref_start_index = indices[1] - 16

            if indices[2] != 1:
                ref_end_index = indices[2] + 18
            elif indices[1] != 1:
                ref_end_index = indices[1] + 34

            # If the reference string is longer than 50, we're not dealing with SNPs so it's not a potential match
            if ref_end_index - ref_start_index > 50:
                return match_found

        # Iterate through the reference string (potential match) and read, and find any mismatching characters at the respective index (possible SNP)
        # Store them in potential_SNP; we only allow up 2 SNPs for it to be a match, so if len(potential_SNP) > 2, we reject this and say no match found
        potential_match = ref[ref_start_index:ref_end_index]
        potential_SNP = []
        for index in range(0, len(potential_match)):
            if seq_read[index] != potential_match[index]:
                if len(potential_SNP) == 2:
                    return match_found
                else:
                    potential_SNP.append(index)

        # Put the SNPs in our mismatch index and then return index in the full reference genome for the match we found
        for SNP_index in potential_SNP:
            if mismatches.get((potential_match[SNP_index], seq_read[SNP_index], ref_start_index + SNP_index), "DNE") == "DNE":
                mismatches[(potential_match[SNP_index], seq_read[SNP_index], ref_start_index + SNP_index)] = 1
            else: 
                mismatches[(potential_match[SNP_index], seq_read[SNP_index], ref_start_index + SNP_index)] += 1
        match_found = ref_start_index
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

    # Remake reference dictionary with substrings that are approximately 1/3 the size of the reads
    # This is about ~16.666 since 50/3 = ~16.666; we'll round down and use 16
    # We'll match each read by each of its thirds to see if we get any partial matches i.e. parts of the string match but are offset
    substr_ref_index = dict()
    final_index = len(reference) - 16
    for index in range(0, final_index + 1):
        substr_ref_index[reference[index:index+16]] = index

    # TARGET DICTIONARIES

    # Preliminary SNPs
    mismatches = dict()
    # Preliminary insertions
    inserts = dict()
    # Preliminary deletes
    deletes = dict()

    # Array to store potential indels (reads that we did not find a perfect match for or only one SNP)
    oriented_potential_indels = []

    # HANDLE SNPS

    # Create dictionary of all mutations (key) and the number of times they are observed (value)
    for pair in input_reads:
        # Case 1: check if left end is in reference dictionary
        first_pos = test_read(pair[0], ref_index)
        # If found a match for left end in Case 1
        if first_pos != -1:
            # If the right end of the read also matches, this read is a perfect match, continue
            if test_read(pair[1][::-1], ref_index) != -1:
                continue
            # Find the mutations in the right end of the read
            else:
                right_match_found = find_match(reference, pair[1][::-1], substr_ref_index, mismatches)
                if right_match_found != -1:
                    continue
                else:
                    oriented_potential_indels.append([pair[1][::-1], first_pos + 50])
                    continue
        else:
            left_match_found = find_match(reference, pair[0], substr_ref_index, mismatches)
            
            if left_match_found != -1:
                # If the right end of the read also matches, this read is a perfect match, continue
                if test_read(pair[1][::-1], ref_index) != -1:
                    continue
                # Find the mutations in the right end of the read
                else:
                    right_match_found = find_match(reference, pair[1][::-1], substr_ref_index, mismatches)
                    if right_match_found != -1:
                        continue
                    else:
                        oriented_potential_indels.append([pair[1][::-1], left_match_found + 50])
                        continue
            # If no close matches for found this read at all, it could either be garbage or Indels, so we try other combinations of
            # else:

        # Case 2: check if reverse of left end is in reference dictionary
        first_rev_pos = test_read(pair[0][::-1], ref_index)
        # If found a match for left end in Case 2
        if first_rev_pos != -1:
            # If the right end of the read also matches, this read is a perfect match, continue
            if test_read(pair[1], ref_index) != -1:
                continue
            # Find the mutations in the right end of the read
            else:
                right_match_found = find_match(reference, pair[1], substr_ref_index, mismatches)
                if right_match_found != -1:
                    continue
                else:
                    oriented_potential_indels.append([pair[1], first_rev_pos + 50])
                    continue
        else:
            left_match_found = find_match(reference, pair[0][::-1], substr_ref_index, mismatches)
            
            if left_match_found != -1:
                # If the right end of the read also matches, this read is a perfect match, continue
                if test_read(pair[1], ref_index) != -1:
                    continue
                # Find the mutations in the right end of the read
                else:
                    right_match_found = find_match(reference, pair[1], substr_ref_index, mismatches)
                    if right_match_found != -1:
                        continue
                    else:
                        oriented_potential_indels.append([pair[1], left_match_found + 50])
                        continue
            # If no close matches for found this read at all, it could either be garbage or Indels
            # else:

        # Case 3: check if right end is in reference dictionary
        second_pos = test_read(pair[1], ref_index)
        # If found a match for left end in Case 3
        if second_pos != -1:
            # If the left end of the read also matches, this read is a perfect match, continue
            if test_read(pair[0][::-1], ref_index) != -1:
                continue
            # Find the mutations in the leftt end of the read 
            else:
                left_match_found = find_match(reference, pair[0][::-1], substr_ref_index, mismatches)
                if left_match_found != -1:
                    continue
                else:
                    oriented_potential_indels.append([pair[0][::-1], second_pos + 50])
                    continue
        else:
            right_match_found = find_match(reference, pair[1], substr_ref_index, mismatches)
            
            if right_match_found != -1:
                # If the left end of the read also matches, this read is a perfect match, continue
                if test_read(pair[0][::-1], ref_index) != -1:
                    continue
                # Find the mutations in the left end of the read
                else:
                    left_match_found = find_match(reference, pair[0][::-1], substr_ref_index, mismatches)
                    if left_match_found != -1:
                        continue
                    else:
                        oriented_potential_indels.append([pair[0][::-1], right_match_found + 50])
                        continue
            # If no close matches for found this read at all, it could either be garbage or Indels
            # else:

        # Case 4: check if reverse of right end is in reference dictionary
        second_rev_pos = test_read(pair[1][::-1], ref_index)
        # If found a match for left end in Case 4
        if second_rev_pos != -1:
            # If the right end of the read also matches, this read is a perfect match, continue
            if test_read(pair[0], ref_index) != -1:
                continue
            # Find the mutations in the right end of the read
            else:
                left_match_found = find_match(reference, pair[0], substr_ref_index, mismatches)
                if left_match_found != -1:
                    continue
                else:
                    oriented_potential_indels.append([pair[0], second_rev_pos + 50])
                    continue
        else:
            right_match_found = find_match(reference, pair[1][::-1], substr_ref_index, mismatches)
            
            if right_match_found != -1:
                # If the right end of the read also matches, this read is a perfect match, continue
                if test_read(pair[0], ref_index) != -1:
                    continue
                # Find the mutations in the right end of the read
                else:
                    left_match_found = find_match(reference, pair[0], substr_ref_index, mismatches)
                    if left_match_found != -1:
                        continue
                    else:
                        oriented_potential_indels.append([pair[0], right_match_found + 50])
                        continue
            # If no close matches for found this read at all, it could either be garbage or Indels
            # else:

        # If we reach the end of the loop (not a single continue was hit), this read has absolutely no matches; treat as garbage

    # HANDLE INDELS

    # We'll test each of these strings that we are confident are valid reads for the paired-end reads (because the other end of this read has a match)
    # However, we'll only test for indels if 2 out of 3 of the thirds of the read match; we'll do read[0:16], read[16:32], read[32:48]
    SNP_iterator = []
    insert_iterator = []
    delete_iterator = []
    for read in oriented_potential_indels:
        match_count = 0
        # Get the indices of anywhere we found matches for the thirds of potential indel strings
        # We do this because indels will be strings where parts of them will align perfectly, but it's offset due to insertion or deletion
        indices = [-1, -1, -1]
        if substr_ref_index.get(read[0][0:16], "DNE") != "DNE":
            indices[0] = substr_ref_index[read[0][0:16]]
            match_count += 1
        if substr_ref_index.get(read[0][16:32], "DNE") != "DNE":
            indices[1] = substr_ref_index[read[0][16:32]]
            match_count += 1
        if substr_ref_index.get(read[0][32:48], "DNE") != "DNE":
            indices[2] = substr_ref_index[read[0][32:48]]
            match_count += 1

        # If we found more than 1 match, we found that there's to parts of read
        if match_count > 1:
            # Calculate the starting and ending indices based on the matching thirds to get the entire portion of the reference that contains the indel
            # See find_match() function for more details about a similar process to this
            ref_start_index = 0
            read_start_index = 0
            if indices[0] != -1:
                ref_start_index = indices[0]
                read_start_index = 0
            elif indices[1] != -1:
                ref_start_index = indices[1]
                read_start_index = 16
            
            ref_end_index = 0
            read_end_index = 0
            if indices[2] != 1:
                ref_end_index = indices[2] + 16
                read_end_index = 48
            elif indices[1] != 1:
                ref_end_index = indices[1] + 16
                read_end_index = 32

            # We limit the reference string size that we suspect to contain the indel to a length of 60, effectively limiting the size of our potential indel
            if ref_end_index - ref_start_index > 60:
                continue

            # Run needleman-wunsch on this potential indel-containing read; not we may not pass in the whole read, but just the part of the read (the thirds)
            # where the indel may lie
            temp_score, SNP_iter, insert_iter, delete_iter = needleman_wunsch(reference[ref_start_index:ref_end_index], read[0][read_start_index:read_end_index], ref_start_index)

            # Append the indels to their respective lists so that we can iterate through and count them later
            # We found SNPs earlier so we don't use the SNPs returned (would probably require tuning in order for the SNPs to be viable to increase my sensitivity)
            # SNP_iterator = SNP_iterator + SNP_iter
            insert_iterator = insert_iterator + insert_iter
            delete_iterator = delete_iterator + delete_iter

    # GET FINAL TARGETS

    # for SNP in SNP_iterator:
    #     if mismatches.get((SNP[0], SNP[1], SNP[2]), "DNE") == "DNE":
    #         mismatches[(SNP[0], SNP[1], SNP[2])] = 1
    #     else: 
    #         mismatches[(SNP[0], SNP[1], SNP[2])] += 1

    # Tally the indels we found
    for insertion in insert_iterator:
        if inserts.get((insertion[0], insertion[1]), "DNE") == "DNE":
            inserts[(insertion[0], insertion[1])] = 1
        else: 
            inserts[(insertion[0], insertion[1])] += 1

    for deletion in delete_iterator:
        if deletes.get((deletion[0], deletion[1]), "DNE") == "DNE":
            deletes[(deletion[0], deletion[1])] = 1
        else: 
            deletes[(deletion[0], deletion[1])] += 1

    # Write dictionaries to pickle files so that I can load them in a testing program in order to tune my voting rules on what passes as an SNP or indel
    # This way, we avoid having to run the program again and again
    # with open('obj/'+ 'mismatches' + '.pkl', 'wb') as f:
    #     pickle.dump(mismatches, f, pickle.HIGHEST_PROTOCOL)

    # with open('obj/'+ 'inserts' + '.pkl', 'wb') as f:
    #     pickle.dump(inserts, f, pickle.HIGHEST_PROTOCOL)

    # with open('obj/'+ 'deletes' + '.pkl', 'wb') as f:
    #     pickle.dump(deletes, f, pickle.HIGHEST_PROTOCOL)

    # Iterate through the dictionary of mutations and pick only the ones that are observed above a certain threshold
    true_snps = []
    for k, v in mismatches.items():
        # Adjust what v > ? for changing confidence level of SNP (v is the number of occurrences we observed this SNP)
        if v > 2:
            true_snps.append([k[0], k[1], k[2]])

    true_insertions = []
    for k, v in inserts.items():
        # Adjust what v > ? for changing confidence level of SNP (v is the number of occurrences we observed this insertion)
        if v > 1:
            true_insertions.append([k[0], k[1]])

    true_deletions = []
    for k, v in deletes.items():
        # Adjust what v > ? for changing confidence level of SNP (v is the number of occurrences we observed this deletion)
        if v > 1:
            true_deletions.append([k[0], k[1]])

    # Sort the final SNPs/indels by index so that we can compare the printout with the answer key easily in a diffchecker
    true_snps = sorted(true_snps, key=lambda x: x[2])
    true_insertions = sorted(true_insertions, key=lambda x: x[1])
    true_deletions = sorted(true_deletions, key=lambda x: x[1])
    # print(true_snps)
    # print(true_insertions)
    # print(true_deletions)


    # Set the final lists for output to the lists we found
    snps = true_snps
    insertions = true_insertions
    deletions = true_deletions

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
