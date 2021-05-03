import sys
import argparse
import time
import zipfile


def parse_reads_file(reads_fn):
    """
    :param reads_fn: the file containing all of the reads
    :return: outputs a list of all paired-end reads
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
# Return position of matched 50 bp end in reference genome, otherwise return -1
def test_read(seq_read, ref_index):
    if ref_index.get(seq_read, "DNE") != "DNE":
        return ref_index[seq_read]
    else:
        return -1

# Take the unmatched end of the read, the chunk of reference we will match it to, the original position in the reference genome of the chunk, and the mismatch dictionary
def test_gap(seq_read, reference, original_pos, mismatches):
    for gap in range(0, 20):
        mismatch_count = 0
        mismatch_pos = 0
        for index in range(0, len(seq_read)):
            if seq_read[index] != reference[index + gap]:
                if mismatch_count == 0:
                    mismatch_pos = index
                mismatch_count += 1
            if mismatch_count >= 2:
                break
        if mismatch_count == 0:
            break
        if mismatch_count == 1:
            if mismatches.get((reference[mismatch_pos + gap], seq_read[mismatch_pos], original_pos + mismatch_pos + gap), "DNE") == "DNE":
                mismatches[(reference[mismatch_pos + gap], seq_read[mismatch_pos], original_pos + mismatch_pos + gap)] = 1
            else: 
                mismatches[(reference[mismatch_pos + gap], seq_read[mismatch_pos], original_pos + mismatch_pos + gap)] += 1
            break

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
    parser = argparse.ArgumentParser(description='basic_aligner.py takes in data for homework assignment 1 consisting '
                                     'of a genome and a set of reads and aligns the reads to the reference genome, '
                                     'then calls SNPs based on this alignment')
    parser.add_argument('-g', '--referenceGenome', required=True, dest='reference_file',
                        help='File containing a reference genome.')
    parser.add_argument('-r', '--reads', required=True, dest='reads_file',
                        help='File containg sequencing reads.')
    parser.add_argument('-o', '--outputFile', required=True, dest='output_file',
                        help='Output file name.')
    parser.add_argument('-t', '--outputHeader', required=True, dest='output_header',
                        help='String that needs to be outputted on the first line of the output file so that the '
                             'online submission system recognizes which leaderboard this file should be submitted to.'
                             'This HAS to be practice_W_1_chr_1 for the practice data and hw1_W_2_chr_1 for the '
                             'for-credit assignment!')
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
    # flatten the list of paired-end reads, so we have list of just reads
    # reads = [item for sublist in input_reads for item in sublist]

    # Create dictionary of all substrings of length 50 (key) and their positions (value)
    ref_index = dict()
    final_index = len(reference) - len(input_reads[0][0])
    for index in range(0, final_index + 1):
        ref_index[reference[index:index+50]] = index

    # Create dictionary of all mutations (key) and the number of times they are observed (value)
    mismatches = dict()
    for pair in input_reads:
        # Case 1: check if left end is in reference dictionary
        first_pos = test_read(pair[0], ref_index)
        # Case 2: check if reverse of left end is in reference dictionary
        first_rev_pos = test_read(pair[0][::-1], ref_index)

        # If Case 1
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

        # If first_normal, second_reversed didn't find a match
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

    snps = true_snps

    output_fn = args.output_file
    zip_fn = output_fn + '.zip'
    with open(output_fn, 'w') as output_file:
        header = '>' + args.output_header + '\n>SNP\n'
        output_file.write(header)
        for x in snps:
            line = ','.join([str(u) for u in x]) + '\n'
            output_file.write(line)

        tails = ('>' + x for x in ('STR', 'CNV', 'ALU', 'INV', 'INS', 'DEL'))
        output_file.write('\n'.join(tails))

    with zipfile.ZipFile(zip_fn, 'w') as myzip:
        myzip.write(output_fn)


"""
    # brute force
    mismatches = []
    for read in reads:
        rev_read = read[::-1]
        index = 0
        for index in range(0, final_index + 1):
            ref_sequence = reference[index:index+50]
            mismatch_count = 0
            rev_mismatch_count = 0
            mismatch_pos = 0
            rev_mismatch_pos = 0
            for pos in range(0, len(read)):
                if read[pos] != ref_sequence[pos]:
                    if mismatch_count == 0:
                        mismatch_pos = pos
                    mismatch_count += 1
                if rev_read[pos] != ref_sequence[pos]:
                    if mismatch_count == 0:
                        rev_mismatch_pos = pos
                    rev_mismatch_count += 1
                if mismatch_count >= 2 and rev_mismatch_count >= 2:
                    break
            if mismatch_count == 1:
                mismatches.append([ref_sequence[mismatch_pos], read[mismatch_pos], mismatch_pos + index])
            elif rev_mismatch_count == 1:
                mismatches.append([ref_sequence[rev_mismatch_pos], read[rev_mismatch_pos], rev_mismatch_pos + index])
    print(mismatches)
    """