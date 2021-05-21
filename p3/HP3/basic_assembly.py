from os.path import join
import sys
import time
from collections import defaultdict, Counter
import sys
import os
import zipfile
import argparse
sys.path.insert(0, os.path.abspath(".."))
sys.path.insert(0, os.path.abspath("../.."))


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


"""
    TODO: Use this space to implement any additional functions you might need

"""
# To keep track where forking branches are
def get_kmers(text_arr, k, kmer_counts, coverage):
    size = len(text_arr[0]) - k + 1
    for text in text_arr:
        for index in range(0, size):
            if kmer_counts.get(text[index:index+k], "DNE") != "DNE":
                kmer_counts[text[index:index+k]] += ((kmer_counts[text[index:index+k]] * coverage) + 1) / coverage
            else:
                kmer_counts[text[index:index+k]] = 1 / coverage
    return

def filter_kmers(kmer_counts, error_threshold):
    to_delete = []
    for k, v in kmer_counts.items():
        if round(v) < error_threshold:
            to_delete.append(k)
        else:
            kmer_counts[k] = round(v)
    for key in to_delete:
        del kmer_counts[key]

def get_graph_dict(kmer_dict, k):
    graph_dict = dict()
    node_degrees = dict()
    size = k - 1
    for kmer, count in kmer_dict.items():
        left_end = kmer[:size]
        right_end = kmer[1:]
        if left_end in graph_dict:
            graph_dict[left_end].append(right_end)
        else:
            graph_dict[left_end] = [right_end]
            
        if left_end in node_degrees:
            node_degrees[left_end][1] += 1
        else:
            node_degrees[left_end] = [0, 1]
        if right_end in node_degrees:
            node_degrees[right_end][0] += 1
        else:
            node_degrees[right_end] = [1, 0]
    return graph_dict, node_degrees

def get_graph_nodes(node_degrees):
    nodes = []
    for key in node_degrees.keys():
        nodes.append(key)
    return nodes

def get_cycle(graph_dict, node, start_node, curr_path):
    if node in graph_dict:
        # Given that this is an isolated cycle where nodes are 1-in-1-out
        curr_path += node + ' -> '
        next_node = graph_dict[node][0]
        del graph_dict[node]
        return get_cycle(graph_dict, next_node, start_node, curr_path)
    else:
        # We deleted the start node to this cycle so if the current_node isn't here, then it must be the start
        if node == start_node:
            curr_path += node
            return [curr_path]
        else:
            return []

def get_maximal_non_branching(graph_dict, node_degrees, nodes):
    path_arrays = []
    for node in nodes:
        degrees = node_degrees[node]
        # Has to not be a 1 in 1 out
        if degrees[0] != 1 or degrees[1] != 1:
            if degrees[1] > 0:
                for next_node in graph_dict[node]:
                    path = [node, next_node]
                    current_node = next_node
                    # Has to be a 1 in 1 out
                    while node_degrees[current_node][0] == 1 and node_degrees[current_node][1] == 1:
                        path.append(graph_dict[current_node][0])
                        current_node = graph_dict[current_node][0]
                    path_arrays.append(path)
                    
    for path in path_arrays:
        for node in path:
            if node in graph_dict:
                del graph_dict[node]
                
    cyclic_nodes = []
    cycles = []
    for key in graph_dict.keys():
        cyclic_nodes.append(key)
    for node in cyclic_nodes:
        if node in graph_dict:
            cycle = get_cycle(graph_dict, node, node, "")
            cycles += cycle
                
    non_branching_paths = []
    for path in path_arrays:
        current_path = ""
        for index, node in enumerate(path):
            if index == 0:
                current_path += node
            else:
                current_path += node[-1]
        non_branching_paths.append(current_path)
    return non_branching_paths + cycles


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='basic_assembly.py takes in data for homework assignment 3 consisting '
                                                 'of a set of reads and aligns the reads to the reference genome.')
    parser.add_argument('-r', '--reads', required=True, dest='reads_file',
                        help='File containg sequencing reads.')
    parser.add_argument('-o', '--outputFile', required=True, dest='output_file',
                        help='Output file name.')
    parser.add_argument('-t', '--outputHeader', required=True, dest='output_header',
                        help='String that needs to be outputted on the first line of the output file so that the\n'
                             'online submission system recognizes which leaderboard this file should be submitted to.\n'
                             'This HAS to be one of the following:\n'
                             '1) spectrum_A_1_chr_1 for 10K spectrum reads practice data\n'
                             '2) practice_A_2_chr_1 for 10k normal reads practice data\n'
                             '3) hw3all_A_3_chr_1 for project 3 for-credit data\n')
    args = parser.parse_args()
    reads_fn = args.reads_file

    input_reads = parse_reads_file(reads_fn)
    if input_reads is None:
        sys.exit(1)

    """
            TODO: Call functions to do the actual assembly here

    """
    # VARIABLES
    # Tune these to get better score
    k = 25
    coverage = 10
    error_threshold = 2

    # Unwrap lists inside input_reads, we'll just assume left end is forward and reverse the right end of the read
    reads = []
    for pair in input_reads:
        reads.append(pair[0])
        reads.append(pair[1][::-1])
    
    # GET THE KMERS
    # First get the kmers from the reads
    kmer_counts = dict()
    get_kmers(reads, k, kmer_counts, coverage)
    # Remove any kmers that don't meet our error threshold
    filter_kmers(kmer_counts, error_threshold)

    # GET DE BRUIJN GRAPH
    graph_dict, node_degrees = get_graph_dict(kmer_counts, k)
    nodes = get_graph_nodes(node_degrees)

    # GET CONTIGS
    all_contigs = get_maximal_non_branching(graph_dict, node_degrees, nodes)


    contigs = all_contigs

    output_fn = args.output_file
    zip_fn = output_fn + '.zip'
    with open(output_fn, 'w') as output_file:
        output_file.write('>' + args.output_header + '\n')
        output_file.write('>ASSEMBLY\n')
        output_file.write('\n'.join(contigs))
    with zipfile.ZipFile(zip_fn, 'w') as myzip:
        myzip.write(output_fn)
