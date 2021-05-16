import sys
import argparse
import itertools
import numpy as np
import time
import zipfile

import pickle

from numpy.lib import insert

def load_obj(name):
    with open('obj/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)

if __name__ == "__main__":
    mismatches = load_obj('mismatches')
    inserts = load_obj('inserts')
    deletes = load_obj('deletes')

    print(mismatches)
    print(inserts)
    print(deletes)

    # Iterate through the dictionary of mutations and pick only the ones that are observed above a certain threshold
    true_snps = []
    for k, v in mismatches.items():
        # Adjust what v > ? for changing confidence level of SNP (v is the number of occurrences we observed this SNP)
        if v > 1:
            true_snps.append([k[0], k[1], k[2]])

    true_insertions = []
    for k, v in inserts.items():
        # Adjust what v > ? for changing confidence level of SNP (v is the number of occurrences we observed this SNP)
        if v > 1:
            true_insertions.append([k[0], k[1]])

    true_deletions = []
    for k, v in deletes.items():
        # Adjust what v > ? for changing confidence level of SNP (v is the number of occurrences we observed this SNP)
        if v > 1:
            true_deletions.append([k[0], k[1]])

    true_snps = sorted(true_snps, key=lambda x: x[2])
    true_insertions = sorted(true_insertions, key=lambda x: x[1])
    true_deletions = sorted(true_deletions, key=lambda x: x[1])
    print(true_snps)
    print(true_insertions)
    print(true_deletions)



    snps = true_snps
    insertions = true_insertions
    deletions = true_deletions

    output_fn = 'test_output.txt'
    zip_fn = output_fn + '.zip'
    with open(output_fn, 'w') as output_file:
        output_file.write('>' + 'practice_W_3_chr_1' + '\n>SNP\n')
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