#! /usr/bin/python3

import argparse
import sys
from itertools import groupby

usage = 'This is a toolbox for handling fastq or fasta files'

parser = argparse.ArgumentParser(description=usage)

parser.add_argument('-v', '--version', action = 'version', version = '%(prog)s 1.0')
parser.add_argument('-i', dest = 'infile', metavar = 'INFILE', type = argparse.FileType('r'), required=True, help = 'define an input file')
parser.add_argument('-o', dest = 'outfile', metavar = 'OUTFILE', type = argparse.FileType('w'), default = sys.stdout, help = 'define an output file, if not defined print to stdout')

parser.add_argument('-convert', dest = 'convert', action = 'store_true', help = 'convert fastq to fasta format',)
parser.add_argument('-gc', dest = 'gc', action = 'store_true', help = 'calculate overal GC percentage')
parser.add_argument('-gc_each', dest = 'gc_each', action = 'store_true', help = 'calculate GC content in each sequence, print sequence ID and GC percentage')
parser.add_argument('-gc_positions', dest = 'gc_positions', action = 'store_true', help = 'calculate GC content in each position and output together with the sequence ID')
parser.add_argument('-reverse_comp', dest = 'reverse_comp', action = 'store_true', help = 'return reverse complements in fastq format')
parser.add_argument('-average_len', dest = 'average_len', action = 'store_true', help = 'return average sequence length of the whole file')
parser.add_argument('-min_max', dest = 'min_max', action = 'store_true', help = 'return the length of shortest and longest sequences')
parser.add_argument('-translate', dest = 'translate', action = 'store_true', help = 'translate dna to protein, output in fasta format')
parser.add_argument('-aa_frequency', dest = 'aa_frequency', action = 'store_true', help = 'calculate')
parser.add_argument('-codon_frequency', dest = 'codon_frequency', action = 'store_true', help = 'calculate codon frequency in fasta file')

args = parser.parse_args()

'''Function for calculating GC content in a sequence'''
def get_gc(seq):
    gc = 0
    at = 0
    for letter in seq:
        if letter in 'GC':
            gc += 1
        elif letter in 'AT':
            at += 1
    return round(100 * gc / len(seq), 2)

'''Function for reverse complement'''
def reverse(sequence):
    return sequence[::-1].translate(str.maketrans('ACGT', 'TGCA'))

'''Function for translating dna to rna'''
def transcribe(sequence):
    return sequence.replace('T', 'U')
'''Function for translating codons to amino acids'''
def translate_rna(s):
    protein = {"TTT" : "F", "CTT" : "L", "ATT" : "I", "GTT" : "V",
               "TTC" : "F", "CTC" : "L", "ATC" : "I", "GTC" : "V",
               "TTA" : "L", "CTA" : "L", "ATA" : "I", "GTA" : "V",
               "TTG" : "L", "CTG" : "L", "ATG" : "M", "GTG" : "V",
               "TCT" : "S", "CCT" : "P", "ACT" : "T", "GCT" : "A",
               "TCC" : "S", "CCC" : "P", "ACC" : "T", "GCC" : "A",
               "TCA" : "S", "CCA" : "P", "ACA" : "T", "GCA" : "A",
               "TCG" : "S", "CCG" : "P", "ACG" : "T", "GCG" : "A",
               "TAT" : "Y", "CAT" : "H", "AAT" : "N", "GAT" : "D",
               "TAC" : "Y", "CAC" : "H", "AAC" : "N", "GAC" : "D",
               "TAA" : "STOP", "CAA" : "Q", "AAA" : "K", "GAA" : "E",
               "TAG" : "STOP", "CAG" : "Q", "AAG" : "K", "GAG" : "E",
               "TGT" : "C", "CGT" : "R", "AGT" : "S", "GGT" : "G",
               "TGC" : "C", "CGC" : "R", "AGC" : "S", "GGC" : "G",
               "TGA" : "STOP", "CGA" : "R", "AGA" : "R", "GGA" : "G",
               "TGG" : "W", "CGG" : "R", "AGG" : "R", "GGG" : "G"
               }
    l = [protein.get(s[n:n+3], 'X') for n in range(0, len(s), 3)]
    return "".join(l)

'''Convert fastq file to fasta'''
if args.convert:
    sequence_list = []
    id_list = []
    for line_number, line in enumerate(args.infile):
        if line_number % 4 == 1: #sequence line
            sequence = line.rstrip().upper()
            sequence_list.append(sequence)
        if line_number % 4 == 0: #id line
            id_line = line[1:].rstrip()
            id_list.append(id_line)

    for id, sequence in zip(id_list, sequence_list):
        print('>{}\n{}'.format(id, sequence), file = args.outfile)

'''Count the overall GC percentage in a fastq file'''
if args.gc:
    total_gc_count = 0
    total_length = 0
    for line_number, line in enumerate(args.infile):
        if line_number % 4 == 1:
            sequence = line.rstrip().upper()
            total_gc_count += sequence.count('C') + sequence.count('G')
            total_length += len(sequence)
    gc_percentage = round((100 * total_gc_count/total_length), 2)
    print ('Total GC content: {}%'.format(gc_percentage), file = args.outfile)

''' Calculate each sequence GC content and output the ID with GC content'''
if args.gc_each:
    sequence_list = []
    id_list = []
    gc_content_list = []
    for line_number, line in enumerate(args.infile):
        if line_number % 4 == 1:
            sequence = line.rstrip()
            sequence = line.upper()
            sequence_list.append(sequence)
        if line_number % 4 == 0:
            id = line.rstrip()
            id_list.append(id)
    for sequence in sequence_list:
        gc_content = get_gc(sequence)
        gc_content_list.append(gc_content)
    for id, gc_content in zip(id_list, gc_content_list):
        print ('{}: {}%'.format(id, gc_content), file = args.outfile)

'''Calculate GC1, GC2, GC3 and output together with the ID'''
if args.gc_positions:
    id_list = []
    gc1_list = []
    gc2_list = []
    gc3_list = []
    for line_number, line in enumerate(args.infile):
        if line_number % 4 == 1:
            sequence = line.rstrip()
            first_position = line[0::3]
            gc1 = get_gc(first_position)
            gc1_list.append(gc1)
            second_position = line [1::3]
            gc2 = get_gc(second_position)
            gc2_list.append(gc2)
            third_position = line [2::3]
            gc3 = get_gc(third_position)
            gc3_list.append(gc3)
        if line_number % 4 == 0:
            id = line.rstrip()
            id_list.append(id)
    for id, gc1, gc2, gc3 in zip(id_list, gc1_list, gc2_list, gc3_list):
        print('{}: GC1 {}%, GC2 {}%, GC3 {}%'.format(id, gc1, gc2, gc3), file = args.outfile)

'''Return reverse complements in a fastq format'''
if args.reverse_comp:
    for line_number, line in enumerate(args.infile):
        if line_number % 4 == 0:
            id = line
        elif line_number % 4 == 1:
            sequence = line.rstrip().upper()
            reversed_sequence = reverse(sequence)
        elif line_number % 4 == 3:
            quality = line
            print ('{}{}+\n{}'.format(id, reversed_sequence, quality), end = '', file = args.outfile)

'''Calculate average sequence length'''
if args.average_len:
    total_sequences = 0
    total_nucleotides = 0
    for line_number, line in enumerate(args.infile):
        if line_number % 4 == 1:
            sequence_line = line
            total_sequences += 1
            total_nucleotides += len(sequence_line)
    average = round((total_nucleotides/total_sequences),2)
    print('Average sequence length is {} nucleotides'.format(average), file = args.outfile)


'''Calculate minimum and maximum sequence lengths'''
if args.min_max:
    sequence_list = []
    for line_number, line in enumerate(args.infile):
        if line_number % 4 == 1:
            sequence = line.rstrip().upper()
            sequence_list.append(sequence)
    minimum_length = min(sequence_list, key = len)
    maximum_length = max(sequence_list, key = len)
    print('The longest sequence is {} nucleotides long'.format(len(maximum_length)))
    print('The shortest sequence is {} nucleotides long'.format(len(minimum_length)), file = args.outfile)


'''Translate dna sequence to protein'''
if args.translate:
    faiter = (x[1] for x in groupby(args.infile, lambda l: l[0] == ">"))
    for header in faiter:
        header = next(header)[1:].strip()
        seq = "".join(s.strip() for s in next(faiter))
        seq1 = translate_rna(transcribe(seq))
        print('>{}\n{}'.format(header,seq1), file=args.outfile)

'''Calculate amino acid frequencies'''
if args.aa_frequency:
    dictionary = {}
    amino_acids = 'ACDEFGHIKLMNPQRSTVWYX'
    unknown_letters = ['B', 'J', 'O', 'U', 'Z']
    for letter in amino_acids:
        dictionary[letter] = 0

    for line in args.infile:
            if line[0] != '>':
                dictionary[key] += line.count(key)
    for letter1 in unknown_letters:
            if line[0] != '>':
                dictionary['X'] += line.count(letter1)

    for element in dictionary:
        print ('{} {}'.format(element,dictionary[element]), file = args.outfile)

'''Calculate codon frequency'''
if args.codon_frequency:
    dictionary = {}
    codons = ['UUU', 'UUC', 'UUA', 'UUG', 'UCU', 'UCC', 'UCA', 'UCG', 'UAU', 'UAC', 'UAA', 'UAG', 'UGU', 'UGC', 'UGA', 'UGG', 'CUU', 'CUC', 'CUA', 'CUG', 'CCU', 'CCC', 'CCA', 'CCG', 'CAU', 'CAC', 'CAA', 'CAG', 'CGU', 'CGC', 'CGA', 'CGG', 'AUU', 'AUC', 'AUA', 'AUG', 'ACU', 'ACC', 'ACA', 'ACG', 'AAU', 'AAC', 'AAA', 'AAG', 'AGU', 'AGC', 'AGA', 'AGG', 'GUU', 'GUC', 'GUA', 'GUG', 'GCU', 'GCC', 'GCA', 'GCG', 'GAU', 'GAC', 'GAA', 'GAG', 'GGU', 'GGC', 'GGA', 'GGG']
    for codon in codons:
        dictionary[codon] = 0

    for line in args.infile:
        for key in dictionary:
            if line[0] != '>':
                line = transcribe(line)
                split_sequence = [line[i:i+3] for i in range(0, len(line), 3)]
                dictionary[key] += split_sequence.count(key)

    for element in dictionary:
        print ('{} {}'.format(element,dictionary[element]))
