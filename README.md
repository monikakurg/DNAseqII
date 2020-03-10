# MKTool
> A toolbox designed to handle fasta or fastq files.

## Description

A lot of bioinformatitians and people who have just started their work in this field handle genomic information on daily basis. Some projects require handling some commands repeatedly on different fasta or fastq files. This toolbox handles some basic and some more advanced commands that can help ease the work. The input files are required to be in fasta or fastq formats. The program prints to standard output, unless defined otherwise. This toolbox contains options such as GC content calculations, convertion, length calculations, translation, amino acid and codon frequrncies calculation.

![](header.png)

## Version

v1.0

## Installation

No installation is required. If the program is being used for the first time, make it executable:

```sh
chmod +x MKTool.py
```

## Usage example

Convert fastq file format in to fasta format

```sh
./MKTool.py -i input.fastq -o output.fasta -convert
```
Translate DNA sequence to protein sequence

```sh
./MKTool.py -i input.fasta -o output.fasta -translate
```
## Available options

```sh
-convert          convert fastq format to fasta
-gc               calculate overal GC percentage (fastq)
-gc_each          calculate GC content in each sequence, print sequence ID and GC percentage (fastq)
-gc_positions     calculate GC content in each position and output together with the sequence ID (fastq)
-reverse_comp     return reverse complements in fastq format (fastq)
-average_len      return average sequence length of the whole file (fastq)
-min_max          return the length of shortest and longest sequences (fastq)
-translate        translate dna to protein, output in fasta format (fasta)
-aa_frequency     calculate the abundance of amino acids (fasta)
-codon_frequency  calculate codon frequency (fasta)
```

## Limitations

Some options do not recognise both fasta and fastq format. In the help section of the program, there are explanations which format is suitable for a desired option.
When using the '-translate' option, the sequences have to be coding sequences with a start codon

## Creator and contact

Monika Kurgonaite – monika.kurg@gmail.com
