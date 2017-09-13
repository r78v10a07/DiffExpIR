DiffExpIR
===

This program execute a differential expression analysis for introns using as
input the output provided by the TPMCalculator.

## Output file

The output file includes these columns:

1. GeneId
2. Chr
3. Start (Gene start position)
4. End (Gene end position)
5. Intron_Start
6. Intron_End
7. Log2TPMRatio (Log2 TPM_sample1/TPM_sample_2)
8. TPM_1 (TPM group of samples 1)
9. TPM_2 (TPM group of samples 2)
10. minusLog10PValue (-1.0 * log10(P-Value))
11. PValue (P-Value)
12. RValue_1 (TPM_Intron_i/(TPM_Exon_i-1 + TPM_Exon_i+1) for group of samples 1)
13. RValue_2 (TPM_Intron_i/(TPM_Exon_i-1 + TPM_Exon_i+1) for group of samples 2)

### BAMTools

Clone the BAMTools repository from GitHub: https://github.com/pezmaster31/bamtools

Compile it on this way and set the environment variables for DiffExpIR:

    cd bamtools
    mkdir build
    cd build
    cmake ..
    make
    cd ..
    export BAMTOOLS_DIR=`pwd`
    export CPPFLAGS="-I $BAMTOOLS_DIR/include"
    export LDFLAGS="-L $BAMTOOLS_DIR/lib -Wl,-rpath,$BAMTOOLS_DIR/lib"

That's it. BAMTools was compiled and the env variables were set for compiling
DiffExpIR.

## Installation

After the instalation of BAMTools go to the DiffExpIR folder and do make:

    make

A bin folder will be created with the DiffExpIR executable.

## Usage

Usage: ./bin/diffexpIR -g GTF_file -o output_file_name -d TPMCalculator_output_dir -p sample_1,sample_2

    ./bin/diffexpIR options:

    -h    Display this usage information.
    -g    GTF file
    -o    Output file name
    -d    Directory with the TPM output files
    -k    Gene key to use from GTF file. Default: gene_id
    -t    Transcript key to use from GTF file. Default: transcript_id
    -c    Smaller size allowed for an intron created for genes. Default: 16. We recommend to use the reads length
    -p    Prefix for grouping samples. (sample_1,sample_2)
    -s    Stat method: ttest (default), wilcox
    -f    Minimum fold change to filter out (default value: 2.0), wilcox
    -v    Minimum P-Value to filter out (default value: 1.0E-6)
    -r    Minimum fold change between intron and neighboring exons (default value: -1.0)

## Credits

Roberto Vera Alvarez, PhD

Emails: veraalva@ncbi.nlm.nih.gov, r78v10a07@gmail.com

# Public Domain notice

## National Center for Biotechnology Information.

This software is a "United States Government Work" under the terms of the United States
Copyright Act. It was written as part of the authors' official duties as United States
Government employees and thus cannot be copyrighted. This software is freely available
to the public for use. The National Library of Medicine and the U.S. Government have not
placed any restriction on its use or reproduction.

Although all reasonable efforts have been taken to ensure the accuracy and reliability
of the software and data, the NLM and the U.S. Government do not and cannot warrant the
performance or results that may be obtained by using this software or data. The NLM and
the U.S. Government disclaim all warranties, express or implied, including warranties
of performance, merchantability or fitness for any particular purpose.

Please cite NCBI in any work or product based on this material.
