# CLOVERS
<br/>

![LOGO](LOGO.png)

<br/>

![C++](https://img.shields.io/badge/language-C%2B%2B-blue)
![Bioinformatics](https://img.shields.io/badge/field-bioinformatics-brightgreen)
![Genomics](https://img.shields.io/badge/domain-genomics-purple)
![AVX](https://img.shields.io/badge/optimized-AVX%2FFMA-red)
![Static](https://img.shields.io/badge/build-static-lightgrey)
![Dependencies](https://img.shields.io/badge/dependencies-minimal-green)
[![GPLv3 License](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![x86_64](https://img.shields.io/badge/arch-x86__64-green)

Accurate *ab initio* prediction of prokaryotic and phagic genes and with CLOVERS !

## Contents
- **[Overview](#overview)** - Project introduction and features
- **[Setup](#setup)** - Installation guide for different operating systems
- **[Usage](#usage)** - Command-line interface and options
- **[Examples](#examples)** - Usage examples with command-line options
- **[Input](#input)** - Description of the input files and formats
- **[Output](#output)** - Description of the output files and formats
- **[Performance](#performance)** - Performance benchmarks and comparisons with other tools
- **[Structure](#structure)** - Project file structure and organization
- **[Server](#server)** - Free available server for running tasks online
- **[Contact](#contact)** - Contact information for questions or support
- **[License](#license)** - GNU General Public License v3.0

## Overview
CLOVERS is an novel ab initio gene finder that utilizes the Z-curve method for robust protein-coding gene prediction on prokaryotic and phage genomes. Its novel technical approach gives it an advantage in detecting more overlapped genes - which are often missed by other state-of-the-art tools - while remaining highly competitive in overall gene recognition.

## Setup

Download the latest precompiled binary file from the [release page](https://tubic.tju.edu.cn/clovers/download), and decompress it to the directory of your choice.

## Usage
### Quick Start
We recommend configure the environment variable `PATH` to include the directory of the executable binary file, such that you can run `clovers` directly in the terminal.

```bash
clovers -i example.fasta -o example.gff -c -f gff
```
### Options

```bash
clovers.exe [OPTION...]
```

#### General Options
* `-h, --help`  
Print help menu and exit. If no option is specified (only `./clovers` was entered or double-clicked the icon), print the help menu and wait for the user to press the Enter key or Ctrl+C to exit.

* `-q, --quiet`  
Run quietly with no stderr output (note that fatal error messages will still be printed under quiet mode).

* `-v, --version`  
Print version info and exit (useful when to quickly check whether it has been configured in the environment).

* `-T, --threads`  
Number of threads to use. Please set as a positive integer. (default: half)

#### Input/Output Options
* `-i, --input`   
Specify FASTA/Genbank/EMBL input file or their compressed versions (gzip). (default: stdin)  

* `-o, --output`  
Specify output file or '-' as standard output (stdout).

* `-f, --format`  
Select output format (gff, gbk, med, gbk-full). (default: gff)

* `-a, --faa`  
Write protein translations of genes to the selected file or '-' as stdout.

* `-d, --fna`  
Write nucleotide sequences of genes to the selected file or '-' as stdout.

#### CLOVERS Options  
* `-g, --table`  
Specify a genetic codon table to use (1, 4, 11, 15, 16, 25 or auto). (default: 11)  

  **Note:** Auto-selection of codon table function can only determine whether TGA or TAG have been reprogrammed into sense codons.

* `-l, --minlen`  
Specify the mininum length (nt) of ORFs. (default: 90)

* `-c, --circ`  
Treat the default topology as circular. Note that the topology for each sequence can be set by the words "circular" or "linear" appear in the header line (FASTA/GenBank/EMBL), and this option will only take effect when the words "circular" or "linear" do not exist.

* `-p, --proc`  
Select prediction procedure (single or meta).

* `-s, --thres`  
Specify putative gene probability score threshold. For chromosomes, the recommended setting is 0.5. For plasmids, viruses, and bacteriophages, the recommended setting is 0.4. (default: 0.5)

* `-t, --train`  
Write (if none exists) or use the specified training file.

#### TriTISA+ Options  
* `-n, --bypass`  
Bypass TriTISA+ and output the longest ORFs (most left 5'-end).

* `-M,--maxiter`  
Max iteration times for RBS revision. (default: 20)

* `-R,--rbs`  
Write (if none exists) or use the specified RBS model file.

#### GOP-Reporter options:  
* `-L, --min-olen`  
Specify the mininum overlapped length between two ORFs. (default: 120)

* `-O, --overlap`  
Write overlapped genes to the selected file (GFF3 format).

* `-A, --amino`  
Write protein translations of overlapped genes to the selected file.

* `-D, --nucl`  
Write nucleotide sequences of overlapped genes to the selected file.

## Examples
```bash
# Predict genes in a FASTA file and output to GFF3
clovers -i example.fasta -o example.gff

# Use auto-detection of translation table
clovers -i input.fasta -o output.gff -g auto

# Generate protein and nucleotide sequences
clovers -i input.fasta -o output.gff -a proteins.faa -d genes.fna

# Detect overlapped genes and their sequences
clovers -i input.fasta -o output.gff -O overlap.gff -A overlap.faa -D overlap.fna

# Predict genes in metagenomes
clovers -i input.fasta -o output.gff -p meta

# Generate training file and RBS model file
clovers -i input.fasta -o output.gff -t training.dat -I rbs_model.bin
```

## Input
CLOVERS supports FASTA format input (stdin or file) or their compressed versions (gzip) with content like follow:
```
>NC_000913.3 circular
AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC
TTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAA
TATAGGCATAGCGCACAGACAGATAAAAATTACAGAGTACACAACATCCATGAAACGCATTAGCACCACC
. . .
```
Or minimal GenBank format input (stdout or file) or their compressed versions (gzip) with content like follow:
```
LOCUS       NC_000913.3
DEFINITION  NC_000913.3
FEATURES             Location/Qualifiers
ORIGIN   
        1 agcttttcat tctgactgca acgggcaata tgtctctgtg tggattaaaa aaagagtgtc
       61 tgatagcagc ttctgaactg gttacctgcc gtgagtaaat taaaatttta ttgacttagg
      121 tcactaaata ctttaaccaa tataggcata gcgcacagac agataaaaat tacagagtac
          . . .
//
```
Or minimal EMBL format input (stdout or file) or their compressed versions (gzip) with content like follow:
```
ID   NC_000913.3
XX
FH   Key             Location/Qualifiers
XX
SQ   
     agcttttcat tctgactgca acgggcaata tgtctctgtg tggattaaaa aaagagtgtc        60
     tgatagcagc ttctgaactg gttacctgcc gtgagtaaat taaaatttta ttgacttagg       120
     tcactaaata ctttaaccaa tataggcata gcgcacagac agataaaaat tacagagtac       180
     . . .
//
```
## Output
- **GFF/GenBank/MED files**: Primary annotation results containing gene locations, scores and other attributes

    **GFF File Example:**
    ```
    ##gff-version 3
    # trans_tbl: 11
    # NC_000913.3	4641652 bp	circular	UNA	19-MAY-2026
    NC_000913.3	CLOVERS_v1.1.0	CDS	337	2799	0.974	+	0	ID=ORF000001
    NC_000913.3	CLOVERS_v1.1.0	CDS	2801	3733	0.968	+	0	ID=ORF000002
    . . .
    ```
    **GenBank File Example:**
    ```
    LOCUS       NC_000913.3          4641652 bp    DNA     circular UNA 19-MAY-2026
    DEFINITION  NC_000913.3
    FEATURES             Location/Qualifiers
         CDS             337..2799
                         /locus_tag=ORF000001
                         /transl_table=11
                         /codon_start=1
                         /note="Derived by protein-coding gene prediction method: CLOVERS_v1.1.0"
                         /translation="MRVLKFGG...LGV*"
         CDS             2801..3733
                         /locus_tag=ORF000002
                         /transl_table=11
                         /codon_start=1
                         /note="Derived by protein-coding gene prediction method: CLOVERS_v1.1.0"
                         /translation="MVKVYAPA...LEN*"
    . . .
    ORIGIN
    //
    ```
    **MED File Example:**
    ```
    ## MED
    # trans_tbl: 11
    # NC_000913.3	4641652 bp	circular	UNA	19-MAY-2026
    337 2799	+
    2801 3733	+
    . . .
    ```
- **Protein sequence files (.faa)**: Translated amino acid sequences of predicted genes
    ```
    >NC_000913.3:2801..3733(+) score=0.968
    MVKVYAPASSANMSVGFDVLGAAVTPVDGALLGDVVTVEAAETFSLNNLGRFADKLPSEPRENIVYQCWERFCQELGKQI
    PVAMTLEKNMPIGSGLGSSACSVVAALMAMNEHCGKPLNDTRLLALMGELEGRISGSIHYDNVAPCFLGGMQLMIEENDI
    ISQQVPGFDEWLWVLAYPGIKVSTAEARAILPAQYRRQDCIAHGRHLAGFIHACYSRQPELAAKLMKDVIAEPYRERLLP
    GFRQARQAVAEIGAVASGISGSGPTLFALCDKPETAQRVADWLGKNYLQNQEGFVHICRLDTAGARVLEN*
    . . .
    ```
- **Nucleotide sequence files (.fna)**: DNA sequences of predicted genes
    ```
    >NC_000913.3:2801..3733(+) score=0.968
    ATGGTTAAAGTTTATGCCCCGGCTTCCAGTGCCAATATGAGCGTCGGGTTTGATGTGCTCGGGGCGGCGGTGACACCTGT
    TGATGGTGCATTGCTCGGAGATGTAGTCACGGTTGAGGCGGCAGAGACATTCAGTCTCAACAACCTCGGACGCTTTGCCG
    ATAAGCTGCCGTCAGAACCACGGGAAAATATCGTTTATCAGTGCTGGGAGCGTTTTTGCCAGGAACTGGGTAAGCAAATT
    CCAGTGGCGATGACCCTGGAAAAGAATATGCCGATCGGTTCGGGCTTAGGCTCCAGTGCCTGTTCGGTGGTCGCGGCGCT
    GATGGCGATGAATGAACACTGCGGCAAGCCGCTTAATGACACTCGTTTGCTGGCTTTGATGGGCGAGCTGGAAGGCCGTA
    TCTCCGGCAGCATTCATTACGACAACGTGGCACCGTGTTTTCTCGGTGGTATGCAGTTGATGATCGAAGAAAACGACATC
    ATCAGCCAGCAAGTGCCAGGGTTTGATGAGTGGCTGTGGGTGCTGGCGTATCCGGGGATTAAAGTCTCGACGGCAGAAGC
    CAGGGCTATTTTACCGGCGCAGTATCGCCGCCAGGATTGCATTGCGCACGGGCGACATCTGGCAGGCTTCATTCACGCCT
    GCTATTCCCGTCAGCCTGAGCTTGCCGCGAAGCTGATGAAAGATGTTATCGCTGAACCCTACCGTGAACGGTTACTGCCA
    GGCTTCCGGCAGGCGCGGCAGGCGGTCGCGGAAATCGGCGCGGTAGCGAGCGGTATCTCCGGCTCCGGCCCGACCTTGTT
    CGCTCTGTGTGACAAGCCGGAAACCGCCCAGCGCGTTGCCGACTGGTTGGGTAAGAACTACCTGCAAAATCAGGAAGGTT
    TTGTTCATATTTGCCGGCTGGATACGGCGGGCGCACGAGTACTGGAAAACTAA
    . . .

## Server
- **CLOVERS Server**: Free available at [https://tubic.tju.edu.cn/clovers/](https://tubic.tju.edu.cn/clovers/).

## Contact
- **Authors**: Zetong Zhang, Yan Lin*, Feng Gao*
- **Emails**: zhangzetong@tju.edu.cn | ylin@tju.edu.cn | fgao@tju.edu.cn

## License
CLOVERS is distributed under the GNU General Public License v3.0. See the [LICENSE](LICENSE) file for details.

