# ERAP tiler

This repository is designed to tile together overlapping ERAP fragments, to provide an assembled CDS/GEN sequence when it is unambiguously possible.

This is not the neatest repository, and probably deserves a cleanup/refactor/implementation into TGSTS once it is stable.

## Requirements

As a jupyter notebook, this must be installed in the python version.
Also, the [TGS Toolset](https://github.com/Anthony-Nolan/TgsToolset) must be installed.
Minimap2 and Bowtie are also requiements.

## Input

An excel sheet as per [250218_assembled.xlsx](https://github.com/richardnatarajanAN/erap_tiler/blob/main/250218_assembled.xlsx), generally created as a dump from FileMakerPro

## Usage
- Open the [tiler.ipynb](https://github.com/richardnatarajanAN/erap_tiler/blob/main/tiler.ipynb).
- In cells 2 and 3, set the input file and output dir.
- Run all cells.


## Methodology

### Data loading
All sequences with the following criteria are removed:
- Analysis code == 0
- Null sequence

If analysis code == 2, this is deemed homozygous, and the fragment is duplicated.

All sample ids are cleaned (replacing '/'s with '-'s)

Following this, each sample is described by pairs of fragments, each containing sequence as well as metadata including fragment name, read count, etc.

### Alignment
All fragments are aligned to the global ERAP reference using minimap2, for speed reasons. Pipes are inserted post alignment, relative to the ERAP reference.

### Valid overlap identification
Between each fragment and it's neighbours, there are 4 possible overlaps. The validity of these overlaps are assessed as per:
- CDS match
- GEN match
- RR-Tolerant GEN match (identical across all sequence ignoring extruncs of repeat regions)
- K-mer distance
An acceptable overlap must be a CDS match and a RR-Tolerant GEN match, in other words, intronic repeat region length mismatches are excusable.

This is performed across all neighbouring fragments, allowing a 'tree' of acceptable paths to be created.

If no possible overlaps are found, assembly fails and alignment information is written to the /fails file.

### Assembly
Based on acceptable overlaps, all possible sequences (ignoring UTR) are constructed. If there are only two possiblilities, the sequence is deemed unambiguously assembleable and taken forward.
Two identical sequences is perfectly acceptable in the homozygous case.
We ignore the UTRs simply because there are many cases where UTR sequencing error caused GEN ambiguity, and in this body of work the phased UTR is not important information. 
We do not trim the UTR as it provides padding needed to identify the exons with our bowtie based exonic sequence identification.

If >2 sequences are possible, this is a GEN ambiguous case. This occurs when there are SNPs outside the overlap which cannot be phased based on the overlaps. In this case, we attempt a CDS assembly.

In a CDS assembly, we verify that there are only 2 possible CDS sequences. If so, this means all assembly ambiguity is intronic. This means we cannot accurately assess the introns, but the CDS
is unambiguous, which is important information within the current body of work. These cases are still written to the 'fails', where extra information describing the nature of the intronic ambiguity
is stored.

If >2 CDS sequences are possible, this is a truly ambiguous case, and alignment information is written to the 'fails'

Once this has been completed, the GEN sequence for each allele (and full GEN sequence if possible) is determined by merging overlaps. The GEN sequence for GEN ambig cases is calculated
for the purposes of typing, but is not written to the results.

### Typing
From this, using the GEN sequence, we calculate the CDS, Exonic and protein sequences. Exonic sequence is calculated using ANTs, CDS and PROT is calculated using SFAT.

SFAT uses the rules that the intronic position 34124 determines splicing (A T at this position results in 19 exons, with a small chunk of intron 19 being included. An A results in the longer case where exon 20 is included).

We then calculate the CDS/Exonic/GEN mismatches, and write that alongside other metadata.
