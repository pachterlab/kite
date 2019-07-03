# _kite_: kallisto indexing and tag extraction

This package enables fast and accurate pre-processing of Feature Barcoding experiments, a common datatype in single-cell genomics. In Feature Barcoding assays, cellular data are recorded as short DNA sequences using procedures adapted from single-cell RNA-seq. 

The __kite ("kallisto indexing and tag extraction__") package prepares input files prior to running the [kallisto | bustools](https://www.kallistobus.tools/getting_started.html) scRNA-seq pipeline. Starting with a .csv file containing Feature Barcode names and Feature Barcode sequences, the program `featuremap.py` generates a "mismatch map" and outputs "mismatch" fasta and "mismatch" transcript-to-gene (t2g) files. The mismatch files, containing the Feature Barcode sequences and their Hamming distance = 1 mismatches, are used to run kallisto | bustools on Feature Barcoding data. 

The mismatch fasta file is used by `kallisto index` with a k-mer length `-k` set to the length of the Feature Barcode. 

The mismatch t2g file is used by `bustools count` to generate a Features x Cells matrix. 

In this way, kallisto | bustools will effectively search the sequencing data for the Feature Barcodes and their Hamming distance = 1 neighbors. We find that for Feature Barcodes of moderate length (6-15bp) pre-processing is remarkably fast and the results equivalent to or better than those from traditional alignment.

A walk-through from the kallisto | bustools [Tutorials](https://www.kallistobus.tools/tutorials) page is reproduced below, and a copmlete Feature Barcode analysis can be found in the [docs](https://github.com/pachterlab/kite/tree/master/docs/) directory of the `kite` GitHub repository.

## kite Installation
Clone the GitHub repo to obtain the core featuremap.py program and some useful accessory files. 
```
$ git clone https://github.com/pachterlab/kite
```
## System Requirements
Feature Barcode pre-processing requires up-to-date versions of `kallisto` and `bustools`
```
Python3
kallisto v0.46 or higher
bustools v0.39.0 or higher
```
For downstream analysis, we use [ScanPy](https://scanpy.readthedocs.io/en/stable/installation.html) (Wolf et. al, Genome Biology 2018) and the [LeidenAlg](https://github.com/vtraag/leidenalg) clustering package (Traag et. al, arXiv 2018).

## kite pre-processing

#### `featuremap.py FeatureBarcodes.csv`
The featuremap.py program is run prior to the standard kallisto | bustools pipeline. It takes a .csv input and outputs "mismatch" transcript-to-gene (t2g) and fasta files that can be used by kallisto | bustools to complete pre-processing (see below and Vignettes). The program takes a single argument, FeatureBarcodes.csv, and outputs mismatch fasta and t2g files to the working directory.

FeatureBarcodes.csv: path to a .csv-formatted file containing Feature Barcode names and sequences (example below).

returns mismatch t2g and fasta files saved to the working directory

### NOTE: Use only odd values for k-mer length during `kallisto index` 
To avoid potential pseudoalignment errors arising from inverted repeats, kallisto only accepts odd values for the k-mer length `-k`. If your Feature Barcodes have an even length, just add an appropriate constant base on one side and follow the protocol as suggested. For example, append an __A__ base to the CD3_TotalSeqB barcode AACAAGACCCTTGAG → AACAAGACCCTTGAGA. Adding constant bases in this way increases specificity and may be useful for experiments with low sequencing quality or very short Feature Barcodes. 

## Brief Example: 1k PBMCs from a Healthy Donor - Gene Expression and Cell Surface Protein

The [docs](https://github.com/pachterlab/kite/tree/master/docs) folder contains a complete analysis ([10x_kiteVignette.ipynb](https://github.com/pachterlab/kite/tree/master/docs/10X_kiteVignette.ipynb)) for a 10x dataset collected on ~730 peripheral blood mononuclear cells (PBMCs) labeled with 17 unique Feature Barcoded antibodies. The dataset can be found [here](https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_1k_protein_v3). 

### 1. Download materials
Prepare a folder:
```
$ mkdir kallisto_bustools_kite/
$ cd kallisto_bustools_kite/
```
Download and unqip the following files.
```
$ wget http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_protein_v3/pbmc_1k_protein_v3_fastqs.tar
$ tar -xvf ./pbmc_1k_protein_v3_fastqs.tar
$ wget http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_protein_v3/pbmc_1k_protein_v3_feature_ref.csv
$ wget http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.tar.gz
$ tar -xvzf ./pbmc_1k_protein_v3_filtered_feature_bc_matrix.tar.gz
$ wget https://github.com/BUStools/getting_started/releases/download/species_mixing/10xv3_whitelist.txt
```
#### 2. Make the mismatch FASTA and t2g files
Start by preparing a csv-formatted matrix of Feature Barcode names and Feaure Barcode sequences, __including a header__, is used as input. Do not include any common or constant sequences. In this case, we parsed the feature_ref.csv file provided by 10x to give a properly formatted csv (below). Example code for this step and a correctly formatted file (FeatureBarcodes.csv) is included in the [kite GitHub repo](https://github.com/pachterlab/kite/docs/).

|Feature Barcode name|Feature Barcode sequence|
| ------------- | ------------- |
|CD3_TotalSeqB|AACAAGACCCTTGAG|
|CD8a_TotalSeqB|TACCCGTAATAGCGT|
|CD14_TotalSeqB|GAAAGTCAAAGCACT|
|CD15_TotalSeqB|ACGAATCAATCTGTG|
|CD16_TotalSeqB|GTCTTTGTCAGTGCA|
|CD56_TotalSeqB|GTTGTCCGACAATAC|
|CD19_TotalSeqB|TCAACGCTTGGCTAG|
|CD25_TotalSeqB|GTGCATTCAACAGTA|
|CD45RA_TotalSeqB|GATGAGAACAGGTTT|
|CD45RO_TotalSeqB|TGCATGTCATCGGTG|
|PD-1_TotalSeqB|AAGTCGTGAGGCATG|
|TIGIT_TotalSeqB|TGAAGGCTCATTTGT|
|CD127_TotalSeqB|ACATTGACGCAACTA|
|IgG2a_control_TotalSeqB|CTCTATTCAGACCAG|
|IgG1_control_TotalSeqB|ACTCACTGGAGTCTC|
|IgG2b_control_TotalSeqB| ATCACATCGTTGCCA|

Now run featuremap.py, which creates a mismatch FASTA file and a mismatch t2g file for the experiment. In this case the mismatch files each have 782 entries. 
```
$ ./kite/featuremap/featuremap.py FeatureBarcodes.csv
```
__Note:__ kallisto only accepts odd values for the k-mer length, so if your Feature Barcodes are even in length, add a constant base on either side before running featuremap.py. For example, append an __A__ base to the CD3_TotalSeqB barcode AACAAGACCCTTGAG &rarr; AACAAGACCCTTGAGA

Processing Feature Barcodes is similar to processing transcripts except instead of looking for transcript fragments of length `-k` (the k-mer length) in the reads, a "mismatch" index is used to search the raw reads for the Feature Barcode whitelist and mismatch sequences. Please refer to the [kallisto documentation](https://www.kallistobus.tools/documentation) for more information on the kallisto | bustools workflow. 

Because Feature Barcodes are typically designed to be robust to some sequencing errors, each Feature Barcode and its mismatches are unique across an experiment, thus each Feature Barcode equivalence class has a one-to-one correspondence to a member of the Feature Barcode whitelist. This is reflected in the t2g file, where each mismatch Feature Barcode points to a unique parent Feature Barcode from the whitelist, analogous to the relationship between genes and transcripts in the case of cDNA processing. 

```
$head -n 4 ./10xFeatures_t2g.txt
CD3_TotalSeqB	CD3_TotalSeqB	CD3_TotalSeqB
CD3_TotalSeqB-0-1	CD3_TotalSeqB	CD3_TotalSeqB
CD3_TotalSeqB-0-2	CD3_TotalSeqB	CD3_TotalSeqB
CD3_TotalSeqB-0-3	CD3_TotalSeqB	CD3_TotalSeqB

$head -n 8 ./10xFeaturesMismatch.fa
>CD3_TotalSeqB
AACAAGACCCTTGAG
>CD3_TotalSeqB-0-1
TACAAGACCCTTGAG
>CD3_TotalSeqB-0-2
GACAAGACCCTTGAG
>CD3_TotalSeqB-0-3
CACAAGACCCTTGAG
```

#### 3. Build Index
Build the kallisto index using the mismatch fasta and a k-mer length `-k` equal to the length of the Feature Barcodes:
```
$ kallisto index -i FeaturesMismatch.idx -k 15 ./FeaturesMismatch.fa
```
#### 4. Run kallisto
Pseudoalign the reads:
```
$ kallisto bus -i FeaturesMismatch.idx -o ./ -x 10xv3 -t 4 \
./pbmc_1k_protein_v3_fastqs/pbmc_1k_protein_v3_antibody_fastqs/pbmc_1k_protein_v3_antibody_S2_L001_R1_001.fastq.gz \
./pbmc_1k_protein_v3_fastqs/pbmc_1k_protein_v3_antibody_fastqs/pbmc_1k_protein_v3_antibody_S2_L001_R2_001.fastq.gz \
./pbmc_1k_protein_v3_fastqs/pbmc_1k_protein_v3_antibody_fastqs/pbmc_1k_protein_v3_antibody_S2_L002_R1_001.fastq.gz \
./pbmc_1k_protein_v3_fastqs/pbmc_1k_protein_v3_antibody_fastqs/pbmc_1k_protein_v3_antibody_S2_L002_R2_001.fastq.gz \
```

#### 5. Run bustools
For `bustools count`, use the mismatch t2g file. 
```
$ bustools correct -w ./10xv3_whitelist.txt ./output.bus -o ./output_corrected.bus
```
```
$ bustools sort -t 4 -o ./output_sorted.bus ./output_corrected.bus
```
```
$ mkdir ./featurecounts/
```
```
$ bustools count -o ./featurecounts/featurecounts --genecounts -g ./Features_t2g.txt -e ./matrix.ec -t ./transcripts.txt ./output_sorted.bus
```
#### 6. Analyze count matrix
`Bustools count` outputs a .mtx-formatted Features x Cells matrix and vectors of gene names and cell barcodes (genes.txt and barcodes.txt). From here, standard analysis packages like ScanPy and Seurat can be used to continue the Feature Barcode analysis. For details, check out the [Jupyter notebook](https://github.com/pachterlab/kite/tree/master/docs/).
