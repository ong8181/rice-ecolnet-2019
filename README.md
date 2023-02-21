# Analysis codes for Ushio et al. (2023) "Detecting and validating influential organisms for rice growth: An ecological network approach"

This repository contains analysis codes to reproduce the results in Ushio et al. (2023) _bioRxiv_ https://doi.org/10.1101/2023.02.19.529115.

# License
See LICENSE.


# Analysis workflow
## Step 1. Detection of causal relationships between ecological community members and rice growth in 2017
- `01_2017_DNAxRice/`: Codes for Unified Information-theoretic Causality (UIC) analysis of quantitative environmental DNA data and rice growth in 2017.<br>

## Step 2. Summary of rice growth data in 2019
- `02_2019_Rice/`: Codes for summarizing and visualizing rice growth in 2019<br>

## Step 3. Sequence data processing of environmental DNA data in 2019
- `03_2019_eDNA/`: Codes for sequence data processing of quantitative environmental DNA data in 2019<br>

## Step 4. Detection of differentially expressed genes in rice in 2019
- `04_2019_riceRNA/`: Codes for RNA sequence data processing and detection of differentially expressed genes in rice leaves in 2019<br>

## Step 5. Format figures
- `FigCode/`: Format figures<br>
- `FigCode/Fig_Combine.Rmd`: Combine figures<br>


# Package versions
see text files in `00_SessionInfo/`


# Data availablity
Sequence data are deposited in DDBJ Sequence Read Archives (DRA). The accession numbers are as follows: DRA009658 (16S), DRA009659 (ITS), DRA009660 (COI), and DRA009661 (18S) for eDNA data of ecological communities in 2017 (Ushio 2022 _Proceedings of the Royal Society B_ https://doi.org/10.1098/rspb.2021.2690), DRA015682 (16S), DRA015683 (ITS), DRA015685 (COI), and DRA015686 (18S) for eDNA data of ecological communities in 2019, and DRA015706 (Rice RNA-seq) for rice RNA expression data.


## Get sequence data
You may download sequence data (`*.fastq.bz2`) from DDBJ Sequence Read Archives (DRA) by executing the following commands. `DRAXXXXXX` should be DRA accession number. Alternatively, you can access https:/ddbj.nig.ac.jp/DRASearch/submission?acc=DRAXXXXXX.

```
# Prepare folders
cd ~/Desktop
mkdir temp
mkdir temp/xml
cd temp

# Specify DRA accession number
DRA_ACCESSION="DRAXXXXXX"

# Download data (downlaod only *.bz2 and *.xml files)
wget -r --no-parent -A "*.bz2","*.xml" "ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/DRA015/${DRA_ACCESSION}/"

# If you need to specify the proxy:
#wget -r --no-parent -A "*.bz2","*.xml" -e FTP_PROXY=proxy.xxx.xx:xxxx "ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/DRA015/${DRA_ACCESSION}/"

# Move files
mv ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/DRA015/${DRA_ACCESSION}/*/*.fastq.bz2 ./
mv ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/DRA015/${DRA_ACCESSION}/*.xml xml

# Delete temporal folder
rm -r ftp.ddbj.nig.ac.jp

# Convert bz2 to gz format (if you wish)
#for f in *.bz2; do
#  bzcat "$f" | gzip -c >"${f%.bz2}.gz"
#  rm $f
#done
```


