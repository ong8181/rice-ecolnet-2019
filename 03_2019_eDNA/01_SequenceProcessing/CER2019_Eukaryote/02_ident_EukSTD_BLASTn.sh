##
## Standard STD identification
##

cd "01_SequenceProcessing/CER2019_Eukaryote"

DBPATH=/home/ushio/Desktop/DADA2_DB/STDseqs/EukaryoteSTD/Euk1391f_STD
QUERYPATH=01_CERrice2019_EukSeqOut/EukASV_seqs.fa
OUTPUT=02_ident_EukSTD_BLASTnOut/EukSTD_out.txt
EVALUE_SET=1e-50

mkdir 02_ident_EukSTD_BLASTnOut
blastn -db ${DBPATH} -query ${QUERYPATH} -evalue ${EVALUE_SET} -outfmt 6 -out ${OUTPUT}

