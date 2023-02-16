##
## Standard STD identification
##

cd "01_SequenceProcessing/CER2019_Prokaryote"

DBPATH=/home/ushio/Desktop/DADA2_DB/STDseqs/ProkaryoteSTD/ProkSTD_515F
QUERYPATH=01_CERrice2019_ProkSeqOut/ProkASV_seqs.fa
OUTPUT=02_ident_ProkSTD_BLASTnOut/ProkSTD_out.txt
EVALUE_SET=1e-100

mkdir 02_ident_ProkSTD_BLASTnOut
blastn -db ${DBPATH} -query ${QUERYPATH} -evalue ${EVALUE_SET} -outfmt 6 -out ${OUTPUT}

