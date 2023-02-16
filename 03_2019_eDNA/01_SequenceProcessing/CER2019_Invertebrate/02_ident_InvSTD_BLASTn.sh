##
## Standard STD identification
##

cd "01_SequenceProcessing/CER2019_Invertebrate"

DBPATH=/home/ushio/Desktop/DADA2_DB/STDseqs/InvertebrateSTD/InvertebrateSTD_mlCOI
QUERYPATH=01_CERrice2019_InvSeqOut/InvASV_seqs.fa
OUTPUT=02_ident_InvSTD_BLASTnOut/InvSTD_out.txt
EVALUE_SET=1e-167

mkdir 02_ident_InvSTD_BLASTnOut
blastn -db ${DBPATH} -query ${QUERYPATH} -evalue ${EVALUE_SET} -outfmt 6 -out ${OUTPUT}

