##
## Standard STD identification
##

cd "01_SequenceProcessing/CER2019_Fungi"

DBPATH=/home/ushio/Desktop/DADA2_DB/STDseqs/FungiSTD/FungiSTD_ITSKYO
QUERYPATH=01_CERrice2019_FungiSeqOut/FungiASV_seqs.fa
OUTPUT=02_ident_FungiSTD_BLASTnOut/FungiSTD_out.txt
EVALUE_SET=1e-118

mkdir 02_ident_FungiSTD_BLASTnOut
blastn -db ${DBPATH} -query ${QUERYPATH} -evalue ${EVALUE_SET} -outfmt 6 -out ${OUTPUT}


