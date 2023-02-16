####
#### No.4: Claident assigntax for EUK
#### (For Euk1391f sequences)
####

#----------- Taxa assignment using claident -----------#
# Taxa assignment using claident (after picking ASV by DADA2)
IN_DIR="01_SequenceProcessing/CER2019_Eukaryote/01_CERrice2019_EukSeqOut"
OUT_DIR="01_SequenceProcessing/CER2019_Eukaryote/03_ClaidentAssignTax_EukOut"
ASV_NAME="EukASV"

# Prepare
mkdir ${OUT_DIR}
cd ${IN_DIR}

# Check overall_class
clmakecachedb --blastdb=overall_class --numthreads=72 ${IN_DIR}/${ASV_NAME}_seqs.fa ${OUT_DIR}/${ASV_NAME}_overall_cache
clidentseq --blastdb=${OUT_DIR}/${ASV_NAME}_overall_cache --numthreads=72 ${IN_DIR}/${ASV_NAME}_seqs.fa ${OUT_DIR}/${ASV_NAME}_overall_clidentseq
classigntax --taxdb=overall_class --maxpopposer=0.05 --minsoratio=19 ${OUT_DIR}/${ASV_NAME}_overall_clidentseq ${OUT_DIR}/${ASV_NAME}_overall_classigntax

# Overall genus
clmakecachedb --blastdb=overall_genus --numthreads=72 ${IN_DIR}/${ASV_NAME}_seqs.fa ${OUT_DIR}/${ASV_NAME}_overallg_cache
clidentseq --blastdb=${OUT_DIR}/${ASV_NAME}_overallg_cache --numthreads=72 ${IN_DIR}/${ASV_NAME}_seqs.fa ${OUT_DIR}/${ASV_NAME}_overallg_clidentseq
classigntax --taxdb=overall_genus --maxpopposer=0.05 --minsoratio=19 ${OUT_DIR}/${ASV_NAME}_overallg_clidentseq ${OUT_DIR}/${ASV_NAME}_overallg_classigntax

# Merge identification results (overall_class + overall_genus)
clmergeassign --priority=descend \
${OUT_DIR}/${ASV_NAME}_overallg_classigntax \
${OUT_DIR}/${ASV_NAME}_overall_classigntax \
${OUT_DIR}/${ASV_NAME}_merge_classigntax
