####
#### Command to demultiplex bcl to fastq
#### (For CMR-011 [600 cycle]; CEReDNA2019 samples)
####

#----------- Demultiplex using claident -----------#
# Set parameters
RUNNAME="CMR-011"
RUN_FOLDER="200320_M00962_0082_000000000-CT7T6"
FASTQ_OUT_FOLDER2="CMR-011_fastq_out"
DEMULTIPLEX_OUT="CMR-011_ClaidentDemultiplexed"
F_PRIMER="../TagPrimerFiles/CMR-011_F_primer.txt"
R_PRIMER="../TagPrimerFiles/CMR-011_R_primer.txt"
I7_INDEX="../TagPrimerFiles/CMR-011_i7_index.txt"
I5_INDEX="../TagPrimerFiles/CMR-011_i5_index.txt"

# Convert Bcl to Fastq (bcl2fastq2 v2.18)
# !!! Rename "SampleSheet.csv" before this command (e.g., "SampleSheet_rename.csv") !!!
mv ${RUN_FOLDER}/SampleSheet.csv ${RUN_FOLDER}/SampleSheet_rename.csv
mv ${RUN_FOLDER}/Data/Intensities/BaseCalls/SampleSheet.csv ${RUN_FOLDER}/Data/Intensities/BaseCalls/SampleSheet_rename.csv
bcl2fastq --processing-threads 72 --use-bases-mask Y300n,I8,I8,Y300n --create-fastq-for-index-reads --runfolder-dir $RUN_FOLDER --output-dir $FASTQ_OUT_FOLDER2

# Demultiplexing
cd $FASTQ_OUT_FOLDER2
clsplitseq --runname=$RUNNAME --index1file=$I7_INDEX --index2file=$I5_INDEX --primerfile=$F_PRIMER --reverseprimerfile=$R_PRIMER --minqualtag=30 --numthreads=72 --truncateN=enable *_R1_001.fastq.gz *_I1_001.fastq.gz *_I2_001.fastq.gz *_R2_001.fastq.gz $DEMULTIPLEX_OUT

# Move undetermined fastq to another folder
cd $DEMULTIPLEX_OUT
mkdir Undetermined
mv *undetermined.*.fastq.gz Undetermined

# Move demultiplex folder
cd ..
mv $DEMULTIPLEX_OUT ..

