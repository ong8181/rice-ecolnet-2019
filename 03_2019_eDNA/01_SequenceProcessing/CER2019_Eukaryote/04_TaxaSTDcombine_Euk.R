####
#### CER rice 2019 eDNA analysis
#### No.4 STD sequences assignment
####

# Create output directory
dir.create("04_TaxaSTDcombine_EukOut")

# Load Claident taxa assignment
claident_tax <- read.delim("03_ClaidentAssignTax_EukOut/EukASV_merge_classigntax")
head(claident_tax)

# Load Blastn taxa assignment for STD sequences
# Select sequences with < 2 mismatches and length > 120
blastn_std0 <- read.table("02_ident_EukSTD_BLASTnOut/EukSTD_out.txt")
cond1 <- blastn_std0$V5 < 2 & blastn_std0$V4 > 120 & blastn_std0$V6 < 1 # Alignment length and No. of mismatch and gap
cond2 <- blastn_std0$V7 < 3 & blastn_std0$V8 > 120 # Start and end points of query
blastn_std <- blastn_std0[cond1 & cond2,]

# Check claident taxa assignments of the potential std sequences
potential_std_id <- match(blastn_std$V1, claident_tax$query)
claident_tax[potential_std_id, "family"] # Not close to any existing tax --> OK
claident_tax[potential_std_id, "species"] # Not close to any existing tax --> OK

# Replace STD taxa names with claident taxa assigment
claident_tax$species <- as.character(claident_tax$species)
claident_tax[potential_std_id, "species"] <- as.character(blastn_std$V2)

# Output new claident tax table
write.csv(claident_tax, "04_TaxaSTDcombine_EukOut/claident_tax_revise.csv", row.names = F)

#### save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/04_SessionInfo_EukSTDident_%s.txt", substr(Sys.time(), 1, 10)))
