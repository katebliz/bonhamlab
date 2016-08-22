#!/bin/bash
# QiimeRunOtus.sh
# Elizabeth Grice Lab
# University of Pennsylvania

# Requires Casey's scripts: Summarize_FASTA_Size.pl, Filter_FASTA_Size.pl
# Assumes R scripts in same folder as this script!
# Written using QIIME 1.8

########################################################################
#                              Set Up                                  #
########################################################################

# Load in the latest version of R
if [ -f /etc/profile.d/modules.sh ]; then
    source /etc/profile.d/modules.sh
fi

# Because denovo is a Boolean variable that will be changed by including the flag, first set the default variable
# Default will use the denovo OTU clustering process
DENOVO="FALSE"
RAREFAC=10100
NUMCORE=1
FACTORS=0 # Set default factors to none

# extract command line options and their arguments into variables based on this while loop taking in conditional position arguments.
usage() { echo -e "Usage: This is the standardized Grice Lab QIIME workflow, with the following options: 
    [-j] Project ID name (required)
    [-m] Mapping file that only conatins wanted samples (required)
    [-i] Input FNA file (required)
    [-w] Working directory (required)
    [-s] Server used: CONSIGN or SCISUB
    [-g] Greengenes reference directory if you want to override
         the default folders based on server (optional)
    [-d] Denovo clustering method. If you want reference based,
         don't include this flag. (optional)
    [-r] Rarefaction sequence upper limit (optional)
    [-c] Number of proc cores (optional)
    [-v] Run version (required) Acceptable values are V1V3, V4, or 
         {Minimum sequence length cutoff},{Maximum sequence length cutoff}
    [-f] List of comma separated factors (optional) to use for
         reads/sample and rarefaction analyses}"; exit 1;}

while getopts ":hj:m:i:w:s:g:d:r:c:v:f:" option; do
    case "$option" in
        h) usage ;;
        j) PROJ="$OPTARG";;
        m) MAP="$OPTARG";;
        i) INPUT_FNA="$OPTARG";;
        w) WORKDIR="$OPTARG";;
        s) SERVER="$OPTARG";;
        g) GGDIR="$OPTARG";;
        d) DENOVO="$OPTARG";; #Right now the denovo is only working for V4 seuquencing runs.
        r) RAREFAC="$OPTARG";;
        c) NUMCORE="$OPTARG";;
        v) RUNVERSION="$OPTARG";;
        f) FACTORS="$OPTARG";;
        :)  echo "Error: -$OPTARG requires an argument" ; exit 1;;
        ?)  echo "Error: unknown option -$OPTARG" ; exit 1;;
    esac
done

# Assign min and max sequence length for size filtering
case "$RUNVERSION" in
    "V1V3") RUNSTAT_MIN=490 ; RUNSTAT_MAX=600;;
    "V4") RUNSTAT_MIN=248 ; RUNSTAT_MAX=255;;
    *) RUNSTAT_MIN=$(echo $RUNVERSION | sed 's/\,.*//') ; RUNSTAT_MAX=$(echo $RUNVERSION | sed 's/.*\,//') ; shift 2 ;;
esac

# Define reference folder depending on server
case "$SERVER" in
    "CONSIGN") GGDIR="/project/egricelab/references/greengenes/" ; source /opt/software/qiime/1.8.0/activate.sh ; module load "R-3.2.2" ;;
    "SCISUB") GGDIR="/project/grice/pub_data/greengenes/" ; module load "R/3.2.2" ; module load "python/2.7.10" ;;
esac

if [[ -z $PROJ ]]; then
    echo "ERROR: ProjectID is required"; exit 1;
elif [[ -z $MAP ]]; then
    echo "ERROR: Map is required"; exit 1;
elif [[ -z $INPUT_FNA ]]; then
    echo "ERROR: Input FNA file is required"; exit 1;
elif [[ -z $WORKDIR ]]; then
    echo "ERROR: Working directory is required"; exit 1;
elif [[ -z $GGDIR ]]; then
    echo "ERROR: Greengenes reference directory is required"; exit 1;
elif [[ -z $RUNSTAT_MIN ]] || [[ -z $RUNSTAT_MAX ]]; then
    echo "ERROR: Run version is required"; exit 1;
fi

echo -e "QIIME Run Starting Variables:
--------------------------------------------------
# ProjectID              : ${PROJ}
# Mapping File           : ${MAP}
# Input FASTA            : ${INPUT_FNA}
# Working Directory      : ${WORKDIR}
# Server Used            : ${SERVER}
# Greengenes Directory   : ${GGDIR}
# Denovo method          : ${DENOVO}
# Rarefaction Upper Limit: ${RAREFAC}
# Cores                  : ${NUMCORE}
# Sequence Size Range    : ${RUNSTAT_MIN} - ${RUNSTAT_MAX}
# Factors for R Scripts  : ${FACTORS}
--------------------------------------------------\n" > qiime_run_setup.log

# Set the working directory
cd ${WORKDIR}
echo -e "\nChanging working directory to ${WORKDIR}...\n"

# Take combined FASTA file and filter out non-project reads with mapping file
echo -e 'Set up: Removing unwanted sample sequences\n'
cut -f1 ${MAP} | tail -n +2 > subject_id.txt
filter_fasta.py -f ${INPUT_FNA} -o ${PROJ}_specific_raw.fna --sample_id_fp subject_id.txt
rm subject_id.txt

########################################################################
#                             Processing                               #
########################################################################

# Summarize sequences based on length
echo -e 'Processing: summarizing read lengths\n' 
# Use Casey's summarizing script - in same folder
perl $(dirname $0)/Summarize_FASTA_Size.pl ${PROJ}_specific_raw.fna

# Perform read length analyses in R: summary file, min, max
Rscript $(dirname $0)/ReadLengthAnalysis.R ${PROJ}_specific_raw_SizeStats.txt RUNSTAT_MIN RUNSTAT_MAX

# Filter sequences based on length
echo -e 'Processing: filtering reads by size\n' 
# Use Casey's size filtering script - in same folder
perl $(dirname $0)/Filter_FASTA_Size.pl ${PROJ}_specific_raw.fna $RUNSTAT_MIN $RUNSTAT_MAX
mv ${PROJ}_specific_raw_kept.fna ${PROJ}_size_filtered.fna

# Cluster sequences into OTUs
echo -e 'Processing: clustering into OTUs\n'
if [[ $DENOVO == "FALSE" ]]; then
    parallel_pick_otus_uclust_ref.py -i ${PROJ}_size_filtered.fna -r ${GGDIR}/gg_13_8_otus/rep_set/97_otus.fasta -o otu_map/ -O $NUMCORE
elif [[ $DENOVO == "swarm" ]]; then
    pick_otus.py -i ${PROJ}_size_filtered.fna -t --otu_picking_method swarm --swarm_resolution 1 -o otu_map --thread $NUMCORE
else
    pick_otus.py -i ${PROJ}_size_filtered.fna -t --otu_picking_method $DENOVO -s 0.97 -o otu_map
fi

#Pick a representative sequence from each OTU
echo -e 'Processing: picking rep set\n'
pick_rep_set.py -f ${PROJ}_size_filtered.fna -i otu_map/${PROJ}_size_filtered_otus.txt --rep_set_picking_method most_abundant -o ${PROJ}_rep_set.fna
echo Number of OTUs: `expr $(wc -l ${PROJ}_rep_set.fna | cut -f 1 -d ' ') / 2`

# Assign taxonomy to each OTU representative
echo -e 'Processing: assign taxonomy\n'
parallel_assign_taxonomy_rdp.py -O $NUMCORE -i ${PROJ}_rep_set.fna -o assigned_taxonomy -c 0.8
# List Unclassified and Cyanobacteria sequences
egrep '(Unclassified|Cyanobacteria)' assigned_taxonomy/${PROJ}_rep_set_tax_assignments.txt | awk '{print$1}' > excluded_otus_taxa.txt

# NOTE: THIS SECTION IS LIKELY TO BE REMOVED DUE TO LOW UTILITY OF CHIMERA DETECTION. ALSO, QIIME SAYS IT SHOULD BE RUN
#       BEFORE YOU RUN THE OTU CHECKING STEP, NOT AFTER.
##Align sequences
#echo -e 'Processing: aligning\n'
#parallel_align_seqs_pynast.py -O $NUMCORE -i ${PROJ}_rep_set.fna -o alignment -e 150 -p 0.75 -t ${REFDIR}/greengenes_alignment/core_set_aligned.fasta.imputed
#
##Identify chimeras
#echo 'Processing: chimera checking'
#mkdir chimeras
#parallel_identify_chimeric_seqs.py -O $NUMCORE -m ChimeraSlayer -i alignment/${PROJ}_rep_set_aligned.fasta -a /project/egricelab/references/greengenes/core_set_aligned.fasta.imputed -o chimeras/chimeric_seqs.txt
#cut -f1 chimeras/chimeric_seqs.txt > excluded_otus_chimera.txt
#Make an OTU table with Chimeras, Unclassifieds and Cyanobacteria removed
#echo 'Processing: make OTU table w/o chimeras and unwanted taxa'
#cat excluded_otus_taxa.txt excluded_otus_chimera.txt | sort | uniq -u > excluded_otus_combined.txt
#make_otu_table.py -i otu_map/${PROJ}_size_filtered_otus.txt -t assigned_taxonomy/${PROJ}_rep_set_tax_assignments.txt -o ${PROJ}_otu_table.biom -e excluded_otus_combined.txt

#Make an OTU table with Chimeras, Unclassifieds and Cyanobacteria removed
echo -e 'Processing: make OTU table w/o unwanted taxa\n'
make_otu_table.py -i otu_map/${PROJ}_size_filtered_otus.txt -t assigned_taxonomy/${PROJ}_rep_set_tax_assignments.txt -o ${PROJ}_otu_table.biom -e excluded_otus_taxa.txt
#Print out stats for each sample in the otu table.
biom summarize-table -i ${PROJ}_otu_table.biom -o ${PROJ}_otu_table_summary.txt

#Remove rare taxa from OTU table
echo -e 'Processing: filter rare OTUs from table\n'
filter_otus_from_otu_table.py -i ${PROJ}_otu_table.biom -o ${PROJ}_otu_table_no_singletons.biom -n 2 -s 2
#Print out stats for each sample in the nonzero category
biom summarize-table -i ${PROJ}_otu_table_no_singletons.biom -o ${PROJ}_otu_table_no_singletons_summary.txt

#Make new rep set with the singletons and unwanted OTUs removed
filter_fasta.py -f ${PROJ}_rep_set.fna -o ${PROJ}_rep_set_no_singletons.fna -b ${PROJ}_otu_table_no_singletons.biom

# Perform a final alignment of the OTUs for phylogenetic purposes
parallel_align_seqs_pynast.py -O $NUMCORE -i ${PROJ}_rep_set_no_singletons.fna -o alignment -e 150 -p 0.75 -t ${GGDIR}/greengenes_alignment/core_set_aligned.fasta.imputed


#######################################################################
#                           Rarefaction                               #
#######################################################################

# Rarefy OTU table to find ideal cutoff for eliminating samples. These rarefied OTU tables
#   will be used to make alpha diversity rarefaction plots.
echo -e 'Rarefaction: Subsampling OTU table\n'
multiple_rarefactions.py -i ${PROJ}_otu_table_no_singletons.biom -o rarefactions -m 100 -x $RAREFAC -s 500 -n 5

# Calculate alpha diversity metrics on rarefied OTU tables. Afterwards, manually remove empty temp
#   folder left behind by QIIME script. This is an error on QIIME's part, but it will stop collate
#   _alphy.py from running.
echo -e 'Rarefaction: Calculating Alpha Diversity\n'
parallel_alpha_diversity.py -i rarefactions -o rarefactions/alpha_metrics -m observed_species,shannon -O $NUMCORE
# This will throw errors but it serves the correct purpose so do not change.
rmdir rarefactions/alpha_metrics/*

# Collate alpha diversity metrics and iterations into one file for plotting.
echo -e 'Rarefaction: Collating Metrics\n'
collate_alpha.py -i rarefactions/alpha_metrics -o rarefactions/collated_alpha

# Perform R analysis of reads per sample and rarefaction curves
if [[ $FACTORS -eq 0 ]]; then
    Rscript $(dirname $0)/SampleQualityandRarefaction.R -p $PROJ -m $MAP -b ${PROJ}_otu_table_no_singletons_summary.txt -a rarefactions/collated_alpha/
else
    Rscript $(dirname $0)/SampleQualityandRarefaction.R -p $PROJ -m $MAP -b ${PROJ}_otu_table_no_singletons_summary.txt -a rarefactions/collated_alpha/ -f FACTORS    
fi

echo 'Rarefaction: Done'

