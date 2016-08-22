#!/bin/bash
# QiimeRunDiversity.sh
# Elizabeth Grice Lab
# University of Pennsylvania
#
# Arguments {project name}, {otu_table}, {# of seqs to keep} 
# Written using QIIME 1.8

########################################################################
#                              Set Up                                  #
########################################################################
# Initialize Flags
WORKDIR=$(pwd)
NUMCORE=1
SUB=0
ALPHA="observed_species,PD_whole_tree,shannon"
BETA="bray_curtis,unweighted_unifrac,weighted_normalized_unifrac"
COUNTS=0

# extract command line options and their arguments into variables based on this while loop taking in conditional position arguments.
usage() { echo -e "\nUsage: This is the second part of standardized Grice Lab QIIME workflow, which generates the diversity data. It has the following options: \n\n
    [-j] Project ID name (required) \n
    [-w] Working directory (required) \n 
    [-o] OTU table to use(required).\n 
    [-r] Aligned rep set to use. Recommend the no_singletons\n
         one (required).\n 
    [-c] Minimum reads per sample value. Will also be used \n
         for subsampling depth if -s flag given. (required) \n 
    [-n] Number of proc cores (optional) \n
    [-s] Subsample to even depth (optional)\n
    [-a] Alpha metrics to use besides default of observed species,\n
         shannon, & pd whole tree. Write 0 if you don't want anything\n
         done. (optional) \n 
    [-b] Beta metrics to use besides default of bray-curtis, weighted,\n
         & unweighted unifrac. Write 0 if you don't want anything \n
         done. (optional)\n
    [-l] Relative abundance tables in counts.\n"; exit 1;}

#
while getopts ":hj:w:n:o:r:c:sa:b:l" option; do
    case "$option" in
        h) usage ;;
        j) PROJ="$OPTARG";;
        w) WORKDIR="$OPTARG";;
        n) NUMCORE="$OPTARG";;
        o) OTU="$OPTARG";;
        r) REP="$OPTARG";;
        c) CUTOFF="$OPTARG";;
        s) SUB=1;;
        a) ALPHA="$OPTARG";;
        b) BETA="$OPTARG";;
        l) COUNTS=1;;
        :) echo "ERROR: -${OPTARG} requires an argument" ; exit 1 ;;
        ?) echo "ERROR: unknown option -${OPTARG}" ; exit 1 ;;
    esac
done

if [[ -z $PROJ ]]; then
    echo "ERROR: ProjectID is required"; exit 1;
elif [[ -z $WORKDIR ]]; then
    echo "ERROR: Working directory is required"; exit 1;
elif [[ -z $OTU ]]; then
    echo "ERROR: OTU table is required"; exit 1;
elif [[ -z $REP ]]; then
    echo "ERROR: Rep set is required"; exit 1;
elif [[ -z $CUTOFF ]]; then
    echo "ERROR: Cutoff is required"; exit 1;
fi

# Print staring variables
echo ProjectID: ${PROJ}
echo Working Directory: ${WORKDIR}
echo Number of Cores: ${NUMCORE}
echo OTU Table: ${OTU}
echo Rep Set: ${REP}
echo Cutoff: ${CUTOFF}
if [ $SUB -eq 1 ]; then
    echo Subsampling: Yes
else
    echo Subsampling: No
fi
echo Alpha Metrics: ${ALPHA}
echo Beta Metrics: ${BETA}
if [ $COUNTS -eq 0 ]; then
    echo Relative Abundance: Proportions
else
    echo Relative Abundance: Counts
fi

# Set the working directory
cd ${WORKDIR}
echo -e 'Changing working directory to' ${WORKDIR} '\n'

# Until QIIME module is fixed...
source /opt/software/qiime/1.8.0/activate.sh

# Depending on subsampling state either remove samples or subsample
if [ $SUB -eq 0 ] && [ $CUTOFF -eq 0 ]; then
    # Don't remove samples if the cutoff is 0
    echo 'Processing: Keeping All Samples'
    OTUNEW=$OTU
elif [ $SUB -eq 0 ]; then
    # Remove samples with <= cutoff sequences. This number is informed by the rarefaction plots 
    # from processing script.
    echo 'Processing: Removing samples with less than' ${CUTOFF} 'sequences'
    filter_samples_from_otu_table.py -i $OTU -o ${PROJ}_otu_table.filtered${CUTOFF}.biom -n $CUTOFF
    OTUNEW=${PROJ}_otu_table.filtered${CUTOFF}.biom
else
    # Subsample the OTU table
    echo 'Processing: Subsampling at' ${CUTOFF} 'per sample'
    single_rarefaction.py -i $OTU -d $CUTOFF -o ${PROJ}_otu_table.subsampled${CUTOFF}.biom
    OTUNEW=${PROJ}_otu_table.subsampled${CUTOFF}.biom
fi

#Print out stats for each sample
biom summarize-table -i $OTUNEW -o ${OTUNEW%.*}.summary.txt

# Subset the repset to the working set if samples were removed
echo 'Processing: Cleaning up rep_set'
if [[ $CUTOFF -eq 0 ]]; then
    REPNEW=$REP
elif [[ $SUB -eq 0 ]]; then
    REPNEW=${REP%.*}_filtered${CUTOFF}.fasta
    filter_fasta.py -f $REP -o $REPNEW -b $OTUNEW
else
    REPNEW=${REP%.*}_subsampled${CUTOFF}.fasta
    filter_fasta.py -f $REP -o $REPNEW -b $OTUNEW
fi

########################################################################
#                          Create Phylogeny                            #
########################################################################

# Remove ambiguously-aligned regions
echo 'Phylogeny: Removing ambiguous alignments'
REPDIR=$(dirname ${REPNEW})
filter_alignment.py -i $REPNEW -o ${REPDIR}/lanemasked -m /project/egricelab/references/greengenes/greengenes_alignment/lanemask_in_1s_and_0s

# Infer a phylogeny of reference sequences
echo 'Phylogeny: Making phylogeny'
REPNAME=${REPNEW##*/}
TREE=${REPNAME%.*}.tre
make_phylogeny.py -i $REPDIR/lanemasked/${REPNAME%.*}_pfiltered.fasta -o ${TREE} -t fasttree

########################################################################
#                         Diversity Analysis                           #
########################################################################

# Make taxonomic summaries of each sample (at various taxonomic levels, each with relative abundances)
echo 'Diversity: Relative Abundance'
if [ $COUNTS -eq 0 ]; then
    summarize_taxa.py -i $OTUNEW -o abundance_relative/ -L 2,3,4,5,6,7
else
    summarize_taxa.py -i $OTUNEW -o abundance_counts/ -a -L 2,3,4,5,6,7
fi

# Calculate within-sample diversity using an OTU table
if [[ $ALPHA != 0 ]]; then
    #statements
    echo 'Diversity: Making alpha tables'
    mkdir -p parallel_alpha_otu_table
    cp $OTUNEW parallel_alpha_otu_table/
    parallel_alpha_diversity.py -O $NUMCORE -i parallel_alpha_otu_table -o alpha -m $ALPHA -t $TREE
    rm -r parallel_alpha_otu_table
fi

#Create distance/dissimilarity matrices based on several metrics
if [[ $BETA != 0 ]]; then
    echo 'Diversity: Making beta tables'
    parallel_beta_diversity.py -O $NUMCORE -i $OTUNEW -o beta -t $TREE -m $BETA
fi

