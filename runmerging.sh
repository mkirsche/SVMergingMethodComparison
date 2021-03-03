# Get source and working directories

if [ "$(uname -s)" = 'Linux' ]; then
    BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
else
    BINDIR=$(dirname "$(readlink "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
fi

WORKINGDIR=`pwd`
echo 'Working directory: '$WORKINGDIR
echo $1 $2 $3

# Locations of necessary software

JASMINE_UTILS_PATH=$BINDIR'/'Jasmine_utils
JASMINE_PATH=$BINDIR'/'../../../Jasmine
SURVIVOR_PATH=$BINDIR'/'../../../SURVIVOR
SVIMMER_PATH=$BINDIR'/'../../figure5/svimmer

javac $JASMINE_PATH/Iris/src/*.java
javac -cp $JASMINE_PATH/Iris/src $JASMINE_PATH/src/*.java
javac -cp $JASMINE_PATH/src $BINDIR/src/*.java

# Command line parameters
FILELIST=$1
OUTPREFIX=$2
DISCSUPPVEC=$3

# Remove extension from filelist for making caller-specific versions of it
filelist_noext=`echo "${FILELIST%.*}"`

#source /home/mkirsche/anaconda3/etc/profile.d/conda.sh
#conda activate py2

if [ ! -r $WORKINGDIR/$OUTPREFIX.svtools_augmented.txt ]
then

    svtoolslist=$filelist_noext'.svtools.txt'
    if [ -r $svtoolslist ]
    then
      rm $svtoolslist
    fi

    echo 'Preprocessing for svtools'
    for filename in `cat $FILELIST`
    do
      noext=`echo "${filename%.*}"`
      outfile=$noext.svtools.vcf 
      echo $outfile >> $svtoolslist
      python $JASMINE_UTILS_PATH/sniffles2svtools.py -o $outfile $filename
    done

    echo 'Running svtools'
    svtools lsort -r -f $svtoolslist > $WORKINGDIR/$OUTPREFIX.svtools.lsort.vcf

    svtools lmerge -i $WORKINGDIR/$OUTPREFIX.svtools.lsort.vcf -f 1000 > $WORKINGDIR/$OUTPREFIX.svtools.lmerge.vcf

    awk '($0 !~ /MATEID=[^;]+_1;/)' $OUTPREFIX.svtools.lmerge.vcf > $WORKINGDIR/$OUTPREFIX.svtools.vcf

    echo 'Tabulating svtools results'
    java -cp $JASMINE_PATH/src:$BINDIR/src BuildMergingTable vcf_file=$WORKINGDIR/$OUTPREFIX.svtools.vcf vcf_filelist=$FILELIST out_file=$WORKINGDIR/$OUTPREFIX.svtools_simple.txt mode=svtools

    java -cp $JASMINE_PATH/src:$BINDIR/src AugmentMergingTable table_file=$WORKINGDIR/$OUTPREFIX.svtools_simple.txt vcf_filelist=$FILELIST out_file=$WORKINGDIR/$OUTPREFIX.svtools_augmented.txt

fi

if [ ! -r $WORKINGDIR/$OUTPREFIX.jasmine_augmented.txt ]
then

    echo 'Running Jasmine'
    java -cp $JASMINE_PATH/Iris/src:$JASMINE_PATH/src Main file_list=$FILELIST out_file=$WORKINGDIR/$OUTPREFIX.jasmine.vcf

    echo 'Tabulating Jasmine results'
    java -cp $JASMINE_PATH/src:$BINDIR/src BuildMergingTable vcf_file=$WORKINGDIR/$OUTPREFIX.jasmine.vcf vcf_filelist=$FILELIST out_file=$WORKINGDIR/$OUTPREFIX.jasmine_simple.txt mode=jasmine

    java -cp $JASMINE_PATH/src:$BINDIR/src AugmentMergingTable table_file=$WORKINGDIR/$OUTPREFIX.jasmine_simple.txt vcf_filelist=$FILELIST out_file=$WORKINGDIR/$OUTPREFIX.jasmine_augmented.txt

fi

if [ ! -r $WORKINGDIR/$OUTPREFIX.jasmineintra_augmented.txt ]
then

    echo 'Running Jasmine (with intrasample)'
    $JASMINE_PATH/jasmine --allow_intrasample file_list=$FILELIST out_file=$WORKINGDIR/$OUTPREFIX.jasmineintra.vcf kd_tree_norm=1 --use_end

    #echo 'Tabulating Jasmine results'
    java -cp $JASMINE_PATH/src:$BINDIR/src BuildMergingTable vcf_file=$WORKINGDIR/$OUTPREFIX.jasmineintra.vcf vcf_filelist=$FILELIST out_file=$WORKINGDIR/$OUTPREFIX.jasmineintra_simple.txt mode=jasmine_intra

    java -cp $JASMINE_PATH/src:$BINDIR/src AugmentMergingTable table_file=$WORKINGDIR/$OUTPREFIX.jasmineintra_simple.txt vcf_filelist=$FILELIST out_file=$WORKINGDIR/$OUTPREFIX.jasmineintra_augmented.txt

fi

if [ ! -r $WORKINGDIR/$OUTPREFIX.survivor_augmented.txt ]
then

    echo 'Running SURVIVOR'
    $SURVIVOR_PATH/Debug/SURVIVOR merge $FILELIST 1000 1 1 1 0 1 $WORKINGDIR/$OUTPREFIX.survivor.vcf

    echo 'Tabulating SURVIVOR results'
    java -cp $JASMINE_PATH/src:$BINDIR/src BuildMergingTable vcf_file=$WORKINGDIR/$OUTPREFIX.survivor.vcf vcf_filelist=$FILELIST out_file=$WORKINGDIR/$OUTPREFIX.survivor_simple.txt mode=survivor

    java -cp $JASMINE_PATH/src:$BINDIR/src AugmentMergingTable table_file=$WORKINGDIR/$OUTPREFIX.survivor_simple.txt vcf_filelist=$FILELIST out_file=$WORKINGDIR/$OUTPREFIX.survivor_augmented.txt

fi

conda activate py3

if [ ! -r $WORKINGDIR/$OUTPREFIX.svimmer_augmented.txt ]
then

    echo 'Preprocessing for svimmer'
    svimmerlist=$filelist_noext'.svimmer.txt'

    if [ -r $svimmerlist ]
    then
      rm $svimmerlist
    fi

    for filename in `cat $FILELIST`
    do
      noext=`echo "${filename%.*}"`
      outfile=$noext.svimmer.vcf.gz 
      echo $outfile >> $svimmerlist
      python $JASMINE_UTILS_PATH/sniffles2svimmer.py -o $outfile $filename
      tabix $outfile
    done
    #chrs='chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM'
    chrs=`java -cp $BINDIR/src GetChromosomeNames $FILELIST`
    echo 'Got chromosomes: '$chrs
    echo 'Running svimmer'
    python $SVIMMER_PATH/svimmer $svimmerlist $chrs --threads 2 --output $WORKINGDIR/$OUTPREFIX.svimmer.vcf --max_distance 1000 --max_size_difference 1000 --ids

    echo 'Tabulating svimmer results'

    java -cp $JASMINE_PATH/src:$BINDIR/src BuildMergingTable vcf_file=$WORKINGDIR/$OUTPREFIX.svimmer.vcf vcf_filelist=$FILELIST out_file=$WORKINGDIR/$OUTPREFIX.svimmer_simple.txt mode=svimmer

    java -cp $JASMINE_PATH/src:$BINDIR/src AugmentMergingTable table_file=$WORKINGDIR/$OUTPREFIX.svimmer_simple.txt vcf_filelist=$FILELIST out_file=$WORKINGDIR/$OUTPREFIX.svimmer_augmented.txt

fi

if [ "$DISCSUPPVEC" != "0" ]
then
  echo 'Counting errors'
  ERRORSFILE=$WORKINGDIR/$OUTPREFIX.errors.txt
  java -cp $JASMINE_PATH/src:$BINDIR/src CountMergingErrors table_file=$WORKINGDIR/$OUTPREFIX.svtools_simple.txt vcf_filelist=$FILELIST out_file=$ERRORSFILE software=svtools disc_supp_vec=$DISCSUPPVEC
  java -cp $JASMINE_PATH/src:$BINDIR/src CountMergingErrors table_file=$WORKINGDIR/$OUTPREFIX.jasmine_simple.txt vcf_filelist=$FILELIST out_file=$ERRORSFILE software=jasmine disc_supp_vec=$DISCSUPPVEC --append
  java -cp $JASMINE_PATH/src:$BINDIR/src CountMergingErrors table_file=$WORKINGDIR/$OUTPREFIX.jasmineintra_simple.txt vcf_filelist=$FILELIST out_file=$ERRORSFILE software=jasmineintra disc_supp_vec=$DISCSUPPVEC --append
  java -cp $JASMINE_PATH/src:$BINDIR/src CountMergingErrors table_file=$WORKINGDIR/$OUTPREFIX.survivor_simple.txt vcf_filelist=$FILELIST out_file=$ERRORSFILE software=survivor disc_supp_vec=$DISCSUPPVEC --append
  java -cp $JASMINE_PATH/src:$BINDIR/src CountMergingErrors table_file=$WORKINGDIR/$OUTPREFIX.svimmer_simple.txt vcf_filelist=$FILELIST out_file=$ERRORSFILE software=svimmer disc_supp_vec=$DISCSUPPVEC --append
fi


