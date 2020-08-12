# Get source and working directories

if [ "$(uname -s)" = 'Linux' ]; then
    BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
else
    BINDIR=$(dirname "$(readlink "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
fi

WORKINGDIR=`pwd`

# Locations of necessary software

JASMINE_UTILS_PATH='/home/mkirsche/git/Jasmine_utils'
JASMINE_PATH='/home/mkirsche/eclipse-workspace/Jasmine'
SURVIVOR_PATH='/home/mkirsche/git/SURVIVOR'
SVIMMER_PATH='/home/mkirsche/git/svimmer'

javac -cp $JASMINE_PATH/src $BINDIR/src/*.java

# Command line parameters
FILELIST=$1
OUTPREFIX=$2

# Remove extension from filelist for lsmaking caller-specific versions of it
filelist_noext=`echo "${FILELIST%.*}"`

source /home/mkirsche/anaconda3/etc/profile.d/conda.sh
conda activate py2

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
    svtools lsort -f $svtoolslist > $WORKINGDIR/$OUTPREFIX.svtools.lsort.vcf

    svtools lmerge -i $WORKINGDIR/$OUTPREFIX.svtools.lsort.vcf -f 1000 > $WORKINGDIR/$OUTPREFIX.svtools.lmerge.vcf

    awk '($0 !~ /MATEID=[^;]+_1;/)' $OUTPREFIX.svtools.lmerge.vcf > $WORKINGDIR/$OUTPREFIX.svtools.vcf

    echo 'Tabulating svtools results'
    java -cp $JASMINE_PATH/src:$BINDIR/src BuildMergingTable vcf_file=$WORKINGDIR/$OUTPREFIX.svtools.vcf vcf_filelist=$FILELIST out_file=$WORKINGDIR/$OUTPREFIX.svtools_simple.txt mode=svtools

    java -cp $JASMINE_PATH/src:$BINDIR/src AugmentMergingTable table_file=$WORKINGDIR/$OUTPREFIX.svtools_simple.txt vcf_filelist=$FILELIST out_file=$WORKINGDIR/$OUTPREFIX.svtools_augmented.txt

fi

if [ ! -r $WORKINGDIR/$OUTPREFIX.jasmine_augmented.txt ]
then

    echo 'Running Jasmine'
    $JASMINE_PATH/jasmine file_list=$FILELIST out_file=$WORKINGDIR/$OUTPREFIX.jasmine.vcf

    echo 'Tabulating Jasmine results'
    java -cp $JASMINE_PATH/src:$BINDIR/src BuildMergingTable vcf_file=$WORKINGDIR/$OUTPREFIX.jasmine.vcf vcf_filelist=$FILELIST out_file=$WORKINGDIR/$OUTPREFIX.jasmine_simple.txt mode=jasmine

    java -cp $JASMINE_PATH/src:$BINDIR/src AugmentMergingTable table_file=$WORKINGDIR/$OUTPREFIX.jasmine_simple.txt vcf_filelist=$FILELIST out_file=$WORKINGDIR/$OUTPREFIX.jasmine_augmented.txt

fi

if [ ! -r $WORKINGDIR/$OUTPREFIX.jasmineintra_augmented.txt ]
then

    echo 'Running Jasmine (with intrasample)'
    $JASMINE_PATH/jasmine --allow_intrasample file_list=$FILELIST out_file=$WORKINGDIR/$OUTPREFIX.jasmineintra.vcf --use_end kd_tree_norm=1

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
    chrs='1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT'

    echo 'Running svimmer'
    python $SVIMMER_PATH/svimmer $svimmerlist $chrs --threads 2 --output $WORKINGDIR/$OUTPREFIX.svimmer.vcf --max_distance 1000 --max_size_difference 1000 --ids

    echo 'Tabulating svimmer results'

    java -cp $JASMINE_PATH/src:$BINDIR/src BuildMergingTable vcf_file=$WORKINGDIR/$OUTPREFIX.svimmer.vcf vcf_filelist=$FILELIST out_file=$WORKINGDIR/$OUTPREFIX.svimmer_simple.txt mode=svimmer

    java -cp $JASMINE_PATH/src:$BINDIR/src AugmentMergingTable table_file=$WORKINGDIR/$OUTPREFIX.svimmer_simple.txt vcf_filelist=$FILELIST out_file=$WORKINGDIR/$OUTPREFIX.svimmer_augmented.txt

fi
