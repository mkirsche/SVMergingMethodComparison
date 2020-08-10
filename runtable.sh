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

svtoolslist=$filelist_noext'.svtools.txt' 
echo 'Tabulating svtools results'
java -cp $JASMINE_PATH/src:$BINDIR/src BuildMergingTable vcf_file=$WORKINGDIR/$OUTPREFIX.svtools.vcf vcf_filelist=$FILELIST out_file=$WORKINGDIR/$OUTPREFIX.svtools_simple.txt mode=svtools
java -cp $JASMINE_PATH/src:$BINDIR/src AugmentMergingTable table_file=$WORKINGDIR/$OUTPREFIX.svtools_simple.txt vcf_filelist=$FILELIST out_file=$WORKINGDIR/$OUTPREFIX.svtools_augmented.txt

echo 'Tabulating Jasmine results'
java -cp $JASMINE_PATH/src:$BINDIR/src BuildMergingTable vcf_file=$WORKINGDIR/$OUTPREFIX.jasmine.vcf vcf_filelist=$FILELIST out_file=$WORKINGDIR/$OUTPREFIX.jasmine_simple.txt mode=jasmine
java -cp $JASMINE_PATH/src:$BINDIR/src AugmentMergingTable table_file=$WORKINGDIR/$OUTPREFIX.jasmine_simple.txt vcf_filelist=$FILELIST out_file=$WORKINGDIR/$OUTPREFIX.jasmine_augmented.txt

echo 'Tabulating SURVIVOR results'
java -cp $JASMINE_PATH/src:$BINDIR/src BuildMergingTable vcf_file=$WORKINGDIR/$OUTPREFIX.survivor.vcf vcf_filelist=$FILELIST out_file=$WORKINGDIR/$OUTPREFIX.survivor_simple.txt mode=survivor
java -cp $JASMINE_PATH/src:$BINDIR/src AugmentMergingTable table_file=$WORKINGDIR/$OUTPREFIX.survivor_simple.txt vcf_filelist=$FILELIST out_file=$WORKINGDIR/$OUTPREFIX.survivor_augmented.txt

svimmerlist=$filelist_noext'.svimmer.txt'
echo 'Tabulating svimmer results'
java -cp $JASMINE_PATH/src:$BINDIR/src BuildMergingTable vcf_file=$WORKINGDIR/$OUTPREFIX.svimmer.vcf vcf_filelist=$FILELIST out_file=$WORKINGDIR/$OUTPREFIX.svimmer_simple.txt mode=svimmer
java -cp $JASMINE_PATH/src:$BINDIR/src AugmentMergingTable table_file=$WORKINGDIR/$OUTPREFIX.svimmer_simple.txt vcf_filelist=$FILELIST out_file=$WORKINGDIR/$OUTPREFIX.svimmer_augmented.txt

