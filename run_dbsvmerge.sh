BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
FILELIST=$1
prefix=$2
DISCSUPPVEC=$3

DBSVMERGEPATH=$BINDIR'/svmerge'
if [ ! -r $DBSVMERGEPATH/build/bin/dbSV_merge ]
then
    curdir=`pwd`
    cd $BINDIR
	git submodule update --init --recursive
	cd svmerge/build
	cmake ..
	make
	cd $curdir
fi
JASMINESRCPATH=$BINDIR'/'../../../Jasmine/src
SVPOPANALYSISPATH=$BINDIR'/src'

javac -cp $JASMINESRCPATH $SVPOPANALYSISPATH/*.java

WORKINGDIR=$4

names=''
tablefilelist=$WORKINGDIR/$prefix'_'tablefilelist.txt
if [ -r $tablefilelist ]
then
  rm $tablefilelist
fi

if [ ! -r $WORKINGDIR/$prefix'.dbsvmerge_augmented.txt' ]
then
	outprefix=$WORKINGDIR/$prefix'_'
	mergedfile=$outprefix'_dbsvmerge_merged.tsv'
	$DBSVMERGEPATH/build/bin/dbSV_merge  -f $FILELIST -o $mergedfile -l 2.0 -r 0.4

	echo 'Building merging table'
	java -cp $JASMINESRCPATH:$SVPOPANALYSISPATH BuildMergingTable vcf_file=$mergedfile out_file=$WORKINGDIR/$prefix'.dbsvmerge_simple.txt' vcf_filelist=$FILELIST mode=dbsvmerge sample_list=$names
	echo 'Augmenting merging table'
	java -cp $JASMINESRCPATH:$SVPOPANALYSISPATH AugmentMergingTable table_file=$WORKINGDIR/$prefix'.dbsvmerge_simple.txt' out_file=$WORKINGDIR/$prefix'.dbsvmerge_augmented.txt' vcf_filelist=$FILELIST
fi

ERRORSFILE=$WORKINGDIR'/'$prefix'.errors.txt'

appendflag=''

if [ -r $ERRORSFILE ]
then
  appendflag='--append'
fi

java -cp $JASMINESRCPATH:$SVPOPANALYSISPATH CountMergingErrors table_file=$WORKINGDIR/$prefix'.dbsvmerge_simple.txt' vcf_filelist=$FILELIST out_file=$ERRORSFILE software=dbsvmerge disc_supp_vec=$DISCSUPPVEC $appendflag
