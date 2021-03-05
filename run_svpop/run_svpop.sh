BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
FILELIST=$1
prefix=$2
DISCSUPPVEC=$3

SVPOPPATH=$BINDIR'/../svpop'
JASMINESRCPATH=$BINDIR'/'../../../../Jasmine/src
SVPOPANALYSISPATH=$BINDIR'/../src'

javac -cp $JASMINESRCPATH $SVPOPANALYSISPATH/*.java

WORKINGDIR=`pwd`

names=''
tablefilelist=$WORKINGDIR/$prefix'_'tablefilelist.txt
if [ -r $tablefilelist ]
then
  rm $tablefilelist
fi

if [ ! -r $WORKINGDIR/$prefix'.svpop_augmented.txt' ]
then
	for i in `cat $FILELIST`
	do 
	  echo $i
	  bn=`basename $i`
	  tablefile=$WORKINGDIR'/'$bn.tab
	  if [ ! -r $tablefile ]
	  then
		echo 'Converting '$i' to tsv'
		vcf2tsv -n NA -g $i > $tablefile
	  fi
	  
	  echo $tablefile >> $tablefilelist
	  echo 'Writing per-type tables for '$i
	  python $SVPOPPATH'/'write_tables.py $tablefile
	  
	  name=${bn%%.*}
	  name=${bn%%vGRCh38*}
	  echo $name
	  
	  if [ ! -z $names ]
	  then
		names=$names','
	  fi
	  names=$names''$name
	done
	cat $tablefilelist
	echo $names

	outprefix=$WORKINGDIR/$prefix'_'
	python $SVPOPPATH/merge.py $tablefilelist $names $outprefix
	mergedfile=$outprefix'merged_all.bed'
	echo 'Building merging table'
	java -cp $JASMINESRCPATH:$SVPOPANALYSISPATH BuildMergingTable vcf_file=$mergedfile out_file=$WORKINGDIR/$prefix'.svpop_simple.txt' vcf_filelist=$FILELIST mode=svpop sample_list=$names
	echo 'Augmenting merging table'
	java -cp $JASMINESRCPATH:$SVPOPANALYSISPATH AugmentMergingTable table_file=$WORKINGDIR/$prefix'.svpop_simple.txt' out_file=$WORKINGDIR/$prefix'.svpop_augmented.txt' vcf_filelist=$FILELIST
fi

ERRORSFILE=$WORKINGDIR'/'$prefix'.errors.txt'

appendflag=''

if [ -r $ERRORSFILE ]
then
  appendflag='--append'
fi

java -cp $JASMINESRCPATH:$SVPOPANALYSISPATH CountMergingErrors table_file=$WORKINGDIR/$prefix'.svpop_simple.txt' vcf_filelist=$FILELIST out_file=$ERRORSFILE software=svpop disc_supp_vec=$DISCSUPPVEC $appendflag
