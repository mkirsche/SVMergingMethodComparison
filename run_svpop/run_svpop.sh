BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")

FILELIST='/home/mkirsche/jasmine_data/figures/figure2/crosstool.filelist.fixed.txt'
OUTFILE=''

SVPOPPATH=$BINDIR'/../svpop'
JASMINESRCPATH='/home/mkirsche/jasmine_data/Jasmine/src'
SVPOPANALYSISPATH=$BINDIR'/../src'
DISCSUPPVEC='100'

javac -cp $JASMINESRCPATH $SVPOPANALYSISPATH/*.java

prefix=$1
WORKINGDIR=`pwd`

names=''
tablefilelist=$BINDIR/$prefix'_'tablefilelist.txt
if [ -r $tablefilelist ]
then
  rm $tablefilelist
fi

for i in `cat $FILELIST`
do 
  echo $i
  bn=`basename $i`
  tablefile=$BINDIR'/'$bn.tab
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

outprefix=$BINDIR/$prefix'_'
python $SVPOPPATH/merge.py $tablefilelist $names $outprefix
mergedfile=$outprefix'merged_all.bed'
echo 'Building merging table'
java -cp $JASMINESRCPATH:$SVPOPANALYSISPATH BuildMergingTable vcf_file=$mergedfile out_file=$outprefix'simple.txt' vcf_filelist=$FILELIST mode=svpop sample_list=$names
echo 'Augmenting merging table'
java -cp $JASMINESRCPATH:$SVPOPANALYSISPATH AugmentMergingTable table_file=$outprefix'simple.txt' out_file=$outprefix'augmented.txt' vcf_filelist=$FILELIST

ERRORSFILE=$WORKINGDIR'/'$prefix'.errors.txt'

appendflag=''

if [ -r $ERRORSFILE ]
then
  appendflag='--append'
fi

java -cp $JASMINESRCPATH:$SVPOPANALYSISPATH CountMergingErrors table_file=$outprefix'simple.txt' vcf_filelist=$FILELIST out_file=$ERRORSFILE software=svpop disc_supp_vec=$DISCSUPPVEC $appendflag
