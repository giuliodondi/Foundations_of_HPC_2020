#!/bin/bash

NN=1000;

MAXN=1000000;

filename="add_these.txt"


if [ -f $filename ];
then
rm $filename
fi

touch $filename

./rand_gen.sh -n $NN -m $MAXN | awk -v fname=$filename 'NR==1 {for ( i=1; i<=NF; i++ ) { if ($i~/^[0-9]+/) {print $i >> fname} } } NR>1 { print $2  >> fname} ' #NR==1 { print $0  >> fname} 
#if ($i~/\d+/) {print $i >> fname;} 
