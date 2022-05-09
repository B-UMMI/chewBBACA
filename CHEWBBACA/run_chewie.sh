#! /bin/bash

WORKDIR=$1
SCHEMA=$2
DATASETS_DIR=$3

for i in {1,2,4,8,16,32,64,128,256,512,1024}
do
	cp -r $SCHEMA $WORKDIR;
	{ /usr/bin/time -f "%e\n%M\n%P" python chewBBACA.py AlleleCall -i $DATASETS_DIR/"$i" -g $WORKDIR/senterica_schema -o $WORKDIR/senterica"$i"_results --cpu 6 1> $WORKDIR/senterica"$i"_runlog ; } 2> $WORKDIR/senterica"$i"_timelog;
	rm -rf $WORKDIR/senterica_schema;
done
