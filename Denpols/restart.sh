#!/bin/bash

HERE=$(pwd)

for i in 0 1
do

	NEWFOLD="RUN_"$i
        echo $NEWFOLD
        #mkdir -p $NEWFOLD

        cd $NEWFOLD
	jid=$(sbatch "dendrongpu.slrm" | cut -d ' ' -f 4)

	for i in {1..5}
	do
        	cp "dendrongpu.slrm" "dendrongpu_1.slrm"
        	sed -i '6s/.*/#SBATCH -d afterany:'$jid'/' "dendrongpu_1.slrm"
		jid=$(sbatch "dendrongpu_1.slrm" | cut -d ' ' -f 4)
	
        done

 	cd $HERE
done

