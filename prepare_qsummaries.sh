#!/bin/bash

cov="4X"

for p in 0.1 0.3 0.5 0.7 0.9
do  
	for div in 5000 10000 25000 50000
	do  
		for ref in S1 S2 S3
		do   
#l=$(./avg_fngs_single.sh $ref.$p.$div.*.4X.ngsadmix_beagleGL_ref*.qopt) # >> q.summary.4X.uncorr.txt
#echo $l 4X uncorr >> q.summary.all.txt
#l=$(./avg_fngs_single.sh $ref.$p.$div.*.4X.ngsadmix_beagleGL_ref*.corrected.qopt) # >> q.summary.4X.corr.txt
#echo $l 4X corr >> q.summary.all.txt
#l=$(./avg_fngs_single.sh $ref.$p.$div.*.ngsadmix_beagleGL_ref*.qopt) # >> q.summary.1X.uncorr.txt
#echo $l 1X uncorr >> q.summary.all.txt
#l=$(./avg_fngs_single.sh $ref.$p.$div.*.ngsadmix_beagleGL_ref*.corrected.qopt) #>> q.summary.1X.corr.txt
#echo $l 1X corr >> q.summary.all.txt

			for i in {1..50}
			do
				rm tmp.q
				for file in $(ls fastngsadmix_results/$ref.$p.$div.$i.ngsadmix_beagleGL_ref*.qopt 2> /dev/null | grep -v 'corr')
				do
					tail -n 1 $file >> tmp.q
				done
				if [ -f tmp.q ]
					then
					S2_avg=$(awk '{ sum += $1 } END { print sum / NR }' tmp.q)
					S3_avg=$(awk '{ sum += $2 } END { print sum / NR }' tmp.q)
					echo $ref $p $div $S2_avg $S3_avg $cov uncorr FNGSA >> q.summary.$cov.txt
				fi
			done


			for i in {1..50}
			do
				rm tmp.q
				for file in $(ls fastngsadmix_results/$ref.$p.$div.$i.ngsadmix_beagleGL_ref*.corrected.qopt 2> /dev/null)
				do
					tail -n 1 $file >> tmp.q
				done
				if [ -f tmp.q ]
					then
					S2_avg=$(awk '{ sum += $1 } END { print sum / NR }' tmp.q)
					S3_avg=$(awk '{ sum += $2 } END { print sum / NR }' tmp.q)
					echo $ref $p $div $S2_avg $S3_avg $cov corr FNGSA >> q.summary.$cov.txt					
				fi 
			done


			for file in $(ls admixture_results/$ref.$p.$div.*.admixture.2.Q 2> /dev/null)
			do
				Rscript summarize_admixture.R $file $ref $p $div $cov >> q.summary.$cov.txt
			done

		done  
	done
done




