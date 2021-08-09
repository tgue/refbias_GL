#!/bin/bash

cov="0.5X.noMAF"

#for p in 0.1 0.3 0.5 0.7 0.9
for p in 0.5
do  
	#for div in 5000 10000 25000 50000
	for div in 50000
	do  
		for ref in S1 S2 S3
		do   

			for file in $(ls admixture_results/$ref.$p.$div.*.supervised.admixture.2.Q 2> /dev/null)
			do
				Rscript summarize_admixture.R $file $ref $p $div $cov uncorr ADMIXTURE pruned haploid single supervised >> q.summary.$cov.txt
			done

			for file in $(ls admixture_results/$ref.$p.$div.*supervised.admixture.nonpruned.2.Q 2> /dev/null)
			do
				Rscript summarize_admixture.R $file $ref $p $div $cov uncorr ADMIXTURE nonpruned haploid single supervised >> q.summary.$cov.txt
			done

			for file in $(ls admixture_results/$ref.$p.$div.*.unsupervised.admixture.2.Q 2> /dev/null)
			do
				Rscript summarize_admixture.R $file $ref $p $div $cov uncorr ADMIXTURE pruned haploid single unsupervised >> q.summary.$cov.txt
			done

			for file in $(ls admixture_results/$ref.$p.$div.*unsupervised.admixture.nonpruned.2.Q 2> /dev/null)
			do
				Rscript summarize_admixture.R $file $ref $p $div $cov uncorr ADMIXTURE nonpruned haploid single unsupervised >> q.summary.$cov.txt
			done

			for file in $(ls admixture_results/$ref.$p.$div.*.supervised.multithread.admixture.2.Q 2> /dev/null)
			do
				Rscript summarize_admixture.R $file $ref $p $div $cov uncorr ADMIXTURE pruned haploid multi supervised >> q.summary.$cov.txt
			done

			for file in $(ls admixture_results/$ref.$p.$div.*supervised.multithread.admixture.nonpruned.2.Q 2> /dev/null)
			do
				Rscript summarize_admixture.R $file $ref $p $div $cov uncorr ADMIXTURE nonpruned haploid multi supervised >> q.summary.$cov.txt
			done

			for file in $(ls admixture_results/$ref.$p.$div.*.supervised.nonhaploid.admixture.2.Q 2> /dev/null)
			do
				Rscript summarize_admixture.R $file $ref $p $div $cov uncorr ADMIXTURE pruned diploid single supervised >> q.summary.$cov.txt
			done

			for file in $(ls admixture_results/$ref.$p.$div.*supervised.nonhaploid.admixture.nonpruned.2.Q 2> /dev/null)
			do
				Rscript summarize_admixture.R $file $ref $p $div $cov uncorr ADMIXTURE nonpruned diploid single supervised >> q.summary.$cov.txt
			done

			for i in {1..50}
			do
				rm tmp.q
				for file in $(ls fastngsadmix_results/$ref.$p.$div.$i.ngsadmix_beagleGL_ref*.qopt 2> /dev/null | grep -v 'corr' | grep -v 'pruned')
				do
					tail -n 1 $file >> tmp.q
				done
				if [ -f tmp.q ]
					then
					S2_avg=$(awk '{ sum += $1 } END { print sum / NR }' tmp.q)
					S3_avg=$(awk '{ sum += $2 } END { print sum / NR }' tmp.q)
					echo $ref $p $div $S2_avg $S3_avg $cov uncorr FNGSA nonpruned >> q.summary.$cov.txt
				fi
			done


			for i in {1..50}
			do
				rm tmp.q
				for file in $(ls fastngsadmix_results/$ref.$p.$div.$i.ngsadmix_beagleGL_ref*.corrected.qopt 2> /dev/null | grep -v "pruned")
				do
					tail -n 1 $file >> tmp.q
				done
				if [ -f tmp.q ]
					then
					S2_avg=$(awk '{ sum += $1 } END { print sum / NR }' tmp.q)
					S3_avg=$(awk '{ sum += $2 } END { print sum / NR }' tmp.q)
					echo $ref $p $div $S2_avg $S3_avg $cov corr FNGSA nonpruned >> q.summary.$cov.txt					
				fi 
			done

			for i in {1..50}
			do
				rm tmp.q

				for file in $(ls fastngsadmix_results/$ref.$p.$div.$i.ngsadmix_beagleGL_ref*.qopt 2> /dev/null | grep -v 'corr' | grep 'pruned')
				do
					tail -n 1 $file >> tmp.q
				done
				if [ -f tmp.q ]
					then
					S2_avg=$(awk '{ sum += $1 } END { print sum / NR }' tmp.q)
					S3_avg=$(awk '{ sum += $2 } END { print sum / NR }' tmp.q)
					echo $ref $p $div $S2_avg $S3_avg $cov uncorr FNGSA pruned >> q.summary.$cov.txt
				fi
			done


			for i in {1..50}
			do
				rm tmp.q
				for file in $(ls fastngsadmix_results/$ref.$p.$div.$i.ngsadmix_beagleGL_ref*.corrected.pruned.qopt 2> /dev/null | grep "pruned")
				do
					tail -n 1 $file >> tmp.q
				done
				if [ -f tmp.q ]
					then
					S2_avg=$(awk '{ sum += $1 } END { print sum / NR }' tmp.q)
					S3_avg=$(awk '{ sum += $2 } END { print sum / NR }' tmp.q)
					echo $ref $p $div $S2_avg $S3_avg $cov corr FNGSA pruned >> q.summary.$cov.txt					
				fi 
			done

			for file in $(ls pcangsd_results/$ref.$p.$div.*.all.pcangsd.corr.pruned.qopt 2> /dev/null)
			do
				Rscript summarize_admixture.R $file $ref $p $div $cov corr PCAngsd pruned >> q.summary.$cov.txt
			done

			for file in $(ls pcangsd_results/$ref.$p.$div.*.all.pcangsd.corr.qopt 2> /dev/null)
			do
				Rscript summarize_admixture.R $file $ref $p $div $cov corr PCAngsd nonpruned >> q.summary.$cov.txt
			done

			for file in $(ls pcangsd_results/$ref.$p.$div.*.all.pcangsd.uncorr.pruned.qopt 2> /dev/null)
			do
				Rscript summarize_admixture.R $file $ref $p $div $cov uncorr PCAngsd pruned >> q.summary.$cov.txt
			done

			for file in $(ls pcangsd_results/$ref.$p.$div.*.all.pcangsd.uncorr.qopt 2> /dev/null)
			do
				Rscript summarize_admixture.R $file $ref $p $div $cov uncorr PCAngsd nonpruned >> q.summary.$cov.txt
			done


			for file in $(ls ngsadmix_results/$ref.$p.$div.*.NGSADM.beagle.pruned.corrected.qopt 2> /dev/null)
			do
				Rscript summarize_admixture.R $file $ref $p $div $cov corr NGSadmix pruned >> q.summary.$cov.txt
			done

			for file in $(ls ngsadmix_results/$ref.$p.$div.*.NGSADM.beagle.corrected.qopt 2> /dev/null)
			do
				Rscript summarize_admixture.R $file $ref $p $div $cov corr NGSadmix nonpruned >> q.summary.$cov.txt
			done

			for file in $(ls ngsadmix_results/$ref.$p.$div.*.NGSADM.beagle.pruned.qopt 2> /dev/null)
			do
				Rscript summarize_admixture.R $file $ref $p $div $cov uncorr NGSadmix pruned >> q.summary.$cov.txt
			done

			for file in $(ls ngsadmix_results/$ref.$p.$div.*.NGSADM.beagle.qopt 2> /dev/null)
			do
				Rscript summarize_admixture.R $file $ref $p $div $cov uncorr NGSadmix nonpruned >> q.summary.$cov.txt
			done


			for file in $(ls qpadm_results/$ref.$p.$div.*.out 2> /dev/null | grep -v "pruned")
			do
				if grep -q "best coeff" $file
					then
					s2p=$(grep "best coeff" $file | awk '{print $3}')
					s3p=$(grep "best coeff" $file | awk '{print $4}')
					echo $ref $p $div $s2p $s3p $cov uncorr qpAdm nonpruned >> q.summary.$cov.txt	
				fi
			done

			for file in $(ls qpadm_results/$ref.$p.$div.*.out 2> /dev/null | grep "pruned")
			do
				if grep -q "best coeff" $file
					then
					s2p=$(grep "best coeff" $file | awk '{print $3}')
					s3p=$(grep "best coeff" $file | awk '{print $4}')
					echo $ref $p $div $s2p $s3p $cov uncorr qpAdm pruned >> q.summary.$cov.txt	
				fi
			done

		done  
	done
done




