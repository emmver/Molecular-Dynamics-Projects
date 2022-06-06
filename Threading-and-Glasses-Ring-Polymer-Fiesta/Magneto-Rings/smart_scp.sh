for i in {1..8}; do sshpass -p 1Q@w3A4s5Z scp -r manos4@vsc4.vsc.ac.at:/gpfs/data/fs71642/manos4/single_ring_dilute_magnetics/RUN_$i/lambda_12/* ./RUN_$i/lambda_12/;echo $i;  done

