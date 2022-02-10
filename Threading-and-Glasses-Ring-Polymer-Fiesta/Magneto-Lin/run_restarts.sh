for i in {1..20}
do 
for j in 1 2  
do 
~/espresso/build/pypresso restart_min_poly.py mycheck_passive $j /gpfs/data/fs71642/manos4/single_chain_dilute_magnetics/RUN_$i/lambda_$j &
done 
done
wait 


for i in {1..20}
do 
for j in  3 4 
do 
~/espresso/build/pypresso restart_min_poly.py mycheck_passive $j /gpfs/data/fs71642/manos4/single_chain_dilute_magnetics/RUN_$i/lambda_$j &
done 
done
wait 





for i in {1..20}
do 
for j in 5 6 
do 
~/espresso/build/pypresso restart_min_poly.py mycheck_passive $j /gpfs/data/fs71642/manos4/single_chain_dilute_magnetics/RUN_$i/labmda_$j &
done
done
