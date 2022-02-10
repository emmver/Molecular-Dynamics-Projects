for i in {1..20} 
do
echo $i  
~/espresso/build/pypresso initial_poly.py /gpfs/data/fs71642/manos4/single_chain_dilute_magnetics/RUN_$i &
done 
