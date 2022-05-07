gen=5
for j in 250 500
do
for i  in 0 1 2 3 4 5 
do 
echo 'G'$gen'_'$j
python3 analysis.py 'G'$gen'_'$j/RUN_$i/ traj.gsd  $gen $j &  
done
done 
wait 

gen=4
for j in 250 500
do
for i  in 0 1 2 3 4 5
do
echo 'G'$gen'_'$j
python3 analysis.py 'G'$gen'_'$j/RUN_$i/ traj.gsd  $gen $j &
done
done
wait



gen=3
for j in 250 500
do
for i  in 0 1 2 3 4 5
do
echo 'G'$gen'_'$j
python3 analysis.py 'G'$gen'_'$j/RUN_$i/ traj.gsd  $gen $j &
done
done
wait


j=1000

python3 analysis.py 'G'3'_'$j/RUN_0/ all_traj.gsd  3 $j & python3 analysis.py 'G'$3'_'$j/RUN_1/ all_traj.gsd  3 $j &  python3 analysis.py 'G'$4'_'$j/RUN_0/ all_traj.gsd  4 $j &  python3 analysis.py 'G'$4'_'$j/RUN_2/ all_traj.gsd  4 $j &  python3 analysis.py 'G'$5'_'$j/RUN_1/ traj_1.gsd  5 $j & python3 analysis.py 'G'$5'_'$j/RUN_0/ traj_1.gsd  5 $j & python3 analysis.py 'G'$5'_'$j/RUN_2/ traj_1.gsd  5 $j



