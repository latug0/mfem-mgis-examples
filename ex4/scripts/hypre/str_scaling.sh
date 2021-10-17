#!/bin/bash

#declare -a nodes=(2  4  8)
#declare -a cores=(32 32 32)
declare -a nodes=(1  1  )
declare -a cores=(16 32 )

### no ./ for executable
ex=Ssna303_3d_hypre 
baserep=hyprer
ndof=""
time_limit="1:00:00"

mkdir ${baserep}_000 2>/dev/null

lastrep=$(ls -rX | grep "^${baserep}" | head -n 1)
echo "lasterp=${lastrep}"
((number=${lastrep#${baserep}_}))
echo "number=${number}"
for i in ${!nodes[@]}; do
   ((number++))
   ((n=${nodes[$i]}))
   ((c=${cores[$i]}))
   ((t=$n*$c))
   if [ $t -lt 10 ]; then
	num='00'$t
   elif [ $t -lt 100 ]; then
	num='0'$t
   else 
	num=$t
   fi  
   newnum=$(printf "%3.3d" $((number)))
   rep="${baserep}_${newnum}"
   echo "rep=$rep"
   mkdir -p $rep
   FILE="slurm${num}.cmd"
   sed "s/NNN/${n}/;s/CCC/${c}/;s/FFF/${num}/;s/EEE/${ex}/;s/JJJ/${rep}/;s/TTT/${ndof}/;s/LLL/${time_limit}/" ./slurmxxx.cmd > $rep/$FILE
   echo "Launching $ex $ndof in $rep/$FILE with $t cores" 
    
   (cd $rep ; 
    ln -s ../Ssna303_3d_hypre
    ln -s ../ssna303_3d.msh 
    ln -s ../src
    ln -s ../Plasticity.mfront 
    ln -s ../include
       sbatch $FILE)
   echo " " $FILE
done

rmdir  ${baserep}_000 2>/dev/null
