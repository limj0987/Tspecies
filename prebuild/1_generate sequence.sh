#!/bin/bash
Ne2_initial=1000  #the theoretical 2*Ne
Ne2_step=4950     
Ne2_end=1000000   
Tspecies_initial=10000     # the theoretical divergence time 
Tspecies_step=165000
Tspecies_end=1000000
miu=0.00000001    #change the miu and generate data

echo "table_fixed_fixed" > result_t_1.csv
echo -e "miu,Tspecies,2Ne,mean,variance" >> result_t_1.csv
Tspecies=$(echo "$Tspecies_initial" | bc -l)
while [ "$(echo "$Tspecies <= $Tspecies_end" | bc -l)" -eq 1 ]
do
   Ne2=$(echo "$Ne2_initial" | bc -l)
   while [ "$(echo "$Ne2 <= $Ne2_end" | bc -l)" -eq 1 ]
   do
       theta=$(echo "scale =10; 2*$Ne2*$miu" | bc -l)
       tx=$(echo "scale=10; $Tspecies/(2*$Ne2)" | bc -l)
       output_tree="treefile"
       output_seq="seqfile"
       out1="output.txt"

       ms 2 10000 -T -t "$theta" -I 2 1 1 -n 1 1 -n 2 1 -ej "$tx" 2 1 > "$out1"
       grep -v "//" "$out1" > "${output_tree}"
       seq-gen -mHKY -l 1000 -s "$theta" <"${output_tree}" > "${output_seq}"

       mean_and_variance=$(Rscript correct.R "${output_seq}")
       last_two_lines=$(echo "$mean_and_variance" | tail -n 2)
       read mean variance <<< "$last_two_lines"
       echo -e "${miu},${Tspecies},${Ne2},${mean},${variance}" >> result_t_1.csv

       Ne2=$(echo "$Ne2 + $Ne2_step" | bc -l)
    done
       Tspecies=$(echo "$Tspecies + $Tspecies_step" | bc -l)
done
