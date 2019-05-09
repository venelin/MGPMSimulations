for id in `seq 257 1 384`
do
mkdir -p Results_t5_NULL_MGPM_A_F_best_clade_2_RR_N318/t5_NULL_MGPM_A_F_best_clade_2_RR_id_$id
cd Results_t5_NULL_MGPM_A_F_best_clade_2_RR_N318/t5_NULL_MGPM_A_F_best_clade_2_RR_id_$id
if [ -f "FinalResult_t5_NULL_MGPM_A_F_best_clade_2_RR_id_"$id".RData" ]
then
rm MPI*.log
rm CurrentResults*.RData
rm *worker*.RData
else
rm MPI*.log
#rm CurrentResults*.RData
rm *worker*.RData
bsub -M 250000 -n 50 -W 23:59 -R "rusage[mem=4096]" -R ib sh R --vanilla --slave -f ../../DetectShifts_t5_NULL_MGPM_A_F_best_clade_2_RR.R --args $id
fi
cd ../..
done


for id in `seq 1 1 128`
do
mkdir -p Results_t5_NULL_MGPM_A_F_best_clade_2_RR_N80/t5_NULL_MGPM_A_F_best_clade_2_RR_id_$id
cd Results_t5_NULL_MGPM_A_F_best_clade_2_RR_N80/t5_NULL_MGPM_A_F_best_clade_2_RR_id_$id
if [ -f "FinalResult_t5_NULL_MGPM_A_F_best_clade_2_RR_id_"$id".RData" ]
then
rm MPI*.log
rm CurrentResults*.RData
rm *worker*.RData
else
rm MPI*.log
rm CurrentResults*.RData
rm *worker*.RData
bsub -M 250000 -n 50 -W 23:59 -R "rusage[mem=4096]" -R ib sh R --vanilla --slave -f ../../DetectShifts_t5_NULL_MGPM_A_F_best_clade_2_RR.R --args $id
fi
cd ../..
done


for id in `seq 129 1 256`
do
mkdir -p Results_t5_NULL_MGPM_A_F_best_clade_2_RR_N159/t5_NULL_MGPM_A_F_best_clade_2_RR_id_$id
cd Results_t5_NULL_MGPM_A_F_best_clade_2_RR_N159/t5_NULL_MGPM_A_F_best_clade_2_RR_id_$id
if [ -f "FinalResult_t5_NULL_MGPM_A_F_best_clade_2_RR_id_"$id".RData" ]
then
rm MPI*.log
rm CurrentResults*.RData
rm *worker*.RData
else
rm MPI*.log
rm CurrentResults*.RData
rm *worker*.RData
bsub -M 250000 -n 40 -W 23:59 -R "rusage[mem=4096]" -R ib sh R --vanilla --slave -f ../../DetectShifts_t5_NULL_MGPM_A_F_best_clade_2_RR.R --args $id
fi
cd ../..
done

