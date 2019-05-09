for id in `seq 385 1 512`
do
mkdir -p Results_t5_NULL_MGPM_A_F_best_clade_2_RR_N638/t5_NULL_MGPM_A_F_best_clade_2_RR_id_$id
cd Results_t5_NULL_MGPM_A_F_best_clade_2_RR_N638/t5_NULL_MGPM_A_F_best_clade_2_RR_id_$id
if [ -f "FinalResult_t5_NULL_MGPM_A_F_best_clade_2_RR_id_"$id".RData" ]
then
rm MPI*.log
rm CurrentResults*.RData
else
bsub -M 400000 -n 200 -W 23:59 -R "rusage[mem=4096]" -R ib sh R --vanilla --slave -f ../../DetectShifts_t5_NULL_MGPM_A_F_best_clade_2_RR.R --args $id
fi
cd ../..
done
