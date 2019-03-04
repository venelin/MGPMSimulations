for id in `seq 1 1 128`
do
mkdir -p Results_t4/MGPM_A_F_all_$id
cd Results_t4/MGPM_A_F_all_$id
if [ -f "FinalResult_MGPM_A_F_all_1_id_"$id"_.RData" ]
then
rm MPI*.log
rm CurrentResults*.RData
else
bsub -M 200000 -n 200 -W 23:59 -R ib sh R --vanilla --slave -f ../../DetectShifts_t4_MGPM_A_F_all.R --args $id
fi
cd ../..
done
