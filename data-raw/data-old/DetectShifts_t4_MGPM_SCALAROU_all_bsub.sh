for id in `seq 1 1 256`
do
mkdir -p Results_t4/MGPM_SCALAROU_all_$id
cd Results_t4/MGPM_SCALAROU_all_$id
if [ -f "FinalResult_MGPM_SCALAROU_all_1_id_"$id"_.RData" ]
then
rm MPI*.log
rm CurrentResults*.RData
else
bsub -M 100000 -n 100 -W 23:59 -R ib sh R --vanilla --slave -f ../../DetectShifts_t4_MGPM_SCALAROU_all.R --args $id
fi
cd ../..
done
