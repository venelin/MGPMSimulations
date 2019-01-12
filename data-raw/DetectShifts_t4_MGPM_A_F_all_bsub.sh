for id in `seq 129 1 130`
do
mkdir -p Results_t4/MGPM_A_F_all_$id
cd Results_t4/MGPM_A_F_all_$id
bsub -M 200000 -n 200 -W 23:59 -R ib sh R --vanilla --slave -f ../../DetectShifts_t4_MGPM_A_F_all.R --args $id
cd ../..
done
