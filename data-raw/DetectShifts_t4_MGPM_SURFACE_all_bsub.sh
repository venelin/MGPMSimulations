for id in `seq 1 1 256`
do
mkdir -p Results_t4/MGPM_SURFACE_all_$id
cd Results_t4/MGPM_SURFACE_all_$id
bsub -M 200000 -n 200 -W 23:59 -R ib sh R --vanilla --slave -f ../../DetectShifts_t4_MGPM_SURFACE_all.R --args $id
cd ../..
done
