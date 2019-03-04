# choose simulations to fit
# Tree−type ultrametric / N=80 / 2 regimes  / Mapping 1. BC 1:32
# Tree−type ultrametric / N=80 / 2 regimes  / Mapping 2. DF 33:64
# Tree−type non−ultrametric / N=80 / 2 regimes  / Mapping 1. FB 257:288
# Tree−type non−ultrametric / N=80 / 2 regimes  / Mapping 2. DA 289:320

# Tree−type ultrametric / N=80 / 3 regimes  / Mapping 1. DAB 129:160
# Tree−type ultrametric / N=80 / 3 regimes  / Mapping 2. BEC 161:192
# Tree−type non−ultrametric / N=80 / 3 regimes  / Mapping 1. FCC 385:416
# Tree−type non−ultrametric / N=80 / 3 regimes  / Mapping 2. DCB 417:448

for id in `seq 1 1 32; seq 33 1 64; seq 257 1 288; seq 289 1 320; seq 129 1 160; seq 161 1 192; seq 385 1 416; seq 417 1 448`
do
mkdir -p Results_t5_N80_SURFACE_mcs10/SURFACE_best_clade_2_mcs10_id_$id
cd Results_t5_N80_SURFACE_mcs10/SURFACE_best_clade_2_mcs10_id_$id
if [ -f "FinalResult_SURFACE_best_clade_2_mcs10_id_"$id".RData" ]
then
rm MPI*.log
rm CurrentResults*.RData
else
bsub -n 10 -W 23:59 -R ib sh R --vanilla --slave -f ../../DetectShifts_t5_SURFACE_best_clade_2_mcs10.R --args $id
fi
cd ../..
done
