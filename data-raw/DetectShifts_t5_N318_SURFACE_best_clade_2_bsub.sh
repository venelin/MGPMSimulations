# choose simulations to fit
# Tree−type ultrametric / N=318 / 2 regimes  / Mapping 3. DC 1089:1120
# Tree−type ultrametric / N=318 / 2 regimes  / Mapping 4. BF 1121:1152
# Tree−type ultrametric / N=318 / 8 regimes  / Mapping 1. DBACFDFE 1153:1184
# Tree−type ultrametric / N=318 / 8 regimes  / Mapping 2. CCAAEACD 1185:1216

# Tree−type non−ultrametric / N=318 / 2 regimes  / Mapping 1. DD 1281:1312
# Tree−type non−ultrametric / N=318 / 2 regimes  / Mapping 3. ED 1345:1376
# Tree−type non−ultrametric / N=318 / 8 regimes  / Mapping 1. ECBEAFDD 1409:1440
# Tree−type non−ultrametric / N=318 / 8 regimes  / Mapping 3. BFCEFCAC 1473:1504

for id in `seq 1089 1 1120; seq 1121 1 1152; seq 1281 1 1312; seq 1345 1 1376; seq 1153 1 1184; seq 1185 1 1216; seq 1409 1 1440; seq 1473 1 1504`
do
mkdir -p Results_t5_N318_SURFACE/SURFACE_best_clade_2_id_$id
cd Results_t5_N318_SURFACE/SURFACE_best_clade_2_id_$id
if [ -f "FinalResult_SURFACE_best_clade_2_id_"$id".RData" ]
then
rm MPI*.log
rm CurrentResults*.RData
else
bsub -n 10 -W 23:59 -R ib sh R --vanilla --slave -f ../../DetectShifts_t5_SURFACE_best_clade_2.R --args $id
fi
cd ../..
done
