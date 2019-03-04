# choose simulations to fit
# Tree−type ultrametric / N=318 / 2 regimes  / Mapping 3. DC 1089:1120
# Tree−type ultrametric / N=318 / 2 regimes  / Mapping 4. BF 1121:1152
# Tree−type ultrametric / N=318 / 8 regimes  / Mapping 1. DBACFDFE 1153:1184
# Tree−type ultrametric / N=318 / 8 regimes  / Mapping 2. CCAAEACD 1185:1216

# Tree−type non−ultrametric / N=318 / 2 regimes  / Mapping 1. DD 1281:1312
# Tree−type non−ultrametric / N=318 / 2 regimes  / Mapping 3. ED 1345:1376
# Tree−type non−ultrametric / N=318 / 8 regimes  / Mapping 1. ECBEAFDD 1409:1440
# Tree−type non−ultrametric / N=318 / 8 regimes  / Mapping 3. BFCEFCAC 1473:1504

#for id in `seq 1537 4 1568; seq 1569 4 1600; seq 1825 4 1856; seq 1889 4 1920; seq 1665 4 1696; seq 1697 4 1728; seq 1921 4 1952; seq 2017 4 2048`
for id in `seq 1538 4 1568; seq 1570 4 1600; seq 1826 4 1856; seq 1890 4 1920; seq 1666 4 1696; seq 1698 4 1728; seq 1922 4 1952; seq 2018 4 2048`
do
mkdir -p Results_t5_N638_SCALAROU_24h/SCALAROU_best_clade_2_id_$id
cd Results_t5_N638_SCALAROU_24h/SCALAROU_best_clade_2_id_$id
if [ -f "FinalResult_SCALAROU_best_clade_2_id_"$id".RData" ]
then
rm MPI*.log
rm CurrentResults*.RData
else
bsub -n 10 -W 23:59 -R ib sh R --vanilla --slave -f ../../DetectShifts_t5_SCALAROU_best_clade_2_24h.R --args $id
fi
cd ../..
done
