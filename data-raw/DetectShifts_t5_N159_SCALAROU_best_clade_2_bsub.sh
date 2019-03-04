# choose simulations to fit
# Tree−type ultrametric / N=159 / 2 regimes  / Mapping 1. ED 513:544
# Tree−type ultrametric / N=159 / 2 regimes  / Mapping 4. AC 609:640
# Tree−type ultrametric / N=159 / 5 regimes  / Mapping 1. EECFC 641:672
# Tree−type ultrametric / N=159 / 5 regimes  / Mapping 3. DCFBC 705:736

# Tree−type non−ultrametric / N=159 / 2 regimes  / Mapping 1. AF 769:800
# Tree−type non−ultrametric / N=159 / 2 regimes  / Mapping 2. CF 801:832
# Tree−type non−ultrametric / N=159 / 5 regimes  / Mapping 1. FCEFC 897:928
# Tree−type non−ultrametric / N=159 / 5 regimes  / Mapping 4. ADFEE 993:1024

for id in `seq 513 1 544; seq 609 1 640; seq 769 1 800; seq 801 1 832; seq 641 1 672; seq 705 1 736; seq 897 1 928; seq 993 1 1024`
do
mkdir -p Results_t5_N159_SCALAROU/SCALAROU_best_clade_2_id_$id
cd Results_t5_N159_SCALAROU/SCALAROU_best_clade_2_id_$id
if [ -f "FinalResult_SCALAROU_best_clade_2_id_"$id".RData" ]
then
rm MPI*.log
rm CurrentResults*.RData
else
bsub -n 10 -W 23:59 -R ib sh R --vanilla --slave -f ../../DetectShifts_t5_SCALAROU_best_clade_2.R --args $id
fi
cd ../..
done
