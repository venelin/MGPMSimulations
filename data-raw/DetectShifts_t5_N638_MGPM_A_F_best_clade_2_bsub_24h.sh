# choose simulations to fit
# Tree−type ultrametric / N=638 / 2 regimes  / Mapping 1. DE            1537:1568
# Tree−type ultrametric / N=638 / 2 regimes  / Mapping 2. DF            1569:1600
# Tree−type ultrametric / N=638 / 8 regimes  / Mapping 1. FBEEFDEC      1665:1696
# Tree−type ultrametric / N=638 / 8 regimes  / Mapping 2. AFBDAFBE      1697:1728

# Tree−type non−ultrametric / N=638 / 2 regimes  / Mapping 2. FC        1825:1856
# Tree−type non−ultrametric / N=638 / 2 regimes  / Mapping 4. BD        1889:1920
# Tree−type non−ultrametric / N=638 / 8 regimes  / Mapping 1. FDBACFCA  1921:1952
# Tree−type non−ultrametric / N=638 / 8 regimes  / Mapping 4. CFEFCDCA  2017:2048

#for id in `seq 1537 1 1568; seq 1569 1 1600; seq 1825 1 1856; seq 1889 1 1920; seq 1665 1 1696; seq 1697 1 1728; seq 1921 1 1952; seq 2017 1 2048`
for id in `seq 1538 4 1568; seq 1570 4 1600; seq 1826 4 1856; seq 1890 4 1920; seq 1666 4 1696; seq 1698 4 1728; seq 1922 4 1952; seq 2018 4 2048`
do
mkdir -p Results_t5_MGPM_A_F_best_clade_2_N638_24h/MGPM_A_F_best_clade_2_id_$id
cd Results_t5_MGPM_A_F_best_clade_2_N638_24h/MGPM_A_F_best_clade_2_id_$id
if [ -f "FinalResult_MGPM_A_F_best_clade_2_id_"$id".RData" ]
then
rm MPI*.log
rm CurrentResults*.RData
else
bsub -M 200000 -n 200 -W 23:59 -R ib sh R --vanilla --slave -f ../../DetectShifts_t5_MGPM_A_F_best_clade_2_24h.R --args $id
fi
cd ../..
done
