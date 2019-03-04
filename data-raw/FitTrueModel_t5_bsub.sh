for id in `seq 1 1 2048`
do
mkdir -p Results_t5/FitTrueModel_t5_2_id_$id
cd Results_t5/FitTrueModel_t5_2_id_$id
if [ -f "FinalResult_FitTrueModel_t5_2_id_"$id".RData" ]
then
echo "Found result for "$id
else
bsub -n 1 -W 4:00 sh R --vanilla --slave -f ../../FitTrueModel_t5.R --args $id
fi
cd ../..
done
