for id in `seq 1 1 2048`
do
mkdir -p Results_FitTrueModel_NoHint_t5/FitTrueModel_NoHint_t5_id_$id
cd Results_FitTrueModel_NoHint_t5/FitTrueModel_NoHint_t5_id_$id
if [ -f "FinalResult_FitTrueModel_NoHint_t5_id_"$id".RData" ]
then
echo "Found result for "$id
else
bsub -n 1 -W 4:00 sh R --vanilla --slave -f ../../FitTrueModel_NoHint_t5.R --args $id
fi
cd ../..
done
