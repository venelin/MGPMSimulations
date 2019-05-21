for id in `seq 1 1 2048`
do
mkdir -p Results_FitModelFWithTrueShiftsOnPRCData_t5_3/FitModelFWithTrueShiftsOnPRCData_t5_id_$id
cd Results_FitModelFWithTrueShiftsOnPRCData_t5_3/FitModelFWithTrueShiftsOnPRCData_t5_id_$id
if [ -f "FinalResult_FitModelFWithTrueShiftsOnPRCData_t5_id_"$id".RData" ]
then
echo "Found result for "$id
else
bsub -n 1 -W 23:59 sh R --vanilla --slave -f ../../FitModelFWithTrueShiftsOnPRCData_t5.R --args $id
fi
cd ../..
done
