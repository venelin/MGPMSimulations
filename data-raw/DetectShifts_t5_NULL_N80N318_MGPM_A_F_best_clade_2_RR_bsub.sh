for id in 1 2 5 6 9 10 13 14 17 18 21 22 25 26 29 30 33 34 37 38 41 42 45 46 49 50 53 54 57 58 61 62 65 66 69 70 73 74 77 78 81 82 85 86 89 90 93 94 97 98 101 102 105 106 109 110 113 114 117 118 121 122 125 126
do
mkdir -p Results_t5_NULL_MGPM_A_F_best_clade_2_RR_N80_2/t5_NULL_MGPM_A_F_best_clade_2_RR_id_$id
cd Results_t5_NULL_MGPM_A_F_best_clade_2_RR_N80_2/t5_NULL_MGPM_A_F_best_clade_2_RR_id_$id
if [ -f "FinalResult_t5_NULL_MGPM_A_F_best_clade_2_RR_id_"$id".RData" ]
then
rm MPI*.log
rm CurrentResults*.RData
rm *worker*.RData
else
rm MPI*.log
bsub -M 250000 -n 40 -W 23:59 -R "rusage[mem=4096]" -R ib sh R --vanilla --slave -f ../../DetectShifts_t5_NULL_MGPM_A_F_best_clade_2_RR.R --args $id
fi
cd ../..
done

for id in 257 258 261 262 265 266 269 270 273 274 277 278 281 282 285 286 289 290 293 294 297 298 301 302 305 306 309 310 313 314 317 318 321 322 325 326 329 330 333 334 337 338 341 342 345 346 349 350 353 354 357 358 361 362 365 366 369 370 373 374 377 378 381 382
do
mkdir -p Results_t5_NULL_MGPM_A_F_best_clade_2_RR_N318_2/t5_NULL_MGPM_A_F_best_clade_2_RR_id_$id
cd Results_t5_NULL_MGPM_A_F_best_clade_2_RR_N318_2/t5_NULL_MGPM_A_F_best_clade_2_RR_id_$id
if [ -f "FinalResult_t5_NULL_MGPM_A_F_best_clade_2_RR_id_"$id".RData" ]
then
rm MPI*.log
rm CurrentResults*.RData
rm *worker*.RData
else
rm MPI*.log
bsub -M 250000 -n 80 -W 23:59 -R "rusage[mem=4096]" -R ib sh R --vanilla --slave -f ../../DetectShifts_t5_NULL_MGPM_A_F_best_clade_2_RR.R --args $id
fi
cd ../..
done
