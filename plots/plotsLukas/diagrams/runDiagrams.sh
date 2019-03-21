#!/bin/bash

#declare -a cards=("TTGamma_SingleLeptFromT_3LineBuggy_test" "TTGamma_SingleLeptFromT_3LinePatched_test" "TTGamma_SingleLeptFromT_1Line_test" )
#declare -a cards=("TTGamma_SingleLeptFromT_3LineBuggy_full" "TTGamma_SingleLeptFromT_3LinePatched_full" "TTGamma_SingleLeptFromT_1Line_full" )
#declare -a cards=("TTGamma_SingleLeptFromT_1Line_full" "TTGamma_SingleLeptFromT_1Line_test")
#declare -a cards=("TTGamma_SingleLeptFromT_1Line")
#declare -a cards=('TTGamma_DiLept_1Line_EFT' 'TTGamma_DiLept_1Line_EFT_test' 'TTGamma_Had_1Line_EFT' 'TTGamma_Had_1Line_EFT_test' 'TTGamma_SingleLeptFromT_1Line_EFT' 'TTGamma_SingleLeptFromT_1Line_EFT_test' 'TTGamma_SingleLeptFromTbar_1Line_EFT' 'TTGamma_SingleLeptFromTbar_1Line_EFT_test')
#declare -a cards=('TTGamma_SingleLeptFromT_1Line' 'TTGamma_SingleLeptFromT_1Line_test' 'TTGamma_SingleLeptFromT_3LineBuggy' 'TTGamma_SingleLeptFromT_3LineBuggy_test' 'TTGamma_SingleLeptFromT_3LinePatched' 'TTGamma_SingleLeptFromT_3LinePatched_test')
#declare -a cards=('TTGamma_DiLept_1Line_EFT_test' 'TTGamma_Had_1Line_EFT_test' 'TTGamma_SingleLeptFromT_1Line_EFT_test' 'TTGamma_SingleLeptFromTbar_1Line_EFT_test')
declare -a cards=('TTGamma_DiLept_1Line_EFT' 'TTGamma_Had_1Line_EFT' 'TTGamma_SingleLeptFromT_1Line_EFT' 'TTGamma_SingleLeptFromTbar_1Line_EFT')

declare -a WC=("ctZ 2" "ctW 2" "ctG 2")
#declare -a WC=()

model="dim6top_LO"

#################################################

for card in "${cards[@]}"
do

    nohup krenew -t -K 10 -- bash -c "getDiagrams.py --process ${card} --model ${model} --overwrite" &

   for param in "${WC[@]}"
   do

    nohup krenew -t -K 10 -- bash -c "getDiagrams.py --process ${card} --model ${model} --couplings ${param} --overwrite" &

   done

done


