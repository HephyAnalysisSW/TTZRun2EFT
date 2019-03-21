#!/bin/bash

#declare -a cards=("TTGamma_SingleLeptFromT_3LineBuggy" "TTGamma_SingleLeptFromTbar_3LineBuggy" "TTGamma_DiLept_3LineBuggy" )
#declare -a cards=("TTGamma_DiLept_3LineBuggy" )
#declare -a cards=("TTGamma_SingleLeptFromT_1Line_SM")
declare -a cards=("TTGamma_Atlas")
#declare -a cards=("TTGamma_SingleLeptFromT_3LineBuggy_test" "TTGamma_SingleLeptFromT_3LinePatched_test" "TTGamma_SingleLeptFromT_1Line_test" )

runCard="/afs/hephy.at/user/l/llechner/public/CMSSW_9_4_10/src/TTGammaEFT/Generation/data/runCards/TTGamma_Atlas.dat"
runTag="Atlas_paper"
#runCard="/afs/hephy.at/user/l/llechner/public/CMSSW_9_4_10/src/TTGammaEFT/Generation/data/runCards/TTGamma.dat"
#runTag="CMS_old"
model="dim6top_LO"
#model="sm-ckm_no_b_mass"

#################################################

for card in "${cards[@]}"
do

#   echo "run.py --calcXSec --runCard ${runCard} --process ${card} --model ${model} --overwrite > logs/${card}.xsec 2>&1" &
#    run.py --calcXSec --runCard ${runCard} --process ${card} --model ${model} --overwrite > logs/${card}.xsec 2>&1 &
   nohup krenew -t -K 10 -- bash -c "run.py --calcXSec --runCard ${runCard} --process ${card} --model ${model} --overwrite" > logs/${card}_${runTag}.xsec 2>&1 &

done


