#!/bin/sh

# ttZ Signal samples with reference point
#python genPostProcessing.py --overwrite all --targetDir TTZRun2EFT_PP_GEN_v1 --logLevel DEBUG --sample ttZ_ll_LO_order2_15weights_ref --addReweights --interpolationOrder 2 #SPLIT200
#python genPostProcessing.py --overwrite all --targetDir TTZRun2EFT_PP_GEN_v1 --logLevel DEBUG --sample ttZ_ll_LO_order4_6weights_ref --addReweights --interpolationOrder 2 #SPLIT200
python genPostProcessing.py --overwrite all --targetDir TTZRun2EFT_PP_GEN_v1 --logLevel DEBUG --sample ttZ_ll_LO_order3_8weights --addReweights --interpolationOrder 2 #SPLIT200
