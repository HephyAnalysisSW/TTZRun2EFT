#ctZ

card16="inputcards/regionsE_2016_xsec_shape_lowUnc_expected_SRandCR_ewkDM_currents_ewkDM_ttZ_ll.txt"
card17="inputcards/regionsE_2017_xsec_shape_lowUnc_expected_SRandCR_ewkDM_currents_ewkDM_ttZ_ll.txt"

# TOP-18-009 caching
python cache_nll.py --cardFileSM16 $card16 --cardFileSM17 $card17 --years 2016 2017 --binning 30 -2  2  30 -2  2  --variables ctZ ctZI --cores 3
python cache_nll.py --cardFileSM16 $card16 --cardFileSM17 $card17 --years 2016 2017 --binning 30 -12 36 30 -28 18 --variables cpQM cpt --cores 3

# full run2 caching
#python cache_nll.py --cardFileSM16 $card16 --cardFileSM17 $card17 --years 2016 2017 2018 --binning 30 -2  2  30 -2  2  --variables ctZ ctZI --cores 6
#python cache_nll.py --cardFileSM16 $card16 --cardFileSM17 $card17 --years 2016 2017 2018 --binning 30 -12 36 30 -28 18 --variables cpQM cpt --cores 6

# 4D NLL caching
#python cache_nll.py --cardFileSM16 $card16 --cardFileSM17 $card17 --years 2016 2017 2018 --binning 30 -2 2 30 -2 2 30 -12 36 30 -28 18 --variables ctZ ctZI cpQM cpt --cores 6
