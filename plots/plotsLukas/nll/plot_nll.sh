#ctZ

#Dileps
#python plot_nll.py $1 --years 2016 --variables ctZ ctZI --selections 2l --contours --smooth --zRange 0 20 --xyRange -1 1 -1 1 --binMultiplier 5
#python plot_nll.py $1 --years 2017 --variables ctZ ctZI --selections 2l --contours --smooth --zRange 0 20 --xyRange -1 1 -1 1 --binMultiplier 5
#python plot_nll.py $1 --years 2018 --variables ctZ ctZI --selections 2l --contours --smooth --zRange 0 20 --xyRange -1 1 -1 1 --binMultiplier 5
#SemiLep
#python plot_nll.py $1 --years 2016 --variables ctZ ctZI --selections 1l --contours --smooth --zRange 0 20 --xyRange -1 1 -1 1 --binMultiplier 5
#python plot_nll.py $1 --years 2017 --variables ctZ ctZI --selections 1l --contours --smooth --zRange 0 20 --xyRange -1 1 -1 1 --binMultiplier 5
#python plot_nll.py $1 --years 2018 --variables ctZ ctZI --selections 1l --contours --smooth --zRange 0 20 --xyRange -1 1 -1 1 --binMultiplier 5
#Combined
#python plot_nll.py $1 --years 2016 2017 2018 --variables ctZ ctZI --selections 2l --contours --smooth --zRange 0 20 --xyRange -1 1 -1 1 --binMultiplier 5
#SemiLep Combined
#python plot_nll.py $1 --years 2016 2017 2018 --variables ctZ ctZI --selections 1l --contours --smooth --zRange 0 20 --xyRange -1 1 -1 1 --binMultiplier 5
#All Combined
#python plot_nll.py $1 --years 2016 2017 2018 --variables ctZ ctZI --selections 1l 2l --contours --smooth --zRange 0 20 --xyRange -0.5 0.5 -0.5 0.5 --binMultiplier 5
# all years
#python plot_nll.py $1 --years 2016 --variables ctZ ctZI --selections 1l 2l --contours --smooth --zRange 0 20 --xyRange -1 1 -1 1 --binMultiplier 5
#python plot_nll.py $1 --years 2017 --variables ctZ ctZI --selections 1l 2l --contours --smooth --zRange 0 20 --xyRange -1 1 -1 1 --binMultiplier 5
#python plot_nll.py $1 --years 2018 --variables ctZ ctZI --selections 1l 2l --contours --smooth --zRange 0 20 --xyRange -1 1 -1 1 --binMultiplier 5
# Unc. studies
#python plot_nll.py $1 --tag high --years 2016 2017 2018 --variables ctZ ctZI --selections 1l 2l --contours --smooth --zRange 0 20 --xyRange -0.5 0.5 -0.5 0.5 --binMultiplier 5 --binning 50 -1 1 50 -1 1
#python plot_nll.py $1 --tag low  --years 2016 2017 2018 --variables ctZ ctZI --selections 1l 2l --contours --smooth --zRange 0 20 --xyRange -0.5 0.5 -0.5 0.5 --binMultiplier 5 --binning 50 -1 1 50 -1 1

#ctW

#Dilep
##python plot_nll.py $1 --years 2016 --variables ctW ctWI --selections 2l --contours --smooth --zRange 0 20 --xyRange -1 1 -1 1 --binMultiplier 5
##python plot_nll.py $1 --years 2017 --variables ctW ctWI --selections 2l --contours --smooth --zRange 0 20 --xyRange -1 1 -1 1 --binMultiplier 5
##python plot_nll.py $1 --years 2018 --variables ctW ctWI --selections 2l --contours --smooth --zRange 0 20 --xyRange -1 1 -1 1 --binMultiplier 5
#SemiLep
##python plot_nll.py $1 --years 2016 --variables ctW ctWI --selections 1l --contours --smooth --zRange 0 20 --xyRange -1 1 -1 1 --binMultiplier 5
##python plot_nll.py $1 --years 2017 --variables ctW ctWI --selections 1l --contours --smooth --zRange 0 20 --xyRange -1 1 -1 1 --binMultiplier 5
##python plot_nll.py $1 --years 2018 --variables ctW ctWI --selections 1l --contours --smooth --zRange 0 20 --xyRange -1 1 -1 1 --binMultiplier 5
#Combined
##python plot_nll.py $1 --years 2016 2017 2018 --variables ctW ctWI --selections 2l --contours --smooth --zRange 0 20 --xyRange -1 1 -1 1 --binMultiplier 5
#SemiLep Combined
##python plot_nll.py $1 --years 2016 2017 2018 --variables ctW ctWI --selections 1l --contours --smooth --zRange 0 20 --xyRange -1 1 -1 1 --binMultiplier 5
#All Combined
##python plot_nll.py $1 --years 2016 2017 2018 --variables ctW ctWI --selections 1l 2l --contours --smooth --zRange 0 20 --xyRange -1 1 -1 1 --binMultiplier 5


# 1D

#Dileps
python plot_nll1D.py $1 --years 2016 --variables ctZ ctZI --selections 2l --contours --smooth --zRange 0 5 --xyRange -0.7 0.7 -0.7 0.7 --binMultiplier 5
python plot_nll1D.py $1 --years 2017 --variables ctZ ctZI --selections 2l --contours --smooth --zRange 0 5 --xyRange -0.7 0.7 -0.7 0.7 --binMultiplier 5
python plot_nll1D.py $1 --years 2018 --variables ctZ ctZI --selections 2l --contours --smooth --zRange 0 5 --xyRange -0.7 0.7 -0.7 0.7 --binMultiplier 5
#SemiLep
python plot_nll1D.py $1 --years 2016 --variables ctZ ctZI --selections 1l --contours --smooth --zRange 0 5 --xyRange -0.5 0.5 -0.5 0.5 --binMultiplier 5
python plot_nll1D.py $1 --years 2017 --variables ctZ ctZI --selections 1l --contours --smooth --zRange 0 5 --xyRange -0.5 0.5 -0.5 0.5 --binMultiplier 5
python plot_nll1D.py $1 --years 2018 --variables ctZ ctZI --selections 1l --contours --smooth --zRange 0 5 --xyRange -0.5 0.5 -0.5 0.5 --binMultiplier 5
#Combined
python plot_nll1D.py $1 --years 2016 2017 2018 --variables ctZ ctZI --selections 2l --contours --smooth --zRange 0 5 --xyRange -0.6 0.6 -0.6 0.6 --binMultiplier 5
#SemiLep Combined
python plot_nll1D.py $1 --years 2016 2017 2018 --variables ctZ ctZI --selections 1l --contours --smooth --zRange 0 5 --xyRange -0.4 0.4 -0.4 0.4 --binMultiplier 5
#All Combined
python plot_nll1D.py $1 --years 2016 2017 2018 --variables ctZ ctZI --selections 1l 2l --contours --smooth --zRange 0 5 --xyRange -0.3 0.3 -0.3 0.3 --binMultiplier 5
#all years
python plot_nll1D.py $1 --years 2016 --variables ctZ ctZI --selections 1l 2l --contours --smooth --zRange 0 5 --xyRange -0.7 0.7 -0.7 0.7 --binMultiplier 5
python plot_nll1D.py $1 --years 2017 --variables ctZ ctZI --selections 1l 2l --contours --smooth --zRange 0 5 --xyRange -0.7 0.7 -0.7 0.7 --binMultiplier 5
python plot_nll1D.py $1 --years 2018 --variables ctZ ctZI --selections 1l 2l --contours --smooth --zRange 0 5 --xyRange -0.7 0.7 -0.7 0.7 --binMultiplier 5
#Unc studies
#python plot_nll1D.py $1 --tag high --years 2016 2017 2018 --variables ctZ ctZI --selections 1l 2l --contours --smooth --zRange 0 5 --xyRange -0.3 0.3 -0.3 0.3 --binMultiplier 5
#python plot_nll1D.py $1 --tag low  --years 2016 2017 2018 --variables ctZ ctZI --selections 1l 2l --contours --smooth --zRange 0 5 --xyRange -0.3 0.3 -0.3 0.3 --binMultiplier 5

#ctW

#Dilep
##python plot_nll1D.py $1 --years 2016 --variables ctW ctWI --selections 2l --contours --smooth --zRange 0 5 --xyRange -0.3 0.3 -0.7 0.7 --binMultiplier 5
##python plot_nll1D.py $1 --years 2017 --variables ctW ctWI --selections 2l --contours --smooth --zRange 0 5 --xyRange -0.3 0.3 -0.7 0.7 --binMultiplier 5
##python plot_nll1D.py $1 --years 2018 --variables ctW ctWI --selections 2l --contours --smooth --zRange 0 5 --xyRange -0.3 0.3 -0.7 0.7 --binMultiplier 5
#SemiLep
##python plot_nll1D.py $1 --years 2016 --variables ctW ctWI --selections 1l --contours --smooth --zRange 0 5 --xyRange -0.25 0.25 -0.5 0.5 --binMultiplier 5
##python plot_nll1D.py $1 --years 2017 --variables ctW ctWI --selections 1l --contours --smooth --zRange 0 5 --xyRange -0.25 0.25 -0.5 0.5 --binMultiplier 5
##python plot_nll1D.py $1 --years 2018 --variables ctW ctWI --selections 1l --contours --smooth --zRange 0 5 --xyRange -0.25 0.25 -0.5 0.5 --binMultiplier 5
#Combined
##python plot_nll1D.py $1 --years 2016 2017 2018 --variables ctW ctWI --selections 2l --contours --smooth --zRange 0 5 --xyRange -0.3 0.3 -0.7 0.7 --binMultiplier 5
#SemiLep Combined
##python plot_nll1D.py $1 --years 2016 2017 2018 --variables ctW ctWI --selections 1l --contours --smooth --zRange 0 5 --xyRange -0.25 0.25 -0.5 0.5 --binMultiplier 5
#All Combined
##python plot_nll1D.py $1 --years 2016 2017 2018 --variables ctW ctWI --selections 1l 2l --contours --smooth --zRange 0 5 --xyRange -0.25 0.25 -0.5 0.5 --binMultiplier 5
