sleep 10
# Data
#python regions.py         --year 2016 --parameters None
#python regions.py         --year 2016 --parameters ctZ 2
#python regions.py         --year 2016 --parameters ctZI 2 ctWI 2 ctZ 2 ctW 2

#python regions_semilep.py --year 2016 --parameters None
#python regions_semilep.py --year 2016 --parameters ctZ 2
#python regions_semilep.py --year 2016 --parameters ctZI 2 ctWI 2 ctZ 2 ctW 2

# Signal Regions
#python regions.py --noData --year 2016 --parameters None
#python regions.py --noData --year 2017 --parameters None
#python regions.py --noData --year 2018 --parameters None

python regions.py --noData --year 2016 --parameters ctZ 2 ctZI 2
python regions.py --noData --year 2017 --parameters ctZ 2 ctZI 2
python regions.py --noData --year 2018 --parameters ctZ 2 ctZI 2

#python regions.py --noData --year 2016 --parameters ctZI 2 ctWI 2 ctZ 2 ctW 2
#python regions.py --noData --year 2017 --parameters ctZI 2 ctWI 2 ctZ 2 ctW 2
#python regions.py --noData --year 2018 --parameters ctZI 2 ctWI 2 ctZ 2 ctW 2

# Signal Regions
#python regions.py --noData --year 2016 --parameters ctZ 2 --normalize
#python regions.py --noData --year 2017 --parameters ctZ 2 --normalize
#python regions.py --noData --year 2018 --parameters ctZ 2 --normalize

#python regions.py --noData --year 2016 --parameters ctZI 2 ctWI 2 ctZ 2 ctW 2 --normalize
#python regions.py --noData --year 2017 --parameters ctZI 2 ctWI 2 ctZ 2 ctW 2 --normalize
#python regions.py --noData --year 2018 --parameters ctZI 2 ctWI 2 ctZ 2 ctW 2 --normalize

# SemiLep
#python regions_semilep.py --year 2016 --parameters None --noData
#python regions_semilep.py --year 2017 --parameters None --noData
#python regions_semilep.py --year 2018 --parameters None --noData

python regions_semilep.py --year 2016 --parameters ctZ 2 ctZI 2 --noData
python regions_semilep.py --year 2017 --parameters ctZ 2 ctZI 2 --noData
python regions_semilep.py --year 2018 --parameters ctZ 2 ctZI 2 --noData

#python regions_semilep.py --year 2016 --parameters ctZI 2 ctWI 2 ctZ 2 ctW 2 --noData
#python regions_semilep.py --year 2017 --parameters ctZI 2 ctWI 2 ctZ 2 ctW 2 --noData
#python regions_semilep.py --year 2018 --parameters ctZI 2 ctWI 2 ctZ 2 ctW 2 --noData

#normalized
#python regions_semilep.py --year 2016 --parameters ctZ 2 --noData --normalize
#python regions_semilep.py --year 2017 --parameters ctZ 2 --noData --normalize
#python regions_semilep.py --year 2018 --parameters ctZ 2 --noData --normalize

#python regions_semilep.py --year 2016 --parameters ctZI 2 ctWI 2 ctZ 2 ctW 2 --noData --normalize
#python regions_semilep.py --year 2017 --parameters ctZI 2 ctWI 2 ctZ 2 ctW 2 --noData --normalize
#python regions_semilep.py --year 2018 --parameters ctZI 2 ctWI 2 ctZ 2 ctW 2 --noData --normalize

