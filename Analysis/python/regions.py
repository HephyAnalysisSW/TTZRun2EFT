from TTZRun2EFT.Analysis.Region import Region
from TTZRun2EFT.Analysis.Region import texString
from TTZRun2EFT.Analysis.Region import allowedVars

from math import pi

def getRegionsFromThresholds(var, vals, gtLastThreshold = True):
    return [Region(var, (vals[i], vals[i+1])) for i in range(len(vals)-1)]

def getRegions2D(varOne, varOneThresholds, varTwo, varTwoThresholds):
    regions_varOne  = getRegionsFromThresholds(varOne,  varOneThresholds)
    regions_varTwo  = getRegionsFromThresholds(varTwo, varTwoThresholds)

    regions2D = []
    for r1 in regions_varOne:
        for r2 in regions_varTwo:
            regions2D.append(r1+r2)

    return regions2D

def simpleStringToDict( simpleString ):

    # replace variables by a string not containing "_"
    for i, var in enumerate(allowedVars):
        simpleString = simpleString.replace(var, "var%i"%i)
    cutList = simpleString.split("_")

    # convert simpleString to threshold tuple, fill in dict
    cutDict = {}
    for cut in cutList:
        for i, var in enumerate(allowedVars):
            if "var"+str(i) in cut:
                cutRange = cut.replace("var%i"%i, "")
                cutRange = cutRange.split("To")
                cutRange = tuple( map( float, cutRange ) )
                if len(cutRange) == 1: cutRange = ( cutRange[0], -1 )
                cutDict.update( {var:cutRange} )

    return cutDict

def dictToCutString( dict ):

    res=[]
    for var in dict.keys():
        svar = var
        s1=svar+">="+str(dict[var][0])
        if dict[var][1]>-1: s1+="&&"+svar+"<"+str(dict[var][1])
        res.append(s1)
    return "&&".join(res)

def simpleStringToCutString( cutString ):
    print cutString
    print simpleStringToDict( cutString )
    print dictToCutString( simpleStringToDict( cutString ) )
    return dictToCutString( simpleStringToDict( cutString ) )

#Put all sets of regions that are used in the analysis, closure, tables, etc.

#differencial
thresholds = [ 20, 120, 220, 320, 420, -999 ]
genTTZRegions  = getRegionsFromThresholds( "GenPhoton_pt[0]", thresholds )
