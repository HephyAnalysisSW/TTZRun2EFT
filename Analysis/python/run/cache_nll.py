''' Plot script WC parameter LogLikelihood
'''

# Standard imports 
import sys, os, pickle, copy, ROOT, time
import numpy as np
import itertools

# Multiprocessing
from multiprocessing import Pool

# RootTools
from RootTools.core.standard   import *

# turn off graphics
ROOT.gROOT.SetBatch( True )

# TTZRun2EFT
from TTZRun2EFT.Tools.user              import combineReleaseLocation, cache_directory, cardfileLocation

from TTZRun2EFT.Tools.genCutInterpreter import cutInterpreter

# get the reweighting function
from Analysis.Tools.WeightInfo          import WeightInfo
from Analysis.Tools.CardFileWriter      import CardFileWriter
from Analysis.Tools.cardFileHelpers     import scaleCardFile, getRegionCuts

from TTZRun2EFT.Tools.Cache             import Cache
from TTZRun2EFT.Analysis.regions        import simpleStringToCutString

# Default Parameter
loggerChoices = ['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET']

# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',           action='store',      default='INFO', nargs='?', choices=loggerChoices,                                help="Log level for logging")
argParser.add_argument('--genSelection',       action='store',      default='onZll')
argParser.add_argument('--keepCards',          action='store_true',                                                                                  help='Keep all cardfiles?', )
argParser.add_argument('--small',              action='store_true',                                                                                  help='Run only on a small subset of the data?', )
argParser.add_argument('--overwrite',          action='store_true',                                                                                  help='Overwrite Database entries?', )
argParser.add_argument('--order',              action='store',      default=2, type=int,                                                             help='Polynomial order of weight string (e.g. 2)')
argParser.add_argument('--years',              action='store',      default=[ 2016, 2017 ], type=int, choices=[2016, 2017, 2018], nargs="*",         help="Which years to combine?")
argParser.add_argument('--cardFileSM16',       action='store',      default=None, type=str,                                                          help="SM cardfile for 2016")
argParser.add_argument('--cardFileSM17',       action='store',      default=None, type=str,                                                          help="SM cardfile for 2017")
argParser.add_argument('--cardFileSM18',       action='store',      default=None, type=str,                                                          help="SM cardfile for 2018")
argParser.add_argument('--variables' ,         action='store',      default = ['ctZ', 'ctZI'],       type=str,   nargs=2,                            help="argument plotting variables")
argParser.add_argument('--binning',            action='store',      default=[30, -2, 2, 30, -2, 2 ], type=float, nargs=6,                            help="argument parameters")
argParser.add_argument('--cores',              action='store',      default=1,                       type=int,                                       help='number of cpu cores for multicore processing')
argParser.add_argument('--checkOnly',          action='store_true',                                                                                  help='Just check if yield is already calculated', )
args = argParser.parse_args()

# Logger
import Analysis.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(    args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger( args.logLevel, logFile = None)

# Gen Samples
from TTZRun2EFT.Samples.genTuples_TTZ_postProcessed import *

combinedAnalysis = len(args.years) > 1
lumi16 = 35.92
lumi17 = 41.86
lumi18 = 58.83

if not args.cardFileSM16 and not args.cardFileSM17 and not args.cardFileSM18:
    raise Exception("What do you want me to do? Give me a cardfile to scale...")
elif (not args.cardFileSM16 and 2016 in args.years):
    raise Exception("Sorry, currently only scaling of at least 2016 cardfiles implemented")

# create missing SM cardfiles in sclaling previous years
if args.cardFileSM16 and (not args.cardFileSM17 and 2017 in args.years):
    if "2016" in args.cardFileSM16:
        args.cardFileSM17 = args.cardFileSM16.replace("2016","2017")
    else:
        args.cardFileSM17 = "2017_" + args.cardFileSM16
    scaleCardFile( args.cardFileSM16, args.cardFileSM17, scale=lumi17/lumi16, copyUncertainties=True)
    logger.info("Scaled 2016 cardfile to 2017")
if args.cardFileSM17 and (not args.cardFileSM18 and 2018 in args.years):
    # prefer a 2017->2018 scaling over 2016->2018 scaling if possible
    if "2017" in args.cardFileSM17:
        args.cardFileSM18 = args.cardFileSM17.replace("2017","2018")
    else:
        args.cardFileSM18 = "2018_" + args.cardFileSM17
    scaleCardFile( args.cardFileSM17, args.cardFileSM18, scale=lumi18/lumi17, copyUncertainties=True)
    logger.info("Scaled 2017 cardfile to 2018")

if args.cardFileSM16 and (not args.cardFileSM18 and 2018 in args.years):
    # only needed for [2016, 2018] combined analyses, which is stupid
    if "2016" in args.cardFileSM16:
        args.cardFileSM18 = args.cardFileSM16.replace("2016","2018")
    else:
        args.cardFileSM18 = "2018_" + args.cardFileSM16
    scaleCardFile( args.cardFileSM16, args.cardFileSM18, scale=lumi18/lumi16, copyUncertainties=True)
    logger.info("Scaled 2016 cardfile to 2018")

cardsSM    = { 2016:args.cardFileSM16 if args.cardFileSM16 else "none", 2017:args.cardFileSM17 if args.cardFileSM17 else "none", 2018:args.cardFileSM18 if args.cardFileSM18 else "none" }
regionCuts = { year:getRegionCuts( cardsSM[year] ) for year in args.years }

tableName = "nllcache"
cache_dir_nll    = os.path.join(cache_directory, "nll")
if not os.path.isdir(cache_dir_nll):
    os.makedirs(cache_dir_nll)
dbFile = "_".join( ["NLLcache"] + map( str, args.years ) ) + ".sql"
dbPath = os.path.join(cache_dir_nll, dbFile)

nllCache  = Cache( dbPath, tableName, ["cardname16", "cardname17", "cardname18", "cardname", "WC", "WC_val", "nll_prefit", "nll_postfit" ] )
if nllCache is None: raise

#signalSample = TTG_SingleLeptFromT_1L_test_EFT
genSignalSample = TTZ_EFT

cardname  = [ genSignalSample.name ]
cardname += map( str, args.years )
cardname += args.genSelection
for i, var in enumerate(args.variables):
    cardname += [ var, "var%i"%i ]
cardname += [ 'small' if args.small else 'full' ]
cardname  = '_'.join( cardname )

logger.info( "General card name: %s" %cardname )

def isInDatabase( point ):
    card = cardname
    for i, var in enumerate(point):
        card = card.replace("var%i", str(var))
    res  = {"cardname16":cardsSM[2016], "cardname17":cardsSM[2017], "cardname18":cardsSM[2018], "cardname":card, "WC":"_".join(args.variables), "WC_val":"_".join(map(str,point))}
    nCacheFiles = nllCache.contains( res )
    return bool(nCacheFiles)

def chunks( l, n ):
    # split list in chuncs
    for i in range( 0, len(l), n ):
        yield l[i:i + n]

binDict = dict(zip( args.variables, [dict(zip( ["nBins", "min", "max"], var )) for var in chunks(args.binning, 3) ] ))

for var in binDict.values():
    if var["nBins"] > 1:
        xRange       = np.linspace( var["min"], var["max"], var["nBins"], endpoint=False)
        halfstepsize = 0.5 * ( xRange[1] - xRange[0] )
        var["range"] = [0] + [ round(el + halfstepsize, 3) for el in xRange ]
    else:
        var["range"] = [0] + [ 0.5 * ( var["min"] + var["max"] ) ]

# get all possible combinations of points for i variables
points2D = list( itertools.product( *(binDict[var]["range"] for var in args.variables) ) ) #2D plots

# check if already in database before calculations
if not args.overwrite:
    logger.info( "Checking %i points" %(len(points2D)) )
    points2D = [ point for point in points2D if not isInDatabase( point ) ]

    if args.checkOnly:
        print "Limit not calculated for %i points" %len(points2D)
        logger.info( "Limit not calculated for %i points" %(len(points2D)) )
        sys.exit(0)

    if not points2D: sys.exit(0) #nothing to calculate

logger.info( "%i points to calculate" %(len(points2D)) )

w = WeightInfo( genSignalSample.reweight_pkl )
w.set_order( args.order )

coeffList = {}
signal_genRateSM = {}
for year in args.years:
    coeffList[year] = {}
    signal_genRateSM[year] = {}

def calculation( point ):
    
    kwargs = dict( zip( args.variables, point ) )

    cacheCard = cardname
    for i, var in enumerate( point ):
        cacheCard = cacheCard.replace("var%i", str(var))

    card     = {}
    cardpath = {}
    #caching results
    for year in args.years:
        logger.info( "Creating cardfile for year %i, %s"%( year, ", ".join( [ str(key) + " = " + str(val) for key, val in kwargs.iteritems() ] ) ) )

        card[year] = cacheCard.replace( "_".join(map(str,args.years)), str(year) )
        cardpath[year] = os.path.join( cardfileLocation, card[year] + '.txt' )

        scale = {}
        for region in regionCuts[year].keys():

            logger.info( "At EFT region %s", region )

            signal_genRateEFT = w.get_weight_yield( coeffList[year][region], **kwargs )
            scale[region]     = signal_genRateEFT / signal_genRateSM[year][region]

        scaleCardFile( cardsSM[year], cardpath[year], scale=scale, scaledProcesses="signal", copyUncertainties=True)

    if combinedAnalysis:
        # combinin years, selections or both
        combinedCard = c.combineCards( cardpath )
        logger.info( "Using combined card: %s"%combinedCard )
    else:
        combinedCard = cardpath[year]
        logger.info( "Using card: %s"%combinedCard )

    nll          = c.calcNLL( combinedCard )
    nll_prefit   = nll['nll0']
    nll_postfit  = nll['nll_abs']
    
    if nll_prefit  is None or abs(nll_prefit) > 10000 or abs(nll_prefit) < 1e-5:   nll_prefit  = 999
    if nll_postfit is None or abs(nll_postfit) > 10000 or abs(nll_postfit) < 1e-5: nll_postfit = 999

    if not args.keepCards and not (var1==0 and var2==0):
        os.remove( combinedCard )
        if combinedAnalysis:
            for sel in args.selections:
                for year in args.years:
                    year = str(year)+str(sel)
                    os.remove( cardpath[year] )

    res  = {"cardname16":cardsSM[2016], "cardname17":cardsSM[2017], "cardname18":cardsSM[2018], "cardname":cacheCard, "WC":"_".join(args.variables), "WC_val":"_".join(point), "nll_prefit":nll_prefit, "nll_postfit":nll_postfit}
    logger.info( "NLL limit for %s: nll_prefit = %f, nll_postfit = %f"%( ", ".join( [" = ".join( [var, point[i]] ) for i,var in enumerate(args.variables) ] ), nll_prefit, nll_postfit) )
    nllCache.add( res, nll_prefit, overwrite=True )

    del c

def setup():
    # preparing gen-sample reweighting
    logger.info( "Preparing reweighting setup" )

    for year in args.years:
        for region, cut in regionCuts[year].iteritems():

            logger.info( "At region %s", region )
            print cut
            print simpleStringToCutString( cut )

            sel = "&&".join( [ cutInterpreter.cutString( args.genSelection ), simpleStringToCutString( cut ) ] )
            print sel
            coeffList[year][region]        = w.getCoeffListFromDraw( genSignalSample, selectionString=sel )
            signal_genRateSM[year][region] = float( w.get_weight_yield( coeffList[year][region] ) )

# run initial setup
setup()

# run calculations
pool = Pool( processes=args.cores )
pool.map( calculation, points2D )
pool.close()

