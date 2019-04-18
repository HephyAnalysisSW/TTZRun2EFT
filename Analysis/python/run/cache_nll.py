''' Plot script WC parameter LogLikelihood
'''

# Standard imports 
import sys, os, pickle, copy, ROOT, time
import numpy as np
import itertools, uuid
from shutil import copyfile, rmtree

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
from TTZRun2EFT.Analysis.Region         import aliases

# Default Parameter
loggerChoices = ['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET']

# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',           action='store',      default='INFO', nargs='?', choices=loggerChoices,                                help="Log level for logging")
argParser.add_argument('--genSelection',       action='store',      default='onZll-nJet3p') #should remove nJet requirement, but have to be equal to TOP-18-009
argParser.add_argument('--keepCards',          action='store_true',                                                                                  help='Keep all cardfiles?', )
argParser.add_argument('--small',              action='store_true',                                                                                  help='Run only on a small subset of the data?', )
argParser.add_argument('--overwrite',          action='store_true',                                                                                  help='Overwrite Database entries?', )
argParser.add_argument('--order',              action='store',      default=2, type=int,                                                             help='Polynomial order of weight string (e.g. 2)')
argParser.add_argument('--years',              action='store',      default=[ 2016, 2017 ], type=int, choices=[2016, 2017, 2018], nargs="*",         help="Which years to combine?")
argParser.add_argument('--cardFileSM16',       action='store',      default=None, type=str,                                                          help="SM cardfile for 2016")
argParser.add_argument('--cardFileSM17',       action='store',      default=None, type=str,                                                          help="SM cardfile for 2017")
argParser.add_argument('--cardFileSM18',       action='store',      default=None, type=str,                                                          help="SM cardfile for 2018")
argParser.add_argument('--variables' ,         action='store',      default = ['ctZ', 'ctZI'],       type=str,   nargs="*",                          help="argument plotting variables")
argParser.add_argument('--binning',            action='store',      default=[30, -2, 2, 30, -2, 2 ], type=float, nargs="*",                          help="argument parameters")
argParser.add_argument('--cores',              action='store',      default=1,                       type=int,                                       help='number of cpu cores for multicore processing')
argParser.add_argument('--mergeAndCheck',      action='store_true',                                                                                  help='Merge databases and check if yield is already calculated', )
argParser.add_argument('--nJobs',              action='store',      nargs='?', type=int, default=1,                                                  help="Maximum number of simultaneous jobs.")
argParser.add_argument('--job',                action='store',      nargs='?', type=int, default=0,                                                  help="Run only job i")
argParser.add_argument('--forceWriting',       action='store_true',                                                                                  help='Retry writing to DB with recovering file, DANGEROUS', )
args = argParser.parse_args()

# Logger
import Analysis.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(    args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger( args.logLevel, logFile = None)

# Gen Samples
logger.info("Loading gen-samples")
from TTZRun2EFT.Samples.genTuples_TTZ_postProcessed import *

combinedAnalysis = len(args.years) > 1
lumi16 = 35.92
lumi17 = 41.86
lumi18 = 58.83

if args.nJobs != 1:
    logger.info("Running batch mode! Setting args.cores to 1")
    args.cores = 1

if len(args.binning) != 3*len(args.variables):
    raise Exception("Binning arguments must be 3 times the number of variables (dimensions)")

if not args.cardFileSM16 and not args.cardFileSM17 and not args.cardFileSM18:
    raise Exception("What do you want me to do? Give me a cardfile to scale...")
elif (not args.cardFileSM16 and 2016 in args.years):
    raise Exception("Sorry, currently only scaling of at least 2016 cardfiles implemented")
regionCuts = { 2016:getRegionCuts( args.cardFileSM16 ) }
# create missing SM cardfiles in sclaling previous years
if args.cardFileSM16 and (not args.cardFileSM17 and 2017 in args.years):
    if "2016" in args.cardFileSM16:
        args.cardFileSM17 = args.cardFileSM16.replace("2016","scaled_2017")
    else:
        args.cardFileSM17 = args.cardFileSM16.replace(".","_scaled_2017.")
    if not args.mergeAndCheck:
        scaleCardFile( args.cardFileSM16, args.cardFileSM17, scale=lumi17/lumi16, copyUncertainties=True, keepObservation=False )
        logger.info("Scaled 2016 cardfile to 2017")
        regionCuts[2017] = getRegionCuts( args.cardFileSM16 )

if args.cardFileSM17 and (not args.cardFileSM18 and 2018 in args.years):
    # prefer a 2017->2018 scaling over 2016->2018 scaling if possible
    if "2017" in args.cardFileSM17:
        args.cardFileSM18 = args.cardFileSM17.replace("2017","scaled_2018")
    else:
        args.cardFileSM18 = args.cardFileSM17.replace(".","_scaled_2018.")
    if not args.mergeAndCheck:
        scaleCardFile( args.cardFileSM17, args.cardFileSM18, scale=lumi18/lumi17, copyUncertainties=True, keepObservation=False )
        logger.info("Scaled 2017 cardfile to 2018")
        regionCuts[2018] = getRegionCuts( args.cardFileSM17 )

if args.cardFileSM16 and (not args.cardFileSM18 and 2018 in args.years):
    # only needed for [2016, 2018] combined analyses, which is stupid
    if "2016" in args.cardFileSM16:
        args.cardFileSM18 = args.cardFileSM16.replace("2016","scaled_2018")
    else:
        args.cardFileSM18 = args.cardFileSM16.replace(".","_scaled_2018.")
    if not args.mergeAndCheck:
        scaleCardFile( args.cardFileSM16, args.cardFileSM18, scale=lumi18/lumi16, copyUncertainties=True, keepObservation=False )
        logger.info("Scaled 2016 cardfile to 2018")
        regionCuts[2018] = getRegionCuts( args.cardFileSM16 )

logger.info("Preparing SM cardfiles done")
cardsSM    = { 2016:args.cardFileSM16 if args.cardFileSM16 else "none", 2017:args.cardFileSM17 if args.cardFileSM17 else "none", 2018:args.cardFileSM18 if args.cardFileSM18 else "none" }
regionCuts.update( { year:getRegionCuts( cardsSM[year] ) for year in args.years if year not in regionCuts.keys() } )

dbFile = "_".join( ["NLLcache"] + map( str, args.years ) + args.variables )
dbFile += ".sql"

tableName     = "nllcache"
cache_dir_nll = os.path.join(cache_directory, "nll")
if not os.path.isdir(cache_dir_nll): os.makedirs(cache_dir_nll)
dbPath        = os.path.join( cache_dir_nll, dbFile )
if args.nJobs>1:
    dbFileWrite   = dbFile.replace(".sql", "_%s"%str(uuid.uuid4()) + ".sql")
    dbPathWrite   = os.path.join(cache_dir_nll, "tmp" )
    if not os.path.isdir(dbPathWrite): os.makedirs( dbPathWrite )
    dbPathWrite   = os.path.join(dbPathWrite, dbFileWrite )
else:
    dbFileWrite   = dbFile
    dbPathWrite   = dbPath

allWC = list(set(["cpt", "cpQM", "ctp", "ctpI", "cpQM", "cpQ3", "cpt", "cpb", "cptb", "cptbI", "ctW", "ctZ", "ctWI", "ctZI", "cbW", "cbWI", "ctG", "ctGI", "cQlM1", "cQlM2", "cQlM3", "cQl31", "cQl32", "cQl33", "cQe1", "cQe2", "cQe3", "ctl1", "ctl2", "ctl3", "cte1", "cte2", "cte3", "ctlS1", "ctlSI1", "ctlS2", "ctlSI2", "ctlS3", "ctlSI3", "ctlT1", "ctlTI1", "ctlT2", "ctlTI2", "ctlT3", "ctlTI3", "cblS1", "cblSI1", "cblS2", "cblSI2", "cblS3", "cblSI3", "cQq83", "cQq81", "cQu8", "cQd8", "ctq8", "ctu8", "ctd8", "cQq13", "cQq11", "cQu1", "cQd1", "ctq1", "ctu1", "ctd1", "cQQ1", "cQQ8", "cQt1", "cQb1", "ctt1", "ctb1", "cQt8", "cQb8", "ctb8"]))
allWCDict = { var:0 for var in allWC }
tablesEntries  = [ "process", "years" ] + allWC
nllCache       = Cache( dbPath,      tableName, tablesEntries )
nllCacheWrite  = Cache( dbPathWrite, tableName, tablesEntries )
if nllCache      is None: raise
if nllCacheWrite is None: raise

logger.info("Loading gen-sample")
genSignalSample = ttZ_ll_LO_order2_15weights_ref

cardname  = [ genSignalSample.name ]
cardname += map( str, args.years )
cardname += [ args.genSelection ]
for i, var in enumerate(args.variables):
    cardname += [ var, "%sVal"%var ]
cardname += [ 'small' if args.small else 'full' ]
cardname  = '_'.join( cardname )

logger.info( "General card name: %s" %cardname )

def recoverDB( path, file ):
    command = os.path.expandvars("$CMSSW_BASE/src/TTZRun2EFT/Tools/scripts/recoverDB.sh")
    os.system("%s %s %s"%(command, file, path))

def isInDatabase( pointDict, database ):
    varDict = allWCDict
    for key, val in pointDict.iteritems():
        varDict[key] = val
    res = { "process":"ttZ", "years":"_".join( map( str, args.years ) ) }
    res.update( varDict )
    nCacheFiles = database.contains( res )
    return bool(nCacheFiles)

def chunks( l, n ):
    # split list in chuncs
    for i in range( 0, len(l), n ):
        yield l[i:i + n]

def splitList( l, n ):
    k, m = divmod(len(l), n)
    return list(l[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in xrange(n))

binDict = dict(zip( args.variables, [dict(zip( ["nBins", "min", "max"], var )) for var in chunks(args.binning, 3) ] ))

for var in binDict.values():
    var["nBins"] = int(var["nBins"])
    if var["nBins"] > 1:
        xRange       = np.linspace( var["min"], var["max"], var["nBins"], endpoint=False)
        halfstepsize = 0.5 * ( xRange[1] - xRange[0] )
        var["range"] = [0] + [ round(el + halfstepsize, 3) for el in xRange ]
    else:
        var["range"] = [0] + [ 0.5 * ( var["min"] + var["max"] ) ]

# get all possible combinations of points for i variables
points2D = list( itertools.product( *(binDict[var]["range"] for var in args.variables) ) ) #2D plots
points2D.sort( key = lambda point: tuple([abs(point[i]) for i in range(len(args.variables))]) )
points2D = [ dict(zip( args.variables, point ) ) for point in points2D ]

# check if already in database before calculations
if not args.overwrite:
    logger.info( "Checking %i points" %(len(points2D)) )
    points2D = [ point for point in points2D if not isInDatabase( point, nllCache ) ]
    logger.info( "Limit not calculated for %i points" %(len(points2D)) )

    if not points2D or args.mergeAndCheck:
        sys.exit(0)

if args.nJobs > len(points2D):
    logger.info("Batch mode job splitting larger than number of points. Setting nJobs to %i"%len(points2D))
    args.nJobs = len(points2D)

if args.job >= args.nJobs:
    logger.info("Job number > total number of jobs. Exiting")
    sys.exit(0)

if args.nJobs != 1:
    logger.info("Running batch mode job %i of total %i jobs"%(args.job,args.nJobs))
    total = len(points2D)
    points2D = splitList( points2D, args.nJobs)[args.job]
    logger.info("Calculating %i of total %i points"%(len(points2D), total))

logger.info( "%i points to calculate" %(len(points2D)) )

w = WeightInfo( genSignalSample.reweight_pkl )
w.set_order( args.order )

coeffList = {}
signal_genRateSM = {}
for year in args.years:
    coeffList[year] = {}
    signal_genRateSM[year] = {}

def calculation( pointDict ):
    
    uniqueDir = "/tmp/" + str(uuid.uuid4())
    if not os.path.isdir( uniqueDir ):
        os.makedirs( uniqueDir )

    c = CardFileWriter()
    c.reset()
    c.releaseLocation = combineReleaseLocation

    kwargs = pointDict

    cacheCard = cardname
    for i, var in enumerate(args.variables):
        cacheCard = cacheCard.replace("%sVal"%var, str(pointDict[var]))

    card     = {}
    cardpath = {}
    #caching results
    for year in args.years:
        logger.info( "Creating cardfile for year %i, %s"%( year, ", ".join( [ str(key) + " = " + str(val) for key, val in kwargs.iteritems() ] ) ) )

        card[year] = cacheCard.replace( "_".join(map(str,args.years)), str(year) )
        cardpath[year] = os.path.join( uniqueDir, card[year] + '.txt' )

        scale = {}
        for region in regionCuts[year].keys():

            logger.info( "At EFT region %s", region )

            signal_genRateEFT = w.get_weight_yield( coeffList[year][region], **kwargs )
            scale[region]     = signal_genRateEFT / signal_genRateSM[year][region]

        scaleCardFile( cardsSM[year], cardpath[year], scale=scale, scaledProcesses=["signal"], copyUncertainties=True, keepObservation=True )

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

    if not args.keepCards and not all( [ p==0 for p in pointDict.values() ] ):
        os.remove( combinedCard )
        if combinedAnalysis:
            for year in args.years:
                os.remove( cardpath[year] )
    else:
        source = combinedCard
        target = os.path.join( cardfileLocation, combinedCard.split("/")[-1] )
        copyfile( source, target )

    rmtree(uniqueDir)

    varDict = allWCDict
    for key, val in pointDict.iteritems():
        varDict[key] = val
    res = { "process":"ttZ", "years":"_".join( map( str, args.years ) ) }
    res.update( varDict )

    logger.info( "NLL limit for %s: nll_prefit = %f, nll_postfit = %f"%( ", ".join( ["%s = %s"%(key, val) for key, val in pointDict.iteritems() ] ), nll_prefit, nll_postfit) )
    nllCacheWrite.add( res, nll_prefit, overwrite=True )
    del c

def replaceAliases( cutString ):
    cut = cutString
    for key, val in aliases.iteritems():
        cut = cut.replace( key, val )
    logger.info( "Replacing variable names: old cut: %s, new cut: %s" %(cutString, cut) )
    return cut

def setup():
    # preparing gen-sample reweighting
    logger.info( "Preparing reweighting setup" )

    sel = {}
    for i, year in enumerate(args.years):

        sel[year] = {}
        logger.info( "At year %i", year )

        for region, cut in regionCuts[year].iteritems():

            logger.info( "At region %s", region )

            regionCut = replaceAliases( simpleStringToCutString( cut ) ) 
            sel[year][region] = "&&".join( [ cutInterpreter.cutString( args.genSelection ), regionCut ] )
            # safe some time
            if i > 0 and sel[year][region] == sel[args.years[0]][region]:
                coeffList[year][region]        = coeffList[args.years[0]][region]
                signal_genRateSM[year][region] = signal_genRateSM[args.years[0]][region]
            else:
                coeffList[year][region]        = w.getCoeffListFromDraw( genSignalSample, selectionString=sel[year][region] )
                signal_genRateSM[year][region] = float( w.get_weight_yield( coeffList[year][region] ) )
            logger.info( "Calculated SM gen-sample signal rate for region %s and year %i: %f" %(region, year, signal_genRateSM[year][region] ) )

# run initial setup
setup()

for i in range(10):

    if not points2D: break

    # run calculations
    pool = Pool( processes=args.cores )
    pool.map( calculation, points2D )
    pool.close()

    # recalculate missing points
    points2D = [ point for point in points2D if not isInDatabase( point, nllCacheWrite ) ]
