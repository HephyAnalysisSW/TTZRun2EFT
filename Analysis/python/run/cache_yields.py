''' Plot script pTG shape + WC
'''

# Standard imports 
import sys, os, pickle, copy, ROOT

# RootTools
from RootTools.core.standard   import *

# turn off graphics
ROOT.gROOT.SetBatch( True )

# TTGammaEFT
from TTGammaEFT.Tools.user              import cache_directory
from TTGammaEFT.Tools.TriggerSelector   import TriggerSelector
from TTGammaEFT.Tools.cutInterpreter    import cutInterpreter
from TTGammaEFT.Tools.Cache             import Cache
from TTGammaEFT.Analysis.regions        import recoTTGammaRegions as regions
from Analysis.Tools.metFilters          import getFilterCut

# Default Parameter
loggerChoices = ['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET']

# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',           action='store',      default='INFO', nargs='?', choices=loggerChoices,                                help="Log level for logging")
argParser.add_argument('--selection',          action='store',      default='dilepOS-nLepVeto2-pTG20-nPhoton1p-offZSFllg-offZSFll-mll40-nJet2p-nBTag1p')
argParser.add_argument('--small',              action='store_true',                                                                                  help='Run only on a small subset of the data?', )
argParser.add_argument('--year',               action='store',      default=None,   type=int,  choices=[2016,2017,2018],                             help="Which year to plot?")
argParser.add_argument('--noData',             action='store_true', default=False,                                                                   help='also plot data?')
argParser.add_argument('--overwrite',          action='store_true',                                                                                  help='Overwrite Database entries?', )
argParser.add_argument('--inclusive',          action='store_true',                                                                                  help='run inclusive regions', )
argParser.add_argument('--checkOnly',          action='store_true',                                                                                  help='Just check if yield is already calculated', )
argParser.add_argument('--cores',              action='store',      default=1,                       type=int,                                       help='number of cpu cores for multicore processing')
args = argParser.parse_args()


if args.inclusive:
    from TTGammaEFT.Analysis.regions        import recoTTGammaRegionsIncl as regions
else:
    from TTGammaEFT.Analysis.regions        import recoTTGammaRegions     as regions


# Samples
if args.year == 2016:
    from TTGammaEFT.Samples.nanoTuples_Summer16_private_postProcessed      import *
    if not args.noData:
        from TTGammaEFT.Samples.nanoTuples_Run2016_14Dec2018_postProcessed import *

elif args.year == 2017:
    from TTGammaEFT.Samples.nanoTuples_Fall17_private_postProcessed        import *
    if not args.noData:
        from TTGammaEFT.Samples.nanoTuples_Run2017_14Dec2018_postProcessed import *

elif args.year == 2018:
    from TTGammaEFT.Samples.nanoTuples_Autumn18_private_postProcessed      import *
    if not args.noData:
        from TTGammaEFT.Samples.nanoTuples_Run2018_14Dec2018_postProcessed import *

# Logger
import Analysis.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

# Sample definition
if args.year == 2016:
    mc = [ DY_LO_16, TT_pow_16, ZG_16, singleTop_16, other_16, TTG_16 ]

elif args.year == 2017:
    mc = [ DY_LO_17, TT_pow_17, singleTop_17, other_17, TTG_17 ]

elif args.year == 2018:
    mc = [ DY_LO_18, TT_pow_18, TTG_18, singleTop_18, other_18 ]

if args.noData:
    if args.year == 2016:   lumi_scale = 35.92
    elif args.year == 2017: lumi_scale = 41.86
    elif args.year == 2018: lumi_scale = 58.83
    allSamples = mc
else:
    if args.year == 2016:   data_sample = Run2016
    elif args.year == 2017: data_sample = Run2017
    elif args.year == 2018: data_sample = Run2018
    data_sample.texName        = "data (legacy)"
    data_sample.name           = "data"
    data_sample.style          = styles.errorStyle( ROOT.kBlack )
    lumi_scale                 = data_sample.lumi * 0.001
    data_sample.setWeightString("weight")
    data_sample.setSelectionString( [ getFilterCut( args.year, isData=True ), cutInterpreter.cutString( args.selection ) ] )
    allSamples = [data_sample] + mc

tr = TriggerSelector( args.year )

for sample in mc:
    sample.style          = styles.fillStyle( sample.color )
    sample.setWeightString("weight*%f*reweightDilepTriggerBackup*reweightPU*reweightLeptonSF*reweightLeptonTrackingSF*reweightPhotonSF*reweightPhotonElectronVetoSF*reweightBTag_SF"%lumi_scale)
    sample.setSelectionString( [ getFilterCut( args.year, isData=False ), tr.getSelection( "MC" ), cutInterpreter.cutString( args.selection ) ] )

if args.year == 2016:
    TTG_16.addSelectionString( "isTTGamma==1" )
    TT_pow_16.addSelectionString( "isTTGamma==0" )
    ZG_16.addSelectionString( "isZWGamma==1" )
    DY_LO_16.addSelectionString( "isZWGamma==0" )
if args.year == 2017:
    TTG_17.addSelectionString( "isTTGamma==1" )
    TT_pow_17.addSelectionString( "isTTGamma==0" )
if args.year == 2018:
    TTG_18.addSelectionString( "isTTGamma==1" )
    TT_pow_18.addSelectionString( "isTTGamma==0" )

if args.small:

    for sample in allSamples:
        sample.normalization=1.
        sample.reduceFiles( factor=15 )
        sample.addWeightString("%f"%(1./sample.normalization))

cache_dir = os.path.join(cache_directory, "yields")
if not os.path.isdir(cache_dir):
    os.mkdir(cache_dir)

def cacheYields( (i_region, i_sample) ):

    region = regions[i_region]
    sample = allSamples[i_sample]

    logger.info( "At region %s for sample %s"%(region, sample.name) )

    dbFilename = "yields_%s.sql"%sample.name
    yieldDB = Cache( os.path.join( cache_dir, dbFilename ), "yields", [ "selection", "year", "small", "region", sample.name] )
    if not yieldDB: raise

    res = {"selection":args.selection, "year":args.year, "small":args.small, "region":str(region), sample.name:sample.name }

    if yieldDB.contains( res ) and not args.overwrite:
        if not args.checkOnly:
            logger.info( "Yield for sample %s at region %s already in database: %s" %(sample.name, region, yieldDB.getDicts( res )[0]['value']) )
        return

    if args.checkOnly:
        logger.info( "Yield for sample %s at region %s not processed" %(sample.name, region) )
        return

    rate = sample.getYieldFromDraw( selectionString=region.cutString() )['val']
    logger.info("Adding yield for sample %s at region %s: %s"%(sample.name, str(region), str(rate)) )
    yieldDB.add( res, str(rate), overwrite=True )

# Multiprocessing
from multiprocessing import Pool

for region in regions:
    region.name = str(region) # needed for Pool

# Objects can't be pickled, so args are indices of lists
input = [ (i_region, i_sample) for i_region, _ in enumerate(regions) for i_sample, _ in enumerate(allSamples) ]

logger.info( "Calculating %i yields"%len(input) )

pool = Pool( processes=args.cores )
_    = pool.map( cacheYields, input )
pool.close()

