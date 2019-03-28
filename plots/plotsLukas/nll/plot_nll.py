''' Plot script WC parameter LogLikelihood
'''

# Standard imports 
import sys, os, pickle, copy, ROOT
import ctypes
import numpy as np

# Multiprocessing
from multiprocessing import Pool

# RootTools
from RootTools.core.standard   import *

# turn off graphics
ROOT.gROOT.SetBatch( True )

# TTZRun2EFT
from TTZRun2EFT.Tools.user              import plot_directory, cache_directory
from TTZRun2EFT.Tools.cutInterpreter    import cutInterpreter
from TTZRun2EFT.Tools.genCutInterpreter import cutInterpreter as genCutInterpreter

# get the reweighting function
from Analysis.Tools.WeightInfo          import WeightInfo
from Analysis.Tools.metFilters          import getFilterCut

from TTZRun2EFT.Tools.Cache             import Cache

# Default Parameter
loggerChoices = ['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET']

# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',           action='store',      default='INFO', nargs='?', choices=loggerChoices,                                help="Log level for logging")
argParser.add_argument('--genSelection',       action='store',      default='onZll')
argParser.add_argument('--small',              action='store_true',                                                                                  help='Run only on a small subset of the data?', )
argParser.add_argument('--profiled',           action='store_true',                                                                                  help='Profile datapoints?', )
argParser.add_argument('--years',              action='store',      default=[ 2016, 2017 ], type=int, choices=[2016, 2017, 2018], nargs="*",         help="Which years to combine?")
argParser.add_argument('--cardFileSM16',       action='store',      default=None, type=str,                                                          help="SM cardfile for 2016")
argParser.add_argument('--cardFileSM17',       action='store',      default=None, type=str,                                                          help="SM cardfile for 2017")
argParser.add_argument('--cardFileSM18',       action='store',      default=None, type=str,                                                          help="SM cardfile for 2018")
argParser.add_argument('--variables' ,         action='store',      default = ['ctZ', 'ctZI'],       type=str,   nargs=*,                            help="argument variables")
argParser.add_argument('--plotVariables' ,     action='store',      default = ['ctZ', 'ctZI'],       type=str,   nargs=2,                            help="argument plotting variables")
argParser.add_argument('--binning',            action='store',      default=[30, -2, 2, 30, -2, 2 ], type=float, nargs=6,                            help="argument parameters")
argParser.add_argument('--checkOnly',          action='store_true',                                                                                  help='Just check if yield is already calculated', )
argParser.add_argument('--contours',           action='store_true',                                                                                  help='draw 1sigma and 2sigma contour line?')
argParser.add_argument('--smooth',             action='store_true',                                                                                  help='smooth histogram?')
argParser.add_argument('--zRange',             action='store',      default=[None, None],      type=float, nargs=2,                                  help="argument parameters")
argParser.add_argument('--xyRange',            action='store',      default=[None, None, None, None],  type=float, nargs=4,                          help="argument parameters")
argParser.add_argument('--binMultiplier',      action='store',      default=3,                 type=int,                                             help='bin multiplication factor')
argParser.add_argument('--skipMissingPoints',  action='store_true',                                                                                  help='Set missing NLL points to 999?')
argParser.add_argument('--tag',                action='store',      default="combined",        type=str,                                             help='bin multiplication factor')
args = argParser.parse_args()

args = argParser.parse_args()

# Logger
import Analysis.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(    args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger( args.logLevel, logFile = None)

if args.cardFileSM16 and (not args.cardFileSM17 and 2017 in args.years):
    if "2016" in args.cardFileSM16:
        args.cardFileSM17 = args.cardFileSM16.replace("2016","scaled_2017")
    else:
        args.cardFileSM17 = args.cardFileSM16.replace(".","_scaled_2017.")

if args.cardFileSM17 and (not args.cardFileSM18 and 2018 in args.years):
    # prefer a 2017->2018 scaling over 2016->2018 scaling if possible
    if "2017" in args.cardFileSM17:
        args.cardFileSM18 = args.cardFileSM17.replace("2017","scaled_2018")
    else:
        args.cardFileSM18 = args.cardFileSM17.replace(".","_scaled_2018.")

if args.cardFileSM16 and (not args.cardFileSM18 and 2018 in args.years):
    # only needed for [2016, 2018] combined analyses, which is stupid
    if "2016" in args.cardFileSM16:
        args.cardFileSM18 = args.cardFileSM16.replace("2016","scaled_2018")
    else:
        args.cardFileSM18 = args.cardFileSM16.replace(".","_scaled_2018.")

tableName = "nllcache"
cache_dir_nll    = os.path.join(cache_directory, "nll")
if not os.path.isdir(cache_dir_nll):
    os.makedirs(cache_dir_nll)
dbFile = "_".join( ["NLLcache"] + map( str, args.years ) ) + ".sql"
dbPath = os.path.join(cache_dir_nll, dbFile)

nllCache  = Cache( dbPath, tableName, ["cardname16", "cardname17", "cardname18", "cardname", "WC", "WC_val", "nll_prefit", "nll_postfit" ] )
if nllCache is None: raise

if set( args.plotVariables ) - set( args.variables ):
    raise Exception( "PlotVariables not in cached data!" )

# sort in same order as args.variables
args.plotVariables = [ var for var in args.variables if var in args.plotVariables ]
if len( args.plotVariables ) != 2:
    raise Exception( "Need two plotVariables!" )

# Gen Samples
from TTZRun2EFT.Samples.genTuples_TTZ_postProcessed import *
genSignalSample = ttZ_ll_LO_order2_15weights_ref

cardname  = [ genSignalSample.name ]
cardname += map( str, args.years )
cardname += [ args.genSelection ]
for i, var in enumerate(args.variables):
    cardname += [ var, "var%i"%i ]
cardname += [ 'small' if args.small else 'full' ]
cardname  = '_'.join( cardname )

logger.info( "General card name: %s" %cardname )

lumi = {}
lumi[2016] = 35.92
lumi[2017] = 41.86
lumi[2018] = 58.83
lumi_scale = sum( [ lumi[year] for year in args.years ])

def getNllData( point ):
    card = cardname
    for i, var in enumerate(point):
        card = card.replace("var%i"%i, str(var))

    res = {"cardname16":cardsSM[2016], "cardname17":cardsSM[2017], "cardname18":cardsSM[2018], "cardname":card, "WC":"_".join(args.variables), "WC_val":"_".join(map(str,point))}
    nCacheFiles = nllCache.contains( res )

    if nCacheFiles:
        cachedDict = nllCache.getDicts( res )[0]
        nll = cachedDict["nll_prefit"]
    else: 
        logger.info("Data for %s=%s and %s=%s not in cache!"%( args.variables[0], str(var1), args.variables[1], str(var2) ) )
        if args.skipMissingPoints: nll = 999
        else: sys.exit(1)
    return float(nll)


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

# get all possible combinations of points for i variables, set non-plot variables tto 0 except for profiled plot
points2D = list( itertools.product( *(binDict[var]["range"] if var in args.plotVariables or args.profiled else [0]*len(binDict[var]["range"]) for var in args.variables) ) ) #2D plots

logger.info("Loading cache data" )

sm_nll  = getNllData( tuple( [0]*len(args.variables) ) )
nllData = [ tuple( list( point ) + [ 2*(getNllData( point )-sm_nll) ] ) for point in points2D ]

# now reduce dimensions for plot, either profiled or just remove other dimensions
if args.profiled:
    profNll = []
    varIndexX = args.variables.index(args.plotVariables[0])
    varIndexY = args.variables.index(args.plotVariables[1])
    for x in binDict[args.plotVariables[0]]["range"]:
        for y in binDict[args.plotVariables[1]]["range"]:
            nll = filter( lambda point: point[varIndexX]==x and point[varIndexY]==y, nllData )
            if not nll: continue
            nll.sort( key = lambda point: point[-1] )
            profNll.append( (x, y, nll[0]) )
    nllData = profNll
else:
    nllData = [ tuple( [ p for i, p in enumerate(point[:-1]) if args.variables[i] in args.plotVariables] + [ point[-1] ] ) for point in nllData ]


def toGraph2D( name, title, data ):
    result = ROOT.TGraph2D( len(data) )
    debug = ROOT.TGraph()
    result.SetName( name )
    result.SetTitle( title )
    for i, datapoint in enumerate(data):
        x, y, val = datapoint
        result.SetPoint(i, x, y, val)
        debug.SetPoint(i, x, y)
    c = ROOT.TCanvas()
    result.Draw()
    debug.Draw()
    del c
    return result, debug

#get TGraph2D from results list
title = "TTZ_%s_%s"%(args.plotVariables[0], args.plotVariables[1])
a, debug = toGraph2D( title, title, nllData )
nxbins   = max( 1, min( 500, xbins*args.binMultiplier ) )
nybins   = max( 1, min( 500, ybins*args.binMultiplier ) )

#re-bin
hist = a.GetHistogram().Clone()
a.SetNpx( nxbins )
a.SetNpy( nybins )
hist = a.GetHistogram().Clone()

#smoothing
if args.smooth: hist.Smooth()

cans = ROOT.TCanvas("can_%s"%title,"",500,500)

contours = [2.28, 5.99]# (68%, 95%) for 2D
if args.contours:
    histsForCont = hist.Clone()
    c_contlist = ((ctypes.c_double)*(len(contours)))(*contours)
    histsForCont.SetContour(len(c_contlist),c_contlist)
    histsForCont.Draw("contzlist")
    cans.Update()
    conts = ROOT.gROOT.GetListOfSpecials().FindObject("contours")
    cont_p1 = conts.At(0).Clone()
    cont_p2 = conts.At(1).Clone()

pads = ROOT.TPad("pad_%s"%title,"",0.,0.,1.,1.)
pads.SetRightMargin(0.20)
pads.SetLeftMargin(0.14)
pads.SetTopMargin(0.11)
pads.Draw()
pads.cd()

hist.Draw("colz")

#draw contour lines
if args.contours:
    for conts in [cont_p2]:
        for cont in conts:
            cont.SetLineColor(ROOT.kOrange+7)
            cont.SetLineWidth(3)
#            cont.SetLineStyle(7)
            cont.Draw("same")
    for conts in [cont_p1]:
        for cont in conts:
            cont.SetLineColor(ROOT.kSpring-1)
            cont.SetLineWidth(3)
#            cont.SetLineStyle(7)
            cont.Draw("same")


hist.GetZaxis().SetTitle("-2 #Delta ln L")

if not None in args.zRange:
    hist.GetZaxis().SetRangeUser( args.zRange[0], args.zRange[1] )
if not None in args.xyRange[:2]:
    hist.GetXaxis().SetRangeUser( args.xyRange[0], args.xyRange[1] )
if not None in args.zRange[2:]:
    hist.GetYaxis().SetRangeUser( args.xyRange[2], args.xyRange[3] )

xTitle = args.plotVariables[0].replace("c", "C_{").replace("I", "}^{[Im]").replace('p','#phi') + '}'
yTitle = args.plotVariables[1].replace("c", "C_{").replace("I", "}^{[Im]").replace('p','#phi') + '}'
hist.GetXaxis().SetTitle( xTitle + ' (#Lambda/TeV)^{2}' )
hist.GetYaxis().SetTitle( yTitle + ' (#Lambda/TeV)^{2}' )

hist.GetXaxis().SetTitleFont(42)
hist.GetYaxis().SetTitleFont(42)
hist.GetZaxis().SetTitleFont(42)
hist.GetXaxis().SetLabelFont(42)
hist.GetYaxis().SetLabelFont(42)
hist.GetZaxis().SetLabelFont(42)

hist.GetXaxis().SetTitleOffset(1.15)
hist.GetYaxis().SetTitleOffset(1.25)

hist.GetXaxis().SetTitleSize(0.045)
hist.GetYaxis().SetTitleSize(0.045)
hist.GetZaxis().SetTitleSize(0.042)
hist.GetXaxis().SetLabelSize(0.04)
hist.GetYaxis().SetLabelSize(0.04)
hist.GetZaxis().SetLabelSize(0.04)

latex1 = ROOT.TLatex()
latex1.SetNDC()
latex1.SetTextSize(0.035)
latex1.SetTextFont(42)
latex1.SetTextAlign(11)

addon = " (profiled)" if args.profiled else ""
latex1.DrawLatex(0.12, 0.92, '#bf{CMS} #it{Simulation Preliminary}' + addon),
latex1.DrawLatex(0.63, 0.92, '#bf{%3.1f fb{}^{-1} (13 TeV)}' % lumi_scale)

latex2 = ROOT.TLatex()
latex2.SetNDC()
latex2.SetTextSize(0.04)
latex2.SetTextFont(42)
latex2.SetTextAlign(11)

y   = str(args.years[0]) if len(args.years)==1 else "combined"
sel = selection[args.selections[0]] if len(args.selections)==1 else "combined"
plot_directory_ = os.path.join( plot_directory, "NLLPlots", y, "_".join(args.variables) )

if not os.path.isdir( plot_directory_ ):
    try: os.makedirs( plot_directory_ )
    except: pass

plotname = "%s%s"%("_".join(args.plotVariables), "_profiled" if args.profiled else "")
for e in [".png",".pdf",".root"]:
    cans.Print( plot_directory_ + "/%s%s"%(plotname, e) )

