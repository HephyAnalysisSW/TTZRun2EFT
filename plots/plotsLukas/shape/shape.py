#!/usr/bin/env python
''' Analysis script for standard plots
'''

# Standard imports
import ROOT, os, imp, sys, copy
ROOT.gROOT.SetBatch(True)
import itertools
import pickle
from math                                import isnan, ceil, pi

# RootTools
from RootTools.core.standard             import *

# Internal Imports
from TTZRun2EFT.Tools.user               import plot_directory
from TTZRun2EFT.Tools.genCutInterpreter  import cutInterpreter

from Analysis.Tools.WeightInfo           import WeightInfo
from Analysis.Tools.helpers           import getCollection, deltaR

# Default Parameter
loggerChoices = ['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET']

# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',           action='store',      default='INFO', nargs='?', choices=loggerChoices,                                help="Log level for logging")
argParser.add_argument('--plotFile',           action='store',      default='all')
argParser.add_argument('--selection',          action='store',      default='all')
argParser.add_argument('--small',              action='store_true',                                                                                  help='Run only on a small subset of the data?', )
argParser.add_argument('--normalize',          action='store_true', default=False,                                                                   help="Normalize yields" )
args = argParser.parse_args()

# Logger
import Analysis.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

# Samples
from TopEFT.samples.gen_fwlite_benchmarks import *
#params = []
#params.append( {
#            'legendText': text,
#            'WC' : kwargs,
#            'name' : "_".join(kwargs.keys()),
#            'color' : ROOT.kRed,
#            } )

from Analysis.Tools.WeightInfo                      import WeightInfo
from TTZRun2EFT.Samples.genTuples_TTZ_postProcessed import *

color = [ ROOT.kBlack, ROOT.kRed ]

genSignalSample = ttZ_ll_LO_order2_15weights_ref
genSignalSample.reduceFiles( factor=9 )
w = WeightInfo( genSignalSample.reweight_pkl )
w.set_order( 2 )
genSignalSample.name = "ttZ"
genSignalSample.read_variables = [
                                  "ref_lumiweight1fb/F",
                                  VectorTreeVariable.fromString('p[C/F]', nMax=2000),
                                  "genLepZ_mass/F",
                                  "genLepZ_cosThetaStar/F",
                                  "genLepZ_pt/F",
                                  "genZ_mass/F",
                                  "genZ_cosThetaStar/F",
                                  "genZ_pt/F",
                                  "genLep[pt/F,eta/F]",
                                 ]

genSignalSample.setSelectionString( "abs(genZ_mass-91.2)<10&&Sum$(genLep_pt>15&&abs(genLep_eta)<2.4)>=2" )
#genSelection = "onZll-nJet3p"
#regionCut    = "(1)" #replaceAliases( simpleStringToCutString( cut ) ) 

#sel   = "&&".join( [ cutInterpreter.cutString( genSelection ), regionCut ] )
#coeffList = w.getCoeffListFromDraw( genSignalSample, selectionString=sel )
#signal_genRateSM  = float( w.get_weight_yield( coeffList ) )



def drawObjects( lumi_scale ):
    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextSize(0.04)
    tex.SetTextAlign(11) # align right
    lines = [
      (0.15, 0.95, 'CMS #bf{#it{Simulation Preliminary}}'), 
      (0.65, 0.95, '%3.1f fb{}^{-1} (13 TeV)' % lumi_scale)
    ]
    return [tex.DrawLatex(*l) for l in lines] 

if args.normalize:
    scaling = { 1:0 }

# Plotting
def drawPlots( plots, text ):
    for log in [False, True]:
        plot_directory_ = os.path.join( plot_directory, 'comparisonPlots', text, args.selection, "log" if log else "lin" )

        for plot in plots:
            if not max(l[0].GetMaximum() for l in plot.histos): 
                continue # Empty plot
            postFix = " (legacy)"
            extensions_ = ["pdf", "png", "root"]

            plotting.draw( plot,
	                       plot_directory = plot_directory_,
                           extensions = extensions_,
#                           ratio = {'yRange': (0.7, 1.3), 'histos':[(1,0)], 'texY':'Ratio'},
#	                       ratio = None,
	                       logX = False, logY = log, sorting = True,
	                       yRange = (0.03, "auto") if log else (0.001, "auto"),
	                       scaling = scaling if args.normalize else {},
	                       legend = [ (0.18,0.85-0.03*sum(map(len, plot.histos)),0.9,0.88), 2],
	                       drawObjects = drawObjects( lumi_scale ) if not args.normalize else drawObjects( lumi_scale ),
                           copyIndexPHP = True,
                         )

# Read variables and sequences
read_variables = []
def changeVarName( event, sample ):

    if sample.name != "ttZ":

        event.genLepZ_mass = event.Z_mass
        event.genLepZ_pt = event.Z_pt
        event.genLepZ_cosThetaStar = event.Z_cosThetaStar

        event.genZ_mass = event.Z_mass
        event.genZ_pt = event.Z_pt
        event.genZ_cosThetaStar = event.Z_cosThetaStar

# Sequence
sequence = [changeVarName]

lumi_scale = 136.6

weight_ = lambda event, sample: 1

for i, TopEFTSample in enumerate( [dim6top_central] + dim6top_currents ):

    try:
        cpQM = float(TopEFTSample.name.split("cpQM")[1].split("_")[1].replace("p",".").replace("m","-")) if i != 0 else 0
        cpt = float(TopEFTSample.name.split("cpt")[1].split("_")[1].replace("p",".").replace("m","-")) if  i != 0 else 0
    except:
        continue
    kwargs = {"cpQM":cpQM, "cpt":cpt}

    text = "cpQM_%.2f_cpt_%.2f"%(round(cpQM,2), round(cpt,2))
    if args.small:     text += "_small"
    if args.normalize: text += "_normalize"


    TopEFTSample.name = "TOP-18-009"
    TopEFTSample.read_variables = [
                                   "Z_mass/F",
                                   "Z_cosThetaStar/F",
                                   "Z_pt/F",
                                   "Z_daughterPdg/F",
                                   "GenLep[pt/F,eta/F]",
                                  ]

    TopEFTSample.setSelectionString("abs(Z_mass-91.2)<10&&(abs(Z_daughterPdg)==11||abs(Z_daughterPdg)==13)&&Sum$(GenLep_pt>15&&abs(GenLep_eta)<2.4)>=2")
    genSignalSample.weight = lambda event, sample:  w.get_weight_func( **kwargs )(event, sample)*event.ref_lumiweight1fb

    comparisonSamples = [ [genSignalSample], [TopEFTSample] ]

    stack      = Stack( *comparisonSamples )

    # Use some defaults (set defaults before you create/import list of Plots!!)
    Plot.setDefaults( stack=stack, weight=staticmethod( weight_ ), selectionString=cutInterpreter.cutString( args.selection ) )#, addOverFlowBin='upper' )

    # plotList
    addPlots = []

    addPlots.append( Plot(
        name      = 'Z_pt',
        texX      = 'p_{T}(ll)',
        texY      = 'a.u.',
        attribute = lambda event, sample: event.genLepZ_pt,
        binning   = [ 40, 0, 400 ],
    ))

    addPlots.append( Plot(
        name      = 'Z_pt_coarse',
        texX      = 'p_{T}(ll)',
        texY      = 'a.u.',
        attribute = lambda event, sample: event.genLepZ_pt,
        binning   = [ 6, 0, 600 ],
    ))

    addPlots.append( Plot(
        name      = 'Z_mass',
        texX      = 'm(ll)',
        texY      = 'a.u.',
        attribute = lambda event, sample: event.genLepZ_mass,
        binning   = [ 40, 60, 120 ],
    ))

    addPlots.append( Plot(
        name      = 'Z_cosThetaStar',
        texX      = 'cos(#theta*)(ll)',
        texY      = 'a.u.',
        attribute = lambda event, sample: event.genLepZ_cosThetaStar,
        binning   = [ 20, -1, 1 ],
    ))

    addPlots.append( Plot(
        name      = 'Z_cosThetaStar_coarse',
        texX      = 'cos(#theta*)(ll)',
        texY      = 'a.u.',
        attribute = lambda event, sample: event.genLepZ_cosThetaStar,
        binning   = [ 5, -1, 1 ],
    ))

    addPlots.append( Plot(
        name      = 'genZ_pt',
        texX      = 'p_{T}(Z)',
        texY      = 'a.u.',
        attribute = lambda event, sample: event.genZ_pt,
        binning   = [ 40, 0, 400 ],
    ))

    addPlots.append( Plot(
        name      = 'genZ_pt_coarse',
        texX      = 'p_{T}(Z)',
        texY      = 'a.u.',
        attribute = lambda event, sample: event.genZ_pt,
        binning   = [ 6, 0, 600 ],
    ))

    addPlots.append( Plot(
        name      = 'genZ_mass',
        texX      = 'm(Z)',
        texY      = 'a.u.',
        attribute = lambda event, sample: event.genZ_mass,
        binning   = [ 40, 60, 120 ],
    ))

    addPlots.append( Plot(
        name      = 'genZ_cosThetaStar',
        texX      = 'cos(#theta*)(Z)',
        texY      = 'a.u.',
        attribute = lambda event, sample: event.genZ_cosThetaStar,
        binning   = [ 20, -1, 1 ],
    ))

    addPlots.append( Plot(
        name      = 'genZ_cosThetaStar_coarse',
        texX      = 'cos(#theta*)(Z)',
        texY      = 'a.u.',
        attribute = lambda event, sample: event.genZ_cosThetaStar,
        binning   = [ 5, -1, 1 ],
    ))

    for j, sample in enumerate(stack.samples):
        sample.style = styles.lineStyle( color[j], width=2  )
        sample.scale = lumi_scale

    if args.small:
        for sample in stack.samples:
            sample.normalization=1.
            sample.reduceFiles( factor=10 )
            sample.scale /= sample.normalization

    # always initialize with [], elso you get in trouble with pythons references!
    plots  = []
    plots += addPlots

    plotting.fill( plots, read_variables=read_variables, sequence=sequence )

    drawPlots( plots, text )

