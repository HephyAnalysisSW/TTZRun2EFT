#!/usr/bin/env python
import os, copy

from TTZRun2EFT.Tools.infoFromCards import *
from Analysis.Tools.CardFileWriter  import CardFileWriter

# Default Parameter
loggerChoices = ['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET']

# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',           action='store',      default='INFO', nargs='?', choices=loggerChoices, help="Log level for logging")
argParser.add_argument('--inFile',             action='store',      default="cardFile.dat",                           help="input args.inFilepath")
argParser.add_argument('--outFile',            action='store',      default=None,                                     help="output args.inFilepath")
argParser.add_argument('--lumiScale',          action='store',      default=1.,     type=float,                       help='Factor for luminosity scale')
argParser.add_argument('--copyUncertainties',  action='store_true',                                                   help='Copy all uncertainties from inputfile')
args = argParser.parse_args()

# Logger
import Analysis.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(    args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger( args.logLevel, logFile = None)

if not args.outFile: args.outFile = args.inFile.split(".")[0] + "_out" + args.inFile.split(".")[1]
if not args.inFile.startswith("/"):    args.inFile    = os.path.join( os.getcwd(), args.inFile )
if not args.outFile.startswith("/"): args.outFile = os.path.join( os.getcwd(), args.outFile )

regions     = getAllBinNames(args.inFile)
processes   = getAllProcesses(args.inFile)
bgProcesses = copy.copy(processes)
bgProcesses.remove("signal")

c = CardFileWriter()

for i_region, region in enumerate(regions):

    totYield = 0
    for proc in processes:
        rate = round( getEstimateFromCard(args.inFile, proc, region).val * args.lumiScale, 5 )
        c.specifyExpectation( region, proc, rate )
        totYield += rate

    c.addBin( region, bgProcesses, region)
    c.specifyObservation( region, int(round(totYield)) )

c.writeToFile( args.outFile )

# read uncertainty lines and simply append it to the new cardfile
flag = False
if args.copyUncertainties:

    with open(args.inFile, "r") as f:
        lines = f.readlines()

    writeLines = []
    for line in lines:
        if flag: writeLines.append(line)
        elif line.startswith("rate"): flag = True

    with open(args.outFile, "a") as f:
        for line in writeLines[1:]:
            f.write(line)

for bin in regions:
    obsIn  = int(getObservationFromCard(args.inFile, bin).val * args.lumiScale)
    obsOut = getObservationFromCard(args.outFile, bin).val
    if abs(obsIn - obsOut) > 1:
        raise Exception("Error in observation in bin %s! In: %i - Out: %i"%(bin, obsIn, obsOut ))
    for proc in processes:
        rateIn  = getEstimateFromCard(args.inFile, proc, bin).val * args.lumiScale
        rateOut = getEstimateFromCard(args.outFile, proc, bin).val
        if abs(rateIn - rateOut) > 0.1:
            raise Exception("Error in rate in process %s in bin %s! In: %f - Out: %f"%(proc, bin, rateIn, rateOut))

logger.info("Output cardfile checks OK!")
