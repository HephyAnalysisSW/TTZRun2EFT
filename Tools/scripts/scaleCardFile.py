#!/usr/bin/env python
import os, copy

from Analysis.Tools.cardFileHelpers import scaleCardFile
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
argParser.add_argument('--scaledProcesses',    action='store',      default=None,   nargs="*", type=str,              help='Which processes for luminosity scale?')
argParser.add_argument('--copyUncertainties',  action='store_true',                                                   help='Copy all uncertainties from inputfile')
args = argParser.parse_args()

# Logger
import Analysis.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(    args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger( args.logLevel, logFile = None)

scaleCardFile(args.inFile, args.outFile, scale=args.lumiScale, scaledProcesses=args.scaledProcesses, copyUncertainties=args.copyUncertainties)

