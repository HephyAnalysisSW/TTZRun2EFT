#!/usr/bin/env python

# Standard imports
import argparse
import os

# TTZRun2EFT imports
from TTZRun2EFT.Generation.Configuration import Configuration
from TTZRun2EFT.Generation.Process       import Process
from TTZRun2EFT.Tools.user               import plot_directory
from Analysis.Tools.u_float              import u_float

# Logging
if __name__=="__main__":
    import Analysis.Tools.logger as logger
    logger = logger.get_logger("INFO", logFile = None )
    import RootTools.core.logger as logger_rt
    logger_rt = logger_rt.get_logger("INFO", logFile = None )

else:
    import logging
    logger = logging.getLogger(__name__)

# Find all processes
process_path     = os.path.expandvars( "$CMSSW_BASE/src/TTZRun2EFT/Generation/data/processCards" )
defaultProcesses = [ os.path.splitext(f)[0] for f in os.listdir( process_path ) if os.path.isfile( os.path.join( process_path, f ) ) and f.endswith( '.dat' ) ]
defaultModels    = [ 'dim6top_LO' ]
defaultLogger    = [ 'CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET' ]

argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--process',     action='store',      default='TTGamma_SingleLeptFromT_3LineBuggy', choices=defaultProcesses, help="Which process?")
argParser.add_argument('--model',       action='store',      default='dim6top_LO',                         choices=defaultModels,    help="Which madgraph model?")
argParser.add_argument('--couplings',   action='store',      default=[],         nargs='*',  type = str,                             help="Give a list of the non-zero couplings with values, e.g. NAME1 VALUE1 NAME2 VALUE2")
argParser.add_argument('--overwrite',   action='store_true',                                                                         help="Overwrite exisiting x-sec calculation and gridpack")
argParser.add_argument('--logLevel',    action='store',      default='INFO',     nargs='?',                choices=defaultLogger,    help="Log level for logging" )

args = argParser.parse_args()

logger = logger.get_logger( args.logLevel, logFile = None )

# Check that we have an even number of arguments
if not len(args.couplings)%2==0:
    logger.error("Need an even number of coupling arguments of the format coupling1, value1, coupling2, value2, ... . Got %r", args.couplings )

# Interpret coupling argument list
coupling_names  = args.couplings[::2]
coupling_values = map( float, args.couplings[1::2] )

modified_couplings = { c:v for c,v in zip( coupling_names, coupling_values ) }

# Let's not leave the user in the dark
logger.info( "Model:        %s", args.model )
logger.info( "Process:      %s", args.process )
logger.info( "Couplings:    %s", ", ".join( [ "%s=%5.4f" % c for c in modified_couplings.items() ] ) )

# Create configuration class
config = Configuration( model_name = args.model )

p        = Process( process=args.process, nEvents=1, config=config )
xsec_val = p.diagrams( plot_dir=os.path.join( plot_directory, "diagrams" ), modified_couplings=modified_couplings )

config.cleanup()

logger.info( "Done! Calculated xsec: %s ", repr(xsec_val) )
