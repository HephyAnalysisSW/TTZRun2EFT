# Standard
import os

# RootTools
from RootTools.core.standard import *

# TTZRun2EFT
from TTZRun2EFT.Tools.user   import gridpack_directory, cache_directory

def get_parser():
    import argparse
    argParser = argparse.ArgumentParser(description = "Argument parser for samples file")
    argParser.add_argument( '--overwrite',   action='store_true',    help="Overwrite current entry in db?" )
    return argParser

# Logging
if __name__=="__main__":
    import Analysis.Tools.logger as logger
    logger = logger.get_logger("DEBUG", logFile = None )
    import RootTools.core.logger as logger_rt
    logger_rt = logger_rt.get_logger("DEBUG", logFile = None )
    options = get_parser().parse_args()
    ov = options.overwrite
else:
    import logging
    logger = logging.getLogger(__name__)
    ov = False

dbFile = cache_directory + "/samples/DB_TTZ_GEN.sql"

logger.info( "Using db file: %s", dbFile )

ttZ_ll_LO_order2_15weights_ref               = FWLiteSample.fromDAS("ttZ_ll_LO_order2_15weights_ref", "/ttZ0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball/llechner-ttZ0j_order2_15weights_18052018_ref-7a5fde3f5bf89006ee3acec926ca87d8/USER", "phys03", dbFile = dbFile, overwrite=ov, prefix='root://hephyse.oeaw.ac.at/' )
ttZ_ll_LO_order2_15weights_ref.reweight_pkl  = os.path.join( gridpack_directory, "EFT/ttZ_ll_LO_order2_15weights_ref/", "ttZ0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl" )
ttZ_ll_LO_order2_15weights_ref.nEvents       = 1000000
ttZ_ll_LO_order2_15weights_ref.xsec          = 0.5205 * 0.0915 / 0.0565 #pb ttZ, Z->ll, ttZ gridpack * ttZ NLO Daniel / ttZ LO run.py UFO

SM = [
]

EFT = [
    ttZ_ll_LO_order2_15weights_ref,
]

allSamples = SM + EFT

for s in allSamples:
    s.isData = False
    print s.reweight_pkl
