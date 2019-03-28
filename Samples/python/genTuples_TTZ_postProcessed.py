# Standard Imports
import os, sys
import ROOT

# RootTools Imports
from RootTools.core.Sample import Sample

# Logging
if __name__=="__main__":
    import Analysis.Tools.logger as logger
    logger = logger.get_logger("INFO", logFile = None )
    import RootTools.core.logger as logger_rt
    logger_rt = logger_rt.get_logger("INFO", logFile = None )
else:
    import logging
    logger = logging.getLogger(__name__)

# Colors
from TTZRun2EFT.Samples.color import color

# Data directory
from TTZRun2EFT.Tools.user import gridpack_directory
from TTZRun2EFT.Tools.user import data_directory1             as data_directory
from TTZRun2EFT.Tools.user import postprocessing_directoryGEN as postprocessing_directory

logger.info( "Loading MC samples from directory %s", os.path.join( data_directory, postprocessing_directory ) )

# Directories
dirs = {}

dirs['TTZ_EFT'] = [ "ttZ_ll_LO_order2_15weights_ref" ]

directories = { key : [ os.path.join( data_directory, postprocessing_directory, dir) for dir in dirs[key] ] for key in dirs.keys() }

# Samples
ttZ_ll_LO_order2_15weights_ref               = Sample.fromDirectory( name="ttZ_ll_LO_order2_15weights_ref", treeName="Events", isData=False, color=color.TTZ, texName="ttZ", directory=directories['TTZ_EFT'] )
ttZ_ll_LO_order2_15weights_ref.reweight_pkl  = os.path.join( gridpack_directory, "EFT/ttZ_ll_LO_order2_15weights_ref/", "ttZ0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl" )

