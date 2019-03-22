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

dbFile = cache_directory + "/samples/DB_TTGamma_GEN.sql"

logger.info( "Using db file: %s", dbFile )

# SM point, no weights
# test samples
TTGamma_SingleLeptFromT_SM_1Line_test                     = FWLiteSample.fromDAS("TTGamma_SingleLeptFromT_SM_1Line_test", "/TTGamma_SingleLeptFromT_1Line_test_GENSIM/llechner-TTGamma_SingleLeptFromT_1Line_test_GENSIM-b3122b77720ca80362007ff086008b03/USER", "phys03", dbFile=dbFile, overwrite=ov, prefix='root://hephyse.oeaw.ac.at/')
TTGamma_SingleLeptFromT_SM_1Line_test.reweight_pkl        = os.path.join( gridpack_directory, "SM/TTGamma_SingleLeptFromT_1Line_test/", "TTGamma_SingleLeptFromT_1Line_test_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl" )
TTGamma_SingleLeptFromT_SM_1Line_test.xSection            = 0.1267
TTGamma_SingleLeptFromT_SM_1Line_test.nEvents             = 100000

TTGamma_SingleLeptFromT_SM_3LineBuggy_test                = FWLiteSample.fromDAS("TTGamma_SingleLeptFromT_SM_3LineBuggy_test", "/TTGamma_SingleLeptFromT_3LineBuggy_test_GENSIM/llechner-TTGamma_SingleLeptFromT_3LineBuggy_test_GENSIM-f3ebe776d7309479b1058aff34e27155/USER", "phys03", dbFile=dbFile, overwrite=ov, prefix='root://hephyse.oeaw.ac.at/')
TTGamma_SingleLeptFromT_SM_3LineBuggy_test.reweight_pkl   = os.path.join( gridpack_directory, "SM/TTGamma_SingleLeptFromT_3LineBuggy_test/", "TTGamma_SingleLeptFromT_3LineBuggy_test_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl" )
TTGamma_SingleLeptFromT_SM_3LineBuggy_test.xSection       = 0.1266
TTGamma_SingleLeptFromT_SM_3LineBuggy_test.nEvents        = 100000

TTGamma_SingleLeptFromT_SM_3LinePatched_test              = FWLiteSample.fromDAS("TTGamma_SingleLeptFromT_SM_3LinePatched_test", "/TTGamma_SingleLeptFromT_3LinePatched_test_GENSIM/llechner-TTGamma_SingleLeptFromT_3LinePatched_test_GENSIM-e94234d067bc2f081eff982fc497c503/USER", "phys03", dbFile=dbFile, overwrite=ov, prefix='root://hephyse.oeaw.ac.at/')
TTGamma_SingleLeptFromT_SM_3LinePatched_test.reweight_pkl = os.path.join( gridpack_directory, "SM/TTGamma_SingleLeptFromT_3LinePatched_test/", "TTGamma_SingleLeptFromT_3LinePatched_test_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl" )
TTGamma_SingleLeptFromT_SM_3LinePatched_test.xSection     = 0.1293
TTGamma_SingleLeptFromT_SM_3LinePatched_test.nEvents      = 100000

# full proc card
TTGamma_SingleLeptFromT_SM_1Line                          = FWLiteSample.fromDAS("TTGamma_SingleLeptFromT_SM_1Line", "/TTGamma_SingleLeptFromT_1Line_SM_GENSIM_v2/llechner-TTGamma_SingleLeptFromT_1Line_SM_GENSIM_v2-75ed994ddf291590a34f99ee2b8b5398/USER", "phys03", dbFile=dbFile, overwrite=ov, prefix='root://hephyse.oeaw.ac.at/')
TTGamma_SingleLeptFromT_SM_1Line.reweight_pkl             = os.path.join( gridpack_directory, "SM/TTGamma_SingleLeptFromT_1Line/", "TTGamma_SingleLeptFromT_1Line_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl" )
TTGamma_SingleLeptFromT_SM_1Line.xSection                 = 0.7491
TTGamma_SingleLeptFromT_SM_1Line.nEvents                  = 100000

TTGamma_SingleLeptFromT_SM_3LinePatched                   = FWLiteSample.fromDAS("TTGamma_SingleLeptFromT_SM_3LinePatched", "/TTGamma_SingleLeptFromT_3LinePatched_SM_GENSIM/llechner-TTGamma_SingleLeptFromT_3LinePatched_SM_GENSIM-cf918753c14f0b139153df060b08cdc2/USER", "phys03", dbFile=dbFile, overwrite=ov, prefix='root://hephyse.oeaw.ac.at/')
TTGamma_SingleLeptFromT_SM_3LinePatched.reweight_pkl      = os.path.join( gridpack_directory, "SM/TTGamma_SingleLeptFromT_3LinePatched/", "TTGamma_SingleLeptFromT_3LinePatched_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl" )
TTGamma_SingleLeptFromT_SM_3LinePatched.xSection          = 0.7664
TTGamma_SingleLeptFromT_SM_3LinePatched.nEvents           = 100000

# EFT weights included + reference point
TTGamma_SingleLeptFromT_EFT_1Line_test                    = FWLiteSample.fromDAS("TTGamma_SingleLeptFromT_EFT_1Line_test", "/TTGamma_SingleLeptFromT_1Line_EFT_test_GENSIM/llechner-TTGamma_SingleLeptFromT_1Line_EFT_test_GENSIM-c9d6edd965db7012e7f57ed560889b60/USER", "phys03", dbFile=dbFile, overwrite=ov, prefix='root://hephyse.oeaw.ac.at/')
TTGamma_SingleLeptFromT_EFT_1Line_test.reweight_pkl       = os.path.join( gridpack_directory, "EFT/TTGamma_SingleLeptFromT_1Line_EFT_test/", "TTGamma_SingleLeptFromT_1Line_EFT_test_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl" )
TTGamma_SingleLeptFromT_EFT_1Line_test.xSection           = 0.1724
TTGamma_SingleLeptFromT_EFT_1Line_test.nEvents            = 100000

TTGamma_DiLept_EFT_1Line_small                            = FWLiteSample.fromDAS("TTGamma_DiLept_EFT_1Line_small", "/TTGamma_DiLept_1Line_EFT_GENSIM_small_v2/llechner-TTGamma_DiLept_1Line_EFT_GENSIM_small_v2-f2b65dd54e97daa1934ff1f579c2bf2e/USER", "phys03", dbFile=dbFile, overwrite=ov, prefix='root://hephyse.oeaw.ac.at/')
TTGamma_DiLept_EFT_1Line_small.reweight_pkl               = os.path.join( gridpack_directory, "EFT/TTGamma_DiLept_1Line_EFT/", "TTGamma_DiLept_1Line_EFT_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl" )
TTGamma_DiLept_EFT_1Line_small.xSection                   = 0.3948
TTGamma_DiLept_EFT_1Line_small.nEvents                    = 100000

TTGamma_DiLept_EFT_1Line                                  = FWLiteSample.fromDAS("TTGamma_DiLept_EFT_1Line", "/TTGamma_DiLept_1Line_EFT_GENSIM_v3/llechner-TTGamma_DiLept_1Line_EFT_GENSIM_v3-f2b65dd54e97daa1934ff1f579c2bf2e/USER", "phys03", dbFile=dbFile, overwrite=ov, prefix='root://hephyse.oeaw.ac.at/')
TTGamma_DiLept_EFT_1Line.reweight_pkl                     = os.path.join( gridpack_directory, "EFT/TTGamma_DiLept_1Line_EFT/", "TTGamma_DiLept_1Line_EFT_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl" )
TTGamma_DiLept_EFT_1Line.xSection                         = 0.3948
TTGamma_DiLept_EFT_1Line.nEvents                          = 1000000

TTGamma_SingleLeptFromT_EFT_1Line                         = FWLiteSample.fromDAS("TTGamma_SingleLeptFromT_EFT_1Line", "/TTGamma_SingleLeptFromT_1Line_EFT_GENSIM_v3/llechner-TTGamma_SingleLeptFromT_1Line_EFT_GENSIM_v3-727ec4c92973027c1a2b5c77f6a1ad97/USER", "phys03", dbFile=dbFile, overwrite=ov, prefix='root://hephyse.oeaw.ac.at/')
TTGamma_SingleLeptFromT_EFT_1Line.reweight_pkl            = os.path.join( gridpack_directory, "newEFT/TTGamma_SingleLeptFromT_1Line_EFT/", "TTGamma_SingleLeptFromT_1Line_EFT_reweight_card.pkl" )
TTGamma_SingleLeptFromT_EFT_1Line.xSection                = 0.664
TTGamma_SingleLeptFromT_EFT_1Line.nEvents                 = 1000000

TTGamma_SingleLeptFromTbar_EFT_1Line                      = FWLiteSample.fromDAS("TTGamma_SingleLeptFromTbar_EFT_1Line", "/TTGamma_SingleLeptFromTbar_1Line_EFT_GENSIM_v3/llechner-TTGamma_SingleLeptFromTbar_1Line_EFT_GENSIM_v3-8f460aa69f1845634346ff2facddea4c/USER", "phys03", dbFile=dbFile, overwrite=ov, prefix='root://hephyse.oeaw.ac.at/')
TTGamma_SingleLeptFromTbar_EFT_1Line.reweight_pkl         = os.path.join( gridpack_directory, "newEFT/TTGamma_SingleLeptFromTbar_1Line_EFT/", "TTGamma_SingleLeptFromTbar_1Line_EFT_reweight_card.pkl" )
TTGamma_SingleLeptFromTbar_EFT_1Line.xSection             = 0.6611
TTGamma_SingleLeptFromTbar_EFT_1Line.nEvents              = 1000000

#################################
### New Generation of Samples ###
#################################

# Comparing SM 1Line to central samples
# central samples are produced with central gridpack + processed with modified 71X gen cfg
# 1Line samples are produced with central gridpack production repo 93X + usual 93X cfg

# Private 1Line Production
ttGamma_SingleLeptFromT_SM_1Line                 = FWLiteSample.fromDAS("ttGamma_SingleLeptFromT_SM_1Line", "/ttGamma_SingleLeptFromT_1line_5f_ckm_LO_v2/llechner-ttGamma_SingleLeptFromT_1line_5f_ckm_LO_v2-3da16231862a69c7a8ab2bd0d654d28a/USER", "phys03", dbFile=dbFile, overwrite=ov, prefix='root://hephyse.oeaw.ac.at/')
ttGamma_SingleLeptFromT_SM_1Line.xSection        = 0.7092
ttGamma_SingleLeptFromT_SM_1Line.nEvents         = 100000

ttGamma_SingleLeptFromTbar_SM_1Line              = FWLiteSample.fromDAS("ttGamma_SingleLeptFromTbar_SM_1Line", "/ttGamma_SingleLeptFromTbar_1line_5f_ckm_LO_v2/llechner-ttGamma_SingleLeptFromTbar_1line_5f_ckm_LO_v2-d188e1ef76b4c681de518734d4a7e505/USER", "phys03", dbFile=dbFile, overwrite=ov, prefix='root://hephyse.oeaw.ac.at/')
ttGamma_SingleLeptFromTbar_SM_1Line.xSection     = 0.7068
ttGamma_SingleLeptFromTbar_SM_1Line.nEvents      = 100000

# Central Samples
ttGamma_SingleLeptFromT_SM_central               = FWLiteSample.fromDAS("ttGamma_SingleLeptFromT_SM_central", "/ttGamma_SingleLeptFromT_5f_ckm_LO_central_v3/llechner-ttGamma_SingleLeptFromT_5f_ckm_LO_central_v3-f6bd36bbeee1fb1571f3e5f5514cd426/USER", "phys03", dbFile=dbFile, overwrite=ov, prefix='root://hephyse.oeaw.ac.at/')
ttGamma_SingleLeptFromT_SM_central.xSection      = 0.7772
ttGamma_SingleLeptFromT_SM_central.nEvents       = 100000

ttGamma_SingleLeptFromTbar_SM_central            = FWLiteSample.fromDAS("ttGamma_SingleLeptFromTbar_SM_central", "/ttGamma_SingleLeptFromTbar_5f_ckm_LO_central_v3/llechner-ttGamma_SingleLeptFromTbar_5f_ckm_LO_central_v3-dda6ab68465cb0862bad77d9f61ac2c5/USER", "phys03", dbFile=dbFile, overwrite=ov, prefix='root://hephyse.oeaw.ac.at/')
ttGamma_SingleLeptFromTbar_SM_central.xSection   = 0.7749
ttGamma_SingleLeptFromTbar_SM_central.nEvents    = 100000

#TTGamma_Dilept_SM_central                        = FWLiteSample.fromDAS("ttGamma_Dilept_SM_central", "", "phys03", dbFile=dbFile, overwrite=ov, prefix='root://hephyse.oeaw.ac.at/')
#TTGamma_Dilept_SM_central.xSection               = 0.6353
#TTGamma_Dilept_SM_central.nEvents                = 100000

ttGamma_Hadronic_SM_central                      = FWLiteSample.fromDAS("ttGamma_Hadronic_SM_central", "/ttGamma_Hadronic_5f_ckm_LO_central_v3/llechner-ttGamma_Hadronic_5f_ckm_LO_central_v3-a0b185cb632020d2d97d17a4c2a505b7/USER", "phys03", dbFile=dbFile, overwrite=ov, prefix='root://hephyse.oeaw.ac.at/')
ttGamma_Hadronic_SM_central.xSection             = 0.7998
ttGamma_Hadronic_SM_central.nEvents              = 100000

SM = [
    TTGamma_SingleLeptFromT_SM_1Line_test,
    TTGamma_SingleLeptFromT_SM_3LineBuggy_test,
    TTGamma_SingleLeptFromT_SM_3LinePatched_test,
    TTGamma_SingleLeptFromT_SM_1Line,
    TTGamma_SingleLeptFromT_SM_3LinePatched,
    ttGamma_SingleLeptFromT_SM_1Line,
    ttGamma_SingleLeptFromTbar_SM_1Line,
    ttGamma_SingleLeptFromT_SM_central,
    ttGamma_SingleLeptFromTbar_SM_central,
#    ttGamma_Dilept_SM_central,
    ttGamma_Hadronic_SM_central,
]

EFT = [
    TTGamma_SingleLeptFromT_EFT_1Line_test,
    TTGamma_DiLept_EFT_1Line_small,
    TTGamma_DiLept_EFT_1Line,
    TTGamma_SingleLeptFromT_EFT_1Line,
    TTGamma_SingleLeptFromTbar_EFT_1Line,
]

allSamples = SM + EFT

for s in allSamples:
    s.isData = False




