import os

if os.environ['USER'] in ['schoef', 'rschoefbeck', 'schoefbeck']:
    results_directory               = "/afs/hephy.at/data/rschoefbeck02/TTZRun2EFT/results/"
    skim_output_directory           = "/afs/hephy.at/data/rschoefbeck02/TTZRun2EFT/skims/"
    skim_directory                  = "/afs/hephy.at/data/dspitzbart01/TTZRun2EFT/skims/"
    tmp_directory                   = "/afs/hephy.at/data/rschoefbeck02/TTZRun2EFT_tmp/"
    plot_directory                  = "/afs/hephy.at/user/r/rschoefbeck/www/TTZRun2EFT/"
    data_directory                  = "/afs/hephy.at/data/rschoefbeck01/cmgTuples/"
    postprocessing_directory        = "TTZRun2EFT_PP_2016_TTG_v5/inclusive/"
    postprocessing_output_directory = "/afs/hephy.at/data/rschoefbeck02/cmgTuples/"
    analysis_results                = results_directory

if os.environ['USER'] in ['llechner']:
    tmp_directory                       = "/afs/hephy.at/data/llechner01/tmp/"
    results_directory                   = "/afs/hephy.at/data/llechner01/TTZRun2EFT/results/"

    plot_directory                      = "/afs/hephy.at/user/l/llechner/www/TTZRun2EFT/"
    postprocessing_output_directory     = "/afs/hephy.at/data/llechner03/TTZRun2EFT/nanoTuples/"

    cache_directory                     = "/afs/hephy.at/data/llechner01/TTZRun2EFT/cache/"
    combineReleaseLocation              = '/afs/hephy.at/user/l/llechner/public/CMSSW_8_1_0/src'
    cardfileLocation                    = '/afs/hephy.at/data/llechner01/TTZRun2EFT/results/cardfiles/'

    dpm_directory                       = '/dpm/oeaw.ac.at/home/cms/store/user/llechner/'

    data_directory1                     = "/afs/hephy.at/data/llechner01/TTGammaEFT/nanoTuples/"
    data_directory2                     = "/afs/hephy.at/data/llechner02/TTGammaEFT/nanoTuples/"
    data_directory3                     = "/afs/hephy.at/data/llechner03/TTGammaEFT/nanoTuples/"

    postprocessing_directoryGEN         = "TTGammaEFT_PP_GEN_TTG_v3/gen/"

    gridpack_directory                  = "/afs/hephy.at/data/llechner01/TTGammaEFT/gridpacks/"
