''' Plot script WC parameter LogLikelihood
'''

from TTZRun2EFT.Tools.genCutInterpreter             import cutInterpreter
from Analysis.Tools.WeightInfo                      import WeightInfo
from TTZRun2EFT.Samples.genTuples_TTZ_postProcessed import *

genSignalSample = ttZ_ll_LO_order2_15weights_ref
#genSignalSample = ttZ_ll_LO_order3_8weights
genSignalSample.setWeightString( "ref_lumiweight1fb" )

w = WeightInfo( genSignalSample.reweight_pkl )
w.set_order( 2 )

kwargsList = [
              { "cpt":0, "cpQM":0, "res":0.07082 },

              { "cpt":-24.5, "cpQM":-8, "res":0.1497 },
              { "cpt":-24.5, "cpQM":4, "res":0.1093 },
              { "cpt":-24.5, "cpQM":8, "res":0.1104 },

              { "cpt":-14,  "cpQM":4, "res":0.05106 },
              { "cpt":-3.5, "cpQM":4, "res":0.04295 },

              { "cpt":7, "cpQM":-4, "res":0.1386 },
              { "cpt":7, "cpQM":8, "res":0.06936 },

              { "cpt":-24.5, "cpQM":0, "res":0.1155 },
              { "cpt":-14, "cpQM":0, "res":0.06286 },
              { "cpt":-8, "cpQM":0, "res":0.05531 },
              { "cpt":-7, "cpQM":0, "res":0.05563 },
              { "cpt":-6, "cpQM":0, "res":0.05647 },
              { "cpt":-5, "cpQM":0, "res":0.05775 },
              { "cpt":-4, "cpQM":0, "res":0.05947 },
              { "cpt":-3, "cpQM":0, "res":0.06165 },
              { "cpt":-2, "cpQM":0, "res":0.06736 },
              { "cpt":-1, "cpQM":0, "res":0.06425 },
              { "cpt":1, "cpQM":0, "res":0.07475 },
              { "cpt":2, "cpQM":0, "res":0.07922 },
              { "cpt":3, "cpQM":0, "res":0.08411 },
              { "cpt":4, "cpQM":0, "res":0.08935 },
              { "cpt":5, "cpQM":0, "res":0.09526 },
              { "cpt":6, "cpQM":0, "res":0.1016 },
              { "cpt":7, "cpQM":0, "res":0.1083 },
              { "cpt":8, "cpQM":0, "res":0.1154 },

              { "cpt":0, "cpQM":-6, "res":0.1135 },
              { "cpt":0, "cpQM":-5, "res":0.1054 },
              { "cpt":0, "cpQM":-4, "res":0.09751 },
              { "cpt":0, "cpQM":-3, "res":0.09008 },
              { "cpt":0, "cpQM":-2, "res":0.08319 },
              { "cpt":0, "cpQM":-1, "res": 0.07679 },
              { "cpt":0, "cpQM":1, "res":0.06538 },
              { "cpt":0, "cpQM":2, "res":0.06019 },
              { "cpt":0, "cpQM":3, "res":0.05566 },
              { "cpt":0, "cpQM":4, "res":0.05146 },
              { "cpt":0, "cpQM":5, "res":0.04776 },
              { "cpt":0, "cpQM":6, "res":0.04453 },
              { "cpt":0, "cpQM":7, "res":0.04173 },
              { "cpt":0, "cpQM":8, "res":0.03942 },
             ]
genSelection = "abs(genZ_mass-91.2)<10&&Sum$(genLep_pt>15&&abs(genLep_eta)<2.4)>=2"
regionCut    = "(1)" #replaceAliases( simpleStringToCutString( cut ) ) 

sel       = genSelection #"&&".join( [ cutInterpreter.cutString( genSelection ), regionCut ] )
coeffList = w.getCoeffListFromDraw( genSignalSample, selectionString=sel )

signal_genRateSM  = float( w.get_weight_yield( coeffList ) )
dSignal_genRateSM  = 0.07082

for kwargs in kwargsList:
    dRes = kwargs["res"]
    del kwargs["res"]
    signal_genRateEFT = w.get_weight_yield( coeffList, **kwargs )

    print "BSM point:", kwargs
    print "SM:", signal_genRateSM, "Daniel:", dSignal_genRateSM
    print "BSM:", signal_genRateEFT, "Daniel:", dRes
    print "ratio:", signal_genRateEFT / signal_genRateSM, "Daniel:", dRes / dSignal_genRateSM, "difference:", (signal_genRateEFT / signal_genRateSM) / (dRes / dSignal_genRateSM)
    print

