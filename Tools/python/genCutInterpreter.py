''' Class to interpret string based cuts
'''

# TTGamma Imports
from Analysis.Tools.CutInterpreter import CutInterpreter

mZ = 91.1876

special_cuts = {
    "onZll":             "abs(mll-%s)<=15"%(mZ),
  }

continous_variables = [ ("met", "GenMET_pt"), ("pTG","GenPhoton_pt[0]"), ("mll", "mll"), ("mllgamma", "mllgamma") ]
discrete_variables  = [ ("nJet", "nGenJet"), ("nBTag", "nGenBJet"), ("nLep","nGenLepton"), ("nPhoton","nGenPhoton") ]

cutInterpreter = CutInterpreter( continous_variables, discrete_variables, special_cuts)

if __name__ == "__main__":
    print cutInterpreter.cutString("dilepOS-pTG20-nPhoton1p-offZSFll-offZSFllg-mll40")

# Bin0: all Z_pt0To100_cosThetaStar-1
# Bin1: all Z_pt0To100_cosThetaStar-0.6To0.6
# Bin2: all Z_pt0To100_cosThetaStar0.6To1
# Bin3: all Z_pt100To200_cosThetaStar-1
# Bin4: all Z_pt100To200_cosThetaStar-0.6To0.6
# Bin5: all Z_pt100To200_cosThetaStar0.6To1
# Bin6: all Z_pt200To400_cosThetaStar-1
# Bin7: all Z_pt200To400_cosThetaStar-0.6To0.6
# Bin8: all Z_pt200To400_cosThetaStar0.6To1
# Bin9: all Z_pt400_cosThetaStar-1
# Bin10: all Z_pt400_cosThetaStar-0.6To0.6
# Bin11: all Z_pt400_cosThetaStar0.6To1
# Bin12: all Z1_pt_4l0To100_Z1_cosThetaStar_4l-1To1
# Bin13: all Z1_pt_4l100To200_Z1_cosThetaStar_4l-1To1
# Bin14: all Z1_pt_4l200_Z1_cosThetaStar_4l-1To1
# Bin15: all Z_pt0To100_cosThetaStar-1
# Bin16: all Z_pt0To100_cosThetaStar-0.6To0.6
# Bin17: all Z_pt0To100_cosThetaStar0.6To1
# Bin18: all Z_pt100To200_cosThetaStar-1
# Bin19: all Z_pt100To200_cosThetaStar-0.6To0.6
# Bin20: all Z_pt100To200_cosThetaStar0.6To1
# Bin21: all Z_pt200To400_cosThetaStar-1
# Bin22: all Z_pt200To400_cosThetaStar-0.6To0.6
# Bin23: all Z_pt200To400_cosThetaStar0.6To1
# Bin24: all Z_pt400_cosThetaStar-1
# Bin25: all Z_pt400_cosThetaStar-0.6To0.6
# Bin26: all Z_pt400_cosThetaStar0.6To1
# Bin27: all Z1_pt_4l0To100_Z1_cosThetaStar_4l-1To1
# Bin28: all Z1_pt_4l100To200_Z1_cosThetaStar_4l-1To1
# Bin29: all Z1_pt_4l200_Z1_cosThetaStar_4l-1To1

