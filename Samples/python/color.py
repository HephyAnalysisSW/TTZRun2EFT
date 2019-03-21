import ROOT

from TTGammaEFT.Samples.helpers import singleton as singleton

@singleton
class color():
  pass

color.data           = ROOT.kBlack
color.DY             = ROOT.kCyan+2
color.ZGamma         = ROOT.kBlue+2
color.WGamma         = ROOT.kBlue
color.TTG            = ROOT.kOrange
color.Other          = ROOT.kViolet
color.TT             = ROOT.kRed+1
color.T              = ROOT.kOrange+1
color.TGamma         = ROOT.kGray
color.W              = ROOT.kCyan+1

color.TTG1L          = ROOT.kBlue+2
color.TTG3LBuggy     = ROOT.kRed+1
color.TTG3LPatched   = ROOT.kGreen+3
