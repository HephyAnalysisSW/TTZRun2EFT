#cardfile="test"
#cardfile="TTG_DiLept_1L_EFT_COMBINED_1l_ctZ_0_ctZI_0_nLepTight1-nLepVeto1-nJet4p-nBTag1p-pTG20-nPhoton1p_full"
#cardfile="TTG_DiLept_1L_EFT_TTG_DiLept_1L_EFT_COMBINED_2l_ctZ_0_ctZI_0_nLepTight1-nLepVeto1-nJet4p-nBTag1p-pTG20-nPhoton1p_dilepOS-nLepVeto2-pTG20-nPhoton1p-offZSFllg-offZSFll-mll40-nJet2p-nBTag1p_full"
cardfile="TTG_DiLept_1L_EFT_TTG_DiLept_1L_EFT_COMBINED_2l_ctZ_0_ctZI_0_nLepTight1-nLepVeto1-nJet4p-nBTag1p-pTG20-nPhoton1p_dilepOS-nLepVeto2-pTG20-nPhoton1p-offZSFllg-offZSFll-mll40-nJet2p-nBTag1p_full_low"
#cardfile="TTG_DiLept_1L_EFT_TTG_DiLept_1L_EFT_COMBINED_2l_ctZ_0_ctZI_0_nLepTight1-nLepVeto1-nJet4p-nBTag1p-pTG20-nPhoton1p_dilepOS-nLepVeto2-pTG20-nPhoton1p-offZSFllg-offZSFll-mll40-nJet2p-nBTag1p_full_high"
#TTG_DiLept_1L_EFT_2016_ctZ_0_ctZI_0_dilepOS-nLepVeto2-pTG20-nPhoton1p-offZSFllg-offZSFll-mll40-nJet2p-nBTag1p_small"
#cardfile=""
#cardfile="TTG_DiLept_1L_EFT_2016_ctZ_0_ctZI_0_dilepOS-nLepVeto2-pTG20-nPhoton1p-offZSFllg-offZSFll-mll40_full"

python impactPlot.py --removeDir --cardfile ${cardfile} --cores 6 --year combined
