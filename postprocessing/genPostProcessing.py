#!/usr/bin/env python
''' Make flat ntuple from GEN data tier 
'''
# Standard imports and batch mode
import ROOT
import os, sys
ROOT.gROOT.SetBatch(True)
import itertools
from math                                import sqrt, cos, sin, pi, acos
import imp

#RootTools
from RootTools.core.standard             import *

#TTZRun2EFT
from TTZRun2EFT.Tools.user               import postprocessing_output_directory
from TTZRun2EFT.Tools.objectSelection    import isGoodGenJet, isGoodGenLepton, isGoodGenPhoton, genJetId

from Analysis.Tools.GenSearch            import GenSearch
from Analysis.Tools.helpers              import deltaPhi, deltaR, deltaR2, cosThetaStar, closestOSDLMassToMZ
from Analysis.Tools.HyperPoly            import HyperPoly
from Analysis.Tools.WeightInfo           import WeightInfo

# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',           action='store',      default='INFO',          nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--small',              action='store_true', help='Run only on a small subset of the data?')#, default = True)
argParser.add_argument('--overwrite',          action='store',      nargs='?', choices = ['none', 'all', 'target'], default = 'none', help='Overwrite?')#, default = True)
argParser.add_argument('--targetDir',          action='store',      default='v5')
argParser.add_argument('--sample',             action='store',      default='fwlite_ttZ_ll_LO_scan', help="Name of the sample loaded from fwlite_benchmarks. Only if no inputFiles are specified")
argParser.add_argument('--inputFiles',         action='store',      nargs = '*', default=[])
argParser.add_argument('--nJobs',              action='store',      nargs='?', type=int, default=1,  help="Maximum number of simultaneous jobs.")
argParser.add_argument('--job',                action='store',      nargs='?', type=int, default=0,  help="Run only job i")
argParser.add_argument('--addReweights',       action='store_true',   help="Add reweights?")
argParser.add_argument('--interpolationOrder', action='store',      nargs='?', type=int, default=3,  help="Interpolation order for EFT weights.")
args = argParser.parse_args()

# Logger
import Analysis.Tools.logger as _logger
import RootTools.core.logger as _logger_rt
logger    = _logger.get_logger(   args.logLevel, logFile = None)
logger_rt = _logger_rt.get_logger(args.logLevel, logFile = None)

sample_file = "$CMSSW_BASE/python/TTZRun2EFT/Samples/genTuples_TTZ.py"
samples = imp.load_source( "samples", os.path.expandvars( sample_file ) )
sample = getattr( samples, args.sample )
logger.debug( 'Loaded sample %s with %i files.', sample.name, len(sample.files) )

maxEvents = -1
if args.small: 
    args.targetDir += "_small"
    maxEvents=100 
    sample.files=sample.files[:1]

xsec          = sample.xsec
nEvents       = sample.nEvents
lumiweight1fb = xsec * 1000. / nEvents

# output directory
output_directory = os.path.join(postprocessing_output_directory, args.targetDir, "gen",  sample.name) 

if not os.path.exists( output_directory ): 
    os.makedirs( output_directory )
    logger.info( "Created output directory %s", output_directory )

# Load reweight pickle file if supposed to keep weights. 
extra_variables = []
if args.addReweights:
    weightInfo = WeightInfo( sample.reweight_pkl )
    # Determine coefficients for storing in vector
    # Sort Ids wrt to their position in the card file

    # weights for the ntuple
    rw_vector       = TreeVariable.fromString( "rw[w/F,"+",".join(w+'/F' for w in weightInfo.variables)+"]")
    rw_vector.nMax  = weightInfo.nid
    extra_variables.append(rw_vector)

    # coefficients for the weight parametrization
    param_vector      = TreeVariable.fromString( "p[C/F]" )
    param_vector.nMax = HyperPoly.get_ndof(weightInfo.nvar, args.interpolationOrder)
    hyperPoly         = HyperPoly( args.interpolationOrder )
    extra_variables.append(param_vector)
    extra_variables.append(TreeVariable.fromString( "chi2_ndof/F"))

def interpret_weight(weight_id):
    str_s = weight_id.split('_')
    res={}
    for i in range(len(str_s)/2):
        res[str_s[2*i]] = float(str_s[2*i+1].replace('m','-').replace('p','.'))
    return res

# Run only job number "args.job" from total of "args.nJobs"
if args.nJobs>1:
    n_files_before = len(sample.files)
    sample = sample.split(args.nJobs)[args.job]
    n_files_after  = len(sample.files)
    logger.info( "Running job %i/%i over %i files from a total of %i.", args.job, args.nJobs, n_files_after, n_files_before)

products = {
    'lhe':{'type':'LHEEventProduct', 'label':("externalLHEProducer")},
    'gp':{'type':'vector<reco::GenParticle>', 'label':("genParticles")},
    'genJets':{'type':'vector<reco::GenJet>', 'label':("ak4GenJets")},
    'genMET':{'type':'vector<reco::GenMET>',  'label':("genMetTrue")},
}

def varnames( vec_vars ):
    return [v.split('/')[0] for v in vec_vars.split(',')]

def vecSumPt(*args):
    return sqrt( sum([o['pt']*cos(o['phi']) for o in args],0.)**2 + sum([o['pt']*sin(o['phi']) for o in args],0.)**2 )

def addIndex( collection ):
    for i  in range(len(collection)):
        collection[i]['index'] = i

# standard variables
variables  = ["run/I", "lumi/I", "evt/l"]

# MET
variables += ["genMet_pt/F", "genMet_phi/F"]
# jet vector
jet_read_vars       =  "pt/F,eta/F,phi/F,isMuon/I,isElectron/I,isPhoton/I"
jet_read_varnames   =  varnames( jet_read_vars )
jet_write_vars      = jet_read_vars+',matchBParton/I' 
jet_write_varnames  =  varnames( jet_write_vars )
variables += ["genJet[%s]"%jet_write_vars]
variables += ["genBj0_%s"%var for var in jet_write_vars.split(',')]
variables += ["genBj1_%s"%var for var in jet_write_vars.split(',')]
# lepton vector 
lep_vars       =  "pt/F,eta/F,phi/F,pdgId/I"
lep_extra_vars =  "motherPdgId/I,grandmotherPdgId/I"
lep_varnames   =  varnames( lep_vars ) 
lep_all_varnames = lep_varnames + varnames(lep_extra_vars)
variables     += ["genLep[%s]"%(','.join([lep_vars, lep_extra_vars]))]
# associated jet indices
variables += [ "genBjLeadlep_index/I", "genBjLeadhad_index/I" ]
variables += [ "genBjNonZlep_index/I", "genBjNonZhad_index/I" ]
# top vector
top_vars       =  "pt/F,eta/F,phi/F,pdgId/I,mass/F"
top_varnames   =  varnames( top_vars ) 
variables     += ["genTop[%s]"%top_vars]

# to be stored for each boson
boson_read_varnames= [ 'pt', 'phi', 'eta', 'mass']
# Z vector from gen collection
variables     += ["genZ_pt/F", "genZ_phi/F", "genZ_eta/F", "genZ_mass/F", "genZ_cosThetaStar/F", "genZ_daughterPdg/I"]
# Z vector from genleps
variables     += ["genLepZ_pt/F", "genLepZ_phi/F", "genLepZ_eta/F", "genLepZ_mass/F", "genLepZ_lldPhi/F", "genLepZ_lldR/F","genLepZ_cosThetaStar/F", "genLepZ_daughterPdg/I", "genLepNonZ_l1_index/I"]
variables     += ["genLepZ_l1_index/I", "genLepZ_l2_index/I", "genLepNonZ_l1_index/I", "genLepNonZ_l2_index/I"]
# W vector
variables     += ["genW_pt/F", "genW_phi/F", "genW_eta/F", "genW_mass/F", "genW_daughterPdg/I"]
# gamma vector
gen_photon_vars = "pt/F,phi/F,eta/F,mass/F,motherPdgId/I,relIso04/F,minLeptonDR/F,minJetDR/F"
variables     += ["genPhoton[%s]"%gen_photon_vars]
gen_photon_varnames = varnames( gen_photon_vars )

# Lumi weight 1fb
variables += ["lumiweight1fb/F"]

variables += ["signalZ/I"] #ttZ Signal

variables += ["signalPhoton/I"] #cat a1
variables += ["isrPhoton/I"] #cat a2
variables += ["lepPhoton/I"] #cat b
variables += ["nonIsoPhoton/I"] #cat c1
variables += ["jetPhoton/I"] #cat c2
variables += ["fakePhoton/I"] #cat d

if args.addReweights:
    variables.append('rw_nominal/F')
    # Lumi weight 1fb / w_0
    variables.append("ref_lumiweight1fb/F")


def fill_vector_collection( event, collection_name, collection_varnames, objects):
    setattr( event, "n"+collection_name, len(objects) )
    for i_obj, obj in enumerate(objects):
        for var in collection_varnames:
            getattr(event, collection_name+"_"+var)[i_obj] = obj[var]
def fill_vector( event, collection_name, collection_varnames, obj):
    for var in collection_varnames:
        try:
            setattr(event, collection_name+"_"+var, obj[var] )
        except TypeError as e:
            logger.error( "collection_name %s var %s obj[var] %r", collection_name, var,  obj[var] )
            raise e
        except KeyError as e:
            logger.error( "collection_name %s var %s obj[var] %r", collection_name, var,  obj[var] )
            raise e

reader = sample.fwliteReader( products = products )

def filler( event ):

    event.run, event.lumi, event.evt = reader.evt
    event.lumiweight1fb = lumiweight1fb

    if reader.position % 100==0: logger.info("At event %i/%i", reader.position, reader.nEvents)

    if args.addReweights:
        event.nrw = weightInfo.nid
        lhe_weights = reader.products['lhe'].weights()
        weights      = []
        param_points = []

        for weight in lhe_weights:
            # Store nominal weight (First position!) 
            if weight.id=='rwgt_1': event.rw_nominal = weight.wgt
            if not weight.id in weightInfo.id: continue
            pos = weightInfo.data[weight.id]
            event.rw_w[pos] = weight.wgt
            weights.append( weight.wgt )
            interpreted_weight = interpret_weight(weight.id) 
            for var in weightInfo.variables:
                getattr( event, "rw_"+var )[pos] = interpreted_weight[var]
            # weight data for interpolation
            if not hyperPoly.initialized: param_points.append( tuple(interpreted_weight[var] for var in weightInfo.variables) )

        # get list of values of ref point in specific order
        ref_point_coordinates = [weightInfo.ref_point_coordinates[var] for var in weightInfo.variables]

        # Initialize with Reference Point
        if not hyperPoly.initialized: hyperPoly.initialize( param_points, ref_point_coordinates )
        coeff = hyperPoly.get_parametrization( weights )

        # = HyperPoly(weight_data, args.interpolationOrder)
        event.np = hyperPoly.ndof
        event.chi2_ndof = hyperPoly.chi2_ndof(coeff, weights)
        #logger.debug( "chi2_ndof %f coeff %r", event.chi2_ndof, coeff )
        if event.chi2_ndof>10**-6: logger.warning( "chi2_ndof is large: %f", event.chi2_ndof )
        for n in xrange(hyperPoly.ndof):
            event.p_C[n] = coeff[n]

        # lumi weight / w0
        event.ref_lumiweight1fb = event.lumiweight1fb / coeff[0]

    # All gen particles
    gp      = reader.products['gp']

    # for searching
    search  = GenSearch( gp )

    # find heavy objects before they decay
    genTops = map( lambda t:{var: getattr(t, var)() for var in top_varnames}, filter( lambda p:abs(p.pdgId())==6 and search.isLast(p),  gp) )

    genTops.sort( key = lambda p:-p['pt'] )
    fill_vector_collection( event, "genTop", top_varnames, genTops ) 

    # generated Z's
    genZs = filter( lambda p:abs(p.pdgId())==23 and search.isLast(p) and abs(p.daughter(0).pdgId()) in [11, 13], gp)
    genZs.sort( key = lambda p: -p.pt() )
    if len(genZs)>0: 
        genZ = genZs[0]
        for var in boson_read_varnames:
           setattr( event, "genZ_"+var,  getattr(genZ, var)() )
    else:
        genZ = None
    
    event.signalZ = 0    #ttZ with Z from gluon or top
    if genZ is not None:

        if abs(search.ascend(genZ).mother(0).pdgId()) in [ 6, 21 ]:
            event.signalZ = 1    #ttZ with Z from gluon or top

        d1, d2 = genZ.daughter(0), genZ.daughter(1)
        if d1.pdgId()>0: 
            lm, lp = d1, d2
        else:
            lm, lp = d2, d1
        event.genZ_daughterPdg = lm.pdgId()
        event.genZ_cosThetaStar = cosThetaStar(genZ.mass(), genZ.pt(), genZ.eta(), genZ.phi(), lm.pt(), lm.eta(), lm.phi())

    # generated W's
    genWs = filter( lambda p:abs(p.pdgId())==24 and search.isLast(p), gp)
    genWs.sort( key = lambda p: -p.pt() )
    # W can't have a top-mother - We're looking for the extra boson (there is always an extra boson)
    genWs = filter( lambda p: abs(search.ascend(p).mother(0).pdgId())!=6, genWs )
    if len(genWs)>0: 
        genW = genWs[0]
        for var in boson_read_varnames:
           setattr( event, "genW_"+var,  getattr(genW, var)() )
    else:
        genW = None
    
    if genW is not None:
        d1, d2 = genW.daughter(0), genW.daughter(1)
        if abs(d1.pdgId()) in [11, 13, 15]: 
            lep, neu = d1, d2
        else:
            lep, neu = d2, d1

        event.genW_daughterPdg = lep.pdgId()


    def printTopMothers():
            genTopCheck = [ (search.ascend(l), l, search.ancestry( search.ascend(l) )) for l in filter( lambda p:abs(p.pdgId())==6,  gp) ]

            print genTopCheck
            genTopCheck.sort( key = lambda p: -p[1].pt() )
            if len(genTopCheck) > 0:
                topPdg_ids = filter( lambda p:p!=2212, [abs(particle.pdgId()) for particle in genTopCheck[0][2]])
                print 'TOP 1'
                print topPdg_ids
                previous = genTopCheck[0][0].mother(1)
                print 'first', genTopCheck[0][0].pdgId()
                for i, pdg in enumerate(topPdg_ids):
                    try:
                        print 'mother', i, previous.pdgId()
                        print 'mothers', i, [p.pdgId() for p in search.ancestry( previous )]
                        previous = previous.mother(0)
                    except Exception as val:
                        print val
                        break
            if len(genTopCheck) > 1:
                topPdg_ids = filter( lambda p:p!=2212, [abs(particle.pdgId()) for particle in genTopCheck[1][2]])
                print 'TOP 2'
                print topPdg_ids
                previous = genTopCheck[1][0].mother(1)
                print 'first', genTopCheck[1][0].pdgId()
                for i, pdg in enumerate(topPdg_ids):
                    try:
                        print 'mother', i, previous.pdgId()
                        print 'mothers', i, [p.pdgId() for p in search.ancestry( previous )]
                        previous = previous.mother(0)
                    except Exception as val:
                        print val
                        break

    # MET
    genMet = {'pt':reader.products['genMET'][0].pt(), 'phi':reader.products['genMET'][0].phi()}
    event.genMet_pt  = genMet['pt']
    event.genMet_phi = genMet['phi'] 


#    for ttgamma and tt events, categorize events:
#        a1) ttgamma vertex, gamma from t/g, isolated, must be in ttgamma sample   = ttgamma signal (photonSignal == 1)
#            -> gamma isolated
#            -> gamma with only tops/gluons in ancestry
#        a2) ttgamma, gamma from ISR, must be in ttgamma sample                    = tt ISR bg (photonISR == 1)
#            -> gamma isolated
#            -> gamma with no top in ancestry
#            -> gamma with only gluons and light quarks in ancestry
#            -> direct gamma mother != gluon!!! (this is a1)
#        b) tt + isolated gamma from W, l, tau, must be in tt sample               = ttgamma (W,l,tau) bg (photonLep == 1)
#            -> gamma isolated
#            -> gamma with direct abs mother being 11, 13, 15 or 24 
#        c1) tt + non isolated gamma from anything including mesons                = tt bg (ttBg == 1)
#            -> gamma non isolated or meson in ancestry
#        c2) tt + isolated gamma from bottom quark (from t decay) or jets (from W) = tt bottom bg (ttBg == 1)
#            -> gamma isolated from b or j (from t/W decay)
#            ATTENTION: gammas from bottoms with off-shell tops where the
#                       top decay is done by pythia are currently labeled as ISR!
#                       However this case is currently not simulated in any sample 
#        d) tt + gamma fake                                                        = ttgamma fake (photonFake == 1)
#            -> everything else does not contain a photon
#            -> if it still passes selection: it is a photon fake

    genPhotonsSignalCheck = [ (search.ascend(l), l, search.ancestry( search.ascend(l) )) for l in filter( lambda p:abs(p.pdgId())==22 and p.pt()>10 and abs(p.eta())<2.5 and search.isLast(p) and p.status()==1, gp) ]
    genPhotonsSignalCheck.sort( key = lambda p: -p[1].pt() )

    photonSignal = 0    #a1
    photonISR = 0       #a2
    photonLep = 0       #b
    ttBg = 0            #c1
    photonJets = 0      #c2
    photonFake = 0      #d


    if len( genPhotonsSignalCheck ) > 0:     
        # check hardest photon with pT>13 and abs(eta)<2.5
        first, last, ancestry = genPhotonsSignalCheck[0]
        # get abs pdgIDs of ancestry
        pdg_ids = filter( lambda p:p!=2212, [abs(particle.pdgId()) for particle in ancestry])
        # check if particles are close by
        close_particles = filter( lambda p: p!=last and p.pt()>5 and deltaR2( {'phi':last.phi(), 'eta':last.eta()}, {'phi':p.phi(), 'eta':p.eta()} )<0.2**2 , search.final_state_particles_no_neutrinos )
    
#        print 'mothers pdg', pdg_ids
#        print 'close', [p.pdgId() for p in close_particles]
#        print 'first mother pdg', first.mother(0).pdgId()

#        if len(pdg_ids) < 999:
#            previous = first.mother(0)
#            for i, pdg in enumerate(pdg_ids):
#                try:
#                    print 'mother', i, previous.pdgId()
#                    previous = previous.mother(0)
#                except Exception as val:
#                    print val
#                    break

        # deside the categories
        if max(pdg_ids)>100 or len(close_particles) != 0:
            #photon with meson in ancestry or non isolated -> cat c1)
            ttBg = 1
        elif abs(first.mother(0).pdgId()) in [ 11, 13, 15, 24 ]:
            #isolated photon with W, l or tau direct mother -> cat b)
            photonLep = 1
#        elif all( [ p in [ 6, 21 ] for p in pdg_ids ] ) or abs(first.mother(0).pdgId()) == 21: # not working for photons from top, as the gluons can come from light quarks
        elif abs(first.mother(0).pdgId()) in [ 6, 21 ]:
            #isolated photon with photon from top or gluon -> cat a1)
            photonSignal = 1
#            printTopMothers()
        elif all( [ p in [ 1, 2, 3, 4, 5, 21 ] for p in pdg_ids ] ):
            #isolated photon with photon ancestry only containing light quarks or gluons (ISR) -> cat a1)
            photonISR = 1
        else:
            #isolated gammas from bottoms originating from the top decay or jets from W -> cat c2) 
            photonJets = 1

    else:
        # if events with photonFake == 1 pass selection: fake gamma -> cat d)
        photonFake = 1

    # if all flags are 0, it is an isolated gamma from a process I havn't thought of!
    # should not be there! - check!
    event.signalPhoton = photonSignal
    event.isrPhoton    = photonISR
    event.lepPhoton    = photonLep
    event.nonIsoPhoton = ttBg
    event.jetPhoton    = photonJets
    event.fakePhoton   = photonFake

    # gen photons: particle-level isolated gen photons
    genPhotons = [ (search.ascend(l), l) for l in filter( lambda p:abs(p.pdgId())==22 and p.pt()>15 and search.isLast(p) and p.status()==1, gp) ]
    genPhotons.sort( key = lambda p: -p[1].pt() )
    genPhotons_   = []

    for first, last in genPhotons[:100]: 
        mother_pdgId = first.mother(0).pdgId() if first.numberOfMothers()>0 else 0
        genPhoton_ = {var: getattr(last, var)() for var in boson_read_varnames}
        # kinematic photon selection
        if not isGoodGenPhoton( genPhoton_): continue
        genPhoton_['motherPdgId'] = mother_pdgId
        genPhoton_['status']      = last.status()
        
        close_particles = filter( lambda p: p!=last and deltaR2( {'phi':last.phi(), 'eta':last.eta()}, {'phi':p.phi(), 'eta':p.eta()} )<0.4**2 , search.final_state_particles_no_neutrinos )
        genPhoton_['relIso04'] = sum( [p.pt() for p in close_particles], 0) / last.pt()
        # require isolation
        if genPhoton_['relIso04']<0.4:
            genPhotons_.append( genPhoton_ )
    
    # genLeptons: prompt gen-leptons 
    genLeptons = [ (search.ascend(l), l) for l in filter( lambda p:abs(p.pdgId()) in [11, 13] and search.isLast(p) and p.pt()>=0 and p.status()==1,  gp) ]
    promptGenLeps    = []
    allGenLeps    = []
    for first, last in genLeptons:
        mother = first.mother(0) if first.numberOfMothers()>0 else None
        if mother is not None:
            mother_pdgId      = mother.pdgId()
            mother_ascend     = search.ascend(mother)
            grandmother       = mother_ascend.mother(0) if mother.numberOfMothers()>0 else None
            grandmother_pdgId = grandmother.pdgId() if grandmother is not None else 0
        else:
            mother_pdgId = 0
            grandmother_pdgId = 0 
        genLep = {var: getattr(last, var)() for var in lep_varnames}
        genLep['motherPdgId']      = mother_pdgId
        genLep['grandmotherPdgId'] = grandmother_pdgId
        allGenLeps.append( genLep )
        if abs(genLep['motherPdgId']) in [ 11, 13, 15, 23, 24, 25 ]:
            promptGenLeps.append(genLep )

    # filter gen leptons
    promptGenLeps =  list( filter( lambda l:isGoodGenLepton( l ), promptGenLeps ) )
    promptGenLeps.sort( key = lambda p:-p['pt'] )
    addIndex( promptGenLeps )

    ## removing photons in dR cone leptons (radiation photons)
    for genPhoton in genPhotons_:
        genPhoton['minLeptonDR'] =  min([999]+[deltaR(genPhoton, l) for l in allGenLeps])
    genPhotons_ = list(filter( lambda g: g['minLeptonDR']>0.4, genPhotons_))
    addIndex( genPhotons_ )

    # jets
    fwlite_genJets = filter( genJetId, reader.products['genJets'] )
    genJets = map( lambda t:{var: getattr(t, var)() for var in jet_read_varnames}, filter( lambda j:j.pt()>30, fwlite_genJets) )
    # filter genJets
    genJets = list( filter( lambda j:isGoodGenJet( j ), genJets ) )
    # cleaning of jets with isolated photons
    genJets = list( filter( lambda j:min([999]+[deltaR2(j, p) for p in genPhotons_ ])>0.4**2, genJets))

    # store minimum DR to jets
    for genPhoton in genPhotons_:
        genPhoton['minJetDR'] =  min([999]+[deltaR(genPhoton, j) for j in genJets])

    # find b's from tops:
    b_partons = [ b for b in filter( lambda p:abs(p.pdgId())==5 and p.numberOfMothers()==1 and abs(p.mother(0).pdgId())==6,  gp) ]

    # store if gen-jet is DR matched to a B parton
    for genJet in genJets:
        genJet['matchBParton'] = ( min([999]+[deltaR2(genJet, {'eta':b.eta(), 'phi':b.phi()}) for b in b_partons]) < 0.2**2 )

    genJets = filter( lambda j: (min([999]+[deltaR2(j, l) for l in promptGenLeps if l['pt']>10]) > 0.3**2 ), genJets )
    genJets.sort( key = lambda p:-p['pt'] )
    addIndex( genJets )

    # gen b jets
    trueBjets = list( filter( lambda j: j['matchBParton'], genJets ) )
    trueNonBjets = list( filter( lambda j: not j['matchBParton'], genJets ) )

    # Mimick b reconstruction ( if the trailing b fails acceptance, we supplement with the leading non-b jet ) 
    genBj0, genBj1 = ( trueBjets + trueNonBjets + [None, None] )[:2]
    if genBj0: fill_vector( event, "genBj0", jet_write_varnames, genBj0) 
    if genBj1: fill_vector( event, "genBj1", jet_write_varnames, genBj1) 

    # reco-bjet/leading lepton association
    if len(promptGenLeps)>0 and genBj0 and genBj1:
        if vecSumPt( genBj0, promptGenLeps[0], genMet ) > vecSumPt( genBj1, promptGenLeps[0], genMet ):
            event.genBjLeadlep_index, event.genBjLeadhad_index = genBj0['index'], genBj1['index']
        else:
            event.genBjLeadlep_index, event.genBjLeadhad_index = genBj1['index'], genBj0['index']

    # find Z in genLep
    (event.genLepZ_mass, genLepZ_l1_index, genLepZ_l2_index) = closestOSDLMassToMZ(promptGenLeps)
    genLepNonZ_indices = [ i for i in range(len(promptGenLeps)) if i not in [genLepZ_l1_index, genLepZ_l2_index] ]
    event.genLepZ_l1_index    = promptGenLeps[genLepZ_l1_index]['index'] if genLepZ_l1_index>=0 else -1
    event.genLepZ_l2_index    = promptGenLeps[genLepZ_l2_index]['index'] if genLepZ_l2_index>=0 else -1
    event.genLepNonZ_l1_index = promptGenLeps[genLepNonZ_indices[0]]['index'] if len(genLepNonZ_indices)>0 else -1
    event.genLepNonZ_l2_index = promptGenLeps[genLepNonZ_indices[1]]['index'] if len(genLepNonZ_indices)>1 else -1
    # store genLepZ stuff
    if event.genLepZ_mass>0:
        genLepZ_l1 = ROOT.TLorentzVector()
        genLepZ_l1.SetPtEtaPhiM(promptGenLeps[event.genLepZ_l1_index]['pt'], promptGenLeps[event.genLepZ_l1_index]['eta'], promptGenLeps[event.genLepZ_l1_index]['phi'], 0 )
        genLepZ_l2 = ROOT.TLorentzVector()
        genLepZ_l2.SetPtEtaPhiM(promptGenLeps[event.genLepZ_l2_index]['pt'], promptGenLeps[event.genLepZ_l2_index]['eta'], promptGenLeps[event.genLepZ_l2_index]['phi'], 0 )
        genLepZ = genLepZ_l1 + genLepZ_l2
        event.genLepZ_pt   = genLepZ.Pt()
        event.genLepZ_eta  = genLepZ.Eta()
        event.genLepZ_phi  = genLepZ.Phi()
        event.genLepZ_lldPhi = deltaPhi(promptGenLeps[event.genLepZ_l1_index]['phi'], promptGenLeps[event.genLepZ_l2_index]['phi'])
        event.genLepZ_lldR   = deltaR(promptGenLeps[event.genLepZ_l1_index], promptGenLeps[event.genLepZ_l2_index])
        genLepMinus_index = event.genLepZ_l1_index if promptGenLeps[event.genLepZ_l1_index]['pdgId'] > 0 else event.genLepZ_l2_index
        event.genLepZ_cosThetaStar = cosThetaStar(event.genLepZ_mass, event.genLepZ_pt, event.genLepZ_eta, event.genLepZ_phi, promptGenLeps[genLepMinus_index]['pt'], promptGenLeps[genLepMinus_index]['eta'], promptGenLeps[genLepMinus_index]['phi'] )

    # reco-bjet/nonZ lepton association
    if event.genLepNonZ_l1_index>=0 and genBj0 and genBj1:
        if vecSumPt( genBj0, promptGenLeps[event.genLepNonZ_l1_index], genMet ) > vecSumPt( genBj1, promptGenLeps[event.genLepNonZ_l1_index], genMet ):
            event.genBjNonZlep_index, event.genBjNonZhad_index = genBj0['index'], genBj1['index']
        else:
            event.genBjNonZlep_index, event.genBjNonZhad_index = genBj1['index'], genBj0['index']

    #for jet in genJets:
    #    print jet['isMuon'], jet['isElectron'], jet['isPhoton'], min([999]+[deltaR2(jet, l) for l in promptGenLeps if l['pt']>10]), jet

    # jet/lepton disambiguation -> remove jets, because gen-jets cluster all leptons
    #if args.logLevel == 'DEBUG':
    #    for jet in filter( lambda j: not (min([999]+[deltaR2(j, l) for l in promptGenLeps if l['pt']>10]) > 0.3**2 ), genJets ):
    #        logger.debug( "Filtered gen %f jet %r lep %r", sqrt((min([999]+[deltaR2(jet, l) for l in promptGenLeps if l['pt']>10]))), jet, [ (l['eta'], jet['pt']/l['pt']) for l in promptGenLeps] )
    #        assert False, ""


    fill_vector_collection( event, "genPhoton", gen_photon_varnames, genPhotons_ ) 
    fill_vector_collection( event, "genLep", lep_all_varnames, promptGenLeps)
    fill_vector_collection( event, "genJet", jet_write_varnames, genJets)

tmp_dir     = ROOT.gDirectory
#post_fix = '_%i'%args.job if args.nJobs > 1 else ''
output_filename =  os.path.join(output_directory, sample.name + '.root')

_logger.   add_fileHandler( output_filename.replace('.root', '.log'), args.logLevel )
_logger_rt.add_fileHandler( output_filename.replace('.root', '_rt.log'), args.logLevel )

if os.path.exists( output_filename ) and args.overwrite =='none' :
    logger.info( "File %s found. Quit.", output_filename )
    sys.exit(0)

output_file = ROOT.TFile( output_filename, 'recreate')
output_file.cd()
maker = TreeMaker(
    sequence  = [ filler ],
    variables = [ TreeVariable.fromString(x) for x in variables ] + extra_variables,
    treeName = "Events"
    )

tmp_dir.cd()

counter = 0
reader.start()
maker.start()

while reader.run( ):
    #if abs(map( lambda p: p.daughter(0).pdgId(), filter( lambda p: p.pdgId()==23 and p.numberOfDaughters()==2, reader.products['gp']))[0])==13: 
    #    maker.run()
    #    break
    maker.run()

         
    counter += 1
    if counter == maxEvents:  break

logger.info( "Done with running over %i events.", reader.nEvents )

output_file.cd()
maker.tree.Write()
output_file.Close()

logger.info( "Written output file %s", output_filename )

