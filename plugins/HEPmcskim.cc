// -*- C++ -*-
//
// Package:    Test/HEPmcskim
// Class:      HEPmcskim
// 
/**\class HEPmcskim HEPmcskim.cc HEPmcskim/plugins/HEPmcskim.cc
 Description: [one line class summary]
 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Wei Yang
//         Created:  Sun, 30 Aug 2016 11:46:30 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"
#include "DataFormats/Luminosity/interface/LumiDetails.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include <string>
#include <iostream>
#include <vector>

#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TLorentzVector.h"

class HEPmcskim : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit HEPmcskim(const edm::ParameterSet&);
  ~HEPmcskim();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  TTree *tree;
  std::vector<reco::Candidate::LorentzVector> recoTracksp4;
  std::vector<double> recoTrackspt, recoTrackseta, recoTracksphi, recoTracksdz, recoTracksd0, recoTracksdzErr, recoTracksd0Err, recoTracksvx, recoTracksvy, recoTracksvz, recoTrackschi2n, recoTracksptErr;
  std::vector<double> recoTracksdxyvtxBS, recoTracksdxyBS, recoTracksdzvtxBS, recoTracksdzBS;
  std::vector<reco::Candidate::LorentzVector> genParticlesp4;
  std::vector<double> genParticlespt, genParticleseta, genParticlesphi, genParticlesmass, genParticlesvx, genParticlesvy, genParticlesvz;
  std::vector<double> vtxx, vtxy, vtxz, vtxxErr, vtxyErr, vtxzErr, vtxchi2;
  std::vector<double> vtxxBS, vtxyBS, vtxzBS, vtxxErrBS, vtxyErrBS, vtxzErrBS, vtxchi2BS;
  double BSx, BSy, BSz, BSzerr, BSxwidth, BSywidth; //no BSyerr, BSxerr, use BSxwidth, BSywidth. // BSdxz removed.
  std::vector<int> recoTrackssize, recoTrackshighPurity, recoTracksalgo, recoTracksnValidHits, recoTracksnLostHits, recoTrackscharge, trigger, triggerMB, triggerHM60, triggerHM85, triggerHM110, triggerHM135, triggerHM160;
  std::vector<int> genParticlessize, genParticlespdgid, genParticlesst, genParticlescharge;
  std::vector<int> vtxisValid, vtxisFake, vtxndof, vtxnTracks;
  std::vector<int> vtxisValidBS, vtxisFakeBS, vtxndofBS, vtxnTracksBS;
  int run, lumi, event, vtx, vtxBS, ibestvtxBS;
  edm::EDGetTokenT<edm::TriggerResults> trigbit;
  //edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> trigObj; 
  edm::EDGetTokenT<std::vector<reco::Track> > genTrk;
  edm::EDGetTokenT<reco::GenParticleCollection> genPart;
  edm::EDGetTokenT<std::vector<reco::Vertex> > hVtcess;
  edm::EDGetTokenT<std::vector<reco::Vertex> > hVtxx;
  edm::EDGetTokenT<reco::BeamSpot> BS;
      // ----------member data --------------------------
  int countEvt, countZB, countMB, countHM60, countHM85, countHM110, countHM135, countHM160;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
HEPmcskim::HEPmcskim(const edm::ParameterSet& iConfig){
  //usesResource("TFileService");
   edm::Service<TFileService> fs;
   tree = fs->make<TTree>("tree","tree");
   
   
   //----------these in root tree--------------
   
   tree -> Branch("Run", &run); 
   tree -> Branch("Lumi", &lumi);
   tree -> Branch("Evt", &event);

   tree -> Branch("trP4", &recoTracksp4);
   tree -> Branch("trSize", &recoTrackssize);
   tree -> Branch("trPt", &recoTrackspt);
   tree -> Branch("trEta", &recoTrackseta);
   tree -> Branch("trPhi", &recoTracksphi);
   tree -> Branch("dz", &recoTracksdz);
   tree -> Branch("dzvtxBS", &recoTracksdzvtxBS);
   tree -> Branch("dzBS", &recoTracksdzBS);
   tree -> Branch("d0", &recoTracksd0);
   tree -> Branch("dxyvtxBS", &recoTracksdxyvtxBS);
   tree -> Branch("dxyBS", &recoTracksdxyBS);
   tree -> Branch("dzerr", &recoTracksdzErr);
   tree -> Branch("d0err", &recoTracksd0Err);
   tree -> Branch("vx", &recoTracksvx);
   tree -> Branch("vy", &recoTracksvy);
   tree -> Branch("vz", &recoTracksvz);
   tree -> Branch("chi2n", &recoTrackschi2n);
   tree -> Branch("ptErr", &recoTracksptErr);
   tree -> Branch("highPurity", &recoTrackshighPurity);
   tree -> Branch("algo", &recoTracksalgo);
   tree -> Branch("nValidHits", &recoTracksnValidHits);
   tree -> Branch("nLostHits", &recoTracksnLostHits);
   tree -> Branch("charge", &recoTrackscharge);

   tree -> Branch("partP4", &genParticlesp4);
   tree -> Branch("partSize", &genParticlessize);
   tree -> Branch("partPt", &genParticlespt);
   tree -> Branch("partEta", &genParticleseta);
   tree -> Branch("partPhi", &genParticlesphi);
   tree -> Branch("partvx", &genParticlesvx);
   tree -> Branch("partvy", &genParticlesvy);
   tree -> Branch("partvz", &genParticlesvz);
   tree -> Branch("partpdgId", &genParticlespdgid);
   tree -> Branch("partStatus", &genParticlesst);
   tree -> Branch("partCharge", &genParticlescharge);
   tree -> Branch("partMass", &genParticlesmass);

   tree -> Branch("vtx", &vtx);
   tree -> Branch("vtxx", &vtxx); 
   tree -> Branch("vtxy", &vtxy); 
   tree -> Branch("vtxz", &vtxz); 
   tree -> Branch("vtxxErr", &vtxxErr); 
   tree -> Branch("vtxyErr", &vtxyErr); 
   tree -> Branch("vtxzErr", &vtxzErr); 
   tree -> Branch("vtxchi2", &vtxchi2); 
   tree -> Branch("vtxisValid", &vtxisValid); 
   tree -> Branch("vtxisFake", &vtxisFake); 
   tree -> Branch("vtxndof", &vtxndof); 
   tree -> Branch("vtxnTracks", &vtxnTracks); 

   tree -> Branch("ibestvtxBS", &ibestvtxBS);
   tree -> Branch("vtxBS", &vtxBS);
   tree -> Branch("vtxxBS", &vtxxBS); 
   tree -> Branch("vtxyBS", &vtxyBS); 
   tree -> Branch("vtxzBS", &vtxzBS); 
   tree -> Branch("vtxxErrBS", &vtxxErrBS); 
   tree -> Branch("vtxyErrBS", &vtxyErrBS); 
   tree -> Branch("vtxzErrBS", &vtxzErrBS); 
   tree -> Branch("vtxchi2BS", &vtxchi2BS); 
   tree -> Branch("vtxisValidBS", &vtxisValidBS); 
   tree -> Branch("vtxisFakeBS", &vtxisFakeBS); 
   tree -> Branch("vtxndofBS", &vtxndofBS); 
   tree -> Branch("vtxnTracksBS", &vtxnTracksBS); 

   tree -> Branch("BSx", &BSx);
   tree -> Branch("BSy", &BSy);
   tree -> Branch("BSz", &BSz);
   //tree -> Branch("BSxerr", &BSxerr);
   //tree -> Branch("BSyerr", &BSyerr);
   tree -> Branch("BSzerr", &BSzerr);
   //tree -> Branch("BSdxz", &BSdxz);
   tree -> Branch("BSxwidth", &BSxwidth);   //corresponds to beamspot x error
   tree -> Branch("BSywidth", &BSywidth);   //corresponds to beamspot y error

   tree -> Branch("trigger", &trigger);
   tree -> Branch("triggerMB", &triggerMB);
   tree -> Branch("triggerHM60", &triggerHM60);
   tree -> Branch("triggerHM85", &triggerHM85);
   tree -> Branch("triggerHM110", &triggerHM110);
   tree -> Branch("triggerHM135", &triggerHM135);
   tree -> Branch("triggerHM160", &triggerHM160);

   trigbit = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults","","HLT"));
   //trigObj = consumes<pat::TriggerObjectStandAloneCollection>(edm::InputTag("selectedPatTrigger"));
   genTrk  = consumes<std::vector<reco::Track> >(edm::InputTag("generalTracks"));
   genPart = consumes<reco::GenParticleCollection>(edm::InputTag("genParticles"));
   hVtcess = consumes<std::vector<reco::Vertex> >(edm::InputTag("offlinePrimaryVertices"));
   hVtxx = consumes<std::vector<reco::Vertex> >(edm::InputTag("offlinePrimaryVerticesWithBS"));
   BS = consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"));
}


HEPmcskim::~HEPmcskim(){}


//
// member functions
//

// ------------ method called for each event  ------------
void HEPmcskim::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  run = iEvent.id().run();
  lumi = iEvent.luminosityBlock();
  event = iEvent.id().event();
  
  using namespace edm;
  using namespace std;
  using namespace reco;

  recoTracksp4.clear();
  recoTrackspt.clear();
  recoTrackssize.clear();
  recoTrackseta.clear();
  recoTracksphi.clear();
  recoTracksdz.clear();
  recoTracksdzvtxBS.clear();
  recoTracksdzBS.clear();
  recoTracksd0.clear();
  recoTracksdxyvtxBS.clear();
  recoTracksdxyBS.clear();
  recoTracksdzErr.clear();
  recoTracksd0Err.clear();
  recoTracksvx.clear();
  recoTracksvy.clear();
  recoTracksvz.clear();
  recoTrackschi2n.clear();
  recoTracksptErr.clear();
  recoTrackshighPurity.clear();
  recoTracksalgo.clear();
  recoTracksnValidHits.clear();
  recoTracksnLostHits.clear();
  recoTrackscharge.clear();

  genParticlesp4.clear();
  genParticlessize.clear();
  genParticlespt.clear();
  genParticleseta.clear();
  genParticlesphi.clear();
  genParticlesvx.clear();
  genParticlesvy.clear();
  genParticlesvz.clear();
  genParticlespdgid.clear();
  genParticlesst.clear();
  genParticlescharge.clear();
  genParticlesmass.clear();

  //vtx.clear(); //integer
  //vtxBS.clear(); //integer
  vtxx.clear();
  vtxy.clear();
  vtxz.clear();
  vtxxErr.clear();
  vtxyErr.clear();
  vtxzErr.clear();
  vtxchi2.clear();
  vtxisValid.clear();
  vtxisFake.clear();
  vtxndof.clear();
  vtxnTracks.clear();

  vtxxBS.clear();
  vtxyBS.clear();
  vtxzBS.clear();
  vtxxErrBS.clear();
  vtxyErrBS.clear();
  vtxzErrBS.clear();
  vtxchi2BS.clear();
  vtxisValidBS.clear();
  vtxisFakeBS.clear();
  vtxndofBS.clear();
  vtxnTracksBS.clear();

  trigger.clear();
  triggerMB.clear();
  triggerHM60.clear();
  triggerHM85.clear();
  triggerHM110.clear();
  triggerHM135.clear();
  triggerHM160.clear();

  edm::Handle<edm::TriggerResults> triggerBits;
  //edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  //InputTag trigbit("TriggerResults","","HLT");
  
  //iEvent.getByToken(trigggg, triggerBits);
  if(!iEvent.getByToken(trigbit, triggerBits)){// || !iEvent.getByToken(trigObj, triggerObjects)){
    cout<<"Some problem getting token.\n\n";
    return;
  }

  countEvt++;
  //iEvent.getByLabel("selectedPatTrigger", triggerObjects);
  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  int trgflag=0, trgMBflag=0, trgHM60flag=0, trgHM85flag=0, trgHM110flag=0, trgHM135flag=0, trgHM160flag=0;
  for(unsigned int itt = 0; itt<names.size(); itt++)
  {
    if((names.triggerName(itt).find("HLT_ZeroBias_") != string::npos)) trgflag += (triggerBits->accept(itt) ? 1 : 0);
    if((names.triggerName(itt).find("HLT_L1MinimumBiasHF1OR_part") != string::npos)) trgMBflag += (triggerBits->accept(itt) ? 1 : 0);
    //if((names.triggerName(itt).find("HLT_PixelTracks_Multiplicity60ForEndOfFill_v1") != string::npos)) trgHM60flag += (triggerBits->accept(itt) ? 1 : 0);
    //if((names.triggerName(itt).find("HLT_PixelTracks_Multiplicity85ForEndOfFill_v1") != string::npos)) trgHM85flag += (triggerBits->accept(itt) ? 1 : 0);
    //if((names.triggerName(itt).find("HLT_PixelTracks_Multiplicity110ForEndOfFill_v1") != string::npos)) trgHM110flag += (triggerBits->accept(itt) ? 1 : 0);
    //if((names.triggerName(itt).find("HLT_PixelTracks_Multiplicity135ForEndOfFill_v1") != string::npos)) trgHM135flag += (triggerBits->accept(itt) ? 1 : 0);
    if((names.triggerName(itt).find("HLT_PixelTracks_Multiplicity60_v2") != string::npos)) trgHM60flag += (triggerBits->accept(itt) ? 1 : 0);
    if((names.triggerName(itt).find("HLT_PixelTracks_Multiplicity85_v2") != string::npos)) trgHM85flag += (triggerBits->accept(itt) ? 1 : 0);
    if((names.triggerName(itt).find("HLT_PixelTracks_Multiplicity110_v2") != string::npos)) trgHM110flag += (triggerBits->accept(itt) ? 1 : 0);
    if((names.triggerName(itt).find("HLT_PixelTracks_Multiplicity135_v2") != string::npos)) trgHM135flag += (triggerBits->accept(itt) ? 1 : 0);
    if((names.triggerName(itt).find("HLT_PixelTracks_Multiplicity160_v2") != string::npos)) trgHM160flag += (triggerBits->accept(itt) ? 1 : 0);
  }
  trigger.push_back(trgflag);
  triggerMB.push_back(trgMBflag);
  triggerHM60.push_back(trgHM60flag);
  triggerHM85.push_back(trgHM85flag);
  triggerHM110.push_back(trgHM110flag);
  triggerHM135.push_back(trgHM135flag);
  triggerHM160.push_back(trgHM160flag);
  countZB += trgflag;
  countMB += trgMBflag;
  countHM60 += trgHM60flag;
  countHM85 += trgHM85flag;
  countHM110 += trgHM110flag;
  countHM135 += trgHM135flag;
  countHM160 += trgHM160flag;
//HM trigger names: HLT_PixelTracks_Multiplicity60ForEndOfFill_v1, HLT_PixelTracks_Multiplicity85ForEndOfFill_v1, HLT_PixelTracks_Multiplicity110ForEndOfFill_v1, HLT_PixelTracks_Multiplicity135ForEndOfFill_v1

  edm::Handle<std::vector<reco::Vertex> > hVtces;
  iEvent.getByToken(hVtcess, hVtces);
  vtx = hVtces->size();
  for(reco::VertexCollection::const_iterator ivtc = hVtces->begin(); ivtc != hVtces->end(); ++ivtc){
    vtxx.push_back(ivtc->x());
    vtxy.push_back(ivtc->y());
    vtxz.push_back(ivtc->z());
    vtxxErr.push_back(ivtc->xError());
    vtxyErr.push_back(ivtc->yError());
    vtxzErr.push_back(ivtc->zError());
    vtxisValid.push_back(ivtc->isValid());
    vtxisFake.push_back(ivtc->isFake());
    vtxchi2.push_back(ivtc->chi2());
    vtxndof.push_back(ivtc->ndof());
    vtxnTracks.push_back(ivtc->nTracks());
  }

  edm::Handle<std::vector<reco::Vertex> > hVtx;
  iEvent.getByToken(hVtxx, hVtx);
  vtxBS = hVtx->size();
  for(reco::VertexCollection::const_iterator ivtx = hVtx->begin(); ivtx != hVtx->end(); ++ivtx){
    vtxxBS.push_back(ivtx->x());
    vtxyBS.push_back(ivtx->y());
    vtxzBS.push_back(ivtx->z());
    vtxxErrBS.push_back(ivtx->xError());
    vtxyErrBS.push_back(ivtx->yError());
    vtxzErrBS.push_back(ivtx->zError());
    vtxisValidBS.push_back(ivtx->isValid());
    vtxisFakeBS.push_back(ivtx->isFake());
    vtxchi2BS.push_back(ivtx->chi2());
    vtxndofBS.push_back(ivtx->ndof());
    vtxnTracksBS.push_back(ivtx->nTracks());
  }  

  reco::BeamSpot beamSpot;
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  //iEvent.getByLabel("offlineBeamSpot", beamSpotHandle);
  iEvent.getByToken(BS, beamSpotHandle);

  if ( beamSpotHandle.isValid() )
  {
      beamSpot = *beamSpotHandle;

  } else
  {
      edm::LogInfo("MyAnalyzer")
        << "No beam spot available from EventSetup \n";
  }

  BSx = beamSpot.x0();
  BSy = beamSpot.y0();
  BSz = beamSpot.z0();
  BSzerr = beamSpot.sigmaZ();
  //BSyerr = beamSpot.sigmaY();
  //BSzerr = beamSpot.sigmaX();
  //BSdxz = beamSpot.dxdz();
  BSxwidth = beamSpot.BeamWidthX();  //beamspot error in x
  BSywidth = beamSpot.BeamWidthY();  //beamspot error in y

  float bestSum = 0;
  int bestVtx = -1;
  for (unsigned int ivtx = 0; ivtx< hVtx->size();++ivtx){
    reco::Vertex::trackRef_iterator it = hVtx->at(ivtx).tracks_begin();
    reco::Vertex::trackRef_iterator itE = hVtx->at(ivtx).tracks_end();
    float sum= 0;
    for (;it!=itE;++it ){
      sum+=(*it)->pt();
    }
    if (bestSum < sum){
      bestSum = sum;
      bestVtx = ivtx;
    }
  }
  ibestvtxBS = bestVtx;
  if(bestVtx < 0) return;

  edm::Handle<reco::GenParticleCollection> ps;
  iEvent.getByToken(genPart, ps);
  genParticlessize.push_back(ps->size());
  for(reco::GenParticleCollection::const_iterator itTrack = ps->begin(); itTrack != ps->end(); ++itTrack) {
    //if(itTrack->pt()<=0.5 || fabs(itTrack->eta())>=3.0 ) continue;
    double px = itTrack->px();
    double py = itTrack->py();
    double pz = itTrack->pz();
    double E  = px*px + py*py + pz*pz;
    reco::Candidate::LorentzVector pp;
    pp.SetPxPyPzE(px,py,pz,E);
    genParticlesp4.push_back(pp);
    genParticlespt.push_back(itTrack->pt());
    genParticleseta.push_back(itTrack->eta());
    genParticlesphi.push_back(itTrack->phi());
    genParticlesmass.push_back(itTrack->mass());
    genParticlesvx.push_back(itTrack->vx());
    genParticlesvy.push_back(itTrack->vy());    
    genParticlesvz.push_back(itTrack->vz());
    genParticlespdgid.push_back(itTrack->pdgId());
    genParticlesst.push_back(itTrack->status());
    genParticlescharge.push_back(itTrack->charge());
  }


  edm::Handle<vector<reco::Track> > trs;
  iEvent.getByToken(genTrk, trs);
  recoTrackssize.push_back(trs->size());
  for (reco::TrackCollection::const_iterator itTrack = trs->begin(); itTrack != trs->end(); ++itTrack) { 
    //if(itTrack->pt()<=0.5 || fabs(itTrack->eta())>=3.0 ) continue;
    double px = itTrack->px();
    double py = itTrack->py();
    double pz = itTrack->pz();
    double E  = px*px + py*py + pz*pz;
    //TLorentzVector pp; 
    reco::Candidate::LorentzVector pp;
    pp.SetPxPyPzE(px,py,pz,E);
    //cout<<"\npx = "<<pp.Px();
    recoTracksp4.push_back(pp);
    recoTrackspt.push_back(itTrack->pt());
    recoTrackseta.push_back(itTrack->eta());
    recoTracksphi.push_back(itTrack->phi());
    recoTracksdz.push_back(itTrack->dz());
    recoTracksdzvtxBS.push_back(itTrack->dz(hVtx->at(ibestvtxBS).position()));
    recoTracksdzBS.push_back(itTrack->dz(beamSpot.position()));
    recoTracksd0.push_back(itTrack->d0());
    recoTracksdxyvtxBS.push_back(itTrack->dxy(hVtx->at(ibestvtxBS).position()));
    recoTracksdxyBS.push_back(itTrack->dxy(beamSpot.position()));
    recoTracksdzErr.push_back(itTrack->dzError());
    recoTracksd0Err.push_back(itTrack->d0Error());
    recoTracksvx.push_back(itTrack->vx());
    recoTracksvy.push_back(itTrack->vy());
    recoTracksvz.push_back(itTrack->vz());
    int highpurity = 1;
    if(!itTrack->quality(reco::TrackBase::highPurity)) highpurity = 0;
    recoTrackschi2n.push_back(itTrack->normalizedChi2() );
    recoTracksptErr.push_back(itTrack->ptError() );
    recoTrackshighPurity.push_back(highpurity);
    recoTracksalgo.push_back(itTrack->algo() );
    recoTracksnValidHits.push_back(itTrack->numberOfValidHits() );
    recoTracksnLostHits.push_back(itTrack->numberOfLostHits() );
    recoTrackscharge.push_back(itTrack->charge());

  }
  
  tree->Fill();

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void HEPmcskim::beginJob(){
  countEvt = 0;
  countZB = 0;
  countMB = 0;
  countHM60 = 0;
  countHM85 = 0;
  countHM110 = 0;
  countHM135 = 0;
  countHM160 = 0;
}

// ------------ method called once each job just after ending the event loop  ------------
void HEPmcskim::endJob(){
  std::cout<< countEvt << " " << countZB << " " <<  countMB << "  " << countHM60 << " " << countHM85 << " " <<  countHM110 << " " <<  countHM135 << " " <<  countHM160 << std::endl;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HEPmcskim::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HEPmcskim);
