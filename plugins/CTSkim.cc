// -*- C++ -*-
//
// Package:    Test/CTSkim
// Class:      CTSkim
// 
/**\class CTSkim CTSkim.cc Test/CTSkim/plugins/CTSkim.cc

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

#include <string>
#include <iostream>
#include <vector>

#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TLorentzVector.h"

class CTSkim : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit CTSkim(const edm::ParameterSet&);
  ~CTSkim();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  TTree *tree;
  std::vector<reco::Candidate::LorentzVector> trp4;
  std::vector<double> trpt, treta, trphi, dz, d0, dzerr, d0err, vx, vy, vz, chi2n, ptErr;
  std::vector<double> vtxx, vtxy, vtxz, vtxxErr, vtxyErr, vtxzErr, vtxchi2;
  std::vector<double> vtxxBS, vtxyBS, vtxzBS, vtxxErrBS, vtxyErrBS, vtxzErrBS, vtxchi2BS;
  std::vector<int> trsize, highPurity, algo, nValidHits, nLostHits, charge, trigger;
  std::vector<int> vtxisValid, vtxisFake, vtxndof, vtxnTracks;
  std::vector<int> vtxisValidBS, vtxisFakeBS, vtxndofBS, vtxnTracksBS;
  int irun, ilumi, ievt, vtx, vtxBS;
  edm::EDGetTokenT<edm::TriggerResults> trigbit;
  edm::EDGetTokenT<std::vector<reco::Track> > genTrk;
  edm::EDGetTokenT<std::vector<reco::Vertex> > hVtcess;
  edm::EDGetTokenT<std::vector<reco::Vertex> > hVtxx;
      // ----------member data ---------------------------
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
CTSkim::CTSkim(const edm::ParameterSet& iConfig){
  //usesResource("TFileService");
   edm::Service<TFileService> fs;
   tree = fs->make<TTree>("tree","mytree");
   tree -> Branch("Run", &irun); 
   tree -> Branch("Lumi", &ilumi);
   tree -> Branch("Evt", &ievt);

   tree -> Branch("trP4", &trp4);
   tree -> Branch("trPt", &trpt);
   tree -> Branch("trSize", &trsize);
   tree -> Branch("trEta", &treta);
   tree -> Branch("trPhi", &trphi);
   tree -> Branch("dz", &dz);
   tree -> Branch("d0", &d0);
   tree -> Branch("dzerr", &dzerr);
   tree -> Branch("d0err", &d0err);
   tree -> Branch("vx", &vx);
   tree -> Branch("vy", &vy);
   tree -> Branch("vz", &vz);
   tree -> Branch("chi2n", &chi2n);
   tree -> Branch("ptErr", &ptErr);
   tree -> Branch("highPurity", &highPurity);
   tree -> Branch("algo", &algo);
   tree -> Branch("nValidHits", &nValidHits);
   tree -> Branch("nLostHits", &nLostHits);
   tree -> Branch("charge", &charge);

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

   tree -> Branch("trigger", &trigger);

   trigbit = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults","","HLT"));
   genTrk  = consumes<std::vector<reco::Track> >(edm::InputTag("generalTracks"));
   hVtcess = consumes<std::vector<reco::Vertex> >(edm::InputTag("offlinePrimaryVertices"));
   hVtxx = consumes<std::vector<reco::Vertex> >(edm::InputTag("offlinePrimaryVerticesWithBS"));
}


CTSkim::~CTSkim(){}


//
// member functions
//

// ------------ method called for each event  ------------
void CTSkim::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  irun = iEvent.id().run();
  ilumi = iEvent.luminosityBlock();
  ievt = iEvent.id().event();
  
  using namespace edm;
  using namespace std;
  using namespace reco;

  trp4.clear();
  trpt.clear();
  trsize.clear();
  treta.clear();
  trphi.clear();
  dz.clear();
  d0.clear();
  dzerr.clear();
  d0err.clear();
  vx.clear();
  vy.clear();
  vz.clear();
  chi2n.clear();
  ptErr.clear();
  highPurity.clear();
  algo.clear();
  nValidHits.clear();
  nLostHits.clear();
  charge.clear();

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

  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  //InputTag trigbit("TriggerResults","","HLT");
  
  //iEvent.getByToken(trigggg, triggerBits);
  if(!iEvent.getByToken(trigbit, triggerBits)){
    cout<<"Some problem is there\n\n";
    return;
  }
  iEvent.getByLabel("selectedPatTrigger", triggerObjects);
  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  int trgflag=0;
  for(unsigned int itt = 0; itt<names.size(); itt++){
    if((names.triggerName(itt).find("HLT_ZeroBias_") != string::npos))trgflag = (triggerBits->accept(itt)? 1 : 0);
  }
  trigger.push_back(trgflag);

  edm::Handle<std::vector<reco::Vertex> > hVtces;
  iEvent.getByToken(hVtcess, hVtces);
  //vtx = hVtces->size();
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
  if(bestVtx < 0) return;

  edm::Handle<vector<reco::Track> > trs;
  iEvent.getByToken(genTrk, trs);
  trsize.push_back(trs->size());
  for (reco::TrackCollection::const_iterator itTrack = trs->begin(); itTrack != trs->end(); ++itTrack) { 
    if(itTrack->pt()<=0.5 || fabs(itTrack->eta())>=3.0 ) continue;
    double px = itTrack->px();
    double py = itTrack->py();
    double pz = itTrack->pz();
    double E  = px*px + py*py + pz*pz;
    //TLorentzVector pp; 
    reco::Candidate::LorentzVector pp;
    pp.SetPxPyPzE(px,py,pz,E);
    //cout<<"\npx = "<<pp.Px();
    trp4.push_back(pp);
    trpt.push_back(itTrack->pt());
    treta.push_back(itTrack->eta());
    trphi.push_back(itTrack->phi());
    dz.push_back(itTrack->dz());
    d0.push_back(itTrack->d0());
    dzerr.push_back(itTrack->dzError());
    d0err.push_back(itTrack->d0Error());
    vx.push_back(itTrack->vx());
    vy.push_back(itTrack->vy());
    vz.push_back(itTrack->vz());
    int highpurity = 1;
    if(!itTrack->quality(reco::TrackBase::highPurity)) highpurity = 0;
    chi2n.push_back(itTrack->normalizedChi2() );
    ptErr.push_back(itTrack->ptError() );
    highPurity.push_back(highpurity);
    algo.push_back(itTrack->algo() );
    nValidHits.push_back(itTrack->numberOfValidHits() );
    nLostHits.push_back(itTrack->numberOfLostHits() );
    charge.push_back(itTrack->charge());

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
void CTSkim::beginJob(){}

// ------------ method called once each job just after ending the event loop  ------------
void 
CTSkim::endJob(){}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
CTSkim::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(CTSkim);
