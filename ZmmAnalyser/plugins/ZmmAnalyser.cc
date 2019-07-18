// -*- C++ -*-
//
// Package:    ZmmAnalyser/ZmmAnalyser
// Class:      ZmmAnalyser
//
/**\class ZmmAnalyser ZmmAnalyser.cc ZmmAnalyser/ZmmAnalyser/plugins/ZmmAnalyser.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Silvia Taroni
//         Created:  Mon, 18 Feb 2019 11:20:44 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/UserData.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"

#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include <DataFormats/Math/interface/deltaR.h>

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include<vector>
#include <math.h>

#include<TTree.h>
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


using reco::TrackCollection;

class ZmmAnalyser : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit ZmmAnalyser(const edm::ParameterSet&);
      ~ZmmAnalyser();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
  edm::EDGetTokenT<reco::MuonCollection> muonsToken_;  //used to select what tracks to read from configuration file
  edm::EDGetTokenT<EcalRecHitCollection> rechits_EB_;
  edm::EDGetTokenT<reco::PFMETCollection> tok_PFMET_;
  TFile *tree_file;
  TTree * tree; 

  UInt_t     	runNumber;		///< run number
  UShort_t      lumiBlock;		///< lumi section
  Long64_t    eventNumber;	///< event number
  UInt_t	eventTime;		///< unix time of the event
  UShort_t	      nBX;			///< bunch crossing
  int    	nTruePU;		///< pu

  float mu1Pt, mu1Eta, mu1Phi,mu2Pt, mu2Eta, mu2Phi, invMass, met;

  std::vector<float> v_recovRecHits,  v_recovIeta, v_recovIphi, v_recovDphimu1, v_recovDphimu2, v_recovDetamu1, v_recovDetamu2, v_recovDRmu1, v_recovDRmu2, v_sum8, v_aveNeigh;
  std::vector< std::vector<float> > v_neighRecHit, v_neighIeta, v_neighIphi, v_neighDetamu1, v_neighDphimu1, v_neighDetamu2, v_neighDphimu2;
 std::vector<float>  vXtalEn,vNeigh0, vNeigh1, vNeigh2, vNeigh3, vNeigh5, vNeigh6, vNeigh7, vNeigh8;
  std::vector<int> vIc, vIeta, vIphi, vIsm;
  std::vector<unsigned long int> vRawId;

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
ZmmAnalyser::ZmmAnalyser(const edm::ParameterSet& iConfig)
 :
  muonsToken_(consumes<reco::MuonCollection>(iConfig.getUntrackedParameter<edm::InputTag>("muons"))),
  rechits_EB_(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("inputRecHitsEB"))),
  tok_PFMET_(consumes<reco::PFMETCollection>(iConfig.getParameter<edm::InputTag>("METInput")))
{
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("selected", "selected");
   //now do what ever initialization is needed
  tree->Branch("runNumber",     &runNumber,   "runNumber/i");
  tree->Branch("lumiBlock",     &lumiBlock,   "lumiBlock/s");
  tree->Branch("eventNumber",   &eventNumber, "eventNumber/l");
  tree->Branch("eventTime",     &eventTime,   "eventTime/i");
  tree->Branch("nBX",           &nBX,         "nBX/s");

  tree->Branch("mu1Pt" , &mu1Pt , "mu1Pt/F");
  tree->Branch("mu1Eta", &mu1Eta, "mu1Eta/F");
  tree->Branch("mu1Phi", &mu1Phi, "mu1Phi/F");
  tree->Branch("mu2Pt" , &mu2Pt , "mu2Pt/F");
  tree->Branch("mu2Eta", &mu2Eta, "mu2Eta/F");
  tree->Branch("mu2Phi", &mu2Phi, "mu2Phi/F");

  tree->Branch("invMass", &invMass, "invMass/F");
  tree->Branch("met", &met, "met/F");

  tree->Branch("sum8", &v_sum8);
  tree->Branch("aveNeigh", &v_aveNeigh);

  tree->Branch("recovRecHits",&v_recovRecHits);
  tree->Branch("neighRecHit", &v_neighRecHit);
  tree->Branch("recovIeta",   &v_recovIeta);
  tree->Branch("recovIphi",   &v_recovIphi);
  tree->Branch("neighIeta",   &v_neighIeta);
  tree->Branch("neighIphi",   &v_neighIphi);

  tree->Branch("recovDphimu1", &v_recovDphimu1);
  tree->Branch("recovDetamu1", &v_recovDetamu1);
  tree->Branch("recovDphimu2", &v_recovDphimu2);
  tree->Branch("recovDetamu2", &v_recovDetamu2);
  tree->Branch("recovDRmu1", &v_recovDRmu1);
  tree->Branch("recovDRmu2", &v_recovDRmu2);


  tree->Branch("neighDetamu1",&v_neighDetamu1);
  tree->Branch("neighDphimu1",&v_neighDphimu1);
  tree->Branch("neighDetamu2",&v_neighDetamu2);
  tree->Branch("neighDphimu2",&v_neighDphimu2);

  tree->Branch("vNeigh0"   , &vNeigh0);
  tree->Branch("vNeigh1"   , &vNeigh1);
  tree->Branch("vNeigh2"   , &vNeigh2);
  tree->Branch("vNeigh3"   , &vNeigh3);
  tree->Branch("vNeigh5"   , &vNeigh5);
  tree->Branch("vNeigh6"   , &vNeigh6);
  tree->Branch("vNeigh7"   , &vNeigh7);
  tree->Branch("vNeigh8"   , &vNeigh8);

  tree->Branch("xtalIeta",&vIeta );  
  tree->Branch("xtalIphi",&vIphi );  
  tree->Branch("xtalIsm",&vIsm );  
  tree->Branch("xtalIc",&vIc );  
  tree->Branch("xtalRawId",&vRawId );  
  tree->Branch("xtalEn",&vXtalEn);


}


ZmmAnalyser::~ZmmAnalyser()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
ZmmAnalyser::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  using namespace reco;
  using namespace cms;
 std::vector<float> v_recovRecHits_,  v_recovIeta_, v_recovIphi_, v_recovDphimu1_, v_recovDphimu2_, v_recovDetamu1_, v_recovDetamu2_, v_recovDRmu1_, v_recovDRmu2_, v_sum8_, v_aveNeigh_;
  std::vector< std::vector<float> > v_neighRecHit_, v_neighIeta_, v_neighIphi_, v_neighDetamu1_, v_neighDphimu1_, v_neighDetamu2_, v_neighDphimu2_;
  std::vector<float> vNeigh0_,vNeigh1_,vNeigh2_,vNeigh3_,vNeigh5_,vNeigh6_,vNeigh7_,vNeigh8_;
  std::vector<float> vXtalEn_;
  std::vector<int >  vIeta_, vIphi_, vIsm_, vIc_;
  std::vector<unsigned long int> vRawId_;


  runNumber=iEvent.id().run();
  lumiBlock=(UShort_t)iEvent.id().luminosityBlock();
  eventNumber=(Long64_t) iEvent.id().event();
  edm::Timestamp time = iEvent.eventAuxiliary().time();
  eventTime= (UInt_t) time.unixTime();
  nBX=(UShort_t) iEvent.bunchCrossing();	
  
  //get the muons
  Handle<reco::MuonCollection> muons;
  iEvent.getByToken(muonsToken_, muons);

  // get the geometry
  edm::ESHandle<CaloSubdetectorGeometry> theBarrelGeometry_handle;
  iSetup.get<EcalBarrelGeometryRecord>().get("EcalBarrel",theBarrelGeometry_handle);
  const CaloSubdetectorGeometry *theBarrelGeometry;
  theBarrelGeometry = &(*theBarrelGeometry_handle);
  
  // get the topology
  edm::ESHandle<CaloTopology> theCaloTopo;
  iSetup.get<CaloTopologyRecord>().get(theCaloTopo);
  const CaloTopology *topology = theCaloTopo.product();

  // get the  RecHits
  edm::Handle<EcalRecHitCollection> rechit_EB_col;
  iEvent.getByToken(rechits_EB_,rechit_EB_col);
  //const EcalRecHitCollection& unfilteredRecHitsEB = *(rechit_EB_col.product());
  //EcalRecHitCollection recHitsEB;

  //get the Met
  edm::Handle<reco::PFMETCollection> pfmt;
  iEvent.getByToken(tok_PFMET_,pfmt);
  const reco::PFMETCollection pfmet = *(pfmt.product());

  for(reco::MuonCollection::const_iterator itMuon = muons->begin();
      itMuon != muons->end();
      ++itMuon) {
    if (itMuon->isGlobalMuon()==false) continue;
    if (itMuon->pt()<2.5) continue;
    if (itMuon->eta()> 2.4)  continue; 
    
    if (itMuon->passed(reco::Muon::CutBasedIdLoose)==false) continue;
    if (itMuon->passed(reco::Muon::PFIsoVeryLoose) == false) continue;
    
    for(reco::MuonCollection::const_iterator jtMuon = itMuon;
	jtMuon != muons->end();
	++jtMuon) {
      
      if (jtMuon->isGlobalMuon()==false) continue;
      if (jtMuon->pt()<2.5) continue;
      if (jtMuon->eta()> 2.4)  continue; 
      
      if (jtMuon->passed(reco::Muon::CutBasedIdLoose)==false) continue;
      if (jtMuon->passed(reco::Muon::PFIsoVeryLoose) == false) continue;
      
      reco::Candidate::LorentzVector zp4 = itMuon->p4() + jtMuon->p4();
      
      invMass= zp4.M();
      if (abs(91.19 - invMass)>30.) continue;
      std::cout << "Z mass " << invMass<< std::endl;

      // now we can fill the tree
      
      mu1Pt =itMuon->pt();
      mu2Pt =jtMuon->pt();
      mu1Eta=itMuon->eta();
      mu2Eta=jtMuon->eta();
      mu1Phi=itMuon->phi();
      mu2Phi=jtMuon->phi();
      
      met=pfmet.begin()->et();; //add the met here
      //clean the vectors
      vIphi_.clear();
      vIsm_.clear();
      vIc_.clear();
      vRawId_.clear();
      vXtalEn_.clear();

      v_aveNeigh_.clear(); 
      v_sum8_.clear();
      v_recovRecHits_.clear();
      v_recovIeta_.clear();
      v_recovIphi_.clear();
      v_recovDphimu1_.clear();
      v_recovDetamu1_.clear();
      v_recovDphimu2_.clear();
      v_recovDetamu2_.clear();
      v_recovDRmu1_.clear();      
      v_recovDRmu2_.clear();      

      v_neighDetamu1_.clear();
      v_neighDphimu1_.clear();
      v_neighDetamu2_.clear();
      v_neighDphimu2_.clear();

      v_neighRecHit_.clear();
      v_neighIeta_.clear();
      v_neighIphi_.clear();
      
      vNeigh0_.clear();
      vNeigh1_.clear();
      vNeigh2_.clear();
      vNeigh3_.clear();
      vNeigh5_.clear();
      vNeigh6_.clear();
      vNeigh7_.clear();
      vNeigh8_.clear();


      std::vector<DetId> v_id;
      for (edm::SortedCollection<EcalRecHit>::const_iterator hit=rechit_EB_col->begin(); hit!=rechit_EB_col->end(); hit++){
	if (hit->checkFlag(EcalRecHit::kNeighboursRecovered)==false) continue; 
	EBDetId idCurrent= hit->id() ;
	float recIeta=idCurrent.ieta();
	float recIphi=idCurrent.iphi();
	v_recovRecHits_.push_back(hit->energy());
	v_recovIeta_.push_back(recIeta);
	v_recovIphi_.push_back(recIphi);

	float eta =  theBarrelGeometry->getGeometry(idCurrent)->getPosition().eta();
	float phi =  theBarrelGeometry->getGeometry(idCurrent)->getPosition().phi();
	
	float dphimu1= itMuon->phi()-phi;
	//std::cout << __LINE__ << " " << M_PI << " deltaPhimu1 " << dphimu1<<  std::endl;
	if (dphimu1 > M_PI)  dphimu1-=2*M_PI;
	if (dphimu1 < -M_PI) dphimu1+=2*M_PI;
	float dphimu2= jtMuon->phi()-phi;
	if (dphimu2 > M_PI ) dphimu2-=2*M_PI;
	if (dphimu2 <-M_PI ) dphimu2+=2*M_PI;
	//std::cout << __LINE__ << " deltaPhimu1 after changes " <<  dphimu1<<  std::endl;
	v_recovDphimu1_.push_back(dphimu1);
	v_recovDphimu2_.push_back(dphimu2);
	v_recovDetamu1_.push_back(itMuon->eta()-eta);
	v_recovDetamu2_.push_back(jtMuon->eta()-eta);
	
	float dp1= std::abs(itMuon->phi()-phi);
	if (dp1>float(M_PI)) dp1-=float(2*M_PI);
	float dp2= std::abs(jtMuon->phi()-phi);
	if (dp2>float(M_PI)) dp2-=float(2*M_PI);


	float dRmu1=std::sqrt((itMuon->eta()-eta)*(itMuon->eta()-eta)+dp1*dp1);
	float dRmu2=std::sqrt((jtMuon->eta()-eta)*(jtMuon->eta()-eta)+dp2*dp2);
	v_recovDRmu1_.push_back(dRmu1);
	v_recovDRmu2_.push_back(dRmu2);
	std::cout << __LINE__ << " dR1 "<< dRmu1 <<" dR2 "<< dRmu2 << " "<<  (itMuon->eta()-eta)*(itMuon->eta()-eta)<< " "<< dp1*dp1<< " "  << (jtMuon->eta()-eta)*(jtMuon->eta()-eta) << " "<< dp2*dp2<< std::endl;
	
	
	// loop on the neighbours
	float neigh0=0.;
	float neigh1=0.;
	float neigh2=0.;
	float neigh3=0.;
	float neigh5=0.;
	float neigh6=0.;
	float neigh7=0.;
	float neigh8=0.;

	//take the 3x3 matrix around the recovered hit
	v_id = EcalClusterTools::matrixDetId( topology, idCurrent, 1 );
	std::vector<float> vEn_, vIeta_, vIphi_, vDphimu1_, vDphimu2_, vDetamu1_, vDetamu2_;
	float sum8=0; //initialize the sum8

	for ( size_t i = 0; i < v_id.size(); ++i ) {
	  edm::SortedCollection<EcalRecHit>::const_iterator hitNeigh = rechit_EB_col->find( v_id[i] );
	  if (hitNeigh ==  rechit_EB_col->end()) continue;
	  if (hitNeigh == hit) continue;
	  
	  if (hitNeigh->checkFlag(EcalRecHit::kL1SpikeFlag)==true) std::cout<< "SPIKE!" << std::endl;
	  if (hitNeigh->checkFlag(EcalRecHit::kSaturated)==true)continue;
	  sum8+=hitNeigh->energy();
	  std::cout << "neighbour energy :"<<  i << " " <<  hitNeigh->energy() << ", sum8 " << sum8 << ", ave "<< sum8/8. << std::endl;
	  vEn_.push_back(hitNeigh->energy());
	  if (i==0)  neigh0=hitNeigh->energy();
	  if (i==1)  neigh1=hitNeigh->energy();
	  if (i==2)  neigh2=hitNeigh->energy();
	  if (i==3)  neigh3=hitNeigh->energy();
	  if (i==5)  neigh5=hitNeigh->energy();
	  if (i==6)  neigh6=hitNeigh->energy();
	  if (i==7)  neigh7=hitNeigh->energy();
	  if (i==8)  neigh8=hitNeigh->energy();

	  EBDetId neighId=v_id[i];
	  vIeta_.push_back(neighId.ieta()); 
	  vIphi_.push_back(neighId.iphi()); 

	  float etaNeig =  theBarrelGeometry->getGeometry(neighId)->getPosition().eta();
	  float phiNeig =  theBarrelGeometry->getGeometry(neighId)->getPosition().phi();
	  
	  dphimu1= itMuon->phi()-phiNeig;
	  if (dphimu1 > M_PI)  dphimu1-=2*M_PI;
	  if (dphimu1 < -M_PI) dphimu1+=2*M_PI;
	  dphimu2= jtMuon->phi()-phiNeig;
	  if (dphimu2 > M_PI ) dphimu2-=2*M_PI;
	  if (dphimu2 <-M_PI ) dphimu2+=2*M_PI;

	  vDphimu1_.push_back(dphimu1);
	  vDphimu2_.push_back(dphimu2);
	  vDetamu1_.push_back(itMuon->eta()-etaNeig);
	  vDetamu2_.push_back(jtMuon->eta()-etaNeig);
  

	}//loop on the matrix ends
	if (vIeta_.size()!=8) std::cout << "PROBLEM: less than 8 neighbours "<< vIeta_.size() << std::endl;
	vNeigh0_.push_back(neigh0);
	vNeigh1_.push_back(neigh1);
	vNeigh2_.push_back(neigh2);
	vNeigh3_.push_back(neigh3);
	vNeigh5_.push_back(neigh5);
	vNeigh6_.push_back(neigh6);
	vNeigh7_.push_back(neigh7);
	vNeigh8_.push_back(neigh8);

	vIsm_.push_back(idCurrent.ism());
	vIc_.push_back(idCurrent.ic());
	vRawId_.push_back(idCurrent.rawId());
	vXtalEn_.push_back(hit->energy());


	v_sum8_.push_back(sum8);
	v_aveNeigh_.push_back(sum8/float(vIeta_.size()));
	v_neighRecHit_.push_back(vEn_);
	v_neighIeta_.push_back(vIeta_);
	v_neighIphi_.push_back(vIphi_); 

	v_neighDetamu1_.push_back(vDphimu1_);
	v_neighDphimu1_.push_back(vDphimu2_); 
	v_neighDetamu2_.push_back(vDetamu1_);
	v_neighDphimu2_.push_back(vDetamu2_); 

	

      }
      

      v_aveNeigh=v_aveNeigh_;
      v_sum8=v_sum8_;
      v_recovRecHits=v_recovRecHits_;
      v_recovIeta=v_recovIeta_;
      v_recovIphi=v_recovIphi_;
      v_recovDphimu1=v_recovDphimu1_;
      v_recovDetamu1=v_recovDetamu1_;
      v_recovDphimu2=v_recovDphimu2_;
      v_recovDetamu2=v_recovDetamu2_;
      v_recovDRmu1=v_recovDRmu1_;
      v_recovDRmu2=v_recovDRmu2_;

      v_neighDetamu1=v_neighDetamu1_;
      v_neighDphimu1=v_neighDphimu1_;
      v_neighDetamu2=v_neighDetamu2_;
      v_neighDphimu2=v_neighDphimu2_;

      v_neighRecHit=v_neighRecHit;
      v_neighIeta=v_neighIeta;
      v_neighIphi=v_neighIphi;

      vNeigh0=vNeigh0_;
      vNeigh1=vNeigh1_;
      vNeigh2=vNeigh2_;
      vNeigh3=vNeigh3_;
      vNeigh5=vNeigh5_;
      vNeigh6=vNeigh6_;
      vNeigh7=vNeigh7_;
      vNeigh8=vNeigh8_;
	
      vIsm=vIsm_;
      vIc=vIc_;
      vRawId=vRawId_;
      vXtalEn=vXtalEn_;


      tree->Fill();
      
      
    }
  }
  
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
void
ZmmAnalyser::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
ZmmAnalyser::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ZmmAnalyser::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ZmmAnalyser);
