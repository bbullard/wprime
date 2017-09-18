#include <iostream>
#include <fstream>
#include <TTree.h>
#include <TCut.h>
using namespace std;

Double_t deltaPhi(Double_t a, Double_t b)
{
	Double_t dPhi = a-b;
	if(dPhi > TMath::Pi()) dPhi -= TMath::Pi()*2.0;
	if(dPhi < -1.0*TMath::Pi()) dPhi += TMath::Pi()*2.0;
	return dPhi;
}

void makeEventTreeQCD(const char* origTreeFile, string outFileName, string source)
{	
	string path = "/export/home/bbullard/thesis/Files/ntuples/selection/";
	path.append(outFileName);
	TFile outF(path.c_str(),"recreate");
	TDirectory *dir = gDirectory;
    TTree *T = new TTree("tree", "New Tree for analysis");
	
    Double_t mT;
	Double_t wmet_et, muon_pt, muon_phi, muon_eta, wmet_phi, muon_d0sig;
	Double_t jet_pt, jet_phi, jet_eta, dR_jetmuon;
	Double_t		KfactorWeight = 0;
	Float_t			mcCrossSection = 0;
	Float_t			mcFilterEfficiency = 0;
	Float_t			mceventweight = 0;
	Double_t 		puweight = 0;
	Float_t totalnAcceptedEvents, totalsumOfEventWeights;
	Bool_t 			isTight;
	
	T->Branch("mT", &mT,"mT/D");
	T->Branch("wmet_et", &wmet_et,"wmet_et/D");
	T->Branch("wmet_phi", &wmet_phi,"wmet_phi/D");
	T->Branch("muon_pt", &muon_pt,"muon_pt/D");
	T->Branch("muon_phi", &muon_phi,"muon_phi/D");
	T->Branch("muon_eta", &muon_eta,"muon_eta/D");
	T->Branch("muon_d0sig", &muon_d0sig,"muon_d0sig/D");
	T->Branch("jet_pt", &jet_pt,"jet_pt/D");
	T->Branch("jet_phi", &jet_phi,"jet_phi/D");
	T->Branch("jet_eta", &jet_eta,"jet_eta/D");
	T->Branch("dR_jetmuon", &dR_jetmuon,"dR_jetmuon/D");
	T->Branch("KfactorWeight", &KfactorWeight,"KfactorWeight/D");
	T->Branch("mcCrossSection", &mcCrossSection,"mcCrossSection/F");
	T->Branch("mcFilterEfficiency", &mcFilterEfficiency,"mcFilterEfficiency/F");
	T->Branch("puweight", &puweight,"puweight/D");
	T->Branch("mceventweight",&mceventweight,"mceventweight/F");
	T->Branch("totalnAcceptedEvents", &totalnAcceptedEvents,"totalnAcceptedEvents/F");
	T->Branch("totalsumOfEventWeights",&totalsumOfEventWeights, "totalsumOfEventWeights/F");
	T->Branch("isTight",&isTight,"isTight/B");
	
	
	ifstream inF(origTreeFile);
	TFile *f, *cf;
	TLeaf *nAcc, *sumW;
	string Rpath, treeFile, treeCountFile;
	TTree *eventTree, *eventTreeCount;
	Bool_t isEventClean, HLT_mu50;
	Int_t 			nmuon_nselec;
	vector<double>  *vmuon_d0sig = 0;
	vector<double>  *vmuon_z0_vrtPVx = 0;
	vector<double>  *vmuon_sintheta = 0;
	vector<float>  	*vmuon_eta = 0;
	vector<float>  	*vmuon_phi = 0;
	vector<float>  	*vmuon_pt = 0;
	Int_t 			njet_nselec;
	vector<float>  	*vjet_eta = 0;
	vector<float>  	*vjet_phi = 0;
	vector<float>  	*vjet_pt = 0;
	vector<bool>  	*vjet_isClean = 0;
	vector<bool>  	*vmuon_isoLooseTrack = 0;
	vector<bool>  	*vmuon_isHighPt = 0;
	vector<bool>  	*vmuon_isBadMuon = 0;
	vector<bool>  	*vmuon_isMedium = 0;
	Int_t 			nelec_nselec;
	vector<float>  	*velec_pt = 0;
	vector<float>  	*velec_eta = 0;
	vector<float>  	*velec_phi = 0;
	vector<bool>  	*velec_isCR = 0;
	Int_t 			goodMuonFlag = -1;
	Int_t			goodMuonFlag2 = -1;
	Int_t			goodJetFlag = -1;
	
	Int_t GRL = 0, clean = 0, mu50 = 0, muon = 0, MV = 0, EV = 0, MET = 0;
	if(source.compare("dataQCD")!=0)
	{
		cout<<"Analyzing MC"<<endl;
		//Loop through samples	
		if(inF.is_open()){
			while(!inF.eof()){
				getline(inF,Rpath);
				if (Rpath.compare("") == 0) continue;
				treeFile = Rpath+"File.root";
				f = TFile::Open(treeFile.c_str(),"read");
				eventTree = (TTree*)f->Get("tree");
				eventTree->SetBranchAddress("isEventClean",&isEventClean);
				eventTree->SetBranchAddress("HLT_mu50",&HLT_mu50);
				eventTree->SetBranchAddress("wmet_et",&wmet_et);
				eventTree->SetBranchAddress("wmet_phi",&wmet_phi);
				eventTree->SetBranchAddress("nmuon_nselec",&nmuon_nselec);
				eventTree->SetBranchAddress("vmuon_d0sig",&vmuon_d0sig);
				eventTree->SetBranchAddress("vmuon_z0_vrtPVx",&vmuon_z0_vrtPVx);
				eventTree->SetBranchAddress("vmuon_sintheta",&vmuon_sintheta);
				eventTree->SetBranchAddress("vmuon_pt",&vmuon_pt);
				eventTree->SetBranchAddress("vmuon_eta",&vmuon_eta);
				eventTree->SetBranchAddress("vmuon_phi",&vmuon_phi);
				eventTree->SetBranchAddress("njet_nselec",&njet_nselec);
				eventTree->SetBranchAddress("vjet_pt",&vjet_pt);
				eventTree->SetBranchAddress("vjet_eta",&vjet_eta);
				eventTree->SetBranchAddress("vjet_phi",&vjet_phi);
				eventTree->SetBranchAddress("vjet_isClean",&vjet_isClean);
				eventTree->SetBranchAddress("vmuon_isoLooseTrack",&vmuon_isoLooseTrack);
				eventTree->SetBranchAddress("vmuon_isHighPt",&vmuon_isHighPt);
				eventTree->SetBranchAddress("vmuon_isBadMuon",&vmuon_isBadMuon);
				eventTree->SetBranchAddress("vmuon_isMedium",&vmuon_isMedium);
				eventTree->SetBranchAddress("nelec_nselec",&nelec_nselec);
				eventTree->SetBranchAddress("velec_pt",&velec_pt);
				eventTree->SetBranchAddress("velec_eta",&velec_eta);
				eventTree->SetBranchAddress("velec_phi",&velec_phi);
				eventTree->SetBranchAddress("velec_isCR",&velec_isCR);
				eventTree->SetBranchAddress("KfactorWeight",&KfactorWeight);
				eventTree->SetBranchAddress("mcCrossSection",&mcCrossSection);
				eventTree->SetBranchAddress("mcFilterEfficiency",&mcFilterEfficiency);
				eventTree->SetBranchAddress("puweight",&puweight);
				eventTree->SetBranchAddress("mceventweight",&mceventweight);

				treeCountFile = Rpath+"CountFile.root";
				cf = TFile::Open(treeCountFile.c_str(),"read");
				eventTreeCount = (TTree*)cf->Get("mcCount");			
				eventTreeCount->SetBranchAddress("totalnAcceptedEvents",&totalnAcceptedEvents);
				eventTreeCount->SetBranchAddress("totalsumOfEventWeights",&totalsumOfEventWeights);
				nAcc = eventTreeCount->GetLeaf("totalnAcceptedEvents");
				sumW = eventTreeCount->GetLeaf("totalsumOfEventWeights");
				nAcc->GetBranch()->GetEntry(0);
				nAcc->GetValue();
				sumW->GetBranch()->GetEntry(0);
				sumW->GetValue();
		
				Long64_t nentries = eventTree->GetEntries();
				for (Long64_t i=0;i<nentries;i++) 
				{
   	  				eventTree->GetEntry(i);
					goodMuonFlag = -1;
					goodMuonFlag2 = -1;
					goodJetFlag = -1;
					vetoElec = 0;
					vetoMuon = 0;
					
					//Find good muon
					for(Int_t j=0;j<nmuon_nselec;j++)
					{
						if(!vmuon_isHighPt->at(j) || vmuon_isBadMuon->at(j)) continue;
						if(vmuon_pt->at(j) < 55000) continue;
						//Is combined?
						//Take higher pt muon
						if(goodMuonFlag >= 0 && vmuon_pt->at(goodMuonFlag) < vmuon_pt->at(j)) {goodMuonFlag2 = goodMuonFlag; goodMuonFlag = j;}
						if(goodMuonFlag == -1) {goodMuonFlag2 = goodMuonFlag; goodMuonFlag = j;}
					}
					if(goodMuonFlag == -1) continue;
					
					isTight = vmuon_isoLooseTrack->at(goodMuonFlag);
										
					//Fake Muon Control Region
					if(wmet_et > 55000) continue;
					
					for(Int_t j=0;j<njet_nselec;j++)
					{
						if(!vjet_isClean->at(j)) continue;
						dR_jetmuon = sqrt(pow(deltaPhi(vmuon_phi->at(goodMuonFlag),vjet_phi->at(j)),2)+pow(vmuon_eta->at(goodMuonFlag)-vjet_eta->at(j),2));
						if(vjet_pt->at(j)>40000 && dR_jetmuon > 0.2) goodJetFlag = j;
					}
					if(goodJetFlag==-1) continue;
					
					if(fabs(deltaPhi(vmuon_phi->at(goodMuonFlag),wmet_phi)) > 0.5) continue;
					
					muon_d0sig = vmuon_d0sig->at(goodMuonFlag);
					
					//Using mT here as the invariant mass of highest energy muons to avoid creating another variable
					if(goodMuonFlag>=0 && goodMuonFlag2 >=0)
					{
						mT = sqrt(2*vmuon_pt->at(goodMuonFlag)*vmuon_pt->at(goodMuonFlag2)*(cosh(vmuon_eta->at(goodMuonFlag)-vmuon_eta->at(goodMuonFlag2))-cos(deltaPhi(vmuon_phi->at(goodMuonFlag),vmuon_phi->at(goodMuonFlag2)))));
						if(mT>80&&mT<100) continue;
					}
					
					mT = sqrt(2*wmet_et*vmuon_pt->at(goodMuonFlag)*(1-cos(deltaPhi(wmet_phi,vmuon_phi->at(goodMuonFlag)))))/1000;
					muon_pt = vmuon_pt->at(goodMuonFlag)/1000;
					wmet_et = wmet_et/1000;
					muon_eta = vmuon_eta->at(goodMuonFlag);
					muon_phi = vmuon_phi->at(goodMuonFlag);
					jet_pt = vjet_pt->at(goodJetFlag)/1000;
					jet_eta = vjet_eta->at(goodJetFlag);
					jet_phi = vjet_phi->at(goodJetFlag);
					T->Fill();
  				}
				f->Close();
				cf->Close();
			}	
			inF.close();
		}
		dir->cd();
    	T->Write();
    	outF.Close();
	}
	else
	{
		cout<<"Analyzing data"<<endl;
		//Loop through samples	
		if(inF.is_open()){
			while(!inF.eof()){
				getline(inF,Rpath);
				if (Rpath.compare("") == 0) continue;
				f = TFile::Open(Rpath.c_str(),"read");
				eventTree = (TTree*)f->Get("tree");
				eventTree->SetBranchAddress("isEventClean",&isEventClean);
				eventTree->SetBranchAddress("HLT_mu50",&HLT_mu50);
				eventTree->SetBranchAddress("wmet_et",&wmet_et);
				eventTree->SetBranchAddress("wmet_phi",&wmet_phi);
				eventTree->SetBranchAddress("nmuon_nselec",&nmuon_nselec);
				eventTree->SetBranchAddress("vmuon_d0sig",&vmuon_d0sig);
				eventTree->SetBranchAddress("vmuon_z0_vrtPVx",&vmuon_z0_vrtPVx);
				eventTree->SetBranchAddress("vmuon_sintheta",&vmuon_sintheta);
				eventTree->SetBranchAddress("vmuon_pt",&vmuon_pt);
				eventTree->SetBranchAddress("vmuon_eta",&vmuon_eta);
				eventTree->SetBranchAddress("vmuon_phi",&vmuon_phi);
				eventTree->SetBranchAddress("njet_nselec",&njet_nselec);
				eventTree->SetBranchAddress("vjet_isClean",&vjet_isClean);
				eventTree->SetBranchAddress("vjet_pt",&vjet_pt);
				eventTree->SetBranchAddress("vjet_eta",&vjet_eta);
				eventTree->SetBranchAddress("vjet_phi",&vjet_phi);
				eventTree->SetBranchAddress("vmuon_isoLooseTrack",&vmuon_isoLooseTrack);
				eventTree->SetBranchAddress("vmuon_isHighPt",&vmuon_isHighPt);
				eventTree->SetBranchAddress("vmuon_isBadMuon",&vmuon_isBadMuon);
				eventTree->SetBranchAddress("vmuon_isMedium",&vmuon_isMedium);
				eventTree->SetBranchAddress("nelec_nselec",&nelec_nselec);
				eventTree->SetBranchAddress("velec_pt",&velec_pt);
				eventTree->SetBranchAddress("velec_eta",&velec_eta);
				eventTree->SetBranchAddress("velec_phi",&velec_phi);
				eventTree->SetBranchAddress("velec_isCR",&velec_isCR);
				eventTree->SetBranchAddress("KfactorWeight",&KfactorWeight);
				eventTree->SetBranchAddress("mcCrossSection",&mcCrossSection);
				eventTree->SetBranchAddress("mcFilterEfficiency",&mcFilterEfficiency);
				eventTree->SetBranchAddress("puweight",&puweight);
				eventTree->SetBranchAddress("mceventweight",&mceventweight);
				
				Long64_t nentries = eventTree->GetEntries();
				for (Long64_t i=0;i<nentries;i++) 
				{
   	  				eventTree->GetEntry(i);
					goodMuonFlag = -1;
					goodMuonFlag2 = -1;
					goodJetFlag = -1;
					vetoElec = 0;
					vetoMuon = 0;
					
					//Find good muon
					for(Int_t j=0;j<nmuon_nselec;j++)
					{
						if(!vmuon_isHighPt->at(j) || vmuon_isBadMuon->at(j)) continue;
						if(vmuon_pt->at(j) < 55000) continue;
						//Is combined?
						//Take higher pt muon
						if(goodMuonFlag >= 0 && vmuon_pt->at(goodMuonFlag) < vmuon_pt->at(j)) {goodMuonFlag2 = goodMuonFlag; goodMuonFlag = j;}
						if(goodMuonFlag == -1) {goodMuonFlag2 = goodMuonFlag; goodMuonFlag = j;}
					}
					if(goodMuonFlag == -1) continue;
					
					isTight = vmuon_isoLooseTrack->at(goodMuonFlag);
										
					//Fake Muon Control Region
					if(wmet_et > 55000) continue;
					
					for(Int_t j=0;j<njet_nselec;j++)
					{
						if(!vjet_isClean->at(j)) continue;
						dR_jetmuon = sqrt(pow(deltaPhi(vmuon_phi->at(goodMuonFlag),vjet_phi->at(j)),2)+pow(vmuon_eta->at(goodMuonFlag)-vjet_eta->at(j),2));
						if(vjet_pt->at(j)>40000 && dR_jetmuon > 0.2) goodJetFlag = j;
					}
					if(goodJetFlag==-1) continue;
					
					if(fabs(deltaPhi(vmuon_phi->at(goodMuonFlag),wmet_phi)) > 0.5) continue;
					
					muon_d0sig = vmuon_d0sig->at(goodMuonFlag);
					
					//Using mT here as the invariant mass of highest energy muons to avoid creating another variable
					if(goodMuonFlag>=0 && goodMuonFlag2 >=0)
					{
						mT = sqrt(2*vmuon_pt->at(goodMuonFlag)*vmuon_pt->at(goodMuonFlag2)*(cosh(vmuon_eta->at(goodMuonFlag)-vmuon_eta->at(goodMuonFlag2))-cos(deltaPhi(vmuon_phi->at(goodMuonFlag),vmuon_phi->at(goodMuonFlag2)))));
						if(mT>80&&mT<100) continue;
					}
					
					mT = sqrt(2*wmet_et*vmuon_pt->at(goodMuonFlag)*(1-cos(deltaPhi(wmet_phi,vmuon_phi->at(goodMuonFlag)))))/1000;
					muon_pt = vmuon_pt->at(goodMuonFlag)/1000;
					wmet_et = wmet_et/1000;
					muon_eta = vmuon_eta->at(goodMuonFlag);
					muon_phi = vmuon_phi->at(goodMuonFlag);
					jet_pt = vjet_pt->at(goodJetFlag)/1000;
					jet_eta = vjet_eta->at(goodJetFlag);
					jet_phi = vjet_phi->at(goodJetFlag);
					T->Fill();
  				}
				f->Close();
			}
			inF.close();
		}
		dir->cd();
    	T->Write();
    	outF.Close();
	}
	
	
	ofstream cutflow;
  	cutflow.open("/export/home/bbullard/thesis/Files/ntuples/selection/"+source+"cutflow.txt", std::ofstream::app);
  	cutflow <<"\nIs on GRL:\t\t"<<GRL;
	cutflow <<"\n+Is clean:\t\t"<<clean;
	cutflow <<"\n+passes HLT_mu50:\t"<<mu50;
	cutflow <<"\n+has good muon:\t\t"<<muon;
	cutflow <<"\n+passes muon veto:\t"<<MV;
	cutflow <<"\n+passes electron veto:\t"<<EV;
	cutflow <<"\n+MET > 55 GeV:\t\t"<<MET<<endl;
  	cutflow.close();
}
