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

void makeEventTree(const char* origTreeFile, string outFileName, string source)
{	
	string path = "/export/home/bbullard/thesis/Files/ntuples/selection/";
	path.append(outFileName);
	TFile outF(path.c_str(),"recreate");
	TDirectory *dir = gDirectory;
    TTree *T = new TTree("tree", "New Tree for analysis");
	
    Double_t mT, dR;
	Double_t wmet_et, muon_pt, muon_phi, muon_eta, wmet_phi, dPhi_metmuon;
	vector<float>  	*vjet_eta = 0;
	vector<float>  	*vjet_phi = 0;
	vector<float>  	*vjet_pt = 0;
	vector<float>  	*vdR_jetmuon = 0;
	Double_t		KfactorWeight = 0;
	Float_t			mcCrossSection = 0;
	Float_t			mcFilterEfficiency = 0;
	Float_t			mceventweight = 0;
	Double_t 		puweight = 0;
	Float_t totalnAcceptedEvents, totalsumOfEventWeights;
	Bool_t 			isTight;
	Double_t 		muon_HighPtSF;
	Double_t 		SF_trigmu50_highptID;
	Bool_t 			muon_HLTtrigmatch;
	
	T->Branch("mT", &mT,"mT/D");
	T->Branch("wmet_et", &wmet_et,"wmet_et/D");
	T->Branch("wmet_phi", &wmet_phi,"wmet_phi/D");
	T->Branch("dPhi_metmuon",&dPhi_metmuon,"dPhi_metmuon/D");
	T->Branch("muon_pt", &muon_pt,"muon_pt/D");
	T->Branch("muon_phi", &muon_phi,"muon_phi/D");
	T->Branch("muon_eta", &muon_eta,"muon_eta/D");
	T->Branch("vjet_pt", &vjet_pt);
	T->Branch("vjet_phi", &vjet_phi);
	T->Branch("vjet_eta", &vjet_eta);
	T->Branch("vdR_jetmuon", &vdR_jetmuon);
	T->Branch("KfactorWeight", &KfactorWeight,"KfactorWeight/D");
	T->Branch("mcCrossSection", &mcCrossSection,"mcCrossSection/F");
	T->Branch("mcFilterEfficiency", &mcFilterEfficiency,"mcFilterEfficiency/F");
	T->Branch("puweight", &puweight,"puweight/D");
	T->Branch("mceventweight",&mceventweight,"mceventweight/F");
	T->Branch("totalnAcceptedEvents", &totalnAcceptedEvents,"totalnAcceptedEvents/F");
	T->Branch("totalsumOfEventWeights",&totalsumOfEventWeights, "totalsumOfEventWeights/F");
	T->Branch("isTight",&isTight,"isTight/B");
	T->Branch("muon_HighPtSF",&muon_HighPtSF,"muon_HighPtSF/D");
	T->Branch("SF_trigmu50_highptID",&SF_trigmu50_highptID,"SF_trigmu50_highptID/D");
	T->Branch("muon_HLTtrigmatch",&muon_HLTtrigmatch,"muon_HLTtrigmatch/B");
	
	ifstream inF(origTreeFile);
	TFile *f, *cf;
	TLeaf *nAcc, *sumW;
	string Rpath, treeFile, treeCountFile;
	TTree *eventTree;
	Bool_t isEventClean, HLT_mu50;
	Int_t 			nmuon_nselec, njet_nselec;
	vector<double>  *vmuon_d0sig = 0;
	vector<double>  *vmuon_z0_vrtPVx = 0;
	vector<double>  *vmuon_sintheta = 0;
	vector<float>  	*vmuon_eta = 0;
	vector<float>  	*vmuon_phi = 0;
	vector<float>  	*vmuon_pt = 0;
	vector<bool>  	*vmuon_isoLooseTrack = 0;
	vector<bool>  	*vmuon_isHighPt = 0;
	vector<bool>  	*vmuon_isBadMuon = 0;
	vector<bool>  	*vmuon_isMedium = 0;
	vector<bool>  	*vmuon_isCombined = 0;
	Int_t 			nelec_nselec;
	vector<float>  	*velec_pt = 0;
	vector<float>  	*velec_eta = 0;
	vector<float>  	*velec_phi = 0;
	vector<bool>  	*velec_isCR = 0;
	vector<double>  *velec_d0sig = 0;
	vector<double>  *velec_z0_vrtPVx = 0;
	vector<double>  *velec_sintheta = 0;
	vector<bool>  	*velec_isoLoose = 0;
	vector<bool>  	*velec_isLLHMedium = 0;
	vector<bool>  	*velec_OQ = 0;
	Int_t 			goodMuonFlag = -1;
	Bool_t			vetoElec = 0;
	Bool_t			vetoMuon = 0;
	vector<bool>  	*vmuon_HLTtrigmatch = 0;
	vector<double>  *vmuon_HighPtSF = 0;
	
	Double_t muPtUSF = 1.0;
	Double_t pm = -1;
	
	Int_t GRL = 0, clean = 0, mu50 = 0, muon = 0, MV = 0, EV = 0, MET = 0;
	if(source.substr(0,4).compare("data")!=0 && source.compare("ZWmunu_v9") != 0)
	{
		TTree *eventTreeCount;
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
				eventTree->SetBranchAddress("vmuon_isoLooseTrack",&vmuon_isoLooseTrack);
				eventTree->SetBranchAddress("vmuon_isHighPt",&vmuon_isHighPt);
				eventTree->SetBranchAddress("vmuon_isBadMuon",&vmuon_isBadMuon);
				eventTree->SetBranchAddress("vmuon_isMedium",&vmuon_isMedium);
				eventTree->SetBranchAddress("nelec_nselec",&nelec_nselec);
				eventTree->SetBranchAddress("velec_pt",&velec_pt);
				eventTree->SetBranchAddress("velec_eta",&velec_eta);
				eventTree->SetBranchAddress("velec_phi",&velec_phi);
				eventTree->SetBranchAddress("velec_isCR",&velec_isCR);
				eventTree->SetBranchAddress("velec_d0sig",&velec_d0sig);
				eventTree->SetBranchAddress("velec_z0_vrtPVx",&velec_z0_vrtPVx);
				eventTree->SetBranchAddress("velec_sintheta",&velec_sintheta);
				eventTree->SetBranchAddress("velec_isoLoose",&velec_isoLoose);
				eventTree->SetBranchAddress("velec_isLLHMedium",&velec_isLLHMedium);
				eventTree->SetBranchAddress("velec_OQ",&velec_OQ);
				eventTree->SetBranchAddress("KfactorWeight",&KfactorWeight);
				eventTree->SetBranchAddress("mcCrossSection",&mcCrossSection);
				eventTree->SetBranchAddress("mcFilterEfficiency",&mcFilterEfficiency);
				eventTree->SetBranchAddress("puweight",&puweight);
				eventTree->SetBranchAddress("mceventweight",&mceventweight);
				eventTree->SetBranchAddress("vmuon_HighPtSF",&vmuon_HighPtSF);
				eventTree->SetBranchAddress("SF_trigmu50_highptID",&SF_trigmu50_highptID);
				eventTree->SetBranchAddress("vmuon_isCombined",&vmuon_isCombined);
				eventTree->SetBranchAddress("vmuon_HLTtrigmatch",&vmuon_HLTtrigmatch);

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
					vdR_jetmuon->clear();
					goodMuonFlag = -1;
					vetoElec = 0;
					vetoMuon = 0;
					GRL++;
     				if(!isEventClean) continue;
					clean++;
					if(!HLT_mu50) continue;
					mu50++;
					
					if(muon_pt<100000) muPtUSF = 1.0;
					else if(muon_pt>=100000 && muon_pt<200000) muPtUSF = sqrt(pow(0.02,2)+pow(0.01,2));
					else if(muon_pt>=200000 && muon_pt<300000) muPtUSF = sqrt(pow(0.02,2)+pow(0.01,2)+pow(0.02,2));
					else if(muon_pt>=300000 && muon_pt<400000) muPtUSF = sqrt(pow(0.02,2)+pow(0.01,2)+pow(0.03,2));
					else if(muon_pt>=400000 && muon_pt<500000) muPtUSF = sqrt(pow(0.02,2)+pow(0.01,2)+pow(0.04,2));
					else if(muon_pt>=500000 && muon_pt<600000) muPtUSF = sqrt(pow(0.02,2)+pow(0.01,2)+pow(0.05,2));
					else if(muon_pt>=600000 && muon_pt<700000) muPtUSF = sqrt(pow(0.02,2)+pow(0.01,2)+pow(0.06,2));
					else if(muon_pt>=700000 && muon_pt<800000) muPtUSF = sqrt(pow(0.02,2)+pow(0.01,2)+pow(0.07,2));
					else if(muon_pt>=800000 && muon_pt<900000) muPtUSF = sqrt(pow(0.02,2)+pow(0.01,2)+pow(0.08,2));
					else if(muon_pt>=900000 && muon_pt<1000000) muPtUSF = sqrt(pow(0.02,2)+pow(0.01,2)+pow(0.09,2));
					else muPtUSF = sqrt(pow(0.02,2)+pow(0.01,2)+pow(0.1,2));
					
					muon_pt = muon_pt*(1+pm*muPtUSF);
					
					//Find good muon
					for(Int_t j=0;j<nmuon_nselec;j++)
					{
						if(vmuon_pt->at(j) <= 55000) continue;
						if(fabs(vmuon_d0sig->at(j)) >= 3) continue;
						if(fabs(vmuon_z0_vrtPVx->at(j)*vmuon_sintheta->at(j))>=0.5) continue;
						if(!vmuon_isHighPt->at(j) || vmuon_isBadMuon->at(j)) continue;
						if(fabs(vmuon_eta->at(j)) >= 2.4) continue;
						if(!vmuon_isCombined->at(j)) continue;
						//Take higher pt muon
						if(goodMuonFlag >= 0 && vmuon_pt->at(goodMuonFlag) < vmuon_pt->at(j)) goodMuonFlag = j;
						if(goodMuonFlag == -1) goodMuonFlag = j;
					}
					if(goodMuonFlag == -1) continue;
					muon++;
				
					isTight = vmuon_isoLooseTrack->at(goodMuonFlag);
										
					//Muon veto
					for(Int_t j=0;j<nmuon_nselec;j++)
					{
						if(!((vmuon_isMedium->at(j)||vmuon_isHighPt->at(j)) && vmuon_pt->at(j) > 20000 && j!=goodMuonFlag)) continue;
						if(fabs(vmuon_eta->at(j)) >= 2.4) continue;
						if(fabs(vmuon_d0sig->at(j)) >= 3) continue;
						if(fabs(vmuon_z0_vrtPVx->at(j)*vmuon_sintheta->at(j))>=0.5) continue;
						if(!vmuon_isoLooseTrack->at(j)) continue;
						vetoMuon = 1;
					}
					if(vetoMuon == 1) continue;
					MV++;
				
					for(Int_t j=0;j<nelec_nselec;j++)
					{
						if(velec_pt->at(j)<=20000) continue;
						if(velec_isCR->at(j)) continue;
						if(fabs(velec_eta->at(j)) >= 2.47) continue;
						dR = sqrt(pow(deltaPhi(vmuon_phi->at(goodMuonFlag),velec_phi->at(j)),2)+pow(vmuon_eta->at(goodMuonFlag)-velec_eta->at(j),2));
						if(dR < 0.1) continue;
						if(!velec_isLLHMedium->at(j)) continue;
						if(!velec_isoLoose->at(j)) continue;
						if(!velec_OQ->at(j)) continue;
						vetoElec = 1;
					}
					if(vetoElec == 1) continue;
					EV++;
				
					if(wmet_et <= 55000) continue;
					MET++;
				
					dPhi_metmuon = deltaPhi(wmet_phi,vmuon_phi->at(goodMuonFlag));
					mT = sqrt(2*wmet_et*vmuon_pt->at(goodMuonFlag)*(1-cos(dPhi_metmuon)))/1000;
					muon_pt = vmuon_pt->at(goodMuonFlag)/1000;
					wmet_et = wmet_et/1000;
					muon_eta = vmuon_eta->at(goodMuonFlag);
					muon_phi = vmuon_phi->at(goodMuonFlag);
					muon_HighPtSF = vmuon_HighPtSF->at(goodMuonFlag);
					muon_HLTtrigmatch = vmuon_HLTtrigmatch->at(goodMuonFlag);
					eventWeight = KfactorWeight*mcCrossSection*mcFilterEfficiency*puweight;
					for(Int_t j = 0; j < njet_nselec; j++)
					{
						vjet_pt->at(j) = vjet_pt->at(j)/1000;
						vdR_jetmuon->push_back(sqrt(pow(deltaPhi(vmuon_phi->at(goodMuonFlag),vjet_phi->at(j)),2)+pow(vmuon_eta->at(goodMuonFlag)-vjet_eta->at(j),2)));
					}
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
				eventTree->SetBranchAddress("velec_d0sig",&velec_d0sig);
				eventTree->SetBranchAddress("velec_z0_vrtPVx",&velec_z0_vrtPVx);
				eventTree->SetBranchAddress("velec_sintheta",&velec_sintheta);
				eventTree->SetBranchAddress("velec_isoLoose",&velec_isoLoose);
				eventTree->SetBranchAddress("velec_isLLHMedium",&velec_isLLHMedium);
				eventTree->SetBranchAddress("velec_OQ",&velec_OQ);
				eventTree->SetBranchAddress("KfactorWeight",&KfactorWeight);
				eventTree->SetBranchAddress("mcCrossSection",&mcCrossSection);
				eventTree->SetBranchAddress("mcFilterEfficiency",&mcFilterEfficiency);
				eventTree->SetBranchAddress("puweight",&puweight);
				eventTree->SetBranchAddress("mceventweight",&mceventweight);
				eventTree->SetBranchAddress("vmuon_HighPtSF",&vmuon_HighPtSF);
				eventTree->SetBranchAddress("SF_trigmu50_highptID",&SF_trigmu50_highptID);
				eventTree->SetBranchAddress("vmuon_isCombined",&vmuon_isCombined);
				eventTree->SetBranchAddress("vmuon_HLTtrigmatch",&vmuon_HLTtrigmatch);

				Long64_t nentries = eventTree->GetEntries();
				for (Long64_t i=0;i<nentries;i++) 
				{
     				eventTree->GetEntry(i);
					vdR_jetmuon->clear();
					goodMuonFlag = -1;
					vetoElec = 0;
					vetoMuon = 0;
					GRL++;
     				if(!isEventClean) continue;
					clean++;
					if(!HLT_mu50) continue;
					mu50++;
					
					//Find good muon
					for(Int_t j=0;j<nmuon_nselec;j++)
					{
						if(vmuon_pt->at(j) <= 55000) continue;
						if(fabs(vmuon_d0sig->at(j)) >= 3) continue;
						if(fabs(vmuon_z0_vrtPVx->at(j)*vmuon_sintheta->at(j))>=0.5) continue;
						if(!vmuon_isHighPt->at(j) || vmuon_isBadMuon->at(j)) continue;
						if(fabs(vmuon_eta->at(j)) >= 2.4) continue;
						if(!vmuon_isCombined->at(j)) continue;
						//Take higher pt muon
						if(goodMuonFlag >= 0 && vmuon_pt->at(goodMuonFlag) < vmuon_pt->at(j)) goodMuonFlag = j;
						if(goodMuonFlag == -1) goodMuonFlag = j;
					}
					if(goodMuonFlag == -1) continue;
					muon++;
				
					isTight = vmuon_isoLooseTrack->at(goodMuonFlag);
										
					//Muon veto
					for(Int_t j=0;j<nmuon_nselec;j++)
					{
						if(!((vmuon_isMedium->at(j)||vmuon_isHighPt->at(j)) && vmuon_pt->at(j) > 20000 && j!=goodMuonFlag)) continue;
						if(fabs(vmuon_eta->at(j)) >= 2.4) continue;
						if(fabs(vmuon_d0sig->at(j)) >= 3) continue;
						if(fabs(vmuon_z0_vrtPVx->at(j)*vmuon_sintheta->at(j))>=0.5) continue;
						if(!vmuon_isoLooseTrack->at(j)) continue;
						vetoMuon = 1;
					}
					if(vetoMuon == 1) continue;
					MV++;
				
					for(Int_t j=0;j<nelec_nselec;j++)
					{
						if(velec_pt->at(j)<=20000) continue;
						if(velec_isCR->at(j)) continue;
						if(fabs(velec_eta->at(j)) >= 2.47) continue;
						dR = sqrt(pow(deltaPhi(vmuon_phi->at(goodMuonFlag),velec_phi->at(j)),2)+pow(vmuon_eta->at(goodMuonFlag)-velec_eta->at(j),2));
						if(dR < 0.1) continue;
						if(!velec_isLLHMedium->at(j)) continue;
						if(!velec_isoLoose->at(j)) continue;
						if(!velec_OQ->at(j)) continue;
						vetoElec = 1;
					}
					if(vetoElec == 1) continue;
					EV++;
				
					if(wmet_et <= 55000) continue;
					MET++;
				
					dPhi_metmuon = deltaPhi(wmet_phi,vmuon_phi->at(goodMuonFlag));
					mT = sqrt(2*wmet_et*vmuon_pt->at(goodMuonFlag)*(1-cos(dPhi_metmuon)))/1000;
					muon_pt = vmuon_pt->at(goodMuonFlag)/1000;
					wmet_et = wmet_et/1000;
					muon_eta = vmuon_eta->at(goodMuonFlag);
					muon_phi = vmuon_phi->at(goodMuonFlag);
					muon_HighPtSF = vmuon_HighPtSF->at(goodMuonFlag);
					muon_HLTtrigmatch = vmuon_HLTtrigmatch->at(goodMuonFlag);
					eventWeight = KfactorWeight*mcCrossSection*mcFilterEfficiency*puweight;
					for(Int_t j = 0; j < njet_nselec; j++)
					{
						vjet_pt->at(j) = vjet_pt->at(j)/1000;
						vdR_jetmuon->push_back(sqrt(pow(deltaPhi(vmuon_phi->at(goodMuonFlag),vjet_phi->at(j)),2)+pow(vmuon_eta->at(goodMuonFlag)-vjet_eta->at(j),2)));
					}
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
