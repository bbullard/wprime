#ifndef Plotter_h
#define Plotter_h

class Plotter {
public :
   	Plotter(const string plotVar = "mT", Double_t L = 20, Double_t xmin = 0, Double_t xmax = 2500, Bool_t useLogX = 0);
	
	TCanvas *splitCan;
	TDirectory *dir;
	string plotVar;
	Double_t Lumi;
	Bool_t useLogX;
	Double_t QCDd0sig_min;
	TCut d0sigCut;
	string fakeEffParamName;
	
	Int_t nxbins;
	Double_t *xbins;
	
	TTree *WTree;
	TTree *TopTree;
	TTree *DYTree;
	TTree *DibTree;
	TTree *dataTree;
	
	TTree *W2000Tree;
	TTree *W3000Tree;
	TTree *W4000Tree;
	TTree *W5000Tree;
	
	TTree *WQCDTree;
	TTree *TopQCDTree;
	TTree *DYQCDTree;
	TTree *DibQCDTree;
	TTree *dataQCDTree;
	
	//Intitialize Histograms
	TH1F *W;
	TH1F *Top;
	TH1F *DY;
	TH1F *Dib;
	TH1F *QCD;
	TH1F *data;
	TH1F *ratio;
	
	TH1F *dataLoose;
	
	TH1F *W2000;
	TH1F *W3000;
	TH1F *W4000;
	TH1F *W5000;
	
	TH1F *eff;
	
	TH1F *WQCDt;
	TH1F *TopQCDt;
	TH1F *DYQCDt;
	TH1F *DibQCDt;
	TH1F *dataQCDt;
	TH1F *QCDt;
	
	TH1F *WQCDl;
	TH1F *TopQCDl;
	TH1F *DYQCDl;
	TH1F *DibQCDl;
	TH1F *dataQCDl;
	TH1F *QCDl;
	
	void FullPlot();
	
	void MakeFakeEff(Double_t xmin = 0, Double_t xmax = 100);
	
	void PlotIntermediates(const string xTitle = "P_{T} [GeV]", const string saveVar = "muon_pt");
	
	void MakeQCDPlot(TH1F *hist, string var);
	
	void MakeRatioPlot();
	
	void PlotWPrimeKinematics(Double_t eta_c = 0);
	
	void OutputBinContents(TH1F *hist, string mod = "");
	
	void MakeHistVector(vector<TH1F*> *histVector, string name, Int_t bins, Double_t xmin, Double_t xmax);
	
	void PlotSignalmT();
	
	void MakeRatioBounds(TH1F *ratioLB, TH1F *ratioUB);
	
};

#endif
