#include "Plotter.h"

#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TDirectory.h"
#include "TH2F.h"
#include "TH1F.h"
#include "THStack.h"
#include "TLegend.h"

#include <iostream>

using namespace std;

Plotter::Plotter(const string plotVar, Double_t L, Double_t xmin, Double_t xmax, Bool_t useLogX)
{
	dir = gDirectory;   
	
	QCDd0sig_min = 1.5;
	d0sigCut = TCut(Form("muon_d0sig>%f",QCDd0sig_min));
	
	Lumi = L;
	
	this->useLogX = useLogX;
	this->plotVar = plotVar;
	
	nxbins = 50;
	if(useLogX) nxbins = 50;
	xbins = new Double_t[nxbins+1];
	Double_t binwidth;
	if(useLogX)
	{
		Double_t logxmin = TMath::Log10(xmin);
		Double_t logxmax = TMath::Log10(xmax);
		cout<<logxmin<<" "<<logxmax<<endl;
		binwidth = (logxmax-logxmin)/Double_t(nxbins);
		xbins[0] = xmin;
		for(Int_t i = 1; i <= nxbins; i++) xbins[i] = TMath::Power(10, logxmin+i*binwidth);
	}
	else
	{
		binwidth = (xmax-xmin)/Double_t(nxbins);
		xbins[0] = xmin;
		for(Int_t i = 1; i <= nxbins; i++) xbins[i] = xmin+i*binwidth;
	}
	
	W = new TH1F("W","",nxbins,xbins);
	Top = new TH1F("Top","",nxbins,xbins);
	DY = new TH1F("DY","",nxbins,xbins);
	Dib = new TH1F("Dib","",nxbins,xbins);
	QCD = new TH1F("QCD","",nxbins,xbins);
	data = new TH1F("data","",nxbins,xbins);
	ratio = new TH1F("ratio","",nxbins,xbins);
	
	dataLoose = new TH1F("dataLoose","",nxbins,xbins);
	
	W2000 = new TH1F("W2000","",nxbins,xbins);
	W3000 = new TH1F("W3000","",nxbins,xbins);
	W4000 = new TH1F("W4000","",nxbins,xbins);
	W5000 = new TH1F("W5000","",nxbins,xbins);
	
	TFile *f;	
	TF1 *f1 = new TF1("f1","1",xmin,xmax);
	string input;
	
	
	f = TFile::Open("Files/ntuples/selection/WFullSelection.root","read");
	dir->cd();
	WTree = (TTree*)f->Get("tree");
	f = TFile::Open("Files/ntuples/selection/TopFullSelection.root","read");
	dir->cd();
	TopTree = (TTree*)f->Get("tree");
	f = TFile::Open("Files/ntuples/selection/DYFullSelection.root","read");
	dir->cd();
	DYTree = (TTree*)f->Get("tree");
	f = TFile::Open("Files/ntuples/selection/DibFullSelection.root","read");
	dir->cd();
	DibTree = (TTree*)f->Get("tree");
	f = TFile::Open("Files/ntuples/selection/dataFullSelection.root","read");
	dir->cd();
	dataTree = (TTree*)f->Get("tree");

	f = TFile::Open("Files/ntuples/selection/signal2FullSelection.root","read");
	dir->cd();
	W2000Tree = (TTree*)f->Get("tree");
	f = TFile::Open("Files/ntuples/selection/signal3FullSelection.root","read");
	dir->cd();
	W3000Tree = (TTree*)f->Get("tree");
	f = TFile::Open("Files/ntuples/selection/signal4FullSelection.root","read");
	dir->cd();
	W4000Tree = (TTree*)f->Get("tree");
	f = TFile::Open("Files/ntuples/selection/signal5FullSelection.root","read");
	dir->cd();
	W5000Tree = (TTree*)f->Get("tree");
	
	f = TFile::Open("Files/ntuples/selection/WQCD.root","read");
	dir->cd();
	WQCDTree = (TTree*)f->Get("tree");
	f = TFile::Open("Files/ntuples/selection/TopQCD.root","read");
	dir->cd();
	TopQCDTree = (TTree*)f->Get("tree");
	f = TFile::Open("Files/ntuples/selection/DYQCD.root","read");
	dir->cd();
	DYQCDTree = (TTree*)f->Get("tree");
	f = TFile::Open("Files/ntuples/selection/DibQCD.root","read");
	dir->cd();
	DibQCDTree = (TTree*)f->Get("tree");
	f = TFile::Open("Files/ntuples/selection/dataQCD.root","read");
	dir->cd();
	dataQCDTree = (TTree*)f->Get("tree");
	
	input=plotVar+">>W";
	WTree->Draw(input.c_str(),"SF_trigmu50_highptID*muon_HighPtSF*mcCrossSection*mcFilterEfficiency*mceventweight*puweight*KfactorWeight*isTight/totalsumOfEventWeights");
	W->Multiply(f1,Lumi*1e6);
	
	input=plotVar+">>Top";	
	TopTree->Draw(input.c_str(),"SF_trigmu50_highptID*muon_HighPtSF*mcCrossSection*mcFilterEfficiency*mceventweight*puweight*KfactorWeight*isTight/totalsumOfEventWeights");
	Top->Multiply(f1,Lumi*1e6);

	input=plotVar+">>DY";
	DYTree->Draw(input.c_str(),"SF_trigmu50_highptID*muon_HighPtSF*mcCrossSection*mcFilterEfficiency*mceventweight*puweight*KfactorWeight*isTight/totalsumOfEventWeights");
	DY->Multiply(f1,Lumi*1e6);
	
	input=plotVar+">>Dib";
	DibTree->Draw(input.c_str(),"SF_trigmu50_highptID*muon_HighPtSF*mcCrossSection*mcFilterEfficiency*mceventweight*puweight*KfactorWeight*isTight/totalsumOfEventWeights");
	Dib->Multiply(f1,Lumi*1e6);
	

	input=plotVar+">>data";
	dataTree->Draw(input.c_str(),"isTight");
	
	input=plotVar+">>dataLoose";	
	dataTree->Draw(input.c_str(),"");
	
	input=plotVar+">>W2000";
	W2000Tree->Draw(input.c_str(),"SF_trigmu50_highptID*muon_HighPtSF*mcCrossSection*mcFilterEfficiency*mceventweight*puweight*KfactorWeight*isTight/totalsumOfEventWeights");
	W2000->Multiply(f1,Lumi*1e6);
	
	input=plotVar+">>W3000";
	W3000Tree->Draw(input.c_str(),"SF_trigmu50_highptID*muon_HighPtSF*mcCrossSection*mcFilterEfficiency*mceventweight*puweight*KfactorWeight*isTight/totalsumOfEventWeights");
	W3000->Multiply(f1,Lumi*1e6);
	
	input=plotVar+">>W4000";
	W4000Tree->Draw(input.c_str(),"SF_trigmu50_highptID*muon_HighPtSF*mcCrossSection*mcFilterEfficiency*mceventweight*puweight*KfactorWeight*isTight/totalsumOfEventWeights");
	W4000->Multiply(f1,Lumi*1e6);
	
	input=plotVar+">>W5000";
	W5000Tree->Draw(input.c_str(),"SF_trigmu50_highptID*muon_HighPtSF*mcCrossSection*mcFilterEfficiency*mceventweight*puweight*KfactorWeight*isTight/totalsumOfEventWeights");
	W5000->Multiply(f1,Lumi*1e6);
	
	
	splitCan = new TCanvas("splitCan","splitCan",800,800);
	gStyle->SetOptStat(0);
	gStyle->SetLegendBorderSize(0);
}

void Plotter::FullPlot()
{
	fakeEffParamName  = "muon_pt";
	MakeFakeEff(50,350);
	//string fakeEffParamName = "muon_eta";
	//MakeFakeEff(-2.7,2.7);
	MakeQCDPlot(QCD,plotVar);
	MakeRatioPlot();
	TH1F *ratioUB = (TH1F*)ratio->Clone("ratioUB");
	TH1F *ratioLB = (TH1F*)ratio->Clone("ratioLB");
	MakeRatioBounds(ratioLB,ratioUB);
	ratioUB->SetFillColor(1);
	ratioLB->SetFillColor(10);
	ratioUB->SetLineColor(0);
	ratioLB->SetLineColor(0);
	ratioUB->SetFillStyle(3002);
	ratioLB->SetFillStyle(3002);
	
	W2000->SetLineColor(kAzure+10);
	W3000->SetLineColor(kAzure);
	W4000->SetLineColor(kViolet);
	W5000->SetLineColor(kPink);
	
	W2000->SetLineWidth(2);
	W3000->SetLineWidth(2);
	W4000->SetLineWidth(2);
	W5000->SetLineWidth(2);
	
	W->SetLineColor(kBlack);
	Top->SetLineColor(kBlack);
	DY->SetLineColor(kBlack);
	Dib->SetLineColor(kBlack);
	QCD->SetLineColor(kBlack);
	data->SetLineColor(kBlack);
	
	W->SetLineWidth(2);
	Top->SetLineWidth(2);
	DY->SetLineWidth(2);
	Dib->SetLineWidth(2);
	QCD->SetLineWidth(2);
	
	W->SetFillColor(kWhite);
	Top->SetFillColor(kRed+1);
	DY->SetFillColor(kAzure-9);
	Dib->SetFillColor(kOrange);
	QCD->SetFillColor(kYellow-9);
	
	THStack *stack = new THStack("stack","");
	stack->Add(QCD,"hist");
	stack->Add(Dib,"hist");
	stack->Add(DY,"hist");
	stack->Add(Top,"hist");
	stack->Add(W,"hist");
	
	stack->SetMaximum(Lumi*5*1e6);
	stack->SetMinimum(Lumi*5e-4);
	
	if(plotVar.compare("mT")==0) stack->SetTitle(";M_{T} [GeV];Events");
	if(plotVar.compare("muon_pt")==0) stack->SetTitle(";Muon P_{T} [GeV];Events");
	if(plotVar.compare("muon_eta")==0) stack->SetTitle(";Muon #eta;Events");
	if(plotVar.compare("wmet_et")==0) stack->SetTitle(";E_{T}^{miss} [GeV];Events");
	if(plotVar.compare("abs(dPhi_metmuon)")==0) {stack->SetTitle(";#Delta#phi_{#mu,E_{T}^{miss}};Events"); stack->SetMaximum(Lumi*3e8);}
	if(plotVar.compare("cos(dPhi_metmuon)")==0) {stack->SetTitle(";#cos(#Delta#phi_{#mu,E_{T}^{miss}});Events"); stack->SetMaximum(Lumi*1e8);}
	
	
	TLegend *leg = new TLegend(0.75,0.6,0.92,0.9);
	leg->AddEntry(data,"Data","pe");
	leg->AddEntry(W,"W","f");
	leg->AddEntry(Top,"Top","f");
	leg->AddEntry(DY,"Z/#gamma*","f");
	leg->AddEntry(Dib,"Diboson","f");
	leg->AddEntry(QCD,"QCD","f");
	leg->SetTextFont(42);
	
	TLegend *sleg = new TLegend(0.58,0.7,0.75,0.9);
	sleg->AddEntry(W2000,"W' (2 TeV)","l");
	sleg->AddEntry(W3000,"W' (3 TeV)","l");
	sleg->AddEntry(W4000,"W' (4 TeV)","l");
	sleg->AddEntry(W5000,"W' (5 TeV)","l");
	sleg->SetTextFont(42);
	
	TPaveText *glum = new TPaveText(.35,.75,.55,.8,"brNDC");
   	glum->AddText(Form("#int Ldt = %.1f fb^{-1}",Lumi));
   	glum->SetBorderSize(0);
	glum->SetFillColor(0);
	
	TPaveText *ATLASinternal = new TPaveText(.35,.85,.55,.9,"brNDC");
   	ATLASinternal->AddText("#bf{#it{ATLAS}} Internal");
   	ATLASinternal->SetBorderSize(0);
	ATLASinternal->SetFillColor(0);
	
	splitCan->cd();
	TPad plotPad("plotPad", "plotPad", 0, 0.3, 1, 1.0);
	if(useLogX) plotPad.SetLogx();
   	plotPad.SetBottomMargin(0); 
	plotPad.SetLogy();
	plotPad.SetTickx();
	plotPad.SetTicky();                 
	plotPad.Draw();
	plotPad.cd();
	stack->Draw();
	data->Draw("same pe");
	W2000->Draw("same hist");
	W3000->Draw("same hist");
	W4000->Draw("same hist");
	W5000->Draw("same hist");
	leg->Draw();
	sleg->Draw();
	glum->Draw();
	ATLASinternal->Draw();
	plotPad.RedrawAxis();
	
	splitCan->cd();
	TPad ratioPad("ratioPad", "ratioPad", 0, 0.05, 1, 0.3);
	if(useLogX) ratioPad.SetLogx();
   	ratioPad.SetTopMargin(0);
	ratioPad.SetBottomMargin(0.4);  
	ratioPad.SetTickx();
	ratioPad.SetTicky(); 
	ratioPad.SetGridy();
	ratioPad.Draw();
	ratioPad.cd();
	ratioUB->Draw("hist");
	ratioLB->Draw("hist same");
	ratio->Draw("pe same");
	ratioPad.RedrawAxis();
	
	splitCan->Print(Form("/export/home/bbullard/thesis/plots/%.1f/%s_%.1f_%.0f-%.0f_%i_%s.png",Lumi,plotVar.c_str(),Lumi,xbins[0],xbins[nxbins],useLogX,fakeEffParamName.c_str()));
	splitCan->Print(Form("/export/home/bbullard/thesis/plots/%.1f/%s_%.1f_%.0f-%.0f_%i_%s.eps",Lumi,plotVar.c_str(),Lumi,xbins[0],xbins[nxbins],useLogX,fakeEffParamName.c_str()));
	cout<<endl;
	
	string csvOutMod = "d0Psig" + to_string(QCDd0sig_min);
	OutputBinContents(QCD,csvOutMod);
}

void Plotter::MakeFakeEff(Double_t xmin, Double_t xmax)
{
	TCanvas can("can","can",800,600);
	can.cd();
	string var = fakeEffParamName;
	Int_t nbins = 8;
	Double_t bins[nbins+1];
	//Double_t binwidth = (xmax-xmin)/Double_t(nbins);
	bins[0] = 55;
	bins[1] = 60;
	bins[2] = 70;
	bins[3] = 85;
	bins[4] = 100;
	bins[5] = 125;
	bins[6] = 150;
	bins[7] = 175;
	bins[8] = 1500;
	//for(Int_t i = 1; i <= nbins; i++) bins[i] = xmin+i*binwidth;
	
	WQCDt = new TH1F("WQCDt","",nbins,bins);
	TopQCDt = new TH1F("TopQCDt","",nbins,bins);
	DYQCDt = new TH1F("DYQCDt","",nbins,bins);
	DibQCDt = new TH1F("DibQCDt","",nbins,bins);
	dataQCDt = new TH1F("dataQCDt","",nbins,bins);
	QCDt = new TH1F("QCDt","",nbins,bins);
	
	WQCDl = new TH1F("WQCDl","",nbins,bins);
	TopQCDl = new TH1F("TopQCDl","",nbins,bins);
	DYQCDl = new TH1F("DYQCDl","",nbins,bins);
	DibQCDl = new TH1F("DibQCDl","",nbins,bins);
	dataQCDl = new TH1F("dataQCDl","",nbins,bins);
	QCDl = new TH1F("QCDl","",nbins,bins);
	
	eff = new TH1F("eff","",nbins,bins);
	
	TF1 *f1 = new TF1("f1","1",bins[0],bins[nbins]);
	string input;

	input=var+">>WQCDt";
	WQCDTree->Draw(input.c_str(),d0sigCut*"mcCrossSection*mcFilterEfficiency*mceventweight*puweight*KfactorWeight*isTight/totalsumOfEventWeights");
	WQCDt->Multiply(f1,Lumi*1e6);
	
	input=var+">>TopQCDt";
	TopQCDTree->Draw(input.c_str(),d0sigCut*"mcCrossSection*mcFilterEfficiency*mceventweight*puweight*KfactorWeight*isTight/totalsumOfEventWeights");
	TopQCDt->Multiply(f1,Lumi*1e6);
	
	input=var+">>DYQCDt";
	DYQCDTree->Draw(input.c_str(),d0sigCut*"mcCrossSection*mcFilterEfficiency*mceventweight*puweight*KfactorWeight*isTight/totalsumOfEventWeights");
	DYQCDt->Multiply(f1,Lumi*1e6);
	
	input=var+">>DibQCDt";
	DibQCDTree->Draw(input.c_str(),d0sigCut*"mcCrossSection*mcFilterEfficiency*mceventweight*puweight*KfactorWeight*isTight/totalsumOfEventWeights");
	DibQCDt->Multiply(f1,Lumi*1e6);

	input=var+">>dataQCDt";		
	dataQCDTree->Draw(input.c_str(),d0sigCut*"isTight");
	
	
	input=var+">>WQCDl";
	WQCDTree->Draw(input.c_str(),d0sigCut*"mcCrossSection*mcFilterEfficiency*mceventweight*puweight*KfactorWeight/totalsumOfEventWeights");
	WQCDl->Multiply(f1,Lumi*1e6);
	
	input=var+">>TopQCDl";
	TopQCDTree->Draw(input.c_str(),d0sigCut*"mcCrossSection*mcFilterEfficiency*mceventweight*puweight*KfactorWeight/totalsumOfEventWeights");
	TopQCDl->Multiply(f1,Lumi*1e6);

	input=var+">>DYQCDl";
	DYQCDTree->Draw(input.c_str(),d0sigCut*"mcCrossSection*mcFilterEfficiency*mceventweight*puweight*KfactorWeight/totalsumOfEventWeights");
	DYQCDl->Multiply(f1,Lumi*1e6);

	input=var+">>DibQCDl";
	DibQCDTree->Draw(input.c_str(),d0sigCut*"mcCrossSection*mcFilterEfficiency*mceventweight*puweight*KfactorWeight/totalsumOfEventWeights");
	DibQCDl->Multiply(f1,Lumi*1e6);

	input=var+">>dataQCDl";	
	dataQCDTree->Draw(input.c_str(),d0sigCut);
	
	WQCDt->SetLineWidth(2);
	TopQCDt->SetLineWidth(2);
	DYQCDt->SetLineWidth(2);
	DibQCDt->SetLineWidth(2);
	dataQCDt->SetLineWidth(2);
	QCDt->SetLineWidth(2);
	
	WQCDl->SetLineWidth(2);
	TopQCDl->SetLineWidth(2);
	DYQCDl->SetLineWidth(2);
	DibQCDl->SetLineWidth(2);
	dataQCDl->SetLineWidth(2);
	QCDl->SetLineWidth(2);
	
	string xtitle;
	if(var.compare("muon_pt")==0) xtitle = "P_{T} [GeV]";
	if(var.compare("muon_eta")==0) xtitle = "#eta";
	
	Double_t MC = 0.0;
	for(Int_t bin=1; bin <= nbins; bin++)
	{
		MC += WQCDt->GetBinContent(bin);
		MC += TopQCDt->GetBinContent(bin);
		MC += DYQCDt->GetBinContent(bin);
		MC += DibQCDt->GetBinContent(bin);
		QCDt->SetBinContent(bin,dataQCDt->GetBinContent(bin)-MC);
		MC = 0.0;
	}

	for(Int_t bin=1; bin <= nbins; bin++)
	{
		MC += WQCDl->GetBinContent(bin);
		MC += TopQCDl->GetBinContent(bin);
		MC += DYQCDl->GetBinContent(bin);
		MC += DibQCDl->GetBinContent(bin);
		QCDl->SetBinContent(bin,dataQCDl->GetBinContent(bin)-MC);
		MC = 0.0;
	}
	
	eff->Sumw2();
	eff->SetLineWidth(2);
	for(Int_t bin=1; bin <= nbins; bin++)
	{
		if(QCDl->GetBinContent(bin) == 0) continue;
		eff->SetBinContent(bin,QCDt->GetBinContent(bin)/QCDl->GetBinContent(bin));
	}
	
	PlotIntermediates(xtitle,var);
}

void Plotter::PlotIntermediates(const string xTitle, const string saveVar)
{
	TCanvas can("can","can",1600,600);
	can.Divide(2,1);
	can.GetPad(1)->SetGrid();
	can.GetPad(1)->SetLogy();
	can.GetPad(1)->SetLogx();
	can.GetPad(2)->SetGrid();
	can.GetPad(2)->SetLogy();
	can.GetPad(2)->SetLogx();
	string input;
	
	WQCDl->SetLineColor(kBlack);
	TopQCDl->SetLineColor(kBlack);
	DYQCDl->SetLineColor(kBlack);
	DibQCDl->SetLineColor(kBlack);
	dataQCDl->SetLineColor(kBlack);
	
	WQCDl->SetFillColor(kWhite);
	TopQCDl->SetFillColor(kRed+1);
	DYQCDl->SetFillColor(kAzure-9);
	DibQCDl->SetFillColor(kOrange);
	dataQCDl->SetFillColor(kWhite);
	
	WQCDt->SetLineColor(kBlack);
	TopQCDt->SetLineColor(kBlack);
	DYQCDt->SetLineColor(kBlack);
	DibQCDt->SetLineColor(kBlack);
	dataQCDt->SetLineColor(kBlack);
	
	WQCDt->SetFillColor(kWhite);
	TopQCDt->SetFillColor(kRed+1);
	DYQCDt->SetFillColor(kAzure-9);
	DibQCDt->SetFillColor(kOrange);
	dataQCDt->SetFillColor(kWhite);
	
	dataQCDl->SetTitle("Loose QCD Control Region");
	dataQCDt->SetTitle("Tight QCD Control Region");
	dataQCDl->SetYTitle("Events");
	dataQCDt->SetYTitle("Events");
	dataQCDl->SetXTitle(xTitle.c_str());
	dataQCDt->SetXTitle(xTitle.c_str());
	dataQCDl->SetMaximum(5e7);
	dataQCDt->SetMaximum(5e7);
	if(saveVar.compare("muon_eta") == 0) 
	{
		dataQCDl->SetMaximum(1e7);
		dataQCDt->SetMaximum(1e7);
	}
	dataQCDl->SetMinimum(1);
	dataQCDt->SetMinimum(1);
	
	THStack *stackl = new THStack("stackl","");
	stackl->Add(DibQCDl,"hist");
	stackl->Add(DYQCDl,"hist");
	stackl->Add(TopQCDl,"hist");
	stackl->Add(WQCDl,"hist");
	
	THStack *stackt = new THStack("stackt","");
	stackt->Add(DibQCDt,"hist");
	stackt->Add(DYQCDt,"hist");
	stackt->Add(TopQCDt,"hist");
	stackt->Add(WQCDt,"hist");
	
	TLegend *leg = new TLegend(0.75,0.6,0.92,0.9);
	leg->AddEntry(dataQCDl,"Data","p");
	leg->AddEntry(WQCDl,"W","f");
	leg->AddEntry(TopQCDl,"Top","f");
	leg->AddEntry(DYQCDl,"Z/#gamma*","f");
	leg->AddEntry(DibQCDl,"Diboson","f");
	leg->SetTextFont(42);
	
	TPaveText *ATLASinternal = new TPaveText(.425,.83,.675,.9,"brNDC");
   	ATLASinternal->AddText("#bf{#it{ATLAS}} Internal");
   	ATLASinternal->SetBorderSize(0);
	ATLASinternal->SetFillColor(0);
	
	can.cd(1);
	dataQCDl->Draw("pe");
	stackl->Draw("same");
	leg->Draw();
	ATLASinternal->Draw();
	can.cd(2);
	dataQCDt->Draw("pe");
	stackt->Draw("same");
	leg->Draw("same");
	ATLASinternal->Draw();
	input = Form("/export/home/bbullard/thesis/plots/%.1f/",Lumi)+saveVar+"_QCDcontrolRegion.png";
	can.Print(input.c_str());
	
	
	TCanvas can2("can2","can2",800,900);
	QCDt->SetLineColor(kRed+1);
	QCDl->SetLineColor(kBlack);
	QCDt->SetMarkerColor(kRed+1);
	QCDl->SetMarkerColor(kBlack);
	QCDl->SetYTitle("Events");
	
	if(saveVar.compare("muon_eta") == 0) 
	{
		QCDl->SetMaximum(1e6);
		QCDl->SetMinimum(1e2);
	}
	
	QCDl->SetMinimum(1e3);
	
	TLegend *leg2 = new TLegend(0.75,.65,0.92,0.9);
	leg2->AddEntry(QCDl,"Loose","l");
	leg2->AddEntry(QCDt,"Tight","l");
	
	TPaveText *ATLASinternal2 = new TPaveText(.425,.75,.675,.85,"brNDC");
   	ATLASinternal2->AddText("#bf{#it{ATLAS}} Internal");
   	ATLASinternal2->SetBorderSize(0);
	ATLASinternal2->SetFillColor(0);
	
	eff->SetTitle(";;#epsilon_{fake}");
	eff->SetXTitle(xTitle.c_str());
	eff->SetMinimum(0);
	eff->SetMaximum(1);
	eff->GetXaxis()->SetTitleSize(0.13);
	eff->GetXaxis()->SetLabelSize(0.1);
	eff->GetXaxis()->SetTitleOffset(0.95);
	eff->GetYaxis()->SetTitleSize(0.15);
	eff->GetYaxis()->SetLabelSize(0.1);
	eff->GetYaxis()->SetTitleOffset(0.45);
	eff->GetYaxis()->SetNdivisions(505);
	eff->SetLineColor(kBlack);
	if(saveVar.compare("muon_eta") == 0)
	{ 
		eff->SetMaximum(0.6);
		eff->GetYaxis()->SetNdivisions(503);
	}
	
	can2.cd();
	TPad plotPad("plotPad", "plotPad", 0, 0.4, 1, 1.0);
	plotPad.SetLogy();
	plotPad.SetLogx();
	plotPad.SetBottomMargin(0);
	plotPad.SetTickx();
	plotPad.SetTicky(); 
	plotPad.Draw();
	plotPad.cd();
	QCDl->Draw("pe");
	QCDt->Draw("pe same");
	leg2->Draw();
	ATLASinternal2->Draw();
	
	can2.cd();
	TPad effPad("effPad","effPad", 0, 0.05, 1, 0.4);
	effPad.SetLogx();
	effPad.SetTopMargin(0);
	effPad.SetBottomMargin(0.4);
	effPad.SetTickx();
	effPad.SetTicky();
	effPad.SetGridy();
	effPad.Draw();
	effPad.cd();
	//eff->Draw("pe");
	eff->Draw("hist same");
	input = Form("/export/home/bbullard/thesis/plots/%.1f/",Lumi)+saveVar+"_fakeEfficiency.png";
	can2.Print(input.c_str());
	
	TCanvas can3("can3","can3",800,600);
	can3.SetLogy();
	data->SetLineColor(kRed+1);
	dataLoose->SetTitle(";M_{T} [GeV];Events");
	can3.cd();
	dataLoose->Draw("hist");
	data->Draw("same hist");
	leg2->Draw();
	input = Form("/export/home/bbullard/thesis/plots/%.1f/%s_%.1f_%.0f-%.0f_%i_looseAndTight.png",Lumi,plotVar.c_str(),Lumi,xbins[0],xbins[nxbins],useLogX);
	can3.Print(input.c_str());
	
	string csvOutMod = "_"+saveVar;//+"_noTop";
	
	/*OutputBinContents(data);
	OutputBinContents(dataLoose);
	OutputBinContents(W);
	OutputBinContents(Top);
	OutputBinContents(DY);
	OutputBinContents(Dib);
	OutputBinContents(QCD,csvOutMod);
	
	OutputBinContents(eff,csvOutMod);
	
	OutputBinContents(dataQCDt,csvOutMod);
	OutputBinContents(WQCDt,csvOutMod);
	OutputBinContents(TopQCDt,csvOutMod);
	OutputBinContents(DYQCDt,csvOutMod);
	OutputBinContents(DibQCDt,csvOutMod);
	OutputBinContents(QCDt,csvOutMod);
	
	OutputBinContents(dataQCDl,csvOutMod);
	OutputBinContents(WQCDl,csvOutMod);
	OutputBinContents(TopQCDl,csvOutMod);
	OutputBinContents(DYQCDl,csvOutMod);
	OutputBinContents(DibQCDl,csvOutMod);
	OutputBinContents(QCDl,csvOutMod);
	*/
}

void Plotter::MakeQCDPlot(TH1F *hist, string var)
{
	Double_t v, param;
	Char_t isTight;
	Double_t e;
	
	if(var.compare("muon_pt") == 0)
	{
		dataTree->SetBranchAddress(var.c_str(),&v);
		dataTree->SetBranchAddress("isTight",&isTight);
		
		Long64_t nentries = dataTree->GetEntries();
		for (Long64_t i=0;i<nentries;i++) 
		{
    	 	dataTree->GetEntry(i);
			e = eff->GetBinContent(eff->GetXaxis()->FindBin(v));
			hist->Fill(v,e/(1-e)*(1-isTight));
  		} 
	}
	else if(var.compare("cos(dPhi_metmuon)") == 0)
	{
		var="dPhi_metmuon";
		dataTree->SetBranchAddress(var.c_str(),&v);
		dataTree->SetBranchAddress(fakeEffParamName.c_str(),&param);
		dataTree->SetBranchAddress("isTight",&isTight);
		
		Long64_t nentries = dataTree->GetEntries();
		for (Long64_t i=0;i<nentries;i++) 
		{
    	 	dataTree->GetEntry(i);
			e = eff->GetBinContent(eff->GetXaxis()->FindBin(param));
			hist->Fill(cos(v),e/(1-e)*(1-isTight));
  		} 
	}
	else if(var.compare("abs(dPhi_metmuon)") == 0)
	{
		var="dPhi_metmuon";
		dataTree->SetBranchAddress(var.c_str(),&v);
		dataTree->SetBranchAddress(fakeEffParamName.c_str(),&param);
		dataTree->SetBranchAddress("isTight",&isTight);
		
		Long64_t nentries = dataTree->GetEntries();
		for (Long64_t i=0;i<nentries;i++) 
		{
    	 	dataTree->GetEntry(i);
			e = eff->GetBinContent(eff->GetXaxis()->FindBin(param));
			hist->Fill(abs(v),e/(1-e)*(1-isTight));
  		} 
	}
	else
	{
		dataTree->SetBranchAddress(var.c_str(),&v);
		dataTree->SetBranchAddress(fakeEffParamName.c_str(),&param);
		dataTree->SetBranchAddress("isTight",&isTight);
		
		Long64_t nentries = dataTree->GetEntries();
		for (Long64_t i=0;i<nentries;i++) 
		{
    	 	dataTree->GetEntry(i);
			e = eff->GetBinContent(eff->GetXaxis()->FindBin(param));
			hist->Fill(v,e/(1-e)*(1-isTight));
  		} 
	}
	TCanvas can("can","can",800,600);
	can.cd();
	hist->Draw();
	string outstring = var+"testQCD.png";
	can.Print(outstring.c_str());
}

void Plotter::MakeRatioPlot()
{
	TH1F *Wclone = (TH1F*)W->Clone("Wclone");
	TH1F *Topclone = (TH1F*)Top->Clone("Topclone");
	TH1F *DYclone = (TH1F*)DY->Clone("DYclone");
	TH1F *Dibclone = (TH1F*)Dib->Clone("Dibclone");
	TH1F *QCDclone = (TH1F*)QCD->Clone("QCDclone");
	TH1F *dataclone = (TH1F*)data->Clone("dataclone");
	
	Wclone->Add(Topclone);
	Wclone->Add(DYclone);
	Wclone->Add(Dibclone);
	Wclone->Add(QCDclone);
	dataclone->Divide(Wclone);
	ratio->Add(dataclone);
	
	if(plotVar.compare("mT")==0) ratio->SetTitle(";M_{T} [GeV];Data/MC");
	if(plotVar.compare("muon_pt")==0) ratio->SetTitle(";Muon P_{T} [GeV];Data/MC");
	if(plotVar.compare("muon_eta")==0) ratio->SetTitle(";Muon #eta;Data/MC");
	if(plotVar.compare("wmet_et")==0) ratio->SetTitle(";E_{T}^{miss} [GeV];Data/MC");
	if(plotVar.compare("abs(dPhi_metmuon)")==0) ratio->SetTitle(";#Delta#phi_{#mu,E_{T}^{miss}};Data/MC");
	if(plotVar.compare("cos(dPhi_metmuon)")==0) ratio->SetTitle(";cos[#Delta#phi_{#mu,E_{T}^{miss}}];Data/MC");

	
	ratio->SetLineColor(kBlack);
   	ratio->SetMinimum(0.25);  // Define Y ..
   	ratio->SetMaximum(1.75); // .. range
   	ratio->Sumw2();
   	ratio->SetStats(0);      // No statistics on lower plot
   	ratio->SetMarkerStyle(20);
	ratio->GetXaxis()->SetTitleSize(0.2);
	ratio->GetXaxis()->SetLabelSize(0.15);
	ratio->GetXaxis()->SetTitleOffset(0.95);
	ratio->GetYaxis()->SetTitleSize(0.15);
	ratio->GetYaxis()->SetLabelSize(0.15);
	ratio->GetYaxis()->SetTitleOffset(0.45);
	ratio->GetYaxis()->SetNdivisions(503);

	
	/*Double_t MC = 0.0;
	for(Int_t bin=1; bin <= nxbins; bin++)
	{
		MC += W->GetBinContent(bin);
		MC += Top->GetBinContent(bin);
		MC += DY->GetBinContent(bin);
		MC += Dib->GetBinContent(bin);
		MC += QCD->GetBinContent(bin);
		if(MC==0.0) continue;
		ratio->SetBinContent(bin,data->GetBinContent(bin)/MC);
		MC = 0.0;
	}*/
}

void Plotter::PlotWPrimeKinematics(Double_t eta_c)
{
	TCanvas can("can","can",800,600);
	can.SetLogy();
	can.cd();
	//Double_t pzbound =2000+(eta_c-1)*750;
	//TF1 *f1 = new TF1("f1","1",-pzbound,pzbound);
	Double_t pzbound =500;
	TF1 *f1 = new TF1("f1","1",0,pzbound);
	
	TCut minJetPTdR("Max$((vjet_pt>30)&&(vdR_jetmuon>2))");
	
	TPaveText *etaTxt = new TPaveText(.45,.725,.7,.9,"brNDC");
   	etaTxt->AddText(Form("#eta_{E_{T}^{miss}} = %.0f#times#eta_{#mu}",eta_c));
   	etaTxt->SetBorderSize(0);
	etaTxt->SetFillColor(0);
	
	vector<TH1F*> *pz = new vector<TH1F*>(6);
	MakeHistVector(pz,"pz",50,0,pzbound);
	
	vector<TH1F*> *dR = new vector<TH1F*>(6);
	MakeHistVector(dR,"dR",50,0,8);
	
	TLorentzVector v_mu, v_nu, v_wp;
	Double_t muon_pt, muon_eta, muon_phi, wmet, wmet_phi;
	Char_t isTight;
	vector<float>  	*vjet_pt = 0;
	vector<float>  	*vdR_jetmuon = 0;
	Int_t jet_i;
	
//----------------- Fill data plot ---------------------------
	dataTree->SetBranchAddress("muon_pt",&muon_pt);
	dataTree->SetBranchAddress("muon_eta",&muon_eta);
	dataTree->SetBranchAddress("muon_phi",&muon_phi);
	dataTree->SetBranchAddress("wmet_et",&wmet);
	dataTree->SetBranchAddress("wmet_phi",&wmet_phi);
	dataTree->SetBranchAddress("isTight",&isTight);
	dataTree->SetBranchAddress("vjet_pt",&vjet_pt);
	dataTree->SetBranchAddress("vdR_jetmuon",&vdR_jetmuon);

	Long64_t nentries = dataTree->GetEntries();
	for (Long64_t i=0;i<nentries;i++) 
	{
     	dataTree->GetEntry(i);
		jet_i = -1;
		for(Int_t j = 0; j < vjet_pt->size(); j++)
		{
			if(jet_i < 0) jet_i = j;
			else if(jet_i>=0)
			{
				if(vjet_pt->at(j) > vjet_pt->at(jet_i)) jet_i = j;
			}	
		}
		if(jet_i == -1) continue;
		dR->at(0)->Fill(vdR_jetmuon->at(jet_i),isTight);
		
		v_mu.SetPtEtaPhiE(muon_pt,muon_eta,muon_phi,muon_pt*cosh(muon_eta));
		v_nu.SetPtEtaPhiE(wmet,eta_c*muon_eta,wmet_phi,wmet*cosh(eta_c*muon_eta));
		v_wp = v_mu+v_nu;
		pz->at(0)->Fill(v_wp.Pt(),isTight);
  	} 
	
	Float_t mcCrossSection, mcFilterEfficiency, mceventweight, totalsumOfEventWeights;
	Double_t puweight, KfactorWeight;
	
//----------------- Fill W plot ---------------------------
	WTree->SetBranchAddress("muon_pt",&muon_pt);
	WTree->SetBranchAddress("muon_eta",&muon_eta);
	WTree->SetBranchAddress("muon_phi",&muon_phi);
	WTree->SetBranchAddress("wmet_et",&wmet);
	WTree->SetBranchAddress("wmet_phi",&wmet_phi);
	WTree->SetBranchAddress("isTight",&isTight);
	WTree->SetBranchAddress("mcCrossSection",&mcCrossSection);
	WTree->SetBranchAddress("mcFilterEfficiency",&mcFilterEfficiency);
	WTree->SetBranchAddress("mceventweight",&mceventweight);
	WTree->SetBranchAddress("puweight",&puweight);
	WTree->SetBranchAddress("KfactorWeight",&KfactorWeight);
	WTree->SetBranchAddress("totalsumOfEventWeights",&totalsumOfEventWeights);
	WTree->SetBranchAddress("vjet_pt",&vjet_pt);
	WTree->SetBranchAddress("vdR_jetmuon",&vdR_jetmuon);
	
	nentries = WTree->GetEntries();
	for (Long64_t i=0;i<nentries;i++) 
	{
     	WTree->GetEntry(i);
		jet_i = -1;
		for(Int_t j = 0; j < vjet_pt->size(); j++)
		{
			if(jet_i < 0) jet_i = j;
			else if(jet_i>=0)
			{
				if(vjet_pt->at(j) > vjet_pt->at(jet_i)) jet_i = j;
			}	
		}
		if(jet_i == -1) continue;
		dR->at(1)->Fill(vdR_jetmuon->at(jet_i),isTight*mcCrossSection*mcFilterEfficiency*mceventweight*puweight*KfactorWeight/totalsumOfEventWeights);
		
		v_mu.SetPtEtaPhiE(muon_pt,muon_eta,muon_phi,muon_pt*cosh(muon_eta));
		v_nu.SetPtEtaPhiE(wmet,eta_c*muon_eta,wmet_phi,wmet*cosh(eta_c*muon_eta));
		v_wp = v_mu+v_nu;
		pz->at(1)->Fill(v_wp.Pt(),isTight*mcCrossSection*mcFilterEfficiency*mceventweight*puweight*KfactorWeight/totalsumOfEventWeights);
  	} 

//----------------- Fill Top plot ---------------------------
	TopTree->SetBranchAddress("muon_pt",&muon_pt);
	TopTree->SetBranchAddress("muon_eta",&muon_eta);
	TopTree->SetBranchAddress("muon_phi",&muon_phi);
	TopTree->SetBranchAddress("wmet_et",&wmet);
	TopTree->SetBranchAddress("wmet_phi",&wmet_phi);
	TopTree->SetBranchAddress("isTight",&isTight);
	TopTree->SetBranchAddress("mcCrossSection",&mcCrossSection);
	TopTree->SetBranchAddress("mcFilterEfficiency",&mcFilterEfficiency);
	TopTree->SetBranchAddress("mceventweight",&mceventweight);
	TopTree->SetBranchAddress("puweight",&puweight);
	TopTree->SetBranchAddress("KfactorWeight",&KfactorWeight);
	TopTree->SetBranchAddress("totalsumOfEventWeights",&totalsumOfEventWeights);
	TopTree->SetBranchAddress("vjet_pt",&vjet_pt);
	TopTree->SetBranchAddress("vdR_jetmuon",&vdR_jetmuon);
	
	nentries = TopTree->GetEntries();
	for (Long64_t i=0;i<nentries;i++) 
	{
     	TopTree->GetEntry(i);
		jet_i = -1;
		for(Int_t j = 0; j < vjet_pt->size(); j++)
		{
			if(jet_i < 0) jet_i = j;
			else if(jet_i>=0)
			{
				if(vjet_pt->at(j) > vjet_pt->at(jet_i)) jet_i = j;
			}	
		}
		if(jet_i == -1) continue;
		dR->at(2)->Fill(vdR_jetmuon->at(jet_i),isTight*mcCrossSection*mcFilterEfficiency*mceventweight*puweight*KfactorWeight/totalsumOfEventWeights);
		
		v_mu.SetPtEtaPhiE(muon_pt,muon_eta,muon_phi,muon_pt*cosh(muon_eta));
		v_nu.SetPtEtaPhiE(wmet,eta_c*muon_eta,wmet_phi,wmet*cosh(eta_c*muon_eta));
		v_wp = v_mu+v_nu;
		pz->at(2)->Fill(v_wp.Pt(),isTight*mcCrossSection*mcFilterEfficiency*mceventweight*puweight*KfactorWeight/totalsumOfEventWeights);
  	} 
		
//----------------- Fill DY plot ---------------------------
	DYTree->SetBranchAddress("muon_pt",&muon_pt);
	DYTree->SetBranchAddress("muon_eta",&muon_eta);
	DYTree->SetBranchAddress("muon_phi",&muon_phi);
	DYTree->SetBranchAddress("wmet_et",&wmet);
	DYTree->SetBranchAddress("wmet_phi",&wmet_phi);
	DYTree->SetBranchAddress("isTight",&isTight);
	DYTree->SetBranchAddress("mcCrossSection",&mcCrossSection);
	DYTree->SetBranchAddress("mcFilterEfficiency",&mcFilterEfficiency);
	DYTree->SetBranchAddress("mceventweight",&mceventweight);
	DYTree->SetBranchAddress("puweight",&puweight);
	DYTree->SetBranchAddress("KfactorWeight",&KfactorWeight);
	DYTree->SetBranchAddress("totalsumOfEventWeights",&totalsumOfEventWeights);
	DYTree->SetBranchAddress("vjet_pt",&vjet_pt);
	DYTree->SetBranchAddress("vdR_jetmuon",&vdR_jetmuon);

	nentries = DYTree->GetEntries();
	for (Long64_t i=0;i<nentries;i++) 
	{
     	DYTree->GetEntry(i);
		jet_i = -1;
		for(Int_t j = 0; j < vjet_pt->size(); j++)
		{
			if(jet_i < 0) jet_i = j;
			else if(jet_i>=0)
			{
				if(vjet_pt->at(j) > vjet_pt->at(jet_i)) jet_i = j;
			}	
		}
		if(jet_i == -1) continue;
		dR->at(3)->Fill(vdR_jetmuon->at(jet_i),isTight*mcCrossSection*mcFilterEfficiency*mceventweight*puweight*KfactorWeight/totalsumOfEventWeights);
		
		v_mu.SetPtEtaPhiE(muon_pt,muon_eta,muon_phi,muon_pt*cosh(muon_eta));
		v_nu.SetPtEtaPhiE(wmet,eta_c*muon_eta,wmet_phi,wmet*cosh(eta_c*muon_eta));
		v_wp = v_mu+v_nu;
		pz->at(3)->Fill(v_wp.Pt(),isTight*mcCrossSection*mcFilterEfficiency*mceventweight*puweight*KfactorWeight/totalsumOfEventWeights);
  	} 
	
//----------------- Fill Dib plot ---------------------------
	DibTree->SetBranchAddress("muon_pt",&muon_pt);
	DibTree->SetBranchAddress("muon_eta",&muon_eta);
	DibTree->SetBranchAddress("muon_phi",&muon_phi);
	DibTree->SetBranchAddress("wmet_et",&wmet);
	DibTree->SetBranchAddress("wmet_phi",&wmet_phi);
	DibTree->SetBranchAddress("isTight",&isTight);
	DibTree->SetBranchAddress("mcCrossSection",&mcCrossSection);
	DibTree->SetBranchAddress("mcFilterEfficiency",&mcFilterEfficiency);
	DibTree->SetBranchAddress("mceventweight",&mceventweight);
	DibTree->SetBranchAddress("puweight",&puweight);
	DibTree->SetBranchAddress("KfactorWeight",&KfactorWeight);
	DibTree->SetBranchAddress("totalsumOfEventWeights",&totalsumOfEventWeights);
	DibTree->SetBranchAddress("vjet_pt",&vjet_pt);
	DibTree->SetBranchAddress("vdR_jetmuon",&vdR_jetmuon);

	nentries = DibTree->GetEntries();
	for (Long64_t i=0;i<nentries;i++) 
	{
     	DibTree->GetEntry(i);
		jet_i = -1;
		for(Int_t j = 0; j < vjet_pt->size(); j++)
		{
			if( jet_i < 0) jet_i = j;
			else if(jet_i>=0)
			{
				if(vjet_pt->at(j) > vjet_pt->at(jet_i)) jet_i = j;
			}	
		}
		if(jet_i == -1) continue;
		dR->at(4)->Fill(vdR_jetmuon->at(jet_i),isTight*mcCrossSection*mcFilterEfficiency*mceventweight*puweight*KfactorWeight/totalsumOfEventWeights);
		
		v_mu.SetPtEtaPhiE(muon_pt,muon_eta,muon_phi,muon_pt*cosh(muon_eta));
		v_nu.SetPtEtaPhiE(wmet,eta_c*muon_eta,wmet_phi,wmet*cosh(eta_c*muon_eta));
		v_wp = v_mu+v_nu;
		pz->at(4)->Fill(v_wp.Pt(),isTight*mcCrossSection*mcFilterEfficiency*mceventweight*puweight*KfactorWeight/totalsumOfEventWeights);
  	} 
	
//------------------- Fill QCD plot ----------------------
	Double_t e;
	
	nentries = dataTree->GetEntries();
	for (Long64_t i=0;i<nentries;i++) 
	{
     	dataTree->GetEntry(i);
		if(fakeEffParamName.compare("muon_pt")==0) e = eff->GetBinContent(eff->GetBin(muon_pt));
		if(fakeEffParamName.compare("muon_eta")==0) e = eff->GetBinContent(eff->GetBin(muon_eta));
		jet_i = -1;
		for(Int_t j = 0; j < vjet_pt->size(); j++)
		{
			if(jet_i < 0) jet_i = j;
			else if(jet_i>=0)
			{
				if(vjet_pt->at(j) > vjet_pt->at(jet_i)) jet_i = j;
			}	
		}
		if(jet_i == -1) continue;
		dR->at(5)->Fill(vdR_jetmuon->at(jet_i),e/(1-e)*(1-isTight));
		
		v_mu.SetPtEtaPhiE(muon_pt,muon_eta,muon_phi,muon_pt*cosh(muon_eta));
		v_nu.SetPtEtaPhiE(wmet,eta_c*muon_eta,wmet_phi,wmet*cosh(eta_c*muon_eta));
		v_wp = v_mu+v_nu;
		pz->at(5)->Fill(v_wp.Pt(),e/(1-e)*(1-isTight));
  	} 
	
//-------- scale MC plots -----------
	for(Int_t i = 1; i<5; i++)
	{
		pz->at(i)->Multiply(f1,Lumi*1e6);
		dR->at(i)->Multiply(f1,Lumi*1e6);
	}
	pz->at(0)->SetTitle(";W' Object P_{T} [GeV]; Events");
	dR->at(0)->SetTitle(";dR(jet,#mu); Events");
	
	THStack *stack = new THStack("stack","");
	stack->Add(pz->at(5),"hist");
	stack->Add(pz->at(4),"hist");
	stack->Add(pz->at(3),"hist");
	stack->Add(pz->at(2),"hist");
	stack->Add(pz->at(1),"hist");
	
	THStack *stack2 = new THStack("stack2","");
	stack2->Add(dR->at(5),"hist");
	stack2->Add(dR->at(4),"hist");
	stack2->Add(dR->at(3),"hist");
	stack2->Add(dR->at(2),"hist");
	stack2->Add(dR->at(1),"hist");
	
	TLegend *leg = new TLegend(0.75,0.6,0.92,0.9);
	leg->AddEntry(pz->at(0),"Data","pe");
	leg->AddEntry(pz->at(1),"W","f");
	leg->AddEntry(pz->at(2),"Top","f");
	leg->AddEntry(pz->at(3),"Z/#gamma*","f");
	leg->AddEntry(pz->at(4),"Diboson","f");
	leg->AddEntry(pz->at(5),"QCD","f");
	leg->SetTextFont(42);
	
	pz->at(0)->SetMinimum(1);
	pz->at(0)->SetMaximum(5e6);
	pz->at(0)->Draw("p");
	stack->Draw("same");
	etaTxt->Draw();
	leg->Draw();
	can.Print(Form("/export/home/bbullard/thesis/plots/%.1f/WPrimePT-nuEta_%.0f.png",Lumi,eta_c));
	
	dR->at(0)->SetMinimum(1);
	dR->at(0)->SetMaximum(5e8);
	dR->at(0)->Draw("p");
	stack2->Draw("same");
	etaTxt->Draw();
	leg->Draw();
	//can.Print(Form("/export/home/bbullard/thesis/plots/%.1f/dRjetmuon-nuEta_%.0f_mindR.png",Lumi,eta_c));
	
	for(Int_t i = 0; i<6; i++)
	{
		delete pz->at(i);
		delete dR->at(i);
	}
	delete pz;
	delete dR;
}

void Plotter::OutputBinContents(TH1F *hist, string mod = "")
{
	Int_t N = hist->GetNbinsX();
	saveAs = Form("/export/home/bbullard/thesis/plots/%.1f/%s%s.csv",Lumi,hist->GetName(),mod.c_str());
	
	ofstream binsF;
  	binsF.open(saveAs, std::ofstream::out);
	for(Int_t b = 1; b <= N; b++)
  		binsF<<b<<","<<hist->GetBinLowEdge(b)<<"-"<<hist->GetBinLowEdge(b)+hist->GetBinWidth(b)<<","<<hist->GetBinContent(b)<<endl;

  	binsF.close();
}

void Plotter::PlotSignalmT()
{
	TCanvas can("can","can",1200,800);
	
	W2000->SetLineColor(kAzure+10);
	W3000->SetLineColor(kAzure);
	W4000->SetLineColor(kViolet);
	W5000->SetLineColor(kPink);
	
	W2000->SetLineWidth(2);
	W3000->SetLineWidth(2);
	W4000->SetLineWidth(2);
	W5000->SetLineWidth(2);
	
	TLegend *sleg = new TLegend(0.75,0.6,0.92,0.9);
	sleg->AddEntry(W2000,"2 TeV","l");
	sleg->AddEntry(W3000,"3 TeV","l");
	sleg->AddEntry(W4000,"4 TeV","l");
	sleg->AddEntry(W5000,"5 TeV","l");
	sleg->SetTextFont(42);
	
	W2000->SetTitle(";mT [TeV];Entries");
	W2000->Draw();
	W3000->Draw("same");
	W4000->Draw("same");
	W5000->Draw("same");
	sleg->Draw();
	
	can.Print(Form("/export/home/bbullard/thesis/plots/%.1f/signalMC.png",Lumi));
}

void Plotter::MakeHistVector(vector<TH1F*> *histVector, string name, Int_t bins, Double_t xmin, Double_t xmax)
{
	string title;
	
	title = name+"data";
	histVector->at(0) = new TH1F(title.c_str(),"",bins,xmin,xmax);
	title = name+"W";
	histVector->at(1) = new TH1F(title.c_str(),"",bins,xmin,xmax);
	title = name+"Top";
	histVector->at(2) = new TH1F(title.c_str(),"",bins,xmin,xmax);
	title = name+"DY";
	histVector->at(3) = new TH1F(title.c_str(),"",bins,xmin,xmax);
	title = name+"Dib";
	histVector->at(4) = new TH1F(title.c_str(),"",bins,xmin,xmax);
	title = name+"QCD";
	histVector->at(5) = new TH1F(title.c_str(),"",bins,xmin,xmax);
	
	histVector->at(1)->SetLineWidth(2);
	histVector->at(2)->SetLineWidth(2);
	histVector->at(3)->SetLineWidth(2);
	histVector->at(4)->SetLineWidth(2);
	histVector->at(5)->SetLineWidth(2);
	
	histVector->at(1)->SetFillColor(kWhite);
	histVector->at(2)->SetFillColor(kRed+1);
	histVector->at(3)->SetFillColor(kAzure-9);
	histVector->at(4)->SetFillColor(kOrange);
	histVector->at(5)->SetFillColor(kYellow-9);
}

void Plotter::MakeRatioBounds(TH1F *ratioLB, TH1F *ratioUB)
{
	Double_t LumiSFup = 1+sqrt(pow(0.032,2)+pow(0.15,2));
	Double_t LumiSFdown = 1-sqrt(pow(0.032,2)+pow(0.15,2));
	
	TFile *f;	
	TF1 *f1 = new TF1("f1","1",xbins[0],xbins[nxbins-1]);
	string input;
	TDirectory *dir1 = gDirectory;   
	
	f = TFile::Open("Files/ntuples/selection/WFullSelection_UB.root","read");
	dir1->cd();
	TTree *WUBTree = (TTree*)f->Get("tree");
	f = TFile::Open("Files/ntuples/selection/TopFullSelection_UB.root","read");
	dir1->cd();
	TTree *TopUBTree = (TTree*)f->Get("tree");
	f = TFile::Open("Files/ntuples/selection/DYFullSelection_UB.root","read");
	dir1->cd();
	TTree *DYUBTree = (TTree*)f->Get("tree");
	f = TFile::Open("Files/ntuples/selection/DibFullSelection_UB.root","read");
	dir1->cd();
	TTree *DibUBTree = (TTree*)f->Get("tree");
	
	TH1F *WUB = new TH1F("WUB","",nxbins,xbins);
	TH1F *TopUB = new TH1F("TopUB","",nxbins,xbins);
	TH1F *DYUB = new TH1F("DYUB","",nxbins,xbins);
	TH1F *DibUB = new TH1F("DibUB","",nxbins,xbins);
	
	input=plotVar+">>WUB";
	WUBTree->Draw(input.c_str(),"SF_trigmu50_highptID*muon_HighPtSF*mcCrossSection*mcFilterEfficiency*mceventweight*puweight*KfactorWeight*isTight/totalsumOfEventWeights");
	
	input=plotVar+">>TopUB";	
	TopUBTree->Draw(input.c_str(),"SF_trigmu50_highptID*muon_HighPtSF*mcCrossSection*mcFilterEfficiency*mceventweight*puweight*KfactorWeight*isTight/totalsumOfEventWeights");

	input=plotVar+">>DYUB";
	DYUBTree->Draw(input.c_str(),"SF_trigmu50_highptID*muon_HighPtSF*mcCrossSection*mcFilterEfficiency*mceventweight*puweight*KfactorWeight*isTight/totalsumOfEventWeights");
	
	input=plotVar+">>DibUB";
	DibUBTree->Draw(input.c_str(),"SF_trigmu50_highptID*muon_HighPtSF*mcCrossSection*mcFilterEfficiency*mceventweight*puweight*KfactorWeight*isTight/totalsumOfEventWeights");
	
	
	f = TFile::Open("Files/ntuples/selection/WFullSelection_LB.root","read");
	dir1->cd();
	TTree *WLBTree = (TTree*)f->Get("tree");
	f = TFile::Open("Files/ntuples/selection/TopFullSelection_LB.root","read");
	dir1->cd();
	TTree *TopLBTree = (TTree*)f->Get("tree");
	f = TFile::Open("Files/ntuples/selection/DYFullSelection_LB.root","read");
	dir1->cd();
	TTree *DYLBTree = (TTree*)f->Get("tree");
	f = TFile::Open("Files/ntuples/selection/DibFullSelection_LB.root","read");
	dir1->cd();
	TTree *DibLBTree = (TTree*)f->Get("tree");
	
	TH1F *WLB = new TH1F("WLB","",nxbins,xbins);
	TH1F *TopLB = new TH1F("TopLB","",nxbins,xbins);
	TH1F *DYLB = new TH1F("DYLB","",nxbins,xbins);
	TH1F *DibLB = new TH1F("DibLB","",nxbins,xbins);
	
	input=plotVar+">>WLB";
	WLBTree->Draw(input.c_str(),"SF_trigmu50_highptID*muon_HighPtSF*mcCrossSection*mcFilterEfficiency*mceventweight*puweight*KfactorWeight*isTight/totalsumOfEventWeights");
	
	input=plotVar+">>TopLB";	
	TopLBTree->Draw(input.c_str(),"SF_trigmu50_highptID*muon_HighPtSF*mcCrossSection*mcFilterEfficiency*mceventweight*puweight*KfactorWeight*isTight/totalsumOfEventWeights");

	input=plotVar+">>DYLB";
	DYLBTree->Draw(input.c_str(),"SF_trigmu50_highptID*muon_HighPtSF*mcCrossSection*mcFilterEfficiency*mceventweight*puweight*KfactorWeight*isTight/totalsumOfEventWeights");
	
	input=plotVar+">>DibLB";
	DibLBTree->Draw(input.c_str(),"SF_trigmu50_highptID*muon_HighPtSF*mcCrossSection*mcFilterEfficiency*mceventweight*puweight*KfactorWeight*isTight/totalsumOfEventWeights");
	
	TCanvas can("can","can",800,600);
	can.SetLogy();
	if(useLogX) can.SetLogx();
	
	WUB->Multiply(f1,Lumi*1e6);
	TopUB->Multiply(f1,Lumi*1e6);
	DYUB->Multiply(f1,Lumi*1e6);
	DibUB->Multiply(f1,Lumi*1e6);
	
	WLB->Multiply(f1,Lumi*1e6);
	TopLB->Multiply(f1,Lumi*1e6);
	DYLB->Multiply(f1,Lumi*1e6);
	DibLB->Multiply(f1,Lumi*1e6);
	
	WLB->SetLineColor(1);
	WUB->SetLineColor(1);
	WLB->SetLineStyle(3);
	WUB->SetLineStyle(2);
	TopLB->SetLineColor(kRed+1);
	TopUB->SetLineColor(kRed+1);
	TopLB->SetLineStyle(3);
	TopUB->SetLineStyle(2);
	DYLB->SetLineColor(kAzure-9);
	DYUB->SetLineColor(kAzure-9);
	DYLB->SetLineStyle(3);
	DYUB->SetLineStyle(2);
	DibLB->SetLineColor(kOrange);
	DibUB->SetLineColor(kOrange);
	DibLB->SetLineStyle(3);
	DibUB->SetLineStyle(2);
	
	THStack *stack = new THStack("stack","");
	//stack->Add(QCD,"hist");
	//stack->Add(DibUB,"hist");
	//stack->Add(DYUB,"hist");
	stack->Add(TopUB,"hist");
	stack->Add(WUB,"hist");
	
	THStack *stack2 = new THStack("stack2","");
	//stack->Add(QCD,"hist");
	//stack->Add(DibLB,"hist");
	//stack->Add(DYLB,"hist");
	stack->Add(TopLB,"hist");
	stack->Add(WLB,"hist");
	
	stack->SetMaximum(5*1e6);
	stack->SetMinimum(5e-4);
	//stack->Draw();
	//stack2->Draw("same");
	//dataLB->Draw("same pe");
	//dataUB->Draw("same pe");
	//can.Print("mT_ranges_beforeLumi.png");
	
	Double_t BC;
	for(Int_t bin=1; bin <= nxbins; bin++)
	{
		if(WUB->GetBinContent(bin)<WLB->GetBinContent(bin))
		{
			cout<<"Switching W bin content for bin "<<bin<<endl;
			BC = WUB->GetBinContent(bin);
			WUB->SetBinContent(bin,WLB->GetBinContent(bin));
			WLB->SetBinContent(bin,BC);
		}
		if(TopUB->GetBinContent(bin)<TopLB->GetBinContent(bin))
		{
			cout<<"Switching Top bin content for bin "<<bin<<endl;
			BC = TopUB->GetBinContent(bin);
			TopUB->SetBinContent(bin,TopLB->GetBinContent(bin));
			TopLB->SetBinContent(bin,BC);
		}
		if(DYUB->GetBinContent(bin)<DYLB->GetBinContent(bin))
		{
			cout<<"Switching DY bin content for bin "<<bin<<endl;
			BC = DYUB->GetBinContent(bin);
			DYUB->SetBinContent(bin,DYLB->GetBinContent(bin));
			DYLB->SetBinContent(bin,BC);
		}
		if(DibUB->GetBinContent(bin)<DibLB->GetBinContent(bin))
		{
			cout<<"Switching Dib bin content for bin "<<bin<<endl;
			BC = DibUB->GetBinContent(bin);
			DibUB->SetBinContent(bin,DibLB->GetBinContent(bin));
			DibLB->SetBinContent(bin,BC);
		}
	}
	
	WUB->Multiply(f1,LumiSFup);
	TopUB->Multiply(f1,LumiSFup);
	DYUB->Multiply(f1,LumiSFup);
	DibUB->Multiply(f1,LumiSFup);
	
	WLB->Multiply(f1,LumiSFdown);
	TopLB->Multiply(f1,LumiSFdown);
	DYLB->Multiply(f1,LumiSFdown);
	DibLB->Multiply(f1,LumiSFdown);
	
	Double_t MC = 0.0;
	for(Int_t bin=1; bin <= nxbins; bin++)
	{
		MC += WUB->GetBinContent(bin);
		MC += TopUB->GetBinContent(bin);
		MC += DYUB->GetBinContent(bin);
		MC += DibUB->GetBinContent(bin);
		MC += QCD->GetBinContent(bin);
		if(MC==0.0) continue;
		ratioLB->SetBinContent(bin,data->GetBinContent(bin)/MC);
		MC = 0.0;
	}
	
	MC = 0.0;
	for(Int_t bin=1; bin <= nxbins; bin++)
	{
		MC += WLB->GetBinContent(bin);
		MC += TopLB->GetBinContent(bin);
		MC += DYLB->GetBinContent(bin);
		MC += DibLB->GetBinContent(bin);
		MC += QCD->GetBinContent(bin);
		if(MC==0.0) continue;
		ratioUB->SetBinContent(bin,data->GetBinContent(bin)/MC);
		MC = 0.0;
	}
	
	OutputBinContents(WLB,plotVar);
	OutputBinContents(TopLB,plotVar);
	OutputBinContents(DYLB,plotVar);
	OutputBinContents(DibLB,plotVar);
	
	OutputBinContents(W,plotVar);
	OutputBinContents(Top,plotVar);
	OutputBinContents(DY,plotVar);
	OutputBinContents(Dib,plotVar);
	OutputBinContents(QCD,plotVar);
	OutputBinContents(data,plotVar);
	
	OutputBinContents(WUB,plotVar);
	OutputBinContents(TopUB,plotVar);
	OutputBinContents(DYUB,plotVar);
	OutputBinContents(DibUB,plotVar);

	
	f = TFile::Open("Files/ntuples/selection/signal2FullSelection_UB.root","read");
	dir1->cd();
	TTree *W2000UBTree = (TTree*)f->Get("tree");
	f = TFile::Open("Files/ntuples/selection/signal3FullSelection_UB.root","read");
	dir1->cd();
	TTree *W3000UBTree = (TTree*)f->Get("tree");
	f = TFile::Open("Files/ntuples/selection/signal4FullSelection_UB.root","read");
	dir1->cd();
	TTree *W4000UBTree = (TTree*)f->Get("tree");
	f = TFile::Open("Files/ntuples/selection/signal4FullSelection_UB.root","read");
	dir1->cd();
	TTree *W5000UBTree = (TTree*)f->Get("tree");
	
	TH1F *W2000UB = new TH1F("W2000UB","",nxbins,xbins);
	TH1F *W3000UB = new TH1F("W3000UB","",nxbins,xbins);
	TH1F *W4000UB = new TH1F("W4000UB","",nxbins,xbins);
	TH1F *W5000UB = new TH1F("W5000UB","",nxbins,xbins);
	
	input=plotVar+">>W2000UB";
	W2000UBTree->Draw(input.c_str(),"SF_trigmu50_highptID*muon_HighPtSF*mcCrossSection*mcFilterEfficiency*mceventweight*puweight*KfactorWeight*isTight/totalsumOfEventWeights");
	W2000UB->Multiply(f1,Lumi*LumiSFup*1e6);
	
	input=plotVar+">>W3000UB";	
	W3000UBTree->Draw(input.c_str(),"SF_trigmu50_highptID*muon_HighPtSF*mcCrossSection*mcFilterEfficiency*mceventweight*puweight*KfactorWeight*isTight/totalsumOfEventWeights");
	W3000UB->Multiply(f1,Lumi*LumiSFup*1e6);

	input=plotVar+">>W4000UB";
	W4000UBTree->Draw(input.c_str(),"SF_trigmu50_highptID*muon_HighPtSF*mcCrossSection*mcFilterEfficiency*mceventweight*puweight*KfactorWeight*isTight/totalsumOfEventWeights");
	W2000UB->Multiply(f1,Lumi*LumiSFup*1e6);
	
	input=plotVar+">>W5000UB";
	W5000UBTree->Draw(input.c_str(),"SF_trigmu50_highptID*muon_HighPtSF*mcCrossSection*mcFilterEfficiency*mceventweight*puweight*KfactorWeight*isTight/totalsumOfEventWeights");
	W5000UB->Multiply(f1,Lumi*LumiSFup*1e6);
	
	
	f = TFile::Open("Files/ntuples/selection/signal2FullSelection_LB.root","read");
	dir1->cd();
	TTree *W2000LBTree = (TTree*)f->Get("tree");
	f = TFile::Open("Files/ntuples/selection/signal3FullSelection_LB.root","read");
	dir1->cd();
	TTree *W3000LBTree = (TTree*)f->Get("tree");
	f = TFile::Open("Files/ntuples/selection/signal4FullSelection_LB.root","read");
	dir1->cd();
	TTree *W4000LBTree = (TTree*)f->Get("tree");
	f = TFile::Open("Files/ntuples/selection/signal4FullSelection_LB.root","read");
	dir1->cd();
	TTree *W5000LBTree = (TTree*)f->Get("tree");
	
	TH1F *W2000LB = new TH1F("W2000LB","",nxbins,xbins);
	TH1F *W3000LB = new TH1F("W3000LB","",nxbins,xbins);
	TH1F *W4000LB = new TH1F("W4000LB","",nxbins,xbins);
	TH1F *W5000LB = new TH1F("W5000LB","",nxbins,xbins);
	
	input=plotVar+">>W2000LB";
	W2000LBTree->Draw(input.c_str(),"SF_trigmu50_highptID*muon_HighPtSF*mcCrossSection*mcFilterEfficiency*mceventweight*puweight*KfactorWeight*isTight/totalsumOfEventWeights");
	W2000LB->Multiply(f1,Lumi*LumiSFdown*1e6);
	
	input=plotVar+">>W3000LB";	
	W3000LBTree->Draw(input.c_str(),"SF_trigmu50_highptID*muon_HighPtSF*mcCrossSection*mcFilterEfficiency*mceventweight*puweight*KfactorWeight*isTight/totalsumOfEventWeights");
	W3000LB->Multiply(f1,Lumi*LumiSFdown*1e6);

	input=plotVar+">>W4000LB";
	W4000LBTree->Draw(input.c_str(),"SF_trigmu50_highptID*muon_HighPtSF*mcCrossSection*mcFilterEfficiency*mceventweight*puweight*KfactorWeight*isTight/totalsumOfEventWeights");
	W2000LB->Multiply(f1,Lumi*LumiSFdown*1e6);
	
	input=plotVar+">>W5000LB";
	W5000LBTree->Draw(input.c_str(),"SF_trigmu50_highptID*muon_HighPtSF*mcCrossSection*mcFilterEfficiency*mceventweight*puweight*KfactorWeight*isTight/totalsumOfEventWeights");
	W5000LB->Multiply(f1,Lumi*LumiSFdown*1e6);
	
	OutputBinContents(W2000UB,plotVar);
	OutputBinContents(W3000UB,plotVar);
	OutputBinContents(W4000UB,plotVar);
	OutputBinContents(W5000UB,plotVar);
	
	OutputBinContents(W2000,plotVar);
	OutputBinContents(W3000,plotVar);
	OutputBinContents(W4000,plotVar);
	OutputBinContents(W5000,plotVar);
	
	OutputBinContents(W2000LB,plotVar);
	OutputBinContents(W3000LB,plotVar);
	OutputBinContents(W4000LB,plotVar);
	OutputBinContents(W5000LB,plotVar);
	
	W5000UB->Draw();
	W5000LB->Draw("same");
	W4000->Draw("same");
	//W3000LB->Draw("same");
	//W3000->Draw("same");
	//W3000UB->Draw("same");

	can.Print(Form("w4000Test_%s.png",plotVar.c_str()));
	
}
