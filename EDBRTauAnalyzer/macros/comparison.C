{
  gROOT->Reset();
  gStyle->SetCanvasColor(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1,0);
  gStyle->SetTitleX(0.5); //title X location 
  gStyle->SetTitleY(0.96); //title Y location 
  gStyle->SetPaintTextFormat(".2f");

  using namespace std;

  vector<TString> name;                                 vector<int> RebinFactor; 
  name.push_back("metPt");                              RebinFactor.push_back(10);
  name.push_back("metPtSelected");			RebinFactor.push_back(10);
  name.push_back("NJetSelected");			RebinFactor.push_back(1);
  name.push_back("NEleSelected");			RebinFactor.push_back(1);
  name.push_back("NMuoSelected");                       RebinFactor.push_back(1);
  name.push_back("NbtagsLSelected");			RebinFactor.push_back(1);
  name.push_back("NbtagsMSelected");			RebinFactor.push_back(1);
  name.push_back("NbtagsTSelected");			RebinFactor.push_back(1);
  name.push_back("NMuoSelected");                       RebinFactor.push_back(1);
  name.push_back("JetZdRSelected");			RebinFactor.push_back(1);
  name.push_back("EleMuoDRSelected");			RebinFactor.push_back(4);
  name.push_back("tauDecay");				RebinFactor.push_back(1);
  name.push_back("tauDecaySelected");			RebinFactor.push_back(1);
  name.push_back("MassVisEleEle");			RebinFactor.push_back(10); //13
  name.push_back("MassEffEleEle");			RebinFactor.push_back(10);
  name.push_back("MassSvfitEleEle");                    RebinFactor.push_back(10);
  name.push_back("MassCAEleEle");			RebinFactor.push_back(10);
  name.push_back("MassVisMuoMuo");			RebinFactor.push_back(10);
  name.push_back("MassSvfitMuoMuo");			RebinFactor.push_back(10);
  name.push_back("MassEffMuoMuo");			RebinFactor.push_back(10);
  name.push_back("MassCAMuoMuo");			RebinFactor.push_back(10);
  name.push_back("MassVisEleMuo");			RebinFactor.push_back(10);
  name.push_back("MassEffEleMuo");			RebinFactor.push_back(10);
  name.push_back("MassSvfitEleMuo");                    RebinFactor.push_back(10);
  name.push_back("MassCAEleMuo");			RebinFactor.push_back(10); //24
  name.push_back("XMassVisEleEle");			RebinFactor.push_back(50); //25
  name.push_back("XMassVisMuoMuo");			RebinFactor.push_back(50);
  name.push_back("XMassVisEleMuo");			RebinFactor.push_back(50);
  name.push_back("XMassEffEleEle");			RebinFactor.push_back(50);
  name.push_back("XMassEffMuoMuo");			RebinFactor.push_back(50);
  name.push_back("XMassEffEleMuo");			RebinFactor.push_back(50);
  name.push_back("XMassSVFitEleEle");			RebinFactor.push_back(50);
  name.push_back("XMassSVFitMuoMuo");			RebinFactor.push_back(50);
  name.push_back("XMassSVFitEleMuo");			RebinFactor.push_back(50);
  name.push_back("XMassCAEleEle");			RebinFactor.push_back(50);
  name.push_back("XMassCAMuoMuo");			RebinFactor.push_back(50);
  name.push_back("XMassCAEleMuo");			RebinFactor.push_back(50); //36
  name.push_back("jetPtSelected");			RebinFactor.push_back(10);
  name.push_back("jetEtaSelected");			RebinFactor.push_back(10);
  name.push_back("jetMassSelected");			RebinFactor.push_back(1);
  name.push_back("jetSubjettinessSelected");		RebinFactor.push_back(1);
  name.push_back("jetDRLepSelected");			RebinFactor.push_back(1);
  name.push_back("jetDRMetSelected");			RebinFactor.push_back(1);
  name.push_back("jetDRVisZSelected");			RebinFactor.push_back(1);
  name.push_back("jetDPhiMetSelected");			RebinFactor.push_back(1);
  name.push_back("muonDR");				RebinFactor.push_back(10);
  name.push_back("muonPtLead");				RebinFactor.push_back(10);
  name.push_back("muonEtaLead");			RebinFactor.push_back(10);
  name.push_back("muonChi2Lead");			RebinFactor.push_back(1);
  name.push_back("muonValidMuonHitsLead");		RebinFactor.push_back(1);
  name.push_back("muonMatchesLead");                    RebinFactor.push_back(1);
  name.push_back("muondBLead");				RebinFactor.push_back(1);
  name.push_back("muondZLead");				RebinFactor.push_back(1);
  name.push_back("muonPixelHitsLead");			RebinFactor.push_back(1);
  name.push_back("muonLayersLead");			RebinFactor.push_back(1);
  name.push_back("muonDRMetLead");			RebinFactor.push_back(1);
  name.push_back("muonDPhiMetLead");			RebinFactor.push_back(1);
  name.push_back("muonPtSublead");			RebinFactor.push_back(10);
  name.push_back("muonEtaSublead");			RebinFactor.push_back(10);
  name.push_back("muonChi2Sublead");			RebinFactor.push_back(1);
  name.push_back("muonValidMuonHitsSublead");		RebinFactor.push_back(1);
  name.push_back("muonMatchesSublead");			RebinFactor.push_back(1);
  name.push_back("muondBSublead");			RebinFactor.push_back(1);
  name.push_back("muondZSublead");			RebinFactor.push_back(1);
  name.push_back("muonPixelHitsSublead");		RebinFactor.push_back(1);
  name.push_back("muonLayersSublead");			RebinFactor.push_back(1);
  name.push_back("muonDRMetSublead");			RebinFactor.push_back(1);
  name.push_back("muonDPhiMetSublead");			RebinFactor.push_back(1);
  name.push_back("muonPtEMU");				RebinFactor.push_back(10);
  name.push_back("muonEtaEMU");				RebinFactor.push_back(10);
  name.push_back("muonChi2EMU");			RebinFactor.push_back(1);
  name.push_back("muonValidMuonHitsEMU");               RebinFactor.push_back(1);
  name.push_back("muonMatchesEMU");			RebinFactor.push_back(1); //70
  name.push_back("muondBEMU");				RebinFactor.push_back(1);
  name.push_back("muondZEMU");				RebinFactor.push_back(1);
  name.push_back("muonPixelHitsEMU");			RebinFactor.push_back(1);
  name.push_back("muonLayersEMU");			RebinFactor.push_back(1);
  name.push_back("muonDRMetEMU");			RebinFactor.push_back(1);
  name.push_back("muonDPhiMetEMU");			RebinFactor.push_back(1);
  name.push_back("electronDR");				RebinFactor.push_back(4);
  name.push_back("electronPtLead");			RebinFactor.push_back(10);
  name.push_back("electronEtaLead");			RebinFactor.push_back(10);
  name.push_back("electronDRMetLead");			RebinFactor.push_back(1);
  name.push_back("electronDPhiMetLead");		RebinFactor.push_back(1);
  name.push_back("electronPtSublead");			RebinFactor.push_back(10);
  name.push_back("electronEtaSublead");			RebinFactor.push_back(10);
  name.push_back("electronDRMetSublead");		RebinFactor.push_back(1);
  name.push_back("electronDPhiMetSublead");		RebinFactor.push_back(1);
  name.push_back("electronPtEMU");			RebinFactor.push_back(10);
  name.push_back("electronEtaEMU");			RebinFactor.push_back(10);
  name.push_back("electronDRMetEMU");			RebinFactor.push_back(1);
  name.push_back("electronDPhiMetEMU");			RebinFactor.push_back(1);
  name.push_back("electronDetIsoEMU");			RebinFactor.push_back(10); //90
  name.push_back("electronPFIsoEMU");			RebinFactor.push_back(10);
  name.push_back("electronDetIsoSublead");		RebinFactor.push_back(10);
  name.push_back("electronDetIsoLead");			RebinFactor.push_back(10);
  name.push_back("electronPFIsoSublead");		RebinFactor.push_back(10);
  name.push_back("electronPFIsoLead");			RebinFactor.push_back(10);
  name.push_back("muonDetIsoEMU");			RebinFactor.push_back(10);
  name.push_back("muonPFIsoEMU");			RebinFactor.push_back(10);
  name.push_back("muonPFIsoSublead");			RebinFactor.push_back(10);
  name.push_back("muonDetIsoLead");			RebinFactor.push_back(10);
  name.push_back("muonDetIsoSublead");			RebinFactor.push_back(10);
  name.push_back("muonPFIsoLead");			RebinFactor.push_back(10);


  int rebin = 10;
  
  TFile *f11= new TFile("../ROOT/BulkGM1000_v9.root");
  TFile *f12= new TFile("../ROOT/BulkGM1500_v9.root");
  TFile *f13= new TFile("../ROOT/BulkGM2000_v9.root");
  TFile *f2 = new TFile("../ROOT/DYJetsToLL_PtZ-100_v9.root");
  TFile *f3 = new TFile("../ROOT/DYJetsToLL_PtZ-70To100_v9.root");
  TFile *f4 = new TFile("../ROOT/TT_v9.root");
  TFile *f5 = new TFile("../ROOT/WW_v9.root");
  TFile *f6 = new TFile("../ROOT/WZ_v9.root");
  TFile *f7 = new TFile("../ROOT/ZZ_v9.root");
  TFile *f8 = new TFile("../ROOT/WJetsPt50To70_v8.root");
  TFile *f9 = new TFile("../ROOT/WJetsPt70To100_v8.root");
  TFile *f0 = new TFile("../ROOT/WJetsPt100_v8.root");
  
  for(int i=0; i<name.size(); i++){
    TCanvas * c1 = new TCanvas("c1", "c1", 800, 600);
    
    rebin=RebinFactor[i];
    float ScaleSignal = 0.1;
    
    TH1F *BulkGM1000;              BulkGM1000 = (TH1F*) f11->Get("demo/"+name[i]);
    TH1F *BulkGM1500;              BulkGM1500 = (TH1F*) f12->Get("demo/"+name[i]);
    TH1F *BulkGM2000;              BulkGM2000 = (TH1F*) f13->Get("demo/"+name[i]);
    TH1F *DYJetsToLL_PtZ100;       DYJetsToLL_PtZ100 = (TH1F*) f2->Get("demo/"+name[i]);
    TH1F *DYJetsToLL_PtZ70To100;   DYJetsToLL_PtZ70To100 = (TH1F*) f3->Get("demo/"+name[i]);
    TH1F *TT;                      TT = (TH1F*) f4->Get("demo/"+name[i]);
    TH1F *WW;                      WW = (TH1F*) f5->Get("demo/"+name[i]);
    TH1F *WZ;                      WZ = (TH1F*) f6->Get("demo/"+name[i]);
    TH1F *ZZ;                      ZZ = (TH1F*) f7->Get("demo/"+name[i]);
    TH1F *WJetsPt50To70;           WJetsPt50To70 = (TH1F*) f8->Get("demo/"+name[i]);
    TH1F *WJetsPt70To100;          WJetsPt70To100 = (TH1F*) f9->Get("demo/"+name[i]);
    TH1F *WJetsPt100;              WJetsPt100 = (TH1F*) f0->Get("demo/"+name[i]);
    
    BulkGM1000->Scale(19.6*0.0000851*1000*1000/(2*131716));
    BulkGM1500->Scale(19.6*0.0000044*1000*20000/(2*124273));
    BulkGM2000->Scale(19.6*0.0000004*1000*200000/(2*121429));
    DYJetsToLL_PtZ70To100->Scale(19.6*62.9*1000/11764538);
    DYJetsToLL_PtZ100->Scale(19.6*39.1*1000/12511326);
    TT->Scale(19.6*225.197*1000/21675970.);
    WW->Scale(19.6*57.1097*1000/10000431);
    WZ->Scale(19.6*33.21*1000/9955839);
    ZZ->Scale(19.6*8.059*1000/9799908);
    WJetsPt50To70->Scale(19.6*1001.0*1000/24950166);
    WJetsPt70To100->Scale(19.6*529.3*1000/20916010);
    WJetsPt100->Scale(19.6*282.5*1000/12106534);

    WJetsPt50To70->Add(WJetsPt70To100);
    WJetsPt50To70->Add(WJetsPt100);

    cout<<"Integral for "<<name[i]<<" (range 500-2200) "<<DYJetsToLL_PtZ70To100->Integral(501,2200)+DYJetsToLL_PtZ100->Integral(501,2200)+TT->Integral(501,2200)+
      WW->Integral(501,2200)+WZ->Integral(501,2200)+ZZ->Integral(501,2200)+WJetsPt50To70->Integral(501,2200)<<endl; 
    cout<<"Integral for "<<name[i]<<" (range 900-1100) "<<DYJetsToLL_PtZ70To100->Integral(901,1100)+DYJetsToLL_PtZ100->Integral(901,1100)+TT->Integral(901,1100)+
      WW->Integral(901,1100)+WZ->Integral(901,1100)+ZZ->Integral(901,1100)+WJetsPt50To70->Integral(901,1100)<<endl; 
    cout<<"Integral for "<<name[i]<<" (range 1400-1600) "<<DYJetsToLL_PtZ70To100->Integral(1401,1600)+DYJetsToLL_PtZ100->Integral(1401,1600)+TT->Integral(1401,1600)+
      WW->Integral(1401,1600)+WZ->Integral(1401,1600)+ZZ->Integral(1401,1600)+WJetsPt50To70->Integral(1401,1600)<<endl; 
    cout<<"Integral for "<<name[i]<<" (range 1900-2100) "<<DYJetsToLL_PtZ70To100->Integral(1901,2100)+DYJetsToLL_PtZ100->Integral(1901,2100)+TT->Integral(1901,2100)+
      WW->Integral(1901,2100)+WZ->Integral(1901,2100)+ZZ->Integral(1901,2100)+WJetsPt50To70->Integral(1901,2100)<<endl; 
    
    BulkGM1000->Rebin(rebin);
    BulkGM1500->Rebin(rebin);
    BulkGM2000->Rebin(rebin);
    DYJetsToLL_PtZ100->Rebin(rebin);
    DYJetsToLL_PtZ70To100->Rebin(rebin);
    TT->Rebin(rebin);
    WW->Rebin(rebin);
    WZ->Rebin(rebin);
    ZZ->Rebin(rebin);
    WJetsPt50To70->Rebin(rebin);

    DYJetsToLL_PtZ100->SetFillColor(kBlue);
    DYJetsToLL_PtZ70To100->SetFillColor(kGray);
    TT->SetFillColor(kGreen+1);
    WW->SetFillColor(kBlue-5);
    WZ->SetFillColor(kOrange+1);
    ZZ->SetFillColor(kAzure-9);
    WJetsPt50To70->SetFillColor(kMagenta);
    
    BulkGM1000->SetLineColor(kRed);
    BulkGM1000->SetLineWidth(3);
    BulkGM1500->SetLineColor(kGreen+3);
    BulkGM1500->SetLineWidth(3);
    BulkGM2000->SetLineColor(1);
    BulkGM2000->SetLineWidth(3);
    
    
    THStack *hs = new THStack("","");
    hs->SetTitle();
    hs->Add(DYJetsToLL_PtZ100);
    hs->Add(DYJetsToLL_PtZ70To100);
    hs->Add(TT);
    hs->Add(WW);
    hs->Add(WZ);
    hs->Add(ZZ);
    hs->Add(WJetsPt50To70);
    
    hs->Draw();
    hs->GetYaxis()->SetTitleOffset(0.80);
    hs->GetYaxis()->SetTitleSize(0.045);
    hs->GetXaxis()->SetTitleSize(0.045);
    hs->GetYaxis()->SetLabelSize(0.045);
    hs->GetXaxis()->SetLabelSize(0.045);
    hs->GetYaxis()->SetTitleOffset(1.1);   
    BulkGM1000->Draw("same");   
    BulkGM1500->Draw("same");   
    BulkGM2000->Draw("same");

    if(name[i]=="jetMassSelected") hs->GetXaxis()->SetRangeUser(60,120);

    TLatex latexLabel2;
    latexLabel2.SetTextSize(0.04);
    latexLabel2.SetTextFont(32);
    latexLabel2.SetNDC();
    latexLabel2.DrawLatex(0.6, 0.86, "L = 19.6 fb^{-1} at #sqrt{s} = 8 TeV");
    
    TLegend *pl = new TLegend(0.49,0.58,0.89,0.84);
    pl->SetTextSize(0.025); 
    pl->SetFillColor(0);
    TLegendEntry *ple = pl->AddEntry(BulkGM1000, "Bulk G (M=1000 GeV, #tilde{k}=0.2 (x 1000))",  "L");
    TLegendEntry *ple = pl->AddEntry(BulkGM1500, "Bulk G (M=1500 GeV, #tilde{k}=0.2 (x 20000))",  "L");
    TLegendEntry *ple = pl->AddEntry(BulkGM2000, "Bulk G (M=2000 GeV, #tilde{k}=0.2 (x 200000))",  "L");
    ple = pl->AddEntry(DYJetsToLL_PtZ100, "DYJetsToLL_PtZ100",  "F"); 
    ple = pl->AddEntry(DYJetsToLL_PtZ70To100, "DYJetsToLL_PtZ70To100",  "F"); 
    ple = pl->AddEntry(TT, "TT",  "F"); 
    ple = pl->AddEntry(WW, "WW",  "F"); 
    ple = pl->AddEntry(WZ, "WZ",  "F"); 
    ple = pl->AddEntry(ZZ, "ZZ",  "F"); 
    ple = pl->AddEntry(WJetsPt50To70, "WJet",  "F");
    pl->Draw();

    
    if(i>=13 && i<=16) hs->SetMaximum(3.0);
    if(i>=17 && i<=20) hs->SetMaximum(6.0);
    if(i>=21 && i<=24) hs->SetMaximum(5.0);
    if(i>=25 && i<=36) hs->SetMaximum(5.0);

    if(i<90) hs->SetMinimum(0);
    if(i>=90) c1->SetLogy();
    c1->SaveAs(name[i]+".png");
  }
    
    TCanvas * c1 = new TCanvas("c1", "c1", 800, 600);
    TH1F *BulkGM1000Den;              BulkGM1000Den = (TH1F*) f11->Get("demo/tauDecay");
    TH1F *BulkGM1500Den;              BulkGM1500Den = (TH1F*) f12->Get("demo/tauDecay");
    TH1F *BulkGM2000Den;              BulkGM2000Den = (TH1F*) f13->Get("demo/tauDecay");
    TH1F *BulkGM1000NumSel;              BulkGM1000NumSel = (TH1F*) f11->Get("demo/tauDecaySelected");
    TH1F *BulkGM1500NumSel;              BulkGM1500NumSel = (TH1F*) f12->Get("demo/tauDecaySelected");
    TH1F *BulkGM2000NumSel;              BulkGM2000NumSel = (TH1F*) f13->Get("demo/tauDecaySelected");

    BulkGM1000NumSel->Divide(BulkGM1000Den);
    BulkGM1500NumSel->Divide(BulkGM1500Den);
    BulkGM2000NumSel->Divide(BulkGM2000Den);

    BulkGM1000NumSel->SetLineColor(kRed);
    BulkGM1000NumSel->SetLineWidth(2);
    BulkGM1500NumSel->SetLineColor(kGreen+3);
    BulkGM1500NumSel->SetLineWidth(2);
    BulkGM2000NumSel->SetLineColor(1);
    BulkGM2000NumSel->SetLineWidth(2);
    
    TString channels[6];
    channels[0]="ee events"; channels[1]="e#mu events"; channels[2]="e#tau_{h} events"; 
    channels[3]="#mu#mu events"; channels[4]="#mu#tau_{h} events"; channels[5]="#tau_{h}#tau_{h} events"; //channels[6]="no tau";
    for (int ich=0;ich<6;ich++){
      BulkGM2000NumSel->GetXaxis()->SetBinLabel(ich+2,channels[ich]);
    }
    BulkGM2000NumSel->GetYaxis()->SetTitle("Efficiency");
    BulkGM2000NumSel->GetXaxis()->SetTitle("#tau#tau decay");
    BulkGM2000NumSel->GetYaxis()->SetTitleSize(0.045);
    BulkGM2000NumSel->GetXaxis()->SetTitleSize(0.045);
    BulkGM2000NumSel->GetXaxis()->SetLabelSize(0.040);
    BulkGM2000NumSel->GetYaxis()->SetLabelSize(0.045);
    BulkGM2000NumSel->GetYaxis()->SetTitleOffset(1.05);
    BulkGM2000NumSel->GetXaxis()->SetRangeUser(0,7);
    BulkGM2000NumSel->SetMinimum(0);

    BulkGM1000NumSel->SetMarkerStyle(21);
    BulkGM1000NumSel->SetMarkerSize(1.7);
    BulkGM1000NumSel->SetMarkerColor(kRed);
    BulkGM1500NumSel->SetMarkerStyle(21);
    BulkGM1500NumSel->SetMarkerSize(1.7);
    BulkGM1500NumSel->SetMarkerColor(kGreen+3);
    BulkGM2000NumSel->SetMarkerStyle(21);
    BulkGM2000NumSel->SetMarkerSize(1.7);
    BulkGM2000NumSel->SetMarkerColor(1);

    BulkGM1000NumSel->SetBinContent(1,-10000);
    BulkGM1500NumSel->SetBinContent(1,-10000);
    BulkGM2000NumSel->SetBinContent(1,-10000);
    BulkGM1000NumSel->SetBinContent(8,-10000);
    BulkGM1500NumSel->SetBinContent(8,-10000);
    BulkGM2000NumSel->SetBinContent(8,-10000);
    BulkGM1000NumSel->SetBinContent(9,-10000);
    BulkGM1500NumSel->SetBinContent(9,-10000);
    BulkGM2000NumSel->SetBinContent(9,-10000);
    BulkGM2000NumSel->SetMaximum(0.6);
    BulkGM2000NumSel->Draw();
    BulkGM1000NumSel->Draw("same histo text1");
    BulkGM1500NumSel->Draw("same histo text1");
    BulkGM2000NumSel->Draw("same histo text1");

    TLatex latexLabel2;
    latexLabel2.SetTextSize(0.04);
    latexLabel2.SetTextFont(32);
    latexLabel2.SetNDC();
    latexLabel2.DrawLatex(0.6, 0.86, "L = 19.6 fb^{-1} at #sqrt{s} = 8 TeV");
    
    TLegend *pl = new TLegend(0.605,0.70,0.895,0.84);
    pl->SetTextSize(0.025); 
    pl->SetFillColor(0);
    TLegendEntry *ple = pl->AddEntry(BulkGM1000NumSel, "Bulk G (M=1000 GeV, #tilde{k}=0.2)",  "L");
    ple = pl->AddEntry(BulkGM1500NumSel, "Bulk G (M=1500 GeV, #tilde{k}=0.2)",  "L");
    ple = pl->AddEntry(BulkGM2000NumSel, "Bulk G (M=2000 GeV, #tilde{k}=0.2)",  "L"); 
    pl->Draw();

    c1->SaveAs("tauDecaySelected.png");


}

