#include <TString.h>
#include <TH2.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TStyle.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>

void palettea();

/* Make a TGraph from some txt file of phase shifts vs energy
 * where the phase shifts are caculated from the NP potential
 * or the EFT potential at either LO, NLO, NNLO.
*/
TGraph *eftGraph(TString name,TString filename,TGraph *error){

  TGraph *g = new TGraph();

  ifstream input; input.open(filename.Data());
  if(input.fail()){
    cout << "Error opening " << filename.Data() << endl;
    return g;
  }
  else{
    const Double_t conv = 180./M_PI;
    Int_t nlines = 0;
    Double_t eLab,NPphaseShift,EFTphaseShift;
    vector<double> x;
    vector<double> yNP;
    vector<double> yEFT;
    vector<double> ydiff;
    while(1){
      input >> eLab >> NPphaseShift >> EFTphaseShift;
      if(!input.good()){break;}
      nlines++;
      x.push_back(eLab); yNP.push_back(NPphaseShift); yEFT.push_back(EFTphaseShift);
      ydiff.push_back(abs((NPphaseShift - EFTphaseShift)/NPphaseShift));
    }
    input.close();
    g->~TGraph();
    TGraph *g;
    Int_t color=1; Int_t markerStyle = 5;
    if(name.CompareTo("NP")==0){
      g = new TGraph(x.size(),&(x[0]),&(yNP[0]));
      markerStyle = 5;
//      color = 1;
    }
    else if((name.CompareTo("eftLO")==0)||(name.CompareTo("eftNLO")==0)||(name.CompareTo("eftNNLO")==0)){
      g = new TGraph(x.size(),&(x[0]),&(yEFT[0]));
      color = 4; markerStyle = 20;
    }
    else{cout << "Plot type not NP, eftLO, eftNLO, eftNNLO." << endl;g = new TGraph();return g;}

    g->SetName(name);
    g->SetMarkerColor(color); g->SetMarkerStyle(markerStyle); g->SetDrawOption("P"); g->SetFillStyle(0);
    g->SetTitle(Form("#delta_{0%s}",name.Data()));

    error->~TGraph();
    error = new TGraph(x.size(),&(x[0]),&(ydiff[0]));
    error->SetName(Form("%serror",name.Data()));
    error->SetMarkerColor(color); error->SetMarkerStyle(); error->SetDrawOption("P"); error->SetFillStyle(0);
    error->SetTitle(Form("| (#delta_{0NP} - #delta_{0%s}) / #delta_{0NP}|",name.Data()));

    return g;

  }
}


TGraph *errGraph(TString name,TString filename){

  TGraph *error = new TGraph();

  ifstream input; input.open(filename.Data());
  if(input.fail()){
    cout << "Error opening " << filename.Data() << endl;
    return error;
  }
  else{
    const Double_t conv = 180./M_PI;
    Int_t nlines = 0;
    Double_t eLab,NPphase,EFTphase;
    vector<double> x;
    vector<double> y;
    while(1){
      input >> eLab >> NPphase >> EFTphase;
//      cout << eLab << "   " << NPphase << "   " << EFTphase << endl;
      if(!input.good()){break;}
      nlines++;
      x.push_back(eLab); y.push_back(abs((NPphase - EFTphase)/NPphase));
    }
    input.close();
    error->~TGraph();
    TGraph *error;
    Int_t color=1;
    error = new TGraph(x.size(),&(x[0]),&(y[0]));
    error->SetName(name);
    error->SetMarkerColor(color); error->SetMarkerStyle(21); error->SetDrawOption("p"); error->SetFillStyle(0);
    error->SetTitle(Form("#delta_{0%s}",name.Data()));
    return error;
  }
}


void Lepage(){

  TString filename1 = "../output/Lepage/LO.txt";
  TString filename2 = "../output/Lepage/NLO.txt";
  TString filename3 = "../output/Lepage/NNLO.txt";
/*
  TString filename1 = "../output/Lepage/cutoffMpion/LO.txt";
  TString filename2 = "../output/Lepage/cutoffMpion/NLO.txt";
  TString filename3 = "../output/Lepage/cutoffMpion/NNLO.txt";
*/
  TGraph *eLO,*eNLO,*eNNLO;
  eLO = errGraph("eftLO",filename1);
  Int_t color=1;Int_t linestyle=2;
  eLO->SetLineColor(color);eLO->SetMarkerColor(color);
  eLO->SetLineStyle(linestyle);
  eNLO= errGraph("eftNLO",filename2);
  color=2;linestyle=9;
  eNLO->SetLineColor(color);eNLO->SetMarkerColor(color);
  eNLO->SetLineStyle(linestyle);
  eNNLO = errGraph("eftNNLO",filename3);
  color=4;linestyle=1;
  eNNLO->SetLineColor(color);eNNLO->SetMarkerColor(color);
  eNNLO->SetLineStyle(linestyle);

  TMultiGraph *lepage = new TMultiGraph();
  lepage->Add(eLO,"l"); lepage->Add(eNLO,"l"); lepage->Add(eNNLO,"l");

  TCanvas *c0 = new TCanvas("c0","c0",500,400);
  c0->SetLogx(1); c0->SetLogy(1);
/**/
  lepage->Draw("A");
  lepage->GetXaxis()->SetTitle("log E _{lab} [MeV]");
  lepage->GetYaxis()->SetTitle("log | (#delta_{NP} - #delta_{EFT}) / #delta_{NP} |");
//  lepage->GetXaxis()->SetRangeUser(1,10);
  gPad->BuildLegend(0.6,0.66,0.88,0.88);

}


void cutoff(){
/*
  TString filename1 = "../output/cutoff/NLO.txt";
  TString filename2 = "../output/cutoff/NLOcSmallest.txt";
  TString filename3 = "../output/cutoff/NLOmid.txt";
  TString filename4 = "../output/cutoff/NLOlarge.txt";
*/
  TString filename1 = "../output/Lepage/NLO.txt";
  TString filename2 = "../output/Lepage/co1NLO.txt";
  TString filename3 = "../output/Lepage/co2NLO.txt";
  TString filename4 = "../output/Lepage/co3NLO.txt";
  TString filename5 = "../output/Lepage/co4NLO.txt";
  TString filename6 = "../output/Lepage/co5NLO.txt";

  Int_t width = 4;

  TGraph *eNLO1,*eNLO2,*eNLO3,*eNLO4,*eNLO5,*eNLO6;
  eNLO1 = errGraph("eftNLO",filename1);
  eNLO1->SetTitle("#Lambda_{c} = 0.71");
  Int_t color=2;Int_t linestyle=1;
  eNLO1->SetLineColor(color);eNLO1->SetMarkerColor(color);
  eNLO1->SetLineStyle(linestyle);
  eNLO1->SetLineWidth(width);

  eNLO2= errGraph("eftNLO",filename2);
  eNLO2->SetTitle("#Lambda_{c} = 0.07");
  color=14;linestyle=9;
  eNLO2->SetLineColor(color);eNLO2->SetMarkerColor(color);
  eNLO2->SetLineStyle(linestyle);
  eNLO2->SetLineWidth(width);

  eNLO3 = errGraph("eftNLO",filename3);
  eNLO3->SetTitle("#Lambda_{c} = 7.1");
  color=14;linestyle=2;
  eNLO3->SetLineColor(color);eNLO3->SetMarkerColor(color);
  eNLO3->SetLineStyle(linestyle);
  eNLO3->SetLineWidth(width);

  eNLO4 = errGraph("eftNLO",filename4);
  eNLO4->SetTitle("#Lambda_{c} = 2.1");
  color=8;linestyle=8;
  eNLO4->SetLineColor(color);eNLO4->SetMarkerColor(color);
  eNLO4->SetLineStyle(linestyle);
  eNLO4->SetLineWidth(width);

  eNLO5 = errGraph("eftNLO",filename5);
  eNLO5->SetTitle("#Lambda_{c} = 1.4");
  color=4;linestyle=2;
  eNLO5->SetLineColor(color);eNLO5->SetMarkerColor(color);
  eNLO5->SetLineStyle(linestyle);
  eNLO5->SetLineWidth(width);

  eNLO6 = errGraph("eftNLO",filename6);
  eNLO6->SetTitle("#Lambda_{c} = 0.35");
  color=41;linestyle=7;
  eNLO6->SetLineColor(color);eNLO6->SetMarkerColor(color);
  eNLO6->SetLineStyle(linestyle);
  eNLO6->SetLineWidth(width);

  TMultiGraph *lepage = new TMultiGraph();
  lepage->Add(eNLO2,"l");
  lepage->Add(eNLO6,"l");
  lepage->Add(eNLO1,"l");
  lepage->Add(eNLO5,"l"); 
  lepage->Add(eNLO4,"l");
  lepage->Add(eNLO3,"l");

  TCanvas *c0 = new TCanvas("c0","c0",500,400);
  c0->SetLogx(1); c0->SetLogy(1);
/**/
  lepage->Draw("A");
  lepage->GetXaxis()->SetTitle("log E _{lab} [MeV]");
  lepage->GetYaxis()->SetTitle("log | (#delta_{NP} - #delta_{EFT}) / #delta_{NP} |");
  gPad->BuildLegend(0.6,0.66,0.88,0.88);

}


void OnePion(){
  const char *order[3] = {"eftLO","eftNLO","eftNNLO"};
  TString file[3];

  file[0] = "../output/OnePion/co1/LO.txt";
  file[1] = "../output/OnePion/co1/NLO.txt";
  file[2] = "../output/OnePion/co1/NNLO.txt";

  TGraph *g0,*g[3];
  TGraph *e = new TGraph();

  TMultiGraph *mg[3];

  TCanvas *c20 = new TCanvas("c20","c20",500,1200);c20->Divide(1,3);

  g0 = eftGraph("NP", file[0], e);

  for(int i=0;i<3;i++){
    g[i] = eftGraph(order[i],file[i],e);
    mg[i] = new TMultiGraph();
    mg[i]->Add(g0,"lp"); mg[i]->Add(g[i],"lp");
    c20->cd(i+1); mg[i]->Draw("A");
    mg[i]->GetXaxis()->SetTitle("E _{lab} [MeV]");
    mg[i]->GetYaxis()->SetTitle("#delta _{0} [rad]");
    gPad->BuildLegend(0.6,0.20,0.88,0.40);
  }

}


void OnePionLepage(){
  TString file[4];
  Int_t color[4] = {1,2,4,8};
  Int_t lineStyle[4] = {1,9,2,7};
  const char *name[4] = {"#Lambda_{c} = 1.4","#Lambda_{c} = 2.8","#Lambda_{c} = 5.6","#Lambda_{c} = 10."};

  file[0] = "../output/OnePion/co2NLO.txt";
  file[1] = "../output/OnePion/co1NLO.txt";
  file[2] = "../output/OnePion/co3NLO.txt";
  file[3] = "../output/OnePion/co4NLO.txt";

  TGraph *g0,*g[4];
  TGraph *e = new TGraph();


  TCanvas *c21 = new TCanvas("c21","c21",500,400);
  c21->SetLogx(1); c21->SetLogy(1);

//  g0 = eftGraph("NP", file[0], e);

  TMultiGraph *mg = new TMultiGraph();

  for(int i=0;i<4;i++){
    g[i] = errGraph("eftNLO",file[i]);
    g[i]->SetTitle(name[i]);
    g[i]->SetLineColor(color[i]);
    g[i]->SetMarkerColor(color[i]);
    g[i]->SetLineStyle(lineStyle[i]);
    mg->Add(g[i],"l");
  }
  
  mg->Draw("A");
  mg->GetXaxis()->SetTitle("log E _{lab} [MeV]");
  mg->GetYaxis()->SetTitle("log | (#delta_{NP} - #delta_{EFT}) / #delta_{NP} |");
  mg->GetXaxis()->SetTitleOffset(1.4);
  mg->GetYaxis()->SetTitleOffset(1.4);
  gPad->BuildLegend(0.6,0.20,0.88,0.40);

}



void plotPionless(){
  const char *order[3] = {"eftLO","eftNLO","eftNNLO"};

  TString file[3];
  file[0] = "../output/Lepage/LO.txt";
  file[1] = "../output/Lepage/NLO.txt";
  file[2] = "../output/Lepage/NNLO.txt";

  TGraph *g0,*g[3];
  TGraph *e = new TGraph();

  TMultiGraph *mg[3];

  TCanvas *c10 = new TCanvas("c10","c10",500,1200);c10->Divide(1,3);

  g0 = eftGraph("NP", file[0], e);
  for(int i=0;i<3;i++){
    g[i] = eftGraph(order[i],file[i],e);
    mg[i] = new TMultiGraph();
    mg[i]->Add(g0,"lp"); mg[i]->Add(g[i],"lp");
    c10->cd(i+1); mg[i]->Draw("A");
    mg[i]->GetXaxis()->SetTitle("E _{lab} [MeV]");
    mg[i]->GetYaxis()->SetTitle("#delta _{0} [rad]");
    gPad->BuildLegend(0.6,0.20,0.88,0.40);
  }

}


void plotLO(TString filename1,Int_t color){
//  TString filename1 = "../output/LO.txt";

  TGraph *g0,*g1;
  TGraph *e2 = new TGraph();
  g0 = eftGraph("NP",filename1,e2);
  g1 = eftGraph("eftLO",filename1,e2);
  g1->SetMarkerColor(color);


  TMultiGraph *mg1 = new TMultiGraph();
  mg1->Add(g0,"lp"); mg1->Add(g1,"lp");
/**/
  TCanvas *c0 = new TCanvas("c0","c0",500,400);

  mg1->Draw("A");
  mg1->GetXaxis()->SetTitle("E _{lab} [MeV]"); mg1->GetYaxis()->SetTitle("#delta _{0} [rad]");
  gPad->BuildLegend(0.6,0.66,0.88,0.88);

}


void plotOnePionLO(){
  TString filename1 = "../output/LO.txt";

  TGraph *g0,*g1;
  TGraph *e2 = new TGraph();
  g0 = eftGraph("NP",filename1,e2);
  g1 = eftGraph("eftLO",filename1,e2);


  TMultiGraph *mg1 = new TMultiGraph();
  mg1->Add(g0,"lp"); mg1->Add(g1,"lp");
/**/
  TCanvas *c0 = new TCanvas("c0","c0",500,400);

  mg1->Draw("A");
  mg1->GetXaxis()->SetTitle("E _{lab} [MeV]"); mg1->GetYaxis()->SetTitle("#delta _{0} [rad]");
  gPad->BuildLegend(0.6,0.66,0.88,0.88);
//  gPad->SetLogx(1);

}



void plotNLO(TString filename2,Int_t color){
//  TString filename2 = "../output/Lepage/NLO.txt";

  TGraph *g0,*g2;
  TGraph *e2 = new TGraph();
  g0 = eftGraph("NP",filename2,e2);
  g2 = eftGraph("eftNLO",filename2,e2);
  g2->SetMarkerColor(color);


  TMultiGraph *mg2 = new TMultiGraph();
  mg2->Add(g0,"lp"); mg2->Add(g2,"lp");
/**/
  TCanvas *c0 = new TCanvas("c0","c0",500,400);

  mg2->Draw("A");
  mg2->GetXaxis()->SetTitle("E _{lab} [MeV]"); mg2->GetYaxis()->SetTitle("#delta _{0} [rad]");
  gPad->BuildLegend(0.6,0.66,0.88,0.88);
}


void plotNNLO(TString filename2,Int_t color){
//  TString filename2 = "../output/Lepage/NNLO.txt";

  TGraph *g0,*g2;
  TGraph *e2 = new TGraph();
  g0 = eftGraph("NP",filename2,e2);
  g2 = eftGraph("eftNNLO",filename2,e2);
  g2->SetMarkerColor(color);


  TMultiGraph *mg2 = new TMultiGraph();
  mg2->Add(g0,"lp"); mg2->Add(g2,"lp");
/**/
  TCanvas *c0 = new TCanvas("c0","c0",500,400);

  mg2->Draw("A");
  mg2->GetXaxis()->SetTitle("E _{lab} [MeV]"); mg2->GetYaxis()->SetTitle("#delta _{0} [rad]");
  gPad->BuildLegend(0.6,0.66,0.88,0.88);
}



void plotEFT(){
  TString filename1 = "../output/eftLO.txt";
  TString filename2 = "../output/eftNLO.txt";
  TString filename3 = "../output/eftNNLO.txt";

  TGraph *g0,*g1,*g2,*g3;
  TGraph *e0 = new TGraph();
  TGraph *e1 = new TGraph();
  TGraph *e2 = new TGraph();
  TGraph *e3 = new TGraph();
  g0 = eftGraph("NP",filename1,e0);
  g1 = eftGraph("eftLO",filename1,e1);
  g2 = eftGraph("eftNLO",filename2,e2);
  g3 = eftGraph("eftNNLO",filename1,e3);


  TMultiGraph *mg1 = new TMultiGraph();
  mg1->Add(g0,"lp"); mg1->Add(g1,"lp");
  TMultiGraph *mg2 = new TMultiGraph();
  mg2->Add(g0,"lp"); mg2->Add(g2,"lp");
  TMultiGraph *mg3 = new TMultiGraph();
  mg3->Add(g0,"lp"); mg3->Add(g3,"lp");
/**/
  TCanvas *c0 = new TCanvas("c0","c0",1500,800);c0->Divide(3,2);
  c0->cd(1); mg1->Draw("A");
  mg1->GetXaxis()->SetTitle("E _{lab} [MeV]"); mg1->GetYaxis()->SetTitle("#delta _{0} [rad]");
  gPad->BuildLegend(0.6,0.66,0.88,0.88);
  c0->cd(4); e1->Draw("AP");//e1->GetYaxis()->SetTitleOffset(1.4);

  c0->cd(2); mg2->Draw("A");
  mg2->GetXaxis()->SetTitle("E _{lab} [MeV]"); mg2->GetYaxis()->SetTitle("#delta _{0} [rad]");
  gPad->BuildLegend(0.6,0.66,0.88,0.88);
  c0->cd(5);// e2->Draw("AP");e2->GetYaxis()->SetTitleOffset(1.4);

  c0->cd(3); mg3->Draw("A");
  mg3->GetXaxis()->SetTitle("E _{lab} [MeV]"); mg3->GetYaxis()->SetTitle("#delta _{0} [rad]");
  gPad->BuildLegend(0.6,0.66,0.88,0.88);
  c0->cd(6);// e3->Draw("AP");e3->GetYaxis()->SetTitleOffset(1.4);

}


void plotLOfit(){
  TString fileC0 = "../output/eftLOfit.txt";
//  TString fileC0 = "../output/eftLOfit_equalWeights.txt";
  ifstream input; input.open(fileC0.Data());
  if(input.fail()){
    cout << "Error opening " << fileC0.Data() << endl;
    return;
  }
  else{
    Int_t nlines=0;
    Double_t C0=0.; Double_t X2=0.;
    vector<double> x;
    vector<double> y;
    while(1){
      input >> C0 >> X2;
      if(!input.good()){break;}
      nlines++;
      x.push_back(C0); y.push_back(X2);
    }
    TGraph *g0 = new TGraph(x.size(),&(x[0]),&(y[0]));
    g0->SetMarkerColor(434); g0->SetMarkerStyle(22);
    TCanvas *c0 = new TCanvas("c0","c0",500,400);
    g0->Draw("AP");
    g0->SetTitle("#chi ^{2} : C_{0};C_{0};#chi^{2}");
    g0->Draw("AP");
//    g0->Fit("pol2");
//    TF1 *f = (TF1*)g0->GetFunction("pol2");
//    Double_t min = f->GetMinimumX(-1.5,-0.5);
//    cout << "x at min: " << min << endl;

  }
}





void plotEFTfit(){
  TString fileC0 = "../output/eftNLOfit.txt";
//  TString fileC0 = "../output/eftLOfit_equalWeights.txt";
//  TString fileC0 = "../output/eftLOfit_Weighted.txt";
  ifstream input; input.open(fileC0.Data());
  if(input.fail()){
    cout << "Error opening " << fileC0.Data() << endl;
    return;
  }
  else{
    Int_t nlines=0;
    Double_t C0=0.; Double_t X2=0.;
    vector<double> x;
    vector<double> y;
    while(1){
      input >> C0 >> X2;
      if(!input.good()){break;}
      nlines++;
      x.push_back(C0); y.push_back(X2);
    }
    TGraph *g0 = new TGraph(x.size(),&(x[0]),&(y[0]));
    Double_t xLo = x[0]; Double_t xHi = x[x.size()-1];
    g0->SetMarkerColor(434); g0->SetMarkerStyle(22);
    TCanvas *c0 = new TCanvas("c0","c0",500,400);
    g0->Draw("AP");
    g0->SetTitle("#chi ^{2} : C;C;#chi^{2}");
    g0->Draw("AP");
    g0->Fit("pol2");
    TF1 *f = (TF1*)g0->GetFunction("pol2");
    Double_t min = f->GetMinimumX(xLo,xHi);
    cout << "x at min: " << min << endl;

  }
  input.close();
}


void plotNLOfit2D(Double_t c0start,Double_t c0inc,Double_t c2start,Double_t c2inc){
  Double_t c0end = c0start + (49. * c0inc);
  Double_t c2end = c2start + (49. * c2inc);
  TString fileC0C2 = "../output/eftNLOfit.txt";
  ifstream input; input.open(fileC0C2.Data());
  if(input.fail()){
    cout << "Error opening " << fileC0C2.Data() << endl;
    return;
  }
  else{
    Int_t nlines=0;
    Double_t C0=0.; Double_t C2=0.; Double_t X2=0.;
//    TH2F *hh = new TH2F("hh","#chi^{2}:C2:C0;C0;C2;#chi^{2}",25,-25,24,25,-25,24);// [-25,24]x[-25,24]
//    TH2F *hh = new TH2F("hh","#chi^{2}:C2:C0;C0;C2;#chi^{2}",50,-3.,0.92,50,-20,0.96);// [-5,4.6]x[0,98]
    TH2F *hh = new TH2F("hh","#chi^{2}:C2:C0;C0;C2;#chi^{2}",50,c0start,c0end,50,c2start,c2end);// [-3,0.92]x[-9,-6.45]
    while(1){
      input >> C0 >> C2 >> X2;
      Int_t xBin = hh->GetXaxis()->FindBin(C0);
      Int_t yBin = hh->GetYaxis()->FindBin(C2);

      hh->SetBinContent(xBin,yBin,X2);
      if(!input.good()){break;}
      nlines++;
    }

    cout << nlines << " points read." << endl;

    input.close();

    TCanvas *c0 = new TCanvas("c0","c0",500,400);
    palettea();
//    hh->Draw("LEGO2");
    hh->SetStats(0);
    hh->Draw("COLZ");

  }

  input.close();
}


/* Function to parse the output for eftSearch4 and write out
 * the C0,C2,C4,C4prime parameter values corresponding to the
 * minimum chi2 value.
*/
void FindChi2Min(){
  ifstream input; input.open("../output/eftNNLOfit.txt");
  if(input.fail()){
    cout << "Error opening ../output/eftNNLOfit.txt" << endl;
  }
  else{
    Int_t nlines=0;
    Int_t minCount=0;
    Double_t C0=0.; Double_t C2=0.; Double_t C4=0.; Double_t C4p=0.;
    Double_t X2=999.; Double_t minX2=999.;
    Double_t foundC0=0.; Double_t foundC2=0.; Double_t foundC4=0.; Double_t foundC4p=0.;
    while(1){
      input >> C0 >> C2 >> C4 >> C4p >> X2;
      if(X2<minX2){
        minX2 = X2;
	foundC0  = C0;
	foundC2  = C2;
	foundC4  = C4;
	foundC4p = C4p;
      }

      if(!input.good()){break;}
      nlines++;

    }

    cout << nlines << " points read." << endl;

    printf("minX2 = %g\n",minX2);
    printf("C0 = %g\nC2 = %g\nC4 = %g\n C4p = %g\n",foundC0,foundC2,foundC4,foundC4p);

  }

  input.close();

}



void palettea()
{
  const Int_t NRGBs = 6;
  const Int_t NCont = 999;

  Double_t stops[NRGBs] = { 0.00, 0.1, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.99, 0.0, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.0, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.99, 0.0, 1.00, 0.12, 0.00, 0.00 };

  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);

  const Double_t min=-5.0;
//  const Double_t max=55.;// allO, F-27 in Si 0
//  const Double_t max=50.;// allO, F-27 in Si 1
  const Double_t max=5.0;// allO, F-27 in Si 2 && Si [i] vs. IC
//  const Double_t max=50.;// allO_Oic3, F-27 in Si 0
//  const Double_t max=25.;// allO_Oic3, F-27 in Si 2
//  const Double_t max=8.;// nvelRecal_27F_O, gated on O-24
//  const Double_t max=7.;// nvelRecal_27F_O, gated on O-23
//  const Double_t max=6.;// nvelRecal_27F_O, gated on O-22

  const Int_t nLevels = 999;
  Double_t levels[nLevels];

  for(int i=1;i<nLevels;i++)
  {  levels[i] = min + (max-min) / (nLevels-1) * (i);  }
  levels[0] = -1.;
}


