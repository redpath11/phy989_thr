#include <TString.h>
#include <TH2.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>

void plotNN(){
  TString filename = "../output/l0ps.txt";
  ifstream input; input.open(filename.Data());
  if(input.fail()){
    cout << "Error opening " << filename.Data() << endl;
    return;
  }
  else{
    /* Set experimental points and errors */
    Double_t Nijmegen[11] = {62.068,63.63,59.96,
    			50.90,40.54,26.78,
			16.94,8.94,1.96,
			-4.46,-10.59};
    Double_t NijmegenError[11] = {0.03,0.08,0.11,
    				0.19,0.28,0.38,
				0.41,0.39,0.37,
				0.43,0.62};
    const Double_t conv = 180./M_PI;
    Int_t nlines = 0;
    Double_t eLab,phaseShift;
    vector<double> x;
    vector<double> y;
    vector<double> ydiff;
    while(1){
      input >> eLab >> phaseShift;
      if(!input.good()){break;}
      nlines++;
/*      if(nlines <=3){cout << setw(12) << eLab
      			<< setw(12) << phaseShift << endl;
      }*/
      x.push_back(eLab); y.push_back(phaseShift*conv);
      ydiff.push_back(abs((phaseShift*conv - Nijmegen[nlines-1])/Nijmegen[nlines-1]));
    }


    TGraph *g = new TGraph(x.size(),&(x[0]),&(y[0]));
    g->SetMarkerColor(4); g->SetMarkerStyle(20); g->SetDrawOption("P"); g->SetFillStyle(0);
    g->SetTitle("#delta_{0}");

    TGraph *gd = new TGraph(x.size(),&(x[0]),&(ydiff[0]));
    gd->SetMarkerColor(1); gd->SetMarkerStyle(20); gd->SetDrawOption("P"); gd->SetFillStyle(0);
    gd->SetTitle("#delta_{0}");
    gd->SetTitle("Error(#delta_{0});T_{lab} [MeV];|#delta_{0} - #delta_{exp}| / #delta _{exp}");


    Double_t xError[11];
    for(int i=0;i<11;i++){xError[11]=0.;}
    TGraphErrors *n = new TGraphErrors(x.size(),&(x[0]),Nijmegen,xError,NijmegenError);
    n->SetMarkerColor(1);n->SetDrawOption("P");n->SetFillStyle(0);// n->SetMarkerStyle(20);
    n->SetTitle("#delta _{exp}");

    TMultiGraph *mg = new TMultiGraph();
    mg->Add(g,"lp"); mg->Add(n,"lp");
/**/
    TCanvas *c0 = new TCanvas("c0","c0",500,800);c0->Divide(1,2);
    c0->cd(1); mg->Draw("A");
    mg->GetXaxis()->SetTitle("E _{lab} [MeV]"); mg->GetYaxis()->SetTitle("#delta _{0} [deg]");
    gPad->BuildLegend(0.6,0.66,0.88,0.88);
    c0->cd(2); gd->Draw("AP");gd->GetYaxis()->SetTitleOffset(1.4);

  }
}


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
    Int_t color=1;
    if(name.CompareTo("NP")==0){
      g = new TGraph(x.size(),&(x[0]),&(yNP[0]));
//      color = 1;
    }
    else if((name.CompareTo("eftLO")==0)||(name.CompareTo("eftNLO")==0)||(name.CompareTo("eftNNLO")==0)){
      g = new TGraph(x.size(),&(x[0]),&(yEFT[0]));
      color = 4;
    }
    else{cout << "Plot type not NP, eftLO, eftNLO, eftNNLO." << endl;g = new TGraph();return g;}

    g->SetName(name);
    g->SetMarkerColor(color); g->SetMarkerStyle(20); g->SetDrawOption("P"); g->SetFillStyle(0);
    g->SetTitle(Form("#delta_{0%s}",name.Data()));

    error->~TGraph();
    error = new TGraph(x.size(),&(x[0]),&(ydiff[0]));
    error->SetName(Form("%serror",name.Data()));
    error->SetMarkerColor(color); error->SetMarkerStyle(); error->SetDrawOption("P"); error->SetFillStyle(0);
    error->SetTitle(Form("| (#delta_{0NP} - #delta_{0%s}) / #delta_{0NP}|",name.Data()));

    return g;

  }
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



TCanvas * plotFS(TString filename){

//  TString filename = "../output/FiniteSpherel0ps.txt";
  ifstream input; input.open(filename.Data());
  TCanvas *c1 = new TCanvas("c1","c1",500,800);c1->Divide(1,2);
  if(input.fail()){
    cout << "Error opening " << filename.Data() << endl;
    return c1;
  }
  else{
    const Double_t conv = 180./M_PI;
    Int_t nlines = 0;
    Double_t eLab,phaseShift,AphaseShift;
    vector<double> x;
    vector<double> y1;
    vector<double> y2;
    vector<double> ydiff;
    while(1){
      input >> eLab >> phaseShift >> AphaseShift;
      if(!input.good()){break;}
      nlines++;
      if(nlines <=3){cout << setw(12) << eLab
      			<< setw(12) << phaseShift
			<< setw(12) << AphaseShift << endl;
      }
      x.push_back(197.*eLab); y1.push_back(phaseShift*conv);
      y2.push_back(AphaseShift*conv);
      ydiff.push_back(abs(phaseShift - AphaseShift) / AphaseShift);
    }

    TGraph *g1 = new TGraph(x.size(),&(x[0]),&(y1[0]));
    TGraph *g2 = new TGraph(x.size(),&(x[0]),&(y2[0]));
    TGraph *g3 = new TGraph(x.size(),&(x[0]),&(ydiff[0]));
    g1->SetMarkerColor(4); g1->SetMarkerStyle(20); g1->SetDrawOption("P"); g1->SetFillStyle(0);
    g1->SetTitle("#delta_{0}");
    g1->SetName("#delta_{0}");
    g2->SetMarkerColor(2); g2->SetMarkerStyle(20); g2->SetDrawOption("P"); g2->SetFillStyle(0);
    g2->SetTitle("Analytic #delta_{0}");
    g2->SetName("#delta_{0} ^{A}");
    g3->SetMarkerColor(1); g3->SetMarkerStyle(20); g3->SetDrawOption("P"); g3->SetFillStyle(0);
    g3->SetTitle("Error(#delta_{0});T_{lab} [MeV];|#delta_{0} - #delta_{A}| / #delta_{A}");
    g3->SetName("Error(#delta_{0})");
    TMultiGraph *mg = new TMultiGraph();
    mg->Add(g1,"lp");
    mg->Add(g2,"lp");
//    mg->Add(g3,"lp");
    c1->cd(1);
    mg->Draw("AP");
    mg->GetXaxis()->SetTitle("E _{lab} [MeV]"); mg->GetYaxis()->SetTitle("#delta _{0} [deg]");
    gPad->BuildLegend(0.6,0.66,0.88,0.88);
    c1->cd(2); g3->Draw("AP"); g3->GetYaxis()->SetTitleOffset(1.4);
    }

    return c1;

}


void plotV(){

  TString filename = "../output/checkPotential.txt";
  ifstream input; input.open(filename.Data());
  if(input.fail()){
    cout << "Error opening " << filename.Data() << endl;
    return;
  }
  else{
    Int_t nlines=0;
    Int_t xi=0;
    Int_t yj=0;
    Double_t potential = 0.;
    TH2F *h = new TH2F("h","potential",101,1,102,101,1,102);

    while(1){
      input >> xi >> yj >> potential;
      h->SetBinContent(xi+1,yj+1,potential);
      if(!input.good()){break;}
      nlines++;
    }

    cout << nlines <<  " points read." << endl;

    input.close();

  TCanvas *c2 = new TCanvas("c2","c2",500,400);
  h->Draw("LEGO2");
//  h->Draw("COLZ");
  }
}


void plotA(){

  TString filename = "../output/checkAbegin.txt";
  ifstream input; input.open(filename.Data());
  if(input.fail()){
    cout << "Error opening " << filename.Data() << endl;
    return;
  }
  else{
    Int_t nlines=0;
    Int_t xi=0;
    Int_t yj=0;
    Double_t potential = 0.;
    TH2F *h = new TH2F("h","A",101,1,102,101,1,102);

    while(1){
      input >> xi >> yj >> potential;
      h->SetBinContent(xi+1,yj+1,potential);
      if(!input.good()){break;}
      nlines++;
    }

    cout << nlines <<  " points read." << endl;

    input.close();

  TCanvas *c3 = new TCanvas("c3","c3",500,400);
  h->Draw("LEGO2");
//  h->Draw("COLZ");
  }
}


void plotAinv(){

  TString filename = "../output/checkAinv.txt";
  ifstream input; input.open(filename.Data());
  if(input.fail()){
    cout << "Error opening " << filename.Data() << endl;
    return;
  }
  else{
    Int_t nlines=0;
    Int_t xi=0;
    Int_t yj=0;
    Double_t potential = 0.;
    TH2F *h = new TH2F("h","A",101,1,102,101,1,102);

    while(1){
      input >> xi >> yj >> potential;
      h->SetBinContent(xi+1,yj+1,potential);
      if(!input.good()){break;}
      nlines++;
    }

    cout << nlines <<  " points read." << endl;

    input.close();

  TCanvas *c4 = new TCanvas("c4","c4",500,400);
  h->Draw("LEGO2");
//  h->Draw("COLZ");
  }
}


void plotAtimeAinv(){
  TString filename    = "../output/checkIdentity.txt";
  ifstream input;    input.open(filename.Data());
  if(input.fail()){
    cout << "Error opening " << filename.Data() << endl;
    return;
  }
  else{
    Int_t nlines=0;
    Int_t xi=0;
    Int_t yj=0;
    Double_t potential = 0.;
    TH2F *h = new TH2F("h","A",101,1,102,101,1,102);

    while(1){
      input >> xi >> yj >> potential;
      h->SetBinContent(xi+1,yj+1,potential);
      if(!input.good()){break;}
      nlines++;
    }

    cout << nlines <<  " points read." << endl;

    input.close();

  TCanvas *c5 = new TCanvas("c5","c5",500,400);
  h->Draw("LEGO2");
//  h->Draw("COLZ");
  }
}


void plotNN2(){
  TString filename = "../output/NPd0.txt";
  ifstream input; input.open(filename.Data());
  if(input.fail()){
    cout << "Error opening " << filename.Data() << endl;
    return;
  }
  else{
    const Double_t conv = 180./M_PI;
    Int_t nlines = 0;
    Double_t eLab,phaseShift;
    vector<double> x;
    vector<double> y;
    while(1){
      input >> eLab >> phaseShift;
      if(!input.good()){break;}
      nlines++;
      if(nlines <=3){cout << setw(12) << eLab
      			<< setw(12) << phaseShift << endl;
      }
      x.push_back(eLab); y.push_back(phaseShift);
    }

    TGraph *g = new TGraph(x.size(),&(x[0]),&(y[0]));
    g->SetMarkerColor(4); g->SetMarkerStyle(20);
    g->SetTitle("#delta_{0} : T_{lab};T_{lab} [MeV];#delta_{0}");
    g->Draw("AP");


  }
}


