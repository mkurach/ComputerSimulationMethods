#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TLine.h>
#include <TCanvas.h>
#include <math.h>

Double_t fitFunction(Double_t *x, Double_t *par) { //  0 - gamma, 1 - omega, 2 - normalizacja
    return par[3]+par[2]*(par[0]/2/M_PI/((x[0]-par[1])*(x[0]-par[1])+par[0]*par[0]/4));

}


void plotEn() {

    ifstream file;
    file.open("wyniki.txt");

    Float_t omega,en;
    TGraph *gr = new TGraph();

    gr->GetXaxis()->SetTitle("#frac{#omega}{#omega_{0}}");
    gr->GetYaxis()->SetTitle("#epsilon_{max}");

    int i = 0;
    while(file>>omega && file>>en){
        gr->SetPoint(i,omega,en);
        i++;
    }

    TF1 *fun = new TF1("fun", fitFunction,0.9,1.1,4);
    fun->SetParameter(0,1);
    fun->SetParameter(1,1);
    fun->SetParameter(2,1);

    gr->Fit(fun);

    TCanvas *can = new TCanvas("can","can",800,800);
    gr->Draw();
    can->SaveAs("en.jpg");
    file.close();

}