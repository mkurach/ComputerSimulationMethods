#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TLine.h>
#include <TCanvas.h>
void wykres() {
    TGraph* gr = new TGraph();
    TCanvas* can = new TCanvas("niewiem","niewiem",800,800);
    can->SetLeftMargin(0.15);
    gr->GetXaxis()->SetTitle("n");
    gr->GetYaxis()->SetTitle("T_{top} (K)");
    gr->SetTitle("");
    gr-> SetMarkerStyle(20);
    gr -> SetMarkerSize(2);
    gr -> SetMarkerColor(kBlue);

    gr->SetPoint(0,6,120);
    gr->SetPoint(1,7,140);
    gr->SetPoint(2,8,150);


    can->cd();
    gr->Draw("AP");
    can->SaveAs("topnienie.png");
  



}