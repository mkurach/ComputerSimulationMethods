#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TLine.h>
#include <TCanvas.h>

void ploty() {

    ifstream file;
    file.open("omega1.txt");

    Float_t t,norm,en,sr;
    TGraph *grNorm = new TGraph();
    TGraph *grEn = new TGraph();
    TGraph *grSr = new TGraph();

    grNorm->GetXaxis()->SetTitle("t");
    grEn->GetXaxis()->SetTitle("t");
    grSr->GetXaxis()->SetTitle("t");

    grNorm->GetYaxis()->SetTitle("|#psi|^{2}");
    grEn->GetYaxis()->SetTitle("#epsilon");
    grSr->GetYaxis()->SetTitle("#bar{x}");


    grNorm->GetXaxis()->SetTitleSize(0.08);
    grEn->GetXaxis()->SetTitleSize(0.08);
    grSr->GetXaxis()->SetTitleSize(0.04);

    grNorm->GetYaxis()->SetTitleSize(0.08);
    grEn->GetYaxis()->SetTitleSize(0.08);
    grSr->GetYaxis()->SetTitleSize(0.08);

    grNorm->GetXaxis()->SetTitleOffset(0.3);
    grEn->GetXaxis()->SetTitleOffset(0.3);
    grSr->GetXaxis()->SetTitleOffset(0.3);

    grNorm->GetYaxis()->SetTitleOffset(0.3);
    grEn->GetYaxis()->SetTitleOffset(0.3);
    grSr->GetYaxis()->SetTitleOffset(0.3);
    int i = 0;
    while(file>>t && file>>norm && file>>en && file>>sr){
        grNorm->SetPoint(i,t,norm);
        grEn->SetPoint(i,t,en);
        grSr->SetPoint(i,t,sr);
        i++;
    }

    
    TCanvas *can = new TCanvas("can","can",800,800);
    can->Divide(1,3);
    can->cd(1);
    grNorm->Draw("AP");
    can->cd(2);
    grEn->Draw("AP");
    can->cd(3);
    grSr->Draw("AP");

    can->SaveAs("omega1.jpg");

    file.close();

}