#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TLine.h>
#include <TCanvas.h>
void wykresy() {
    //FILLING FROM TXT FILE
    ifstream file;
    file.open("outputMain.txt");

    Float_t t, Ecal, Epot, Ekin, T, P;
    TGraph *gr[5];
    TString names[5] = {"Ecal","Epot", "Ekin", "temp", "p"};
    TString titles[5] = {"E_{całkowita} (t)","E_{potencjalna} (t)", "E_{kinetyczna} (t)", "Temperatura w funkcji czasu", "Ciśnienie w funkcji czasu"};
    TString xLabel = "t (ps)";
    TString yLabels[5] = {"E_{cal} (kJ/mol)", "E_{pot} (kJ/mol)", "E_{kin} (kJ/mol)", "T (K)", "p (#times 16.6 atm)"};
    TCanvas* can[5];
    for (int i = 0; i < 5; i ++) {
        gr[i] = new TGraph();
        can[i] = new TCanvas(names[i].Data(),names[i].Data(),800,800);
        can[i]->SetLeftMargin(0.15);
        gr[i]->GetXaxis()->SetTitle(xLabel.Data());
        gr[i]->GetYaxis()->SetTitle(yLabels[i].Data());
        gr[i]->SetTitle("");
        gr[i]->SetLineColor(kBlue);
        gr[i]->SetLineWidth(2);
    }

    int i = 0;
    while(file>>t && file>>Ecal && file>>Epot && file>>Ekin && file>>T && file>>P){
        gr[0]->SetPoint(i,t,Ecal);
        gr[1]->SetPoint(i,t,Epot);
        gr[2]->SetPoint(i,t,Ekin);
        gr[3]->SetPoint(i,t,T);
        gr[4]->SetPoint(i,t,P);
        i++;
    }
    

    file.close();

    for (int i = 0; i < 5; i ++) {
        can[i]->cd();
        gr[i]->Draw();
        can[i]->SaveAs(Form("%s.png",names[i].Data()));
    }



}