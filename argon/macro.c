
#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TLine.h>
#include <TCanvas.h>
void macro() {
    //FILLING FROM TXT FILE
    ifstream file;
    file.open("pedy.txt");

    Float_t px,py,pz;
    TH1D* histx = new TH1D("x","x",150,-7,7);
    TH1D* histy = new TH1D("y","y",50,-50,50);
    TH1D* histz = new TH1D("z","z",50,-50,50);

    while(file>>px && file>>py && file>>pz){
        histx->Fill(px);
        histy->Fill(py);
        histz->Fill(pz);
    }

    histx->Draw();

}