#include <TCanvas.h>
#include <TLegend.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TStyle.h>

#include <iostream>
#include <fstream>
#include <string>

double fitconst(double *x, double *par){
    return par[0];
}

double fitnorm(double *x, double *par){
    return 1/par[1]/sqrt(2*TMath::Pi())*exp(-1/2*((x[0]-par[0])/par[1])*((x[0]-par[0])/par[1]));
}

void readfromfile(std::string filename, TH1D *hist)
{
    std::ifstream file(filename, std::ios::in);
    double x;
    while (file >> x) {hist->Fill(x);}
}

void hists() {
    TH1D* hist1 = new TH1D("h1", "h1", 100, 500, 600);
    TH1D* hist2 = new TH1D("h2", "h2", 100, 500, 600);

    readfromfile("./data_1.dat", hist1);
    readfromfile("./data_2.dat", hist2);

    TF1* normcurve = new TF1("normcurve", "fitnorm", 500, 600, 2);
    TF1* line = new TF1("line", "fitconst", 500, 600, 1);

    TCanvas *canvas = new TCanvas("canvas", "title.png", 1200, 900);
    canvas->Divide(1, 2);

    canvas->cd(1);

    hist1->Fit("normcurve");
    hist1->Draw("e");

    canvas->cd(2);

    hist2->Fit("line");
    hist2->Draw("e");
}