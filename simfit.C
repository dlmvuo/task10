#include <fstream>
#include <iostream>

#include <TH1D.h>
#include <TCanvas.h>
#include <TMinuit.h>
#include <TF1.h>
#include <TFile.h>
#include <TArrayD.h>

double h1[100];
double h2[100];
double err1[100];
double err2[100];
double x[100];
double chi2;

double func(double x, double *par) {
    return par[0] * TMath::Gaus(x, par[1], par[2]);
}

void fcn(int &npar, double *gin, double &f, double *par, int iflag) {
    const int nbins = 100;
    double chisq = 0;
    chi2 = 0;
    for (int i = 0; i < nbins; i++) {
        double delta1  = (h1[i] - func(x[i], par) - par[3]) / err1[i];
        double delta2 = (h2[i] - par[3]) / err2[i];
        chisq += delta1 * delta1 + delta2 * delta2;
        chi2 += chisq;
    }
    chi2 /= 100.;
    f = chisq;
}

void readfromfile(std::string filename, TH1D *hist)
{
    std::ifstream file(filename, std::ios::in);
    double x;
    while (file >> x) {hist->Fill(x);}
}

void simfit() {
    TH1D* hist1 = new TH1D("h1", "Signal", 100, 500, 600);
    TH1D* hist2 = new TH1D("h2", "Noise", 100, 500, 600);

    readfromfile("./data_1.dat", hist1);
    readfromfile("./data_2.dat", hist2);
    
    for (int i = 0; i < 100; i++) {
        h1[i] = hist1->GetBinContent(i);
        h2[i] = hist2->GetBinContent(i);
    }
    
    for (int i = 0; i < 100; i++) x[i] = i + 500;
    
    //The errors values
    for (int i = 0; i < 100; i++) {
        if (h1[i] == 0) err1[i] = sqrt(3.09);
        else err1[i] = sqrt(h1[i]);
    }

    for (int i = 0; i < 100; i++) {
        if (h2[i] == 0) err2[i] = sqrt(3.09);
        else err2[i] = sqrt(h2[i]);
    }
    
    TMinuit *gMinuit = new TMinuit(4);
    gMinuit->SetFCN(fcn);

    double arglist[10];
    int ierflg = 0;

    arglist[0] = 1;
    gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);

    // Set starting values and step sizes for parameters
    static double vstart[4] = {1., 530., 10., 1.};
    static double step[4] = {0.1, 0.1 , 0.1 , 0.1};
    gMinuit->mnparm(0, "ampl", vstart[0], step[0], 0, 0, ierflg);
    gMinuit->mnparm(1, "mean gaus", vstart[1], step[1], 0, 0, ierflg);
    gMinuit->mnparm(2, "sigma gaus", vstart[2], step[2], 0, 0, ierflg);
    gMinuit->mnparm(3, "const", vstart[3], step[3], 0, 0, ierflg);

    // Now ready for minimization step
    arglist[0] = 500;
    arglist[1] = 0.1;
    gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);

    // Print results
    double amin, edm, errdef;
    int nvpar, nparx, icstat;
    gMinuit->mnstat(amin, edm, errdef, nvpar, nparx, icstat);
    gMinuit->mnprin(3, amin);

    double ampl, mean, sigma, Const;
    double ampl_err, mean_err, sigma_err, Const_err;
    gMinuit->GetParameter(0, ampl, ampl_err);
    gMinuit->GetParameter(1, mean, mean_err);
    gMinuit->GetParameter(2, sigma, sigma_err);
    gMinuit->GetParameter(3, Const, Const_err);

    auto f1 = new TF1("f1", "[0] * TMath::Gaus(x, [1], [2]) + [3]", 500, 600);
    f1->SetParameter(0, ampl);
    f1->SetParameter(1, mean);
    f1->SetParameter(2, sigma);
    f1->SetParameter(3, Const);

    auto f2 = new TF1("f2", "[0]", 500, 600);
    f2->SetParameter(0, Const);

    std::cout << "chi2: " << chi2 <<std::endl;
    std::cout << "number of events = " << -(int)(ampl * (sqrt(2*TMath::Pi()) * sigma)) << std::endl;


    TCanvas *canvas = new TCanvas("canvas", "title.png", 1200, 900);

    canvas->Divide(1, 2);
    canvas->cd(1);
    hist1->SetLineWidth(2);
    hist1->Draw("e");
    f1->Draw("same");
    canvas->cd(2);
    hist2->SetLineWidth(2);
    hist2->Draw("e");
    f2->Draw("same");

    TFile* ofile = new TFile("simfit.root", "recreate");
    hist1->Write();
    hist2->Write();
    f1->Write();
    f2->Write();
    ofile->Close();   
}