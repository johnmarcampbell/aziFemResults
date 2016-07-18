#include "TString.h"
void correctHistogram(TH3D** exp, TH3D** corr, const Double_t res);

void doCorrectHistogram(const TString inFileName = "19Fit_FixedPID_bin2.root"
                        , const TString outFileName = "test.root"
                        , const Double_t res = 0.5)
{

    TFile* fIn = new TFile(inFileName.Data(),"READ");
    TFile* fOut = new TFile(outFileName.Data(),"RECREATE");

    const Int_t nPhiBins = 8;

    Int_t phiLabels[8] = {0, 22, 45, 67, 90, 112, 135, 157};

    TH3D* hExp[6][nPhiBins];
    TH3D* hCorr[6][nPhiBins];

    TString nameBase[6] = {"NumPiPlus_phi","DenPiPlus_phi","CoulPiPlus_phi","NumPiMinus_phi","DenPiMinus_phi","CoulPiMinus_phi"};
    TString histNames[6][8];

    for(Int_t i = 0; i <= 5; i++)
    {
        for(Int_t j = 0; j <= (nPhiBins - 1); j++)
        {
            histNames[i][j] += nameBase[i];
            histNames[i][j] += phiLabels[j];
            histNames[i][j] += "_kt4";
            fIn->cd();
            hExp[i][j] = (TH3D*)fIn->Get(histNames[i][j].Data());
        }

        correctHistogram(hExp[i], hCorr[i], res);

        for(Int_t j = 0; j <= (nPhiBins - 1); j++)
        {
            fOut->cd();
            hCorr[i][j]->Write(hCorr[i][j]->GetName(),2);

        }
    }

    fOut->Close();

}

void correctHistogram(TH3D** exp, TH3D** corr, const Double_t res)
{

    const Int_t nPhiBins = 8;
    const Double_t pi = TMath::Pi();
    const Double_t arc = pi/8.;
    const Double_t phi[8] = {0, 1*arc, 2*arc, 3*arc, 4*arc, 5*arc, 6*arc, 7*arc};
    const Int_t nqBins = exp[0]->GetNbinsX();
    const Double_t qLow = exp[0]->GetXaxis()->GetXmin(), qHigh = exp[0]->GetXaxis()->GetXmax();

    // Clone the histograms and reset them, to get empty histo's with the same binning scheme
    for(Int_t i = 0; i <= (nPhiBins - 1); i++)
    {
        TString corrName(exp[i]->GetName());
        TString newExpName = TString::Format("%s_%i", exp[i]->GetName(), i);
        exp[i]->SetName(newExpName.Data());
        corr[i] = new TH3D(corrName.Data(),corrName.Data(),nqBins,qLow,qHigh,nqBins,qLow,qHigh,nqBins,qLow,qHigh);
        corr[i]->Sumw2();
    }

    TH3D* hcExp = new TH3D("hcExp","hcExp",nqBins,qLow,qHigh,nqBins,qLow,qHigh,nqBins,qLow,qHigh); 
    TH3D* hsExp = new TH3D("hsExp","hsExp",nqBins,qLow,qHigh,nqBins,qLow,qHigh,nqBins,qLow,qHigh); 
    hcExp->Sumw2();

    // Calculate the Fourier components
    for(Int_t i = 0; i <= (nPhiBins - 1); i++)
    {
        hcExp->Add(exp[i],cos(2.*phi[i]));
        hsExp->Add(exp[i],sin(2.*phi[i]));
    }

    hsExp->Scale((1./nPhiBins));
    hcExp->Scale((1./nPhiBins));

    // Calculate zeta parameter
    const Double_t arg = (pi / nPhiBins);
    Float_t zeta = arg / (sin(arg) * res) - 1; 

    // Put it all together 
    for(Int_t i = 0; i <= (nPhiBins - 1); i++)
    {
        corr[i]->Add(hcExp, cos(2 * phi[i]));
        corr[i]->Add(hsExp, sin(2 * phi[i]));
        corr[i]->Scale(2 * zeta);
        corr[i]->Add(exp[i]);

    }

    delete hcExp; 
    delete hsExp; 


}
