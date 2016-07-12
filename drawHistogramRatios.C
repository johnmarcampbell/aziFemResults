Bool_t isFirstIteration(Int_t iq2OrMult, Int_t iMult, Int_t iPhi);
Bool_t isLastIteration(Int_t iq2OrMult, Int_t deltaMult, Int_t deltaPhi);
void doProjections(TFile* lowFile, TFile* highFile, TH1D** ratio, TString phiLabel);
void projectHistogram(TH3D* fullHist, TH1D** projHist, Float_t projRange);
TCanvas* makeCFCanvas(TH1D** ratios, TString canvasLabel);

void drawHistogramRatios()
{

    // TString baseDirectory = "193GeVData";
    // TString species = "UUFemto193";
    // TString outFileNameBase = "UUFemto193CFRatios.pdf";
    TString baseDirectory = "200GeVData";
    TString species = "AuAuFemto200";
    TString outFileNameBase = "AuAuFemto200CFRatios.pdf";
    TString outFileName;
    TString fileNameBaseZdcLow = TString::Format("%s/0to0x25PercZdc", baseDirectory.Data());
    TString fileNameBaseZdcHigh = TString::Format("%s/0to0x5PercZdc", baseDirectory.Data());
    Int_t nMultBins = 5, nPhiBins = 8;
    // Int_t nMultBins = 1, nPhiBins = 1;
    Bool_t firstIteration = kTRUE, lastIteration = kFALSE;
    TString phiLabels[8] = {"0", "22", "45", "67", "90", "112", "135", "157"};
    TH1D* ratios[3] = {0};

    for (Int_t iq2OrMult = 0; iq2OrMult <= 1; ++iq2OrMult)
    {
        TString q2OrMultString;
        if (!iq2OrMult) {
            q2OrMultString = "q2";
        } else {
            q2OrMultString = "mult";
        }
        for (Int_t iMult = 0; iMult <= (nMultBins - 1); iMult++)
        {
            // TString fileNameZdcLow = TString::Format("%s/%s_%s_%i_Corrected.root", fileNameBaseZdcLow.Data(), species.Data(), q2OrMultString.Data(), iMult);
            // TString fileNameZdcHigh = TString::Format("%s/%s_%s_%i_Corrected.root", fileNameBaseZdcHigh.Data(), species.Data(), q2OrMultString.Data(), iMult);
            TString fileNameZdcLow = TString::Format("%s/%s_%s_%i_Corrected.root", fileNameBaseZdcLow.Data(), species.Data(), q2OrMultString.Data(), iMult);
            TString fileNameZdcHigh = TString::Format("%s/%s_%s_%i_Corrected.root", fileNameBaseZdcHigh.Data(), species.Data(), q2OrMultString.Data(), iMult);

            // --- Open Files --- //
            TFile* zdcLowFile = new TFile(fileNameZdcLow.Data(), "READ");
            if (zdcLowFile->IsZombie()) { 
                cout << "Unable to open " << fileNameZdcLow.Data() <<". It is a zombie.\n";
                exit(-1);
            }

            TFile* zdcHighFile = new TFile(fileNameZdcHigh.Data(), "READ");
            if (zdcHighFile->IsZombie()) { 
                cout << "Unable to open " << fileNameZdcHigh.Data() <<". It is a zombie.\n";
                exit(-1);
            }
            
            for (Int_t iPhi = 0; iPhi <= (nPhiBins - 1); iPhi++)
            {
                TString canvasLabel = TString::Format("%s - %s: %i - phi: %s", species.Data(), q2OrMultString.Data(), iMult, phiLabels[iPhi].Data());
                doProjections(zdcLowFile, zdcHighFile, ratios, phiLabels[iPhi]);
                TCanvas* canvas = makeCFCanvas(ratios, canvasLabel.Data());

                // Check for first/last iteration and make pdf page
                firstIteration = isFirstIteration(iq2OrMult, iMult, iPhi);
                lastIteration = isLastIteration(iq2OrMult, (nMultBins - iMult), (nPhiBins - iPhi));

                if( firstIteration && !lastIteration ) {
                    outFileName = (outFileNameBase + "(");
                } else if ( !firstIteration && lastIteration) {
                    outFileName = (outFileNameBase + ")");
                } else {
                    outFileName = outFileNameBase;
                }

                canvas->Print(outFileName.Data());
                delete canvas;

            } // phi bins

            zdcLowFile->Close();
            zdcHighFile->Close();

        } // mult/q2 bins
    } // multOrq2


}

Bool_t isFirstIteration(Int_t iq2OrMult, Int_t iMult, Int_t iPhi)
    { return  !iq2OrMult && !iMult && !iPhi; }

Bool_t isLastIteration(Int_t iq2OrMult, Int_t deltaMult, Int_t deltaPhi)
    { return (iq2OrMult == 1) && (deltaMult == 1) && (deltaPhi == 1); }

void doProjections(TFile* lowFile, TFile* highFile, TH1D** ratio, TString phiLabel)
{
    TString histoBaseLabels[4] = {"NumPiMinus_phi", "NumPiPlus_phi", "DenPiMinus_phi", "DenPiPlus_phi"};
    TString histoNewNames[4] = {"numLow", "denLow", "numHigh", "denHigh"};
    TH3D* numLow;
    TH3D* numHigh;
    TH3D* denLow;
    TH3D* denHigh;
    TString histName;
    TH1D* tempProjections[3];
    Float_t projRange = 0.011;
    // TH1::SetDefaultSumw2;


    // --- Get Histos from zdcLow --- //
    histName = TString::Format("%s%s_kt4",histoBaseLabels[0].Data(), phiLabel.Data());
    numLow = (TH3D*)lowFile->Get(histName.Data());
    numLow->SetNameTitle(histoNewNames[0], histoNewNames[0]);

    histName = TString::Format("%s%s_kt4",histoBaseLabels[1].Data(), phiLabel.Data());
    numLow->Add( (TH3D*)lowFile->Get(histName.Data()) );

    histName = TString::Format("%s%s_kt4",histoBaseLabels[2].Data(), phiLabel.Data());
    denLow = (TH3D*)lowFile->Get(histName.Data());
    denLow->SetNameTitle(histoNewNames[1], histoNewNames[1]);

    histName = TString::Format("%s%s_kt4",histoBaseLabels[3].Data(), phiLabel.Data());
    denLow->Add( (TH3D*)lowFile->Get(histName.Data()) );

    // --- Get Histos from zdcHigh --- //
    histName = TString::Format("%s%s_kt4",histoBaseLabels[0].Data(), phiLabel.Data());
    numHigh = (TH3D*)highFile->Get(histName.Data());
    numHigh->SetNameTitle(histoNewNames[2], histoNewNames[2]);

    histName = TString::Format("%s%s_kt4",histoBaseLabels[1].Data(), phiLabel.Data());
    numHigh->Add( (TH3D*)highFile->Get(histName.Data()) );

    histName = TString::Format("%s%s_kt4",histoBaseLabels[2].Data(), phiLabel.Data());
    denHigh = (TH3D*)highFile->Get(histName.Data());
    denHigh->SetNameTitle(histoNewNames[3], histoNewNames[3]);

    histName = TString::Format("%s%s_kt4",histoBaseLabels[3].Data(), phiLabel.Data());
    denHigh->Add( (TH3D*)highFile->Get(histName.Data()) );

    // --- Noralize Histograms --- //
    numLow->Scale((Float_t) (1./numLow->GetEntries() ));
    numHigh->Scale((Float_t) (1./numHigh->GetEntries() ));
    denLow->Scale((Float_t) (1./denLow->GetEntries() ));
    denHigh->Scale((Float_t) (1./denHigh->GetEntries() ));

    // // --- Set ratio equal to numLow projections --- //
    projectHistogram(numLow, ratio, projRange);

    // --- Divide by numHigh --- //
    projectHistogram(numHigh, tempProjections, projRange);
    for (Int_t i = 0; i <= 2; i++)
    { ratio[i]->Divide(tempProjections[i]); }

    // --- Divide by denLow --- //
    projectHistogram(denLow, tempProjections, projRange);
    for (Int_t i = 0; i <= 2; i++)
    { ratio[i]->Divide(tempProjections[i]); }

    // --- Multiply by denHigh --- //
    projectHistogram(denHigh, tempProjections, projRange);
    for (Int_t i = 0; i <= 2; i++)
    { ratio[i]->Multiply(tempProjections[i]); }

    // --- For making denominator ratios --- //
    // projectHistogram(denLow, ratio, projRange);
    // projectHistogram(denHigh, tempProjections, projRange);
    // for (Int_t i = 0; i <= 2; i++)
    // { ratio[i]->Divide(tempProjections[i]); }
}

void projectHistogram(TH3D* fullHist, TH1D** projHist, Float_t projRange)
{
    Double_t fullMin = fullHist->GetXaxis()->GetXmin();
    Double_t fullMax = fullHist->GetXaxis()->GetXmax();

    // Do x-projection
    fullHist->SetAxisRange(fullMin, fullMax, "x");
    fullHist->SetAxisRange(-1.*projRange, projRange, "y");
    fullHist->SetAxisRange(-1.*projRange, projRange, "z");
    projHist[0] = (TH1D*)fullHist->Project3D("xeNOFNUF");

    // Do y-projection
    fullHist->SetAxisRange(fullMin, fullMax, "y");
    fullHist->SetAxisRange(-1.*projRange, projRange, "x");
    fullHist->SetAxisRange(-1.*projRange, projRange, "z");
    projHist[1] = (TH1D*)fullHist->Project3D("yeNOFNUF");

    // Do z-projection
    fullHist->SetAxisRange(fullMin, fullMax, "z");
    fullHist->SetAxisRange(-1.*projRange, projRange, "x");
    fullHist->SetAxisRange(-1.*projRange, projRange, "y");
    projHist[2] = (TH1D*)fullHist->Project3D("zeNOFNUF");

}

TCanvas* makeCFCanvas(TH1D** ratios, TString canvasLabel)
{

    TPaveText* paveText = new TPaveText(0.1,0.9,0.9,0.95,"NDC");
    paveText->AddText(canvasLabel.Data());
    paveText->SetTextSize(.04);

    Float_t xLow = -.18, xHigh = .18;
    Float_t yLow = 0.95, yHigh = 1.05;
    Int_t canvasSize[2] = {1000, 700};


    const TString paramNames[3] = {"qOut", "qSide", "qLong"};

    TCanvas* canvas = new TCanvas(canvasLabel.Data(), canvasLabel.Data(), canvasSize[0], canvasSize[1]);
    gStyle->SetOptStat("");
    canvas->Divide(3,1,0,0);
    canvas->cd(3)->SetRightMargin(0.01);

    // Read the TH1D's 
    for (Int_t i = 0; i <= 2; ++i)
    {
        ratios[i]->SetMarkerStyle(kFullStar);
        ratios[i]->SetMarkerSize(2);
        ratios[i]->GetXaxis()->SetRangeUser(xLow, xHigh);
        ratios[i]->GetYaxis()->SetRangeUser(yLow, yHigh);
        ratios[i]->GetXaxis()->SetNdivisions(105,1);
        ratios[i]->GetYaxis()->SetNdivisions(105,1);
        ratios[i]->SetTitle(paramNames[i]);
    }

    // Put the TH1D's on the canvas
    for (Int_t i = 0; i <= 2; i++)
    {
        canvas->cd(i+1);
        ratios[i]->Draw();

    }

    canvas->cd(1);
    paveText->Draw();

    return canvas;

}
