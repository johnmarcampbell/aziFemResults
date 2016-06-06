void  makePdf(TString inFileName, TString outFileNameBase, TString directoryFormat);
TCanvas* makeCFCanvas(TFile* inFile, TString dirName, TString paveText);
Bool_t isFirstIteration(Int_t iZdc, Int_t iMult, Int_t iPhi);
Bool_t isLastIteration(Int_t deltaZdc, Int_t deltaMult, Int_t deltaPhi);

void drawCFs()
{

    TString inFileName[2] = {"./193GeVData/UU193GeVFitResults.root", "./200GeVData/AuAu200GeVFitResults.root"};
    TString outFileNameBase[4] = {"CFUUMultOut.pdf", "CFUUq2Out.pdf", "CFAuAuMultOut.pdf", "CFAuAuq2Out.pdf"};
    TString directoryFormat[2] = {"zdc_%d/mult_%d/piCombined/phi_%d", "zdc_%d/q2_%d/piCombined/phi_%d"};

    makePdf(inFileName[0], outFileNameBase[0], directoryFormat[0]);
    makePdf(inFileName[0], outFileNameBase[1], directoryFormat[1]);
    makePdf(inFileName[1], outFileNameBase[2], directoryFormat[0]);
    makePdf(inFileName[1], outFileNameBase[3], directoryFormat[1]);
}

void  makePdf(TString inFileName, TString outFileNameBase, TString directoryFormat)
{

    TFile* inFile = new TFile(inFileName.Data(), "READ");
    if (inFile->IsZombie()) { 
        cout << "Unable to open " << inFileName.Data() <<". It is a zombie.\n";
        exit(-1);
    }

    TCanvas* canvas = 0;

    TString outFileName;
    // Int_t nZdcBins = 1, nMultBins = 1, nPhiBins = 1;
    Int_t nZdcBins = 2, nMultBins = 5, nPhiBins = 8;
    Bool_t firstIteration = kTRUE, lastIteration = kFALSE;

    for (Int_t iZdc = 0; iZdc <= (nZdcBins - 1); ++iZdc)
    {
        for (Int_t iMult = 0; iMult <= (nMultBins - 1); iMult++)
        {
            for (Int_t iPhi = 0; iPhi <= (nPhiBins - 1); iPhi++)
            {
                // TString dirName = TString::Format("zdc_%d/mult_%d/piCombined/phi_%d", iZdc, iMult, iPhi);
                // TString paveText = TString::Format("Zdc: %d - Mult: %d - Phi: %d", iZdc, iMult, iPhi);
                TString dirName = TString::Format(directoryFormat.Data(), iZdc, iMult, iPhi);
                TString paveText = TString::Format("Zdc: %d - q2: %d - Phi: %d", iZdc, iMult, iPhi);
                if(canvas) {delete canvas;}
                canvas = makeCFCanvas(inFile, dirName.Data(), paveText);

                // Check for first/last iteration and make pdf page
                firstIteration = isFirstIteration(iZdc, iMult, iPhi);
                lastIteration = isLastIteration((nZdcBins - iZdc), (nMultBins - iMult), (nPhiBins - iPhi));

                if( firstIteration && !lastIteration ) {
                    outFileName = (outFileNameBase + "(");
                } else if ( !firstIteration && lastIteration) {
                    outFileName = (outFileNameBase + ")");
                } else {
                    outFileName = outFileNameBase;
                }

                canvas->Print(outFileName.Data());
                
            }
        }
    }

    if(canvas) {delete canvas;}

}
TCanvas* makeCFCanvas(TFile* inFile, TString dirName, TString paveText)
{

    inFile->cd(dirName.Data());

    TPaveText* pt = new TPaveText(0,0.9,0.4,1,"NDC");
    pt->AddText(paveText.Data());
    pt->SetTextSize(.04);

    Float_t xLow = -.18, xHigh = .18;
    Float_t yLow = 0.95, yHigh = 1.5;
    Float_t yLowRatio = 0.95, yHighRatio = 1.1;
    Int_t canvasSize[2] = {1000, 700};


    TH1D* params[6];
    TH1D* paramRatios[3];
    const TString paramNames[6] = {"qOut", "qSide", "qLong", "qOutFit", "qSideFit", "qLongFit"};
    const TString paramRatioNames[3] = {"qOutRatio", "qSideRatio", "qLongRatio"};

    TCanvas* canvas = new TCanvas(paveText.Data(), paveText.Data(), canvasSize[0], canvasSize[1]);
    canvas->Divide(3,2,0,0);
    canvas->cd(3)->SetRightMargin(0.01);
    canvas->cd(6)->SetRightMargin(0.01);

    // Read the TH1D's 
    for (Int_t i = 0; i <= 5; ++i)
    {
        params[i] = (TH1D*)gDirectory->Get(paramNames[i].Data());
        params[i]->SetMarkerStyle(kFullStar);
        params[i]->SetMarkerSize(2);
        if(i >= 3) { params[i]->SetLineWidth(4); }
        params[i]->GetXaxis()->SetRangeUser(xLow, xHigh);
        params[i]->GetXaxis()->SetNdivisions(105,1);
        params[i]->GetYaxis()->SetNdivisions(105,1);
    }

    for (Int_t i = 0; i <= 2; i++)
    {
        paramRatios[i] = (TH1D*)params[i+3]->Clone(paramRatioNames[i]);
        paramRatios[i]->SetTitle(paramRatioNames[i]);
        paramRatios[i]->SetMarkerColor(kBlue);
        paramRatios[i]->SetLineColor(kBlue);
        paramRatios[i]->SetLineWidth(2);
        paramRatios[i]->Divide(params[i]);
        paramRatios[i]->GetXaxis()->SetRangeUser(xLow, xHigh);
        paramRatios[i]->GetYaxis()->SetRangeUser(yLowRatio, yHighRatio);
    }

    // Put the TH1D's on the canvas
    for (Int_t i = 0; i <= 2; i++)
    {
        canvas->cd(i+1);
        params[i]->Draw();
        params[i+3]->Draw("lhistsame");

        canvas->cd(i+4);
        paramRatios[i]->Draw();

    }

    canvas->cd(1);
    pt->Draw();

    return canvas;

}

Bool_t isFirstIteration(Int_t iZdc, Int_t iMult, Int_t iPhi)
    { return  !iZdc && !iMult && !iPhi; }

Bool_t isLastIteration(Int_t deltaZdc, Int_t deltaMult, Int_t deltaPhi)
    { return (deltaZdc == 1) && (deltaMult == 1) && (deltaPhi == 1); }

