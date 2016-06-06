#include "fit.C"
TH3D* histCopy(TH3* hist, TString nameTitle);
void doProjections(TH3D* num, TH3D* fitNum, TH3D* den, TH1D** projCF, TH1D** projFitCF, Double_t normalization, const Double_t projRange);
void projectHistogram(TH3D* fullHist, TH1D** projHist, Float_t projRange);
void readData(TH3D* &num, TH3D* &den, TH3D* &coul, const TString inFile, const Int_t pm, const Int_t ktBin, const Int_t phiBin);
TDirectory*  makeDirectory(TFile* theFile, const Int_t zdcBin, const Int_t q2MultBin, const Int_t q2Mult, const Int_t pm, const Int_t phiBin);
void writeTGraphs(TMinuit* minuit, const Int_t binNumber, const Int_t xValue);

void fitManager(
    const TString inputfile = "bin0.root",  //input file
    const TString outFileName = "testOut.root",  //output file
    const Int_t zdcBin = 0,             // zdc bin
    const Int_t pm = 0,                 // -1 = piMinus, 0 = piCombined, 1 = piPlus
    const Int_t q2Mult = 0,             // 0 = q2, 1 = mult
    const Int_t q2MultBin = 0,          // q2MultBin
    const Int_t phiBin = 0,             // phiBin
    const Int_t ktBin = 4,              // 4 = Integrated kT
    const Double_t fitRange = 0.149,    // Do the fit from -fitRange to +fitRange
    const Double_t projRange = 0.03,    // Project from -projRange to +projRange
    const Double_t initNorm = 0.16,   // Initial Values for the fit
    const Double_t initLambda = 0.45, 
    const Double_t initRo = 27, 
    const Double_t initRs = 22, 
    const Double_t initRl = 32, 
    const Double_t initRos = 0, 
    const Double_t initRol = 0, 
    const Double_t initRsl = 0 
    )
{

    TH3D* num;
    TH3D* den;
    TH3D* coul;
    TH3D* fitNum;
    TGraph* contour = 0;

    readData(num,den,coul,inputfile.Data(),pm,ktBin,phiBin);
    coul->Divide(den);

    Int_t nFitPars = 8;
    Double_t fitPars[8] = {0.};
    Double_t fitParErrors[8] = {0.};
	TMinuit* minuit = new TMinuit(nFitPars);

    //Do the fit
    TStopwatch* fitTimer = new TStopwatch();
    fit(num, den, coul, nFitPars, fitRange, minuit, initNorm, initLambda, initRo, initRs, initRl, initRos, initRol, initRsl, &contour);
    cout << "\nFit took " << fitTimer->RealTime() << " seconds to finish.\n\n";

    for (Int_t i = 0; i <= 7; ++i) { minuit->GetParameter(i, fitPars[i], fitParErrors[i]); }

    //Project histograms
    fitNum = histCopy(den, "Theoretical Numerator");
    makeFitNumerator(fitNum, den, coul, fitPars);

    TH1D* projectedCF[3];
    TH1D* projectedFitCF[3];

    // Save correlation function projections
    TFile* outFile = new TFile(outFileName.Data(),"UPDATE");
    TDirectory* fitDirectory = makeDirectory(outFile, zdcBin, q2MultBin, q2Mult, pm, phiBin);
    fitDirectory->cd();
    doProjections(num, fitNum, den, projectedCF, projectedFitCF, fitPars[0], projRange);
    minuit->Write("Minuit", 2);

    // Save fit parameters to TGraphs
    Int_t phiLabels[8] = {0,22,45,67,90,112,135,157};
    gDirectory->cd("..");
    writeTGraphs(minuit, phiBin, phiLabels[phiBin]);
    cout << contour << endl;
    TCanvas* test = new TCanvas();
    contour->Draw();
    test->Print("test.pdf");

    outFile->Close();
    delete minuit;
    delete fitTimer;
}

TH3D* histCopy(TH3* hist, TString nameTitle)
{

    Int_t nBins = hist->GetNbinsX();   
    Double_t lo = hist->GetXaxis()->GetXmin();
    Double_t hi = hist->GetXaxis()->GetXmax();

    TH3D* rHist = new TH3D(nameTitle.Data(),nameTitle.Data(),nBins,lo,hi,nBins,lo,hi,nBins,lo,hi);

    return rHist;

}

void readData(TH3D* &num, TH3D* &den, TH3D* &coul, const TString inFile, const Int_t pm, const Int_t ktBin, const Int_t phiBin)
{

    // Create histogram names
    Int_t phiLabels[8] = {0,22,45,67,90,112,135,157};
//    Int_t phiLabels[4] = {0,45,90,135};
    Float_t ktLabels[4] = {0.22,0.33,0.42,0.52};
    TString histNames[6] = { "NumPiPlus" , "DenPiPlus" , "CoulPiPlus" , "NumPiMinus" , "DenPiMinus" , "CoulPiMinus" };


    for(Int_t i = 0; i <= 5; i++)
    {
        if(phiBin >= 0) {histNames[i] += "_phi"; histNames[i] += phiLabels[phiBin];}

        histNames[i] += "_kt";
        histNames[i] += ktBin;

    }

    // Get histograms from input file
	TFile* file = new TFile(inFile.Data());
    TH3S* pNum = (TH3S*)file->Get(histNames[0].Data());
    TH3S* pDen = (TH3S*)file->Get(histNames[1].Data());
    TH3D* pCoul = (TH3D*)file->Get(histNames[2].Data());
    TH3S* mNum = (TH3S*)file->Get(histNames[3].Data());
    TH3S* mDen = (TH3S*)file->Get(histNames[4].Data());
    TH3D* mCoul = (TH3D*)file->Get(histNames[5].Data());

    num = histCopy(pNum, "Numerator");
    den = histCopy(pDen, "Denominator");
    coul = histCopy(pCoul, "Coulomb");

	if( pm != -1 )
	{
		num->Add(pNum);
		den->Add(pDen);
		coul->Add(pCoul);
	}

	if( pm != 1 )
	{
		num->Add(mNum);
		den->Add(mDen);
		coul->Add(mCoul);
	}

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

void doProjections(TH3D* num, TH3D* fitNum, TH3D* den, TH1D** projCF, TH1D** projFitCF, Double_t normalization, const Double_t projRange)
{
    TH1D* projectedNum[3];
    TH1D* projectedFitNum[3];
    TH1D* projectedDen[3];

    projectHistogram(num, projectedNum, projRange);
    projectHistogram(fitNum, projectedFitNum, projRange);
    projectHistogram(den, projectedDen, projRange);

    TString histNames[3] = {"qOut", "qSide", "qLong"};

    // check for zero normalization
    if( normalization == 0 ){
        cout << "doProjections(): Received a value of zero for the normalization." << endl;
        cout << "doProjections(): Setting normilzation = num->GetEntries() / den->GetEntries()" << endl;
        normalization = (Double_t) num->GetEntries() / den->GetEntries();
    }

    for (Int_t i = 0; i <= 2; ++i)
    {
        projCF[i] = new TH1D(*projectedNum[i]);
        projCF[i]->Divide(projectedDen[i]);
        projCF[i]->Scale(1.0/normalization);
        projCF[i]->SetNameTitle(histNames[i].Data(), histNames[i].Data());
        
        projFitCF[i] = new TH1D(*projectedFitNum[i]);
        projFitCF[i]->Divide(projectedDen[i]);
        histNames[i] += "Fit";
        projFitCF[i]->SetNameTitle(histNames[i].Data(), histNames[i].Data());

        // Some basic graphics options
        projCF[i]->SetMarkerStyle(kFullTriangleUp);
        projCF[i]->SetMarkerColor(kBlue);
        projCF[i]->SetMaximum(1.4);
        projCF[i]->SetMinimum(0.9);
        
        projFitCF[i]->SetMarkerColor(kRed);
        projFitCF[i]->SetLineColor(kRed);
        projFitCF[i]->SetLineWidth(2);
        projFitCF[i]->SetMaximum(1.4);
        projFitCF[i]->SetMinimum(0.9);

        //Write the histograms
        projCF[i]->Write(projCF[i]->GetName(), 2);
        projFitCF[i]->Write(projFitCF[i]->GetName(), 2);

        delete projCF[i];
        delete projFitCF[i];
        
    }

}

TDirectory*  makeDirectory(TFile* theFile, const Int_t zdcBin, const Int_t q2MultBin, const Int_t q2Mult, const Int_t pm, const Int_t phiBin)
{

    TString zdcDir = "zdc_";
    zdcDir += zdcBin;
    if(!gDirectory->GetDirectory(zdcDir.Data())) {gDirectory->mkdir(zdcDir.Data());}
    gDirectory->cd(zdcDir.Data());

    TString q2MultDir = "";
    if(!q2Mult) { 
        q2MultDir += "q2_";
    } else {
        q2MultDir += "mult_";
    }
    q2MultDir += q2MultBin;
    if(!gDirectory->GetDirectory(q2MultDir.Data())) {gDirectory->mkdir(q2MultDir.Data());}
    gDirectory->cd(q2MultDir.Data());

    TString pmDir = ""; 
    if( pm == -1 ) {pmDir += "piMinus"; }
    else if( pm == 0 ) {pmDir += "piCombined"; }
    else if( pm == 1 ) {pmDir += "piPlus"; }
    if(!gDirectory->GetDirectory(pmDir.Data())) {gDirectory->mkdir(pmDir.Data());}
    gDirectory->cd(pmDir.Data());

    TString phiDir = "phi_";
    phiDir += phiBin;
    if(!gDirectory->GetDirectory(phiDir.Data())) {gDirectory->mkdir(phiDir.Data());}

    return gDirectory->GetDirectory(phiDir.Data());

}

void writeTGraphs(TMinuit* minuit, const Int_t binNumber, const Int_t xValue)
{

    const Int_t nPhiBins = 8;
    Double_t fitParameters[8] = {0.};
    Double_t fitParameterErrors[8] = {0.};

    for (Int_t i = 0; i <= 7; ++i) { minuit->GetParameter(i, fitParameters[i], fitParameterErrors[i]); }

    TGraphErrors* graphs[13];
    TString graphNames[13] = {"norm", "lam", "Ro^2", "Rs^2", "Rl^2", "Ros^2", "Rol^2", "Rsl^2", "Ro", "Rs", "Rl", "RoOverRs", "Ro2MinusRs2"};
    TString graphTitles[13] = {"Normalization", "#lambda", "R_{o}^{2}", "R_{s}^{2}", "R_{l}^{2}", "R_{os}^{2}", "R_{ol}^{2}", "R_{sl}^{2}", "R_{o}", "R_{s}", "R_{l}", "R_{o} / R_{s}", "R_{o}^{2} - R_{s}^{2}"};
    for (Int_t i = 0; i <= 12; ++i)
    {

        graphs[i] = (TGraphErrors*)gDirectory->Get(graphNames[i].Data());

        // If the TGraphErrors's don't already exist, create them
        if(!graphs[i]) 
        {
            graphs[i] = new TGraphErrors(nPhiBins);
            graphs[i]->SetNameTitle(graphNames[i].Data(), graphTitles[i].Data());

        } else { // If they *do* exist, delete the old copies so we don't get a bunch of duplicates
            TString nameCycle = graphNames[i];
            nameCycle += ";1";
            gDirectory->Delete(nameCycle.Data());
        }
        
    }

    // The first 8 TGraphErrors correspond directly to fit parameters. Just copy them over
    for (Int_t i = 0; i <= 7; ++i)
    {
        graphs[i]->SetPoint(binNumber, xValue, fitParameters[i]);
        graphs[i]->SetPointError(binNumber, 0, fitParameterErrors[i]);
    }

    // For Ro/Rs/Rl, we just take the square roots of the values
    for (Int_t i = 8; i <= 10; ++i)
    {
        Double_t yValue = sqrt(fitParameters[i-6]);
        Double_t yValueError = fitParameterErrors[i-6] / (2*yValue);
        graphs[i]->SetPoint(binNumber, xValue, yValue);
        graphs[i]->SetPointError(binNumber, 0, yValueError);
    }

    // (Ro^2 - Rs^2) and (Ro/Rs) are a little trickier
    Double_t Ro2MinusRs2 = fitParameters[2] - fitParameters[3];
    Double_t Ro2MinusRs2Error = sqrt(fitParameters[2]*fitParameters[2] + fitParameters[3]*fitParameters[3]);
    graphs[11]->SetPoint(binNumber, xValue, Ro2MinusRs2);
    graphs[11]->SetPointError(binNumber, 0, Ro2MinusRs2Error);

    Double_t Ro = sqrt(fitParameters[2]);
    Double_t RoErr = fitParameterErrors[2] / (2*Ro);
    Double_t Rs = sqrt(fitParameters[3]);
    Double_t RsErr = fitParameterErrors[3] / (2*Rs);
    Double_t RoOverRs = Ro / Rs;
    Double_t RoOverRsErr = (Ro/Rs) * sqrt( (RoErr/Ro)*(RoErr/Ro) + (RsErr/Rs)*(RsErr/Rs) );
    graphs[12]->SetPoint(binNumber, xValue, RoOverRs);
    graphs[12]->SetPointError(binNumber, 0, RoOverRsErr);


    for (Int_t i = 0; i <= 12; ++i)
    {
        graphs[i]->Write(graphs[i]->GetName(), 2);
        delete graphs[i];
        
    }

}
