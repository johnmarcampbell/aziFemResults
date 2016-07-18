void getEpsilon(const TGraphErrors* Rside, Double_t* epsilon, Double_t* errEpsilon);
TCanvas* makeRadiiCanvas(TDirectory* llDirectory, TDirectory* csDirectory, TGraphErrors** params, TString canvasName, TString paveText);
TCanvas* makeF0Canvas(TGraphErrors** graph, TString canvasName, TString paveText);
TCanvas* makeF2Canvas(TGraphErrors** graph, TString canvasName, TString paveText);
void getFComp(TGraphErrors* graph, TGraphErrors* compSin, TGraphErrors* compCos, const Double_t offset = 0);
void getFCompFit(TGraphErrors* graph, TGraphErrors* compSin, TGraphErrors* compCos, const Double_t offset = 0);

void compareRadii(
    const TString inFileNameLL = "./200GeVData/AuAu200GeVFitResults_WithRes.root",
    const TString inFileNameCS = "./200GeVData/AuAu200GeVFitResults_chiSquare.root",
    const TString outFileName = "./200GeVData/AuAu200Eps_WithRes.root"
)
{
    TFile* inFileLL = new TFile(inFileNameLL.Data(), "READ");
    TFile* inFileCS = new TFile(inFileNameCS.Data(), "READ");
    TFile* outFile = new TFile(outFileName.Data(), "RECREATE");
    TCanvas* radiiCanvas;
    TCanvas* f0Canvas;
    TCanvas* f2Canvas;
    TGraphErrors* params[5][6];
    TGraphErrors* fCompsCos[5][6];
    TGraphErrors* fCompsSin[5][6];
    TGraphErrors* f0[6];
    TGraphErrors* f2[6];
    TGraphErrors* eps[2];
    TGraphErrors* Rl[2];

    TString base = "auau";
    // TString base = "uu";

    for (Int_t i = 0; i <= 5; ++i)
    {
       f0[i] = new TGraphErrors(5);
       f2[i] = new TGraphErrors(5);
       eps[0] = new TGraphErrors(5);
       eps[1] = new TGraphErrors(5);
       Rl[0] = new TGraphErrors(5);
       Rl[1] = new TGraphErrors(5);

       for (Int_t j = 0; j <= 4; ++j)
       {
          fCompsSin[j][i] = new TGraphErrors(2); 
          fCompsCos[j][i] = new TGraphErrors(2); 
       } 
    }

    // const Int_t nZdcBins = 1;
    const Int_t nZdcBins = 2;
    const Int_t nQ2MultBins = 5;

    for (Int_t iZdc = 0; iZdc <= (nZdcBins - 1); iZdc++)
    {
        
        for (Int_t iq2Mult = 1; iq2Mult <= 1; iq2Mult++)
        // for (Int_t iq2Mult = 0; iq2Mult <= 1; iq2Mult++)
        {
            TString epsName, RlName;

            for (Int_t i = 0; i <= (nQ2MultBins - 1); ++i)
            {
                TString dirName, paveText;
                if(!iq2Mult) {
                    epsName = base.Data() + TString::Format("EpsF_Zdc_%d_q2",iZdc);
                    RlName = base.Data() + TString::Format("RLong_Zdc_%d_q2",iZdc);
                    dirName = TString::Format("zdc_%d/q2_%d/piCombined", iZdc, i);
                    paveText = TString::Format("Zdc - %d, q2: %d", iZdc, i+1);
                } else {
                    epsName = base.Data() + TString::Format("EpsF_Zdc_%d_mult",iZdc);
                    RlName = base.Data() + TString::Format("RLong_Zdc_%d_mult",iZdc);
                    dirName = TString::Format("zdc_%d/mult_%d/piCombined", iZdc, i);
                    paveText = TString::Format("Zdc - %d, Mult: %d", iZdc, i+1);
                }

                TDirectory* llDirectory = inFileLL->GetDirectory(dirName.Data());
                TDirectory* csDirectory = inFileCS->GetDirectory(dirName.Data());

                if(llDirectory && csDirectory)
                {
                    radiiCanvas = makeRadiiCanvas(llDirectory, csDirectory, params[i], dirName.Data(), paveText);

                    for (Int_t iParams = 0; iParams <= 5; ++iParams)
                    {
                        // getFCompFit(params[i][iParams],fCompsSin[i][iParams],fCompsCos[i][iParams]);
                        getFComp(params[i][iParams],fCompsSin[i][iParams],fCompsCos[i][iParams]);
                        Double_t* cos = fCompsCos[i][iParams]->GetY();
                        Double_t* sin = fCompsSin[i][iParams]->GetY();
                        f0[iParams]->SetPoint(i, i+1, cos[0]);
                        f0[iParams]->SetPointError(i, 0, fCompsCos[i][iParams]->GetErrorY(0));
                        f0[iParams]->SetTitle(params[i][iParams]->GetTitle());
                        f2[iParams]->SetPoint(i, i+1, cos[1]);
                        f2[iParams]->SetPointError(i, 0, fCompsCos[i][iParams]->GetErrorY(1));
                        f2[iParams]->SetTitle(params[i][iParams]->GetTitle());

                        // Do epsF
                        if(iParams == 1)
                        {
                            Double_t epsValue = cos[1]/cos[0];
                            Double_t meanErr = fCompsCos[i][0]->GetErrorY(0);
                            Double_t ampErr = fCompsCos[i][1]->GetErrorY(1);
                            Double_t epsErr = epsValue * sqrt((meanErr/cos[0])*(meanErr/cos[0]) + (ampErr/cos[1])*(ampErr/cos[1]));
                            eps[iq2Mult]->SetPoint(i,i+1,epsValue);
                            eps[iq2Mult]->SetPointError(i,0,fabs(epsErr));
                            eps[iq2Mult]->SetNameTitle(epsName.Data(),epsName.Data());
                        }

                        // Do Rl
                        if(iParams == 2)
                        {
                            Rl[iq2Mult]->SetPoint(i,i+1,cos[0]);
                            Double_t meanErr = fCompsCos[i][0]->GetErrorY(0);
                            Rl[iq2Mult]->SetPointError(i,0,fabs(meanErr));
                            Rl[iq2Mult]->SetNameTitle(RlName.Data(),RlName.Data());
                        }
                        
                    }

                    eps[iq2Mult]->SetNameTitle(epsName.Data(),epsName.Data());

                    TString radiiCanvasName;
                    if( (iZdc==0) && (iq2Mult==1) && (i==0) ) {
                        radiiCanvasName = TString::Format("%sRadii.pdf(", base.Data());
                    } else if ( (iZdc==(nZdcBins - 1)) && (iq2Mult==1) && (i==(nQ2MultBins - 1)) ) {
                        radiiCanvasName = TString::Format("%sRadii.pdf)", base.Data());
                    } else {
                        radiiCanvasName = TString::Format("%sRadii.pdf", base.Data());
                    }

                    radiiCanvas->Print(radiiCanvasName.Data());
                    delete radiiCanvas;
                } // if llDirectory exists
            } // loop over q2/mult bins

            outFile->cd();
            TString f2CanvasName, f0CanvasName;
            if(!iq2Mult) {
                f0CanvasName = TString::Format("%s - f0 - Zdc: %d - q2",base.Data(), iZdc);
                f2CanvasName = TString::Format("%s - f2 - Zdc: %d - q2",base.Data(), iZdc);
            } else {
                f0CanvasName = TString::Format("%s - f0 - Zdc: %d - mult",base.Data(), iZdc);
                f2CanvasName = TString::Format("%s - f2 - Zdc: %d - mult",base.Data(), iZdc);
            }

            f0Canvas = makeF0Canvas(f0,f0CanvasName, f0CanvasName);
            f2Canvas = makeF2Canvas(f2,f2CanvasName, f2CanvasName);
            TString f0Name, f2Name;
            if( (iZdc==0) && (iq2Mult==0) ) {
                f0Name = TString::Format("%sf0.pdf(",base.Data());
                f2Name = TString::Format("%sf2.pdf(",base.Data());
            } else if ( (iZdc==(nZdcBins - 1)) && (iq2Mult==1) ) {
                f0Name = TString::Format("%sf0.pdf)",base.Data());
                f2Name = TString::Format("%sf2.pdf)",base.Data());
            } else { 
                f0Name = TString::Format("%sf0.pdf",base.Data());
                f2Name = TString::Format("%sf2.pdf",base.Data());
            }
            f0Canvas->Print(f0Name.Data());
            f2Canvas->Print(f2Name.Data());

            eps[iq2Mult]->Write();
            f0[2]->SetNameTitle(RlName.Data(),RlName.Data());
            f0[2]->Write();
            delete f0Canvas;
            delete f2Canvas;
        } // loop over q2 or mult
    } // loop over zdc Bins
}

TCanvas* makeRadiiCanvas(TDirectory* llDirectory, TDirectory* csDirectory, TGraphErrors** params, TString canvasName, TString paveText)
{

    llDirectory->cd();

    TPaveText *pt = new TPaveText(0,0.9,0.4,1,"NDC");
    pt->AddText(paveText.Data());
    pt->SetTextSize(.04);

    const TString paramNames[6] = {"Ro^2", "Rs^2", "Rl^2", "Ros^2", "lam", "norm"};
    const TString llNames[6] = {"Ro^2ll", "Rs^2ll", "Rl^2ll", "Ros^2ll", "lamll", "normll"};
    const TString csNames[6] = {"Ro^2cs", "Rs^2cs", "Rl^2cs", "Ros^2cs", "lamcs", "normcs"};
    TCanvas* canvas = new TCanvas(canvasName.Data(),canvasName.Data());
    canvas->Divide(3,2);

    Float_t yAxisRange[6][2] = {{24, 36}, {17, 29}, {27, 39}, {-4, 4}, {0.4, 0.5}, {0.10, .14}};
    // Float_t yAxisRange[6][2] = {{24, 36}, {17, 29}, {27, 39}, {-4, 4}, {0.4, 0.5}, {.115, .122}};

    // --- Get/Draw log-likelihood fits --- //
    for (Int_t i = 0; i <= 5; ++i)
    {
        params[i] = (TGraphErrors*)gDirectory->Get(paramNames[i].Data());
        params[i]->SetName(llNames[i].Data());
        params[i]->SetMarkerStyle(kFullStar);
        params[i]->SetMarkerSize(2);
        params[i]->GetXaxis()->SetRangeUser(-20,200);
        params[i]->GetYaxis()->SetRangeUser(yAxisRange[i][0], yAxisRange[i][1]);

        canvas->cd(i+1);
        params[i]->Draw("alp");
    }

    // --- Get/Draw chi-square fits --- //
    csDirectory->cd();
    for (Int_t i = 0; i <= 5; ++i)
    {
        params[i] = (TGraphErrors*)gDirectory->Get(paramNames[i].Data());
        params[i]->SetName(csNames[i].Data());
        params[i]->SetMarkerStyle(kFullStar);
        params[i]->SetMarkerColor(kBlue);
        params[i]->SetLineColor(kBlue);
        params[i]->SetMarkerSize(2);

        canvas->cd(i+1);
        params[i]->Draw("same");
    }

    canvas->cd(1);
    pt->Draw();

    return canvas;

}

TCanvas* makeF0Canvas(TGraphErrors** graph, TString canvasName, TString paveText)
{

    TCanvas* canvas = new TCanvas(canvasName.Data(),canvasName.Data());
    canvas->Divide(3,2);

    TPaveText *pt = new TPaveText(0,0.9,0.4,1,"NDC");
    pt->AddText(paveText.Data());
    pt->SetTextSize(.04);


    Float_t yAxisRange[6][2] = {{26, 31}, {21, 26}, {29, 34}, {-0.5, 0.5}, {0.42, 0.45}, {.117, .122}};

    for (Int_t i = 0; i <= 5; ++i)
    {
        canvas->cd(i+1);
        graph[i]->Draw("ap");
        graph[i]->SetMarkerStyle(kFullStar);
        graph[i]->SetMarkerSize(2);
//        graph[i]->GetYaxis()->SetRangeUser(yAxisRange[i][0], yAxisRange[i][1]);

    }

    canvas->cd(1);
    pt->Draw();

    return canvas;

}

TCanvas* makeF2Canvas(TGraphErrors** graph, TString canvasName, TString paveText)
{

    TCanvas* canvas = new TCanvas(canvasName.Data(),canvasName.Data());
    canvas->Divide(3,2);
    TPaveText *pt = new TPaveText(0,0.9,0.4,1,"NDC");
    pt->AddText(paveText.Data());
    pt->SetTextSize(.04);


    Float_t yAxisRange[6][2] = {{0, 4.5}, {0, 4.5}, {0, 4.5}, {0, 1.3}, {0, 0.02}, {2e-4, 7e-4}};

    for (Int_t i = 0; i <= 5; ++i)
    {
        canvas->cd(i+1);
        graph[i]->Draw("ap");
        graph[i]->SetMarkerStyle(kFullStar);
        graph[i]->SetMarkerSize(2);
//        graph[i]->GetYaxis()->SetRangeUser(yAxisRange[i][0], yAxisRange[i][1]);

    }

    canvas->cd(1);
    pt->Draw();

    return canvas;

}

void getFComp(TGraphErrors* graph, TGraphErrors* compSin, TGraphErrors* compCos, const Double_t offset)
{
	
	Double_t* RX = graph->GetX();
	Double_t* RY = graph->GetY();
	Double_t tempSin = 0, tempCos = 0, tempSinErr = 0, tempCosErr = 0, tempExp = 0, tempExpErr = 0;
	
	for(int i = 0; i <= 1; i++) //Loop over harmonic orders, i
	{
		tempSin = 0, tempCos = 0, tempSinErr = 0, tempCosErr = 0;

		for(int j = 0; j <= 7; j++) // Loop over azimuthal bins, j
		{

			tempSin += RY[j] * sin(2 * i * RX[j] / 57.30);
			tempSinErr += graph->GetErrorY(j)*graph->GetErrorY(j) * sin(2 * i * RX[j] / 57.30)*sin(2 * i * RX[j] / 57.30);
			tempCos += RY[j] * cos(2 * i * RX[j] / 57.30);			
			tempCosErr += graph->GetErrorY(j)*graph->GetErrorY(j) * cos(2 * i * RX[j] / 57.30)*cos(2 * i * RX[j] / 57.30);
			//if(j==0){tempSin += RY[j] * sin(2 * i * RX[j] / 57.30); tempCos += RY[j] * cos(2 * i * RX[j] / 57.30);}

		}
		
		tempSinErr = sqrt(tempSinErr);
		tempCosErr = sqrt(tempCosErr);
		
		if(i==0){tempCos /= 8.; tempSin /= 8.;} 
		else {tempCos /= 4.; tempSin /= 4.;}
		
		tempSinErr /=8.;
		tempCosErr /=8.;

		tempExp = sqrt(tempSin*tempSin + tempCos*tempCos);
		tempExpErr = pow(tempSin * tempSinErr, 2) + pow(tempCos * tempCosErr, 2);
		tempExpErr /= tempSin*tempSin + tempCos*tempCos;
		tempExpErr = sqrt(tempExpErr);

		compSin->SetPoint(i, 2*i + offset, tempSin);
		compCos->SetPoint(i, 2*i + offset, tempCos);
		compSin->SetPointError(i, 0, tempSinErr);
		compCos->SetPointError(i, 0, tempCosErr);
         
	}
}

void getEpsilon(const TGraphErrors* Rside, Double_t* epsilon, Double_t* errEpsilon)
{

	Double_t temp1 = 0;
	Double_t Rs_0 = 0, Rs_2 = 0, errRs_0 = 0, errRs_2 = 0;
	
	Rside->GetPoint(0,temp1,Rs_0);
	Rside->GetPoint(1,temp1,Rs_2);
	errRs_0 = Rside->GetErrorY(0);
	errRs_2 = Rside->GetErrorY(1);
	
	*epsilon = Rs_2 / Rs_0;
	*errEpsilon = *epsilon * sqrt(pow(errRs_0/Rs_0,2) + pow(errRs_2/Rs_2,2));

    cout << *errEpsilon << endl;

}

void getFCompFit(TGraphErrors* graph, TGraphErrors* compSin, TGraphErrors* compCos, const Double_t offset)
{
	
	Double_t* RX = graph->GetX();
	Double_t* RY = graph->GetY();
	Double_t tempSin = 0, tempCos = 0, tempSinErr = 0, tempCosErr = 0, tempExp = 0, tempExpErr = 0;
    TF1* cosFit = new TF1("cosFit","[0] + [1]*cos(2 * x / 57.30)", 0, 180);
    TF1* sinFit = new TF1("sinFit","[0] + [1]*sin(2 * x / 57.30)", 0, 180);

    graph->Fit(sinFit, "Q");
    graph->Fit(cosFit, "Q");
    
	
	for(int i = 0; i <= 1; i++) //Loop over harmonic orders, i
	{

		compSin->SetPoint(i, 2*i + offset, sinFit->GetParameter(i));
		compCos->SetPoint(i, 2*i + offset, cosFit->GetParameter(i));
		compSin->SetPointError(i, 0, sinFit->GetParError(i));
		compCos->SetPointError(i, 0, cosFit->GetParError(i));
         
	}
}

