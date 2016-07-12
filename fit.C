const TH3D* hNum;
const TH3D* hDen;
const TH3D* hCoul;
TMinuit* tmFit;
Int_t loBin;
Int_t hiBin;

Double_t cfTheory(const Double_t* fitPars, Double_t qx, Double_t qy, Double_t qz, Double_t coul);
void logLikelihoodWrapper(int &npar, double *gin, double &f, double *par, int iflag);
void logLikelihood(double &f, double *par);
void makeFitNumerator(TH3D* fitNum, const TH3D* den, const TH3D* coul, const Double_t* fitPars);
TGraph* scanDistribution(Double_t* pars, const Double_t scanFraction, const Int_t nPar);
TGraph* scanNewDistribution(Double_t* pars, const Double_t scanFraction, const Int_t nPar);

void fit(
    const TH3D* pNum,
    const TH3D* pDen,
    const TH3D* pCoul,
    const Int_t nFitPars,
    const Double_t fitRange,
    TMinuit* pFit,
    const Double_t initNorm = 0.122, 
    const Double_t initLambda = 0.511, 
    const Double_t initRo = 38.1, 
    const Double_t initRs = 27.0, 
    const Double_t initRl = 49.9, 
    const Double_t initRos = -0.259, 
    const Double_t initRol = 0, 
    const Double_t initRsl = 0, 
    TGraph** contour = 0, 
    TGraph** dist = 0,
    TGraph** newDist = 0
)
{
    hNum = pNum;
    hDen = pDen;
    hCoul = pCoul;
    tmFit = pFit;

    //Determine loBin and hiBin
	Int_t nBin = hNum->GetNbinsX();
	loBin = hNum->GetXaxis()->FindBin(-1.*fitRange);
	hiBin = hNum->GetXaxis()->FindBin(fitRange);

    // Set up fit parameters
    Double_t initPars[8] = {initNorm, initLambda, initRo, initRs, initRl, initRos, initRol, initRsl};
    // Double_t initRange[8] = {0.1, 0, 3, 3, 3, 1, 0, 0};
    Double_t initRange[8] = {0.1, 0.1, 3, 3, 3, 1, 0, 0};
    TString parNames[8] = {"Normalization" , "Lambda", "RoutSquared", "RsideSquared", "RlongSquared", "RoutsideSquared", "RoutlongSquared", "RsidelongSquared"};
    tmFit->SetFCN(logLikelihoodWrapper);
    for (Int_t i = 0; i <= 7; ++i) { tmFit->DefineParameter(i, parNames[i].Data(), initPars[i], initRange[i], 0, 0); }

    Double_t arglist[8];
    Int_t ierflg = 1;
    Int_t fitLoop = 0;

    // Do the fit
	do
	{
        fitLoop++;
        cout << "***** Starting Fit loop " << fitLoop << endl;
        tmFit->mnexcm("MIGRAD",arglist,0,ierflg);
	}
	while ( (ierflg != 0) && (fitLoop != 50) );

    tmFit->mnexcm("MINOS",arglist,0,ierflg);

    // ---- Find Error Contours ---- //
    Int_t n = 0;
    Int_t nPoints = 10;
    for (Int_t i = 0; i <= 6; i++)
    {
        for (Int_t j = i+1; j <= 7; j++)
        {
            TString cloneTitle = TString::Format("%s_%s", parNames[i].Data(), parNames[j].Data());

            // Only make the contours if we actually tried to fit the values
            if(initRange[i] && initRange[j]) {
                // contour[n] = (TGraph*)(tmFit->Contour(nPoints,i,j))->Clone(cloneTitle.Data());
            }
            n++;
        }
    }

    // ---- Find f-distribution ---- //
    Double_t tempPars[8] = {0};
    Double_t tempParErrors[8] = {0};
    Double_t scanFraction = 0.20;

    for (Int_t i = 0; i <= 7; ++i) {
        // Only make the dist if we actually tried to fit the value
        if(initRange[i]) {
            for (Int_t i = 0; i <= 7; ++i) { tmFit->GetParameter(i, tempPars[i], tempParErrors[i]); }
            dist[i] = scanDistribution(tempPars, scanFraction, i);
            TString title = TString::Format("%s_Dist", parNames[i].Data());
            dist[i]->SetNameTitle(title.Data(), title.Data());

            for (Int_t i = 0; i <= 7; ++i) { tmFit->GetParameter(i, tempPars[i], tempParErrors[i]); }
            newDist[i] = scanNewDistribution(tempPars, scanFraction, i);
            TString newTitle = TString::Format("%s_DistNew", parNames[i].Data());
            newDist[i]->SetNameTitle(newTitle.Data(), newTitle.Data());
        }
    }
}


Double_t cfTheory(const Double_t* fitPars, Double_t qx, Double_t qy, Double_t qz, Double_t coul) 
{

    Double_t hBar = 0.1973269718; // in GeV*fm, from Particle Data Group
    Double_t arg = 0, theory = 0;

    // construct the argument of the exponential
    arg = qx*qx*fitPars[2]; // q^2_Out*R^2_Out
    arg += qy*qy*fitPars[3]; // q^2_Side*R^2_Side
    arg += qz*qz*fitPars[4]; // q^2_Long*R^2_Long
    arg += 2*qx*qy*fitPars[5]; // q_Out*q_Side*R^2_os
    arg += 2*qx*qz*fitPars[6]; // q_Out*q_Side*R^2_ol
    arg += 2*qy*qz*fitPars[7]; // q_Side*q_Long*R^2_sl
    arg /= -1.0*hBar*hBar;

    theory = (1 - fitPars[1]) + fitPars[1]*coul*(1 + exp(arg)); //UnNormalized CF

    return theory;

}

void logLikelihoodWrapper(int &npar, double *gin, double &f, double *par, int iflag) 
{
    logLikelihood(f, par);
}

void newLogLikelihood(double &f, double *par) 
{
	Float_t qx, qy, qz;
    Double_t n, d, c, t; // (n)umerator, (d)enominator, (c)oulomb, (t)heory

	f = 0.;
	
	for (Int_t x = loBin; x <= hiBin; x++) 
	{
        for (Int_t y = loBin; y <= hiBin; y++) 
		{
            for (Int_t z = loBin; z <= hiBin; z++) 
			{
                qx = hNum->GetXaxis()->GetBinCenter(x);
                qy = hNum->GetYaxis()->GetBinCenter(y);
				qz = hNum->GetZaxis()->GetBinCenter(z);

				n =  hNum->GetBinContent(x,y,z);
				d =  hDen->GetBinContent(x,y,z);
				c =  hCoul->GetBinContent(x,y,z);

				if((d > 0.0001) && (n > 0.0001)) 
				{
                    t = par[0] * cfTheory(par, qx, qy, qz, c);
					f  += 2. * ( (d+1)*log(t+1) - n*log( t/(t+1) ));
				}

			} // z bins
		} // y bins
	} // x bins

}

void logLikelihood(double &f, double *par) 
{
	Float_t qx, qy, qz;
    Double_t n, d, c, t; // (n)umerator, (d)enominator, (c)oulomb, (t)heory
    Double_t tempF = 0; // placeholder variable for f

	f = 0.;
	
	for (Int_t x = loBin; x <= hiBin; x++) 
	{
        for (Int_t y = loBin; y <= hiBin; y++) 
		{
            for (Int_t z = loBin; z <= hiBin; z++) 
			{
                qx = hNum->GetXaxis()->GetBinCenter(x);
                qy = hNum->GetYaxis()->GetBinCenter(y);
				qz = hNum->GetZaxis()->GetBinCenter(z);

				n =  hNum->GetBinContent(x,y,z);
				d =  hDen->GetBinContent(x,y,z);
				c =  hCoul->GetBinContent(x,y,z);

				if((d > 0.0001) && (n > 0.0001)) 
				{
                    t = par[0] * cfTheory(par, qx, qy, qz, c);
					tempF  += (n*log( (t/n)*((n+d) / (t+1)) ) + d*log( (1.0/d) * ((n+d) / (t+1))));
					// f  += -2. * (n*log( (t / (t+1)) ) + d*log( 1. / (t+1)));
				}

			} // z bins
		} // y bins
	} // x bins

    f = -2 * tempF;
}

void deltaF(double &f, double *par) 
{
	Float_t qx, qy, qz;
    Double_t n, d, c, t; // (n)umerator, (d)enominator, (c)oulomb, (t)heory

	f = 0.;
	
	for (Int_t x = loBin; x <= hiBin; x++) 
	{
        for (Int_t y = loBin; y <= hiBin; y++) 
		{
            for (Int_t z = loBin; z <= hiBin; z++) 
			{
                qx = hNum->GetXaxis()->GetBinCenter(x);
                qy = hNum->GetYaxis()->GetBinCenter(y);
				qz = hNum->GetZaxis()->GetBinCenter(z);

				n =  hNum->GetBinContent(x,y,z);
				d =  hDen->GetBinContent(x,y,z);
				c =  hCoul->GetBinContent(x,y,z);

				if((d > 0.0001) && (n > 0.0001)) 
				{
                    t = par[0] * cfTheory(par, qx, qy, qz, c);
					// f  += -2. * (n*log( (t/n)*((n+d) / (t+1)) ) + d*log( (1.0/d) * ((n+d) / (t+1))));
					f  += 4. * log(t+1); 
				}

			} // z bins
		} // y bins
	} // x bins

}

void makeFitNumerator(TH3D* fitNum, const TH3D* den, const TH3D* coul, const Double_t* fitPars)
{
    const Int_t xMax = den->GetNbinsX();
    const Int_t yMax = den->GetNbinsY();
    const Int_t zMax = den->GetNbinsZ();

	Double_t qx = 0, qy = 0, qz = 0;
    Double_t theoryCFValue = 0, theoryNumValue = 0;

    for (Int_t x = 1; x <= xMax; ++x)
    {
        for (Int_t y = 1; y <= yMax; ++y)
        {
            for (Int_t z = 1; z <= zMax; ++z)
            {
                qx = fitNum->GetXaxis()->GetBinCenter(x);
                qy = fitNum->GetYaxis()->GetBinCenter(y);
                qz = fitNum->GetZaxis()->GetBinCenter(z);

                theoryCFValue = cfTheory(fitPars,qx,qy,qz, coul->GetBinContent(x,y,z));
                theoryNumValue = theoryCFValue * den->GetBinContent(x,y,z); 
                fitNum->Fill(qx,qy,qz,theoryNumValue);
                
            } // z bins
        } // y bins
    } // x bins
}

TGraph* scanNewDistribution(Double_t* pars, const Double_t scanFraction, const Int_t nPar)
{

    Double_t tempValue = 0;
    Int_t nPoints = 200;
    Double_t distPointsX[200] = {0};
    Double_t distPointsY[200] = {0};
    Double_t min = pars[nPar] * (1 - scanFraction / 2.); 
    Double_t stepSize = pars[nPar] * scanFraction / nPoints;

    Double_t offset = 0;
    newLogLikelihood(offset, pars);


    for (Int_t i = 0; i <= (nPoints - 1); i++)
    {
        pars[nPar] = min + i*stepSize;
        newLogLikelihood(tempValue, pars);
        distPointsX[i] = min + i*stepSize;
        distPointsY[i] = tempValue - offset;
    }

    TGraph* gr = new TGraph(nPoints, distPointsX, distPointsY);

    return gr;

}

TGraph* scanDistribution(Double_t* pars, const Double_t scanFraction, const Int_t nPar)
{

    Double_t tempValue = 0;
    Int_t nPoints = 200;
    Double_t distPointsX[200] = {0};
    Double_t distPointsY[200] = {0};
    Double_t min = pars[nPar] * (1 - scanFraction / 2.); 
    Double_t stepSize = pars[nPar] * scanFraction / nPoints;

    Double_t offset = 0;
    logLikelihood(offset, pars);


    for (Int_t i = 0; i <= (nPoints - 1); i++)
    {
        pars[nPar] = min + i*stepSize;
        logLikelihood(tempValue, pars);
        distPointsX[i] = min + i*stepSize;
        distPointsY[i] = tempValue - offset;
    }

    TGraph* gr = new TGraph(nPoints, distPointsX, distPointsY);

    return gr;

}
