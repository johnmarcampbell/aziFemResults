const TH3D* hNum;
const TH3D* hDen;
const TH3D* hCoul;
TMinuit* tmFit;
Int_t loBin;
Int_t hiBin;

Double_t cfTheory(const Double_t* fitPars, Double_t qx, Double_t qy, Double_t qz, Double_t coul);
void logLikelihood(int &npar, double *gin, double &f, double *par, int iflag);
void makeFitNumerator(TH3D* fitNum, const TH3D* den, const TH3D* coul, const Double_t* fitPars);

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
    TGraph** contour = 0
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
    Double_t initRange[8] = {0.1, 0.1, 3, 3, 3, 1, 0, 0};
    TString parNames[8] = {"Normalization" , "Lambda" , "Rout Squared" , "Rside Squared" , "Rlong Squared" , "Routside Squared" , "Routlong Squared" , "Rsidelong Squared"};
    tmFit->SetFCN(logLikelihood);
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
    TCanvas* test = new TCanvas();
    *contour = (TGraph*)(tmFit->Contour(5,0,1))->Clone("test");
    // contour->Draw();
    // test->Print("test.pdf");


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

// This 
void logLikelihood(int &npar, double *gin, double &f, double *par, int iflag) 
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
					f  += -2. * (n*log( (t/n)*((n+d) / (t+1)) ) + d*log( (1.0/d) * ((n+d) / (t+1))));
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
