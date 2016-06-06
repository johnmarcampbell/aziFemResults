Double_t getChi(const Double_t res);
Double_t getRealRes(const Double_t chi);

void resCalc(
    const TString AuAuFileName = "./data/AuAu200GeVResolution.root",
    const TString UUFileName = "./data/UU193GeVResolution.root"
    
    )
{

	TFile* fInAuAu = new TFile(AuAuFileName.Data());
	TFile* fInUU = new TFile(UUFileName.Data());
    TMultiGraph* auMultiGraph = makeMultiGraph(fInAuAu);
    TMultiGraph* uMultiGraph = makeMultiGraph(fInUU);

    TCanvas* canvas = new TCanvas("resolution","Event Plane Resolution");
    canvas->Divide(2,1);

    canvas->cd(1);
    auMultiGraph->Draw("ap");
    auMultiGraph->SetMinimum(0.01);
    auMultiGraph->SetMaximum(0.5);

    canvas->cd(2);
    uMultiGraph->Draw("ap");
    uMultiGraph->SetMinimum(0.01);
    uMultiGraph->SetMaximum(0.5);

}

TMultiGraph* makeMultiGraph(TFile* fIn)
{

	TProfile2D* subResSquaredMult = (TProfile2D*)fIn->Get("subResSquaredRefmultInclusiveZdcCut");
	TProfile2D* subAuq2 = (TProfile2D*)fIn->Get("subResSquaredq2InclusiveZdcCut");

    TMultiGraph* multiGraph = new TMultiGraph();
    Int_t nBins = subResSquaredMult->GetNbinsX();
	TGraphErrors* multRes[2];
	TGraphErrors* q2Res[2];

    for (Int_t i = 0; i <= 1; ++i)
    {
        multRes[i] = new TGraphErrors();
        q2Res[i] = new TGraphErrors(nBins);
        makeResGraph(subResSquaredMult, multRes[i], i+1);
        makeResGraph(subAuq2, q2Res[i], i+1);

        Int_t lineStyle = i ? 1 : 3;
        Int_t markerStyle = i ? 29 : 30;
        multRes[i]->SetLineStyle(lineStyle);
        multRes[i]->SetLineWidth(3);
        multRes[i]->SetMarkerStyle(markerStyle);
        multRes[i]->SetMarkerSize(2);
        q2Res[i]->SetLineStyle(lineStyle);
        q2Res[i]->SetLineWidth(3);
        q2Res[i]->SetMarkerStyle(markerStyle);
        q2Res[i]->SetMarkerSize(2);
        q2Res[i]->SetMarkerColor(kRed);
        q2Res[i]->SetLineColor(kRed);
        
        multiGraph->Add(multRes[i]);
        multiGraph->Add(q2Res[i]);
    }

    return multiGraph;
}

void makeResGraph(TProfile2D* subResSquared, TGraphErrors* fullResGraph, const Int_t yBin)
{
	Double_t chiSub = 0;
    Double_t subRes = 0;
    Double_t fullRes = 0;
    Int_t nBins = subResSquared->GetNbinsX();

	cout << "{";
	for (int i = 1; i <= nBins; ++i)
	{

        subRes = sqrt(subResSquared->GetBinContent(i,yBin));
        chiSub = getChi(subRes);
        fullRes = getRealRes(sqrt(2) * chiSub);

        cout << fullRes;
        if(i != nBins) {cout << ", ";}
        fullResGraph->SetPoint(i,i,fullRes);
        fullResGraph->SetPointError(i,0.5,0);

	}

	cout << "}" << endl;


}

Double_t getChi(const Double_t res)
{
	Double_t chi;
	TF1 Subf( "Subf_som", "sqrt(pi)*x*exp(-x*x/2)*(1/2)*( TMath::BesselI(1/2,x*x/2) + TMath::BesselI(3/2,x*x/2) ) - [0]", 0, 10.0 ); 
	Subf.SetParameters(res, 0.0);
	ROOT::Math::WrappedTF1 wf1(Subf);
	ROOT::Math::BrentRootFinder brf;
	brf.SetFunction( wf1, 0, 10.0 );
	brf.Solve();
	chi = brf.Root();

	return chi;

}

Double_t getRealRes(const Double_t chi)
{
	Double_t pi = 3.14159;
	Double_t res = 0;

	res = TMath::BesselI(1/2,chi*chi/2.0);
	res += TMath::BesselI(3/2,chi*chi/2.0);
	res *= exp(-chi*chi/2.0);
	res *= sqrt(pi) * chi / 2.0;
	return res;

}

Double_t getResFromHistogram(TH2D* hist)
{

    Int_t nBins = hist->GetNbinsX();
    Float_t subResSquared = 0;
    Int_t tally = 0;

    for (Int_t i = 1; i <= nBins; ++i)
    {
        for (Int_t j = 1; j <= nBins; ++j)
        {
            Float_t psiA = hist->GetXaxis()->GetBinCenter(i);
            Float_t psiB = hist->GetYaxis()->GetBinCenter(j);

            subResSquared += cos(2*(psiA-psiB)) * hist->GetBinContent(i,j);
            tally += hist->GetBinContent(i,j);

        }
            
    }

    subResSquared /= tally;
    cout << subResSquared << endl;
    Float_t subRes = sqrt(subResSquared);
    Float_t chiSub = getChi(subRes);
    Float_t fullRes = getRealRes(sqrt(2) * chiSub);

    cout << fullRes << endl;

    return fullRes;

}
