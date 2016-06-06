void q2TheoryPlots()
{
	TF1* dPdq2[10];
    TH1F* histDPdq2[10];
    TH1F* histRatios[9];
    Float_t nSample = 1000000;
    Float_t v2Low = 0.01;
    Int_t mult = 600;

    for (Int_t i = 0; i <= 9; ++i)
    {
        
        TString fTitle = TString::Format("func_%d",i);
        TString hTitle = TString::Format("hist_%d",i);
        TString rTitle = TString::Format("ratio_%d",i);
        dPdq2[i] = new TF1(fTitle.Data(), "(x/[2])*exp( -([0]*[0]*[1] + x*x)/(2*[2]) )* TMath::BesselI(0, (x*[0]*sqrt([1]))/([2]) )", 0, 5); 
        dPdq2[i]->SetParameters(v2Low + 0.001*i, mult, 0.5);

        histDPdq2[i] = new TH1F(hTitle.Data(), hTitle.Data(), 100, 0, 5);
        histDPdq2[i]->FillRandom(fTitle.Data(), nSample);
        dPdq2[i]->SetLineColor(i);
        histDPdq2[i]->SetLineColor(i);

        cout << "v2: " << v2Low + 0.001*i << "\t <q2>: " << dPdq2[i]->Mean(0,5) << endl;

        // if(!i) { histDPdq2[i]->Draw();}
        // else { histDPdq2[i]->Draw("same");}

    }

    for (Int_t i = 0; i <= 8; ++i)
    {

        TString rTitle = TString::Format("ratio_%d",i);
        histRatios[i] = (TH1F*)histDPdq2[i+1]->Clone(rTitle.Data());
        histRatios[i]->Divide(histDPdq2[i]);

        if(i==0) { histRatios[i]->Draw();}
        else { histRatios[i]->Draw("same");}
    }
        

}
