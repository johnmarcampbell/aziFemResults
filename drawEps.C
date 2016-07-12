void formatTMGraph(TMultiGraph* tmGraph);
void drawEps()
{

    TFile* auauFile = new TFile("193GeVData/UU193Eps_smallFitRange.root","READ");
    // TFile* auauFile = new TFile("200GeVData/AuAu200Eps_smallFitRange.root","READ");
    // TFile* uuFile = new TFile("193GeVData/UU193Eps_noRes.root","READ");
    // TFile* auauFile = new TFile("200GeVData/AuAu200Eps_noRes.root","READ");
    TFile* uuFile = new TFile("193GeVData/UU193Eps_WithRes.root","READ");
    // TFile* auauFile = new TFile("200GeVData/AuAu200Eps_WithRes.root","READ");
    TGraphErrors* tgUU[2][2][2]; //  [epsF/Rl][q2/mult][zdc0/zdc1]
    TGraphErrors* tgAuAu[2][2][2]; //[epsF/Rl][q2/mult][zdc0/zdc1]
    TMultiGraph* tmg[2][2][2];

    for (Int_t i = 0; i <= 1; ++i) { // epsF/Rl
        for (Int_t j = 0; j <= 1; j++) { //q2/mult
           for (Int_t k = 0; k <= 1; ++k) { //zdc0/zdc1
                TString uuName = "uu";
                TString auauName = "uu";
                // TString auauName = "auau";
                if( i == 0 ) {
                    uuName += "EpsF_Zdc_";
                    auauName += "EpsF_Zdc_";
                } else {
                    uuName += "RLong_Zdc_";
                    auauName += "RLong_Zdc_";
                }

                uuName += k;
                auauName += k;
                if( j == 0 ) {
                    uuName += "_q2";
                    auauName += "_q2";
                } else {
                    uuName += "_mult";
                    auauName += "_mult";
                }

                tgUU[i][j][k] = (TGraphErrors*)uuFile->Get(uuName.Data());
                tgUU[i][j][k]->SetMarkerStyle(kFullStar);
                tgUU[i][j][k]->SetMarkerColor(kBlue);
                tgUU[i][j][k]->SetMarkerSize(2);
                tgUU[i][j][k]->SetLineColor(kBlue);

                tgAuAu[i][j][k] = (TGraphErrors*)auauFile->Get(auauName.Data());
                tgAuAu[i][j][k]->SetMarkerStyle(kFullStar);
                tgAuAu[i][j][k]->SetMarkerColor(kRed);
                tgAuAu[i][j][k]->SetMarkerSize(2);
                tgAuAu[i][j][k]->SetLineColor(kRed);

                tmg[i][j][k] = new TMultiGraph();
                tmg[i][j][k]->Add(tgUU[i][j][k]);
                tmg[i][j][k]->Add(tgAuAu[i][j][k]);
                tmg[0][j][k]->SetMinimum(-0.12);
                tmg[0][j][k]->SetMaximum(0.12);
            } 

            // tgUU[i][j][1]->SetLineStyle(2);
            // tgAuAu[i][j][1]->SetLineStyle(2);
        }
    }

    TCanvas* epsFCanvas = new TCanvas("epsFCanvas","epsFCanvas",1000,1000);
    epsFCanvas->Divide(2,2);

    epsFCanvas->cd(1);
    tmg[0][0][0]->Draw("alp");
    tmg[0][0][0]->SetTitle("0-0.25% Zdc -- q2");

    epsFCanvas->cd(2);
    tmg[0][1][0]->Draw("alp");
    tmg[0][1][0]->SetTitle("0-0.25% Zdc -- mult");

    epsFCanvas->cd(3);
    tmg[0][0][1]->Draw("alp");
    tmg[0][0][1]->SetTitle("0-0.5% Zdc -- q2");

    epsFCanvas->cd(4);
    tmg[0][1][1]->Draw("alp");
    tmg[0][1][1]->SetTitle("0-0.5% Zdc -- mult");

    TCanvas* RlCanvas = new TCanvas("RlCanvas","RlCanvas",1000,1000);
    RlCanvas->Divide(2,2);

    RlCanvas->cd(1);
    tmg[1][0][0]->Draw("alp");
    tmg[1][0][0]->SetTitle("0-0.25% Zdc -- q2");

    RlCanvas->cd(2);
    tmg[1][1][0]->Draw("alp");
    tmg[1][1][0]->SetTitle("0-0.25% Zdc -- mult");

    RlCanvas->cd(3);
    tmg[1][0][1]->Draw("alp");
    tmg[1][0][1]->SetTitle("0-0.5% Zdc -- q2");

    RlCanvas->cd(4);
    tmg[1][1][1]->Draw("alp");
    tmg[1][1][1]->SetTitle("0-0.5% Zdc -- mult");

}

void formatTMGraph(TMultiGraph* tmGraph)
{

    tmGraph->GetXaxis()->SetNdivisions(105);

}
