#include "TCanvas.h"
#include "TH2.h"
#include "THStack.h"
#include "TRandom.h"

void multicolor(Int_t stack=0) {
   TCanvas *c1 = new TCanvas;
   c1->cd();
   c1->Draw();
   Int_t nbins = 20;
   TH2F *h1 = new TH2F("h1","h1",nbins,-4,4,nbins,-4,4);
   h1->SetFillColor(kBlue);
//   h1->SetMarkerColor(kBlue);
   TH2F *h2 = new TH2F("h2","h2",nbins,-4,4,nbins,-4,4);
   h2->SetFillColor(kRed);
//   h2->SetMarkerColor(kRed);
   TH2F *h3 = new TH2F("h3","h3",nbins,-4,4,nbins,-4,4);
   h3->SetFillColor(kYellow);
//   h3->SetMarkerColor(kYellow);
   THStack *hs = new THStack("hs","three plots");
   hs->Add(h1);
   hs->Add(h2);
   hs->Add(h3);
   TRandom r;
   Int_t i;
   for (i=0;i<20000;i++) h1->Fill(r.Gaus(),r.Gaus());
   for (i=0;i<200;i++) {
      Int_t ix = (Int_t)r.Uniform(0,nbins);
      Int_t iy = (Int_t)r.Uniform(0,nbins);
      Int_t bin = h1->GetBin(ix,iy);
      Double_t val = h1->GetBinContent(bin);
      if (val <= 0) continue;
      if (!stack) h1->SetBinContent(bin,0);
      if (r.Rndm() > 0.5) {
         if (!stack) h2->SetBinContent(bin,0);
         h3->SetBinContent(bin,val);
      } else {
         if (!stack) h3->SetBinContent(bin,0);
         h2->SetBinContent(bin,val);
      }
   }

//   hs->Draw("cont1");
//   hs->Draw("text,same");
   hs->Draw("textColZ");
}
