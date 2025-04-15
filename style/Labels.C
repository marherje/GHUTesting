#include "Labels.h"
#include "TLatex.h"
#include "TLine.h"
#include "TPave.h"
#include "TMarker.h"

void QQBARLabel(Double_t x,Double_t y,TString text,Color_t color, Double_t size)
{
  TLatex l;
  //l.SetTextAlign(12);
  l.SetNDC();
  l.SetTextFont(42);
  l.SetTextSize(size);

  if(color>0) l.SetTextColor(color);
  else {
    l.SetTextFont(62);
    Int_t ci = TColor::GetFreeColorIndex();
    TColor *color = new TColor(ci, 0, 0.259, 0.549);
    l.SetTextColor(ci);
  }
  double delx = 1200*gPad->GetWh()/(1000*gPad->GetWw());
  
  l.DrawLatex(x,y,"ILD");
  if (text) {
    TLatex p;
    p.SetNDC();
    p.SetTextSize(0.03);
    p.SetTextFont(50);
    p.SetTextColor(kRed-3);
    p.DrawLatex(x-0.132,y-0.03,text);
  }
}

void QQBARLabel2(Double_t x,Double_t y,TString text,Color_t color,Double_t textsize,Double_t angle)
{
  
  TLatex p;
  p.SetNDC();
  p.SetTextSize(textsize);
  p.SetTextFont(52);
  p.SetTextColor(color);
  p.SetTextAngle(angle);
  p.DrawLatex(x,y,text);
  
}

void QQBARLabel3(Double_t x,Double_t y,TString text,Color_t color,Double_t textsize)
{

  text = text.Strip(TString::kBoth); // Remove leading and trailing spaces
  TLatex p;
  p.SetNDC();
  p.SetTextSize(textsize);
  p.SetTextFont(42);
  p.SetTextColor(color);
  p.DrawLatex(x,y,text);

}

void QQBARLabel4(Double_t x,Double_t y,TString text,Color_t color,Double_t textsize,Double_t angle)
{

  TLatex p;
  p.SetNDC();
  p.SetTextSize(textsize);
  p.SetTextFont(42);
  p.SetTextColor(color);
  p.SetTextAngle(angle);
  p.DrawLatex(x,y,text);

}

