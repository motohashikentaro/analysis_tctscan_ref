#include "plot_supporter.h"
#include "analysis_support.h"

#include <array>

#include <TAxis.h>
#include <TColor.h>
#include <TH1.h>
#include <TLatex.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>

namespace tctscan{

void PlotSupporter::Initialize(){
    gROOT->SetBatch(true);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    TH1::AddDirectory(false);
}

int PlotSupporter::ChannelColor(int channel){
    Require(channel >= 0 && channel < 4, "Invalid plot channel");
    static const std::array<int, 4> kColors = {
        TColor::GetColor("#0072B2"), TColor::GetColor("#D55E00"),
        TColor::GetColor("#009E73"), TColor::GetColor("#CC79A7")
    };
    return kColors[channel];
}

void PlotSupporter::SetPadStyle(TVirtualPad* pad){
    pad->SetFillColor(kWhite);
    pad->SetFrameFillColor(kWhite);
    pad->SetBorderMode(0);
    pad->SetGrid(0, 0);
    pad->SetTicks(1, 1);
    pad->SetLeftMargin(0.13);
    pad->SetRightMargin(0.03);
    pad->SetBottomMargin(0.12);
    pad->SetTopMargin(0.03);
}

void PlotSupporter::SetAxesStyle(TAxis* x_axis, TAxis* y_axis){
    for(auto* axis: {x_axis, y_axis}){
        axis->SetTitleFont(42);
        axis->SetLabelFont(42);
        axis->SetTitleSize(0.045);
        axis->SetLabelSize(0.04);
        axis->SetNdivisions(510);
    }
    x_axis->SetTitleOffset(1.0);
    y_axis->SetTitleOffset(1.1);
}

void PlotSupporter::DrawLabel(double x, double y, const std::string& text, double size, int color){
    TLatex label;
    label.SetNDC();
    label.SetTextFont(42);
    label.SetTextSize(size);
    label.SetTextColor(color);
    label.DrawLatex(x, y, text.c_str());
}

} // namespace tctscan
