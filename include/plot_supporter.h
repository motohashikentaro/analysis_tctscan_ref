#ifndef ANALYSIS_TCTSCAN_PLOT_SUPPORTER_H_
#define ANALYSIS_TCTSCAN_PLOT_SUPPORTER_H_

#include <string>

class TVirtualPad;
class TAxis;

namespace tctscan{

class PlotSupporter{
    public:
        static void Initialize();
        static int ChannelColor(int channel);
        static void SetPadStyle(TVirtualPad* pad);
        static void SetAxesStyle(TAxis* x_axis, TAxis* y_axis);
        static void DrawLabel(double x, double y, const std::string& text, double size = 0.030, int color = 1);
};

} // namespace tctscan

#endif
