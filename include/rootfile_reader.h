#ifndef ANALYSIS_TCTSCAN_ROOTFILE_READER_H_
#define ANALYSIS_TCTSCAN_ROOTFILE_READER_H_

#include "waveform_data.h"

#include <filesystem>
#include <map>
#include <memory>
#include <string>

#include <RtypesCore.h>

class TFile;
class TTree;

namespace tctscan{

class RootfileReader{
    public:
        explicit RootfileReader(const std::filesystem::path& input_path);
        ~RootfileReader();
        RootfileReader(const RootfileReader&) = delete;
        RootfileReader& operator=(const RootfileReader&) = delete;

        const std::vector<ScanPoint>& GetScanPoints() const{return points_;}
        const std::string& GetSourceDirectory() const{return source_directory_;}
        const ScanPoint& GetScanPoint(long long event_number) const;
        void IndexEvents(const std::vector<long long>& event_numbers);
        ShotData ReadShot(long long event_number, int segment_number);
        std::vector<ShotData> ReadShots(long long event_number, const std::vector<int>& segments);

    private:
        void ReadScanInfo();
        void BindWaveforms();

        std::unique_ptr<TFile> file_;
        TTree* waveforms_ = nullptr;
        std::string source_directory_;
        std::vector<ScanPoint> points_;
        std::map<long long, size_t> point_indices_;
        std::map<long long, std::vector<Long64_t>> entries_;

        Long64_t event_number_ = 0;
        Int_t segment_number_ = 0;
        Int_t nsamples_[kNChannels]{};
        Double_t dt_s_[kNChannels]{};
        Double_t t0_s_[kNChannels]{};
        Double_t x_um_ = 0.0;
        Double_t y_um_ = 0.0;
        Double_t z_um_ = 0.0;
        Bool_t clipped_[kNChannels]{};
        std::array<std::vector<float>, kNChannels> voltage_v_;
        std::array<std::vector<float>*, kNChannels> voltage_addresses_{};
};

std::vector<AveragedEvent> ReadAverageRootFile(const std::filesystem::path& input_path);

} // namespace tctscan

#endif
