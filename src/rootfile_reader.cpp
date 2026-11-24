#include "rootfile_reader.h"
#include "analysis_support.h"
#include "waveform_analyzer.h"

#include <algorithm>
#include <cmath>
#include <set>
#include <type_traits>

#include <TBranch.h>
#include <TFile.h>
#include <TNamed.h>
#include <TTree.h>
#include <TLeaf.h>

namespace tctscan{
namespace{

template<typename T>
void BindBranch(TTree& tree, const char* name, T* address, int length = 1){
    auto* branch = tree.GetBranch(name);
    Require(branch != nullptr, std::string("Missing branch: ") + name);
    if constexpr(!std::is_pointer_v<T>){
        auto* leaves = branch->GetListOfLeaves();
        Require(leaves && leaves->GetEntriesFast() == 1, std::string("Unexpected branch shape: ") + name);
        auto* leaf = static_cast<TLeaf*>(leaves->At(0));
        Require(leaf->GetLeafCount() == nullptr && leaf->GetLenStatic() == length,
                std::string("Unexpected branch length: ") + name);
    }
    tree.SetBranchStatus(name, true);
    Require(tree.SetBranchAddress(name, address) == 0, std::string("Unexpected branch type: ") + name);
}

} // namespace

RootfileReader::RootfileReader(const std::filesystem::path& input_path){
    file_.reset(TFile::Open(input_path.c_str(), "READ"));
    Require(file_ && !file_->IsZombie() && !file_->TestBit(TFile::kRecovered),
            "Invalid ROOT file: " + input_path.string());
    auto* source = file_->Get<TNamed>("source_trace_directory");
    if(source) source_directory_ = source->GetTitle();
    waveforms_ = file_->Get<TTree>("waveforms");
    Require(waveforms_ != nullptr, "Missing waveforms TTree");
    ReadScanInfo();
    BindWaveforms();
}

RootfileReader::~RootfileReader(){
    if(waveforms_) waveforms_->ResetBranchAddresses();
}

void RootfileReader::ReadScanInfo(){
    auto* tree = file_->Get<TTree>("scan_info");
    Require(tree != nullptr, "Missing scan_info TTree");
    Long64_t event_number = 0;
    Double_t x_um = 0.0;
    Double_t y_um = 0.0;
    Double_t z_um = 0.0;
    Int_t nsegments = 0;
    tree->SetBranchStatus("*", false);
    BindBranch(*tree, "event_number", &event_number);
    BindBranch(*tree, "x_um", &x_um);
    BindBranch(*tree, "y_um", &y_um);
    BindBranch(*tree, "z_um", &z_um);
    BindBranch(*tree, "nsegments", &nsegments);
    for(Long64_t entry=0; entry<tree->GetEntries(); entry++){
        Require(tree->GetEntry(entry) > 0, "Failed reading scan_info");
        Require(event_number >= 0 && point_indices_.count(event_number) == 0,
                "Invalid or duplicate scan_info event_number");
        Require(nsegments > 0 && nsegments <= waveforms_->GetEntries() &&
                std::isfinite(x_um) && std::isfinite(y_um) && std::isfinite(z_um),
                "Invalid scan metadata");
        point_indices_[event_number] = points_.size();
        points_.push_back({event_number, {x_um, y_um, z_um}, nsegments});
    }
    tree->ResetBranchAddresses();
    Require(!points_.empty(), "scan_info contains no measurement points");
}

const ScanPoint& RootfileReader::GetScanPoint(long long event_number) const{
    auto found = point_indices_.find(event_number);
    Require(found != point_indices_.end(), "Unknown event_number: " + std::to_string(event_number));
    return points_[found->second];
}

void RootfileReader::IndexEvents(const std::vector<long long>& event_numbers){
    std::map<long long, std::vector<Long64_t>> pending;
    for(long long event: event_numbers){
        if(entries_.count(event) > 0) continue;
        const auto& point = GetScanPoint(event);
        pending[event] = std::vector<Long64_t>(point.nsegments, -1);
    }
    if(pending.empty()) return;

    waveforms_->SetBranchStatus("*", false);
    waveforms_->SetBranchStatus("event_number", true);
    waveforms_->SetBranchStatus("segment_number", true);
    for(Long64_t entry=0; entry<waveforms_->GetEntries(); entry++){
        Require(waveforms_->GetEntry(entry) > 0, "Failed indexing waveforms");
        auto found = pending.find(event_number_);
        if(found == pending.end()) continue;
        auto& event_entries = found->second;
        Require(segment_number_ >= 0 && size_t(segment_number_) < event_entries.size(),
                "Invalid segment_number at event " + std::to_string(event_number_));
        Require(event_entries[segment_number_] == -1,
                "Duplicate segment_number at event " + std::to_string(event_number_));
        event_entries[segment_number_] = entry;
    }
    BindWaveforms();
    for(auto& [event, event_entries]: pending){
        Require(std::none_of(event_entries.begin(), event_entries.end(), [](Long64_t entry){return entry < 0;}),
                "Missing waveform at event " + std::to_string(event));
        entries_[event] = std::move(event_entries);
    }
}

void RootfileReader::BindWaveforms(){
    waveforms_->SetBranchStatus("*", false);
    BindBranch(*waveforms_, "event_number", &event_number_);
    BindBranch(*waveforms_, "segment_number", &segment_number_);
    BindBranch(*waveforms_, "nsamples", nsamples_, kNChannels);
    BindBranch(*waveforms_, "dt_s", dt_s_, kNChannels);
    BindBranch(*waveforms_, "t0_s", t0_s_, kNChannels);
    BindBranch(*waveforms_, "x_um", &x_um_);
    BindBranch(*waveforms_, "y_um", &y_um_);
    BindBranch(*waveforms_, "z_um", &z_um_);
    BindBranch(*waveforms_, "clipped", clipped_, kNChannels);
    for(int ch=0; ch<kNChannels; ch++){
        voltage_addresses_[ch] = &voltage_v_[ch];
        std::string branch = "voltage_c" + std::to_string(ch + 1);
        BindBranch(*waveforms_, branch.c_str(), &voltage_addresses_[ch]);
        waveforms_->GetBranch(branch.c_str())->SetAutoDelete(false);
    }
}

ShotData RootfileReader::ReadShot(long long event_number, int segment_number){
    const auto& point = GetScanPoint(event_number);
    Require(segment_number >= 0 && segment_number < point.nsegments, "Requested segment does not exist");
    IndexEvents({event_number});
    Require(waveforms_->GetEntry(entries_.at(event_number)[segment_number]) > 0, "Failed reading selected shot");
    Require(event_number_ == event_number && segment_number_ == segment_number,
            "Indexed waveform identity changed");
    Require(x_um_ == point.position.x_um && y_um_ == point.position.y_um && z_um_ == point.position.z_um,
            "Coordinates differ between waveforms and scan_info");

    ShotData shot;
    shot.event_number = event_number;
    shot.segment_number = segment_number;
    shot.position = point.position;
    for(int ch=0; ch<kNChannels; ch++){
        Require(voltage_addresses_[ch] == &voltage_v_[ch] && nsamples_[ch] >= 2 &&
                size_t(nsamples_[ch]) == voltage_v_[ch].size(), "Waveform vector length mismatch");
        Require(std::isfinite(t0_s_[ch]) && std::isfinite(dt_s_[ch]) && dt_s_[ch] > 0.0,
                "Invalid waveform time grid");
        std::vector<double> values(voltage_v_[ch].begin(), voltage_v_[ch].end());
        for(double value: values) Require(std::isfinite(value), "Nonfinite voltage sample");
        shot.waveforms[ch] = MakeWaveform(t0_s_[ch], dt_s_[ch], values, clipped_[ch]);
    }
    return shot;
}

std::vector<ShotData> RootfileReader::ReadShots(long long event_number, const std::vector<int>& segments){
    Require(!segments.empty(), "No shots selected");
    std::set<int> used;
    std::vector<ShotData> shots;
    shots.reserve(segments.size());
    for(int segment: segments){
        Require(used.insert(segment).second, "Duplicate segment in ReadShots");
        shots.push_back(ReadShot(event_number, segment));
    }
    return shots;
}

std::vector<AveragedEvent> ReadAverageRootFile(const std::filesystem::path& input_path){
    TFile file(input_path.c_str(), "READ");
    Require(!file.IsZombie() && !file.TestBit(TFile::kRecovered), "Invalid average ROOT file");
    auto* tree = file.Get<TTree>("averages");
    Require(tree != nullptr, "Missing averages TTree");
    Long64_t event = 0;
    Int_t channel = 0;
    Int_t total = 0;
    Double_t x = 0.0, y = 0.0, z = 0.0, dt = 0.0, t0_min = 0.0, t0_max = 0.0;
    ULong64_t n_used = 0, n_clipped = 0;
    std::vector<int>* segments = nullptr;
    std::vector<double>* times = nullptr;
    std::vector<double>* values = nullptr;
    BindBranch(*tree, "event_number", &event);
    BindBranch(*tree, "channel", &channel);
    BindBranch(*tree, "nsegments", &total);
    BindBranch(*tree, "x_um", &x);
    BindBranch(*tree, "y_um", &y);
    BindBranch(*tree, "z_um", &z);
    BindBranch(*tree, "dt_s", &dt);
    BindBranch(*tree, "t0_min_s", &t0_min);
    BindBranch(*tree, "t0_max_s", &t0_max);
    BindBranch(*tree, "n_used", &n_used);
    BindBranch(*tree, "n_clipped", &n_clipped);
    BindBranch(*tree, "segments", &segments);
    BindBranch(*tree, "time_s", &times);
    BindBranch(*tree, "voltage_v", &values);
    std::vector<AveragedEvent> events;
    std::map<long long, size_t> indices;
    std::set<std::pair<long long, int>> seen;
    try{
        for(Long64_t entry=0; entry<tree->GetEntries(); entry++){
            Require(tree->GetEntry(entry) > 0 && channel >= 1 && channel <= kNChannels &&
                    segments && times && values && n_used == segments->size() && n_used > 0 &&
                    n_clipped <= n_used && times->size() == values->size() && times->size() >= 2,
                    "Invalid saved average");
            Require(seen.insert({event, channel}).second, "Duplicate average channel");
            if(indices.count(event) == 0){
                indices[event] = events.size();
                AveragedEvent result;
                result.point = {event, {x, y, z}, total};
                result.segments = *segments;
                events.push_back(std::move(result));
            }
            auto& result = events[indices.at(event)];
            Require(result.segments == *segments && result.point.position.x_um == x &&
                    result.point.position.y_um == y && result.point.position.z_um == z &&
                    result.point.nsegments == total, "Inconsistent saved average metadata");
            result.waveforms[channel - 1] = {{*times, *values, dt, n_clipped > 0}, t0_min, t0_max,
                                            size_t(n_used), size_t(n_clipped)};
        }
        Require(!events.empty() && seen.size() == events.size() * kNChannels,
                "Incomplete saved average channels");
    }catch(...){
        tree->ResetBranchAddresses();
        delete segments;
        delete times;
        delete values;
        throw;
    }
    tree->ResetBranchAddresses();
    delete segments;
    delete times;
    delete values;
    return events;
}

} // namespace tctscan
