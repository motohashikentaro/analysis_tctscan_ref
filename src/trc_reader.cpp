#include "trc_reader.h"
#include "analysis_support.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <regex>
#include <stdexcept>

namespace tctscan{
namespace{
namespace fs = std::filesystem;
using Bytes = std::vector<unsigned char>;

struct BinaryReader{
    Bytes bytes;
    size_t descriptor_start = 0;
    bool little_endian = true;
    explicit BinaryReader(const fs::path& path){
        std::ifstream f(path, std::ios::binary | std::ios::ate);
        Require(bool(f), "Cannot open trace");
        auto size = f.tellg();
        Require(size > 0 && size < 2147483648LL, "Only nonempty traces <2 GiB supported");
        bytes.resize(static_cast<size_t>(size));
        f.seekg(0);
        f.read(reinterpret_cast<char*>(bytes.data()), size);
        Require(bool(f), "Trace read failed");
        const std::string tag = "WAVEDESC";
        auto it = std::search(bytes.begin(), bytes.end(), tag.begin(), tag.end());
        Require(it != bytes.end(), "WAVEDESC missing");
        descriptor_start = static_cast<size_t>(it - bytes.begin());
        if(bytes[0] == '#'){
            Require(bytes.size() > 2 && bytes[1] >= '1' && bytes[1] <= '9', "Unsupported IEEE block header");
            size_t n = bytes[1] - '0';
            Require(bytes.size() >= 2 + n && descriptor_start == 2 + n, "Invalid IEEE block prefix");
            std::string digits(bytes.begin() + 2, bytes.begin() + 2 + n);
            Require(std::all_of(digits.begin(), digits.end(), [](char c){
                return c >= '0' && c <= '9';
            }), "Invalid block length");
            Require(std::stoull(digits) == bytes.size() - descriptor_start, "Declared block length != file length");
        }else Require(descriptor_start == 0, "Unsupported trace prefix");
        Require(bytes.size() - descriptor_start >= 346, "Truncated descriptor");
        Require(ReadString(16, 16) == "LECROY_2_3", "Only verified LECROY_2_3 template supported");
        if(bytes[descriptor_start + 34] == 1 && bytes[descriptor_start + 35] == 0) little_endian = true;
        else if(bytes[descriptor_start + 34] == 0 && bytes[descriptor_start + 35] == 0) little_endian = false;
        else throw std::runtime_error("Invalid COMM_ORDER");
    }
    uint64_t ReadUnsigned(size_t off, size_t n) const{
        Require(off <= bytes.size() - descriptor_start && n <= bytes.size() - descriptor_start - off, "Out-of-bounds binary field");
        uint64_t v = 0;
        for(size_t i = 0; i < n; ++i) v |= uint64_t(bytes[descriptor_start + off + i]) << (8 * (little_endian ? i : n - 1 - i));
        return v;
    }
    int32_t ReadInt32(size_t o) const{
        auto x = ReadUnsigned(o, 4);
        return static_cast<int32_t>(x <= 0x7fffffff ? int64_t(x) : int64_t(x) - 0x100000000LL);
    }
    int ReadInt16(size_t o) const{
        auto x = ReadUnsigned(o, 2);
        return x < 32768 ? int(x) : int(x) - 65536;
    }
    float ReadFloat(size_t o) const{
        uint32_t x = static_cast<uint32_t>(ReadUnsigned(o, 4));
        float v;
        std::memcpy(&v, &x, 4);
        return v;
    }
    double ReadDouble(size_t o) const{
        uint64_t x = ReadUnsigned(o, 8);
        double v;
        std::memcpy(&v, &x, 8);
        return v;
    }
    std::string ReadString(size_t o, size_t n) const{
        Require(o + n <= bytes.size() - descriptor_start, "Invalid string field");
        std::string s;
        for(size_t i = 0; i < n && bytes[descriptor_start + o + i]; ++i) s += char(bytes[descriptor_start + o + i]);
        return s;
    }
};

} // namespace

TraceData ReadTrc(const fs::path& path, int ch){
    try{
        BinaryReader reader(path);
        TraceData trace;
        auto length = [&](size_t off){
            int v = reader.ReadInt32(off);
            Require(v >= 0, "Negative block length");
            return size_t(v);
        };
        size_t descriptor_length = length(36), user_length = length(40), trigger_bytes = length(48), waveform_bytes = length(60);
        Require(descriptor_length >= 346, "Descriptor shorter than required fields");
        for(size_t off :{
            44, 52, 56, 64, 68, 72
        }) Require(length(off) == 0, "Unsupported auxiliary/RIS/second data block");
        Require(reader.ReadInt32(132) == 0 && reader.ReadInt32(136) == 1, "Partial/sparsed waveforms unsupported");
        Require(reader.ReadInt32(148) == 1 && reader.ReadInt16(318) == 0, "Averaged/processed waveforms unsupported");
        int count = reader.ReadInt32(116);
        trace.nsegments = reader.ReadInt32(144);
        trace.comm_type = reader.ReadInt16(32);
        trace.nominal_bits = reader.ReadInt16(172);
        Require(trace.comm_type == 0 || trace.comm_type == 1, "Unknown COMM_TYPE");
        Require(trace.nsegments > 1 && count > 0 && count % trace.nsegments == 0, "Invalid sequence shape");
        Require(reader.ReadInt32(124) == 0 && reader.ReadInt32(128) == count - 1, "Partial/invalid sample range");
        Require(reader.ReadInt16(344) == ch - 1, "Filename channel != WAVE_SOURCE");
        Require(waveform_bytes == size_t(count) * (trace.comm_type ? 2 : 1), "Waveform byte count mismatch");
        Require(trigger_bytes == size_t(trace.nsegments) * 16, "Trigger byte count mismatch");
        size_t trigger_start = descriptor_length + user_length, waveform_start = trigger_start + trigger_bytes;
        Require(waveform_start <= reader.bytes.size() - reader.descriptor_start && waveform_bytes == reader.bytes.size() - reader.descriptor_start - waveform_start, "Block lengths do not match file size");
        trace.nsamples = count / trace.nsegments;
        trace.gain_v_per_count = reader.ReadFloat(156);
        trace.vertical_offset_v = reader.ReadFloat(160);
        trace.dt_s = reader.ReadFloat(176);
        Require(std::isfinite(trace.gain_v_per_count) && trace.gain_v_per_count > 0 && std::isfinite(trace.vertical_offset_v) && std::isfinite(trace.dt_s) && trace.dt_s > 0, "Invalid scaling");
        Require(reader.ReadString(196, 48) == "V" && reader.ReadString(244, 48) == "S", "Expected voltage [V] and time [s]");
        trace.instrument = reader.ReadString(76, 16);
        std::copy_n(reader.bytes.begin() + reader.descriptor_start + 296, 16, trace.scope_timestamp_bytes.begin());
        trace.trigger_elapsed_s.resize(trace.nsegments);
        trace.t0_s.resize(trace.nsegments);
        for(int j = 0; j < trace.nsegments; ++j){
            trace.trigger_elapsed_s[j] = reader.ReadDouble(trigger_start + 16 * size_t(j));
            trace.t0_s[j] = reader.ReadDouble(trigger_start + 16 * size_t(j) + 8);
            Require(std::isfinite(trace.trigger_elapsed_s[j]) && std::isfinite(trace.t0_s[j]), "Nonfinite trigger data");
            Require(j == 0 ? trace.trigger_elapsed_s[j] == 0 : trace.trigger_elapsed_s[j] > trace.trigger_elapsed_s[j - 1], "Invalid trigger sequence");
        }
        Require(trace.t0_s[0] == reader.ReadDouble(180), "First trigger offset inconsistent with descriptor");
        trace.raw.resize(count);
        for(int k = 0; k < count; ++k){
            int x;
            if(trace.comm_type) x = reader.ReadInt16(waveform_start + 2 * size_t(k));
            else{
                auto v = reader.ReadUnsigned(waveform_start + size_t(k), 1);
                x = v < 128 ? int(v) : int(v) - 256;
            }
            trace.raw[k] = static_cast<short>(x);
        }
        return trace;
    }catch(const std::exception & e){
        throw std::runtime_error(path.string() + ": " + e.what());
    }
}

TraceFiles FindTraceFiles(const fs::path& dir){
    TraceFiles files;
    std::regex re(R"(^C([1-4])--.+--([0-9]+)\.trc$)", std::regex::icase);
    for(const auto& e : fs::directory_iterator(dir)){
        if(!e.is_regular_file()) continue;
        std::smatch m;
        auto name = e.path().filename().string();
        if(!std::regex_match(name, m, re)) continue;
        int ch = std::stoi(m[1]) - 1;
        auto id = std::stoll(m[2]);
        Require(files[id][ch].empty(), "Duplicate trace channel/index: " + name);
        files[id][ch] = e.path();
    }
    return files;
}


} // namespace tctscan
