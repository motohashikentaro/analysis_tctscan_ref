#include "analysis_support.h"

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <set>
#include <sstream>
#include <stdexcept>

namespace tctscan{

void Require(bool condition, const std::string& message){
    if(!condition) throw std::runtime_error(message);
}

std::string AlignmentName(Alignment alignment){
    return alignment == Alignment::kTrigger ? "trigger" : "sample-index";
}

std::string PaddedNumber(long long number, int width){
    std::ostringstream text;
    text << std::setw(width) << std::setfill('0') << number;
    return text.str();
}

std::string NumberText(double value){
    if(!std::isfinite(value)) return "N/A";
    std::ostringstream text;
    text << std::setprecision(5) << value;
    return text.str();
}

std::string SelectionText(const std::vector<int>& numbers){
    std::ostringstream text;
    for(size_t first=0; first<numbers.size();){
        size_t last = first;
        while(last + 1 < numbers.size() && numbers[last + 1] - numbers[last] == 1) last++;
        if(first > 0) text << ',';
        text << numbers[first];
        if(last > first) text << ':' << numbers[last];
        first = last + 1;
    }
    return text.str();
}

std::vector<long long> ParseNumbers(const std::string& text, long long limit){
    Require(!text.empty(), "Empty number selection");
    auto parse_number = [](const std::string& word){
        long long number = 0;
        auto parsed = std::from_chars(word.data(), word.data() + word.size(), number);
        Require(!word.empty() && parsed.ec == std::errc{} &&
                parsed.ptr == word.data() + word.size() && number >= 0,
                "Invalid nonnegative integer: " + word);
        return number;
    };

    std::vector<long long> numbers;
    std::set<long long> used;
    size_t start = 0;
    while(start < text.size()){
        size_t stop = text.find(',', start);
        std::string token = text.substr(start, stop == std::string::npos ? stop : stop - start);
        size_t colon = token.find(':');
        long long first = parse_number(token.substr(0, colon));
        long long last = colon == std::string::npos ? first : parse_number(token.substr(colon + 1));
        Require(first <= last && last <= limit, "Invalid selection range: " + token);
        Require(last - first <= 1000000 && numbers.size() + size_t(last - first) + 1 <= 1000000,
                "Selection exceeds one million entries");
        for(long long number=first;; number++){
            Require(used.insert(number).second, "Duplicate number: " + std::to_string(number));
            numbers.push_back(number);
            if(number == last) break;
        }
        if(stop == std::string::npos) break;
        start = stop + 1;
        Require(start < text.size(), "Trailing comma in selection");
    }
    return numbers;
}

std::pair<double, double> ParseRange(const std::string& text){
    size_t colon = text.find(':');
    Require(colon != std::string::npos, "Time window requires BEGIN:END");
    size_t nbegin = 0;
    size_t nend = 0;
    double begin = std::stod(text.substr(0, colon), &nbegin);
    double end = std::stod(text.substr(colon + 1), &nend);
    Require(nbegin == colon && nend == text.size() - colon - 1 &&
            std::isfinite(begin) && std::isfinite(end) && begin <= end,
            "Invalid time window: " + text);
    return {begin, end};
}

TimeWindow ParseTimeWindow(const std::string& text){
    auto [begin, end] = ParseRange(text);
    return {begin, end};
}

std::vector<int> SelectSegments(int nsegments, const SegmentSelection& selection){
    Require(nsegments > 0 && selection.start >= 0 && selection.count >= 0,
            "Invalid segment selection");
    if(!selection.numbers.empty()){
        std::vector<int> numbers = selection.numbers;
        std::sort(numbers.begin(), numbers.end());
        Require(numbers.front() >= 0 && numbers.back() < nsegments,
                "Requested segment does not exist");
        Require(std::adjacent_find(numbers.begin(), numbers.end()) == numbers.end(),
                "Duplicate segment number");
        return numbers;
    }
    Require(selection.start < nsegments, "Segment start is outside measurement point");
    int count = selection.count == 0 ? nsegments - selection.start : selection.count;
    Require(count <= nsegments - selection.start, "Not enough shots at measurement point");
    std::vector<int> numbers;
    for(int i=0; i<count; i++) numbers.push_back(selection.start + i);
    return numbers;
}

std::string CsvQuote(const std::string& value){
    std::string quoted = "\"";
    for(char c: value){
        quoted += c;
        if(c == '"') quoted += '"';
    }
    return quoted + '"';
}

std::string JsonQuote(const std::string& value){
    std::ostringstream quoted;
    quoted << '"';
    for(unsigned char c: value){
        if(c == '"' || c == '\\') quoted << '\\' << c;
        else if(c < 0x20) quoted << "\\u" << std::hex << std::setw(4) << std::setfill('0') << int(c) << std::dec;
        else quoted << c;
    }
    quoted << '"';
    return quoted.str();
}

void CheckPng(const std::filesystem::path& path){
    std::ifstream input(path, std::ios::binary);
    std::array<unsigned char, 8> signature{};
    input.read(reinterpret_cast<char*>(signature.data()), signature.size());
    Require(input.gcount() == 8 && signature == std::array<unsigned char, 8>{137,80,78,71,13,10,26,10},
            "PNG output failed: " + path.string());
}

} // namespace tctscan
