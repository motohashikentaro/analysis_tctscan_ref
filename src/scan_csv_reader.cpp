#include "scan_csv_reader.h"
#include "analysis_support.h"

#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>

namespace tctscan{
namespace{
namespace fs = std::filesystem;

std::string Trim(std::string s){
    auto a = s.find_first_not_of(" \t\r\n");
    return a == std::string::npos ? "" : s.substr(a, s.find_last_not_of(" \t\r\n") - a + 1);
}

std::vector<std::string> SplitCsvFields(const std::string& line){
    std::vector<std::string> out;
    std::string cell;
    bool quoted = false;
    for(size_t i = 0; i < line.size(); ++i){
        char c = line[i];
        if(c == '"'){
            if(quoted && i + 1 < line.size() && line[i + 1] == '"'){
                cell += '"';
                ++i;
            }else quoted = !quoted;
        }else if(c == ',' && !quoted){
            out.push_back(Trim(cell));
            cell.clear();
        }else cell += c;
    }
    Require(!quoted, "Unclosed CSV quote (multiline cells unsupported)");
    out.push_back(Trim(cell));
    return out;
}

long long ReadInteger(const std::string& s){
    size_t n = 0;
    auto v = std::stoll(s, &n);
    Require(n == s.size(), "Invalid integer: " + s);
    return v;
}

double ReadReal(const std::string& s){
    size_t n = 0;
    auto v = std::stod(s, &n);
    Require(n == s.size() && std::isfinite(v), "Invalid number: " + s);
    return v;
}


} // namespace

std::vector<ScanCsvRow> ReadScanCsv(const fs::path& path){
    std::ifstream input(path);
    Require(bool(input), "Cannot open CSV: " + path.string());
    std::map<std::string, size_t> columns;
    std::vector<ScanCsvRow> rows;
    std::set<long long> event_numbers;
    std::string line;
    size_t lineno = 0;
    while(std::getline(input, line)){
        ++lineno;
        if(lineno == 1 && line.compare(0, 3, "\xEF\xBB\xBF") == 0) line.erase(0, 3);
        line = Trim(line);
        if(line.empty()) continue;
        auto cells = SplitCsvFields(line);
        if(cells.empty() || cells[0].empty() || cells[0][0] == '#') continue;
        if(columns.empty()){
            for(size_t i = 0; i < cells.size(); ++i) Require(columns.emplace(cells[i], i).second, "Duplicate CSV column");
            for(auto key :{
                "index", "x_um", "y_um", "z_um"
            }) Require(columns.count(key), std::string("Missing CSV column: ") + key);
            continue;
        }
        auto get = [&](const std::string& key)->std::string{
            auto it = columns.find(key);
            if(it == columns.end()) return "";
            Require(it->second < cells.size(), "CSV row too short at line " + std::to_string(lineno));
            return cells[it->second];
        };
        ScanCsvRow row;
        row.event_number = ReadInteger(get("index"));
        Require(row.event_number >= 0 && event_numbers.insert(row.event_number).second, "Invalid/duplicate CSV index");
        row.position.x_um = ReadReal(get("x_um"));
        row.position.y_um = ReadReal(get("y_um"));
        row.position.z_um = ReadReal(get("z_um"));
        auto tune = get("tune_pct");
        row.tune_pct = tune.empty() ? std::numeric_limits<double>::quiet_NaN() : ReadReal(tune);
        row.status = get("status");
        row.timestamp = get("timestamp");
        Require(row.status.empty() || row.status == "ok", "CSV status is not ok at index " + std::to_string(row.event_number));
        rows.push_back(row);
    }
    Require(!rows.empty(), "No CSV data rows");
    return rows;
}


} // namespace tctscan
