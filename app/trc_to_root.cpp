#include "rootfile_writer.h"

#include <iostream>

int main(int argc, char* argv[]){
    if(argc == 2 && std::string(argv[1]) == "--help"){
        std::cout << "Usage: " << argv[0] << " INPUT_BASE [OUTPUT.root]\n";
        return 0;
    }
    if(argc < 2 || argc > 3){
        std::cerr << "Usage: " << argv[0] << " INPUT_BASE [OUTPUT.root]\n";
        return 1;
    }
    try{
        std::filesystem::path output = argc == 3 ? argv[2] : "";
        tctscan::ConvertScanToRoot(argv[1], output);
        return 0;
    }catch(const std::exception& error){
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}
