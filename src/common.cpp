
#include "common.h"
#include <string>

std::string to_upper(std::string &str) {
    std::string upper = str;
    for (auto & c: upper) 
        c = toupper(c);
    return upper;
}
