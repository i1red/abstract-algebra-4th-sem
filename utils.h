#ifndef ABSTRACT_ALGEBRA2_UTILS_H
#define ABSTRACT_ALGEBRA2_UTILS_H

#include <string>
#include <sstream>
#include <algorithm>
#include <vector>
#include <map>
#include "BigInt/BigInt.h"


const char DEFAULT_VAR = 'x';


std::vector<std::string> split(const std::string& s, char delimiter) {
    std::vector<std::string> res;
    std::stringstream ss(s);
    std::string token;

    while (std::getline(ss, token, delimiter)) {
        res.push_back(token);
    }

    return res;
}

std::string removeSpaces(std::string s) {
    s.erase(remove_if(s.begin(), s.end(), isspace), s.end());
    return s;
}

template <typename T>
T fromString(const std::string& num) {
    T res = 0;
    std::stringstream converter(num);
    converter >> res;

    return res;
}

template <>
inline BigInt fromString(const std::string& num) {
    return BigInt(num);
}

template <typename T>
std::map<size_t, T> toPolynomial(const std::string& strRepresentation, char variable) {
    std::vector<std::string> strMonoms = split(removeSpaces(strRepresentation), '+');
    std::map<size_t, T> monoms;

    for (std::string& sMonom: strMonoms) {
        std::vector<std::string> tmp = split(sMonom, variable);

        std::string coefString = tmp[0].empty() ? "1" : tmp[0];
        std::string degreeString = tmp.size() == 2 ? tmp[1].substr(1) : (sMonom[sMonom.size() - 1] == variable ? "1" : "0");

        auto coef = fromString<T>(coefString);

        if (coef != 0) {
            auto degree = fromString<size_t>(degreeString);
            monoms.emplace(degree, coef);
        }

    }

    return monoms;
}

std::vector<int> VecAdd( std::vector<int>& v1,  std::vector<int>& v2){
std::vector<int> temp(v1.size(),0);
for(int iter=0; iter<v1.size(); iter++ ){
    temp[iter] = v1[iter] + v2[iter];
}
return temp;
}

#endif //ABSTRACT_ALGEBRA2_UTILS_H
