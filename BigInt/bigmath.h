#ifndef ABSTRACT_ALGEBRA2_BIGMATH_H
#define ABSTRACT_ALGEBRA2_BIGMATH_H


#include "BigInt.h"


template <typename T>
T mod(const T& lt, const T& rt) {
    T sign = 1, left = lt, right = rt;
    if (rt < 0) {
        sign = -1;
        right = -rt;
        left = -lt;
    }


    while (left < 0) {
        T lenDif = std::log10(std::abs(left)) - std::log10(right);
        left = left + right * (lenDif > 0 ? std::pow(10, lenDif) : 1);
    }

    return sign * (left % right);
}

template <>
inline BigInt mod(const BigInt& lt, const BigInt& rt) {
    return lt.mod(rt);
}

template <typename T>
T intPow(const T& lt, const T& rt) {
    if (rt < 0) {
        throw std::invalid_argument("Degree must be greater or equal 0");
    }

    return (T)pow(lt, rt);
}

template <>
inline BigInt intPow(const BigInt& lt, const BigInt& rt) {
    return lt.power(rt);
}


#endif //ABSTRACT_ALGEBRA2_BIGMATH_H
