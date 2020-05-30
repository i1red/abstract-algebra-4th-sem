#pragma once
#include <iostream>
#include "GaloisField/Polynomial.h"
#include "GaloisField/Vector.h"
#include "ModArithmetic/ModArithmetic.h"
#include "utils.h"
#include "string.h"
template <typename T>
gf::Polynomial<T> MakeMonic(gf::Polynomial<T> p) {
	std::map<size_t, T> pol = toPolynomial<T>(p.toString() , 'x');
	T max_coef = pol[p.n()];
	T q = p.q();
	if (max_coef == 1) {
		return p;
	}
	else {
		for (auto& item : pol) {
				item.second = ModArithmetic<T>::divide(item.second, max_coef, q);
		}
		return gf::Polynomial<T>(pol, q, p.n());
	}
}

/* template <typename T>
gf::Polynomial<T> MakeMonic(gf::Polynomial<T> p) {
	std::map<size_t, gf::Vector<T>> pol = toPolynomial<T>(p.toString(), 'x');
	gf::Vector<T> max_coef = pol[p.n()];
	T q = p.q();
	for (auto& item : pol) {
		item.second = gf::vecPolMod(item.second, max_coef);
	}
	return gf::Polynomial<T>(pol, q, p.n());
}
*/

template <typename T>
gf::Polynomial<T> Derivative(gf::Polynomial<T> p) {
	std::map<size_t, T> pol = toPolynomial<T>(p.toString(), 'x');
	std::map<size_t, T> derivative;
	T q = p.q();
	for (auto& item : pol) {
		if (item.first == 0) continue;
		derivative[item.first - 1] = ModArithmetic<T>::multiply(item.first, item.second, q);
	}
	return gf::Polynomial<T>(derivative, q, p.n()-1);
}

template <typename T>
T PointValue(gf::Polynomial<T> p, T x) {
	std::map<size_t, T> pol = toPolynomial<T>(p.toString(), 'x');
	T result=0;
	T q = p.q();
	for (auto& item : pol) {
		result = ModArithmetic<T>::add(result, ModArithmetic<T>::multiply(item.second, ModularPow(x, item.first, q),q), q);
	}
	return result;
}

template <typename T>
long int ModularPow(T base, long int degree,T mod) {
	if (degree == 0) return 1;
	T res = base;
	for (long int i = 1; i < degree; i++) {
		res = ModArithmetic<T>::multiply(res, base, mod);
	}
	return res;
}