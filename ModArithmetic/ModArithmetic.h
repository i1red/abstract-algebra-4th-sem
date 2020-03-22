#ifndef ABSTRACT_ALGEBRA_MODARITHMETIC_H
#define ABSTRACT_ALGEBRA_MODARITHMETIC_H


#include <stdexcept>
#include "BigInt/BigInt.h"
#include "BigInt/bigmath.h"


template <typename T>
class ModArithmetic {
    T modulo_;
public:
    static T add(const T&, const T&, const T&);
    static T subtract(const T&, const T&, const T&);
    static T multiply(const T&, const T&, const T&);
    static T divide(const T&, const T&, const T&);
    static T inverseEl(const T&, const T&);
    static T gcdex(const T&, const T&, T&, T&, const T&);
    explicit ModArithmetic(const T&);
    T modulo() const;
    T ringAdd(const T&, const T&) const;
    T ringSubtract(const T&, const T&) const;
    T ringMultiply(const T&, const T&) const;
    T ringInverseEl(const T&) const;
    T ringDivide(const T&, const T&) const;
};


template <typename T>
T ModArithmetic<T>::add(const T& lt, const T& rt, const T& modul) {
    return mod(lt + rt, modul);
}

template <typename T>
T ModArithmetic<T>::subtract(const T& lt, const T& rt, const T& modul) {
    return mod(lt - rt, modul);
}

template <typename T>
T ModArithmetic<T>::multiply(const T& lt, const T& rt, const T& modul) {
    return mod(lt * rt, modul);
}

template <typename T>
T ModArithmetic<T>::divide(const T& lt, const T& rt, const T& modul) {
    T inverseRt = inverseEl(rt, modul);

    if (inverseRt == -1) {
        throw std::logic_error("Division error. Inverse element for right operand doesn't exist");
    }

    return multiply(lt, inverseRt, modul);
}

template <typename T>
T ModArithmetic<T>::inverseEl(const T& el, const T& modul) {
    T x, y, res;
    T g = gcdex (el, modul, x, y, modul);

    if (!(g == 1)) {
        res = -1;
    }

    else {
        res = (x % modul + modul) % modul;
    }

    return res;
}


template <typename T>
inline ModArithmetic<T>::ModArithmetic(const T& n) {
    this->modulo_ = n;
}

template <typename T>
inline T ModArithmetic<T>::modulo() const {
    return this->modulo_;
}

template <typename T>
inline T ModArithmetic<T>::ringAdd(const T& lt, const T& rt) const {
    return add(lt, rt, this->modulo_);
}

template <typename T>
inline T ModArithmetic<T>::ringSubtract(const T& lt, const T& rt) const {
    return subtract(lt, rt, this->modulo_);
}

template <typename T>
inline T ModArithmetic<T>::ringMultiply(const T& lt, const T& rt) const {
    return multiply(lt, rt, this->modulo_);
}

template <typename T>
inline T ModArithmetic<T>::ringDivide(const T& lt, const T& rt) const {
    return divide(lt, rt, this->modulo_);
}

template <typename T>
inline T ModArithmetic<T>::ringInverseEl(const T & el) const {
    return inverseEl(el, this->modulo_);
}

template <typename T>
T ModArithmetic<T>::gcdex(const T &lt, const T &rt, T &x, T &y, const T& modul) {

    if (lt == 0) {
        x = 0;
        y = 1;

        return rt;
    }

    T x1, y1, mod_lt;

    mod_lt = mod(lt, modul);
    T d = gcdex(rt % mod_lt, mod_lt, x1, y1, modul);

    x = y1 - (rt / mod_lt) * x1;
    y = x1;

    return d;
}

#endif //ABSTRACT_ALGEBRA_MODARITHMETIC_H