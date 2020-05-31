#ifndef ABSTRACT_ALGEBRA2_POLYNOMIAL_H
#define ABSTRACT_ALGEBRA2_POLYNOMIAL_H


#include <vector>
#include <algorithm>
#include <ostream>
#include <stdexcept>
#include <string>
#include <map>
#include "ModArithmetic/ModArithmetic.h"
#include "utils.h"


namespace gf {
    template <typename T>
    class Polynomial {
        std::vector<T> values;
        T p_;
    public:
        Polynomial(const std::initializer_list<T>&, const T&);
        Polynomial(const std::vector<T>&, const T&);
        Polynomial(T, const T&, size_t);
        Polynomial(const std::string&, const T&, char variable=DEFAULT_VAR);
        size_t n() const;
        T p() const;
        T toInt() const;
        std::string toString(char variable=DEFAULT_VAR) const;

        template <typename X>
        friend std::ostream& operator<<(std::ostream&, const Polynomial<X>&);
        template <typename X>
        friend Polynomial<X> add(const Polynomial<X>&, const Polynomial<X>&);
        template <typename X>
        friend Polynomial<X> subtract(const Polynomial<X>&, const Polynomial<X>&);
        template <typename X>
        friend Polynomial<X> multiply(const Polynomial<X>&, const Polynomial<X>&, const Polynomial<X>&);
        template <typename X>
        friend Polynomial<X> modDivide(const Polynomial<X>&, const Polynomial<X>&);
    };

    template<typename T>
    Polynomial<T>::Polynomial(const std::initializer_list<T>& initValues, const T& p) : values(initValues), p_(p) {}

    template <typename T>
    Polynomial<T>::Polynomial(const std::vector<T>& values, const T& p) : values(values), p_(p) {}

    template <typename T>
    Polynomial<T>::Polynomial(T number, const T& mod, size_t len) : p_(mod) {
        while (number > 0) {
            T tmp = number % mod;
            number = number / mod;
            this->values.push_back(tmp);
        }

        if (values.size() > len) {
            throw std::invalid_argument("Number is to large");
        }

        for (size_t i = values.size(); i < len; ++i) {
            this->values.push_back(0);
        }

        std::reverse(this->values.begin(), this->values.end());
    }

    template <typename T>
    Polynomial<T>::Polynomial(const std::string& polyForm, const T& mod, char variable) : p_(mod) {
        std::map<size_t, T> monoms = toPolynomial<T>(polyForm, variable);

        this->values = std::vector<T>((monoms.size() > 0 ? monoms.rbegin()->first : 0) + 1, 0);

        for (auto& monom: monoms) {
            this->values[this->values.size() - monom.first - 1] = monom.second;
        }
    }

    template <typename T>
    size_t Polynomial<T>::n() const {
        return this->values.size();
    }

    template <typename T>
    T Polynomial<T>::p() const {
        return this->p_;
    }

    template <typename T>
    T Polynomial<T>::toInt() const {
        T res = 0, tmp = 1;

        for (int i = this->values.size() - 1; i >= 0; --i) {
            res = res + tmp * this->values[i];
            tmp = tmp * this->p();
        }

        return res;
    }

    template<typename T>
    std::string Polynomial<T>::toString(char variable) const {
        std::stringstream res;

        for (int i = 0; i + 2 < this->values.size(); ++i) {
            const T& coef = this->values[i];

            if (coef > 0) {
                if (res.rdbuf()->in_avail() > 0) {
                    res << " + ";
                }

                if (coef != 1) {
                    res << coef;
                }
                res << variable << '^' << this->values.size() - i - 1;
            }
        }

        const T& coefPow1 = this->values.size() >= 2 ? this->values[this->values.size() - 2] : 0;
        const T& coefPow0 = this->values.size() >= 1 ? this->values[this->values.size() - 1] : 0;

        if (coefPow1 > 0 ) {
            if (res.rdbuf()->in_avail() > 0) {
                res << " + ";
            }

            if (coefPow1 != 1) {
                res << coefPow1;
            }

            res << variable;
        }

        if (coefPow0 > 0) {
            if (res.rdbuf()->in_avail() > 0) {
                res << " + ";
            }
            res << coefPow0;
        }

        return res.str();
    }

    template <typename X>
    std::ostream &operator<<(std::ostream &os, const Polynomial<X>& polynomial) {
        os << "Polynomial({";

        for (size_t i = 0; i < polynomial.values.size(); ++i) {
            os << polynomial.values[i];

            if (i + 1 < polynomial.values.size()) {
                os << ", ";
            }
        }

        os << "}, " << polynomial.p() << ")";
        return os;
    }

    template <typename X>
    Polynomial<X> add(const Polynomial<X>& lt, const Polynomial<X>& rt) {
        if (lt.n() != rt.n() || lt.p() != rt.p()) {
            throw std::invalid_argument("Polynomials must be in the same field");
        }

        std::vector<X> values;
        for (size_t i = 0; i < lt.n(); ++i) {
            values.push_back(ModArithmetic<X>::add(lt.values[i], rt.values[i], lt.p()));
        }

        return Polynomial<X>(values, lt.p());
    }

    template <typename X>
    Polynomial<X> subtract(const Polynomial<X>& lt, const Polynomial<X>& rt) {
        if (lt.n() != rt.n() || lt.p() != rt.p()) {
            throw std::invalid_argument("Polynomials must be in the same field");
        }

        std::vector<X> values;
        for (size_t i = 0; i < lt.n(); ++i) {
            values.push_back(ModArithmetic<X>::subtract(lt.values[i], rt.values[i], lt.p()));
        }

        return Polynomial<X>(values, lt.p());
    }

    template<typename X>
    Polynomial<X> multiply(const Polynomial<X>& lt, const Polynomial<X>& rt, const Polynomial<X>& primitive) {
        if (lt.n() != rt.n() || lt.p() != rt.p()) {
            throw std::invalid_argument("Polynomials must be in the same field");
        }

        std::vector<X> values(2 * lt.n() - 1, 0);

        for (size_t i = 0; i < lt.n(); ++i) {
            for (size_t j = 0; j < rt.n(); ++j) {
                values[i + j] = values[i + j] + lt.values[i] * rt.values[j];
            }
        }

        for (X& value: values) {
            value = value % lt.p();
        }

        return modDivide(Polynomial<X>(values, lt.p()), primitive);
    }

    template<typename X>
    Polynomial<X> modDivide(const Polynomial<X>& lt, const Polynomial<X>& rt) {
        if (lt.p() != rt.p()) {
            throw std::invalid_argument("Mods are not equal");
        }

        auto tmp = lt;

        while (tmp.n() >= rt.n()) {
            while (tmp.values[0] != 0) {
                auto multiplier = tmp.values[0] / rt.values[0];

                if (tmp.values[0] % rt.values[0] != 0) {
                    multiplier = multiplier + 1;
                }

                for (size_t i = 0; i < rt.n(); ++i) {
                    tmp.values[i] = ModArithmetic<X>::subtract(tmp.values[i], multiplier * rt.values[i], lt.p());
                }
            }

            tmp.values.erase(tmp.values.begin());
        }

        return tmp;
    }
}



#endif //ABSTRACT_ALGEBRA2_POLYNOMIAL_H
