#ifndef ABSTRACT_ALGEBRA2_VECTOR_H
#define ABSTRACT_ALGEBRA2_VECTOR_H


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
    class Vector {
        std::vector<T> values;
        T mod_;
    public:
        Vector(const std::initializer_list<T>&, const T&);
        Vector(const std::vector<T>&, const T&);
        Vector(T, const T&, size_t);
        Vector(const std::string&, const T&, char variable=DEFAULT_VAR);
        size_t length() const;
        T mod() const;
        T toInt() const;
        std::string toPolString(char variable=DEFAULT_VAR) const;

        template <typename X>
        friend std::ostream& operator<<(std::ostream& os, const Vector<X>&);
        template <typename X>
        friend Vector<X> vecAdd(const Vector<X> &lt, const Vector<X> &rt);
        template <typename X>
        friend Vector<X> vecSubtract(const Vector<X> &lt, const Vector<X> &rt);
        template <typename X>
        friend Vector<X> vecPolMul(const Vector<X> &lt, const Vector<X> &rt);
        template <typename X>
        friend Vector<X> vecPolMod(const Vector<X> &lt, const Vector<X> &rt);
    };

    template<typename T>
    Vector<T>::Vector(const std::initializer_list<T>& initValues, const T& mod) : values(initValues),  mod_(mod) {}

    template <typename T>
    Vector<T>::Vector(const std::vector<T>& values, const T& mod) : values(values), mod_(mod) {}

    template <typename T>
    Vector<T>::Vector(T number, const T& mod, size_t len) : mod_(mod) {
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
    Vector<T>::Vector(const std::string& polyForm, const T& mod, char variable) : mod_(mod) {
        std::map<size_t, T> monoms = toPolynomial<T>(polyForm, variable);

        this->values = std::vector<T>((monoms.size() > 0 ? monoms.rbegin()->first : 0) + 1, 0);

        for (auto& monom: monoms) {
            this->values[this->values.size() - monom.first - 1] = monom.second;
        }
    }

    template <typename T>
    size_t Vector<T>::length() const {
        return this->values.size();
    }

    template <typename T>
    T Vector<T>::mod() const {
        return this->mod_;
    }

    template <typename T>
    T Vector<T>::toInt() const {
        T res = 0, tmp = 1;

        for (int i = this->values.size() - 1; i >= 0; --i) {
            res = res + tmp * this->values[i];
            tmp = tmp * this->mod();
        }

        return res;
    }

    template<typename T>
    std::string Vector<T>::toPolString(char variable) const {
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
    std::ostream &operator<<(std::ostream &os, const Vector<X>& vec) {
        os << "Vector({";

        for (size_t i = 0; i < vec.values.size(); ++i) {
            os << vec.values[i];

            if (i + 1 < vec.values.size()) {
                os << ", ";
            }
        }

        os << "}, " << vec.mod() << ")";
        return os;
    }

    template <typename X>
    Vector<X> vecAdd(const Vector<X>& lt, const Vector<X>& rt) {
        if (lt.length() != rt.length() || lt.mod() != rt.mod()) {
            throw std::invalid_argument("Vectors must be in the same field");
        }

        std::vector<X> values;
        for (size_t i = 0; i < lt.length(); ++i) {
            values.push_back(ModArithmetic<X>::add(lt.values[i], rt.values[i], lt.mod()));
        }

        return Vector<X>(values, lt.mod());
    }

    template <typename X>
    Vector<X> vecSubtract(const Vector<X>& lt, const Vector<X>& rt) {
        if (lt.length() != rt.length() || lt.mod() != rt.mod()) {
            throw std::invalid_argument("Vectors must be in the same field");
        }

        std::vector<X> values;
        for (size_t i = 0; i < lt.length(); ++i) {
            values.push_back(ModArithmetic<X>::subtract(lt.values[i], rt.values[i], lt.mod()));
        }

        return Vector<X>(values, lt.mod());
    }

    template<typename X>
    Vector<X> vecPolMul(const Vector<X>& lt, const Vector<X>& rt) {
        if (lt.length() != rt.length() || lt.mod() != rt.mod()) {
            throw std::invalid_argument("Vectors must be in the same field");
        }

        std::vector<X> values(2 * lt.length() - 1, 0);

        for (size_t i = 0; i < lt.length(); ++i) {
            for (size_t j = 0; j < rt.length(); ++j) {
                values[i + j] = values[i + j] + lt.values[i] * rt.values[j];
            }
        }

        for (X& value: values) {
            value = value % lt.mod();
        }

        return Vector<X>(values, lt.mod());
    }

    template<typename X>
    Vector<X> vecPolMod(const Vector<X>& lt, const Vector<X>& rt) {
        if (lt.mod() != rt.mod()) {
            throw std::invalid_argument("Mods are not equal");
        }

        auto tmp = lt;

        while (tmp.length() >= rt.length()) {
            while (tmp.values[0] != 0) {
                auto multiplier = tmp.values[0] / rt.values[0];

                if (tmp.values[0] % rt.values[0] != 0) {
                    multiplier = multiplier + 1;
                }

                for (size_t i = 0; i < rt.length(); ++i) {
                    tmp.values[i] = ModArithmetic<X>::subtract(tmp.values[i], multiplier * rt.values[i], lt.mod());
                }
            }

            tmp.values.erase(tmp.values.begin());
        }

        return tmp;
    }

    template <typename T>
    class VectorFactory {
        T mod;
        size_t len;
    public:
        VectorFactory(const T&, size_t);
        Vector<T> fromInt(T);
    };

    template <typename T>
    VectorFactory<T>::VectorFactory(const T& mod, size_t len) : mod(mod), len(len) {}

    template <typename T>
    Vector<T> VectorFactory<T>::fromInt(T number) {
        return Vector<T>(number, this->mod, this->len);
    }
}



#endif //ABSTRACT_ALGEBRA2_VECTOR_H
