#ifndef ABSTRACT_ALGEBRA2_POLYNOMIAL_H
#define ABSTRACT_ALGEBRA2_POLYNOMIAL_H


#include <map>
#include "Vector.h"
#include "utils.h"


namespace gf {
    template <typename T>
    class Polynomial {
        std::map<size_t, Vector<T>> monoms;
        T q_;
        size_t n_;
    public:
        Polynomial(const std::map<size_t, T>&, const T&, size_t);
        Polynomial(const std::map<size_t, Vector<T>>&, const T&, size_t);
        Polynomial(const std::string& strRepr, const T& q, size_t n, char variable=DEFAULT_VAR);

        std::string toString(char variable=DEFAULT_VAR) const;
        T q() const;
        size_t n() const;

        template <typename X>
        friend std::ostream& operator<<(std::ostream& os, const Polynomial<X>&);
        template <typename X>
        friend Polynomial<X> polAdd(const Polynomial<X>&, const Polynomial<X>&);
        template <typename X>
        friend Polynomial<X> polSubtract(const Polynomial<X>&, const Polynomial<X>&);
        template <typename X>
        friend Polynomial<X> polMultiply(const Polynomial<X>&, const Polynomial<X>&, const Vector<X>&);
    };

    template<typename T>
    Polynomial<T>::Polynomial(const std::map<size_t, T>& monoms, const T& q, size_t n) : q_(q), n_(n) {
        for(auto& item: monoms) {
            this->monoms.emplace(item.first, Vector<T>(item.second, q, n));
        }
    }

    template <typename T>
    Polynomial<T>::Polynomial(const std::map<size_t, Vector<T>>& monoms, const T& q, size_t n) : monoms(monoms), q_(q), n_(n) {}

    template <typename T>
    Polynomial<T>::Polynomial(const std::string &strRepr, const T &q, size_t n, char variable) :
            Polynomial(toPolynomial<T>(strRepr, variable), q, n) {}

    template<typename T>
    std::string Polynomial<T>::toString(char variable) const {
        std::stringstream res;

        for (auto it = this->monoms.rbegin(); it != this->monoms.rend(); it++) {
            T coef = it->second.toInt();
            size_t degree = it->first;

            bool printDegree = degree > 1, printVar = degree != 0;

            if (!printVar || coef > 1) {
                res << coef;
            }
            if (printDegree) {
                res << variable << '^' << degree;
            }
            else if (printVar) {
                res << variable;
            }

            auto next = it;
            next++;

            if (next != this->monoms.rend()) {
                res << " + ";
            }
        }

        return res.str();
    }

    template<typename T>
    T Polynomial<T>::q() const {
        return this->q_;
    }

    template<typename T>
    size_t Polynomial<T>::n() const {
        return this->n_;
    }

    template<typename X>
    std::ostream &operator<<(std::ostream &os, const Polynomial<X>& pol) {
        os << "Polynomial(" << pol.toString()  << ", " << pol.q() << ", " << pol.n() << ", x)";
        return os;
    }

    template <typename X>
    Polynomial<X> polAdd(const Polynomial<X>& lt, const Polynomial<X>& rt) {
        if (lt.q() != rt.q() || rt.n() != lt.n()) {
            throw std::invalid_argument("Polynomials must be in the same field");
        }

        std::map<size_t, Vector<X>> init_map = lt.monoms;

        for (auto& monom: rt.monoms) {
            if (init_map.count(monom.first) > 0) {
                init_map.at(monom.first) = vecAdd(init_map.at(monom.first), monom.second);;
            }
            else {
                init_map.insert(monom);
            }

            if (init_map.at(monom.first).toInt() == 0) {
                init_map.erase(monom.first);
            }
        }

        return Polynomial<X>(init_map, lt.q(), lt.n());
    }

    template <typename X>
    Polynomial<X> polSubtract(const Polynomial<X>& lt, const Polynomial<X>& rt) {
        if (lt.q() != rt.q() || rt.n() != lt.n()) {
            throw std::invalid_argument("Polynomials must be in the same field");
        }

        std::map<size_t, Vector<X>> init_map = lt.monoms;

        for (auto& monom: rt.monoms) {
            if (init_map.count(monom.first) > 0) {
                init_map.at(monom.first) = vecSubtract(init_map.at(monom.first), monom.second);
            }
            else {
                init_map.emplace(monom.first, vecSubtract(Vector<X>(0, monom.second.mod(), monom.second.length()), monom.second));
            }

            if (init_map.at(monom.first).toInt() == 0) {
                init_map.erase(monom.first);
            }
        }

        return Polynomial<X>(init_map, lt.q(), lt.n());
    }

    template <typename X>
    Polynomial<X> polMultiply(const Polynomial<X>& lt, const Polynomial<X>& rt, const Vector<X>& primitive) {
        if (lt.q() != rt.q() || rt.n() != lt.n()) {
            throw std::invalid_argument("Polynomials must be in the same field");
        }

        if (primitive.length() != lt.q() + 2) {
            std::stringstream msg;

            msg << "Primitive degree must be " << lt.q() + 1 << " for field with prime " << lt.q()
                << ". Got " << primitive.length() - 1 << std::endl;

            throw std::logic_error(msg.str());
        }

        std::map<size_t, Vector<X>> init_map;

        for (const auto &ltMonom: lt.monoms) {
            for (const auto &rtMonom: rt.monoms) {
                size_t degree = ltMonom.first + rtMonom.first;
                Vector<X> coef = vecPolMod(vecPolMul(ltMonom.second, rtMonom.second), primitive);

                if (init_map.count(degree) > 0) {
                    init_map.at(degree) = vecAdd(init_map.at(degree), coef);
                } else {
                    init_map.emplace(degree, coef);
                }

                if (init_map.at(degree).toInt() == 0) {
                    init_map.erase(degree);
                }
            }
        }

        return Polynomial<X>(init_map, lt.q(), lt.n());
    }

}


#endif //ABSTRACT_ALGEBRA2_POLYNOMIAL_H
