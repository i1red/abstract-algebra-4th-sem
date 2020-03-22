#ifndef ABSTRACT_ALGEBRA2_GALOISFIELD_H
#define ABSTRACT_ALGEBRA2_GALOISFIELD_H


#include "Polynomial.h"
#include "Vector.h"


namespace gf {
    template <typename T>
    class GaloisField {
        T prime_;
        size_t n_;
        char variable_;
        Vector<T> vecPrimitive_;
    public:
        GaloisField(const T&, size_t, char, const std::string& primitive="");
        T prime() const;
        size_t n() const;
        char variable() const;
        std::string primitive() const;
        void setPrimitive(std::string primitive);
        std::string add(const std::string&, const std::string&) const;
        std::string subtract(const std::string&, const std::string&) const;
        std::string multiply(const std::string&, const std::string&) const;
    };

    template <typename T>
    GaloisField<T>::GaloisField(const T& prime, size_t n, char variable, const std::string& primitive) :
    vecPrimitive_(primitive, prime, variable), prime_(prime), n_(n), variable_(variable) {}

    template<typename T>
    T GaloisField<T>::prime() const {
        return this->prime_;
    }

    template<typename T>
    size_t GaloisField<T>::n() const {
        return this->n_;
    }

    template<typename T>
    char GaloisField<T>::variable() const {
        return this->variable_;
    }

    template<typename T>
    std::string GaloisField<T>::primitive() const {
        return this->vecPrimitive_.toPolString(this->variable());
    }

    template <typename T>
    void GaloisField<T>::setPrimitive(std::string primitive) {
        this->vecPrimitive_ = Vector<T>(primitive, this->prime(), this->variable());
    }

    template<typename T>
    std::string GaloisField<T>::add(const std::string& lt, const std::string& rt) const {
        return polAdd(Polynomial<T>(lt, this->prime(), this->n(), this->variable()),
                      Polynomial<T>(rt, this->prime(), this->n(), this->variable())).toString();
    }

    template<typename T>
    std::string GaloisField<T>::subtract(const std::string& lt, const std::string& rt) const {
        return polSubtract(Polynomial<T>(lt, this->prime(), this->n(), this->variable()),
                           Polynomial<T>(rt, this->prime(), this->n(), this->variable())).toString();
    }

    template<typename T>
    std::string GaloisField<T>::multiply(const std::string& lt, const std::string& rt) const {
        return polMultiply(Polynomial<T>(lt, this->prime(), this->n(), this->variable()),
                           Polynomial<T>(rt, this->prime(), this->n(), this->variable()),
                           this->vecPrimitive_).toString();
    }

}


#endif //ABSTRACT_ALGEBRA2_GALOISFIELD_H
