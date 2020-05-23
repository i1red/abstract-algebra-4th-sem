//
// Created by maxim on 23.05.2020.
//

#include "CyclicPolynomial.h"

CyclicPolynomial::CyclicPolynomial() {}

std::vector<int> CyclicPolynomial::calculatePolynomial(int n) {
    std::vector <int> divisors = this->findDivisors(n);
    std::vector <int> numerator(divisors[0] + 1, 0);
    std::vector <int> denominator(divisors[0] + 1, 0);

    numerator[0] = -1;
    numerator[numerator.size() - 1] = 1;

    denominator[0] = -1;
    denominator[denominator.size() - 1] = 1;

    for (auto d: divisors) {
        short int mobius = this->calculateMobiusFunction(n / d);

        if (mobius == 1) {
            std::vector <int> factor(d + 1, 0);

            factor[0] = -1;
            factor[factor.size() - 1] = 1;

            numerator = this->multiply(factor, numerator);

        } else if (mobius == -1) {

            std::vector <int> factor(d + 1, 0);

            factor[0] = -1;
            factor[factor.size() - 1] = 1;

            denominator = this->multiply(factor, denominator);
        }
    }

    return this->divisionPolynomials(numerator, denominator);
}

void CyclicPolynomial::printPoly(std::vector<int> poly) {
    int n = poly.size();

    for (int i = 0; i < n; i++) {
        if (poly[i] != 0) {
            std::cout << poly[i];

            if (i != 0) {
                std::cout << "x^" << i << " ";
            } else {
                std::cout << " ";
            }
        }
    }

    std::cout << std::endl;
}

std::vector<int> CyclicPolynomial::findDivisors(int n) {
    std::vector<int> divisors;

    for (int i = 1; i <= n; i++) {
        if (n % i == 0) {
            divisors.push_back(i);
        }
    }

    return divisors;
}

bool CyclicPolynomial::isPrime(int num) {
    bool flag = true;

    for (size_t i = 2; i <= num / 2; i++) {
        if (num % i == 0) {
            flag = false;
            break;
        }
    }

    return flag;
}

unsigned int CyclicPolynomial::splitToPrimes(unsigned int n) {
    unsigned int counter = 0;

    for(size_t i = 2; i < n; i++) {
        if (n % i == 0) {
            counter++;
            n = n / i;
            i--;
        }
    }

    return ++counter;
}

short int CyclicPolynomial::calculateMobiusFunction(unsigned int n) {
    if (n == 1) {
        return 1;
    }

    for (size_t i = 2; i <= n; i++) {
        if (isPrime(i)) {
            if (n % (int) pow(i, 2) == 0) {
                return 0;
            }
        }
    }

    return pow(-1, this->splitToPrimes(n));
}

std::vector<int> CyclicPolynomial::multiply(std::vector<int> A, std::vector<int> B) {
    int m = A.size();
    int n = B.size();
    std::vector<int> prod;

    for (int i = 0; i < m + n - 1; i++) {
        prod.push_back(0);
    }


    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            prod[i+j] += A[i] * B[j];
        }
    }

    return prod;
}

std::vector<int> CyclicPolynomial::substract(const std::vector<int> &A, const std::vector<int> &B) {
    std::vector <int> sub;

    for (int val : A) {
        sub.push_back(val);
    }

    for (int i = 0; i < B.size(); i++) {
        sub[i] -= B[i];
    }

    return sub;
}

int CyclicPolynomial::maxDegree(std::vector<int> poly) {
    for (size_t i = poly.size() - 1; i > 0; i--) {
        if (poly[i] != 0) {
            return i;
        }
    }
}

int CyclicPolynomial::koefOfMaxDegree(std::vector<int> poly) {
    for (size_t i = poly.size() - 1; i > 0; i--) {
        if (poly[i] != 0) {
            return poly[i];
        }
    }
}

std::vector<int> CyclicPolynomial::divisionPolynomials(std::vector<int> numerator, std::vector<int> denominator) {
    size_t len = numerator.size() - denominator.size();
    std::vector <int> result(len + 1, 0);

    result[result.size() - 1] = (numerator[numerator.size() - 1] / denominator[denominator.size() - 1]);

    numerator = this->substract(numerator, this->multiply(result, denominator));

    while (this->maxDegree(numerator) >= this->maxDegree(denominator)) {
        result[this->maxDegree(numerator) - this->maxDegree(denominator)] = this->koefOfMaxDegree(numerator) / this->koefOfMaxDegree(denominator);

        size_t len = this->maxDegree(numerator) - this->maxDegree(denominator);

        std::vector <int> temp(len + 1, 0);

        temp[temp.size() - 1] = (this->koefOfMaxDegree(numerator) / this->koefOfMaxDegree(denominator));

        numerator = this->substract(numerator, this->multiply(temp, denominator));
    }

    return result;
}

CyclicPolynomial::~CyclicPolynomial() = default;
