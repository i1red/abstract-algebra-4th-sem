//
// Created by maxim on 23.05.2020.
//

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

#ifndef ABSTRACT_ALGEBRA2_CYCLICPOLYNOMIAL_H
#define ABSTRACT_ALGEBRA2_CYCLICPOLYNOMIAL_H


class CyclicPolynomial {
public:
    CyclicPolynomial();
    ~CyclicPolynomial();
    std::vector <int> calculatePolynomial(int n); /* Main method (one parameter: degree of cyclic polynomial;
        * return array : index - degree, value[i] - koef) */
    void printPoly(std:: vector<int> poly); /* Helper function to print polynom in convenient format */

private:
    std::vector <int> findDivisors(int n); /* Find all positive divisors of number */
    bool isPrime(int num); /* Check if the current number is prime */
    unsigned int splitToPrimes(unsigned int n); /* Calculates the number of prime divisors of a number */
    short int calculateMobiusFunction(unsigned int n); /* Calculates value (-1, 0, 1) according to Mobius function */
    std::vector <int> multiply(std::vector <int> A, std::vector <int> B); /* Multiply two polynomials */
    std::vector <int> substract(const std::vector <int>& A, const std::vector <int>& B); /* Substract two polynomials */
    int maxDegree(std::vector <int> poly); /* Find max degree in polynomial */
    int koefOfMaxDegree(std::vector <int> poly); /* Find koef near with max degree in polynomial */
    std::vector <int> divisionPolynomials(std::vector <int> numerator, std::vector <int> denominator); /* Division two polynomials */
};


#endif //ABSTRACT_ALGEBRA2_CYCLICPOLYNOMIAL_H
