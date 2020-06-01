#include <iostream>
#include "GaloisField/Polynomial.h"
#include "CyclicPolynomial/CyclicPolynomial.h"

namespace test {
    std::vector<std::vector<int>> f9mul = {{0, 0, 0, 0, 0, 0, 0, 0, 0},
                                           {0, 1, 2, 3, 4, 5, 6, 7, 8},
                                           {0, 2, 1, 6, 8, 7, 3, 5, 4},
                                           {0, 3, 6, 2, 5, 8, 1, 4, 7},
                                           {0, 4, 8, 5, 6, 1, 7, 2, 3},
                                           {0, 5, 7, 8, 1, 3, 4, 6, 2},
                                           {0, 6, 3, 1, 7, 4, 2, 8, 5},
                                           {0, 7, 5, 4, 2, 6, 8, 3, 1},
                                           {0, 8, 4, 7, 3, 2, 5, 1, 6}};

    gf::Polynomial<int> f9primitive("x^2 + 1", 3, 3);

    std::vector<std::vector<int>> f8mul = {{0, 0, 0, 0, 0, 0, 0, 0},
                                           {0, 1, 2, 3, 4, 5, 6, 7},
                                           {0, 2, 4, 6, 3, 1, 7, 5},
                                           {0, 3, 6, 5, 7, 4, 1, 2},
                                           {0, 4, 3, 7, 6, 2, 5, 1},
                                           {0, 5, 1, 4, 2, 7, 3, 6},
                                           {0, 6, 7, 1, 5, 3, 2, 4},
                                           {0, 7, 5, 2, 1, 6, 4, 3}};

    gf::Polynomial<int> f8primitive("x^3 + x + 1", 2, 4);


    void testMul(int p, int n, std::vector<std::vector<int>> table, const gf::Polynomial<int> primitive) {
        std::cout << "TEST multiplication p=" << p << ", n=" << n << std::endl;

        int bound = intPow(p, n);
        int successfulTest = 0;

        for (int i = 0; i < bound; ++i) {
            for (int j = 0; j < bound; ++j) {
                int mulRes = gf::multiply(gf::Polynomial<int>(i, p, n), gf::Polynomial<int>(j, p, n), primitive).toInt();
                if (mulRes == table[i][j]) {
                    successfulTest += 1;
                }
                else {
                    std::cout << "FAILED multiplication: " << i << " * " << j << " returned " << mulRes << ". Should return " << table[i][j] << std::endl;
                }
            }
        }

        std::cout << "PASSED " << successfulTest << "/" << bound * bound << " TESTS" << std::endl;
    }

    void run() {
        testMul(3, 2, f9mul, f9primitive);
        testMul(2, 3, f8mul, f8primitive);
    }

    void testCalculatingCyclicPolynomial() { /* You can compare these results with results on Wikipedia (cyclic polynomial) */
        auto cyclicPolinomial = new CyclicPolynomial();

        cyclicPolinomial->printPoly(cyclicPolinomial->calculatePolynomial(3));
        cyclicPolinomial->printPoly(cyclicPolinomial->calculatePolynomial(4));
        cyclicPolinomial->printPoly(cyclicPolinomial->calculatePolynomial(5));
        cyclicPolinomial->printPoly(cyclicPolinomial->calculatePolynomial(6));
        cyclicPolinomial->printPoly(cyclicPolinomial->calculatePolynomial(7));
        cyclicPolinomial->printPoly(cyclicPolinomial->calculatePolynomial(8));
        cyclicPolinomial->printPoly(cyclicPolinomial->calculatePolynomial(9));
        cyclicPolinomial->printPoly(cyclicPolinomial->calculatePolynomial(10));
        cyclicPolinomial->printPoly(cyclicPolinomial->calculatePolynomial(11));
        cyclicPolinomial->printPoly(cyclicPolinomial->calculatePolynomial(12));
    }
    void DividePoly(){
        gf:: Polynomial<int> first("6x^5 + 3x^4 + 6x^3 + 5x^2 + 6x + 3" , 7, 6);
        gf::Polynomial<int>second("4x^3 + x^2 + 3x + 1", 7,6);
        gf::Polynomial<int> f8primitive("x^7 + x^4 + 3x^3 + x + 5", 7, 7);

        gf:: Polynomial<int> third("3x^2 + 3x + 2" , 5, 3);
        gf::Polynomial<int>fourth("2x + 3", 5,3);
        gf::Polynomial<int> f9primitive("x^3 + 3x + 2", 5, 4);

        gf:: Polynomial<int> sixth("x^5+x^3+2x^2+x+2" , 3, 6);
        gf::Polynomial<int>seventh("x^3+2x^2+x+1", 3,6);
        gf::Polynomial<int> f7primitive("x^8+x^4+2x^3+x+2", 3, 7);

        std::pair<gf::Polynomial<int>,gf::Polynomial<int>>res1 = gf::divide(third,fourth,f9primitive);
        std::cout<<"Quotient: " << res1.first  << " Remainder: " << res1.second << std::endl;

        std::pair<gf::Polynomial<int>,gf::Polynomial<int>>res2 = gf::divide(first,second,f8primitive);
        std::cout<<"Quotient: " << res2.first  << " Remainder: " << res2.second << std::endl;

        std::pair<gf::Polynomial<int>,gf::Polynomial<int>>res3 = gf::divide(sixth,seventh,f7primitive);
        std::cout<<"Quotient: " << res3.first  << " Remainder: " << res3.second << std::endl;

    }

}


int main() {

    //test::testCalculatingCyclicPolynomial();
    //test::run();
    test::DividePoly();



    return 0;
}
