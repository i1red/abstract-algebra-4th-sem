#include <iostream>
#include "GaloisField/GaloisField.h"
#include "CyclicPolynomial/CyclicPolynomial.h"
#include "CountRoots/CountRoots.h"

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

    gf::Vector<int> f9primitive("x^2 + 1", 3);

    std::vector<std::vector<int>> f8mul = {{0, 0, 0, 0, 0, 0, 0, 0},
                                           {0, 1, 2, 3, 4, 5, 6, 7},
                                           {0, 2, 4, 6, 3, 1, 7, 5},
                                           {0, 3, 6, 5, 7, 4, 1, 2},
                                           {0, 4, 3, 7, 6, 2, 5, 1},
                                           {0, 5, 1, 4, 2, 7, 3, 6},
                                           {0, 6, 7, 1, 5, 3, 2, 4},
                                           {0, 7, 5, 2, 1, 6, 4, 3}};

    gf::Vector<int> f8primitive("x^3 + x + 1", 2);


    void testMul(int q, int n, std::vector<std::vector<int>> table, const gf::Vector<int> primitive) {
        std::cout << "TEST multiplication q=" << q << ", n=" << n << std::endl;

        int bound = intPow(q, n);
        int successfulTest = 0;

        for (int i = 0; i < bound; ++i) {
            for (int j = 0; j < bound; ++j) {
                int mulRes = gf::vecPolMod(gf::vecPolMul(gf::Vector<int>(i, q, n), gf::Vector<int>(j, q, n)), primitive).toInt();
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

}

namespace ui {
    void help() {
        std::cout << "Enter command with the needed arguments from the list below\n";
        std::cout << "\t->smth ARG1 ARG2\n";
        std::cout << "\t->smth ARG1 ARG2\n";
        std::cout << "\t->quit\n";
        std::cout << "\t\tfor exit\n";
    }

    void execute(std::string command) {
        double arg1, arg2;
        std::cin >> arg1 >> arg2;
        if (command == "add") std::cout << arg1 + arg2;
        else std::cout << "error";
    }

    void interactionLoop() {
        std::string command;
        ui::help();

        while (true) {
            std::cout << "->";
            std::cin >> command;
            if (command == "quit") break;
            else execute(command);
        }

        std::cout << "bye";
    };
}


int main() {

//    test::testCalculatingCyclicPolynomial();
    //test::run();
    ui::interactionLoop();
    return 0;
}
