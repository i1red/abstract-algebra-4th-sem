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

}

namespace ui {

    int p = 11, n = 20;

    void newLine() {
        std::cout << "->";
    }

    void help() {
        std::cout << "Quick guide for using our software:\n";
        std::cout << "\tFirst of all, you need to enter value of field in which you will work,\n";
        std::cout << "\tthen order of cyclic polynomial, decomposition of which will extend field,\n";
        std::cout << "\tafter this steps you must choose one of provided polynomials who will extend\n";
        std::cout << "\tgiven field. Finally, you can use commands from the list below\n";
        std::cout << "\tfor executing needed operations\n\n";
        std::cout << "Existing commands:\n";
        std::cout << "\t->enter\n";
        std::cout << "\t->add\n";
        std::cout << "\t->subtract\n";
        std::cout << "\t->multiply\n";
        std::cout << "\t->derivative\n";
        std::cout << "\t->monic\n";
        std::cout << "\t->point value\n";
        std::cout << "\t->roots\n";
        std::cout << "\t->amount of roots\n";
        std::cout << "\t->quit\n";
        std::cout << "\t\tfor exit\n";
    }

    int readIntValue(const std::string& textToShow) {
        int buffer;
        std::cout << textToShow << std::endl;
        newLine();
        std::cin >> buffer;

        while (std::cin.fail()) {
            std::cin.clear(); // clear input buffer to restore cin to a usable state
            std::cin.ignore(INT_MAX, '\n'); // ignore last input
            std::cout << "You can only enter numbers." << std::endl;
            std::cout << textToShow << std::endl;
            std::cin >> buffer;
        }

        return buffer;
    }

    std::string readStringValue(const std::string& textToShow) {
        std::string buffer;
        std::cout << textToShow << std::endl;
        newLine();
        std::cin >> buffer;

        while (std::cin.fail()) {
            std::cin.clear(); // clear input buffer to restore cin to a usable state
            std::cin.ignore(INT_MAX, '\n'); // ignore last input
            std::cout << "You can only enter numbers." << std::endl;
            std::cout << textToShow << std::endl;
            std::cin >> buffer;
        }

        return buffer;
    }

    gf::Polynomial<int> getPolynomialByInput(std::string text) {
        return gf::Polynomial<int>(readStringValue(text), p, n);
    }

    void execute(const std::string& command) {
        if (command == "enter") {
            int arg1, arg2;
            arg1 = readIntValue("Enter the field");
            arg2 = readIntValue("Enter the n");
            std::cout << arg1 + arg2 << std::endl;

        } else if (command == "add") {
            gf::Polynomial<int> value1 = getPolynomialByInput("Enter first polynomial for addition");
            gf::Polynomial<int> value2 = getPolynomialByInput("Enter second polynomial for addition");
            value1 = gf::add(value1, value2);
            std::cout << value1.toString() << std::endl;

        } else if (command == "subtract") {
            gf::Polynomial<int> value1 = getPolynomialByInput("Enter first polynomial for subtracting");
            gf::Polynomial<int> value2 = getPolynomialByInput("Enter second polynomial for subtracting");
            value1 = gf::subtract(value1, value2);
            std::cout << value1.toString() << std::endl;

        } else if (command == "multiply") {
            gf::Polynomial<int> value1 = getPolynomialByInput("Enter first polynomial for multiplying");
            gf::Polynomial<int> value2 = getPolynomialByInput("Enter second polynomial for multiplying");
//            value1 = gf::multiply(value1, value2, );
            std::cout << value1.toString() << std::endl;

        } else if (command == "derivative") {
            gf::Polynomial<int> value = getPolynomialByInput("Enter polynomial for getting derivative");
            value = gf::Derivative(value);
            std::cout << value.toString() << std::endl;

        } else if (command == "monic") {
            gf::Polynomial<int> value = getPolynomialByInput("Enter polynomial for getting monic");
            value = gf::MakeMonic(value);
            std::cout << value.toString() << std::endl;

        } else if (command == "point value") {
            gf::Polynomial<int> value = getPolynomialByInput("Enter polynomial for getting point value of this polynomial");
//            value = gf::PointValue(value, );
            std::cout << value.toString() << std::endl;

        } else if (command == "roots") {
            std::string pol;
//            pol = readStringValue("Enter polynomial for getting point value of this polynomial");
//            gf::Polynomial<int> value = gf::Polynomial<int>(pol, p, n);
//            value = gf::PointValue(value, );
//            std::cout << value.toString() << std::endl;

        } else if (command == "amount of roots") {
            std::string pol;
//            pol = readStringValue("Enter polynomial for getting point value of this polynomial");
//            gf::Polynomial<int> value = gf::Polynomial<int>(pol, p, n);
//            value = gf::PointValue(value, );
//            std::cout << value.toString() << std::endl;

        } else {
            std::cout << "error" << std::endl;
        }
    }

    void interactionLoop() {
        std::string command;

        while (true) {
            std::cout << "Enter the command or help for getting information" << std::endl;
            std::cout << "->";
            std::cin >> command;
            if (command == "quit") break;
            else if (command == "help") ui::help();
            else execute(command);
        }

        std::cout << "Goodbye";
    };
}


int main() {

//    test::testCalculatingCyclicPolynomial();
    //test::run();
    ui::interactionLoop();
    return 0;
}
