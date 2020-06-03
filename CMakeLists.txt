cmake_minimum_required(VERSION 3.15)
project(abstract_algebra2)

set(CMAKE_CXX_STANDARD 17)


include_directories(
        ${PROJECT_SOURCE_DIR}
        ${BigInt}
        ${ModArithmetic}
        ${MontgomeryReduction}
)

add_executable(abstract_algebra2 main.cpp
               BigInt/BigInt.cpp BigInt/BigInt.h
               ModArithmetic/ModArithmetic.h
               BigInt/bigmath.h
        GaloisField/Polynomial.h
        utils.h CyclicPolynomial/CyclicPolynomial.cpp CyclicPolynomial/CyclicPolynomial.h)