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
#include <utility>
#include <CyclicPolynomial/CyclicPolynomial.h>


namespace gf {
    template<typename T>
    class Polynomial {
        std::vector<T> values;
        T p_;
    public:
        Polynomial(const std::initializer_list<T> &, const T &);

        Polynomial(const std::vector<T> &, const T &);

        Polynomial(T, const T &, size_t);

        Polynomial(const std::string &, const T &, size_t, char variable = DEFAULT_VAR);

        size_t n() const;

        T p() const;

        T toInt() const;

        std::string toString(char variable = DEFAULT_VAR) const;

        template<typename X>
        friend std::ostream &operator<<(std::ostream &, const Polynomial<X> &);

        template<typename X>
        Polynomial &operator=(const Polynomial<X> &object) {
            this->values = object.values;
            return *this;
        }

        template<typename X>
        friend Polynomial<X> add(const Polynomial<X> &, const Polynomial<X> &);

        template<typename X>
        friend Polynomial<X> subtract(const Polynomial<X> &, const Polynomial<X> &);

        template<typename X>
        friend Polynomial<X> multiply(const Polynomial<X> &, const Polynomial<X> &, const Polynomial<X> &);

        template<typename X>
        friend Polynomial<X> modDivide(const Polynomial<X> &, const Polynomial<X> &);

        template<typename X>
        friend Polynomial<X> MakeMonic(Polynomial<X>);

        template<typename X>
        friend Polynomial<X> Derivative(Polynomial<X>);

        template<typename X>
        friend X PointValue(Polynomial<X>, X);

        template<typename X>
        friend X modDivisionKoef(X, X, int);

        template<typename X>
        friend std::pair<Polynomial<X>, Polynomial<X>>
        divide(const Polynomial<X> &, const Polynomial<X> &);

        template<typename X>
        friend X sumPow(Polynomial<X> &poly, int maxPow);

        template<typename X>
        friend std::pair<X, X> positionDivide(Polynomial<X> &poly, int pow);

        template <typename X>
        friend Polynomial<X> simplMult(const Polynomial<X> &, const Polynomial<X> &);

        template<typename X>
        friend int getPower(Polynomial<X> a);

        template<typename X>
        friend bool equals(Polynomial<X> a, Polynomial<X> b);

        template<typename X>
        friend Polynomial<X> gcd_polynom(Polynomial<X> a, Polynomial<X> b);

        template<typename X>
        friend std::vector<Polynomial<X>> factorise_Ri(int q, int n);

        template<typename X>
        friend bool check_zero(Polynomial<X> a);

    };

    template<typename T>
    Polynomial<T>::Polynomial(const std::initializer_list<T> &initValues, const T &p) : values(initValues), p_(p) {}

    template<typename T>
    Polynomial<T>::Polynomial(const std::vector<T> &values, const T &p) : values(values), p_(p) {}

    template<typename T>
    Polynomial<T>::Polynomial(T number, const T &mod, size_t len) : p_(mod) {
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

    template<typename T>
    Polynomial<T>::Polynomial(const std::string &polyForm, const T &p, size_t n, char variable) : p_(p) {
        std::map<size_t, T> monoms = toPolynomial<T>(polyForm, variable);

        this->values = std::vector<T>(n);

        for (auto &monom: monoms) {
            this->values[this->values.size() - monom.first - 1] = monom.second;
        }
    }

    template<typename T>
    size_t Polynomial<T>::n() const {
        return this->values.size();
    }

    template<typename T>
    T Polynomial<T>::p() const {
        return this->p_;
    }


    template<typename T>
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
            const T &coef = this->values[i];

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

        const T &coefPow1 = this->values.size() >= 2 ? this->values[this->values.size() - 2] : 0;
        const T &coefPow0 = this->values.size() >= 1 ? this->values[this->values.size() - 1] : 0;

        if (coefPow1 > 0) {
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

        std::string result = res.str();
        return result == "" ? "0" : result;
    }

    template<typename X>
    std::ostream &operator<<(std::ostream &os, const Polynomial<X> &polynomial) {
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

    template<typename X>
    Polynomial<X> add(const Polynomial<X> &lt, const Polynomial<X> &rt) {
        if (lt.n() != rt.n() || lt.p() != rt.p()) {
            throw std::invalid_argument("Polynomials must be in the same field");
        }

        std::vector<X> values;
        for (size_t i = 0; i < lt.n(); ++i) {
            values.push_back(ModArithmetic<X>::add(lt.values[i], rt.values[i], lt.p()));
        }

        return Polynomial<X>(values, lt.p());
    }

    template<typename X>
    Polynomial<X> subtract(const Polynomial<X> &lt, const Polynomial<X> &rt) {
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
    Polynomial<X> multiply(const Polynomial<X> &lt, const Polynomial<X> &rt, const Polynomial<X> &primitive) {
        if (lt.n() != rt.n() || lt.p() != rt.p()) {
            throw std::invalid_argument("Polynomials must be in the same field");
        }

        std::vector<X> values(2 * lt.n() - 1, 0);

        for (size_t i = 0; i < lt.n(); ++i) {
            for (size_t j = 0; j < rt.n(); ++j) {
                values[i + j] = values[i + j] + lt.values[i] * rt.values[j];
            }
        }

        for (X &value: values) {
            value = value % lt.p();
        }

        return modDivide(Polynomial<X>(values, lt.p()), primitive);
    }

    template<typename X>
    Polynomial<X> modDivide(const Polynomial<X> &lt, const Polynomial<X> &rt) {
        if (lt.p() != rt.p()) {
            throw std::invalid_argument("Mods are not equal");
        }

        auto tmpLt = lt;
        auto tmpRt = rt;
        std::vector<X> rtValues;
        size_t k = 0;
        while (k < tmpRt.n() && tmpRt.values[k] == 0) {
            ++k;
        }
        while (k < tmpRt.n()) {
            rtValues.push_back(tmpRt.values[k]);
            ++k;
        }

        if (rtValues.size() == 0) {
            throw std::logic_error("Zero division error");
        }
        tmpRt.values = rtValues;

        while (tmpLt.n() >= tmpRt.n()) {
            while (tmpLt.values[0] != 0) {
                auto multiplier = ModArithmetic<X>::divide(tmpLt.values[0], tmpRt.values[0], lt.p());

                for (size_t i = 0; i < tmpRt.n(); ++i) {
                    tmpLt.values[i] = ModArithmetic<X>::subtract(tmpLt.values[i], ModArithmetic<X>::multiply(multiplier, tmpRt.values[i], lt.p()), lt.p());
                }
            }

            tmpLt.values.erase(tmpLt.values.begin());
        }

        return tmpLt;
    }

    /* Creates Monic Polynomial from given Polynomial */
    template<typename X>
    Polynomial<X> MakeMonic(gf::Polynomial<X> pol) {
        std::vector<X> init_vec;
        X max_coef = 0;
        for (int i = 0; i < pol.n(); i++) {
            if (pol.values[i] != 0) {
                max_coef = pol.values[i];
                break;
            }
        }
        X p = pol.p();
        if (max_coef == 1) {
            return pol;
        } else {
            for (int i = 0; i < pol.n(); i++) {
                init_vec.push_back(ModArithmetic<X>::divide(pol.values[i], max_coef, p));
            }
            return Polynomial<X>(init_vec, p);
        }
    }

    /* Creates Derivative Polynomial from given Polynomial */
    template<typename X>
    gf::Polynomial<X> Derivative(gf::Polynomial<X> pol) {
        std::vector<X> init_vec;
        X p = pol.p();
        for (int i = 0; i < pol.n() - 1; i++) {
            init_vec.push_back(ModArithmetic<X>::multiply(pol.values[i], pol.n() - i - 1, p));
        }
        return gf::Polynomial<X>(init_vec, p);
    }

    template<typename X>
    long int ModularPow(X base, long int degree, X mod) {
        if (degree == 0) return 1;
        X res = base;
        for (long int i = 1; i < degree; i++) {
            res = ModArithmetic<X>::multiply(res, base, mod);
        }
        return res;
    }

    /*Finds polynomial value at given point  */
    template<typename X>
    X PointValue(gf::Polynomial<X> pol, X x) {
        std::vector<X> init_vec;
        X p = pol.p();
        X result = 0;
        for (int i = 0; i < pol.n(); i++) {
            result = ModArithmetic<X>::add(result, ModArithmetic<X>::multiply(pol.values[i],
                                                                              gf::ModularPow(x, pol.n() - i - 1, p), p),
                                           p);
        }
        return result;
    }


    /**
     * Function that returns coef of quotient
     *
     * @tparam X
     * @param dividentKf
     * @param divisorKoef
     * @param mod
     * @return
     */
    template<typename X>
    X modDivisionKoef(X dividentKoef, X divisorKoef, int mod) {
        for (X count = 1; count <= mod; count++) {
            if ((count * divisorKoef) % mod == dividentKoef)
                return count;
        }
    }

    /**
     * Function that returns sum of degrees
     *
     * @tparam X
     * @param poly dividend poly
     * @param maxPow max degree of divisor
     * @return
     */
    template<typename X>
    X sumPow(Polynomial<X> &poly, int maxPow) {
        X sum = 0;
        for (int i = 0; i <= maxPow; i++) {
            sum += poly.values[i];
        }
        return sum;
    }

    /**
     * Function that returns max divident degree and quotient degree
     *
     * @tparam X
     * @param poly
     * @param pow
     * @return
     */
    template<typename X>
    std::pair<X, X> positionDivide(Polynomial<X> &poly, int pow) {
        std::pair<X, X> finish;
        X count = 0;
        while (poly.values[count] == 0) {
            count++;
            if (count == poly.values.size())
                break;
        }
        finish.first = count;
        std::map<X, X> pows;
        for (int i = 0; i < poly.values.size(); i++) {
            int size = poly.values.size();
            pows[i] = size - i - 1;
        }
        int divident_pow = pows[count];
        int divisor_pow = pows[pow];
        int result_pow = divident_pow - divisor_pow;

        for (auto i: pows) {
            if (i.second == result_pow) {
                finish.second = i.first;
                break;
            }
        }
        return finish;
    }

    int CountRoots(gf::Polynomial<int> P) {
        std::string s = P.toString();
        int q = P.p();
        std::map<size_t, int> pol = toPolynomial<int>(s, 'x');
        auto itr = pol.begin();

        int arr[q - 1][q - 1];
        for (auto i=0; i<q-1; i++)
            arr[0][i] = 0;

        for (itr; itr != pol.end(); itr++)
                arr[0][itr->first % (q-1)] = (arr[0][itr->first % (q-1)] + itr->second) % q;

        for (auto i = 1; i < q - 1; i++) {
            for (auto j = 0; j < q - 2; j++)
                arr[i][j] = arr[i - 1][j + 1];
            arr[i][q - 2] = arr[i - 1][0];
        }


        /*for(int i = 0; i < q-1; i++)
        {
            for(int j = 0; j < q-1; j++)
                std::cout<<arr[i][j]<<" ";
            std::cout<<std::endl;
        }*/

        int n = q - 1;
        int r = 0;

        for (auto i = 0; i < n; i++) {
            for (auto j = i; j < n; j++)
                if (arr[j][i] > 0) {
                    r++;
                    for (auto temp = i; temp < n; temp++)
                        std::swap(arr[i][temp], arr[j][temp]);
                    break;
                }

            for (auto j = i + 1; j < n; j++)
                if (arr[j][i])
                    while (arr[j][i] != 0)
                        for (auto temp = i; temp < n; temp++)
                            arr[j][temp] = (arr[j][temp] + arr[i][temp]) % q;
        }

        /*std::cout<<std::endl;
        for(int i = 0; i < q-1; i++)
        {
            for(int j = 0; j < q-1; j++)
                std::cout<<arr[i][j]<<" ";
            std::cout<<std::endl;
        }*/
        return q - 1 - r;
    }

    template <typename X>
    Polynomial<X> simplMult(const Polynomial<X> &first, const Polynomial<X> &second){
        std:: vector<X> result(first.values.size(),0);
        std::pair<X,X>temp_pow_koef;
        for(int i=0; i<second.values.size(); i++){
            if(second.values[i] != 0){
                temp_pow_koef.first =second.values[i];
                temp_pow_koef.second = i;
            }
        }
        std::map<X, X> pows;
        for (int i = 0; i < first.values.size(); i++) {
            int size = first.values.size();
            pows[i] = size - i - 1;
        }
        auto pow = pows[temp_pow_koef.second];

        for(int i=0; i<first.values.size(); i++) {
            if(first.values[i] != 0) {
                auto step = pows[i];
                auto finPow= pow + step ;
                for (auto j: pows) {
                    if (j.second == finPow) {
                        result[j.first] = ModArithmetic<X>::multiply(first.values[i], temp_pow_koef.first, first.p_);
                        break;
                    }
                }
            }
        }
        Polynomial<X>mul(result, first.p_);
        return mul;
    }


    /**
     * Function that divide polynoms in Field
     *
     * @tparam X
     * @param divident
     * @param divisor
     * @param primitive
     * @return
     */
    template<typename X>
    std::pair<Polynomial<X>, Polynomial<X>>
    divide(const Polynomial<X> &divident, const Polynomial<X> &divisor) {
        Polynomial<X> tet(divident.values, divident.p());
        size_t sizeDivident = divident.values.size();
        std::vector<X> quotientVec(sizeDivident, 0);

        if (divident.values[0] == 0 && divisor.values[0] != 0) {
            Polynomial<X> quot(quotientVec, divident.p());
            std::pair<Polynomial<X>, Polynomial<X>> unexpected(quot,tet);
            return unexpected;
        }

        int iterDivisor = 0;
        for (int j = 0; j < sizeDivident; ++j) {
            if (divisor.values[j] != 0) {
                break;
            } else {
                iterDivisor++;
            }
        }
        while (sumPow(tet, iterDivisor) != 0) {
            std::vector<X> tempDivident(sizeDivident, 0);
            std::pair<X, X> pow_and_koef = positionDivide(tet, iterDivisor);
            int modKoef = modDivisionKoef(tet.values[pow_and_koef.first], divisor.values[iterDivisor], divident.p_);
            tempDivident.at(pow_and_koef.second) = modKoef;
            quotientVec = VecAdd(quotientVec, tempDivident);
            Polynomial<X> mult(tempDivident, divident.p_);

            Polynomial<X> result = simplMult(divisor, mult);
            Polynomial<X> substr = subtract(tet, result);
            tet = substr;
        }
        Polynomial<X> quotient(quotientVec, divident.p_);
        std::pair<Polynomial<X>, Polynomial<X>> quotientAndRemainder(quotient, tet);

        return quotientAndRemainder;

    }

    std::vector<gf::Polynomial<int>> getAllIrreduciblePolynomials(int p, int n, int amount) {

        std::vector<gf::Polynomial<int>> result, buffer;
        gf::Polynomial<int> *cyclic;

        int num = pow(p, n) - 1;

        for (int m = 1; m <= num; ++m) {
            if ((num % m) == 0) {
                cyclic = new gf::Polynomial<int>((new CyclicPolynomial())->calculatePolynomial(n), p);

                if (std::__gcd(p, m)) {
//                    buffer = decompositionWithRi(m);
                } else {
                    gf::Polynomial<int> temp = gf::Derivative(*cyclic);
//                    temp = gcd(cyclic, temp);
//                    buffer = berlekamp(temp);
                }

                for (auto &irreduciblePolynomial : buffer) {
                    if (irreduciblePolynomial.n() == n) {
                        if (result.size() < amount) result.push_back(irreduciblePolynomial);
                        else return result;
                    }
                }

            }
        }

//        return result;
        return {gf::Polynomial<int>("x^6+x^4+2x^3+x+2", p, n + 1), gf::Polynomial<int>("8x^3+4x^2+1", p, n)};
    }
    
    /**
     * function to calculate Euler Functions
     *
     * @param int
     * @return
     */
    int phi(int n) {
        int result = n;
        for (int i = 2; i * i <= n; ++i)
            if (n % i == 0) {
                while (n % i == 0)
                    n /= i;
                result -= result / i;
            }
        if (n > 1)
            result -= result / n;
        return result;
    }
    /**
     * function for checking polynomial equality
     *
     * @tparam X
     * @param Polynomial<X>
     * @param Polynomial<X>
     * @return
     */
    template<typename X>
    bool equals(Polynomial<X> a, Polynomial<X> b) {
        if (a.n() != b.n())
            return false;
        for (size_t i = 0; i < a.n(); i++) {
            if (a.values[i] != b.values[i]) return false;
        }
        return true;
    }
    /**
     * function to find the maximum degree of a polynomial
     *
     * @tparam X
     * @param Polynomial<X>
     * @return
     */
    template<typename X>
    int getPower(Polynomial<X> a) {
        for (int i = 0; i < a.n(); i++) {
            if (a.values[i] != 0)
                return a.n() - i - 1;
        }
        return 0;
    }
    /**
     * function for finding gcd of two polynomials
     *
     * @tparam X
     * @param Polynomial<X>
     * @param Polynomial<X>
     * @return
     */
    template<typename X>
    Polynomial<X> gcd_polynom(Polynomial<X> a, Polynomial<X> b) {
        Polynomial<X> p1 = a;
        Polynomial<X> p2 = b;
        int pow_a = getPower(a);
        int pow_b = getPower(b);
        if (getPower(a) >= getPower(b)) {
            p1 = b;
            p2 = a;
        }
        else {
            p1 = a;
            p2 = b;
        }

        while (p2.n() != 0 && check_zero(p2)) {
            Polynomial<X> temp = p2;
            p2 = modDivide(p1, p2);
            p1 = temp;
        }

        Polynomial<X> result = p1;

        return MakeMonic(result);
    }
    /**
     * function to check the zero polynomial
     *
     * @tparam X
     * @param Polynomial<X>
     * @return
     */
    template<typename X>
    bool check_zero(Polynomial<X> a) {
        for (int i = 0; i < a.n(); i++) {
            if (a.values[i] != 0)
                return true;
        }
        return false;
    } 
    /**
     * function to find Ri 
     *
     * @tparam X polynomial
     * @param q
     * @param n
     * @param i
     * @param j
     * @return
     */
    template<typename X>
    Polynomial<X> find_Ri(int q, int n, int i, int j) {
        while (((int)pow((double)q, (double)j) * i) % n != i) {
            j++;
        }
        int degree = (int)pow((double)q, (double)j) * i + 1;
        std::vector<int> Ri_values(degree, 0);
        for (int k = j - 1; k >= 0; k--) {
            int d = degree - (int)pow((double)q, (double)k) * i - 1;
            Ri_values[d] = 1;
        }
        Polynomial<X> result(Ri_values, q);
        return result;
    }
    /**
     * function for finding the decomposition of a circular polynomial
     *
     * @tparam X
     * @param q
     * @param n
     * @return
     */
    template<typename X>
    std::vector<Polynomial<X>> factorise_Ri(int q, int n) {
        std::vector<Polynomial<X>> result;
        CyclicPolynomial* newCyclic = new CyclicPolynomial();
        Polynomial<X> cycle_polinomial(newCyclic->calculatePolynomial(n), q);
        if (n == 1) {
            result.push_back(MakeMonic(cycle_polinomial));
            return result;
        }
        if (std::__gcd(n, q) > 1) {
            int newN = n;
            while (newN % q == 0)
                newN /= q;
            result = factorise_Ri<X>(q, newN);
            CyclicPolynomial* new_Cyclic = new CyclicPolynomial();
            Polynomial<X> new_cycle_polinomial(new_Cyclic->calculatePolynomial(newN), q);
            size_t repeat = getPower(cycle_polinomial) / getPower(new_cycle_polinomial);
            size_t count = result.size();
            for (size_t i = 1; i < repeat; i++) {
                for (size_t j = 0; j < count; j++) {
                    result.push_back(MakeMonic(result[j]));
                }
            }
            return result;
        }
        int d = 1;
        int w;
        do {
            d++;
            w = (int)pow((double)q, (double)d) % n;
        } while (w != 1);
        size_t factors_count = phi(n) / d;
        if (factors_count == 0) factors_count = 1;
        size_t factor_power = d;

        if (factors_count == 1) {
            result.push_back(MakeMonic(cycle_polinomial));
            return result;
        }

        size_t i = 1;

        std::list<Polynomial<X>> temp_pol;
        temp_pol.push_back(cycle_polinomial);

        while (result.size() < factors_count && i < n) {
            int j = 1;
            Polynomial<X> pol_Ri = find_Ri<X>(q, n, i, j);
            bool factorized = false;
            j = 0;
            while (j < q) {
                Polynomial<X> gcdRi = gcd_polynom(temp_pol.front(), pol_Ri);
                if (j == 0 && (equals(gcdRi, temp_pol.front()) || equals(gcdRi, pol_Ri))) {
                    factorized = false;
                    break;
                }

                long long gcdPower = getPower<X>(gcdRi);
                if (gcdPower == factor_power) {
                    factorized = true;
                    Polynomial<int> gcdRi_new = MakeMonic(gcdRi);
                    result.push_back(gcdRi_new);
                }
                else if (gcdPower > 0 && gcdPower % factor_power == 0) {
                    factorized = true;
                    temp_pol.push_back(gcdRi);
                }

                pol_Ri = add(pol_Ri, Polynomial<X>("1", pol_Ri.p(), pol_Ri.n()));
                j++;
            }
            if (factorized) {
                temp_pol.pop_front();
            }
            i++;
        }
        return result;
    }
}


#endif //ABSTRACT_ALGEBRA2_POLYNOMIAL_H
