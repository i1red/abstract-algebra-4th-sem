
#include <iostream>
#include <map>
#include <math.h>

//все функции
void factorize(int);
int pol_degree(int*, int);
int* polynomialDivision(int*, int*, int);
bool isDivisible(int*, int*, int);
int inverseElement(int);
void shiftPolyLeft(int* tmp_residue, int* divisor, int tmp_residue_degree, int divisor_degree);
bool isIrredusible();
bool polynomesAreEqual(int* pol1, int* pol2);

int q, m;//q - характеристика поля, m - степень расширения поля(степень полинома+1)
int* polynomial;//полином
std::map<int, int> factors;//разложение на множители q^m-1

//считаем степень полинома 
int pol_degree(int* poly, int size = m) {
	for (int i = size - 1; i >= 0; i--) {
		if (poly[i] != 0)
			return i;
	}

	return -1;
}

int* polynomialDivision(int* divident, int* divisor, int divident_size = m) {
	int* tmp_residue = new int[divident_size];

	int divisor_degree = pol_degree(divisor);

	//создаем временную переменную для делителя
	int* tmp_divisor = new int[divident_size];

	memcpy(tmp_residue, divident, divident_size * sizeof(int));

	do {
		int tmp_residue_degree = pol_degree(tmp_residue, divident_size);

		memset(tmp_divisor, 0, divident_size * sizeof(int));
		memcpy(tmp_divisor, divisor, m * sizeof(int)); //Обновляем содержимое временного делителя

		//сдвиг влево на n = умножение полинома на х^n
		shiftPolyLeft(tmp_residue, tmp_divisor, tmp_residue_degree, divisor_degree);

		for (int i = 0; i < divident_size; i++) {
			tmp_residue[i] = (tmp_residue[i] - tmp_divisor[i]) % q;
			if (tmp_residue[i] < 0)
				tmp_residue[i] += q;
		}

	} while (pol_degree(tmp_residue, divident_size) >= m - 1);

	return tmp_residue;
}

// x^(-1)(mod q)
int inverseElement(int x) {
	for (int i = 0; i < q; i++) {
		if ((x * i) % q == 1)
			return i;
	}
	return 0;
}

//Сдвинуть делитель влево(умножить на x^(...)) и если надо, то умножить на коэффициент (changes content of tmp_divisor)
void shiftPolyLeft(int* tmp_residue, int* tmp_divisor, int tmp_residue_degree, int divisor_degree) {
	int shift = tmp_residue_degree - divisor_degree;
	int div_multiplicator = 1;

	if (tmp_residue[tmp_residue_degree] != tmp_divisor[divisor_degree]) {
		div_multiplicator = inverseElement(tmp_divisor[divisor_degree]) * tmp_residue[tmp_residue_degree];
	}
	for (int i = divisor_degree + shift; i >= shift; i--) {
		tmp_divisor[i] = (div_multiplicator * tmp_divisor[i - shift]) % q;
	}
	for (int i = shift - 1; i >= 0; i--) {
		tmp_divisor[i] = 0;
	}

}

//Проверяем есть ли у полинома делители
bool isDivisible(int* divident, int* divisor, int divident_size = m) {

	//если степень делителя больше чем делимого, то не делится
	if (pol_degree(divident, divident_size) < pol_degree(divisor)) {
		return false;
	}

	int* result = polynomialDivision(divident, divisor, divident_size);
	if (pol_degree(result, divident_size) == -1)// если остаток = 0
	{

		return true;
	}

	return false;
}

bool isIrredusible() {
	int* polynomial_divisor = new int[m];
	for (int i = 0; ++i < pow<int, int>(q, m);) {
		for (int j = 0; j < m; j++) {
			polynomial_divisor[j] = (int)(i / pow<int, int>(q, j)) % q;
		}
		if (pol_degree(polynomial_divisor) > 0 && !polynomesAreEqual(polynomial, polynomial_divisor)) {//проверяем делиться ли он на многочлены степени > 1
			if (isDivisible(polynomial, polynomial_divisor))
				return false;
		}
	}
	return true;
}

bool polynomesAreEqual(int* pol1, int* pol2) {
	for (int j = 1; j < q; j++) {
		bool areEqual = true;
		int i = 0;
		for (; i < m; i++) {
			if (pol1[i] != ((j * pol2[i]) % q)) {
				areEqual = false;
				break;
			}
		}
		if (i == m && areEqual == true)
			return true;
	}
	return false;
}

//passing factors map reference
void factorize(int number) {
	for (int i = 2; i <= (int)(sqrt(number)); i++) {
		if (number % i == 0) {
			number /= i;
			factors.insert(std::pair<int, int>(i, 1));
			while (number % i == 0) {
				number /= i;
				factors.at(i) += 1;
			}
		}
	}
	if (number != 1) {
		factors.insert(std::pair<int, int>(number, 1));
	}
}

//Непосредственно находим порядок по алгоритму
int polynomialOrder() {
	int e = 1;
	int common_power = pow<int, int>(q, m) - 1;
	for (std::map<int, int>::iterator itr = factors.begin(); itr != factors.end(); itr++) {
		int tmp_power = common_power / itr->first;

		//устанавливаем коєффициенты полинома x^((q^m-1)/pj)-1
		int* tmp_poly_divident = new int[tmp_power + 1];
		memset(tmp_poly_divident, 0, (tmp_power + 1) * sizeof(int));
		tmp_poly_divident[tmp_power] = 1;
		tmp_poly_divident[0] = q - 1;

		if (isDivisible(tmp_poly_divident, polynomial, tmp_power + 1 >= m ? tmp_power + 1 : m)) {

			printf("x^((q^m-1)/%d) - 1 делится на f(x) -> e не делиться на %d^%d\nпроверяем младшие степени:\n", itr->first, itr->first, itr->second);

			//проверяем делиться ли E на степени pj^(rj-1),..., pj
			for (int i = itr->second - 1; i > 0; i--) {

				//выделяем память для временного полинома x^( (q^m-1) / (p^k) ), k = 2,3,4...,rj 
				int tmp_power2 = common_power / (int)pow<int, int>(itr->first, itr->second - i + 1);
				int* tmp_poly_divident2 = new int[tmp_power2 + 1];
				memset(tmp_poly_divident2, 0, (tmp_power2 + 1) * sizeof(int));
				tmp_poly_divident2[tmp_power2] = 1;
				tmp_poly_divident2[0] = q - 1;

				if (isDivisible(tmp_poly_divident2, polynomial, tmp_power2 + 1 >= m ? tmp_power2 + 1 : m)) {
					printf("е не делиться на %d^%d\n", itr->first, i);
				}
				else {
					e *= (int)pow(itr->first, i);
					printf("е делиться на %d^%d\n", itr->first, i);
					break;
				}
			}
		}
		else {
			printf("x^%d-1 не делиться на f(x) -> e делиться на %d^%d\n", itr->first, itr->first, itr->second);
			e *= (int)pow<int, int>(itr->first, itr->second);
		}
	}

	return e;
}

int main()
{
	setlocale(LC_ALL, "Russian");

	std::cout << "Введите характеристику поля(простое число)\n\t";
	std::cin >> q;

	std::cout << "Введите степень полинома\n\t";
	std::cin >> m;
	m++;

	polynomial = new int[m];

	std::cout << "Введите коэффициенты полинома\n";
	for (int i = m - 1; i >= 0; i--) {
		std::cout << "\t_*x^" << i << " : ";
		std::cin >> polynomial[i];
		std::cout << "\n";
	}

	if (!isIrredusible()) {
		std::cout << "Введенный полином приводим завершение программы\n\n";
		return 0;
	}
	else {
		std::cout << "Введенный полином неприводим\n\n";
	}


	//разкладываем на множители и записываем результаты в map
	factorize((int)pow<int, int>(q, m - 1) - 1);

	std::cout << "Делители q^m-1 и их степени:\n";
	for (std::map<int, int>::iterator itr = factors.begin(); itr != factors.end(); itr++) {
		printf("\t%d -> %d\n", itr->first, itr->second);
	}
	printf("\n");

	int poly_order = polynomialOrder();

	std::cout << "\nпорядок введенного неприводимого полинома = " << poly_order << "\n";


}
