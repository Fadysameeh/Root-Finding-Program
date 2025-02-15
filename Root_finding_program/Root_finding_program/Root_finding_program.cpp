#include "stdlib.h"
#include <iostream>
#include <string>
#include <cstdlib>
#include <cstring>
#include <complex>
using namespace std;
using namespace std::complex_literals;
float const EPSILON = 0.000000001;
float* coefficients;
// counting terms
int countOperator(string str)
{
	// the count value increases by 1 whenever a sign is read by the for loop
	int count = 1;
	for (int i = 0; i < str.size(); i++)
	{
		if (str[i] == '+' || str[i] == '-')
			count++;
	}
	return count;
}
//ERRORS
bool isNumPowDigit(string str)
{
	int i = 0;
	if (str[i] == '+' || str[i] == '-') i++;
	while ((str[i] >= '0' && str[i] <= '9') || str[i] == '.') {
		i++;
	}
	if (!(str[i] == 'x' || str[i] == 'X'))
		return false;
	i++;
	if (!(str[i] == '^')) {
		return false;
	}
	i++;
	if (!(str[i] >= '0' && str[i] <= '5')) {
		return false;
	}
	i++;
	if (!(str[i] == '\0'))
		return false;
	return true;
}
//extract exponents, get degree of poly
int getdegree(string str)
{
	// initializing array to save all exponent values in equation
	int exp[10] = { 0 };
	// for loop to read all values in string
	for (int i = 0, j = 0; i <= str.size(); i++)
	{
		// also checking en el value da raqam msh sign
		if ((str[i] == 'x' || str[i] == 'X') && (str[i + 1] == '^') && (str[i] != '+' || str[i] != '-'))
		{
			exp[j] = str[i + 2] - '0';
			j++;
		}
		else if ((str[i] == 'x' || str[i] == 'X') && (str[i] != '+' || str[i] != '-'))
		{
			exp[j] = 1;
			j++;
		}
	}
	int max = 1;
	for (int i = 0; i < 10; i++)
	{
		if (exp[i] >= max)
			max = exp[i];
	}
	return max;
}
string* extractTerms(string poly, int termsCount) {
	// extracting the polynomial's terms
	int k = 0;
	int j = 0;
	// create and initialize a dynamic array to save terms
	string* terms = new string[100];
	for (int i = 1; i < poly.size(); i++)
		if (poly[i] == '+' || poly[i] == '-')
		{
			terms[k++] = poly.substr(j, i - j);
			j = i;
		}
	if (j != poly.size()) {
		terms[k] = poly.substr(j, poly.size() - j);
	}
	for (int i = 0; i < termsCount; i++)
	{
		// if x is with coeff 1
		if (terms[i][0] == 'x' && i == 0) {
			terms[i].insert(0, "1");
		}
		if (terms[i][1] == 'x' && i != 0) {
			terms[i].insert(1, "1");
		}
		bool xPresent = false;
		for (int j = 1; j < terms[i].size(); j++)
			// if x is with exponent 1
			if (terms[i][j] == 'x')
			{
				xPresent = true;
				if (j == terms[i].size() - 1)
					terms[i].append("^1");
				break;
			}
		// in case of constant
		if (!xPresent)
			terms[i].append("x^0");
	}
	return terms;
}
bool checkError(string* terms, int termsCount) {
	for (int i = 0; i < termsCount; i++)
	{
		if (!isNumPowDigit(terms[i]))
		{
			return false;
		}
	}
	return true;
}
// extract coeff
float* getcoeff(string* terms, int termsCount)
{
	// create , initialize array for exp
	int expon[100] = { 0 };
	// create, initialize dynamic array for coeff with 0 3shan coeff le exp 0 yekon b 0
	float* coeff = new float[100] { 0 };
	// for loop le terms array el kbeera
	for (int i = 0; i < termsCount; ++i)
	{
		// for loop lel array el soghyra bta3t kol term
		for (int j = 0; j < terms[i].size(); ++j)
		{
			if (terms[i][j] == '^')
			{
				expon[i] = stoi(terms[i].substr(j + 1, terms[i].size() - j));
				coeff[expon[i]] = stof(terms[i].substr(0, j));
				break;
			}
		}
	}
	return coeff;
}
// calculate roots of poly
void Linear(float* values)
{
	float res = (-values[0]) / values[1];
	cout << "Roots = " << res << endl;
}
void SecondDegree(float* values)
{
	float x1, x2, discr, real, imaginary;
	discr = pow(values[1], 2) - 4 * values[2] * values[0];
	if (discr > 0) {
		cout << "Roots are real." << endl;
		x1 = (-values[1] + sqrt(discr)) / (2 * values[2]);
		x2 = (-values[1] - sqrt(discr)) / (2 * values[2]);
		cout << "x1 = " << x1 << " ,x2 = " << x2 << endl;
	}
	else if (discr == 0) {
		cout << "Roots are real." << endl;
		x1 = -values[1] / (2 * values[2]);
		cout << "x1 = x2 =" << x1 << endl;
	}
	else {
		real = -values[1] / (2 * values[2]);
		imaginary = sqrt(-discr) / (2 * values[2]);
		cout << "Roots are complex." << endl;
		cout << "x1 = " << real << "+" << imaginary << "i" << endl;
		cout << "x2 = " << real << "-" << imaginary << "i" << endl;
	}
}
void ThirdDegree(float* values)
{
	const double PI = 4.0 * atan(1.0);
	// Reduced equation: X^3 - 3zX - 2y = 0, where X = x-b/(3a)
	double z = (pow(values[2], 2) - 3.0 * values[3] * values[1]) / (9.0 * pow(values[3], 2));
	double y = (9.0 * values[3] * values[2] * values[1] - 27.0 * pow(values[3], 2) * values[0] - 2.0
		* pow(values[2], 3))
		/ (54.0 * pow(values[3], 3));
	double offset = values[2] / (3.0 * values[3]);
	// Discr
	double discr = pow(z, 3) - pow(y, 2);
	cout << "Roots: " << endl;
	if (discr > 0)
		// set X = 2 sqrt(z) cos(theta) and compare 4 cos^3(theta)-3 cos(theta) = cos(3 theta)
	{
		double theta = acos(y / (z * sqrt(z)));
		double r = 2.0 * sqrt(z);
		for (int n = 0; n < 3; n++)
		{
			cout << r * cos((theta + 2.0 * n * PI) / 3.0) - offset << endl;
		}
	}
	else
	{
		double J = cbrt(y + sqrt(-discr));
		double K = cbrt(y - sqrt(-discr));
		cout << J + K - offset << '\n';
		double Real = -0.5 * (J + K) - offset;
		double Imaginary = (J - K) * sqrt(3.0) / 2.0;
		if (discr == 0.0)
			// Equal roots
		{
			cout << Real << endl;
			cout << Real << endl;
		}
		else
		{
			cout << Real << " + " << Imaginary << " i" << endl;
			cout << Real << " - " << Imaginary << " i" << endl;
		}
	}
}
std::complex<double> cuberoot(std::complex<double> z) {
	if (z.real() < 0) {
		return -pow(-z, 1.0 / 3.0);
	}
	else {
		return pow(z, 1.0 / 3.0);
	}
}
void FourthDegree(float* values) {
	double a = values[4];
	double b = values[3];
	double c = values[2];
	double d = values[1];
	double e = values[0];
	complex<double> p1 = 2 * c * c * c - 9 * b * c * d + 27 * a * d * d + 27 * b * b * e - 72
		* a * c * e;
	complex<double> p2 = p1 + sqrt(-4 * pow((c * c - 3 * b * d + 12 * a * e), 3) + pow(p1, 2));
	complex<double> p3 = (c * c - 3 * b * d + 12 * a * e) / (3 * a * cuberoot(p2 / 2.0)) +
		(cuberoot(p2 / 2.0) / (3.0 * a));
	complex<double> p4 = sqrt(b * b / (4 * a * a) - 2 * c / (3 * a) + p3);
	complex<double> p5 = b * b / (2 * a * a) - 4 * c / (3 * a) - p3;
	complex<double> p6 = (-b * b * b / (a * a * a) + 4 * b * c / (a * a) - 8 * d / a) / (4.0 *
		p4);
	complex<double> x1 = -b / (4 * a) - p4 / 2.0 - sqrt(p5 - p6) / 2.0;
	complex<double> x2 = -b / (4 * a) - p4 / 2.0 + sqrt(p5 - p6) / 2.0;
	complex<double> x3 = -b / (4 * a) + p4 / 2.0 - sqrt(p5 + p6) / 2.0;
	complex<double> x4 = -b / (4 * a) + p4 / 2.0 + sqrt(p5 + p6) / 2.0;
	if (x1.imag() == 0)
		cout << "First root: " << x1.real() << endl;
	else cout << "First root: " << x1.real() << (x1.imag() > 0 ? "+" : "") << x1.imag() << "i"
		<< endl;
	if (x2.imag() == 0)
		cout << "Second root: " << x2.real() << endl;
	else cout << "Second root: " << x2.real() << (x2.imag() > 0 ? "+" : "") << x2.imag() <<
		"i" << endl;
	if (x3.imag() == 0)
		cout << "Third root: " << x3.real() << endl;
	else cout << "Third root: " << x3.real() << (x3.imag() > 0 ? "+" : "") << x3.imag() << "i"
		<< endl;
	if (x4.imag() == 0)
		cout << "Fourth root: " << x4.real() << endl;
	else cout << "Fourth root: " << x4.real() << (x4.imag() > 0 ? "+" : "") << x4.imag() <<
		"i" << endl;
}
double func(double x) {
	double res = 0;
	for (int i = 0; i < 6; i++)
		res += coefficients[i] * pow(x, 5 - i);
	return res;
}
double derivFunc(double x) {
	double res = 0;
	for (int i = 0; i < 5; i++)
		res += coefficients[i] * (5 - i) * pow(x, 5 - i - 1);
	return res;
}
double newtonRaphson(double x)
{
	double h = func(x) / derivFunc(x);
	while (abs(h) >= EPSILON)
	{
		h = func(x) / derivFunc(x);
		x = x - h;
	}
	return x;
}
void FifthDegree(float* values)
{
	double* N = new double[6];
	N[5] = coefficients[0];
	N[4] = coefficients[1];
	N[3] = coefficients[2];
	N[2] = coefficients[3];
	N[1] = coefficients[4];
	N[0] = coefficients[5];
	double* D = new double[2];
	double* q;
	double* r;
	double* d;
	double root1 = newtonRaphson(0);
	cout << "root" << root1 << endl;
	D[0] = 1;
	D[1] = -root1;
	int dN = 5, dD = 1, dd = 0, dq = 4, dr = 4;
	d = new double[dN + 1];
	q = new double[dq + 1];
	r = new double[dr + 1];
	while (dN >= dD) {
		// d equals D shifted right
		d = new double[2];
		fill(d, d + 2, 0);
		for (int i = 0; i <= dD; i++)
			d[i + dN - dD] = D[i];
		dd = dN;
		// calculating one element of q
		q[dN - dD] = N[dN] / d[dd];
		// d equals d * q[dN-dD]
		for (int i = 0; i < dq + 1; i++)
			d[i] = d[i] * q[dN - dD];
		// N equals N - d
		for (int i = 0; i < dN + 1; i++)
			N[i] = N[i] - d[i];
		dN--;
	}
	coefficients[0] = q[4];
	coefficients[1] = q[3];
	coefficients[2] = q[2];
	coefficients[3] = q[1];
	coefficients[4] = q[0];
}
int main()
{
	// create and initialize char for poly
	string poly;
	cout << "enter equation:";
	getline(cin, poly);
	// calculate number of terms in given equation
	int termsCount = countOperator(poly);
	cout << termsCount << " terms" << endl;
	// Extracting the terms and reformatting how it's written
	string* terms = extractTerms(poly, termsCount);
	// Checking for errors
	if (!checkError(terms, termsCount)) {
		cout << "Error, please re-try" << endl;
		// deallocate the dynamic array
		delete[] terms;
		system("pause");
		return 0;
	}
	// Calculating the degree of equation
	int degree = getdegree(poly);
	cout << "degree of poly equation is: " << degree << endl;
	// Creating and initializing the coefficients array and calling the function that
	// extracts the coeff values from the equation
	coefficients = getcoeff(terms, termsCount);
	// for loop for printing the coeff values
	cout << "coefficients are from exp 0 till max: " << endl;
	for (int i = 0; i <= 5; i++) {
		cout << coefficients[i] << endl;
	}
	// Calling the root functions depending on degree value
	switch (degree)
	{
	case 1:
		Linear(coefficients);
		break;
	case 2:
		SecondDegree(coefficients);
		break;
	case 3:
		ThirdDegree(coefficients);
		break;
	case 4:
		FourthDegree(coefficients);
		break;
	case 5:
		FifthDegree(coefficients);
		break;
	}
	// deallocating dynamic array
	delete[] terms;
	system("pause");
	return 0;
}