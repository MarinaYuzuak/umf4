#ifndef _SOLVER_HPP
#define _SOLVER_HPP

#include <filesystem>
#include <ctime>

#include "lightweight.hpp"
#include "Vector.hpp"
using std::filesystem::path;

class Solver {
public:
	struct Matrix {
		Vector di;
		Vector ggl;
		Vector ggu;
		std::vector<size_t> ig;
		std::vector<size_t> jg;
		Vector operator*(const Vector& v) const;
	};

	Solver(const path& path) {
		readParams(path);
		clock_t start = clock();

		//size_t iter = losLU();
		size_t iter = BСGSTABLU();

		clock_t end = clock();
		double seconds = double(end - start) / CLOCKS_PER_SEC;
		std::cout << "Number of iterations: " << iter << "\nTime: " << seconds << " c.\n" << std::endl;
		printVector(x);
	};

private:
	void readParams(const path& path);
	void readSparseMatrix(const path& path);       // чтение матрицы, если она в разреженном формате.
	void readDenseMatrix(std::ifstream& fin);      // чтение матрицы в плотном формате.
	void readGg(Vector& gg, std::ifstream& fin);   // чтение массивов ggl и ggu.


	void convertDenseMatrix();                     // из плотного формата в разреженый строчно-стоблцовый. 

	void incompleteLUDecomposition(Vector& L, Vector& U, Vector& D);
	void straightMove(const Vector& L, const Vector& D, const Vector& v1, Vector& v2);
	void backMove(const Vector& U, const Vector& v1, Vector& v2);
	size_t BСGSTABLU();
	size_t losLU();


private:
	size_t n;                                      // размерность матрицы системы.
	size_t max_iter;                               // макисмальное количество итераций.
	size_t is_sparse;                              // нужно или не нужно конвертировать плотный формат в разреженный.

	double eps;                                    // точность вычислений.

	Matrix A_sparse;                               // матрица системы в разреженном строчно-столбцовом формате.
	Vector b;
	Vector x;

	std::vector<std::vector<double>> A;            // матрица в плотном формате (только для матриц небольших размерностей!)
};

size_t Solver::BСGSTABLU() {
	size_t n = A_sparse.di.size();

	Vector r(n), r_0(n), p(n), s(n), LAUp(n), LAUs(n);

	Vector L(A_sparse.ig[n]), U(A_sparse.ig[n]), D(n);

	double alpha, omega, betta, dotPrr_0r;
	double normb = sqrt(b * b);
	double residual;
	incompleteLUDecomposition(L, U, D);

	straightMove(L, D, b, r_0);
	r = r_0;
	backMove(U, r_0, p);

	size_t k = 1;
	for (; k < max_iter && (residual = sqrt(r * r) / normb) > eps; k++) {

		backMove(U, p, LAUp);
		LAUp = A_sparse * LAUp;
		straightMove(L, D, LAUp, LAUp);
		dotPrr_0r = r_0 * r;
		alpha = dotPrr_0r / (LAUp * r_0);

		s = r - alpha * LAUp;

		backMove(U, s, LAUs);
		LAUs = A_sparse * LAUs;
		straightMove(L, D, LAUs, LAUs);
		omega = (LAUs * s) / (LAUs * LAUs);

		x += alpha * p + omega * s;

		r = s - omega * LAUs;

		betta = alpha * (r * r_0) / (omega * dotPrr_0r);

		p = r + betta * (p - omega * LAUp);

	}
	backMove(U, x, x);

	return k;
}

size_t Solver::losLU() {
	size_t n = A_sparse.di.size();

	Vector r(n), z(n), p(n), LAUr(n), Ur(n);

	Vector L(A_sparse.ig[n]);
	Vector U(A_sparse.ig[n]);
	Vector D(n);

	double alpha, betta, dotPrpp;
	double normb = sqrt(b * b);
	double residual;
	incompleteLUDecomposition(L, U, D);

	straightMove(L, D, b, r);
	backMove(U, r, z);
	straightMove(L, D, A_sparse * z, p);

	size_t k = 1;
	for (; k < max_iter && (residual = sqrt(r * r) / normb) > eps; k++) {
		dotPrpp = p * p;
		alpha = (p * r) / dotPrpp;

		x += alpha * z;
		r -= alpha * p;

		backMove(U, r, Ur);
		LAUr = A_sparse * Ur;
		straightMove(L, D, LAUr, LAUr);

		betta = (-p * LAUr) / dotPrpp;
		z = Ur + betta * z;
		p = LAUr + betta * p;

	}
	return k;
}

void Solver::incompleteLUDecomposition(Vector& L, Vector& U, Vector& D) {
	size_t n = A_sparse.di.size();

	D = A_sparse.di;
	L = A_sparse.ggl;
	U = A_sparse.ggu;

	for (size_t i = 0; i < n; i++) {
		double d = 0;
		size_t temp = A_sparse.ig[i];
		for (size_t j = A_sparse.ig[i]; j < A_sparse.ig[i + 1]; j++) {
			double ls = 0;
			double us = 0;
			for (size_t h = A_sparse.ig[A_sparse.jg[j]], k = temp; h < A_sparse.ig[A_sparse.jg[j] + 1] && k < j;)
				if (A_sparse.jg[k] == A_sparse.jg[h]) {
					ls += L[k] * U[h];
					us += L[h++] * U[k++];
				}
				else (A_sparse.jg[k] < A_sparse.jg[h]) ? k++ : h++;

			L[j] -= ls;
			U[j] = (U[j] - us) / D[A_sparse.jg[j]];
			d += L[j] * U[j];
		}
		D[i] -= d;
	}
}

void Solver::straightMove(const Vector& L, const Vector& D, const Vector& v1, Vector& v2) {
	size_t n = A_sparse.di.size();
	for (size_t i = 0; i < n; i++) {
		double sum = 0;
		for (size_t j = A_sparse.ig[i]; j < A_sparse.ig[i + 1]; j++)
			sum += v2[A_sparse.jg[j]] * L[j];

		v2[i] = (v1[i] - sum) / D[i];
	}
}

void Solver::backMove(const Vector& U, const Vector& v1, Vector& v2) {
	size_t n = A_sparse.di.size();
	for (size_t i = 0; i < n; i++)
		v2[i] = v1[i];

	for (int i = n - 1; i >= 0; i--)
		for (size_t j = A_sparse.ig[i]; j < A_sparse.ig[i + 1]; j++)
			v2[A_sparse.jg[j]] -= v2[i] * U[j];

}

void Solver::readParams(const path& path) {
	std::ifstream fin(path / "isSparse.txt");
	if (fin.is_open()) { fin >> is_sparse; fin.close(); }
	else printErrorFile("isSparse.txt");

	fin.open(path / "params.txt");
	if (fin.is_open()) { fin >> eps >> max_iter; fin.close(); }
	else printErrorFile("params.txt");

	fin.open(path / "A_size.txt");
	if (fin.is_open()) { fin >> n; fin.close(); }
	else printErrorFile("A_size.txt");

	b.resize(n);
	fin.open(path / "b.txt");
	if (fin.is_open()) { readVector(b, fin); fin.close(); }
	else printErrorFile("b.txt");

	x.resize(n);
	fin.open(path / "init_x.txt");
	if (fin.is_open()) { readVector(x, fin); fin.close(); }
	else printErrorFile("init_x.txt");

	switch (is_sparse) {
	case 0:
		A.resize(n);
		fin.open(path / "A.txt");
		if (fin.is_open()) { readDenseMatrix(fin); fin.close(); }
		else printErrorFile("A.txt");

		convertDenseMatrix();
		break;
	case 1:
		A_sparse.di.resize(n);
		A_sparse.ig.resize(n + 1);
		readSparseMatrix(path);
		break;
	default:
		std::cerr << "There is no such matrix format. 0 - sparse format, 1 - dense format." << std::endl;
		std::exit(1);
		break;
	}
}

void Solver::convertDenseMatrix() {

}

void Solver::readSparseMatrix(const path& path) {
	std::ifstream fin(path / "di.txt");
	if (fin.is_open()) { readVector(A_sparse.di, fin); fin.close(); }
	else printErrorFile("di.txt");

	fin.open(path / "ig.txt");
	if (fin.is_open()) { readVector(A_sparse.ig, fin); fin.close(); }
	else printErrorFile("ig.txt");

	fin.open(path / "ggu.txt");
	if (fin.is_open()) { readGg(A_sparse.ggu, fin); fin.close(); }
	else printErrorFile("ggu.txt");

	fin.open(path / "ggl.txt");
	if (fin.is_open()) { readGg(A_sparse.ggl, fin); fin.close(); }
	else printErrorFile("ggl.txt");

	A_sparse.jg.resize(A_sparse.ggu.size());
	fin.open(path / "jg.txt");
	if (fin.is_open()) { readVector(A_sparse.jg, fin); fin.close(); }
	else printErrorFile("jg.txt");
}

void Solver::readDenseMatrix(std::ifstream& fin) {
	for (size_t i = 0; i < A.size(); i++) {
		for (size_t j = 0; j < A.size(); j++)
			fin >> A[i][j];
	}
}

void Solver::readGg(Vector& gg, std::ifstream& fin) {
	double value;
	while (fin >> value)
		gg.push_back(value);
}

Vector Solver::Matrix::operator * (const Vector& v) const {
	Vector result(v.size());

	for (size_t i = 0; i < v.size(); i++) {
		result[i] = di[i] * v[i];
		for (size_t k = ig[i]; k < ig[i + 1]; k++) {
			size_t j = jg[k];
			result[i] += ggl[k] * v[j];
			result[j] += ggu[k] * v[i];
		}
	}
	return result;
}
#endif // !_SOLVER_HPP
