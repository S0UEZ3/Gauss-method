#include "matrix.h"

Matrix::Matrix(const std::size_t n, const std::size_t m, const double default_value) :
	n(n), m(m), mtx(n, std::vector<double>(m, default_value)) {}

Matrix::Matrix(const std::size_t n, const std::size_t m) :
	n(n), m(m), mtx(n, std::vector<double> (m)) {}

Matrix::Matrix(const std::size_t n) :
	n(n), m(n), mtx(n, std::vector<double>(n,0)) {
	for (std::size_t i = 0; i < n; ++i) {
		mtx[i][i] = 1;
	}
}

double& Matrix::operator()(const std::size_t x, const std::size_t y) {
	if (x >= n || y >= m || x < 0 || y < 0) {
		throw std::out_of_range("Out of bounds");
	}
	return mtx[x][y];
}

Matrix& Matrix::operator+= (const Matrix& other) {
	if (n != other.n || m != other.m) {
		throw std::logic_error("Shapes not match");
	}
	for (std::size_t i = 0; i < n; ++i) {
		for (std::size_t j = 0; j < m; ++j) {
			mtx[i][j] += other.mtx[i][j];
		}
	}
	return *this;
};

Matrix Matrix::operator+ (const Matrix& other) const {
	Matrix summ = *this;
	summ += other;
	return summ;
};

Matrix& Matrix::operator++ () {
	if (n != m) {
		throw std::logic_error("Not square matrix");
	}
	for (std::size_t i = 0; i < n; ++i) {
		mtx[i][i]++;
	}
	return *this;
};

Matrix Matrix::operator++ (int) {
	Matrix copy = *this;
	++(*this);
	return copy;
};
 
bool Matrix::operator== (const Matrix& other) const {
	if (n != other.n || m != other.m) {
		return false;
	}
	for (std::size_t i = 0; i < n; ++i) {
		for (std::size_t j = 0; j < n; ++j) {
			if (mtx[i][j] != other.mtx[i][j]) {
				return false;
			}
		}
	}
	return true;
};

bool Matrix::operator!= (const Matrix & other) const {
	return !((*this) == other);
};


Matrix Matrix::operator- (const Matrix& other) const {
	Matrix diff = *this;
	diff -= other;
	return diff;
}


Matrix& Matrix::operator-= (const Matrix& other) {
	if (n != other.n || m != other.m) {
		throw std::logic_error("Shapes not match");
	}
	for (std::size_t i = 0; i < n; ++i) {
		for (std::size_t j = 0; j < m; ++j) {
			mtx[i][j] -= other.mtx[i][j];
		}
	}
	return *this;
}

Matrix Matrix::operator* (const Matrix& other) const {
	Matrix product = *this;
	product *= other;
	return product;
}

Matrix& Matrix::operator*= (const Matrix& other) {
	if (m != other.n) {
		throw std::logic_error("Shapes not match");
	}
	Matrix result(n, other.m, 0);
	for (std::size_t i = 0; i < n; ++i) {
		for (std::size_t j = 0; j < other.m; ++j) {
			for (std::size_t k = 0; k < m; ++k) {
				result.mtx[i][j] += mtx[i][k] * other.mtx[k][j];
			}
		}
	}
	*this = result;
	return *this;
}

std::pair<std::size_t, std::size_t> Matrix::size() const {
	return std::make_pair( n, m );
};


Matrix Matrix::pow(int s) const {
	if (n != m) {
		throw std::logic_error("Not square matrix");
	}
	Matrix result(n); 
	Matrix base = *this;
	while (s > 0) {
		if (s % 2 == 1) {
			result *= base;
		}
		base *= base;
		s /= 2;
	}
	return result;
}

Matrix Matrix::inverse() const {
	if (n != m) {
		throw std::logic_error("Not square matrix");
	}

	Matrix inv(n);
	Matrix gs(*this);

	// Initialize inv as identity matrix
	for (std::size_t i = 0; i < n; ++i) {
		inv(i, i) = 1.0;
	}

	for (std::size_t i = 0; i < n; ++i) {
		if (gs(i, i) == 0) {
			for (std::size_t j = i + 1; j < n; ++j) {
				if (gs(j, i) != 0) {
					for (std::size_t k = 0; k < n; ++k) {
						std::swap(gs(i, k), gs(j, k));
						std::swap(inv(i, k), inv(j, k));
					}
					break;
				}
			}
			if (gs(i, i) == 0) {
				throw std::logic_error("Singular matrix");
			}
		}

		double diag_val = gs(i, i);
		for (std::size_t j = 0; j < n; ++j) {
			gs(i, j) /= diag_val;
			inv(i, j) /= diag_val;
		}

		for (std::size_t j = 0; j < n; ++j) {
			if (i != j) {
				double f = gs(j, i);
				for (std::size_t k = 0; k < n; ++k) {
					gs(j, k) -= f * gs(i, k);
					inv(j, k) -= f * inv(i, k);
				}
			}
		}
	}
	return inv;
}

double Matrix::det() const {
	if (n != m) {
		throw std::logic_error("Not square matrix");
	}
	Matrix temp(*this);
	double determinant = 1.0;
	for (std::size_t i = 0; i < n; ++i) {
		if (temp(i, i) == 0) {
			for (std::size_t j = i + 1; j < n; ++j) {
				if (temp(j, i) != 0) {
					for (std::size_t k = 0; k < n; ++k) {
						std::swap(temp(i, k), temp(j, k));
					}
					determinant = -determinant;
					break;
				}
			}
			if (temp(i, i) == 0) {
				return 0;
			}
		}
		determinant *= temp(i, i);
		for (std::size_t j = i + 1; j < n; ++j) {
			double factor = temp(j, i) / temp(i, i);
			for (std::size_t k = i; k < n; ++k) {
				temp(j, k) -= factor * temp(i, k);
			}
		}
	}
	return determinant;
}

Matrix Matrix::t() const {
	Matrix tr(m, n);
	for (std::size_t i = 0; i < n; ++i) {
		for (std::size_t j = 0; j < m; ++j) {
			tr.mtx[j][i] = mtx[i][j];
		}
	}
	return tr;
};

std::ostream& operator << (std::ostream& os, const Matrix& matrix) {
	os << matrix.size().first << " " << matrix.size().second << '\n';
	for (std::size_t i = 0; i < matrix.size().first; ++i) {
		for (std::size_t j = 0; j < matrix.size().second; ++j) {
			os << matrix.mtx[i][j] << " ";
		}
		os << '\n';
	}
	return os;
}