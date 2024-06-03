#pragma once
#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <vector>

class Matrix {
public:
	Matrix() = delete;
	Matrix(const std::size_t n, const std::size_t m, const double default_value);
	Matrix(const std::size_t n, const std::size_t m);
	Matrix(const std::size_t n);
	Matrix(const Matrix&) = default;

	double& operator()(const std::size_t x, const std::size_t y);
	Matrix operator+ (const Matrix& other) const;
	Matrix& operator++ ();
	Matrix operator++ (int);
	bool operator== (const Matrix& other) const;
	bool operator!= (const Matrix& other) const;
	Matrix& operator+= (const Matrix& other);
	Matrix operator- (const Matrix& other) const;
	Matrix& operator-= (const Matrix& other);
	Matrix operator* (const Matrix& other) const;
	Matrix& operator*= (const Matrix& other);

	std::pair<std::size_t, std::size_t> size() const;
	Matrix pow(int exp) const;
	Matrix inverse() const;
	double det() const;
	Matrix t() const;



private:
	std::size_t n, m;
	std::vector<std::vector<double>> mtx;

	friend std::ostream& operator << (std::ostream&, const Matrix&); 
};

std::ostream& operator << (std::ostream& os, const Matrix& matrix);