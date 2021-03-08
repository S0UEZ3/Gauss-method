#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <vector>

class Rational {
private:
    int num = 0;
    int den = 1;

    int gcd() const {
        int x = std::abs(num);
        int y = std::abs(den);
        while (x > 0 && y > 0) {
            if (x > y) {
                x %= y;
            }
            else {
                y %= x;
            }
        }
        return x + y;
    }

    Rational& make_correct() {
        int numGcd = gcd();
        num /= numGcd;
        den /= numGcd;
        if (den < 0) {
            num *= -1;
            den *= -1;
        }
        return *this;
    }

public:
    Rational(int _num = 0, int _den = 1) : num(_num), den(_den) {
        int numGcd = gcd();
        num /= numGcd;
        den /= numGcd;
        if (den < 0) {
            num *= -1;
            den *= -1;
        }
    }

    Rational operator+(const Rational other) const {
        return {
            num * other.den + other.num * den,
            den * other.den
        };
    }

    Rational operator-(const Rational other) const {
        return {
            num * other.den - other.num * den,
            den * other.den
        };
    }

    Rational operator+() const {
        return { +num, den };
    }

    Rational operator-() const {
        return { -num, den };
    }

    Rational operator*(const Rational other) const {
        return {
            num * other.num,
            den * other.den
        };
    }

    Rational operator/(const Rational other) const {
        return {
            num * other.den,
            den * other.num
        };
    }

    Rational& operator+=(const Rational other) {
        num = num * other.den + other.num * den;
        den *= other.den;
        return make_correct();
    }

    Rational& operator-=(const Rational other) {
        num = num * other.den - other.num * den;
        den *= other.den;
        return make_correct();
    }

    Rational& operator*=(const Rational other) {
        num *= other.num;
        den *= other.den;
        return make_correct();
    }

    Rational& operator/=(const Rational other) {
        num *= other.den;
        den *= other.num;
        return make_correct();
    }

    int numerator() const {
        return num;
    }

    int denominator() const {
        return den;
    }

    bool operator==(const Rational other) const {
        return (num == other.num && den == other.den);
    }

    bool operator!=(const Rational other) const {
        return !(*this == other);
    }

    Rational& operator++() {
        num += den;
        return make_correct();
    }

    Rational operator++(int) {
        num += den;
        make_correct();
        return { num - den, den };
    }

    Rational& operator--() {
        num -= den;
        return make_correct();
    }

    Rational operator--(int) {
        num -= den;
        make_correct();
        return { num + den, den };
    }
};

Rational operator+(int num, Rational fraction) {
    return Rational(num) + fraction;
}

Rational operator-(int num, Rational fraction) {
    return Rational(num) + fraction;
}

Rational operator*(int num, Rational fraction) {
    return Rational(num) * fraction;
}

Rational operator/(int num, Rational fraction) {
    return Rational(num) / fraction;
}

bool operator==(int num, Rational fraction) {
    return Rational(num) == fraction;
}

bool operator!=(int num, Rational fraction) {
    return Rational(num) != fraction;
}

std::ostream& operator<<(std::ostream& out, const Rational& number) {
    if (number.denominator() == 1) {
        out << number.numerator();
    }
    else if (number.denominator() == 2) {
        out << (double)number.numerator() / (double)number.denominator();
    }
    else {
        out << number.numerator() << "/" << number.denominator();
    }
    return out;
}

std::istream& operator>>(std::istream& in, Rational& number) {
    int x;
    in >> x;
    number = Rational(x);
    return in;
}

using namespace std;

static const double ep = 1e-4;


size_t n, m, k;
vector<vector<Rational>> matrix(100, vector<Rational>(100));
const vector<string> RomeNum{ "I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X" };
bool was_smth = 1;

void devide(size_t i, Rational x) {
    if (x == 0) {
        return;
    }
    for (size_t j = 0; j < m + k; ++j) {
        matrix[i][j] /= x;
    }
    cout << RomeNum[i] + " = " + RomeNum[i] << " : " << x << endl;
    was_smth = 1;
}

void out() {
    was_smth = 0;
    cout << endl;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j) {
            cout << matrix[i][j] << ' ';
        }
        cout << "| ";
        for (size_t j = m; j < m + k; ++j) {
            cout << matrix[i][j] << ' ';
        }
        cout << endl;
    }
    cout << endl;
}

void subtract(size_t i, size_t i2, Rational x) {
    for (size_t j = 0; j < m + k; ++j) {
        matrix[i][j] -= x * matrix[i2][j];
    }
    if (x.numerator() < 0) {
        cout << RomeNum[i] + " = " + RomeNum[i] << " + " << -x << " * " << RomeNum[i2] << endl;
    }
    else {
        cout << RomeNum[i] + " = " + RomeNum[i] << " - " << x << " * " << RomeNum[i2] << endl;
    }
    was_smth = 1;
}

void swap_rows(size_t i, size_t i2) {
    if (i != i2) {
        cout << RomeNum[i] << " <=> " << RomeNum[i2] << endl;
        was_smth = 1;
    }
    for (size_t j = 0; j < m + k; ++j) {
        swap(matrix[i][j], matrix[i2][j]);
    }

}

bool find_and_swap(size_t i) {
    for (size_t j = i; j < n; ++j) {
        if (matrix[j][i] != 0) {
            swap_rows(i, j);
            return true;
        }
    }
    return false;
}

void out_answer() {
    cout << "answer:\n\n";
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < k; ++j) {
            if (matrix[i][i] == 0) {
                cout << "x" << i + 1 << j + 1 << "\t";
                continue;
            }
            if (matrix[i][j + m] == 0) {
                cout << "0";
            }
            else cout << matrix[i][j + m];
            for (size_t j2 = i + 1; j2 < m; ++j2) {
                if (matrix[i][j2] == 0) {
                    continue;
                }
                if (matrix[i][j2].numerator() < 0) {
                    cout << '+';
                }
                cout << -matrix[i][j2] << "*x" << j2 + 1 << j + 1;
            }

            cout << "\t";
        }
        cout << endl;
    }
}

void check() {
    for (size_t i = 0; i < n; ++i) {
        bool exist = 0;
        for (size_t j = 0; j < m; ++j) {
            if (matrix[i][j] != 0) {
                exist = 1;
            }
        }
        if (!exist) {
            for (size_t j = m; j < m + k; ++j) {
                if (matrix[i][j] != 0) {
                    cout << "No solution";
                    exit(0);
                }
            }
        }
    }
}

int main() {
    FILE* stream;
    freopen_s(&stream, "input.txt", "r", stdin);
    freopen_s(&stream, "output.txt", "w", stdout);
    cin >> n >> m >> k;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j) {
            string num;
            cin >> num;
            int f = 0, s = 0, t = 0, k = 1, zn = 1;
            if (num[0] == '-') zn = -1, num.erase(0, 1);
            for (auto x : num) {
                if (x == '/') t = 1;
                else if (x == '.') t = 2;
                else if (t == 0) f *= 10, f += x - '0';
                else s *= 10, s += x - '0', k *= 10;
            }
            if (t == 0) matrix[i][j] = Rational(f, 1);
            if (t == 1) matrix[i][j] = Rational(f, s);
            if (t == 2) matrix[i][j] = Rational(f * k + s, k);
            matrix[i][j] *= zn;
        }
    }
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = m; j < m + k; ++j) {
            string num;
            cin >> num;
            int f = 0, s = 0, t = 0, k = 1, zn = 1;
            if (num[0] == '-') zn = -1, num.erase(0, 1);
            for (auto x : num) {
                if (x == '/') t = 1;
                else if (x == '.') t = 2;
                else if (t == 0) f *= 10, f += x - '0';
                else s *= 10, s += x - '0', k *= 10;
            }
            if (t == 0) matrix[i][j] = Rational(f, 1);
            if (t == 1) matrix[i][j] = Rational(f, s);
            if (t == 2) matrix[i][j] = Rational(f * k + s, k);
            matrix[i][j] *= zn;
        }
    }
    out();
    check();
    for (size_t j = 0; j < m; ++j) {
        if (!find_and_swap(j)) {
            continue;
        }
        if (matrix[j][j] != 1) {
            devide(j, matrix[j][j]);
        }
        for (size_t i = j + 1; i < n; ++i) {
            if (matrix[i][j] != 0) {
                subtract(i, j, matrix[i][j]);
            }
        }
        if (was_smth) {
            out();
            check();
        }
    }
    cout << "--------------------\n\n";
    for (int j = m - 1; j >= 0; --j) {
        if (matrix[j][j] == 0) {
            continue;
        }
        for (int i = j - 1; i >= 0; --i) {
            if (matrix[i][j] != 0) {
                subtract(i, j, matrix[i][j]);
            }
        }
        if (was_smth) {
            out();
            check();
        }
    }

    out_answer();
}
