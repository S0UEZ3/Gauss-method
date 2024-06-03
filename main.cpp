#include <assert.h>
#include <iostream>
#include <vector>
#include "matrix.h"
using namespace std;

int main() {
    double eps = 1e-5;
    {
        Matrix m1(3, 4, 1.0);
        assert(m1.size() == make_pair(3, 4));
        assert(abs(m1(0, 0) - 1.0) < eps);
        assert(abs(m1(2, 3) - 1.0) < eps);
        cout << "Constructor tests passed" << endl;
    }

    {
        Matrix m2(2, 3);
        assert(m2.size() == make_pair(2, 3));
        assert(abs(m2(0, 0)) < eps);
        assert(abs(m2(1, 2)) < eps);
        cout << "Constructor 2 tests passed" << endl;
    }

    {
        Matrix m3(3);
        assert(m3.size() == make_pair(3, 3));
        assert(abs(m3(0, 0) - 1.0) < eps);
        assert(abs(m3(1, 1) - 1.0) < eps);
        assert(abs(m3(2, 2) - 1.0) < eps);
        assert(abs(m3(0, 1)) < eps);
        assert(abs(m3(1, 0)) < eps);
        cout << "Identity matrix constructor tests passed" << endl;
    }

    {
        Matrix m4(2, 2, 2.0);
        m4(0, 0) = 5.0;
        assert(abs(m4(0, 0) - 5.0) < eps);
        cout << "Operator() tests passed" << endl;
    }

    {
        Matrix m5(2, 2, 1.0);
        Matrix m6(2, 2, 2.0);
        m5 += m6;
        assert(abs(m5(0, 0) - 3.0) < eps);
        cout << "Operator+= tests passed" << endl;
    }

    {
        Matrix m5(2, 2, 3.0);
        Matrix m6(2, 2, 2.0);
        Matrix m7 = m5 + m6;
        assert(abs(m7(0, 0) - 5.0) < eps);
        cout << "Operator+ tests passed" << endl;
    }

    {
        Matrix m8(2);
        ++m8;
        assert(abs(m8(0, 0) - 2.0) < eps);
        Matrix m9 = m8++;
        assert(abs(m9(0, 0) - 2.0) < eps);
        assert(abs(m8(0, 0) - 3.0) < eps);
        cout << "Operator++ tests passed" << endl;
    }

    {
        Matrix m10(2, 2, 1.0);
        Matrix m11(2, 2, 1.0);
        assert(m10 == m11);
        cout << "Operator== tests passed" << endl;
    }

    {
        Matrix m10(2, 2, 1.0);
        Matrix m11(2, 2, 2.0);
        assert(m10 != m11);
        cout << "Operator!= tests passed" << endl;
    }

    {
        Matrix m12(2, 2, 3.0);
        Matrix m13(2, 2, 1.0);
        Matrix m14 = m12 - m13;
        assert(abs(m14(0, 0) - 2.0) < eps);
        cout << "Operator- tests passed" << endl;
    }

    {
        Matrix m12(2, 2, 3.0);
        Matrix m13(2, 2, 1.0);
        m12 -= m13;
        assert(abs(m12(0, 0) - 2.0) < eps);
        cout << "Operator-= tests passed" << endl;
    }

    {
        Matrix m15(2, 3, 2.0);
        Matrix m16(3, 2, 3.0);
        Matrix m17 = m15 * m16;
        assert(m17.size() == make_pair(2, 2));
        assert(abs(m17(0, 0) - 18.0) < eps);
        cout << "Operator* tests passed" << endl;
    }

    {
        Matrix m15(2, 3, 2.0);
        Matrix m16(3, 2, 3.0);
        m15 *= m16;
        assert(abs(m15(0, 0) - 18.0) < eps);
        cout << "Operator*= tests passed" << endl;
    }

    {
        Matrix m18(2, 2, 2.0);
        Matrix m19 = m18.pow(3);
        assert(abs(m19(0, 0) - 32.0) < eps);
        cout << "Pow tests passed" << endl;
    }

    //{
    //    Matrix m20(2, 2);
    //    m20(0, 0) = 4; m20(0, 1) = 7;
    //    m20(1, 0) = 2; m20(1, 1) = 6;
    //    Matrix m21 = m20.inverse();
    //    assert(abs(m21(0, 0) - 0.6) < eps);
    //    assert(abs(m21(0, 1) + 0.7) < eps);
    //    assert(abs(m21(1, 0) - 0.2) < eps);
    //    assert(abs(m21(1, 1) - 0.4) < eps);
    //    cout << "Inverse tests passed" << endl;
    //}

    {
        Matrix m22(2, 2);
        m22(0, 0) = 4; m22(0, 1) = 7;
        m22(1, 0) = 2; m22(1, 1) = 6;
        double det = m22.det();
        assert(abs(det - 10.0) < eps);
        cout << "Determinant tests passed" << endl;
    }

    {
        Matrix m23(2, 3);
        m23(0, 0) = 1; m23(0, 1) = 2; m23(0, 2) = 3;
        m23(1, 0) = 4; m23(1, 1) = 5; m23(1, 2) = 6;
        Matrix m24 = m23.t();
        assert(m24.size() == make_pair(3, 2));
        assert(abs(m24(0, 0) - 1.0) < eps);
        assert(abs(m24(1, 0) - 2.0) < eps);
        assert(abs(m24(2, 0) - 3.0) < eps);
        assert(abs(m24(0, 1) - 4.0) < eps);
        assert(abs(m24(1, 1) - 5.0) < eps);
        assert(abs(m24(2, 1) - 6.0) < eps);
        cout << "Transpose tests passed" << endl;
    }

    cout << "All tests passed!" << endl;
}


