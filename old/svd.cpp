/*
 * This is free and unencumbered software released into the public domain.
 *
 * Anyone is free to copy, modify, publish, use, compile, sell, or
 * distribute this software, either in source code form or as a compiled
 * binary, for any purpose, commercial or non-commercial, and by any
 * means.
 *
 * In jurisdictions that recognize copyright laws, the author or authors
 * of this software dedicate any and all copyright interest in the
 * software to the public domain. We make this dedication for the benefit
 * of the public at large and to the detriment of our heirs and
 * successors. We intend this dedication to be an overt act of
 * relinquishment in perpetuity of all present and future rights to this
 * software under copyright law.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
 * OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 *
 * For more information, please refer to <http://unlicense.org/>
 */
#include "svd.h"
#include <math.h>

namespace svd
{

    Mat3::Mat3()
    {
        this->clear();
    }

    Mat3::Mat3(const Real m00, const Real m01, const Real m02,
               const Real m10, const Real m11, const Real m12,
               const Real m20, const Real m21, const Real m22)
    {
        this->set(m00, m01, m02, m10, m11, m12, m20, m21, m22);
    }

    Mat3::Mat3(const Mat3 &rhs)
    {
        this->set(rhs);
    }

    void Mat3::clear()
    {
        this->set(0, 0, 0, 0, 0, 0, 0, 0, 0);
    }

    void Mat3::set(const Real m00, const Real m01, const Real m02,
                   const Real m10, const Real m11, const Real m12,
                   const Real m20, const Real m21, const Real m22)
    {
        this->m00 = m00;
        this->m01 = m01;
        this->m02 = m02;
        this->m10 = m10;
        this->m11 = m11;
        this->m12 = m12;
        this->m20 = m20;
        this->m21 = m21;
        this->m22 = m22;
    }

    void Mat3::set(const Mat3 &rhs)
    {
        this->set(rhs.m00, rhs.m01, rhs.m02, rhs.m10, rhs.m11, rhs.m12, rhs.m20,
                  rhs.m21, rhs.m22);
    }

    void Mat3::setSymmetric(const SMat3 &rhs)
    {
        this->setSymmetric(rhs.m00, rhs.m01, rhs.m02, rhs.m11, rhs.m12, rhs.m22);
    }

    void Mat3::setSymmetric(const Real a00, const Real a01, const Real a02,
                            const Real a11, const Real a12, const Real a22)
    {
        this->set(a00, a01, a02, a01, a11, a12, a02, a12, a22);
    }

    SMat3::SMat3()
    {
        this->clear();
    }

    SMat3::SMat3(const Real m00, const Real m01, const Real m02,
                 const Real m11, const Real m12, const Real m22)
    {
        this->setSymmetric(m00, m01, m02, m11, m12, m22);
    }

    SMat3::SMat3(const SMat3 &rhs)
    {
        this->setSymmetric(rhs);
    }

    void SMat3::clear()
    {
        this->setSymmetric(0, 0, 0, 0, 0, 0);
    }

    void SMat3::setSymmetric(const SMat3 &rhs)
    {
        this->setSymmetric(rhs.m00, rhs.m01, rhs.m02, rhs.m11, rhs.m12, rhs.m22);
    }

    void SMat3::setSymmetric(const Real a00, const Real a01, const Real a02,
                             const Real a11, const Real a12, const Real a22)
    {
        this->m00 = a00;
        this->m01 = a01;
        this->m02 = a02;
        this->m11 = a11;
        this->m12 = a12;
        this->m22 = a22;
    }

    Vec3::Vec3() : x(0), y(0), z(0) { }

    Vec3::Vec3(const Vec3 &rhs)// : Vec3()
    {
        this->set(rhs);
    }

    Vec3::Vec3(const Real x, const Real y, const Real z)// : Vec3()
    {
        this->set(x, y, z);
    }

    void Vec3::clear()
    {
        this->set(0, 0, 0);
    }

    void Vec3::set(const Real x, const Real y, const Real z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    }

    void Vec3::set(const Vec3 &rhs)
    {
        this->set(rhs.x, rhs.y, rhs.z);
    }

#ifndef NO_OSTREAM
    std::ostream &operator<<(std::ostream &os, const Mat3 &m)
    {
        os << "[[" << m.m00 << ", " << m.m01 << ", " << m.m02 << "]" <<
        std::endl;
        os << " [" << m.m10 << ", " << m.m11 << ", " << m.m12 << "]" <<
        std::endl;
        os << " [" << m.m20 << ", " << m.m21 << ", " << m.m22 << "]]";
        return os;
    }

    std::ostream &operator<<(std::ostream &os, const SMat3 &m)
    {
        os << "[[" << m.m00 << ", " << m.m01 << ", " << m.m02 << "]" <<
        std::endl;
        os << " [" << m.m01 << ", " << m.m11 << ", " << m.m12 << "]" <<
        std::endl;
        os << " [" << m.m02 << ", " << m.m12 << ", " << m.m22 << "]]";
        return os;
    }

    std::ostream &operator<<(std::ostream &os, const Vec3 &v)
    {
        os << "[" << v.x << ", " << v.y << ", " << v.z << "]";
        return os;
    }
#endif

    Real MatUtils::fnorm(const Mat3 &a)
    {
        return sqrt((a.m00 * a.m00) + (a.m01 * a.m01) + (a.m02 * a.m02)
                    + (a.m10 * a.m10) + (a.m11 * a.m11) + (a.m12 * a.m12)
                    + (a.m20 * a.m20) + (a.m21 * a.m21) + (a.m22 * a.m22));
    }

    Real MatUtils::fnorm(const SMat3 &a)
    {
        return sqrt((a.m00 * a.m00) + (a.m01 * a.m01) + (a.m02 * a.m02)
                    + (a.m01 * a.m01) + (a.m11 * a.m11) + (a.m12 * a.m12)
                    + (a.m02 * a.m02) + (a.m12 * a.m12) + (a.m22 * a.m22));
    }

    Real MatUtils::off(const Mat3 &a)
    {
        return sqrt((a.m01 * a.m01) + (a.m02 * a.m02) + (a.m10 * a.m10)
                    + (a.m12 * a.m12) + (a.m20 * a.m20) + (a.m21 * a.m21));
    }

    Real MatUtils::off(const SMat3 &a)
    {
        return sqrt(2 * ((a.m01 * a.m01) + (a.m02 * a.m02) + (a.m12 * a.m12)));
    }

    void MatUtils::mmul(Mat3 &out, const Mat3 &a, const Mat3 &b)
    {
        out.set(a.m00 * b.m00 + a.m01 * b.m10 + a.m02 * b.m20,
                a.m00 * b.m01 + a.m01 * b.m11 + a.m02 * b.m21,
                a.m00 * b.m02 + a.m01 * b.m12 + a.m02 * b.m22,
                a.m10 * b.m00 + a.m11 * b.m10 + a.m12 * b.m20,
                a.m10 * b.m01 + a.m11 * b.m11 + a.m12 * b.m21,
                a.m10 * b.m02 + a.m11 * b.m12 + a.m12 * b.m22,
                a.m20 * b.m00 + a.m21 * b.m10 + a.m22 * b.m20,
                a.m20 * b.m01 + a.m21 * b.m11 + a.m22 * b.m21,
                a.m20 * b.m02 + a.m21 * b.m12 + a.m22 * b.m22);
    }

    void MatUtils::mmul_ata(SMat3 &out, const Mat3 &a)
    {
        out.setSymmetric(a.m00 * a.m00 + a.m10 * a.m10 + a.m20 * a.m20,
                         a.m00 * a.m01 + a.m10 * a.m11 + a.m20 * a.m21,
                         a.m00 * a.m02 + a.m10 * a.m12 + a.m20 * a.m22,
                         a.m01 * a.m01 + a.m11 * a.m11 + a.m21 * a.m21,
                         a.m01 * a.m02 + a.m11 * a.m12 + a.m21 * a.m22,
                         a.m02 * a.m02 + a.m12 * a.m12 + a.m22 * a.m22);
    }

    void MatUtils::transpose(Mat3 &out, const Mat3 &a)
    {
        out.set(a.m00, a.m10, a.m20, a.m01, a.m11, a.m21, a.m02, a.m12, a.m22);
    }

    void MatUtils::vmul(Vec3 &out, const Mat3 &a, const Vec3 &v)
    {
        out.x = (a.m00 * v.x) + (a.m01 * v.y) + (a.m02 * v.z);
        out.y = (a.m10 * v.x) + (a.m11 * v.y) + (a.m12 * v.z);
        out.z = (a.m20 * v.x) + (a.m21 * v.y) + (a.m22 * v.z);
    }

    void MatUtils::vmul_symmetric(Vec3 &out, const SMat3 &a, const Vec3 &v)
    {
        out.x = (a.m00 * v.x) + (a.m01 * v.y) + (a.m02 * v.z);
        out.y = (a.m01 * v.x) + (a.m11 * v.y) + (a.m12 * v.z);
        out.z = (a.m02 * v.x) + (a.m12 * v.y) + (a.m22 * v.z);
    }

    void VecUtils::addScaled(Vec3 &v, const Real s, const Vec3 &rhs)
    {
        v.x += s * rhs.x;
        v.y += s * rhs.y;
        v.z += s * rhs.z;
    }

    Real VecUtils::dot(const Vec3 &a, const Vec3 &b)
    {
        return a.x * b.x + a.y * b.y + a.z * b.z;
    }

    void VecUtils::normalize(Vec3 &v)
    {
        const Real len2 = VecUtils::dot(v, v);

        if (fabs(len2) < 1e-12) {
            v.clear();
        } else {
            VecUtils::scale(v, 1 / sqrt(len2));
        }
    }

    void VecUtils::scale(Vec3 &v, const Real s)
    {
        v.x *= s;
        v.y *= s;
        v.z *= s;
    }

    void VecUtils::sub(Vec3 &c, const Vec3 &a, const Vec3 &b)
    {
        const Real v0 = a.x - b.x;
        const Real v1 = a.y - b.y;
        const Real v2 = a.z - b.z;
        c.x = v0;
        c.y = v1;
        c.z = v2;
    }

    void Givens::rot01_post(Mat3 &m, const Real c, const Real s)
    {
        const Real m00 = m.m00, m01 = m.m01, m10 = m.m10, m11 = m.m11, m20 = m.m20,
                m21 = m.m21;
        m.set(c * m00 - s * m01, s * m00 + c * m01, m.m02, c * m10 - s * m11,
              s * m10 + c * m11, m.m12, c * m20 - s * m21, s * m20 + c * m21, m.m22);
    }

    void Givens::rot02_post(Mat3 &m, const Real c, const Real s)
    {
        const Real m00  = m.m00, m02  = m.m02, m10  = m.m10, m12  = m.m12,
                m20 = m.m20, m22  = m.m22 ;
        m.set(c * m00 - s * m02, m.m01, s * m00 + c * m02, c * m10 - s * m12, m.m11,
              s * m10 + c * m12, c * m20 - s * m22, m.m21, s * m20 + c * m22);
    }

    void Givens::rot12_post(Mat3 &m, const Real c, const Real s)
    {
        const Real m01  = m.m01, m02  = m.m02, m11  = m.m11, m12  = m.m12,
                m21  = m.m21, m22  = m.m22;
        m.set(m.m00, c * m01 - s * m02, s * m01 + c * m02, m.m10, c * m11 - s * m12,
              s * m11 + c * m12, m.m20, c * m21 - s * m22, s * m21 + c * m22);
    }

    static void calcSymmetricGivensCoefficients(const Real a_pp,
                                                const Real a_pq, const Real a_qq, Real &c, Real &s)
    {
        if (a_pq == 0) {
            c = 1;
            s = 0;
            return;
        }

        const Real tau = (a_qq - a_pp) / (2 * a_pq);
        const Real stt = sqrt(1.0f + tau * tau);
        const Real tan = 1.0f / ((tau >= 0) ? (tau + stt) : (tau - stt));
        c = 1.0f / sqrt(1.f + tan * tan);
        s = tan * c;
    }

    void Schur2::rot01(SMat3 &m, Real &c, Real &s)
    {
        svd::calcSymmetricGivensCoefficients(m.m00, m.m01, m.m11, c, s);
        const Real cc = c * c;
        const Real ss = s * s;
        const Real mix = 2 * c * s * m.m01;
        m.setSymmetric(cc * m.m00 - mix + ss * m.m11, 0, c * m.m02 - s * m.m12,
                       ss * m.m00 + mix + cc * m.m11, s * m.m02 + c * m.m12, m.m22);
    }

    void Schur2::rot02(SMat3 &m, Real &c, Real &s)
    {
        svd::calcSymmetricGivensCoefficients(m.m00, m.m02, m.m22, c, s);
        const Real cc = c * c;
        const Real ss = s * s;
        const Real mix = 2 * c * s * m.m02;
        m.setSymmetric(cc * m.m00 - mix + ss * m.m22, c * m.m01 - s * m.m12, 0,
                       m.m11, s * m.m01 + c * m.m12, ss * m.m00 + mix + cc * m.m22);
    }

    void Schur2::rot12(SMat3 &m, Real &c, Real &s)
    {
        svd::calcSymmetricGivensCoefficients(m.m11, m.m12, m.m22, c, s);
        const Real cc = c * c;
        const Real ss = s * s;
        const Real mix = 2 * c * s * m.m12;
        m.setSymmetric(m.m00, c * m.m01 - s * m.m02, s * m.m01 + c * m.m02,
                       cc * m.m11 - mix + ss * m.m22, 0, ss * m.m11 + mix + cc * m.m22);
    }

    static void rotate01(SMat3 &vtav, Mat3 &v)
    {
        if (vtav.m01 == 0) {
            return;
        }

        Real c, s;
        Schur2::rot01(vtav, c, s);
        Givens::rot01_post(v, c, s);
    }

    static void rotate02(SMat3 &vtav, Mat3 &v)
    {
        if (vtav.m02 == 0) {
            return;
        }

        Real c, s;
        Schur2::rot02(vtav, c, s);
        Givens::rot02_post(v, c, s);
    }

    static void rotate12(SMat3 &vtav, Mat3 &v)
    {
        if (vtav.m12 == 0) {
            return;
        }

        Real c, s;
        Schur2::rot12(vtav, c, s);
        Givens::rot12_post(v, c, s);
    }

    void Svd::getSymmetricSvd(const SMat3 &a, SMat3 &vtav, Mat3 &v,
                              const Real tol,
                              const int max_sweeps)
    {
        vtav.setSymmetric(a);
        v.set(1, 0, 0, 0, 1, 0, 0, 0, 1);
        const Real delta = tol * MatUtils::fnorm(vtav);

        for (int i = 0; i < max_sweeps
                        && MatUtils::off(vtav) > delta; ++i) {
            rotate01(vtav, v);
            rotate02(vtav, v);
            rotate12(vtav, v);
        }
    }

    static Real calcError(const Mat3 &A, const Vec3 &x,
                           const Vec3 &b)
    {
        Vec3 vtmp;
        MatUtils::vmul(vtmp, A, x);
        VecUtils::sub(vtmp, b, vtmp);
        return VecUtils::dot(vtmp, vtmp);
    }

    static Real calcError(const SMat3 &origA, const Vec3 &x,
                           const Vec3 &b)
    {
        Mat3 A;
        Vec3 vtmp;
        A.setSymmetric(origA);
        MatUtils::vmul(vtmp, A, x);
        VecUtils::sub(vtmp, b, vtmp);
        return VecUtils::dot(vtmp, vtmp);
    }

    static Real pinv(const Real x, const Real tol)
    {
        return (fabs(x) < tol || fabs(1 / x) < tol) ? 0 : (1 / x);
    }

    void Svd::pseudoinverse(Mat3 &out, const SMat3 &d, const Mat3 &v,
                            const Real tol)
    {
        const Real d0 = pinv(d.m00, tol), d1 = pinv(d.m11, tol), d2 = pinv(d.m22,
                                                                            tol);
        out.set(v.m00 * d0 * v.m00 + v.m01 * d1 * v.m01 + v.m02 * d2 * v.m02,
                v.m00 * d0 * v.m10 + v.m01 * d1 * v.m11 + v.m02 * d2 * v.m12,
                v.m00 * d0 * v.m20 + v.m01 * d1 * v.m21 + v.m02 * d2 * v.m22,
                v.m10 * d0 * v.m00 + v.m11 * d1 * v.m01 + v.m12 * d2 * v.m02,
                v.m10 * d0 * v.m10 + v.m11 * d1 * v.m11 + v.m12 * d2 * v.m12,
                v.m10 * d0 * v.m20 + v.m11 * d1 * v.m21 + v.m12 * d2 * v.m22,
                v.m20 * d0 * v.m00 + v.m21 * d1 * v.m01 + v.m22 * d2 * v.m02,
                v.m20 * d0 * v.m10 + v.m21 * d1 * v.m11 + v.m22 * d2 * v.m12,
                v.m20 * d0 * v.m20 + v.m21 * d1 * v.m21 + v.m22 * d2 * v.m22);
    }

    Real Svd::solveSymmetric(const SMat3 &A, const Vec3 &b, Vec3 &x,
                              const Real svd_tol, const int svd_sweeps, const Real pinv_tol)
    {
        Mat3 mtmp, pinv, V;
        SMat3 VTAV;
        Svd::getSymmetricSvd(A, VTAV, V, svd_tol, svd_sweeps);
        Svd::pseudoinverse(pinv, VTAV, V, pinv_tol);
        MatUtils::vmul(x, pinv, b);
        return svd::calcError(A, x, b);
    }

    Real
    LeastSquares::solveLeastSquares(const Mat3 &a, const Vec3 &b, Vec3 &x,
                                    const Real svd_tol, const int svd_sweeps, const Real pinv_tol)
    {
        Mat3 at;
        SMat3 ata;
        Vec3 atb;
        MatUtils::transpose(at, a);
        MatUtils::mmul_ata(ata, a);
        MatUtils::vmul(atb, at, b);
        return Svd::solveSymmetric(ata, atb, x, svd_tol, svd_sweeps, pinv_tol);
    }

}