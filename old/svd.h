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
#ifndef SVD_H
#define SVD_H
#ifndef NO_OSTREAM
#include <iostream>
#endif

#include "../Constants.h"

namespace svd
{
    class SMat3
    {
    public:
        Real m00, m01, m02, m11, m12, m22;
    public:
        SMat3();
        SMat3(const Real m00, const Real m01, const Real m02,
              const Real m11, const Real m12, const Real m22);
        void clear() ;
        void setSymmetric(const Real m00, const Real m01, const Real m02,
                          const Real m11,
                          const Real m12, const Real m22) ;
        void setSymmetric(const SMat3 &rhs) ;
    private:
        SMat3(const SMat3 &rhs);
        SMat3 &operator=(const SMat3 &rhs);
    };
    class Mat3
    {
    public:
        Real m00, m01, m02, m10, m11, m12, m20, m21, m22;
    public:
        Mat3();
        Mat3(const Real m00, const Real m01, const Real m02,
             const Real m10, const Real m11, const Real m12,
             const Real m20, const Real m21, const Real m22);
        void clear() ;
        void set(const Real m00, const Real m01, const Real m02,
                 const Real m10, const Real m11, const Real m12,
                 const Real m20, const Real m21, const Real m22) ;
        void set(const Mat3 &rhs) ;
        void setSymmetric(const Real a00, const Real a01, const Real a02,
                          const Real a11, const Real a12, const Real a22) ;
        void setSymmetric(const SMat3 &rhs);
    private:
        Mat3(const Mat3 &rhs);
        Mat3 &operator=(const Mat3 &rhs);
    };
    class Vec3
    {
    public:
        Real x, y, z;
    public:
        Vec3();
        Vec3(const Real x, const Real y, const Real z);
        void clear();
        void set(const Real x, const Real y, const Real z);
        void set(const Vec3 &rhs);
    private:
        Vec3(const Vec3 &rhs);
        Vec3 &operator=(const Vec3 &rhs);
    };
#ifndef NO_OSTREAM
    std::ostream &operator<<(std::ostream &os, const Mat3 &m) ;
    std::ostream &operator<<(std::ostream &os, const SMat3 &m) ;
    std::ostream &operator<<(std::ostream &os, const Vec3 &v) ;
#endif
    class MatUtils
    {
    public:
        static Real fnorm(const Mat3 &a) ;
        static Real fnorm(const SMat3 &a) ;
        static Real off(const Mat3 &a) ;
        static Real off(const SMat3 &a) ;

    public:
        static void mmul(Mat3 &out, const Mat3 &a, const Mat3 &b) ;
        static void mmul_ata(SMat3 &out, const Mat3 &a) ;
        static void transpose(Mat3 &out, const Mat3 &a);
        static void vmul(Vec3 &out, const Mat3 &a, const Vec3 &v) ;
        static void vmul_symmetric(Vec3 &out, const SMat3 &a, const Vec3 &v) ;
    };
    class VecUtils
    {
    public:
        static void addScaled(Vec3 &v, const Real s, const Vec3 &rhs) ;
        static Real dot(const Vec3 &a, const Vec3 &b) ;
        static void normalize(Vec3 &v) ;
        static void scale(Vec3 &v, const Real s) ;
        static void sub(Vec3 &c, const Vec3 &a, const Vec3 &b) ;
    };
    class Givens
    {
    public:
        static void rot01_post(Mat3 &m, const Real c, const Real s);
        static void rot02_post(Mat3 &m, const Real c, const Real s);
        static void rot12_post(Mat3 &m, const Real c, const Real s);
    };
    class Schur2
    {
    public:
        static void rot01(SMat3 &out, Real &c, Real &s) ;
        static void rot02(SMat3 &out, Real &c, Real &s) ;
        static void rot12(SMat3 &out, Real &c, Real &s) ;
    };
    class Svd
    {
    public:
        static void getSymmetricSvd(const SMat3 &a, SMat3 &vtav, Mat3 &v,
                                    const Real tol, const int max_sweeps);
        static void pseudoinverse(Mat3 &out, const SMat3 &d, const Mat3 &v,
                                  const Real tol);
        static Real solveSymmetric(const SMat3 &A, const Vec3 &b, Vec3 &x,
                                    const Real svd_tol, const int svd_sweeps, const Real pinv_tol);
    };
    class LeastSquares
    {
    public:
        static Real solveLeastSquares(const Mat3 &a, const Vec3 &b, Vec3 &x,
                                       const Real svd_tol, const int svd_sweeps, const Real pinv_tol) ;
    };
};
#endif