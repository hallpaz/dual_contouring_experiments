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
#include "qef.h"
#include <stdexcept>
namespace svd
{

    QefData::QefData()
    {
        this->clear();
    }

    QefData::QefData(const Real ata_00, const Real ata_01,
                     const Real ata_02, const Real ata_11, const Real ata_12,
                     const Real ata_22, const Real atb_x, const Real atb_y,
                     const Real atb_z, const Real btb, const Real massPoint_x,
                     const Real massPoint_y, const Real massPoint_z,
                     const int numPoints)
    {
        this->set(ata_00, ata_01, ata_02, ata_11, ata_12, ata_22, atb_x, atb_y,
                  atb_z, btb, massPoint_x, massPoint_y, massPoint_z, numPoints);
    }

    QefData::QefData(const QefData &rhs)
    {
        this->set(rhs);
    }

    QefData& QefData::operator=(const QefData& rhs)
    {
        this->set(rhs);
        return *this;
    }

    void QefData::add(const QefData &rhs)
    {
        this->ata_00 += rhs.ata_00;
        this->ata_01 += rhs.ata_01;
        this->ata_02 += rhs.ata_02;
        this->ata_11 += rhs.ata_11;
        this->ata_12 += rhs.ata_12;
        this->ata_22 += rhs.ata_22;
        this->atb_x += rhs.atb_x;
        this->atb_y += rhs.atb_y;
        this->atb_z += rhs.atb_z;
        this->btb += rhs.btb;
        this->massPoint_x += rhs.massPoint_x;
        this->massPoint_y += rhs.massPoint_y;
        this->massPoint_z += rhs.massPoint_z;
        this->numPoints += rhs.numPoints;
    }

    void QefData::clear()
    {
        this->set(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    }

    void QefData::set(const Real ata_00, const Real ata_01,
                      const Real ata_02, const Real ata_11, const Real ata_12,
                      const Real ata_22, const Real atb_x, const Real atb_y,
                      const Real atb_z, const Real btb, const Real massPoint_x,
                      const Real massPoint_y, const Real massPoint_z,
                      const int numPoints)
    {
        this->ata_00 = ata_00;
        this->ata_01 = ata_01;
        this->ata_02 = ata_02;
        this->ata_11 = ata_11;
        this->ata_12 = ata_12;
        this->ata_22 = ata_22;
        this->atb_x = atb_x;
        this->atb_y = atb_y;
        this->atb_z = atb_z;
        this->btb = btb;
        this->massPoint_x = massPoint_x;
        this->massPoint_y = massPoint_y;
        this->massPoint_z = massPoint_z;
        this->numPoints = numPoints;
    }

    void QefData::set(const QefData &rhs)
    {
        this->set(rhs.ata_00, rhs.ata_01, rhs.ata_02, rhs.ata_11, rhs.ata_12,
                  rhs.ata_22, rhs.atb_x, rhs.atb_y, rhs.atb_z, rhs.btb,
                  rhs.massPoint_x, rhs.massPoint_y, rhs.massPoint_z,
                  rhs.numPoints);
    }

    QefData QefData::operator*(const Real &number)
    {
        return QefData(this->ata_00 * number, this->ata_01 * number, this->ata_02 * number,
                       this->ata_11 * number, this->ata_12 * number, this->ata_22 * number,
                       this->atb_x * number, this->atb_y * number, this->atb_z * number,
                       this->btb * number, this->massPoint_x, this->massPoint_y,
                       this->massPoint_z, this->numPoints);
    }
    QefData operator*(const Real &number, QefData &qefData)
    {
        return qefData * number;
    }

    QefData QefData::operator+(const QefData &rhs)
    {
        QefData newqef;
        newqef.add(*this);
        newqef.add(rhs);
        return newqef;
    }

#ifndef NO_OSTREAM
    std::ostream &operator<<(std::ostream &os, const QefData &qef)
    {
        SMat3 ata;
        Vec3 atb, mp;
        ata.setSymmetric(qef.ata_00, qef.ata_01, qef.ata_02, qef.ata_11,
                         qef.ata_12, qef.ata_22);
        atb.set(qef.atb_x, qef.atb_y, qef.atb_z);
        mp.set(qef.massPoint_x, qef.massPoint_y, qef.massPoint_z);

        if (qef.numPoints > 0) {
            VecUtils::scale(mp, 1.0f / qef.numPoints);
        }

        os << "QefData [ " << std::endl
        << " ata =" << std::endl <<  ata << "," << std::endl
        << " atb = " << atb << "," << std::endl
        << " btb = " << qef.btb << "," << std::endl
        << " massPoint = " << mp << "," << std::endl
        << " numPoints = " << qef.numPoints << "]";
        return os;
    }
#endif

    QefSolver::QefSolver() : data(), ata(), atb(), massPoint(), x(),
                             hasSolution(false) {}

    static void normalize(Real &nx, Real &ny, Real &nz)
    {
        Vec3 tmpv(nx, ny, nz);
        VecUtils::normalize(tmpv);
        nx = tmpv.x;
        ny = tmpv.y;
        nz = tmpv.z;
    }

    void QefSolver::add(const Real px, const Real py, const Real pz,
                        Real nx, Real ny, Real nz)
    {
        this->hasSolution = false;
        normalize(nx, ny, nz);
        this->data.ata_00 += nx * nx;
        this->data.ata_01 += nx * ny;
        this->data.ata_02 += nx * nz;
        this->data.ata_11 += ny * ny;
        this->data.ata_12 += ny * nz;
        this->data.ata_22 += nz * nz;
        const Real dot = nx * px + ny * py + nz * pz;
        this->data.atb_x += dot * nx;
        this->data.atb_y += dot * ny;
        this->data.atb_z += dot * nz;
        this->data.btb += dot * dot;
        this->data.massPoint_x += px;
        this->data.massPoint_y += py;
        this->data.massPoint_z += pz;
        ++this->data.numPoints;
    }

    void QefSolver::add(const Vec3 &p, const Vec3 &n)
    {
        this->add(p.x, p.y, p.z, n.x, n.y, n.z);
    }

    void QefSolver::add(const QefData &rhs)
    {
        this->hasSolution = false;
        this->data.add(rhs);
    }

    QefData QefSolver::getData()
    {
        return data;
    }

    Real QefSolver::getError()
    {
        if (!this->hasSolution) {
            throw std::runtime_error("illegal state");
        }

        return this->getError(this->x);
    }

    Real QefSolver::getError(const Vec3 &pos)
    {
        if (!this->hasSolution) {
            this->setAta();
            this->setAtb();
        }

        Vec3 atax;
        MatUtils::vmul_symmetric(atax, this->ata, pos);
        return VecUtils::dot(pos, atax) - 2 * VecUtils::dot(pos, this->atb)
               + this->data.btb;
    }

    void QefSolver::reset()
    {
        this->hasSolution = false;
        this->data.clear();
    }

    void QefSolver::setAta()
    {
        this->ata.setSymmetric(this->data.ata_00, this->data.ata_01,
                               this->data.ata_02, this->data.ata_11, this->data.ata_12,
                               this->data.ata_22);
    }

    void QefSolver::setAtb()
    {
        this->atb.set(this->data.atb_x, this->data.atb_y, this->data.atb_z);
    }

    Real QefSolver::solve(Vec3 &outx, const Real svd_tol,
                           const int svd_sweeps, const Real pinv_tol)
    {
        if (this->data.numPoints == 0) {
            throw std::invalid_argument("...");
        }

        this->massPoint.set(this->data.massPoint_x, this->data.massPoint_y,
                            this->data.massPoint_z);
        VecUtils::scale(this->massPoint, 1.0f / this->data.numPoints);
        this->setAta();
        this->setAtb();
        Vec3 tmpv;
        MatUtils::vmul_symmetric(tmpv, this->ata, this->massPoint);
        VecUtils::sub(this->atb, this->atb, tmpv);
        this->x.clear();
        const Real result = Svd::solveSymmetric(this->ata, this->atb,
                                                 this->x, svd_tol, svd_sweeps, pinv_tol);
        VecUtils::addScaled(this->x, 1, this->massPoint);
        this->setAtb();
        outx.set(x);
        this->hasSolution = true;
        return result;
    }

    void QefSolver::setData(QefData rhs)
    {
        this->data.set(rhs);
    }
}

