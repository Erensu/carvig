/*------------------------------------------------------------------------------
* math-unit.h : math unit functions
*
* version : $Revision:$ $Date:$
* history : 2018/11/10 1.0  new (first version)
*-----------------------------------------------------------------------------*/
#ifndef MATH_UNIT_H
#define MATH_UNIT_H
#include <vtype.h>

namespace vig {
    template <typename _Scalar>
    inline Matrix3<_Scalar> vectorToSkewSymmetric(const Vector3<_Scalar>& Vec) {

        /* Returns skew-symmetric form of a 3-d vector */
        Matrix3<_Scalar> M;

        M << 0     , -Vec(2),  Vec(1),
             Vec(2),       0, -Vec(0),
            -Vec(1),  Vec(0),       0;
        return M;
    }
    template <typename _Scalar>
    inline Matrix4<_Scalar> omegaMat(const Vector3<_Scalar>& omega) {
        /* Compute the omega-matrix of a 3-d vector omega */
        Matrix4<_Scalar> bigOmega;
        bigOmega.setZero();

        bigOmega.template block<3,3>(0,0)=-vectorToSkewSymmetric(omega);
        bigOmega.template block<3,1>(0,3)= omega;
        bigOmega.template block<1,3>(3,0)=-omega.transpose();
        return bigOmega;
    }
}
#endif