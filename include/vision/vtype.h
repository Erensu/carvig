#ifndef VISION_TYPE_H
#define VISION_TYPE_H

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <vector>

namespace vig {
    template <typename _Scalar>
    using Quaternion=Eigen::Quaternion<_Scalar>;

    template <typename _Scalar>
    using Matrix3=Eigen::Matrix<_Scalar,3,3>;

    template <typename _Scalar>
    using Matrix4=Eigen::Matrix<_Scalar,4,4>;

    template <typename _Scalar>
    using MatrixX=Eigen::Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic>;

    template <typename _Scalar>
    using RowVector3=Eigen::Matrix<_Scalar,1,3>;

    template <typename _Scalar>
    using Vector2=Eigen::Matrix<_Scalar,2,1>;

    template <typename _Scalar>
    using Vector3=Eigen::Matrix<_Scalar,3,1>;

    template <typename _Scalar>
    using Vector4=Eigen::Matrix<_Scalar,4,1>;

    template <typename _Scalar>
    using VectorX=Eigen::Matrix<_Scalar,Eigen::Dynamic,1>;

    template <typename _Scalar>
    using Point=Vector3<_Scalar>;

    template <typename _Scalar>
    using GyroscopeReading=Vector3<_Scalar>;

    template <typename _Scalar>
    using AccelerometerReading=Vector3<_Scalar>;

    template <typename _Scalar>
    using Isometry3=Eigen::Transform<_Scalar,3,Eigen::Isometry>;

    template <typename _Scalar>
    struct camState {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        Point<_Scalar> p_C_G;
        Quaternion<_Scalar> q_CG;
        _Scalar time;
        int state_id;
        int last_correlated_id;
        std::vector<size_t> tracked_feature_ids;
    };

    template <typename _Scalar>
    struct featureTrackToResidualize {
        size_t feature_id;
        std::vector<Vector2<_Scalar>,
                Eigen::aligned_allocator<Vector2<_Scalar>>> observations;

        std::vector<camState<_Scalar>> cam_states;
        std::vector<size_t> cam_state_indices;

        bool initialized;
        Vector3<_Scalar> p_f_G;

        featureTrackToResidualize():initialized(false) {}
    };

    template <typename _Scalar>
    struct featureTrack {
        size_t feature_id;
        std::vector<Vector2<_Scalar>,
                Eigen::aligned_allocator<Vector2<_Scalar>>> observations;

        /* state_ids of cam states corresponding to observations */
        std::vector<size_t> cam_state_indices;

        bool initialized=false;
        Vector3<_Scalar> p_f_G;
    };
}
#endif