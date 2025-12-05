#pragma once
#include <vector>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <array>
#include <fstream>
#include <string>

namespace Data_structs {

    // Forward declaration
    struct Matrix3x3;

    struct Vector{
        double x,y,z;
        Vector(double x=0, double y=0, double z=0):x(x),y(y),z(z){}
        
        Vector operator+(const Vector& v) const {
            return Vector(x + v.x, y + v.y, z + v.z);
        }
        Vector operator-(const Vector& v) const {
            return Vector(x - v.x, y - v.y, z - v.z);
        }
        // V * scalar and scalar * V
        // friend Vector operator*(double scalar, const Vector& v) {
        //    return Vector(v.x * scalar, v.y * scalar, v.z * scalar);
        // }
        Vector operator*(double scalar) const {
            return Vector(x * scalar, y * scalar, z * scalar);
        }
        // friend Vector operator/(const Vector& v, double scalar) {
        //    return Vector(v.x / scalar, v.y / scalar, v.z / scalar);
        // }
        Vector operator/(double scalar) const {
            return Vector(x / scalar, y / scalar, z / scalar);
        }
        double dot(const Vector& v) const {
            return x * v.x + y * v.y + z * v.z;
        }
        double norm() const {
            return std::sqrt(x * x + y * y + z * z);
        }

        Matrix3x3 outer(const Vector& v) const;

        Vector normalize() const {
            double n = norm();
            return n > 0 ? (*this) / n : Vector(0, 0, 0);
        }
        Vector cross(const Vector& v) const {
            return Vector(
                y * v.z - z * v.y,
                z * v.x - x * v.z,
                x * v.y - y * v.x
            );
        }

        Vector& operator+=(const Vector& rhs) {
            x += rhs.x;
            y += rhs.y;
            z += rhs.z;
            return *this;
        }

        // elementwise subtraction
        Vector& operator-=(const Vector& rhs) {
            x -= rhs.x;
            y -= rhs.y;
            z -= rhs.z;
            return *this;
        }

        double& operator[](int i) {
            if(i==0) return x;
            if(i==1) return y;
            if(i==2) return z;
            throw std::out_of_range("Vector index out of range");
        }

        double operator[](int i) const {
            if(i==0) return x;
            if(i==1) return y;
            if(i==2) return z;
            throw std::out_of_range("Vector index out of range");
        }
    };

    // Define operator* (scalar * Vector) outside the class to ensure visibility
    inline Vector operator*(double scalar, const Vector& v) {
        return Vector(v.x * scalar, v.y * scalar, v.z * scalar);
    }

    struct Matrix3x3 {
        double m[3][3];
        
        static Matrix3x3 zero() {
            Matrix3x3 result;
            for(int i=0;i<3;++i)
                for(int j=0;j<3;++j)
                    result.m[i][j] = 0;
            return result;
        }

        Matrix3x3() {
            for(int i=0;i<3;++i)
                for(int j=0;j<3;++j)
                    m[i][j] = 0;
        }
        
        Matrix3x3 operator+(const Matrix3x3& mat) const {
            Matrix3x3 result;
            for(int i=0;i<3;++i)
                for(int j=0;j<3;++j)
                    result.m[i][j] = m[i][j] + mat.m[i][j];
            return result;
        }

        Matrix3x3 operator-(const Matrix3x3& mat) const {
            Matrix3x3 result;
            for(int i=0;i<3;++i)
                for(int j=0;j<3;++j)
                    result.m[i][j] = m[i][j] - mat.m[i][j];
            return result;
        }

        Matrix3x3 operator*(const Matrix3x3& mat) const {
            Matrix3x3 result;
            for(int i=0;i<3;++i)
                for(int j=0;j<3;++j)
                    for(int k=0;k<3;++k)
                        result.m[i][j] += m[i][k] * mat.m[k][j];
            return result;
        }

        Matrix3x3 operator*(double scalar) const {
            Matrix3x3 result;
            for(int i=0;i<3;++i)
                for(int j=0;j<3;++j)
                    result.m[i][j] = m[i][j] * scalar;
            return result;
        }

        Matrix3x3 operator/(double scalar) const {
            Matrix3x3 result;
            for(int i=0;i<3;++i)
                for(int j=0;j<3;++j)
                    result.m[i][j] = m[i][j] / scalar;
            return result;
        }

        Vector operator*(const Vector& v) const {
            return Vector(
                m[0][0]*v.x + m[0][1]*v.y + m[0][2]*v.z,
                m[1][0]*v.x + m[1][1]*v.y + m[1][2]*v.z,
                m[2][0]*v.x + m[2][1]*v.y + m[2][2]*v.z
            );
        }

        double determinant() const {
    double det =
        m[0][0]*(m[1][1]*m[2][2] - m[1][2]*m[2][1]) -
        m[0][1]*(m[1][0]*m[2][2] - m[1][2]*m[2][0]) +
        m[0][2]*(m[1][0]*m[2][1] - m[1][1]*m[2][0]);

    if (!std::isfinite(det)) {
        throw std::runtime_error("det is NaN or INF – upstream matrix corrupted");
    }
    return det;
}

        Matrix3x3 invert() const {
    double det = determinant();

    // Hard threshold gegen unterlaufende Determinanten
    const double det_min = 1e-12;

    if (!std::isfinite(det) || std::fabs(det) < det_min) {
        throw std::runtime_error("Matrix nearly singular");
    }

    Matrix3x3 inv;

    double idet = 1.0 / det;

    inv.m[0][0] =  (m[1][1]*m[2][2] - m[1][2]*m[2][1]) * idet;
    inv.m[0][1] = -(m[0][1]*m[2][2] - m[0][2]*m[2][1]) * idet;
    inv.m[0][2] =  (m[0][1]*m[1][2] - m[0][2]*m[1][1]) * idet;

    inv.m[1][0] = -(m[1][0]*m[2][2] - m[1][2]*m[2][0]) * idet;
    inv.m[1][1] =  (m[0][0]*m[2][2] - m[0][2]*m[2][0]) * idet;
    inv.m[1][2] = -(m[0][0]*m[1][2] - m[1][0]*m[0][2]) * idet;

    inv.m[2][0] =  (m[1][0]*m[2][1] - m[1][1]*m[2][0]) * idet;
    inv.m[2][1] = -(m[0][0]*m[2][1] - m[0][1]*m[2][0]) * idet;
    inv.m[2][2] =  (m[0][0]*m[1][1] - m[0][1]*m[1][0]) * idet;

    return inv;
}

        Matrix3x3 transpose() const {
            Matrix3x3 result;
            for(int i = 0; i < 3; ++i){
                for(int j = 0; j < 3; ++j){
                    result.m[i][j] = m[j][i];
                }
            }
        return result;
        }
    };

    // Implementation of outer() now that Matrix3x3 is defined
    inline Matrix3x3 Vector::outer(const Vector& v) const {
        Matrix3x3 result;
        result.m[0][0] = x * v.x; result.m[0][1] = x * v.y; result.m[0][2] = x * v.z;
        result.m[1][0] = y * v.x; result.m[1][1] = y * v.y; result.m[1][2] = y * v.z;
        result.m[2][0] = z * v.x; result.m[2][1] = z * v.y; result.m[2][2] = z * v.z;
        return result;
    }

    inline Matrix3x3 id() {
        Matrix3x3 result;
        result.m[0][0] = 1.0; result.m[0][1] = 0.0; result.m[0][2] = 0.0;
        result.m[1][0] = 0.0; result.m[1][1] = 1.0; result.m[1][2] = 0.0;
        result.m[2][0] = 0.0; result.m[2][1] = 0.0; result.m[2][2] = 1.0;
        return result;
    }

    inline Matrix3x3 zero() {
        Matrix3x3 result;
        for(int i=0;i<3;++i)
            for(int j=0;j<3;++j)
                result.m[i][j] = 0.0;
        return result;
    }

    struct Particle {
        Vector position;
        Vector velocity;
        Vector acceleration;
        double density;
        double drho_dt;
        double rho_pred;
        double smoothing_length;
        double sound_speed;
        double pressure;
        double mass;
        double shepard;
        Matrix3x3 correction_tensor;
        Matrix3x3 stress;
        Matrix3x3 dstress_dt;

        bool operator==(const Particle& other) const {
            return position.x == other.position.x && position.y == other.position.y && position.z == other.position.z;
        }
    };
}
