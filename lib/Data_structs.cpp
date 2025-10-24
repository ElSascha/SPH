#pragma once
#include <vector>
#include <cmath>

namespace Data_structs {


    struct Vector{
        double x,y,z;
        Vector(double x=0, double y=0, double z=0):x(x),y(y),z(z){}
        Vector operator+(const Vector& v) const {
            return Vector(x + v.x, y + v.y, z + v.z);
        }
        Vector operator-(const Vector& v) const {
            return Vector(x - v.x, y - v.y, z - v.z);
        }
        Vector operator*(double scalar) const {
            return Vector(x * scalar, y * scalar, z * scalar);
        }
        Vector operator/(double scalar) const {
            return Vector(x / scalar, y / scalar, z / scalar);
        }
        double dot(const Vector& v) const {
            return x * v.x + y * v.y + z * v.z;
        }
        double norm() const {
            return std::sqrt(x * x + y * y + z * z);
        }

        Matrix3x3 outer(const Vector& v) const {
        
             Matrix3x3 result;
             result.m[0][0] = x * v.x; result.m[0][1] = x * v.y; result.m[0][2] = x * v.z;
             result.m[1][0] = y * v.x; result.m[1][1] = y * v.y; result.m[1][2] = y * v.z;
             result.m[2][0] = z * v.x; result.m[2][1] = z * v.y; result.m[2][2] = z * v.z;
             return result;
        }

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

    };

    struct Matrix3x3 {
        double m[3][3];
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
            return 
                m[0][0]*(m[1][1]*m[2][2] - m[1][2]*m[2][1]) -
                m[0][1]*(m[1][0]*m[2][2] - m[1][2]*m[2][0]) +
                m[0][2]*(m[1][0]*m[2][1] - m[1][1]*m[2][0]);
        }

        Matrix3x3 invert(){
            if(determinant() == 0){
                throw std::runtime_error("Matrix is singular and cannot be inverted.");
            }
            Matrix3x3 inv; 
            double det = determinant();
            inv.m[0][0] = (m[1][1]*m[2][2] - m[1][2]*m[2][1])/det; inv.m[0][1] = (m[0][2]*m[2][1] - m[0][1]*m[2][2])/det; inv.m[0][2] = (m[0][1]*m[1][2] - m[0][2]*m[1][1])/det;
            inv.m[1][0] = (m[1][2]*m[2][0] - m[1][0]*m[2][2])/det; inv.m[1][1] = (m[0][0]*m[2][2] - m[0][2]*m[2][0])/det; inv.m[1][2] = (m[0][2]*m[1][0] - m[0][0]*m[1][2])/det;
            inv.m[2][0] = (m[1][0]*m[2][1] - m[1][1]*m[2][0])/det; inv.m[2][1] = (m[0][1]*m[2][0] - m[0][0]*m[2][1])/det; inv.m[2][2] = (m[0][0]*m[1][1] - m[0][1]*m[1][0])/det;
        }
    };

    struct Particle {
        Vector position;
        Vector velocity;
        Vector acceleration;
        double density;
        double pressure;
        double mass;
        double speed_of_sound;
        double shepard;
        Matrix3x3 correction_tensor;

        bool operator==(const Particle& other) const {
            return position.x == other.position.x && position.y == other.position.y && position.z == other.position.z;
        }
    };
}