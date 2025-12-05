#include "Kernel.hpp"
#include <Eigen/Dense>


// Gradient of the kernel (stable version)
Vector gradW(Vector pos_i, Vector pos_j, double h, int dim) {
    Vector r_vec = pos_i - pos_j;
    double r_norm = r_vec.norm();
    if (r_norm < 1e-12) return Vector(0, 0, 0);  // avoid division by zero
    double factor = W_deriv(r_norm, h, dim) / r_norm;
    return r_vec * factor;
}


double W(double r, double h, int dim) {
    double sigma = 0.0;
    double q = r / h;

    if (dim == 3) {
        sigma = 495.0 / (256.0 * M_PI * h * h * h);
    } else {
        throw std::invalid_argument("Wendland kernel only implemented for 3D.");
    }

    if (q < 2.0) {
        double term = 1.0 - q / 2.0;
        return sigma * pow(term, 6.0) * (1.0 + 6.0 * (q / 2.0) + 35.0 / 3.0 * (q / 2.0) * (q / 2.0));
    } else {
        return 0.0;
    }
}


double W_deriv(double r, double h, int dim) {
    double sigma = 0.0;
    double q = r / h;

    if (dim == 3) {
        sigma = 495.0 / (256.0 * M_PI * h * h * h);
    } else {
        throw std::invalid_argument("Wendland kernel only implemented for 3D.");
    }

    if (q < 2.0) {
        double term = 1.0 - q / 2.0;
        double poly = 1.0 + 6.0 * (q / 2.0) + 35.0 / 3.0 * (q / 2.0) * (q / 2.0);
        double dpoly = 3.0 + 35.0 / 3.0 * (q / 2.0);
        return sigma * (-3.0 * pow(term, 5.0) * poly + pow(term, 6.0) * dpoly) / h;
    } else {
        return 0.0;
    }
}


void shepard_correction(std::vector<Particle>& particles, int dim) {
    if (dim != 3) {
        throw std::invalid_argument("Shepard correction only implemented for 3D.");
    }
    for (auto& pi : particles) {
        double shepard_sum = 0.0;
        for (auto& pj : particles) {
            Vector r = pi.position - pj.position;
            double r_norm = r.norm();
            shepard_sum += (pj.mass / pj.density) * W(r_norm, pi.smoothing_length, dim);
        }
        if (shepard_sum > 1e-12) {
            pi.shepard = 1.0 / shepard_sum;
        } else {
            pi.shepard = 1.0;
        }
    }
}


void consistent_shepard_interpolation(std::vector<Particle>& particles, int dim)
{
    if (dim != 3)
        throw std::invalid_argument("Shepard correction only implemented for 3D.");

    const int N = static_cast<int>(particles.size());
    if (N == 0) return;

    // Tolerances like in the paper
    const double eps_global = 1e-6;
    const double eps_local  = 1e-6;
    const double tiny = 1e-14;

    // Step 1: Build matrix A (N x N)
    Eigen::MatrixXd A(N, N);
    for (int i = 0; i < N; ++i) {
        double h_i = particles[i].smoothing_length;

        for (int j = 0; j < N; ++j) {
            Vector rij = particles[i].position - particles[j].position;
            double r = rij.norm();
            double w = W(r, h_i, dim);

            A(i, j) = particles[j].mass / std::max(particles[j].density, tiny) * w;
        }
    }

    // Step 2: Initialize c0 = 1
    Eigen::VectorXd c = Eigen::VectorXd::Ones(N);
    Eigen::VectorXd c_new = c;

    // Power iteration
    while (true)
    {
        // Step: compute A * c
        Eigen::VectorXd Ac = A * c;

        // Max-norm of A c  (∞-Norm)
        double norm_inf = Ac.cwiseAbs().maxCoeff();
        if (norm_inf < tiny)
            norm_inf = tiny;

        // Update c^{r+1}_i = (A c)_i / ||A c||_∞
        c_new = Ac / norm_inf;

        // Compute global + local error
        double global_error = (c_new - c).cwiseAbs().sum();
        double local_error  = (c_new - c).cwiseAbs().maxCoeff();

        // Convergence check
        if (global_error < eps_global && local_error < eps_local)
            break;

        // Continue iteration
        c = c_new;
    }

    // Write Shepard factors back
    for (int i = 0; i < N; ++i) {
        particles[i].shepard = 1.0/c_new(i);
    }
}

void tensor_correction(std::vector<Particle>& particles, int dim) {
    if (dim != 3) {
        throw std::invalid_argument("Tensor correction only implemented for 3D.");
    }
// For each particle, compute the correction tensor
    for (auto& pi : particles) {
        Eigen::Matrix3d L = Eigen::Matrix3d::Zero();

        for (auto& pj : particles) {
            if (&pi == &pj) continue;

            Vector r = pj.position - pi.position;
            Vector gradW_ij = gradW(pi.position, pj.position, pi.smoothing_length, dim);

            // (r * (m_j / rho_j)) outer gradW
            double vol_j = pj.mass / pj.density;
            Eigen::Vector3d rv(r.x * vol_j, r.y * vol_j, r.z * vol_j);
            Eigen::Vector3d gv(gradW_ij.x, gradW_ij.y, gradW_ij.z);

            L += rv * gv.transpose();
        }

        // SVD-based pseudo-inverse with thresholding
        Eigen::JacobiSVD<Eigen::Matrix3d> svd(L, Eigen::ComputeFullU | Eigen::ComputeFullV);
        Eigen::Matrix3d U = svd.matrixU();
        Eigen::Matrix3d V = svd.matrixV();
        Eigen::Vector3d S = svd.singularValues();

        double epsilon = 1e-4;  // threshold for singular values

        // Build Σ⁻¹ with thresholding
        Eigen::Matrix3d Sinv = Eigen::Matrix3d::Zero();
        for (int k = 0; k < 3; ++k) {
            if (S(k) > epsilon) {
                Sinv(k, k) = 1.0 / S(k);
            } else {
                Sinv(k, k) = 1.0;  // uncorrected direction
            }
        }

        Eigen::Matrix3d C = V * Sinv * U.transpose();

        // Convert to Matrix3x3
        Matrix3x3 correction;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                correction.m[i][j] = C(i, j);
            }
        }

        pi.correction_tensor = correction;
    }
}