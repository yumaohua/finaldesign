#ifndef BEM_RESONANCE_H
#define BEM_RESONANCE_H

#include "bem_matrix.h"
#include "bem_gauss_integrator.h"
#include <complex>
#include <tuple>
#include "../include/eigen-3.4.0/Eigen/Dense"
#include "../include/eigen-3.4.0/Eigen/SVD"
#include "bem_microcavity.h"
// #include "../include/eigen-3.4.0/spectra/include/Spectra/MatOp/DenseSymMatProd.h"
// #include "../include/eigen-3.4.0/spectra/include/Spectra/HermEigsSolver.h"

namespace bem {

// using Spectra::HermEigsSolver;
// using Spectra::LARGEST_MAGN;   // 按模最大的特征值（我们稍后取倒序）
// using Spectra::SMALLEST_MAGN;  // 按模最小的特征值
// using Spectra::CompInfo;
/**
 * ResonanceResult
 * 存储求得的共振波数 k_res 及对应的边界解 φ, ψ
 */
struct ResonanceResult {
    std::complex<double> k_res;  // 共振波数（复数）
    Eigen::VectorXcd phi;        // 边界法向导数 φ(s_l)
    Eigen::VectorXcd psi;        // 边界波函数 ψ(s_l)
};

inline double complex_compute_C_ll(
    const Microcavity::BEMElement& elem,
    double eps = EPS
) {
    double Δs = elem.length;
    double κ  = elem.curvature;
    return -1.0 + (κ * Δs) / (2.0 * M_PI);
}

inline std::complex<double> complex_compute_B_ll(
    const Microcavity::BEMElement& elem,
    std::complex<double> k,
    double n_j = 1.466,
    double eps = EPS
) {
    double Δs = elem.length;
    std::complex<double> z = n_j * k * (Δs/4.0);
    std::complex<double> lnz = std::log(z);
    return (Δs / M_PI) * (1.0 - lnz+ std::complex<double>(0.0, M_PI/2.0)-EULER_GAMMA);
}

inline
Eigen::MatrixXcd complex_build_B_matrix_block(
    const std::vector<Microcavity::BEMElement>& elems,
    Polarization pol,
    std::complex<double> k,
    double n_inner = 1.466,
    bool out = false
) {
    const int N = static_cast<int>(elems.size());
    Eigen::MatrixXcd B = Eigen::MatrixXcd::Zero(N, N);

    // 对每一对 (i,l) 元
    for (int i = 0; i < N; ++i) {
        const auto& ei = elems[i];
        for (int l = 0; l < N; ++l) {
            const auto& el = elems[l];

            if (i == l) {
                // 自作用项，用解析展开补偿 r→0 奇异
                B(i,l) = complex_compute_B_ll(el, k, /*n_j=*/n_inner);
                if(out) {
                   B(i,l)=-1.0*B(i,l);
                }
            } else {
                // 数值积分
                std::complex<double> I = GaussIntegrator::complex_integrate(
                    el,
                    ei.mid,
                    GaussIntegrator::complex_kernel_tm_B<CurveSegment>,
                    k,
                    n_inner
                );
                B(i,l) = I;
                if(out) {
                   B(i,l)=-1.0*B(i,l);
                }
            }
        }
        
    }
   //  flagB=false;
    return B;
}


inline
Eigen::MatrixXcd complex_build_C_matrix_block(
    const std::vector<Microcavity::BEMElement>& elems,
    Polarization pol,
    std::complex<double> k,
    double n_inner = 1.466,
    bool out=false
) {
    const int N = static_cast<int>(elems.size());
    Eigen::MatrixXcd C = Eigen::MatrixXcd::Zero(N, N);

    for (int i = 0; i < N; ++i) {
        const auto& ei = elems[i];
        for (int l = 0; l < N; ++l) {
            const auto& el = elems[l];

            if (i == l) {
                // 自作用，对角解析修正
                C(i,l) = std::complex<double>( complex_compute_C_ll(el), 0.0 );
                if (out)  {
                   C(i,l)=-1.0*C(i,l)-2.0;
                }
            } else {
                std::complex<double> I = GaussIntegrator::complex_integrate(
                    el,
                    ei.mid,
                    GaussIntegrator::complex_kernel_tm_C<CurveSegment>,
                    k,
                    n_inner
                );
                C(i,l) = I;
                if(out) {
                   C(i,l)=-1.0*C(i,l);
                }
            }
        }
    }
    return C;
}


inline
Eigen::MatrixXcd
complex_build_matrix(
    const std::vector<Microcavity::BEMElement>& all_elems,
    const std::vector<std::vector<Microcavity::BEMElement>>& per_cav,
    Polarization pol,
    std::complex<double> k,
    double n_inner = 1.466,
    double n_outer = 1.0
) {
    // 所有边界元的总数
    size_t N_el = all_elems.size();
    size_t N    = 2 * N_el;

    Eigen::MatrixXcd M  = Eigen::MatrixXcd::Zero(N, N);
    // Eigen::MatrixXcd M0 = Eigen::MatrixXcd::Zero(N, N);

    // 1) J-1 个腔体块
    size_t row = 0;
    for (auto const& elems : per_cav) {
        int n_j = static_cast<int>(elems.size());
        // 构建 B_j, C_j
        auto B_j = complex_build_B_matrix_block(elems, pol, k, n_inner);
        auto C_j = complex_build_C_matrix_block(elems, pol, k, n_inner);

        // 插入 B_j
        M.block(row, row, n_j, n_j) = B_j;
        // 插入 C_j
        M.block(row, N_el + row, n_j, n_j) = C_j;

        row += n_j;
    }

    // 2) 外部区域块
    //    注意折射率互换: 内→外, 外→内
    auto B_ext = complex_build_B_matrix_block(all_elems, pol, k, n_outer, true);
    auto C_ext = complex_build_C_matrix_block(all_elems, pol, k, n_outer, true);

    // 将它们插到底部，同时也写入 M0
    M .block(row,           0, N_el, N_el) = B_ext;
    M .block(row, N_el, N_el, N_el) = C_ext;
    // M0.block(row,           0, N_el, N_el) = B_ext;
    // M0.block(row, N_el, N_el, N_el) = C_ext;

    return M;
}

/**
 * detM
 * ----
 * 计算离散 BIE 矩阵 M(k) 的行列式 g(k)=det M(k)
 * 根据论文方程 (14)：只有 det M(k_res)=0 时才有非平凡解
 */
inline std::complex<double> detM(
    const std::vector<Microcavity::BEMElement>& all_elems,
    const std::vector<std::vector<Microcavity::BEMElement>>& per_cavity,
    Polarization pol,
    std::complex<double> k,
    double n_inner = 1.466,
    double n_outer = 1.0
) {
    Eigen::MatrixXcd M;
    M =
        complex_build_matrix(all_elems, per_cavity, pol, k, n_inner, n_outer);
    return M.determinant();
}

/**
 * detM_derivative
 * ----------------
 * 用论文方程 (36) 的数值近似计算 g'(k)=d/dk detM(k):
 *   g'(k) ≈ [g(k+Δ) - g(k)]/(2Δ) - i [g(k+iΔ) - g(k)]/(2Δ)
 */
inline std::complex<double> detM_derivative(
    const std::vector<Microcavity::BEMElement>& all_elems,
    const std::vector<std::vector<Microcavity::BEMElement>>& per_cavity,
    Polarization pol,
    std::complex<double> k,
    double delta = 1e-6,
    double n_inner = 1.466,
    double n_outer = 1.0
) {
    auto g0   = detM(all_elems, per_cavity, pol, k,             n_inner, n_outer);
    auto gre  = detM(all_elems, per_cavity, pol, k + delta,     n_inner, n_outer);
    auto gim  = detM(all_elems, per_cavity, pol, k + std::complex<double>(0,delta),
                     n_inner, n_outer);
    return (gre - g0)/(2.0*delta)
         - std::complex<double>(0,1)*(gim - g0)/(2.0*delta);
}

/**
 * find_resonance
 * ---------------
 * 牛顿迭代找零点 det M(k)=0，使用改进的“迹法”（论文方程 (37)）：
 *   k_{l+1} = k_l - q / tr[ M^{-1}(k_l) M'(k_l) ]
 *
 * @param all_elems    — 平坦的所有边界元列表
 * @param per_cavity   — 每个腔体对应的边界元列表（用于分块）
 * @param pol          — 偏振 (TM 或 TE)
 * @param k0           — 初始波数近似（α - i γ/2）
 * @param delta        — 数值微分步长 Δ
 * @param tol          — 收敛阈值 |k_{l+1}-k_l| < tol
 * @param max_iter     — 最大迭代次数
 * @param n_inner,n_outer — 内/外折射率
 * @return ResonanceResult 包含 k_res 及其对应的 φ, ψ
 */
inline ResonanceResult find_resonance(
    const std::vector<Microcavity::BEMElement>& all_elems,
    const std::vector<std::vector<Microcavity::BEMElement>>& per_cavity,
    Polarization pol,
    std::complex<double> k0,
    double delta = 1e-7,
    double tol = 1e-8,
    int max_iter = 200000,
    double n_inner = 1.466,
    double n_outer = 1.0
) {
    std::complex<double> k = k0;
    Eigen::MatrixXcd M, M0;

    for (int iter = 0; iter < max_iter; ++iter) {
        std::cout<<"begin build M("<<k<<")----"<<std::endl;
        // 1) 构造 M(k)
        M =
            complex_build_matrix(all_elems, per_cavity, pol, k, n_inner, n_outer);

        // 2) 数值微分 M'(k) 按 eq (36)
        //    M_forward = M(k+Δ), M_imag = M(k + iΔ)
        Eigen::MatrixXcd M_fwd, M_imag, dummy;
        std::cout<<"begin build M_fwd("<<k<<")----"<<std::endl;
        M_fwd=
            complex_build_matrix(all_elems, per_cavity, pol, k + delta, n_inner, n_outer);
        std::cout<<"begin build M_imag("<<k<<")----"<<std::endl;
        M_imag =
            complex_build_matrix(all_elems, per_cavity,
                         pol, k + std::complex<double>(0,delta),
                         n_inner, n_outer);
        std::cout<<"begin compute Mprime("<<k<<")----"<<std::endl;
        Eigen::MatrixXcd Mprime =
            (M_fwd  - M)/(2.0*delta)
          - std::complex<double>(0,1)*(M_imag - M)/(2.0*delta);

        std::cout<<"begin compute inverse----"<<std::endl;

        Eigen::MatrixXcd Minv = M.inverse();
        std::complex<double> trace = (Minv * Mprime).trace();
        std::complex<double> dk = 1.0 / trace;

        // 4) 更新 k
        std::complex<double> k1 = k - dk;

        if (std::abs(k1 - k) < tol) {
            k = k1;
            break;
        }
        std::cout<<"iter: "<<iter<<" , k: "<<k1<<std::endl;
        k = k1;
    }
    std::cout<<"finish"<<" , k: "<<k<<std::endl;

    // k 已收敛到 k_res
    ResonanceResult res;
    res.k_res = k;

    // // 用 SVD 求 M(k_res) 的零特征向量（奇异值最小对应的 V 列）
    // M =
    //     complex_build_matrix(all_elems, per_cavity, pol, k, n_inner, n_outer);
    // Eigen::JacobiSVD<Eigen::MatrixXcd> svd(
    //     M, Eigen::ComputeFullV
    // );
    // // 最小奇异值对应的 V 列在最后一列
    // Eigen::VectorXcd nullvec = svd.matrixV().col(M.cols()-1);

    // // 前 N 分量 = φ，后 N 分量 = ψ
    // int N = static_cast<int>(all_elems.size());
    // res.phi = nullvec.head(N);
    // res.psi = nullvec.tail(N);

    return res;
}
inline ResonanceResult find_resonance_secant(
    const std::vector<Microcavity::BEMElement>& all_elems,
    const std::vector<std::vector<Microcavity::BEMElement>>& per_cavity,
    Polarization pol,
    std::complex<double> k0,
    std::complex<double> delta_k0    = std::complex<double>(1e-3, 1e-3),
    double tol       = 1e-8,
    int    max_iter  = 200
) {
    using CD = std::complex<double>;

    // --- Helper: 构造 M(k) 并快速算 det(M) ---
    auto build_and_det = [&](const CD& kx, CD& det_out) {
        Eigen::MatrixXcd Mkx = complex_build_matrix(all_elems, per_cavity, pol,
                                                    kx, /*n_inner=*/1.466, /*n_outer=*/1.0);
        Eigen::PartialPivLU<Eigen::MatrixXcd> lu(Mkx);
        det_out = lu.determinant();
        return std::move(Mkx);
    };

    // 初始两点
    CD k_prev = k0;
    CD k      = k0 + delta_k0;
    CD f_prev, f;
    std::cout<<"begin build M("<<k_prev<<")"<<std::endl;
    Eigen::MatrixXcd M_prev = build_and_det(k_prev, f_prev);
    std::cout<<"begin build M("<<k<<")"<<std::endl;
    Eigen::MatrixXcd M      = build_and_det(k,      f);

    for (int iter = 0; iter < max_iter; ++iter) {
        // Secant 更新
        CD denom = (f - f_prev);
        if (std::abs(denom) < 1e-20) break;    // 防止除零
        CD k_new = k - f * (k - k_prev) / denom;

        // 判断收敛
        if (std::abs(k_new - k) < tol) {
            k_prev = k;  f_prev = f;
            k      = k_new;
            break;
        }

        // 轮换变量：k_prev←k, f_prev←f； k←k_new, f←det(M(k_new))
        k_prev = k;  f_prev = f;
        k      = k_new;
        std::cout<<"begin build M("<<k<<")"<<std::endl;
        M      = build_and_det(k, f);

        // （可选）调试输出
        std::cout << "iter=" << iter
                  << "  k=" << k
                  << "  |Δk|=" << std::abs(k - k_prev)
                  << "  |f|="   << std::abs(f)
                  << std::endl;
    }

    // 最终 k 上构造一次 M，并取最小奇异向量
    M = complex_build_matrix(all_elems, per_cavity, pol,
                             k, /*n_inner=*/1.466, /*n_outer=*/1.0);
    Eigen::JacobiSVD<Eigen::MatrixXcd> svd(M, Eigen::ComputeFullV);
    Eigen::VectorXcd nullvec = svd.matrixV().col(M.cols()-1);

    // 填充结果
    ResonanceResult res;
    res.k_res = k;
    int N = static_cast<int>(all_elems.size());
    res.phi = nullvec.head(N);
    res.psi = nullvec.tail(N);
    return res;
}
// inline ResonanceResult resonance(
//     const std::vector<Microcavity::BEMElement>& all_elems,
//     const std::vector<std::vector<Microcavity::BEMElement>>& per_cavity,
//     Polarization pol,
//     std::complex<double> k,
//     double n_inner = 1.466,
//     double n_outer = 1.0
// ) {
//     Eigen::MatrixXcd M=complex_build_matrix(all_elems, per_cavity, TM, k);

//     // k 已收敛到 k_res
//     ResonanceResult res;
//     res.k_res = k;

//     Eigen::JacobiSVD<Eigen::MatrixXcd> svd(
//         M, Eigen::ComputeFullV
//     );
//     // 最小奇异值对应的 V 列在最后一列
//     Eigen::VectorXcd nullvec = svd.matrixV().col(M.cols()-1);

//     // 前 N 分量 = φ，后 N 分量 = ψ
//     int N = static_cast<int>(all_elems.size());
//     res.phi = nullvec.head(N);
//     res.psi = nullvec.tail(N);

//     return res;
// }
inline ResonanceResult resonance(
    const std::vector<Microcavity::BEMElement>& all_elems,
    const std::vector<std::vector<Microcavity::BEMElement>>& per_cavity,
    Polarization pol,
    std::complex<double> k,
    double n_inner = 1.466,
    double n_outer = 1.0
) {
    // 1) 装配矩阵 M(k)
    Eigen::MatrixXcd M = complex_build_matrix(all_elems, per_cavity, pol, k, n_inner, n_outer);

    // 2) 填充返回结果的 k_res
    ResonanceResult res;
    res.k_res = k;

    // 3) 用 Divide‐and‐Conquer SVD（BDCSVD）只计算最小奇异向量
    //    需要在编译时定义 EIGEN_USE_LAPACKE，并链接 -llapacke -llapack -lblas
    Eigen::BDCSVD<Eigen::MatrixXcd> svd(
        M,
        Eigen::ComputeThinU | Eigen::ComputeThinV   // 只算 U,V 的薄版本
    );

    // 最小奇异值对应的 V 的最后一列
    Eigen::VectorXcd nullvec = svd.matrixV().col(M.cols() - 1);

    // 4) 拆分 φ 和 ψ
    int N = static_cast<int>(all_elems.size());
    res.phi = nullvec.head(N);
    res.psi = nullvec.tail(N);

    return res;
}
// 自定义 MatHMatOp：实现 y = (M^H * M) * x
struct MatHMatOp
{
    using Scalar     = std::complex<double>;
    using RealScalar = double;

    const Eigen::MatrixXcd& M;
    mutable Eigen::VectorXcd tmp;  // ← 加上 mutable

    MatHMatOp(const Eigen::MatrixXcd& mat)
        : M(mat), tmp(mat.cols()) {}

    int rows() const { return M.cols(); }
    int cols() const { return M.cols(); }

    void perform_op(const Scalar* x_in, Scalar* y_out) const
    {
        Eigen::Map<const Eigen::VectorXcd> x(x_in, cols());
        tmp.noalias() = M * x;                  // 现在允许在 const 函数中修改 tmp
        Eigen::Map<Eigen::VectorXcd> y(y_out, rows());
        y.noalias() = M.adjoint() * tmp;
    }
};

// // 原有的 ResonanceResult 结构体略…

// inline ResonanceResult resonance_spectra(
//     const std::vector<Microcavity::BEMElement>& all_elems,
//     const std::vector<std::vector<Microcavity::BEMElement>>& per_cavity,
//     Polarization pol,
//     std::complex<double> k,
//     double n_inner = 1.466,
//     double n_outer = 1.0
// ) {
//     // 1) 装配 M(k)
//     Eigen::MatrixXcd M = complex_build_matrix(all_elems, per_cavity, pol, k, n_inner, n_outer);

//     // 2) 准备 Spectra 操作对象
//     MatHMatOp op(M);

//     // 3) 构造 Hermitian 特征求解器：只要最小的 1 个特征值／向量
//     const int nev = 1;
//     const int ncv = 4 * nev + 1;  // 或者更大一些：std::min( M.cols(), 10*nev );
//     Spectra::HermEigsSolver<MatHMatOp> eigs(op, nev, ncv);

//     // 4) 初始化并计算，指定查找最小代数值特征值
//     eigs.init();
//     eigs.compute(Spectra::SortRule::SmallestAlge);

//     // 5) 检查收敛
//     if (eigs.info() != Spectra::CompInfo::Successful) {
//         throw std::runtime_error("Spectra 特征值求解未能收敛");
//     }

//     // 6) 提取 null vector
//     Eigen::VectorXcd nullvec = eigs.eigenvectors().col(0);

//     // 7) 拆分 φ 和 ψ
//     int N = static_cast<int>(all_elems.size());
//     ResonanceResult res;
//     res.k_res = k;
//     res.phi = nullvec.head(N);
//     res.psi = nullvec.tail(N);
//     return res;
// }
/// @brief 在单个点 r 处计算 |ψ(r)|²
/// @param res       已归一化的共振结果（包含 k_res, φ, ψ）
/// @param elements  边界元列表
/// @param r         待计算的空间点 (x,y)
/// @return          |ψ(r)|²
double compute_psi2_at_point_out(
    const ResonanceResult&                  res,
    const std::vector<Microcavity::BEMElement>& elements,
    const Point&                            r
) {
    std::complex<double> sum1 = 0.0, sum2 = 0.0;
    // 离散化公式 (38):
    // ψ(r) = ∑ ψ_l ∫_l ds ∂_ν G(s,r)  -  ∑ φ_l ∫_l ds G(s,r)
    for (size_t l = 0; l < elements.size(); ++l) {
        auto const& el = elements[l];
        // 1) ∫_l ∂ν G(s,r) ds
        auto I1 = GaussIntegrator::complex_integrate(
            el, r,
            GaussIntegrator::kernal_partial_G<CurveSegment>,
            res.k_res, /*n=*/1.0
        );
        // 2) ∫_l G(s,r) ds
        auto I2 = GaussIntegrator::complex_integrate(
            el, r,
            GaussIntegrator::kernel_G<CurveSegment>,
            res.k_res, /*n=*/1.0
        );
        sum1 += res.psi[l] * I1;
        sum2 += res.phi[l]   * I2;
    }
    std::complex<double> psi_r = sum1 - sum2;
    return std::norm(psi_r);  // 返回 |ψ(r)|^2
}


/// 2) 在圆周上（固定半径）按等分 θ 采样，写入 CSV
///
/// @param res           共振态结果
/// @param elements      边界元列表
/// @param radius        圆半径
/// @param n_theta       θ 分段数（总采样点 = n_theta）
/// @param theta0        起始角度（弧度）
/// @param theta1        结束角度（弧度）
/// @param out_filename  输出文件名
void compute_psi2_on_circle_txt(
    const ResonanceResult&      res,
    const std::vector<Microcavity::BEMElement>& elements,
    double                      radius,
    int                         n_theta,
    double                      theta0,
    double                      theta1,
    const std::string&          out_filename
) {
    std::ofstream fout(out_filename);
    fout << "theta,psi2\n";
    for (int i = 0; i < n_theta; ++i) {
        double theta = theta0 + (theta1 - theta0) * i / double(n_theta - 1);
        Point r{ radius * std::cos(theta),
                 radius * std::sin(theta) };
        double psi2 = compute_psi2_at_point_out(res, elements, r);
        fout << theta << "," << psi2 << "\n";
    }
    fout.close();
}

/// 将 psi 网格写入 TXT 文件
/// @param all_elems       所有边界元（平坦列表）
/// @param res             已归一化的共振结果（含 k_res, phi, psi）
/// @param cavities        每个腔体的多边形顶点列表（用于点内测试）
/// @param elems_per_cav   每个腔体各有多少个边界元，长度 = cavities.size()
/// @param nx, ny         网格的横纵分辨率
/// @param x0,x1,y0,y1    区域范围
/// @param filename       输出文件名
inline void compute_psi2_grid_txt(
    const std::vector<Microcavity::BEMElement>& all_elems,
    const ResonanceResult& res,
    const std::vector<std::vector<Point>>& cavities,
    const std::vector<int>& elems_per_cav,
    int nx, int ny,
    double x0, double x1,
    double y0, double y1,
    const std::string& filename
) {
    // 先构造每个腔体在 all_elems 中的区间 [start_j, end_j)
    std::vector<int> start_idx(cavities.size(), 0);
    for (size_t j = 1; j < cavities.size(); ++j) {
        start_idx[j] = start_idx[j-1] + elems_per_cav[j-1];
    }

    std::ofstream ofs(filename);
    // 对每个网格点
    for (int ix = 0; ix < nx; ++ix) {
        double x = x0 + (x1 - x0) * ix / (nx - 1);
        for (int iy = 0; iy < ny; ++iy) {
            double y = y0 + (y1 - y0) * iy / (ny - 1);
            Point r{x,y};

            // 1) 判定落入哪个腔体
            int which = -1;
            for (size_t j = 0; j < cavities.size(); ++j) {
                if (point_in_polygon(cavities[j], r)) {
                    which = int(j);
                    break;
                }
            }

            // 2) 对应的边界元区间
            int l0, l1;
            double index;
            if (which >= 0) {
                l0 = start_idx[which];
                l1 = l0 + elems_per_cav[which];
                index=1.466;
            } else {
                l0 = 0;
                l1 = int(all_elems.size());
                index=1.0;
            }

            // 3) 用公式(38)离散化计算 ψ(r)
            std::complex<double> sum1 = 0, sum2 = 0;
            for (int l = l0; l < l1; ++l) {
                const auto &el = all_elems[l];
                // ∫ ds ∂_ν G(s,r)
                auto I1 = GaussIntegrator::complex_integrate(
                    el, r,
                    GaussIntegrator::kernal_partial_G<CurveSegment>,
                    res.k_res, index
                );
                // ∫ ds G(s,r)
                auto I2 = GaussIntegrator::complex_integrate(
                    el, r,
                    GaussIntegrator::kernel_G<CurveSegment>,
                    res.k_res, index
                );
                sum1 += res.psi[l] * I1;
                sum2 += res.phi[l] * I2;
            }
            std::complex<double> psi_r = sum1 - sum2;
            double psi2 = std::norm(psi_r);//正负其实没关系

            ofs << x << "," << y << "," << psi2 << "\n";
        }
    }
    ofs.close();
}
// 解析 "(pr,pi),(sr,si)" 这一行
static bool parse_pair_line(
    const std::string& line,
    std::complex<double>& out_phi,
    std::complex<double>& out_psi
) {
    // 去掉所有 '(' 和 ')'
    std::string s;
    s.reserve(line.size());
    for (char c : line) {
        if (c != '(' && c != ')') s.push_back(c);
    }
    // 现在 s 应该形如 "pr,pi,sr,si"
    std::istringstream iss(s);
    double pr, pi, sr, si;
    char c1, c2, c3;
    if (!(iss >> pr >> c1 >> pi >> c2 >> sr >> c3 >> si)) {
        return false;
    }
    if (c1!=',' || c2!=',' || c3!=',') return false;
    out_phi = {pr, pi};
    out_psi = {sr, si};
    return true;
}

ResonanceResult load_and_normalize_txt(
    const std::string& filename,
    std::complex<double> k_res
) {
    std::ifstream ifs(filename);
    if (!ifs) throw std::runtime_error("无法打开文件 " + filename);

    std::vector<std::complex<double>> phi, psi;
    std::string line;
    while (std::getline(ifs, line)) {
        if (line.empty()) continue;
        std::complex<double> p, s;
        if (!parse_pair_line(line, p, s)) {
            throw std::runtime_error("解析失败，行格式不正确: " + line);
        }
        phi.push_back(p);
        psi.push_back(s);
    }
    if (phi.size() != psi.size()) {
        throw std::runtime_error("读取到的 φ 和 ψ 数量不一致");
    }

    // 找到 ψ 的最大模
    double max_psi = 0;
    for (auto& v : psi) {
        max_psi = std::max(max_psi, std::abs(v));
    }
    // 防止除以零
    if (max_psi < 1e-30) max_psi = 1.0;

    // 填充并归一化
    size_t N = phi.size();
    ResonanceResult res;
    res.k_res = k_res;
    res.phi .resize(N);
    res.psi .resize(N);
    for (size_t i = 0; i < N; ++i) {
        res.phi[i] = phi[i] / max_psi;
        res.psi[i] = psi[i] / max_psi;
    }
    return res;
}

} // namespace bem

#endif // BEM_RESONANCE_H
