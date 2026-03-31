#ifndef BEM_SCATTERING_H
#define BEM_SCATTERING_H

#include "bem_matrix.h"
#include "bem_gauss_integrator.h"
#include <vector>
#include <complex>
#include <cmath>
#include "../include/eigen-3.4.0/Eigen/Sparse"
#include "../include/eigen-3.4.0/Eigen/SVD"
#include "../include/eigen-3.4.0/Eigen/IterativeLinearSolvers"
#include "../include/eigen-3.4.0/unsupported/Eigen/IterativeSolvers"  // for IncompleteLUT
#include "bem_resonance.h"
// #include <mkl_lapacke.h>  // 添加 MKL 头文件

namespace bem {

// 存储散射解（边界上的 φ 和 ψ）
struct ScatteringResult {
    Eigen::VectorXcd phi;  // ∂νψ on boundary
    Eigen::VectorXcd psi;  // ψ on boundary
};

// // 用 long double 提升精度
// using ld      = long double;
// using Cld     = std::complex<ld>;
// using MatrixXcld = Eigen::Matrix<Cld, Eigen::Dynamic, Eigen::Dynamic>;
// using VectorXcld = Eigen::Matrix<Cld, Eigen::Dynamic, 1>;

// inline ScatteringResult solve_scattering(
//     const Eigen::VectorXcd &inc, 
//     const std::vector<Microcavity::BEMElement>& all_elems,
//     const std::vector<std::vector<Microcavity::BEMElement>>& per_cav,
//     Polarization pol, double k,
//     const Point& incident_dir = Point(1.0,0.0),
//     double n_inner = 1.466, double n_outer = 1.0
// ) {
//     // 1) 装配 double 精度矩阵
//     Eigen::MatrixXcd M_d, M0_d;
//     std::tie(M_d, M0_d) = build_matrix(all_elems, per_cav, pol, k, n_inner, n_outer);

//     // // 快速估条件数
//     // Eigen::FullPivLU<Eigen::MatrixXcd> lu_cond(M_d);
//     // double rcond = lu_cond.rcond();
//     // double condM = (rcond>0?1.0/rcond:std::numeric_limits<double>::infinity());
//     // std::cout << "[debug] approx cond(M) = " << condM << std::endl;

//     // 2) 构造 rhs（double）
//     Eigen::VectorXcd rhs_d = M0_d * inc;

//     // 3) 转到高精度 long double
//     int N = M_d.rows();
//     MatrixXcld M_hp(N,N);
//     VectorXcld rhs_hp(N);
//     for(int i=0;i<N;++i){
//         rhs_hp(i) = Cld{ (ld)rhs_d(i).real(), (ld)rhs_d(i).imag() };
//         for(int j=0;j<N;++j){
//             auto v = M_d(i,j);
//             M_hp(i,j) = Cld{ (ld)v.real(), (ld)v.imag() };
//         }
//     }

//     // 4) 行缩放：D_row 使每行范数≈1
//     Eigen::Matrix<ld, Eigen::Dynamic,1> row_norm(N);
//     for(int i=0;i<N;++i){
//         ld s = 0;
//         for(int j=0;j<N;++j) s += std::abs(M_hp(i,j));
//         row_norm(i) = (s>0?s:1.0L);
//     }
//     Eigen::DiagonalMatrix<Cld,Eigen::Dynamic> D_row(N);
//     for(int i=0;i<N;++i) D_row.diagonal()(i) = Cld{ 1.0L/row_norm(i), 0.0L };
//     M_hp   = D_row * M_hp;
//     rhs_hp = D_row * rhs_hp;

//     // 5) 列缩放：D_col 使每列范数≈1
//     Eigen::Matrix<ld, Eigen::Dynamic,1> col_norm(N);
//     for(int j=0;j<N;++j){
//         ld s = 0;
//         for(int i=0;i<N;++i) s += std::abs(M_hp(i,j));
//         col_norm(j) = (s>0?s:1.0L);
//     }
//     Eigen::DiagonalMatrix<Cld,Eigen::Dynamic> D_col(N);
//     for(int j=0;j<N;++j) D_col.diagonal()(j) = Cld{ 1.0L/col_norm(j), 0.0L };
//     M_hp = M_hp * D_col;
//     // 注意：rhs_hp 不做列缩放

//         // ——— 6) 估算缩放后矩阵的条件数 ———
//         {
//             Eigen::FullPivLU<decltype(M_hp)> lu_bal(M_hp);
//             ld rcond_bal = lu_bal.rcond();
//             ld cond_bal  = (rcond_bal > 0 ? 1.0L / rcond_bal : std::numeric_limits<ld>::infinity());
//             std::cout << "[debug] approx cond(balanced M) = " << (double)cond_bal << "\n";
//         }

//     // 6) 用高精度 LU 求解
//     Eigen::PartialPivLU<MatrixXcld> lu_hp(M_hp);
//     VectorXcld sol_hp = lu_hp.solve(rhs_hp);

//     // 7) 恢复原尺度： sol_true = D_col * sol_hp
//     VectorXcld sol_true_hp = D_col * sol_hp;

//     // 8) 转回 double
//     Eigen::VectorXcd sol_d(N);
//     for(int i=0;i<N;++i){
//         sol_d(i) = std::complex<double>(
//             (double)sol_true_hp(i).real(),
//             (double)sol_true_hp(i).imag()
//         );
//     }

//     // 9) 残差检查
//     double res_norm = (M_d * sol_d - rhs_d).norm();
//     std::cout << "[debug] residual ‖M·sol - rhs‖ = " << res_norm
//               << ", ‖rhs‖ = " << rhs_d.norm() << std::endl;

//     // 10) 拆分返回
//     size_t Nel = all_elems.size();
//     ScatteringResult R;
//     R.phi = sol_d.head(Nel);
//     R.psi = sol_d.tail(Nel);
//     return R;
// }


// // 提升精度用的别名
// using ld = long double;
// using Cld = std::complex<ld>;
// using MatrixXcld = Eigen::Matrix<Cld, Eigen::Dynamic, Eigen::Dynamic>;
// using VectorXcld = Eigen::Matrix<Cld, Eigen::Dynamic, 1>;

// inline ScatteringResult solve_scattering(
//     const Eigen::VectorXcd &inc,
//     const std::vector<Microcavity::BEMElement>& all_elems,
//     const std::vector<std::vector<Microcavity::BEMElement>>& per_cav,
//     Polarization pol,
//     double k,
//     const Point& incident_dir = Point(1.0, 0.0),
//     double n_inner = 1.466,
//     double n_outer = 1.0
// ) {
//     // 1) 装配 double 精度的 M, M0
//     Eigen::MatrixXcd M_d, M0_d;
//     std::tie(M_d, M0_d) = build_matrix(all_elems, per_cav, pol, k, n_inner, n_outer);

//     // 1a) 估算条件数（快速近似）
//     Eigen::FullPivLU<Eigen::MatrixXcd> lu_cond(M_d);
//     double rcond = lu_cond.rcond();
//     double condM = (rcond > 0 ? 1.0 / rcond : std::numeric_limits<double>::infinity());
//     std::cout << "[debug] approx cond(M) = " << condM << std::endl;

//     // 2) 构造 rhs（double）
//     Eigen::VectorXcd rhs_d = M0_d * inc;

//     // 3) 转到 long double 版 M_hp, rhs_hp
//     int N = M_d.rows();
//     MatrixXcld M_hp(N, N);
//     VectorXcld rhs_hp(N);
//     for (int i = 0; i < N; ++i) {
//         rhs_hp(i) = Cld((ld)rhs_d(i).real(), (ld)rhs_d(i).imag());
//         for (int j = 0; j < N; ++j) {
//             const auto &v = M_d(i, j);
//             M_hp(i, j) = Cld((ld)v.real(), (ld)v.imag());
//         }
//     }

//     // 4) 行缩放 (row‐scaling)：构造 D，使得新矩阵每行范数≈1
//     Eigen::Matrix<ld, Eigen::Dynamic, 1> row_norm(N);
//     for (int i = 0; i < N; ++i) {
//         ld s = 0;
//         for (int j = 0; j < N; ++j) {
//             s += std::abs(M_hp(i, j));
//         }
//         if (s < 1e-30L) s = 1.0L;
//         row_norm(i) = s;
//     }
//     Eigen::DiagonalMatrix<Cld, Eigen::Dynamic> D(N);
//     for (int i = 0; i < N; ++i) {
//         D.diagonal()(i) = Cld(1.0L / row_norm(i), 0.0L);
//     }
//     M_hp   = D * M_hp;
//     rhs_hp = D * rhs_hp;

//     // 5) long double 精度下列主元 LU 求解
//     Eigen::PartialPivLU<MatrixXcld> lu_hp(M_hp);
//     VectorXcld sol_hp = lu_hp.solve(rhs_hp);

//     // 6) 转回 double 精度
//     Eigen::VectorXcd sol(N);
//     for (int i = 0; i < N; ++i) {
//         sol(i) = std::complex<double>(
//             (double)sol_hp(i).real(),
//             (double)sol_hp(i).imag()
//         );
//     }

//     // 7) 原矩阵残差检查
//     double res_norm = (M_d * sol - rhs_d).norm();
//     std::cout << "[debug] residual ‖M·sol - rhs‖ = " << res_norm
//               << ", ‖rhs‖ = " << rhs_d.norm() << std::endl;

//     // 8) 拆分 φ 和 ψ 并返回
//     size_t Nel = all_elems.size();
//     ScatteringResult R;
//     R.phi = sol.head(Nel);
//     R.psi = sol.tail(Nel);
//     return R;
// }

/**
 * @brief 求解散射问题 (只用实数 k)
 *
 *   M [φ; ψ] = M0 [φ_in; ψ_in]
 */
inline
ScatteringResult solve_scattering(
    const Eigen::VectorXcd & inc,
    const std::vector<Microcavity::BEMElement>& all_elems,
    const std::vector<std::vector<Microcavity::BEMElement>>& per_cav,
    Polarization pol,
    double k,
    const Point& incident_dir = Point(1.0, 0.0),
    double n_inner = 1.466,
    double n_outer = 1.0
) {
    // 1) 装配矩阵
    Eigen::MatrixXcd M, M0;
    std::tie(M, M0) = build_matrix(all_elems, per_cav, pol, k, n_inner, n_outer);
    Eigen::VectorXd row_norms(all_elems.size());

    //     // ——— 更快的条件数估计 ———
    // // FullPivLU 会在分解时计算一个矩阵的近似秩，并给出 rcond() = 1/cond(A)
    // Eigen::FullPivLU<Eigen::MatrixXcd> lu_for_cond(M);
    // double rcond = lu_for_cond.rcond();
    // double condM = (rcond > 0 ? 1.0 / rcond : std::numeric_limits<double>::infinity());
    // std::cout << "[debug] approx cond(M) = " << condM << std::endl;

    // 2) 构造 RHS
    // Eigen::VectorXcd rhs = M0 * inc;
    Eigen::VectorXcd rhs = inc;
    // auto inc_ori=compute_incident_vectors(all_elems, k, incident_dir);
    // auto rhs_ori = M0*inc_ori;
    // {
    //     std::ofstream ofs("./data/rhs_compare.txt");
    //     ofs << "# idx    rhs_re    rhs_ori_re    rhs_im    rhs_ori_im\n";
    //     for (int i = 0; i < (int)rhs.size(); ++i) {
    //         ofs << i << " "
    //             << rhs[i].real()     << " "
    //             << rhs_ori[i].real() << " "
    //             << rhs[i].imag()     << " "
    //             << rhs_ori[i].imag()
    //             << "\n";
    //     }
    // }
    // 3) 列主元 LU 分解求解（PartialPivLU 更快一些）
    Eigen::PartialPivLU<Eigen::MatrixXcd> lu(M);
    Eigen::VectorXcd sol = lu.solve(rhs);
    
    std::cout << "‖rhs‖ = " << rhs.norm() 
              << std::endl;

    // 2) 构造入射向量 [φ_in; ψ_in]
    // Eigen::VectorXcd inc = compute_incident_vectors(all_elems, k, incident_dir);

    // 3) 计算 RHS
    // Eigen::VectorXcd rhs = M0 * inc;

    // 4) 更快的列主元 LU 分解
    // Eigen::PartialPivLU<Eigen::MatrixXcd> lu(M);
    // Eigen::VectorXcd sol = lu.solve(rhs);

    // Eigen::JacobiSVD<Eigen::MatrixXcd> svd(M, Eigen::ComputeThinU | Eigen::ComputeThinV);
    // Eigen::VectorXcd sol = svd.solve(rhs);

    // Eigen::VectorXcd sol = M.colPivHouseholderQr().solve(rhs);
          
    // 检查是否解成功
    // std::cout<<std::endl;
    // std::cout<<"norm: "<<(M * sol - rhs).norm() <<","<< rhs.norm()<<std::endl;
 
     // 5) 拆分结果
     size_t N_el = all_elems.size();
     ScatteringResult res;
     res.phi = sol.head(N_el);
     res.psi = sol.tail(N_el);
     return res;
 }

// /**
//  * @brief 求解散射问题 (只用实数 k)
//  *
//  *   M [φ; ψ] = M0 [φ_in; ψ_in]
//  */
//  inline
//  ScatteringResult solve_scattering(
//      const Eigen::VectorXcd & inc,
//      const std::vector<Microcavity::BEMElement>& all_elems,
//      const std::vector<std::vector<Microcavity::BEMElement>>& per_cav,
//      Polarization pol,
//      double k,
//      const Point& incident_dir = Point(1.0, 0.0),
//      double n_inner = 1.466,
//      double n_outer = 1.0
//  ) {
//      // 1) 装配矩阵
//      Eigen::MatrixXcd M, M0;
//      std::tie(M, M0) = build_matrix(all_elems, per_cav, pol, k, n_inner, n_outer);
 
//      // 2) 构造 RHS
//      Eigen::VectorXcd rhs = M0 * inc;
 
//      // 3) 预处理：对角缩放
//      Eigen::VectorXd row_norms = M.rowwise().norm();
//      Eigen::DiagonalMatrix<double, Eigen::Dynamic> D(row_norms.cwiseInverse());
//      Eigen::MatrixXcd M_pre = D * M;
//      Eigen::VectorXcd rhs_pre = D * rhs;
 
//      // 4) 使用 BiCGSTAB 求解
//      Eigen::BiCGSTAB<Eigen::MatrixXcd> solver;
//      solver.setTolerance(1e-4); // 设置收敛容差
//      solver.setMaxIterations(2000); // 设置最大迭代次数
//      solver.compute(M_pre);
//      std::cout<<std::endl;
//      std::cout<<3<<std::endl;
//      std::cout<<std::endl;

//      Eigen::VectorXcd sol_pre = solver.solve(rhs_pre);
//      Eigen::VectorXcd sol = D.inverse() * sol_pre; // 恢复原始解
//      std::cout<<std::endl;
//      std::cout<<4<<std::endl;
//      std::cout<<std::endl;
//      // 检查求解是否成功
//      if (solver.info() != Eigen::Success) {
//          std::cerr << "BiCGSTAB 求解失败！" << std::endl;
//      }
 
//      // 输出实际迭代次数
//      std::cout << "[debug] BiCGSTAB iterations: " << solver.iterations() << std::endl;
 
//      // 5) 残差检查
//      double res_norm = (M * sol - rhs).norm();
//      std::cout << "[debug] residual ‖M·sol - rhs‖ = " << res_norm
//                << ", ‖rhs‖ = " << rhs.norm() << std::endl;
 
//      // 6) 拆分结果
//      size_t N_el = all_elems.size();
//      ScatteringResult res;
//      res.phi = sol.head(N_el);
//      res.psi = sol.tail(N_el);
//      return res;
//  }     

/**
 * @brief 纯散射振幅积分核
 *
 * 对应文献公式 (19) 中被积分的那一段：
 *    integrand(s) = e^{−i k_f·r(s)} 
 *                   * { i (k_f·ν(s)) [ψ(s)−ψ_in(s)] + [φ(s)−φ_in(s)] }
 *
 * 参数说明：
 *   - kr       : 实数波数 k 的 real 部
 *   - kf_dir   : 出射方向单位向量 k_f / |k_f|
 *   - dpsi     : ψ(s) − ψ_in(s)（在当前元素上视为常数）
 *   - dphi     : φ(s) − φ_in(s)（在当前元素上视为常数）
 *   - s        : 当前积分点的坐标 r(s)
 *   - raw_n    : 曲线段给出的未归一化法向量 n(s)
 *
 * 返回值：
 *   该积分点处的被积函数值：phase * bracket
 */
 inline std::complex<double> scattering_kernel(
    double                  kr,
    const Point&            kf_dir,
    const std::complex<double>& dpsi,
    const std::complex<double>& dphi,
    const Point&            s,
    const Point&            raw_n
) {
    // 1) 单位外法线 ν(s)
    Point nu = -1.0*raw_n.normalized();

    // 2) 相位因子 exp(-i k_f·r)
    //    phase_arg = k * (k_f_unit · r)
    double phase_arg = kr * kf_dir.dot(s);
    std::complex<double> phase = std::exp(std::complex<double>(0.0, -phase_arg));

    // 3) 大括号内的那一项:
    //    i (k_f·ν) * dpsi  +  dphi
    std::complex<double> bracket =
        std::complex<double>(0.0, 1.0) * (kr*kf_dir.dot(nu)) * dpsi
      - dphi;

    // 4) 最终被积函数值
    return phase * bracket;
}


/**
 * @brief 计算纯散射振幅 f_scat(θ)（使用 15 点 Gauss–Legendre 积分）
 *
 * 根据文献公式 (19):
 *   f_scat(θ) = (1+i)/(4 √(π k)) ∮_{\partial Ω_J} ds  scattering_kernel(...)
 *
 * 本实现将边界积分统一使用 15 点 Gauss–Legendre 配权。
 */
 inline std::complex<double> compute_scattering_amplitude(
    const ScatteringResult&                sol,
    const std::vector<Microcavity::BEMElement>& elements,
    const Eigen::VectorXcd&                phi_psi_in,
    std::complex<double>                   k,
    const double&                          theta,
    const Point&                           incident_dir,
    bool                                   filter = false
) {
    int N = int(elements.size());
    // 拆分入射 φ_in, ψ_in
    Eigen::VectorXcd phi_in = phi_psi_in.head(N);
    Eigen::VectorXcd psi_in = phi_psi_in.tail(N);

    // 前置因子 (1+i)/(4 √(π k))
    double kr = k.real();
    std::complex<double> prefac = std::complex<double>(1.0, 1.0)
                               / (4.0 * std::sqrt(M_PI * kr));

    // 出射方向单位向量 k_f / |k_f|
    Point kf_dir{ std::cos(theta), std::sin(theta) };

    // 15 点 Gauss–Legendre 配权
    auto points = GaussIntegrator::get_points(15);

    std::complex<double> sum = 0.0;

    // 对 ∂Ω_J 上的每个元素做积分
    for (int ℓ = 0; ℓ < N; ++ℓ) {
        const auto& elem = elements[ℓ];
        const auto& seg  = *elem.seg_ptr;
        double t0 = elem.seg_start_t;
        double t1 = elem.seg_end_t;
        double dt = t1 - t0;

        for (auto& gp : points) {
            // 参数映射
            double t = t0 + gp.t * dt;
            Point  s        = seg(t);
            Point  nu       = 1.0 * seg.normal(t);     // 保证使用外法线
            double w        = gp.w * seg.tangent(t).norm() * dt;

            // 计算 ψ_in, φ_in
            Point k_dir = incident_dir.normalized();
            double dot_r = k_dir.x_ * s.x_ + k_dir.y_ * s.y_;
            std::complex<double> psi_i = std::exp(std::complex<double>(0.0, kr * dot_r));
            double nu_dot_k = nu.dot(k_dir);
            std::complex<double> phi_i = std::complex<double>(0.0,1.0) * kr * nu_dot_k * psi_i;

            // 散射部分
            std::complex<double> dpsi = sol.psi[ℓ] - psi_i;
            std::complex<double> dphi = sol.phi[ℓ] - phi_i;
            if (filter && std::abs(dphi) > 1e3) {
                continue;
            }

            // 积分核
            auto I = scattering_kernel(kr, kf_dir, dpsi, dphi, s, nu);
            sum += w * I;
        }
    }

    return prefac * sum;
}

/**
 * @brief 计算纯散射振幅 f_scat(θ)（使用 15 点 Gauss–Legendre 积分）
 *
 * 根据文献公式 (19):
 *   f_scat(θ) = (1+i)/(4 √(π k)) ∮_{\partial Ω_J} ds  scattering_kernel(...)
 *
 * 本实现将边界积分统一使用 15 点 Gauss–Legendre 配权。
 */
 inline std::complex<double> out_compute_scattering_amplitude(
    const ResonanceResult&                res,
    const std::vector<Microcavity::BEMElement>& elements,
    const double&                          theta
) {
    int N = int(elements.size());
    std::complex<double> prefac = std::complex<double>(1.0, 1.0)
                               / (4.0 * std::sqrt(M_PI * res.k_res));

    // 出射方向单位向量 k_f / |k_f|
    Point kf_dir{ std::cos(theta), std::sin(theta) };

    // 15 点 Gauss–Legendre 配权
    // auto points = GaussIntegrator::get_points(15);

    std::complex<double> sum = 0.0;

    // 对 ∂Ω_J 上的每个元素做积分
    for (int l = 0; l < N; ++l) {
        const auto& elem = elements[l];
        const auto& seg  = *elem.seg_ptr;
        double t0 = elem.seg_start_t;
        double t1 = elem.seg_end_t;
        double dt = t1 - t0;
        Point r = elem.mid;
        // 1) 单位外法线 ν(s)
        Point nu = -1.0*seg.normal(t0+dt/2.0).normalized();
        std::complex<double> phase = std::exp(std::complex<double>(0.0, -1)
    *res.k_res*kf_dir.dot(r));

        std::complex<double> bracket =
            std::complex<double>(0.0, 1.0) * (res.k_res*kf_dir.dot(nu)) * res.psi[l]
        - res.phi[l];

        // 4) 最终被积函数值
        auto I= phase * bracket;
        sum += I*elem.length;
    }

    return prefac * sum;
}
// /**
//  * @brief 计算纯散射振幅 f_scat(θ)
//  *
//  * 根据文献公式 (19)：  
//  *   f_scat(θ) = (1+i)/(4 √(π k)) ∮_{\partial Ω_J} ds  scattering_kernel(...)
//  *
//  * 在这里我们先对 φ_diff = φ − φ_in 进行“跳变”检测：如果某个元素的 |φ_diff[l]|
//  * 远大于其左右相邻元素的平均值，就用左右两点线性插值代替该点的 φ_diff[l]，
//  * 以去除孤立的尖峰。然后再带入高斯积分。
//  */
//  inline std::complex<double> compute_scattering_amplitude(
//     const ScatteringResult&                sol,
//     const std::vector<Microcavity::BEMElement>& elements,
//     const Eigen::VectorXcd&                phi_psi_in,
//     std::complex<double>                   k,
//     const double&             theta,
//     const Point&                           incident_dir,
//     bool filter = false
// ) {
//     int N = int(elements.size());
//     // 拆分入射 φ_in, ψ_in
//     Eigen::VectorXcd phi_in = phi_psi_in.head(N);
//     Eigen::VectorXcd psi_in = phi_psi_in.tail(N);

//     // 1) 构造原始差分
//     std::vector<std::complex<double>> phi_diff(N), psi_diff(N);
//     for (int l = 0; l < N; ++l) {
//         phi_diff[l] = sol.phi[l] - phi_in[l];
//         psi_diff[l] = sol.psi[l] - psi_in[l];
//     }

//     // 2) 对 φ_diff 做“跳变”检测与平滑
//     //    如果 |φ_diff[l]| > thresh * 0.5*(|φ_diff[l−1]| + |φ_diff[l+1]|)
//     //    就令 φ_diff[l] = 0.5*(φ_diff[l−1] + φ_diff[l+1])
//     const double thresh = 10.0;  // 阈值因子，可根据需要调整
//     std::vector<std::complex<double>> phi_smooth = phi_diff;
//     for (int l = 0; l < N/2; ++l) {
//         int lm = (l == 0   ? N/2-1 : l-1);
//         int lp = (l == N/2-1 ? 0   : l+1);
//         double dm=(elements[lm].length+elements[l].length)/2.0;
//         double dp=(elements[lp].length+elements[l].length)/2.0;
//         double neigh_avg =  (dp*std::abs(phi_diff[lm]) + dm*std::abs(phi_diff[lp]))/(dp+dm);
//         if (neigh_avg > 0.0 && std::abs(phi_diff[l]) > thresh * neigh_avg) {
//             // 发生“跳变”，用两侧线性插值
//             phi_smooth[l] = (dp*phi_diff[lm] + dm*phi_diff[lp])/(dp+dm);
//         }
//     }
//     for (int l = N/2; l < N; ++l) {
//         int lm = (l == N/2   ? N-1 : l-1);
//         int lp = (l == N-1 ? N/2   : l+1);
//         double dm=(elements[lm].length+elements[l].length)/2.0;
//         double dp=(elements[lp].length+elements[l].length)/2.0;
//         double neigh_avg =  (dp*std::abs(phi_diff[lm]) + dm*std::abs(phi_diff[lp]))/(dp+dm);
//         if (neigh_avg > 0.0 && std::abs(phi_diff[l]) > thresh * neigh_avg) {
//             // 发生“跳变”，用两侧线性插值
//             phi_smooth[l] = (dp*phi_diff[lm] + dm*phi_diff[lp])/(dp+dm);
//         }
//     }

//     // 3) 准备前置因子
//     double kr = k.real();
//     std::complex<double> prefac = std::complex<double>(1.0,1.0)
//                                / (4.0 * std::sqrt(M_PI * kr));
//     Point kf_dir{ std::cos(theta), std::sin(theta) };
//     // auto points = GaussIntegrator::get_points(15);

//     std::complex<double> sum = 0.0;
//     for (int l = 0; l < N; ++l) {
//         auto part1 = std::exp(std::complex<double>(0,-1)*k*kf_dir.dot(elements[l].mid));
//         double t0=elements[l].seg_start_t;
//         double t1=elements[l].seg_end_t;
//         auto part2=-std::complex<double>(0,1)*k*kf_dir.dot(elements[l].seg_ptr->normal((t0+t1)/2.0));
//         sum+=part1*(part2*psi_diff[l]-phi_smooth[l])*elements[l].length;
//     }

//     // 5) 乘以前置因子
//     return prefac * sum;
// }

// /**
//  * @brief 计算纯散射振幅 f_scat(θ)
//  *
//  * 根据文献公式 (19)：
//  *  f_scat(θ) = (1+i)/(4 √(π k)) ∮_{\partial Ω_J} ds  scattering_kernel(...)
//  *
//  * @param sol            BEM 求得的散射解（φ, ψ 向量）
//  * @param elements       ∂Ω_J 上所有边界元
//  * @param phi_psi_in     入射场 [φ_in; ψ_in] 向量
//  * @param k              波数（complex，但这里只用 real 部）
//  * @param thetas         观测角度列表 θ
//  * @param incident_dir   入射方向单位矢量（仅用于与 θ=φ 对齐）
//  * @return               对应每个 θ 的散射振幅 f_scat(θ)
//  */
// inline std::complex<double> compute_scattering_amplitude(
//     const ScatteringResult&                sol,
//     const std::vector<Microcavity::BEMElement>& elements,
//     const Eigen::VectorXcd&                phi_psi_in,
//     std::complex<double>                   k,
//     const double&             theta,
//     const Point&                           incident_dir,
//     bool filter = false
// ) {
//     int N = int(elements.size());
//     // 拆分入射 φ_in, ψ_in
//     Eigen::VectorXcd phi_in = phi_psi_in.head(N);
//     Eigen::VectorXcd psi_in = phi_psi_in.tail(N);

//     // 前置因子 (1+i)/(4√(π k))
//     double kr = k.real();
//     std::complex<double> prefac = std::complex<double>(1.0,1.0)
//                                / (4.0 * std::sqrt(M_PI * kr));

//     // std::vector<std::complex<double>> f_scat;

//     // 出射方向单位向量 k_f / |k_f|
//     Point kf_dir{ std::cos(theta), std::sin(theta) };

//     std::complex<double> sum = 0.0;
//     auto points = GaussIntegrator::get_points(15);
//     // 对 ∂Ω_J 上的每个元素做 Gauss 积分
//     for (int l = 0; l < N; ++l) {
//         auto const& elem = elements[l];

//         // 散射部分 ψ−ψ_in, φ−φ_in
//         // std::complex<double> dpsi = sol.psi[ℓ] - psi_in[ℓ];
//         // std::complex<double> dphi = sol.phi[ℓ] - phi_in[ℓ];

//         // if(filter) {
//         //     if (std::abs(dphi)>500) {
//         //         continue;
//         //     }
//         // }
//         const auto& seg = *elem.seg_ptr;
//         double t0 = elem.seg_start_t, t1 = elem.seg_end_t;
//         for (auto& gp : points) {
//             double t = t0 + gp.t*(t1 - t0);
//             Point s   = seg(t);
//             // double mid_t = 0.5 * (elem.seg_start_t + elem.seg_end_t);
//             Point nu = -1.0*elem.seg_ptr->normal(t);
//             // 在该元素上执行数值积分
//             Point r = seg(t);
//             Point k_dir = incident_dir.normalized();
//             double dot_ = k_dir.x_ * r.x_ + k_dir.y_ * r.y_;
//             double nu_dot_k = nu.dot(k_dir);
//             std::complex<double> psi = std::exp(std::complex<double>(0,k.real()*dot_));
//             std::complex<double> phi = std::complex<double>(0,1)*k*nu_dot_k*psi;
//             std::complex<double> dpsi = sol.psi[l]-psi;
//             std::complex<double> dphi = -sol.phi[l]-phi;//注意将法向量方向校准到一致

//             auto I = scattering_kernel(kr, kf_dir, dpsi, dphi, r, nu);
//             // if(abs(I)>100){I=std::complex<double>(0,0);}
//             // sum += I*(elem.seg_ptr->tangent(mid_t).norm())*(elem.seg_end_t - elem.seg_start_t);
//             sum += gp.w *I*(seg.tangent(t).norm())*(t1 - t0);
//             // sum+=elem.length;
//         }
//     }

//     // 最终乘以前置因子
//     return prefac * sum;
//     // return sum;
    
// }

/**
 * compute_total_cross_section
 * ---------------------------
 * 按光学定理 (20)：σ = 2√(π/k) Im[(1−i) f_scat_forward]
 *
 * @param f_scat_forward  θ=φ（入射方向）的散射振幅 f_scat
 * @param k               波数（仅用 real 部）
 * @return                总散射截面 σ
 */
inline double compute_total_cross_section(
    std::complex<double> f_scat_forward,
    std::complex<double> k
) {
    double kr = k.real();
    // (1−i) f_scat_forward
    std::complex<double> expr = (std::complex<double>(1.0,-1.0))
                              * f_scat_forward;
    return 2.0 * std::sqrt(M_PI/kr) * std::imag(expr);
}


// /**
//  * compute_scattering_amplitude（Gauss 版本）
//  * ------------------------------------------
//  * 返回散射部分 f_scat(θ)
//  *
//  * @param sol           BEM 求得的散射解（φ, ψ 向量）
//  * @param elements      ∂Ω_J 上所有边界元（BEMElement 列表）
//  * @param phi_psi_in    入射场 [φ_in; ψ_in] 向量
//  * @param k             波数（复数，但这里只用 real 部）
//  * @param thetas        观测角度列表 θ
//  * @param incident_dir  入射方向单位矢量
//  * @return              对应每个 θ 的散射振幅 f_scat(θ)
//  */
//  inline std::vector<std::complex<double>> compute_scattering_amplitude(
//     const ScatteringResult& sol,
//     const std::vector<Microcavity::BEMElement>& elements,
//     const Eigen::VectorXcd& phi_psi_in,
//     std::complex<double> k,
//     const std::vector<double>& thetas,
//     const Point& incident_dir
//  ) {
//     int N = int(elements.size());
//     Eigen::VectorXcd phi_in = phi_psi_in.head(N);
//     Eigen::VectorXcd psi_in = phi_psi_in.tail(N);
 
//     std::vector<std::complex<double>> f_scat;
//     f_scat.reserve(thetas.size());
 
//     // 公式 (19) 中的前置因子 i/(4√(π k))
//     double kr = k.real();
//     std::complex<double> prefac = std::complex<double>(1,1)
//                                / (4.0 * std::sqrt(M_PI * kr));
 
//     // 观测方向单位矢量函数
//     auto kf = [&](double θ){
//       return Point{std::cos(θ), std::sin(θ)};
//     };
 
//     for(double θ : thetas) {
//       Point kf_dir = kf(θ);
//       std::complex<double> sum = 0.0;
 
//       // 对 ∂Ω_J 做高斯积分
//       for(int ℓ = 0; ℓ < N; ++ℓ) {
//         auto const& elem = elements[ℓ];
//         auto Δψ = sol.psi[ℓ] - psi_in[ℓ];
//         auto Δφ = sol.phi[ℓ] - phi_in[ℓ];
 
//         // 这里传入 k.real()
//         auto I = GaussIntegrator::integrate(
//           elem,
//           elem.mid,
//           [&](auto /*t*/, Point const& s, Point const& raw_n,
//               Point const& /*dsdt*/, Point const& /*s_i*/, double /*k_unused*/) {
//             Point ν = raw_n.normalized();
//             double phase_arg = kr * (kf_dir.dot(s));
//             auto phase = std::exp(std::complex<double>(0, -phase_arg));
//             auto bracket =
//               std::complex<double>(0,1)*(kf_dir.dot(ν))*Δψ
//               + Δφ;
//             return phase * bracket;
//           },
//           /*k=*/ kr
//         );
 
//         sum += I;
//       }
 
//       f_scat.push_back(prefac * sum);
//     }
 
//     return f_scat;
//  }
 
//  /**
//   * compute_total_cross_section
//   * 再代入光学定理 σ = (2√π / k) Im[(1−i) f_tot ]
//   */
//  inline double compute_total_cross_section(
//     std::complex<double> f_scat_forward,
//     std::complex<double> k
//  ) {
//     double kr = k.real();
//     auto expr  = (std::complex<double>(1.0,0.0)
//                 - std::complex<double>(0,1)) * f_scat_forward;
//     return 2.0 * std::sqrt(M_PI*kr) * std::imag(expr);
//  }


// /**
//  * compute_scattering_amplitude（Gauss 版本）
//  * ------------------------------------------
//  * 用 Gauss–Legendre 数值积分精确计算公式 (19)。
//  *
//  * @param sol           BEM 求得的散射解（φ, ψ 向量）
//  * @param elements      ∂Ω_J 上所有边界元（BEMElement 列表）
//  * @param phi_psi_in    入射场 [φ_in; ψ_in] 向量
//  * @param k             波数（复数）
//  * @param thetas        观测角度列表 θ
//  * @param incident_dir  入射方向单位矢量（仅用于 φ_in, ψ_in 已在 phi_psi_in 中）
//  * @return              对应每个 θ 的散射振幅 f(θ,k)
//  */
//  inline std::vector<std::complex<double>> compute_scattering_amplitude(
//     const ScatteringResult& sol,
//     const std::vector<Microcavity::BEMElement>& elements,
//     const Eigen::VectorXcd& phi_psi_in,
//     std::complex<double> k,
//     const std::vector<double>& thetas,
//     const Point& incident_dir = {1.0, 0.0}
// ) {
//     int N = int(elements.size());
//     // 拆分出 φ_in 和 ψ_in
//     Eigen::VectorXcd phi_in = phi_psi_in.head(N);
//     Eigen::VectorXcd psi_in = phi_psi_in.tail(N);

//     std::vector<std::complex<double>> f_vals;
//     f_vals.reserve(thetas.size());

//     // 前置常数 (1+i)/(4 sqrt(pi k))
//     const std::complex<double> prefac =
//       (std::complex<double>(1,1)) / (4.0 * std::sqrt(M_PI * k));

//     // 对每个散射角度 θ 做积分
//     for(double theta : thetas) {
//         // 出射方向波矢单位向量 k_f
//         Point kf_dir{ std::cos(theta), std::sin(theta) };

//         // 积分结果累加
//         std::complex<double> sum = 0.0;

//         // 论文 (19) 中 ∮_∂Ω_J ds { … }
//         // 将 ∂Ω_J 划分为 N 个 BEMElement，逐元素高斯积分
//         for(int l=0; l<N; ++l) {
//             const auto& elem = elements[l];

//             // Δψ 和 Δφ 在整个元素上常数（BEM 近似）
//             std::complex<double> Δψ = sol.psi[l] - psi_in[l];
//             std::complex<double> Δφ = sol.phi[l] - phi_in[l];

//             // 用 GaussIntegrator 在该元素上做积分
//             // kernel 返回被积函数值： phase * { i (k_f·ν) Δψ + Δφ }
//             auto integral = GaussIntegrator::integrate(
//                 elem,
//                 /*s_i=*/elem.mid,  // 不影响 kernel 中 Δψ/Δφ
//                 [&](auto /*t*/, const Point& s, const Point& raw_n,
//                     const Point& /*dsdt*/, const Point& /*s_i*/, std::complex<double> /*k_unused*/) {
//                   // 1) 计算单位外法线 ν(s)
//                   Point nu = raw_n.normalized();
//                   // 2) 相位因子 e^{-i k_f·r(s)}
//                   double phase_arg = k.real() * kf_dir.dot(s);
//                   std::complex<double> phase = std::exp(std::complex<double>(0, -phase_arg));
//                   // 3) 被积函数 { … }
//                   std::complex<double> bracket =
//                     std::complex<double>(0,1) * (kf_dir.dot(nu)) * Δψ
//                     + Δφ;
//                   return phase * bracket;
//                 },
//                 /*k=*/k
//             );

//             sum += integral;
//         }

//         // 最终振幅
//         f_vals.push_back(prefac * sum);
//     }

//     return f_vals;
// }


// /**
//  * compute_differential_cross_section
//  *   dσ/dθ = |f(θ,k)|²
//  */
// inline std::vector<double> compute_differential_cross_section(
//     const std::vector<std::complex<double>>& f_vals
// ) {
//     std::vector<double> ds;
//     ds.reserve(f_vals.size());
//     for (auto& f : f_vals) {
//         ds.push_back(std::norm(f));
//     }
//     return ds;
// }

// /**
//  * compute_total_cross_section
//  * 光学定理：σ(k) = (2√π / k) Im[(1−i) f_forward]
//  */
// inline double compute_total_cross_section(
//     std::complex<double> f_forward,
//     std::complex<double> k
// ) {
//     std::complex<double> expr = (1.0 - std::complex<double>(0,1)) * f_forward;
//     return 2.0 * std::sqrt(M_PI) / k.real() * std::imag(expr);
// }

} // namespace bem

#endif // BEM_SCATTERING_H
