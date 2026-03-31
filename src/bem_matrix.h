#ifndef BEM_MATRIX_H
#define BEM_MATRIX_H
#include "bem_gauss_integrator.h"
#include "../include/eigen-3.4.0/Eigen/Dense"
#include "../include/eigen-3.4.0/Eigen/Sparse"
#include "bem_microcavity.h"
#include <complex>
#include <fstream>
#include <mutex>

namespace bem {

// 定义偏振类型的枚举
enum Polarization {
    TM,
    TE
};

/**
 * @brief 构建 B_j 块 (边界积分方程中 ψ 部分)
 *
 * 对应论文中 B(s',s) = -2 G(s',s) ，
 * G = (i/4) H_0^{(1)}(n k |r - r'|)
 */
 inline
 Eigen::MatrixXcd build_B_matrix_block(
     const std::vector<Microcavity::BEMElement>& elems,
     Polarization pol,
     double k,
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
                 B(i,l) = compute_B_ll(el, k, /*n_j=*/n_inner);
                 if(out) {
                    B(i,l)=-1.0*B(i,l);
                 }
             } else {
                 // 数值积分
                 std::complex<double> I = GaussIntegrator::integrate(
                     el,
                     ei.mid,
                     GaussIntegrator::kernel_tm_B<CurveSegment>,
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
 
 /**
  * @brief 构建 C_j 块 (边界积分方程中 ∂_ν ψ 部分)
  *
  * 对应论文中 C(s',s) = 2 ∂_ν G(s',s) ，
  * ∂_ν G = (i/4) n k cosα H_1^{(1)}(n k r)
  */
 inline
 Eigen::MatrixXcd build_C_matrix_block(
     const std::vector<Microcavity::BEMElement>& elems,
     Polarization pol,
     double k,
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
                 C(i,l) = std::complex<double>( compute_C_ll(el), 0.0 );
                 if (out)  {
                    C(i,l)=-1.0*C(i,l)-2.0;
                 }
             } else {
                 std::complex<double> I = GaussIntegrator::integrate(
                     el,
                     ei.mid,
                     GaussIntegrator::kernel_tm_C<CurveSegment>,
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


/**
 * @brief 构建稠密矩阵 M 和 M0
 *
 * 参照论文式 (17):
 *   M [φ;ψ] = M0 [φ_in; ψ_in]
 *
 * - 对前 J-1 个腔体：B_j,C_j 放在对角块  
 * - 最后外部区域：B_ext,C_ext 放在最下方，且同样填入 M0
 */
 inline
 std::pair<Eigen::MatrixXcd, Eigen::MatrixXcd>
 build_matrix(
     const std::vector<Microcavity::BEMElement>& all_elems,
     const std::vector<std::vector<Microcavity::BEMElement>>& per_cav,
     Polarization pol,
     double k,
     double n_inner = 1.466,
     double n_outer = 1.0
 ) {
     // 所有边界元的总数
     size_t N_el = all_elems.size();
     size_t N    = 2 * N_el;
 
     Eigen::MatrixXcd M  = Eigen::MatrixXcd::Zero(N, N);
     Eigen::MatrixXcd M0 = Eigen::MatrixXcd::Zero(N, N);
 
     // 1) J-1 个腔体块
     size_t row = 0;
     for (auto const& elems : per_cav) {
         int n_j = static_cast<int>(elems.size());
         // 构建 B_j, C_j
         auto B_j = build_B_matrix_block(elems, pol, k, n_inner);
         auto C_j = build_C_matrix_block(elems, pol, k, n_inner);
 
         // 插入 B_j
         M.block(row, row, n_j, n_j) = B_j;
         // 插入 C_j
         M.block(row, N_el + row, n_j, n_j) = C_j;
 
         row += n_j;
     }
 
     // 2) 外部区域块
     //    注意折射率互换: 内→外, 外→内
     auto B_ext = build_B_matrix_block(all_elems, pol, k, n_outer, true);
     auto C_ext = build_C_matrix_block(all_elems, pol, k, n_outer, true);
 
     // 将它们插到底部，同时也写入 M0
     M .block(row,           0, N_el, N_el) = B_ext;
     M .block(row, N_el, N_el, N_el) = C_ext;
     M0.block(row,           0, N_el, N_el) = B_ext;
     M0.block(row, N_el, N_el, N_el) = C_ext;
 
     return {M, M0};
 }
 
inline Eigen::VectorXcd compute_incident_vectors(
    const std::vector<Microcavity::BEMElement>& all_discretized_elements,
    double k,
    const Point& incident_dir
) {
    size_t N = all_discretized_elements.size();
    Eigen::VectorXcd psi_in(N), phi_in(N);

    // 单位入射方向
    Point k_dir = incident_dir.normalized();

    // std::ofstream ofs("./data/incident_vectors.csv", std::ios::out);
    // ofs << "idx,phi_re,phi_im,psi_re,psi_im\n";

    // 对每个元素，用 Simpson 公式（三点二次）来近似 ψ_in, φ_in
    for (size_t i = 0; i < N; ++i) {
        const auto& elem = all_discretized_elements[i];
        auto seg = elem.seg_ptr;

        // 参数 t0, tm, t1
        double t0 = elem.seg_start_t;
        double t1 = elem.seg_end_t;
        double tm = 0.5 * (t0 + t1);

        // 取样三点
        std::array<double,3> ts = { t0, tm, t1 };
        std::array<Point,3>   rs;
        std::array<Point,3>   nus;
        for (int j = 0; j < 3; ++j) {
            double t = ts[j];
            rs[j]  = (*seg)(t);
            // 对所有元素，用外法线方向
            nus[j] = seg->normal(t).normalized();
        }

        // 计算各点的 ψ 和 φ
        std::array<std::complex<double>,3> psis, phis;
        for (int j = 0; j < 3; ++j) {
            auto const& r  = rs[j];
            auto const& nu = nus[j];
            double dot   = k_dir.x_ * r.x_ + k_dir.y_ * r.y_;
            // ψ = e^{i k (k_dir·r)}
            psis[j] = std::exp(std::complex<double>(0, k * dot));
            // φ = ∂_ν ψ = i k (ν·k_dir) ψ
            double nu_dot_k = nu.dot(k_dir);
            phis[j] = std::complex<double>(0,1) * k * nu_dot_k * psis[j];
        }

        // Simpson 加权平均: ∫₀¹ f(t) dt ≈ (f0 + 4 f_m + f1)/6
        std::complex<double> psi_avg = (psis[0] + std::complex<double>(4.0) * psis[1] + psis[2]) / 6.0;
        std::complex<double> phi_avg = (phis[0] + std::complex<double>(4.0) * phis[1] + phis[2]) / 6.0;

        // 存储
        psi_in[i] = psi_avg;
        phi_in[i] = phi_avg;

        // // 写入 CSV
        // ofs << i << ","
        //     << phi_avg.real() << "," << phi_avg.imag() << ","
        //     << psi_avg.real() << "," << psi_avg.imag() << "\n";
    }

    // ofs.close();

    Eigen::VectorXcd result(2*N);
    result << phi_in, psi_in;
    return result;
}
// inline Eigen::VectorXcd compute_incident_vectors(
//     const std::vector<Microcavity::BEMElement>& all_discretized_elements,
//     double k,
//     const Point& incident_dir
// ) {
//     size_t N = all_discretized_elements.size();
//     Eigen::VectorXcd psi_in(N), phi_in(N);

//     // 单位入射方向
//     Point k_dir = incident_dir.normalized();

//     // 对每个元素，用 Boole 公式（五点四次）来近似 ψ_in, φ_in
//     for (size_t i = 0; i < N; ++i) {
//         const auto& elem = all_discretized_elements[i];
//         auto seg = elem.seg_ptr;

//         // 参数 t0, t1
//         double t0 = elem.seg_start_t;
//         double t1 = elem.seg_end_t;
//         double h = (t1 - t0) / 4.0;

//         // 取样五点
//         std::array<double, 5> ts = { t0, t0 + h, t0 + 2 * h, t0 + 3 * h, t1 };
//         std::array<Point, 5> rs;
//         std::array<Point, 5> nus;
//         for (int j = 0; j < 5; ++j) {
//             double t = ts[j];
//             rs[j] = (*seg)(t);
//             // 对所有元素，用外法线方向
//             nus[j] = seg->normal(t).normalized();
//         }

//         // 计算各点的 ψ 和 φ
//         std::array<std::complex<double>, 5> psis, phis;
//         for (int j = 0; j < 5; ++j) {
//             auto const& r = rs[j];
//             auto const& nu = nus[j];
//             double dot = k_dir.x_ * r.x_ + k_dir.y_ * r.y_;
//             // ψ = e^{i k (k_dir·r)}
//             psis[j] = std::exp(std::complex<double>(0, k * dot));
//             // φ = ∂_ν ψ = i k (ν·k_dir) ψ
//             double nu_dot_k = nu.dot(k_dir);
//             phis[j] = std::complex<double>(0, 1) * k * nu_dot_k * psis[j];
//         }

//         // Boole 加权平均: ∫₀¹ f(t) dt ≈ (7f0 + 32f1 + 12f2 + 32f3 + 7f4) / 90
//         std::complex<double> psi_avg = (7.0 * psis[0] + 32.0 * psis[1] + 12.0 * psis[2] + 32.0 * psis[3] + 7.0 * psis[4]) / 90.0;
//         std::complex<double> phi_avg = (7.0 * phis[0] + 32.0 * phis[1] + 12.0 * phis[2] + 32.0 * phis[3] + 7.0 * phis[4]) / 90.0;

//         psi_in[i] = psi_avg;
//         phi_in[i] = phi_avg;
//     }

//     Eigen::VectorXcd result(2 * N);
//     result << phi_in, psi_in;
//     return result;
// }

// inline Eigen::VectorXcd compute_incident_vectors(
//     const std::vector<Microcavity::BEMElement>& all_discretized_elements,
//     double k,
//     const Point& incident_dir
// ) {
//     size_t N = all_discretized_elements.size();
//     Eigen::VectorXcd psi_in(N), phi_in(N);

//     // 单位入射方向
//     Point k_dir = incident_dir.normalized();

//     std::ofstream ofs("./data/incident_vectors.csv", std::ios::out);
//     ofs << "idx,phi_re,phi_im,psi_re,psi_im\n";

//     for (size_t i = 0; i < N; ++i) {
//         const auto& elem = all_discretized_elements[i];
//         // r = 元素中点坐标
//         Point r = elem.mid;
//         // 单位外法线
//         double mid_t = 0.5 * (elem.seg_start_t + elem.seg_end_t);
//         Point nu = 1.0*elem.seg_ptr->normal(mid_t);

//         // 1) 计算相位：phase = k·(k_dir · r)
//         double dot = k_dir.x_ * r.x_ + k_dir.y_ * r.y_;
//         // 2) 构造虚指数 i k dot
//         std::complex<double> ikr = std::complex<double>(0,1) * k * dot;
//         // 3) 入射波 ψ_in = e^{i k·r}
//         psi_in[i] = std::exp(ikr);

//         // 4) 法向导数 φ_in = ∂_ν ψ_in = i k (ν·k_dir) ψ_in
//         double nu_dot_k = nu.dot(k_dir);
//         phi_in[i] = std::complex<double>(0,1) * k * nu_dot_k * psi_in[i];

//                 // 写一行 CSV
//                 ofs << i << ","
//                 << std::fixed << std::setprecision(6)
//                     << phi_in[i].real() << "," << phi_in[i].imag() << ","
//                     << psi_in[i].real() << "," << psi_in[i].imag()
//                 << "\n";
//     }

//     ofs.close();

//     Eigen::VectorXcd result(2*N);
//     result << phi_in, psi_in;
//     return result;
// }

inline Eigen::VectorXcd build_rhs(    
    const std::vector<Microcavity::BEMElement>& all_discretized_elements,
    const double& k,
    const double& index,
    const Point& incident_dir){
        // auto phi_psi_in  = compute_incident_vectors(all_discretized_elements, k, incident_dir);
        size_t N = all_discretized_elements.size();
        Eigen::VectorXcd zero(N), rhs(N);
            // 打开输出文件
        // std::ofstream diff_ofs("./data/dsum_diff.txt");
        // if (!diff_ofs) {
        //     throw std::runtime_error("Cannot open dsum_diff.txt for writing");
        // }
        for (int i=0;i<N;i++) {
            std::complex<double> sum=0;
            for(int l=0;l<N;l++) {
                if(i==l) {
                sum+=compute_B_ll_out_phi_l(all_discretized_elements[i], k, incident_dir)
                +compute_C_ll_out_psi_l(all_discretized_elements[i], k, incident_dir);
                
                // if(i==0){
                // auto dsum=(compute_B_ll_out_phi_l(all_discretized_elements[i], k, incident_dir)
                // +compute_C_ll_out_psi_l(all_discretized_elements[i], k, incident_dir));
                // auto dsum2=(-compute_B_ll(all_discretized_elements[i], k)*phi_psi_in[l]
                // +(-compute_C_ll(all_discretized_elements[i])-2)*phi_psi_in[l+N]);
                // auto diff = dsum - dsum2;
                // diff_ofs << l << " " << diff.real() << " " << diff.imag() << "\n";
                
                // auto dsum=(compute_B_ll_out_phi_l(all_discretized_elements[i], k, incident_dir));
                // auto dsum2=(-compute_B_ll(all_discretized_elements[i], k)*phi_psi_in[l]);
                // auto diff = dsum-dsum2;
                // diff_ofs << l << " " << diff.real() << " " << diff.imag() << "\n";
                
                // auto dsum=(compute_C_ll_out_psi_l(all_discretized_elements[i], k, incident_dir));
                // auto dsum2=(-compute_C_ll(all_discretized_elements[i])-2)*phi_psi_in[l+N];
                // auto diff = dsum-dsum2;
                // diff_ofs << l << " " << diff.real() << " " << diff.imag() << "\n";
                // }
                    
                }
                else {
                    sum+=GaussIntegrator::integrate(
                        all_discretized_elements[l],
                        all_discretized_elements[i].mid,
                        GaussIntegrator::kernel_tm_X_il_out_x_l<CurveSegment>,
                        k,
                        index,
                        incident_dir
                    );

                    
                    // if(i==0){
                    //     auto dsum=GaussIntegrator::integrate(
                    //         all_discretized_elements[l],
                    //         all_discretized_elements[i].mid,
                    //         GaussIntegrator::kernel_tm_X_il_out_x_l<CurveSegment>,
                    //         k,
                    //         index,
                    //         incident_dir
                    //     );
                    //     auto dsum2=-GaussIntegrator::integrate(
                    //         all_discretized_elements[l],
                    //         all_discretized_elements[i].mid,
                    //         GaussIntegrator::kernel_tm_B<CurveSegment>,
                    //         k,
                    //         index,
                    //         incident_dir
                    //     )*phi_psi_in[l]-GaussIntegrator::integrate(
                    //         all_discretized_elements[l],
                    //         all_discretized_elements[i].mid,
                    //         GaussIntegrator::kernel_tm_C<CurveSegment>,
                    //         k,
                    //         index,
                    //         incident_dir
                    //     )*phi_psi_in[l+N];
                    //     auto diff = dsum - dsum2;
                    // diff_ofs << l << " " << diff.real() << " " << diff.imag() << "\n";
                    // }
                }
            }
            
            
            zero(i)=0;
            rhs(i)=sum;
        }
    // diff_ofs.close();
    Eigen::VectorXcd result(2*N);
    result << zero, rhs;
    return result;
}


} // namespace bem
#endif