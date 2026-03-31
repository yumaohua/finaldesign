// main.cpp

#include <boost/math/special_functions/math_fwd.hpp>
#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <iomanip>
#include <mutex>
#include <omp.h>
#include <algorithm>
#include <map>
#include <fstream>
// #include "bem_gauss_integrator.h"
// #include "bem_hankel_manual.h"
#include "bem_curve_segment.h"
#include "bem_hankel_manual.h"
#include "bem_microcavity.h"
#include "bem_matrix.h"
#include "bem_scattering.h"
#include "bem_output.h"
#include "bem_resonance.h"

using namespace bem;
// Filename: compute_resonant_mode.cpp

#include <tuple>

#include "bem_gauss_integrator.h"   // complex_integrate(), kernels
#include "../include/eigen-3.4.0/Eigen/Dense"
#include "../include/eigen-3.4.0/Eigen/SVD"

using namespace bem;

int main(){
    // --- 1) 复波数 k_res
    std::complex<double> k_res( 2.294444 , -0.009696);

    // --- 2) 构造两个微腔并离散化（你原来调用的方式）
    double R = 10.0, corner_radius = 0.06, n_inside = 1.466, n_outer = 1.0;
    Microcavity cav1 = make_hermite_rounded_hexagonal_cavity(R, corner_radius, n_inside);
    Microcavity cav2 = make_hermite_rounded_hexagonal_cavity(R, corner_radius, n_inside);
    cav1.move({-0.9*R, 0.25*R});
    cav2.move({ 0.9*R,-0.25*R});
    auto elems1 = cav1.discretize(240, 24, 10/240.0,18);
    auto elems2 = cav2.discretize(240, 24, 10/240.0,18);

    // 合并
    std::vector<Microcavity::BEMElement> all_elems;
    all_elems.reserve(elems1.size()+elems2.size());
    all_elems.insert(all_elems.end(), elems1.begin(), elems1.end());
    all_elems.insert(all_elems.end(), elems2.begin(), elems2.end());
    std::vector<std::vector<Microcavity::BEMElement>> per_cav{elems1, elems2};

    int N = (int)all_elems.size();
    int NN = 2*N;

    // --- 3) 构建 M(k_res)
    Eigen::MatrixXcd M = complex_build_matrix(all_elems, per_cav, TM, k_res, n_inside, n_outer);

    // --- 4) SVD 求零特征向量
    Eigen::JacobiSVD<Eigen::MatrixXcd> svd(M, Eigen::ComputeFullV);
    Eigen::VectorXcd vec = svd.matrixV().col(NN-1);  // 最后一列
    Eigen::VectorXcd phi = vec.head(N);
    Eigen::VectorXcd psi = vec.tail(N);

    // --- 5) 输出边界元中心点
    std::ofstream felems("./dataf/elements_midpoints.txt");
    for(auto &e : all_elems){
        felems << e.mid.x_ << "," << e.mid.y_ << "\n";
    }
    felems.close();

    // --- 6) 在区域上做灰度网格
    std::ofstream fgrid("./dataf/psi2_grid.txt");
    int nx = 400, ny = 300;
    double x0=-20, x1=20, y0=-15, y1=15;
    for(int ix=0; ix<nx; ++ix){
        double x = x0 + (x1-x0)*ix/(nx-1);
        for(int iy=0; iy<ny; ++iy){
            double y = y0 + (y1-y0)*iy/(ny-1);
            // 对每个边界元 l，计算两个积分
            std::complex<double> sum1=0, sum2=0;
            Point r{ x,y };
            for(int l=0; l<N; ++l){
                auto &el = all_elems[l];
                // ∫_l ds ∂_ν G(s,r)
                auto I1 = GaussIntegrator::complex_integrate(
                    el, r,
                    GaussIntegrator::kernal_partial_G<CurveSegment>,
                    k_res, n_outer
                );
                // ∫_l ds G(s,r)
                auto I2 = GaussIntegrator::complex_integrate(
                    el, r,
                    GaussIntegrator::kernel_G<CurveSegment>,
                    k_res, n_outer
                );
                sum1 += -psi[l] * I1;
                sum2 += -phi[l] * I2;
            }
            std::complex<double> psi_r = sum1 - sum2;
            double psi2 = std::norm(psi_r);
            fgrid << x << "," << y << "," << psi2 << "\n";
        }
    }
    fgrid.close();

    std::cout<<"All data written.\n";
    return 0;
}

// int main() {
//     std::ofstream outFile("./dataf/error_data.txt");
//     if (!outFile.is_open()) {
//         std::cerr << "无法打开文件" << std::endl;
//         return 1;
//     }

//     outFile << "z,error_nu_0,error_nu_1" << std::endl;

//     for (double z_val = 0.0001; z_val <= 100; z_val += 0.0001) {
//         double z=z_val;

//         // 计算阶数为 0 的误差
//         std::complex<double> hankel_1_0_custom = complex_hankel_1(0, z);
//         std::complex<double> hankel_1_0_boost = boost::math::cyl_hankel_1(0, z);
//         double error_0 = std::abs(hankel_1_0_custom - hankel_1_0_boost);

//         // 计算阶数为 1 的误差
//         std::complex<double> hankel_1_1_custom = complex_hankel_1(1, z);
//         std::complex<double> hankel_1_1_boost = boost::math::cyl_hankel_1(1, z);
//         double error_1 = std::abs(hankel_1_1_custom - hankel_1_1_boost);

//         outFile << z_val << "," << error_0 << "," << error_1 << std::endl;
//     }

//     outFile.close();
//     std::cout << "误差数据已输出到 error_data.txt" << std::endl;

//     return 0;
// }    

// int main() {
//     // —— 1) 全局参数 ——
//     // double input=1;
//     // while(input){
//     //     std::cin>>input;
//     //     std::cout<<"0: "<<complex_hankel_1(0, input)<<", "<<cyl_hankel_1(0,input)<<std::endl;
//     //     std::cout<<"1: "<<complex_hankel_1(1, input)<<", "<<cyl_hankel_1(1,input)<<std::endl;
//     // }
//     double R             = 10.0;           // 六边形边长
//     double corner_radius = 0.06;       // 圆角半径
//     // double corner_radius = 0.0044*R;       // 圆角半径
//     double n_inside      = 1.466;          // 腔内折射率
//     double n_outside     = 1.0;            // 腔外折射率
//     Polarization pol     = TM;             // TM 偏振
//     double phi_inc_deg   = 15.0;           // 入射角（度）
//     Point  incident_dir(std::cos(phi_inc_deg*M_PI/180.0),
//                         std::sin(phi_inc_deg*M_PI/180.0));

//         // —— 2) 构造两个圆角六边形腔体 ——
//     std::cout << "[Info] Constructing two rounded hexagonal cavities...\n";
//     Microcavity cav1 = make_hermite_rounded_hexagonal_cavity(R, corner_radius, n_inside);
//     // Microcavity cav1 = make_bezier_rounded_hexagonal_cavity(R, corner_radius, n_inside);
//     // Microcavity cav1 = make_rounded_hexagonal_cavity(R, corner_radius, n_inside);
//     Microcavity cav2 = make_hermite_rounded_hexagonal_cavity(R, corner_radius, n_inside);
//     // Microcavity cav2 = make_bezier_rounded_hexagonal_cavity(R, corner_radius, n_inside);
//     // Microcavity cav2 = make_rounded_hexagonal_cavity(R, corner_radius, n_inside);

//         // —— 4) 平移两腔体到指定中心 —— 
//     Point offset1(-0.9*R,  0.25*R);
//     Point offset2( 0.9*R, -0.25*R);

//     cav1.move(offset1);
//     cav2.move(offset2);

//         // auto elems1 = cav1.discretize(249, 18);
//     auto elems1 = cav1.discretize(240, 24, 10/240.0,18);
//     // auto elems2 = cav2.discretize(249, 18);
//     auto elems2 = cav2.discretize(240, 24,10/240.0,18);
//     std::cout << "[Info]  Cav1 elements: " << elems1.size()
//               << ", Cav2 elements: " << elems2.size() << "\n";

//             // 合并所有元素及分区信息
//     std::vector<Microcavity::BEMElement> all_elems;
//     all_elems.reserve(elems1.size() + elems2.size());
//     all_elems.insert(all_elems.end(), elems1.begin(), elems1.end());
//     all_elems.insert(all_elems.end(), elems2.begin(), elems2.end());
//     std::vector<std::vector<Microcavity::BEMElement>> per_cav = { elems1, elems2 };
//     std::cout<<"begin find resonance-----"<<std::endl;
//     auto res = find_resonance(all_elems, per_cav, pol, std::complex<double>(2.295,-0.0098));
//     // auto res = find_resonance(all_elems, per_cav, pol, std::complex<double>(22.95,-0.098));
//     std::cout<<res.k_res<<std::endl;
//     return 0;
// }

// int main(){
//     // 并行线程
//     omp_set_num_threads(2);
//     std::mutex cout_mutex, file_mutex;

//     // —— 1) 全局参数 ——
//     double R             = 10.0;           // 六边形边长
//     double corner_radius = 0.06;       // 圆角半径
//     // double corner_radius = 0.0044*R;       // 圆角半径
//     double n_inside      = 1.466;          // 腔内折射率
//     double n_outside     = 1.0;            // 腔外折射率
//     Polarization pol     = TM;             // TM 偏振
//     double phi_inc_deg   = 15.0;           // 入射角（度）
//     Point  incident_dir(std::cos(phi_inc_deg*M_PI/180.0),
//                         std::sin(phi_inc_deg*M_PI/180.0));

//         // —— 2) 构造两个圆角六边形腔体 ——
//     std::cout << "[Info] Constructing two rounded hexagonal cavities...\n";
//     Microcavity cav1 = make_hermite_rounded_hexagonal_cavity(R, corner_radius, n_inside);
//     // Microcavity cav1 = make_bezier_rounded_hexagonal_cavity(R, corner_radius, n_inside);
//     // Microcavity cav1 = make_rounded_hexagonal_cavity(R, corner_radius, n_inside);
//     Microcavity cav2 = make_hermite_rounded_hexagonal_cavity(R, corner_radius, n_inside);
//     // Microcavity cav2 = make_bezier_rounded_hexagonal_cavity(R, corner_radius, n_inside);
//     // Microcavity cav2 = make_rounded_hexagonal_cavity(R, corner_radius, n_inside);

//         // —— 4) 平移两腔体到指定中心 —— 
//     Point offset1(-0.9*R,  0.25*R);
//     Point offset2( 0.9*R, -0.25*R);

//     cav1.move(offset1);
//     cav2.move(offset2);

//         // auto elems1 = cav1.discretize(249, 18);
//     auto elems1 = cav1.discretize(240, 24, 10/240.0,18);
//     // auto elems2 = cav2.discretize(249, 18);
//     auto elems2 = cav2.discretize(240, 24,10/240.0,18);
//     std::cout << "[Info]  Cav1 elements: " << elems1.size()
//               << ", Cav2 elements: " << elems2.size() << "\n";
//     // 导出中点
//     Point off1(0,  0),
//         off2( 0, 0);
//     save_midpoints_to_txt(elems1, off1, elems2, off2, "./data/cavities_midpoints.txt");
//     std::cout << "[Info] midpoints written to cavities_midpoints.txt\n";

//         // 合并所有元素及分区信息
//     std::vector<Microcavity::BEMElement> all_elems;
//     all_elems.reserve(elems1.size() + elems2.size());
//     all_elems.insert(all_elems.end(), elems1.begin(), elems1.end());
//     all_elems.insert(all_elems.end(), elems2.begin(), elems2.end());
//     std::vector<std::vector<Microcavity::BEMElement>> per_cav = { elems1, elems2 };

//     int    N_k    = 101;
//     double k_min  = 2.0, k_max = 2.5;
//     std::vector<double> k_vals(N_k), kR_vals(N_k);
//     for(int i=0;i<N_k;++i){
//         k_vals[i]  = k_min + (k_max - k_min)*i/(N_k-1);
//         kR_vals[i] = k_vals[i]*R;
//     }


//     // 3) 读已有结果，初始化 done[] 和 sigma_norm[]
//     const std::string cs_file = "./data/cross_section.txt";
//     std::map<double,double> existing;
//     {
//         std::ifstream ifs(cs_file);
//         double kR, sigma;
//         while(ifs >> kR >> sigma){
//             existing[kR] = sigma;
//         }
//     }
//     int begin = std::max((int)existing.size()-4,0);
//     std::vector<double> sigma_norm(N_k, NAN);
//     std::vector<bool>   done     (N_k, false);
//     for(int i=0; i<N_k; ++i){
//         auto it = existing.find(kR_vals[i]);
//         if(it != existing.end()){
//             sigma_norm[i] = it->second;
//             done[i] = true;

//         }
//     }

//     // 4) 并行扫描
//     #pragma omp parallel for schedule(dynamic)
//     for(int i=begin;i<N_k;++i){
//         double k  = k_vals[i];
//         double kR = kR_vals[i];
//         {
//             std::lock_guard<std::mutex> lk(cout_mutex);
//             std::cout << "[Step " << (i+1)<<"/"<<N_k<<"] k="<<k<<" ... ";
//         }

//         // 如果这一条已经在文件里了，就直接跳过
//         if(done[i]){
//             std::lock_guard<std::mutex> lk(cout_mutex);
//             std::cout << "skipped (σ/R="<<sigma_norm[i]<<")\n";
//             continue;
//         }

//         // —— 否则正常跑：compute_incident_vectors, build_rhs, solve_scattering ……
//         auto phi_psi_in = compute_incident_vectors(all_elems, k, incident_dir);
//         auto rhs        = build_rhs(all_elems, k, n_outside, incident_dir);
//         auto sol        = solve_scattering(rhs, all_elems, per_cav,
//                                            pol, k, incident_dir,
//                                            n_inside, n_outside);

//         double phi_rad = phi_inc_deg * M_PI/180.0;
//         auto f_scat = compute_scattering_amplitude(
//             sol, all_elems, phi_psi_in,
//             std::complex<double>(k,0.0), phi_rad,
//             incident_dir, false
//         );
//         double sigma = compute_total_cross_section(f_scat, k);
//         sigma_norm[i] = sigma / R;

//         // —— 将新结果追加到文件末尾 —— 
//         {
//             std::lock_guard<std::mutex> lk(file_mutex);
//             std::ofstream ofs(cs_file, std::ofstream::app);
//             ofs << std::setprecision(8)
//                 << kR << " " << sigma_norm[i] << "\n";
//         }

//         {
//             std::lock_guard<std::mutex> lk(cout_mutex);
//             std::cout << "done. σ/R="<<sigma_norm[i]<<"\n";
//         }
//     }

//     std::cout << "[Done] results in " << cs_file << "\n";
//     return 0;
// }

// int main() {
//     // —— 并行线程数设置 ——
//     omp_set_num_threads(2);

//     std::cout << "[Info] Setting up parameters...\n";
//     // —— 1) 全局参数 ——
//     double R             = 10.0;           // 六边形边长
//     double corner_radius = 0.06;       // 圆角半径
//     // double corner_radius = 0.0044*R;       // 圆角半径
//     double n_inside      = 1.466;          // 腔内折射率
//     double n_outside     = 1.0;            // 腔外折射率
//     Polarization pol     = TM;             // TM 偏振
//     double phi_inc_deg   = 15.0;           // 入射角（度）
//     Point  incident_dir(std::cos(phi_inc_deg*M_PI/180.0),
//                         std::sin(phi_inc_deg*M_PI/180.0));

//     // —— 2) 构造两个圆角六边形腔体 ——
//     std::cout << "[Info] Constructing two rounded hexagonal cavities...\n";
//     Microcavity cav1 = make_hermite_rounded_hexagonal_cavity(R, corner_radius, n_inside);
//     // Microcavity cav1 = make_bezier_rounded_hexagonal_cavity(R, corner_radius, n_inside);
//     // Microcavity cav1 = make_rounded_hexagonal_cavity(R, corner_radius, n_inside);
//     Microcavity cav2 = make_hermite_rounded_hexagonal_cavity(R, corner_radius, n_inside);
//     // Microcavity cav2 = make_bezier_rounded_hexagonal_cavity(R, corner_radius, n_inside);
//     // Microcavity cav2 = make_rounded_hexagonal_cavity(R, corner_radius, n_inside);

//     // —— 4) 平移两腔体到指定中心 —— 
//     Point offset1(-0.9*R,  0.25*R);
//     Point offset2( 0.9*R, -0.25*R);

//     cav1.move(offset1);
//     cav2.move(offset2);
//     // std::cout<<std::endl;
//     // std::cout<<2<<std::endl;
//     // std::cout<<std::endl;
//     // auto elems1 = cav1.discretize(249, 18);
//     auto elems1 = cav1.discretize(240, 24, 10/240.0,18);
//     // auto elems2 = cav2.discretize(249, 18);
//     auto elems2 = cav2.discretize(240, 24,10/240.0,18);
//     std::cout << "[Info]  Cav1 elements: " << elems1.size()
//               << ", Cav2 elements: " << elems2.size() << "\n";
//     // 导出中点
//     Point off1(0,  0),
//         off2( 0, 0);
//     save_midpoints_to_txt(elems1, off1, elems2, off2, "./data/cavities_midpoints.txt");
//     std::cout << "[Info] midpoints written to cavities_midpoints.txt\n";



//     // 合并所有元素及分区信息
//     std::vector<Microcavity::BEMElement> all_elems;
//     all_elems.reserve(elems1.size() + elems2.size());
//     all_elems.insert(all_elems.end(), elems1.begin(), elems1.end());
//     all_elems.insert(all_elems.end(), elems2.begin(), elems2.end());
//     std::vector<std::vector<Microcavity::BEMElement>> per_cav = { elems1, elems2 };

//     // —— 5) 扫描 k ∈ [2,2.5] —— 
//     int    N_k    = 6;
//     double k_min  = 2.0, k_max = 2.5;
//     std::vector<double> k_vals(N_k), kR_vals(N_k), sigma_vals(N_k);

//     // 预先填好 k 和 kR=k*R
//     for (int i = 0; i < N_k; ++i) {
//         k_vals[i]  = k_min + (k_max - k_min) * i / (N_k - 1);
//         kR_vals[i] = k_vals[i] * R;
//     }

//     std::cout << "[Info] Starting parallel scan over " << N_k << " k values...\n";
//     std::mutex cout_mutex;
//     #pragma omp parallel for schedule(dynamic)
//     for (int i = 0; i < N_k; ++i) {
//         double k = k_vals[i];
//         {
//             std::lock_guard<std::mutex> lock(cout_mutex);
//             std::cout << "[Step " << (i+1) << "/" << N_k << "] "
//                       << "k=" << std::fixed << std::setprecision(4) << k << " ... " << std::flush;
//         }

//         try {


//             // 5.2) 构造入射场向量
//             auto phi_psi_in = compute_incident_vectors(
//                 all_elems,
//                 k,
//                 incident_dir
//             );
//             auto rhs = build_rhs(all_elems, k, n_outside, incident_dir);
//             // 5.1) 求解散射边界值问题
//             auto sol = solve_scattering( rhs,
//                 all_elems, per_cav,
//                 pol, k, incident_dir,
//                 n_inside, n_outside
//             );

//                     // === 输出 Δφ, Δψ 到 CSV ===
//             {
//                 size_t N_el = all_elems.size();
//                 // 拆分出 φ_in, ψ_in
//                 Eigen::VectorXcd phi_in = phi_psi_in.head(N_el);
//                 Eigen::VectorXcd psi_in = phi_psi_in.tail(N_el);
//                 // 创建文件名，方便区分
//                 std::ostringstream fname;
//                 fname << "./data/delta_" << std::setw(2) << std::setfill('0') << (i+1) << ".csv";
//                 std::ofstream ofs(fname.str());
//                 ofs << "idx,delta_phi_re,delta_phi_im,delta_psi_re,delta_psi_im\n";
//                 for (size_t j = 0; j < N_el; ++j) {
//                     std::complex<double> dphi = sol.phi[j] - phi_in[j];
//                     std::complex<double> dpsi = sol.psi[j] - psi_in[j];
//                     ofs 
//                     << j << ","
//                     << dphi.real() << "," << dphi.imag() << ","
//                     << dpsi.real() << "," << dpsi.imag() << "\n";
//                 }
//             }

//             // 5.3) 计算前向散射振幅 f(θ=φ_inc)
//             double phi_rad = phi_inc_deg * M_PI / 180.0;
//             auto f_scat = compute_scattering_amplitude(
//                 sol, all_elems, phi_psi_in,
//                 std::complex<double>(k,0.0),
//                 phi_rad,
//                 incident_dir,
//                 false
//             );

//             // 5.4) 根据光学定理得总散射截面 σ(k)
//             double sigma = compute_total_cross_section(f_scat, k);

//             // 保存并打印
//             {
//                 std::lock_guard<std::mutex> lock(cout_mutex);
//                 sigma_vals[i] = sigma;
//                 std::cout << "done. |f0|=" << std::scientific << std::setprecision(3)
//                           << std::abs(f_scat)
//                           << ", σ=" << std::fixed << std::setprecision(6) << sigma
//                           << "\n";
//             }
//         } catch (const std::exception& e) {
//             std::lock_guard<std::mutex> lock(cout_mutex);
//             std::cerr << "[Error] step " << i << " failed: " << e.what() << "\n";
//             sigma_vals[i] = NAN;
//         }
//     }

//     // —— 6) 写出归一化的 total cross section —— 
//     //      横坐标 kR，纵坐标 σ/R
//     std::vector<double> sigma_norm(N_k);
//     for (int i = 0; i < N_k; ++i) {
//         sigma_norm[i] = sigma_vals[i] / R;
//     }

//     std::cout << "[Info] Writing cross_section.txt ...\n";
//     write_cross_section(kR_vals, sigma_norm, "./data/cross_section.txt");

//     std::cout << "[Done] All outputs written.\n";
//     return 0;
// }

// int main() {
//     // —— 并行线程数设置 ——
//     omp_set_num_threads(2);

//     std::cout << "[Info] Setting up parameters...\n";
//     // —— 1) 全局参数 ——
//     double R             = 10.0;           // 六边形边长
//     // double corner_radius = 0.0276*R;       // 圆角半径
//     double corner_radius = 0.0044*R;       // 圆角半径
//     double n_inside      = 1.466;          // 腔内折射率
//     double n_outside     = 1.0;            // 腔外折射率
//     Polarization pol     = TM;             // TM 偏振
//     double phi_inc_deg   = 15.0;           // 入射角（度）
//     Point  incident_dir(std::cos(phi_inc_deg*M_PI/180.0),
//                         std::sin(phi_inc_deg*M_PI/180.0));

//     // —— 2) 构造两个圆角六边形腔体 ——
//     std::cout << "[Info] Constructing two rounded hexagonal cavities...\n";
//     Microcavity cav1 = make_rounded_hexagonal_cavity(R, corner_radius, n_inside);
//     Microcavity cav2 = make_rounded_hexagonal_cavity(R, corner_radius, n_inside);

//     auto elems1 = cav1.discretize(240, 24);
//     // auto elems1 = cav1.discretize(249, 18);
//     auto elems2 = cav2.discretize(240, 24);
//     // auto elems2 = cav2.discretize(249, 18);
//     std::cout << "[Info]  Cav1 elements: " << elems1.size()
//               << ", Cav2 elements: " << elems2.size() << "\n";
//     // 导出中点
//     Point off1(-0.9*R,  0.25*R),
//         off2( 0.9*R, -0.25*R);
//     save_midpoints_to_txt(elems1, off1, elems2, off2, "./data/cavities_midpoints.txt");
//     std::cout << "[Info] midpoints written to cavities_midpoints.txt\n";

//     // —— 4) 平移两腔体到指定中心 —— 
//     Point offset1(-0.9*R,  0.25*R);
//     Point offset2( 0.9*R, -0.25*R);
//     for (auto &e : elems1) {
//         e.p1  += offset1; e.p2  += offset1; e.mid += offset1;
//     }
//     for (auto &e : elems2) {
//         e.p1  += offset2; e.p2  += offset2; e.mid += offset2;
//     }

//     // 合并所有元素及分区信息
//     std::vector<Microcavity::BEMElement> all_elems;
//     all_elems.reserve(elems1.size() + elems2.size());
//     all_elems.insert(all_elems.end(), elems1.begin(), elems1.end());
//     all_elems.insert(all_elems.end(), elems2.begin(), elems2.end());
//     std::vector<std::vector<Microcavity::BEMElement>> per_cav = { elems1, elems2 };

//     // —— 5) 扫描 k ∈ [2,2.5] —— 
//     int    N_k    = 6;
//     double k_min  = 2.0, k_max = 2.5;
//     std::vector<double> k_vals(N_k), kR_vals(N_k), sigma_vals(N_k);

//     // 预先填好 k 和 kR=k*R
//     for (int i = 0; i < N_k; ++i) {
//         k_vals[i]  = k_min + (k_max - k_min) * i / (N_k - 1);
//         kR_vals[i] = k_vals[i] * R;
//     }

//     std::cout << "[Info] Starting parallel scan over " << N_k << " k values...\n";
//     std::mutex cout_mutex;

//     N_k=1;
//     #pragma omp parallel for schedule(dynamic)
//     for (int i = 0; i < N_k; ++i) {
//         double k = k_vals[i];
//         {
//             std::lock_guard<std::mutex> lock(cout_mutex);
//             std::cout << "[Step " << (i+1) << "/" << N_k << "] "
//                       << "k=" << std::fixed << std::setprecision(4) << k << " ... " << std::flush;
//         }

//         try {


//             // 5.2) 构造入射场向量
//             auto phi_psi_in = compute_incident_vectors(
//                 all_elems,
//                 k,
//                 incident_dir
//             );

//             // 5.1) 求解散射边界值问题
//             auto sol = solve_scattering( phi_psi_in,
//                 all_elems, per_cav,
//                 pol, k, incident_dir,
//                 n_inside, n_outside
//             );

//                     // === 输出 Δφ, Δψ 到 CSV ===
//             {
//                 size_t N_el = all_elems.size();
//                 // 拆分出 φ_in, ψ_in
//                 Eigen::VectorXcd phi_in = phi_psi_in.head(N_el);
//                 Eigen::VectorXcd psi_in = phi_psi_in.tail(N_el);
//                 // 创建文件名，方便区分
//                 std::ostringstream fname;
//                 fname << "./data/delta_" << std::setw(2) << std::setfill('0') << (i+1) << ".csv";
//                 std::ofstream ofs(fname.str());
//                 ofs << "idx,delta_phi_re,delta_phi_im,delta_psi_re,delta_psi_im\n";
//                 for (size_t j = 0; j < N_el; ++j) {
//                     std::complex<double> dphi = sol.phi[j] - phi_in[j];
//                     std::complex<double> dpsi = sol.psi[j] - psi_in[j];
//                     ofs 
//                     << j << ","
//                     << dphi.real() << "," << dphi.imag() << ","
//                     << dpsi.real() << "," << dpsi.imag() << "\n";
//                 }
//             }

//             // 5.3) 计算前向散射振幅 f(θ=φ_inc)
//             double phi_rad = phi_inc_deg * M_PI / 180.0;
//             auto f_scat = compute_scattering_amplitude(
//                 sol, all_elems, phi_psi_in,
//                 std::complex<double>(k,0.0),
//                 phi_rad,
//                 incident_dir
//             );

//             // 5.4) 根据光学定理得总散射截面 σ(k)
//             double sigma = compute_total_cross_section(f_scat, k);

//             // 保存并打印
//             {
//                 std::lock_guard<std::mutex> lock(cout_mutex);
//                 sigma_vals[i] = sigma;
//                 std::cout << "done. |f0|=" << std::scientific << std::setprecision(3)
//                           << std::abs(f_scat)
//                           << ", σ=" << std::fixed << std::setprecision(6) << sigma
//                           << "\n";
//             }
//         } catch (const std::exception& e) {
//             std::lock_guard<std::mutex> lock(cout_mutex);
//             std::cerr << "[Error] step " << i << " failed: " << e.what() << "\n";
//             sigma_vals[i] = NAN;
//         }
//     }

//     // —— 6) 写出归一化的 total cross section —— 
//     //      横坐标 kR，纵坐标 σ/R
//     std::vector<double> sigma_norm(N_k);
//     for (int i = 0; i < N_k; ++i) {
//         sigma_norm[i] = sigma_vals[i] / R;
//     }

//     std::cout << "[Info] Writing cross_section.txt ...\n";
//     write_cross_section(kR_vals, sigma_norm, "./data/cross_section.txt");

//     std::cout << "[Done] All outputs written.\n";
//     return 0;
// }

// #include <iostream>
// #include <vector>
// #include <complex>
// #include <cmath>
// #include <algorithm>
// #include <iomanip>
// #include <mutex>
// #include <omp.h>  // 使用 OpenMP

// #include "bem_microcavity.h"
// #include "bem_matrix.h"
// #include "bem_scattering.h"
// #include "bem_output.h"

// using namespace bem;

// int main() {
//     omp_set_num_threads(2);
//     std::cout << "[Info] Setting up parameters...\n";
//     double R             = 1.0;
//     double corner_radius = 0.02 * R;
//     double n_inside      = 1.466;
//     double n_outside     = 1.0;
//     Polarization pol     = TM;
//     double phi_inc_deg   = 15.0;
//     Point incident_dir(std::cos(phi_inc_deg * M_PI / 180.0),
//                        std::sin(phi_inc_deg * M_PI / 180.0));

//     std::cout << "[Info] Constructing two rounded hexagonal cavities...\n";
//     Microcavity cav1 = make_rounded_hexagonal_cavity(R, corner_radius, n_inside);
//     Microcavity cav2 = make_rounded_hexagonal_cavity(R, corner_radius, n_inside);

//     double approx_kR     = 22.5;
//     double approx_lambda = 2 * M_PI / approx_kR;
//     double target_ds     = approx_lambda / 16.0;
//     std::cout << "[Info] Discretizing boundaries with target ds = "
//               << std::setprecision(4) << target_ds << " ...\n";
//     auto elems1 = cav1.discretize(target_ds);
//     auto elems2 = cav2.discretize(target_ds);
//     std::cout << "[Info]  Cav1 elements: " << elems1.size()
//               << ", Cav2 elements: " << elems2.size() << "\n";

//     Point offset(1.8 * R, 0.5 * R);
//     for (auto &e : elems2) {
//         e.p1 += offset;
//         e.p2 += offset;
//         e.mid += offset;
//     }

//     std::vector<Microcavity::BEMElement> all_elems;
//     all_elems.reserve(elems1.size() + elems2.size());
//     all_elems.insert(all_elems.end(), elems1.begin(), elems1.end());
//     all_elems.insert(all_elems.end(), elems2.begin(), elems2.end());
//     std::vector<std::vector<Microcavity::BEMElement>> per_cav = { elems1, elems2 };

//     int N_k = 5;
//     double kR_min = 20.0, kR_max = 25.0;

//     std::vector<double> kR_vals(N_k), sigma_vals(N_k);
//     for (int i = 0; i < N_k; ++i) {
//         kR_vals[i] = kR_min + (kR_max - kR_min) * i / (N_k - 1);
//     }

//     std::cout << "[Info] Starting parallel scan over " << N_k << " kR values...\n";

//     std::mutex cout_mutex;
// #pragma omp parallel for schedule(dynamic)
//     for (int i = 0; i < N_k; ++i) {
//         double kR = kR_vals[i];
//         double k = kR;

//         try {
//             {
//                 std::lock_guard<std::mutex> lock(cout_mutex);
//                 std::cout << "[Step " << (i+1) << "/" << N_k << "] "
//                           << "Solving scattering (kR=" << std::fixed << std::setprecision(4) << kR << ") ... " << std::flush;
//             }

//             auto sol = solve_scattering(all_elems, per_cav, pol, k, incident_dir, n_inside, n_outside);

//             auto phi_psi_in = compute_incident_vectors(all_elems, k, incident_dir);

//             auto f0 = compute_scattering_amplitude(
//                 sol, all_elems, phi_psi_in,
//                 std::complex<double>(k, 0.0),
//                 std::vector<double>{ phi_inc_deg * M_PI / 180.0 },
//                 incident_dir
//             )[0];

//             double sigmaR = compute_total_cross_section(f0, k);

//             {
//                 std::lock_guard<std::mutex> lock(cout_mutex);
//                 sigma_vals[i] = sigmaR;
//                 std::cout << "done.\n";
//                 std::cout << "   debug: |f0|=" << std::abs(f0)
//                           << ", f0=" << f0 << "\n";
//                 std::cout << "  Total cross section σR = "
//                           << std::setprecision(6) << sigmaR << "\n\n";
//             }

//         } catch (const std::exception& e) {
//             std::lock_guard<std::mutex> lock(cout_mutex);
//             std::cerr << "[Error] Step " << i << ": " << e.what() << "\n";
//             sigma_vals[i] = std::nan("");
//         }
//     }

//     std::cout << "[Info] Writing cross_section.txt ...\n";
//     write_cross_section(kR_vals, sigma_vals, "cross_section.txt");

//     std::cout << "[Done] All outputs written.\n";
//     return 0;
// }
// // main.cpp

// #include <iostream>
// #include <vector>
// #include <complex>
// #include <cmath>
// #include <algorithm>
// #include <iomanip>              // 为了 std::setprecision

// #include "bem_microcavity.h"
// #include "bem_matrix.h"
// #include "bem_scattering.h"
// #include "bem_output.h"         // 包含 write_cross_section、write_far_field、write_near_field

// using namespace bem;

// int main() {
//     // —————————————————————————————————————
//     // 1) 参数设置
//     // —————————————————————————————————————
//     std::cout << "[Info] Setting up parameters...\n";
//     double R             = 1.0;           // 归一化六边形边长 R
//     double corner_radius = 0.02 * R;      // 平滑半径
//     double n_inside      = 1.466;         // 腔体内折射率
//     double n_outside     = 1.0;           // 外部折射率
//     Polarization pol     = TM;            // TM 偏振
//     double phi_inc_deg   = 15.0;          // 入射角 φ = 15°
//     Point  incident_dir(std::cos(phi_inc_deg*M_PI/180.0),
//                         std::sin(phi_inc_deg*M_PI/180.0));

//     // —————————————————————————————————————
//     // 2) 构造两个耦合六边形腔体
//     //    中心分别在 (0,0) 和 (1.8R, 0.5R)
//     // —————————————————————————————————————
//     std::cout << "[Info] Constructing two rounded hexagonal cavities...\n";
//     Microcavity cav1 = make_rounded_hexagonal_cavity(R, corner_radius, n_inside);
//     Microcavity cav2 = make_rounded_hexagonal_cavity(R, corner_radius, n_inside);

//     // —————————————————————————————————————
//     // 3) 离散化：目标弧长 ≈ λ/16
//     // —————————————————————————————————————
//     double approx_kR     = 22.5;
//     double approx_lambda = 2*M_PI / approx_kR;   // λ ≈ 2π/kR
//     double target_ds     = approx_lambda/16.0;  // 每段 ≈ λ/16
//     std::cout << "[Info] Discretizing boundaries with target ds = "
//               << std::setprecision(4) << target_ds << " ...\n";
//     auto elems1 = cav1.discretize(target_ds);
//     auto elems2 = cav2.discretize(target_ds);
//     std::cout << "[Info]  Cav1 elements: " << elems1.size()
//               << ", Cav2 elements: " << elems2.size() << "\n";

//     // 平移第二个腔体的所有 BEM 元素
//     Point offset(1.8*R, 0.5*R);
//     for (auto &e : elems2) {
//         e.p1  += offset;
//         e.p2  += offset;
//         e.mid += offset;
//     }

//     // 合并所有元素
//     std::vector<Microcavity::BEMElement> all_elems;
//     all_elems.reserve(elems1.size() + elems2.size());
//     all_elems.insert(all_elems.end(), elems1.begin(), elems1.end());
//     all_elems.insert(all_elems.end(), elems2.begin(), elems2.end());

//     // 每个腔体的元素分块
//     std::vector<std::vector<Microcavity::BEMElement>> per_cav = { elems1, elems2 };

//     // —————————————————————————————————————
//     // 4) 扫描 kR ∈ [20,25]，计算总散射截面 σ(kR)
//     // —————————————————————————————————————
//     int    N_k     = 5;
//     double kR_min  = 20.0, kR_max = 25.0;
//     std::vector<double> kR_vals, sigma_vals;
//     kR_vals.reserve(N_k);
//     sigma_vals.reserve(N_k);

//     std::cout << "[Info] Starting scan over kR from "
//               << kR_min << " to " << kR_max
//               << " in " << N_k << " steps...\n";

//     for (int i = 0; i < N_k; ++i) {
//         double kR = kR_min + (kR_max - kR_min) * i / (N_k - 1);
//         kR_vals.push_back(kR);

//         // — 将 kR 视作实数波数 k —
//         double k = kR;  

//         std::cout << "[Step " << (i+1) << "/" << N_k << "] "
//                   << "Solving scattering (kR=" << std::fixed << std::setprecision(4) << kR << ") ... " 
//                   << std::flush;

//         // 求解散射边界解，此时将实数 k 转成复数 (k,0)
//         auto sol = solve_scattering(
//             all_elems, per_cav,
//             pol,
//             k,
//             incident_dir,
//             n_inside, n_outside
//         );
//         std::cout << "done.\n";

//         // 计算入射场向量（同样把 k 作为实数传入）
//         std::cout << "  Computing incident vectors... " << std::flush;
//         auto phi_psi_in = compute_incident_vectors(
//             all_elems,
//             k,
//             incident_dir
//         );
//         std::cout << "done.\n";


//         // 计算前向散射振幅
//         std::cout << "  Computing forward scattering amplitude... " << std::flush;
//         auto f0 = compute_scattering_amplitude(
//             sol, all_elems, phi_psi_in,
//             std::complex<double>(k, 0.0),
//             std::vector<double>{ phi_inc_deg * M_PI / 180.0 },
//             incident_dir
//         )[0];
//         std::cout << "done.\n";

//         // 输出一下 |f0| 方便调试
//         std::cout << "   debug: |f0|=" << std::abs(f0)
//                   << ", f0=" << f0 << "\n";

//         // 根据光学定理计算总截面 σ(k)
//         double sigmaR = compute_total_cross_section(f0, k);
//         sigma_vals.push_back(sigmaR);
//         std::cout << "  Total cross section σR = "
//                   << std::setprecision(6) << sigmaR << "\n\n";
//     }

//     // 写出 cross_section.txt
//     std::cout << "[Info] Writing cross_section.txt ...\n";
//     write_cross_section(kR_vals, sigma_vals, "cross_section.txt");

//     std::cout << "[Done] All outputs written.\n";
//     return 0;
// }

                // // —— 在第一次 (i == 0) 时把 phi_psi_in 和 sol.dump 到文件 —— 
                // if (i == 0) {
                //     std::ofstream ofs("debug_solution_kR_" +
                //                       std::to_string(kR) + ".txt");
                //     ofs << "# idx   phi_in.re    phi_in.im    "
                //            "psi_in.re    psi_in.im    "
                //            "phi_sc.re    phi_sc.im    "
                //            "psi_sc.re    psi_sc.im\n";
                //     int N_el = int(all_elems.size());
                //     for (int idx = 0; idx < N_el; ++idx) {
                //         auto phi_in = phi_psi_in[idx];
                //         auto psi_in = phi_psi_in[idx + N_el];
                //         auto phi_sc = sol.phi[idx];
                //         auto psi_sc = sol.psi[idx];
                //         ofs
                //           << std::setw(4) << idx << "   "
                //           << std::fixed << std::setprecision(4)
                //           << std::setw(12) << phi_in.real() << " "
                //           << std::setw(12) << phi_in.imag() << "   "
                //           << std::setw(12) << psi_in.real() << " "
                //           << std::setw(12) << psi_in.imag() << "   "
                //           << std::setw(12) << phi_sc.real() << " "
                //           << std::setw(12) << phi_sc.imag() << "   "
                //           << std::setw(12) << psi_sc.real() << " "
                //           << std::setw(12) << psi_sc.imag()
                //           << "\n";
                //     }
                //     ofs.close();
                //     std::cout << "  [Debug] dumped phi_psi_in & sol to debug_solution_kR_"
                //               << kR << ".txt\n";
                // }



// // main.cpp

// #include <iostream>
// #include <vector>
// #include <complex>
// #include <cmath>
// #include <algorithm>
// #include <iomanip>              // <-- 为了 std::setprecision

// #include "bem_microcavity.h"
// #include "bem_matrix.h"
// #include "bem_scattering.h"
// #include "bem_output.h"         // 包含 write_cross_section、write_far_field、write_near_field

// using namespace bem;

// int main() {
//     // —————————————————————————————————————
//     // 1) 参数设置
//     // —————————————————————————————————————
//     std::cout << "[Info] Setting up parameters...\n";
//     double R             = 1.0;           // 归一化六边形边长 R
//     double corner_radius = 0.02 * R;      // 平滑半径
//     double n_inside      = 1.466;         // 腔体内折射率
//     double n_outside     = 1.0;           // 外部折射率
//     Polarization pol     = TM;            // TM 偏振
//     double phi_inc_deg   = 15.0;          // 入射角 φ=15°
//     Point  incident_dir(std::cos(phi_inc_deg*M_PI/180.0),
//                         std::sin(phi_inc_deg*M_PI/180.0));

//     // —————————————————————————————————————
//     // 2) 构造两个耦合六边形腔体
//     //    中心分别在 (0,0) 和 (1.8R, 0.5R)
//     // —————————————————————————————————————
//     std::cout << "[Info] Constructing two rounded hexagonal cavities...\n";
//     Microcavity cav1 = make_rounded_hexagonal_cavity(R, corner_radius, n_inside);
//     Microcavity cav2 = make_rounded_hexagonal_cavity(R, corner_radius, n_inside);

//     // —————————————————————————————————————
//     // 3) 离散化：目标弧长 ≈ λ/16
//     // —————————————————————————————————————
//     double approx_kR = 22.5;
//     double approx_lambda = 2*M_PI/(approx_kR);        // λ ≈ 2π/kR
//     double target_ds = approx_lambda / 8.0;        // 每段 ≈ λ/16
//     std::cout << "[Info] Discretizing boundaries with target ds = "
//               << std::setprecision(4) << target_ds << " ...\n";
//     auto elems1 = cav1.discretize(target_ds);
//     auto elems2 = cav2.discretize(target_ds);
//     std::cout << "[Info]  Cav1 elements: " << elems1.size()
//               << ", Cav2 elements: " << elems2.size() << "\n";

//     // 平移第二个腔体的所有 BEM 元素
//     Point offset(1.8*R, 0.5*R);
//     for (auto &e : elems2) {
//         e.p1  += offset;
//         e.p2  += offset;
//         e.mid += offset;
//     }

//     // 合并所有元素
//     std::vector<Microcavity::BEMElement> all_elems;
//     all_elems.reserve(elems1.size() + elems2.size());
//     all_elems.insert(all_elems.end(), elems1.begin(), elems1.end());
//     all_elems.insert(all_elems.end(), elems2.begin(), elems2.end());

//     // 每个腔体的元素分块
//     std::vector<std::vector<Microcavity::BEMElement>> per_cav = { elems1, elems2 };

//     // —————————————————————————————————————
//     // 4) 扫描 kR ∈ [20,25]，计算总散射截面 σ(kR)
//     // —————————————————————————————————————
//     int    N_k     = 200;
//     double kR_min  = 20.0, kR_max = 25.0;
//     std::vector<double> kR_vals, sigma_vals;
//     kR_vals.reserve(N_k);
//     sigma_vals.reserve(N_k);

//     std::cout << "[Info] Starting scan over kR from "
//               << kR_min << " to " << kR_max
//               << " in " << N_k << " steps...\n";

//     for (int i = 0; i < N_k; ++i) {
//         double kR = kR_min + (kR_max - kR_min)*i/(N_k-1);
//         kR_vals.push_back(kR);
//         std::complex<double> k = kR;  // R=1 归一化后

//         // — 重计算1：求解散射边界解 —
//         std::cout << "[Step " << (i+1) << "/" << N_k << "] "
//                   << "Solving scattering (kR=" << std::setprecision(4) << kR << ") ... "
//                   << std::flush;
//         auto sol = solve_scattering(
//             all_elems, per_cav,
//             pol, k, incident_dir,
//             n_inside, n_outside
//         );
//         std::cout << "done.\n";

//         // — 重计算2：入射场向量 —
//         std::cout << "  Computing incident vectors... " << std::flush;
//         auto phi_psi_in = compute_incident_vectors(
//             all_elems, k, incident_dir
//         );
//         std::cout << "done.\n";

//         // — 重计算3：散射振幅 —
//         std::cout << "  Computing forward scattering amplitude... " << std::flush;
//         auto f0 = compute_scattering_amplitude(
//             sol, all_elems, phi_psi_in, k,
//             std::vector<double>{ phi_inc_deg*M_PI/180.0 },
//             incident_dir
//         )[0];
        
//         std::cout << "done.\n";

//         std::cout<<"   debug: |f0|="<<std::abs(f0)<<",  f0="<<f0<<"\n";

//         // — 计算总截面 —
//         double sigmaR = compute_total_cross_section(f0, k);
//         sigma_vals.push_back(sigmaR);
//         std::cout << "  Total cross section σR = "
//                   << std::setprecision(6) << sigmaR << "\n\n";
//     }

//     // 写出 cross_section.txt
//     std::cout << "[Info] Writing cross_section.txt ...\n";
//     write_cross_section(kR_vals, sigma_vals, "cross_section.txt");

//     std::cout << "[Done] All outputs written.\n";
//     return 0;
// }


    // —————————————————————————————————————
    // 5) 找到主峰 k1，做牛顿迭代求精确共振 k_res
    //    论文示例：初始猜测 k1R = 22.95 - i·0.098
    // —————————————————————————————————————
    // std::complex<double> k1 = { 22.95, -0.098 };
    // double tol      = 1e-8;
    // double delta    = 1e-3;
    // int    max_iter = 20;
    // std::complex<double> k_res = k1;
    // for (int iter = 0; iter < max_iter; ++iter) {
    //     // 行列式 g(k)
    //     auto [M, M0] = build_matrix(all_elems, per_cav,
    //                                 pol, k_res,
    //                                 n_inside, n_outside);
    //     std::complex<double> g  = M.determinant();
    //     // g(k+Δ)
    //     auto [M1, M01] = build_matrix(all_elems, per_cav,
    //                                   pol, k_res + delta,
    //                                   n_inside, n_outside);
    //     std::complex<double> g1 = M1.determinant();
    //     // g(k + iΔ)
    //     auto [M2, M02] = build_matrix(all_elems, per_cav,
    //                                   pol, k_res + std::complex<double>(0,delta),
    //                                   n_inside, n_outside);
    //     std::complex<double> g2 = M2.determinant();
    //     // 数值微分 g'
    //     std::complex<double> gp = (g1 - g)/(2.0*delta)
    //                               - std::complex<double>(0,1)*(g2 - g)/(2.0*delta);
    //     // Newton 更新
    //     std::complex<double> k_new = k_res - g/gp;
    //     if (std::abs(k_new - k_res) < tol) {
    //         k_res = k_new;
    //         break;
    //     }
    //     k_res = k_new;
    // }
    // std::cout << "Found resonance: k_res R = "
    //           << k_res.real() << "  "
    //           << k_res.imag() << " i\n";

    // —————————————————————————————————————
    // 6) 用 M(k_res) 求零特征向量 [φ_l, ψ_l]
    // —————————————————————————————————————
    // {
    //     auto [M_res, M0_res] = build_matrix(all_elems, per_cav,
    //                                         pol, k_res,
    //                                         n_inside, n_outside);
    //     // SVD 得到 nullspace
    //     Eigen::JacobiSVD<Eigen::MatrixXcd> svd(
    //         M_res, Eigen::ComputeFullV);
    //     auto V      = svd.matrixV();
    //     Eigen::VectorXcd sol0 = V.col(V.cols()-1); // 最小奇异值对应的特征向量
    //     int N = all_elems.size();
    //     Eigen::VectorXcd phi_sol = sol0.head(N);
    //     Eigen::VectorXcd psi_sol = sol0.tail(N);

        // ——————————————————————————————
        // 7) 重构近场 (公式 38)
        //    在 [x,y] 网格上计算 |ψ(r)|^2
        // ——————————————————————————————
        // const int Nx = 101, Ny = 101;
        // std::vector<double> xs(Nx), ys(Ny);
        // for (int i = 0; i < Nx; ++i) xs[i] = -2.0*R + 4.0*R*i/(Nx-1);
        // for (int j = 0; j < Ny; ++j) ys[j] = -2.0*R + 4.0*R*j/(Ny-1);
        // std::vector<std::complex<double>> near_field(Nx*Ny);
        // // 对每个观察点用数值积分重构 ψ(r)
        // for (int ix = 0; ix < Nx; ++ix) {
        //     for (int jy = 0; jy < Ny; ++jy) {
        //         Point r{ xs[ix], ys[jy] };
        //         std::complex<double> val = 0.0;
        //         // ∑ ψ_l ∫ ∂νG  –  ∑ φ_l ∫ G
        //         for (int l = 0; l < N; ++l) {
        //             const auto &e = all_elems[l];
        //             // 2 种积分核
        //             auto Cint = GaussIntegrator::integrate(
        //                 e, r,
        //                 GaussIntegrator::kernel_tm_C<CurveSegment>,
        //                 k_res);
        //             auto Gint = GaussIntegrator::integrate(
        //                 e, r,
        //                 GaussIntegrator::kernel_tm_B<CurveSegment>,
        //                 k_res);
        //             val += psi_sol[l] * Cint - phi_sol[l] * Gint;
        //         }
        //         near_field[ix*Ny + jy] = val;
        //     }
        // }
        // write_near_field(xs, ys, near_field, "near_field.txt");
        
        // ——————————————————————————————
        // 8) 计算远场发射 (公式 19)
        // ——————————————————————————————
        // const int Nθ = 360;
        // std::vector<double> thetas(Nθ);
        // for (int t = 0; t < Nθ; ++t) thetas[t] = 2*M_PI * t/Nθ;
        // auto f_vals = compute_scattering_amplitude(
        //     {phi_sol, psi_sol},
        //     all_elems,
        //     phi_psi_in,   // 重用前面计算的入射向量
        //     k_res,
        //     thetas,
        //     incident_dir
        // );
        // write_far_field(thetas, f_vals, "far_field.txt");
    // }


// // main.cpp
// #include <iostream>
// #include <vector>
// #include <complex>
// #include <cmath>

// #include "bem_microcavity.h"
// #include "bem_matrix.h"
// #include "bem_scattering.h"
// #include "bem_output.h"

// int main() {
//     using namespace bem;

//     // —— 1) 构造几何 & 材料参数 —— 
//     double side_length = 1.0;            // 六边形边长（单位任意）
//     double n_inside   = 1.466;           // 腔体折射率
//     double n_outside  = 1.0;             // 外部折射率
//     Polarization pol  = TM;              // TM 偏振
//     Point incident_dir(1.0, 0.0);        // 入射方向：x 轴正向

//     // 便捷构造：一个六边形微腔
//     Microcavity cav = make_rounded_hexagonal_cavity(side_length, n_inside);
//     std::vector<Microcavity> cavities = { cav };

//     // —— 2) 离散化边界 —— 
//     //    返回 flat_all：所有元素扁平列表
//     //           per_cav：每个腔体对应的元素列表
//     auto [flat_all, per_cav] = build_discretize(cavities, /*n_per_seg=*/32);

//     // —— 3) 扫描 k 区间，计算总散射截面 σ(k) —— 
//     std::vector<double> k_vals, sigma_vals;
//     double k_min = 20.0, k_max = 25.0;
//     int    N_k   = 200;
//     k_vals.reserve(N_k);
//     sigma_vals.reserve(N_k);

//     for(int i = 0; i < N_k; ++i) {
//         double k = k_min + (k_max - k_min) * i / (N_k - 1);
//         k_vals.push_back(k);

//         // 3.1) 求解散射问题，得到边界解 [φ,ψ]
//         auto result = solve_scattering(flat_all, per_cav, pol, k, incident_dir, n_inside, n_outside);

//         // 3.2) 计算前向散射幅度 f(θ=0)
//         auto phi_psi_in = compute_incident_vectors(flat_all, k, incident_dir);
//         std::vector<double> theta0 = { 0.0 };
//         auto f0 = compute_scattering_amplitude(result, flat_all, phi_psi_in, k, theta0, incident_dir)[0];

//         // 3.3) 通过光学定理得到 σ(k)
//         double sigma = compute_total_cross_section(f0, k);
//         sigma_vals.push_back(sigma);
//     }

//     // 写出 σ(k) 到 cross_section.txt
//     write_cross_section(k_vals, sigma_vals, "cross_section.txt");

//     // —— 4) 在总截面峰值附近取一个 k_peak，计算并输出远场 —— 
//     // 简单取中点
//     int    ip = N_k / 2;
//     double k_peak = k_vals[ip];

//     auto sol_peak     = solve_scattering(flat_all, per_cav, pol, k_peak, incident_dir, n_inside, n_outside);
//     auto phi_psi_in_p = compute_incident_vectors(flat_all, k_peak, incident_dir);

//     // 设置远场角度列表
//     std::vector<double> thetas;
//     const int Nθ = 360;
//     thetas.reserve(Nθ);
//     for(int j=0; j<Nθ; ++j) {
//         thetas.push_back(2*M_PI * j / Nθ);
//     }

//     auto f_vals = compute_scattering_amplitude(sol_peak, flat_all, phi_psi_in_p, k_peak, thetas, incident_dir);
//     write_far_field(thetas, f_vals, "far_field.txt");

//     std::cout << "Done. Data written to:\n"
//               << "  cross_section.txt\n"
//               << "  far_field.txt\n";

//     return 0;
// }
