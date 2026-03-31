// #include <boost/math/special_functions/math_fwd.hpp>
// #include <iostream>
// #include <vector>
// #include <complex>
// #include <cmath>
// #include <iomanip>
// #include <mutex>
// #include <omp.h>
// #include <algorithm>
// #include <map>
// #include <fstream>
// // #include "bem_gauss_integrator.h"
// // #include "bem_hankel_manual.h"
// #include "bem_curve_segment.h"
// #include "bem_hankel_manual.h"
// #include "bem_microcavity.h"
// #include "bem_matrix.h"
// #include "bem_scattering.h"
// #include "bem_output.h"
// #include "bem_resonance.h"

// #include <tuple>

// #include "bem_gauss_integrator.h"   // complex_integrate(), kernels
// #include "../include/eigen-3.4.0/Eigen/Dense"
// #include "../include/eigen-3.4.0/Eigen/SVD"

// using namespace bem;

// // main.cpp
#include <cmath>
#include <iostream>
#include <fstream>
#include <complex>
#include <vector>
#include <iomanip>

#include "bem_microcavity.h"
#include "bem_matrix.h"
#include "bem_scattering.h"
#include "bem_resonance.h"
#include "bem_output.h"

using namespace bem;

int main() {


    const double R = 10.0;
    const double corner_radius = 0.001;  // 可按需调整
    const double n_inside  = 1.466, n_outside = 1.0;
    Polarization pol = TM;
    Microcavity cav1 = make_hermite_rounded_hexagonal_cavity(R, corner_radius, n_inside);
    Microcavity cav2 = make_hermite_rounded_hexagonal_cavity(R, corner_radius, n_inside);

    // 平移至指定位置
    cav1.move(Point(-0.9*R,  0.25*R));
    cav2.move(Point( 0.9*R, -0.25*R));

    // 2) 离散化
    auto elems1 = cav1.discretize(/*n_line_per_edge=*/240, /*n_arc_per_corner=*/36,10.0/240,18);
    auto elems2 = cav2.discretize(/*n_line_per_edge=*/240, /*n_arc_per_corner=*/36,10.0/240,18);

    std::vector<Microcavity::BEMElement> all_elems;
    all_elems.reserve(elems1.size()+elems2.size());
    all_elems.insert(all_elems.end(), elems1.begin(), elems1.end());
    all_elems.insert(all_elems.end(), elems2.begin(), elems2.end());
    std::vector<std::vector<Microcavity::BEMElement>> per_cav = { elems1, elems2 };
        // 1) 装配 double 精度矩阵
    auto res = find_resonance(all_elems, per_cav, pol, std::complex<double>(2.295,-0.0098));

    // // 快速估条件数
    // Eigen::FullPivLU<Eigen::MatrixXcd> lu_cond(M_d);
    // double rcond = lu_cond.rcond();
    // double condM = (rcond>0?1.0/rcond:std::numeric_limits<double>::infinity());
    // std::cout << "[debug] approx cond(M) = " << condM << std::endl;


    // auto res = resonance(all_elems, per_cav, TM, std::complex<double>( 2.294444,-0.009696));
    // std::ofstream fres("./dataf/resonancetest.txt");
    // for(int i=0;i<all_elems.size();i++){
    //     fres<<res.phi[i]<<","<<res.psi[i]<<"\n";
    // }
    // fres.close();
    // 3) 求共振态 φ,ψ
    // auto res = load_and_normalize_txt("./dataf/resonance.txt", std::complex<double>( 2.294444,-0.009696));
    
    // // --- 5) 输出边界元中心点
    // std::ofstream felems("./dataf/elements_midpoints.txt");
    // for(auto &e : all_elems){
    //     felems << e.mid.x_ << "," << e.mid.y_ << "\n";
    // }
    // felems.close();
    // // --- 6) 在区域上做灰度网格
    // std::vector<Point> cavity1,cavity2;
    // cavity1.reserve(6);
    // cavity2.reserve(6);
    // for (int i=0;i<6;i++){
    //     double theta=i*M_PI/3;
    //     Point ini = Point::polar(R, theta);
    //     cavity1.push_back(ini+Point(-0.9*R,  0.25*R));
    //     cavity2.push_back(ini+Point(0.9*R,  -0.25*R));
    // }
    // compute_psi2_grid_txt(all_elems, res, {cavity1,cavity2}, {(int)elems1.size(),(int)elems2.size()}, 400, 300, -20.0, 20.0, -15.0, 15.0, "./dataf/compute_psi2_grid.txt");
    // std::cout<<"All data written.\n";
    // 4) 在 θ ∈ [0,π] 上采样
    // compute_psi2_on_circle_txt(res, all_elems, 100*R, 361, -M_PI, M_PI, "./dataf/compute_psi2_circle.txt");
    return 0;
}

// // Filename: compute_resonant_mode.cpp
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
//     std::complex<double> k_res(2.294444, -0.009696);

//     // --- 2) 构造两个微腔并离散化（你原来调用的方式）
//     double R = 10.0, corner_radius = 0.01, n_inside = 1.466, n_outer = 1.0;
//     Microcavity cav1 = make_hermite_rounded_hexagonal_cavity(R, corner_radius, n_inside);
//     Microcavity cav2 = make_hermite_rounded_hexagonal_cavity(R, corner_radius, n_inside);
//     cav1.move({-0.9 * R, 0.25 * R});
//     cav2.move({0.9 * R, -0.25 * R});
//     auto elems1 = cav1.discretize(240, 24, 10 / 240.0, 18);
//     auto elems2 = cav2.discretize(240, 24, 10 / 240.0, 18);

//     // 合并
//     std::vector<Microcavity::BEMElement> all_elems;
//     all_elems.reserve(elems1.size() + elems2.size());
//     all_elems.insert(all_elems.end(), elems1.begin(), elems1.end());
//     all_elems.insert(all_elems.end(), elems2.begin(), elems2.end());
//     std::vector<std::vector<Microcavity::BEMElement>> per_cav{elems1, elems2};

//     // 定义特定区域
//     double x_angle = M_PI / 3; // 60 度转换为弧度
//     double center_x = R * std::cos(x_angle) - 0.9 * R;
//     double center_y = R * std::sin(x_angle) + 0.25 * R;
//     double width = R / 250;
//     double height = R / 250;
//     double x_min = center_x - width ;
//     double x_max = center_x + width ;
//     double y_min = center_y - height;
//     double y_max = center_y + height;

//     // 将离散元中心点存储到 txt 文件
//     std::ofstream outFile("./dataf/center_points.txt");
//     if (outFile.is_open()) {
//         for (const auto& elem : all_elems) {
//             double x = elem.mid.x_;
//             double y = elem.mid.y_;
//             // 检查点是否在特定区域内
//             if (x >= x_min && x <= x_max && y >= y_min && y <= y_max) {
//                 outFile << x << " " << y << std::endl;
//             }
//         }
//         outFile.close();
//         std::cout << "Center points saved to center_points.txt" << std::endl;
//     } else {
//         std::cerr << "Unable to open file" << std::endl;
//     }

//     return 0;
// }