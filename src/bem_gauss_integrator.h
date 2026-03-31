// 文件名：bem_gauss_integrator.h
#ifndef BEM_GAUSS_INTEGRATOR_H
#define BEM_GAUSS_INTEGRATOR_H

#include "bem_constants.h"
#include "bem_curve_segment.h"
#include "bem_microcavity.h"
#include "bem_hankel_manual.h"
#include <boost/math/special_functions/math_fwd.hpp>
#include <vector>
#include <cmath>
#include <complex>
#include <functional>
#include <stdexcept>

#include <boost/math/special_functions/hankel.hpp>
using boost::math::cyl_hankel_1;


namespace bem {

// C_ll 曲率修正 (公式 32)
inline double compute_C_ll(
    const Microcavity::BEMElement& elem,
    double eps = EPS
) {
    double Δs = elem.length;
    double κ  = elem.curvature;
    return -1.0 + (κ * Δs) / (2.0 * M_PI);
}
// inline std::complex<double> compute_C_ll_out_psi_l(
//     const Microcavity::BEMElement& elem,
//     double k,
//     const Point& incident_dir,
//     double eps = EPS
// ) {
//     // 1) 先计算解析 C_ll（实数）
//     double Cll = -compute_C_ll(elem, eps)-2;

//     // 2) 在参数 t 上取三点
//     double t0 = elem.seg_start_t;
//     double t1 = elem.seg_end_t;
//     double tm = 0.5*(t0 + t1);

//     auto seg = elem.seg_ptr;
//     auto sample_psi = [&](double t){
//         Point s = (*seg)(t);
//         Point k_dir = incident_dir.normalized();
//         double dot = k_dir.x_*s.x_ + k_dir.y_*s.y_;
//         return std::exp(std::complex<double>(0, k * dot));
//     };

//     // 3) 三点取样
//     std::complex<double> psi0 = sample_psi(t0);
//     std::complex<double> psim = sample_psi(tm);
//     std::complex<double> psi1 = sample_psi(t1);

//     // 4) Simpson 加权平均
//     std::complex<double> psi_avg = (psi0 + std::complex<double>(4.0)*psim + psi1) / 6.0;

//     // 5) 返回加权后的 ∫C(s_i,s) ψ(s) ds ≈ C_ll * ψ_avg
//     return std::complex<double>(Cll) * psi_avg;
// }


// B_ll 高阶展开 (公式 34)
// 这里仍接受 std::complex<double> k，但在本阶段 k 实为实数
inline std::complex<double> compute_B_ll(
    const Microcavity::BEMElement& elem,
    std::complex<double> k,
    double n_j = 1.466,
    double eps = EPS
) {
    double Δs = elem.length;
    double z = n_j * k.real() * (Δs/4.0);
    double lnz = std::log(z);
    return (Δs / M_PI) * (1.0 - lnz+ std::complex<double>(0.0, M_PI/2.0)-EULER_GAMMA);
}

inline std::complex<double> compute_B_ll_out_phi_l(
    const Microcavity::BEMElement& elem,
    double k,
    const Point& incident_dir,
    double n_j = 1.466,
    double eps = EPS
) {
    // 1) 先计算解析 B_ll
    std::complex<double> Bll = compute_B_ll(elem, std::complex<double>(k,0), n_j, eps);

    // 2) 在参数 t 上取三点
    double t0 = elem.seg_start_t;
    double t1 = elem.seg_end_t;
    double tm = 0.5*(t0 + t1);

    auto seg = elem.seg_ptr;
    auto sample_phi = [&](double t){
        Point s = (*seg)(t);
        // 计算 ψ_in 的法向导数 φ(s) = ∂ν ψ_in = i k (ν·k_dir) e^{i k (k_dir·s)}
        Point nu = -1.0*seg->normal(t);                      // 取 “外” 法线
        Point k_dir = incident_dir.normalized();
        double dot = k_dir.x_*s.x_ + k_dir.y_*s.y_;
        double nu_dot_k = nu.dot(k_dir);
        std::complex<double> phase = std::exp(std::complex<double>(0, k * dot));
        return std::complex<double>(0,1) * k * nu_dot_k * phase;
    };

    // 3) 三点取样
    std::complex<double> phi0 = sample_phi(t0);
    std::complex<double> phim = sample_phi(tm);
    std::complex<double> phi1 = sample_phi(t1);

    // 4) Simpson 加权平均
    std::complex<double> phi_avg = (phi0 + std::complex<double>(4.0)*phim + phi1) / 6.0;

    // 5) 返回加权后的 ∫B(s_i,s) φ(s) ds ≈ B_ll * φ_avg
    return Bll * phi_avg;
}

// inline std::complex<double> compute_C_ll_out_psi_l(
//     const Microcavity::BEMElement& elem,
//     double k,
//     const Point& incident_dir,
//     double eps = EPS
// ) {
//     // 1) 计算 C_ll 并包含 “-2” 修正
//     double Cll = -compute_C_ll(elem, eps) - 2.0;

//     // 2) 参数区间 [t0,t1] 五等分
//     double t0 = elem.seg_start_t;
//     double t1 = elem.seg_end_t;
//     double h  = (t1 - t0) / 4.0;

//     auto seg = elem.seg_ptr;
//     auto sample_psi = [&](double t){
//         Point s = (*seg)(t);
//         Point k_dir = incident_dir.normalized();
//         double dot = k_dir.x_*s.x_ + k_dir.y_*s.y_;
//         return std::exp(std::complex<double>(0, k * dot));
//     };

//     // 3) 五点取样
//     std::complex<double> ψ0 = sample_psi(t0);
//     std::complex<double> ψ1 = sample_psi(t0 + h);
//     std::complex<double> ψ2 = sample_psi(t0 + 2*h);
//     std::complex<double> ψ3 = sample_psi(t0 + 3*h);
//     std::complex<double> ψ4 = sample_psi(t1);

//     // 4) Boole’s 公式加权
//     std::complex<double> ψ_avg =
//         (7.0*ψ0 + 32.0*ψ1 + 12.0*ψ2 + 32.0*ψ3 + 7.0*ψ4) / 90.0;

//     // 5) 返回近似 ∫ C(s_i,s) ψ(s) ds
//     return Cll * ψ_avg;
// }
inline std::complex<double> compute_C_ll_out_psi_l(
    const Microcavity::BEMElement& elem,
    double k,
    const Point& incident_dir,
    double eps = EPS
) {
    auto r = elem.mid;
    auto k_dir = incident_dir.normalized();
    auto i = std::complex<double>(0,1);
    auto alpha=std::exp(i*k*k_dir.dot(r));
    auto mid_t=(elem.seg_start_t+elem.seg_end_t)/2;
    auto tan =(elem.seg_ptr->tangent(mid_t)).normalized();
    auto beta = k*k_dir.dot(tan);
    auto epsilon = elem.length/2;
    auto J=alpha*2.0*std::sin(beta*epsilon)/(beta);
    return -alpha-elem.curvature*J/(2*M_PI);
}

// inline std::complex<double> compute_B_ll_out_phi_l(
//     const Microcavity::BEMElement& elem,
//     double k,
//     const Point& incident_dir,
//     double n_j = 1.466,
//     double eps = EPS
// ) {
//     auto k_dir = incident_dir.normalized();
//     auto mid_t=(elem.seg_start_t+elem.seg_end_t)/2;
//     auto v = elem.seg_ptr->normal(mid_t);
//     auto alpha = k*k_dir.dot(v);
//     auto tan =elem.seg_ptr->tangent(mid_t).normalized();
//     auto beta = k*k_dir.dot(tan);
//     auto epsilon = elem.length/2;
//     auto i=std::complex<double>(0,1);
//     auto I=2*alpha*std::exp(i*k*k_dir.dot(elem.mid))/beta*
//     (i*std::sin(beta*epsilon)*std::log(tan.norm())+std::cos(beta*epsilon)*std::log(epsilon)
// +2*EULER_GAMMA+std::log(epsilon)+2*std::log(std::abs(beta)));
//     auto J = 2.0*i*k*k_dir.dot(v)*std::exp(i*k*k_dir.dot(elem.mid))*epsilon;
//     auto res =  -I/M_PI+(-std::log(n_j*k/2)/M_PI+i/2.0-EULER_GAMMA/M_PI)*J;
//     return -res;
// }


struct GaussPoint {
    double t, w;
};

class GaussIntegrator {
public:
    static std::vector<GaussPoint> get_points(int order) {
        if (order == 5) {
            return {{0.04691007703066798,0.1184634425},
                    {0.2307653449471585, 0.2393143352},
                    {0.5,                0.2844444444},
                    {0.7692346550528415, 0.2393143352},
                    {0.9530899229693320, 0.1184634425}};
        } else if (order == 7) {
            return {        { 0.025446043828620, 0.064742483084435 },
            { 0.129234407200303, 0.139852695744638 },
            { 0.297077424311301, 0.190915025252559 },
            { 0.500000000000000, 0.208979591836735 },
            { 0.702922575688699, 0.190915025252559 },
            { 0.870765592799697, 0.139852695744638 },
            { 0.974553956171380, 0.064742483084435 }};
        } else if (order==15) {
            return {            {0.006003740989257, 0.0153766209980586},
            {0.031363303799647, 0.0351830237440540},
            {0.075896708294212, 0.0535796102335859},
            {0.137791134319915, 0.0697853389630771},
            {0.214513913695521, 0.0831346029084965},
            {0.302924326461218, 0.0930845029494211},
            {0.399402953001282, 0.0992157426635558},
            {0.5,                0.1012891209627806},
            {0.600597046998718, 0.0992157426635558},
            {0.697075673538782, 0.0930845029494211},
            {0.785486086304479, 0.0831346029084965},
            {0.862208865680085, 0.0697853389630771},
            {0.924103291705788, 0.0535796102335859},
            {0.968636696200353, 0.0351830237440540},
            {0.993996259010743, 0.0153766209980586}};
        }
        throw std::invalid_argument("Unsupported integration order");
    }
    template<typename Kernel>
    static std::complex<double> complex_integrate(
        const Microcavity::BEMElement& elem,
        const Point& s_i,
        Kernel&& kernel,
        std::complex<double> k,
        double index,
        const Point& incident_dir={0,0}
    ){
        // 相邻单元很近时用高阶 (7 点)，否则用 5 点
        double dist = elem.mid.distance(s_i);
        int order = (dist < 10*elem.length ? 15 : 15);

        auto points = get_points(order);
        std::complex<double> sum = 0.0;

        const auto& seg = *elem.seg_ptr;
        double t0 = elem.seg_start_t, t1 = elem.seg_end_t;
        for (auto& gp : points) {
            double t = t0 + gp.t*(t1 - t0);
            Point s   = seg(t);
            double distance = (s-s_i).norm();
            Point n   = seg.normal(t);
            double dot = n.dot(s-s_i);
            sum += gp.w *(seg.tangent(t).norm())* kernel(distance,index,k,dot,incident_dir,elem,t);
        }
        return sum*(t1-t0);
    }
    /**
     * @param elem   待积分的边界元
     * @param s_i    源点
     * @param kernel 自定义核函数 (t,s,n,dsdt,s_i,k) -> complex
     * @param k      实数波数
     */
    template<typename Kernel>
    static std::complex<double> integrate(
        const Microcavity::BEMElement& elem,
        const Point& s_i,
        Kernel&& kernel,
        double k,
        double index,
        const Point& incident_dir={0,0}
    ) {
        // 相邻单元很近时用高阶 (7 点)，否则用 5 点
        double dist = elem.mid.distance(s_i);
        int order = (dist < 10*elem.length ? 15 : 15);

        auto points = get_points(order);
        std::complex<double> sum = 0.0;

        const auto& seg = *elem.seg_ptr;
        double t0 = elem.seg_start_t, t1 = elem.seg_end_t;
        for (auto& gp : points) {
            double t = t0 + gp.t*(t1 - t0);
            Point s   = seg(t);
            double distance = (s-s_i).norm();
            Point n   = seg.normal(t);
            double dot = n.dot(s-s_i);
            sum += gp.w *(seg.tangent(t).norm())* kernel(distance,index,k,dot,incident_dir,elem,t);
        }
        return sum*(t1-t0);
    }

    // TM 偏振 B 矩阵积分核（-2G）
    template<typename Curve>
    static std::complex<double> kernel_tm_B(
        const double& distance, 
        const double& index,
        double k,
        double dot,
        const Point& incident_dir={0,0},
        const Microcavity::BEMElement& elem={},
        const double& t=0
    ) {
        return std::complex<double>(0, 0.5) * cyl_hankel_1(0, index*k*distance );
    }
    // TM 偏振 B 矩阵积分核（-2G）
    template<typename Curve>
    static std::complex<double> complex_kernel_tm_B(
        const double& distance, 
        const double& index,
        std::complex<double> k,
        double dot,
        const Point& incident_dir={0,0},
        const Microcavity::BEMElement& elem={},
        const double& t=0
    ) {
        return std::complex<double>(0, 0.5) * complex_hankel_1(0, index*k*distance );
    }
    template<typename Curve>
    static std::complex<double> kernel_tm_B_il_out_phi_l(
        const double& distance, 
        const double& index,
        double k,
        double dot,
        const Point& incident_dir={0,0},
        const Microcavity::BEMElement& elem={},
        const double& t=0
    ) {
        auto seg = elem.seg_ptr;
        Point r = (*seg)(t);
        double mid_t = t;
        Point nu = -1.0*elem.seg_ptr->normal(mid_t);
        Point k_dir = incident_dir.normalized();
        double dot_ = k_dir.x_ * r.x_ + k_dir.y_ * r.y_;
        double nu_dot_k = nu.dot(k_dir);
        std::complex<double> phi = std::complex<double>(0,1)*k*nu_dot_k*std::exp(std::complex<double>(0,k*dot_));
        return std::complex<double>(0, 0.5) * cyl_hankel_1(0, index*k*distance )*phi;
    }

    template<typename Curve>
    static std::complex<double> kernel_tm_X_il_out_x_l(
        const double& distance, 
        const double& index,
        double k,
        double n_dot_r,
        const Point& incident_dir={0,0},
        const Microcavity::BEMElement& elem={},
        const double& t=0
    ) {
        auto seg = elem.seg_ptr;
        Point r = (*seg)(t);
        Point nu = -1.0*elem.seg_ptr->normal(t);
        Point k_dir = incident_dir.normalized();
        double dot_ = k_dir.x_ * r.x_ + k_dir.y_ * r.y_;
        double nu_dot_k = nu.dot(k_dir);
        std::complex<double> psi = std::exp(std::complex<double>(0,k*dot_));
        std::complex<double> phi = std::complex<double>(0,1)*k*nu_dot_k*psi;
        return std::complex<double>(0, 0.5) * cyl_hankel_1(0, index*k*distance )*phi
        +std::complex<double>(0,-0.5)*index*k*n_dot_r/distance*cyl_hankel_1(1,index*k*distance)*psi;
        // return -std::complex<double>(0, 0.5) * cyl_hankel_1(0, index*k*distance )*phi
        // +std::complex<double>(0,0.5)*index*k*dot/distance*cyl_hankel_1(1,index*k*distance)*psi;
    }

    // TM 偏振 C 矩阵积分核（2∂G/∂ν）
    template<typename Curve>
    static std::complex<double> kernel_tm_C(
        const double& distance, 
        const double& index,
        double k,
        double dot,
        const Point& incident_dir={0,0},
        const Microcavity::BEMElement& elem={},
        const double& t=0
    ) {
        return std::complex<double>(0,0.5)*index*k*dot/distance*cyl_hankel_1(1,index*k*distance);
    }
    // TM 偏振 C 矩阵积分核（2∂G/∂ν）
    template<typename Curve>
    static std::complex<double> complex_kernel_tm_C(
        const double& distance, 
        const double& index,
        std::complex<double> k,
        double dot,
        const Point& incident_dir={0,0},
        const Microcavity::BEMElement& elem={},
        const double& t=0
    ) {
        return std::complex<double>(0,0.5)*index*k*dot/distance*complex_hankel_1(1,index*k*distance);
    }

    template<typename Curve>
    static std::complex<double> kernal_partial_G(
        const double& distance, 
        const double& index,
        std::complex<double> k,
        double dot,
        const Point& incident_dir={0,0},
        const Microcavity::BEMElement& elem={},
        const double& t=0
    ) {
        return std::complex<double>(0,0.25)*index*k*dot/distance*complex_hankel_1(1,index*k*distance);
    }

    template<typename Curve>
    static std::complex<double> kernel_G(
        const double& distance, 
        const double& index,
        std::complex<double> k,
        double dot,
        const Point& incident_dir={0,0},
        const Microcavity::BEMElement& elem={},
        const double& t=0
    ) {
        return std::complex<double>(0,-0.25)*complex_hankel_1(0,index*k*distance);
    }

    // TM 偏振 C 矩阵积分核（2∂G/∂ν）
    template<typename Curve>
    static std::complex<double> kernel_tm_C_il_out_psi_l(
        const double& distance, 
        const double& index,
        double k,
        double dot,
        const Point& incident_dir={0,0},
        const Microcavity::BEMElement& elem={},
        const double& t=0
    ) {
        auto seg = elem.seg_ptr;
        Point r = (*seg)(t);
        Point k_dir = incident_dir.normalized();
        double dot_ = k_dir.x_ * r.x_ + k_dir.y_ * r.y_;
        std::complex<double> psi = std::exp(std::complex<double>(0,k*dot_));
        return std::complex<double>(0,-0.5)*index*k*dot/distance*cyl_hankel_1(1,index*k*distance)*psi;
    }

    template<typename Curve>
    static std::complex<double> kernel_test(
        const double& distance, 
        const double& index,
        double k,
        double dot
    ) {
        return 1;
    }
};



} // namespace bem

#endif // BEM_GAUSS_INTEGRATOR_H