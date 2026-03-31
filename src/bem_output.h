// 文件名：bem_output.h
#ifndef BEM_OUTPUT_H
#define BEM_OUTPUT_H

#include "bem_scattering.h"
#include <fstream>
#include <iomanip>

namespace bem {

/**
 * 将总散射截面 σ(k) 输出到 cross_section.txt
 * 格式：每行 "k σ"
 */
inline void write_cross_section(
    const std::vector<double>& k_vals,
    const std::vector<double>& sigma_vals,
    const std::string& filename = "cross_section.txt"
) {
    std::ofstream ofs(filename);
    ofs << std::setprecision(8);
    for (size_t i = 0; i < k_vals.size(); ++i) {
        ofs << k_vals[i] << " " << sigma_vals[i] << "\n";
    }
}

/**
 * 将近场 ψ(r) 输出到 near_field.txt
 * 第一行：Nx Ny （网格尺寸）
 * 其后每行： x y Re(psi) Im(psi)
 */
inline void write_near_field(
    const std::vector<std::pair<Point,std::complex<double>>>& grid, 
    int Nx, int Ny,
    const std::string& filename = "near_field.txt"
) {
    std::ofstream ofs(filename);
    ofs << Nx << " " << Ny << "\n";
    ofs << std::setprecision(8);
    for (auto& p : grid) {
        const auto& pt = p.first;
        auto psi = p.second;
        ofs << pt.x_ << " " << pt.y_ << " "
            << psi.real() << " " << psi.imag() << "\n";
    }
}

/**
 * 将远场 f(θ) 输出到 far_field.txt
 * 格式：每行 "theta Re(f) Im(f)"
 */
inline void write_far_field(
    const std::vector<double>& thetas,
    const std::vector<std::complex<double>>& f_vals,
    const std::string& filename = "far_field.txt"
) {
    std::ofstream ofs(filename);
    ofs << std::setprecision(8);
    for (size_t i = 0; i < thetas.size(); ++i) {
        ofs << thetas[i] << " "
            << f_vals[i].real() << " " << f_vals[i].imag() << "\n";
    }
}

// 导出两个腔的 BEMElement.mid
void save_midpoints_to_txt(
    const std::vector<Microcavity::BEMElement>& elems1,
    const Point& off1,
    const std::vector<Microcavity::BEMElement>& elems2,
    const Point& off2,
    const std::string& filename
) {
    std::ofstream ofs(filename);
    // 第一腔
    for (auto &e : elems1) {
        Point p = e.mid + off1;
        ofs << p.x_ << "," << p.y_ << "\n";
    }
    ofs << "\n";  // 空行分隔
    // 第二腔
    for (auto &e : elems2) {
        Point p = e.mid + off2;
        ofs << p.x_ << "," << p.y_ << "\n";
    }
    ofs.close();
}

}
#endif // BEM_OUTPUT_H
