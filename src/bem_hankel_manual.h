#pragma once
#include <complex>
#include <cmath>
#include <cassert>
#include <vector>
#include <iostream>
#include "bem_constants.h"
// Euler 常数 γ
// static const double EULER_GAMMA = 0.57721566490153286060;


// 用于存储调和数的容器
std::vector<double> harmonic_cache;

// 初始化调和数缓存
void init_harmonic_cache() {
    harmonic_cache.resize(1000001);
    harmonic_cache[0] = 0.0;
    for (int i = 1; i <= 1000000; ++i) {
        harmonic_cache[i] = harmonic_cache[i - 1] + 1.0 / i;
    }
}

// 获取调和数，从缓存中读取
double harmonic_number(int m) {
    static bool initialized = false;
    if (!initialized) {
        init_harmonic_cache();
        initialized = true;
    }
    if (m < 0) {
        return 0.0;
    }
    if (m > 1000000) {
        std::cerr << "Harmonic number requested is out of pre - calculated range." << std::endl;
        return 0.0;
    }
    return harmonic_cache[m];
}


// 第一类贝塞尔函数 J0(z)
std::complex<double> BesselJ0(const std::complex<double>& z) {
    double absz = std::abs(z);
    // 小 |z| 时使用幂级数展开
    if (absz < 30.0) {
        std::complex<double> term = 1.0;
        std::complex<double> sum = term;
        for(int k = 1; k < 5000000; ++k) {
            term *= -z*z / (4.0 * k * k);
            sum += term;
            if (std::abs(term) < 1e-16) break;
        }
        return sum;
    }
    // 大 |z| 时使用渐近展开（参考 Abramowitz-Stegun 9.2.5）
    std::complex<double> w = std::sqrt(2.0/(M_PI * z));
    std::complex<double> phase = z - std::complex<double>(M_PI/4.0, 0.0);
    std::complex<double> C = std::cos(phase);
    std::complex<double> S = std::sin(phase);

    // 系数展开：添加更多项
    std::complex<double> z2 = z*z, z4 = z2*z2, z6 = z4*z2, z8 = z4*z4;
    auto Ccorr = 1.0
               -  (9.0/(128.0 * z2))
               +  (3675.0/(262144.0 * z4))
               -  (2278125.0/(4294967296.0 * z6))
               +  (202169175.0/(68719476736.0 * z8));
    auto Scorr = (1.0/(8.0 * z))
               -  (75.0/(3072.0 * z2 * z))
               +  (59535.0/(983040.0 * z4 * z))
               -  (54916875.0/(31539609600.0 * z6 * z))
               +  (563656335.0/(1009569761280.0 * z8 * z));

    return w * (C * Ccorr - S * Scorr);
}

// 第一类贝塞尔函数 J1(z)
std::complex<double> BesselJ1(const std::complex<double>& z) {
    double absz = std::abs(z);
    if (absz < 25.0) {
        std::complex<double> term = z/2.0;
        std::complex<double> sum = term;
        for(int k = 1; k < 5000000; ++k) {
            term *= -z*z / (4.0 * k * (k+1));
            sum += term;
            if (std::abs(term) < 1e-16) break;
        }
        return sum;
    }
    // 大 z 渐近展开到 O(z^-5)，添加更多项
    auto w = std::sqrt(2.0/(M_PI * z));
    auto phi = z - std::complex<double>(3.0*M_PI/4.0, 0.0);
    auto C = std::cos(phi), S = std::sin(phi);
    auto z2 = z*z, z4 = z2*z2, z6 = z4*z2, z8 = z4*z4;
    // C-part: 1 + 15/(128 z^2) - 14175/(98304 z^4) + ...
    auto Ccorr = 1.0
               + (15.0/(128.0 * z2))
               - (14175.0/(98304.0 * z4))
               + (229725.0/(3145728.0 * z6))
               - (42373725.0/(100663296.0 * z8));
    // S-part: 3/(8z) - 105/(1024 z^3) + 1091475/(3932160 z^5) + ...
    auto Scorr = (3.0/(8.0 * z))
               - (105.0/(1024.0 * z2 * z))
               + (1091475.0/(3932160.0 * z4 * z))
               - (129201375.0/(125829120.0 * z6 * z))
               + (24724348125.0/(4026531840.0 * z8 * z));
    return w * (C * Ccorr - S * Scorr);
}

// 第二类贝塞尔函数 Y0(z)
std::complex<double> BesselY0(const std::complex<double>& z) {
    double absz = std::abs(z);
    if (absz < 0.005) {
        std::complex<double> first_term = (2.0 / M_PI) * (std::log(z / 2.0) + bem::EULER_GAMMA);
        std::complex<double> sum = 0.0;
        const int MAX_TERMS = 2000; // 级数展开的最大项数
        for (int k = 1; k < MAX_TERMS; ++k) {
            std::complex<double> term = std::pow(-1, k) * std::pow(z / 2.0, 2 * k);
            double harmonic_sum = 0.0;
            for (int m = 1; m <= k; ++m) {
                harmonic_sum += 1.0 / m;
            }
            term *= harmonic_sum / (std::tgamma(k + 1) * std::tgamma(k + 1));
            sum += term;
        }
        std::complex<double> second_term = (2.0 / M_PI) * sum;
        return first_term - second_term;
    }
    else if(absz<30.0) {

        const int MAXM = 1000000;
        const double tol = 1e-14;
        const double gamma = 0.5772156649015328606;
        auto J0 = BesselJ0(z);
        auto L  = std::log(z/2.0) + gamma;
        std::complex<double> sum  = 0.0;
        std::complex<double> term = 1.0;            // m=0 时 (z/2)^0/(0!·0!) = 1
        std::complex<double> z2   = (z*z) / 4.0;
        for (int m = 1; m < MAXM; ++m) {
            term *= z2 / double(m*m);               // (z/2)^{2m}/(m!·m!)
            double sign = ((m % 2)==1 ? +1.0 : -1.0);// (-1)^{m+1}
            auto add = sign * harmonic_number(m) * term;
            sum += add;
            if (std::abs(add) < tol * std::abs(sum)) break;
        }
        return (2.0/M_PI) * ( L*J0 + sum );
    }
    auto w = std::sqrt(2.0/(M_PI * z));
    auto phi = z - std::complex<double>(M_PI/4.0, 0.0);
    auto C = std::cos(phi), S = std::sin(phi);
    auto z2 = z*z, z4 = z2*z2, z6 = z4*z2, z8 = z4*z4;
    // S-part: 1 - 9/(128 z^2) + 3675/(262144 z^4) + ...
    auto Scorr = 1.0
               - (9.0/(128.0 * z2))
               + (3675.0/(262144.0 * z4))
               - (2278125.0/(4294967296.0 * z6))
               + (202169175.0/(68719476736.0 * z8));
    // C-part: 1/(8z) - 75/(3072 z^3) + 59535/(983040 z^5) + ...
    auto Ccorr = (1.0/(8.0 * z))
               - (75.0/(3072.0 * z2 * z))
               + (59535.0/(983040.0 * z4 * z))
               - (54916875.0/(31539609600.0 * z6 * z))
               + (563656335.0/(1009569761280.0 * z8 * z));
    return w * (S * Scorr - C * Ccorr);
}

// 第二类贝塞尔函数 Y1(z)
std::complex<double> BesselY1(const std::complex<double>& z) {
    double absz = std::abs(z);
    if (absz < 0.005) {
        std::complex<double> first_term = -2.0 / (M_PI * z);
        std::complex<double> second_term = (2.0 / M_PI) * (std::log(z / 2.0) + bem::EULER_GAMMA) * (z / 2.0);
        std::complex<double> sum = 0.0;
        const int MAX_TERMS = 2000; // 级数展开的最大项数
        for (int k = 1; k < MAX_TERMS; ++k) {
            std::complex<double> term = std::pow(-1, k) * std::pow(z / 2.0, 2 * k + 1);
            double harmonic_sum1 = 0.0;
            double harmonic_sum2 = 0.0;
            for (int m = 1; m <= k; ++m) {
                harmonic_sum1 += 1.0 / m;
            }
            for (int m = 1; m <= k + 1; ++m) {
                harmonic_sum2 += 1.0 / m;
            }
            term *= (harmonic_sum1 + harmonic_sum2) / (std::tgamma(k + 1) * std::tgamma(k + 2));
            sum += term;
        }
        std::complex<double> third_term = sum / (2 * M_PI);
        return first_term + second_term + third_term;
    }
    else if(absz<25.0) {
                // 常量和收敛控制
        constexpr int    MAXM = 1000000;
        constexpr double tol  = 1e-14;
        constexpr double γ    = 0.5772156649015328606;  // 欧拉常数

        // 1) 先用已有的幂级数实现计算 J1(z)
        auto J1 = BesselJ1(z);

        // 2) 构造 L = ln(z/2) + γ
        std::complex<double> L = std::log(z/2.0) + γ;

        // 3) 计算第二项 2/(π z)
        std::complex<double> singular = 2.0 / (M_PI * z);

        // 4) 幂级数求和 ∑_{k=0}∞ (−1)^k (H_k+H_{k+1})/(k!(k+1)!) (z/2)^{2k+1}
        std::complex<double> sum2 = 0.0;
        std::complex<double> term = z / 2.0;               // k=0 项: (z/2)^{1}/(0!·1!)
        auto                  z2   = (z * z) / 4.0;       // (z/2)^2

        for (int k = 0; k < MAXM; ++k) {
            // 调和数 H_k, H_{k+1}
            double Hk   = harmonic_number(k);
            double Hkp1 = harmonic_number(k+1);
            // (-1)^k
            double sign = (k % 2 == 0 ? +1.0 : -1.0);

            // 累加
            sum2 += sign * (Hk + Hkp1) * term;

            // 更新 term -> (z/2)^{2(k+1)+1}/((k+1)!(k+2)!)
            term *= z2 / double((k+1) * (k+2));

            // 收敛判断
            if (std::abs(term) < tol * std::abs(sum2)) break;
        }

        // 5) 最终组合
        //    Y1(z) = 2/π [ L·J1(z) ] - singular - (1/π)·sum2
        return std::complex<double>(2.0/M_PI) * (L * J1)
            - singular
            - sum2 * std::complex<double>(1.0/M_PI);
    }
    auto w = std::sqrt(2.0/(M_PI * z));
    auto phi = z - std::complex<double>(3.0*M_PI/4.0, 0.0);
    auto C = std::cos(phi), S = std::sin(phi);
    auto z2 = z*z, z4 = z2*z2, z6 = z4*z2, z8 = z4*z4;
    // S-part: 1 + 15/(128 z^2) - 14175/(98304 z^4) + ...
    auto Scorr = 1.0
               + (15.0/(128.0 * z2))
               - (14175.0/(98304.0 * z4))
               + (229725.0/(3145728.0 * z6))
               - (42373725.0/(100663296.0 * z8));
    // C-part: 3/(8z) - 105/(1024 z^3) +1091475/(3932160 z^5) + ...
    auto Ccorr = (3.0/(8.0 * z))
               - (105.0/(1024.0 * z2 * z))
               + (1091475.0/(3932160.0 * z4 * z))
               - (129201375.0/(125829120.0 * z6 * z))
               + (24724348125.0/(4026531840.0 * z8 * z));
    return w * (S * Scorr + C * Ccorr);
}

// 第一类汉克尔函数 H0^(1)(z) = J0(z) + i·Y0(z)
std::complex<double> HankelH0(const std::complex<double>& z) {
    return BesselJ0(z) + std::complex<double>(0.0, 1.0) * BesselY0(z);
}

// 第一类汉克尔函数 H1^(1)(z) = J1(z) + i·Y1(z)
std::complex<double> HankelH1(const std::complex<double>& z) {
    return BesselJ1(z) + std::complex<double>(0.0, 1.0) * BesselY1(z);
}

std::complex<double> complex_hankel_1(const int& nu, const std::complex<double>& z) {
    if(nu==1) return BesselJ1(z) + std::complex<double>(0.0, 1.0) * BesselY1(z);
    if(nu==0) return BesselJ0(z) + std::complex<double>(0.0, 1.0) * BesselY0(z);

}    