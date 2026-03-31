// 文件名：bem_microcavity.h
#ifndef BEM_MICROCAVITY_H
#define BEM_MICROCAVITY_H

#include "bem_curve_segment.h"
#include <cmath>
#include <memory>
#include <vector>
#include <stdexcept>
#include <algorithm>

namespace bem {

// 判断两个点是否接近相等（防止浮点误差）
inline bool point_equal(const Point& a, const Point& b, double tol = EPS) {
    return (a - b).norm() < tol;
}

class Microcavity {
public:
    using SegmentPtr = std::shared_ptr<CurveSegment>;
    using SegmentList = std::vector<SegmentPtr>;
    using PointList = std::vector<Point>;
    using BoundaryCondition = enum { Dirichlet, Neumann, Robin };

    // 构造/析构
    Microcavity(double index = 1.466) : index_(index) {
        check_closed_ = true; // 默认开启封闭性检查
    }

    const std::vector<std::shared_ptr<CurveSegment>>& segments() const {
        return segments_;
    }

    void add_segment(SegmentPtr seg, BoundaryCondition bc = Neumann) {
        if (segments_.empty()) {
            start_point_ = seg->start();
        } else {
            auto& last_seg = segments_.back();
            if (!point_equal(last_seg->end(), seg->start()) && check_closed_) {
                throw std::runtime_error("Curve segments are not connected");
            }
        }
        segments_.push_back(seg);
        bcs_.push_back(bc);
    }
    
    void close() {
        if (segments_.empty()) {
            is_closed_ = false;
            return;
          }
        auto first = segments_.front();
        auto last  = segments_.back();
        if (check_closed_ && !point_equal(first->start(), last->end())) {
            is_closed_ = false;
            throw std::runtime_error("Microcavity boundary is not closed");
          } else {
            is_closed_ = true;
          }
    }
    

    // 几何属性
    double perimeter() const {
        double len = 0.0;
        for(const auto& seg : segments_) {
            len += seg->arc_length();
        }
        return len;
    }

    // 边界元离散（生成边界元所需的单元和节点）
    struct BEMElement {
        Point p1, p2;        // 单元端点（物理坐标）
        double seg_start_t, seg_end_t; // 在原始CurveSegment上的参数范围（0≤t≤1）
        Point mid;           // 中点（用于源点近似）
        double length;       // 单元长度（弧长）
        BoundaryCondition bc; // 边界条件
        const CurveSegment* seg_ptr; // 指向原始曲线段的指针
        double curvature;    // 中点曲率（新增，用于对角元素计算）
    };


    // 对称性处理（论文4.1节镜像格林函数）
    void set_symmetry(double angle, bool mirror = true) {
        symmetry_angle_ = angle;
        is_mirror_ = mirror;
    }

    // 辅助检查
    bool is_closed() const { return is_closed_; }
    void enable_close_check(bool enable) { check_closed_ = enable; }

    std::vector<BEMElement> discretize(int n_line_per_edge, int n_arc_per_corner,double ds_ref,
    double rho_ds_ratio) const;
    SegmentList segments_;       // 组成微腔的曲线段
    void move(const Point& p) {
        start_point_+=p;
        for(auto seg: segments_) {
            seg->move(p);
        }
    }
private:
    // SegmentList segments_;       // 组成微腔的曲线段
    std::vector<BoundaryCondition> bcs_; // 各段边界条件
    Point start_point_;          // 起始点（用于封闭检查）
    double index_;               // 微腔折射率（默认1.466）
    bool is_closed_ = false;      // 是否已闭合
    bool check_closed_ = true;    // 是否开启封闭性检查
    double symmetry_angle_ = 0.0; // 对称轴角度（弧度）
    bool is_mirror_ = false;      // 是否镜像对称
};

inline std::vector<Microcavity::BEMElement>
Microcavity::discretize(
    int n_line_per_edge,
    int n_arc_per_corner,
    double ds_ref,
    double rho_ds_ratio
) const {
    std::vector<BEMElement> elems;
    // 先粗略估算总元数
    double perim = perimeter();
    int total_est = std::max(1, int(std::ceil(perim / ds_ref)));
    elems.reserve(total_est);

    for (size_t si = 0; si < segments_.size(); ++si) {
        auto seg = segments_[si];
        // 直线段
        if (auto line = std::dynamic_pointer_cast<LineSegment>(seg)) {
            double L = line->arc_length();
            int m = n_line_per_edge;
            double dt = 1.0 / m;
            double t0 = 0.0;
            for (int k = 0; k < m; ++k) {
                double t1 = t0 + dt;
                double tm = 0.5 * (t0 + t1);
                Point p0 = (*line)(t0), p1 = (*line)(t1), pm = (*line)(tm);
                double len = p0.distance(p1);
                elems.push_back(Microcavity::BEMElement{p0, p1, t0, t1, pm, len,
                                 bcs_[si], seg.get(),
                                 line->curvature(tm)});
                t0 = t1;
            }
        }
        // 圆弧段
        else if (auto arc = std::dynamic_pointer_cast<ArcSegment>(seg)) {
            double L = arc->arc_length();
            int m = n_arc_per_corner > 0 ? n_arc_per_corner : int(std::ceil(L / ds_ref));
            double dt = 1.0 / m;
            double t0 = 0.0;
            for (int k = 0; k < m; ++k) {
                double t1 = t0 + dt;
                double tm = 0.5 * (t0 + t1);
                Point p0 = (*arc)(t0), p1 = (*arc)(t1), pm = (*arc)(tm);
                double len = p0.distance(p1);
                elems.push_back(Microcavity::BEMElement{p0, p1, t0, t1, pm, len,
                                 bcs_[si], seg.get(),
                                 arc->curvature(tm)});
                t0 = t1;
            }
        }
        // 贝塞尔倒角段
        else if (auto bez = std::dynamic_pointer_cast<FilletBezierSegment>(seg)) {
            double t0 = 0.0;
            while (t0 < 1.0) {
                // 1) 先用 ds_ref 在 t0 处估算一个初步 dt
                double speed0 = bez->tangent(t0).norm();
                double dt_ref = ds_ref / speed0;
                double t1_ref = std::min(1.0, t0 + dt_ref);
        
                // 2) 在 [t0, t1_ref] 的端点和中点取样，找出最大曲率
                double tm_ref = 0.5*(t0 + t1_ref);
                double κ0 = std::abs(bez->curvature(t0));
                double κ1 = std::abs(bez->curvature(t1_ref));
                double κm = std::abs(bez->curvature(tm_ref));
                double κ_max = std::max({κ0, κ1, κm});
                double ρ_min = (κ_max>0.0 ? 1.0/κ_max : std::numeric_limits<double>::infinity());
        
                // 3) 允许的弧长：既要 <= ds_ref，又要满足 ρ_min/Δs >= rho_ds_ratio
                double ds_allowed = std::min(ds_ref, ρ_min / rho_ds_ratio);
        
                // 4) 用中点速度更准确地反推出 dt
                double speed_m = bez->tangent(tm_ref).norm();
                double dt = ds_allowed / speed_m;
                if (dt <= 1e-6) dt = 1e-6;
        
                // 5) 最终截断到 [0,1]
                double t1 = std::min(1.0, t0 + dt);
                double tm = 0.5*(t0 + t1);
        
                // 6) 生成这一小段
                Point p0 = (*bez)(t0), p1 = (*bez)(t1), pm = (*bez)(tm);
                double len = p0.distance(p1);
                elems.push_back(Microcavity::BEMElement{
                  p0, p1, t0, t1, pm, len,
                  bcs_[si], seg.get(),
                  bez->curvature(tm)
                });
        
                if (t1 >= 1.0) break;
                t0 = t1;
            }
        }else if (auto bez = dynamic_cast<QuinticHermiteSegment*>(seg.get())) {
            double t0 = 0.0;
            while (t0 < 1.0) {
                // 1) 先用 ds_ref 在 t0 处估算一个初步 dt
                double speed0 = bez->tangent(t0).norm();
                double dt_ref = ds_ref / speed0;
                double t1_ref = std::min(1.0, t0 + dt_ref);
        
                // 2) 在 [t0, t1_ref] 的端点和中点取样，找出最大曲率
                double tm_ref = 0.5*(t0 + t1_ref);
                double κ0 = std::abs(bez->curvature(t0));
                double κ1 = std::abs(bez->curvature(t1_ref));
                double κm = std::abs(bez->curvature(tm_ref));
                double κ_max = std::max({κ0, κ1, κm});
                double ρ_min = (κ_max>0.0 ? 1.0/κ_max : std::numeric_limits<double>::infinity());
        
                // 3) 允许的弧长：既要 <= ds_ref，又要满足 ρ_min/Δs >= rho_ds_ratio
                double ds_allowed = std::min(ds_ref, ρ_min / rho_ds_ratio);
        
                // 4) 用中点速度更准确地反推出 dt
                double speed_m = bez->tangent(tm_ref).norm();
                double dt = ds_allowed / speed_m;
                if (dt <= 1e-6) dt = 1e-6;
        
                // 5) 最终截断到 [0,1]
                double t1 = std::min(1.0, t0 + dt);
                double tm = 0.5*(t0 + t1);
        
                // 6) 生成这一小段
                Point p0 = (*bez)(t0), p1 = (*bez)(t1), pm = (*bez)(tm);
                double len = p0.distance(p1);
                elems.push_back(Microcavity::BEMElement{
                  p0, p1, t0, t1, pm, len,
                  bcs_[si], seg.get(),
                  bez->curvature(tm)
                });
        
                if (t1 >= 1.0) break;
                t0 = t1;
            }
        }
        else {
            throw std::runtime_error("Unknown segment type in discretize");
        }
    }
    return elems;
}


/// @brief 构造带圆角（fillet）的正六边形微腔
/// @param R             正六边形的边长（同时也是顶点到中心的距离）
/// @param corner_radius 圆角半径（fillet radius）
/// @param index         折射率
inline Microcavity make_rounded_hexagonal_cavity(
    double R,
    double corner_radius,
    double index    = 1.466
) {
    Microcavity cav(index);

    const int N = 6;
    // 1) 先算出规则六边形的 6 个顶点（逆时针）
    std::vector<Point> P(N);
    std::vector<double> theta(N);
    double start_ang = M_PI/3.0;
    for (int i = 0; i < N; ++i) {
        double θ = start_ang + M_PI/3.0*i;
        theta[i]=θ;
        P[i] = Point::polar(R, θ);
    }

    // 2) 计算每个顶点处圆弧的起点 S[i] 和终点 E[i]
    std::vector<Point> S(N), E(N);
    for (int i = 0; i < N; ++i) {
        S[i] = P[i]+Point::polar(corner_radius, theta[i]-M_PI/6.0);
        E[i] = P[i]+Point::polar(corner_radius, theta[i]+M_PI/6.0);
    }

    // 3) 按顺序插入：先直线，再圆弧
    for (int i = 0; i < N; ++i) {
        int j = (i+N-1)%N;  // 上一个圆角的索引
        // 3.1) 直线： E[j] --> S[i]
        {
            auto line = std::make_shared<LineSegment>(E[j], S[i]);
            cav.add_segment(line, Microcavity::Neumann);
        }
        // 3.2) 圆弧：S[i] --> E[i]，圆心是 P[i]
        {
            auto arc = std::make_shared<ArcSegment>(
                P[i],        // 圆心
                theta[i]-M_PI/6.0,
                theta[i]+M_PI/6.0,
                corner_radius,
                true         // ccw 方向保持逆时针
            );
            cav.add_segment(arc, Microcavity::Neumann);
        }
    }

    cav.close();
    return cav;
}
inline Microcavity make_bezier_rounded_hexagonal_cavity(
    double R,
    double corner_radius,
    double index = 1.466
) {
    Microcavity cav(index);
    const int N = 6;
    std::vector<Point> P(N), S(N), E(N);
    std::vector<double> theta(N);

    // 1) 计算正六边形 6 个顶点
    double start_ang = M_PI / 3.0;
    for (int i = 0; i < N; ++i) {
        theta[i] = start_ang + i * M_PI/3.0;
        P[i]     = Point::polar(R, theta[i]);
    }

    // 2) 在每条原始边上，距离顶点 corner_radius 处确定倒角起点 S[i] 和终点 E[i]
    for (int i = 0; i < N; ++i) {
        int j = (i + N - 1) % N;  // 上一个顶点
        int k = (i + 1) % N;      // 下一个顶点
        Point v_prev = (P[j] - P[i]).normalized();  // 由 P[i] 指向 P[j]
        Point v_next = (P[k] - P[i]).normalized();  // 由 P[i] 指向 P[k]
        S[i] = P[i] + v_prev * corner_radius;       // 沿前一条边回退
        E[i] = P[i] + v_next * corner_radius;       // 沿下一条边前进
    }

    // 3) 预先构造所有 Bézier 圆角段
    //    选 L 使得 t=0.5 处曲率 = 1/corner_radius
    const double C = (8.0 + 6.0*std::sqrt(3.0)
                      - 4.0*std::sqrt(4.0 + 6.0*std::sqrt(3.0)))
                     / 3.0;
    double L = C * corner_radius;

    std::vector<std::shared_ptr<FilletBezierSegment>> beziers(N);
    for (int i = 0; i < N; ++i) {
        int j = (i + N - 1) % N;
        int k = (i + 1) % N;
        Point v0 = (P[i] - P[j]).normalized();  // 从 S[i] 进入圆角的切线方向
        Point v3 = (P[k] - P[i]).normalized();  // 从圆角退出到 E[i] 的切线方向
        beziers[i] = std::make_shared<FilletBezierSegment>(
            S[i], E[i], v0, v3, L
        );
    }

    // 4) 按顺序插入：先直线段，再 Bézier 圆角
    for (int i = 0; i < N; ++i) {
        int j = (i + N - 1) % N;
        // 4.1) 直线段：从上一个圆角的结束 E[j] 到本圆角的开始 S[i]
        cav.add_segment(
            std::make_shared<LineSegment>( E[j], S[i] ),
            Microcavity::Neumann
        );
        // 4.2) Bézier 圆角段
        cav.add_segment(
            beziers[i],
            Microcavity::Neumann
        );
    }

    cav.close();
    return cav;
}


// 计算圆角切向量长度 L0
double solve_hermite_L0(double corner_radius,
    const Point& S, const Point& E,
    const Point& dir_in, const Point& dir_out)
{
// Hermite 基函数在 t=0.5 处的导数系数（常量）
const double t = 0.5;
const double dH0 = -1.875;    // dH0/dt |_{t=0.5}
const double dH1 = -0.4375;   // dH1/dt |_{t=0.5}
const double dH3 =  1.875;    // dH3/dt |_{t=0.5}
const double dH4 = -0.4375;   // dH4/dt |_{t=0.5}
const double ddH1 = -1.5;     // d^2H1/dt^2 |_{t=0.5}
const double ddH4 =  3.375;   // d^2H4/dt^2 |_{t=0.5}

// 计算给定 L 时曲率与目标曲率之差
auto curvature_diff = [&](double L) {
// 计算切向量 V0, V1
double V0x = dir_in.x_ * L;
double V0y = dir_in.y_ * L;
double V1x = dir_out.x_ * L;
double V1y = dir_out.y_ * L;
// 一阶导 P'(0.5) = S*dH0 + V0*dH1 + E*dH3 + V1*dH4
double dx = S.x_ * dH0 + V0x * dH1 + E.x_ * dH3 + V1x * dH4;
double dy = S.y_ * dH0 + V0y * dH1 + E.y_ * dH3 + V1y * dH4;
// 二阶导 P''(0.5) = S*0 + V0*ddH1 + E*0 + V1*ddH4
double ddx = V0x * ddH1 + V1x * ddH4;
double ddy = V0y * ddH1 + V1y * ddH4;
// 计算曲率
double num = dx * ddy - dy * ddx;  // 分子（向量叉积）
double denom = pow(dx*dx + dy*dy, 1.5);
if (denom == 0.0) return 1e30;    // 防止除零
double curvature = fabs(num) / denom;
return curvature - (1.0 / corner_radius);
};

// 二分法查找根（寻找小的正 L 解）
double Llo = 0.0;
double Lhi = std::max(1e-3, std::hypot(E.x_-S.x_, E.y_-S.y_)); // 初始上界
// 调整上界使 f(Lhi) >= 0（曲率大于等于目标）
while (curvature_diff(Lhi) < 0.0) {
Lhi *= 2.0;
if (Lhi > 1e3) break; // 避免无限扩张
}
// 二分迭代
double Lmid = 0.0;
for (int iter = 0; iter < 60; ++iter) {
Lmid = 0.5 * (Llo + Lhi);
double fmid = curvature_diff(Lmid);
if (fmid > 0) {
Lhi = Lmid;
} else {
Llo = Lmid;
}
}
return 0.5*(Llo + Lhi);
}

/// @brief 正六边形上的 Hermite 圆角微腔（G² 端点平滑）
/// @param R    原六边形边长
/// @param index 折射率
inline Microcavity make_hermite_rounded_hexagonal_cavity(
    double R, double corner_radius, double index=1.466
) {
    Microcavity cav(index);
    const int N = 6;
    std::vector<Point> P(N), S(N), E(N);
    std::vector<double> θ(N);
    double start_ang = M_PI/3.0;
    // 1) 三角形顶点
    for(int i=0;i<N;++i){
        θ[i] = start_ang + i*M_PI/3.0;
        P[i] = Point::polar(R, θ[i]);
    }
    // 2) 倒角起点/终点
    for(int i=0;i<N;++i){
        int j=(i+N-1)%N, k=(i+1)%N;
        Point v_prev = (P[j]-P[i]).normalized();
        Point v_next = (P[k]-P[i]).normalized();
        S[i] = P[i] + v_prev*R/800;
        E[i] = P[i] + v_next*R/800;
    }
    // 3) Hermite 圆角段：端点曲率=0，切向长度 = R/50
    double L0 = solve_hermite_L0(corner_radius, S[0], E[0], (P[0] - P[N-1]).normalized(), (P[1] - P[0]).normalized());
    std::vector<std::shared_ptr<QuinticHermiteSegment>> segs(N);
    for(int i=0;i<N;++i){
        int j=(i+N-1)%N, k=(i+1)%N;
        // 切向向量
        Point m0 = (P[i] - P[j]).normalized()*L0;   // 进
        Point m1 = (P[k] - P[i]).normalized()*L0; // 出 （方向取反）
        // P''(0)=P''(1)=0 保证端点曲率为 0
        segs[i] = std::make_shared<QuinticHermiteSegment>(
            S[i], E[i], m0, m1,
            Point(0,0), Point(0,0)
        );
    }
    // 4) 拼装：Line + Hermite + Line + Hermite …
    for(int i=0;i<N;++i){
        int j=(i+N-1)%N;
        cav.add_segment(
            std::make_shared<LineSegment>( E[j], S[i] ),
            Microcavity::Neumann
        );
        cav.add_segment(segs[i], Microcavity::Neumann);
    }
    cav.close();
    return cav;
}

/// 判断点 q 是否在由有序顶点 poly 构成的多边形内部（射线法）
inline bool point_in_polygon(const std::vector<Point>& poly, const Point& q) {
    bool inside = false;
    int n = (int)poly.size();
    for (int i = 0, j = n-1; i < n; j = i++) {
        const Point &pi = poly[i], &pj = poly[j];
        // 检查边 (pj->pi) 与水平射线 y = q.y 相交
        bool intersect = ((pi.y_ > q.y_) != (pj.y_ > q.y_))
                      && (q.x_ < (pj.x_-pi.x_)*(q.y_-pi.y_)/(pj.y_-pi.y_) + pi.x_);
        if (intersect) inside = !inside;
    }
    return inside;
}

} // namespace bem

#endif // BEM_MICROCAVITY_H