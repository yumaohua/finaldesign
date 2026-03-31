#include <vector>
#include <functional>
#include <cmath>
#include <stdexcept>
#include "../include/eigen-3.4.0/Eigen/Dense"
#include "../include/complex_bessel/include/complex_bessel.h"

using namespace Eigen;
using Complex = std::complex<double>;
constexpr double PI = 3.14159265358979323846;
constexpr double EPSILON = 1e-6;

// 基础数据结构
struct Point {
    double x, y;
    Point(double x_ = 0, double y_ = 0) : x(x_), y(y_) {}
    double norm() const { return std::hypot(x, y); }

    Point operator-(const Point& other) const {
        return Point(x - other.x, y - other.y);
    }
    Point operator+(const Point& other) const {
        return Point(x + other.x, y + other.y);
    }

    Point operator*(double scalar) const {
        return Point(x * scalar, y * scalar);
    }
};

// 标量乘法：double * Point
Point operator*(double scalar, const Point& p) {
    return p * scalar;
}

// 参数化曲线段类
class ParametricSegment {
private:
    std::function<Point(double)> r_func_;
    std::function<Point(double)> dr_dt_func_;
    double t_start_, t_end_;
    mutable double cached_length_ = -1;

public:
    ParametricSegment(
        std::function<Point(double)> r_func,
        std::function<Point(double)> dr_dt_func,
        double t_start, double t_end
    ) : r_func_(r_func), dr_dt_func_(dr_dt_func), t_start_(t_start), t_end_(t_end) {
        if (t_end_ <= t_start_) {
            throw std::invalid_argument("[ParametricSegment] t_end must be > t_start");
        }
    }

    Point eval(double t) const { return r_func_(t); }

    Point derivative(double t) const { return dr_dt_func_(t); }

    double length() const {
        if (cached_length_ < 0) {
            auto f = [&](double t) { return dr_dt_func_(t).norm(); };
            cached_length_ = GaussIntegrateReal(f, t_start_, t_end_, 8);
        }
        return cached_length_;
    }

    std::vector<Point> discretize(int N) const {
        std::vector<Point> points;
        double dt = (t_end_ - t_start_) / (N - 1);
        for (int i = 0; i < N; ++i) {
            double t = t_start_ + i * dt;
            points.push_back(eval(t));
        }
        return points;
    }

    std::pair<double, double> parameter_range() const { return {t_start_, t_end_}; }

private:
    // 辅助函数：实数高斯积分
    static double GaussIntegrateReal(
        std::function<double(double)> f,
        double a, double b, int order
    ) {
        // 简单的梯形积分近似，实际应使用高斯积分公式
        double h = (b - a) / (order - 1);
        double sum = 0.5 * (f(a) + f(b));
        for (int i = 1; i < order - 1; ++i) {
            sum += f(a + i * h);
        }
        return sum * h;
    }
};

// 微腔类（管理闭合边界）
class Microcavity {
private:
    std::vector<ParametricSegment> segments_;
    std::vector<Point> boundary_points_;
    std::vector<Point> normals_;
    std::vector<double> curvatures_;
    double max_segment_length_ = 0;

public:
    void add_segment(const ParametricSegment& seg) {
        if (!segments_.empty()) {
            Point prev_end = segments_.back().eval(segments_.back().parameter_range().second);
            Point new_start = seg.eval(seg.parameter_range().first);
            if ((prev_end - new_start).norm() > EPSILON) {
                throw std::runtime_error("[Microcavity] Segments are not connected!");
            }
        }
        segments_.push_back(seg);
        max_segment_length_ = std::max(max_segment_length_, seg.length());
    }

    const ParametricSegment& get_segment(int index) const {
        if (index < 0 || index >= segments_.size()) {
            throw std::out_of_range("Invalid segment index");
        }
        return segments_[index];
    }

    void close() {
        if (segments_.empty()) return;
        Point first_point = segments_.front().eval(segments_.front().parameter_range().first);
        Point last_point = segments_.back().eval(segments_.back().parameter_range().second);
        if ((first_point - last_point).norm() > EPSILON) {
            auto close_func = [&](double t) {
                return last_point + t * (first_point - last_point);
            };
            auto close_deriv = [&](double t) {
                return first_point - last_point;
            };
            segments_.emplace_back(close_func, close_deriv, 0, 1);
        }
    }

    void generate_boundary(int points_per_segment) {
        boundary_points_.clear();
        for (const auto& seg : segments_) {
            auto points = seg.discretize(points_per_segment);
            boundary_points_.insert(boundary_points_.end(), points.begin(), points.end());
        }
        calculate_normals();
        calculate_curvatures();
    }

    const std::vector<Point>& boundary_points() const { return boundary_points_; }

    const std::vector<Point>& normals() const { return normals_; }

    const std::vector<double>& curvatures() const { return curvatures_; }

    double max_segment_length() const { return max_segment_length_; }

private:
    void calculate_normals() {
        normals_.resize(boundary_points_.size());
        for (size_t i = 0; i < boundary_points_.size(); ++i) {
            Point prev = boundary_points_[(i - 1 + boundary_points_.size()) % boundary_points_.size()];
            Point curr = boundary_points_[i];
            Point next = boundary_points_[(i + 1) % boundary_points_.size()];
            Point tangent = {next.x - prev.x, next.y - prev.y};
            Point normal = {-tangent.y, tangent.x};
            double len = normal.norm();
            normals_[i] = {normal.x / len, normal.y / len};
        }
    }

    void calculate_curvatures() {
        curvatures_.resize(boundary_points_.size());
        for (size_t i = 0; i < boundary_points_.size(); ++i) {
            Point prev = boundary_points_[(i - 1 + boundary_points_.size()) % boundary_points_.size()];
            Point curr = boundary_points_[i];
            Point next = boundary_points_[(i + 1) % boundary_points_.size()];
            double dx1 = curr.x - prev.x, dy1 = curr.y - prev.y;
            double dx2 = next.x - curr.x, dy2 = next.y - curr.y;
            double cross = dx1 * dy2 - dx2 * dy1;
            double denominator = std::pow(dx1 * dx1 + dy1 * dy1 + dx2 * dx2 + dy2 * dy2, 1.5);
            curvatures_[i] = 2 * cross / denominator;
        }
    }
};

// 高斯积分工具类（动态阶数）
class GaussIntegrator {
public:
    static int select_order(double distance, double max_length) {
        double ratio = distance / max_length;
        if (ratio < 0.1)      return 8;
        else if (ratio < 0.3) return 6;
        else                  return 4;
    }

    static Complex integrate(
        std::function<Complex(Point)> f,
        const ParametricSegment& seg,
        int order
    ) {
        const auto& nodes = get_gauss_nodes(order);
        double a = seg.parameter_range().first;
        double b = seg.parameter_range().second;
        Complex sum = 0;
        double scale = (b - a) / 2.0;
        double shift = (a + b) / 2.0;
        for (const auto& [x, w] : nodes) {
            double t = scale * x + shift;
            Point r = seg.eval(t);
            Point dr_dt = seg.derivative(t);
            double ds_dt = dr_dt.norm();
            sum += w * f(r) * ds_dt;
        }
        return sum * scale;
    }

private:
    static const std::vector<std::pair<double, double>>& get_gauss_nodes(int order) {
        static const std::vector<std::pair<double, double>> nodes_4 = {
            {-std::sqrt(1.0 / 3.0), 1.0},
            {std::sqrt(1.0 / 3.0), 1.0}
        };
        static const std::vector<std::pair<double, double>> nodes_6 = {
            {-std::sqrt(3.0 / 5.0), 5.0 / 9.0},
            {0, 8.0 / 9.0},
            {std::sqrt(3.0 / 5.0), 5.0 / 9.0}
        };
        static const std::vector<std::pair<double, double>> nodes_8 = {
            {-std::sqrt(7.0 + 2.0 * std::sqrt(14.0 / 15.0)) / 3.0, (18.0 + std::sqrt(30.0)) / 36.0},
            {-std::sqrt(7.0 - 2.0 * std::sqrt(14.0 / 15.0)) / 3.0, (18.0 - std::sqrt(30.0)) / 36.0},
            {std::sqrt(7.0 - 2.0 * std::sqrt(14.0 / 15.0)) / 3.0, (18.0 - std::sqrt(30.0)) / 36.0},
            {std::sqrt(7.0 + 2.0 * std::sqrt(14.0 / 15.0)) / 3.0, (18.0 + std::sqrt(30.0)) / 36.0}
        };
        if (order == 4) return nodes_4;
        else if (order == 6) return nodes_6;
        else if (order == 8) return nodes_8;
        else throw std::invalid_argument("Unsupported Gauss integration order");
    }
};

// 边界元矩阵组装类
class BEMBuilder {
public:
    static MatrixXcd build_matrix(
        const Microcavity& cavity,
        double k,
        double n_in,
        double n_out,
        int points_per_segment
    ) {
        const auto& points = cavity.boundary_points();
        const auto& normals = cavity.normals();
        const auto& curvatures = cavity.curvatures();
        int N = points.size();
        MatrixXcd B = MatrixXcd::Zero(N, N);
        MatrixXcd C = MatrixXcd::Zero(N, N);
        double max_length = cavity.max_segment_length();

        for (int i = 0; i < N; ++i) {
            Point r_i = points[i];
            Point normal_i = normals[i];
            for (int l = 0; l < N; ++l) {
                Point r_l = points[l];
                double distance = (r_i - r_l).norm();
                int order = GaussIntegrator::select_order(distance, max_length);

                if (i == l) {
                    // 对角项（公式 32, 34）
                    double ds = cavity.get_segment(l).length() / (points_per_segment - 1);
                    B(i, l) = -1.0 + (curvatures[i] * ds) / (2 * PI);
                    C(i, l) = (ds / PI) * (
                        1.0 - std::log(n_in * std::abs(k) * ds / 4.0) +
                        Complex(0, PI / 2) - 0.5772156649
                    );
                } else {
                    // 非对角项（动态积分）
                    auto f_B = [&](Point r) -> Complex {
                        return GreenDerivative(r, r_i, normal_i, k, n_out);
                    };
                    auto f_C = [&](Point r) -> Complex {
                        return GreenFunction(r, r_i, k, n_out);
                    };
                    B(i, l) = GaussIntegrator::integrate(f_B, cavity.get_segment(l), order);
                    C(i, l) = GaussIntegrator::integrate(f_C, cavity.get_segment(l), order);
                }
            }
        }

        MatrixXcd M(2 * N, 2 * N);
        M << B, C,
             MatrixXcd::Identity(N, N), MatrixXcd::Zero(N, N);
        return M;
    }

private:
    static Complex GreenFunction(const Point& r, const Point& r_prime, double k, double n) {
        double rho = (r - r_prime).norm();
        return -Complex(0, 0.25) * sp_bessel::hankelH1(0, n * k * rho);
    }

    static Complex GreenDerivative(const Point& r, const Point& r_prime,
                                  const Point& normal, double k, double n) {
        double dx = r.x - r_prime.x;
        double dy = r.y - r_prime.y;
        double rho = std::hypot(dx, dy);
        if (rho < EPSILON) return Complex(0, 0);
        Complex h1 = sp_bessel::hankelH1(1, n * k * rho);
        double cos_theta = (dx * normal.x + dy * normal.y) / rho;
        return -Complex(0, 0.25) * n * k * h1 * cos_theta;
    }
};

// 共振波数计算相关函数（简单示例，未完全实现论文中的牛顿法等复杂逻辑）
Complex find_resonant_wavenumber(const Microcavity& cavity, double initial_k, double n_in, double n_out, int points_per_segment) {
    MatrixXcd M = BEMBuilder::build_matrix(cavity, initial_k, n_in, n_out, points_per_segment);
    // 简单的行列式计算示例，实际应使用更高效的方法
    Complex det_M = M.determinant();
    // 这里只是简单返回，实际应使用牛顿法迭代
    return initial_k - det_M / (1e-6 + det_M);
}

int main() {
    // 示例：构建一个简单的圆形微腔
    auto circle_r = [](double t) {
        return Point(1.0 * std::cos(t), 1.0 * std::sin(t));
    };
    auto circle_dr_dt = [](double t) {
        return Point(-1.0 * std::sin(t), 1.0 * std::cos(t));
    };
    ParametricSegment circle_segment(circle_r, circle_dr_dt, 0, 2 * PI);

    Microcavity cavity;
    cavity.add_segment(circle_segment);
    cavity.close();
    int points_per_segment = 100;
    cavity.generate_boundary(points_per_segment);

    double k = 1.0;
    double n_in = 1.5;
    double n_out = 1.0;
    MatrixXcd M = BEMBuilder::build_matrix(cavity, k, n_in, n_out, points_per_segment);
    Complex resonant_k = find_resonant_wavenumber(cavity, k, n_in, n_out, points_per_segment);
    std::cout << "Resonant wavenumber: " << resonant_k << std::endl;

    return 0;
}