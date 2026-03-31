#include <vector>
#include <functional>
#include <cmath>
#include <stdexcept>
#include "../include/eigen-3.4.0/Eigen/Dense"
#include "../include/complex_bessel/include/complex_bessel.h"

using namespace Eigen;
using Complex = std::complex<double>;
constexpr double PI = 3.14159265358979323846;
constexpr double EPSILON = 1e-8;

struct Point {
    double x, y;
    Point(double x_ = 0, double y_ = 0) : x(x_), y(y_) {}
    double norm() const { return std::hypot(x, y); }
    Point operator-(const Point& other) const { return {x - other.x, y - other.y}; }
    Point operator+(const Point& other) const { return {x + other.x, y + other.y}; }
    Point operator*(double scalar) const { return {x * scalar, y * scalar}; }
};

Point operator*(double scalar, const Point& p) { return p * scalar; }

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
        if (t_end_ <= t_start_) throw std::invalid_argument("t_end must be > t_start");
    }

    Point eval(double t) const { return r_func_(t); }
    Point derivative(double t) const { return dr_dt_func_(t); }
    double length() const {
        if (cached_length_ < 0) {
            auto f = [&](double t) { return dr_dt_func_(t).norm(); };
            cached_length_ = GaussIntegrator::integrate_real(f, t_start_, t_end_, 8);
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
};

// 微腔类
class Microcavity {
private:
    std::vector<ParametricSegment> segments_;
    std::vector<Point> boundary_points_;
    std::vector<Point> normals_;
    std::vector<double> curvatures_;
    std::vector<std::pair<const ParametricSegment*, double>> boundary_info_; // (段指针, 参数t)
    double max_segment_length_ = 0;

public:
    void add_segment(const ParametricSegment& seg) {
        if (!segments_.empty()) {
            Point prev_end = segments_.back().eval(segments_.back().parameter_range().second);
            Point new_start = seg.eval(seg.parameter_range().first);
            if ((prev_end - new_start).norm() > EPSILON) {
                throw std::runtime_error("Segments are not connected!");
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
        Point first = segments_.front().eval(segments_.front().parameter_range().first);
        Point last = segments_.back().eval(segments_.back().parameter_range().second);
        if ((first - last).norm() > EPSILON) {
            segments_.emplace_back(
                [&](double t) { return last + t * (first - last); },
                [](double) { return first - last; },
                0, 1
            );
        }
    }

    void generate_boundary(int points_per_segment) {
        boundary_points_.clear();
        boundary_info_.clear();

        for (const auto& seg : segments_) {
            auto points = seg.discretize(points_per_segment);
            auto [t_start, t_end] = seg.parameter_range();
            double dt = (t_end - t_start) / (points_per_segment - 1);
            for (int i = 0; i < points_per_segment; ++i) {
                double t = t_start + i * dt;
                boundary_points_.push_back(points[i]);
                boundary_info_.emplace_back(&seg, t);
            }
        }

        calculate_normals();
        calculate_curvatures();
    }

    const std::vector<Point>& boundary_points() const { return boundary_points_; }
    const std::vector<Point>& normals() const { return normals_; }
    const std::vector<double>& curvatures() const { return curvatures_; }
    double max_segment_length() const { return max_segment_length_; }

    std::pair<const ParametricSegment*, double> get_boundary_info(int index) const {
        if (index < 0 || index >= boundary_info_.size()) {
            throw std::out_of_range("Invalid boundary index");
        }
        return boundary_info_[index];
    }

private:
    void calculate_normals() {
        normals_.resize(boundary_points_.size());
        for (size_t i = 0; i < boundary_points_.size(); ++i) {
            size_t prev = (i == 0) ? boundary_points_.size() - 1 : i - 1;
            size_t next = (i == boundary_points_.size() - 1) ? 0 : i + 1;
            Point tangent = boundary_points_[next] - boundary_points_[prev];
            Point normal = {-tangent.y, tangent.x};
            double len = normal.norm();
            normals_[i] = len > 0 ? normal / len : Point(0, 1);
        }
    }

    void calculate_curvatures() {
        curvatures_.resize(boundary_points_.size());
        for (size_t i = 0; i < boundary_points_.size(); ++i) {
            size_t prev = (i == 0) ? boundary_points_.size() - 1 : i - 1;
            size_t next = (i == boundary_points_.size() - 1) ? 0 : i + 1;
            Point p0 = boundary_points_[prev], p1 = boundary_points_[i], p2 = boundary_points_[next];
            double dx1 = p1.x - p0.x, dy1 = p1.y - p0.y;
            double dx2 = p2.x - p1.x, dy2 = p2.y - p1.y;
            double cross = dx1 * dy2 - dx2 * dy1;
            double denom = std::pow((dx1*dx1 + dy1*dy1 + dx2*dx2 + dy2*dy2), 1.5);
            curvatures_[i] = denom > 0 ? 2 * cross / denom : 0;
        }
    }
};

// 高斯积分工具类
class GaussIntegrator {
public:
    static int select_order(double distance, double max_length) {
        if (distance < EPSILON) return 8;
        double ratio = distance / max_length;
        return ratio < 0.1 ? 8 : (ratio < 0.3 ? 6 : 4);
    }

    template<typename F>
    static Complex integrate(F f, const ParametricSegment& seg, int order) {
        auto [a, b] = seg.parameter_range();
        auto nodes = get_gauss_nodes(order);
        Complex sum(0, 0);
        double scale = (b - a) / 2, shift = (a + b) / 2;
        for (const auto& [x, w] : nodes) {
            double t = scale * x + shift;
            Point r = seg.eval(t);
            Point dr_dt = seg.derivative(t);
            sum += w * f(r) * dr_dt.norm();
        }
        return sum * scale;
    }

    static double integrate_real(std::function<double(double)> f, double a, double b, int order) {
        auto nodes = get_gauss_nodes(order);
        double sum = 0.0;
        double scale = (b - a) / 2, shift = (a + b) / 2;
        for (const auto& [x, w] : nodes) {
            double t = scale * x + shift;
            sum += w * f(t);
        }
        return sum * scale;
    }

private:
    static const std::vector<std::pair<double, double>>& get_gauss_nodes(int order) {
        static const std::vector<std::pair<double, double>> nodes_4 = {
            {-std::sqrt(1.0/3.0), 1.0},
            {std::sqrt(1.0/3.0), 1.0}
        };
        static const std::vector<std::pair<double, double>> nodes_6 = {
            {-std::sqrt(3.0/5.0), 5.0/9.0},
            {0, 8.0/9.0},
            {std::sqrt(3.0/5.0), 5.0/9.0}
        };
        static const std::vector<std::pair<double, double>> nodes_8 = {
            {-0.8611363115940526, 0.3478548451374538},
            {-0.3399810435848563, 0.6521451548625461},
            {0.3399810435848563, 0.6521451548625461},
            {0.8611363115940526, 0.3478548451374538}
        };
        switch(order) {
            case 4: return nodes_4;
            case 6: return nodes_6;
            case 8: return nodes_8;
            default: throw std::invalid_argument("Unsupported order");
        }
    }
};

// 边界元矩阵组装类
class BEMBuilder {
public:
    static MatrixXcd build_matrix(
        const Microcavity& cavity,
        double k,
        double n_in,
        double n_out
    ) {
        const auto& points = cavity.boundary_points();
        int N = points.size();
        MatrixXcd B(N, N), C(N, N);
        double max_length = cavity.max_segment_length();

        for (int i = 0; i < N; ++i) {
            auto [seg_i, t_i] = cavity.get_boundary_info(i);
            Point r_i = points[i];
            Point normal_i = cavity	normals()[i];

            for (int l = 0; l < N; ++l) {
                auto [seg_l, t_l] = cavity.get_boundary_info(l);
                double distance = (r_i - points[l]).norm();
                int order = GaussIntegrator::select_order(distance, max_length);

                if (i == l) {
                    double ds = seg_l->length() / (cavity.get_segment(0).discretize(100).size() - 1);
                    double curvature = cavity.curvatures()[i];
                    B(i, l) = -1.0 + (curvature * ds) / (2 * PI);
                    C(i, l) = (ds / PI) * (
                        1.0 - std::log(n_in * std::abs(k) * ds / 4.0) +
                        Complex(0, PI/2) - 0.5772156649
                    );
                } else {
                    auto f_B = [&](Point r) { return GreenDerivative(r, r_i, normal_i, k, n_out); };
                    B(i, l) = GaussIntegrator::integrate(f_B, *seg_l, order);

                    auto f_C = [&](Point r) { return GreenFunction(r, r_i, k, n_out); };
                    C(i, l) = GaussIntegrator::integrate(f_C, *seg_l, order);
                }
            }
        }

        MatrixXcd M(2*N, 2*N);
        M << B, C,
             MatrixXcd::Identity(N, N), MatrixXcd::Zero(N, N);
        return M;
    }

private:
    static Complex GreenFunction(const Point& r, const Point& r_prime, double k, double n) {
        double rho = (r - r_prime).norm();
        return rho < EPSILON ? Complex(0,0) : -Complex(0, 0.25) * sp_bessel::hankelH1(0, n*k*rho);
    }

    static Complex GreenDerivative(const Point& r, const Point& r_prime,
                                   const Point& normal, double k, double n) {
        Point delta = r - r_prime;
        double rho = delta.norm();
        if (rho < EPSILON) return Complex(0,0);
        Complex h1 = sp_bessel::hankelH1(1, n*k*rho);
        double cos_theta = (delta.x*normal.x + delta.y*normal.y) / rho;
        return -Complex(0, 0.25) * n*k * h1 * cos_theta;
    }
};

Complex find_resonant_wavenumber(const Microcavity& cavity, double initial_k, double n_in, double n_out) {
    MatrixXcd M = BEMBuilder::build_matrix(cavity, initial_k, n_in, n_out);
    Complex det = M.determinant();
    return initial_k - det / (1e-6 + det);
}

int main() {
    auto circle_r = [](double t) { return Point(std::cos(t), std::sin(t)); };
    auto circle_dr = [](double t) { return Point(-std::sin(t), std::cos(t)); };
    ParametricSegment circle(circle_r, circle_dr, 0, 2*PI);

    Microcavity cavity;
    cavity.add_segment(circle);
    cavity.close();
    cavity.generate_boundary(100);

    double k = 1.0, n_in = 1.5, n_out = 1.0;
    Complex resonant_k = find_resonant_wavenumber(cavity, k, n_in, n_out);
    std::cout << "Resonant wavenumber: " << resonant_k << std::endl;

    return 0;
}