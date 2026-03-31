// 文件名：bem_point.h
#ifndef BEM_POINT_H
#define BEM_POINT_H

#include <cmath>
#include <iostream>
#include <iomanip>
#include "bem_constants.h"

namespace bem {

// const double EPS = 1e-8;  // 浮点精度控制

class Point {
public:
    // 构造函数
    Point() = default;
    Point(double x, double y) : x_(x), y_(y) {}
    
    // 基础运算
    Point& operator+=(const Point& other) {
        x_ += other.x_;
        y_ += other.y_;
        return *this;
    }
    
    Point& operator-=(const Point& other) {
        x_ -= other.x_;
        y_ -= other.y_;
        return *this;
    }
    
    Point& operator*=(double scalar) {
        x_ *= scalar;
        y_ *= scalar;
        return *this;
    }
    
    Point& operator/=(double scalar) {
        return *this *= 1.0 / scalar;
    }

    // 几何运算
    double norm() const { return std::hypot(x_, y_); }         // 模长
    double distance(const Point& other) const { return (*this - other).norm(); } // 距离
    double dot(const Point& other) const { return x_ * other.x_ + y_ * other.y_; } // 点积
    double cross(const Point& other) const { return x_ * other.y_ - y_ * other.x_; } // 叉积（二维z分量）
    
    // 极坐标转换
    static Point polar(double r, double theta) {  // 静态工厂方法
        return {r * std::cos(theta), r * std::sin(theta)};
    }

    Point normalized() const {
        double n = norm();
        return (n > EPS/100.0) ? Point(x_ / n, y_ / n) : Point(0.0, 0.0);  // 避免零除
    }

    // 边界元专用方法
    // Point tangent() const { return {-y_, x_}; }       // 切线方向（逆时针）
    Point normal() const { return {y_, -x_}; }        // 外法线方向（逆时针边界）
    double arc_length(const Point& prev, const Point& next) const { 
        // 三点弧长近似（用于曲线离散）
        return 0.5 * (prev.distance(*this) + this->distance(next)); 
    }

    // 运算符重载（友元保持对称性）
    friend Point operator+(Point a, const Point& b) { return a += b; }
    friend Point operator-(Point a, const Point& b) { return a -= b; }
    friend Point operator*(Point a, double s) { return a *= s; }
    friend Point operator*(double s, Point a) { return a *= s; }
    friend Point operator/(Point a, double s) { return a /= s; }

    // 比较运算符（带精度控制）
    friend bool operator==(const Point& a, const Point& b) {
        return std::abs(a.x_ - b.x_) < EPS && std::abs(a.y_ - b.y_) < EPS;
    }
    
    friend bool operator!=(const Point& a, const Point& b) {
        return !(a == b);
    }

    // 流运算符（用于调试输出）
    friend std::ostream& operator<<(std::ostream& os, const Point& p) {
        return os << std::fixed << std::setprecision(6) 
                  << "(" << p.x_ << ", " << p.y_ << ")";
    }

    // 数据成员（开放访问，性能优先）
    double x_ = 0.0;
    double y_ = 0.0;
};

} // namespace bem

#endif // BEM_POINT_H