// 文件名：bem_curve_segment.h
#ifndef BEM_CURVE_SEGMENT_H
#define BEM_CURVE_SEGMENT_H

#include "bem_point.h"
#include <vector>
#include <cmath>

namespace bem {

class CurveSegment {
public:
    // 基础类型定义
    using Param = double;  // 参数t ∈ [0, 1]
    using PointList = std::vector<Point>;

    // 构造/析构
    CurveSegment() 
        :is_closed_(false) {}
    virtual ~CurveSegment() = default;

    // 纯虚函数：参数评估
    virtual Point operator()(Param t) const = 0;       // 参数t对应的点
    virtual Point tangent(Param t) const = 0;         // 切线向量（dt方向）
    virtual double curvature(Param t) const = 0;      // 曲率（用于拐角处理）
    virtual double arc_length() const = 0;            // 曲线段总弧长
    virtual void move(const Point& p) = 0;            // 曲线段总弧长
    // virtual double curvature_at_midpoint(Param t) const = 0; // 中点曲率

    // 通用几何运算
    Point normal(Param t) const { 
        return tangent(t).normalized().normal();  // 外法线（逆时针边界约定）
    }

    // 积分支持（参数t → 物理坐标的雅可比）
    double jacobian(Param t) const { 
        return tangent(t).norm();  // dt到ds的转换系数
    }

    // 辅助属性
    virtual Point start() const=0;
    virtual Point end() const=0;
    bool is_closed() const { return is_closed_; }
    void set_closed(bool closed) { is_closed_ = closed; }

protected:
    // Point start_;    // 曲线段起点
    // Point end_;      // 曲线段终点
    bool is_closed_; // 是否闭合（用于环形边界）
};

// 线性曲线段（论文中的分段线性离散）
class LineSegment : public CurveSegment {
private:
    Point start_;    // 曲线段起点
    Point end_;      // 曲线段终点
public:
    LineSegment(const Point& p1, const Point& p2) 
        : CurveSegment(),start_(p1),end_(p2){}

    Point operator()(Param t) const override {
        return start_ + (end_ - start_) * t;
    }
    void move(const Point& p) override {
        start_+=p;
        end_+=p;
    }

    Point tangent(Param t) const override {
        return end_ - start_;  // 切线恒定
    }

    double curvature(Param t) const override { return 0.0; }

    double arc_length() const override {
        return start_.distance(end_);
    }

    Point start() const override {
        return start_;
    }
    Point end() const override {
        return end_;
    }
};

// 圆弧段（用于圆形微腔边界）
class ArcSegment : public CurveSegment {
public:
    // 构造函数：起点、终点、圆心、逆时针方向
    ArcSegment(const Point& center, const double& theta_start,const double& theta_end,const double& radius, bool ccw=true)
        : CurveSegment(),center_(center),theta_start_(theta_start),theta_end_(theta_end),radius_(radius),ccw_(ccw){
            start_=center+Point::polar(radius, theta_start);
            end_=center+Point::polar(radius, theta_end);
    }

    Point operator()(Param t) const override {
        double theta = theta_start_ + (theta_end_ - theta_start_) * t;
        return {center_.x_ + std::cos(theta)*radius_, 
                center_.y_ + std::sin(theta)*radius_};
    }

    Point tangent(Param t) const override {
        double theta = theta_start_ + (theta_end_ - theta_start_) * t;
        double dx = -std::sin(theta)*radius_;  // 逆时针切线方向
        double dy =  std::cos(theta)*radius_;
        return Point(dx, dy);
    }

    double curvature(Param t) const override { 
        return 1.0 / radius_;  // 圆弧曲率恒定
    }

    double arc_length() const override {
        return std::abs(theta_end_ - theta_start_) * radius_;
    }

    Point start() const override {
        return start_;
    }
    Point end() const override {
        return end_;
    }

    void move(const Point& p) override {
        start_+=p;
        end_+=p;
        center_+=p;
    }

private:
    Point center_;      // 圆心
    double radius_;
    double theta_start_ = 0.0;
    double theta_end_ = 0.0;
    bool ccw_ = true;   // 逆时针方向
    Point start_;
    Point end_;
};
class FilletBezierSegment : public CurveSegment {
public:
    FilletBezierSegment(
        const Point& p0, const Point& p3,
        const Point& v0, const Point& v3,
        double L
    )
        : p0_(p0), p3_(p3)
    {
        // 三次 Bézier：P(t)=∑ B_i^3(t) p_i
        // 控制点 p1,p2 取在切线方向上，距离端点 L/3
        p1_ = p0 + v0.normalized() * (L / 3.0);
        p2_ = p3 - v3.normalized() * (L / 3.0);
    }

    Point operator()(Param t) const override {
        double u = 1 - t;
        return u*u*u*p0_
            + 3*u*u*t * p1_
            + 3*u*t*t * p2_
            + t*t*t * p3_;
    }

    Point tangent(Param t) const override {
        double u = 1 - t;
        // P′(t) = 3 [ u*u (p1−p0) + 2u t (p2−p1) + t*t (p3−p2) ]
        return ((p1_ - p0_)*(3*u*u)
            + (p2_ - p1_)*(6*u*t)
            + (p3_ - p2_)*(3*t*t));
    }

    double curvature(Param t) const override {
        // κ = |P′×P″|/|P′|^3
        Point d1 = tangent(t);
        // 二阶导 P″(t) = 6 [ u(p2−2p1+p0) + t(p3−2p2+p1) ]
        double u = 1 - t;
        Point d2 = (p2_ - p1_ - (p1_ - p0_)) * (6*u)
                + (p3_ - p2_ - (p2_ - p1_)) * (6*t);
        double num = std::abs(d1.x_*d2.y_ - d1.y_*d2.x_);
        double den = std::pow(d1.norm(), 3);
        return den>1e-16 ? num/den : 0.0;
    }


    double arc_length() const override {
        // 可用高阶数值积分
        const int M = 200;
        double L = 0;
        for(int i=0;i<M;i++){
        double t1 = double(i)/M, t2=double(i+1)/M;
        Point a = (*this)((t1+t2)/2);
        Point b = operator()(t2);
        L += (b - operator()(t1)).norm();
        }
        return L;
    }


    Point start() const override { return p0_; }
    Point end()   const override { return p3_; }

    void move(const Point& p) override {
        p0_+=p;
        p1_+=p;
        p2_+=p;
        p3_+=p;
    }

private:
    Point p0_, p1_, p2_, p3_;
};


class QuinticHermiteSegment : public CurveSegment {
public:
    /**
        * p0,p1: 端点坐标
        * m0,m1: 端点处的一阶导向量 P'(0), P'(1)
        * a0,a1: 端点处的二阶导向量 P''(0), P''(1)
        */
    QuinticHermiteSegment(
        const Point& p0, const Point& p1,
        const Point& m0, const Point& m1,
        const Point& a0 = {0,0}, const Point& a1 = {0,0}
    ) : p0_(p0), p1_(p1), m0_(m0), m1_(m1), a0_(a0), a1_(a1) {}

    // 位置插值
    Point operator()(Param t) const override {
        double t2 = t*t, t3 = t2*t, t4 = t3*t, t5 = t4*t;
        // Quintic Hermite basis
        double h00 =  1 - 10*t3 + 15*t4 - 6*t5;
        double h10 =      t -  6*t3 +  8*t4 - 3*t5;
        double h20 = 0.5*t2 -1.5*t3 +1.5*t4 -0.5*t5;
        double h01 =     10*t3 - 15*t4 + 6*t5;
        double h11 =    -4*t3 +  7*t4 - 3*t5;
        double h21 = 0.5*t3 -   t4 + 0.5*t5;
        return p0_*h00
                + m0_*h10
                + a0_*h20
                + p1_*h01
                + m1_*h11
                + a1_*h21;
    }

    // 切线
    Point tangent(Param t) const override {
        double t2 = t*t, t3 = t2*t, t4 = t3*t;
        // 导数基函数
        double dh00 = -30*t2 + 60*t3 - 30*t4;
        double dh10 =  1 - 18*t2 + 32*t3 - 15*t4;
        double dh20 =       t  - 4.5*t2 + 6*t3 -2.5*t4;
        double dh01 =  30*t2 - 60*t3 + 30*t4;
        double dh11 = -12*t2 + 28*t3 - 15*t4;
        double dh21 = 1.5*t2 - 4*t3 + 2.5*t4;
        return p0_*dh00
                + m0_*dh10
                + a0_*dh20
                + p1_*dh01
                + m1_*dh11
                + a1_*dh21;
    }

    // 曲率 = |P'×P''|/|P'|^3
    double curvature(Param t) const override {
        Point d1 = tangent(t);
        // 二阶导数
        double t2 = t*t, t3 = t2*t;
        double d2h00 = -60*t + 180*t2 - 120*t3;
        double d2h10 =     -36*t + 96*t2 - 60*t3;
        double d2h20 =      1 - 9*t + 18*t2 - 10*t3;
        double d2h01 =  60*t - 180*t2 + 120*t3;
        double d2h11 =   -24*t + 84*t2 - 60*t3;
        double d2h21 =    3*t - 12*t2 + 12.5*t3;
        Point d2 = p0_*d2h00
                    + m0_*d2h10
                    + a0_*d2h20
                    + p1_*d2h01
                    + m1_*d2h11
                    + a1_*d2h21;
        double num = std::abs(d1.x_*d2.y_ - d1.y_*d2.x_);
        double den = std::pow(d1.norm(), 3);
        return den>1e-16 ? num/den : 0.0;
    }

    double arc_length() const override {
        // 简单切分积分
        const int M = 200;
        double L = 0;
        for(int i=0;i<M;i++){
            double u = i/(double)M, v = (i+1)/(double)M;
            L += ((*this)( (u+v)/2 ) - (*this)(u)).norm();
        }
        return L;
    }

    Point start() const override { return p0_; }
    Point end()   const override { return p1_; }
    void move(const Point& p) override {
        p0_+=p;
        p1_+=p;
    }
private:
    Point p0_, p1_, m0_, m1_, a0_, a1_;
};
    

} // namespace bem

#endif // BEM_CURVE_SEGMENT_H