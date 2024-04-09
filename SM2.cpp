#include "SM2.h"

// 初始化大整数
void init_big_integer(BigInteger* num) {
    *num = NULL;
}
//初始化椭圆曲线上的点
void init_elliptic_curve_ecpoint(ECPoint* point) {
    init_big_integer(&(point->x));
    init_big_integer(&(point->y));
    ((point)->is_infinity) = 0;
}
// 初始化椭圆曲线参数
void init_elliptic_curve_params(EllipticCurveParams* params) {
    init_big_integer(&(params->p));
    init_big_integer(&(params->a));
    init_big_integer(&(params->b));
    init_big_integer(&(params->n));
    init_elliptic_curve_ecpoint(&(params)->G);
}
// 赋值大整数
void assign_big_integer(BigInteger* dest, BigInteger src) {
    *dest = src;
}
// 赋值椭圆曲线上的点
void assign_elliptic_curve_ecpoint(ECPoint* dest, ECPoint src) {
    assign_big_integer(&(dest->x), src.x);
    assign_big_integer(&(dest->y), src.y);
    dest->is_infinity = src.is_infinity;
}
// 赋值椭圆曲线参数
void assign_elliptic_curve_params(EllipticCurveParams* dest, EllipticCurveParams src) {
    assign_big_integer(&(dest->p), src.p);
    assign_big_integer(&(dest->a), src.a);
    assign_big_integer(&(dest->b), src.b);
    assign_big_integer(&(dest->n), src.n);
    assign_elliptic_curve_ecpoint(&(dest->G), src.G);
}
// 打印大整数
void print_big_integer(BigInteger num) {
    cout << num << endl;
}
// 打印椭圆曲线上的点
void print_elliptic_curve_ecpoint(ECPoint point) {
    print_big_integer(point.x);
    print_big_integer(point.y);
}
// 打印椭圆曲线参数
void print_elliptic_curve_params(EllipticCurveParams params) {
    cout<<"p: "; print_big_integer(params.p); 
    cout<<"a: "; print_big_integer(params.a); 
    cout<<"b: "; print_big_integer(params.b); 
    cout<<"n: "; print_big_integer(params.n); 
    cout<<"G: "; print_elliptic_curve_ecpoint(params.G); 
}
// 求大整数取模
BigInteger mod(BigInteger num, BigInteger mod_num) {
    BigInteger result;
    // 在这里实现大整数取模的逻辑
    result = NULL;
    return result;
}

BigInteger mod_inverse(BigInteger a, BigInteger m) {
    // 计算模m下的a的逆元
    BigInteger result;
    result = NULL;
    return result;
}
// 椭圆曲线点相加 P+Q
ECPoint ec_point_add(ECPoint P, ECPoint Q, ECParams params) {
    ECPoint result;
    // 在这里实现椭圆曲线点相加的逻辑
    return result;
}

// 椭圆曲线点倍乘 2P
ECPoint ec_point_double(ECPoint P, ECParams params) {
    ECPoint result;
    // 在这里实现椭圆曲线点倍乘的逻辑
    return result;
}

