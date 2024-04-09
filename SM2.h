/*
推荐使用素数域256位椭圆曲线。
椭圆曲线方程：y2=x3+ax+b。
曲线参数
p=FFFFFFFE FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF 00000000 FFFFFFFF FFFFFFFF
a=FFFFFFFE FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF 00000000 FFFFFFFF FFFFFFFC
b=28E9FA9E 9D9F5E34 4D5A9E4B CF6509A7 F39789F5 15AB8F92 DDBCBD41 4D940E93
n=FFFFFFFE FFFFFFFF FFFFFFFF FFFFFFFF 7203DF6B 21C6052B 53BBF409 39D54123
Gx=32C4AE2C 1F198119 5F990446 6A39C994 8FE30BBF F2660BE1 715A4589 334C74C7
Gy=BC3736A2 F4F6779C 59BDCEE3 6B692153 D0A9877C C62A4740 02DF32E5 2139F0A0
*/
// 定义大整数类型
#pragma once
#include <string>
#include <iostream>
using namespace std;
typedef char* BigInteger;

//定义椭圆曲线上的点
typedef struct
{
    BigInteger x; //x坐标
    BigInteger y; //y坐标
    int is_infinity;  // 用于表示无穷远点
}ECPoint;

// 定义椭圆曲线参数结构体
typedef struct ECParams {
    BigInteger p;  // 有限域的模数
    BigInteger a;  // 椭圆曲线参数a
    BigInteger b;  // 椭圆曲线参数b
    BigInteger n;  // 群的阶
    ECPoint G;     // 群基点
} EllipticCurveParams;

void init_big_integer(BigInteger* num);
void init_elliptic_curve_ecpoint(ECPoint* point);
void init_elliptic_curve_params(EllipticCurveParams* params);

void assign_big_integer(BigInteger* dest, BigInteger src);
void assign_elliptic_curve_ecpoint(ECPoint* dest, ECPoint src);
void assign_elliptic_curve_params(EllipticCurveParams* dest, EllipticCurveParams src);

void print_big_integer(BigInteger num);
void print_elliptic_curve_ecpoint(ECPoint point);
void print_elliptic_curve_params(EllipticCurveParams params);

BigInteger mod(BigInteger num, BigInteger mod_num);
BigInteger mod_inverse(BigInteger a, BigInteger m);
ECPoint ec_point_add(ECPoint P, ECPoint Q, ECParams params);
ECPoint ec_point_double(ECPoint P, ECParams params);

