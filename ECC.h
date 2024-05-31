#pragma once
//这里对实现椭圆曲线仿射坐标，标准射影坐标，Jacobian加重射影坐标进行实现
//利用拓展欧几里得实现求逆运算
//实现三种坐标系下二进制表示，NAF表示，w-NAF表示的点乘运算
//使用蒙哥马利算法对点乘运算进行优化
/*
曲线参数
p=FFFFFFFE FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF 00000000 FFFFFFFF FFFFFFFF
a=FFFFFFFE FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF 00000000 FFFFFFFF FFFFFFFC
b=28E9FA9E 9D9F5E34 4D5A9E4B CF6509A7 F39789F5 15AB8F92 DDBCBD41 4D940E93
n=FFFFFFFE FFFFFFFF FFFFFFFF FFFFFFFF 7203DF6B 21C6052B 53BBF409 39D54123
Gx=32C4AE2C 1F198119 5F990446 6A39C994 8FE30BBF F2660BE1 715A4589 334C74C7
Gy=BC3736A2 F4F6779C 59BDCEE3 6B692153 D0A9877C C62A4740 02DF32E5 2139F0A0
*/
#include"BIGNUM.h"
#include<windows.h>
using namespace std;

//仿射坐标 Affine
class Point {
public:
	BIGNUM x;
	BIGNUM y;
};

//标准射影坐标
class SPoint {
public:
	BIGNUM x;
	BIGNUM y;
	BIGNUM z;
};

//Jacobian加重射影坐标
class JPoint {
public:
	BIGNUM x;
	BIGNUM y;
	BIGNUM z;
};

//定义椭圆曲线参数结构体
class  Params {
public:
	BIGNUM p;	//有限域的模数
	BIGNUM a;	//椭圆曲线参数a
	BIGNUM b;	//椭圆曲线参数b
	BIGNUM n;	//群的阶
	BIGNUM Gx;	//基点坐标
	BIGNUM Gy;
};

void print_Params(Params params);		//打印椭圆曲线的参数
void print_Point(Point point);			//打印仿射坐标
void print_SPoint(SPoint point);		//打印标准射影坐标
void print_JPoint(JPoint point);		//打印Jacobian加重射影坐标

bool is_in_Params(Point point, Params params);		//仿射
bool is_in_Params_S(SPoint point, Params params);	//标准
bool is_in_Params_J(JPoint point, Params params);	//Jacobian

SPoint Affine_To_StandardProjection(Point point); //仿射坐标转换为标准射影坐标
Point StandardProjection_To_Affine(SPoint point, Params params); //标准射影坐标转换为仿射坐标 AffineToStandardProjection
JPoint Affine_To_Jacobian(Point point);//仿射坐标转换为雅可比坐标
Point Jacobian_To_Affine(JPoint point, Params params);//雅可比坐标转换为仿射坐标

//拓展gcd求逆元
BIGNUM EXGCD(BIGNUM a, BIGNUM b, BIGNUM& x, BIGNUM& y);
BIGNUM INVERSE(BIGNUM a, BIGNUM p);

//蒙哥马利模乘
BIGNUM M_Multiply(BIGNUM a, BIGNUM b, BIGNUM N);
//蒙哥马利约简
BIGNUM M_Reduction(BIGNUM X, BIGNUM N);

BIGNUM Random(int d, BIGNUM N); //生成一个随机数 1 < m < n - 1

Point Point_Add(Point P, Point Q, Params C);	//两点加 仿射坐标
Point Point_Add_M(Point PR, Point QR, Params C);	//蒙哥马利
SPoint SPoint_Add(SPoint P, SPoint Q, Params C); //两点加 标准射影坐标
SPoint SPoint_Add_M(SPoint PR, SPoint QR, Params C);//蒙哥马利
JPoint JPoint_Add(JPoint P, JPoint Q, Params C);//两点加 Jacobian加重射影坐标 
JPoint JPoint_Add_M(JPoint PR, JPoint QR, Params C);//蒙哥马利

//点乘方法实现
//二进制表示
Point Point_Mul_Bin(BIGNUM k, Point P, Params C);
Point Point_Mul_Bin_M(BIGNUM k, Point P, Params C);//蒙哥马利
SPoint SPoint_Mul_Bin(BIGNUM k, SPoint P, Params C);
SPoint SPoint_Mul_Bin_M(BIGNUM k, SPoint P, Params C);//蒙哥马利
JPoint JPoint_Mul_Bin(BIGNUM k, JPoint P, Params C);
JPoint JPoint_Mul_Bin_M(BIGNUM k, JPoint P, Params C);//蒙哥马利

//NAF表示
Point Point_Mul_NAF(BIGNUM k, Point P, Params C);
Point Point_Mul_NAF_M(BIGNUM k, Point P, Params C);//蒙哥马利
SPoint SPoint_Mul_NAF(BIGNUM k, SPoint P, Params C);
SPoint SPoint_Mul_NAF_M(BIGNUM k, SPoint P, Params C);//蒙哥马利
JPoint JPoint_Mul_NAF(BIGNUM k, JPoint P, Params C);
JPoint JPoint_Mul_NAF_M(BIGNUM k, JPoint P, Params C);//蒙哥马利

//w-NAF表示
Point Point_Mul_wNAF(BIGNUM k, Point P, Params C, int w);
Point Point_Mul_wNAF_M(BIGNUM k, Point P, Params C, int w);//蒙哥马利
SPoint SPoint_Mul_wNAF(BIGNUM k, SPoint P, Params C, int w);
SPoint SPoint_Mul_wNAF_M(BIGNUM k, SPoint P, Params C, int w);//蒙哥马利
JPoint JPoint_Mul_wNAF(BIGNUM k, JPoint P, Params C, int w);
JPoint JPoint_Mul_wNAF_M(BIGNUM k, JPoint P, Params C, int w);//蒙哥马利