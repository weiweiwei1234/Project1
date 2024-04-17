//椭圆曲线上的运算
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
#include "BIGNUM.h"
#include <windows.h>
using namespace std;
//暂不考虑无穷远点，实际计算基本使用不到
//定义椭圆曲线上的点	仿射坐标 Affine
class EccPoint
{
public:
	BIGNUM x;	//x
	BIGNUM y;	//y
};

class EccPointStandardProjection	//标准射影坐标
{
public:
	BIGNUM x;
	BIGNUM y;
	BIGNUM z;
};

class EccPointJacobian	//Jacobian加重射影坐标
{
public:
	BIGNUM x;
	BIGNUM y;
	BIGNUM z;
};

//定义椭圆曲线参数结构体
typedef struct EccParams {
	BIGNUM p;	//有限域的模数
	BIGNUM a;	//椭圆曲线参数a
	BIGNUM b;	//椭圆曲线参数b
	BIGNUM n;	//群的阶
	BIGNUM Gx;	//基点坐标
	BIGNUM Gy;
}EllipticCurveParams;

void printEccParams(EccParams C);	//打印椭圆曲线的参数
void printEccPoint(EccPoint point);		//打印点坐标
void printEccPointStandarProjection(EccPointStandardProjection point); //打印标准射影坐标
void printEccPointJacobian(EccPointJacobian point);//打印Jacobian加重射影坐标
bool isinEccParams(EccPoint point, EccParams C);	//判断点是否在椭圆曲线上
bool isinEccParamsStandardProjection(EccPointStandardProjection point, EccParams C);//判断点是否在标准射影坐标上的椭圆曲线上
bool isinEccParamsJacobian(EccPointJacobian point, EccParams C);//判断是否在Jacobian的椭圆曲线上
EccPointStandardProjection AffineToStandardProjection(EccPoint P); //仿射坐标转换为标准射影坐标
EccPoint StandardProjectionToAffine(EccPointStandardProjection P,EccParams C); //标准射影坐标转换为仿射坐标 AffineToStandardProjection
EccPointJacobian AffineTOJacobian(EccPoint P);//仿射坐标转换为雅可比坐标
EccPoint JacobianToAffine(EccPointJacobian P, EccParams C);//雅可比坐标转换为仿射坐标
EccPoint EccPointadd(EccPoint P, EccPoint Q, EccParams C);	//两点加 仿射坐标
EccPointStandardProjection EccPointaddStandardProjection(EccPointStandardProjection P, EccPointStandardProjection Q, EccParams C); //两点加 标准射影坐标
EccPointJacobian EccPointaddJacobian(EccPointJacobian P, EccPointJacobian Q, EccParams C);//两点加 雅可比坐标
//点乘算法实现
EccPoint EccPointmul1(BIGNUM k, EccPoint P, EccParams C);	//简单循环
//加减链方法
EccPoint EccPointmulBIN(BIGNUM k, EccPoint P, EccParams C);	//将k二进制表示
EccPoint EccPointmulNAF(BIGNUM k, EccPoint P, EccParams C);	//将k用NAF表示
EccPoint EccPointmulW_NAF(BIGNUM k, EccPoint P,int w,EccParams C);
/*
w-NAF算法 预见计算
输入：k, P, 窗口宽度w, 椭圆曲线参数C
输出：计算结果kP
*/
EccPointStandardProjection EccPointmulNAFStandardProjection(BIGNUM k, EccPointStandardProjection P, EccParams C);	//标准射影坐标下的NAF点乘
EccPointJacobian EccPointmulNAKJacobian(BIGNUM k, EccPointJacobian P, EccParams C);//雅可比坐标下的NAF点乘

EccPoint EccPointmul4(BIGNUM k, EccPoint P, EccParams C);

//拓展gcd求逆元 a * a^-1 = 1 (mod b)
BIGNUM exgcd(BIGNUM a, BIGNUM b, BIGNUM& x, BIGNUM& y);
BIGNUM modinverse(BIGNUM a,	BIGNUM b);  
BIGNUM random(BIGNUM n); //生成一个随机数 ∈[1,n-1]

