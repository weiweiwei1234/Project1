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
#include "BigNumber.h"
#include <windows.h>
using namespace std;
//暂不考虑无穷远点，实际计算基本使用不到
//定义椭圆曲线上的点	仿射坐标 Affine 
typedef struct ECPoint
{
	BigNumber x;	//x
	BigNumber y;	//y
}Point;

typedef struct ECPointStandardProjection	//标准射影坐标
{
	BigNumber x;
	BigNumber y;
	BigNumber z;
}PointSP;

typedef struct ECPointJacobian	//Jacobian加重射影坐标
{
	BigNumber x;
	BigNumber y;
	BigNumber z;
}PointJ;

//定义椭圆曲线参数结构体
typedef struct ECParams {
	BigNumber p;	//有限域的模数
	BigNumber a;	//椭圆曲线参数a
	BigNumber b;	//椭圆曲线参数b
	BigNumber n;	//群的阶
	BigNumber Gx;	//基点坐标
	BigNumber Gy;
}EllipticCurveParams;

void printecparams(ECParams C);	//打印椭圆曲线的参数
void printecpoint(ECPoint point);		//打印点坐标
void printecpointStandarProjection(ECPointStandardProjection point); //打印标准射影坐标
void printecpointJacobian(ECPointJacobian point);//打印Jacobian加重射影坐标
bool isinparams(ECPoint point, ECParams C);	//判断点是否在椭圆曲线上
bool isinparamsStandardProjection(ECPointStandardProjection point, ECParams C);//判断点是否在标准射影坐标上的椭圆曲线上
bool isinparamsJacobian(ECPointJacobian point, ECParams C);//判断是否在Jacobian的椭圆曲线上
ECPointStandardProjection AffineToStandardProjection(ECPoint P); //仿射坐标转换为标准射影坐标
ECPoint StandardProjectionToAffine(ECPointStandardProjection P,ECParams C); //标准射影坐标转换为仿射坐标 AffineToStandardProjection
ECPointJacobian AffineTOJacobian(ECPoint P);//仿射坐标转换为雅可比坐标
ECPoint JacobianToAffine(ECPointJacobian P, ECParams C);//雅可比坐标转换为仿射坐标
ECPoint ecpointadd(ECPoint P, ECPoint Q, ECParams C);	//两点加 仿射坐标
ECPointStandardProjection ecpointaddStandardProjection(ECPointStandardProjection P, ECPointStandardProjection Q, ECParams C); //两点加 标准射影坐标
ECPointJacobian ecpointaddJacobian(ECPointJacobian P, ECPointJacobian Q, ECParams C);//两点加 雅可比坐标
//点乘算法实现
ECPoint ecpointmul1(BigNumber k, ECPoint P, ECParams C);	//简单循环
//加减链方法
ECPoint ecpointmulBIN(BigNumber k, ECPoint P, ECParams C);	//将k二进制表示
ECPoint ecpointmulNAF(BigNumber k, ECPoint P, ECParams C);	//将k用NAF表示
ECPoint ecpointmulW_NAF(BigNumber k, ECPoint P,int w,ECParams C);
/*
w-NAF算法 预见计算
输入：k, P, 窗口宽度w, 椭圆曲线参数C
输出：计算结果kP
*/
ECPointStandardProjection ecpointmulNAFStandardProjection(BigNumber k, ECPointStandardProjection P, ECParams C);	//标准射影坐标下的NAF点乘
ECPointJacobian ecpointmulNAKJacobian(BigNumber k, ECPointJacobian P, ECParams C);//雅可比坐标下的NAF点乘

ECPoint ecpointmul4(BigNumber k, ECPoint P, ECParams C);

//拓展gcd求逆元 a * a^-1 = 1 (mod b)
BigNumber exgcd(BigNumber a, BigNumber b, BigNumber& x, BigNumber& y);
BigNumber modinverse(BigNumber a,	BigNumber b);  
BigNumber random(BigNumber n); //生成一个随机数 ∈[1,n-1]

