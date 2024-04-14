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
using namespace std;
//定义椭圆曲线上的点	仿射坐标 Affine 
typedef struct ECPoint
{
	BigNumber x;	//x
	BigNumber y;	//y
};

typedef struct ECPointStandardProjection	//标准射影坐标
{
	BigNumber x;
	BigNumber y;
	BigNumber z;
};

typedef struct ECPointJacobian	//Jacobian加重射影坐标
{
	BigNumber x;
	BigNumber y;
	BigNumber z;
};

//定义椭圆曲线参数结构体
typedef struct ECParams {
	BigNumber p;	//有限域的模数
	BigNumber a;	//椭圆曲线参数a
	BigNumber b;	//椭圆曲线参数b
	BigNumber n;	//群的阶
	BigNumber Gx;	//基点坐标
	BigNumber Gy;
}EllipticCurveParams;

void printecparams(ECParams params);	//打印椭圆曲线的参数
void printecpoint(ECPoint point);		//打印点坐标
void printecpointStandarProjection(ECPointStandardProjection point); //打印标准射影坐标

bool isinparams(ECPoint point, ECParams params);	//判断点是否在椭圆曲线上
ECPointStandardProjection StandardProjectionToAffine(ECPoint P); //仿射坐标转换为标准射影坐标
ECPoint AffineToStandardProjection(ECPointStandardProjection P,ECParams C); //仿射坐标转换为标准射影坐标

ECPoint ecpointadd(ECPoint P, ECPoint Q, ECParams params);	//两点加 仿射坐标
ECPointStandardProjection ecpointaddStandardProjection(ECPointStandardProjection P, ECPointStandardProjection Q, ECParams C); //两点加 标准射影坐标

//点乘算法实现
ECPoint ecpointmul1(BigNumber k, ECPoint P, ECParams C);		
ECPoint ecpointmulBIN(BigNumber k, ECPoint P, ECParams C);
ECPoint ecpointmulNAF(BigNumber k, ECPoint P, ECParams C);
ECPointStandardProjection ecpointmulNAFStandardProjection(BigNumber k, ECPointStandardProjection P, ECParams C);
ECPoint ecpointmul4(BigNumber k, ECPoint P, ECParams C);

//拓展gcd求逆元 a * a^-1 = 1 (mod b)
BigNumber exgcd(BigNumber a, BigNumber b, BigNumber& x, BigNumber& y);
BigNumber modinverse(BigNumber a,	BigNumber b);  
BigNumber random(BigNumber n); //生成一个随机数 ∈[1,n-1]

