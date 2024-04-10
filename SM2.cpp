#include "SM2.h"

void print_ecparams(ECParams params)
{
	cout << "椭圆曲线的有限域：" << params.p << endl;
	cout << "椭圆曲线的参数a：" << params.a << endl;
	cout << "椭圆曲线的参数b：" << params.b << endl;
	cout << "椭圆曲线的阶n" << params.n << endl;
	cout << "椭圆曲线的基点X坐标：" << params.Gx << endl;
	cout << "椭圆曲线的基点Y坐标" << params.Gy << endl;
}

void print_ecpoint(ECPoint point)
{
	cout << "X坐标：" << point.x << endl;
	cout << "Y坐标：" << point.y << endl;
}

bool is_in_params(ECPoint point, ECParams params)
{
	if ((point.y * point.y) % params.p == (point.x * point.x * point.x + params.a * point.x + params.b) % params.p)
		return true;
	return false;
}

ECPoint ecpoint_mul(BigNumber k, ECPoint P, ECParams params)
{
	return ECPoint();
}

void exgcd(BigNumber a, BigNumber b, BigNumber& x, BigNumber& y)
{

}

BigNumber mod_inverse(BigNumber a, BigNumber m)
{
	return BigNumber();
}
