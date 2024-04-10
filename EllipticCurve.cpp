#include "EllipticCurve.h"

void print_ecparams(ECParams params)
{
	cout << "椭圆曲线的有限域p：" << params.p << endl;
	cout << "椭圆曲线的参数a：" << params.a << endl;
	cout << "椭圆曲线的参数b：" << params.b << endl;
	cout << "椭圆曲线的阶n:" << params.n << endl;
	cout << "椭圆曲线的基点X坐标：" << params.Gx << endl;
	cout << "椭圆曲线的基点Y坐标:" << params.Gy << endl;
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

//求 P + Q
ECPoint ecpoint_add(ECPoint P, ECPoint Q, ECParams params)
{
	ECPoint R; //R=P+Q
	BigNumber x1, x2, x3, y1, y2, y3, a, b, p;
	x1 = P.x;
	y1 = P.y;
	x2 = Q.x;
	y2 = Q.y;
	a = params.a;
	b = params.b;
	p = params.p;
	if (!is_in_params(P, params) || !is_in_params(Q, params))
		cout << "P、Q必须在椭圆曲线上" << endl;
	BigNumber k;	//斜率k
	if (P.x != Q.x) {
		k = ((y2 - y1) * mod_inverse(x2 - x1, p)) % p;
	}
	else {
		k = ((3 * x1 * x1 + a) * mod_inverse(2 * y1 , p)) % p;
	}
	x3 = (k * k - x1 - x2) % p;
	y3 = (k * (x1 - x3) - y1) % p;
	x3 = (x3 % p + p) % p;
	y3 = (y3 % p + p) % p;
	R.x = x3;
	R.y = y3;
	if (is_in_params(R, params) == 1)
		cout << "结果正确!" << endl;
	else
		cout << "结果错误!" << endl;
	return R;
}

//求 kP
ECPoint ecpoint_mul_1(BigNumber k, ECPoint P, ECParams params)
{
	//朴素乘法1-k 非常慢
	ECPoint R = P;
	BigNumber p = params.p;
	k = (k % p + p) % p; //转化k为正值
	k = k - 1;
	while (k > 0) {
		R = ecpoint_add(R, P,params);
		k = k - 1;
	}
	return R;
}

//求最大公约数的同时求 ax + by = gcd(a,b) 的解
//exgcd结束时，x、y就是所求解
BigNumber exgcd(BigNumber a, BigNumber b, BigNumber& x, BigNumber& y)
{
	//递归边界
	if (b == 0) {
		x = 1, y = 0;
		return a;
	}
	//递归计算最大公约数gcd
	BigNumber gcd = exgcd(b, a % b, x, y);
	//递推公式，求解
	BigNumber temp = x;
	x = y;
	y = temp - a / b * y;
	return gcd;
}

BigNumber mod_inverse(BigNumber a, BigNumber b)
{
	BigNumber x, y;
	BigNumber gcd = exgcd(a, b, x, y);
	//cout << gcd << endl;
	//if (gcd != 1) {
	//	// 如果最大公约数不是 1，则没有逆元
	//	std::cerr << "逆元不存在，因为最大公约数不是1。" << std::endl;
	//	return -1;
	//}
	//x= (x % b + b) % b;	//最小整数解，也就是a模b的逆元
	//cout << a * x % b << endl;
	return x * gcd;	//若a、b不是互质的则返回x*gcd以确保计算正确
}

BigNumber random(BigNumber n)
{
	return BigNumber();
}
