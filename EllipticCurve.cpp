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
	if (P.x == 0 && P.y == 0) return Q;
	if (Q.x == 0 && Q.y == 0) return P;
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
//朴素乘法 1-k 非常慢
ECPoint ecpoint_mul_1(BigNumber k, ECPoint P, ECParams C)
{
	if (k == 0) return { BigNumber(0), BigNumber(0) };
	if (k == 1) return P;
	ECPoint R = P;
	while (k > 1) {
		R = ecpoint_add(R, P,C);
		k = k - 1;
	}
	return R;
}
//将k二进制表示
ECPoint ecpoint_mul_2(BigNumber k, ECPoint P, ECParams C)
{
	if (k == 0) return { BigNumber(0), BigNumber(0) };
	if (k == 1) return P;
	ECPoint R = {BigNumber(0),BigNumber(0)};
	ECPoint L = P;
	BigNumber k2 = k;
	while (k2 % 2 == 0) {
		k2 = k2 / 2;
		L = ecpoint_add(L, L, C);
	}
	while (k2 > 0) {
		if (k2 % 2 == 1) {
			R = ecpoint_add(R, L, C);
		}
		L = ecpoint_add(L, L, C);
		k2 = k2 / 2;
	}
	return R;
}
//使用NAF算法(比二进制表示好点)
//标量的非相邻表示称为NAF
//定义：一个正整数k的非相邻表示型是表达式, 其中, 并且没有两个连续的数字ki是非零的, NAF的长度是l.
//NAF的性质：
//1 k有唯一的NAF, 并记作NAF(k)。
//2 NAF(k)在k 的所有带符号表示中具有最少的非零位。
//3 NAF(k)的长度最多比k的二进制表示的长度大1.
//4 若NAF(k)的长度是l, 则2l / 3 < k < 2(l + 1) / 3.
//5 所有长度为l的NAF中非零数字的平均密度约为1 / 3.

ECPoint ecpoint_mul_NAF(BigNumber k, ECPoint P, ECParams C)
{
	if (k == 0) return { BigNumber(0), BigNumber(0) };
	if (k == 1) return P;
	ECPoint R = { BigNumber(0),BigNumber(0) };
	ECPoint L = P;
	BigNumber k2 = k;

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
