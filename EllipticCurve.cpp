#include "EllipticCurve.h"
#include "SM3.h"
void printecparams(ECParams params)
{
	cout << "椭圆曲线参数：" << endl;
	cout << "有限域p：\t" << params.p << endl;
	cout << "系数a：\t\t" << params.a << endl;
	cout << "系数b：\t\t" << params.b << endl;
	cout << "阶n：\t\t" << params.n << endl;
	cout << "基点X坐标：\t" << params.Gx << endl;
	cout << "基点Y坐标：\t" << params.Gy << endl;
}

void printecpoint(ECPoint point)
{
	cout << "仿射坐标表示：" << endl;
	cout << "X坐标：\t" << point.x << endl;
	cout << "Y坐标：\t" << point.y << endl;
}

void printecpointStandarProjection(ECPointStandardProjection point)
{
	cout << "标准射影坐标表示" << endl;
	cout << "X坐标：\t" << point.x << endl;
	cout << "Y坐标：\t" << point.y << endl;
	cout << "Z坐标：\t" << point.z << endl;
}

bool isinparams(ECPoint point, ECParams params)
{
	if ((point.y * point.y) % params.p == (point.x * point.x * point.x + params.a * point.x + params.b) % params.p)
		return true;
	return false;
}

ECPointStandardProjection StandardProjectionToAffine(ECPoint P)
{
	ECPointStandardProjection R = { P.x,P.y,1 };
	return R;
}

ECPoint AffineToStandardProjection(ECPointStandardProjection P,ECParams C)
{
	ECPoint R;
	R.x = P.x * modinverse(P.z, C.p);
	R.y = P.y * modinverse(P.z, C.p);
	R.x = (R.x % C.p + C.p) % C.p;
	R.y = (R.y % C.p + C.p) % C.p;
	return R;
}

//求 P + Q
ECPoint ecpointadd(ECPoint P, ECPoint Q, ECParams params)
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
	if (!isinparams(P, params) || !isinparams(Q, params))
		cout << "P、Q必须在椭圆曲线上" << endl;
	BigNumber k;	//斜率k
	if (P.x != Q.x) {
		k = (y2 - y1) * modinverse(x2 - x1, p);
	}
	else {
		k = (3 * x1 * x1 + a) * modinverse(2 * y1 , p);
	}
	x3 = k * k - x1 - x2;
	y3 = k * (x1 - x3) - y1;
	x3 = (x3 % p + p) % p;
	y3 = (y3 % p + p) % p;
	R = { x3,y3 };
	return R;
}

//标准射影坐标的两点加
ECPointStandardProjection ecpointaddStandardProjection(ECPointStandardProjection P, ECPointStandardProjection Q, ECParams C)
{
	if (P.z == 0) return Q;
	if (Q.z == 0) return P;
	ECPointStandardProjection R;
	BigNumber x1, x2, x3, y1, y2, y3, z1, z2, z3, a, b, p;
	x1 = P.x;
	y1 = P.y;
	z1 = P.z;
	x2 = Q.x;
	y2 = Q.y;
	z2 = Q.z;
	a = C.a;
	b = C.b;
	p = C.p;
	if (P.x != Q.x && P.y != Q.y) {
		BigNumber t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11;
		t1 = x1 * z2;
		t2 = x2 * z1;
		t3 = t1 - t2;
		t4 = y1 * z2;
		t5 = y2 * z1;
		t6 = t4 - t5;
		t7 = t1 + t2;
		t8 = z1 * z2;
		t9 = t3 * t3;
		t10 = t3 * t9;
		t11 = t8 * t6 * t6 - t7 * t9;
		//得到R的射影坐标
		x3 = t3 * t11;
		y3 = t6 * (t9 * t1 - t11) - t4 * t10;
		z3 = t10 * t8;
	}
	else {
		BigNumber t1, t2, t3, t4, t5, t6;
		t1 = 3 * x1 * x1 + a * z1 * z1;
		t2 = 2 * y1 * z1;
		t3 = y1 * y1;
		t4 = t3 * x1 * z1;
		t5 = t2 * t2;
		t6 = t1 * t1 - 8 * t4;
		//得到R的射影坐标
		x3 = t2 * t6;
		y3 = t1 * (4 * t4 - t6) - 2 * t5 * t3;
		z3 = t2 * t5;
	}
	x3 = (x3 % p + p) % p;
	y3 = (y3 % p + p) % p;
	z3 = (z3 % p + p) % p;
	R = { x3,y3,z3 };
	return R;
}

//求 kP
//朴素乘法 1-k 非常慢
ECPoint ecpointmul1(BigNumber k, ECPoint P, ECParams C)
{
	if (k == 0) return { BigNumber(0), BigNumber(0) };
	if (k == 1) return P;
	ECPoint R = P;
	while (k > 1) {
		R = ecpointadd(R, P,C);
		k = k - 1;
	}
	return R;
}
//将k二进制表示
ECPoint ecpointmulBIN(BigNumber k, ECPoint P, ECParams C)
{
	if (k == 0) return { BigNumber(0), BigNumber(0) };
	if (k == 1) return P;
	ECPoint R = {BigNumber(0),BigNumber(0)};
	ECPoint L = P;
	BigNumber k1 = k;
	string k1_str = k1.get_value();
	string k1_str_BIN = HexToBin(k1_str);
	cout << k1 << "的二进制表示为" << k1_str_BIN << endl;
	while (k1 > 0) {
		if (k1 % 2 == 1) {
			R = ecpointadd(R, L, C);
		}
		L = ecpointadd(L, L, C);
		k1 = k1 / 2;
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

ECPoint ecpointmulNAF(BigNumber k, ECPoint P, ECParams C)
{
	if (k == 0) return { BigNumber(0), BigNumber(0) };
	if (k == 1) return P;
	ECPoint R = { BigNumber(0),BigNumber(0) };
	ECPoint _P = {P.x,0-P.y}; //-P
	BigNumber k1 = k;
	//计算k的NAF值
	int i = 0;
	BigNumber NAF_k[1025];
	while (k1 >= 1) {
		if (k1 % 2 == 1) {
			NAF_k[i] = 2 - (k1 % 4);
			k1 = k1 - NAF_k[i];
		}
		else {
			NAF_k[i] = 0;
		}
		k1 = k1 / 2;
		i++;
	}
	cout << k << "的NAF表示为：";
	int j;
	for (j = i - 1; j >= 0; j--) {
		cout << NAF_k[j];
	}
	cout << endl;
	for (j = i - 1; j >= 0; j--) {
		R = ecpointadd(R, R, C);
		if (NAF_k[j] == 1)
			R = ecpointadd(R, P, C);
		if (NAF_k[j] == -1)
			R = ecpointadd(R, _P, C);
	}
	return R;
}

//射影坐标
ECPointStandardProjection ecpointmulNAFStandardProjection(BigNumber k, ECPointStandardProjection P, ECParams C)
{
	if (k == 0) return { BigNumber(0),BigNumber(0),BigNumber(0) };
	if (k == 1) return P;
	ECPointStandardProjection R = { BigNumber(0),BigNumber(0) };
	ECPointStandardProjection _P = { P.x,0 - P.y ,P.z}; //-P
	BigNumber k1 = k;
	//计算k的NAF值
	int i = 0;
	BigNumber NAF_k[1025];
	while (k1 >= 1) {
		if (k1 % 2 == 1) {
			NAF_k[i] = 2 - (k1 % 4);
			k1 = k1 - NAF_k[i];
		}
		else {
			NAF_k[i] = 0;
		}
		k1 = k1 / 2;
		i++;
	}
	cout << k << "的NAF表示为：";
	int j;
	for (j = i - 1; j >= 0; j--) {
		cout << NAF_k[j];
	}
	cout << endl;
	for (j = i - 1; j >= 0; j--) {
		R = ecpointaddStandardProjection(R, R, C);
		if (NAF_k[j] == 1)
			R = ecpointaddStandardProjection(R, P, C);
		if (NAF_k[j] == -1)
			R = ecpointaddStandardProjection(R, _P, C);
	}
	return R;
}

ECPoint ecpointmul4(BigNumber k, ECPoint P, ECParams C)
{
	return ECPoint();
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

BigNumber modinverse(BigNumber a, BigNumber b)
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
