#include "EllipticCurve.h"
#include "SM3.h"
void printEccParams(EccParams C)
{
	cout << "椭圆曲线参数：" << endl;
	cout << "有限域p：\t" << C.p << endl;
	cout << "系数a：\t\t" << C.a << endl;
	cout << "系数b：\t\t" << C.b << endl;
	cout << "阶n：\t\t" << C.n << endl;
	cout << "基点X坐标：\t" << C.Gx << endl;
	cout << "基点Y坐标：\t" << C.Gy << endl;
}

void printEccPoint(EccPoint point)
{
	cout << "仿射坐标表示：" << endl;
	cout << "X坐标：\t" << point.x << endl;
	cout << "Y坐标：\t" << point.y << endl;
}

void printEccPointStandarProjection(EccPointStandardProjection point)
{
	cout << "标准射影坐标表示:" << endl;
	cout << "X坐标：\t" << point.x << endl;
	cout << "Y坐标：\t" << point.y << endl;
	cout << "Z坐标：\t" << point.z << endl;
}

void printEccPointJacobian(EccPointJacobian point)
{
	cout << "雅可比坐标表示:" << endl;
	cout << "X坐标：\t" << point.x << endl;
	cout << "Y坐标：\t" << point.y << endl;
	cout << "Z坐标：\t" << point.z << endl;
}
//Y^2=X^3+aX+b
bool isinEccParams(EccPoint point, EccParams C)
{
	BIGNUM x, y, a, b, p;
	x = point.x;
	y = point.y;
	a = C.a;
	b = C.b;
	p = C.p;
	if ((y * y) % p == (x * x * x + a * x + b) % p)
		return true;
	return false;
}
//Y^2*Z=X^3+aXZ^2+bZ^3
bool isinEccParamsStandardProjection(EccPointStandardProjection point, EccParams C)
{
	BIGNUM x, y, z, a, b, p;
	x = point.x;
	y = point.y;
	z = point.z;
	a = C.a;
	b = C.b;
	p = C.p;
	if ((y * y * z) % p == (x * x * x + a * x * z * z + b * z * z * z) % p) {
		return true;
	}
	return false;
}
//Y^2=X^3+aXZ^4+bZ^6
bool isinEccParamsJacobian(EccPointJacobian point, EccParams C)
{
	BIGNUM x, y, z, a, b, p;
	x = point.x;
	y = point.y;
	z = point.z;
	a = C.a;
	b = C.b;
	p = C.p;
	if ((y * y) % p == (x * x * x + a * x * z * z * z * z + b * z * z * z * z * z * z) % p){
		return true;
	}
	return false;
}

EccPointStandardProjection AffineToStandardProjection(EccPoint P)
{
	EccPointStandardProjection R = { P.x,P.y,1 };
	return R;
}

EccPoint StandardProjectionToAffine(EccPointStandardProjection P,EccParams C)
{
	EccPoint R;
	R.x = P.x * Mod_inverse(P.z, C.p);
	R.y = P.y * Mod_inverse(P.z, C.p);
	R.x = (R.x % C.p + C.p) % C.p;
	R.y = (R.y % C.p + C.p) % C.p;
	return R;
}

EccPointJacobian AffineTOJacobian(EccPoint P)
{
	EccPointJacobian R = { P.x,P.y,1 };
	return R;
}

EccPoint JacobianToAffine(EccPointJacobian P, EccParams C)
{
	EccPoint R;
	R.x = P.x * Mod_inverse(P.z * P.z, C.p);
	R.y = P.y * Mod_inverse(P.z * P.z * P.z, C.p);
	R.x = (R.x % C.p + C.p) % C.p;
	R.y = (R.y % C.p + C.p) % C.p;
	return R;
}

//求 P + Q
EccPoint EccPointAdd(EccPoint P, EccPoint Q, EccParams C)
{
	if (P.x == 0 && P.y == 0) return Q;
	if (Q.x == 0 && Q.y == 0) return P;
	EccPoint R; //R=P+Q
	BIGNUM x1, x2, x3, y1, y2, y3, a, b, p;
	x1 = P.x;
	y1 = P.y;
	x2 = Q.x;
	y2 = Q.y;
	a = C.a;
	b = C.b;
	p = C.p;
	if (!isinEccParams(P, C) || !isinEccParams(Q, C))
		cout << "P、Q必须在椭圆曲线上" << endl;
	BIGNUM k;	//斜率k
	if (P.x != Q.x) {
		k = (y2 - y1) * Mod_inverse(x2 - x1, p);
	}
	else {
		k = (3 * x1 * x1 + a) * Mod_inverse(2 * y1 , p);
	}
	x3 = k * k - x1 - x2;
	y3 = k * (x1 - x3) - y1;
	x3 = (x3 % p + p) % p;
	y3 = (y3 % p + p) % p;
	R = { x3,y3 };
	return R;
}

//标准射影坐标的两点加 一般不使用
EccPointStandardProjection EccPointAddStandardProjection(EccPointStandardProjection P, EccPointStandardProjection Q, EccParams C)
{
	if (P.z == 0) return Q;
	if (Q.z == 0) return P;
	if (!isinEccParamsStandardProjection(P, C) || !isinEccParamsStandardProjection(Q, C))
		cout << "P、Q必须在椭圆曲线上" << endl;
	EccPointStandardProjection R;
	BIGNUM x1, x2, x3, y1, y2, y3, z1, z2, z3, a, b, p;
	x1 = P.x;
	y1 = P.y;
	z1 = P.z;
	x2 = Q.x;
	y2 = Q.y;
	z2 = Q.z;
	a = C.a;
	b = C.b;
	p = C.p;
	if (P.x != Q.x && P.y != Q.y && P.z != Q.z) {
		BIGNUM t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11;
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
		BIGNUM t1, t2, t3, t4, t5, t6;
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
//两点加 雅可比坐标
EccPointJacobian EccPointAddJacobian(EccPointJacobian P, EccPointJacobian Q, EccParams C)
{
	if (P.z == 0) return Q;
	if (Q.z == 0) return P;
	if (!isinEccParamsJacobian(P, C) || !isinEccParamsJacobian(Q, C))
		cout << "P、Q必须在椭圆曲线上" << endl;
	EccPointJacobian R;
	BIGNUM x1, x2, x3, y1, y2, y3, z1, z2, z3, a, b, p;
	x1 = P.x;
	y1 = P.y;
	z1 = P.z;
	x2 = Q.x;
	y2 = Q.y;
	z2 = Q.z;
	a = C.a;
	b = C.b;
	p = C.p;
	if (P.x != Q.x && P.y != Q.y && P.z != Q.z) {
		BIGNUM t1, t2, t3, t4, t5, t6, t7, t8, t9;
		t1 = x1 * z2 * z2;
		t2 = x2 * z1 * z1;
		t3 = t1 - t2;
		t4 = y1 * z2 * z2 * z2;
		t5 = y2 * z1 * z1 * z1;
		t6 = t4 - t5;
		t7 = t1 + t2;
		t8 = t4 + t5;
		//得到R的雅可比坐标
		x3 = t6 * t6 - t7 * t3 * t3;
		t9 = t7 * t3 * t3 - 2 * x3;
		y3 = (t9 * t6 - t8 * t3 * t3 * t3) * Mod_inverse(2, p);
		z3 = z1 * z2 * t3;
	}
	else {
		BIGNUM t1, t2, t3;
		t1 = 3 * x1 * x1 + a * z1 * z1 * z1 * z1;
		t2 = 4 * x1 * y1 * y1;
		t3 = 8 * y1 * y1 * y1 * y1;
		//得到R的雅可比坐标
		x3 = t1 * t1 - 2 * t2;
		y3 = t1 * (t2 - x3) - t3;
		z3 = 2 * y1 * z1;
	}
	x3 = (x3 % p + p) % p;
	y3 = (y3 % p + p) % p;
	z3 = (z3 % p + p) % p;
	R = { x3,y3,z3 };
	return R;
}

//求 kP
//循环 1-k 非常慢
EccPoint EccPointMul1(BIGNUM k, EccPoint P, EccParams C)
{
	if (k == 0) return { BIGNUM(0), BIGNUM(0) };
	if (k == 1) return P;
	EccPoint R = P;
	while (k > 1) {
		R = EccPointAdd(R, P,C);
		k = k - 1;
	}
	return R;
}
//将k二进制表示
EccPoint EccPointMulBIN(BIGNUM k, EccPoint P, EccParams C)
{
	if (k == 0) return { BIGNUM(0), BIGNUM(0) };
	if (k == 1) return P;
	EccPoint R = {BIGNUM(0),BIGNUM(0)};
	EccPoint L = P;
	//string k_str = k.get_value();
	//string k_str_BIN = HexToBin(k_str);
	//cout << k << "的二进制表示为" << k_str_BIN << endl;
	while (k > 0) {
		if (k % 2 == 1) {
			R = EccPointAdd(R, L, C);
		}
		L = EccPointAdd(L, L, C);
		k = k / 2;
	}
	return R;
}
//使用NAF算法
//标量的非相邻表示称为NAF
//定义：一个正整数k的非相邻表示型是表达式, 其中, 并且没有两个连续的数字ki是非零的, NAF的长度是l.
//NAF的性质：
//1 k有唯一的NAF, 并记作NAF(k)。
//2 NAF(k)在k 的所有带符号表示中具有最少的非零位。
//3 NAF(k)的长度最多比k的二进制表示的长度大1.
//4 若NAF(k)的长度是l, 则2l / 3 < k < 2(l + 1) / 3.
//5 所有长度为l的NAF中非零数字的平均密度约为1 / 3.

EccPoint EccPointMulNAF(BIGNUM k, EccPoint P, EccParams C)
{
	if (k == 0) return { BIGNUM(0), BIGNUM(0) };
	if (k == 1) return P;
	EccPoint R = { BIGNUM(0),BIGNUM(0) };
	BIGNUM k1 = k;
	//计算k的NAF值
	int i = 0;
	int NAF_k[1025] = { 0 };
	BIGNUM temp;
	while (k1 >= 1) {
		if (k1 % 2 == 1) {
			temp = 2 - (k1 % 4);
			k1 = k1 - temp;
			NAF_k[i] = temp.to_int();
		}
		else {
			NAF_k[i] = 0;
		}
		k1 = k1 / 2;
		i++;
	}
	for (int j = i - 1; j >= 0; j--) {
		R = EccPointAdd(R, R, C);
		if (NAF_k[j] == 1)
			R = EccPointAdd(R, P, C);
		if (NAF_k[j] == -1)
			R = EccPointAdd(R, { P.x,0 - P.y }, C);
	}
	return R;
}

EccPoint EccPointMulW_NAF(BIGNUM k, EccPoint P, int w, EccParams C)
{
	if (k == 0) return { BIGNUM(0),BIGNUM(0) };
	if (k == 1) return P;
	EccPoint R{ BIGNUM(0), BIGNUM(0) };
	BIGNUM k1 = k;
	//用w计算预计算表 计算iP i=1,3,5,...,2^(w-1)-1
	EccPoint Pi[64];
	EccPoint P_2 = EccPointAdd(P, P, C);
	for (int j = 1; j < (int)pow(2, w); j = j + 2) {
		if (j == 1) Pi[j] = P;
		else
			Pi[j] = EccPointAdd(Pi[j - 2], P_2, C);
	}
	int i = 0;
	int NAFW[1025] = { 0 };
	BIGNUM w2((int)pow(2, w));
	while (k1 >= 1) {
		if (k1 % 2 == 1) {
			NAFW[i] = BIGNUM(k1 % w2).to_int();
			while (NAFW[i] > pow(2, w - 1)) {
				NAFW[i] = NAFW[i] - pow(2, w);
			}
			k1 = k1 - BIGNUM(NAFW[i]);
		}
		else {
			NAFW[i] = 0;
		}
		k1 = k1 / 2;
		i++;
	}
	long t1, t2;//计算运行时间，t1:开始时间,t2:结束时间
	t1 = GetTickCount64();
	for (int j = i - 1; j >= 0; j--) {
		R = EccPointAdd(R, R, C);
		if (NAFW[j] > 0) {
			R = EccPointAdd(R, Pi[NAFW[j]], C);
		}
		if (NAFW[j] < 0) {
			R = EccPointAdd(R, {Pi[-NAFW[j]].x,0-Pi[-NAFW[j]].y}, C);
		}
	}
	t2 = GetTickCount64();
	cout << "执行时间：" << t2 - t1 << endl;  //程序运行的时间得到的时间单位为毫秒 /1000为秒
	return R;
}

//射影坐标 不使用
EccPointStandardProjection EccPointMulNAFStandardProjection(BIGNUM k, EccPointStandardProjection P, EccParams C)
{
	if (k == 0) return { BIGNUM(0),BIGNUM(1),BIGNUM(0) };
	if (k == 1) return P;
	EccPointStandardProjection R = { BIGNUM(0),BIGNUM(1),BIGNUM(0) };
	EccPointStandardProjection _P = { P.x,0 - P.y ,P.z}; //-P
	BIGNUM k1 = k;
	//计算k的NAF值
	int i = 0;
	int NAF_k[1025] = { 0 };
	BIGNUM temp;
	while (k1 >= 1) {
		if (k1 % 2 == 1) {
			temp = 2 - (k1 % 4);
			k1 = k1 - temp;
			NAF_k[i] = temp.to_int();
		}
		else {
			NAF_k[i] = 0;
		}
		k1 = k1 / 2;
		i++;
	}
	for (int j = i - 1; j >= 0; j--) {
		R = EccPointAddStandardProjection(R, R, C);
		if (NAF_k[j] == 1)
			R = EccPointAddStandardProjection(R, P, C);
		if (NAF_k[j] == -1)
			R = EccPointAddStandardProjection(R, _P, C);
	}
	return R;
}

EccPointJacobian EccPointMul_NAF_Jacobian(BIGNUM k, EccPointJacobian P, EccParams C)
{
	if (k == 0) return { BIGNUM(1),BIGNUM(1),BIGNUM(0) };
	if (k == 1) return P;
	EccPointJacobian R = { BIGNUM(1),BIGNUM(1),BIGNUM(0) };
	EccPointJacobian _P = { P.x,0 - P.y ,P.z }; //-P
	BIGNUM k1 = k;
	//计算k的NAF值
	int i = 0;
	int NAF_k[1025] = { 0 };
	BIGNUM temp;
	while (k1 >= 1) {
		if (k1 % 2 == 1) {
			temp = 2 - (k1 % 4);
			k1 = k1 - temp;
			NAF_k[i] = temp.to_int();
		}
		else {
			NAF_k[i] = 0;
		}
		k1 = k1 / 2;
		i++;
	}
	for (int j = i - 1; j >= 0; j--) {
		R = EccPointAddJacobian(R, R, C);
		if (NAF_k[j] == 1)
			R = EccPointAddJacobian(R, P, C);
		if (NAF_k[j] == -1)
			R = EccPointAddJacobian(R, _P, C);
	}
	return R;
}

EccPointJacobian EccPointMul_W_NAF_Jacobian(BIGNUM k, EccPointJacobian P, int w, EccParams C)
{
	if (k == 0) return { BIGNUM(1),BIGNUM(1),BIGNUM(0)};
	if (k == 1) return P;
	EccPointJacobian R{ BIGNUM(1), BIGNUM(1),BIGNUM(0) };
	BIGNUM k1 = k;
	//用w计算预计算表 计算iP i=1,3,5,...,2^(w-1)-1
	EccPointJacobian Pi[64];
	EccPointJacobian P_2 = EccPointAddJacobian(P, P, C);
	for (int j = 1; j < (int)pow(2, w); j = j + 2) {
		if (j == 1) Pi[j] = P;
		else
			Pi[j] = EccPointAddJacobian(Pi[j - 2], P_2, C);
	}
	int i = 0;
	int NAFW[1025] = { 0 };
	BIGNUM w2((int)pow(2, w));
	while (k1 >= 1) {
		if (k1 % 2 == 1) {
			NAFW[i] = BIGNUM(k1 % w2).to_int();
			while (NAFW[i] > pow(2, w - 1)) {
				NAFW[i] = NAFW[i] - pow(2, w);
			}
			k1 = k1 - BIGNUM(NAFW[i]);
		}
		else {
			NAFW[i] = 0;
		}
		k1 = k1 / 2;
		i++;
	}
	long t1, t2;//计算运行时间，t1:开始时间,t2:结束时间
	t1 = GetTickCount64();
	for (int j = i - 1; j >= 0; j--) {
		R = EccPointAddJacobian(R, R, C);
		if (NAFW[j] > 0) {
			R = EccPointAddJacobian(R, Pi[NAFW[j]], C);
		}
		if (NAFW[j] < 0) {
			R = EccPointAddJacobian(R, { Pi[-NAFW[j]].x,0 - Pi[-NAFW[j]].y, Pi[-NAFW[j]].z}, C);
		}
	}
	t2 = GetTickCount64();
	cout << "执行时间：" << t2 - t1 << endl;  //程序运行的时间得到的时间单位为毫秒 /1000为秒
	return R;
}

EccPoint EccPointMul4(BIGNUM k, EccPoint P, EccParams C)
{
	return EccPoint();
}

//求最大公约数的同时求 ax + by = gcd(a,b) 的解
//exgcd结束时，x、y就是所求解
BIGNUM exgcd(BIGNUM a, BIGNUM b, BIGNUM &x, BIGNUM &y)
{
	//递归边界
	if (b == 0) {
		x = 1, y = 0;
		return a;
	}
	//递归计算最大公约数gcd
	BIGNUM d = exgcd(b, a % b, y, x);//a,b交换位置，x,y也交换位置
	//递推公式，求解
	y = y - a / b * x;
	return d;
}

BIGNUM Mod_inverse(BIGNUM a, BIGNUM b)
{
	BIGNUM x, y;
	BIGNUM gcd = exgcd(a, b, x, y); //最大公因子可能为1或者-1
	if (gcd == -1) return 0 - x;
	return x;
}
// 输入：a,b,N
// 输出：y=ab (mod N)
// 1.计算
//		R=2^len(N);
//		a'=aR(mod N),b'=bR(mod N),X=a'b'
// 2.调用蒙哥马利约简算法
//		X1=MD(X,R,N) = X/R=abR(mod N)
// 3.再调用蒙哥马利约简算法
//		y=MD(X1,R,N) = X1/R=ab(mod N)
BIGNUM Montgomery_Multiply(BIGNUM a , BIGNUM b, BIGNUM N)
{
	BIGNUM R;
	BIGNUM a_ = a * R % N;
	BIGNUM b_ = b * R % N;
	BIGNUM X = a_ * b_;
	BIGNUM X1 = Montgomery_Reduction(X, R, N);
	BIGNUM y = Montgomery_Reduction(X1, R, N);
	return y;
}
// 已知 a,b,N
// 输入：X，R，N
// 输出 y
// 1.计算
//		N'=-N^-1 (mod R)
//		m=XN'(mod R);
//	2.计算 y=(X+mN)/R: X+mN >> k;
//  3.若y>N,则y=y-N; 
BIGNUM Montgomery_Reduction(BIGNUM X, BIGNUM R, BIGNUM N)
{
	BIGNUM N_;
	BIGNUM m = X * N_ % R;
	BIGNUM y = (X + m * N) / R;//R=2^k y = X + mN >> k
	if (y > N) y = y - N;
	return y;
}

BIGNUM random(BIGNUM n)
{
	return BIGNUM();
}
