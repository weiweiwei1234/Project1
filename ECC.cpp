#include "ECC.h"

void print_Params(Params params)
{
	cout << "椭圆曲线参数：" << endl;
	cout << "有限域p：\t" << params.p << endl;
	cout << "系数a：\t\t" << params.a << endl;
	cout << "系数b：\t\t" << params.b << endl;
	cout << "阶n：\t\t" << params.n << endl;
	cout << "基点X坐标：\t" << params.Gx << endl;
	cout << "基点Y坐标：\t" << params.Gy << endl;
}

void print_Point(Point point)
{
	cout << "仿射坐标表示：" << endl;
	cout << "X坐标：\t" << point.x << endl;
	cout << "Y坐标：\t" << point.y << endl;
}

void print_SPoint(SPoint point)
{
	cout << "标准射影坐标表示:" << endl;
	cout << "X坐标：\t" << point.x << endl;
	cout << "Y坐标：\t" << point.y << endl;
	cout << "Z坐标：\t" << point.z << endl;
}

void print_JPoint(JPoint point)
{
	cout << "雅可比坐标表示:" << endl;
	cout << "X坐标：\t" << point.x << endl;
	cout << "Y坐标：\t" << point.y << endl;
	cout << "Z坐标：\t" << point.z << endl;
}
//Y^2=X^3+aX+b
bool is_in_Params(Point point, Params params)
{
	BIGNUM x, y, a, b, p;
	x = point.x;
	y = point.y;
	a = params.a;
	b = params.b;
	p = params.p;
	BIGNUM Left = (y * y) %p;
	BIGNUM Right = (x * x * x+ a * x + b) % p;
	Left = (Left + p) % p;
	Right = (Right + p) % p;
	if (Left == Right)
		return true;
	return false;
}
//Y^2*Z=X^3+aXZ^2+bZ^3
bool is_in_Params_S(SPoint point, Params params)
{
	BIGNUM x, y, z, a, b, p;
	x = point.x;
	y = point.y;
	z = point.z;
	a = params.a;
	b = params.b;
	p = params.p;
	BIGNUM Left = (y * y * z) % p;
	BIGNUM Right = (x * x * x + a * x * z * z + b * z * z * z) % p;
	Left = (Left + p) % p;
	Right = (Right + p) % p;
	if (Left == Right)
		return true;
	return false;
}
//Y^2=X^3+aXZ^4+bZ^6
bool is_in_Params_J(JPoint point, Params params)
{
	BIGNUM x, y, z, a, b, p;
	x = point.x;
	y = point.y;
	z = point.z;
	a = params.a;
	b = params.b;
	p = params.p;
	BIGNUM Left = (y * y) % p;
	BIGNUM Right = (x * x * x + a * x * z * z * z * z + b * z * z * z * z * z * z) % p;
	Left = (Left + p) % p;
	Right = (Right + p) % p;
	if (Left == Right)
		return true;
	return false;
}

SPoint Affine_To_StandardProjection(Point point)
{
	SPoint R = { point.x,point.y,1 };//z=1
	return R;
}

Point StandardProjection_To_Affine(SPoint point, Params params)
{
	Point R;
	BIGNUM inv = INVERSE(point.z, params.p);
	R.x = (point.x * inv) % params.p;
	R.y = (point.y * inv) % params.p;
	R.x = (R.x + params.p) % params.p;
	R.y = (R.y + params.p) % params.p;
	return R;
}

JPoint Affine_To_Jacobian(Point point)
{
	JPoint R = { point.x,point.y,1 };
	return R;
}

Point Jacobian_To_Affine(JPoint point, Params params)
{
	Point R;
	BIGNUM inv1 = INVERSE(point.z * point.z % params.p, params.p);
	BIGNUM inv2 = INVERSE(point.z * point.z * point.z % params.p, params.p);
	R.x = (point.x * inv1) % params.p;
	R.y = (point.y * inv2) % params.p;
	R.x = (R.x + params.p) % params.p;
	R.y = (R.y + params.p) % params.p;
	return R;
}

BIGNUM EXGCD(BIGNUM a, BIGNUM b, BIGNUM& x, BIGNUM& y)
{
	if (b == 0) {
		x = 1, y = 0;
		return a;
	}
	BIGNUM d = EXGCD(b, a % b, y, x);
	y = y - a / b * x;
	return d;
}

BIGNUM INVERSE(BIGNUM a, BIGNUM p) //p--b
{
	if (a == 2) return BIGNUM("7FFFFFFF7FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF800000008000000000000000");
	Sleep(100);
	BIGNUM x, y;
	BIGNUM GCD = EXGCD(a, p, x, y);
	if (GCD == -1)  x = 0 - x;
	return (x % p + p) % p;
}

BIGNUM M_Multiply(BIGNUM a, BIGNUM b, BIGNUM N)
{
	//R=2^N.bitlen()
	BIGNUM aR = (a << N.bitlen()) % N;
	BIGNUM bR = (b << N.bitlen()) % N;
	BIGNUM X = aR * bR;
	BIGNUM X1 = M_Reduction(X, N);
	BIGNUM y = M_Reduction(X, N);
	return y;
}

BIGNUM M_Reduction(BIGNUM X, BIGNUM N)
{
		/*
	N=FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFF
	R=10000000000000000000000000000000000000000000000000000000000000000
	N_=-3FFFFFFFE00000001FFFFFFFF00000000FFFFFFFEFFFFFFFFFFFFFFFF
	R_=3FFFFFFFA00000003FFFFFFFD00000001FFFFFFFA00000006FFFFFFFB
	*/
	/*
	* 	BIGNUM N_, R_;
	BIGNUM d = exgcd(R, N, R_, N_);
	if (d == 1) {
		N_ = BIGNUM(0) - N_;
		R_ = BIGNUM(0) - R_;
	}
	*/
	//BIGNUM N("FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFF");
	BIGNUM R("10000000000000000000000000000000000000000000000000000000000000000");
	BIGNUM N_("-3FFFFFFFE00000001FFFFFFFF00000000FFFFFFFEFFFFFFFFFFFFFFFF");
	BIGNUM R_("3FFFFFFFA00000003FFFFFFFD00000001FFFFFFFA00000006FFFFFFFB");
	BIGNUM t = X * N_;
	BIGNUM m = t - ((t >> 256) << 256);
	BIGNUM y = (X + m * N) >> 256;
	if (y > N) y = y - N;
	return y;
}

BIGNUM Random(int d, BIGNUM N)
{
	string res = "";
	srand(time(NULL));
	char hex[16] = { '0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F' };
	for (int i = 1; i <= d; i++) {
		int x = rand() % 16;
		res += hex[x];
	}
	if (BIGNUM(res) > N) return Random(d, N);
	return BIGNUM(res);
}

Point Point_Add(Point P, Point Q, Params C)
{
	if (P.x == 0 && P.y == 0) return Q;
	if (Q.x == 0 && Q.y == 0) return P;
	Point R;
	BIGNUM x1, x2, x3, y1, y2, y3, a, b, p;
	x1 = P.x;
	y1 = P.y;
	x2 = Q.x;
	y2 = Q.y;
	a = C.a;
	b = C.b;
	p = C.p;
	BIGNUM k;
	if (x1 != x2)
		k = (y2 - y1) % p * INVERSE((x2 - x1) % p, p);
	else
		k = ((3 * x1 * x1) % p + a) % p * INVERSE(2 * y1 % p, p);
	x3 = (k * k % p - x1 - x2) % p;
	y3 = (k * (x1 - x3) % p - y1) % p;
	x3 = (x3 + p) % p;
	y3 = (y3 + p) % p;
	R = { x3,y3 };
	return R;
}

Point Point_Add_M(Point PR, Point QR, Params C)
{
	if (PR.x == 0 && PR.y == 0) return QR;
	if (QR.x == 0 && QR.y == 0) return PR;
	Point RR; //RR=PR+QR
	BIGNUM x1R, x2R, x3R, y1R, y2R, y3R, a, b, p;
	x1R = PR.x;
	y1R = PR.y;
	x2R = QR.x;
	y2R = QR.y;
	a = C.a;
	b = C.b;
	p = C.p;
	BIGNUM aRR = (a << 512) % p;
	BIGNUM kR,k2R;
	if (x1R != x2R) {
		kR = (y2R - y1R) * INVERSE(x2R - x1R, p);
		kR = (kR << 256) % p;
		k2R = M_Reduction(kR * kR, p);
	}
	else {
		BIGNUM inv = INVERSE(2 * y1R, p);
		kR = (3 * x1R * x1R + aRR) * inv % p;
		k2R = M_Reduction(kR * kR, p);
	}
	x3R = k2R - x1R - x2R;
	y3R = M_Reduction(kR * (x1R - x3R), p) - y1R;
	x3R = (x3R % p + p) % p;
	y3R = (y3R % p + p) % p;
	RR = { x3R,y3R };
	return RR;
}

SPoint SPoint_Add(SPoint P, SPoint Q, Params C)
{
	if (P.z == 0) return Q;
	if (Q.z == 0) return P;
	SPoint R;
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
	if (x1 != x2 && y1 != y2 && z1 != z2) {
		BIGNUM t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11;
		t1 = x1 * z2 % p;
		t2 = x2 * z1 % p;
		t3 = (t1 - t2) % p;
		t4 = y1 * z2 % p;
		t5 = y2 * z1 % p;
		t6 = (t4 - t5) % p;
		t7 = (t1 + t2) % p;
		t8 = z1 * z2 % p;
		t9 = t3 * t3 % p;
		t10 = t3 * t9 % p;
		t11 = (t8 * t6 % p * t6 % p - t7 * t9 % p) % p;
		//得到R的射影坐标
		x3 = t3 * t11 % p;
		y3 = (t6 * (t9 * t1 % p - t11) % p - t4 * t10 % p) % p;
		z3 = t10 * t8 % p;
	}
	else {
		BIGNUM t1, t2, t3, t4, t5, t6;
		t1 = 3 * x1 * x1 % p + a * z1 % p * z1 % p;
		t2 = 2 * y1 * z1 % p;
		t3 = y1 * y1 % p;
		t4 = t3 * x1 % p * z1 % p;
		t5 = t2 * t2 % p;
		t6 = (t1 * t1 % p - 8 * t4) % p;
		//得到R的射影坐标
		x3 = t2 * t6 % p;
		y3 = t1 * (4 * t4 - t6) % p - 2 * t5 * t3 % p;
		z3 = t2 * t5 % p;
	}
	x3 = (x3 + p) % p;
	y3 = (y3 + p) % p;
	z3 = (z3 + p) % p;
	R = { x3,y3,z3 };
	return R;
}

SPoint SPoint_Add_M(SPoint PR, SPoint QR, Params C)
{
	if (PR.z == 0) return QR;
	if (QR.z == 0) return PR;
	SPoint RR;
	BIGNUM x1R, x2R, x3R, y1R, y2R, y3R, z1R, z2R, z3R, a, b, p;
	x1R = PR.x;
	y1R = PR.y;
	z1R = PR.z;
	x2R = QR.x;
	y2R = QR.y;
	z2R = QR.z;
	a = C.a;
	b = C.b;
	p = C.p;
	if (PR.x != QR.x && PR.y != QR.y && PR.z != QR.z) {
		BIGNUM t1R, t2R, t3R, t4R, t5R, t6R, t7R, t8R, t9R, t10R, t11R;
		t1R = M_Reduction(x1R * z2R, p);
		t2R = M_Reduction(x2R * z1R, p);
		t3R = t1R - t2R;
		t4R = M_Reduction(y1R * z2R, p);
		t5R = M_Reduction(y2R * z1R, p);
		t6R = t4R - t5R;
		t7R = t1R + t2R;
		t8R = M_Reduction(z1R * z2R, p);
		t9R = M_Reduction(t3R * t3R, p);
		t10R = M_Reduction(t3R * t9R, p);
		t11R = M_Reduction(M_Reduction(t8R * t6R, p) * t6R - t7R * t9R, p);
		//得到R的射影坐标
		x3R = M_Reduction(t3R * t11R, p);
		y3R = M_Reduction(t6R * (M_Reduction(t9R * t1R, p) - t11R), p) - M_Reduction(t4R * t10R, p);
		z3R = M_Reduction(t10R * t8R, p);
	}
	else {
		BIGNUM t1R, t2R, t3R, t4R, t5R, t6R;
		t1R = 3 * M_Reduction(x1R * x1R, p) + a * M_Reduction(z1R * z1R, p) % p;
		t2R = 2 * M_Reduction(y1R * z1R, p);
		t3R = M_Reduction(y1R * y1R, p);
		t4R = M_Reduction(M_Reduction(t3R * x1R, p) * z1R, p);
		t5R = M_Reduction(t2R * t2R, p);
		t6R = M_Reduction(t1R * t1R, p) - 8 * t4R;
		//得到R的射影坐标
		x3R = M_Reduction(t2R * t6R, p);
		y3R = M_Reduction(t1R * (4 * t4R - t6R), p) - 2 * M_Reduction(t5R * t3R, p);
		z3R = M_Reduction(t2R * t5R, p);
	}
	x3R = (x3R % p + p) % p;
	y3R = (y3R % p + p) % p;
	z3R = (z3R % p + p) % p;
	RR = { x3R,y3R,z3R };
	return RR;
}

JPoint JPoint_Add(JPoint P, JPoint Q, Params C)
{
	if (P.z == 0) return Q;
	if (Q.z == 0) return P;
	JPoint R;
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
	if (x1 != x2 && y1 != y2 && z1 != z2) {
		BIGNUM t1, t2, t3, t4, t5, t6, t7, t8, t9;
		t1 = x1 * z2 % p * z2 % p;
		t2 = x2 * z1 % p * z1 % p;
		t3 = (t1 - t2) % p;
		t4 = y1 * z2 % p * z2 % p * z2 % p;
		t5 = y2 * z1 % p * z1 % p * z1 % p;
		t6 = (t4 - t5) % p;
		t7 = (t1 + t2) % p;
		t8 = (t4 + t5) % p;
		//得到R的雅可比坐标
		x3 = (t6 * t6 % p - t7 * t3 % p * t3 % p) % p;
		t9 = t7 * t3 * t3 - 2 * x3;
		y3 = (t9 * t6 % p - t8 * t3 % p * t3 % p * t3 % p) * INVERSE(2, p) % p;
		z3 = z1 * z2 % p * t3 % p;
	}
	else {
		BIGNUM t1, t2, t3;
		t1 = 3 * x1 * x1 % p + a * z1 % p * z1 % p * z1 % p * z1 % p;
		t2 = 4 * x1 * y1 % p * y1 % p;
		t3 = 8 * y1 * y1 % p * y1 % p * y1 % p;
		//得到R的雅可比坐标
		x3 = (t1 * t1 % p - 2 * t2) % p;
		y3 = (t1 * (t2 - x3) % p - t3) % p;
		z3 = 2 * y1 % p * z1 % p;
	}
	x3 = (x3 + p) % p;
	y3 = (y3 + p) % p;
	z3 = (z3 + p) % p;
	R = { x3,y3,z3 };
	return R;
}

Point Point_Mul_Bin(BIGNUM k, Point P, Params C)
{
	if (k == 0) return { BIGNUM(0), BIGNUM(0) };
	if (k == 1) return P;
	Point R = { BIGNUM(0),BIGNUM(0) };
	Point L = P;
	while (k > 0) {
		if (k % 2 == 1) {
			R = Point_Add(R, L, C);
		}
		L = Point_Add(L, L, C);
		k = k / 2;
	}
	return R;
}

Point Point_Mul_Bin_M(BIGNUM k, Point P, Params C)
{
	if (k == 0) return { BIGNUM(0), BIGNUM(0) };
	if (k == 1) return P;
	Point RR = { BIGNUM(0),BIGNUM(0) };
	Point LR = { (P.x << 256) % C.p, (P.y << 256) % C.p };
	while (k > 0) {
		if (k % 2 == 1) {
			RR = Point_Add_M(RR, LR, C);
		}
		LR = Point_Add_M(LR, LR, C);
		k = k / 2;
	}
	Point R = { M_Reduction(RR.x, C.p) ,M_Reduction(RR.y, C.p) };
	return R;
}

SPoint SPoint_Mul_Bin(BIGNUM k, SPoint P, Params C)
{
	if (k == 0) return { BIGNUM(0), BIGNUM(0),BIGNUM{0} };
	if (k == 1) return P;
	SPoint R = { BIGNUM(0),BIGNUM(0) ,BIGNUM{0} };
	SPoint L = P;
	while (k > 0) {
		if (k % 2 == 1) {
			R = SPoint_Add(R, L, C);
		}
		L = SPoint_Add(L, L, C);
		k = k / 2;
	}
	return R;
}

SPoint SPoint_Mul_Bin_M(BIGNUM k, SPoint P, Params C)
{
	if (k == 0) return { BIGNUM(0), BIGNUM(0),BIGNUM{0} };
	if (k == 1) return P;
	SPoint RR = { BIGNUM(0),BIGNUM(0) ,BIGNUM{0} };
	SPoint LR = { (P.x << 256) % C.p, (P.y << 256) % C.p,(P.z << 256) % C.p };
	while (k > 0) {
		if (k % 2 == 1) {
			RR = SPoint_Add_M(RR, LR, C);
		}
		LR = SPoint_Add_M(LR, LR, C);
		k = k / 2;
	}
	SPoint R = { M_Reduction(RR.x, C.p) ,M_Reduction(RR.y, C.p) ,M_Reduction(RR.z,C.p) };
	return R;
}

JPoint JPoint_Mul_Bin(BIGNUM k, JPoint P, Params C)
{
	if (k == 0) return { BIGNUM(1), BIGNUM(1),BIGNUM{0} };
	if (k == 1) return P;
	JPoint R = { BIGNUM(1),BIGNUM(1) ,BIGNUM{0} };
	JPoint L = P;
	while (k > 0) {
		if (k % 2 == 1) {
			R = JPoint_Add(R, L, C);
		}
		L = JPoint_Add(L, L, C);
		k = k / 2;
	}
	return R;
}

Point Point_Mul_NAF(BIGNUM k, Point P, Params C)
{
	if (k == 0) return { BIGNUM(0), BIGNUM(0) };
	if (k == 1) return P;
	Point R = { BIGNUM(0),BIGNUM(0) };
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
		R = Point_Add(R, R, C);
		if (NAF_k[j] == 1)
			R = Point_Add(R, P, C);
		if (NAF_k[j] == -1)
			R = Point_Add(R, { P.x,0 - P.y }, C);
	}
	return R;
}

Point Point_Mul_NAF_M(BIGNUM k, Point P, Params C)
{
	if (k == 0) return { BIGNUM(0), BIGNUM(0) };
	if (k == 1) return P;
	Point RR = { BIGNUM(0),BIGNUM(0) };
	Point PR = { (P.x << 256) % C.p, (P.y << 256) % C.p };
	Point _PR = { PR.x,0 - PR.y };
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
		RR = Point_Add_M(RR, RR, C);
		if (NAF_k[j] == 1)
			RR = Point_Add_M(RR, PR, C);
		if (NAF_k[j] == -1)
			RR = Point_Add_M(RR, _PR, C);
	}
	Point R = { M_Reduction(RR.x, C.p) ,M_Reduction(RR.y, C.p) };
	return R;
}

SPoint SPoint_Mul_NAF(BIGNUM k, SPoint P, Params C)
{
	if (k == 0) return { BIGNUM(0),BIGNUM(1),BIGNUM(0) };
	if (k == 1) return P;
	SPoint R = { BIGNUM(0),BIGNUM(1),BIGNUM(0) };
	SPoint _P = { P.x,0 - P.y ,P.z }; //-P
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
		R = SPoint_Add(R, R, C);
		if (NAF_k[j] == 1)
			R = SPoint_Add(R, P, C);
		if (NAF_k[j] == -1)
			R = SPoint_Add(R, _P, C);
	}
	return R;
}

SPoint SPoint_Mul_NAF_M(BIGNUM k, SPoint P, Params C)
{
	if (k == 0) return { BIGNUM(0),BIGNUM(1),BIGNUM(0) };
	if (k == 1) return P;
	SPoint RR = { BIGNUM(0),BIGNUM(1),BIGNUM(0) };
	SPoint PR = { (P.x << 256) % C.p, (P.y << 256) % C.p ,(P.z << 256) % C.p };
	SPoint _PR = { PR.x,0 - PR.y,PR.z };
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
		RR = SPoint_Add_M(RR, RR, C);
		if (NAF_k[j] == 1)
			RR = SPoint_Add_M(RR, P, C);
		if (NAF_k[j] == -1)
			RR = SPoint_Add_M(RR, _PR, C);
	}
	SPoint R = { M_Reduction(RR.x, C.p) ,M_Reduction(RR.y, C.p),M_Reduction(RR.z,C.p) };
	return R;
}

JPoint JPoint_Mul_NAF(BIGNUM k, JPoint P, Params C)
{
	if (k == 0) return { BIGNUM(1),BIGNUM(1),BIGNUM(0) };
	if (k == 1) return P;
	JPoint R = { BIGNUM(1),BIGNUM(1),BIGNUM(0) };
	JPoint _P = { P.x,0 - P.y ,P.z }; //-P
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
		R = JPoint_Add(R, R, C);
		if (NAF_k[j] == 1)
			R = JPoint_Add(R, P, C);
		if (NAF_k[j] == -1)
			R = JPoint_Add(R, _P, C);
	}
	return R;
}

Point Point_Mul_wNAF(BIGNUM k, Point P, Params C, int w)
{
	if (k == 0) return { BIGNUM(0),BIGNUM(0) };
	if (k == 1) return P;
	Point R{ BIGNUM(0), BIGNUM(0) };
	BIGNUM k1 = k;

	long t1, t2;
	t1 = GetTickCount64();
	//用w计算预计算表 计算iP i=1,3,5,...,2^(w-1)-1
	Point Pi[1024];
	Point P_2 = Point_Add(P, P, C);
	for (int j = 1; j < (int)pow(2, w); j = j + 2) {
		if (j == 1) Pi[j] = P;
		else
			Pi[j] = Point_Add(Pi[j - 2], P_2, C);
	}
	t2 = GetTickCount64();
	cout << "计算预计算表的执行时间：" << t2 - t1 << endl;

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

	for (int j = i - 1; j >= 0; j--) {
		R = Point_Add(R, R, C);
		if (NAFW[j] > 0) {
			R = Point_Add(R, Pi[NAFW[j]], C);
		}
		if (NAFW[j] < 0) {
			R = Point_Add(R, { Pi[-NAFW[j]].x,0 - Pi[-NAFW[j]].y }, C);
		}
	}
	return R;
}

Point Point_Mul_wNAF_M(BIGNUM k, Point P, Params C, int w)
{
	if (k == 0) return { BIGNUM(0),BIGNUM(0) };
	if (k == 1) return P;
	Point RR{ BIGNUM(0), BIGNUM(0) };
	Point PR = { (P.x << 256) % C.p, (P.y << 256) % C.p };
	BIGNUM k1 = k;

	long t1, t2;
	t1 = GetTickCount64();
	//用w计算预计算表 计算iP i=1,3,5,...,2^(w-1)-1
	Point PRi[1024];
	Point PR_2 = Point_Add_M(P, P, C);
	for (int j = 1; j < (int)pow(2, w); j = j + 2) {
		if (j == 1) PRi[j] = P;
		else
			PRi[j] = Point_Add_M(PRi[j - 2], PR_2, C);
	}
	t2 = GetTickCount64();
	cout << "计算预计算表的执行时间：" << t2 - t1 << endl;

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

	for (int j = i - 1; j >= 0; j--) {
		RR = Point_Add_M(RR, RR, C);
		if (NAFW[j] > 0) {
			RR = Point_Add_M(RR, PRi[NAFW[j]], C);
		}
		if (NAFW[j] < 0) {
			RR = Point_Add_M(RR, { PRi[-NAFW[j]].x,0 - PRi[-NAFW[j]].y }, C);
		}
	}
	Point R = { M_Reduction(RR.x, C.p) ,M_Reduction(RR.y, C.p) };
	return R;
}

SPoint SPoint_Mul_wNAF(BIGNUM k, SPoint P, Params C, int w)
{
	if (k == 0) return { BIGNUM(0),BIGNUM(0),BIGNUM(0) };
	if (k == 1) return P;
	SPoint R{ BIGNUM(0), BIGNUM(0),BIGNUM(0) };
	BIGNUM k1 = k;

	long t1, t2;
	t1 = GetTickCount64();
	//用w计算预计算表 计算iP i=1,3,5,...,2^(w-1)-1
	SPoint Pi[1024];
	SPoint P_2 = SPoint_Add(P, P, C);
	for (int j = 1; j < (int)pow(2, w); j = j + 2) {
		if (j == 1) Pi[j] = P;
		else
			Pi[j] = SPoint_Add(Pi[j - 2], P_2, C);
	}
	t2 = GetTickCount64();
	cout << "计算预计算表的执行时间：" << t2 - t1 << endl;

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
	for (int j = i - 1; j >= 0; j--) {
		R = SPoint_Add(R, R, C);
		if (NAFW[j] > 0) {
			R = SPoint_Add(R, Pi[NAFW[j]], C);
		}
		if (NAFW[j] < 0) {
			R = SPoint_Add(R, { Pi[-NAFW[j]].x,0 - Pi[-NAFW[j]].y, Pi[-NAFW[j]].z }, C);
		}
	}
	return R;
}

SPoint SPoint_Mul_wNAF_M(BIGNUM k, SPoint P, Params C, int w)
{
	if (k == 0) return { BIGNUM(0),BIGNUM(0),BIGNUM(0) };
	if (k == 1) return P;
	SPoint RR{ BIGNUM(0), BIGNUM(0),BIGNUM(0) };
	SPoint PR = { (P.x << 256) % C.p, (P.y << 256) % C.p ,(P.z << 256) % C.p };
	BIGNUM k1 = k;

	long t1, t2;
	t1 = GetTickCount64();
	//用w计算预计算表 计算iP i=1,3,5,...,2^(w-1)-1
	SPoint PRi[1024];
	SPoint PR_2 = SPoint_Add_M(P, P, C);
	for (int j = 1; j < (int)pow(2, w); j = j + 2) {
		if (j == 1) PRi[j] = P;
		else
			PRi[j] = SPoint_Add_M(PRi[j - 2], PR_2, C);
	}
	t2 = GetTickCount64();
	cout << "计算预计算表的执行时间：" << t2 - t1 << endl;

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
	for (int j = i - 1; j >= 0; j--) {
		RR = SPoint_Add_M(RR, RR, C);
		if (NAFW[j] > 0) {
			RR = SPoint_Add_M(RR, PRi[NAFW[j]], C);
		}
		if (NAFW[j] < 0) {
			RR = SPoint_Add_M(RR, { PRi[-NAFW[j]].x,0 - PRi[-NAFW[j]].y, PRi[-NAFW[j]].z }, C);
		}
	}
	SPoint R = { M_Reduction(RR.x, C.p) ,M_Reduction(RR.y, C.p),M_Reduction(RR.z,C.p) };
	return R;
}

JPoint JPoint_Mul_wNAF(BIGNUM k, JPoint P, Params C, int w)
{
	if (k == 0) return { BIGNUM(1),BIGNUM(1),BIGNUM(0) };
	if (k == 1) return P;
	JPoint R{ BIGNUM(1), BIGNUM(1),BIGNUM(0) };
	BIGNUM k1 = k;

	long t1, t2;
	t1 = GetTickCount64();
	//用w计算预计算表 计算iP i=1,3,5,...,2^(w-1)-1
	JPoint Pi[1024];
	JPoint P_2 = JPoint_Add(P, P, C);
	for (int j = 1; j < (int)pow(2, w); j = j + 2) {
		if (j == 1) Pi[j] = P;
		else
			Pi[j] = JPoint_Add(Pi[j - 2], P_2, C);
	}
	t2 = GetTickCount64();
	cout << "计算预计算表的执行时间：" << t2 - t1 << endl;

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
	for (int j = i - 1; j >= 0; j--) {
		R = JPoint_Add(R, R, C);
		if (NAFW[j] > 0) {
			R = JPoint_Add(R, Pi[NAFW[j]], C);
		}
		if (NAFW[j] < 0) {
			R = JPoint_Add(R, { Pi[-NAFW[j]].x,0 - Pi[-NAFW[j]].y, Pi[-NAFW[j]].z }, C);
		}
	}
	return R;
}

