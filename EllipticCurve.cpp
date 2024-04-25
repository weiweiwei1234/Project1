#include "EllipticCurve.h"
#include "SM3.h"
void printEccParams(EccParams C)
{
	cout << "��Բ���߲�����" << endl;
	cout << "������p��\t" << C.p << endl;
	cout << "ϵ��a��\t\t" << C.a << endl;
	cout << "ϵ��b��\t\t" << C.b << endl;
	cout << "��n��\t\t" << C.n << endl;
	cout << "����X���꣺\t" << C.Gx << endl;
	cout << "����Y���꣺\t" << C.Gy << endl;
}

void printEccPoint(EccPoint point)
{
	cout << "���������ʾ��" << endl;
	cout << "X���꣺\t" << point.x << endl;
	cout << "Y���꣺\t" << point.y << endl;
}

void printEccPointStandarProjection(EccPointStandardProjection point)
{
	cout << "��׼��Ӱ�����ʾ:" << endl;
	cout << "X���꣺\t" << point.x << endl;
	cout << "Y���꣺\t" << point.y << endl;
	cout << "Z���꣺\t" << point.z << endl;
}

void printEccPointJacobian(EccPointJacobian point)
{
	cout << "�ſɱ������ʾ:" << endl;
	cout << "X���꣺\t" << point.x << endl;
	cout << "Y���꣺\t" << point.y << endl;
	cout << "Z���꣺\t" << point.z << endl;
}
//Y^2=X^3+aX+b
bool isinEccParams(EccPoint point, EccParams C)
{
	BIGNUM x, y, a, b, p;
	a = C.a;
	b = C.b;
	p = C.p;
	x = (point.x + p) % p;
	y = (point.y + p) % p;
	BIGNUM Left = y * y % p;
	BIGNUM Right = x * x % p * x % p + a * x % p + b;
	Left = (Left + p) % p;
	Right = (Right + p) % p;
	/*cout << Left << endl;
	cout << Right << endl;*/
	if (Left == Right)
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
	if (y * y % p * z % p == (x * x % p * x % p + a * x % p * z %p * z %p + b * z % p * z % p * z %p) % p) {
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
	if ((y * y) % p == (x * x % p * x + a * x % p * z % p * z % p * z % p * z + b * z %p * z % p  * z % p * z % p * z % p * z % p) % p){
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
	R.x = P.x * Mod_inverse(P.z * P.z % C.p, C.p);
	R.y = P.y * Mod_inverse(P.z * P.z % C.p * P.z % C.p, C.p);
	R.x = (R.x % C.p + C.p) % C.p;
	R.y = (R.y % C.p + C.p) % C.p;
	return R;
}

//�� P + Q
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
		cout << "P��Q��������Բ������" << endl;
	BIGNUM k;	//б��k
	if (P.x != Q.x) {
		k = (y2 - y1) * Mod_inverse(x2 - x1, p) % p;
	}
	else {
		k = (3 * x1 * x1 % p + a) * Mod_inverse(2 * y1 % p, p) % p;
	}
	x3 = k * k % p- x1 - x2;
	y3 = k * (x1 - x3) % p - y1;
	x3 = (x3 % p + p) % p;
	y3 = (y3 % p + p) % p;
	R = { x3,y3 };
	return R;
}

//��׼��Ӱ���������� һ�㲻ʹ��
EccPointStandardProjection EccPointAddStandardProjection(EccPointStandardProjection P, EccPointStandardProjection Q, EccParams C)
{
	if (P.z == 0) return Q;
	if (Q.z == 0) return P;
	if (!isinEccParamsStandardProjection(P, C) || !isinEccParamsStandardProjection(Q, C))
		cout << "P��Q��������Բ������" << endl;
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
		t1 = x1 * z2 % p;
		t2 = x2 * z1 % p;
		t3 = t1 - t2;
		t4 = y1 * z2 % p;
		t5 = y2 * z1 % p;
		t6 = t4 - t5;
		t7 = t1 + t2;
		t8 = z1 * z2 % p;
		t9 = t3 * t3 % p;
		t10 = t3 * t9 % p;
		t11 = t8 * t6 % p * t6 % p - t7 * t9 % p;
		//�õ�R����Ӱ����
		x3 = t3 * t11 % p;
		y3 = t6 * (t9 * t1 % p - t11) % p - t4 * t10 % p;
		z3 = t10 * t8 % p;
	}
	else {
		BIGNUM t1, t2, t3, t4, t5, t6;
		t1 = 3 * x1 % p * x1 % p + a * z1 % p * z1 % p;
		t2 = 2 * y1 % p * z1 % p;
		t3 = y1 * y1 % p;
		t4 = t3 * x1 % p * z1 % p;
		t5 = t2 * t2 % p;
		t6 = t1 * t1 % p - 8 * t4 % p;
		//�õ�R����Ӱ����
		x3 = t2 * t6 % p;
		y3 = t1 * (4 * t4 % p - t6) % p - 2 * t5 % p * t3 % p;
		z3 = t2 * t5 % p;
	}
	x3 = (x3 % p + p) % p;
	y3 = (y3 % p + p) % p;
	z3 = (z3 % p + p) % p;
	R = { x3,y3,z3 };
	return R;
}
//����� �ſɱ�����
EccPointJacobian EccPointAddJacobian(EccPointJacobian P, EccPointJacobian Q, EccParams C)
{
	if (P.z == 0) return Q;
	if (Q.z == 0) return P;
	if (!isinEccParamsJacobian(P, C) || !isinEccParamsJacobian(Q, C))
		cout << "P��Q��������Բ������" << endl;
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
		t1 = x1 * z2 % p * z2 % p;
		t2 = x2 * z1 % p * z1 % p;
		t3 = t1 - t2;
		t4 = y1 * z2 % p * z2 % p * z2 % p;
		t5 = y2 * z1 % p * z1 % p * z1 % p;
		t6 = t4 - t5;
		t7 = t1 + t2;
		t8 = t4 + t5;
		//�õ�R���ſɱ�����
		x3 = t6 * t6 % p - t7 * t3 % p * t3 % p;
		t9 = t7 * t3 % p * t3 % p - 2 * x3 % p;
		y3 = (t9 * t6 % p - t8 * t3 % p * t3 % p * t3 % p) * Mod_inverse(2, p) % p;
		z3 = z1 * z2 % p * t3 % p;
	}
	else {
		BIGNUM t1, t2, t3;
		t1 = 3 * x1 % p * x1 % p + a * z1 % p * z1 % p * z1 % p * z1 % p;
		t2 = 4 * x1 % p * y1 % p * y1 % p;
		t3 = 8 * y1 % p * y1 % p * y1 % p * y1 % p;
		//�õ�R���ſɱ�����
		x3 = t1 * t1 % p - 2 * t2 % p;
		y3 = t1 * (t2 - x3) % p - t3;
		z3 = 2 * y1 % p * z1 % p;
	}
	x3 = (x3 % p + p) % p;
	y3 = (y3 % p + p) % p;
	z3 = (z3 % p + p) % p;
	R = { x3,y3,z3 };
	return R;
}

//�� kP
//ѭ�� 1-k �ǳ���
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
//��k�����Ʊ�ʾ
EccPoint EccPointMulBIN(BIGNUM k, EccPoint P, EccParams C)
{
	if (k == 0) return { BIGNUM(0), BIGNUM(0) };
	if (k == 1) return P;
	EccPoint R = {BIGNUM(0),BIGNUM(0)};
	EccPoint L = P;
	while (k > 0) {
		if (k % 2 == 1) {
			R = EccPointAdd(R, L, C);
		}
		L = EccPointAdd(L, L, C);
		k = k / 2;
	}
	return R;
}
//ʹ��NAF�㷨
//�����ķ����ڱ�ʾ��ΪNAF
//���壺һ��������k�ķ����ڱ�ʾ���Ǳ��ʽ, ����, ����û����������������ki�Ƿ����, NAF�ĳ�����l.
//NAF�����ʣ�
//1 k��Ψһ��NAF, ������NAF(k)��
//2 NAF(k)��k �����д����ű�ʾ�о������ٵķ���λ��
//3 NAF(k)�ĳ�������k�Ķ����Ʊ�ʾ�ĳ��ȴ�1.
//4 ��NAF(k)�ĳ�����l, ��2l / 3 < k < 2(l + 1) / 3.
//5 ���г���Ϊl��NAF�з������ֵ�ƽ���ܶ�ԼΪ1 / 3.

EccPoint EccPointMulNAF(BIGNUM k, EccPoint P, EccParams C)
{
	if (k == 0) return { BIGNUM(0), BIGNUM(0) };
	if (k == 1) return P;
	EccPoint R = { BIGNUM(0),BIGNUM(0) };
	BIGNUM k1 = k;
	//����k��NAFֵ
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

	long t1, t2;
	t1 = GetTickCount64();
	//��w����Ԥ����� ����iP i=1,3,5,...,2^(w-1)-1
	EccPoint Pi[1024];
	EccPoint P_2 = EccPointAdd(P, P, C);
	for (int j = 1; j < (int)pow(2, w); j = j + 2) {
		if (j == 1) Pi[j] = P;
		else
			Pi[j] = EccPointAdd(Pi[j - 2], P_2, C);
	}
	t2 = GetTickCount64();
	cout << "����Ԥ������ִ��ʱ�䣺" << t2 - t1 << endl;  

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
		R = EccPointAdd(R, R, C);
		if (NAFW[j] > 0) {
			R = EccPointAdd(R, Pi[NAFW[j]], C);
		}
		if (NAFW[j] < 0) {
			R = EccPointAdd(R, {Pi[-NAFW[j]].x,0-Pi[-NAFW[j]].y}, C);
		}
	}
	return R;
}
EccPoint Mul_Montgomery_ladder(BIGNUM k, EccPoint P, EccParams C)
{
	EccPoint R = { BIGNUM(0),BIGNUM(0) };
	EccPoint S = P;
	string k_str = k.HexToBin();
	cout << k_str << endl;
	for (int i = k_str.length() - 1; i >= 0; i--) {
		if (k_str[i] = '1') {
			R = EccPointAdd(R, S, C);
			S = EccPointAdd(S, S, C);
		}
		else {
			S = EccPointAdd(R, S, C);
			R = EccPointAdd(R, R, C);
		}
	}
	return R;
}
//��׼��Ӱ�����µĶ����Ʊ�ʾ
EccPointStandardProjection EccPointMulBINStandardProjection(BIGNUM k, EccPointStandardProjection P, EccParams C)
{
	if (k == 0) return { BIGNUM(0), BIGNUM(0),BIGNUM{0} };
	if (k == 1) return P;
	EccPointStandardProjection R = { BIGNUM(0),BIGNUM(0) ,BIGNUM{0} };
	EccPointStandardProjection L = P;
	while (k > 0) {
		if (k % 2 == 1) {
			R = EccPointAddStandardProjection(R, L, C);
		}
		L = EccPointAddStandardProjection(L, L, C);
		k = k / 2;
	}
	return R;
}

//��׼��Ӱ�����µ�NAF��ʾ
EccPointStandardProjection EccPointMulNAFStandardProjection(BIGNUM k, EccPointStandardProjection P, EccParams C)
{
	if (k == 0) return { BIGNUM(0),BIGNUM(1),BIGNUM(0) };
	if (k == 1) return P;
	EccPointStandardProjection R = { BIGNUM(0),BIGNUM(1),BIGNUM(0) };
	EccPointStandardProjection _P = { P.x,0 - P.y ,P.z}; //-P
	BIGNUM k1 = k;
	//����k��NAFֵ
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

EccPointStandardProjection EccPointMul_W_NAF_StandardProjection(BIGNUM k, EccPointStandardProjection P, int w, EccParams C)
{
	if (k == 0) return { BIGNUM(0),BIGNUM(0),BIGNUM(0) };
	if (k == 1) return P;
	EccPointStandardProjection R{ BIGNUM(0), BIGNUM(0),BIGNUM(0) };
	BIGNUM k1 = k;

	long t1, t2;
	t1 = GetTickCount64();
	//��w����Ԥ����� ����iP i=1,3,5,...,2^(w-1)-1
	EccPointStandardProjection Pi[1024];
	EccPointStandardProjection P_2 = EccPointAddStandardProjection(P, P, C);
	for (int j = 1; j < (int)pow(2, w); j = j + 2) {
		if (j == 1) Pi[j] = P;
		else
			Pi[j] = EccPointAddStandardProjection(Pi[j - 2], P_2, C);
	}
	t2 = GetTickCount64();
	cout << "����Ԥ������ִ��ʱ�䣺" << t2 - t1 << endl;

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
		R = EccPointAddStandardProjection(R, R, C);
		if (NAFW[j] > 0) {
			R = EccPointAddStandardProjection(R, Pi[NAFW[j]], C);
		}
		if (NAFW[j] < 0) {
			R = EccPointAddStandardProjection(R, { Pi[-NAFW[j]].x,0 - Pi[-NAFW[j]].y, Pi[-NAFW[j]].z }, C);
		}
	}
	return R;
}

EccPointJacobian EccPointMulBINJacobian(BIGNUM k, EccPointJacobian P, EccParams C)
{
	if (k == 0) return { BIGNUM(1), BIGNUM(1),BIGNUM{0} };
	if (k == 1) return P;
	EccPointJacobian R = { BIGNUM(1),BIGNUM(1) ,BIGNUM{0} };
	EccPointJacobian L = P;
	while (k > 0) {
		if (k % 2 == 1) {
			R = EccPointAddJacobian(R, L, C);
		}
		L = EccPointAddJacobian(L, L, C);
		k = k / 2;
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
	//����k��NAFֵ
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

	long t1, t2;
	t1 = GetTickCount64();
	//��w����Ԥ����� ����iP i=1,3,5,...,2^(w-1)-1
	EccPointJacobian Pi[1024];
	EccPointJacobian P_2 = EccPointAddJacobian(P, P, C);
	for (int j = 1; j < (int)pow(2, w); j = j + 2) {
		if (j == 1) Pi[j] = P;
		else
			Pi[j] = EccPointAddJacobian(Pi[j - 2], P_2, C);
	}
	t2 = GetTickCount64();
	cout << "����Ԥ������ִ��ʱ�䣺" << t2 - t1 << endl; 

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
		R = EccPointAddJacobian(R, R, C);
		if (NAFW[j] > 0) {
			R = EccPointAddJacobian(R, Pi[NAFW[j]], C);
		}
		if (NAFW[j] < 0) {
			R = EccPointAddJacobian(R, { Pi[-NAFW[j]].x,0 - Pi[-NAFW[j]].y, Pi[-NAFW[j]].z}, C);
		}
	}
	return R;
}

EccPoint EccPointMul4(BIGNUM k, EccPoint P, EccParams C)
{
	return EccPoint();
}

//�����Լ����ͬʱ�� ax + by = gcd(a,b) �Ľ�
//exgcd����ʱ��x��y���������
BIGNUM exgcd(BIGNUM a, BIGNUM b, BIGNUM &x, BIGNUM &y)
{
	//�ݹ�߽�
	if (b == 0) {
		x = 1, y = 0;
		return a;
	}
	//�ݹ�������Լ��gcd
	BIGNUM d = exgcd(b, a % b, y, x);//a,b����λ�ã�x,yҲ����λ��
	//���ƹ�ʽ�����
	y = y - a / b * x;
	return d;
}

BIGNUM Mod_inverse(BIGNUM a, BIGNUM b)
{
	if (a == 2) return BIGNUM("7FFFFFFF7FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF800000008000000000000000");
	Sleep(100);
	BIGNUM x, y;
	BIGNUM gcd = exgcd(a, b, x, y); //������ӿ���Ϊ1����-1
	if (gcd == -1) return 0 - x;
	return (x % b + b) % b;
}
// ���룺a,b,N
// �����y=ab (mod N)
// 1.����
//		R=2^len(N);
//		a'=aR(mod N),b'=bR(mod N),X=a'b'
// 2.�����ɸ�����Լ���㷨
//		X1=MD(X,R,N) = X/R=abR(mod N)
// 3.�ٵ����ɸ�����Լ���㷨
//		y=MD(X1,R,N) = X1/R=ab(mod N)
BIGNUM Montgomery_Multiply(BIGNUM a , BIGNUM b, BIGNUM N)
{
	//BIGNUM R = 2 << N.bitlen();
	BIGNUM a_ = (a << N.bitlen()) % N;
	BIGNUM b_ = (b << N.bitlen()) % N;
	BIGNUM X = a_ * b_;
	BIGNUM X1 = Montgomery_Reduction(X, N);
	BIGNUM y = Montgomery_Reduction(X1, N);
	return y;
}
// ��֪ a,b,N
// ���룺X��R��N
// ��� y
// 1.����
//		N'=-N^-1 (mod R)
//		m=XN'(mod R);
//	2.���� y=(X+mN)/R: X+mN >> k;
//  3.��y>N,��y=y-N; 
BIGNUM Montgomery_Reduction(BIGNUM X, BIGNUM N)
{
	BIGNUM R = BIGNUM(2) << N.bitlen();
	BIGNUM N_, R_;
	BIGNUM d = exgcd(R, N, R_, N_);
	if (d == 1) N_ = BIGNUM(0) - N_;
	else N_ = N_;
	BIGNUM m = X * N_ % R;
	BIGNUM y = (X + m * N) >> N.bitlen();//R=2^k y = X + mN >> k
	if (y > N) y = y - N;
	return y;
}

BIGNUM random(int len,BIGNUM N)
{
	string res = "";
	srand(time(NULL));
	char hex[16] = { '0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F' };
	for (int i = 1; i <= len; i++) {
		int x = rand() % 16;
		res += hex[x];
	}
	if (BIGNUM(res) > N) return random(len, N);
	return BIGNUM(res);
}
