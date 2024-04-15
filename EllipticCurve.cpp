#include "EllipticCurve.h"
#include "SM3.h"
void printecparams(ECParams C)
{
	cout << "��Բ���߲�����" << endl;
	cout << "������p��\t" << C.p << endl;
	cout << "ϵ��a��\t\t" << C.a << endl;
	cout << "ϵ��b��\t\t" << C.b << endl;
	cout << "��n��\t\t" << C.n << endl;
	cout << "����X���꣺\t" << C.Gx << endl;
	cout << "����Y���꣺\t" << C.Gy << endl;
}

void printecpoint(ECPoint point)
{
	cout << "���������ʾ��" << endl;
	cout << "X���꣺\t" << point.x << endl;
	cout << "Y���꣺\t" << point.y << endl;
}

void printecpointStandarProjection(ECPointStandardProjection point)
{
	cout << "��׼��Ӱ�����ʾ:" << endl;
	cout << "X���꣺\t" << point.x << endl;
	cout << "Y���꣺\t" << point.y << endl;
	cout << "Z���꣺\t" << point.z << endl;
}

void printecpointJacobian(ECPointJacobian point)
{
	cout << "�ſɱ������ʾ:" << endl;
	cout << "X���꣺\t" << point.x << endl;
	cout << "Y���꣺\t" << point.y << endl;
	cout << "Z���꣺\t" << point.z << endl;
}
//Y^2=X^3+aX+b
bool isinparams(ECPoint point, ECParams C)
{
	BigNumber x, y, a, b, p;
	x = point.x;
	y = point.y;
	a = C.a;
	b = C.b;
	p = C.p;
	if ((y * y) % p == (x * x * x + a * x + b) % p)
		return true;
	return false;
}
bool isinparamsStandardProjection(ECPointStandardProjection point, ECParams C)
{
	BigNumber x, y, z, a, b, p;
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
bool isinparamsJacobian(ECPointJacobian point, ECParams C)
{
	BigNumber x, y, z, a, b, p;
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

ECPointStandardProjection AffineToStandardProjection(ECPoint P)
{
	ECPointStandardProjection R = { P.x,P.y,1 };
	return R;
}

ECPoint StandardProjectionToAffine(ECPointStandardProjection P,ECParams C)
{
	ECPoint R;
	R.x = P.x * modinverse(P.z, C.p);
	R.y = P.y * modinverse(P.z, C.p);
	R.x = (R.x % C.p + C.p) % C.p;
	R.y = (R.y % C.p + C.p) % C.p;
	return R;
}

ECPointJacobian AffineTOJacobian(ECPoint P)
{
	ECPointJacobian R = { P.x,P.y,1 };
	return R;
}

ECPoint JacobianToAffine(ECPointJacobian P, ECParams C)
{
	ECPoint R;
	R.x = P.x * modinverse(P.z * P.z, C.p);
	R.y = P.y * modinverse(P.z * P.z * P.z, C.p);
	R.x = (R.x % C.p + C.p) % C.p;
	R.y = (R.y % C.p + C.p) % C.p;
	return R;
}

//�� P + Q
ECPoint ecpointadd(ECPoint P, ECPoint Q, ECParams C)
{
	if (P.x == 0 && P.y == 0) return Q;
	if (Q.x == 0 && Q.y == 0) return P;
	ECPoint R; //R=P+Q
	BigNumber x1, x2, x3, y1, y2, y3, a, b, p;
	x1 = P.x;
	y1 = P.y;
	x2 = Q.x;
	y2 = Q.y;
	a = C.a;
	b = C.b;
	p = C.p;
	if (!isinparams(P, C) || !isinparams(Q, C))
		cout << "P��Q��������Բ������" << endl;
	BigNumber k;	//б��k
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
	//if (isinparams(R, C) == 1) {
	//	cout << "�����ȷ" << endl;
	//}
	//else {
	//	cout << "�������" << endl;
	//}
	return R;
}

//��׼��Ӱ����������
ECPointStandardProjection ecpointaddStandardProjection(ECPointStandardProjection P, ECPointStandardProjection Q, ECParams C)
{
	if (P.z == 0) return Q;
	if (Q.z == 0) return P;
	if (!isinparamsStandardProjection(P, C) || !isinparamsStandardProjection(Q, C))
		cout << "P��Q��������Բ������" << endl;
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
	if (P.x != Q.x && P.y != Q.y && P.z != Q.z) {
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
		//�õ�R����Ӱ����
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
		//�õ�R����Ӱ����
		x3 = t2 * t6;
		y3 = t1 * (4 * t4 - t6) - 2 * t5 * t3;
		z3 = t2 * t5;
	}
	x3 = (x3 % p + p) % p;
	y3 = (y3 % p + p) % p;
	z3 = (z3 % p + p) % p;
	R = { x3,y3,z3 };
	//if (isinparamsStandardProjection(R, C) == 1) {
	//	cout << "�����ȷ" << endl;
	//}
	//else {
	//	cout << "�������" << endl;
	//}
	return R;
}
//����� �ſɱ�����
ECPointJacobian ecpointaddJacobian(ECPointJacobian P, ECPointJacobian Q, ECParams C)
{
	if (P.z == 0) return Q;
	if (Q.z == 0) return P;
	if (!isinparamsJacobian(P, C) || !isinparamsJacobian(Q, C))
		cout << "P��Q��������Բ������" << endl;
	ECPointJacobian R;
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
	if (P.x != Q.x && P.y != Q.y && P.z != Q.z) {
		BigNumber t1, t2, t3, t4, t5, t6, t7, t8, t9;
		t1 = x1 * z2 * z2;
		t2 = x2 * z1 * z1;
		t3 = t1 - t2;
		t4 = y1 * z2 * z2 * z2;
		t5 = y2 * z1 * z1 * z1;
		t6 = t4 - t5;
		t7 = t1 + t2;
		t8 = t4 + t5;
		//�õ�R���ſɱ�����
		x3 = t6 * t6 - t7 * t3 * t3;
		t9 = t7 * t3 * t3 - 2 * x3;
		y3 = (t9 * t6 - t8 * t3 * t3 * t3) * modinverse(2, p);
		z3 = z1 * z2 * t3;
	}
	else {
		BigNumber t1, t2, t3;
		t1 = 3 * x1 * x1 + a * z1 * z1 * z1 * z1;
		t2 = 4 * x1 * y1 * y1;
		t3 = 8 * y1 * y1 * y1 * y1;
		//�õ�R���ſɱ�����
		x3 = t1 * t1 - 2 * t2;
		y3 = t1 * (t2 - x3) - t3;
		z3 = 2 * y1 * z1;
	}
	x3 = (x3 % p + p) % p;
	y3 = (y3 % p + p) % p;
	z3 = (z3 % p + p) % p;
	R = { x3,y3,z3 };
	//if (isinparamsJacobian(R, C) == 1) {
	//	cout << "�����ȷ" << endl;
	//}
	//else {
	//	cout << "�������" << endl;
	//}
	return R;
}

//�� kP
//���س˷� 1-k �ǳ���
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
//��k�����Ʊ�ʾ
ECPoint ecpointmulBIN(BigNumber k, ECPoint P, ECParams C)
{
	if (k == 0) return { BigNumber(0), BigNumber(0) };
	if (k == 1) return P;
	ECPoint R = {BigNumber(0),BigNumber(0)};
	ECPoint L = P;
	BigNumber k1 = k;
	string k1_str = k1.get_value();
	string k1_str_BIN = HexToBin(k1_str);
	cout << k1 << "�Ķ����Ʊ�ʾΪ" << k1_str_BIN << endl;
	while (k1 > 0) {
		if (k1 % 2 == 1) {
			R = ecpointadd(R, L, C);
		}
		L = ecpointadd(L, L, C);
		k1 = k1 / 2;
	}
	return R;
}
//ʹ��NAF�㷨(�ȶ����Ʊ�ʾ�õ�)
//�����ķ����ڱ�ʾ��ΪNAF
//���壺һ��������k�ķ����ڱ�ʾ���Ǳ��ʽ, ����, ����û����������������ki�Ƿ����, NAF�ĳ�����l.
//NAF�����ʣ�
//1 k��Ψһ��NAF, ������NAF(k)��
//2 NAF(k)��k �����д����ű�ʾ�о������ٵķ���λ��
//3 NAF(k)�ĳ�������k�Ķ����Ʊ�ʾ�ĳ��ȴ�1.
//4 ��NAF(k)�ĳ�����l, ��2l / 3 < k < 2(l + 1) / 3.
//5 ���г���Ϊl��NAF�з������ֵ�ƽ���ܶ�ԼΪ1 / 3.

ECPoint ecpointmulNAF(BigNumber k, ECPoint P, ECParams C)
{
	if (k == 0) return { BigNumber(0), BigNumber(0) };
	if (k == 1) return P;
	ECPoint R = { BigNumber(0),BigNumber(0) };
	ECPoint _P = {P.x,0-P.y}; //-P
	BigNumber k1 = k;
	//����k��NAFֵ
	int i = 0;
	int NAF_k[1025];
	BigNumber temp;
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
	cout << k << "\n��NAF��ʾΪ��\n";
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

//��Ӱ����
ECPointStandardProjection ecpointmulNAFStandardProjection(BigNumber k, ECPointStandardProjection P, ECParams C)
{
	if (k == 0) return { BigNumber(0),BigNumber(1),BigNumber(0) };
	if (k == 1) return P;
	ECPointStandardProjection R = { BigNumber(0),BigNumber(1),BigNumber(0) };
	ECPointStandardProjection _P = { P.x,0 - P.y ,P.z}; //-P
	BigNumber k1 = k;
	//����k��NAFֵ
	int i = 0;
	int NAF_k[1025];
	BigNumber temp;
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
	cout << k << "\n��NAF��ʾΪ��\n";
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

ECPointJacobian ecpointmulNAKJacobian(BigNumber k, ECPointJacobian P, ECParams C)
{
	if (k == 0) return { BigNumber(1),BigNumber(1),BigNumber(0) };
	if (k == 1) return P;
	ECPointJacobian R = { BigNumber(1),BigNumber(1),BigNumber(0) };
	ECPointJacobian _P = { P.x,0 - P.y ,P.z }; //-P
	BigNumber k1 = k;
	//����k��NAFֵ
	int i = 0;
	int NAF_k[1025];
	BigNumber temp;
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
	cout << k << "\n��NAF��ʾΪ��\n";
	int j;
	for (j = i - 1; j >= 0; j--) {
		cout << NAF_k[j];
	}
	cout << endl;
	for (j = i - 1; j >= 0; j--) {
		R = ecpointaddJacobian(R, R, C);
		if (NAF_k[j] == 1)
			R = ecpointaddJacobian(R, P, C);
		if (NAF_k[j] == -1)
			R = ecpointaddJacobian(R, _P, C);
	}
	return R;
}

ECPoint ecpointmul4(BigNumber k, ECPoint P, ECParams C)
{
	return ECPoint();
}

//�����Լ����ͬʱ�� ax + by = gcd(a,b) �Ľ�
//exgcd����ʱ��x��y���������
BigNumber exgcd(BigNumber a, BigNumber b, BigNumber& x, BigNumber& y)
{
	//�ݹ�߽�
	if (b == 0) {
		x = 1, y = 0;
		return a;
	}
	//�ݹ�������Լ��gcd
	BigNumber gcd = exgcd(b, a % b, x, y);
	//���ƹ�ʽ�����
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
	//	// ������Լ������ 1����û����Ԫ
	//	std::cerr << "��Ԫ�����ڣ���Ϊ���Լ������1��" << std::endl;
	//	return -1;
	//}
	//x= (x % b + b) % b;	//��С�����⣬Ҳ����aģb����Ԫ
	//cout << a * x % b << endl;
	return x * gcd;	//��a��b���ǻ��ʵ��򷵻�x*gcd��ȷ��������ȷ
}

BigNumber random(BigNumber n)
{
	return BigNumber();
}
