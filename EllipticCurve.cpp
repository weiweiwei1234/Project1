#include "EllipticCurve.h"
#include "SM3.h"
void print_ecparams(ECParams params)
{
	cout << "��Բ���ߵ�������p��" << params.p << endl;
	cout << "��Բ���ߵĲ���a��" << params.a << endl;
	cout << "��Բ���ߵĲ���b��" << params.b << endl;
	cout << "��Բ���ߵĽ�n:" << params.n << endl;
	cout << "��Բ���ߵĻ���X���꣺" << params.Gx << endl;
	cout << "��Բ���ߵĻ���Y����:" << params.Gy << endl;
}

void print_ecpoint(ECPoint point)
{
	cout << "X���꣺" << point.x << endl;
	cout << "Y���꣺" << point.y << endl;
}

bool is_in_params(ECPoint point, ECParams params)
{
	if ((point.y * point.y) % params.p == (point.x * point.x * point.x + params.a * point.x + params.b) % params.p)
		return true;
	return false;
}

//�� P + Q
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
		cout << "P��Q��������Բ������" << endl;
	BigNumber k;	//б��k
	if (P.x != Q.x) {
		k = (y2 - y1) * mod_inverse(x2 - x1, p);
	}
	else {
		k = (3 * x1 * x1 + a) * mod_inverse(2 * y1 , p);
	}
	x3 = k * k - x1 - x2;
	y3 = k * (x1 - x3) - y1;
	R.x = x3;
	R.y = y3;
	R.x = (x3 % p + p) % p;
	R.y = (y3 % p + p) % p;
	if (is_in_params(R, params) == 1)
		cout << "�����ȷ!" << endl;
	else
		cout << "�������!" << endl;
	return R;
}

//��׼��Ӱ����������
ECPoint ecpoint_add_(ECPoint P, ECPoint Q, ECParams C)
{
	if (P.x == 0 && P.y == 0) return Q;
	if (Q.x == 0 && Q.y == 0) return P;
	ECPoint R;
	BigNumber x1, x2, x3, y1, y2, y3, z1, z2, z3, a, b, p;
	x1 = P.x;
	y1 = P.y;
	z1 = 1;
	x2 = Q.x;
	y2 = Q.y;
	z2 = 1;
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
	//ת��Ϊ��������
	R.x = x3 * mod_inverse(z3, p) % p;
	R.y = y3 * mod_inverse(z3, p) % p;
	R.x = (R.x % p + p) % p;
	R.y = (R.y % p + p) % p;
	if (is_in_params(R, C) == 1)
		cout << "�����ȷ!" << endl;
	else
		cout << "�������!" << endl;
	return R;
}

//�� kP
//���س˷� 1-k �ǳ���
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
//��k�����Ʊ�ʾ
ECPoint ecpoint_mul_BIN(BigNumber k, ECPoint P, ECParams C)
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
			R = ecpoint_add(R, L, C);
		}
		L = ecpoint_add(L, L, C);
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

ECPoint ecpoint_mul_NAF(BigNumber k, ECPoint P, ECParams C)
{
	if (k == 0) return { BigNumber(0), BigNumber(0) };
	if (k == 1) return P;
	ECPoint R = { BigNumber(0),BigNumber(0) };
	ECPoint _P = {P.x,0-P.y}; //-P
	BigNumber k1 = k;
	//����k��NAFֵ
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
	cout << k << "��NAF��ʾΪ��";
	int j;
	for (j = i - 1; j >= 0; j--) {
		cout << NAF_k[j];
	}
	cout << endl;
	for (j = i - 1; j >= 0; j--) {
		R = ecpoint_add(R, R, C);
		if (NAF_k[j] == 1)
			R = ecpoint_add(R, P, C);
		if (NAF_k[j] == -1)
			R = ecpoint_add(R, _P, C);
	}
	return R;
}

//��Ӱ����
ECPoint ecpoint_mul_NAF_(BigNumber k, ECPoint P, ECParams C)
{
	if (k == 0) return { BigNumber(0), BigNumber(0) };
	if (k == 1) return P;
	ECPoint R = { BigNumber(0),BigNumber(0) };
	ECPoint _P = { P.x,0 - P.y }; //-P
	BigNumber k1 = k;
	//����k��NAFֵ
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
	cout << k << "��NAF��ʾΪ��";
	int j;
	for (j = i - 1; j >= 0; j--) {
		cout << NAF_k[j];
	}
	cout << endl;
	for (j = i - 1; j >= 0; j--) {
		R = ecpoint_add(R, R, C);
		if (NAF_k[j] == 1)
			R = ecpoint_add_(R, P, C);
		if (NAF_k[j] == -1)
			R = ecpoint_add_(R, _P, C);
	}
	return R;
}

ECPoint ecpoint_mul_4(BigNumber k, ECPoint P, ECParams C)
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

BigNumber mod_inverse(BigNumber a, BigNumber b)
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
