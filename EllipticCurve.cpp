#include "EllipticCurve.h"

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
		cout << "�����ȷ!" << endl;
	else
		cout << "�������!" << endl;
	return R;
}

//�� kP
ECPoint ecpoint_mul_1(BigNumber k, ECPoint P, ECParams params)
{
	//���س˷�1-k �ǳ���
	ECPoint R = P;
	BigNumber p = params.p;
	k = (k % p + p) % p; //ת��kΪ��ֵ
	k = k - 1;
	while (k > 0) {
		R = ecpoint_add(R, P,params);
		k = k - 1;
	}
	return R;
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
