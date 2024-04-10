#include "SM2.h"

void print_ecparams(ECParams params)
{
	cout << "��Բ���ߵ�������" << params.p << endl;
	cout << "��Բ���ߵĲ���a��" << params.a << endl;
	cout << "��Բ���ߵĲ���b��" << params.b << endl;
	cout << "��Բ���ߵĽ�n" << params.n << endl;
	cout << "��Բ���ߵĻ���X���꣺" << params.Gx << endl;
	cout << "��Բ���ߵĻ���Y����" << params.Gy << endl;
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
