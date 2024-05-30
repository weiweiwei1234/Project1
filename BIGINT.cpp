#include"BIGINT.h"

BIGINT::BIGINT()
{
	sgn = true;
	for (int i = 0; i < 64; i++) {
		data[i] = 0;
		p[i] = 0;
	}
}

BIGINT::BIGINT(long num,long mod)
{
	if (mod < 0) {
		cout << "错误，模数应为正整数！" << endl;
	}
	while (num<0)
	{
		num += mod;
	}
	sgn = true;
	m = 0;
	while (pow(2, m) < num) m++;
	k = m / 8 + (m % 8) > 0;
	int i = 0;
	for (int i = 0; i < k; i++) {
		data[i] = num % 256;
		num /= 256;
	}
	int m_p = 0;
	while (pow(2, m_p) < mod) m_p++;
	int k_p = m_p / 8 + (m_p % 8) > 0;

}
