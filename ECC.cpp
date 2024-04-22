#include "SM3.h"
#include "EllipticCurve.h"
#include "BIGNUM.h"
#include "test1.h"
#include <windows.h>
#include <random>
#include <sstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <cctype>
//椭圆曲线方程：y2 = x3 + ax + b。
int main() {
	BIGNUM p("FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFF");
	BIGNUM a("FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFC");
	BIGNUM b("28E9FA9E9D9F5E344D5A9E4BCF6509A7F39789F515AB8F92DDBCBD414D940E93");
	BIGNUM n("FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFF7203DF6B21C6052B53BBF40939D54123");
	BIGNUM Gx("32C4AE2C1F1981195F9904466A39C9948FE30BBFF2660BE1715A4589334C74C7");
	BIGNUM Gy("BC3736A2F4F6779C59BDCEE36B692153D0A9877CC62A474002DF32E52139F0A0");
	EccPoint G = { Gx,Gy };
	EccParams C = { p,a,b,n,Gx,Gy };
	EccPointJacobian G_Jacobian = AffineTOJacobian(G);
	long t1, t2;//计算运行时间，t1:开始时间,t2:结束时间
	printEccParams(C);
	printEccPoint(G);
	printEccPointJacobian(G_Jacobian);

	/*cout << "a = " << a << endl;
	cout << "b = " << b << endl;
	cout << "p = " << p << endl;
	cout << "a + b = " << a + b << endl;
	cout << "a - b = " << a - b << endl;
	cout << "a * b = " << a * b << endl;
	cout << "a * b % p = " << a * b % p << endl;
	cout << "a << 64 = " << (a << 64) << endl;
	cout << "a >> 64 = " << (a >> 64) << endl;*/
	test1();
	/*cout << random(64) << endl;
    BIGNUM a1("ABCD");
	cout << "2 mod p的逆元:" << Mod_inverse(2, p) << endl;
	cout << "2 * Mod_inverse(2,p) mod p = " << 2 * Mod_inverse(2, p) % p << endl;
	BIGNUM x = Mod_inverse(a, p);
	cout << "a mod p的逆元x：" << x << endl;
	cout << "a * x mod p = " << a * x % p << endl;
	cout << a.bitlen() << endl;
	cout << (a1 << 24 ) << endl;
	cout << a * b << endl;
	cout << a * b % p << endl;
	cout << Montgomery_Multiply(a, a, p) << endl;*/

	//string str;
	//cout << "请输入字符串：" << endl;
	//getline(cin, str);
	//string hash = SM3Hash(str);
	//string hex_str;
	////将含有中文的字符串编码为16进制字符串形式
	//stringstream ss;
	//for (char c : str) {
	//	ss << std::hex << (int)(unsigned char)c;
	//}
	//hex_str = ss.str();
	//cout << hex_str << endl;
	//transform(hex_str.begin(), hex_str.end(), hex_str.begin(), ::toupper);
	//cout << hex_str << endl;
	//cout << hash << endl;
	return 0;
}