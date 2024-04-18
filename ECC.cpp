#include "SM3.h"
#include "EllipticCurve.h"
#include "BIGNUM.h"
#include "test1.h"
#include <windows.h>
#include <random>
#include <sstream>
#include <iomanip>
#include <limits>
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

	BIGNUM a1("ABCD");
	cout << a.bitlen() << endl;
	cout << (a1 << 24 ) << endl;
	cout << a * a << endl;
	cout << a * a % p << endl;
	cout << Montgomery_Multiply(a, a, p) << endl;
	//test1();
	return 0;
}