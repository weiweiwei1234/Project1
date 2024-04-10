#include "SM3.h"
#include "EllipticCurve.h"
#include "BigNumber.h"
//椭圆曲线方程：y2 = x3 + ax + b。
int main() {
	BigNumber p("FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFF");
	BigNumber a("FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFC");
	BigNumber b("28E9FA9E9D9F5E344D5A9E4BCF6509A7F39789F515AB8F92DDBCBD414D940E93");
	BigNumber n("FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFF7203DF6B21C6052B53BBF40939D54123");
	BigNumber Gx("32C4AE2C1F1981195F9904466A39C9948FE30BBFF2660BE1715A4589334C74C7");
	BigNumber Gy("BC3736A2F4F6779C59BDCEE36B692153D0A9877CC62A474002DF32E52139F0A0");
	BigNumber o("0000000000000000000000000000000000000000000000000000000000000000");
	ECPoint G = { Gx,Gy };
	ECParams C = { p,a,b,n,Gx,Gy };
	cout << "以下为测试内容:" << endl;
	print_ecparams(C);
	print_ecpoint(G);
	if (is_in_params(G, C) == 1) {
		cout << "G在椭圆曲线上" << endl;
	}
	else {
		cout << "G不在椭圆曲线上" << endl;
	}
	BigNumber a_1 = mod_inverse(a, p);
	BigNumber b_1 = mod_inverse(b, p);
	cout << a_1 << endl;
	cout << b_1 << endl;
	cout << a * a_1 % p << endl;
	cout << b * b_1 % p << endl;
	ECPoint P2 = ecpoint_add(G, G, C);
	print_ecpoint(P2);
	if (is_in_params(P2, C) == 1) {
		cout << "P2在椭圆曲线上" << endl;
	}
	else {
		cout << "P2不在椭圆曲线上" << endl;
	}
	ECPoint P3 = ecpoint_add(G, P2, C);
	print_ecpoint(P3);
	if (is_in_params(P3, C) == 1) {
		cout << "P3在椭圆曲线上" << endl;
	}
	else {
		cout << "P3不在椭圆曲线上" << endl;
	}
	ECPoint P4 = ecpoint_add(P2, P2, C);
	print_ecpoint(P4);
	if (is_in_params(P4, C) == 1) {
		cout << "P4在椭圆曲线上" << endl;
	}
	else {
		cout << "P4不在椭圆曲线上" << endl;
	}
	ECPoint P5 = ecpoint_add(P4, G, C);
	print_ecpoint(P5);
	if (is_in_params(P5, C) == 1) {
		cout << "P5在椭圆曲线上" << endl;
	}
	else {
		cout << "P5不在椭圆曲线上" << endl;
	}
	ECPoint P5_ = ecpoint_mul_1(5, G, C);
	print_ecpoint(P5_);
	BigNumber k("F");
	ECPoint P20 = ecpoint_mul_1(k, G, C);
	print_ecpoint(P20);
	cout << "以上为测试内容。" << endl;
	return 0;
}