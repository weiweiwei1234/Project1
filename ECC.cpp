#include "SM3.h"
#include "EllipticCurve.h"
#include "BigNumber.h"
#include <windows.h>
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
	printecparams(C);
	printecpoint(G);
	if (isinparams(G, C) == 1) {
		cout << "G在椭圆曲线上\n" << endl;
	}
	else {
		cout << "G不在椭圆曲线上\n" << endl;
	}
	ECPoint P2 = ecpointadd(G, G, C);
	printecpoint(P2);
	if (isinparams(P2, C) == 1) {
		cout << "P2在椭圆曲线上\n" << endl;
	}
	else {
		cout << "P2不在椭圆曲线上\n" << endl;
	}
	BigNumber k("FFFFFFE");

	/*cout << "使用朴素点乘:" << endl;
	ECPoint Pk1 = ecpoint_mul_1(k, G, C);
	print_ecpoint(Pk1);*/

	long t1,t2;//计算运行时间，t1:开始时间,t2:结束时间

	cout << "将k二进制表示:" << endl;
	t1 = GetTickCount64();
	ECPoint Pk2 = ecpointmulBIN(k, G, C);
	t2 = GetTickCount64();
	cout << "执行时间：" << t2 - t1 << endl;  //程序运行的时间得到的时间单位为毫秒 /1000为秒
	printecpoint(Pk2);

	cout << "使用NAF乘法(仿射坐标）:" << endl;
	t1 = GetTickCount64();
	ECPoint Pk3 = ecpointmulNAF(k, G, C);
	t2 = GetTickCount64();
	cout << "执行时间：" << t2 - t1 << endl;  //程序运行的时间得到的时间单位为毫秒 /1000为秒
	printecpoint(Pk3);

	cout << "使用NAF乘法(射影坐标）:" << endl;
	ECPointStandardProjection G_ECPointStandardProjection = StandardProjectionToAffine(G);
	t1 = GetTickCount64();
	ECPointStandardProjection Pk4 = ecpointmulNAFStandardProjection(k, G_ECPointStandardProjection, C);
	t2 = GetTickCount64();
	cout << "执行时间：" << t2 - t1 << endl;  //程序运行的时间得到的时间单位为毫秒 /1000为秒
	ECPoint Pk4_ = AffineToStandardProjection(Pk4, C);
	printecpoint(Pk4_);
	printecpointStandarProjection(Pk4);

	cout << "以上为测试内容。" << endl;
	return 0;
}