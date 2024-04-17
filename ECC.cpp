#include "SM3.h"
#include "EllipticCurve.h"
#include "BIGNUM.h"
#include <windows.h>
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
	cout << "以下为测试内容:" << endl;
	printEccParams(C);
	printEccPoint(G);
	printEccPointJacobian(G_Jacobian);
	BIGNUM k("ABCDEF");

	long t1,t2;//计算运行时间，t1:开始时间,t2:结束时间

	cout << "\n二进制表示：" << endl;
	t1 = GetTickCount64();
	EccPoint Pk2 = EccPointmulBIN(k, G, C);
	t2 = GetTickCount64();
	cout << "执行时间：" << t2 - t1 << endl;  //程序运行的时间得到的时间单位为毫秒 /1000为秒
	printEccPoint(Pk2);

	cout << "\n使用NAF方法（仿射坐标）：" << endl;
	t1 = GetTickCount64();
	EccPoint Pk3 = EccPointmulNAF(k, G, C);
	t2 = GetTickCount64();
	cout << "执行时间：" << t2 - t1 << endl;  //程序运行的时间得到的时间单位为毫秒 /1000为秒
	printEccPoint(Pk3);

	cout << "\n使用w-NAF方法（仿射坐标）:" << endl;
	EccPoint Pk3_ = EccPointmulW_NAF(k, G, 4, C);
	printEccPoint(Pk3_);

	cout << "\n使用NAF方法（Jacobian加重射影坐标）：" << endl;
	t1 = GetTickCount64();
	EccPointJacobian Pk5 = EccPointmulNAKJacobian(k, G_Jacobian, C);
	t2 = GetTickCount64();
	cout << "执行时间：" << t2 - t1 << endl;  //程序运行的时间得到的时间单位为毫秒 /1000为秒
	EccPoint Pk5_ = JacobianToAffine(Pk5, C);
	printEccPoint(Pk5_);
	printEccPointJacobian(Pk5);

	cout << "以上为测试内容。" << endl;
	return 0;
}