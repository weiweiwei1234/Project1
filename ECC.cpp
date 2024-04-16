#include "SM3.h"
#include "EllipticCurve.h"
#include "BigNumber.h"
#include <windows.h>
//椭圆曲线方程：y2 = x3 + ax + b。
int main() {
	BigNumber p("115792089237316195423570985008687907853269984665640564039457584007913129639932");
	BigNumber a("FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFC");
	BigNumber b("28E9FA9E9D9F5E344D5A9E4BCF6509A7F39789F515AB8F92DDBCBD414D940E93");
	BigNumber n("FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFF7203DF6B21C6052B53BBF40939D54123");
	BigNumber Gx("32C4AE2C1F1981195F9904466A39C9948FE30BBFF2660BE1715A4589334C74C7");
	BigNumber Gy("BC3736A2F4F6779C59BDCEE36B692153D0A9877CC62A474002DF32E52139F0A0");
	BigNumber o("0000000000000000000000000000000000000000000000000000000000000000");
	ECPoint G = { Gx,Gy };
	ECParams C = { p,a,b,n,Gx,Gy };
	ECPointStandardProjection G_StandardProjection = AffineToStandardProjection(G);
	ECPointJacobian G_Jacobian = AffineTOJacobian(G);
	cout << "以下为测试内容:" << endl;
	printecparams(C);
	printecpoint(G);
	printecpointStandarProjection(G_StandardProjection);
	printecpointJacobian(G_Jacobian);
	//cout << isinparams(G, C) << endl;
	//cout << isinparamsStandardProjection(G_StandardProjection, C) << endl;
	//cout << isinparamsJacobian(G_Jacobian, C)<< endl;
	BigNumber k("28E9FA9E9D9F5E344D5A9E4BCF6509A7F39789F515AB8F92DDBCBD414D940E93");
	//cout << HexToBin(k.get_value());
	/*cout << "使用朴素点乘:" << endl;
	ECPoint Pk1 = ecpoint_mul_1(k, G, C);
	print_ecpoint(Pk1);*/

	long t1,t2;//计算运行时间，t1:开始时间,t2:结束时间

	cout << "\n加减链方法（二进制）表示：" << endl;
	t1 = GetTickCount64();
	ECPoint Pk2 = ecpointmulBIN(k, G, C);
	t2 = GetTickCount64();
	cout << "执行时间：" << t2 - t1 << endl;  //程序运行的时间得到的时间单位为毫秒 /1000为秒
	printecpoint(Pk2);

	cout << "\n使用NAF方法（仿射坐标）：" << endl;
	t1 = GetTickCount64();
	ECPoint Pk3 = ecpointmulNAF(k, G, C);
	t2 = GetTickCount64();
	cout << "执行时间：" << t2 - t1 << endl;  //程序运行的时间得到的时间单位为毫秒 /1000为秒
	printecpoint(Pk3);

	cout << "\n使用w-NAF方法（仿射坐标）:" << endl;
	ECPoint Pk3_ = ecpointmulW_NAF(k, G, 4, C);
	printecpoint(Pk3_);

	cout << "\n使用NAF方法（标准射影坐标）：" << endl;
	t1 = GetTickCount64();
	ECPointStandardProjection Pk4 = ecpointmulNAFStandardProjection(k, G_StandardProjection, C);
	t2 = GetTickCount64();
	cout << "执行时间：" << t2 - t1 << endl;  //程序运行的时间得到的时间单位为毫秒 /1000为秒
	ECPoint Pk4_ = StandardProjectionToAffine(Pk4, C);
	printecpoint(Pk4_);
	printecpointStandarProjection(Pk4);

	cout << "\n使用NAF方法（Jacobian加重射影坐标）：" << endl;
	t1 = GetTickCount64();
	ECPointJacobian Pk5 = ecpointmulNAKJacobian(k, G_Jacobian, C);
	t2 = GetTickCount64();
	cout << "执行时间：" << t2 - t1 << endl;  //程序运行的时间得到的时间单位为毫秒 /1000为秒
	ECPoint Pk5_ = JacobianToAffine(Pk5, C);
	printecpoint(Pk5_);
	printecpointJacobian(Pk5);

	cout << "以上为测试内容。" << endl;
	return 0;
}