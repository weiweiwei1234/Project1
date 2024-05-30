#include "SM3.h"
#include "EllipticCurve.h"
#include "BIGNUM.h"
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
	EccPointStandardProjection G_SP = AffineToStandardProjection(G);
	EccPointJacobian G_Jacobian = AffineTOJacobian(G);
	long t1, t2;//计算运行时间，t1:开始时间,t2:结束时间
	printEccPoint(EccPointAdd(G, G, C));
	cout << Mod_inverse(Gx, p) << endl;
	cout << (Gx + Gy) % p << endl;
	cout << (Gx - Gy + p) % p << endl;
	cout << Gx * Gy % p << endl;
	cout << Montgomery_Multiply(a, b, p) << endl;
	cout << a * b % p << endl;
	BIGNUM a_R = (a << 256) % p;
	BIGNUM b_R = (b << 256) % p;
	BIGNUM N_R = (p << 256) % p;
	cout << "a的蒙哥马利表示：" << a_R << endl;
	cout << "b的蒙哥马利表示：" << b_R << endl;
	cout << "p的蒙哥马利表示：" << N_R << endl;
	BIGNUM _a = Mod_inverse(a, p);
	BIGNUM _a_R = Mod_inverse(a_R, p);
	cout << _a << endl;
	cout << (_a_R << 256) % p << endl;

	
	BIGNUM k = random(64, n);
	cout << k << endl;

	cout << "\n二进制表示（仿射坐标）：" << endl;
	t1 = GetTickCount64();
	EccPointMulBIN(k, G, C);
	t2 = GetTickCount64();
	cout << "执行时间：" << t2 - t1 << endl;  //程序运行的时间得到的时间单位为毫秒 /1000为秒

	cout << "\n二进制表示（仿射坐标）：" << endl;
	t1 = GetTickCount64();
	EccPointMul_Montomery(k, G, C);
	t2 = GetTickCount64();
	cout << "执行时间：" << t2 - t1 << endl;  //程序运行的时间得到的时间单位为毫秒 /1000为秒
	//cout << a * b << endl;
	
	/*
	cout << Gx * Gy << endl;
	cout << Gx * Gy % p << endl;
	cout << Gx + Gy % p << endl;
	cout << (Gx - Gy + p)%p << endl;
	*/









	/*
	//计算整数基本运算时间
	long time_random = 0;	//生成随机数
	long time_add = 0;		//加法
	long time_sub = 0;		//减法
	long time_mul = 0;		//乘法
	long time_mod = 0;		//求余
	long time_weiyi = 0;	//位移
	long time_inverse = 0;	//求逆元
	BIGNUM A, B, D;
	for (int i = 0; i < 1000; i++) {

		t1 = GetTickCount64();
		A = random(64, n);
		B = random(64, n);
		t2 = GetTickCount64();
		time_random += t2 - t1;

		t1 = GetTickCount64();
		A + B;
		t2 = GetTickCount64();
		time_add += t2 - t1;

		t1 = GetTickCount64();
		A - B;
		t2 = GetTickCount64();
		time_sub += t2 - t1;

		t1 = GetTickCount64();
		D = A * B;
		t2 = GetTickCount64();
		time_mul += t2 - t1;

		t1 = GetTickCount64();
		D % p;
		t2 = GetTickCount64();
		time_mod += t2 - t1;

		t1 = GetTickCount64();
		D >> 512;
		D << 512;
		t2 = GetTickCount64();
		time_weiyi += (t2 - t1) / 2;

		t1 = GetTickCount64();
		//Mod_inverse(A, B);
		t2 = GetTickCount64();
		time_inverse += t2 - t1;
	}
	cout << "1000次随机数生成时间:" << time_random << endl;
	cout << "1000次加法操作时间:" << time_add << endl;
	cout << "1000次减法操作时间:" << time_sub << endl;
	cout << "1000次乘法操作时间:" << time_mul << endl;
	cout << "1000次求余法操作时间:" << time_mod << endl;
	cout << "1000次位移法操作时间:" << time_weiyi << endl;
	cout << "1000次求逆元操作时间:" << time_inverse << endl;
	
	






















	cout << "以下为椭圆曲线点乘运算测试内容:" << endl;

	BIGNUM k = random(64, n);
	cout << k << endl;
	//cout << "montgomery_ladder算法：" << endl;
	//t1 = GetTickCount64();
	//Mul_Montgomery_ladder(k, G, C);
	//t2 = GetTickCount64();
	//cout << "执行时间：" << t2 - t1 << endl;

	cout << "\n二进制表示（仿射坐标）：" << endl;
	t1 = GetTickCount64();
	EccPointMulBIN(k, G, C);
	t2 = GetTickCount64();
	cout << "执行时间：" << t2 - t1 << endl;  //程序运行的时间得到的时间单位为毫秒 /1000为秒

	cout << "\n蒙哥马利优化后的二进制表示（仿射坐标）：" << endl;
	t1 = GetTickCount64();

	t2 = GetTickCount64();
	cout << "执行时间：" << t2 - t1 << endl;  //程序运行的时间得到的时间单位为毫秒 /1000为秒

	cout << "\n使用NAF方法（仿射坐标）：" << endl;
	t1 = GetTickCount64();
	EccPointMulNAF(k, G, C);
	t2 = GetTickCount64();
	cout << "执行时间：" << t2 - t1 << endl;  

	cout << "\n使用w-NAF方法（仿射坐标）:" << endl;
	t1 = GetTickCount64();
	EccPointMulW_NAF(k, G, 8, C);
	t2 = GetTickCount64();
	cout << "执行时间：" << t2 - t1 << endl;  

	cout << "\n使用二进制表示（标准射影坐标）：" << endl;
	t1 = GetTickCount64();
	EccPointMulBINStandardProjection(k, G_SP, C);
	t2 = GetTickCount64();
	cout << "执行时间：" << t2 - t1 << endl;

	cout << "\n使用NAF方法（标准射影坐标）：" << endl;
	t1 = GetTickCount64();
	EccPointMulNAFStandardProjection(k, G_SP, C);
	t2 = GetTickCount64();
	cout << "执行时间：" << t2 - t1 << endl;  

	cout << "\n使用w-NAF方法（标准射影坐标）：" << endl;
	t1 = GetTickCount64();
	EccPointMul_W_NAF_StandardProjection(k, G_SP, 8, C);
	t2 = GetTickCount64();
	cout << "执行时间：" << t2 - t1 << endl;  

	cout << "\n使用二进制表示（Jacobian加重射影坐标）：" << endl;
	t1 = GetTickCount64();
	EccPointMulBINJacobian(k, G_Jacobian, C);
	t2 = GetTickCount64();
	cout << "执行时间：" << t2 - t1 << endl;

	cout << "\n使用NAF方法（Jacobian加重射影坐标）：" << endl;
	t1 = GetTickCount64();
	EccPointMul_NAF_Jacobian(k, G_Jacobian, C);
	t2 = GetTickCount64();
	cout << "执行时间：" << t2 - t1 << endl;  

	cout << "\n使用w-NAF方法（Jacobian坐标）:" << endl;
	t1 = GetTickCount64();
	EccPointMul_W_NAF_Jacobian(k, G_Jacobian, 8, C);
	t2 = GetTickCount64();
	cout << "执行时间：" << t2 - t1 << endl;  
	*/
	cout << "以上为测试内容。" << endl;
	return 0;
}