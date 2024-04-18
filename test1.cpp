#include "SM3.h"
#include "EllipticCurve.h"
#include "BIGNUM.h"
#include "test1.h"
#include <windows.h>
#include <random>
#include <sstream>
#include <iomanip>
#include <limits>
using namespace std;
static void test1() {
	BIGNUM p("FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFF");
	BIGNUM a("FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFC");
	BIGNUM b("28E9FA9E9D9F5E344D5A9E4BCF6509A7F39789F515AB8F92DDBCBD414D940E93");
	BIGNUM n("FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFF7203DF6B21C6052B53BBF40939D54123");
	BIGNUM Gx("32C4AE2C1F1981195F9904466A39C9948FE30BBFF2660BE1715A4589334C74C7");
	BIGNUM Gy("BC3736A2F4F6779C59BDCEE36B692153D0A9877CC62A474002DF32E52139F0A0");
	EccPoint G = { Gx,Gy };
	EccParams C = { p,a,b,n,Gx,Gy };
	EccPointJacobian G_Jacobian = AffineTOJacobian(G);
	long t1, t2;//��������ʱ�䣬t1:��ʼʱ��,t2:����ʱ��

	cout << Mod_inverse(a, p) << endl;

	cout << "����Ϊ��Բ�����������������:" << endl;
	printEccParams(C);
	printEccPoint(G);
	printEccPointJacobian(G_Jacobian);
	BIGNUM k("ABEEEEEEEEEEEEEEEEEEECDEF");

	cout << "\n�����Ʊ�ʾ��" << endl;
	t1 = GetTickCount64();
	EccPoint Pk2 = EccPointMulBIN(k, G, C);
	t2 = GetTickCount64();
	cout << "ִ��ʱ�䣺" << t2 - t1 << endl;  //�������е�ʱ��õ���ʱ�䵥λΪ���� /1000Ϊ��
	printEccPoint(Pk2);

	cout << "\nʹ��NAF�������������꣩��" << endl;
	t1 = GetTickCount64();
	EccPoint Pk3 = EccPointMulNAF(k, G, C);
	t2 = GetTickCount64();
	cout << "ִ��ʱ�䣺" << t2 - t1 << endl;  //�������е�ʱ��õ���ʱ�䵥λΪ���� /1000Ϊ��
	printEccPoint(Pk3);

	cout << "\nʹ��w-NAF�������������꣩:" << endl;
	EccPoint Pk3_ = EccPointMulW_NAF(k, G, 6, C);
	printEccPoint(Pk3_);

	cout << "\nʹ��NAF������Jacobian������Ӱ���꣩��" << endl;
	t1 = GetTickCount64();
	EccPointJacobian Pk5 = EccPointMul_NAF_Jacobian(k, G_Jacobian, C);
	t2 = GetTickCount64();
	cout << "ִ��ʱ�䣺" << t2 - t1 << endl;  //�������е�ʱ��õ���ʱ�䵥λΪ���� /1000Ϊ��
	printEccPointJacobian(Pk5);

	cout << "\nʹ��w-NAF������Jacobian���꣩:" << endl;
	EccPointJacobian Pk6 = EccPointMul_W_NAF_Jacobian(k, G_Jacobian, 6, C);
	printEccPoint(JacobianToAffine(Pk6, C));
	printEccPointJacobian(Pk6);

	cout << "����Ϊ�������ݡ�" << endl;
}