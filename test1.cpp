#include "test1.h"
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
void test1() {
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
	long t1, t2;//��������ʱ�䣬t1:��ʼʱ��,t2:����ʱ��

	cout << Mod_inverse(a, p) << endl;

	cout << "����Ϊ��Բ�����������������:" << endl;
	printEccParams(C);
	printEccPoint(G);
	printEccPointJacobian(G_Jacobian);
	BIGNUM k("ABCDEF");

	cout << "\n�����Ʊ�ʾ���������꣩��" << endl;
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

	cout << "\nʹ�ö����Ʊ�ʾ����׼��Ӱ���꣩��" << endl;
	t1 = GetTickCount64();
	EccPointStandardProjection SP1 = EccPointMulBINStandardProjection(k,G_SP,C);
	t2 = GetTickCount64();
	cout << "ִ��ʱ�䣺" << t2 - t1 << endl;

	cout << "\nʹ��NAF��������׼��Ӱ���꣩��" << endl;
	t1 = GetTickCount64();
	EccPointStandardProjection SP2 = EccPointMulNAFStandardProjection(k, G_SP, C);
	t2 = GetTickCount64();
	cout << "ִ��ʱ�䣺" << t2 - t1 << endl;  //�������е�ʱ��õ���ʱ�䵥λΪ���� /1000Ϊ��

	cout << "\nʹ��w-NAF������" << endl;
	t1 = GetTickCount64();
	EccPointStandardProjection SP2 = EccPointMul_W_NAF_StandardProjection(k, G_SP, 6, C);
	t2 = GetTickCount64();
	cout << "ִ��ʱ�䣺" << t2 - t1 << endl;  //�������е�ʱ��õ���ʱ�䵥λΪ���� /1000Ϊ��

	cout << "\nʹ�ö����Ʊ�ʾ��Jacobian������Ӱ���꣩��" << endl;
	t1 = GetTickCount64();
	EccPointJacobian SP1 = EccPointMulBINJacobian(k, G_Jacobian, C);
	t2 = GetTickCount64();
	cout << "ִ��ʱ�䣺" << t2 - t1 << endl;

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

void test2()
{
/*
1���û�Aѡ��һ���ʺϼ��ܵ���Բ����Ep(a,b)(��:y2=x3+ax+b)����ȡ��Բ������һ�㣬��Ϊ����G��
2���û�Aѡ��һ��˽����Կk�������ɹ�����ԿK=kG��
3���û�A��Ep(a,b)�͵�K��G�����û�B��
4���û�B�ӵ���Ϣ�� ��������������ı��뵽Ep(a,b)��һ��M��������һ���������r��r<n����
5���û�B�����C1=M+rK��C2=rG��
6���û�B��C1��C2�����û�A��
7���û�A�ӵ���Ϣ�󣬼���C1-kC2��������ǵ�M����Ϊ
  C1-kC2=M+rK-k(rG)=M+rK-r(kG)=M
�ٶԵ�M���н���Ϳ��Եõ����ġ�
*/
}
