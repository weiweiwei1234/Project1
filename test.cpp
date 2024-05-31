#include "SM3.h"
#include "ECC.h"
#include "BIGNUM.h"
#include <windows.h>
#include <random>
#include <sstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <cctype>
//��Բ���߷��̣�y2 = x3 + ax + b��
int main() {
	BIGNUM p("FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFF");
	BIGNUM a("FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFC");
	BIGNUM b("28E9FA9E9D9F5E344D5A9E4BCF6509A7F39789F515AB8F92DDBCBD414D940E93");
	BIGNUM n("FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFF7203DF6B21C6052B53BBF40939D54123");
	BIGNUM Gx("32C4AE2C1F1981195F9904466A39C9948FE30BBFF2660BE1715A4589334C74C7");
	BIGNUM Gy("BC3736A2F4F6779C59BDCEE36B692153D0A9877CC62A474002DF32E52139F0A0");
	Point G = { Gx,Gy };
	Params C = { p,a,b,n,Gx,Gy };
	SPoint G_S = Affine_To_StandardProjection(G);
	JPoint G_J = Affine_To_Jacobian(G);
	long t1, t2;//��������ʱ�䣬t1:��ʼʱ��,t2:����ʱ��
	
	BIGNUM k = Random(64, n);
	cout << k << endl;
	
	cout << "\n�����Ʊ�ʾ���������꣩��" << endl;
	t1 = GetTickCount64();
	Point result1=Point_Mul_Bin(k, G, C);
	print_Point(result1);
	cout << is_in_Params(result1, C) << endl;
	t2 = GetTickCount64();
	cout << "ִ��ʱ�䣺" << t2 - t1 << endl;  

	cout << "\n�����Ʊ�ʾ���������꣩(�ɸ������Ż�����" << endl;
	t1 = GetTickCount64();
	Point result2 = Point_Mul_Bin_M(k, G, C);
	print_Point(result2);
	cout << is_in_Params(result2, C) << endl;
	t2 = GetTickCount64();
	cout << "ִ��ʱ�䣺" << t2 - t1 << endl;

	cout << "\n�����Ʊ�ʾ����׼��Ӱ���꣩��" << endl;
	t1 = GetTickCount64();
	SPoint result3 = SPoint_Mul_Bin(k, G_S, C);
	print_SPoint(result3);
	cout << is_in_Params_S(result3, C) << endl;
	t2 = GetTickCount64();
	cout << "ִ��ʱ�䣺" << t2 - t1 << endl;

	cout << "\n�����Ʊ�ʾ����׼��Ӱ���꣩(�ɸ������Ż�����" << endl;
	t1 = GetTickCount64();
	SPoint result4 = SPoint_Mul_Bin_M(k, G_S, C);
	print_SPoint(result4);
	cout << is_in_Params_S(result4, C) << endl;
	t2 = GetTickCount64();
	cout << "ִ��ʱ�䣺" << t2 - t1 << endl;

	cout << "\n�����Ʊ�ʾ��Jacobian������Ӱ���꣩��" << endl;
	t1 = GetTickCount64();
	JPoint result5 = JPoint_Mul_Bin(k, G_J, C);
	print_JPoint(result5);
	cout << is_in_Params_J(result5, C) << endl;
	t2 = GetTickCount64();
	cout << "ִ��ʱ�䣺" << t2 - t1 << endl;

	cout << "\nNAF��ʾ���������꣩��" << endl;
	t1 = GetTickCount64();
	Point result6 = Point_Mul_NAF(k, G, C);
	print_Point(result6);
	cout << is_in_Params(result6, C) << endl;
	t2 = GetTickCount64();
	cout << "ִ��ʱ�䣺" << t2 - t1 << endl;

	cout << "\nNAF��ʾ���������꣩(�ɸ������Ż�����" << endl;
	t1 = GetTickCount64();
	Point result7 = Point_Mul_NAF_M(k, G, C);
	print_Point(result7);
	cout << is_in_Params(result7, C) << endl;
	t2 = GetTickCount64();
	cout << "ִ��ʱ�䣺" << t2 - t1 << endl;

	cout << "\nNAF��ʾ����׼��Ӱ���꣩��" << endl;
	t1 = GetTickCount64();
	SPoint result8 = SPoint_Mul_NAF(k, G_S, C);
	print_SPoint(result8);
	cout << is_in_Params_S(result8, C) << endl;
	t2 = GetTickCount64();
	cout << "ִ��ʱ�䣺" << t2 - t1 << endl;

	cout << "\nNAF��ʾ����׼��Ӱ���꣩(�ɸ������Ż�����" << endl;
	t1 = GetTickCount64();
	SPoint result9 = SPoint_Mul_NAF_M(k, G_S, C);
	print_SPoint(result9);
	cout << is_in_Params_S(result9, C) << endl;
	t2 = GetTickCount64();
	cout << "ִ��ʱ�䣺" << t2 - t1 << endl;

	cout << "\nNAF��ʾ��Jacobian������Ӱ���꣩��" << endl;
	t1 = GetTickCount64();
	JPoint result10 = JPoint_Mul_NAF(k, G_J, C);
	print_JPoint(result10);
	cout << is_in_Params_J(result10, C) << endl;
	t2 = GetTickCount64();
	cout << "ִ��ʱ�䣺" << t2 - t1 << endl;

	cout << "\nwNAF��ʾ���������꣩��" << endl;
	t1 = GetTickCount64();
	Point result11 = Point_Mul_wNAF(k, G, C,8);
	print_Point(result11);
	cout << is_in_Params(result11, C) << endl;
	t2 = GetTickCount64();
	cout << "ִ��ʱ�䣺" << t2 - t1 << endl;

	cout << "\nwNAF��ʾ���������꣩(�ɸ������Ż�����" << endl;
	t1 = GetTickCount64();
	Point result12 = Point_Mul_wNAF_M(k, G, C,8);
	print_Point(result12);
	cout << is_in_Params(result12, C) << endl;
	t2 = GetTickCount64();
	cout << "ִ��ʱ�䣺" << t2 - t1 << endl;

	cout << "\nwNAF��ʾ����׼��Ӱ���꣩��" << endl;
	t1 = GetTickCount64();
	SPoint result13 = SPoint_Mul_wNAF(k, G_S, C,8);
	print_SPoint(result13);
	cout << is_in_Params_S(result13, C) << endl;
	t2 = GetTickCount64();
	cout << "ִ��ʱ�䣺" << t2 - t1 << endl;

	cout << "\nwNAF��ʾ����׼��Ӱ���꣩(�ɸ������Ż�����" << endl;
	t1 = GetTickCount64();
	SPoint result14 = SPoint_Mul_wNAF_M(k, G_S, C,8);
	print_SPoint(result14);
	cout << is_in_Params_S(result14, C) << endl;
	t2 = GetTickCount64();
	cout << "ִ��ʱ�䣺" << t2 - t1 << endl;

	cout << "\nwNAF��ʾ��Jacobian������Ӱ���꣩��" << endl;
	t1 = GetTickCount64();
	JPoint result15 = JPoint_Mul_wNAF(k, G_J, C,8);
	cout << is_in_Params_J(result15, C) << endl;
	print_JPoint(result15);
	t2 = GetTickCount64();
	cout << "ִ��ʱ�䣺" << t2 - t1 << endl;
	/*
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
	cout << "a���ɸ�������ʾ��" << a_R << endl;
	cout << "b���ɸ�������ʾ��" << b_R << endl;
	cout << "p���ɸ�������ʾ��" << N_R << endl;
	BIGNUM _a = Mod_inverse(a, p);
	BIGNUM _a_R = Mod_inverse(a_R, p);
	cout << _a << endl;
	cout << (_a_R << 256) % p << endl;

	
	BIGNUM k = random(64, n);
	cout << k << endl;

	cout << "\n�����Ʊ�ʾ���������꣩��" << endl;
	t1 = GetTickCount64();
	EccPointMulBIN(k, G, C);
	t2 = GetTickCount64();
	cout << "ִ��ʱ�䣺" << t2 - t1 << endl;  //�������е�ʱ��õ���ʱ�䵥λΪ���� /1000Ϊ��

	cout << "\n�����Ʊ�ʾ���������꣩��" << endl;
	t1 = GetTickCount64();
	EccPointMul_Montomery(k, G, C);
	t2 = GetTickCount64();
	cout << "ִ��ʱ�䣺" << t2 - t1 << endl;  //�������е�ʱ��õ���ʱ�䵥λΪ���� /1000Ϊ��
	//cout << a * b << endl;
	*/
	/*
	cout << Gx * Gy << endl;
	cout << Gx * Gy % p << endl;
	cout << Gx + Gy % p << endl;
	cout << (Gx - Gy + p)%p << endl;
	*/









	/*
	//����������������ʱ��
	long time_random = 0;	//���������
	long time_add = 0;		//�ӷ�
	long time_sub = 0;		//����
	long time_mul = 0;		//�˷�
	long time_mod = 0;		//����
	long time_weiyi = 0;	//λ��
	long time_inverse = 0;	//����Ԫ
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
	cout << "1000�����������ʱ��:" << time_random << endl;
	cout << "1000�μӷ�����ʱ��:" << time_add << endl;
	cout << "1000�μ�������ʱ��:" << time_sub << endl;
	cout << "1000�γ˷�����ʱ��:" << time_mul << endl;
	cout << "1000�����෨����ʱ��:" << time_mod << endl;
	cout << "1000��λ�Ʒ�����ʱ��:" << time_weiyi << endl;
	cout << "1000������Ԫ����ʱ��:" << time_inverse << endl;
	
	






















	cout << "����Ϊ��Բ���ߵ�������������:" << endl;

	BIGNUM k = random(64, n);
	cout << k << endl;
	//cout << "montgomery_ladder�㷨��" << endl;
	//t1 = GetTickCount64();
	//Mul_Montgomery_ladder(k, G, C);
	//t2 = GetTickCount64();
	//cout << "ִ��ʱ�䣺" << t2 - t1 << endl;

	cout << "\n�����Ʊ�ʾ���������꣩��" << endl;
	t1 = GetTickCount64();
	EccPointMulBIN(k, G, C);
	t2 = GetTickCount64();
	cout << "ִ��ʱ�䣺" << t2 - t1 << endl;  //�������е�ʱ��õ���ʱ�䵥λΪ���� /1000Ϊ��

	cout << "\n�ɸ������Ż���Ķ����Ʊ�ʾ���������꣩��" << endl;
	t1 = GetTickCount64();

	t2 = GetTickCount64();
	cout << "ִ��ʱ�䣺" << t2 - t1 << endl;  //�������е�ʱ��õ���ʱ�䵥λΪ���� /1000Ϊ��

	cout << "\nʹ��NAF�������������꣩��" << endl;
	t1 = GetTickCount64();
	EccPointMulNAF(k, G, C);
	t2 = GetTickCount64();
	cout << "ִ��ʱ�䣺" << t2 - t1 << endl;  

	cout << "\nʹ��w-NAF�������������꣩:" << endl;
	t1 = GetTickCount64();
	EccPointMulW_NAF(k, G, 8, C);
	t2 = GetTickCount64();
	cout << "ִ��ʱ�䣺" << t2 - t1 << endl;  

	cout << "\nʹ�ö����Ʊ�ʾ����׼��Ӱ���꣩��" << endl;
	t1 = GetTickCount64();
	EccPointMulBINStandardProjection(k, G_SP, C);
	t2 = GetTickCount64();
	cout << "ִ��ʱ�䣺" << t2 - t1 << endl;

	cout << "\nʹ��NAF��������׼��Ӱ���꣩��" << endl;
	t1 = GetTickCount64();
	EccPointMulNAFStandardProjection(k, G_SP, C);
	t2 = GetTickCount64();
	cout << "ִ��ʱ�䣺" << t2 - t1 << endl;  

	cout << "\nʹ��w-NAF��������׼��Ӱ���꣩��" << endl;
	t1 = GetTickCount64();
	EccPointMul_W_NAF_StandardProjection(k, G_SP, 8, C);
	t2 = GetTickCount64();
	cout << "ִ��ʱ�䣺" << t2 - t1 << endl;  

	cout << "\nʹ�ö����Ʊ�ʾ��Jacobian������Ӱ���꣩��" << endl;
	t1 = GetTickCount64();
	EccPointMulBINJacobian(k, G_Jacobian, C);
	t2 = GetTickCount64();
	cout << "ִ��ʱ�䣺" << t2 - t1 << endl;

	cout << "\nʹ��NAF������Jacobian������Ӱ���꣩��" << endl;
	t1 = GetTickCount64();
	EccPointMul_NAF_Jacobian(k, G_Jacobian, C);
	t2 = GetTickCount64();
	cout << "ִ��ʱ�䣺" << t2 - t1 << endl;  

	cout << "\nʹ��w-NAF������Jacobian���꣩:" << endl;
	t1 = GetTickCount64();
	EccPointMul_W_NAF_Jacobian(k, G_Jacobian, 8, C);
	t2 = GetTickCount64();
	cout << "ִ��ʱ�䣺" << t2 - t1 << endl;  
	*/
	cout << "����Ϊ�������ݡ�" << endl;
	return 0;
}