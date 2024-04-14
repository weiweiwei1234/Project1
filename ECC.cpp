#include "SM3.h"
#include "EllipticCurve.h"
#include "BigNumber.h"
#include <windows.h>
//��Բ���߷��̣�y2 = x3 + ax + b��
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
	cout << "����Ϊ��������:" << endl;
	printecparams(C);
	printecpoint(G);
	if (isinparams(G, C) == 1) {
		cout << "G����Բ������\n" << endl;
	}
	else {
		cout << "G������Բ������\n" << endl;
	}
	ECPoint P2 = ecpointadd(G, G, C);
	printecpoint(P2);
	if (isinparams(P2, C) == 1) {
		cout << "P2����Բ������\n" << endl;
	}
	else {
		cout << "P2������Բ������\n" << endl;
	}
	BigNumber k("FFFFFFE");

	/*cout << "ʹ�����ص��:" << endl;
	ECPoint Pk1 = ecpoint_mul_1(k, G, C);
	print_ecpoint(Pk1);*/

	long t1,t2;//��������ʱ�䣬t1:��ʼʱ��,t2:����ʱ��

	cout << "��k�����Ʊ�ʾ:" << endl;
	t1 = GetTickCount64();
	ECPoint Pk2 = ecpointmulBIN(k, G, C);
	t2 = GetTickCount64();
	cout << "ִ��ʱ�䣺" << t2 - t1 << endl;  //�������е�ʱ��õ���ʱ�䵥λΪ���� /1000Ϊ��
	printecpoint(Pk2);

	cout << "ʹ��NAF�˷�(�������꣩:" << endl;
	t1 = GetTickCount64();
	ECPoint Pk3 = ecpointmulNAF(k, G, C);
	t2 = GetTickCount64();
	cout << "ִ��ʱ�䣺" << t2 - t1 << endl;  //�������е�ʱ��õ���ʱ�䵥λΪ���� /1000Ϊ��
	printecpoint(Pk3);

	cout << "ʹ��NAF�˷�(��Ӱ���꣩:" << endl;
	ECPointStandardProjection G_ECPointStandardProjection = StandardProjectionToAffine(G);
	t1 = GetTickCount64();
	ECPointStandardProjection Pk4 = ecpointmulNAFStandardProjection(k, G_ECPointStandardProjection, C);
	t2 = GetTickCount64();
	cout << "ִ��ʱ�䣺" << t2 - t1 << endl;  //�������е�ʱ��õ���ʱ�䵥λΪ���� /1000Ϊ��
	ECPoint Pk4_ = AffineToStandardProjection(Pk4, C);
	printecpoint(Pk4_);
	printecpointStandarProjection(Pk4);

	cout << "����Ϊ�������ݡ�" << endl;
	return 0;
}