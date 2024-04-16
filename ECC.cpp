#include "SM3.h"
#include "EllipticCurve.h"
#include "BigNumber.h"
#include <windows.h>
//��Բ���߷��̣�y2 = x3 + ax + b��
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
	cout << "����Ϊ��������:" << endl;
	printecparams(C);
	printecpoint(G);
	printecpointStandarProjection(G_StandardProjection);
	printecpointJacobian(G_Jacobian);
	//cout << isinparams(G, C) << endl;
	//cout << isinparamsStandardProjection(G_StandardProjection, C) << endl;
	//cout << isinparamsJacobian(G_Jacobian, C)<< endl;
	BigNumber k("28E9FA9E9D9F5E344D5A9E4BCF6509A7F39789F515AB8F92DDBCBD414D940E93");
	//cout << HexToBin(k.get_value());
	/*cout << "ʹ�����ص��:" << endl;
	ECPoint Pk1 = ecpoint_mul_1(k, G, C);
	print_ecpoint(Pk1);*/

	long t1,t2;//��������ʱ�䣬t1:��ʼʱ��,t2:����ʱ��

	cout << "\n�Ӽ��������������ƣ���ʾ��" << endl;
	t1 = GetTickCount64();
	ECPoint Pk2 = ecpointmulBIN(k, G, C);
	t2 = GetTickCount64();
	cout << "ִ��ʱ�䣺" << t2 - t1 << endl;  //�������е�ʱ��õ���ʱ�䵥λΪ���� /1000Ϊ��
	printecpoint(Pk2);

	cout << "\nʹ��NAF�������������꣩��" << endl;
	t1 = GetTickCount64();
	ECPoint Pk3 = ecpointmulNAF(k, G, C);
	t2 = GetTickCount64();
	cout << "ִ��ʱ�䣺" << t2 - t1 << endl;  //�������е�ʱ��õ���ʱ�䵥λΪ���� /1000Ϊ��
	printecpoint(Pk3);

	cout << "\nʹ��w-NAF�������������꣩:" << endl;
	ECPoint Pk3_ = ecpointmulW_NAF(k, G, 4, C);
	printecpoint(Pk3_);

	cout << "\nʹ��NAF��������׼��Ӱ���꣩��" << endl;
	t1 = GetTickCount64();
	ECPointStandardProjection Pk4 = ecpointmulNAFStandardProjection(k, G_StandardProjection, C);
	t2 = GetTickCount64();
	cout << "ִ��ʱ�䣺" << t2 - t1 << endl;  //�������е�ʱ��õ���ʱ�䵥λΪ���� /1000Ϊ��
	ECPoint Pk4_ = StandardProjectionToAffine(Pk4, C);
	printecpoint(Pk4_);
	printecpointStandarProjection(Pk4);

	cout << "\nʹ��NAF������Jacobian������Ӱ���꣩��" << endl;
	t1 = GetTickCount64();
	ECPointJacobian Pk5 = ecpointmulNAKJacobian(k, G_Jacobian, C);
	t2 = GetTickCount64();
	cout << "ִ��ʱ�䣺" << t2 - t1 << endl;  //�������е�ʱ��õ���ʱ�䵥λΪ���� /1000Ϊ��
	ECPoint Pk5_ = JacobianToAffine(Pk5, C);
	printecpoint(Pk5_);
	printecpointJacobian(Pk5);

	cout << "����Ϊ�������ݡ�" << endl;
	return 0;
}