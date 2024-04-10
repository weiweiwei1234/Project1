//��Բ�����ϵ�����
/*
�Ƽ�ʹ��������256λ��Բ���ߡ�
��Բ���߷��̣�y2=x3+ax+b��
���߲���
p=FFFFFFFE FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF 00000000 FFFFFFFF FFFFFFFF
a=FFFFFFFE FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF 00000000 FFFFFFFF FFFFFFFC
b=28E9FA9E 9D9F5E34 4D5A9E4B CF6509A7 F39789F5 15AB8F92 DDBCBD41 4D940E93
n=FFFFFFFE FFFFFFFF FFFFFFFF FFFFFFFF 7203DF6B 21C6052B 53BBF409 39D54123
Gx=32C4AE2C 1F198119 5F990446 6A39C994 8FE30BBF F2660BE1 715A4589 334C74C7
Gy=BC3736A2 F4F6779C 59BDCEE3 6B692153 D0A9877C C62A4740 02DF32E5 2139F0A0
*/
#include "BigNumber.h"
using namespace std;
//������Բ�����ϵĵ�
typedef struct ECPoint
{
	BigNumber x;	//x
	BigNumber y;	//y
};

//������Բ���߲����ṹ��
typedef struct ECParams {
	BigNumber p;	//�������ģ��
	BigNumber a;	//��Բ���߲���a
	BigNumber b;	//��Բ���߲���b
	BigNumber n;	//Ⱥ�Ľ�
	BigNumber Gx;	//��������
	BigNumber Gy;
}EllipticCurveParams;
void print_ecparams(ECParams params);	//��ӡ��Բ���ߵĲ���
void print_ecpoint(ECPoint point);		//��ӡ������
bool is_in_params(ECPoint point, ECParams params);	//�жϵ��Ƿ�����Բ������
ECPoint ecpoint_add(ECPoint P, ECPoint Q, ECParams params);	//�����

//����㷨ʵ��
ECPoint ecpoint_mul_1(BigNumber k, ECPoint P, ECParams params);	
ECPoint ecpoint_mul_2(BigNumber k, ECPoint P, ECParams params);	

//��չgcd����Ԫ a * a^-1 = 1 (mod b)
BigNumber exgcd(BigNumber a, BigNumber b, BigNumber& x, BigNumber& y);
BigNumber mod_inverse(BigNumber a,	BigNumber b);  
BigNumber random(BigNumber n); //����һ������� ��[1,n-1]
