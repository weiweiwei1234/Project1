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
#include <windows.h>
using namespace std;
//�ݲ���������Զ�㣬ʵ�ʼ������ʹ�ò���
//������Բ�����ϵĵ�	�������� Affine 
typedef struct ECPoint
{
	BigNumber x;	//x
	BigNumber y;	//y
}Point;

typedef struct ECPointStandardProjection	//��׼��Ӱ����
{
	BigNumber x;
	BigNumber y;
	BigNumber z;
}PointSP;

typedef struct ECPointJacobian	//Jacobian������Ӱ����
{
	BigNumber x;
	BigNumber y;
	BigNumber z;
}PointJ;

//������Բ���߲����ṹ��
typedef struct ECParams {
	BigNumber p;	//�������ģ��
	BigNumber a;	//��Բ���߲���a
	BigNumber b;	//��Բ���߲���b
	BigNumber n;	//Ⱥ�Ľ�
	BigNumber Gx;	//��������
	BigNumber Gy;
}EllipticCurveParams;

void printecparams(ECParams C);	//��ӡ��Բ���ߵĲ���
void printecpoint(ECPoint point);		//��ӡ������
void printecpointStandarProjection(ECPointStandardProjection point); //��ӡ��׼��Ӱ����
void printecpointJacobian(ECPointJacobian point);//��ӡJacobian������Ӱ����
bool isinparams(ECPoint point, ECParams C);	//�жϵ��Ƿ�����Բ������
bool isinparamsStandardProjection(ECPointStandardProjection point, ECParams C);//�жϵ��Ƿ��ڱ�׼��Ӱ�����ϵ���Բ������
bool isinparamsJacobian(ECPointJacobian point, ECParams C);//�ж��Ƿ���Jacobian����Բ������
ECPointStandardProjection AffineToStandardProjection(ECPoint P); //��������ת��Ϊ��׼��Ӱ����
ECPoint StandardProjectionToAffine(ECPointStandardProjection P,ECParams C); //��׼��Ӱ����ת��Ϊ�������� AffineToStandardProjection
ECPointJacobian AffineTOJacobian(ECPoint P);//��������ת��Ϊ�ſɱ�����
ECPoint JacobianToAffine(ECPointJacobian P, ECParams C);//�ſɱ�����ת��Ϊ��������
ECPoint ecpointadd(ECPoint P, ECPoint Q, ECParams C);	//����� ��������
ECPointStandardProjection ecpointaddStandardProjection(ECPointStandardProjection P, ECPointStandardProjection Q, ECParams C); //����� ��׼��Ӱ����
ECPointJacobian ecpointaddJacobian(ECPointJacobian P, ECPointJacobian Q, ECParams C);//����� �ſɱ�����
//����㷨ʵ��
ECPoint ecpointmul1(BigNumber k, ECPoint P, ECParams C);	//��ѭ��
//�Ӽ�������
ECPoint ecpointmulBIN(BigNumber k, ECPoint P, ECParams C);	//��k�����Ʊ�ʾ
ECPoint ecpointmulNAF(BigNumber k, ECPoint P, ECParams C);	//��k��NAF��ʾ
ECPoint ecpointmulW_NAF(BigNumber k, ECPoint P,int w,ECParams C);
/*
w-NAF�㷨 Ԥ������
���룺k, P, ���ڿ��w, ��Բ���߲���C
�����������kP
*/
ECPointStandardProjection ecpointmulNAFStandardProjection(BigNumber k, ECPointStandardProjection P, ECParams C);	//��׼��Ӱ�����µ�NAF���
ECPointJacobian ecpointmulNAKJacobian(BigNumber k, ECPointJacobian P, ECParams C);//�ſɱ������µ�NAF���

ECPoint ecpointmul4(BigNumber k, ECPoint P, ECParams C);

//��չgcd����Ԫ a * a^-1 = 1 (mod b)
BigNumber exgcd(BigNumber a, BigNumber b, BigNumber& x, BigNumber& y);
BigNumber modinverse(BigNumber a,	BigNumber b);  
BigNumber random(BigNumber n); //����һ������� ��[1,n-1]

