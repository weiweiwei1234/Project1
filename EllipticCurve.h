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
#include "BIGNUM.h"
#include <windows.h>
using namespace std;
//�ݲ���������Զ�㣬ʵ�ʼ������ʹ�ò���
//������Բ�����ϵĵ�	�������� Affine
class EccPoint{
public:
	//bool inf;	//����Զ��
	BIGNUM x, y;
};

class EccPointStandardProjection	//��׼��Ӱ����
{
public:
	BIGNUM x;
	BIGNUM y;
	BIGNUM z;
};

class EccPointJacobian	//Jacobian������Ӱ����
{
public:
	BIGNUM x;
	BIGNUM y;
	BIGNUM z;
};

//������Բ���߲����ṹ��
class  EccParams {
public:
	BIGNUM p;	//�������ģ��
	BIGNUM a;	//��Բ���߲���a
	BIGNUM b;	//��Բ���߲���b
	BIGNUM n;	//Ⱥ�Ľ�
	BIGNUM Gx;	//��������
	BIGNUM Gy;
};

void printEccParams(EccParams);	//��ӡ��Բ���ߵĲ���
void printEccPoint(EccPoint);		//��ӡ������
void printEccPointStandarProjection(EccPointStandardProjection); //��ӡ��׼��Ӱ����
void printEccPointJacobian(EccPointJacobian);//��ӡJacobian������Ӱ����

bool isinEccParams(EccPoint, EccParams);	//�жϵ��Ƿ�����Բ������
bool isinEccParamsStandardProjection(EccPointStandardProjection, EccParams);//�жϵ��Ƿ��ڱ�׼��Ӱ�����ϵ���Բ������
bool isinEccParamsJacobian(EccPointJacobian, EccParams);//�ж��Ƿ���Jacobian

EccPointStandardProjection AffineToStandardProjection(EccPoint); //��������ת��Ϊ��׼��Ӱ����
EccPoint StandardProjectionToAffine(EccPointStandardProjection,EccParams); //��׼��Ӱ����ת��Ϊ�������� AffineToStandardProjection
EccPointJacobian AffineTOJacobian(EccPoint);//��������ת��Ϊ�ſɱ�����
EccPoint JacobianToAffine(EccPointJacobian, EccParams);//�ſɱ�����ת��Ϊ��������

EccPoint EccPointAdd(EccPoint, EccPoint, EccParams);	//����� ��������
EccPointStandardProjection EccPointAddStandardProjection(EccPointStandardProjection, EccPointStandardProjection, EccParams); //����� ��׼��Ӱ����
EccPointJacobian EccPointAddJacobian(EccPointJacobian, EccPointJacobian, EccParams);//����� 

//����㷨ʵ��
EccPoint EccPointMul1(BIGNUM, EccPoint, EccParams);	//��ѭ��
//�Ӽ�������
EccPoint EccPointMulBIN(BIGNUM, EccPoint, EccParams);	//��k�����Ʊ�ʾ
EccPoint EccPointMulNAF(BIGNUM, EccPoint, EccParams);	//��k��NAF��ʾ
EccPoint EccPointMulW_NAF(BIGNUM, EccPoint,int,EccParams); //w-NAF�㷨
/*
w-NAF�㷨 Ԥ�����
���룺k, P, ���ڿ��w, ��Բ���߲���C
�����������kP
*/
EccPointStandardProjection EccPointMulNAFStandardProjection(BIGNUM, EccPointStandardProjection, EccParams);	//��׼��Ӱ�����µ�NAF���

EccPointJacobian EccPointMul_NAF_Jacobian(BIGNUM, EccPointJacobian, EccParams);//�ſɱ������µ�NAF���
EccPointJacobian EccPointMul_W_NAF_Jacobian(BIGNUM, EccPointJacobian, int,EccParams);//�ſɱ������µ�w-NAF���


EccPoint EccPointMul4(BIGNUM, EccPoint, EccParams);

//��չgcd����Ԫ
BIGNUM exgcd(BIGNUM, BIGNUM, BIGNUM&, BIGNUM&);
BIGNUM Mod_inverse(BIGNUM,BIGNUM);  

//�ɸ�����ģ��
BIGNUM Montgomery_Multiply(BIGNUM, BIGNUM, BIGNUM);
//�ɸ�����Լ��
BIGNUM Montgomery_Reduction(BIGNUM, BIGNUM);

BIGNUM random(BIGNUM); //����һ������� 1 < m < n - 1

