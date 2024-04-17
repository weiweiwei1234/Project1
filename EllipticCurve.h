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
class EccPoint
{
public:
	BIGNUM x;	//x
	BIGNUM y;	//y
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
typedef struct EccParams {
	BIGNUM p;	//�������ģ��
	BIGNUM a;	//��Բ���߲���a
	BIGNUM b;	//��Բ���߲���b
	BIGNUM n;	//Ⱥ�Ľ�
	BIGNUM Gx;	//��������
	BIGNUM Gy;
}EllipticCurveParams;

void printEccParams(EccParams C);	//��ӡ��Բ���ߵĲ���
void printEccPoint(EccPoint point);		//��ӡ������
void printEccPointStandarProjection(EccPointStandardProjection point); //��ӡ��׼��Ӱ����
void printEccPointJacobian(EccPointJacobian point);//��ӡJacobian������Ӱ����
bool isinEccParams(EccPoint point, EccParams C);	//�жϵ��Ƿ�����Բ������
bool isinEccParamsStandardProjection(EccPointStandardProjection point, EccParams C);//�жϵ��Ƿ��ڱ�׼��Ӱ�����ϵ���Բ������
bool isinEccParamsJacobian(EccPointJacobian point, EccParams C);//�ж��Ƿ���Jacobian����Բ������
EccPointStandardProjection AffineToStandardProjection(EccPoint P); //��������ת��Ϊ��׼��Ӱ����
EccPoint StandardProjectionToAffine(EccPointStandardProjection P,EccParams C); //��׼��Ӱ����ת��Ϊ�������� AffineToStandardProjection
EccPointJacobian AffineTOJacobian(EccPoint P);//��������ת��Ϊ�ſɱ�����
EccPoint JacobianToAffine(EccPointJacobian P, EccParams C);//�ſɱ�����ת��Ϊ��������
EccPoint EccPointadd(EccPoint P, EccPoint Q, EccParams C);	//����� ��������
EccPointStandardProjection EccPointaddStandardProjection(EccPointStandardProjection P, EccPointStandardProjection Q, EccParams C); //����� ��׼��Ӱ����
EccPointJacobian EccPointaddJacobian(EccPointJacobian P, EccPointJacobian Q, EccParams C);//����� �ſɱ�����
//����㷨ʵ��
EccPoint EccPointmul1(BIGNUM k, EccPoint P, EccParams C);	//��ѭ��
//�Ӽ�������
EccPoint EccPointmulBIN(BIGNUM k, EccPoint P, EccParams C);	//��k�����Ʊ�ʾ
EccPoint EccPointmulNAF(BIGNUM k, EccPoint P, EccParams C);	//��k��NAF��ʾ
EccPoint EccPointmulW_NAF(BIGNUM k, EccPoint P,int w,EccParams C);
/*
w-NAF�㷨 Ԥ������
���룺k, P, ���ڿ��w, ��Բ���߲���C
�����������kP
*/
EccPointStandardProjection EccPointmulNAFStandardProjection(BIGNUM k, EccPointStandardProjection P, EccParams C);	//��׼��Ӱ�����µ�NAF���
EccPointJacobian EccPointmulNAKJacobian(BIGNUM k, EccPointJacobian P, EccParams C);//�ſɱ������µ�NAF���

EccPoint EccPointmul4(BIGNUM k, EccPoint P, EccParams C);

//��չgcd����Ԫ a * a^-1 = 1 (mod b)
BIGNUM exgcd(BIGNUM a, BIGNUM b, BIGNUM& x, BIGNUM& y);
BIGNUM modinverse(BIGNUM a,	BIGNUM b);  
BIGNUM random(BIGNUM n); //����һ������� ��[1,n-1]

