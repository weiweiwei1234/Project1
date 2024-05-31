#pragma once
//�����ʵ����Բ���߷������꣬��׼��Ӱ���꣬Jacobian������Ӱ�������ʵ��
//������չŷ�����ʵ����������
//ʵ����������ϵ�¶����Ʊ�ʾ��NAF��ʾ��w-NAF��ʾ�ĵ������
//ʹ���ɸ������㷨�Ե����������Ż�
/*
���߲���
p=FFFFFFFE FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF 00000000 FFFFFFFF FFFFFFFF
a=FFFFFFFE FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF 00000000 FFFFFFFF FFFFFFFC
b=28E9FA9E 9D9F5E34 4D5A9E4B CF6509A7 F39789F5 15AB8F92 DDBCBD41 4D940E93
n=FFFFFFFE FFFFFFFF FFFFFFFF FFFFFFFF 7203DF6B 21C6052B 53BBF409 39D54123
Gx=32C4AE2C 1F198119 5F990446 6A39C994 8FE30BBF F2660BE1 715A4589 334C74C7
Gy=BC3736A2 F4F6779C 59BDCEE3 6B692153 D0A9877C C62A4740 02DF32E5 2139F0A0
*/
#include"BIGNUM.h"
#include<windows.h>
using namespace std;

//�������� Affine
class Point {
public:
	BIGNUM x;
	BIGNUM y;
};

//��׼��Ӱ����
class SPoint {
public:
	BIGNUM x;
	BIGNUM y;
	BIGNUM z;
};

//Jacobian������Ӱ����
class JPoint {
public:
	BIGNUM x;
	BIGNUM y;
	BIGNUM z;
};

//������Բ���߲����ṹ��
class  Params {
public:
	BIGNUM p;	//�������ģ��
	BIGNUM a;	//��Բ���߲���a
	BIGNUM b;	//��Բ���߲���b
	BIGNUM n;	//Ⱥ�Ľ�
	BIGNUM Gx;	//��������
	BIGNUM Gy;
};

void print_Params(Params params);		//��ӡ��Բ���ߵĲ���
void print_Point(Point point);			//��ӡ��������
void print_SPoint(SPoint point);		//��ӡ��׼��Ӱ����
void print_JPoint(JPoint point);		//��ӡJacobian������Ӱ����

bool is_in_Params(Point point, Params params);		//����
bool is_in_Params_S(SPoint point, Params params);	//��׼
bool is_in_Params_J(JPoint point, Params params);	//Jacobian

SPoint Affine_To_StandardProjection(Point point); //��������ת��Ϊ��׼��Ӱ����
Point StandardProjection_To_Affine(SPoint point, Params params); //��׼��Ӱ����ת��Ϊ�������� AffineToStandardProjection
JPoint Affine_To_Jacobian(Point point);//��������ת��Ϊ�ſɱ�����
Point Jacobian_To_Affine(JPoint point, Params params);//�ſɱ�����ת��Ϊ��������

//��չgcd����Ԫ
BIGNUM EXGCD(BIGNUM a, BIGNUM b, BIGNUM& x, BIGNUM& y);
BIGNUM INVERSE(BIGNUM a, BIGNUM p);

//�ɸ�����ģ��
BIGNUM M_Multiply(BIGNUM a, BIGNUM b, BIGNUM N);
//�ɸ�����Լ��
BIGNUM M_Reduction(BIGNUM X, BIGNUM N);

BIGNUM Random(int d, BIGNUM N); //����һ������� 1 < m < n - 1

Point Point_Add(Point P, Point Q, Params C);	//����� ��������
Point Point_Add_M(Point PR, Point QR, Params C);	//�ɸ�����
SPoint SPoint_Add(SPoint P, SPoint Q, Params C); //����� ��׼��Ӱ����
SPoint SPoint_Add_M(SPoint PR, SPoint QR, Params C);//�ɸ�����
JPoint JPoint_Add(JPoint P, JPoint Q, Params C);//����� Jacobian������Ӱ���� 
JPoint JPoint_Add_M(JPoint PR, JPoint QR, Params C);//�ɸ�����

//��˷���ʵ��
//�����Ʊ�ʾ
Point Point_Mul_Bin(BIGNUM k, Point P, Params C);
Point Point_Mul_Bin_M(BIGNUM k, Point P, Params C);//�ɸ�����
SPoint SPoint_Mul_Bin(BIGNUM k, SPoint P, Params C);
SPoint SPoint_Mul_Bin_M(BIGNUM k, SPoint P, Params C);//�ɸ�����
JPoint JPoint_Mul_Bin(BIGNUM k, JPoint P, Params C);
JPoint JPoint_Mul_Bin_M(BIGNUM k, JPoint P, Params C);//�ɸ�����

//NAF��ʾ
Point Point_Mul_NAF(BIGNUM k, Point P, Params C);
Point Point_Mul_NAF_M(BIGNUM k, Point P, Params C);//�ɸ�����
SPoint SPoint_Mul_NAF(BIGNUM k, SPoint P, Params C);
SPoint SPoint_Mul_NAF_M(BIGNUM k, SPoint P, Params C);//�ɸ�����
JPoint JPoint_Mul_NAF(BIGNUM k, JPoint P, Params C);
JPoint JPoint_Mul_NAF_M(BIGNUM k, JPoint P, Params C);//�ɸ�����

//w-NAF��ʾ
Point Point_Mul_wNAF(BIGNUM k, Point P, Params C, int w);
Point Point_Mul_wNAF_M(BIGNUM k, Point P, Params C, int w);//�ɸ�����
SPoint SPoint_Mul_wNAF(BIGNUM k, SPoint P, Params C, int w);
SPoint SPoint_Mul_wNAF_M(BIGNUM k, SPoint P, Params C, int w);//�ɸ�����
JPoint JPoint_Mul_wNAF(BIGNUM k, JPoint P, Params C, int w);
JPoint JPoint_Mul_wNAF_M(BIGNUM k, JPoint P, Params C, int w);//�ɸ�����