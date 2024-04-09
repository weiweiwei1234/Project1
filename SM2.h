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
// �������������
#pragma once
#include <string>
#include <iostream>
using namespace std;
typedef char* BigInteger;

//������Բ�����ϵĵ�
typedef struct
{
    BigInteger x; //x����
    BigInteger y; //y����
    int is_infinity;  // ���ڱ�ʾ����Զ��
}ECPoint;

// ������Բ���߲����ṹ��
typedef struct ECParams {
    BigInteger p;  // �������ģ��
    BigInteger a;  // ��Բ���߲���a
    BigInteger b;  // ��Բ���߲���b
    BigInteger n;  // Ⱥ�Ľ�
    ECPoint G;     // Ⱥ����
} EllipticCurveParams;

void init_big_integer(BigInteger* num);
void init_elliptic_curve_ecpoint(ECPoint* point);
void init_elliptic_curve_params(EllipticCurveParams* params);

void assign_big_integer(BigInteger* dest, BigInteger src);
void assign_elliptic_curve_ecpoint(ECPoint* dest, ECPoint src);
void assign_elliptic_curve_params(EllipticCurveParams* dest, EllipticCurveParams src);

void print_big_integer(BigInteger num);
void print_elliptic_curve_ecpoint(ECPoint point);
void print_elliptic_curve_params(EllipticCurveParams params);

BigInteger mod(BigInteger num, BigInteger mod_num);
BigInteger mod_inverse(BigInteger a, BigInteger m);
ECPoint ec_point_add(ECPoint P, ECPoint Q, ECParams params);
ECPoint ec_point_double(ECPoint P, ECParams params);

