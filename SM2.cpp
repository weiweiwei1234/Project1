#include "SM2.h"

// ��ʼ��������
void init_big_integer(BigInteger* num) {
    *num = NULL;
}
//��ʼ����Բ�����ϵĵ�
void init_elliptic_curve_ecpoint(ECPoint* point) {
    init_big_integer(&(point->x));
    init_big_integer(&(point->y));
    ((point)->is_infinity) = 0;
}
// ��ʼ����Բ���߲���
void init_elliptic_curve_params(EllipticCurveParams* params) {
    init_big_integer(&(params->p));
    init_big_integer(&(params->a));
    init_big_integer(&(params->b));
    init_big_integer(&(params->n));
    init_elliptic_curve_ecpoint(&(params)->G);
}
// ��ֵ������
void assign_big_integer(BigInteger* dest, BigInteger src) {
    *dest = src;
}
// ��ֵ��Բ�����ϵĵ�
void assign_elliptic_curve_ecpoint(ECPoint* dest, ECPoint src) {
    assign_big_integer(&(dest->x), src.x);
    assign_big_integer(&(dest->y), src.y);
    dest->is_infinity = src.is_infinity;
}
// ��ֵ��Բ���߲���
void assign_elliptic_curve_params(EllipticCurveParams* dest, EllipticCurveParams src) {
    assign_big_integer(&(dest->p), src.p);
    assign_big_integer(&(dest->a), src.a);
    assign_big_integer(&(dest->b), src.b);
    assign_big_integer(&(dest->n), src.n);
    assign_elliptic_curve_ecpoint(&(dest->G), src.G);
}
// ��ӡ������
void print_big_integer(BigInteger num) {
    cout << num << endl;
}
// ��ӡ��Բ�����ϵĵ�
void print_elliptic_curve_ecpoint(ECPoint point) {
    print_big_integer(point.x);
    print_big_integer(point.y);
}
// ��ӡ��Բ���߲���
void print_elliptic_curve_params(EllipticCurveParams params) {
    cout<<"p: "; print_big_integer(params.p); 
    cout<<"a: "; print_big_integer(params.a); 
    cout<<"b: "; print_big_integer(params.b); 
    cout<<"n: "; print_big_integer(params.n); 
    cout<<"G: "; print_elliptic_curve_ecpoint(params.G); 
}
// �������ȡģ
BigInteger mod(BigInteger num, BigInteger mod_num) {
    BigInteger result;
    // ������ʵ�ִ�����ȡģ���߼�
    result = NULL;
    return result;
}

BigInteger mod_inverse(BigInteger a, BigInteger m) {
    // ����ģm�µ�a����Ԫ
    BigInteger result;
    result = NULL;
    return result;
}
// ��Բ���ߵ���� P+Q
ECPoint ec_point_add(ECPoint P, ECPoint Q, ECParams params) {
    ECPoint result;
    // ������ʵ����Բ���ߵ���ӵ��߼�
    return result;
}

// ��Բ���ߵ㱶�� 2P
ECPoint ec_point_double(ECPoint P, ECParams params) {
    ECPoint result;
    // ������ʵ����Բ���ߵ㱶�˵��߼�
    return result;
}

