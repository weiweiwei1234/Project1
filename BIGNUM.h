//���ݵĻ������㣬�ӡ������ˡ�����mod���Ƚ�
#pragma once
#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdint>
#include <string>

#define BIGGER 1
#define SMALLER -1
#define EQUAL 0

using namespace std;
// Big number class definition
class BIGNUM {
private:
    // ��Ա����
    bool sgn;   //trueΪ+��falseΪ-
    std::vector<int8_t> data;     //int8_t:signed char ���� �����ݵľ���ֵ���Ŵ洢
    // �ж���� ����ǰ��0
    static int abs_compare(const BIGNUM&, const BIGNUM&);
    static void discard_leading_zero(std::vector<int8_t>&);

public:
    BIGNUM();
    BIGNUM(long long); // ֱ��ת��int
    BIGNUM(const std::string&);
    BIGNUM(bool, const std::vector<int8_t>&);
    BIGNUM(bool, std::vector<int8_t>&&);
    std::string get_value(); //����ֵ
    int to_int(); //С��תint
    int bitlen(); //�����Ƴ���
    std::string HexToBin();//ʮ������ת������

    // �����߼������
    friend bool operator==(const BIGNUM&, const BIGNUM&); 
    friend bool operator!=(const BIGNUM&, const BIGNUM&);
    friend bool operator>(const BIGNUM&, const BIGNUM&);
    friend bool operator<(const BIGNUM&, const BIGNUM&);
    friend bool operator>=(const BIGNUM&, const BIGNUM&);
    friend bool operator<=(const BIGNUM&, const BIGNUM&);

    // �������������
    friend const BIGNUM operator+(const BIGNUM&, const BIGNUM&);
    friend const BIGNUM operator-(const BIGNUM&, const BIGNUM&);
    friend const BIGNUM operator*(const BIGNUM&, const BIGNUM&);
    friend const BIGNUM operator/(const BIGNUM&, const BIGNUM&); //����
    friend const BIGNUM operator%(const BIGNUM&, const BIGNUM&); //���� mod
    friend const BIGNUM operator>>(const BIGNUM&, const int&);   //����
    friend const BIGNUM operator<<(const BIGNUM&, const int&);   //����

    // �������
    friend std::ostream& operator<<(std::ostream&, const BIGNUM&);


};//���ŷ��룬data����ֻ��ʾ��ֵ����ͨ��sgn����ʾ����