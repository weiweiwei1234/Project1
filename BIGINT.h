#pragma once
#include<iostream>
#include<string>
#include<cmath>

#define BIGGER 1
#define SMALLER -1
#define EQUAL 0

using namespace std;

//�������������ͣ���ģ�� ��data,p)
class BIGINT {
private:
	bool sgn;					//����	ʵ���о�ȡ��ֵ 
	uint8_t data[32] = { 0 };	//32λ�ֽ� ���򴢴��������� uint8_t 0-256��00-FF)
	uint8_t p[32] = { 0 };		//������p
	int k;						//�ֽڴ���С
	int m;						//���ش���С
	static int compare(const BIGINT&, const BIGINT&);//�ж��������Ƿ����
public:
	BIGINT();				 //Ĭ�Ϲ��캯��
	BIGINT(long,long);		 //ʹ�ý�С��������������� ������ʹ��
	BIGINT(bool, BIGINT);	 //һ������ȡ��ֵ
	BIGINT(bool, string);	 //ʹ���ַ������������ ȷ��������ֵ
	BIGINT(string);			 //��ȷ������ֵ
	int get_byte_length();	 //��ȡ�ֽڴ���С
	int get_bite_length();	 //��ȡ���ش���С
	BIGINT inverse(BIGINT m);//���� a mod m �ĳ˷���Ԫ 
	int to_int();			 //ת10������������ֵ��С��ʱ��
	string to_bite();		 //���ر�ʾ �����ַ���

	//�����߼������
	friend bool operator==(const BIGINT&, const BIGINT&);
	friend bool operator!=(const BIGINT&, const BIGINT&);
	friend bool operator>(const BIGINT&, const BIGINT&);
	friend bool operator<(const BIGINT&, const BIGINT&);
	friend bool operator>=(const BIGINT&, const BIGINT&);
	friend bool operator<=(const BIGINT&, const BIGINT&);

	// �������������
	friend const BIGINT operator+(const BIGINT&, const BIGINT&);
	friend const BIGINT operator-(const BIGINT&, const BIGINT&);
	friend const BIGINT operator*(const BIGINT&, const BIGINT&);
	friend const BIGINT operator/(const BIGINT&, const BIGINT&); //����
	friend const BIGINT operator%(const BIGINT&, const BIGINT&); //���� mod
	friend const BIGINT operator>>(const BIGINT&, const int&);   //����
	friend const BIGINT operator<<(const BIGINT&, const int&);   //����

	// �������
	friend std::ostream& operator<<(std::ostream&, const BIGINT&);
};