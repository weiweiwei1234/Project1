#pragma once
#include<iostream>
#include<string>

#define BIGGER 1
#define SMALLER -1
#define EQUAL 0

using namespace std;
class BIGINT {
private:
	bool sgn;			//����	ʵ���о�ȡ��ֵ
	int8_t data[64];	//����64λ16������ �൱��256λ������	���򴢴���������
	int8_t p[64];		//�������ϵ����� ģ��
	int size;			//�����Ĵ�С 
	static int compare(const BIGINT&, const BIGINT&);//�ж��������Ƿ����
public:

};