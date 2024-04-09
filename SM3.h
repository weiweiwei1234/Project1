#pragma once
#include <iostream>
#include <string>
#include <cmath>
using namespace std;
string BinToHex(string str);							//������ת��Ϊʮ�����ƺ���ʵ��
string HexToBin(string str);							//ʮ������ת��Ϊ�����ƺ���ʵ��
int BinToDec(string str);								//������ת��Ϊʮ���Ƶĺ���ʵ��
string DecToBin(int str);								//ʮ����ת��Ϊ�����Ƶĺ���ʵ��
int HexToDec(string str);								//ʮ������ת��Ϊʮ���Ƶĺ���ʵ��
string DecToHex(int str);								//ʮ����ת��Ϊʮ�����Ƶĺ���ʵ��
string padding(string str);								//�����ݽ������ 
string LeftShift(string str, int len);					//ʵ��ѭ������lenλ����
string XOR(string str1, string str2);					//ʵ��������
string AND(string str1, string str2);					//ʵ�������
string OR(string str1, string str2);					//ʵ�ֻ����
string NOT(string str);									//ʵ�ַǲ���
char binXor(char str1, char str2);						//ʵ�ֵ����ص�������
char binAnd(char str1, char str2);						//ʵ�ֵ����ص������
string ModAdd(string str1, string str2);				//mod 2^32����ĺ���ʵ��
string P1(string str);									//ʵ���û�����P1��X��
string P0(string str);									//ʵ���û�����P0��X��
string T(int j);										//����Tj����ֵ�ĺ���ʵ��
string FF(string str1, string str2, string str3, int j);//ʵ�ֲ�������FF����
string GG(string str1, string str2, string str3, int j);//ʵ�ֲ�������GG����
string extension(string str);							//��Ϣ��չ����
string compress(string str1, string str2);				//��Ϣѹ������
string iteration(string str);							//����ѹ������ʵ��
string SM3(string str);									//�����Ӵ�ֵ