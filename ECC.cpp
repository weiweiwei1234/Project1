#include "SM3.h"
#include "SM2.h"
int main() {//������
	string str;
	cout << "������Ϣ:" << endl;
	cin >> str;
	cout << "�Ӵ�ֵ:\n" << SM3(str) << endl;
	return 0;
}