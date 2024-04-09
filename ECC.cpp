#include "SM3.h"
#include "SM2.h"
int main() {//主函数
	string str;
	cout << "输入消息:" << endl;
	cin >> str;
	cout << "杂凑值:\n" << SM3(str) << endl;
	return 0;
}