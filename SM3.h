#pragma once
#include <iostream>
#include <string>
#include <cmath>
using namespace std;
string BinToHex(string str);							//二进制转换为十六进制函数实现
string HexToBin(string str);							//十六进制转换为二进制函数实现
int BinToDec(string str);								//二进制转换为十进制的函数实现
string DecToBin(int str);								//十进制转换为二进制的函数实现
int HexToDec(string str);								//十六进制转换为十进制的函数实现
string DecToHex(int str);								//十进制转换为十六进制的函数实现
string padding(string str);								//对数据进行填充 
string LeftShift(string str, int len);					//实现循环左移len位功能
string XOR(string str1, string str2);					//实现异或操作
string AND(string str1, string str2);					//实现与操作
string OR(string str1, string str2);					//实现或操作
string NOT(string str);									//实现非操作
char binXor(char str1, char str2);						//实现单比特的异或操作
char binAnd(char str1, char str2);						//实现单比特的与操作
string ModAdd(string str1, string str2);				//mod 2^32运算的函数实现
string P1(string str);									//实现置换功能P1（X）
string P0(string str);									//实现置换功能P0（X）
string T(int j);										//返回Tj常量值的函数实现
string FF(string str1, string str2, string str3, int j);//实现布尔函数FF功能
string GG(string str1, string str2, string str3, int j);//实现布尔函数GG功能
string extension(string str);							//消息扩展函数
string compress(string str1, string str2);				//消息压缩函数
string iteration(string str);							//迭代压缩函数实现
string SM3(string str);									//生成杂凑值