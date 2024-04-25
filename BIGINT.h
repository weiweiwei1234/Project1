#pragma once
#include<iostream>
#include<string>

#define BIGGER 1
#define SMALLER -1
#define EQUAL 0

using namespace std;
class BIGINT {
private:
	bool sgn;			//符号	实际中均取正值
	int8_t data[64];	//储存64位16进制数 相当于256位大整数	逆序储存在数组中
	int8_t p[64];		//整数域上的运算 模数
	int size;			//整数的大小 
	static int compare(const BIGINT&, const BIGINT&);//判断两个数是否相等
public:

};