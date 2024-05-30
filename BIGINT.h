#pragma once
#include<iostream>
#include<string>
#include<cmath>

#define BIGGER 1
#define SMALLER -1
#define EQUAL 0

using namespace std;

//大整数数据类型，带模数 （data,p)
class BIGINT {
private:
	bool sgn;					//符号	实际中均取正值 
	uint8_t data[32] = { 0 };	//32位字节 逆序储存在数组中 uint8_t 0-256（00-FF)
	uint8_t p[32] = { 0 };		//有限域p
	int k;						//字节串大小
	int m;						//比特串大小
	static int compare(const BIGINT&, const BIGINT&);//判断两个数是否相等
public:
	BIGINT();				 //默认构造函数
	BIGINT(long,long);		 //使用较小的整数构造大整数 几乎不使用
	BIGINT(bool, BIGINT);	 //一般用于取负值
	BIGINT(bool, string);	 //使用字符串构造大整数 确定的正负值
	BIGINT(string);			 //不确定正负值
	int get_byte_length();	 //获取字节串大小
	int get_bite_length();	 //获取比特串大小
	BIGINT inverse(BIGINT m);//求逆 a mod m 的乘法逆元 
	int to_int();			 //转10进制整数，数值很小的时候
	string to_bite();		 //比特表示 返回字符串

	//重载逻辑运算符
	friend bool operator==(const BIGINT&, const BIGINT&);
	friend bool operator!=(const BIGINT&, const BIGINT&);
	friend bool operator>(const BIGINT&, const BIGINT&);
	friend bool operator<(const BIGINT&, const BIGINT&);
	friend bool operator>=(const BIGINT&, const BIGINT&);
	friend bool operator<=(const BIGINT&, const BIGINT&);

	// 重载算术运算符
	friend const BIGINT operator+(const BIGINT&, const BIGINT&);
	friend const BIGINT operator-(const BIGINT&, const BIGINT&);
	friend const BIGINT operator*(const BIGINT&, const BIGINT&);
	friend const BIGINT operator/(const BIGINT&, const BIGINT&); //整除
	friend const BIGINT operator%(const BIGINT&, const BIGINT&); //求余 mod
	friend const BIGINT operator>>(const BIGINT&, const int&);   //右移
	friend const BIGINT operator<<(const BIGINT&, const int&);   //左移

	// 重载输出
	friend std::ostream& operator<<(std::ostream&, const BIGINT&);
};