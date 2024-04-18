//数据的基本运算，加、减、乘、除、mod、比较
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
    // 成员数据
    bool sgn;   //true为+，false为-
    std::vector<int8_t> data;     //int8_t:signed char 类型 将数据的绝对值倒着存储
    // 判断相等 消除前置0
    static int abs_compare(const BIGNUM&, const BIGNUM&);
    static void discard_leading_zero(std::vector<int8_t>&);

public:
    BIGNUM();
    BIGNUM(long long); // 直接转换int
    BIGNUM(const std::string&);
    BIGNUM(bool, const std::vector<int8_t>&);
    BIGNUM(bool, std::vector<int8_t>&&);
    std::string get_value(); //整数值
    int to_int(); //小数转int
    int bitlen(); //二进制长度
    std::string HexToBin();//十六进制转二进制

    // 重载逻辑运算符
    friend bool operator==(const BIGNUM&, const BIGNUM&); 
    friend bool operator!=(const BIGNUM&, const BIGNUM&);
    friend bool operator>(const BIGNUM&, const BIGNUM&);
    friend bool operator<(const BIGNUM&, const BIGNUM&);
    friend bool operator>=(const BIGNUM&, const BIGNUM&);
    friend bool operator<=(const BIGNUM&, const BIGNUM&);

    // 重载算术运算符
    friend const BIGNUM operator+(const BIGNUM&, const BIGNUM&);
    friend const BIGNUM operator-(const BIGNUM&, const BIGNUM&);
    friend const BIGNUM operator*(const BIGNUM&, const BIGNUM&);
    friend const BIGNUM operator/(const BIGNUM&, const BIGNUM&); //整除
    friend const BIGNUM operator%(const BIGNUM&, const BIGNUM&); //求余 mod
    friend const BIGNUM operator>>(const BIGNUM&, const int&);   //右移
    friend const BIGNUM operator<<(const BIGNUM&, const int&);   //左移

    // 重载输出
    friend std::ostream& operator<<(std::ostream&, const BIGNUM&);


};//负号分离，data部分只表示数值部分通过sgn来表示正负