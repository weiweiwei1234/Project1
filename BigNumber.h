//数据的基本运算，加减乘除、比较
#pragma once
#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdint>


#define BIGGER 1
#define SMALLER -1
#define EQUAL 0


// Big number class definition
class BigNumber {
private:
    // member
    bool sgn;
    std::vector<int8_t> data;     //int8_t:signed char 类型
    // private handle method
    static int abs_compare(const BigNumber&, const BigNumber&);
    static void discard_leading_zero(std::vector<int8_t>&);

public:
    // constructors
    BigNumber();
    BigNumber(long long); // directly convert from an int
    BigNumber(const std::string&);
    BigNumber(bool, const std::vector<int8_t>&);
    BigNumber(bool, std::vector<int8_t>&&);
    std::string get_value();
    int to_int();

    // overloaded logical operators as member functions
    friend bool operator==(const BigNumber&, const BigNumber&);
    friend bool operator!=(const BigNumber&, const BigNumber&);
    friend bool operator>(const BigNumber&, const BigNumber&);
    friend bool operator<(const BigNumber&, const BigNumber&);
    friend bool operator>=(const BigNumber&, const BigNumber&);
    friend bool operator<=(const BigNumber&, const BigNumber&);

    // overloaded arithmetic operators as member functions
    friend const BigNumber operator+(const BigNumber&, const BigNumber&);
    friend const BigNumber operator+=(const BigNumber&, const BigNumber&);
    friend const BigNumber operator-(const BigNumber&, const BigNumber&);
    friend const BigNumber operator*(const BigNumber&, const BigNumber&);
    friend const BigNumber operator/(const BigNumber&, const BigNumber&);
    friend const BigNumber operator%(const BigNumber&, const BigNumber&);

    // ouput format for BigNumber
    friend std::ostream& operator<<(std::ostream&, const BigNumber&);


};//负号分离，data部分只表示数值部分通过sgn来表示正负