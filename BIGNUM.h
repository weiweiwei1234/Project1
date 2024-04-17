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
class BIGNUM {
private:
    // member
    bool sgn;
    std::vector<int8_t> data;     //int8_t:signed char 类型
    // private handle method
    static int abs_compare(const BIGNUM&, const BIGNUM&);
    static void discard_leading_zero(std::vector<int8_t>&);

public:
    // constructors
    BIGNUM();
    BIGNUM(long long); // directly convert from an int
    BIGNUM(const std::string&);
    BIGNUM(bool, const std::vector<int8_t>&);
    BIGNUM(bool, std::vector<int8_t>&&);
    std::string get_value();
    int to_int();

    // overloaded logical operators as member functions
    friend bool operator==(const BIGNUM&, const BIGNUM&);
    friend bool operator!=(const BIGNUM&, const BIGNUM&);
    friend bool operator>(const BIGNUM&, const BIGNUM&);
    friend bool operator<(const BIGNUM&, const BIGNUM&);
    friend bool operator>=(const BIGNUM&, const BIGNUM&);
    friend bool operator<=(const BIGNUM&, const BIGNUM&);

    // overloaded arithmetic operators as member functions
    friend const BIGNUM operator+(const BIGNUM&, const BIGNUM&);
    friend const BIGNUM operator-(const BIGNUM&, const BIGNUM&);
    friend const BIGNUM operator*(const BIGNUM&, const BIGNUM&);
    friend const BIGNUM operator/(const BIGNUM&, const BIGNUM&);
    friend const BIGNUM operator%(const BIGNUM&, const BIGNUM&);

    // ouput format for BIGNUM
    friend std::ostream& operator<<(std::ostream&, const BIGNUM&);


};//负号分离，data部分只表示数值部分通过sgn来表示正负