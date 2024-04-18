#include "BIGNUM.h"
// constructor
BIGNUM::BIGNUM() {
    BIGNUM(0);
}

BIGNUM::BIGNUM(long long input_number) {
    long long unsign_number;

    // 判断正负
    sgn = !(input_number < 0);

    // 取数据绝对值
    unsign_number = (input_number < 0) ? -input_number : input_number;

    // 转换为16进制字符串
    while (unsign_number >= 16) {
        data.push_back(unsign_number & 15); // mod 16
        unsign_number = unsign_number >> 4; // div 16
    }
    data.push_back(unsign_number);
}

BIGNUM::BIGNUM(const std::string& input_string) {
    sgn = !(input_string.front() == '-');

    for (auto i = input_string.rbegin(), end = input_string.rend(); i != end; ++i) {
        if (*i >= '0' && *i <= '9') {
            data.push_back(*i - '0');
        }
        else if (*i >= 'A' && *i <= 'F') {
            data.push_back(*i - 'A' + 10);
        }
        else if (*i >= 'a' && *i <= 'f') {
            data.push_back(*i - 'a' + 10);
        }
    }
}

BIGNUM::BIGNUM(bool input_sgn, const std::vector<int8_t>& input_data) {
    sgn = input_sgn;
    data = input_data;
}

BIGNUM::BIGNUM(bool input_sgn, std::vector<int8_t>&& input_data) {
    sgn = input_sgn;
    data = std::move(input_data);
}


// logical operators
bool operator==(const BIGNUM& lhs, const BIGNUM& rhs) {
    return (lhs.sgn == rhs.sgn) && (BIGNUM::abs_compare(lhs, rhs) == EQUAL);
}

bool operator!=(const BIGNUM& lhs, const BIGNUM& rhs) {
    return !(lhs == rhs);
}

bool operator>(const BIGNUM& lhs, const BIGNUM& rhs) {
    int abs_cmp;

    if (lhs.sgn == rhs.sgn) {
        abs_cmp = BIGNUM::abs_compare(lhs, rhs);
        return ((lhs.sgn && abs_cmp == BIGGER) || (!lhs.sgn && abs_cmp == SMALLER));
    }
    else {
        return lhs.sgn;
    }
}

bool operator<(const BIGNUM& lhs, const BIGNUM& rhs) {
    return rhs > lhs;
}

bool operator>=(const BIGNUM& lhs, const BIGNUM& rhs) {
    return !(lhs < rhs);
}

bool operator<=(const BIGNUM& lhs, const BIGNUM& rhs) {
    return !(lhs > rhs);
}


// 判断相等
int BIGNUM::abs_compare(const BIGNUM& lhs, const BIGNUM& rhs) {
    if (lhs.data.size() > rhs.data.size()) {
        return BIGGER;
    }
    else if (lhs.data.size() < rhs.data.size()) {
        return SMALLER;
    }

    // same size
    for (auto i = lhs.data.rbegin(), j = rhs.data.rbegin(), end = lhs.data.rend(); i != end; ++i, ++j) {
        if (*i > *j) {
            return BIGGER;
        }
        else if (*i < *j) {
            return SMALLER;
        }
    }

    return EQUAL;
}
//消除前置0
void BIGNUM::discard_leading_zero(std::vector<int8_t>& input) {
    while (input.back() == 0 && input.size() != 1) {
        input.pop_back();
    }
}

// arithmetic operators
const BIGNUM operator+(const BIGNUM& lhs, const BIGNUM& rhs) {
    bool sgn;
    std::vector<int8_t> abs_result;
    unsigned long min_size;
    int8_t carry, sum;
    //min_size取数据长度较小的那个
    min_size = (lhs.data.size() < rhs.data.size()) ? lhs.data.size() : rhs.data.size();

    // 同号相加，异号相减
    if (lhs.sgn == rhs.sgn) {
        // 同号
        sgn = lhs.sgn;
        carry = 0;

        // 做加法
        for (unsigned long i = 0; i < min_size; i++) {
            sum = lhs.data[i] + rhs.data[i];
            abs_result.push_back(sum);
        }

        // 插入剩余的数
        if (lhs.data.size() > rhs.data.size()) {
            for (unsigned long i = min_size; i < lhs.data.size(); i++) {
                abs_result.push_back(lhs.data[i]);
            }
        }
        else {
            for (unsigned long i = min_size; i < rhs.data.size(); i++) {
                abs_result.push_back(rhs.data[i]);
            }
        }

        // 处理进位
        carry = 0;
        for (unsigned long i = 0; i < abs_result.size(); i++) {
            sum = abs_result[i] + carry;
            abs_result[i] = sum & 15;
            carry = sum >> 4;
        }
        if (carry > 0) {
            abs_result.push_back(carry);
        }
    }
    else {
        //异号 相当于做减法
        if (!lhs.sgn) {
            return rhs - BIGNUM(!lhs.sgn, lhs.data);
        }
        else {
            return lhs - BIGNUM(!rhs.sgn, rhs.data);
        }
    }

    return BIGNUM(sgn, std::move(abs_result));
}

const BIGNUM operator-(const BIGNUM& lhs, const BIGNUM& rhs) {
    bool sgn;
    std::vector<int8_t> abs_result;
    unsigned long min_size;

    // 同号相加，异号相减
    if (lhs.sgn != rhs.sgn) {
        // 异号 相当于做加法
        return lhs + BIGNUM(!rhs.sgn, rhs.data);
    }
    else {
        if (BIGNUM::abs_compare(lhs, rhs) == EQUAL) {
            return BIGNUM(0);
        }
        else if (BIGNUM::abs_compare(lhs, rhs) == BIGGER) {
            sgn = lhs.sgn;
            // 做减法
            min_size = rhs.data.size();
            abs_result.resize(lhs.data.size());
            std::transform(lhs.data.begin(), lhs.data.begin() + min_size, rhs.data.begin(), abs_result.begin(), std::minus<int8_t>());
            // 插入剩余的数
            std::copy(lhs.data.begin() + min_size, lhs.data.end(), abs_result.begin() + min_size);
        }
        else {
            sgn = !rhs.sgn;
            // 做减法
            min_size = lhs.data.size();
            abs_result.resize(rhs.data.size());
            std::transform(rhs.data.begin(), rhs.data.begin() + min_size, lhs.data.begin(), abs_result.begin(), std::minus<int8_t>());
            // 插入剩余的数
            std::copy(rhs.data.begin() + min_size, rhs.data.end(), abs_result.begin() + min_size);
        }
        // 处理借位
        for (unsigned long i = 0; i < abs_result.size(); i++) {
            if (abs_result[i] < 0) {
                abs_result[i] += 16;
                abs_result[i + 1]--;
            }
        }

        //消除前置0
        BIGNUM::discard_leading_zero(abs_result);
    }

    return BIGNUM(sgn, std::move(abs_result));
}

const BIGNUM operator*(const BIGNUM& lhs, const BIGNUM& rhs) {
    bool sgn;
    std::vector<int8_t> abs_result;

    if (lhs == 0 || rhs == 0) {
        return BIGNUM(0);
    }

    sgn = (lhs.sgn == rhs.sgn);

    abs_result.resize(lhs.data.size() + rhs.data.size(), 0);
    int16_t sum = 0;
    for (unsigned long i = 0; i < lhs.data.size(); i++) {
        for (unsigned long j = 0; j < rhs.data.size(); j++) {
            sum = abs_result[i + j] + lhs.data[i] * rhs.data[j];
            abs_result[i + j] = sum % 16; // sum
            abs_result[i + j + 1] += sum / 16; // carry
        }
    }
    //消除前置0
    BIGNUM::discard_leading_zero(abs_result);

    return BIGNUM(sgn, abs_result);
}

const BIGNUM operator/(const BIGNUM& lhs, const BIGNUM& rhs) {
    BIGNUM temp(0);
    BIGNUM remainder(true, lhs.data);
    BIGNUM divisor(true, rhs.data);
    BIGNUM quotient(0);
    quotient.sgn = (lhs.sgn == rhs.sgn);

    int8_t count = 0;
    while (remainder.data.size() != 0) {
        while (temp < divisor && remainder.data.size() != 0) {
            temp.data.insert(temp.data.begin(), remainder.data.back());
            remainder.data.pop_back();
            quotient.data.insert(quotient.data.begin(), 0);
        }
        BIGNUM::discard_leading_zero(temp.data);

        count = 0;
        while (temp >= divisor) {
            count++;
            temp = temp - divisor;
        }
        quotient.data.front() = count;
    }
    BIGNUM::discard_leading_zero(quotient.data);

    // -0 -> +0
    if (quotient.data.size() == 1 && quotient.data.back() == 0) {
        quotient.sgn = true;
    }

    return quotient;
}

const BIGNUM operator%(const BIGNUM& lhs, const BIGNUM& rhs) {
    BIGNUM temp(0);
    BIGNUM remainder(true, lhs.data);
    BIGNUM divisor(true, rhs.data);

    if (BIGNUM::abs_compare(lhs, rhs) == EQUAL) {
        return BIGNUM(0);
    }
    if (BIGNUM::abs_compare(lhs, rhs) == SMALLER) {
        return lhs;
    }
    while (remainder >= divisor) {
        temp.data.assign(remainder.data.end() - divisor.data.size(), remainder.data.end());
        remainder.data.erase(remainder.data.end() - divisor.data.size(), remainder.data.end());
        if (temp < divisor) {
            temp.data.insert(temp.data.begin(), remainder.data.back());
            remainder.data.pop_back();
        }

        BIGNUM::discard_leading_zero(temp.data);

        while (temp >= divisor) {
            temp = temp - divisor;
        }
        remainder.data.insert(remainder.data.end(), temp.data.begin(), temp.data.end());
        temp.data.clear();
    }
    BIGNUM::discard_leading_zero(remainder.data);
    remainder.sgn = lhs.sgn;
    // -0 -> +0
    if (remainder.data.size() == 1 && remainder.data.back() == 0) {
        remainder.sgn = true;
    }
    return remainder;
}

//const BIGNUM operator^(const BIGNUM& lhs, const BIGNUM& rhs)
//{
//   //rhs>1
//    BIGNUM base = lhs, exponent = rhs;
//    BIGNUM result = base;
//    if (exponent == 1) return result;
//    while (exponent > 1) {
//        result = result * base;
//    }
//    return result;
//}
// 
//右移 k 位 相当于除以2^k
const BIGNUM operator>>(const BIGNUM& lhs, const int& rhs) {
    BIGNUM temp(lhs.sgn, lhs.data);
    std::string str = temp.HexToBin();
    if (str.length() <= rhs)
        return BIGNUM(0);
    str = str.substr(0, str.length() - rhs);

    std::string hex = "";//用来存储最后生成的十六进制数
    int t = 0;//用来存储每次四位二进制数的十进制值
    while (str.size() % 4 != 0) {//因为每四位二进制数就能够成为一个十六进制数，所以将二进制数长度转换为4的倍数
        str = "0" + str;//最高位添0直到长度为4的倍数即可
    }
    for (int i = 0; i < str.size(); i += 4) {
        t = (str[i] - '0') * 8 + (str[i + 1] - '0') * 4 + (str[i + 2] - '0') * 2 + (str[i + 3] - '0') * 1;//判断出4位二进制数的十进制大小为多少
        if (t < 10) {//当得到的值小于10时，可以直接用0-9来代替
            hex += to_string(t);
        }
        else {//当得到的值大于10时，需要进行A-F的转换
            hex += 'A' + (t - 10);
        }
    }
    if (temp.sgn == false)
        hex = "-" + hex;
    return BIGNUM(hex);
}
//左移 k 位 相当于乘以2^k
const BIGNUM operator<<(const BIGNUM& lhs, const int& rhs) {
    BIGNUM temp(lhs.sgn, lhs.data);
    std::string str = temp.HexToBin();
    for (int i = 0; i < rhs; i++) {
        str = str + "0";
    }

    std::string hex = "";//用来存储最后生成的十六进制数
    int t = 0;//用来存储每次四位二进制数的十进制值
    while (str.size() % 4 != 0) {//因为每四位二进制数就能够成为一个十六进制数，所以将二进制数长度转换为4的倍数
        str = "0" + str;//最高位添0直到长度为4的倍数即可
    }
    for (int i = 0; i < str.size(); i += 4) {
        t = (str[i] - '0') * 8 + (str[i + 1] - '0') * 4 + (str[i + 2] - '0') * 2 + (str[i + 3] - '0') * 1;//判断出4位二进制数的十进制大小为多少
        if (t < 10) {//当得到的值小于10时，可以直接用0-9来代替
            hex += to_string(t);
        }
        else {//当得到的值大于10时，需要进行A-F的转换
            hex += 'A' + (t - 10);
        }
    }
    if (temp.sgn == false)
        hex = "-" + hex;
    return BIGNUM(hex);
}

// output format
std::ostream& operator<<(std::ostream& os, const BIGNUM& rhs) {
    if (!rhs.sgn) {
        os << "-";
    }

    for (auto i = rhs.data.rbegin(); i != rhs.data.rend(); ++i) {
        // i is a pointer, point to a certain position in the vector (rhs.data)
        // *i is the value store in THAT position
        if (*i >= 10) {
            os << static_cast<char>(*i + 'A' - 10); // 10 -> 'a'
        }
        else {
            os << static_cast<char>(*i + '0'); // 1 -> '1'
        }
    }

    return os;
}

std::string BIGNUM::get_value() {
    std::string A = "";

    if (!sgn) {
        A[0] = '-';
    }
    for (auto i = data.rbegin(); i != data.rend(); ++i) {
        // i is a pointer, point to a certain position in the vector (rhs.data)
        // *i is the value store in THAT position
        if (*i >= 10) {
            A += static_cast<char>(*i + 'A' - 10); // 10 -> 'a'
        }
        else {
            A += static_cast<char>(*i + '0'); // 1 -> '1'
        }
    }
    return A;
}//十六进制数字符串输出

int BIGNUM::to_int() {
    int sum = 0, n = 0;
    for (auto i = data.begin(); i != data.end(); ++i) {
        // i is a pointer, point to a certain position in the vector (rhs.data)
        // *i is the value store in THAT position
        if (*i >= 10) {
            sum += *i * pow(16, n); // 10 -> 'a'
        }
        else {
            sum += *i * pow(16, n); // 1 -> '1'
        }
        n++;
    }
    if (sgn == 0) {
        return -sum;
    }
    else {
        return sum;
    }
}
int BIGNUM::bitlen()
{
    int length = HexToBin().length();
    return length;
}
std::string BIGNUM::HexToBin()
{
    std::string str = get_value();
    std::string bin = "";
    std::string table[16] = { "0000","0001","0010","0011","0100","0101","0110","0111","1000","1001","1010","1011","1100","1101","1110","1111" };
    for (int i = 0; i < str.size(); i++) {
        if (str[i] >= 'A' && str[i] <= 'F') {
            bin += table[str[i] - 'A' + 10];
        }
        else {
            bin += table[str[i] - '0'];
        }
    }
    while (bin[0]=='0') {
            bin.erase(0, 1);
    }
    return bin;
}
//十六进制数转十进制