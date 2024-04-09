#include "SM3.h"
#include "SM2.h"
//全局变量
extern BigInteger p = (char*)"115792089237316195423570985008687907853269984665640564039457584007913129639935";
extern BigInteger a = (char*)"115792089237316195423570985008687907853269984665640564039457584007913129639932";
extern BigInteger b = (char*)"41058363725152142129326129780047268409114441015993725554835256314039467401291";
extern BigInteger n = (char*)"115792089237316195423570985008687907853269984665640564039457584007913129639923";
extern BigInteger Gx = (char*)"602046282375688656758213480587526111916698976636884684818";
extern BigInteger Gy = (char*)"171943586981501236713555768347108872879680030001983424494933422011205721715898";
int main() {
	string str;  //禁止输入汉字
	cin >> str;
	cout << SM3(str) << endl;
    //椭圆曲线的参数
    EllipticCurveParams curve_params;
    init_elliptic_curve_params(&curve_params);
    assign_elliptic_curve_params(&curve_params, { p,a,b,n,{Gx,Gy} });
    //基点G
    ECPoint G;
    init_elliptic_curve_ecpoint(&G);
    assign_elliptic_curve_ecpoint(&G, curve_params.G);
	return 0;
}
/*
int main() {

    return 0;
}
*/