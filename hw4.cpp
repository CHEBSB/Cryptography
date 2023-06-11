#include <stdio.h>
#include <stdlib.h>
#define L 84

class BigNum;
class Elem;

class BigNum {              // big number
    friend class Elem;
    friend class point;
private:
    int sign;               // its sign (+1 / -1)
    int length;             // length of the number in hex
    int S[L];               // sequence of number
public:
    BigNum();
    BigNum(const int&);     // directly from a integer
    BigNum(const Elem&);    // from a field element
    BigNum(const BigNum&);  // copy another big number
    BigNum(int, int, int*);
    // arithmetic: assume the inputs are positive
    BigNum operator+(const BigNum&);
    BigNum operator+(const int&);
    BigNum operator-(const BigNum&);
    BigNum operator-(const int&);
    BigNum operator*(const BigNum&);
    BigNum operator*(const int&);
    BigNum operator/(BigNum&);
    BigNum operator/(const int&);
    BigNum operator%(const BigNum&);
    bool operator==(const BigNum&) const;      
    bool operator>(const BigNum&) const;
    // compare abs value (true for >=)
    bool abscmp(const BigNum&) const; 
    bool geqP() const;      // whether this >= p
    // modular exponentiation
    BigNum modExp(const BigNum& b, BigNum& n);
    BigNum Lshift(int z) const; // left-shift z positions
    void lenCheck();        // renew one's length
    void print();           // print the number in hex
    void printPlain();      // print w/out padded bits
};
class Elem {                // finite field element
    friend class BigNum;
    friend class point;
private:
    int S[40];              // sequence of number
public:
    Elem();
    Elem(const int&);       // directly from a integer
    Elem(const Elem&);      // copy another element
    Elem(const BigNum&);    // from a big number
    Elem(int, int*);
    // arithmetic: assume the inputs are positive
    Elem operator+(const Elem&);
    Elem operator-(const Elem&);
    Elem operator*(const Elem&);
    Elem operator/(const Elem&);  
    int len() const;        // length count in hex
    // modular exponentiation (mod p)
    Elem modExp(const Elem& b);
    void print();           // print out in hex
};
class point {                   // point on elliptic curve
private:
    Elem x;                  // x coordinate
    Elem y;                  // y coordinate
public:
    point();   
    point(const Elem&, const Elem&);    // given x, y
    point(const point&);        // copy constructor
    point operator+(point&);
    point mulBy2();             // multiply the point by 2
    point operator*(const int&);
    void print();               // print out x and y in hex
};
// elliptic curve parameters
int P[] = {0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF,
0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF,
0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF,
0xF, 0xF, 0x7, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF};
int A[] = {0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF,
0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF,
0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF,
0xF, 0xF, 0x7, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xC};
int B[] = {0x1, 0xC, 0x9, 0x7, 0xB, 0xE, 0xF, 0xC,
0x5, 0x4, 0xB, 0xD, 0x7, 0xA, 0x8, 0xB, 0x6, 0x5, 0xA,
0xC, 0xF, 0x8, 0x9, 0xF, 0x8, 0x1, 0xD, 0x4, 0xD, 0x4,
0xA, 0xD, 0xC, 0x5, 0x6, 0x5, 0xF, 0xA, 0x4, 0x5};
int Gx[] = {0x4, 0xA, 0x9, 0x6, 0xB, 0x5, 0x6, 0x8,
0x8, 0xE, 0xF, 0x5, 0x7, 0x3, 0x2, 0x8, 0x4, 0x6, 0x6,
0x4, 0x6, 0x9, 0x8, 0x9, 0x6, 0x8, 0xC, 0x3, 0x8, 0xB,
0xB, 0x9, 0x1, 0x3, 0xC, 0xB, 0xF, 0xC, 0x8, 0x2};
int Gy[] = {0x2, 0x3, 0xA, 0x6, 0x2, 0x8, 0x5, 0x5,
0x3, 0x1, 0x6, 0x8, 0x9, 0x4, 0x7, 0xD, 0x5, 0x9, 0xD,
0xC, 0xC, 0x9, 0x1, 0x2, 0x0, 0x4, 0x2, 0x3, 0x5, 0x1,
0x3, 0x7, 0x7, 0xA, 0xC, 0x5, 0xF, 0xB, 0x3, 0x2};
int N[] = {0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0xF, 4, 0xC, 8, 0xF, 9,
2, 7, 0xA, 0xE, 0xD, 3, 0xC, 0xA, 7, 5, 2, 2, 5, 7};
const BigNum p(1, 40, P);          // p == 3 mod 4
const BigNum a(1, 40, A);
const BigNum b(1, 40, B);
const BigNum gx(1, 40, Gx);
const BigNum gy(1, 40, Gy);
const BigNum n(1, 42, N);

// for Euclid's algorithm
BigNum uj[3];                 
BigNum rj[3];
BigNum vj[3];
BigNum qj[3];

// Extended Euclid's Algorithm to find inverse in GF(p)
int EuclidAlg(const BigNum&);

int main(void)
{
    BigNum a = 5;
    (a.Lshift(3)).print();

    return 0;
}

// default constructor: init as 0
BigNum::BigNum(void)
{          
    int i;                      // looping index

    sign = 1;
    length = 1;
    for (i = 0; i < L; i++) 
        S[i] = 0;
}
// copy constructor
BigNum::BigNum(const BigNum& another)
{
    int i;                      // looping index

    sign = another.sign;
    length = another.length;
    for (i = 0; i < L; i++) 
        S[i] = another.S[i];
}
// from a finite field element
BigNum::BigNum(const Elem& another)
{
    int i;                      // looping index

    sign = 1;
    length = 40;
    for (i = 0; i < 40; i++) 
        S[i] = another.S[i];
    for (; i < L; i++)
        S[i] = 0;
}
// convert an integer into BigNum
BigNum::BigNum(const int& n)
{  
    int i = 0;                  // looping index
    int K;

    sign = (n >= 0)? 1: (-1);
    K = n * sign;               // K = |n|
    while (K != 0 && i < L) {   // starting from LSB
        S[i] = K % 16;
        K = (K - S[i]) / 16;
        i += 1;
    }
    if (i >= L) printf("An integer too large to convert.\n");
    length = i;
    for (; i < L; i++)
        S[i] = 0;
    if (n == 0) length = 1;     // handle the case n = 0
}
// detail-specify constructor
BigNum::BigNum(int sgn, int l , int* seq)
{
    int i;                      // looping index

    sign = sgn;
    length = l;
    for (i = 0; i < l; i++) 
        S[i] = seq[l - 1 - i];
    for (i = l; i < L; i++) 
        S[i] = 0;
}
void BigNum::lenCheck(void)
{
    int i;                      // looping index

    for (i = L - 1; i > 0; i--) // find largest nonzero digit
        if (this->S[i] != 0) {
            this->length = i + 1;
            return;
        }
    this->length = 1;
}
// modular exponentiation: compute (a^b) mod n
BigNum BigNum::modExp(const BigNum& b, BigNum& n)
{
    int* b_bit;             // b in binary representation
    int k = 4 * (b.length);
    int i, j;               // looping indices
    int temp;               // temp storage
    BigNum d(1);            // d = 1

    b_bit = (int*)calloc(k, sizeof(int));
    // convert b to bit array
    for (i = 0, j = 0; i < b.length; i++, j += 4) {
        temp = b.S[i]; 
        b_bit[j] = temp % 2;
        temp = temp / 2;
        b_bit[j + 1] = temp % 2;
        temp = temp / 2;
        b_bit[j + 2] = temp % 2;
        temp = temp / 2;
        b_bit[j + 3] = temp % 2;
    }
    // modular exponentiation
    for (i = k - 1, j = 0; i >= 0; i--) {
        j = 2 * j;
        d = (d * d) % n;
        if (b_bit[i] == 1) {
            j += 1;
            d = (d * (*this)) % n;
        }
    }
    free(b_bit);
    return d;
}
// multiply the number by 16 to the power of z
BigNum BigNum::Lshift(int z) const
{
    int i;                      // looping index
    BigNum newB;

    newB.length = this->length + z;
    newB.sign = this->sign;
    for (i = 0; i < z; i++) 
        newB.S[i] = 0;
    for (; i < this->length + z; i++)
        newB.S[i] = this->S[i - z];
    for (; i < L; i++)
        newB.S[i] = 0;
    return newB;
}
// compare absolute value (return true when >=)
bool BigNum::abscmp(const BigNum& another) const
{
    int i;                      // looping index

    if (this->length > another.length)  return true;
    if (this->length < another.length)  return false;
    // having same length                  
    for (i = this->length; i >= 0; i--) {   // digit-by-digit
        if (this->S[i] > another.S[i])  return true;
        if (this->S[i] < another.S[i])  return false;
    }
    return true;                // they are equal
}
// whether this number >= p
bool BigNum::geqP() const            
{
    int i;                  // loopiog index

    if (this->sign == -1) printf("Negative comparison!\n");
    if (this->length > 40) return true;
    if (this->length < 40) return false;
    for (i = 39; i >= 0; i--) 
        if (i != 7) 
            if (this->S[i] < 0xf) 
                return false;
        else {
            if (this->S[i] < 7) return false;
            if (this->S[i] > 7) return true;
        }
    return true;        // they are equal
}
bool BigNum::operator>(const BigNum& another) const {
    int i;                      // looping index

    if (this->sign == 1 && another.sign == -1) return true;
    if (this->sign == -1 && another.sign == 1) return false;
    if (this->sign == 1) {      // both are positive
        if (this->length > another.length)  return true;
        if (this->length < another.length)  return false;
        // having same bytes                  
        for (i = this->length; i >= 0; i--) {
            if (this->S[i] > another.S[i])  return true;
            if (this->S[i] < another.S[i])  return false;
        }
        return false;           // they are equal
    }
    // both negative
    if (this->length > another.length)  return false;
    if (this->length < another.length)  return true;             
    for (i = this->length; i >= 0; i--) {
        if (this->S[i] > another.S[i])  return false;
        if (this->S[i] < another.S[i])  return true;
    }
    return false;               // they are equal
}
bool BigNum::operator==(const BigNum& another) const
{
    int i;                      // looping index

    if (this->sign != another.sign)
        return false;
    if (this->length != another.length)
        return false;                 
    for (i = this->length; i >= 0; i--) {
        if (this->S[i] != another.S[i])
            return false;
    }
    return true;               // they are equal
}
BigNum BigNum::operator+(const BigNum& another)
{
    int i;                      // looping index
    int maxL;                   // larger length of the 2
    BigNum sum;                 // sum of the 2 numebrs

    maxL = (this->length > another.length)? this->length: another.length;
    if (this->sign == another.sign) {
        sum.sign = this->sign;
        for (i = 0; i < maxL; i++) {
            sum.S[i] += (this->S[i] + another.S[i]);
            if (sum.S[i] > 15) {
                sum.S[i] -= 16;
                sum.S[i + 1] += 1;
            }
        }  
    } else if (this->sign == 1) {   // this > 0, another < 0
        if (this->abscmp(another)) {        // result > 0
            sum.sign = 1;
            for (i = 0; i < maxL; i++) {
                sum.S[i] += (this->S[i] - another.S[i]);
                if (sum.S[i] < 0) {
                    sum.S[i] += 16;
                    sum.S[i + 1] -= 1;
                }
            }
        } else {                            // result < 0
            sum.sign = -1;
            for (i = 0; i < maxL; i++) {
                sum.S[i] += (another.S[i] - this->S[i]);
                if (sum.S[i] < 0) {
                    sum.S[i] += 16;
                    sum.S[i + 1] -= 1;
                }
            }
        }
    } else {                        // this < 0, another > 0
        if (this->abscmp(another)) {        // result < 0
            sum.sign = -1;
            for (i = 0; i < maxL; i++) {
                sum.S[i] += (this->S[i] - another.S[i]);
                if (sum.S[i] < 0) {
                    sum.S[i] += 16;
                    sum.S[i + 1] -= 1;
                }
            }
        } else {                            // result > 0
            sum.sign = 1;
            for (i = 0; i < maxL; i++) {
                sum.S[i] += (another.S[i] - this->S[i]);
                if (sum.S[i] < 0) {
                    sum.S[i] += 16;
                    sum.S[i + 1] -= 1;
                }
            }
        }
    }
    sum.lenCheck();
    return sum;
}
BigNum BigNum::operator+(const int& b)
{
    BigNum B(b);
    return ((*this) + B);
}
BigNum BigNum::operator-(const BigNum& another)
{
    int i;                      // looping index
    int maxL;                   // = max{Length1, Length2}
    BigNum result;              // = (this - another)

    maxL = (this->length > another.length)? this->length: another.length;
    if ((*this) > another) {    // result is positive
        result.sign = 1;
        // this > 0, another < 0
        if (this->sign != another.sign) {
            for (i = 0; i < maxL; i++) {
                result.S[i] += (this->S[i] + another.S[i]);
                if (result.S[i] > 15) {
                    result.S[i] -= 16;
                    result.S[i + 1] += 1;
                }
            }
        } else if (this->sign == 1) {   // both > 0
            for (i = 0; i < maxL; i++) {
                result.S[i] += (this->S[i] - another.S[i]);
                if (result.S[i] < 0) {
                    result.S[i] += 16;
                    result.S[i + 1] -= 1;
                }
            }
        } else {                        // both < 0
            for (i = 0; i < maxL; i++) {
                result.S[i] += (another.S[i] - this->S[i]);
                if (result.S[i] < 0) {
                    result.S[i] += 16;
                    result.S[i + 1] -= 1;
                }
            }
        }
    } else {            // result is negative
        result.sign = -1;
        // this < 0, another > 0
        if (this->sign != another.sign) {
            for (i = 0; i < maxL; i++) {
                result.S[i] += (another.S[i] + this->S[i]);
                if (result.S[i] > 15) {
                    result.S[i] -= 16;
                    result.S[i + 1] += 1;
                }
            }
        } else if (this->sign == 1) {   // both > 0
            for (i = 0; i < maxL; i++) {
                result.S[i] += (another.S[i] - this->S[i]);
                if (result.S[i] < 0) {
                    result.S[i] += 16;
                    result.S[i + 1] -= 1;
                }
            }
        } else {                        // both < 0
            for (i = 0; i < maxL; i++) {
                result.S[i] += (this->S[i] - another.S[i]);
                if (result.S[i] < 0) {
                    result.S[i] += 16;
                    result.S[i + 1] -= 1;
                }
            }
        }
    }
    result.lenCheck();
    return result;
}
BigNum BigNum::operator-(const int& b)
{
    BigNum B(b);
    return ((*this) - B);
}
BigNum BigNum::operator*(const BigNum& another)
{
    int i, j;                   // looping indices
    BigNum prod;                // product of the 2

    prod.sign = (this->sign) * (another.sign);
    // shift-and-add multiplication
    for (i = 0; i < another.length; i++) {
        for (j = 0; j < this->length; j++) {
            prod.S[j + i] += ((another.S[i]) * (this->S[j]));
            if (prod.S[j + i] > 15) {
                prod.S[j + i + 1] += (prod.S[j + i] / 16);
                prod.S[j + i] = prod.S[j + i] % 16;
            }
        }
    }
    prod.lenCheck();
    return prod;
}
BigNum BigNum::operator*(const int& b)
{
    BigNum B(b);
    return ((*this) * B);
}
BigNum BigNum::operator/(BigNum& another)
{
    int i, j;                   // looping index
    BigNum r(*this);            // remainder
    BigNum temp(another);       // temporary variable
    BigNum q;                   // quotient

    if (another == BigNum(0)) {
        printf("Division by zero!\n");
        return 0;
    }
    i = -1;
    for (i = -1; this->abscmp(temp); i++)
        temp = temp.Lshift(1);  // go to next digit
    if (i == -1)                // if this < another         
        return q;               // return 0
    // if i = k, then this >= (temp left shifted by k)
    q.length = i + 1;
    for (j = i; j >= 0; j--) {
        temp = another.Lshift(j);
        while (r.abscmp(temp)) {
            r = r - temp;
            q.S[j] += 1;
        }
    }
    return q;
}   
BigNum BigNum::operator/(const int& b)
{
    BigNum B(b);
    return ((*this) / B);
}
BigNum BigNum::operator%(const BigNum& another)
{
    int i, j;                   // looping index
    BigNum r(*this);            // remainder
    BigNum temp(another);       // temporary variable

    if (this->sign == 1) {
        for (j = this->length - another.length; j >= 0; j--) {
            temp = another.Lshift(j);
            while (r.abscmp(temp)) 
                r = r - temp;
        }
    } else {                    // if this < 0
        r.sign = 1;
        for (j = this->length - another.length; j >= 0; j--) {
            temp = another.Lshift(j);
            while (r.abscmp(temp)) 
                r = r - temp;
        }
        r = temp - r;
    }
    r.lenCheck();
    return r;
}
void BigNum::print(void)        // print the number in hex
{      
    int i;                      // looping index
    int k;                      // counter

    if (sign == -1) printf("-");
    for (i = length - 1, k = 0; i >= 0; i--, k++) {
        if (k == 8) {
            printf(" ");
            k = 0;
        }
        printf("%x", S[i]);     // starting from MSB
    }
    printf("\n");
}
// Extended Euclid's Algorithm stop at rj = 1
int EuclidAlg(const BigNum& q)
{
    int i;          // looping index
    BigNum uni(1);  // 1

    // u[-1] = 1    v[-1] = 0   r[-1] = p
    // u[0] = 0     v[0] = 1    r[0] = q
    uj[0] = 1;
    vj[0] = 0;
    rj[0] = p;
    uj[1] = 0;
    vj[1] = 1;
    rj[1] = q;
    i = 0;
    while (!(rj[(i + 1) % 3] == uni)) {
        // find qj
        qj[(i + 2) % 3] = rj[i] / rj[(i + 1) % 3];
        // update rj, uj, and vj
        rj[(i + 2) % 3] = rj[i] - qj[(i + 2) % 3] * rj[(i + 1) % 3];
        uj[(i + 2) % 3] = uj[i] - qj[(i + 2) % 3] * uj[(i + 1) % 3]; 
        vj[(i + 2) % 3] = vj[i] - qj[(i + 2) % 3] * vj[(i + 1) % 3];
        // index increment
        i = (i + 1) % 3;
    } 
    i = (i + 1) % 3;
    return i;
}
// default constructor: init as 0
Elem::Elem(void)
{          
    int i;                      // looping index

    for (i = 0; i < 40; i++) 
        S[i] = 0;
}
// copy constructor
Elem::Elem(const Elem& another)
{
    int i;                      // looping index

    for (i = 0; i < 40; i++) 
        S[i] = another.S[i];
}
// construct from a big number
Elem::Elem(const BigNum& another)
{
    int i;                      // looping index

    if (another.sign == -1)
        printf("Negative input!\n");
    if (another.length > 40)
        printf("A big number too long!\n");
    for (i = 0; i < 40; i++) 
        S[i] = another.S[i];
}
// convert an integer into Elem
Elem::Elem(const int& n)
{  
    int i = 0;                  // looping index
    int K = n;

    while (K != 0 && i < 40) {   // starting from LSB
        S[i] = K % 16;
        K = (K - S[i]) / 16;
        i += 1;
    }
    if (i >= 40) printf("An integer too large to convert.\n");
    for (; i < 40; i++)
        S[i] = 0;
}
// detail-specify constructor
Elem::Elem(int l , int* seq)
{
    int i;                      // looping index

    for (i = 0; i < l; i++) 
        S[i] = seq[l - 1 - i];
    for (i = l; i < 40; i++) 
        S[i] = 0;
}
// arithmetics in prime field
Elem Elem::operator+(const Elem& another)
{
    int i;                      // looping index
    BigNum sum = BigNum(*this) + BigNum(another);

    if (sum.geqP()) sum = sum - p;
    return Elem(sum);
}
Elem Elem::operator-(const Elem& another)
{
    int i;                      // looping index
    BigNum temp;        // additive inverse of another

    // temp = p - another
    for (i = 0; i < 40; i++) {
        temp.S[i] = p.S[i] - another.S[i];
        if (temp.S[i] < 0) {
            temp.S[i] += 16;
            temp.S[i + 1] -= 1;
        }
    }
    temp.lenCheck();
    return Elem(BigNum(*this) + temp);
}
Elem Elem::operator*(const Elem& another)
{
    int i, j;                   // looping indices
    BigNum prod = (BigNum(*this) * BigNum(another)) % p;

    return Elem(prod);
}
Elem Elem::operator/(const Elem& another)
{
    int i;                  // index

    // find mul inverse of another
    i = EuclidAlg(BigNum(another));
    return ((*this) * Elem(vj[i] % p));
}
// length of the element in hex
int Elem::len(void) const
{
    int i;                      // looping index

    for (i = 39; i >= 0; i--)   // find largest nonzero digit
        if (this->S[i] != 0) 
            return (i + 1);
    return 0;
}
// modular exponentiation: compute (a^b) mod p
Elem Elem::modExp(const Elem& b)
{
    int* b_bit;             // b in binary representation
    int l = b.len();        // length of b in hex
    int k = 4 * l;          // length of b in binary
    int i, j;               // looping indices
    int temp;               // temp storage
    Elem d(1);              // d = 1

    b_bit = (int*)calloc(k, sizeof(int));
    // convert b to bit array
    for (i = 0, j = 0; i < l; i++, j += 4) {
        temp = b.S[i]; 
        b_bit[j] = temp % 2;
        temp = temp / 2;
        b_bit[j + 1] = temp % 2;
        temp = temp / 2;
        b_bit[j + 2] = temp % 2;
        temp = temp / 2;
        b_bit[j + 3] = temp % 2;
    }
    // modular exponentiation
    for (i = k - 1, j = 0; i >= 0; i--) {
        j = 2 * j;
        d = d * d;
        if (b_bit[i] == 1) {
            j += 1;
            d = d * (*this);
        }
    }
    free(b_bit);
    return d;
}
void Elem::print(void)        // print the number in hex
{      
    int i;                      // looping index
    int k;                      // counter

    for (i = 39, k = 0; i >= 0; i--, k++) {
        if (k == 8) {
            printf(" ");
            k = 0;
        }
        printf("%x", S[i]);     // starting from MSB
    }
    printf("\n");
}
// default constructor
point::point(void)
{          
    x = 0;      // point at infinity
    y = 0;
}
// construct given x and y coordinates
point::point(const Elem& X, const Elem& Y)
{
    x = X;
    y = Y;
}
point::point(const point& another)
{
    x = another.x;
    y = another.y;
}
point point::operator+(point& another)
{
    point sum;
    Elem lmb, X, Y;

    lmb = (another.y - this->y) / (another.x - this->x);
    X = lmb * lmb - this->x - another.x;
    Y = lmb * (this->x - X) - this->y;
    sum.x = X;
    sum.y = Y;
    return sum;
}
void point::print(void)
{
    this->x.print();
    this->y.print();
}
