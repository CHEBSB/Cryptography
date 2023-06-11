#include <stdio.h>
#include <stdlib.h>
#define L 42

class BigNum {                  // big finite field element
    friend class point;
private:
    int length;                 // length of the number in hex
    int S[L];                   // sequence of number
public:
    BigNum();
    BigNum(const int&);         // directly from a integer
    BigNum(const BigNum&);      // copy another big number
    BigNum(int, int*);
    // arithmetic: assume the inputs are positive
    BigNum operator+(const BigNum&);
    BigNum operator-(const BigNum&);
    BigNum operator*(const BigNum&);
    BigNum operator/(BigNum&);  
    BigNum Lshift(int p);       // left-shift p positions
    void lenCheck();            // renew one's length
    void print();               // print out the number in hex
};
class point {                   // point on elliptic curve
private:
    BigNum x;                  // x coordinate
    BigNum y;                  // y coordinate
public:
    point();   
    point(const BigNum&, const BigNum&);    // given x, y
    point(const point&);        // copy constructor
};
// elliptic curve parameters
const int P[] = {0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF,
0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF,
0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF,
0xF, 0xF, 0x7, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF};
const int A[] = {0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF,
0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF,
0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF,
0xF, 0xF, 0x7, 0xF, 0xF, 0xF, 0xF, 0xF, 0xF, 0xC};
const int B[] = {0x1, 0xC, 0x9, 0x7, 0xB, 0xE, 0xF, 0xC,
0x5, 0x4, 0xB, 0xD, 0x7, 0xA, 0x8, 0xB, 0x6, 0x5, 0xA,
0xC, 0xF, 0x8, 0x9, 0xF, 0x8, 0x1, 0xD, 0x4, 0xD, 0x4,
0xA, 0xD, 0xC, 0x5, 0x6, 0x5, 0xF, 0xA, 0x4, 0x5};
const int Gx[] = {0x4, 0xA, 0x9, 0x6, 0xB, 0x5, 0x6, 0x8,
0x8, 0xE, 0xF, 0x5, 0x7, 0x3, 0x2, 0x8, 0x4, 0x6, 0x6,
0x4, 0x6, 0x9, 0x8, 0x9, 0x6, 0x8, 0xC, 0x3, 0x8, 0xB,
0xB, 0x9, 0x1, 0x3, 0xC, 0xB, 0xF, 0xC, 0x8, 0x2};
const int Gy[] = {0x2, 0x3, 0xA, 0x6, 0x2, 0x8, 0x5, 0x5,
0x3, 0x1, 0x6, 0x8, 0x9, 0x4, 0x7, 0xD, 0x5, 0x9, 0xD,
0xC, 0xC, 0x9, 0x1, 0x2, 0x0, 0x4, 0x2, 0x3, 0x5, 0x1,
0x3, 0x7, 0x7, 0xA, 0xC, 0x5, 0xF, 0xB, 0x3, 0x2};
const int N[] = {0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0xF, 4, 0xC, 8, 0xF, 9,
2, 7, 0xA, 0xE, 0xD, 3, 0xC, 0xA, 7, 5, 2, 2, 5, 7};

int main(void)
{
    point Test;

    return 0;
}

// default constructor: init as 0
BigNum::BigNum(void)
{          
    int i;                      // looping index

    length = 1;
    for (i = 0; i < L; i++) 
        S[i] = 0;
}
// copy constructor
BigNum::BigNum(const BigNum& another)
{
    int i;                      // looping index

    length = another.length;
    for (i = 0; i < L; i++) 
        S[i] = another.S[i];
}
// convert an integer into BigNum
BigNum::BigNum(const int& n)
{  
    int i = 0;                  // looping index
    int K;

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
BigNum::BigNum(int l , int* seq)
{
    int i;                      // looping index

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
// multiply the number by 16 to the power of p
BigNum BigNum::Lshift(int p)
{
    int i;                      // looping index
    BigNum newB;

    newB.length = this->length + p;
    for (i = 0; i < p; i++) 
        newB.S[i] = 0;
    for (; i < this->length + p; i++)
        newB.S[i] = this->S[i - p];
    for (; i < L; i++)
        newB.S[i] = 0;
    return newB;
}
BigNum BigNum::operator+(const BigNum& another)
{
    int i;                      // looping index
    int maxL;                   // larger length of the 2
    BigNum sum;                 // sum of the 2 numebrs

    maxL = (this->length > another.length)? this->length: another.length;
    sum.length = maxL;
    for (i = 0; i < maxL; i++) {
        sum.S[i] += (this->S[i] + another.S[i]);
        if (sum.S[i] > 15) {
            sum.S[i] -= 16;
            sum.S[i + 1] += 1;
        }
    }
    if (sum.S[maxL] != 0) sum.length += 1;
    return sum;
}
BigNum BigNum::operator-(const BigNum& another)
{
    int i;                      // looping index
    BigNum result;              // = (this - another)

    // return [this + (-another]
    result.lenCheck();
    return result;
}
BigNum BigNum::operator*(const BigNum& another)
{
    int i, j;                   // looping indices
    BigNum prod;                // product of the 2

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
BigNum BigNum::operator/(BigNum& another)
{
    int i, j;                   // looping index
    BigNum r(*this);            // remainder
    BigNum temp(another);       // temporary variable
    BigNum q;                   // quotient

    // return [this * (another^-1)]
    return q;
}
void BigNum::print(void)        // print the number in hex
{      
    int i;                      // looping index
    int k;                      // counter

    for (i = length - 1, k = 0; i >= 0; i--, k++) {
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
    x = 0;      // note that this is not a point on curve
    y = 0;
}
// construct given x and y coordinates
point::point(const BigNum& X, const BigNum& Y)
{
    x = X;
    y = Y;
}
point::point(const point& another)
{
    x = another.x;
    y = another.y;
}
