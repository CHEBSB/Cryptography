// COM5335 Homework 3: Miller-Rabin Primality Test
// 111064502, 曾雋卿
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define L 150
class BigNum {
private:
    int sign;                   // its sign (+1 / -1)
    int length;                 // length of the number in hex
    int S[L];                   // sequence of number
public:
    BigNum();
    BigNum(const int&);         // directly from a integer
    BigNum(const BigNum&);      // copy another big number
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
    BigNum operator%(BigNum&);
    BigNum operator%(const int&);
    bool operator==(const BigNum&);      
    bool operator>(const BigNum&);
    bool abscmp(const BigNum&); // compare abs value (true for >=)
    bool IsEven();              // check if this is even
    bool leq1();                // check if this <= 1
    bool eq1();                 // check if this == 1
    bool eq0();                 // check if this == 0
    // modular exponentiation
    BigNum modExp(const BigNum& b, BigNum& n);
    BigNum Lshift(int p);       // left-shift p positions
    bool TrialDiv();            // Trial Division
    // assign an array to be S[]
    void assg(const int*, const int&);                
    void lenCheck();            // renew one's length
    void print();               // print the number in hex
};
// Miller-Rabin Primality Test
bool MillerRabin(BigNum&, const int&);

// 168 small prime numbers for trial division
const int smallP[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151,
157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211,
223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271,
277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347,
349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409,
419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467,
479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557,
563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617,
619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683,
691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761,
769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839,
853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919,
929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997};
BigNum SmallP[168];

int main(void)
{
    int i;                      // looping index
    int R[64];                  // array of 64 random numbers
    bool flag;                  // whether we got a prime
    BigNum r;                   // 256-bit number
    int j = 0;                  // for debug


    // convert prime number table to big number
    for (i = 0; i < 168; i++)
        SmallP[i] = smallP[i];
    flag = false;
    srand(time(NULL));          // set seed
    // loop until a prime is generated
    while (flag == false) {
        // generate 63 integers within [0, 15]
        for (i = 1; i < 64; i++)
            R[i] = rand() % 16;
        // least significant byte = odd
        R[0] = 2 * (rand() % 8) + 1;
        r.assg(R, 64);              // 256-bit number  
        // primality test
        if (!(r.leq1())) {
            if (r.TrialDiv()) {   
                if (MillerRabin(r, 8)) {
                    printf("r = ");
                    r.print();
                    flag = true;
                }
            }
        }
        j += 1;
        if (j % 50000 == 0) {
            printf("50000 iterations done.\n");
            j = 0;
        }
    }
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
// multiply the number by 16 to the power of p
BigNum BigNum::Lshift(int p)
{
    int i;                      // looping index
    BigNum newB;

    newB.length = this->length + p;
    newB.sign = this->sign;
    for (i = 0; i < p; i++) 
        newB.S[i] = 0;
    for (; i < this->length + p; i++)
        newB.S[i] = this->S[i - p];
    for (; i < L; i++)
        newB.S[i] = 0;
    return newB;
}
// compare absolute value (return true when >=)
bool BigNum::abscmp(const BigNum& another)
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
bool BigNum::operator>(const BigNum& another)
{
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
bool BigNum::operator==(const BigNum& another)
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
BigNum BigNum::operator%(BigNum& another)
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
BigNum BigNum::operator%(const int& b)
{
    BigNum B(b);
    return ((*this) % B);
}   
bool BigNum::IsEven()           // check if this is even
{          
    if ((this->S[0]) % 2 == 0) return true;
    return false;
}            
void BigNum::print(void)        // print the number in hexadecimal
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
// assign an array to be its value
void BigNum::assg(const int* A, const int& l) 
{
    int i;              // looping index
    
    this->length = 1;
    for (i = 0; i < l; i++) {
        this->S[i] = A[i];
        if (A[i] != 0) 
            this->length = i + 1;
    }
    for (i = l; i < L; i++) 
        this->S[i] = 0;
}
// check if this number <= 1
bool BigNum::leq1()
{
    if (this->length > 1) return false;
    if (this->S[0] < 2) return true;
    return false;
}
// check if this number == 1
bool BigNum::eq1()
{
    if (this->length > 1) return false;
    if (this->S[0] == 1) return true;
    return false;
}
// check if this number == 0
bool BigNum::eq0()
{
    int i;      // looping index

    for (i = 0; i < this->length; i++)
        if (this->S[0] != 0) return false;
    return true;
}
// test primality using test division
bool BigNum::TrialDiv()
{
    int i;              // looping index
    BigNum rem;         // remainder

    for (i = 0; i < 168; i++) {
        rem = (*this) % (SmallP[i]);
        if (rem.eq0())
            return false;
    }
    return true;
}
// Test if a number is prime with high probability
bool MillerRabin(BigNum& n, const int& T)
{
    int i, j;           // looping indices
    int k;              // n - 1 = m * (2^k)
    BigNum m;
    BigNum n_1;         // represent (n - 1)
    BigNum a;           // random number in iteration
    BigNum y;           // y = a^m mod n

    // compute k and m
    m = n - 1;
    n_1 = m;            // record value of (n - 1)
    k = 0;
    while (m.IsEven()) {
        m = m / 2;
        k += 1;
    }
    // repeat T times
    for (i = 0; i < T; i++) {
        // generate a
        a = 0;
        while (a.leq1()) {
            a = rand();
            a = a % n_1;
        }
        // compute y
        y = a.modExp(m, n);
        if (!(y.eq1()) && !(y == n_1)) {
            j = 1;
            while (j < k && !(y == n_1)) {
                y = (y * y) % n;
                if (y.eq1()) return false;  // composite
                j += 1;
            }
            if (!(y == n_1)) return false;
        }
    } 
    return true;            // high probability it is prime
}
