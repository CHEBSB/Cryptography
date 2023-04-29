// COM5335 Homework 3: Rabin Decryption
// 111064502, 曾雋卿
#include <stdio.h>
#include <stdlib.h>
#define L 200

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
    bool operator==(const BigNum&);      
    bool operator>(const BigNum&);
    bool abscmp(const BigNum&); // compare abs value (true for >=)
    bool Is3mod4();             // check if this = 3 mod 4
    bool LegalPlain();          // check for padding
    // modular exponentiation
    BigNum modExp(const BigNum& b, BigNum& n);
    BigNum Lshift(int p);       // left-shift p positions
    void lenCheck();            // renew one's length
    void print();               // print the number in hex
    void printPlain();          // print w/out padded bits
};

int pow(int a, int x);              // a to the power of x
int charToInt(char*, int, int*);    // convert string to integer array
// find quadratic residue of a mod p
BigNum Sqrt3mod4(BigNum& p, BigNum& a); // p = 3 mod 4
BigNum Sqrt5mod8(BigNum& p, BigNum& a); // p = 5 mod 8
// extended Euclid algorithm on 2 primes p and q
int EuclidAlg(const BigNum&, const BigNum&);
// for Euclid's Algorithm
BigNum uj[3];                 
BigNum rj[3];
BigNum vj[3];
BigNum qj[3];

int main(void)
{
    char p_s[36];               // input prime numbers 
    char q_s[36];
    char ciphtxt[72];           // input plain text
    int ciph[72];               // plain text in hex 
    int P[32];                  // prime numbers in hex
    int Q[32];
    int l;                      // len of hex array
    int i;                      // looping index
    char temp;                  // temp storage
    BigNum r;                   // %p quadratic residue 
    BigNum s;                   // %q quadratic residue 
    BigNum N;                   // N = p * q
    BigNum Temp1, Temp2;        // Temp1 = rdq, Temp2 = scp
    BigNum x1, x2, y1, y2;      // 4 sqroot of ciphertext

    // read 64-byte cipher text
    printf("Ciphertext = ");
    i = 0;
    while (i < 64) {
        scanf("%c", &temp);
        if ((('0' <= temp) && (temp <= '9'))
        || (('a' <= temp) && (temp <= 'f'))
        || (('A' <= temp) && (temp <= 'F'))) {
            ciphtxt[i] = temp;
            i += 1;
        }
    }
    ciphtxt[i] = '\0';
    // convert char to integer
    l = charToInt(ciphtxt, 64, ciph);
    // read p
    printf("Private Key:\n");
    printf("p = ");
    i = 0;
    while (i < 32) {
        scanf("%c", &temp);
        if ((('0' <= temp) && (temp <= '9'))
        || (('a' <= temp) && (temp <= 'f'))
        || (('A' <= temp) && (temp <= 'F'))) {
            p_s[i] = temp;
            i += 1;
        }
    }
    p_s[i] = '\0';
    charToInt(p_s, 36, P);
    // read q
    printf("q = ");
    i = 0;
    while (i < 32) {
        scanf("%c", &temp);
        if ((('0' <= temp) && (temp <= '9'))
        || (('a' <= temp) && (temp <= 'f'))
        || (('A' <= temp) && (temp <= 'F'))) {
            q_s[i] = temp;
            i += 1;
        }
    }
    q_s[i] = '\0';
    charToInt(q_s, 36, Q);
    // Rabin Decryption
    // convert to big number
    BigNum C(1, 64, ciph);
    BigNum p(1, 32, P);
    BigNum q(1, 32, Q);
    // compute N = p * q
    N = p * q;
    // find quadratic residue of C mod p
    Temp1 = C % p;          //  C mod p before finding sqrt
    if (p.Is3mod4()) 
        r = Sqrt3mod4(p, Temp1);
    else
        r = Sqrt5mod8(p, Temp1);
    // find quadratic residue of C mod q
    Temp2 = C % q;          //  C mod q before finding sqrt
    if (q.Is3mod4()) 
        s = Sqrt3mod4(q, Temp2);
    else
        s = Sqrt5mod8(q, Temp2);
    // perform extended Euclid Algorithm
    i = EuclidAlg(p, q);
    // c = uj[i], d = vj[1]
    Temp1 = r * (vj[i] * q);
    Temp2 = s * (uj[i] * p);
    // test the 4 square roots one by one
    x1 = (Temp1 + Temp2) % N;
    x2 = N - x1;
    y1 = (Temp1 - Temp2) % N;
    y2 = N - y1;
    printf("Plaintext = ");
    if (x1.LegalPlain()) 
        x1.printPlain();
    else if (x2.LegalPlain())
        x2.printPlain();
    else if (y1.LegalPlain())
        y1.printPlain();
    else if (y2.LegalPlain())
        y2.printPlain();
    else 
        printf("No legal plaintext available.");
    return 0;
}

// compute (a to the power of x). integers only
int pow(int a, int x) {
    int y = 1;                  // output

    while (x > 0) {
        y = y * a;
        x -= 1;
    }
    return y;
}       
// convert a string of length l to hex integer array
int charToInt(char* s, int l, int* A) {
    int i, j;           // looping indices
    int k = 0;          // counter of digit
    int temp;           // temp storage
 
    for (i = 0; i < l; i++) {
        if (('0' <= s[i]) && (s[i] <= '9')) {
            temp = s[i] - '0';
            A[k] = temp;
            k += 1;
        }
        else if (('a' <= s[i]) && (s[i] <= 'f')) {
            temp = s[i] - 'a' + 0x0a;
            A[k] = temp;
            k += 1;
        }
        else if (('A' <= s[i]) && (s[i] <= 'F')) {
            temp = s[i] - 'A' + 0x0a;
            A[k] = temp;
            k += 1;
        }
        else if (s[i] == '\0')  // reach the end
            return k;
        // else, skip whitespace
    }
    return k;
}
// default constructor: init as 0
BigNum::BigNum(void) {          
    int i;                      // looping index

    sign = 1;
    length = 1;
    for (i = 0; i < L; i++) 
        S[i] = 0;
}
// copy constructor
BigNum::BigNum(const BigNum& another) {
    int i;                      // looping index

    sign = another.sign;
    length = another.length;
    for (i = 0; i < L; i++) 
        S[i] = another.S[i];
}
// convert an integer into BigNum
BigNum::BigNum(const int& n) {  
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
BigNum::BigNum(int sgn, int l , int* seq) {
    int i;                      // looping index

    sign = sgn;
    length = l;
    for (i = 0; i < l; i++) 
        S[i] = seq[l - 1 - i];
    for (i = l; i < L; i++) 
        S[i] = 0;
}
void BigNum::lenCheck(void) {
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
    return d;
}
// multiply the number by 16 to the power of p
BigNum BigNum::Lshift(int p) {
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
bool BigNum::abscmp(const BigNum& another) {
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
bool BigNum::operator>(const BigNum& another) {
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
bool BigNum::operator==(const BigNum& another) {
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
BigNum BigNum::operator+(const BigNum& another) {
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
BigNum BigNum::operator+(const int& b) {
    BigNum B(b);
    return ((*this) + B);
}
BigNum BigNum::operator-(const BigNum& another) {
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
BigNum BigNum::operator-(const int& b) {
    BigNum B(b);
    return ((*this) - B);
}
BigNum BigNum::operator*(const BigNum& another) {
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
BigNum BigNum::operator*(const int& b) {
    BigNum B(b);
    return ((*this) * B);
}
BigNum BigNum::operator/(BigNum& another) {
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
BigNum BigNum::operator/(const int& b) {
    BigNum B(b);
    return ((*this) / B);
}
BigNum BigNum::operator%(BigNum& another) {
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
// check the last 16 bits to see if this is plaintext
bool BigNum::LegalPlain() {
    if (this->S[0] != this->S[4])
        return false;
    if (this->S[1] != this->S[5])
        return false;
    if (this->S[2] != this->S[6])
        return false;
    if (this->S[3] != this->S[7])
        return false;
    return true;
}        
bool BigNum::Is3mod4() {        // check if this = 3 mod 4
    if ((this->S[0]) % 4 == 3) return true;
    return false;
}            
void BigNum::print(void) {      // print the number in hexadecimal
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
void BigNum::printPlain() {     // omit 16 padded bits
    int i;                      // looping index
    int k;                      // counter

    if (sign == -1) printf("-");
    for (i = length - 1, k = 0; i >= 4; i--, k++) {
        if (k == 8) {
            printf(" ");
            k = 0;
        }
        printf("%x", S[i]);     // starting from MSB
    }
    printf("\n");
}
// find 2 sqrt of a mod p, where p = 3 mod 4
BigNum Sqrt3mod4(BigNum& p, BigNum& a) {
    BigNum k;       // k = (p + 1) / 4
    BigNum r;       // one sqare root

    k = (p + 1) / 4;
    r = a.modExp(k, p);
    return r;
}
// find 2 sqrt of a mod p, where p = 5 mod 8
BigNum Sqrt5mod8(BigNum& p, BigNum& a) {
    BigNum k, temp; // temp storage
    BigNum d;       // the first check
    BigNum r;       // one sqare root

    k = (p - 1) / 4;
    d = a.modExp(k, p);
    if (d == p - 1) {
        k = (p - 5) / 8;
        temp = a * 4;
        r = temp.modExp(k, p);
        r = (r * (a * 2)) % p;
    } else {
        k = (p + 3) / 8;
        r = a.modExp(k, p);
    }
    return r;
}
// Extended Euclid's Algorithm stop at rj = 1
int EuclidAlg(const BigNum& p, const BigNum& q)
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
