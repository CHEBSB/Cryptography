// COM5335 Homework 3: Rabin Encryption
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
    BigNum operator-(const BigNum&);
    BigNum operator*(const BigNum&);
    BigNum operator/(BigNum&);
    BigNum operator%(BigNum&);      
    bool operator>(const BigNum&);
    bool abscmp(const BigNum&); // compare abs value (true for >=)
    BigNum Lshift(int p);       // left-shift p positions
    void lenCheck();            // renew one's length
    void print();               // print out the number in hex
};

int pow(int a, int x);              // a to the power of x
int charToInt(char*, int, int*);    // convert string to integer array

int main(void)
{
    char p_s[36];               // input prime numbers 
    char q_s[36];
    char plaintxt[64];          // input plain text
    int plain[64];              // plain text in hex 
    int P[32];                  // prime numbers in hex
    int Q[32];
    BigNum ciph, N;             // ciphertext; N = p * q
    int l;                      // len of hex array
    int i;                      // looping index
    char temp;                  // temp storage
    
    // read p
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
    // convert to big number
    BigNum p(1, 32, P);
    BigNum q(1, 32, Q);
    // compute and print N = p * q
    N = p * q;
    printf("n = pq = ");
    N.print();
    // read 56-byte plain text
    printf("Plaintext: ");
    i = 0;
    while (i < 56) {
        scanf("%c", &temp);
        if ((('0' <= temp) && (temp <= '9'))
        || (('a' <= temp) && (temp <= 'f'))
        || (('A' <= temp) && (temp <= 'F'))) {
            plaintxt[i] = temp;
            i += 1;
        }
    }
    plaintxt[i] = '\0';
    // convert char to integer
    l = charToInt(plaintxt, 64, plain);
    // repetitively pad the last 16 bits
    plain[l] = plain[l - 4];
    plain[l + 1] = plain[l - 3];
    plain[l + 2] = plain[l - 2];
    plain[l + 3] = plain[l - 1];
    l += 4;
    // Rabin Encryption = Big Number Arithmetic
    BigNum Plain(1, l, plain);
    ciph = (Plain * Plain) % N;
    printf("Ciphertext = ");
    ciph.print();

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
BigNum BigNum::operator+(const BigNum& another) {
    int i;                      // looping index
    int maxL;                   // larger length of the 2
    BigNum sum;                 // sum of the 2 numebrs

    maxL = (this->length > another.length)? this->length: another.length;
    sum.sign = 1;
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
BigNum BigNum::operator-(const BigNum& another) {
    int i;                      // looping index
    BigNum result;              // = (this - another)

    if (this->abscmp(another)) {            // result is nonnegative
        result.sign = 1;
        for (i = 0; i < this->length; i++) {
            result.S[i] += (this->S[i] - another.S[i]);
            if (result.S[i] < 0) {
                result.S[i] += 16;
                result.S[i + 1] -= 1;
            }
        }
    }
    else {                      // negative result => -(B - A)
        result.sign = -1;           
        for (i = 0; i < another.length; i++) {
            result.S[i] += (another.S[i] - this->S[i]);
            if (result.S[i] < 0) {
                result.S[i] += 16;
                result.S[i + 1] -= 1;
            }
        }
    }
    result.lenCheck();
    return result;
}
BigNum BigNum::operator*(const BigNum& another) {
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
BigNum BigNum::operator%(BigNum& another) {
    int i, j;                   // looping index
    BigNum r(*this);            // remainder
    BigNum temp(another);       // temporary variable

    i = -1;
    for (i = -1; this->abscmp(temp); i++)
        temp = temp.Lshift(1);  // go to next digit
    if (i == -1)                // if this < another         
        return r;               // return the original number
    // if i = k, then this >= (temp left shifted by k)
    for (j = i; j >= 0; j--) {
        temp = another.Lshift(j);
        while (r.abscmp(temp)) 
            r = r - temp;
    }
    r.lenCheck();
    return r;
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
