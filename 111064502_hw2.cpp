#include <stdio.h>
#include <stdlib.h>

// element of GF(256): polynomial over GF(2)
class poly {
public:
    int c;                          // the 8-bit integer
    poly();                         // init as 0 by default
    poly(int);                      // directly from a integer
    poly(const poly&);              // copy another element
    void print(int = 0);            // print out the value in hex
    int deg();                      // degree of the polynomial
    poly mul_x();                   // multiply  by x
    poly operator+(const poly&);    // addition (bitwise XOR)
    poly operator*(const poly&);    // multiplication mod 0x11B
    poly operator*(const int&);     // multipli with input as int
    // long division: this = Q(x) * another + R(x)
    poly longDiv(const poly&, poly*);
    int mulInv();                   // multiplicative inverse
    void bytesub();                 // perform byte substitution 
};


void round(poly** S);           // 1 round in AES except AddRoundKey
void finalRound(poly** S);              // the last round
void AddRoundKey(poly** S, poly** key); // add round key to state
void keyExpansion(poly** key);          // key Expansion
void printState(poly** S);              // print current text
void AES_Encrypt(poly** S, poly** key); // Encryption of AES

// for Euclid's algorithm
poly uj[3];                 
poly rj[3];
poly vj[3];
poly qj[3];
// round constant
poly RC(0x01);

int main(void)
{
    int i, j, k;        // looping indices
    int temp;           // temp storage for reading input
    poly** state;       // 16-byte text
    poly** key;         // 16-byte round key
    char plain[36];     // plaintext as string
    char KEY[36];       // initial key

    // allocate memory
    state = (poly**)calloc(4, sizeof(poly*));
    key = (poly**)calloc(4, sizeof(poly*));
    for (i = 0; i < 4; i++) {
        state[i] = (poly*)calloc(4, sizeof(poly));
        key[i] = (poly*)calloc(4, sizeof(poly));
    } 
    // read input
    printf("Plaintext: ");
    scanf("%s", plain);
    printf("Key: ");
    scanf("%s", KEY);
    // convert char to integer
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            k = 8 * i + 2 * j;  // the higher byte
            if (('0' <= plain[k]) && (plain[k] <= '9')) 
                temp = plain[k] - '0';
            else if (('a' <= plain[k]) && (plain[k] <= 'f'))  
                temp = plain[k] - 'a' + 0x0a;

            else if (('A' <= plain[k]) && (plain[k] <= 'F')) 
                temp = plain[k] - 'A' + 0x0a;
            else 
                printf("invalid: %c\n", plain[k]);
            state[j][i].c = temp << 4;
            // the lower byte
            if (('0' <= plain[k + 1]) && (plain[k + 1] <= '9')) 
                temp = plain[k + 1] - '0';
            else if (('a' <= plain[k + 1]) && (plain[k + 1] <= 'f')) 
                temp = plain[k + 1] - 'a' + 0x0a;
            else if (('A' <= plain[k + 1]) && (plain[k + 1] <= 'F')) 
                temp = plain[k + 1] - 'A' + 0x0a; 
            else 
                printf("invalid: %c\n", plain[k + 1]);
            state[j][i].c += temp;
        }
    }
    // same procedure for key
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            k = 8 * i + 2 * j;  // the higher byte
            if (('0' <= KEY[k]) && (KEY[k] <= '9')) 
                temp = KEY[k] - '0';
            else if (('a' <= KEY[k]) && (KEY[k] <= 'f')) 
                temp = KEY[k] - 'a' + 0x0a;
            else if (('A' <= KEY[k]) && (KEY[k] <= 'F')) 
                temp = KEY[k] - 'A' + 0x0a;
            else 
                printf("invalid: %c\n", KEY[k]);
            key[j][i].c = temp << 4;
            // the lower byte
            if (('0' <= KEY[k + 1]) && (KEY[k + 1] <= '9')) 
                temp = KEY[k + 1] - '0';
            else if (('a' <= KEY[k + 1]) && (KEY[k + 1] <= 'f'))  
                temp = KEY[k + 1] - 'a' + 0xa;
            else if (('A' <= KEY[k + 1]) && (KEY[k + 1] <= 'F')) 
                temp = KEY[k + 1] - 'A' + 0xa;
            else 
                printf("invalid: %c\n", KEY[k + 1]);
            key[j][i].c += temp;
        }
    }
    // perform AES encryption and print out the result
    AES_Encrypt(state, key);
    return 0;
}

// default constructor: init as 0
poly::poly(void) {
    c = 0;
}
// constructor from an integer in [0, 255]
poly::poly(int n) {
    if (n < 0 || n > 255)
        printf("integer out of range!\n");
    c = n;
}
// copy constructor
poly::poly(const poly& another) {
    c = another.c;
}
// print out the value in hexadecimal
void poly::print(int flag) {
    if (this->c < 0x10) printf("0");    // if the left byte is zero
    printf("%x", this->c);
    if (flag == 0) printf("\n");
    else if (flag == 2) printf(" ");
    return;
}
// return the degree of the polynomial
int poly::deg(void) {
    if (this->c >= 256) // for 0x11B
        return 8;
    if (this->c >= 128)
        return 7;
    if (this->c >= 64)
        return 6;
    if (this->c >= 32)
        return 5;
    if (this->c >= 16)
        return 4;
    if (this->c >= 8)
        return 3;
    if (this->c >= 4)
        return 2;
    if (this->c >= 2)
        return 1;
    return 0;
}
// addition in GF(256): bitwise XOR
poly poly::operator+(const poly& another) {
    poly sum;       // the sum to be returned

    sum.c = (this->c) ^ (another.c);
    return sum;
}
// multiply the element by x in GF(256)
poly poly::mul_x(void) {
    poly result;        // the result to be returned

    if (this->c >= 128) {       // highest bit = 1
        result.c = this->c - 128;
        result.c = result.c << 1;
        result.c = result.c ^ 0x1B;     // modulo
    }
    else                        // highest bit = 0
        result.c = this->c << 1;
    return result;
}
// multiplication in GF(256)
poly poly::operator*(const poly& another) {
    int i, j;           // looping indices
    poly prod;          // the prod to be returned
    int y;              // a copy of another that can be modified
    poly temP;          // temp storage of a polynomial

    // do multiplication using addition
    y = another.c;
    for (i = 0; i < 8; i++) {
        if (y % 2 == 1) {
            temP.c = this->c;           // temP = this
            for (j = 0; j < i; j++)     // to the power of i
                temP = temP.mul_x();
            prod = prod + temP;
        }
        y = y >> 1;                     // go to the next bit
    }
    return prod;
}
// multiplication with input = int
poly poly::operator*(const int& a) {
    int i, j;           // looping indices
    poly prod;          // the prod to be returned
    int y;              // a copy of another that can be modified
    poly temP;          // temp storage of a polynomial

    // do multiplication using addition
    y = a;
    for (i = 0; i < 8; i++) {
        if (y % 2 == 1) {
            temP.c = this->c;           // temP = this
            for (j = 0; j < i; j++)     // to the power of i
                temP = temP.mul_x();
            prod = prod + temP;
        }
        y = y >> 1;                     // go to the next bit
    }
    return prod;
}
// long division of polynomial over GF(2)
// return the remainder while modifying the quotient
poly poly::longDiv(const poly& another, poly* Q) {
    int i, j;           // degree of remainder and divisor
    int q = 0;          // quotient as integer in Z(256)
    poly R(*this);      // the remainder, init as dividand
    poly D(another);    // a copy of divisor

    if (another.c == 0) {
        printf("Division by 0!\n");
        return D;
    }
    R.c = this->c;
    i = R.deg();
    j = D.deg();
    while (i >= j) {
        q += (1 << (i - j));
        R.c = R.c ^ (D.c << (i - j));
        i = R.deg();
    }
    Q->c = q;           // copy q to the value of Q
    return R;           // return the remainder
}
/* 
Euclid's Algorithm with 11B = a(x) > b(x) = this polynomial.
stop when remainder = 1.
multiplicative inverse of this = vj[i + 1] mod (0x11B).
*/
int poly::mulInv(void)
{
    int i;          // looping index
    poly inv;       // the inverse to be returned

    if (this->c == 0) return 0;     // special case: 0 -> 0
    // u[-1] = 1    v[-1] = 0   r[-1] = a(x)
    // u[0] = 0     v[0] = 1    r[0] = b(x)
    uj[0].c = 1;
    vj[0].c = 0;
    rj[0].c = 0x11B;
    uj[1].c = 0;
    vj[1].c = 1;
    rj[1].c = this->c;
    i = 0;
    while (rj[(i + 1) % 3].c != 1) {
        // find qj(x)
        rj[(i + 2) % 3] =
        rj[i].longDiv(rj[(i + 1) % 3], &(qj[(i + 2) % 3]));
        // update uj(x) and vj(x)
        uj[(i + 2) % 3] = uj[i] + qj[(i + 2) % 3] * uj[(i + 1) % 3]; 
        vj[(i + 2) % 3] = vj[i] + qj[(i + 2) % 3] * vj[(i + 1) % 3];
        // index increment
        i = (i + 1) % 3;
    } 
    i = (i + 1) % 3;
    // if vj[i] has degree too large
    if (vj[i].c > 255) {
        uj[0].c = 0x11B;        // modulo 0x11B
        vj[i] = vj[i].longDiv(uj[0], &(qj[0]));
    }
    return vj[i].c;
}
// byte substitution of 1 byte
void poly::bytesub(void) {
    int y;          // the result as integer
    int inv;        // multipli inverse of a

    // multiplication and addition over GF(2)
    y = 0x63;       // init as 01100011
    inv = this->mulInv();
    if ((inv & 1) == 1)
        y = y ^ 0x1f;
    if ((inv & 2) == 2)
        y = y ^ 0x3e;
    if ((inv & 4) == 4)
        y = y ^ 0x7c;
    if ((inv & 8) == 8)
        y = y ^ 0xf8;
    if ((inv & 16) == 16)
        y = y ^ 0xf1;
    if ((inv & 32) == 32)
        y = y ^ 0xe3;
    if ((inv & 64) == 64)
        y = y ^ 0xc7;
    if ((inv & 128) == 128)
        y = y ^ 0x8f;
    this->c = y;
    return;
}
// Add round key to current state
void AddRoundKey(poly** S, poly** key) {
    int i;          // looping index

    for (i = 0; i < 4; i++) {
        S[i][0] = S[i][0] + key[i][0];
        S[i][1] = S[i][1] + key[i][1];
        S[i][2] = S[i][2] + key[i][2];
        S[i][3] = S[i][3] + key[i][3];
    }
}
// each round (except the last one) in AES
// S: current state, key: current round key
void round(poly** S) {
    int i;              // looping index
    int temp[4];        // temp storage for shifting
    poly temP[4][4];    // temp storage during mix column

    // ByteSub
    for (i = 0; i < 4; i++) {
        S[i][0].bytesub();
        S[i][1].bytesub();
        S[i][2].bytesub();
        S[i][3].bytesub();
    }
    // ShiftRow
    // row 1: cyclic left shift by 1 byte
    for (i = 0; i < 4; i++)
        temp[i] = S[1][i].c;
    S[1][0].c = temp[1];
    S[1][1].c = temp[2];
    S[1][2].c = temp[3];
    S[1][3].c = temp[0];
    // row 2: cyclic left shift by 2 bytes
    for (i = 0; i < 4; i++)
        temp[i] = S[2][i].c;
    S[2][0].c = temp[2];
    S[2][1].c = temp[3];
    S[2][2].c = temp[0];
    S[2][3].c = temp[1];
    // row 3: cyclic left shift by 3 bytes
    for (i = 0; i < 4; i++)
        temp[i] = S[3][i].c;
    S[3][0].c = temp[3];
    S[3][1].c = temp[0];
    S[3][2].c = temp[1];
    S[3][3].c = temp[2];
    // MixColumn
    // matrix multiplication over GF(256)
    for (i = 0; i < 4; i++) {
        temP[0][i] =
    S[0][i].mul_x() + S[1][i].mul_x() + S[1][i] + S[2][i] + S[3][i];
        temP[1][i] = 
    S[0][i] + S[1][i].mul_x() + S[2][i].mul_x() + S[2][i] + S[3][i];
        temP[2][i] = 
    S[0][i] + S[1][i] + S[2][i].mul_x() + S[3][i].mul_x() + S[3][i];
        temP[3][i] = 
    S[0][i].mul_x() + S[0][i] + S[1][i] + S[2][i] + S[3][i].mul_x();
    }
    // copy temP back to S
    for (i = 0; i < 4; i++) {
        S[i][0].c = temP[i][0].c;
        S[i][1].c = temP[i][1].c;
        S[i][2].c = temP[i][2].c;
        S[i][3].c = temP[i][3].c;
    }
    return;
}
// the final round of AES encryption
void finalRound(poly** S) {
    int i;              // looping index
    int temp[4];        // temp storage for shifting

    // ByteSub
    for (i = 0; i < 4; i++) {
        S[i][0].bytesub();
        S[i][1].bytesub();
        S[i][2].bytesub();
        S[i][3].bytesub();
    }
    // ShiftRow
    // row 1: cyclic left shift by 1 byte
    for (i = 0; i < 4; i++)
        temp[i] = S[1][i].c;
    S[1][0].c = temp[1];
    S[1][1].c = temp[2];
    S[1][2].c = temp[3];
    S[1][3].c = temp[0];
    // row 2: cyclic left shift by 2 bytes
    for (i = 0; i < 4; i++)
        temp[i] = S[2][i].c;
    S[2][0].c = temp[2];
    S[2][1].c = temp[3];
    S[2][2].c = temp[0];
    S[2][3].c = temp[1];
    // row 3: cyclic left shift by 3 bytes
    for (i = 0; i < 4; i++)
        temp[i] = S[3][i].c;
    S[3][0].c = temp[3];
    S[3][1].c = temp[0];
    S[3][2].c = temp[1];
    S[3][3].c = temp[2];
    return;
}
// key Expansion after each round
void keyExpansion(poly** key) {
    poly temP[4];       // temp storage of a word

    // RotWord
    // copy the rightmost word but shifted
    temP[0].c = key[1][3].c;
    temP[1].c = key[2][3].c;
    temP[2].c = key[3][3].c;
    temP[3].c = key[0][3].c;
    // SubWord
    temP[0].bytesub();
    temP[1].bytesub();
    temP[2].bytesub();
    temP[3].bytesub();
    // XOR with Rcon
    temP[0] = temP[0] + RC;
    // add temP to the 0th word
    key[0][0] = key[0][0] + temP[0];
    key[1][0] = key[1][0] + temP[1];
    key[2][0] = key[2][0] + temP[2];
    key[3][0] = key[3][0] + temP[3];
    // add the 0th word to the 1st word
    key[0][1] = key[0][1] + key[0][0];
    key[1][1] = key[1][1] + key[1][0];
    key[2][1] = key[2][1] + key[2][0];
    key[3][1] = key[3][1] + key[3][0];
    // add the 1st word to the 2nd word
    key[0][2] = key[0][2] + key[0][1];
    key[1][2] = key[1][2] + key[1][1];
    key[2][2] = key[2][2] + key[2][1];
    key[3][2] = key[3][2] + key[3][1];
    // add the 2nd word to the 3rd word
    key[0][3] = key[0][3] + key[0][2];
    key[1][3] = key[1][3] + key[1][2];
    key[2][3] = key[2][3] + key[2][2];
    key[3][3] = key[3][3] + key[3][2];
    // double Rc
    RC = RC.mul_x();
    return;
}
// print current text
void printState(poly** S) {
    S[0][0].print(2);
    S[1][0].print(2);
    S[2][0].print(2);
    S[3][0].print(2);
    S[0][1].print(2);
    S[1][1].print(2);
    S[2][1].print(2);
    S[3][1].print(2);
    S[0][2].print(2);
    S[1][2].print(2);
    S[2][2].print(2);
    S[3][2].print(2);
    S[0][3].print(2);
    S[1][3].print(2);
    S[2][3].print(2);
    S[3][3].print(1);
    return;
}  
// Encryption of AES
void AES_Encrypt(poly** S, poly** key) {
    int i;          // looping index

    // add round key at the beginning
    AddRoundKey(S, key);
    // the first 9 rounds
    for (i = 0; i < 9; i++) {
        keyExpansion(key);
        round(S);
        AddRoundKey(S, key);
        printf("Round %d:\n", i + 1);
        printState(S);
        printf("\n");
    }
    // the last round
    keyExpansion(key);
    finalRound(S);
    AddRoundKey(S, key);
    printf("Round 10:\n");
    printState(S);
    return;
}