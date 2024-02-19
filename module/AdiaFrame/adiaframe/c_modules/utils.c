
#define XZ_COEF_POW(x, z) bit_count(x&z)
#define ENCODE_XZCODE(x, z)  (insert_zeros(x) << 1) + insert_zeros(z)
#include <stdio.h>

unsigned int  bit_count(unsigned int);
unsigned int  insert_zeros(unsigned int);
unsigned int xz_coef_pow(unsigned int, unsigned int);
void decode_pcode(unsigned int, size_t, unsigned int *, unsigned int * );


unsigned int bit_count(unsigned int n){
    //Brian Kernighan’s Algorithm
    int num = 0;
    while (n){n = n&(n-1);num++;}
    return num;
}

unsigned int xz_coef_pow(unsigned int x, unsigned int z){
    return XZ_COEF_POW(x, z);
}

unsigned int insert_zeros(unsigned int n){ 
    // Insert zeros to nested gaps 
    // 110 -> 010100 = (01)(01)(00)
    int result = 0;
    int bit_position = 0;
    int r_most_bit = 0;
    while (n>0){
        r_most_bit = n&1; //  (...101)&(...001) = (...001), 1 or 0
        result = result | (r_most_bit << (bit_position <<1));
        n = n>>1;
        bit_position++;
    }
    return result;
}
unsigned int encode_xzcode(unsigned int x, unsigned int z){
    return ENCODE_XZCODE(x, z);
}

void decode_pcode(
    unsigned int p_int, 
    size_t p_len, 
    unsigned int * x, 
    unsigned int * z){
    *x = 0;
    *z = 0;

    unsigned int mask =1;
    int bit_position =0;
    while (p_int > 0) {
        if (bit_position % 2 == 0) {
            // Even bit position: Add to x
            *x |= (p_int & 1) << (bit_position/ 2);
        } else {
            // Odd bit position: Add to z
            *z |= (p_int & 1) << (bit_position/ 2);
        }
        p_int >>= 1; // Shift right to get the next bit
        bit_position++;
    }
}