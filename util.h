#ifndef _JPEG_HEADER
#define _JPEG_HEADER
#include <stdint.h>
#include <inttypes.h>

typedef uint8_t byte;

typedef struct {
    //big endian
    byte higher_byte;
    byte lowwer_byte;
}word_unit;

typedef uint16_t word;

uint16_t read_word_to_bigendian(FILE* fp) {
    word c;
    fscanf(fp,"%"SCNd16,&c);
    // byte high = c.higher_byte;
    byte high = c >> 8;
    // byte low = c.lowwer_byte;
    byte low = c & 0x00ff;
    high = high << 8;
    return high|low;
}

#define START 0xFF
#define SOI 0xD8
#define APP0 0xE0
#define APP1 0xE1
#define APP2 0xE2
#define APP3 0xE3
#define APP4 0xE4
#define APP5 0xE5
#define APP6 0xE6
#define APP7 0xE7
#define APP8 0xE8
#define APP9 0xE9
#define APP10 0xEA
#define APP11 0xEB
#define APP12 0xEC
#define APP13 0xED
#define APP14 0xEE
#define APP15 0xEF
#define SOF0 0xC0
#define DQT 0xDB
#define DHT 0xC4
#define SOS 0xDA
#define EOI 0xD9

#endif
