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

uint16_t read_word_to_bigendian(FILE* fp) {
    word_unit c;
    fscanf(fp,"%"SCNd16,&c);
    byte high = c.higher_byte;
    byte low = c.lowwer_byte;
    high = high << 8;
    return high|low;
}

#define START 0xFF
#define SOI 0xD8
#define APP0 0xE0
#define APP1 0xE1
#define APP15 0xEF
#define SOF0 0xC0
#define DQT 0xDB
#define DHT 0xC4
#define SOS 0xDA
#define EOI 0xD9

#endif
