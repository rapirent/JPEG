#ifndef _JPEG_HEADER
#define _JPEG_HEADER
#include <stdint.h>
#include <inttypes.h>

typedef uint8_t byte;

// typedef struct {
//     //big endian
//     byte higher_byte;
//     byte lowwer_byte;
// }word_unit;

typedef uint16_t word;

const int zigzag_index[8][8] = {
    0, 1, 5, 6,14,15,27,28,
    2, 4, 7,13,16,26,29,42,
    3, 8,12,17,25,30,41,43,
    9,11,18,24,31,40,44,53,
    10,19,23,32,39,45,52,54,
    20,22,33,38,46,51,55,60,
    21,34,37,47,50,56,59,61,
    35,36,48,49,57,58,62,63
};

const int quantize_table[8][8] = {
    16, 11, 10, 16, 24,  40, 51, 61,
    12, 12, 14, 19, 26,  58, 60, 55,
    14, 13, 16, 24, 40,  57, 69, 56,
    14, 17, 22, 29, 51,  87, 80, 62,
    18, 22, 37, 56, 68,  109, 103, 77,
    24, 35, 55, 64, 81,  104, 113, 92,
    49, 64, 78, 87, 103, 121, 120, 101,
    72, 92, 95, 98, 112, 100, 103, 99
};

typedef struct {
    byte r;
    byte g;
    byte b;
}rgb_element;
// 我們完全可以每個點保存一個 8bit 的亮度值, 每 2x2 個點保存一個 Cr Cb 值, 而圖像在肉眼中的感覺不會起太大的變化.
// 原來用 RGB 模型, 4 個點需要 4x3=12 字節. 而現在僅需要 4+2=6 字節
typedef struct {
    byte y[4];
    byte cb;
    byte cr;
}yuv12_element;

typedef struct {
    byte y;
    byte cb;
    byte cr;
}yuv_element;

#define BYTE_TO_BINARY_PATTERN "%c%c%c%c%c%c%c%c"
#define BYTE_TO_BINARY(byte)  \
  (byte & 0x0080 ? '1' : '0'), \
  (byte & 0x0040 ? '1' : '0'), \
  (byte & 0x0020 ? '1' : '0'), \
  (byte & 0x0010 ? '1' : '0'), \
  (byte & 0x0008 ? '1' : '0'), \
  (byte & 0x0004 ? '1' : '0'), \
  (byte & 0x0002 ? '1' : '0'), \
  (byte & 0x0001 ? '1' : '0')

typedef double mcu_small_block[8][8];

uint16_t read_word_to_bigendian(FILE* fp)
{
    byte high,low;
    fread(&high, 1, 1, fp);
    fread(&low, 1, 1, fp);
    // printf("get a bigendian word: %.2x%.2x\n", high, low);
    // 0000 0000 0100 0011
    //
    // byte high = c.higher_byte;
    //exchange the byte: Motorola format
    word word_value = 0x0000;
    word_value = ((word_value | high) << 8) &0xff00;
    word_value = word_value | low;
    // 0000 0000 0000 0000
    // byte low = c.lowwer_byte;
    // word low = (c & 0x00ff << 8);
    // 0000 0000 0100 0011
    return word_value;
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
