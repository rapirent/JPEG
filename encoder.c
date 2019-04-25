#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include "util.h"
#include "./bitmap/bitmap_image.hpp"


// Implementation of LLM DCT.
void llm_dct(const double in[8], double out[8]) {
    // Constants:
    const double s1 = sin(1. * M_PI / 16.);
    const double c1 = cos(1. * M_PI / 16.);
    const double s3 = sin(3. * M_PI / 16.);
    const double c3 = cos(3. * M_PI / 16.);
    const double r2s6 = sqrt(2.) * sin(6. * M_PI / 16.);
    const double r2c6 = sqrt(2.) * cos(6. * M_PI / 16.);

    // After stage 1:
    const double s1_0 =  in[0] + in[7];
    const double s1_1 =  in[1] + in[6];
    const double s1_2 =  in[2] + in[5];
    const double s1_3 =  in[3] + in[4];
    const double s1_4 = -in[4] + in[3];
    const double s1_5 = -in[5] + in[2];
    const double s1_6 = -in[6] + in[1];
    const double s1_7 = -in[7] + in[0];

    // After stage 2:
    const double s2_0 =  s1_0 + s1_3;
    const double s2_1 =  s1_1 + s1_2;
    const double s2_2 = -s1_2 + s1_1;
    const double s2_3 = -s1_3 + s1_0;

    const double z1 = c3 * (s1_7 + s1_4);
    const double s2_4 = ( s3-c3) * s1_7 + z1;
    const double s2_7 = (-s3-c3) * s1_4 + z1;

    const double z2 = c1 * (s1_6 + s1_5);
    const double s2_5 = ( s1-c1) * s1_6 + z2;
    const double s2_6 = (-s1-c1) * s1_5 + z2;

    // After stage 3:
    const double s3_0 =  s2_0 + s2_1;
    const double s3_1 = -s2_1 + s2_0;

    const double z3 = r2c6 * (s2_3 + s2_2);
    const double s3_2 = ( r2s6-r2c6) * s2_3 + z3;
    const double s3_3 = (-r2s6-r2c6) * s2_2 + z3;

    const double s3_4 =  s2_4 + s2_6;
    const double s3_5 = -s2_5 + s2_7;
    const double s3_6 = -s2_6 + s2_4;
    const double s3_7 =  s2_7 + s2_5;

    // After stage 4:
    const double s4_4 = -s3_4 + s3_7;
    const double s4_5 =  s3_5 * sqrt(2.);
    const double s4_6 =  s3_6 * sqrt(2.);
    const double s4_7 =  s3_7 + s3_4;

    // Shuffle and scaling:
    out[0] = s3_0 / sqrt(8.);
    out[4] = s3_1 / sqrt(8.);
    out[2] = s3_2 / sqrt(8.);
    out[6] = s3_3 / sqrt(8.);
    out[7] = s4_4 / sqrt(8.);
    out[3] = s4_5 / sqrt(8.);  // Alternative: s3_5 / 2
    out[5] = s4_6 / sqrt(8.);
    out[1] = s4_7 / sqrt(8.);
}

double cos_cache[200];

void init_cos_cache()
{
    for (int i = 0; i < 200; i++) {
        cos_cache[i] = cos(i * M_PI / 16.0);
    }
}

double c(int i)
{
    if (i == 0) {
        return sqrt(1.0/2.0);
    } else {
        return 1.0;
    }
}

void zero_rle() {

}
void canonical_huffman() {

}
void encode() {

}


int image_height,image_width;

void calculate_mcu_block(int block[][8]) {
    double tmp[8][8],s[8][8];
    printf("Start fdct!!!\n");
    memset(tmp,0,sizeof(tmp));
    memset(s,0,sizeof(s));

    for (int jj = 0; jj < 8; jj++) {
        for (int x = 0; x < 8; x++) {
            for (int y = 0; y < 8; y++) {
                s[jj][x] += c (y) * tmp[x][y] * cos_cache[(jj*2 + 1) * y];
            }
            s[jj][x] = s[jj][x] / 2.0;
        }
    }

    for (int ii = 0; ii < 8; ii++) {
        for (int jj = 0; jj < 8; jj++) {
            for (int x = 0; x < 8; x++) {
                block[ii][jj] += c(x) * s[jj][x] * cos_cache[(ii*2 + 1) * x];
            }
            tmp[ii][jj] = tmp[ii][jj] / 2.0;
            tmp[ii][jj]-=128;
        }
    }
    printf("zigzag arrange start!!!\n");
    for (int x = 0; x < 8; x++) {
        for (int y = 0; y < 8; y++) {
            s[x][y] = tmp[zigzag_index[x][y]][zigzag_index[x][y]];
        }
    }

    // 4. 量化
    //     對於前面得到的 64 個空間頻率振幅值, 我們將對它們作幅度分層量化操作.方
    // 法就是分別除以量化表裡對應值並四捨五入.

    for (int i = 0 ; i<=63; i++ ) {
        block[i/8][i%8] = (int) (s[i/8][i%8] / quantize_table[i/8][i%8] + 0.5);
    }

}

void calculate_mcu(int block_x, int block_y, yuv_element* yuv_image) {

    // JPEG 裡是對每 8x8 個點為一個單位處理的.
    // 所以如果原始圖片的長寬不是 8 的倍數, 都需要先補成 8 的倍數, 好一塊塊的處理.
    // 另外, 記得剛才我說的 Cr Cb 都是 2x2 記錄一次嗎?
    // 所以大多數情況, 是要補成 16x16 的整數塊.按從左到右, 從上到下的次序排列
    //  (和我們寫字的次序一樣). JPEG 裡是對 Y Cr Cb 分別做 DCT 變換的. 這裡進行 DCT 變換
    // 的 Y, Cr, Cb 值的範圍都是 -128~127. (Y 被減去 128)
    int mcu[8][8];
    for (int yuv_id = 0;yuv_id<3;yuv_id++) {
        for (int y=0;y<8;y++) {
            for (int x=0;x<8;x++) {
                // JPEG 裡是對每 8x8 個點為一個單位處理的.
                // 所以如果原始圖片的長寬不是 8 的倍數, 都需要先補成 8 的倍數, 好一塊塊的處理.
                if (block_y - y>=0 && block_x + x < image_width) {
                    switch(yuv_id) {
                        case 0:
                            mcu[y][x] = yuv_image[(block_y - y) * image_width + block_x + x].y;
                            break;
                        case 1:
                            mcu[y][x] = yuv_image[(block_y - y) * image_width + block_x + x].cb;
                            break;
                        case 2:
                            mcu[y][x] = yuv_image[(block_y - y) * image_width + block_x + x].cr;
                            break;
                    }
                }
                else {
                    mcu[y][x] = 0.0;
                }
            }
        }
        calculate_mcu_block(mcu);
        zero_rle();
        canonical_huffman();
        encode();
    }
}



int main(int argc, char* argv[]) {
    if (argc != 2  && argc != 3) {
        printf("[ERROR]:\nusage: ./encoder <input_file> [<output_file_name>]");
        exit(1);
    }
    bitmap_image image(argv[1]);
    if (!image) {
        printf("%s can't be opened\n", argv[1]);
        exit(1);
    }
    FILE* fp;
    if (argc==3 && (fp = fopen(argv[2], "w")) == NULL) {
        printf("%s can't be opened\n", argv[2]);
        exit(1);
    }
    else if ((fp = fopen("image.jpeg", "w")) == NULL) {
        printf("image.jpeg can't be opened\n");
        exit(1);
    }
    byte c;
    c = 0xff;
    fwrite(&c, 1, 1, fp);
    c = SOI;
    fwrite(&c, 1, 1, fp);
    c = 0xff;
    fwrite(&c, 1, 1, fp);
    c = APP0;
    fwrite(&c, 1, 1, fp);

    image_height = image.height();
    image_width  = image.width();

    yuv_element* yuv_image = (yuv_element*) malloc(image_height*image_width*sizeof(yuv_element));

    word R,G,B;
    //color model
    rgb_t colour;
    for (int y = 0; y < image_height; y++)
    {
        for (int x = 0; x < image_width; x++)
        {

            // Y = 0.299*R + 0.587*G + 0.114*B  (亮度)
            // Cb =  - 0.1687*R - 0.3313*G + 0.5   *B + 128
            // Cr =    0.5   *R - 0.4187*G - 0.0813*B + 128
            // 因為人眼對圖片上的亮度 Y 的變化遠比色度 C 的變化敏感. 我們完全可以每個點保存一個 8bit 的亮
            // 度值, 每 2x2 個點保存一個 Cr Cb 值, 而圖像在肉眼中的感覺不會起太大的變化.
            // 所以, 原來用 RGB 模型, 4 個點需要 4x3=12 字節. 而現在僅需要 4+2=6 字節; 平
            // 均每個點佔 12bit. 當然 JPEG 格式裡允許每個點的 C 值都記錄下來; 不過 MPEG 裡
            // 都是按 12bit 一個點來存放的, 我們簡寫為 YUV12.
            image.get_pixel(x, y, colour);
            yuv_image[x*image_width+y].y = 0.299*colour.red + 0.587*colour.green + 0.114*colour.blue;
            yuv_image[x*image_width+y].cb =  - 0.1687*colour.red - 0.3313*colour.green + 0.5   *colour.blue + 128;
            yuv_image[x*image_width+y].cr =    0.5   *colour.red - 0.4187*colour.green - 0.0813*colour.blue + 128;
        }
    }
    // JPEG 裡, 要對數據壓縮, 先要做一次 DCT 變換. DCT 變換的原理, 涉及到數學
    // 知識, 這裡我們不必深究. 反正和傅立葉變換(學過高數的都知道) 是差不多了. 經過
    // 這個變換, 就把圖片裡點和點間的規律呈現出來了, 更方便壓縮.JPEG 裡是對每 8x8
    // 個點為一個單位處理的. 所以如果原始圖片的長寬不是 8 的倍數, 都需要先補成 8
    // 的倍數, 好一塊塊的處理. 另外, 記得剛才我說的 Cr Cb 都是 2x2 記錄一次嗎? 所
    // 以大多數情況, 是要補成 16x16 的整數塊.按從左到右, 從上到下的次序排列 (和我
    // 們寫字的次序一樣). JPEG 裡是對 Y Cr Cb 分別做 DCT 變換的. 這裡進行 DCT 變換
    // 的 Y, Cr, Cb 值的範圍都是 -128~127. (Y 被減去 128)
    init_cos_cache();
    for (int y = image_height - 1; y>=0;y-=8) {
        for (int x = 0;x<image_width;x+=8) {
            calculate_mcu(x,y,yuv_image);
        }
    }

    // 壓縮算法簡介
    // 1. 色彩模型
    // 2. DCT (離散餘弦變換)
    // 3. 重排列 DCT 結果
    // 4. 量化
    // 5. 0 RLE 編碼
    // 6. 範式 Huffman 編碼
    // 7. DC 的編碼
    return 0;
}
