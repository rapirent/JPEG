#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include "util.h"
#include "bitmap_image.hpp"

int image_height,image_width;
byte category[3] = {0,1,1};
byte zigzag_quantize_table[2][64];
bitmap_image rgb_image;

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

void write_one_bit(FILE *fp, int value, bool mode)
{
    static byte buffer = 0x00;
    static byte count = 0;
    if(mode==true) {
        for (int i = count; i<7; i++) {
            buffer = (buffer << 1);
        }
        fwrite(&buffer,1,1,fp);
        return;
    }
    if (value != 1 && value != 0) {
        printf("bit number should be a binary\n");
        exit(1);
    }
    buffer = (buffer << 1) | value;
    if(count==7) {
        fwrite(&buffer,1,1,fp);
        if(buffer == 0xff) {
            byte c = 0x00;
            fwrite(&c,1,1,fp);
        }
        buffer = 0x00;
    }
    count = (count == 7 ? 0 : count + 1);
}

void codeword_encode (FILE *fp, int value, byte length)
{
    for (int i = length-1; i>=0 ; i--) {
        write_one_bit(fp,(value >> i) & 0x01,false);
    }
}

byte find_length (int* value)
{
    byte need_write_length = 0;
    int tmp = (*value);
    if (*value < 0) {
        tmp = -(*value);
    }
    for (; tmp ; tmp = tmp >> 1) {
        need_write_length++;
    }
    *value = (*value) > 0 ? (*value) : ((0x1 << need_write_length) + (*value) -1);
    return need_write_length;
}
void calculate_mcu_block(FILE *fp,char block[][8],int yuv_id)
{
    double tmp[8][8], s[8][8];
    short quantized_block[64];
    memset(tmp,0,sizeof(tmp));
    memset(s,0,sizeof(s));
    for (int ii = 0; ii < 8; ii++) {
        for (int x = 0; x < 8; x++) {
            for (int y = 0; y < 8; y++) {
                s[ii][x] += (double)block[ii][y] * cos_cache[(y*2 + 1) * x];
            }
            s[ii][x] = c(x) *s[ii][x] / 2.0;
        }
    }
    for (int ii = 0; ii < 8; ii++) {
        for (int jj = 0; jj < 8; jj++) {
            for (int x = 0; x < 8; x++) {
                tmp[ii][jj] += s[x][jj] * cos_cache[(x*2 + 1) * ii];
            }
            tmp[ii][jj] = c(ii) * tmp[ii][jj] / 2.0;
        }
    }
    for (int x = 0; x < 8; x++) {
        for (int y = 0; y < 8; y++) {
            quantized_block[zigzag_index[x][y]] = (short)round(tmp[x][y]/(double)quantize_table[category[yuv_id]][x*8+y]);
        }
    }
    static short dc_block[5] = {0,0,0,0,0};
    int diff = (int)(quantized_block[0] - dc_block[yuv_id]);
    int need_write_length;
    int index;
    dc_block[yuv_id] = quantized_block[0];
    if (!diff) {
        for (int i = 0; i<strlen((char*)dc_diff[category[yuv_id]][0]); i++) {
            write_one_bit(fp, dc_diff[category[yuv_id]][0][i] - '0',false);
        }
    } else {
        need_write_length = find_length(&diff);
        for (int i = 0; i<strlen((char*)dc_diff[category[yuv_id]][need_write_length]); i++) {
            write_one_bit(fp, dc_diff[category[yuv_id]][need_write_length][i] - '0',false);
        }
        codeword_encode(fp,diff,need_write_length);
    }
    int count_zero = 0;
    int value;
    byte cof_index;
    int padding_start_index=63;
    while (padding_start_index>1&&!quantized_block[padding_start_index]) {
        padding_start_index--;
    }
    for (int i = 1; i <= padding_start_index; i++) {
        if (!quantized_block[i]) {
            count_zero++;
            continue;
        }
        if (count_zero>=16) {
            for (int j = 0; j<count_zero/16; j++) {
                for (int k = 0; k<strlen((char*)ac_cof[category[yuv_id]][0xf0]); k++) {
                    write_one_bit(fp,ac_cof[category[yuv_id]][0xf0][k] - '0',false);
                }
            }
            count_zero = count_zero % 16;
        }
        value = quantized_block[i];
        need_write_length = find_length(&value);
        for (int j = 0; j<strlen((char*)ac_cof[category[yuv_id]][((count_zero << 4)&0xf0)|need_write_length]); j++) {
            write_one_bit(fp,ac_cof[category[yuv_id]][((count_zero << 4)&0xf0)|need_write_length][j] - '0',false);
        }
        codeword_encode(fp,value,need_write_length);
        count_zero = 0;
    }
    if (padding_start_index!=63) {
        for (int i = 0; i<strlen((char*)ac_cof[category[yuv_id]][0x00]); i++) {
            write_one_bit(fp,ac_cof[category[yuv_id]][0x00][i] - '0',false);
        }
    }
}

void calculate_mcu(FILE *fp,int block_x, int block_y)
{

    rgb_t colour;
    char mcu[8][8];
    byte R,G,B;
    for (int yuv_id = 0; yuv_id<3; yuv_id++) {
        for (int y=0; y<8; y++) {
            for (int x=0; x<8; x++) {
                if (image_width - block_x < 7 && image_height - block_y < 7) {
                    rgb_image.get_pixel((x + block_x - 8), (y + block_y - 8), colour);
                } else if (image_width - block_x < 7) {
                    rgb_image.get_pixel((x + block_x - 8), y, colour);
                } else if(image_height - block_y < 7) {
                    rgb_image.get_pixel(x, (y + block_y - 8), colour);
                } else {
                    rgb_image.get_pixel(x+block_x, y+block_y, colour);
                }
                switch(yuv_id) {
                    case 0:
                        mcu[y][x] = (byte)(0.299 * colour.red + 0.587 * colour.green + 0.114 * colour.blue - 128);
                        break;
                    case 1:
                        mcu[y][x] = (byte)(-0.1687 * colour.red - 0.3313 * colour.green + 0.5 * colour.blue);
                        break;
                    case 2:
                        mcu[y][x] = (byte)(0.5 * colour.red - 0.4187 * colour.green - 0.0813 * colour.blue);
                        break;
                }
            }
        }
        calculate_mcu_block(fp,mcu,yuv_id);
    }

}

int main(int argc, char* argv[])
{
    if (argc != 2  && argc != 3) {
        printf("[ERROR]:\nusage: ./encoder <input_file> [<output_file_name>]");
        exit(1);
    }


    FILE* fp;
    if (argc==3 && (fp = fopen(argv[2], "w")) == NULL) {
        printf("%s can't be opened\n", argv[2]);
        exit(1);
    } else if (argc==2 && (fp = fopen("image.jpg", "w")) == NULL) {
        printf("image.jpeg can't be opened\n");
        exit(1);
    }

    bitmap_image image(argv[1]);
    if (!image) {
        printf("%s can't be opened\n", argv[1]);
        exit(1);
    }

    image_height = image.height();
    image_width  = image.width();
    rgb_image = image;

    byte c;
    c = 0xff;
    fwrite(&c, 1, 1, fp);
    c = SOI;
    fwrite(&c, 1, 1, fp);
    c = 0xff;
    fwrite(&c, 1, 1, fp);
    c = APP0;
    fwrite(&c, 1, 1, fp);
    write_word_to_littleendian(fp, 16);
    char comment[5] = "JFIF";
    fwrite(comment, 1, 5, fp);
    c = 0x01;
    fwrite(&c, 1, 1, fp);
    c = 0x01;
    fwrite(&c, 1, 1, fp);
    c = 0x00;
    fwrite(&c, 1, 1, fp);
    c = 0x00;
    fwrite(&c, 1, 1, fp);
    c = 0x01;
    fwrite(&c, 1, 1, fp);
    c = 0x00;
    fwrite(&c, 1, 1, fp);
    c = 0x01;
    fwrite(&c, 1, 1, fp);
    c = 0x00;
    fwrite(&c, 1, 1, fp);
    c = 0x00;
    fwrite(&c, 1, 1, fp);

    c = 0xff;
    fwrite(&c, 1, 1, fp);
    c = DQT;
    fwrite(&c, 1, 1, fp);
    write_word_to_littleendian(fp, 132);
    c = 0x00;
    fwrite(&c, 1, 1, fp);
    for (int i = 0; i<8; i++) {
        for (int j = 0; j<8; j++) {
            zigzag_quantize_table[0][zigzag_index[i][j]] = quantize_table[0][i*8+j];
            zigzag_quantize_table[1][zigzag_index[i][j]] = quantize_table[1][i*8+j];
        }
    }
    for (int i = 0; i<8; i++) {
        for (int j = 0; j<8; j++) {
            fwrite(&(zigzag_quantize_table[0][i*8+j]),1,1,fp);
        }
    }
    c = 0x01;
    fwrite(&c, 1, 1, fp);
    for (int i = 0; i<8; i++) {
        for (int j = 0; j<8; j++) {
            fwrite(&(zigzag_quantize_table[1][i*8+j]),1,1,fp);
        }
    }
    c = 0xff;
    fwrite(&c, 1, 1, fp);
    c = SOF0;
    fwrite(&c, 1, 1, fp);
    write_word_to_littleendian(fp, 17);
    c = 8;
    fwrite(&c, 1, 1, fp);
    write_word_to_littleendian(fp,(word)image_height);
    write_word_to_littleendian(fp,(word)image_width);
    c = 3;
    fwrite(&c, 1, 1, fp);
    c = 1;
    fwrite(&c, 1, 1, fp);
    c = 0x11;
    fwrite(&c, 1, 1, fp);
    c = 0x00;
    fwrite(&c, 1, 1, fp);
    c = 2;
    fwrite(&c, 1, 1, fp);
    c = 0x11;
    fwrite(&c, 1, 1, fp);
    c = 1;
    fwrite(&c, 1, 1, fp);
    c = 3;
    fwrite(&c, 1, 1, fp);
    c = 0x11;
    fwrite(&c, 1, 1, fp);
    c = 1;
    fwrite(&c, 1, 1, fp);

    c = 0xff;
    fwrite(&c, 1, 1, fp);
    c = DHT;
    fwrite(&c, 1, 1, fp);
    write_word_to_littleendian(fp, 0x01a2);
    c = 0x00;
    fwrite(&c, 1, 1, fp);
    fwrite(dc_cof_code_length[0], 1, 16, fp);
    fwrite(dc_cof_value[0], 1, 12, fp);
    c = 0x01;
    fwrite(&c, 1, 1, fp);
    fwrite(dc_cof_code_length[1], 1, 16, fp);
    fwrite(dc_cof_value[1], 1, 12, fp);
    c = 0x10;
    fwrite(&c, 1, 1, fp);
    fwrite(ac_cof_code_length[0], 1, 16, fp);
    fwrite(ac_cof_value[0], 1, 162, fp);
    c = 0x11;
    fwrite(&c, 1, 1, fp);
    fwrite(ac_cof_code_length[1], 1, 16, fp);
    fwrite(ac_cof_value[1], 1, 162, fp);

    c = 0xff;
    fwrite(&c, 1, 1, fp);
    c = SOS;
    fwrite(&c, 1, 1, fp);
    write_word_to_littleendian(fp, 12);
    c = 3;
    fwrite(&c, 1, 1, fp);
    c = 1;
    fwrite(&c, 1, 1, fp);
    c = 0x00;
    fwrite(&c, 1, 1, fp);
    c = 2;
    fwrite(&c, 1, 1, fp);
    c = 0x11;
    fwrite(&c, 1, 1, fp);
    c = 3;
    fwrite(&c, 1, 1, fp);
    c = 0x11;
    fwrite(&c, 1, 1, fp);
    c = 0x00;
    fwrite(&c, 1, 1, fp);
    c = 0x3f;
    fwrite(&c, 1, 1, fp);
    c = 0x00;
    fwrite(&c, 1, 1, fp);

    init_cos_cache();
    for (int y = 0; y<image_height; y+=8) {
        for (int x = 0; x<image_width; x+=8) {
            calculate_mcu(fp,x,y);
        }
    }
    write_one_bit(fp,0,true);
    c = 0xff;
    fwrite(&c, 1, 1, fp);
    c = EOI;
    fwrite(&c, 1, 1, fp);
    fclose(fp);
    return 0;
}
