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

typedef struct {
    byte horizontal_sample;
    byte vertical_sample;
    byte qantize_table_id;
} frame_component;

typedef struct {
    word height;
    word width;
    byte vmax;
    byte hmax;
    byte components_num;
    byte precision;
    frame_component frame_components[5];
} frame_data;

typedef struct huffman_leaf {
    word codeword;
    int codeword_len;
    byte value;
} huffman_leaf;

typedef struct {
    huffman_leaf* start;
    int ht_node_num;
} huffman_table_list;

int huffman_table_index (byte index)
{
    switch(index) {
        case 0x10:
            return 2;
        case 0x11:
            return 3;
        case 0x00:
            return 0;
        case 0x01:
            return 1;
        default:
            printf("huffman table index error");
            exit(1);
    }
}

double cos_cache[200];

void init_cos_cache()
{
    for (int i = 0; i < 200; i++) {
        cos_cache[i] = cos(i * M_PI / 16.0);
    }
}

int quantize_table_list[4][64];
frame_data f0;
huffman_table_list huffman_tables[4];
byte component_mapping_huffman[5];
double block[8][8];

void read_qt(FILE* fp)
{
    byte len = read_word_to_bigendian(fp);
    byte id_prec;
    len = len - 2;
    byte pq,id;
    byte c;
    word tmp;
    while (len > 0) {
        fread(&id_prec,1,1,fp);
        len-=1;
        pq = (id_prec >> 4) ? 1: 0;
        id = id_prec & 0x0f;
        for(int i = 0; i<64; i++) {
            if (pq) {
                tmp = read_word_to_bigendian(fp);
                len = len - 2;
            } else {
                fread(&c, 1, 1, fp);
                len--;
                tmp = c;
            }
            quantize_table_list[id][i] = tmp;
        }
    }
}

void read_frame(FILE *fp)
{
    int len = read_word_to_bigendian(fp);

    fread(&(f0.precision),1,1,fp);
    f0.height = read_word_to_bigendian(fp);
    f0.width = read_word_to_bigendian(fp);
    fread(&(f0.components_num),1,1,fp);

    byte component_id;
    byte sample;
    f0.vmax = 0x00;
    f0.hmax = 0x00;
    for (int i = 0; i < f0.components_num; i++) {
        fread(&component_id,1,1,fp);
        fread(&sample,1,1,fp);
        (f0.frame_components[component_id-1]).horizontal_sample = (sample >> 4) & 0x0f;
        (f0.frame_components[component_id-1]).vertical_sample = sample & 0x0f;
        if( (f0.frame_components[component_id-1]).horizontal_sample > f0.hmax) {
            f0.hmax = (f0.frame_components[component_id-1]).horizontal_sample;
        }
        if( (f0.frame_components[component_id-1]).vertical_sample > f0.vmax) {
            f0.vmax = (f0.frame_components[component_id-1]).vertical_sample;
        }
        fread(&((f0.frame_components[component_id-1]).qantize_table_id),1,1,fp);
    }
}

void read_ht(FILE* fp)
{
    byte huffman_length[16];
    int len = read_word_to_bigendian(fp);
    len = len - 2;
    int ht_node_num;
    byte ht_length;
    byte class_id;
    byte ht_value;
    while (len > 0) {
        fread(&class_id,1,1,fp);
        len--;
        ht_node_num = 0;
        for(int i = 0; i<16; i++) {
            fread(&(ht_length),1,1,fp);
            huffman_length[i] = ht_length;
            len--;
            ht_node_num += huffman_length[i];
        }

        huffman_tables[huffman_table_index(class_id)].ht_node_num = ht_node_num;
        huffman_leaf* ht_leaf = (huffman_leaf*)malloc((ht_node_num)*sizeof(
                                    huffman_leaf));
        word codeword = 0;
        for(int height=0, leaf_index = 0; height<16; height++) {
            for(int j = 0; j<huffman_length[height]; j++) {

                fread(&(ht_value),1,1,fp);
                len--;
                ht_leaf[leaf_index].value = ht_value;
                ht_leaf[leaf_index].codeword = (codeword++);
                ht_leaf[leaf_index].codeword_len = height+1;
                leaf_index++;
            }
            codeword = codeword<<1;
        }
        huffman_tables[huffman_table_index(class_id)].start = ht_leaf;
    }
}

void read_sos(FILE* fp)
{
    int len = read_word_to_bigendian(fp);
    byte component_num;
    fread(&component_num,1,1,fp);
    assert(len == 6+2*component_num);
    len = len - 3;
    assert(component_num==3);//JFIF defined
    byte component_id,destination;
    for (int i =0; i<component_num; i++) {

        fread(&component_id,1,1,fp);
        fread(&destination,1,1,fp);
        len = len-2;
        component_mapping_huffman[component_id-1] = destination;
    }
    fread(&destination,1,1,fp);
    fread(&destination,1,1,fp);
    fread(&destination,1,1,fp);
}

byte read_one_bit(FILE *fp)
{
    static byte buffer;
    static byte count = 0;
    byte check_ff00;
    if (!count) {
        fread(&buffer,1,1,fp);
        if (buffer == 0xff) {
            byte check_ff00;
            fread(&check_ff00,1,1,fp);
            if (check_ff00 != 0x00) {
                printf("missing 0xff00 sequence!");
                exit(1);
            }
        }
    }
    byte bit = (buffer >> (7 - count))&0x01;
    count = (count == 7 ? 0 : count + 1);
    return bit;
}


int codeword_decode (FILE *fp, byte need_read_length)
{
    byte leading = read_one_bit(fp);
    int decoding_code = 1;
    byte c;
    for (int i = 1; i < need_read_length; i++) {
        c = read_one_bit(fp);
        decoding_code = decoding_code << 1;
        decoding_code += leading ? c : !c;
    }
    decoding_code = leading ? decoding_code: -decoding_code;
    return decoding_code;
}

byte decode_huffman(FILE *fp, int huffman_table_id)
{
    huffman_leaf* huffman_table = huffman_tables[huffman_table_id].start;
    int ht_node_num = huffman_tables[huffman_table_id].ht_node_num;
    word codeword = 0x0000;
    for(int i = 0; i<16; i++) {
        codeword =codeword<<1;
        codeword |= (word)read_one_bit(fp);
        for(int j = 0; j<ht_node_num; j++) {
            if(huffman_table[j].codeword==codeword&&huffman_table[j].codeword_len==(i+1)) {
                return huffman_table[j].value;
            }
        }
    }
    printf("do not find\n");
    exit(1);
}
double c(int i)
{
    if (i == 0) {
        return sqrt(1.0/2.0);
    } else {
        return 1.0;
    }
}
void calcualte_mcu_block(FILE *fp, byte component_id)
{
    static int dc_block[5] = {0,0,0,0,0};
    memset(block,0,sizeof(block));
    int dc_table_id = huffman_table_index((component_mapping_huffman[component_id]>> 4)&0x0f);
    byte need_read_bit = decode_huffman(fp,dc_table_id);
    byte leading_bit;
    if (need_read_bit) {
        dc_block[component_id] += codeword_decode(fp,need_read_bit);
    }
    block[0][0] = dc_block[component_id];
    int ac_table_id = huffman_table_index((component_mapping_huffman[component_id]
                                           &0x0f)|0x10);
    byte zerosnum_needread,zerosnum;
    for (int i = 1; i< 64;) {
        zerosnum_needread = decode_huffman(fp,ac_table_id);
        if(zerosnum_needread == 0x00) {
            break;
        }
        if(zerosnum_needread == 0xf0) {
            i+=16;
            continue;
        }
        need_read_bit = zerosnum_needread & 0x0f;
        zerosnum = (zerosnum_needread >> 4) & 0x0f;
        for (byte j = 0; j< zerosnum; j++) {
            block[i/8][i%8] = 0.0;
            i++;
        }
        if (need_read_bit) {
            block[i/8][i%8] = codeword_decode(fp,need_read_bit);
            i++;
        }
    }
    for (int quantize_i = 0; quantize_i<64; quantize_i++) {
        block[quantize_i/8][quantize_i%8] *=
            quantize_table_list[(f0.frame_components[component_id]).qantize_table_id][quantize_i];
    }
    double tmp[8][8];
    for (int x = 0; x < 8; x++) {
        for (int y = 0; y < 8; y++) {
            tmp[x][y] = block[zigzag_index[x][y]/8][zigzag_index[x][y]%8];
        }
    }
    double s[8][8];
    memset(s,0,sizeof(s));
    for (int jj = 0; jj < 8; jj++) {
        for (int x = 0; x < 8; x++) {
            for (int y = 0; y < 8; y++) {
                s[jj][x] += c (y) * tmp[x][y] * cos_cache[(jj + jj + 1) * y];
            }
            s[jj][x] = s[jj][x] / 2.0;
        }
    }
    memset(block,0,sizeof(block));
    for (int ii = 0; ii < 8; ii++) {
        for (int jj = 0; jj < 8; jj++) {
            for (int x = 0; x < 8; x++) {
                block[ii][jj] += c(x) * s[jj][x] * cos_cache[(ii + ii + 1) * x];
            }
            block[ii][jj] = block[ii][jj] / 2.0;
            block[ii][jj] += 128.0;
        }
    }
}
void calculate_mcu(FILE* fp, const char* fwrite)
{
    int mcu_width = 8 * f0.hmax;
    int mcu_height = 8 * f0.vmax;
    int mcus_on_x = (f0.width - 1) / mcu_width + 1;
    int mcus_on_y = (f0.height - 1) / mcu_height + 1;
    mcu_small_block* data_unit = (mcu_small_block*) malloc(mcus_on_y*mcus_on_x*5*mcu_width*mcu_height*sizeof(
                                     mcu_small_block));

    for (int i  = 0; i<mcus_on_y; i++) {
        for (int j = 0; j<mcus_on_x; j++) {
            for (int component_id = 0; component_id<f0.components_num; component_id++)  {
                for (int a =0; a<f0.frame_components[component_id].vertical_sample; a++) {
                    for (int b = 0; b<f0.frame_components[component_id].horizontal_sample; b++) {
                        calcualte_mcu_block(fp, component_id);
                        for (int x = 0; x < 8; x++) {
                            for (int y = 0; y < 8; y++) {
                                data_unit[i*mcus_on_x*5*mcu_width*mcu_height
                                          + j*5*mcu_width*mcu_height
                                          + component_id*mcu_width*mcu_height
                                          + a*mcu_height
                                          + b][x][y] = block[x][y];
                            }
                        }
                    }
                }
            }
        }
    }
    rgb_element mcu_rgb[mcus_on_x*mcus_on_y][mcu_height][mcu_width];
    rgb_element rgb_image[f0.width][f0.height];
    int yV = f0.vmax/f0.frame_components[0].vertical_sample,
        yH = f0.hmax/f0.frame_components[0].horizontal_sample,
        cbV = f0.vmax/f0.frame_components[1].vertical_sample,
        cbH =f0.hmax/ f0.frame_components[1].horizontal_sample,
        crV = f0.vmax/f0.frame_components[2].vertical_sample,
        crH =f0.hmax/ f0.frame_components[2].horizontal_sample;
    for (int i  = 0; i<mcus_on_y; i++) {
        for (int j = 0; j<mcus_on_x; j++) {
            //1 = Y, 2 = Cb, 3 = Cr, 4 = I, 5 = Q
            for (int a = 0; a < mcu_height; a++) {
                for (int b = 0; b < mcu_width; b++) {
                    double Y = data_unit[i*mcus_on_x*5*mcu_width*mcu_height
                                         + j*5*mcu_width*mcu_height
                                         + 0*mcu_width*mcu_height
                                         + (a/(8*yV))*mcu_height
                                         + (b/(8*yH))][(a % (8*yV)) /yV][(b % (8*yH)) /yH],
                               Cb = data_unit[i*mcus_on_x*5*mcu_width*mcu_height
                                              + j*5*mcu_width*mcu_height
                                              + 1*mcu_width*mcu_height
                                              + (a/(8*cbV))*mcu_height
                                              + b/(8*cbH)][(a % (8*cbV)) /cbV][(b % (8*cbH)) /cbH],
                                    Cr = data_unit[i*mcus_on_x*5*mcu_width*mcu_height
                                                   + j*5*mcu_width*mcu_height
                                                   + 2*mcu_width*mcu_height
                                                   + (a/(8*crV))*mcu_height
                                                   + b/(8*crH)][(a % (8*crV)) /crV][(b % (8*crH)) /crH];
                    double R = Y + 1.402*(Cr - 128),
                           G = Y - 0.34414*(Cb - 128) - 0.71414*(Cr -  128),
                           B = Y + 1.772*(Cb - 128);
                    mcu_rgb[i*mcus_on_x + j][a][b].r = (byte)(R > 255.0 ? 255.0 : (R < 0.0 ? 0.0 : R));
                    mcu_rgb[i*mcus_on_x + j][a][b].g = (byte)(G > 255.0 ? 255.0 : (G < 0.0 ? 0.0 : G));
                    mcu_rgb[i*mcus_on_x + j][a][b].b = (byte)(B > 255.0 ? 255.0 : (B < 0.0 ? 0.0 : B));
                }
            }
        }
    }
    free(data_unit);
    bitmap_image outimg(f0.width,f0.height);

    for (int y = 0; y < f0.height; y++) {
        for (int x = 0; x < f0.width; x++) {
            outimg.set_pixel(x,y,mcu_rgb[y/mcu_height * mcus_on_x + x/mcu_width][y%mcu_height][x%mcu_width].r
                             ,mcu_rgb[y/mcu_height * mcus_on_x + x/mcu_width][y%mcu_height][x%mcu_width].g
                             ,mcu_rgb[y/mcu_height * mcus_on_x + x/mcu_width][y%mcu_height][x%mcu_width].b);
        }
    }
    outimg.save_image(fwrite);
}

int main(int argc,char* argv[])
{
    if (argc != 2  && argc != 3) {
        printf("[ERROR]:\nusage: ./decoder <input_file> [<output_file_name>]");
        exit(1);
    }
    FILE* fp;
    if ((fp = fopen(argv[1], "r")) == NULL) {
        printf("%s can't be opened\n", argv[1]);
        exit(1);
    }
    const char* fwrite = argc ==3 ? argv[2] : "image.bmp";

    byte h;
    byte l;
    bool b_SOI = false;
    bool b_EOI = false;
    init_cos_cache();
    while (fread(&h,1, 1,fp)) {
        fread(&l,1, 1,fp);
        if (h == 0xff) {
            if (!b_SOI && l == SOI) {
                b_SOI = true;
            }
            assert(b_SOI == true);
            switch(l) {
            case DQT:
                read_qt(fp);
                break;
            case SOF0:
                read_frame(fp);
                break;
            case DHT:
                read_ht(fp);
                break;
            case SOS:
                read_sos(fp);
                calculate_mcu(fp,fwrite);
                break;
            case APP0:
            case APP1:
            case APP2:
            case APP3:
            case APP4:
            case APP5:
            case APP6:
            case APP7:
            case APP8:
            case APP9:
            case APP10:
            case APP11:
            case APP12:
            case APP13:
            case APP14:
            case APP15:
                break;
            case EOI:
                b_EOI = true;
                break;
            }
        }
    }
    assert(b_EOI==true);
    return 0;
}
