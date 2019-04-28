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

int image_height,image_width;
byte category[3] = {0,1,1}; // table for luminance or chrominance
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
    // printf("receive write request for %d\n",value);
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
    // assert(length==0);
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
    //add leading 1 if negative
    *value = (*value) > 0 ? (*value) : ((0x1 << need_write_length) + (*value) -1);
    return need_write_length;
}
void calculate_mcu_block(FILE *fp,char block[][8],int yuv_id)
{
    double tmp[8][8], s[8][8];
    // printf("Start fdct!!!\n");
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
    // printf("\n---------\n");
    // printf("\n---------\n");
    // printf("zigzag arrange start!!!\n");
    // for (int x = 0; x < 8; x++) {
    //     for (int y = 0; y < 8; y++) {
    //         // printf("%d %d\n",x,y);
    //         quantized_block[zigzag_index[x][y]] = (int)block[x][y];
    //     }
    // }

    // 4. 量化
    //     對於前面得到的 64 個空間頻率振幅值, 我們將對它們作幅度分層量化操作.方
    // 法就是分別除以量化表裡對應值並四捨五入.
    // printf("quantization start!!!\n");
    for (int x = 0; x < 8; x++) {
        for (int y = 0; y < 8; y++) {
            quantized_block[zigzag_index[x][y]] = (short)round(tmp[x][y]/(double)zigzag_quantize_table[category[yuv_id]][x*8+y]);
        }
    }
    // printf("gogogo!\n");
    // DC即一塊圖像樣本的平均值. 就是說, 它包含了原始 8x8 圖像塊裡的很多能量. (通常
    // 會得到一個很大的數值)
    static short dc_block[5] = {0,0,0,0,0};
    // JPEG 的作者指出連續塊的 DC 率之間有很緊密的聯繫,  因此他們決定對 8x8 塊的
    // DC 值的差別進行編碼. (Y, Cb, Cr 分別有自己的 DC)
    // Diff = DC(i)  - DC(i-1)
    // 所以這一塊的 DC(i) 就是:  DC(i)  = DC(i-1)  + Diff
    int diff = (int)(quantized_block[0] - dc_block[yuv_id]);
    // printf("dc encode\n");
    int need_write_length;
    int index;
    dc_block[yuv_id] = quantized_block[0];//dc_block[yuv_id] + diff;
    if (!diff) {
        for (int i = 0; i<strlen((char*)dc_diff[category[yuv_id]][0]); i++) {
            write_one_bit(fp, dc_diff[category[yuv_id]][0][i] - '0',false);
        }
    } else {
        need_write_length = find_length(&diff);
        // printf("length = %d, number = %d\n",need_write_length,diff);
        // printf("%s\n----\n",dc_diff[category[yuv_id]][need_write_length]);
        for (int i = 0; i<strlen((char*)dc_diff[category[yuv_id]][need_write_length]); i++) {
            write_one_bit(fp, dc_diff[category[yuv_id]][need_write_length][i] - '0',false);
        }
        codeword_encode(fp,diff,need_write_length);
    }
    // 例如上面例子中, Diff 是 -511, 就編碼成 ((0,9), 000000000)
    // 如果 9 的 Huffman 編碼是 1111110
    //(在 JPG 文件中, 一般有兩個 Huffman 表, 一個是 DC 用, 一個是 AC 用)
    //那麼在 JPG 文件中, DC 的 2 進製表示為 1111110 000000000
    // (0,57) ; (0,45) ; (4,23) ; (1,-30) ; (0,-8) ; (2,1) ; (0,0)
    // 只處理每對數右邊的那個:
    // 57 是第 6 組的, 實際保存值為 111001 , 所以被編碼為 (6,111001)
    // 45 , 同樣的操作, 編碼為 (6,101101)
    //     23  ->  (5,10111)
    //    -30  ->  (5,00001)
    //     -8  ->  (4,0111)
    //      1  ->  (1,1)

    // 前面的那串數字就變成了:
    // (0,6), 111001 ; (0,6), 101101 ; (4,5), 10111; (1,5), 00001; (0,4) , 0111 ; (2,1), 1 ; (0,0)

    // 括號裡的數值正好合成一個word. 後面被編碼的數字表示範圍是  -32767..32767.
    // 合成的字節裡, 高 4 位是前續 0 的個數, 低 4 位描述了後面數字的位數.
    // for (int i = 0;!(value >> 15*i) ; value = value >> (15*i)) {


    //AC
    // 現在我們矢量中有許多連續的 0. 我們可以使用 RLE 來壓縮掉這些 0. 這裡我們
    // 將跳過第一個矢量 (後面將解釋為什麼) 因為它的編碼比較特別. 假設有一組矢量
    // (64 個的後 63 個) 是
    //     57,45,0,0,0,0,23,0,-30,-16,0,0,1,0,0,0, 0 , 0 ,0 , 0,..,0
    // 經過 RLE 壓縮後就是
    //     (0,57) ; (0,45) ; (4,23) ; (1,-30) ; (0,-16) ; (2,1) ; EOB
    // 我們用 (0,0) 表示 EOB
    // 但是, 如果這組數字不以 0 結束,  那麼就不需要 EOB.
    // 另外需要注意的是, 由於後面 huffman 編碼的要求, 每組數字前一個表示 0 的
    // 數量的必須是 4 bit, 就是說, 只能是 0~15, 所以, 如果有這麼一組數字:
    //     57, 十八個0, 3, 0, 0, 0, 0, 2, 三十三個0, 895, EOB
    // 我們實際這樣編碼:
    //     (0,57) ; (15,0) (2,3) ; (4,2) ; (15,0) (15,0) (1,895) , (0,0)
    // 注意 (15,0) 表示了 16 個連續的 0.
    // printf("ac encode\n");
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
                // printf("???");
                // printf("%s\n----\n",ac_cof[category[yuv_id]][0xf0]);
                for (int k = 0; k<strlen((char*)ac_cof[category[yuv_id]][0xf0]); k++) {
                    write_one_bit(fp,ac_cof[category[yuv_id]][0xf0][k] - '0',false);
                }
                // fwrite(fp,1,16,ac_cof[category[yuv_id]][0xF][0x0]);
            }
            count_zero = count_zero % 16;
        }
        value = quantized_block[i];
        need_write_length = find_length(&value);
        // c = (count_zero << 4) & 0xf0;
        // 括號裡的數值正好合成一個word. 後面被編碼的數字表示範圍是  -32767..32767.
        // 合成的字節裡, 高 4 位是前續 0 的個數, 低 4 位描述了後面數字的位數.
        // c = c | ((length = find_length(block[i/8][i%8]) >> 4) & 0x0f);
        // printf("need_length = %d\n",need_write_length);
        // printf("count_zero is %u, match word = %d\n", (count_zero << 4)&0xf0,need_write_length);
        // printf("%s\n----\n",ac_cof[category[yuv_id]][(count_zero << 4)&0xf0|need_write_length]);
        for (int j = 0; j<strlen((char*)ac_cof[category[yuv_id]][((count_zero << 4)&0xf0)|need_write_length]); j++) {
            // printf("writing...%d\n",ac_cof[category[yuv_id]][count_zero][index][j] - '0');
            write_one_bit(fp,ac_cof[category[yuv_id]][((count_zero << 4)&0xf0)|need_write_length][j] - '0',false);
        }
        codeword_encode(fp,value,need_write_length);
        count_zero = 0;
    }
    if (padding_start_index!=63) {
        // EOB 0x00
        // fwrite(&c,1,1,fp);
        // printf("eob%s\n----\n",ac_cof[category[yuv_id]][0x00]);
        for (int i = 0; i<strlen((char*)ac_cof[category[yuv_id]][0x00]); i++) {
            write_one_bit(fp,ac_cof[category[yuv_id]][0x00][i] - '0',false);
        }
    }
}

void calculate_mcu(FILE *fp,int block_x, int block_y)
{

    rgb_t colour;

    // JPEG 裡是對每 8x8 個點為一個單位處理的.
    // 所以如果原始圖片的長寬不是 8 的倍數, 都需要先補成 8 的倍數, 好一塊塊的處理.
    // 另外, 記得剛才我說的 Cr Cb 都是 2x2 記錄一次嗎?
    // 所以大多數情況, 是要補成 16x16 的整數塊.按從左到右, 從上到下的次序排列
    //  (和我們寫字的次序一樣). JPEG 裡是對 Y Cr Cb 分別做 DCT 變換的. 這裡進行 DCT 變換
    // 的 Y, Cr, Cb 值的範圍都是 -128~127. (Y 被減去 128)
    //要char才可保留正負= =...
    char mcu[8][8];
    byte R,G,B;
    for (int yuv_id = 0; yuv_id<3; yuv_id++) {
        for (int y=0; y<8; y++) {
            for (int x=0; x<8; x++) {
                // JPEG 裡是對每 8x8 個點為一個單位處理的.
                // 所以如果原始圖片的長寬不是 8 的倍數, 都需要先補成 8 的倍數, 好一塊塊的處理.
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
                        // mcu[y][x] = yuv_image[(block_x + x) * image_height + (block_y + y)].y;
                        // printf("for %d %d %d %d choose value = %d %g\n",x,y,block_x,block_y,mcu[y][x],0.299 * colour.red + 0.587 * colour.green + 0.114 * colour.blue);
                        break;
                    case 1:
                        mcu[y][x] = (byte)(-0.1687 * colour.red - 0.3313 * colour.green + 0.5 * colour.blue);

                        // mcu[y][x] = yuv_image[(block_x + x) * image_height + (block_y + y)].cb;
                        break;
                    case 2:
                        mcu[y][x] = (byte)(0.5 * colour.red - 0.4187 * colour.green - 0.0813 * colour.blue);

                        // mcu[y][x] = yuv_image[(block_x + x) * image_height + (block_y + y)].cr;
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

    // yuv_element* yuv_image = (yuv_element*) malloc(image_height*image_width*sizeof(yuv_element));

    byte c;
    //Start of image
    c = 0xff;
    fwrite(&c, 1, 1, fp);
    c = SOI;
    fwrite(&c, 1, 1, fp);
    //APP0 commment
    c = 0xff;
    fwrite(&c, 1, 1, fp);
    c = APP0;
    fwrite(&c, 1, 1, fp);
    //APP0 length
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

    //DQT1
    c = 0xff;
    fwrite(&c, 1, 1, fp);
    c = DQT;
    fwrite(&c, 1, 1, fp);
    //length
    write_word_to_littleendian(fp, 132);
    //id & precision 0x00
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
            // printf("%u ",zigzag_quantize_table[0][i*8+j]);
            fwrite(&(zigzag_quantize_table[0][i*8+j]),1,1,fp);
        }
        // printf("\n");
    }
    //id & precision 0x01
    c = 0x01;
    fwrite(&c, 1, 1, fp);
    for (int i = 0; i<8; i++) {
        for (int j = 0; j<8; j++) {
            fwrite(&(zigzag_quantize_table[1][i*8+j]),1,1,fp);
        }
    }
    //SOF
    c = 0xff;
    fwrite(&c, 1, 1, fp);
    c = SOF0;
    fwrite(&c, 1, 1, fp);
    //length
    write_word_to_littleendian(fp, 17);
    //   - 數據精度 (1 byte) 每個樣本位數, 通常是 8 (大多數軟件不支持 12 和 16)
    c = 8;
    fwrite(&c, 1, 1, fp);
    write_word_to_littleendian(fp,(word)image_height);
    write_word_to_littleendian(fp,(word)image_width);
    // # of component = 3
    c = 3;
    fwrite(&c, 1, 1, fp);
    // 1 = Y, 2 = Cr, 3 = Cb
    c = 1;
    fwrite(&c, 1, 1, fp);
    //      - 採樣係數 (bit 0-3 vert., 4-7 hor.)
    c = 0x11;
    fwrite(&c, 1, 1, fp);
    c = 0x00;
    fwrite(&c, 1, 1, fp);
    // 1 = Y, 2 = Cr, 3 = Cb
    c = 2;
    fwrite(&c, 1, 1, fp);
    //      - 採樣係數 (bit 0-3 vert., 4-7 hor.)
    c = 0x11;
    fwrite(&c, 1, 1, fp);
    c = 1;
    fwrite(&c, 1, 1, fp);
    // 1 = Y, 2 = Cr, 3 = Cb
    c = 3;
    fwrite(&c, 1, 1, fp);
    //      - 採樣係數 (bit 0-3 vert., 4-7 hor.)
    c = 0x11;
    fwrite(&c, 1, 1, fp);
    c = 1;
    fwrite(&c, 1, 1, fp);

    //DHT
    c = 0xff;
    fwrite(&c, 1, 1, fp);
    c = DHT;
    fwrite(&c, 1, 1, fp);
    //DHT length
    write_word_to_littleendian(fp, 0x01a2);
    // 寫入spec指定的表(ref K.3.3)
    // 0x00表示DC直流0號表；
    c = 0x00;
    fwrite(&c, 1, 1, fp);
    fwrite(dc_cof_code_length[0], 1, 16, fp);
    fwrite(dc_cof_value[0], 1, 12, fp);
    // 0x01表示DC直流1號表；
    c = 0x01;
    fwrite(&c, 1, 1, fp);
    fwrite(dc_cof_code_length[1], 1, 16, fp);
    fwrite(dc_cof_value[1], 1, 12, fp);
    // 0x10表示AC交流0號表；
    c = 0x10;
    fwrite(&c, 1, 1, fp);
    fwrite(ac_cof_code_length[0], 1, 16, fp);
    fwrite(ac_cof_value[0], 1, 162, fp);
    // 0x11表示AC交流1號表。
    c = 0x11;
    fwrite(&c, 1, 1, fp);
    fwrite(ac_cof_code_length[1], 1, 16, fp);
    fwrite(ac_cof_value[1], 1, 162, fp);


    //SOS
    c = 0xff;
    fwrite(&c, 1, 1, fp);
    c = SOS;
    fwrite(&c, 1, 1, fp);
    //SOS length
    write_word_to_littleendian(fp, 12);
    //component num should be 3
    c = 3;
    fwrite(&c, 1, 1, fp);
    //component_id 1 = Y, 2 = Cr, 3 = Cb
    c = 1;
    fwrite(&c, 1, 1, fp);
    //0 ac & 0 dc table
    c = 0x00;
    fwrite(&c, 1, 1, fp);
    //component_id 1 = Y, 2 = Cr, 3 = Cb
    c = 2;
    fwrite(&c, 1, 1, fp);
    //1 ac & 1 dc table
    c = 0x11;
    fwrite(&c, 1, 1, fp);
    //component_id 1 = Y, 2 = Cr, 3 = Cb
    c = 3;
    fwrite(&c, 1, 1, fp);
    //1 ac & 1 dc table
    c = 0x11;
    fwrite(&c, 1, 1, fp);
    //忽略三個byte (???) 超謎= =
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
    //記得清掃
    write_one_bit(fp,0,true);
    c = 0xff;
    fwrite(&c, 1, 1, fp);
    c = EOI;
    fwrite(&c, 1, 1, fp);
    // 壓縮算法簡介
    // 1. 色彩模型
    // 2. DCT (離散餘弦變換)
    // 3. 重排列 DCT 結果
    // 4. 量化
    // 5. 0 RLE 編碼
    // 6. 範式 Huffman 編碼
    // 7. DC 的編碼
    // free(yuv_image);
    fclose(fp);
    return 0;
}
