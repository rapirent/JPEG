#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include "util.h"
// #include "./bitmap/bitmap_image.hpp"

#define min(a,b) (a<b?a:b)
#define max(a,b) (a>b?a:b)
#define abs(a) (a<0?-a:a)
#define clip(a) max(min(a,255),0)



int image_height,image_width;
byte category[3] = {0,1,1}; // table for luminance or chrominance
byte zigzag_quantize_table[2][64];

// Implementation of LLM DCT.
// void llm_dct(const double in[8], double out[8]) {
//     // Constants:
//     const double s1 = sin(1. * M_PI / 16.);
//     const double c1 = cos(1. * M_PI / 16.);
//     const double s3 = sin(3. * M_PI / 16.);
//     const double c3 = cos(3. * M_PI / 16.);
//     const double r2s6 = sqrt(2.) * sin(6. * M_PI / 16.);
//     const double r2c6 = sqrt(2.) * cos(6. * M_PI / 16.);

//     // After stage 1:
//     const double s1_0 =  in[0] + in[7];
//     const double s1_1 =  in[1] + in[6];
//     const double s1_2 =  in[2] + in[5];
//     const double s1_3 =  in[3] + in[4];
//     const double s1_4 = -in[4] + in[3];
//     const double s1_5 = -in[5] + in[2];
//     const double s1_6 = -in[6] + in[1];
//     const double s1_7 = -in[7] + in[0];

//     // After stage 2:
//     const double s2_0 =  s1_0 + s1_3;
//     const double s2_1 =  s1_1 + s1_2;
//     const double s2_2 = -s1_2 + s1_1;
//     const double s2_3 = -s1_3 + s1_0;

//     const double z1 = c3 * (s1_7 + s1_4);
//     const double s2_4 = ( s3-c3) * s1_7 + z1;
//     const double s2_7 = (-s3-c3) * s1_4 + z1;

//     const double z2 = c1 * (s1_6 + s1_5);
//     const double s2_5 = ( s1-c1) * s1_6 + z2;
//     const double s2_6 = (-s1-c1) * s1_5 + z2;

//     // After stage 3:
//     const double s3_0 =  s2_0 + s2_1;
//     const double s3_1 = -s2_1 + s2_0;

//     const double z3 = r2c6 * (s2_3 + s2_2);
//     const double s3_2 = ( r2s6-r2c6) * s2_3 + z3;
//     const double s3_3 = (-r2s6-r2c6) * s2_2 + z3;

//     const double s3_4 =  s2_4 + s2_6;
//     const double s3_5 = -s2_5 + s2_7;
//     const double s3_6 = -s2_6 + s2_4;
//     const double s3_7 =  s2_7 + s2_5;

//     // After stage 4:
//     const double s4_4 = -s3_4 + s3_7;
//     const double s4_5 =  s3_5 * sqrt(2.);
//     const double s4_6 =  s3_6 * sqrt(2.);
//     const double s4_7 =  s3_7 + s3_4;

//     // Shuffle and scaling:
//     out[0] = s3_0 / sqrt(8.);
//     out[4] = s3_1 / sqrt(8.);
//     out[2] = s3_2 / sqrt(8.);
//     out[6] = s3_3 / sqrt(8.);
//     out[7] = s4_4 / sqrt(8.);
//     out[3] = s4_5 / sqrt(8.);  // Alternative: s3_5 / 2
//     out[5] = s4_6 / sqrt(8.);
//     out[1] = s4_7 / sqrt(8.);
// }

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
        for (int i = count; i<7;i++) {
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
    for (;tmp ;tmp = tmp >> 1) {
        need_write_length++;
    }
    //add leading 1 if negative
    *value = (*value) > 0 ? (*value) : ((0x1 << need_write_length) + (*value) -1);
    return need_write_length;
}
double PI = acos(-1);
void calculate_mcu_block(FILE *fp,int block[][8],int yuv_id) {
    //TODO
    double tmp[8][8], s[8][8];
    // printf("Start fdct!!!\n");
    int quantized_block[64];
    memset(tmp,0,sizeof(tmp));
    memset(s,0,sizeof(s));
    // for (int ii = 0; ii < 8; ii++) {
    //     for (int x = 0; x < 8; x++) {
    //         for (int y = 0; y < 8; y++) {
    //             s[ii][x] += block[ii][y] * cos_cache[(y*2 + 1) * x];
    //         }
    //         s[ii][x] = c(x) *s[ii][x] / 2.0;
    //     }
    // }
    // for (int ii = 0; ii < 8; ii++) {
    //     for (int jj = 0; jj < 8; jj++) {
    //         for (int x = 0; x < 8; x++) {
    //             tmp[ii][jj] += s[x][jj] * cos_cache[(x*2 + 1) * ii];
    //         }
    //         block[ii][jj] = c(ii) * tmp[ii][jj] / 2.0;
    //         // printf("%d ",(int)block[ii][jj]);
    //     }
    //     // printf("\n");
    // }
    // // printf("\n---------\n");
	for(int i=0;i<8;i++)
		for(int j=0;j<8;j++){
			double sum=0;
			for(int k=0;k<8;k++)
				sum+=(double)block[i][k]*cos((2*k+1)*j*PI/16);
			tmp[i][j] = (j==0 ? (1.0/sqrt(2)):1.0)*sum/2;
		}
	for(int i=0;i<8;i++) {
		for(int j=0;j<8;j++){
			double sum=0;
			for(int k=0;k<8;k++)
				sum+=tmp[k][j]*cos((2*k+1)*i*PI/16);
			block[i][j] = round((i==0 ? (1.0/sqrt(2)):1.0)*sum/2);
            // printf("%d ",j);
            // printf("%d ",block[i][j]);
		}
        // printf("\n");
    }
    // printf("\n---------\n");
    // printf("zigzag arrange start!!!\n");
    for (int x = 0; x < 8; x++) {
        for (int y = 0; y < 8; y++) {
            // printf("%d %d\n",x,y);
            quantized_block[zigzag_index[x][y]] = block[x][y];
        }
    }

    // 4. 量化
    //     對於前面得到的 64 個空間頻率振幅值, 我們將對它們作幅度分層量化操作.方
    // 法就是分別除以量化表裡對應值並四捨五入.
    // printf("quantization start!!!\n");
    for (int i = 0 ; i<64; i++ ) {
        quantized_block[i] = round((double)quantized_block[i]/zigzag_quantize_table[category[yuv_id]][i]);
    }
    // printf("gogogo!\n");
    // DC即一塊圖像樣本的平均值. 就是說, 它包含了原始 8x8 圖像塊裡的很多能量. (通常
    // 會得到一個很大的數值)
    static int dc_block[5] = {0,0,0,0,0};
    // JPEG 的作者指出連續塊的 DC 率之間有很緊密的聯繫,  因此他們決定對 8x8 塊的
    // DC 值的差別進行編碼. (Y, Cb, Cr 分別有自己的 DC)
    // Diff = DC(i)  - DC(i-1)
    // 所以這一塊的 DC(i) 就是:  DC(i)  = DC(i-1)  + Diff
    int diff = quantized_block[0] - dc_block[yuv_id];
    // printf("dc encode\n");
    int need_write_length;
    int index;
    dc_block[yuv_id] = dc_block[yuv_id] + diff;
    need_write_length = find_length(&diff);
    printf("length = %d, number = %d\n",need_write_length,diff);
    printf("%s\n----\n",dc_diff[category[yuv_id]][need_write_length]);
    for (int i = 0; i<strlen((char*)dc_diff[category[yuv_id]][need_write_length]);i++) {
        write_one_bit(fp, dc_diff[category[yuv_id]][need_write_length][i] - '0',false);
    }
    codeword_encode(fp,diff,need_write_length);
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
            for (int j = 0;j<count_zero/16;j++) {
                // printf("???");
                printf("%s\n----\n",ac_cof[category[yuv_id]][0xf0]);
                for (int k = 0; k<strlen((char*)ac_cof[category[yuv_id]][0xf0]);k++) {
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
        printf("%s\n----\n",ac_cof[category[yuv_id]][(count_zero << 4)&0xf0|need_write_length]);
        for (int j = 0; j<strlen((char*)ac_cof[category[yuv_id]][(count_zero << 4)&0xf0|need_write_length]);j++) {
            // printf("writing...%d\n",ac_cof[category[yuv_id]][count_zero][index][j] - '0');
            write_one_bit(fp,ac_cof[category[yuv_id]][(count_zero << 4)&0xf0|need_write_length][j] - '0',false);
        }
        codeword_encode(fp,value,need_write_length);
        count_zero = 0;
    }
    if (padding_start_index!=63) {
        // EOB 0x00
        // fwrite(&c,1,1,fp);
        printf("eob%s\n----\n",ac_cof[category[yuv_id]][0x00]);
        for (int i = 0; i<strlen((char*)ac_cof[category[yuv_id]][0x00]);i++) {
            write_one_bit(fp,ac_cof[category[yuv_id]][0x00][i] - '0',false);
        }
    }
}

void calculate_mcu(FILE *fp,int block_x, int block_y, yuv_element* yuv_image) {

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
                    // printf("for %d %d %d %d choose value = %d\n",x,y,block_x,block_y,yuv_image[(block_y - y) * image_width + (block_x + x)].y);
                    switch(yuv_id) {
                        case 0:
                            mcu[y][x] = yuv_image[(block_y - y) * image_width + (block_x + x)].y;
                            break;
                        case 1:
                            mcu[y][x] = yuv_image[(block_y - y) * image_width + (block_x + x)].cr;
                            break;
                        case 2:
                            mcu[y][x] = yuv_image[(block_y - y) * image_width + (block_x + x)].cb;
                            break;
                    }
                }
                else {
                    mcu[y][x] = 0;
                }
                mcu[y][x]-=128;
            }
        }
        calculate_mcu_block(fp,mcu,yuv_id); 
    }
    // for (int y=0;y<8;y++) {
    //     for (int x=0;x<8;x++) {
    //         // JPEG 裡是對每 8x8 個點為一個單位處理的.
    //         // 所以如果原始圖片的長寬不是 8 的倍數, 都需要先補成 8 的倍數, 好一塊塊的處理.
    //         if (block_y - y>=0 && block_x + x < image_width) {
    //             // printf("for %d %d %d %d choose value = %d\n",x,y,block_x,block_y,yuv_image[(block_y - y) * image_width + (block_x + x)].cr);
    //             mcu[y][x] = yuv_image[(block_y - y) * image_width + (block_x + x)].cr;
    //         }
    //         else {
    //             mcu[y][x] = 0;
    //         }
    //         mcu[y][x]-=128;
    //     }
    // }
    // calculate_mcu_block(fp,mcu,1);

    // for (int y=0;y<8;y++) {
    //     for (int x=0;x<8;x++) {
    //         // JPEG 裡是對每 8x8 個點為一個單位處理的.
    //         // 所以如果原始圖片的長寬不是 8 的倍數, 都需要先補成 8 的倍數, 好一塊塊的處理.
    //         if (block_y - y>=0 && block_x + x < image_width) {
    //             // printf("for %d %d %d %d choose value = %d\n",x,y,block_x,block_y,yuv_image[(block_y - y) * image_width + (block_x + x)].cb);
    //             mcu[y][x] = yuv_image[(block_y - y) * image_width + (block_x + x)].cb;
    //         }
    //         else {
    //             mcu[y][x] = 0;
    //         }
    //         mcu[y][x]-=128;
    //     }
    // }
    // calculate_mcu_block(fp,mcu,2);   
}

int main(int argc, char* argv[]) {
    if (argc != 2  && argc != 3) {
        printf("[ERROR]:\nusage: ./encoder <input_file> [<output_file_name>]");
        exit(1);
    }
    // bitmap_image image(argv[1]);
    FILE *bitmap_f;
    if ((bitmap_f=fopen(argv[1],"r"))==NULL) {
        printf("%s can't be opened\n", argv[1]);
        exit(1);
    }
    fseek(bitmap_f, 18, SEEK_CUR);
    fread(&image_width, 1, sizeof(int), bitmap_f);
    fread(&image_height, 1, sizeof(int), bitmap_f);
    if((image_width&7)!=0 || (image_height&7)!=0) {
        exit(1);
    }
	fseek(bitmap_f, 28, SEEK_CUR);
    yuv_element* yuv_image = (yuv_element*) malloc(image_height*image_width*sizeof(yuv_element));

    // word R,G,B;
    byte discard[3] = {0,0,0};
    byte R,G,B;
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
            // image.get_pixel(x, y, colour);
            fread(&R, 1, 1, bitmap_f);
			fread(&G, 1, 1, bitmap_f);
			fread(&B, 1, 1, bitmap_f);
            // yuv_image[y*image_width+x].y = 0.299*colour.red + 0.587*colour.green + 0.114*colour.blue;
            // yuv_image[y*image_width+x].cb =  - 0.1687*colour.red - 0.3313*colour.green + 0.5   *colour.blue + 128;
            // yuv_image[y*image_width+x].cr =    0.5   *colour.red - 0.4187*colour.green - 0.0813*colour.blue + 128;
            // yuv_image[y*image_width+x].y = 0.299*R + 0.587*G + 0.114*B;
            // yuv_image[y*image_width+x].cb =  - 0.1687*R - 0.3313*G + 0.5   *B + 128;
            // yuv_image[y*image_width+x].cr =    0.5   *R - 0.4187*G - 0.0813*B + 128;

            yuv_image[y*image_width+x].y = clip(round(0.2126*R + 0.7152*G + 0.0722*B));
            // printf("fucking value is %u\n",yuv_image[y*image_width+x].y);
            yuv_image[y*image_width+x].cb = clip(round(-0.09991*R - 0.33609*G + 0.436*B + 128));
            yuv_image[y*image_width+x].cr = clip(round(0.615*R -0.55861*G -0.05639*B + 128));
        }
        fread(discard,1,(4-(image_width*3)%4)%4,bitmap_f);
    }

    FILE* fp;
    if (argc==3 && (fp = fopen(argv[2], "w")) == NULL) {
        printf("%s can't be opened\n", argv[2]);
        exit(1);
    }
    else if ((fp = fopen("image.jpg", "w")) == NULL) {
        printf("image.jpeg can't be opened\n");
        exit(1);
    }

    // image_height = image.height();
    // image_width  = image.width();

    // JPEG 裡, 要對數據壓縮, 先要做一次 DCT 變換. DCT 變換的原理, 涉及到數學
    // 知識, 這裡我們不必深究. 反正和傅立葉變換(學過高數的都知道) 是差不多了. 經過
    // 這個變換, 就把圖片裡點和點間的規律呈現出來了, 更方便壓縮.JPEG 裡是對每 8x8
    // 個點為一個單位處理的. 所以如果原始圖片的長寬不是 8 的倍數, 都需要先補成 8
    // 的倍數, 好一塊塊的處理. 另外, 記得剛才我說的 Cr Cb 都是 2x2 記錄一次嗎? 所
    // 以大多數情況, 是要補成 16x16 的整數塊.按從左到右, 從上到下的次序排列 (和我
    // 們寫字的次序一樣). JPEG 裡是對 Y Cr Cb 分別做 DCT 變換的. 這裡進行 DCT 變換
    // 的 Y, Cr, Cb 值的範圍都是 -128~127. (Y 被減去 128)

    byte c;
    //Start of image
    c = 0xff;fwrite(&c, 1, 1, fp);
    c = SOI;fwrite(&c, 1, 1, fp);
    //APP0 commment
    c = 0xff;fwrite(&c, 1, 1, fp);
    c = APP0;fwrite(&c, 1, 1, fp);
    //APP0 length
    write_word_to_littleendian(fp, 16);
    char comment[5] = "JFIF";
    fwrite(comment, 1, 5, fp);
    c = 0x01;fwrite(&c, 1, 1, fp);
    c = 0x01;fwrite(&c, 1, 1, fp);
    c = 0x00;fwrite(&c, 1, 1, fp);
    c = 0x00;fwrite(&c, 1, 1, fp);
    c = 0x01;fwrite(&c, 1, 1, fp);
    c = 0x00;fwrite(&c, 1, 1, fp);
    c = 0x01;fwrite(&c, 1, 1, fp);
    c = 0x00;fwrite(&c, 1, 1, fp);
    c = 0x00;fwrite(&c, 1, 1, fp);

    //DQT1
    c = 0xff;fwrite(&c, 1, 1, fp);
    c = DQT;fwrite(&c, 1, 1, fp);
    //length
    write_word_to_littleendian(fp, 132);
    //id & precision 0x00
    c = 0x00;fwrite(&c, 1, 1, fp);
    for (int i = 0;i<8; i++) {
        for (int j = 0;j<8; j++) {
            zigzag_quantize_table[0][zigzag_index[i][j]] = quantize_table[0][i*8+j];
            zigzag_quantize_table[1][zigzag_index[i][j]] = quantize_table[1][i*8+j];
        }
    }
    for (int i = 0;i<8; i++) {
        for (int j = 0;j<8; j++) {
            // printf("%u ",zigzag_quantize_table[0][i*8+j]);
            fwrite(&(zigzag_quantize_table[0][i*8+j]),1,1,fp);
        }
        // printf("\n");
    }
    //id & precision 0x01
    c = 0x01;fwrite(&c, 1, 1, fp);
    for (int i = 0;i<8; i++) {
        for (int j = 0;j<8; j++) {
            fwrite(&(zigzag_quantize_table[1][i*8+j]),1,1,fp);
        }
    }
    //SOF
    c = 0xff;fwrite(&c, 1, 1, fp);
    c = SOF0;fwrite(&c, 1, 1, fp);
    //length
    write_word_to_littleendian(fp, 17);
    //   - 數據精度 (1 byte) 每個樣本位數, 通常是 8 (大多數軟件不支持 12 和 16)
    c = 8;fwrite(&c, 1, 1, fp);
    write_word_to_littleendian(fp,(word)image_height);
    write_word_to_littleendian(fp,(word)image_width);
    // # of component = 3
    c = 3;fwrite(&c, 1, 1, fp);
    // 1 = Y, 2 = Cr, 3 = Cb
    c = 1;fwrite(&c, 1, 1, fp);
    //      - 採樣係數 (bit 0-3 vert., 4-7 hor.)
    c = 0x11;fwrite(&c, 1, 1, fp);
    c = 0x00;fwrite(&c, 1, 1, fp);
    // 1 = Y, 2 = Cr, 3 = Cb
    c = 2;fwrite(&c, 1, 1, fp);
    //      - 採樣係數 (bit 0-3 vert., 4-7 hor.)
    c = 0x11;fwrite(&c, 1, 1, fp);
    c = 1;fwrite(&c, 1, 1, fp);
    // 1 = Y, 2 = Cr, 3 = Cb
    c = 3;fwrite(&c, 1, 1, fp);
    //      - 採樣係數 (bit 0-3 vert., 4-7 hor.)
    c = 0x11;fwrite(&c, 1, 1, fp);
    c = 1;fwrite(&c, 1, 1, fp);

    //DHT
    c = 0xff;fwrite(&c, 1, 1, fp);
    c = DHT;fwrite(&c, 1, 1, fp);
    //DHT length
    write_word_to_littleendian(fp, 0x01a2);
    // 寫入spec指定的表(ref K.3.3)
    // 0x00表示DC直流0號表；
    c = 0x00;fwrite(&c, 1, 1, fp);
    fwrite(dc_cof_code_length[0], 1, 16, fp);
    fwrite(dc_cof_value[0], 1, 12, fp);
    // 0x01表示DC直流1號表；
    c = 0x01;fwrite(&c, 1, 1, fp);
    fwrite(dc_cof_code_length[1], 1, 16, fp);
    fwrite(dc_cof_value[1], 1, 12, fp);
    // 0x10表示AC交流0號表；
    c = 0x10;fwrite(&c, 1, 1, fp);
    fwrite(ac_cof_code_length[0], 1, 16, fp);
    fwrite(ac_cof_value[0], 1, 162, fp);
    // 0x11表示AC交流1號表。
    c = 0x11;fwrite(&c, 1, 1, fp);
    fwrite(ac_cof_code_length[1], 1, 16, fp);
    fwrite(ac_cof_value[1], 1, 162, fp);


    //SOS
    c = 0xff;fwrite(&c, 1, 1, fp);
    c = SOS;fwrite(&c, 1, 1, fp);
    //SOS length
    write_word_to_littleendian(fp, 12);
    //component num should be 3
    c = 3;fwrite(&c, 1, 1, fp);
    //component_id 1 = Y, 2 = Cr, 3 = Cb
    c = 1;fwrite(&c, 1, 1, fp);
    //0 ac & 0 dc table
    c = 0x00;fwrite(&c, 1, 1, fp);
    //component_id 1 = Y, 2 = Cr, 3 = Cb
    c = 2;fwrite(&c, 1, 1, fp);
    //1 ac & 1 dc table
    c = 0x11;fwrite(&c, 1, 1, fp);
    //component_id 1 = Y, 2 = Cr, 3 = Cb
    c = 3;fwrite(&c, 1, 1, fp);
    //1 ac & 1 dc table
    c = 0x11;fwrite(&c, 1, 1, fp);
    //忽略三個byte (???) 超謎= =
    c = 0x00;fwrite(&c, 1, 1, fp);
    c = 0x3f;fwrite(&c, 1, 1, fp);
    c = 0x00;fwrite(&c, 1, 1, fp);

    init_cos_cache();
    for (int y = image_height - 1; y>=0;y-=8) {
        for (int x = 0;x<image_width;x+=8) {
            calculate_mcu(fp,x,y,yuv_image);
        }
    }
    //記得清掃
    write_one_bit(fp,0,true);
    c = 0xff;fwrite(&c, 1, 1, fp);
    c = EOI;fwrite(&c, 1, 1, fp);
    // 壓縮算法簡介
    // 1. 色彩模型
    // 2. DCT (離散餘弦變換)
    // 3. 重排列 DCT 結果
    // 4. 量化
    // 5. 0 RLE 編碼
    // 6. 範式 Huffman 編碼
    // 7. DC 的編碼
    free(yuv_image);
    return 0;
}
