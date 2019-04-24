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
            //AC 0
            return 2;
        case 0x11:
            //AC 1
            return 3;
        case 0x00:
            //DC 0
            return 0;
        case 0x01:
            //DC 1
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

int quantize_table_list[4][64];//取值範圍0~3
frame_data f0;
huffman_table_list huffman_tables[4];
byte component_mapping_huffman[5];
double block[8][8];

void read_qt(FILE* fp)
{
    byte len = read_word_to_bigendian(fp);
    byte id_prec;
    printf("section length = %x (hex) %u (dec)\n",len,len);
    len = len - 2;
    byte pq,id;
    byte c;
    word tmp;
    while (len > 0) {
        fread(&id_prec,1,1,fp);
        len-=1;
        //       - QT 信息 (1 byte):
        //  bit 4..7: QT 精度, 0 = 8 bit, 否则 16 bit
        pq = (id_prec >> 4) ? 1: 0;
        // 低4位：量化表ID，取值範圍為0～3
        //tq
        id = id_prec & 0x0f;
        printf("the qt id is %u, its precision is (%u + 1) * 8bits\n",id,pq);
        for(int i = 0; i<64; i++) {
            //  - n 字节的 QT (8x8 martrix), n = 64*(精度+1)
            if (pq) {
                tmp = read_word_to_bigendian(fp);
                len = len - 2;
            } else {
                fread(&c, 1, 1, fp);
                len--;
                tmp = c;
            }
            quantize_table_list[id][i] = tmp;
            printf("%u ",tmp);
            if((i+1)%8==0) {
                printf("\n");
            }
        }
        // len-=64*(precision+1);
    }
}

void read_frame(FILE *fp)
{
    //长度 (高字节, 低字节), 8+components*3
    int len = read_word_to_bigendian(fp);
    printf("section length = %x (hex) %u (dec)\n",len,len);

    //   - 数据精度 (1 byte) 每个样本位数, 通常是 8 (大多数软件不支持 12 和 16)
    fread(&(f0.precision),1,1,fp);
    printf("precision is %d(dec) \n",f0.precision);
    // Number of lines – Specifies the maximum number of lines in the source image.
    //   - 图片高度 (高字节, 低字节), 如果不支持 DNL 就必须 >0
    f0.height = read_word_to_bigendian(fp);
    // Number of samples per line
    //   - 图片宽度 (高字节, 低字节), 如果不支持 DNL 就必须 >0
    f0.width = read_word_to_bigendian(fp);
    printf("image height %d(dec), image widht %d(dec)\n",f0.height,f0.width);
    //   - components 数量(1 byte), 灰度图是 1, YCbCr/YIQ 彩色图是 3, CMYK 彩色图是 4
    fread(&(f0.components_num),1,1,fp);
    printf("# of image components is %d (dec)\n",f0.components_num);

    byte component_id;
    byte sample;
    f0.vmax = 0x00;
    f0.hmax = 0x00;
    for (int i = 0; i < f0.components_num; i++) {
        //   - 每个 component: 3 bytes
        //      - component id (1 = Y, 2 = Cb, 3 = Cr, 4 = I, 5 = Q)
        fread(&component_id,1,1,fp);
        printf("------------\nnow is reading component whose id=%d(dec)\n",
               component_id);
        //      - 采样系数 (bit 0-3 vert., 4-7 hor.)
        fread(&sample,1,1,fp);
        // 先Horizontal sampling factor 再Vertical sampling factor
        (f0.frame_components[component_id-1]).horizontal_sample = (sample >> 4) & 0x0f;
        (f0.frame_components[component_id-1]).vertical_sample = sample & 0x0f;
        printf("horizontal sample factor=%d,\n",
               (f0.frame_components[component_id-1]).horizontal_sample);
        printf("horizontal sample factor=%d\n",
               (f0.frame_components[component_id-1]).vertical_sample);
        //MCU 的寬 = 8 * 最高水平採樣率
        // MCU 的高 = 8 * 最高垂直採樣率
        if( (f0.frame_components[component_id-1]).horizontal_sample > f0.hmax) {
            f0.hmax = (f0.frame_components[component_id-1]).horizontal_sample;
        }
        if( (f0.frame_components[component_id-1]).vertical_sample > f0.vmax) {
            f0.vmax = (f0.frame_components[component_id-1]).vertical_sample;
        }
        //Quantization table destination selector
        fread(&((f0.frame_components[component_id-1]).qantize_table_id),1,1,fp);
        printf("qt destination = %d\n",
               (f0.frame_components[component_id-1]).qantize_table_id);
    }
}

void read_ht(FILE* fp)
{
    byte huffman_length[16];
    int len = read_word_to_bigendian(fp);
    printf("section length = %x (hex) %u (dec)\n",len,len);
    len = len - 2;
    int ht_node_num;
    byte ht_length;
    byte class_id;
    byte ht_value;
    while (len > 0) {
        fread(&class_id,1,1,fp);
        len--;
        //      bit 0..3: HT 號 (0..3, 否則錯誤)
        //      bit 4   : HT 類型, 0 = DC table, 1 = AC table
        //      bit 5..7: 必須是 0
        // 0x00表示DC直流0號表；
        // 0x01表示DC直流1號表；
        // 0x10表示AC交流0號表；
        // 0x11表示AC交流1號表。
        printf("Huffman Table class & destination %.2x\n",class_id);
        // （2~17字節）為不同位數的碼字的數量。
        //這16個數值實際意義為：沒有1位和4位的哈夫曼碼字；2位和3位的碼字各有2個；5位碼字有5個；6位和8位碼字各有1個；7位碼字各有6個；沒有9位或以上的碼字。
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
        printf("----------\ndecoding (total have %d node)...\n",ht_node_num);

        //綠色部分（18~34字節）為編碼內容。
        //由藍色部分數據知道，此哈夫曼樹有0+2+2+0+5+1+6+1=17個葉子結點，即本字段應該有17個字節。
        //這段數據表示17個葉子結點按從小到大排列，其權值依次為0、1、11、2、21、3、31、41……
        //Value associated with each Huffman code –
        //Specifies, for each i, the value associated with each Huffman code of length i. The meaning of each value is determined by the Huffman coding model. The Vi,j’s are the elements of the list HUFFVAL.
        for(int height=0, leaf_index = 0; height<16; height++) {
            //huffman code length = huffman tree height
            printf("reading the %d bits length codeword (%.3d total)\n", height+1,
                   huffman_length[height]);
            for(int j = 0; j<huffman_length[height]; j++) {

                fread(&(ht_value),1,1,fp);
                len--;
                ht_leaf[leaf_index].value = ht_value;
                // 若高度相等：leaf[n] = leaf[n — 1] + 1
                // 若高度差 1 ：leaf[n] = (leaf[n — 1] + 1) * 2
                // 若高度差 k ：leaf[n] = (leaf[n — 1] + 1) * 2^k
                //int foo = 0b1010;
                //https://medium.com/@yc1043/%E8%B7%9F%E6%88%91%E5%AF%AB-jpeg-%E8%A7%A3%E7%A2%BC%E5%99%A8-%E4%B8%89-%E8%AE%80%E5%8F%96%E9%87%8F%E5%8C%96%E8%A1%A8-%E9%9C%8D%E5%A4%AB%E6%9B%BC%E8%A1%A8-585f2cf4c494
                // ht_leaf[leaf_index].codeword = ((int)pow(2,
                // height - ht_leaf[leaf_index-1].codeword_len)) * (codeword);
                ht_leaf[leaf_index].codeword = (codeword++);
                ht_leaf[leaf_index].codeword_len = height+1;
                printf("leaf #%d | codword %.2x | source value %.2x | length=%d \n",leaf_index,
                       ht_leaf[leaf_index].codeword, ht_value,ht_leaf[leaf_index].codeword_len);
                leaf_index++;
            }
            //不斷往左移一位
            codeword = codeword<<1;
        }
        huffman_tables[huffman_table_index(class_id)].start = ht_leaf;
    }
}

void read_sos(FILE* fp)
{
    //     SOS: Start Of Scan:
    //   - 長度 (高字節, 低字節), 必須是 6+2*(掃瞄行內組件的數量)
    int len = read_word_to_bigendian(fp);
    printf("section length = %x (hex) %u (dec)\n",len,len);
    //   - 掃瞄行內組件的數量 (1 byte), 必須 >= 1 , <=4 (否則是錯的) 通常是 3
    byte component_num;
    fread(&component_num,1,1,fp);
    printf(" Number of img components = %d (dec)\n",component_num);
    assert(len == 6+2*component_num);
    len = len - 3;
    assert(component_num==3);//JFIF defined
    //   - 每個組件: 2 bytes
    //      - component id (1 = Y, 2 = Cb, 3 = Cr, 4 = I, 5 = Q), 見 SOF0
    byte component_id,destination;
    for (int i =0; i<component_num; i++) {

        fread(&component_id,1,1,fp);
        fread(&destination,1,1,fp);
        len = len-2;
        printf("# %d component select %.2x (DC) and %.2x(AC)  destination\n",component_id,destination>>4,
               (destination&0x0f)|0x10);
        component_mapping_huffman[component_id-1] = destination;
    }
    //   - 忽略 3 bytes (???)
    fread(&destination,1,1,fp);
    printf("discrad byte = %.2x\n",destination);
    fread(&destination,1,1,fp);
    printf("discrad byte = %.2x\n",destination);
    fread(&destination,1,1,fp);
    printf("discrad byte = %.2x\n",destination);
}
//only read one bit
byte read_one_bit(FILE *fp)
{
    static byte buffer;
    static byte count = 0;
    byte check_ff00;
    //每八一次
    if (!count) {
        fread(&buffer,1,1,fp);
        // printf("-------\nread this buffer=%.2x\n",buffer);
        //https://www.jianshu.com/p/ccb52e9cd2e4
        //由於JPEG中以0XFF來做為特殊標記符，
        //因此，如果某個像素的取值為0XFF，
        //那麼實際在保存的時候，是以0XFF00來保存的，
        //從而避免其跟特殊標記符0XFF之間產生混淆。
        //在讀取文件信息的時候，如果遇0XFF00，就必須去除後面的00；即，將0XFF00當做0XFF；
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
    //每八一次
    // printf("------\nread bit %d\n",bit);
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
        //正數就是照一般二進位計算，負數就是對正數的碼字取反而已
    }
    decoding_code = leading ? decoding_code: -decoding_code;
    printf("decoding = %d\n",decoding_code);
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
        // printf("now codeword is %.2x(hex) |"BYTE_TO_BINARY_PATTERN,codeword,
        //        BYTE_TO_BINARY(codeword >> 8));
        // printf(BYTE_TO_BINARY_PATTERN"(binary) with length = %d \n",
        //        BYTE_TO_BINARY(codeword),i+1);
        for(int j = 0; j<ht_node_num; j++) {
            // printf("the huffman=%.2x(hex) len=%d \n",huffman_table[j].codeword,huffman_table[j].codeword_len);
            //i+1 = 位移幾次 + 1 = 有幾位
            if(huffman_table[j].codeword==codeword&&huffman_table[j].codeword_len==(i+1)) {
                // printf("find codeword= %.2x (len = %d) in huffman tree\n",codeword,i+1);
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
// void calcualte_mcu_block(FILE *fp,byte component_id) {
    //https://github.com/MROS/jpeg_tutorial/blob/master/doc/%E8%B7%9F%E6%88%91%E5%AF%ABjpeg%E8%A7%A3%E7%A2%BC%E5%99%A8%EF%BC%88%E5%9B%9B%EF%BC%89%E8%AE%80%E5%8F%96%E5%A3%93%E7%B8%AE%E5%9C%96%E5%83%8F%E6%95%B8%E6%93%9A.md
    // 每個 block 都是 8 * 8 ，最左上角的數值就是直流變量，要使用直流霍夫曼表來解碼；
    //而餘下的 63 個數值是交流變量，使用交流霍夫曼表來解碼。
    // DC:不斷從數據流中讀取一個 bit，直到可以對上直流霍夫曼表中的一個碼字，
    // 取出的對應信源編碼代表的是接下來還要讀入幾位，
    // 假如是 n 就繼續讀取 n bits ，以下表解碼後，就是直流變量。
    //  - bit 0..3: AC table (0..3)
    //  - bit 4..7: DC table (0..3)
    // 前一個代表直流(0)或交流(1)
    // 0000 0000 | 0000 0001 | 0001 0000 | 0001 0001
    //在整個圖片解碼的開始, 你需要先初始化 DC 值為 0.
    static int dc_block[5] = {0,0,0,0,0};
    memset(block,0,sizeof(block));
    //DC在高位，必須轉為0x00或0x01
    printf("pickup DC table %.2x\n",
           (component_mapping_huffman[component_id] >> 4)&0x0f);
    int dc_table_id = huffman_table_index((component_mapping_huffman[component_id]>> 4)&0x0f);
    //從此顏色份量單元數據流的起點開始一位一位的讀入，直到讀入的編碼與該份量直流哈夫曼樹的某個碼字（葉子結點）一致，
    //然後用直流哈夫曼樹查得該碼字對應的權值。
    //權值（共8位）表示該直流份量數值的二進制位數，也就是接下來需要讀入的位數
    byte need_read_bit = decode_huffman(fp,dc_table_id);
    byte leading_bit;
    if (need_read_bit) {
        //譯碼
        //讀入數據流並對照直流哈夫曼樹，第一個哈夫曼編碼為110，其權值為6，
        //所以往後讀入6位數據“1001101”，
        //譯碼成數值為77。因為每個顏色份量單元只有一個直流份量，所以下一個就是第一個交流份量了。
        //c) 取得 N 位, 計算 Diff 值
        //  d) DC + = Diff
        dc_block[component_id] += codeword_decode(fp,need_read_bit);
    }
    printf("DC completed!\n");
    //  e) 寫入 DC 值:      " vector[0]=DC "
    block[0][0] = dc_block[component_id];
    printf("now block[0][0] = %lf\n",block[0][0]);
    printf("pickup AC table %.2x\n",(component_mapping_huffman[component_id]&0x0f)|0x10);
    int ac_table_id = huffman_table_index((component_mapping_huffman[component_id]
                                           &0x0f)|0x10);
    //------- 循環處理每個 AC 直到 EOB 或者處理到 64 個 AC
    byte zerosnum_needread,zerosnum;
    for (int i = 1; i< 64;) {
        printf("this is # %d AC block\n",i);
        zerosnum_needread = decode_huffman(fp,ac_table_id);
        //權值的高4位表示當前數值前面有多少個連續的零，低4位表示該交流份量數值的二進制位數，也就是接下來需要讀入的位數。
        // b) Huffman 解碼, 得到 (前面 0 數量, 組號)
        // [記住: 如果是(0,0) 就是 EOB 了]
        if(zerosnum_needread == 0x00) {
            printf("0x00 EOB!\n");
            break;
        }
        if(zerosnum_needread == 0xf0) {
            printf("0xf0 16bits");
            i+=16;
            continue;
        }
        need_read_bit = zerosnum_needread & 0x0f;
        zerosnum = (zerosnum_needread >> 4) & 0x0f;
        printf("there are %u leading zero, %u need to read\n",need_read_bit,zerosnum);
        for (byte j = 0; j< zerosnum; j++) {
            block[i/8][i%8] = 0.0;
            i++;
        }
        if (need_read_bit) {
            //c) 取得 N 位(組號) 計算 AC
            //d) 寫入相應數量的 0
            block[i/8][i%8] = codeword_decode(fp,need_read_bit);
            i++;
        }
    }
    for (int x = 0; x < 8; x++) {
        for (int y = 0; y < 8; y++) {
            printf("%g ",block[x][y]);
        }
        printf("\n");
    }
    printf("\n");
    // 1) 反量化 64 個矢量 : "for (i=0;i<=63;i++) vector[i]*=quant[i]" (注意防止溢出)
    printf("dequantize start!\n");
    printf("use quantize table id = %d\n quantize =\n",(f0.frame_components[component_id]).qantize_table_id);
    for (int quantize_i = 0; quantize_i<64; quantize_i++) {
        block[quantize_i/8][quantize_i%8] *=
            quantize_table_list[(f0.frame_components[component_id]).qantize_table_id][quantize_i];
    }
    //4) 將所有的 8bit 數加上 128
    //5) 轉換 YCbCr 到 RGB
    // printf("\nlocation = %8X\n", (unsigned int) &block);
    // printf("\nlocation2 = %8X\n", (unsigned int) &(block[5][4]));
    // 2) 重排列 64 個矢量到 8x8 的塊中
    printf("zigzag rearrange start!\n");
    double tmp[8][8];
    for (int x = 0; x < 8; x++) {
        for (int y = 0; y < 8; y++) {
            tmp[x][y] = block[zigzag_index[x][y]/8][zigzag_index[x][y]%8];
        }
    }
    // 3) 對 8x8 的塊作 IDCT
    //TODO: find AA&N?????
    double s[8][8];
    memset(s,0,sizeof(s));
    printf("IDCT start!!!!\n");
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
void calculate_mcu(FILE* fp)
{
    //一個 (Hmax*8,Vmax*8) 的塊被稱作 MCU (Minimun Coded Unix) 前面例子中一個
    // MCU = YDU,YDU,YDU,YDU,CbDU,CrDU

    // 如果  HY =1, VY=1
    //       HCb=1, VCb=1
    //       HCr=1, VCr=1
    // 這樣 (Hmax=1,Vmax=1), MCU 只有 8x8 大, MCU = YDU,CbDU,CrDU
    //MCU 的寬 = 8 * 最高水平採樣率
    int mcu_width = 8 * f0.hmax;
    // MCU 的高 = 8 * 最高垂直採樣率
    int mcu_height = 8 * f0.vmax;
    //算出有幾個MCU
    int mcus_on_x = (f0.width - 1) / mcu_width + 1;
    int mcus_on_y = (f0.height - 1) / mcu_height + 1;
    // Add access it by arr[i*M + j]
    mcu_small_block* data_unit = (mcu_small_block*) malloc(mcus_on_y*mcus_on_x*5*mcu_width*mcu_height*sizeof(
                                     mcu_small_block));

    // double (******data_unit) = (double*******) malloc(sizeof(double******)*mcus_on_y+)
    //MCU 的順序:由左到右，再由上到下
    for (int i  = 0; i<mcus_on_y; i++) {
        for (int j = 0; j<mcus_on_x; j++) {
            //MCU[i][j].Cb = Cb[new_i / 8][new_j / 8][new_i % 8][new_j % 8]
            for (int component_id = 0; component_id<f0.components_num; component_id++)  {
                // double** mcu = (double**) malloc(f0.hmax*f0.vmax*sizeof(double**));
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
                    //5) 轉換 YCbCr 到 RGB
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
    outimg.save_image("image.bmp");
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
                printf("\nread DQT\n");
                read_qt(fp);
                break;
            case SOF0:
                printf("\nread SOF0\n");
                read_frame(fp);
                break;
            case DHT:
                printf("\nread DHT\n");
                read_ht(fp);
                break;
            case SOS:
                printf("\nread SOS\n");
                read_sos(fp);
                calculate_mcu(fp);
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
                printf("\nread APP: discard\n");
                break;
            case EOI:
                printf("End of Image\n");
                b_EOI = true;
                break;
            }
        }
    }
    assert(b_EOI==true);
    return 0;
}
