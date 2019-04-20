#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include "util.h"

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
        exit(1);
        printf("huffman table index error");
    }
}

double cos_cache[200];

void init_cos_cache()
{
    for (int i = 0; i < 200; i++) {
        cos_cache[i] = cos(i * M_PI / 16.0);
    }
}

// 查阅标记SOF0，可以得到图像不同颜色分量的采样因子，
// 即Y、Cr、Cb三个分量各自的水平采样因子和垂直采样因子。
// 大多图片的采样因子为4：1：1或1：1：1。
// 其中，4：1：1即（2*2）：（1*1）：（1*1））；1：1：1即（1*1）：（1*1）：（1*1）。
// 记三个分量中水平采样因子最大值为Hmax，垂直采样因子最大值为Vmax，
// 那么单个MCU矩阵的宽就是Hmax*8像素，高就是Vmax*8像素。

int quantize_table_list[4][64];//取值範圍0~3
frame_data f0;
huffman_table_list huffman_tables[4];
byte component_mapping_huffman[5];

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
        // fscanf(fp,"%"SCNd8,&id_prec);
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
    // fscanf(fp,"%"SCNd8,&f0.precision);
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
    // fscanf(fp,"%"SCNd8,&f0.components_num);
    fread(&(f0.components_num),1,1,fp);
    printf("# of image components is %d (dec)\n",f0.components_num);

    byte component_id;
    byte sample;
    f0.vmax = 0x00;
    f0.hmax = 0x00;
    for (int i = 0; i < f0.components_num; i++) {
        //   - 每个 component: 3 bytes
        //      - component id (1 = Y, 2 = Cb, 3 = Cr, 4 = I, 5 = Q)
        // fscanf(fp,"%"SCNd8,&component_id);
        fread(&component_id,1,1,fp);
        printf("------------\nnow is reading component whose id=%d(dec)\n",
               component_id);
        //      - 采样系数 (bit 0-3 vert., 4-7 hor.)
        // fscanf(fp,"%"SCNd8,&sample);
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
        // fscanf(fp,"%"SCNd8,&((f0.frame_components[component_id-1]).qantize_table_id));
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
        // fscanf(fp,"%"SCNd8,&class_id);
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
            //huffman tree height = codeword length
            // if(fscanf(fp,"%"SCNd8,&huffman_length[i])!=1) {
            //     printf("read HT height info error");
            //     exit(1);
            // }
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
                // fscanf(fp,"%"SCNd8,&(ht_leaf[leaf_index].content));

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
    // fscanf(fp,"%"SCNd8,&component_num);
    fread(&component_num,1,1,fp);
    printf(" Number of img components = %d (dec)\n",component_num);
    assert(len == 6+2*component_num);
    len = len - 3;
    assert(component_num==3);//JFIF defined
    //   - 每個組件: 2 bytes
    //      - component id (1 = Y, 2 = Cb, 3 = Cr, 4 = I, 5 = Q), 見 SOF0
    byte component_id,destination;
    for (int i =0; i<component_num; i++) {

        // fscanf(fp,"%"SCNd8,&component_id);
        fread(&component_id,1,1,fp);
        // fscanf(fp,"%"SCNd8,&component_mapping_huffman[component_id]);
        fread(&destination,1,1,fp);
        len = len-2;
        printf("# %d component select %.2x (DC) and %.2x(AC)  destination\n",component_id,destination>>4,(destination&0x0f)|0x10);
        component_mapping_huffman[component_id-0x01] = destination;
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
        // fscanf(fp,"%"SCNd8,&buffer);
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
            // fscanf(fp,"%"SCNd8,&check_ff00);
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

byte decode_huffman(FILE *fp, byte huffman_table_id)
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
                printf("find codeword= %.2x (len = %d) in huffman tree\n",codeword,i+1);
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
mcu_small_block* calcualte_mcu_block(FILE *fp, byte component_id)
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
    double block[8][8];
    memset(block,0,sizeof(block));
    //DC在高位，必須轉為0x00或0x01
    printf("pickup DC table %.2x\n",
           (component_mapping_huffman[component_id] >> 4)&0x0f);
    byte dc_table_id = huffman_table_index((component_mapping_huffman[component_id]>> 4)&0x0f);
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
    byte ac_table_id = huffman_table_index((component_mapping_huffman[component_id]
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
    printf("gogogogo\n");
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            printf("%g ",(block[i][j]));
        }
        printf("\n");
    }
    // 1) 反量化 64 個矢量 : "for (i=0;i<=63;i++) vector[i]*=quant[i]" (注意防止溢出)
    printf("dequantize start!\n");
    for (int i = 0; i<64; i++) {
        block[i/8][i%8] *= quantize_table_list[component_id][i];
    }
    // 2) 重排列 64 個矢量到 8x8 的塊中
    printf("zigzag rearrange start!\n");
    double tmp[8][8];
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            tmp[i][j] = block[zigzag_index[i][j]/8][zigzag_index[i][j]%8];
        }
    }
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            block[i][j] = tmp[i][j];
        }
    }

    // 3) 對 8x8 的塊作 IDCT
    //TODO: find AA&N?????
    double s[8][8];
    memset(tmp,0,sizeof(tmp));
    memset(s,0,sizeof(s));
    printf("IDCT start!!!!\n");
    for (int j = 0; j < 8; j++) {
        for (int x = 0; x < 8; x++) {
            for (int y = 0; y < 8; y++) {
                s[j][x] += c (y) * block[x][y] * cos_cache[(j + j + 1) * y];
            }
            s[j][x] = s[j][x] / 2.0;
        }
    }
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            for (int x = 0; x < 8; x++) {
                tmp[i][j] += c(x) * s[j][x] * cos_cache[(i + i + 1) * x];
            }
            tmp[i][j] = tmp[i][j] / 2.0;
        }
    }
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            block[i][j] = tmp[i][j];
            // printf("%8x ",&(block[i][j]));
        }
        // printf("\n");
    }
    //4) 將所有的 8bit 數加上 128
    //5) 轉換 YCbCr 到 RGB
    // printf("\nlocation = %8X\n", (unsigned int) &block);
    // printf("\nlocation2 = %8X\n", (unsigned int) &(block[5][4]));
    return &block;
}

byte chomp(double x)
{
    if (x > 255.0) {
        return 255;
    } else if (x < 0) {
        return 0;
    } else {
        return (byte) x;
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
    int mcu_number_col = ceil((f0.height) / mcu_height);
    int mcu_number_row = ceil((f0.width) / mcu_width);
    // double image[mcu_height][mcu_width];

    rgb_element rbg_image[mcu_number_col][mcu_number_row];
    // double** mcu_block = (double**) malloc(8*8*sizeof(double**));
    // double** data_unit[5][mcu_width][mcu_height]; //Y Cb Cr
    //1 = Y, 2 = Cb, 3 = Cr, 4 = I, 5 = Q
    int block_index_i_Y = f0.frame_components[0].horizontal_sample / f0.hmax,
        block_index_j_Y = f0.frame_components[0].vertical_sample / f0.vmax,
        block_index_i_Cb = f0.frame_components[1].horizontal_sample / f0.hmax,
        block_index_j_Cb = f0.frame_components[1].vertical_sample / f0.vmax,
        block_index_i_Cr = f0.frame_components[2].horizontal_sample / f0.hmax,
        block_index_j_Cr = f0.frame_components[2].vertical_sample / f0.vmax;
    double data_unit[5][mcu_width][mcu_height][8][8];
    mcu_small_block* pseudo_block;
    //MCU 的順序:由左到右，再由上到下
    for (int i  = 0; i<mcu_number_col; i++) {
        for (int j = 0; j<mcu_number_row; j++) {
            //一個顏色分量內部各個 block 的順序:由左到右，再由上到下
            //# Cb 是一個四階陣列
            //# 前兩階描述 block 的位置，後兩階描述要擷取的是這 8*8 中的哪一個點
            //MCU[i][j].Cb = Cb[new_i / 8][new_j / 8][new_i % 8][new_j % 8]
            for (int qt_id = 0; qt_id<f0.components_num; qt_id++)  {
                // double** mcu = (double**) malloc(f0.hmax*f0.vmax*sizeof(double**));
                for (int a =0; a<f0.frame_components[qt_id].horizontal_sample; a++) {
                    for (int b = 0; b<f0.frame_components[qt_id].vertical_sample; b++) {
                        pseudo_block = calcualte_mcu_block(fp, qt_id);
                        for (int x = 0; x < 8; x++) {
                            for (int y = 0; y < 8; y++) {
                                data_unit[qt_id][a][b][x][y] = (*pseudo_block)[x][y];
                            }
                        }
                        // printf("\nlocation = %8X\n", (unsigned int) &(*(data_unit[qt_id][a][b])));
                        // printf("\nlocation2 = %8X\n", (unsigned int) &((*data_unit[qt_id][a][b])[5][4]));
                        // printf("\nlocation2 = %8X\n", (unsigned int) &((*test)[5][4]));
                        for (int x = 0; x < 8; x++) {
                            for (int y = 0; y < 8; y++) {
                                printf("%g ",data_unit[qt_id][a][b][x][y]);
                            }
                            printf("\n");
                        }
                    }
                }
            }
            //1 = Y, 2 = Cb, 3 = Cr, 4 = I, 5 = Q
            //4) 將所有的 8bit 數加上 128
            //另外，由於離散餘弦變化要求定義域的對稱，所以在編碼時把RGB的數值範圍從[0，255]統一減去128偏移成[-128，127]。
            //因此解碼時必須為每個份量加上128。具體公式如下：
            // R=Y                       +1.402*Cb     +128;
            // G=Y-0.34414*Cr    -0.71414*Cb   +128;
            // B=Y                       +1.772*Cb     +128;
            //5) 轉換 YCbCr 到 RGB
            //雖然我們在本章還不用知道這個問題的答案，但是看到每個分量的大小不同，很自然就會產生這個疑問吧？
            // 所有採樣率都爲 1 的時候，MCU 跟各個顏色分量都是 8 * 8 的正方形，直接一一對應即可。
            // 但是在第一張圖中的例子，MCU 有 16 * 16 ，Cb, Cr 都只有 8 * 8 ，那要如何對應呢？
            // 按照採樣率來對應，例如下圖， Cb, Cr 的最左上角對應到了 MCU 的左上四個 px ：
            // 用虛擬碼表示會更清晰：
            // # 欲計算 MCU[i][j].Cb ，也就是 MCU 的第 i, j 個像素的 Cb 分量數值。

            // # Cb 儲存 Cb 顏色分量的數值
            // # horizontal sampling（水平採樣率） 縮寫爲 hs
            // # vertical sampling（垂直採樣率） 縮寫爲 vs
            // max_hs  = max(Y.hs, Cb.hs, Cr.hs)
            // max_vs = max(Y.vs, Cb.vs, Cr.vs)

            // # 計算 i, j 縮小到 Cb 的大小時，應該索引爲多少
            // # 注意此處的除號 "/" 都會捨去到整數
            // new_i = i * Cb.hs / max_hs
            // new_j = j * Cb.vs / max_vs

            // # Cb 是一個四階陣列
            // # 前兩階描述 block 的位置，後兩階描述要擷取的是這 8*8 中的哪一個點
            // MCU[i][j].Cb = Cb[new_i / 8][new_j / 8][new_i % 8][new_j % 8]

            // # Y, Cr 的算法跟 Cb 完全相同，省略之
            printf("????\n");
            rbg_image[i][j].r = chomp(data_unit[0][i*block_index_i_Y/8][j*block_index_j_Y/8][i*block_index_i_Y%8][j*block_index_j_Y%8]
                                      + 1.402*data_unit[2][i*block_index_i_Cr/8][j*block_index_j_Cr/8][i*block_index_i_Cr%8][j*block_index_j_Cr%8]
                                      + 128);
            rbg_image[i][j].g = chomp(data_unit[0][i*block_index_i_Y/8][j*block_index_j_Y/8][i*block_index_i_Y%8][j*block_index_j_Y%8]
                                      - 0.34414*data_unit[1][i*block_index_i_Cb/8][j*block_index_j_Cb/8][i*block_index_i_Cb%8][j*block_index_j_Cb%8]
                                      - 0.71414*data_unit[2][i*block_index_i_Cr/8][j*block_index_j_Cr/8][i*block_index_i_Cr%8][j*block_index_j_Cr%8]
                                      + 128);
            rbg_image[i][j].b = chomp(data_unit[0][i*block_index_i_Y/8][j*block_index_j_Y/8][i*block_index_i_Y%8][j*block_index_j_Y%8]
                                      + 1.772*data_unit[1][i*block_index_i_Cb/8][j*block_index_j_Cb/8][i*block_index_i_Cb%8][j*block_index_j_Cb%8]
                                      + 128);

        }
    }
    // double** mcu[mcu_number_row][mcu_number_col][3];
    // 每個 component: 3 bytes
    //  - component id (1 = Y, 2 = Cb, 3 = Cr, 4 = I, 5 = Q)
    //  - 採樣係數 (bit 0-3 vert., 4-7 hor.)
    //  - quantization table 號
    // return mcu;
}

int main(int argc,char* argv[])
{
    if (argc != 2) {
        printf("[ERROR]:\nusage: ");
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
    double*** mcu;
    init_cos_cache();
    while (fread(&h,1, 1,fp)) {
        // printf("fuck");
        fread(&l,1, 1,fp);
        if (h == 0xff) {
            if (!b_SOI && l == SOI) {
                b_SOI = true;
            }
            assert(b_SOI == true);
            // if (c.lowwer_byte >= APP0 && c.lowwer_byte <= APP15) {
            //     continue;
            // }
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
                // b_EOI = true;
                break;
            }
        }
    }
    // assert(b_EOI==true);
    return 0;
}
