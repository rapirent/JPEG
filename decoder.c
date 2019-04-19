#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>
#include "util.h"

typedef struct {
    byte horizontal_sample;
    byte vertical_sample;
    byte qantize_table_id;
} frame_component;

typedef struct {
    int height;
    int width;
    int vmax;
    int hmax;
    byte components_num;
    byte precision;
    frame_component frame_components[5];
} frame_data;

typedef struct huffman_leaf {
    int codeword;
    int codeword_len;
    byte content;
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

// 查阅标记SOF0，可以得到图像不同颜色分量的采样因子，
// 即Y、Cr、Cb三个分量各自的水平采样因子和垂直采样因子。
// 大多图片的采样因子为4：1：1或1：1：1。
// 其中，4：1：1即（2*2）：（1*1）：（1*1））；1：1：1即（1*1）：（1*1）：（1*1）。
// 记三个分量中水平采样因子最大值为Hmax，垂直采样因子最大值为Vmax，
// 那么单个MCU矩阵的宽就是Hmax*8像素，高就是Vmax*8像素。

int quantize_table_list[4][64];//取值範圍0~3
frame_data f0;
huffman_table_list huffman_table[4];
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
            }
            else {
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
    //   - 数据精度 (1 byte) 每个样本位数, 通常是 8 (大多数软件不支持 12 和 16)
    // fscanf(fp,"%"SCNd8,&f0.precision);
    fread(&(f0.precision),1,1,fp);
    //   - 图片高度 (高字节, 低字节), 如果不支持 DNL 就必须 >0
    f0.height = read_word_to_bigendian(fp);
    //   - 图片宽度 (高字节, 低字节), 如果不支持 DNL 就必须 >0
    f0.width = read_word_to_bigendian(fp);
    //   - components 数量(1 byte), 灰度图是 1, YCbCr/YIQ 彩色图是 3, CMYK 彩色图是 4
    // fscanf(fp,"%"SCNd8,&f0.components_num);
    fread(&(f0.components_num),1,1,fp);
    byte component_id;
    byte sample_cof;
    f0.vmax = -1;
    f0.hmax = -1;
    for (int i = 0; i < f0.components_num; i++) {
        //   - 每个 component: 3 bytes
        //      - component id (1 = Y, 2 = Cb, 3 = Cr, 4 = I, 5 = Q)
        // fscanf(fp,"%"SCNd8,&component_id);
        fread(&component_id,1,1,fp);
        //MCU 的寬 = 8 * 最高水平採樣率
        // MCU 的高 = 8 * 最高垂直採樣率

        //      - 采样系数 (bit 0-3 vert., 4-7 hor.)
        // fscanf(fp,"%"SCNd8,&sample_cof);
        fread(&sample_cof,1,1,fp);
        (f0.frame_components[component_id-1]).horizontal_sample = sample_cof >> 4;
        (f0.frame_components[component_id-1]).vertical_sample = sample_cof & 0x0f;
        if( (f0.frame_components[component_id-1]).horizontal_sample > f0.hmax) {
            f0.hmax = (f0.frame_components[component_id-1]).horizontal_sample;
        }
        if( (f0.frame_components[component_id-1]).vertical_sample > f0.vmax) {
            f0.vmax = (f0.frame_components[component_id-1]).vertical_sample;
        }
        // fscanf(fp,"%"SCNd8,&((f0.frame_components[component_id-1]).qantize_table_id));
        fread(&((f0.frame_components[component_id-1]).qantize_table_id),1,1,fp);
    }
    // return f0;
}

void read_ht(FILE* fp)
{
    //      bit 0..3: HT 號 (0..3, 否則錯誤)
    //      bit 4   : HT 類型, 0 = DC table, 1 = AC table
    //      bit 5..7: 必須是 0
    //   - 16 bytes: 長度是 1..16 代碼的符號數. 這 16 個數的和應該 <=256
    //   - n bytes: 一個包含了按遞增次序代碼長度排列的符號表
    //     (n = 代碼總數)
    byte height_info[16];
    int len = read_word_to_bigendian(fp);
    int ht_node_num;
    byte content;
    byte type_id;
    while (len > 0) {
        // fscanf(fp,"%"SCNd8,&type_id);
        fread(&type_id,1,1,fp);
        // (type_id >> 4) & 0x0f
        // type_id & 0x0f
        // 0x0F bit-pattern 0000 1111
        //      bit 0..3: HT 號 (0..3, 否則錯誤)
        //      bit 4   : HT 類型, 0 = DC table, 1 = AC table
        //      bit 5..7: 必須是 0
        // 這個字節的值為一般只有四個0x00、0x01、0x10、0x11。
        // 0000 0000 | 0000 0001 | 0001 0000 | 0001 0001
        // 0x00表示DC直流0號表；
        // 0x01表示DC直流1號表；
        // 0x10表示AC交流0號表；
        // 0x11表示AC交流1號表。
        len-=1;
        memset(height_info, 0, sizeof(height_info));
        ht_node_num = 0;
        for(int i = 0; i<16; i++) {
            // if(fscanf(fp,"%"SCNd8,&height_info[i])!=1) {
            //     printf("read HT height info error");
            //     exit(1);
            // }
            fread(&height_info[i],1,1,fp);
            ht_node_num += height_info[i];
        }
        len-=16;
        //array implement HT, root is HT[1]
        huffman_table[huffman_table_index(type_id)].ht_node_num = ht_node_num;
        huffman_leaf* ht_leaf = (huffman_leaf*)malloc((ht_node_num)*sizeof(
                                    huffman_leaf));
        int leaf_index = 0;
        unsigned int codeword = 0;
        for(int height=1; height<=16; height++) {
            for(int j = 0; j<height_info[height]; j++) {
                // fscanf(fp,"%"SCNd8,&(ht_leaf[leaf_index].content));
                fread(&(ht_leaf[leaf_index].content),1,1,fp);
                ht_leaf[leaf_index].codeword_len = height;
                // 若高度相等：leaf[n] = leaf[n — 1] + 1
                // 若高度差 1 ：leaf[n] = (leaf[n — 1] + 1) * 2
                // 若高度差 k ：leaf[n] = (leaf[n — 1] + 1) * 2^k
                //int foo = 0b1010;
                //https://medium.com/@yc1043/%E8%B7%9F%E6%88%91%E5%AF%AB-jpeg-%E8%A7%A3%E7%A2%BC%E5%99%A8-%E4%B8%89-%E8%AE%80%E5%8F%96%E9%87%8F%E5%8C%96%E8%A1%A8-%E9%9C%8D%E5%A4%AB%E6%9B%BC%E8%A1%A8-585f2cf4c494
                ht_leaf[leaf_index].codeword = ((int)pow(2,
                                                height - ht_leaf[leaf_index-1].codeword_len)) * (codeword);
                leaf_index++;
                codeword++;
            }
            codeword=codeword<<1;
        }
        len -= ht_node_num;
        huffman_table[huffman_table_index(type_id)].start = ht_leaf;
    }
    // return huffman_table;
}

void read_sos(FILE* fp)
{
    //     SOS: Start Of Scan:
    //   - 長度 (高字節, 低字節), 必須是 6+2*(掃瞄行內組件的數量)
    int len = read_word_to_bigendian(fp);
    //   - 掃瞄行內組件的數量 (1 byte), 必須 >= 1 , <=4 (否則是錯的) 通常是 3
    byte component_num;
    // fscanf(fp,"%"SCNd8,&component_num);
    fread(&component_num,1,1,fp);
    assert(len == 6+2*component_num);
    assert(component_num==3);//JFIF defined
    //   - 每個組件: 2 bytes
    //      - component id (1 = Y, 2 = Cb, 3 = Cr, 4 = I, 5 = Q), 見 SOF0
    byte component_mapping_huffman[5];
    byte component_id,discard;
    // for(int i =0;i<component_num;i++) {
    for (int i =0; i<3; i++) {

        // fscanf(fp,"%"SCNd8,&component_id);
        fread(&component_id,1,1,fp);
        //      - 使用的 Huffman 表:
        //  - bit 0..3: AC table (0..3)
        //  - bit 4..7: DC table (0..3)
        // fscanf(fp,"%"SCNd8,&component_mapping_huffman[component_id]);
        fread(&component_mapping_huffman[component_id],1,1,fp);
    }
    //   - 忽略 3 bytes (???)
    // fscanf(fp,"%"SCNd8,&discard);
    // fscanf(fp,"%"SCNd8,&discard);
    // fscanf(fp,"%"SCNd8,&discard);
    fread(&discard,1,1,fp);
    fread(&discard,1,1,fp);
    fread(&discard,1,1,fp);
    // return component_mapping_huffman;
}

bool fscanf_bit(FILE *fp)
{
    static byte buffer;
    //use count to do 8 times after read a byte
    static byte count = 0;
    byte check_ff00;
    //每八一次
    if (!count) {
        // fscanf(fp,"%"SCNd8,&buffer);
        fread(&buffer,1,1,fp);
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
    bool bit = buffer & (1 << (7 - count));
    count = (count == 7 ? 0 : count + 1);
    return bit;
}


int codeword_decode (FILE *fp, byte code_length)
{
    byte leading;
    leading = fscanf_bit(fp);
    int decoding_code = 1;
    byte c;
    for (int i = 1; i < code_length; i++) {
        c = fscanf_bit(fp);
        decoding_code = decoding_code << 1;
        decoding_code += leading ? c : !c;
        //正數就是照一般二進位計算，負數就是對正數的碼字取反而已
    }
    decoding_code = leading ? decoding_code: -decoding_code;
    return decoding_code;
}
double (*calcualte_mcu_block(FILE *fp, byte component_id))[8] {
    //https://github.com/MROS/jpeg_tutorial/blob/d90271bf96da4f0ea772597aa2be74cb83e09296/doc/%E8%B7%9F%E6%88%91%E5%AF%ABjpeg%E8%A7%A3%E7%A2%BC%E5%99%A8%EF%BC%88%E5%9B%9B%EF%BC%89%E8%AE%80%E5%8F%96%E5%A3%93%E7%B8%AE%E5%9C%96%E5%83%8F%E6%95%B8%E6%93%9A.md
    static int dc_block[5] = {0, 0, 0, 0, 0};
    // 直流變量：
    // 不斷從數據流中讀取一個 bit，直到可以對上直流霍夫曼表中的一個碼字，
    // 取出的對應信源編碼代表的是接下來還要讀入幾位，
    // 假如是 n 就繼續讀取 n bits ，以下表解碼後，就是直流變量。
    //  - bit 0..3: AC table (0..3)
    //  - bit 4..7: DC table (0..3)
    // 前一個代表直流(0)或交流(1)
    // 0000 0000 | 0000 0001 | 0001 0000 | 0001 0001
    double block[8][8];
    memset(block,0,sizeof(block));
    //DC在高位，必須轉為0x00或0x01
    huffman_leaf* dc_table = huffman_table[huffman_table_index(component_mapping_huffman[component_id] >> 4)].start;
    byte code_length;
    unsigned int codeword = 0;
    bool find_it = false;
    for (int i = 0; i< huffman_table[huffman_table_index(component_mapping_huffman[component_id] >> 4)].ht_node_num; i++)
    {
        codeword = codeword << 1;
        codeword+=fscanf_bit(fp);
        if(dc_table[i].codeword == codeword) {
            //找到對應碼字
            code_length = dc_table[i].content;
            find_it = true;
        }
    }
    assert(find_it==true);
    dc_block[component_id] += !code_length ? 0 : codeword_decode(fp, code_length);
    block[0][0] = dc_block[component_id];
    //     讀取交流變量
    // 再接着取出 bit 直到對上交流霍夫曼表的一個碼字，取出對應信源編碼。
    // 這個信源編碼代表的意義爲
    // AC在低位，只需要在高位加入0x1
    huffman_leaf* ac_table = huffman_table[huffman_table_index((component_mapping_huffman[component_id] & 0x0f) | 0x10)].start;
    for(int i=0; i<63;)
    {
        codeword = 0;
        find_it = false;
        for (int j = 0;
             j< huffman_table[huffman_table_index((component_mapping_huffman[component_id] &
                                                                                           0x0f) | 0x10)].ht_node_num; j++) {
            codeword = codeword << 1;
            codeword+=fscanf_bit(fp);
            if(dc_table[j].codeword == codeword) {
                //找到對應碼字
                code_length = dc_table[j].content;
                find_it = true;
            }
        }
        assert(find_it==true);
        // 高 4 bits 代表接下來的數值連續有幾個 0
        // 低 4 bits 代表這些 0 之後跟著的數值的位數，假如是 n 就繼續讀取 n bits ，以下表解碼後，就是這些 0 之後跟着的數值。
        // 有兩個特殊的信源編碼，含義與上述不同
        // 0x00 代表接下來所有的交流變量全爲 0
        // 0xF0 代表接下來有 16 個 0
        // 當以下兩種狀況發生時，代表交流變量已經讀取完畢
        // 1.63 個交流編碼都已經讀完
        // 2.讀到 0x00 ，直到剩下的交流變量全爲 0
        switch(code_length) {
        case 0x00:
            for (int k = 0; k<16; k++) {
                block[i/8][i%8] = 0.0;
                i++;
            }
            break;
        case 0xf0:
            while(i<63) {
                block[i/8][i%8] = 0.0;
                i++;
            }
            break;
        default:
            for (int k=0; k<(code_length >> 4); k++) {
                block[i/8][i%8] = 0.0;
                i++;
            }
            block[i/8][i%8] = codeword_decode(fp, code_length & 0x0f);
            i++;
            break;
        }
    }
    return block;
}

double**** mcu_component(FILE *fp, byte id)
{
    //一個顏色分量內部各個 block 的順序:由左到右，再由上到下
    double** mcu_block[f0.frame_components[id].horizontal_sample][f0.frame_components[id].vertical_sample];
    for (int i =1; i<f0.frame_components[id].horizontal_sample; i ++) {
        for (int j = 1; j<f0.frame_components[id].vertical_sample; j++) {
            mcu_block[i][j] = calcualte_mcu_block(fp, id);
        }
    }
    return mcu_block;
}
//TODO
void dequantize() {}
void dezigzag() {}
void inversedct() {}
void upsample() {}

void calculate_mcu(FILE* fp)
{
    //MCU 的寬 = 8 * 最高水平採樣率
    int mcu_width = 8 * f0.hmax;
    // MCU 的高 = 8 * 最高垂直採樣率
    int mcu_height = 8 * f0.vmax;
    //TODO
    int mcu_number_col = ((f0.height + mcu_height - 1) / mcu_height);
    int mcu_number_row = ((f0.width + mcu_width - 1) / mcu_width);
    //MCU 的順序:由左到右，再由上到下
    double** mcu[mcu_number_col][mcu_number_row][3];
    for (int i  = 0; i<mcu_number_col; i++) {
        for (int j = 0; j<mcu_number_row; j++) {
            for (int id = 0 ; id<3; id++) {
                mcu[i][j][id] = mcu_component(fp,id);
            }
            //decode
            dequantize();
            dezigzag();
            inversedct();
            upsample();
        }
    }
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
                printf("read DQT\n");
                read_qt(fp);
                break;
            case SOF0:
                printf("read SOF0\n");
                // read_frame(fp);
                break;
            case DHT:
                printf("read DHT\n");
                // read_ht(fp);
                break;
            case SOS:
                printf("read SOS\n");
                // read_sos(fp);
                // calculate_mcu(fp);
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
                printf("APP section: discard\n");
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
