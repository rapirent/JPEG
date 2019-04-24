#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include "util.h"
#include "./bmitmap/bitmap_image.hpp"




int main(int argc, char* argv[]) {
    if (argc != 2  && argc != 3) {
        printf("[ERROR]:\nusage: ./encoder <input_file> [<output_file_name>]");
        exit(1);
    }
    FILE* fp;
    if ((fp = fopen(argv[1], "r")) == NULL) {
        printf("%s can't be opened\n", argv[1]);
        exit(1
    }
    fwrite(0xff, 1, 1, fp);
    fwrite(SOI, 1, 1, fp);
    fwrite(0xff, 1, 1, fp);
    fwrite(APP0, 1, ,1, fp);
    color_model();
    //DQT
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
