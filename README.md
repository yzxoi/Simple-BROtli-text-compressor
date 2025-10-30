# Simple BROtli 文件压缩器

Simple-BROtli-text-compressor is a lightweight C++ implementation of a Brotli-inspired lossless compressor. It uses a 32KB LZ77 front-end, 4-context canonical Huffman coding, and bucket-based length/distance encoding to achieve good compression on structured text and log files, while keeping it easy to study.

这是我的同济大学 2025 秋《高级语言程序设计（进阶）》课程“大作业题组（二）文件压缩”作业。本项目实现了一个类 Brotli 思路的、但极为简化的通用无损压缩。

本项目的目标是：
 - 不调用现成压缩库，完整实现无损压缩/解压缩；
 - 结构化设计并在压缩的过程中自校验；
 - 压缩率优秀；
 - 便于后续中加入“分块压缩 / 静态字典”等扩展。

在测试实例 `ser.log` 中，压缩率达 $4.12\%$，压缩耗时 136 ms，解压缩耗时 33 ms。

## Quick Start

### 1. Compile 
We recommend you to use -O2 command to get better performance, although we still can pass the test without it.
```shell
g++ sbro.cpp -o sbro -O2
```

For users using Visual Studio with msvc compiler, please select Release x64 mode and use -O2 command for better performance.

### 2. Run
```shell
# Compress
./sbro ser.log ser.log.sbro zip
# Decompress
./sbro ser.log.sbro recover.log unzip

For Visual Studio Users:
# Compress
.\Simple-BROtli-text-compressor.exe ..\..\..\ser.log ..\..\..\ser.log.sbro zip
# Decompress
.\Simple-BROtli-text-compressor.exe ..\..\..\ser.log.sbro ..\..\..\recover.log unzip
```

The program will output the time and compression ratio automatically.

## 算法介绍

本项目融合了 LZ77 匹配、4 路上下文 Huffman 编码、分桶编码。

### 2.1 容器格式

压缩后的 `.sbro` 文件整体结构如下：
1.	`'S' 'B' 'R' 'O'`
2.	`<1 byte: version>` 版本号：1 字节
3.	`<uint32_t: raw size>` 原始数据长度（4B, LE）
4.	`<uint16_t: ins/cop/dst alphabet size>` 3 个 Huffman 字典的大小：插入长度表大小、拷贝长度表大小、距离表大小（各 2B, LE）
5.	`<256 bytes * 4: literal code lengths for 4 contexts>` 4 组 Huffman 码长表（每组 256B，一共 1024B）
6.	`<insA bytes: insert length code lengths>` 插入长度 Huffman 码长表
7.	`<copA bytes: copy length code lengths>` 拷贝长度 Huffman 码长表
8.	`<dstA bytes: distance code lengths>` 距离 Huffman 码长表
9.	`<bitstream>` bit 流（命令序列）

解码器 decompress_sbro 就是按这个顺序把码长表读出来，重建 Huffman，解码命令流。

### 2.2 LZ77（32KiB 窗口）

我们设定窗口大小：32 KiB（`WND = 32768`），与较多生产环境下的压缩算法一致；最小匹配长度：3 字节（`MIN_MATCH = 3`）；同时，为了控制复杂度，限制每个哈希槽最多追溯 64 个候选位置（`MAX_CANDS = 64`）。

算法流程：
1.	以前 3 字节作为 key，去哈希表里找历史上出现过的相同 3 字节的位置；
2.	对这些候选位置，回溯比对，找出距离不超过 32 KiB 的最长的匹配；
3.	如果匹配长度 ≥ 3，就发出一条 “有匹配的命令”，并打包之前的 literals；
4.	否则就把当前字节塞进 literals；
5.	整个输入结束后，把末尾的 literals 也打包。

LZ77 Command 格式如下：
```cpp
struct Command {
    std::vector<uint8_t> literals; // 先要输出的字面量
    bool hasMatch;                 // 有没有匹配
    uint32_t matchLen;             // 匹配长度
    uint32_t distance;             // 回溯距离
};
```

### 2.3 4 路上下文建模

在很多文本中，前一个字符的类型会影响下一个字符的分布：比如字母之后大概率还是字母、数字之后大概率可能是分隔符、空白后大概率可能是大写字母......

我们利用这个简单的性质简单做了一个 4 分类（见 `charContext(...)`），这其实就是一种非常轻量化的 context modeling：
1.	0：A–Z, a–z
2.	1：0–9
3.	2：空白字符（空格、制表、换行等）
4.	3：其他

于是建立了 4 张独立的 256 符号 Huffman 表：
```cpp
Huffman lit[4];
```
编码时，对当前要写入的字节，先看前一个解码出的字节是什么类型，选择对应的那一张 Huffman 表去编码。

### 2.4 BucketCoder：长度 / 距离（自然数）的分桶编码

匹配长度 matchLen 和距离 distance 等这类值（自然数）有一个共同特点：小值多，大值偶尔出现。直接做 Huffman 编码可能不太经济，于是我们实现了一个非常简单的分桶编码：
- 把一个非负整数 $\text{value}$ 拆成两部分：
- 符号位 Sym：表示它落在哪个 $2^k$ 桶里，其中 $k$ 为 $\text{value}$ 的二进制下最高位
- 附加位 Ex：表示它在这个桶里的偏移，只需记 $\text{value} - 2^k$ 即可。

这样做的优点是：
- 常见的小值都会变成很小的 Sym，对应的 Huffman 码也会特别短
- 不常见的大值只是在后面多写几比特，并不会污染整个 Huffman 结构
- 解压缩算法清晰

在实际实现中，我们把 3 类量都用这套方式处理：
1.	插入字节数 literals.size()
2.	匹配长度偏移 matchLen - 3
3.	距离偏移 distance - 1

### 2.5 Huffman

本项目中的 Huffman 不只是做编码，还会从码长表中重建树（解码）。

值得注意的是，为了方便存储读写，我们将所有码值都用反转码存储。

## 3. 理论压缩率分析

注意：这里的“理论”是指基于算法结构的上、下界与主要影响因素，不是指对所有数据都成立的固定数值。压缩比跟数据的冗余度、字符分布、是否有长重复段强相关。

### 3.1 固定开销

头部至少包括：
1. 标识符 SBRO + 版本号：5 B
2. 原始大小：4 B
3. 3 个 alphabet 大小：6 B
4. 4 × 256 B 的字面量码长表：1024 B

所以无论输入如何，都会有 ~1039 B 左右的固定开销。

### 3.2 LZ77

LZ77 的核心作用是把很多重复的字节串，变成“(len, dist)”的命令。
经典场景：
- 日志里大量重复的时间戳/路径前缀/IP 段；
- 日志、HTML、JSON 等结构化文本（有固定关键词，比如 GET, HTTP 等）；
- 多行之间结构一致、只有少量字段不同的文本。

### 3.3 4 路 Huffman 的贡献

我们需要为此付出 1024 B 的码长存储，所以它的贡献需要中等以上体量的文本里才能显现。

## 4. 理论时空复杂度分析

### 4.1 时间复杂度
压缩算法的复杂度主要在于 LZ77 解析，以及 4 路 Huffman 的构建。

由于我们限制 LZ77 的窗口大小为 32 KiB，且设定每个位置最多回溯 64 个候选，因此实际开销可以写成：$\mathcal O(n \cdot C)$，其中 $C=64\times 256$（实际更小）。

采用优先队列构建 Huffman 树，其复杂度为 $\mathcal O(256\log 256)$

所以压缩端总体时间复杂度可以写成：

$$T_{\text{compress}}(n) = \mathcal O(n\cdot C)$$

其中，$C$ 代表一个小因子（哈希查重 + 候选匹配）。

另一方面，由于每个输出字节恰好只被“写一次”，所以解压就是严格的 $\mathcal O(n)$。

即，

$$T_{\text{decompress}}(n) = \mathcal O(n)$$

### 4.2 空间复杂度
- LZ77 要维护一个哈希表，大小是输入长度的常数倍；
- LZ77 的命令序列要完整存一份，最坏的空间复杂度为 $\mathcal O(n)$；
- 自检过程占用 $\mathcal O(n)$ 空间。

所以总空间复杂度是 $\mathcal O(n)$。

## TODO

- [ ] 把 charContext 换成基于字频/字符集的自适应划分
- [ ] 对特定日志字段做结构化压缩（日期、IP、方法、路径）
- [ ] 添加静态字典
- [ ] 分块压缩

## LICENSE
[MIT](LICENSE)