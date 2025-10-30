#include <cstdint>
#include <vector>
#include <string>
#include <array>
#include <queue>
#include <unordered_map>
#include <stdexcept>
#include <algorithm>
#include <functional>
#include <fstream>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <sys/stat.h>

namespace edu {

// ========== Bit I/O (LSB-first) ==========
struct BitWriter {
    std::vector<uint8_t> out;
    uint64_t buf = 0;
    int bitcnt = 0;

    void writeBits(uint32_t v, int n) {
        if (n <= 0) return;
        uint64_t mask = (n==32)? 0xFFFFFFFFull : ((1ull<<n)-1ull);
        buf |= (uint64_t(v) & mask) << bitcnt;
        bitcnt += n;
        while (bitcnt >= 8) {
            out.push_back(uint8_t(buf & 0xFF));
            buf >>= 8;
            bitcnt -= 8;
        }
    }
    void writeBit(uint32_t b){ writeBits(b&1u, 1); }

    void flushTo(std::vector<uint8_t>& dst){
        if (bitcnt > 0) {
            out.push_back(uint8_t(buf & 0xFF));
            buf = 0; bitcnt = 0;
        }
        dst.insert(dst.end(), out.begin(), out.end());
        out.clear();
    }
};

struct BitReader {
    const uint8_t* p; size_t n;
    size_t idx = 0;
    uint64_t buf = 0;
    int bitcnt = 0;

    BitReader(const uint8_t* p_, size_t n_) : p(p_), n(n_) {}

    uint32_t readBits(int nbits){
        while (bitcnt < nbits){
            if (idx >= n) throw std::runtime_error("BitReader: out of bytes");
            buf |= (uint64_t)p[idx++] << bitcnt;
            bitcnt += 8;
        }
        uint32_t v = (nbits==32)? (uint32_t)buf : (uint32_t)(buf & ((1ull<<nbits)-1ull));
        buf >>= nbits; bitcnt -= nbits;
        return v;
    }
    uint32_t readBit(){ return readBits(1); }
};

// ========== Huffman (canonical) ==========
struct Huffman {
    std::vector<uint8_t> codeLen;    // per symbol; 0 => not used
    std::vector<uint32_t> code;      // MSB-first canonical code
    std::vector<uint32_t> codeRev;   // LSB-first to write
    int alphabet = 0;

    struct Node { int l=-1, r=-1, sym=-1; };
    std::vector<Node> nodes;

    static uint32_t reverseBits(uint32_t x, int len){
        uint32_t r=0;
        for(int i=0;i<len;i++){ r = (r<<1) | ((x>>i)&1u); }
        return r;
    }

    void buildFromFreq(const std::vector<uint64_t>& freq){
        alphabet = (int)freq.size();
        codeLen.assign(alphabet, 0);
        code.assign(alphabet, 0);
        codeRev.assign(alphabet, 0);

        std::vector<int> syms; syms.reserve(alphabet);
        for(int i=0;i<alphabet;i++) if (freq[i]>0) syms.push_back(i);

        if (syms.empty()){
            if (alphabet==0) { buildDecTree(); return; }
            codeLen[0]=1; code[0]=0; codeRev[0]=0; buildDecTree(); return;
        }
        if ((int)syms.size()==1){
            codeLen[syms[0]] = 1;
            code[syms[0]] = 0;
            codeRev[syms[0]] = 0;
            buildDecTree(); return;
        }

        struct T{ uint64_t f; int id; int l=-1, r=-1; int sym=-1; };
        std::vector<T> tn;
        auto mkLeaf=[&](int s, uint64_t f){ tn.push_back({f,(int)tn.size(),-1,-1,s}); return tn.back().id; };
        auto mkParen=[&](int L,int R, uint64_t f){ tn.push_back({f,(int)tn.size(),L,R,-1}); return tn.back().id; };

        struct Item{ uint64_t f; int id; bool operator>(Item const& o)const{return f>o.f;} };
        std::priority_queue<Item, std::vector<Item>, std::greater<Item>> pq;
        for(int s: syms) pq.push({freq[s], mkLeaf(s, freq[s])});

        int root=-1;
        while(pq.size()>=2){
            Item a=pq.top(); pq.pop();
            Item b=pq.top(); pq.pop();
            int id = mkParen(a.id, b.id, a.f + b.f);
            pq.push({tn[id].f, id});
            root = id;
        }
        if (root==-1) root = pq.top().id;

        std::function<void(int,int)> dfs = [&](int u,int d){
            if (tn[u].sym!=-1){
                codeLen[tn[u].sym] = (d==0?1:d);
                return;
            }
            dfs(tn[u].l, d+1);
            dfs(tn[u].r, d+1);
        };
        dfs(root, 0);

        int maxL=0; for (auto L: codeLen) maxL = std::max(maxL, (int)L);
        std::vector<int> bl_count(maxL+1,0);
        for(int i=0;i<alphabet;i++) if (codeLen[i]) bl_count[codeLen[i]]++;

        std::vector<uint32_t> next_code(maxL+1,0);
        uint32_t c=0;
        for(int bits=1; bits<=maxL; ++bits){
            c = (c + bl_count[bits-1]) << 1;
            next_code[bits] = c;
        }
        for(int i=0;i<alphabet;i++){
            int len = codeLen[i];
            if (!len) continue;
            code[i] = next_code[len]++;
            codeRev[i] = reverseBits(code[i], len);
        }
        buildDecTree();
    }

    void buildDecTree(){
        nodes.clear(); nodes.push_back(Node());
        for(int s=0;s<alphabet;s++){
            int len = codeLen[s];
            if (!len) continue;
            uint32_t rev = codeRev[s];
            int cur=0;
            for(int i=0;i<len;i++){
                int b = (rev>>i)&1;
                int nxt = (b==0? nodes[cur].l : nodes[cur].r);
                if (nxt==-1){
                    nxt = (int)nodes.size();
                    nodes.push_back(Node());
                    if (b==0) nodes[cur].l = nxt; else nodes[cur].r = nxt;
                }
                cur = nxt;
            }
            nodes[cur].sym = s;
        }
    }

    void encSymbol(BitWriter& bw, int s) const {
        int len = codeLen[s];
        bw.writeBits(codeRev[s], len);
    }
    int decSymbol(BitReader& br) const {
        int cur=0;
        while(true){
            if (nodes[cur].sym!=-1) return nodes[cur].sym;
            int bit = br.readBit();
            cur = (bit==0? nodes[cur].l : nodes[cur].r);
            if (cur==-1) throw std::runtime_error("Huffman decode: invalid path");
        }
    }
};

// ========== Bucket coding ==========
struct BucketCoder {
    static inline int ilog2_u32(uint32_t v){
#if defined(_MSC_VER)
        #include <intrin.h>
        unsigned long idx;
        _BitScanReverse(&idx, v);
        return (int)idx;
#else
        return 31 - __builtin_clz(v);
#endif
    }
    struct Enc { uint32_t sym; int exBits; uint32_t exVal; };
    static inline Enc encode(uint32_t value){
        if (value==0) return {0,0,0};
        int k = ilog2_u32(value);
        uint32_t base = 1u<<k;
        return { (uint32_t)(k+1), k, value - base };
    }
    static inline uint32_t decode(uint32_t sym, uint32_t exVal){
        if (sym==0) return 0;
        uint32_t base = 1u<<(sym-1);
        return base + exVal;
    }
	static inline uint32_t decodeFromStream(const Huffman& H, BitReader& br) {
		int Sym = H.decSymbol(br);
        if (Sym==0) return 0;
        else {
            int eb = Sym - 1;
            uint32_t ev = eb? br.readBits(eb) : 0;
            return decode((uint32_t)Sym, ev);
        }
	}
};

// ========== LZ77 command ==========
struct Command {
    std::vector<uint8_t> literals;
    bool hasMatch = false;
    uint32_t matchLen = 0;
    uint32_t distance = 0;
};

inline uint8_t charContext(uint8_t prev){
    if ((prev>='A'&&prev<='Z')||(prev>='a'&&prev<='z')) return 0;
    if (prev>='0'&&prev<='9') return 1;
    if (prev==' '||prev=='\t'||prev=='\n'||prev=='\r'||prev=='\f'||prev=='\v') return 2;
    return 3;
}

// ========== LZ77 (32KiB window) ==========
struct LZ77 {
    static constexpr int WND = 32768;
    static constexpr int MIN_MATCH = 3;
    static constexpr int MAX_CANDS = 64;

    static std::vector<Command> parse(const std::vector<uint8_t>& in){
        std::vector<Command> cmds;
        std::vector<int> litbuf; litbuf.reserve(256);

        std::unordered_map<uint32_t, std::vector<int>> ht;
        ht.reserve(in.size()/4+1);

        auto push_key = [&](int pos){
            if (pos+2 >= (int)in.size()) return;
            uint32_t key = (uint32_t(in[pos])<<16) | (uint32_t(in[pos+1])<<8) | uint32_t(in[pos+2]);
            std::vector<int>& v = ht[key];
            v.push_back(pos);
            if ((int)v.size()>MAX_CANDS*4) v.erase(v.begin(), v.begin()+ (int)v.size()-MAX_CANDS*2);
        };

        int n = (int)in.size();
        push_key(0); push_key(1);

        int i=0;
        while(i<n){
            int bestLen=0, bestDist=0;
            if (i+MIN_MATCH-1 < n){
                uint32_t key = (i+2<n)? ((uint32_t(in[i])<<16)|(uint32_t(in[i+1])<<8)|in[i+2]) : 0xFFFFFFu;
                auto it = ht.find(key);
                if (it!=ht.end()){
                    const std::vector<int>& cand = it->second;
                    int checked=0;
                    for (int idx=(int)cand.size()-1; idx>=0 && checked<MAX_CANDS; --idx, ++checked){
                        int p = cand[idx];
                        int dist = i - p;
                        if (dist<=0 || dist>WND) continue;
                        int L=0;
                        int maxL = std::min(n - i, WND);
                        while (L<maxL && in[p+L]==in[i+L]) ++L;
                        if (L>=MIN_MATCH && L>bestLen){
                            bestLen=L; bestDist=dist;
                            if (L>=258) break;
                        }
                    }
                }
            }

            if (bestLen>=MIN_MATCH){
                Command cmd;
                if (!litbuf.empty()){
                    cmd.literals.reserve(litbuf.size());
                    for(int v: litbuf) cmd.literals.push_back((uint8_t)v);
                    litbuf.clear();
                }
                cmd.hasMatch = true;
                cmd.matchLen = (uint32_t)bestLen;
                cmd.distance = (uint32_t)bestDist;
                cmds.push_back(std::move(cmd)); // 以 std::move 方式放入，避免拷贝大数组，原内容被直接转移

                for(int j=0;j<bestLen;j++) push_key(i+j);
                i += bestLen;
            }else{
                litbuf.push_back(in[i]);
                push_key(i);
                ++i;
            }
        }
        if (!litbuf.empty()){
            Command cmd;
            cmd.literals.reserve(litbuf.size());
            for(int v: litbuf) cmd.literals.push_back((uint8_t)v);
            litbuf.clear();
            cmd.hasMatch=false;
            cmds.push_back(std::move(cmd));
        }
        return cmds;
    }
};
// ---- class-out definitions to avoid ODR linker issues on some toolchains
constexpr int LZ77::WND;
constexpr int LZ77::MIN_MATCH;

// ========== Codebooks ==========
struct Codebooks {
    Huffman lit[4];
    Huffman insLen, copLen, dist;

    std::array<std::vector<uint8_t>,4> litCodeLen;
    std::vector<uint8_t> insCodeLen, copCodeLen, distCodeLen;

    void build(const std::vector<Command>& cmds, const std::vector<uint8_t>& original){
        std::array<std::vector<uint64_t>,4> litFreq;
        for(int c=0;c<4;c++) litFreq[c].assign(256,0);

        std::vector<uint64_t> insFreq(1,0), copFreq(1,0), distFreq(1,0);
        std::vector<uint8_t> out; out.reserve(original.size());

        auto bumpLenBucket = [&](std::vector<uint64_t>& F, uint32_t val){
            auto e = BucketCoder::encode(val);
            int need = (int)e.sym + 1;
            if ((int)F.size() < need) F.resize(need,0);
            F[e.sym]++;
        };

        for (const auto& cmd : cmds){
            for (uint8_t b : cmd.literals){
                uint8_t ctx = out.empty()? 3 : charContext(out.back());
                litFreq[ctx][b]++;
                out.push_back(b);
            }
            bumpLenBucket(insFreq, (uint32_t)cmd.literals.size());

            if (cmd.hasMatch){
                bumpLenBucket(copFreq, cmd.matchLen - 3);
                bumpLenBucket(distFreq, cmd.distance - 1);
                if (cmd.distance==0 || cmd.distance > out.size())
                    throw std::runtime_error("Invalid distance while simulating");
                size_t start = out.size() - cmd.distance;
                for (uint32_t k=0;k<cmd.matchLen;k++) out.push_back(out[start + k]);
            }
        }
        if (out != original) throw std::runtime_error("Command stream does not reconstruct input.");

        for(int c=0;c<4;c++){
            lit[c].buildFromFreq(litFreq[c]);
            litCodeLen[c].assign(256,0);
            for(int s=0;s<256;s++) litCodeLen[c][s]=lit[c].codeLen[s];
        }
        insLen.buildFromFreq(insFreq);
        copLen.buildFromFreq(copFreq);
        dist.buildFromFreq(distFreq);

        insCodeLen.assign(insLen.codeLen.begin(), insLen.codeLen.end()); if (insCodeLen.empty()) insCodeLen.resize(1,1);
        copCodeLen.assign(copLen.codeLen.begin(), copLen.codeLen.end()); if (copCodeLen.empty()) copCodeLen.resize(1,1);
        distCodeLen.assign(dist.codeLen.begin(), dist.codeLen.end());   if (distCodeLen.empty()) distCodeLen.resize(1,1);
    }
};

// ========== Little-endian helpers ==========
static void write_u32_le(std::vector<uint8_t>& buf, uint32_t v){
    buf.push_back(uint8_t(v & 0xFF));
    buf.push_back(uint8_t((v>>8) & 0xFF));
    buf.push_back(uint8_t((v>>16)& 0xFF));
    buf.push_back(uint8_t((v>>24)& 0xFF));
}
static void write_u16_le(std::vector<uint8_t>& buf, uint16_t v){
    buf.push_back(uint8_t(v & 0xFF));
    buf.push_back(uint8_t((v>>8) & 0xFF));
}
static uint32_t read_u32_le(const uint8_t* p){
    return uint32_t(p[0]) | (uint32_t(p[1])<<8) | (uint32_t(p[2])<<16) | (uint32_t(p[3])<<24);
}
static uint16_t read_u16_le(const uint8_t* p){
    return uint16_t(p[0]) | (uint16_t(p[1])<<8);
}

// ========== Encoder ==========
static std::vector<uint8_t> compress_sbro(const std::vector<uint8_t>& input){
    auto cmds = LZ77::parse(input);

    Codebooks cb; cb.build(cmds, input);

    std::vector<uint8_t> out;
    out.insert(out.end(), {'S','B','R','O'});
    out.push_back(1);
    write_u32_le(out, (uint32_t)input.size());

    uint16_t insA = (uint16_t)std::max<size_t>(cb.insCodeLen.size(), 1);
    uint16_t copA = (uint16_t)std::max<size_t>(cb.copCodeLen.size(), 1);
    uint16_t dstA = (uint16_t)std::max<size_t>(cb.distCodeLen.size(),1);
    write_u16_le(out, insA);
    write_u16_le(out, copA);
    write_u16_le(out, dstA);

    for(int c=0;c<4;c++) for(int s=0;s<256;s++) out.push_back(cb.litCodeLen[c][s]);
    for(size_t i=0;i<insA;i++) out.push_back( i<cb.insCodeLen.size()? cb.insCodeLen[i] : 0 );
    for(size_t i=0;i<copA;i++) out.push_back( i<cb.copCodeLen.size()? cb.copCodeLen[i] : 0 );
    for(size_t i=0;i<dstA;i++) out.push_back( i<cb.distCodeLen.size()? cb.distCodeLen[i] : 0 );

    BitWriter bw;
    std::vector<uint8_t> recon; recon.reserve(input.size());

    for (const auto& cmd : cmds){
        auto encIns = BucketCoder::encode((uint32_t)cmd.literals.size());
        cb.insLen.encSymbol(bw, (int)encIns.sym);
        if (encIns.exBits>0) bw.writeBits(encIns.exVal, encIns.exBits);

        for (uint8_t b : cmd.literals){
            uint8_t ctx = recon.empty()? 3 : charContext(recon.back());
            cb.lit[ctx].encSymbol(bw, (int)b);
            recon.push_back(b);
        }

        bw.writeBit(cmd.hasMatch?1u:0u);

        if (cmd.hasMatch){
            auto encLen = BucketCoder::encode(cmd.matchLen - 3);
            cb.copLen.encSymbol(bw, (int)encLen.sym);
            if (encLen.exBits>0) bw.writeBits(encLen.exVal, encLen.exBits);

            auto encDst = BucketCoder::encode(cmd.distance - 1);
            cb.dist.encSymbol(bw, (int)encDst.sym);
            if (encDst.exBits>0) bw.writeBits(encDst.exVal, encDst.exBits);

            if (cmd.distance==0 || cmd.distance > recon.size())
                throw std::runtime_error("Encoder: bad distance");
            size_t start = recon.size() - cmd.distance;
            for (uint32_t k=0;k<cmd.matchLen;k++) recon.push_back(recon[start + k]);
        }
    }
    if (recon != input) throw std::runtime_error("Encoder self-check failed.");

    bw.flushTo(out);
    return out;
}
/*
SBRO
1
<uint32_t: raw size>
<uint16_t: ins alphabet size>
<uint16_t: cop alphabet size>
<uint16_t: dst alphabet size>
<256 bytes * 4: literal code lengths for 4 contexts>
<insA bytes: insert length code lengths>
<copA bytes: copy length code lengths>
<dstA bytes: distance code lengths>
<bitstream>
*/

// ========== Decoder ==========
static std::vector<uint8_t> decompress_sbro(const std::vector<uint8_t>& in){
    if (in.size()<4+1+4+6+4*256) throw std::runtime_error("Input too small");
    if (!(in[0]=='S'&&in[1]=='B'&&in[2]=='R'&&in[3]=='O')) throw std::runtime_error("Bad magic");
    uint8_t ver = in[4]; if (ver!=1) throw std::runtime_error("Unsupported version");
    size_t off = 5;
    uint32_t rawSize = read_u32_le(&in[off]); off+=4;

    uint16_t insA = read_u16_le(&in[off]); off+=2;
    uint16_t copA = read_u16_le(&in[off]); off+=2;
    uint16_t dstA = read_u16_le(&in[off]); off+=2;

    std::array<std::vector<uint8_t>,4> litCL;
    for(int c=0;c<4;c++){ litCL[c].resize(256); for(int s=0;s<256;s++) litCL[c][s]=in[off++]; }
    std::vector<uint8_t> insCL(insA), copCL(copA), dstCL(dstA);
    for(int i=0;i<insA;i++) insCL[i]=in[off++];
    for(int i=0;i<copA;i++) copCL[i]=in[off++];
    for(int i=0;i<dstA;i++) dstCL[i]=in[off++];

    auto buildFromCL = [](Huffman& h, const std::vector<uint8_t>& cl){
        h.alphabet = (int)cl.size();
        h.codeLen.assign(cl.begin(), cl.end());
        h.code.assign(h.alphabet,0);
        h.codeRev.assign(h.alphabet,0);
        if (h.alphabet==0){ h.buildDecTree(); return; }

        int maxL=0; for(auto L: h.codeLen) maxL=std::max(maxL,(int)L);
        if (maxL==0){ h.codeLen.assign(std::max(1,h.alphabet),0); h.codeLen[0]=1; maxL=1; }

        std::vector<int> bl_count(maxL+1,0);
        for(auto L: h.codeLen) if (L) bl_count[L]++;

        std::vector<uint32_t> next_code(maxL+1,0);
        uint32_t codev=0;
        for(int bits=1; bits<=maxL; ++bits){
            codev = (codev + bl_count[bits-1]) << 1;
            next_code[bits] = codev;
        }
        for(int s=0;s<h.alphabet;s++){
            int len=h.codeLen[s];
            if (!len) continue;
            h.code[s]=next_code[len]++;
            h.codeRev[s]=Huffman::reverseBits(h.code[s], len);
        }
        h.buildDecTree();
    };

    Huffman lit[4], insH, copH, dstH;
    for(int c=0;c<4;c++) buildFromCL(lit[c], litCL[c]);
    buildFromCL(insH, insCL);
    buildFromCL(copH, copCL);
    buildFromCL(dstH, dstCL);

    BitReader br(in.data()+off, in.size()-off);

    std::vector<uint8_t> out; out.reserve(rawSize);
    while (out.size() < rawSize){
		uint32_t insVal = BucketCoder::decodeFromStream(insH, br);
        for(uint32_t i=0;i<insVal;i++){
            uint8_t ctx = out.empty()? 3 : charContext(out.back());
            int litSym = lit[ctx].decSymbol(br);
            out.push_back((uint8_t)litSym);
            if (out.size() > rawSize) throw std::runtime_error("Decoded beyond raw size (literals).");
        }
        if (out.size() >= rawSize) break;

        uint32_t hasM = br.readBit();
        if (!hasM) continue;

		uint32_t lenVal = BucketCoder::decodeFromStream(copH, br);
        uint32_t matchLen = lenVal + 3;
		uint32_t dstVal = BucketCoder::decodeFromStream(dstH, br);
        uint32_t dist = dstVal + 1;
        if (dist==0 || dist>out.size()) throw std::runtime_error("Bad distance while decoding");
        size_t start = out.size() - dist;
        for(uint32_t k=0;k<matchLen;k++){
            out.push_back(out[start + k]);
            if (out.size() > rawSize) throw std::runtime_error("Decoded beyond raw size (match).");
        }
    }
    if (out.size()!=rawSize) throw std::runtime_error("Decoded size mismatch");
    return out;
}

// ========== File I/O ==========
// static std::vector<uint8_t> readAll(const std::string& path){
//     FILE* f = std::fopen(path.c_str(), "rb");
//     if (!f) throw std::runtime_error("Cannot open input: "+path);
//     std::fseek(f,0,SEEK_END); long sz = std::ftell(f); std::fseek(f,0,SEEK_SET);
//     std::vector<uint8_t> buf; buf.resize(std::max<long>(sz,0));
//     if (sz>0) std::fread(buf.data(),1,(size_t)sz,f);
//     std::fclose(f);
//     return buf;
// }
// static void writeAll(const std::string& path, const std::vector<uint8_t>& data){
//     FILE* f = std::fopen(path.c_str(), "wb");
//     if (!f) throw std::runtime_error("Cannot open output: "+path);
//     if (!data.empty()) std::fwrite(data.data(),1,data.size(),f);
//     std::fclose(f);
// }
// Using fstream to read/write file instead of C-style FILE due to the fucking problem requirements
static std::vector<uint8_t> readAll(const std::string& path){
    std::ifstream fin(path, std::ios::binary);
    if (!fin.is_open()) throw std::runtime_error("Cannot open input: " + path);
    fin.seekg(0, std::ios::end);
    std::streampos sz = fin.tellg();
    if (sz < 0) throw std::runtime_error("Failed to get size of: " + path);
    std::vector<uint8_t> buf; buf.resize(static_cast<size_t>(sz));
	fin.seekg(0, std::ios::beg);
    if (sz > 0) {
        fin.read((char *)buf.data(), sz);
        if (!fin) throw std::runtime_error("Failed to read file: " + path);
    }
    return buf;
}
static void writeAll(const std::string& path, const std::vector<uint8_t>& data){
    std::ofstream fout(path, std::ios::binary);
    if (!fout.is_open()) throw std::runtime_error("Cannot open output: " + path);
    if (!data.empty()) {
        fout.write((char *)data.data(), (std::streampos)data.size());
        if (!fout) throw std::runtime_error("Failed to write file: " + path);
    }
	// No need to close, destructor will handle it
}

} // namespace edu

// ========== CLI ==========
int main(int argc, char** argv){
    if (argc!=4){
        std::cerr << "Usage:\n"
                  << "  " << argv[0] << " <input> <output> zip\n"
                  << "  " << argv[0] << " <input> <output> unzip\n"
				  << "Example:\n"
				  << "  " << argv[0] << " ser.log ser.log.sbro zip\n"
				  << "  " << argv[0] << " ser.log.sbro ser_rec.log unzip\n";
        return 1;
    }
    std::string inPath=argv[1], outPath=argv[2], mode = argv[3];
    try{
        auto data = edu::readAll(inPath);
        if (mode=="zip"){
			auto start_time = std::chrono::high_resolution_clock::now();
            auto enc = edu::compress_sbro(data);
            edu::writeAll(outPath, enc);
			auto end_time = std::chrono::high_resolution_clock::now();
			auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

			// 获取文件大小用于计算压缩率
			struct stat input_stat, output_stat;
			if (stat(inPath.c_str(), &input_stat) == 0 && stat(outPath.c_str(), &output_stat) == 0)
			{
				double compression_ratio = (double)output_stat.st_size / input_stat.st_size * 100;
				std::cout << "Compression completed in " << duration.count() << " ms\n";
				std::cout << "Original size: " << input_stat.st_size << " bytes\n";
				std::cout << "Compressed size: " << output_stat.st_size << " bytes\n";
				std::cout << "Compression ratio: " << std::fixed << std::setprecision(2) << compression_ratio << "%\n";
			}
        }else if (mode=="unzip"){
			auto start_time = std::chrono::high_resolution_clock::now();
            auto dec = edu::decompress_sbro(data);
            edu::writeAll(outPath, dec);
			auto end_time = std::chrono::high_resolution_clock::now();
			auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

			// 获取文件大小用于计算压缩率
			struct stat input_stat, output_stat;
			if (stat(inPath.c_str(), &input_stat) == 0 && stat(outPath.c_str(), &output_stat) == 0)
			{
				double compression_ratio = (double)input_stat.st_size / output_stat.st_size * 100;
				std::cout << "Decompression completed in " << duration.count() << " ms\n";
				std::cout << "Compressed size: " << input_stat.st_size << " bytes\n";
				std::cout << "Decompressed size: " << output_stat.st_size << " bytes\n";
				std::cout << "Decompression ratio: " << std::fixed << std::setprecision(2) << compression_ratio << "%\n";
			}
        }else{
            throw std::runtime_error("Unknown mode (use zip or unzip)");
        }
    }catch(const std::exception& e){
        std::cerr << "[ERROR] " << e.what() << "\n";
        return 2;
    }
    return 0;
}