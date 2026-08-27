/*
 * Protein-coding Gene Annoation System of CLOVERS
 * 
 * @copyright (C)2026 TUBIC, Tianjin University
 * @authors:  Zetong Z, You Z, Lin Y*, Gao F*
 * @version   1.1.0
 * @date      2026-05-01
 * @modified  2026-08-31
 * @license:  GNU-GPLv3
 */
#pragma once
/* C/C++ Standard Libraries */
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <chrono>
#include <cmath>
#include <assert.h>
#include <array>
#include <sys/stat.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <iomanip>
#include <unordered_set>
#include <iterator>
#include <numeric>
/* Third-party Libraries */
#include "cxxopts.hpp"
#include "svm.hpp"
// Optional
#if defined(__AVX__)
    #include <immintrin.h>
#endif
#if defined(_OPENMP)
    #include <omp.h>
#endif
#if defined(ZLIB)
    #include <zlib.h>
#endif
/* Metagenomic model parameters */
extern const float META[];

/* embedded MLP weights are IEEE-754 binary32 (see src/Meta.cpp) */
static_assert(sizeof(float) == 4, "embedded MLP weights require 32-bit float");
/* system clock */
#define SYSTEM_CLOCK std::chrono::system_clock::
/* to_string */
#define Str std::to_string
/* software version */
#define VERSION "CLOVERS_v1.1.0"
/* no throw memory error new */
#define NEW new (std::nothrow)
/* The buffer size for file reading. */
#define BUFF_SIZE 65536

/* Number of heuristic models */
#define N_MODELS  61
/* Mininum length of ORFs */
#define MIN_LEN   75

/* error codes */
#define NON_ERR exit (0x0000)
#define MEM_ERR exit (0x0001)
#define CFG_ERR exit (0x0002)
#define SUB_ERR exit (0x0003)
#define ARG_ERR exit (0x0004)
#define IOP_ERR exit (0x0005)

/* memory error info */
#define MEM_ERR_INFO "\nError: memory error\n"

/* dimension of all-in features. */
#define DIM_A 190
/* dimension of Z-curve params.  */
#define DIM_S 189

/* abbreviations */
typedef std::string                 str;
typedef std::vector<char *>     str_arr;
typedef std::vector<size_t>     idx_arr;
typedef std::vector<int>        int_arr;
typedef std::unordered_set<int> int_set;
typedef std::vector<float>      flt_arr;
typedef SYSTEM_CLOCK time_point _time_t;

/* debug log stream (inline: one shared instance across translation units) */
inline std::ostream* debug = &std::cerr;
inline std::ofstream DEVNULL;
inline std::ostream& Debug() { return *debug; }

/* translation table */
extern char trans_tbl[64];
/* base-to-address map */
// T/U = 0  C = 1  A = 2  G = 3
const char BASE2ADDR[] = {
    -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1,
    -1,
  /* A   B   C   D   E   F   G   H */
    +2, -1, +1, -1, -1, -1, +3, -1, 
  /* I   J   K   L   M   N   O   P */
    -1, -1, -1, -1, -1, -1, -1, -1, 
  /* Q   R   S   T   U   V   W   X */ 
    -1, -1, -1, +0, +0, -1, -1, -1,
  /* Y   Z */
    -1, -1, -1, -1, -1, -1, -1, -1, 
  /* a   b   c   d   e   f   g   h */
    +2, -1, +1, -1, -1, -1, +3, -1, 
  /* i   j   k   l   m   n   o   p  */
    -1, -1, -1, -1, -1, -1, -1, -1, 
  /* q   r   s   t   u   v   w   x  */
    -1, -1, -1, +0, +0, -1, -1, -1, 
  /* y   z */
    -1, -1
};
// offset constants
const int   X = 0, Y = 1, Z = 2, PHASE = 3;
// float constants
const float V = 1.0F/2, W=1.0F/3, Q=1.0F/4;
/* complement base table */
constexpr char COMP[] = {
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 
    'N',
  /* A    B    C    D    E    F    G    H  */
    'T', 'V', 'G', 'H', 'N', 'N', 'C', 'D', 
  /* I    J    K    L    M    N    O    P  */
    'N', 'N', 'M', 'N', 'K', 'N', 'N', 'N',
  /* Q    R    S    T    U    V    W    X  */
    'N', 'Y', 'W', 'A', 'A', 'B', 'S', 'N', 
  /* Y    Z */
    'R', 'T', 'N', 'N', 'N', 'N', 'N', 'N',
  /* a    b    c    d    e    f    g    h  */
    't', 'v', 'g', 'h', 'n', 'n', 'c', 'd', 
  /* i    j    k    l    m    n    o    p  */
    'n', 'n', 'm', 'n', 'k', 'n', 'n', 'n',
  /* q    r    s    t    u    v    w    x  */
    'n', 'y', 'w', 'a', 'a', 'b', 's', 'n', 
  /* y    z */
    'r', 't'
};
/* Map of one-hot encoding for DNA sequences */
const float ONE_HOT[][4] = 
{   
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0},
    {0, 0, 0, 0},
    {1, 0, 0, 0}, {0, W, W, W}, {0, 0, 1, 0}, {W, W, 0, W},
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 1, 0, 0}, {W, 0, W, W}, 
    {Q, Q, Q, Q}, {0, 0, 0, 0}, {0, V, 0, V}, {0, 0, 0, 0},
    {V, 0, V, 0}, {Q, Q, Q, Q}, {0, 0, 0, 0}, {0, 0, 0, 0},
    {0, 0, 0, 0}, {V, V, 0, 0}, {0, V, V, 0}, {0, 0, 0, 1},
    {0, 0, 0, 1}, {W, W, W, 0}, {V, 0, 0, V}, {0, 0, 0, 0}, 
    {0, 0, V, V}, {1, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0},
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {1, 0, 0, 0}, {0, W, W, W}, {0, 0, 1, 0}, {W, W, 0, W},
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 1, 0, 0}, {W, 0, W, W}, 
    {Q, Q, Q, Q}, {0, 0, 0, 0}, {0, V, 0, V}, {0, 0, 0, 0}, 
    {V, 0, V, 0}, {Q, Q, Q, Q}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {V, V, 0, 0}, {0, V, V, 0}, {0, 0, 0, 1},
    {0, 0, 0, 1}, {W, W, W, 0}, {V, 0, 0, V}, {0, 0, 0, 0},
    {0, 0, V, V}, {1, 0, 0, 0}
};
/* Map of Z-curve encoding for DNA sequences */
const float Z_COORD[][3] = {
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0},
    {+1, +1, +1}, {-W, -W, -W}, {-1, +1, -1}, {+W, -W, +W},
    {+0, +0, +0}, {+0, +0, +0}, {+1, -1, -1}, {-W, +W, +W},
    {+0, +0, +0}, {+0, +0, +0}, {+0, -1, +0}, {+0, +0, +0},
    {+0, +1, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+1, +0, +0}, {+0, +0, -1}, {-1, -1, +1},
    {-1, -1, +1}, {+W, +W, -W}, {+0, +0, +1}, {+0, +0, +0},
    {+1, +0, +0}, {+1, +1, +1}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+1, +1, +1}, {-W, -W, -W}, {-1, +1, -1}, {+W, -W, +W},
    {+0, +0, +0}, {+0, +0, +0}, {+1, -1, -1}, {-W, +W, +W},
    {+0, +0, +0}, {+0, +0, +0}, {+0, -1, +0}, {+0, +0, +0},
    {+0, +1, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+1, +0, +0}, {+0, +0, -1}, {-1, -1, +1},
    {-1, -1, +1}, {+W, +W, -W}, {+0, +0, +1}, {+0, +0, +0},
    {+1, +0, +0}, {+1, +1, +1}
};
/*
 * format debug info to [key:[pad]value]-like logs
 * @param retr  should print a <return> inplace or not
 */
inline str format_log(str key, str value, bool retr=false) {
    static const int TXTLEN = 49;
    int vlen = value.length(), ilen = key.length();
    int pad = std::max(TXTLEN, ilen+vlen+1)-ilen-vlen;
    key += str(pad, ' ') + value;
    if (retr) key.insert(0, "\r"); else key += '\n';
    return key;
}
/* determine whether the file is a Gzip compressed file */
inline bool is_gzip(const str &filepath) {
    std::ifstream file(filepath, std::ios::binary);
    if (!file.is_open()) return false;
    unsigned char magic[2] = {0, 0};
    file.read(reinterpret_cast<char*>(magic), 2);
    // 1st and 2nd magic of gzip are 0x1F and 0x8B
    return file.gcount() == 2 && magic[0] == 0x1F && magic[1] == 0x8B;
}
/* convert the time into the standard date format of GenBank. */
inline str format_date(const _time_t& tp) {
    static const char* months[] = 
        {"JAN", "FEB", "MAR", "APR", "MAY", "JUN", 
        "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"};
    auto time_t_val = SYSTEM_CLOCK to_time_t(tp);
    std::tm* tm_ptr = std::localtime(&time_t_val);
    if (!tm_ptr) return "01-JAN-1970";  // clock conversion failed: fixed fallback
    std::ostringstream oss;
    oss << std::setw(2) << std::setfill('0') << tm_ptr->tm_mday << '-';
    // Month abbreviations
    oss << months[tm_ptr->tm_mon] << '-';
    oss << (tm_ptr->tm_year + 1900);
    return oss.str();
}
/* convert DNA sequences into the standard format of GenBank */
inline str format_gbk(const str& seq, size_t line_len = 60, size_t group = 10) {
    std::ostringstream out;
    size_t n = seq.size();
    int width = std::max((int)std::to_string(n).size(), 9);

    for (size_t i = 0; i < n; i += line_len) {
        size_t remain = std::min(line_len, n - i);
        str line = seq.substr(i, remain);
        // insert a space every `group` bases (step +1 accounts for inserted spaces)
        for (size_t j = group; j < line.size(); j += group + 1)
            line.insert(j, " ");
        out << std::setw(width) << (i + 1) << ' ' << line << '\n';
    }
    return out.str();
}
/* 
 * Perform vector dot product operation 
 * (AVX and FMA will be invoked when applicable) 
 */
inline float dot_ff(const float *x, const float *y, int dim) {
    #ifdef __AVX__
        __m256 sum = _mm256_setzero_ps();
        int i = 0;

    #ifdef __FMA__
        for (; i + 7 < dim; i += 8) {
            __m256 vx = _mm256_loadu_ps(x + i);
            __m256 vy = _mm256_loadu_ps(y + i);
            sum = _mm256_fmadd_ps(vx, vy, sum);
        }
    #else
        for (; i + 7 < dim; i += 8) {
            __m256 vx = _mm256_loadu_ps(x + i);
            __m256 vy = _mm256_loadu_ps(y + i);
            sum = _mm256_add_ps(sum, _mm256_mul_ps(vx, vy));
        }
    #endif
        // horizontal reduction of the 8 SIMD lanes
        __m128 lo = _mm256_castps256_ps128(sum);
        __m128 hi = _mm256_extractf128_ps(sum, 1);
        __m128 s2 = _mm_add_ps(lo, hi);
        s2 = _mm_add_ps(s2, _mm_movehl_ps(s2, s2));
        s2 = _mm_add_ss(s2, _mm_shuffle_ps(s2, s2, 1));
        float result = _mm_cvtss_f32(s2);
        for (; i < dim; ++i)
            result += x[i] * y[i];

        return result;
    #else
        float result = 0.0f;
        for (int i = 0; i < dim; ++i)
            result += x[i] * y[i];
        return result;
    #endif
}
/* common structs and classes for bioinfomatics */
namespace bioinfo {
    /* a contig, scaffold or replicon in an assembly */
    class scaffold {
      public:
        str                  name;  // name
        std::array<str,2> strands;  // two strands
        size_t                len;  // length
        int            gc_count=0;  // G+C count
        bool                 circ;  // circular
    /*  flt_arr             curve;  // GC-profile */
        scaffold() {};
        scaffold(str &header, str &&origin, bool circ) {
            // FASTA
            if (header[0] == '>') {
                size_t sep = header.find(' ');
                name = sep != str::npos ? header.substr(1, sep-1):header.substr(1);
                if (name.empty()) {
                    std::cerr << "\nError: invalid FASTA header line\n";
                    IOP_ERR;
                }
            // GenBank/EMBL
            } else {
                std::istringstream header_s(header);
                header_s >> name;
                if (!(header_s >> name)) {
                    std::cerr << "\nError: invalid LOCUS/ID line\n";
                    IOP_ERR;
                }
            }
            // Topology
            if (header.find(" circular") != str::npos) this->circ = true;
            else if (header.find(" linear") != str::npos) this->circ = false;
            else this->circ = circ;
            strands[0] = std::move(origin);
            len = strands[0].length();
            strands[1].resize(len);
            // build the reverse complement strand and count G+C
            for (size_t i = 0; i < len; i ++) {
                int c = BASE2ADDR[strands[0][i]];
                strands[1][i] = COMP[strands[0][len-i-1]];
                gc_count += (c == 1 || c == 3);
            }
        }
        /* safe circular sequence indexing function */
        char at(int idx, bool strand) const noexcept {
            idx += (idx < 0) * len;
            idx -= (idx >= (int)len) * len;
            return strands[strand][idx];
        }
        /* convert base to offset address */
        int addr(int idx, char strand) {
            return BASE2ADDR[at(idx, strand)];
        }
        /* convert kmer to a continuous integer encoding */
        int addr(int idx, char strand, int order) {
            static const int digits[] = { 1, 4, 16, 64 };
            if (order == 0) return addr(idx, strand);
            int addr = 0;
            for (int d = order; d >= 0; d --) {
                int a = BASE2ADDR[at(idx+d, strand)];
                if (a < 0) return -1;
                addr += digits[order-d] * a;
            }
            return addr;
        }
        /* determine whether the triplet at this position belongs to a set of codons */
        int codon_type(size_t idx, char strand, str_arr &codons) {
            char codon[3] = { 
                BASE2ADDR[at(idx + 0, strand)], 
                BASE2ADDR[at(idx + 1, strand)], 
                BASE2ADDR[at(idx + 2, strand)]
            };
            for (int i = 0; i < codons.size(); i ++) {
                if (!std::memcmp(codon, codons[i], 3)) return i;
            }
            return -1;
        }
        /* convert to GenBank header */
        str to_header(str &date) {
            std::stringstream oss;
            oss << name << ' ' << std::setw(27-(int)name.length()) << len << " bp"
                << "    DNA" << std::setw(13) << (circ?"circular":"  linear") << " UNA " << date;
            return oss.str();
        }
    };
    /* a genome consist of a group of scaffolds */
    typedef std::vector<scaffold> genome;
    /* open reading frame */
    class orf {
      private:
        str        prot;  // protein sequence
        str         seq;  // nucleotide sequence
      public:
        int         idx;  // unique index
        scaffold  *host;  // host scaffold
        size_t    start;  // start
        idx_arr  starts;  // alternative starts
        int_arr   types;  // codon types of starts
        size_t      end;  // end
        size_t      len;  // length
        int      rstart;  // relative start
        int        rend;  // relative end
        char     strand;  // strand

        float   gc_cont;  // G+C content
        float   G2_cont;  // G2 bias

        float i_score=0;  // init score
        float c_score=0;  // CSVM score
        float r_score=0;  // RBS score
        float t_score=0;  // total score

        int edge_type=0;  // edge (partial 5'-end/3'-end)
        size_t   minlen;  // mininum length
        orf(): host(nullptr) {}
        orf(scaffold *host, idx_arr &&starts, int_arr &&types, size_t end, char strand, size_t minlen): 
        host(host),starts(std::move(starts)),types(std::move(types)),end(end),strand(strand),minlen(minlen)
        {   
            // refine start
            start = this->starts[0];
            if (this->types[0] != 0) for (int i = 1; i < this->types.size(); i ++) {
                size_t new_start = this->starts[i];
                if (this->types[i] == 0 && (new_start-start) <= MIN_LEN && (end-new_start) >= minlen) {
                    start = new_start;
                    break;
                }
            }
            len = end - start;
        }
        /* safe circular sequence indexing function */
        char operator[](size_t idx){
            return host->at(start+idx, strand);
        }
        /* overlapped length between two ORFs */
        int operator-(orf &other) {
            if (other.host == this->host) {
                int t_rst = rstart, t_rnd = rend;
                int o_rst = other.rstart, o_rnd = other.rend;
                // unwrap coordinates that cross the origin of a circular scaffold
                if (t_rst > t_rnd) t_rst -= host->len;
                if (o_rst > o_rnd) o_rst -= host->len;
                if (t_rnd > o_rst && o_rnd > t_rst) 
                    return std::min(t_rnd, o_rnd) - std::max(t_rst, o_rst);
            }
            return 0;
        }
        /* set alternative start as true start */
        int set_start(int idx) {
            int plen = end - starts[idx];
            if (plen >= minlen && starts[idx] != start) {
                start = starts[idx]; 
                len = plen;
                return 0;
            }
            return 1;
        }
        /* initialize edge type, GC content and G2 bias */
        void init(bool partial3=false) {
            // codon type 5 marks a partial 5'-end (no in-frame start codon)
            if (!types.empty() && types[0] == 5) edge_type = 5;
            if (partial3) edge_type += 3;  // mark a partial 3'-end
            // calculate gc content and g2 bias
            int gc = 0, G2 = 0;
            for (int ps = start; ps < end; ps ++) {
                int c = BASE2ADDR[host->at(ps, strand)];
                gc += (int) (c == 1 || c == 3);
                G2 += (int) (c == 3 && (ps-start)%3 == 1);
            }
            gc_cont = (float) gc / len;
            G2_cont = (float) G2 * 3 / len;
        }
        /* convert base to offset address */
        int addr(int idx, int offset, int order) {
            int off = (idx>=0?starts[idx]:start)+offset;
            return host->addr(off, strand, order);
        }
        /* freeze 1-based relative coordinates, unwrapping origin-spanning ORFs */
        orf &froze() {
            const long long L = (long long) host->len;
            long long begin = strand ? -(long long) end : (long long) start;
            begin = ((begin % L) + L) % L;
            rstart = (int) begin + 1;
            rend = (int) ((begin + (long long) len - 1) % L) + 1;  // len == end - start
            return *this;
        }
        /* convert to FASTA header */
        str to_header(int &id) {
            std::ostringstream oss;
            oss << host->name << ':' << rstart << ".." << rend << '(' << (strand?'-':'+') << ')'
                << " ID=ORF" << std::setw(6) << std::setfill('0') << id+1
                << ";score=" << std::fixed << std::setprecision(3) << t_score
                << ";paritial=" << (edge_type>=5) << (edge_type==3 || edge_type==8);
            return oss.str();
        }
        /* get the corresponding protein sequence */
        str translate() {
            if (prot.empty()) {
                int prolen = (int)len/3;
                prot.resize(prolen);
                for (int i = 0; i < prolen; i ++) {
                    // a complete ORF always starts with Met
                    if (!i && edge_type < 5) {
                        prot[0] = 'M';
                        continue;
                    }
                    int off = addr(-1, i*3, 2);
                    prot[i] = (off >= 0 ? trans_tbl[off] : 'X');
                }
            }
            return prot;
        }
        /* get the nucleotide sequence */
        str sequence() {
            if (seq.empty()) {
                str &chr = host->strands[strand];
                if (end <= host->len) seq = chr.substr(start, end-start);
                // wrap around the origin of a circular scaffold
                else {
                    seq = chr.substr(start, host->len - start);
                    seq += chr.substr(0, end - host->len);
                }
            }
            return seq;
        }
        /* convert to GFF line */
        str to_gff() {
            str type = (i_score >= 0 ? "CDS" : "ORF");
            std::ostringstream oss;
            oss << host->name << '\t' << VERSION << '\t' << type << '\t' << rstart << '\t' << rend << '\t'
                << std::fixed << std::setprecision(3) << t_score << '\t' << (strand?'-':'+')
                << "\t0\tID=ORF" << std::setw(6) << std::setfill('0') << idx+1
                << ";partial=" << (edge_type>=5) << (edge_type==3 || edge_type==8);
            return oss.str();
        }
        /* convert to BED line (BED6: chrom, chromStart(0-based), chromEnd(half-open), name, score, strand) */
        str to_bed() {
            int score = (int)(t_score * 1000);
            if (score < 0) score = 0;
            if (score > 1000) score = 1000;
            std::ostringstream oss;
            oss << host->name << '\t' << rstart-1 << '\t' << rend << "\tORF"
                << std::setw(6) << std::setfill('0') << idx+1 << '\t' << score
                << '\t' << (strand?'-':'+') << '\n';
            return oss.str();
        }
        /* convert to GenBank feature */
        str to_gbk(str &table) {
            std::stringstream oss;
            oss << "     " << (i_score >= 0 ? "CDS" : "ORF") << "             ";
            if (strand) oss << "complement(";
            if (rstart > rend) {
                oss << "join(" << rstart << ".." << host->len << ",1.." << rend << ")";
            } else {
                // '<' and '>' mark partial 5'-end and 3'-end respectively
                if (strand?(edge_type==3||edge_type==8):(edge_type>=5)) oss << '<'; 
                oss << rstart << "..";
                if (strand?(edge_type>=5):(edge_type==3||edge_type==8)) oss << '>'; 
                oss << rend;
            }
            oss << (strand ? ")\n" : "\n");
            oss << "                     /locus_tag=ORF" << std::setw(6) << std::setfill('0') << (idx+1) << "\n"
                   "                     /transl_table=" << table << "\n"
                   "                     /codon_start=1\n"
                   "                     /note=\"Derived by protein-coding gene prediction method: " << VERSION << "\"\n"
                   "                     /translation=\"" << translate() << "\"\n";
            return oss.str();
        }
    };
    
    /* orf set manager of a genome */
    class orfs {
      private: 
        const int DIM;
        std::vector<orf> arr;   // inner array
        float *data = nullptr;  // zcurve data

      public:

        orfs(int dim = DIM_S): DIM(dim) {}

        /* override methods of vector */
        size_t size() { return arr.size(); }
        orf &operator[](size_t idx) { return arr.at(idx); }
        void clear() { arr.clear(); delete[] data; data = nullptr;}
        void push_back(orf &item) { arr.push_back(item); }

        /* sort ORFs by GC content and allocate the feature matrix */
        void init() {
            std::sort(arr.begin(), arr.end(), [](orf& a, orf& b) 
                { return a.gc_cont < b.gc_cont; });
            data = NEW float[size()*(DIM+1)];
            if (data == nullptr) {
                std::cerr << MEM_ERR_INFO;
                MEM_ERR;
            }
        }
        /* get zcurve features of ith orf */
        float *row(int i = 0) { return data + i * DIM_A; }
        /* filter out orfs length > minlen and G2 < thres */
        int_arr filter(int min_len, float G2) {
            int_arr indices;
            for (int i = 0; i < size(); i ++) {
                if (arr[i].len >= min_len && arr[i].G2_cont < G2) 
                    indices.push_back(i);
            }
            return indices;
        };
        /* filter out orfs length > minlen and score > min_score */
        int_arr filter(int min_len, float min_score, bool not_edge) {
            int_arr indices;
            for (int i = 0; i < size(); i ++) {
                if (arr[i].len >= min_len && arr[i].i_score > min_score) {
                    if (not_edge && arr[i].edge_type) continue;
                    indices.push_back(i);
                }
            }
            return indices;
        }
        /* get positive genes */
        void yield(str proc, float thres, orfs &genes) {
            for (int i = 0, j = 0; i < size(); i ++) {
                float score = std::max(arr[i].i_score, arr[i].c_score);
                // empirical score calibration for metagenomic mode
                if (proc == "Metagenome") score = std::pow(score, 0.64F);
                if ((arr[i].t_score = score) >= thres) 
                    genes.arr.push_back(arr[i].froze());
            }
            std::sort(genes.arr.begin(), genes.arr.end(), [](orf& a, orf& b) {
                if (a.host < b.host) return true;
                else if (a.host > b.host) return false;
                else return a.rend < b.rend;
            });
            for (int i = 0; i < genes.size(); i ++) genes[i].idx = i;
        }
        /* concat RBS score to Z-curve features */
        orfs &concat() {
            // store the RBS score in the last column (index DIM) of each row
            for (int i = 0; i < size(); i ++)
                data[i*DIM+DIM+i] = arr[i].r_score;
            return *this;
        }
        /* seek all the orfs in the sequence */
        void locate_orfs(scaffold &seq, str_arr &start_types, str_arr &stop_types, size_t minlen) {
            std::vector<orf> cache;
            for (int strand = 0; strand < 2; strand ++) {
                size_t count = 0;
                for (int phase = 0; phase < 3; phase ++) {
                    idx_arr starts; int_arr types;
                    size_t pcount = 0, send = seq.len;
                    // circular: scan twice around to catch origin-spanning ORFs
                    if (seq.circ) send += seq.len;
                    else {
                        // non-circular: seed a partial 5' start (type 5) at frame origin
                        starts.push_back(phase);
                        types.push_back(5);
                    }
                    
                    for (size_t ps = phase; ps < send; ps += 3) {
                        int match = seq.codon_type(ps, strand, start_types);
                        if (match > -1) {
                            if (starts.empty()) { if (ps >= seq.len) break; }
                            else if (ps == phase) {
                                starts.pop_back();
                                types.pop_back();
                            }
                            starts.push_back(ps);
                            types.push_back(match);
                        } else if (starts.size() && seq.codon_type(ps, strand, stop_types) > -1) {
                            arr.emplace_back(&seq, std::move(starts), std::move(types), ps+3, strand, minlen);
                            if (arr.back().len < minlen || arr.back().len > 50000) arr.pop_back();
                            else { arr.back().init(); pcount ++; }
                            starts.clear(); types.clear();
                        }
                    }
                    if (seq.circ) {
                        if (pcount > 0) {
                            orf &back = arr.back();
                            if (back.end > seq.len) cache.push_back(back);
                        }
                        count += pcount;
                    } else if (starts.size()) {
                        arr.emplace_back(&seq, std::move(starts), std::move(types), seq.len, strand, minlen);
                        if (arr.back().len < minlen || arr.back().len > 50000) arr.pop_back();
                        else arr.back().init(true);
                    }
                }
                
                // remove duplicate ORFs found twice around a circular scaffold
                if (seq.circ) for (const auto& cached : cache) {
                    for (int i = size()-count; i < size(); i ++) {
                        if (cached.end - (*this)[i].end == seq.len) {
                            arr.erase(arr.begin() + i);
                            --count;
                            break;
                        }
                    }
                }
            }
        }

        ~orfs() { delete[] data; }

        /* collect all pairs of genes overlapping by at least min_olen */
        static orfs overlap(int min_olen, orfs &genes) {
            orfs pairs;
            for (int i = 0; i < genes.size(); i ++) {
                for (int j = i + 1; j < genes.size(); j ++) {
                    int olen = genes[i] - genes[j];
                    if (olen >= min_olen) {
                        pairs.arr.push_back(genes[i]);
                        pairs.arr.push_back(genes[j]);
                    }
                }
            }
            return pairs;
        }
    };
}

namespace zcurve {
    /* mono-nucleotide Z-curve transform per codon position (9 features) */
    inline void mono_trans(bioinfo::orf &orf, float *params) noexcept {
        float counts[PHASE][3] = {{0.0f}};
        int len = orf.len / 3 * 3;  // truncate to complete codons
            
        for (int i = 0; i < len; i ++) {
            int p = i % PHASE;

            counts[p][X] += Z_COORD[orf[i]][X];
            counts[p][Y] += Z_COORD[orf[i]][Y];
            counts[p][Z] += Z_COORD[orf[i]][Z];
        }
            
        for (int p = 0; p < PHASE; p ++, params += 3) {
            params[X] = counts[p][X] / len * PHASE;
            params[Y] = counts[p][Y] / len * PHASE;
            params[Z] = counts[p][Z] / len * PHASE;
        }
    }

    /* di-nucleotide Z-curve transform per codon position (36 features) */
    inline void di_trans(bioinfo::orf &orf, float *params) noexcept {
        float counts[PHASE][4][3] = {{{0.0f}}};
        int len = orf.len / 3 * 3;

        for (int i = 0; i < len-1; i++) {
            int p = i % PHASE;
            for (int b = 0; b < 4; b ++) {
                counts[p][b][X] += ONE_HOT[orf[i]][b] * Z_COORD[orf[i + 1]][X];
                counts[p][b][Y] += ONE_HOT[orf[i]][b] * Z_COORD[orf[i + 1]][Y];
                counts[p][b][Z] += ONE_HOT[orf[i]][b] * Z_COORD[orf[i + 1]][Z];
            }
        }
            
        for (int p = 0; p < PHASE; p ++)
            for (int b = 0; b < 4; b ++, params += 3) {
                params[X] = counts[p][b][X] / len * PHASE;
                params[Y] = counts[p][b][Y] / len * PHASE;
                params[Z] = counts[p][b][Z] / len * PHASE;
            }
    }

    /* tri-nucleotide Z-curve transform per codon position (144 features) */
    inline void tri_trans(bioinfo::orf &orf, float *params) noexcept {
        float counts[PHASE][4][4][3] = {{{{0.0f}}}};
        int len = orf.len / 3 * 3;

        for (int i = 0; i < len-2; i ++) {
            int p = i % PHASE;
            for (int s = 0; s < 4; s ++)
            for (int b = 0; b < 4; b ++) {
                counts[p][s][b][X] += ONE_HOT[orf[i]][s] * ONE_HOT[orf[i + 1]][b] * 
                                    Z_COORD[orf[i + 2]][X];
                counts[p][s][b][Y] += ONE_HOT[orf[i]][s] * ONE_HOT[orf[i + 1]][b] * 
                                    Z_COORD[orf[i + 2]][Y];
                counts[p][s][b][Z] += ONE_HOT[orf[i]][s] * ONE_HOT[orf[i + 1]][b] * 
                                    Z_COORD[orf[i + 2]][Z];
            }
        }

        for (int p = 0; p < PHASE; p ++)
        for (int s = 0; s < 4; s ++)
        for (int b = 0; b < 4; b ++, params += 3) {
            params[X] = counts[p][s][b][X] / len * PHASE;
            params[Y] = counts[p][s][b][Y] / len * PHASE;
            params[Z] = counts[p][s][b][Z] / len * PHASE;
        }
    }
    
    /* compute Z-curve features of all orders for every ORF (OpenMP parallelized) */
    inline void encode(bioinfo::orfs &orfs, int max_order=3) {
        static void (*k_trans[3])(bioinfo::orf &orf, float *params) = {
            mono_trans, di_trans, tri_trans
        };
        static int dims[] = { 9, 36, 144 };
        const int count = (int) orfs.size();
    #ifdef _OPENMP
        #pragma omp parallel for
    #endif
        for (int i = 0; i < count; i ++) {
            float *head = orfs.row(i);
            for (int j = 0; j < max_order; j ++) {
                (*k_trans[j])(orfs[i], head);
                head += dims[j];
            }
        }
    }
/*
    inline void gc_profile(bioinfo::genome &scaffolds, float gc_cont) {
        float k = 1.0F - 2 * gc_cont;
        for (auto &scaffold : scaffolds) {
            scaffold.curve.resize(scaffold.len);
            float cnt = 0.0F;
            for (int i = 0; i < scaffold.len; i ++) {
                cnt += Z_COORD[scaffold.at(i, 0)][Z];
                scaffold.curve[i] = cnt - k * i;
            }
        }
    }
*/
    /* dump Z-curve features and labels of selected ORFs to a binary file */
    inline void save_bin(str &file, bioinfo::orfs &orfs, int_arr &indices) {
        std::ofstream handle(file, std::ios::binary);
        if (!handle.is_open()) {
            std::cerr << "\nError: failed to open output file for debug\n";
            IOP_ERR;
        }
        int n_rep = (int) indices.size();
        handle.write((char*)&n_rep, sizeof(int));

        /* write the full 189D Z-curve features (mono 9D + di 36D + tri 144D) */
        for (int i = 0; i < n_rep; i ++) {
            float *row = orfs.row(indices[i]);
            handle.write((char*) row, sizeof(float)*DIM_S);
        }

        for (int i = 0; i < n_rep; i ++) {
            int label = (int) (orfs[indices[i]].t_score >= 0.5);
            handle.write((char*) &label, sizeof(int));
        }
    }
}

namespace model {
    /* single hidden-layer MLP for initial scoring (pre-trained, embedded in META) */
    class mlp {
      private:
        // per-model param count: 2*189 scaler + 100*189+100 hidden + 100+1 output
        static const int N_PARAMS=19479, N_HIDDEN=100;
        const float *MODELS[N_MODELS];

        static float sigmoid(float v) { 
            return 1.0F/(1.0F+std::exp(-v)); 
        }

        /* standardize features with the pre-computed mean and std */
        static float *std_scale(float *data, int size, const float *scale) {
            float *cache = NEW float[size*DIM_S];
            const float *stds = scale + DIM_S;
            if (cache == nullptr) {
                std::cerr << MEM_ERR_INFO;
                MEM_ERR;
            }
            for (int i = 0; i < size; i++) {
                int j;
                float *p = cache + i*DIM_S;
                for (j = 0; j < DIM_S; j++, p ++) {
                    *p = (data[i*DIM_A+j]-scale[j])/stds[j];
                }
            }
            return cache;
        }
        mlp() {
            for (int i = 0; i < N_MODELS; ++i)
                MODELS[i] = META + i * N_PARAMS;
        }
      public:
        /* run the index-th heuristic model on orfs[offset, offset+size) */
        static void predict(int index, bioinfo::orfs &orfs, int offset, int size) {
            // lazily-initialized singleton holding all model parameters
            static mlp *heuristic = NEW mlp();
            if (!size) return;
            
            const float *scale = heuristic->MODELS[index];
            float *scaled = std_scale(orfs.row(offset), size, scale);

            const float *hid_ws = scale + 2*DIM_S,          // hidden w
                        *hid_bs = hid_ws + N_HIDDEN*DIM_S,  // hidden b
                        *out_ws = hid_bs + N_HIDDEN;        // output w
            const float out_b = *(out_ws + N_HIDDEN);       // output b
        #ifdef _OPENMP
            #pragma omp parallel for
        #endif
            for (int i = 0; i < size; i ++) {
                float *x = scaled + i*DIM_S;
                // Hidden Layer (189, ) -> (100, )
                float hid_out[N_HIDDEN] = { 0.0F };
                for (int j = 0; j < N_HIDDEN; j ++) {
                    const float *hid_w = hid_ws + j*DIM_S;
                    // hid_out = hid_w * x + hid_b
                    hid_out[j] = dot_ff(x, hid_w, DIM_S)+hid_bs[j];
                    // ReLU activation function
                    hid_out[j] = std::max(0.0F, hid_out[j]);
                }
                // out = out_w * hid_out + out_b
                float proba = dot_ff(hid_out, out_ws, N_HIDDEN) + out_b;
                // Sigmoid activation function (divide by 2 as temperature scaling)
                orfs[i+offset].i_score = sigmoid(proba / 2.0F);
            }
            delete[] scaled;
        }
    };

    /* positional Markov model (TriTISA+) for start-codon/RBS scoring */
    class pmm {
      private:
        bioinfo::orfs &orfs;
        // scoring window: 50 bp upstream to 15 bp downstream of the start codon
        static const int UPSREAM = 50, DWSTREAM = 15;
        int order, dim, pwm_s, max_cand = 6;
        float *pwm = nullptr, uf, df;
        str table;
        /* normalize each PWM column and convert to log-probabilities */
        void norm_log() {
            for (int j = 0; j < pwm_s*3; j += 4) {
                float sum = pwm[j] + pwm[j+1] + pwm[j+2] + pwm[j+3];
                for (int k = 0; k < 4; k ++)
                    pwm[j+k] = log(pwm[j+k] / sum);
            }
        }
      public:
        pmm(bioinfo::orfs &orfs, str &table): orfs(orfs), table(table) {}
        void reset(int order) {
            this->order = order;
            dim = 1<<(2*(order+1));  // 4^(order+1), number of distinct k-mers
            pwm_s = dim*(UPSREAM+DWSTREAM);
            delete[] pwm;
            pwm = NEW float[pwm_s*3];
            if (pwm == nullptr) {
                std::cerr << MEM_ERR_INFO;
                MEM_ERR;
            }
        }
        /* train the PWMs from the start-codon regions of seed genes */
        void train(int_arr &seeds, str_arr &starts) {
            flt_arr bkg(4, 0);
            std::fill(pwm, pwm+pwm_s*3, 1.0F);

            const int n_seeds = (int) seeds.size();
            for (size_t i = 0; i < n_seeds; i ++) {
                bioinfo::orf &gene = orfs[seeds[i]];
                if (gene.len < 300) continue;
                int t_start = gene.start, start;

                const int stop = std::min((int)gene.starts.size(), max_cand);
                for (int j = 0; j < stop; j ++) {
                    float *wm = nullptr, *ip = nullptr;
                    if ((start = gene.starts[j]) < t_start) continue;

                    // pwm layout: [background | true starts | alternative starts]
                    wm = pwm + ((start!=t_start)+1) * pwm_s;
                    for (int k = -UPSREAM; k < DWSTREAM - order; k ++) {
                        int addr = gene.addr(j, k, order);
                        if (addr >= 0) wm[(k+UPSREAM)*dim+addr] ++;
                    }
                }
                
                // estimate background frequencies around start-like codons in the
                // upstream non-coding region (150..3 bp before the true start)
                for (int j = t_start - 150; j < t_start - 3; j ++) {
                    if (gene.host->codon_type(j, gene.strand, starts) >= 0) {
                        for (int k = -UPSREAM; k < DWSTREAM - order; k ++) {
                            int addr = gene.host->addr(j+k, gene.strand, order);
                            if (addr >= 0) pwm[(k+UPSREAM)*dim+addr]++;
                            if ((addr = gene.host->addr(j+k, gene.strand)) >= 0) bkg[addr]++;
                        }
                    }
                }
            }
            
            float C = bkg[0] + bkg[1] + bkg[2] + bkg[3];
            if (C > 0.0F) for (int i = 0; i < 4; i ++) bkg[i] /= C;
            
            float t = bkg[0], c = bkg[1], a = bkg[2], g = bkg[3];

            float s = a*t*g;
            // support for standard
            if (table != "1") s += g*t*g+t*t*g;
            float e = t*a*a;
            // support for mycoplasma, etc.
            if (table != "4" && table != "25") e += t*g*a;
            if (table != "15" && table != "16") e += t*a*g;

            // p/q: probability that a codon is/isn't a valid start under this
            // genetic code table; uf/df: expected numbers of up/downstream
            // candidate starts from the geometric distribution
            float p = s / (s + e), q = e / (s + e), P = 0.0F;

            int i = 0; float qpi = q; uf = 0;
            while (true) {
                P += qpi;
                if (P > 0.999 || i > 100) {
                    max_cand = i + 1;
                    break;
                } else {
                    uf += i * qpi;
                    qpi *= p;
                    i ++;
                }
            }
            df = (float) max_cand - 1 - uf;

            this->norm_log();
        }
        /* re-select start codons by PWM scores; return the fraction unchanged */
        float revise(int_arr &seeds) {
            const int n = (int) seeds.size();
            float *up = pwm, *ts = pwm + pwm_s, *dw = pwm+2*pwm_s;
            int unc = 0;
            for (size_t i = 0; i < n; i ++) {
                bioinfo::orf &gene = orfs[seeds[i]];
                if (gene.edge_type != 0) continue;
                int n_alter = (int)gene.starts.size();
                float max_score = 0.0F, max_idx = -1;
                for (int j = 0; j < std::min(n_alter, max_cand); j ++) {
                    float pu = 0.0F, pt = 0.0F, pd = 0.0F;
                    for (int k = -UPSREAM; k < DWSTREAM - order; k ++) {
                        int addr = gene.addr(j, k, order);
                        if (addr < 0) continue;  // skip ambiguous bases (e.g. N)
                        pu += up[(k+UPSREAM)*dim+addr];
                        pt += ts[(k+UPSREAM)*dim+addr];
                        pd += dw[(k+UPSREAM)*dim+addr];
                    }
                    // posterior-like score combining upstream, true-start and
                    // downstream (alternative-start) evidence
                    float score = 1.0F / (1 + exp(pu-pt) * uf + exp(pd-pt) * df);
                    if (score > max_score) { max_score = score; max_idx = j; }
                }
                if (max_idx >= 0) {
                    unc += gene.set_start(max_idx);
                    gene.r_score = max_score;
                }
            }
            return (float) unc / n;
        }
        void revise_all() {
            int_arr indices(orfs.size());
            std::iota(indices.begin(), indices.end(), 0);
            revise(indices);
        }
        /* save the TriTISA+ model to a binary file */
        void dump(str &file) {
            if (file.empty()) return;
            std::ofstream outfile(file, std::ios::binary);
            if (!outfile.is_open()) {
                std::cerr << "\nError: failed to open " << file << '\n';
                IOP_ERR;
            }
            outfile.write((char*) &order, sizeof(int));
            outfile.write((char*) &uf, sizeof(float));
            outfile.write((char*) &df, sizeof(float));
            outfile.write((char*) pwm, pwm_s*3*sizeof(float));
            if (!outfile) {
                std::cerr << "\nError: failed to write TriTISA+ model to " << file << '\n';
                IOP_ERR;
            }
        }
        /* load the TriTISA+ model from a binary file */
        void load(str &file) {
            if (file.empty()) return;
            std::ifstream infile(file, std::ios::binary);
            if (!infile.is_open()) {
                std::cerr << "\nError: failed to open " << file << '\n';
                IOP_ERR;
            }
            infile.read((char*) &order, sizeof(int));
            // supported Markov orders are 0..3 (bounded by the 4^4 k-mer space)
            if (!infile || order < 0 || order > 3) {
                std::cerr << "\nError: invalid TriTISA+ model file " << file << '\n';
                IOP_ERR;
            }
            this->reset(order);
            infile.read((char*) &uf, sizeof(float));
            infile.read((char*) &df, sizeof(float));
            infile.read((char*) pwm, pwm_s*3*sizeof(float));
            if (!infile) {
                std::cerr << "\nError: truncated TriTISA+ model file " << file << '\n';
                IOP_ERR;
            }
            max_cand = (int) (uf + df + 1);
            if (max_cand < 1) max_cand = 1;
        }
        ~pmm() { delete[] pwm; }
    };

    /* RBF-SVM classifier for coding-potential scoring (libsvm-based) */
    class svm {
      private:
        svm_parameter param = {
            C_SVC,    /* svm_type     */
            RBF,      /* kernel_type  */
            3,        /* degree       */
            1.0/DIM_S,/* gamma        */
            0.0,      /* coef0        */
            200,      /* cache_size   */
            1e-3,     /* eps          */
            1.0,      /* C            */
            0,        /* nr_weight    */
            nullptr,  /* weight_label */
            nullptr,  /* weight       */
            0.0,      /* nu           */
            0.0,      /* p            */
            true,     /* shrinking    */
            false     /* probability  */
        };
        svm_model *model;    // svm model
        svm_problem prob;    // svm problem

        bioinfo::orfs &orfs; // total orfs

        flt_arr mins, maxs;  // min-max scaler
        flt_arr data, y;     // dataset

        /* apply min-max scaler to current dataset */
        void transform(float *p_data, const int n) {
            for (int i = 0; i < n; ++ i) {
                float *row = p_data + i * DIM_A;
                for (int j = 0; j < DIM_A; j ++) {
                    float intv = maxs[j] - mins[j];
                    if (intv > 0.0F) row[j] = (row[j]-mins[j]) / intv;
                }
            }
        }

      public:
        svm(bioinfo::orfs &orfs): orfs(orfs), model(nullptr) {
            data.reserve(orfs.size()*DIM_A);
            mins.resize(DIM_A), maxs.resize(DIM_A);

        }

        svm &set_C(float val) { this->param.C = val; return *this; }

        /* pick out training set of svm */
        svm &sample(float up, float down, float eps=1E-5) {
            y.clear();
            
            std::fill(mins.begin(), mins.end(), 1.0F);
            std::fill(maxs.begin(), maxs.end(), 0.0F);

            if (model) svm_free_model_content(model);

            int_arr X;
            for (int i = 0; i < orfs.size(); i ++) {
                float s = std::max(orfs[i].i_score, orfs[i].c_score);
                if (s > up || (s <= down && s > eps)) {
                    X.push_back(i);
                    y.push_back(s > up ? 1.0F : -1.0F);
                }
            }

            /* construct svm problem */
            prob.l = (int) X.size();
            prob.x = NEW float*[prob.l];
            prob.y = y.data();

            data.resize(prob.l*DIM_A);

            for (int i = 0; i < prob.l; i ++) {
                float *row = data.data() + i * DIM_A;
                std::memcpy(row, orfs.row(X[i]), DIM_A*sizeof(float));
                for (int j = 0; j < DIM_A; j ++) {
                    mins[j] = std::min(row[j], mins[j]);
                    maxs[j] = std::max(row[j], maxs[j]);
                }
                prob.x[i] = row;
            }

            transform(data.data(), prob.l);
            return *this;
        }

        /* train svm model */
        svm &train() {
            if (prob.l == 0) return *this;
            
            /* calculate gamma */
            const int total = prob.l * DIM_A;
            float sum = 0.0F, sum2 = 0.0F;
            
            for (int i = 0; i < prob.l; i++) {
                for (int j = 0; j < DIM_A; ++j) {
                    float v = prob.x[i][j];
                    sum += v;
                    sum2 += v * v;
                }
            }

            float mean = sum / total;
            float mean_sq = mean * mean;
            float var = (sum2 / total) - mean_sq;
            if (var <= 0) var = 1E-12;
            param.gamma = 1.0F / (DIM_A * var);

            /* train model */
            model = svm_train(&prob, &param, DIM_A);

            delete[] prob.x; prob.x = nullptr;

            if (model == nullptr) {
                std::cerr << "Error: failed to build RBF-SVM classifier\n";
                SUB_ERR;
            }
            return *this;
        }

        /* prediction on all the orfs */
        void predict() {
            const int n = (int) orfs.size();
            if (n == 0 || !model) return;

            flt_arr samples(n*DIM_A);
            std::memcpy(samples.data(), orfs.row(), n*DIM_A*sizeof(float));
            transform(samples.data(), n);

            #ifdef _OPENMP
                #pragma omp parallel for
            #endif
            for (int i = 0; i < n; i ++) {
                float *row = samples.data() + i * DIM_A;
                orfs[i].c_score = svm_predict_score(model, row, DIM_A);
            }
        }

        /* save the SVM model (scaler, support vectors, coefficients, intercept) */
        void dump(str &file) {
            if (file.empty()) return;
            if (model == nullptr) {
                Debug() << "Warning: no SVM model was trained; skip dumping\n";
                return;
            }
            std::ofstream outfile(file, std::ios::binary);
            if (!outfile.is_open()) {
                std::cerr << "\nError: failed to open " << file << '\n';
                IOP_ERR;
            }
            outfile.write((char *) mins.data(), DIM_A * sizeof(float));
            outfile.write((char *) maxs.data(), DIM_A * sizeof(float));

            float gamma = model->param.gamma;
            outfile.write((char*) &gamma, sizeof(float));

            int n_sv = model->l;
            outfile.write((char*) &n_sv, sizeof(int));

            for (int i = 0; i < n_sv; i ++) {
                float *sv = model->SV[i];
                outfile.write((char*) sv, DIM_A * sizeof(float));
            }

            for (int i = 0; i < n_sv; i ++) {
                float coef = model->sv_coef[0][i];
                outfile.write((char*) &coef, sizeof(float));
            }

            float intercept = model->rho[0];
            outfile.write((char*) &intercept, sizeof(float));
            if (!outfile) {
                std::cerr << "\nError: failed to write svm model to " << file << '\n';
                IOP_ERR;
            }
        }
        /* load the SVM model from a binary file */
        void load(str &file) {
            std::ifstream infile(file, std::ios::binary);
            if (!infile.is_open()) {
                std::cerr << "\nError: failed to open " << file << '\n';
                IOP_ERR;
            }
            model = NEW svm_model();  // value-initialized: all fields zeroed
            if (model == nullptr) {
                std::cerr << MEM_ERR_INFO;
                MEM_ERR;
            }
            infile.read((char *) mins.data(), DIM_A * sizeof(float));
            infile.read((char *) maxs.data(), DIM_A * sizeof(float));

            float gamma;
            infile.read((char*) &gamma, sizeof(float));

            int n_sv;
            infile.read((char*) &n_sv, sizeof(int));
            // reject truncated headers and implausible support-vector counts
            if (!infile || n_sv <= 0 || n_sv > 10000000) {
                std::cerr << "\nError: invalid SVM model file " << file << '\n';
                IOP_ERR;
            }
            param.gamma = gamma;
            param.kernel_type = RBF;
            model->param = param;
            model->l = n_sv;
            model->nr_class = 2;  // binary classifier (coding / non-coding)

            data.resize(n_sv * DIM_A);

            infile.read((char*) data.data(), n_sv * DIM_A * sizeof(float));

            // svm_free_model_content() releases these with free(); allocate with malloc
            model->SV = (float**) std::malloc(n_sv * sizeof(float*));
            model->sv_coef = (float**) std::malloc(sizeof(float*));
            model->rho = (float*) std::malloc(sizeof(float));
            if (model->SV && model->sv_coef && model->rho)
                model->sv_coef[0] = (float*) std::malloc(n_sv * sizeof(float));
            if (!model->SV || !model->sv_coef || !model->rho || !model->sv_coef[0]) {
                std::cerr << MEM_ERR_INFO;
                MEM_ERR;
            }
            for (int i = 0; i < n_sv; i ++) {
                model->SV[i] = data.data() + i * DIM_A;
            }

            infile.read((char*) model->sv_coef[0], n_sv * sizeof(float));
            infile.read((char*) model->rho, sizeof(float));
            if (!infile) {
                std::cerr << "\nError: truncated SVM model file " << file << '\n';
                IOP_ERR;
            }
        }
        ~svm() { svm_free_model_content(model); }
    };
}

namespace utils {
    inline bool file_exists(str &path) {
        if (path.empty()) return false;
        struct stat buffer;
        return (stat(path.c_str(), &buffer) == 0);
    }

    inline str filename(str path) {
        size_t pos = path.find_last_of("/\\");
        if (pos == str::npos) return path;
        else return path.substr(pos + 1);
    }

    /* whether the line starts with the given keyword (structural lines only) */
    inline bool starts_with(const str &line, const char *prefix) {
        size_t n = std::strlen(prefix);
        return line.size() >= n && line.compare(0, n, prefix) == 0;
    }

    /* parse FASTA/GenBank/EMBL content from a stream into a genome */
    inline void read_stream(std::istream &handle, bioinfo::genome &genome, bool circ) {
        str header;
        // skip empty lines
        while(std::getline(handle, header))
            if (!header.empty()) break;
        if (header.empty()) {
            std::cerr << "\nError: empty content \n";
            IOP_ERR;
        }
        str line, buffer;
        // FASTA files should start with '>'
        if (header.at(0) == '>') {
            while (std::getline(handle, line)) {
                if (line.empty()) continue;
                else if (line.at(0) == '>') {
                    genome.emplace_back(header, std::move(buffer), circ);
                    header = std::move(line);
                    buffer.clear();
                } else {
                    if (!std::isalpha((unsigned char) line.back())) line.pop_back();
                    buffer.append(std::move(line));
                }
                line.clear();
            }
        // GenBank files should start with 'LOCUS ... '
        } else if (starts_with(header, "LOCUS")) {
            bool in_origin = false;
            while (std::getline(handle, line)) {
                if (line.empty()) continue;
                else if (starts_with(line, "ORIGIN")) in_origin = true;
                else if (starts_with(line, "//")) {
                    if (!buffer.empty()) {
                        genome.emplace_back(header, std::move(buffer), circ);
                        buffer.clear();
                    }
                    in_origin = false;
                } else if (starts_with(line, "LOCUS")) {
                    if (!buffer.empty()) {
                        genome.emplace_back(header, std::move(buffer), circ);
                        buffer.clear();
                    }
                    in_origin = false; 
                    header = std::move(line);
                } else if (in_origin) {
                    for (char c : line) if (std::isalpha((unsigned char) c)) buffer += (char) std::toupper((unsigned char) c);
                }
            }
        // EMBL files should start with 'ID   '
        } else if (starts_with(header, "ID   ")) {
            bool in_sequence = false;
            while (std::getline(handle, line)) {
                if (line.empty()) continue;
                else if (starts_with(line, "SQ   ")) in_sequence = true;
                else if (starts_with(line, "//")) {
                    if (!buffer.empty()) {
                        genome.emplace_back(header, std::move(buffer), circ);
                        buffer.clear();
                    }
                    in_sequence = false;
                } else if (starts_with(line, "ID   ")) {
                    if (!buffer.empty()) {
                        genome.emplace_back(header, std::move(buffer), circ);
                        buffer.clear();
                    }
                    in_sequence = false;
                    header = std::move(line);
                } else if (in_sequence) {
                    for (char c : line) if (std::isalpha((unsigned char) c)) buffer += (char) std::toupper((unsigned char) c);
                }
            }
        } else {
            std::cerr << "\nError: unsupported file type (options: GenBank, EMBL, FASTA)\n";
            IOP_ERR;
        }
        if (!buffer.empty()) genome.emplace_back(header, std::move(buffer), circ);
    }

    /* write predicted genes to a stream in the requested output format */
    inline void write_stream(std::ostream &handle, bioinfo::orfs &genes, str date, str table, str format) {
        bioinfo::scaffold *last = nullptr;
        if (format == "gff") {
            handle << "##gff-version 3\n# trans_tbl=" << table << '\n';
            for (int idx = 0; idx < genes.size(); idx ++) {
                if (last != genes[idx].host) {
                    last = genes[idx].host;
                    handle << "# " << last->to_header(date) << '\n';
                }
                handle << genes[idx].to_gff() << '\n';
            }
        } else if (format == "faa") {
            for (int idx = 0; idx < genes.size(); idx ++) {
                handle << '>' << genes[idx].to_header(idx) << '\n';
                handle << genes[idx].translate() << '\n';
            }
        } else if (format == "faa-tmp") {
            for (int idx = 0; idx < genes.size(); idx ++) {
                handle << '>' << idx << '\n' << 
                genes[idx].translate() << '\n';
            }
        } else if (format == "fna") {
            for (int idx = 0; idx < genes.size(); idx ++) {
                handle << '>' << genes[idx].to_header(idx) << ";gc_cont=" << std::fixed 
                       << std::setprecision(3) << genes[idx].gc_cont << '\n';
                handle << genes[idx].sequence() << '\n';
            }
        } else if (format == "gbk" || format == "gbk-full") {
            for (int idx = 0; idx < genes.size(); idx ++) {
                if (last != genes[idx].host) {
                    if (last != nullptr) {
                        handle << "ORIGIN\n";
                        if (format == "gbk-full") {
                            handle << format_gbk(last->strands[0]);
                        }
                        handle << "//\n";
                    }
                    last = genes[idx].host;
                    handle << "LOCUS       " << last->to_header(date) << '\n'
                        << "DEFINITION  " << last->name << "\n"
                        << "FEATURES             Location/Qualifiers\n";
                }
                handle << genes[idx].to_gbk(table);
            }
            if (last != nullptr) {  // no genes predicted -> nothing to close
                handle << "ORIGIN\n";
                if (format == "gbk-full") {
                    handle << format_gbk(last->strands[0]);
                }
                handle << "//\n";
            }
        } else if (format == "bed") {
            for (int idx = 0; idx < genes.size(); idx ++) {
                if (last != genes[idx].host) {
                    last = genes[idx].host;
                    handle << "# " << last->to_header(date) << '\n';
                }
                handle << genes[idx].to_bed();
            }
        } else {
            std::cerr << "Error: unknown format '" << format << "'\n"
                    << "Supported output formats are gff, bed, gbk, gbk-full\n";
            ARG_ERR;
        }
    }

    /* read sequences from stdin ("-"), a gzip file, or a plain file */
    #ifdef ZLIB
    /* streambuf that decompresses a gzip file chunk by chunk (no full-file buffering) */
    class gzinbuf : public std::streambuf {
      public:
        gzinbuf(const str &path) {
            file = gzopen(path.c_str(), "rb");
            if (file == nullptr) {
                std::cerr << "\nError: cannot open " << path << "\n";
                IOP_ERR;
            }
        }
        ~gzinbuf() { if (file) gzclose(file); }
      protected:
        int_type underflow() override {
            int n = gzread(file, buff, BUFF_SIZE);
            if (n <= 0) {
                int err_code = Z_OK;
                gzerror(file, &err_code);
                if (n < 0 || (err_code != Z_OK && err_code != Z_STREAM_END)) {
                    std::cerr << "\nError: failed to decompress the input file\n";
                    IOP_ERR;
                }
                return traits_type::eof();
            }
            setg(buff, buff, buff + n);
            return traits_type::to_int_type(*gptr());
        }
      private:
        gzFile file;
        char buff[BUFF_SIZE];
    };
    #endif

    inline void read_source(str &filename, bioinfo::genome &genome, bool circ) {
        if (filename == "-") {
            Debug() << "(note: reading content from the standard input :)\n";
            read_stream(std::cin, genome, circ);
        }
    #ifdef ZLIB
        else if (is_gzip(filename)) {
            gzinbuf buf(filename);
            std::istream handle(&buf);
            read_stream(handle, genome, circ);
        }
    #endif
        else {
            std::ifstream handle(filename);
            if (!handle.is_open()) {
                std::cerr << "\nError: cannot open " << filename << "\n";
                IOP_ERR;
            }
            read_stream(handle, genome, circ);
        }
    }

    /* write results to stdout ("-") or to an output file */
    inline void write_output(
        str &filename, bioinfo::orfs &genes, 
        _time_t start_t, str table, str format
    ) {
        str date = format_date(start_t);
        if (filename == "-") write_stream(std::cout, genes, date, table, format);
        else {
            std::ofstream handle(filename);
            if (!handle.is_open()) {
                std::cerr << "\nError: failed to open " << filename << '\n';
                IOP_ERR;
            }
            return write_stream(handle, genes, date, table, format);
        }
    }
}