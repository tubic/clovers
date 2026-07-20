#include "Common.hpp"
#include "svm.hpp"
#include <array>
#include <random>

#define DIM_K  144
#define S_LEN  550

/* ignore */
const unsigned char META[] = {0};

char trans_tbl[] = {
    'F', 'F', 'L', 'L', 'S', 'S', 'S', 'S',
    'Y', 'Y', '*', '*', 'C', 'C', '*', 'W',
    'L', 'L', 'L', 'L', 'P', 'P', 'P', 'P',
    'H', 'H', 'Q', 'Q', 'R', 'R', 'R', 'R', 
    'I', 'I', 'I', 'M', 'T', 'T', 'T', 'T',
    'N', 'N', 'K', 'K', 'S', 'S', 'R', 'R', 
    'V', 'V', 'V', 'V', 'A', 'A', 'A', 'A',
    'D', 'D', 'E', 'E', 'G', 'G', 'G', 'G'
};
/* start codons */
char ATG[]{2,0,3}, GTG[]{3,0,3}, TTG[]{0,0,3};
/* stop codons  */
char TAA[]{0,2,2}, TAG[]{0,2,3}, TGA[]{0,3,2};

class Kmeans {
    bioinfo::orfs &orfs;    // all orfs from the genome
    flt_arr c;              // k-means centroids
    flt_arr icache, jcache; // square norm cache

    /* sample pointer */
    float *irow(int idx) {
        static int off = DIM_S - DIM_K;
        return orfs.row(seeds[idx]) + off;
    }
    /* centroid pointer */
    float *jrow(int idx) {
        return c.data() + idx * DIM_K;
    }
    /* calculate distance between sample and centroid */
    float dist_sq(int i, int j) {
        return icache[i] + jcache[j] - 2 * dot_ff(irow(i), jrow(j), DIM_K);
    }
    /* add sample to centroid (vector) */
    void add(int i, int j) {
        float *i_row = irow(i), *j_row = jrow(j);
        for (int m = 0; m < DIM_K; m ++) {
            j_row[m] += i_row[m];
        }
    }

  public:
    /* n_clusters */
    int k, n;
    /* indices of samples and their labels */
    int_arr inits, seeds, labels;

    Kmeans(bioinfo::orfs &orfs, int k): orfs(orfs), k(k) {
        c = flt_arr(DIM_K * k, 0.0F);
        jcache.resize(k);
        zcurve::encode(orfs);
    }

    /* pick out initial seed orfs */
    void init(float gc_cont) {
        // get long orfs with expected G2 bias
        inits = orfs.filter(S_LEN, 0.2F * gc_cont + 0.12F);
        n = (int) inits.size();
        // remove overlapped orfs
        std::vector<bool> overlap(n, false);
        for (int i = 0; i < n; i ++) {
            for (int j = i + 1; j < n; j ++) {
                int olen = orfs[inits[i]] - orfs[inits[j]];
                if (olen > MIN_LEN) overlap[i] = overlap[j] = true;
            }
        }
        for (int i = n-1; i > -1; i --) {
            if (overlap[i]) inits.erase(inits.begin()+i);
        }
        // init first cluster center
        seeds = inits;
        n = (int) seeds.size();
        for (int i = 0; i < n; i ++) add(i, 0);
        for (int j = 0; j < DIM_K; j ++) c[j] /= n;
        jcache[0] = dot_ff(c.data(), c.data(), DIM_K);

        Debug() << format_log("- Initial Seed ORFs:", Str(n));
    }

    /* sample long orfs as seed ORFs */
    void seed(int LONG = 300) {
        seeds = orfs.filter(LONG, 1.0F);
        n = (int) seeds.size();

        labels.resize(n), icache.resize(n);
        std::fill(labels.begin(), labels.end(), -1);
        for (int i = 0; i < seeds.size(); i ++)
            icache[i] = dot_ff(irow(i), irow(i), DIM_K);

        Debug() << format_log("- Seed ORFs:", Str(n));
    }

    /* init 2~k cluster centers using kmeans ++ */
    Kmeans &kmeans_pp(int seed) {
        std::default_random_engine engine(seed);
        flt_arr min_d(n, INFINITY);
        std::fill(c.begin() + DIM_K, c.end(), 0.0F);

        for (int j = 1; j < k; j ++) {
            float sum = 0.0F;

            for (int i = 0; i < n; i ++) {
                float dk = dist_sq(i, j-1);
                min_d[i] = std::min(dk, min_d[i]);
                sum += min_d[i];
            }

            std::uniform_real_distribution<float> udist(0.0F, sum);
            float pin = udist(engine);
            
            sum = 0.0F;
            for (int i = 0; i < n; i ++) {
                if (pin < (sum += min_d[i])) {
                    std::memcpy(jrow(j), irow(i), DIM_K*sizeof(float));
                    jcache[j] = dot_ff(jrow(j), jrow(j), DIM_K);
                    break;
                }
            }
        }

        return *this;
    }

    /* predict k-means labels */
    float predict() {
        int unc = 0;

        for (int i = 0; i < n; i ++) {
            float min_d = INFINITY;
            int label = -1;

            for (int j = 0; j < k; j ++) {
                float d = dist_sq(i, j);
                if (d < min_d) min_d = d, label = j;
            }

            unc += (int) (label == labels[i]);
            labels[i] = label;
        }

        return (float)unc/n;
    }

    /* calculate sum of squared errors */
    float SSE() {
        float sum = 0.0F;
        for (int i = 0; i < n; i ++) {
            sum += dist_sq(i, labels[i]);
        }
        return sum;
    }

    void update() {
        std::fill(c.begin(), c.end(), 0.0F);
        int_arr cnt(k, 0);

        for (int i = 0; i < n; i ++) {
            cnt[labels[i]] ++;
            add(i, labels[i]);
        }
        
        for (int j = 0; j < k; j ++) {
            float *j_row = jrow(j);
            for (int m = 0; m < DIM_K; m ++) 
                j_row[m] /= cnt[j];
            jcache[j] = dot_ff(jrow(j), jrow(j), DIM_K);
        }
    }

    int positive() {
        int_arr cnt(k, 0);
        int_set iset(inits.begin(), inits.end());

        for (int i = 0; i < n; i ++) {
            if (iset.find(seeds[i]) != iset.end())
                cnt[labels[i]] ++;
        }

        return argmax(cnt);
    }

    int train(float gc_cont, int LONG) {
        init(gc_cont); seed(LONG);
        float best_SSE = INFINITY;
        flt_arr best_c;

        for (int seed = 0; seed < 10; seed ++) {
            kmeans_pp(seed);
            while(true) {
                if (predict() > 0.99F) break;
                update();
            }

            float sse = SSE();
            if (sse > best_SSE) {
                best_SSE = sse;
                best_c = c;
            }
        }

        c = best_c;
        predict();

        return positive();
    }
};

class Labeler {
   /* ab initio label task */
  private:
    bioinfo::genome    genome;  // genome read from input
    _time_t           start_t;  // program start time
    cxxopts::ParseResult args;  // program arguments
    bioinfo::orfs        orfs;  // open reading frames
    str                 table;  // translation table
    bioinfo::orfs       genes;  // putative genes
    float             gc_cont;  // GC content of the genome
    /* ATG GTG TTG (AUG GUG UUG) */
    str_arr STARTS { ATG, GTG, TTG };
    /* TAA TAG TGA (UAA UAG UGA) */
    str_arr STOPS  { TAA, TAG, TGA };
    
    static cxxopts::ParseResult parse_args(int argc, char* argv[]) {
        // executable file name
        str exe_name = utils::filename(argv[0]);
    
        /* set and parse parameters */
        cxxopts::Options options(exe_name);

        /* general parameters */
        options.add_options("General")
            ("h,help",     "Print help menu and exit.")

            ("q,quiet",    "Run quietly with no stderr output.")

            ("v,version",  "Print version info and exit.")

            ("T,threads",  "Number of threads to use. (default: half)",
             cxxopts::value<uint32_t>());
        
        /* input/output parameters */
        options.add_options("Input/Output")
            ("i,input",    "Specify FASTA/Genbank/EMBL input file or their compressed versions (gzip). (default: stdin)",
            cxxopts::value<str>())
            
            ("o,output",   "Write results to the output file or '-' as stdout. (default: none)",
            cxxopts::value<str>())
            
            ("a,faa",      "Write protein translations to the selected file or '-' as stdout.",
            cxxopts::value<str>())
    
            ("d,fna",      "Write nucleotide sequences of genes to the selected file or '-' as stdout.",
            cxxopts::value<str>());
        
        /* orf-finder parameters */
        options.add_options("ORF-finder")
            ("g,table",    "Specify a translation table to use (1, 4, 11, 15, 16, 25).",
            cxxopts::value<str>()->default_value("11"))

            ("l,minlen",   "Specify the mininum length of putative genes.",
            cxxopts::value<uint32_t>()->default_value("90"))

            ("c,circ",     "Treat default topology as circular.");

        /* labeler parameters */
        options.add_options("Labeler")
            ("L,long",     "Specify the mininum length of long orfs.",
            cxxopts::value<uint32_t>()->default_value("300"))

            ("n,clusters", "Specify number of clusters for K-means.",
            cxxopts::value<uint32_t>()->default_value("4"))

            ("z,zcurve",   "Output the 189-digit Z-curve parameters.",
            cxxopts::value<str>());
        
        cxxopts::ParseResult args;
        
        try { 
            args = options.parse(argc, argv); 
        } catch (const std::exception& e) {
            std::cerr << "Error: " << e.what() 
                      << "\nUse '-h' or '--help' to show help info\n";
            ARG_ERR;
        }

        if (argc <= 1 || args.count("help")) {
            Debug() << options.help() << "\nExample: " << options.program() << ' '
                    << "-i example.fa -o example.gff -c -f gff                 \n";
            if (argc <= 1) {
                Debug() << "\nPress Enter or Ctrl+C to exit the help ... ";
                std::cin.get();
            }
            NON_ERR;
        } else if (args.count("version")) { Debug() << VERSION << "\n"; NON_ERR; 
        } else if (args.count("quiet"))   { debug = &DEVNULL; };

        return args;
    }

    void locate_orfs(str_arr &starts, str_arr &stops, int minlen) {
        orfs.clear();
        for (auto &scaffold : genome) {
            if (scaffold.len >= minlen) {
                orfs.locate_orfs(scaffold, starts, stops, minlen);
                if (orfs.size() > 1000000) {
                    std::cerr << "\nError: too many candidate ORFs ( > 1000000)\n";
                    SUB_ERR;
                }
            } else Debug() << "Warning: skip " << scaffold.name << " (too short)\n";
        }
        
        if (orfs.size() <= 0) {
            std::cerr << "\nError: no candidate ORF was found.\n";
            SUB_ERR;
        } else { orfs.init(); }
    }

  public:
    Labeler(int argc, char* argv[]) {
        start_t = SYSTEM_CLOCK now();
        args = parse_args(argc, argv);
        Debug() << "\nAB INITIO LABELLING ORFS USING THE Z-CURVE METHOD\n\n";
        auto start_time = SYSTEM_CLOCK to_time_t(start_t);
        std::ostringstream oss;
        oss << std::put_time(std::localtime(&start_time), "%Y-%m-%d %H:%M:%S")
            << format_log("Program Start Time:", oss.str());

        str extensions;
        #ifdef _OPENMP
        #ifdef ZLIB
        extensions += " zlib";
        #endif
        extensions += " OpenMP";
        auto n_threads = omp_get_max_threads() / 2;
        if(args.count("threads")) n_threads = args["threads"].as<uint32_t>();
        omp_set_num_threads(n_threads);
        Debug() << format_log("Threads Used: ", Str(n_threads));
        #endif
        #ifdef __AVX__
        extensions += " AVX/FMA";
        #endif
        if (extensions.empty()) extensions = "None";

        Debug() << format_log("Extensions Enabled:", extensions);
    }

    Labeler &Load_Genome() {
        bool circ = (bool) args.count("circ");
        Debug() << format_log("Default Topology:", circ ? "Circular" : "Linear");
        str filename = args["input"].as<str>();
        if (filename.empty()) filename = "-";
        utils::read_source(filename, genome, circ);
        size_t size = 0; int gc_count = 0;
        for (auto &scaffold : genome) {
            size += scaffold.len;
            gc_count += scaffold.gc_count;
        }
        gc_cont = (float) gc_count / size;
        Debug() << format_log("Number of Scaffolds:", Str(genome.size()))
                << format_log("Genome Size:", Str(size) + " bp")
                << "GC Content:" << std::setw(36) << std::fixed 
                << std::setprecision(2) << gc_cont*100 << " %\n";
        return *this;
    }

    Labeler &Check_Table() {
        table = args["table"].as<str>();
        if (table == "1") {
            STARTS.assign({ ATG });
        } else if (table == "4") {
            STOPS. assign({ TAA, TAG });
            trans_tbl[14] = 'W';
        } else if (table == "15") {
            STOPS. assign({ TAA, TGA });
            trans_tbl[11] = 'Q';
        } else if (table == "16") {
            STOPS. assign({ TAA, TGA });
            trans_tbl[11] = 'L';
        } else if (table == "25") {
            STOPS. assign({ TAA, TAG });
            trans_tbl[14] = 'G';
        } else if (table != "11") {
            std::cerr << "\nError: unsupported translation table " << table << ".\n\n"
                         "Genetic code options: \n\n 1  - Standard\n"
                         " 4  - Mycoplasma, Spiroplasma\n"
                         " 11 - Bacteria, Archaea\n"
                         " 15 - Blepharisma Nuclear\n"
                         " 16 - Chlorophycean Mitochondria\n"
                         " 25 - Candidate Division SR1, Gracilibacteria\n\n"
                         "Please see https://www.ddbj.nig.ac.jp/ddbj/geneticcode-e.html.\n";
            ARG_ERR;
        }
        Debug() << format_log("Translation Table:", table);
        return *this;
    }

    Labeler &Locate_Orfs() {
        int minlen = (int)args["minlen"].as<uint32_t>();
        if (minlen < MIN_LEN) {
            std::cerr << "Error: mininum gene length should not be less than 75 nt\n";
            ARG_ERR;
        }
        Debug() << format_log("Mininum Gene Length:", Str(minlen) + " nt");
        locate_orfs(STARTS, STOPS, minlen); 
        for (int i = 0; i < orfs.size(); i ++) orfs[i].froze();
        Debug() << format_log("Number of Candidate ORFs:", Str(orfs.size()));
        return *this;
    }

    Labeler &Cluster_Orf() {
        static int seed = 42;
        
        int K = (int)args["clusters"].as<uint32_t>();
        int LONG = (int)args["long"].as<uint32_t>();

        Kmeans kms = Kmeans(orfs, K);

        int pos_label = kms.train(gc_cont, LONG);
        
        int n_genes = 0;
        for (int i = 0; i < kms.n; i ++) {
            if (kms.labels[i] == pos_label) {
                orfs[kms.seeds[i]].i_score = +1.0F;
                n_genes ++;
            } else {
                orfs[kms.seeds[i]].i_score = -1.0F;
            }
            genes.push_back(orfs[kms.seeds[i]]);
        }

        Debug() << format_log("Number of Putative Genes:", Str(n_genes));

        return *this;
    }

    Labeler &Save_Result() {
        if (args.count("output")) {
            str outfile = args["output"].as<str>();
            utils::write_output(outfile, genes, start_t, table, "gff");
        }

        if (args.count("faa")) {
            auto faaout = args["faa"].as<str>();
            utils::write_output(faaout, genes, start_t, table, "faa");
        }

        if (args.count("fna")) {
            auto fnaout = args["fna"].as<str>();
            utils::write_output(fnaout, genes, start_t, table, "fna");
        }

        return *this;
    }

    ~Labeler() {
        using namespace std::chrono;
        auto end_t = SYSTEM_CLOCK now();
        auto duration = duration_cast<milliseconds>(end_t - start_t);
        auto seconds = duration.count() / 1000.0;
        Debug() << "\nFinished in " << std::fixed << std::setprecision(3) << seconds << " s\n";
    }
};

int main(int argc, char *argv[]) {
    std::ios::sync_with_stdio(false);
    
    Labeler(argc, argv)
        . Load_Genome()
        . Check_Table()
        . Locate_Orfs()
        . Cluster_Orf()
        . Save_Result();
}
