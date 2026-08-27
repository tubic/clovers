#include "Common.hpp"
#include <random>

/* feature window length used by k-means (last DIM_K dims of the Z-curve vector) */
#define DIM_K     144
/* offset of the feature window within the full DIM_S-length feature vector */
#define FEAT_OFF  (DIM_S - DIM_K)
/* min ORF length (nt) for trusted initial seeds */
#define S_LEN     550

/* ignore */
const float META[] = {0};

/* codon -> amino acid table; index = 3-digit base-4 codon (T=0, C=1, A=2, G=3) */
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

/* k-means clustering of ORFs over standardized Z-curve features */
class Kmeans {
    bioinfo::orfs &orfs;    // all orfs from the genome
    flt_arr c;              // k-means centroids
    flt_arr icache, jcache; // square norm cache

    /* sample pointer (feature window starts at FEAT_OFF) */
    float *irow(int idx) {
        static int off = FEAT_OFF;
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
    /* calculate distance between two samples */
    float sample_d2(int a, int b) {
        return icache[a] + icache[b] - 2 * dot_ff(irow(a), irow(b), DIM_K);
    }
    /* add weighted sample to centroid (vector) */
    void add(int i, int j) {
        float w = weights[i];
        float *i_row = irow(i), *j_row = jrow(j);
        for (int m = 0; m < DIM_K; m ++) {
            j_row[m] += w * i_row[m];
        }
    }

  public:
    /* n_clusters */
    int k, n;
    /* max iterations for k-means convergence */
    static constexpr int MAX_ITER = 100;
    /* centroid drift threshold for convergence */
    static constexpr float DRIFT_EPS = 1E-4F;
    /* indices of samples and their labels */
    int_arr inits, seeds, labels;
    /* per-sample weight (uniform for now; kept so init and update agree) */
    flt_arr weights;

    /* encode Z-curve features, then z-score standardize the feature window in place */
    Kmeans(bioinfo::orfs &orfs, int k): orfs(orfs), k(k) {
        c = flt_arr(DIM_K * k, 0.0F);
        jcache.resize(k);
        zcurve::encode(orfs);
        /* z-score standardize the feature window in place */
        const int n_all = (int) orfs.size();
        for (int m = 0; m < DIM_K; m ++) {
            float mean = 0.0F;
            for (int i = 0; i < n_all; i ++) mean += orfs.row(i)[FEAT_OFF+m];
            mean /= n_all;
            float var = 0.0F;
            for (int i = 0; i < n_all; i ++) {
                float d = orfs.row(i)[FEAT_OFF+m] - mean;
                var += d * d;
            }
            var /= n_all;
            float sd = sqrtf(var);
            if (sd < 1E-9F) sd = 1.0F;   /* constant dimension */
            for (int i = 0; i < n_all; i ++)
                orfs.row(i)[FEAT_OFF+m] = (orfs.row(i)[FEAT_OFF+m] - mean) / sd;
        }
    }

    /* pick out initial seed orfs */
    void init(float gc_cont) {
        // get long orfs with expected G2 bias
        float G2_bias = 0.2F * gc_cont + 0.12F;
        inits = orfs.filter(S_LEN, gc_cont > 0.57F ? G2_bias : 1.0F);
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
        // init first cluster center (weighted mean of inits)
        seeds = inits;
        n = (int) seeds.size();
        weights.assign(n, 1.0F);
        if (n > 0) {
            /* compute the geometric center (weighted mean) of the init orfs */
            for (int i = 0; i < n; i ++) add(i, 0);
            for (int j = 0; j < DIM_K; j ++) c[j] /= n;
            /* use the actual sample closest to the center as the first centroid */
            float best_d = INFINITY;
            int best_i = 0;
            for (int i = 0; i < n; i ++) {
                float *s = irow(i);
                float d = 0.0F;
                for (int m = 0; m < DIM_K; m ++) {
                    float t = s[m] - c[m];
                    d += t * t;
                }
                if (d < best_d) { best_d = d; best_i = i; }
            }
            std::memcpy(c.data(), irow(best_i), DIM_K*sizeof(float));
        }
        jcache[0] = dot_ff(c.data(), c.data(), DIM_K);

        Debug() << format_log("- Indicator ORFs:", Str(n));
    }

    /* sample long orfs as seed ORFs (hard length threshold) */
    void seed(int LONG = 240) {
        seeds.clear(); weights.clear();
        for (int i = 0; i < (int)orfs.size(); i ++) {
            int len = orfs[i].len;
            if (len >= LONG) {
                seeds.push_back(i);
                weights.push_back(1.0F);
            }
        }
        n = (int) seeds.size();

        if (n <= 0) {
            std::cerr << "\nError: no long ORFs found (length >= " << LONG << ")\n";
            SUB_ERR;
        }

        labels.resize(n), icache.resize(n);
        std::fill(labels.begin(), labels.end(), -1);
        for (int i = 0; i < seeds.size(); i ++)
            icache[i] = dot_ff(irow(i), irow(i), DIM_K);

        Debug() << format_log("- Long ORFs:", Str(n));
    }

    /* seed the remaining centers (2~k) by Arthur/Vassilvitskii k-means++ */
    Kmeans &kmeans_pp(int seed) {
        std::default_random_engine engine(seed);
        std::fill(c.begin() + DIM_K, c.end(), 0.0F);

        int n_local_trials = 2 + (int)std::log((double)k);

        flt_arr closest_d2(n);
        float current_pot = 0.0F;
        for (int i = 0; i < n; i ++) {
            closest_d2[i] = dist_sq(i, 0);
            current_pot += closest_d2[i] * weights[i];
        }

        for (int j = 1; j < k; j ++) {
            int best_pick = -1;
            float best_pot = INFINITY;
            flt_arr cand_d2(n), best_d2(n);

            for (int t = 0; t < n_local_trials; t ++) {
                int pick = -1;
                if (current_pot > 0.0F) {
                    std::uniform_real_distribution<float> udist(0.0F, current_pot);
                    float pin = udist(engine);
                    float acc = 0.0F;
                    for (int i = 0; i < n; i ++) {
                        if ((acc += closest_d2[i] * weights[i]) >= pin) { pick = i; break; }
                    }
                }
                if (pick < 0) {
                    std::uniform_int_distribution<int> udist(0, n-1);
                    pick = udist(engine);
                }

                /* candidate potential: closest distances become
                   min(D(x)^2, d2 to the candidate), summed with weights */
                float pot = 0.0F;
                for (int i = 0; i < n; i ++) {
                    float d = std::min(closest_d2[i], sample_d2(i, pick));
                    cand_d2[i] = d;
                    pot += d * weights[i];
                }

                /* greedily keep the candidate that reduces the potential most */
                if (pot < best_pot) {
                    best_pot = pot;
                    best_pick = pick;
                    best_d2.swap(cand_d2);
                }
            }

            /* permanently add the winning candidate as center j */
            std::memcpy(jrow(j), irow(best_pick), DIM_K*sizeof(float));
            jcache[j] = dot_ff(jrow(j), jrow(j), DIM_K);
            closest_d2.swap(best_d2);
            current_pot = best_pot;
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

            unc += (int) (label == labels[i]);   /* count samples whose label stayed the same */
            labels[i] = label;
        }

        return (float)unc/n;   /* fraction of unchanged labels (convergence ratio) */
    }

    /* calculate sum of squared errors */
    float SSE() {
        float sum = 0.0F;
        for (int i = 0; i < n; i ++) {
            sum += dist_sq(i, labels[i]);
        }
        return sum;
    }

    /* recompute each centroid as the weighted mean of its assigned samples */
    void update() {
        flt_arr old_c = c;
        std::fill(c.begin(), c.end(), 0.0F);
        flt_arr wsum(k, 0.0F);

        for (int i = 0; i < n; i ++) {
            wsum[labels[i]] += weights[i];
            add(i, labels[i]);
        }
        
        for (int j = 0; j < k; j ++) {
            float *j_row = jrow(j);
            if (wsum[j] > 0.0F) {
                for (int m = 0; m < DIM_K; m ++) 
                    j_row[m] /= wsum[j];
            } else {
                /* keep the previous centroid to avoid empty-cluster NaN */
                std::memcpy(j_row, old_c.data() + j*DIM_K, DIM_K*sizeof(float));
            }
            jcache[j] = dot_ff(jrow(j), jrow(j), DIM_K);
        }
    }

    /* max squared centroid displacement between two snapshots */
    float drift_sq(const flt_arr &old_c) {
        float max_d = 0.0F;
        for (int j = 0; j < k; j ++) {
            const float *a = old_c.data() + j*DIM_K;
            const float *b = jrow(j);
            float d2 = 0.0F;
            for (int m = 0; m < DIM_K; m ++) {
                float e = a[m] - b[m];
                d2 += e * e;
            }
            if (d2 > max_d) max_d = d2;
        }
        return max_d;
    }

    /* pick the cluster most enriched in trusted init ORFs (the putative gene cluster);
       near-ties broken by closeness to centroid 0 */
    int positive() {
        int_arr cnt(k, 0);
        flt_arr size(k, 0.0F);
        int_set iset(inits.begin(), inits.end());

        for (int i = 0; i < n; i ++) {
            size[labels[i]] += weights[i];
            if (iset.find(seeds[i]) != iset.end())
                cnt[labels[i]] ++;
        }

        /* soft decision: score by inits density (cnt^2/size) */
        float best_score = -1.0F, best_d = INFINITY;
        int best = 0;
        for (int j = 0; j < k; j ++) {
            if (size[j] <= 0.0F) continue;
            float score = (float)cnt[j] * cnt[j] / size[j];
            float d = sqrtf(std::max(jcache[j] + jcache[0] - 2 * dot_ff(jrow(j), jrow(0), DIM_K), 0.0F));
            if (score > best_score * 1.05F || (score > best_score * 0.95F && d < best_d)) {
                best_score = score;
                best_d = d;
                best = j;
            }
        }

        return best;
    }

    /* inits-density score of the best cluster (gene enrichment) */
    float pos_density() {
        int_arr cnt(k, 0);
        flt_arr size(k, 0.0F);
        int_set iset(inits.begin(), inits.end());
        for (int i = 0; i < n; i ++) {
            size[labels[i]] += weights[i];
            if (iset.find(seeds[i]) != iset.end()) cnt[labels[i]] ++;
        }
        float best = -1.0F;
        for (int j = 0; j < k; j ++)
            if (size[j] > 0.0F) best = std::max(best, (float)cnt[j] * cnt[j] / size[j]);
        return best;
    }

    /* run k-means with 50 random restarts, keep the best run, return its positive cluster */
    int train(float gc_cont, int LONG) {
        init(gc_cont); seed(LONG);
        float best_SSE = INFINITY, best_den = -1.0F;
        flt_arr best_c;
        bool have_best = false;

        for (int seed = 0; seed < 50; seed ++) {
            kmeans_pp(seed);
            for (int iter = 0; iter < MAX_ITER; iter ++) {
                predict();
                flt_arr prev_c = c;
                update();
                if (drift_sq(prev_c) < DRIFT_EPS * DRIFT_EPS) break;
            }

            float sse = SSE();
            /* reject degenerate runs that collapsed to a single cluster */
            int nonempty = 0;
            int_arr sz(k, 0);
            for (int i = 0; i < n; i ++) sz[labels[i]] ++;
            for (int j = 0; j < k; j ++) nonempty += (sz[j] > 0);

            /* pick the restart whose best cluster has the highest inits density */
            float den = pos_density();
            bool take = nonempty >= 2 && (!have_best || den > best_den * 1.01F);
            if (take) best_den = den;
            if (take) {
                have_best = true;
                best_SSE = sse;
                best_c = c;
            }
        }

        c = best_c;
        /* sync jcache with the restored centroids (predict() depends on it) */
        for (int j = 0; j < k; j ++)
            jcache[j] = dot_ff(jrow(j), jrow(j), DIM_K);
        if (!have_best) {
            /* all runs collapsed (shouldn't happen): keep the last centroids */
            Debug() << "Warning: all k-means runs collapsed; using last run\n";
            kmeans_pp(0);
            for (int iter = 0; iter < MAX_ITER; iter ++) {
                predict();
                flt_arr prev_c = c;
                update();
                if (drift_sq(prev_c) < DRIFT_EPS * DRIFT_EPS) break;
            }
        }
        predict();

        {
            int_arr sz(k, 0);
            for (int i = 0; i < n; i ++) sz[labels[i]] ++;
            std::ostringstream oss;
            for (int j = 0; j < k; j ++) oss << " C" << j << "=" << sz[j];
            Debug() << format_log("DBG Clusters:", oss.str());
        }

        return positive();
    }
};

/* ab initio gene-calling pipeline: genome -> ORFs -> k-means -> RBS -> SVM -> output */
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
    /* ---- automatic (L, n) selection -----------------------------------
         - L: 180 for GC < 45% or GC >= 63%, otherwise 240;
         - n grows with GC content (re-optimized boundaries):
              GC < 35%  -> 2
              GC < 40%  -> 3
              GC < 57%  -> 4
              GC < 62.6%-> 5
              else      -> 6
       Used only when the user does not pass -L / -n explicitly. */
    static int auto_long(float gc_cont) {
        float gc = gc_cont * 100.0F;   /* fraction -> percent */
        return (gc < 45.0F || gc >= 63.0F) ? 180 : 240;
    }
    static int auto_clusters(float gc_cont, str tbl) {
        float gc = gc_cont * 100.0F;   /* fraction -> percent */
        if (gc < 35.0F)   return (tbl != "11" ? 3 : 2);
        if (gc < 40.0F)   return 3;
        if (gc < 57.0F)   return 4;
        if (gc < 62.6F)   return 5;
        return 6;
    }
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
            ("L,long",     "Specify the mininum length of long orfs (default: auto by GC).",
            cxxopts::value<uint32_t>())

            ("n,clusters", "Specify number of clusters for K-means (default: auto by GC%).",
            cxxopts::value<uint32_t>());
        
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
                    << "-i example.fa -o example.gff -c                    \n";
            NON_ERR;
        } else if (args.count("version")) { Debug() << VERSION << "\n"; NON_ERR; 
        } else if (args.count("quiet"))   { debug = &DEVNULL; };

        return args;
    }

    /* scan every scaffold for ORFs bounded by the given start/stop codons */
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
    /* parse args, print banner, configure threads */
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
        auto n_threads = std::max(1, omp_get_max_threads() / 2);
        if(args.count("threads")) {
            n_threads = (int) args["threads"].as<uint32_t>();
            if (n_threads < 1) {
                std::cerr << "\nError: number of threads should be at least 1\n";
                ARG_ERR;
            }
        }
        omp_set_num_threads(n_threads);
        Debug() << format_log("Threads Used: ", Str(n_threads));
        #endif
        #ifdef __AVX__
        extensions += " AVX/FMA";
        #endif
        if (extensions.empty()) extensions = "None";

        Debug() << format_log("Extensions Enabled:", extensions);
    }

    /* read the genome (optionally circular) and compute its GC content */
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

    /* adjust start/stop codons and the translation table for the chosen genetic code */
    Labeler &Check_Table() {
        table = args["table"].as<str>();
        if (table == "1") {
            STARTS.assign({ ATG });
        } else if (table == "4") {
            STOPS. assign({ TAA, TAG });
            trans_tbl[14] = 'W';   /* TGA reassigned as Trp */
        } else if (table == "15") {
            STOPS. assign({ TAA, TGA });
            trans_tbl[11] = 'Q';   /* TAG reassigned as Gln */
        } else if (table == "16") {
            STOPS. assign({ TAA, TGA });
            trans_tbl[11] = 'L';   /* TAG reassigned as Leu */
        } else if (table == "25") {
            STOPS. assign({ TAA, TAG });
            trans_tbl[14] = 'G';   /* TGA reassigned as Gly */
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

    /* find candidate ORFs and freeze their boundaries */
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

    /* cluster ORFs by Z-curve features and give seeds a prior score by cluster membership */
    Labeler &Cluster_Orf() {
        /* use the learned automatic rule unless the user overrides */
        int K = args.count("clusters")
            ? (int)args["clusters"].as<uint32_t>()
            : auto_clusters(gc_cont, table);
        int LONG = args.count("long")
            ? (int)args["long"].as<uint32_t>()
            : auto_long(gc_cont);
        if (K < 2 || K > 64) {
            std::cerr << "\nError: number of clusters should be between 2 and 64\n";
            ARG_ERR;
        }
        if (!args.count("clusters"))
            Debug() << format_log("Auto Clusters (By GC%):", Str(K));
        if (!args.count("long"))
            Debug() << format_log("Auto Long Threshold:", Str(LONG) + " nt");

        Kmeans kms = Kmeans(orfs, K);

        int pos_label = kms.train(gc_cont, LONG);
        
        int n_pos = 0, n_neg = 0;
        for (int i = 0; i < kms.n; i ++) {
            if (kms.labels[i] == pos_label) {
                orfs[kms.seeds[i]].i_score = 0.75;   /* in the positive (gene) cluster */
                n_pos ++;
            } else {
                orfs[kms.seeds[i]].i_score = 0.25;   /* in a negative cluster */
                n_neg ++;
            }
        }

        Debug() << format_log("Number of Seed ORFs:", Str(n_pos) + ":" + Str(n_neg));

        return *this;
    }

    /* iteratively refine the RBS (Shine-Dalgarno) model over Markov orders 0..2 */
    Labeler &Revise_Orfs() {
        Debug() << format_log("Revising RBS Model:", "Round #0", true);
        model::pmm model(orfs, table);

        int_arr seeds = orfs.filter(300, 0.5, true);
        int round = 0;
        if (seeds.size() >= 0) for (int order = 0; order <= 2; order ++) {   /* train each Markov order, up to 20 rounds */
            model.reset(order);
            for (int i = 0; i < 20; i ++) {
                Debug() << format_log("Revising RBS Model:", str("Round #")+Str(++round), true);
                model.train(seeds, STARTS);
                if (model.revise(seeds) > 0.99) break;
            }
        }
        
        model.revise_all(); 
        Debug() << '\n';
        return *this;
    }

    /* train a genome-specific RBF-SVM on the clustered seeds and call genes by its score */
    Labeler &Csvm_Refine() {
        zcurve::encode(orfs.concat());
        auto model = model::svm(orfs);
        Debug() << "Constructing Genome-Specific Gene Model (RBF-SVM)\n";
        model.sample(0.7, 0.3).set_C(0.15F).train().predict();
        int cnt = 0;
        for (int i = 0; i < orfs.size(); i ++) {
            orfs[i].i_score = orfs[i].c_score - 0.5F;
            cnt += (orfs[i].i_score > 0);
        }
        orfs.yield("ignore", 0.0F, genes);
        Debug() << format_log("Putative Genes:", Str(cnt));
        return *this;
    }

    /* write GFF and the optional protein/nucleotide outputs */
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

    /* print total runtime */
    ~Labeler() {
        using namespace std::chrono;
        auto end_t = SYSTEM_CLOCK now();
        auto duration = duration_cast<milliseconds>(end_t - start_t);
        auto seconds = duration.count() / 1000.0;
        Debug() << "\nFinished in " << std::fixed << std::setprecision(3) << seconds << " s\n";
    }
};

/* pipeline entry point: chain the Labeler stages, destruct at end of full expression */
int main(int argc, char *argv[]) {
    std::ios::sync_with_stdio(false);
    
    try {
        Labeler(argc, argv)
            . Load_Genome()
            . Check_Table()
            . Locate_Orfs()
            . Cluster_Orf()
            . Revise_Orfs()
            . Csvm_Refine()
            . Save_Result();
    } catch (const std::exception &e) {
        std::cerr << "\nError: " << e.what() << '\n';
        SUB_ERR;
    }
}
