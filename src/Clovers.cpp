#include "Common.hpp"
#include "Meta.cpp"

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

class Clovers {
  /* clovers prediction task */
  private:
    bioinfo::genome    genome;  // genome read from input
    _time_t           start_t;  // program start time
    cxxopts::ParseResult args;  // program arguments
    bioinfo::orfs        orfs;  // open reading frames
    bioinfo::orfs       genes;  // putative genes
    str                 table;  // translation table
    str                  proc;  // procedure
    bool        revised=false;  // revise gene start site
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
            
            ("f,format",   "Select output format (gff, med, gbk, gbk-full).",
            cxxopts::value<str>()->default_value("gff"))
            
            ("a,faa",      "Write protein translations to the selected file or '-' as stdout.",
            cxxopts::value<str>())
    
            ("d,fna",      "Write nucleotide sequences of genes to the selected file or '-' as stdout.",
            cxxopts::value<str>());
        
        /* clovers parameters */
        options.add_options("CLOVERS")
            ("g,table",    "Specify a translation table to use (1, 4, 11, 15, 16, 25, auto).",
            cxxopts::value<str>()->default_value("11"))

            ("l,minlen",   "Specify the mininum length of putative genes.",
            cxxopts::value<uint32_t>()->default_value("90"))

            ("c,circ",     "Treat default topology as circular.")
        
            ("p,proc",     "Select procedure (single or meta).",
            cxxopts::value<str>()->default_value("single"))

            ("s,thres",    "Specify putative gene score threshold.",
            cxxopts::value<float>()->default_value("0.499"))
            
            ("t,train",    "Write (if none exists) or use the specified SVM training file.",
            cxxopts::value<str>());
        
        /* tri-tisa+ parameters */
        options.add_options("TriTISA+")
            ("n,bypass",   "Bypass TriTISA+ and output longest ORFs.")

            ("M,maxiter",  "Max iteration times for gene start revision.",
            cxxopts::value<uint32_t>()->default_value("20"))
            
            ("R,rbs",      "Write (if none exists) or use the RBS training file.",
            cxxopts::value<str>());

        /* gol-reporter parameters */
        options.add_options("GOL-Reporter")
            ("L,min-olen", "Specify the mininum overlapped length between two genes.",
            cxxopts::value<uint32_t>()->default_value("120"))

            ("O,overlap",  "Write overlapped genes to the selected file (GFF3 format).",
            cxxopts::value<str>())

            ("A,amino",    "Write protein translations of overlapped genes to the selected file.",
            cxxopts::value<str>())

            ("D,nucl",     "Write nucleotide sequences of overlapped genes to the selected file.",
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
            Debug() << "- - - - - - - - - - - - - - - - - - - - - - - - - - - -\n"
                       "PROTEIN-CODING GENE RECOGNITION SYSTEM OF CLOVERS 1.1.0\n\n"
                       "Copyright:  (C) 2003-2026 TUBIC, Tianjin University    \n"
                       "Authors:    Zetong Zhang, Yan Lin*, Feng Gao*          \n"
                       "Date:       March 31, 2026                             \n"
                       "Contact:    ylin@tju.edu.cn | fgao@tju.edu.cn          \n"
                       "- - - - - - - - - - - - - - - - - - - - - - - - - - - -\n"
                    << options.help() << "\nExample: " << options.program() << ' '
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
        for (auto &seq : genome) {
            if (seq.len >= minlen) {
                orfs.locate_orfs(seq, starts, stops, minlen);
                if (orfs.size() > 1000000) {
                    std::cerr << "\nError: too many candidate ORFs ( > 1000000)\n";
                    SUB_ERR;
                }
            } else Debug() << "Warning: skip " << seq.name << " (too short)\n";
        }
        
        if (orfs.size() <= 0) {
            std::cerr << "\nError: no candidate ORF was found.\n";
            SUB_ERR;
        } else { orfs.init(); }
    }

    str check_code() {
        static str_arr starts { ATG };
        static std::map<str, str_arr> stops = {
            {"11", { TAA, TAG, TGA }},
            {"15", { TAA, TGA }     },
            {"4",  { TAA, TAG }     },
        };
        char buff[128];
        str best_table; float max_score = -1;
        for (auto &it : stops) {
            int poslen = 0, sumlen = 0;
            locate_orfs(starts, it.second, MIN_LEN);
            Init_Scores(false);
            for (int i = 0; i < orfs.size(); i ++) {
                if (orfs[i].i_score > 0.5) poslen += orfs[i].len;
                sumlen += orfs[i].len;
            }
            float score = (float) poslen / sumlen;
            sprintf(buff, "%.3f", score);
            Debug() << format_log(str("Scoring Translation Table #")+it.first+":", buff);
            if (score > max_score) {
                best_table = it.first;
                max_score = score;
            }
        }
        return best_table;
    }

  public:
    Clovers(int argc, char* argv[]) {
        start_t = SYSTEM_CLOCK now();
        auto start_time = SYSTEM_CLOCK to_time_t(start_t);
        std::ostringstream oss;
        oss << std::put_time(std::localtime(&start_time), "%Y-%m-%d %H:%M:%S");

        args = parse_args(argc, argv);

        Debug() << "\nPROTEIN-CODING GENE RECOGNITION SYSTEM OF CLOVERS\n\n"
                << format_log("Program Start Time:", oss.str())
                << format_log("Loaded Model:", "Prokaryote/Phage");

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

        proc = std::tolower(args["proc"].as<str>()[0])=='m' ? "Metagenome" : "Single";

        Debug() << format_log("Extensions Enabled:", extensions)
                << format_log("Procedure:", proc);
    }

    Clovers &Load_Genome() {
        bool circ = (bool) args.count("circ");
        Debug() << format_log("Default Topology:", circ ? "Circular" : "Linear");
        str filename = args["input"].as<str>();
        if (filename.empty()) filename = "-";
        utils::read_source(filename, genome, circ);
        size_t size = 0; int gc_count = 0;
        for (auto &seq : genome) {
            size += seq.len;
            gc_count += seq.gc_count;
        }
        float gc_cont = (float) gc_count / size * 100;
        Debug() << format_log("Number of Scaffolds:", Str(genome.size()))
                << format_log("Genome Size:", Str(size) + " bp")
                << "GC Content:" << std::setw(36) << std::fixed << std::setprecision(2) << gc_cont << " %\n";
        return *this;
    }

    Clovers &Init_Scores(bool debug=true) {
        zcurve::encode(orfs);
        int gc_intv[N_MODELS] = {0};
        {
            float max_gc = 0.20;
            for (int i = 0, j = 0; i < orfs.size(); i ++) {
                if (orfs[i].gc_cont > max_gc) {
                    max_gc += 0.01;
                    if ((++j) >= (N_MODELS)) break;
                }
                gc_intv[j] ++;
            }
        }
        int offset = 0;
        if (debug) Debug() << format_log("Initialization:", "0 %", true);
        for (int i = 0; i < N_MODELS; i ++) {
            model::mlp::predict(i, orfs, offset, gc_intv[i]);
            if (debug) {
                Debug() << format_log("Initialization:", Str((int)(i*1.67))+" %", true);
                if (i == N_MODELS-1) Debug() << '\n';
            }
            offset += gc_intv[i];
        }
        orfs.froze();
        
        return *this;
    }

    Clovers &Check_Table() {
        table = args["table"].as<str>();
        if (table == "auto") table = check_code();
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

    Clovers &Locate_Orfs() {
        int minlen = (int)args["minlen"].as<uint32_t>();
        if (minlen < MIN_LEN) {
            std::cerr << "Error: mininum gene length should not be less than 75 nt\n";
            ARG_ERR;
        }
        Debug() << format_log("Mininum Gene Length:", Str(minlen) + " nt");
        locate_orfs(STARTS, STOPS, minlen); 
        Debug() << format_log("Number of Candidate ORFs:", Str(orfs.size()));
        return *this;
    }

    Clovers &Revise_Orfs() {
        if (proc == "Metagenome" || args.count("bypass")) return *this;
        Debug() << format_log("Revising RBS Model:", "Round #0", true);
        auto maxiter = args["maxiter"].as<uint32_t>();
        model::pmm model(orfs, table);

        str mfile;
        if (args.count("rbs")) mfile = args["rbs"].as<str>();
        if (utils::file_exists(mfile)) model.load(mfile);
        else {
            int_arr seeds = orfs.filter(300, 0.5, true);
            int round = 0;
            if (seeds.size() >= 0) for (int order = 0; order <= 2; order ++) {
                model.reset(order);
                for (int i = 0; i < maxiter; i ++) {
                    Debug() << format_log("Revising RBS Model:", str("Round #")+Str(++round), true);
                    model.train(seeds, STARTS);
                    if (model.revise(seeds) > 0.99) break;
                }
            }
        }
        model.revise_all(); 
        Debug() << '\n'; revised = true;
        if (!utils::file_exists(mfile)) model.dump(mfile);
        return *this;
    }

    Clovers &Csvm_Scores() {
        if (proc != "Metagenome" && orfs.weight >= 0.8) {
            if (revised) zcurve::encode(orfs.concat());

            model::svm model(orfs);
            
            str mfile;
            if (args.count("train")) mfile = args["train"].as<str>();
            if (utils::file_exists(mfile)) model.load(mfile);
            else {
                Debug() << "Constructing Genome-Specific Gene Model (RBF-SVM)\n";
                model.sample(0.5, 0.05).train();
            }
            if (!utils::file_exists(mfile)) model.dump(mfile);
            model.predict();
            
        }
        orfs.yield(proc, args["thres"].as<float>(), genes);

        Debug() << format_log("Number of Putative Genes:", Str(genes.size()));
        return *this;
    }

    Clovers &Save_Result() {
        str format = args["format"].as<str>();

        if (args.count("output")) {
            str outfile = args["output"].as<str>();
            utils::write_output(outfile, genes, start_t, table, format);
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

    Clovers &Report_Gols() {
        if (!args.count("overlap") && !args.count("amino") && !args.count("nucl"))
            return *this;

        auto pairs = bioinfo::orfs::overlap(args["min-olen"].as<uint32_t>(), genes);

        Debug() << format_log("Number of Overlapped Pairs:", Str(pairs.size()/2));

        if (args.count("overlap")) {
            str outfile = args["overlap"].as<str>();
            utils::write_output(outfile, pairs, start_t, table, "gff");
        }

        if (args.count("amino")) {
            auto faaout = args["amino"].as<str>();
            utils::write_output(faaout, pairs, start_t, table, "faa");
        }

        if (args.count("nucl")) {
            auto fnaout = args["nucl"].as<str>();
            utils::write_output(fnaout, pairs, start_t, table, "fna");
        }
        
        return *this;
    }

    ~Clovers() { 
        using namespace std::chrono;
        auto end_t = SYSTEM_CLOCK now();
        auto duration = duration_cast<milliseconds>(end_t - start_t);
        auto seconds = duration.count() / 1000.0;
        Debug() << "\nFinished in " << std::fixed << std::setprecision(3) << seconds << " s\n";
    }
};

int main(int argc, char *argv[]) {
    std::ios::sync_with_stdio(false);
    
    Clovers(argc, argv)
        . Load_Genome()
        . Check_Table()
        . Locate_Orfs()
        . Init_Scores()
        . Revise_Orfs()
        . Csvm_Scores()
        . Save_Result()
        . Report_Gols()
    ;
}
