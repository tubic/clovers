#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <filesystem>

namespace fs = std::filesystem;

class Orf {
public:
    std::string host;
    int start;
    int end;
    char strand;

    Orf(const std::string& h, int s, int e, char st)
        : start(s), end(e), strand(st)
    {
        size_t sep = h.rfind('.');
        if (sep != std::string::npos)
            host = h.substr(0, sep);
        else
            host = h;
    }

    inline bool eqstrand(const Orf& query) const {
        return host == query.host && strand == query.strand;
    }

    inline bool match5end(const Orf& query) const {
        if (strand == '+')
            return start == query.start;
        else
            return end == query.end;
    }

    inline bool match3end(const Orf& query) const {
        if (strand == '+')
            return end == query.end;
        else
            return start == query.start;
    }

    inline bool operator==(const Orf& query) const {
        return eqstrand(query) &&
               (match5end(query) || match3end(query));
    }
};

typedef std::vector<Orf> orf_array;

void parse(const fs::path file_path, orf_array &genes, orf_array &ctrls) {
    bool flag = file_path.extension() == ".predict";

    std::ifstream fin(file_path);
    if (!fin) {
        std::cerr << "Error: cannot open file "
                  << file_path << ".\n";
        return;
    }

    std::string line;
    std::string cur_host;

    while (std::getline(fin, line)) {
        bool orf = false;

        if (line.empty() || line[0] == '#' || line[0] == '\n')
            continue;

        std::string host;
        int cds_start;
        int cds_end;
        char strand;

        if (flag) {
            if (line[0] == '>') {
                int sep = line.find(' ');
                if (sep == std::string::npos) {
                    cur_host = line.substr(1);
                } else {
                    cur_host = line.substr(1, sep-1);
                }
                continue;
            }

            std::istringstream iss(line);
            std::string name, start_str, end_str, frame_str;
            if (!(iss >> name >> start_str >> end_str >> frame_str))
                continue;

            host = cur_host;
            cds_start = std::stoi(start_str);
            cds_end   = std::stoi(end_str);
            strand    = frame_str[0];
        } else {
            std::vector<std::string> fields;
            fields.reserve(9);

            size_t start_pos = 0;
            size_t pos;
            while ((pos = line.find('\t', start_pos)) != std::string::npos) {
                fields.emplace_back(line.substr(start_pos, pos - start_pos));
                start_pos = pos + 1;
            }
            fields.emplace_back(line.substr(start_pos));

            if (fields.size() < 8)
                continue;
            if (fields[2] != "CDS") {
                if (fields[2] == "ORF") orf = true;
                else continue;
            }

            host = fields[0];
            cds_start = std::stoi(fields[3]);
            cds_end   = std::stoi(fields[4]);
            strand    = fields[6][0];
        }

        if (orf) ctrls.emplace_back(host, cds_start, cds_end, strand);
        else genes.emplace_back(host, cds_start, cds_end, strand);
    }
}

int inter_count(const orf_array& pred_orfs,
                const orf_array& true_orfs)
{
    std::vector<bool> matched(true_orfs.size(), false);

    int intersection = 0;

    for (const auto& p : pred_orfs) {
        for (size_t i = 0; i < true_orfs.size(); ++i) {
            if (!matched[i] && p == true_orfs[i]) {
                matched[i] = true;
                ++intersection;
                break;
            }
        }
    }

    return intersection;
}

void print_usage(const char* prog) {
    std::string path(prog);
    size_t pos = path.find_last_of("/\\");
    if (pos != std::string::npos) {
        path = path.substr(pos + 1);
    }
    std::cerr << "Usage:\n" << path
              << " -p <prediction_gff_dir> -t <truth_gff_dir>\n";
}

int main(int argc, char* argv[]) {
    std::string pred_dir;
    std::string true_dir;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if ((arg == "-p" || arg == "--pred") && i + 1 < argc) {
            pred_dir = argv[++i];
        }
        else if ((arg == "-t" || arg == "--true") && i + 1 < argc) {
            true_dir = argv[++i];
        }
    }

    if (pred_dir.empty() || true_dir.empty()) {
        print_usage(argv[0]);
        return 1;
    }

    if (!fs::is_directory(pred_dir)) {
        std::cerr << "Error: prediction directory '"
                  << pred_dir
                  << "' does not exist.\n";
        return 1;
    }

    if (!fs::is_directory(true_dir)) {
        std::cerr << "Error: truth directory '"
                  << true_dir
                  << "' does not exist.\n";
        return 1;
    }

    std::vector<std::string> pred_files;

    for (const auto& entry : fs::directory_iterator(pred_dir)) {
        if (entry.is_regular_file()) {
            std::string ext = entry.path().extension().string();
            if (ext == ".predict" || ext == ".gff") {
                pred_files.emplace_back(entry.path().filename().string());
            }
        }
    }

    if (pred_files.empty()) {
        std::cerr << "Error: no .gff files found in prediction directory.\n";
        return 1;
    }

    std::cout << "Assembly\tPos\tTrue\tTrue_Pos\tNeg\tFalse_Neg\n";

    const size_t total_files = pred_files.size();

    for (size_t idx = 0; idx < total_files; ++idx) {
        const std::string& fname = pred_files[idx];

        std::cerr << "Processing "
                  << fname
                  << " (" << idx
                  << "/" << total_files
                  << ")\n";

        std::string acc = fname;
        if (acc.size() >= 4 &&
            acc.substr(acc.size() - 4) == ".gff")
        {
            acc.resize(acc.size() - 4);
        } else if (acc.size() >= 8 &&
            acc.substr(acc.size() - 8) == ".predict")
        {
            acc.resize(acc.size() - 8);
        }

        fs::path pred_path = fs::path(pred_dir) / fname;
        fs::path true_path = fs::path(true_dir) / (acc + ".gff");

        if (!fs::exists(true_path)) {
            std::cerr << "Truth file '"
                      << true_path.string()
                      << "' not found for prediction '"
                      << pred_path.string()
                      << "'\n";
        }

        orf_array pred_genes, pred_ctrls, true_genes, ignore;

        parse(pred_path.string(), pred_genes, pred_ctrls);
        parse(true_path.string(), true_genes, ignore);

        int pred_count = static_cast<int>(pred_genes.size());
        int true_count = static_cast<int>(true_genes.size());
        int nega_count = static_cast<int>(pred_ctrls.size());
        int tp = inter_count(pred_genes, true_genes);
        int fn = inter_count(pred_ctrls, true_genes);

        std::cout << acc << '\t' << pred_count << '\t' << true_count 
                  << '\t' << tp << '\t' << nega_count << '\t' << fn << '\n';
    }

    return 0;
}