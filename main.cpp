#include <iostream>
#include <fstream>
#include <chrono>
#include <algorithm>
#include <vector>
#include <set>
#include "Pagerank.h"
#include <cstring>

using namespace std;
using namespace chrono;

void usage() {
    cerr << "pagerank [-t] [-a alpha ] [-c convergence] [-b blocknums]"
         << "[-m max_iterations] <graph_file>" << endl
         << " -t enable tracing " << endl
         << "integer vertex names" << endl
         << " -a alpha" << endl
         << "    default is 0.85 " << endl
         << " -c convergence" << endl
         << "    the convergence criterion default is 1e-5" << endl
         << " -m max_iterations" << endl
         << "    maximum number of iterations to perform" << endl
         << " -b blocknums" << endl
         << "    the number of blocknums, for some reason, the number shouldn't be larger than 100" << endl
         << " the output file named out.txt" << endl;
}

int main(int argc, char **argv) {
    string input = "";
    Pagerank s;
    s.setMaxIterations(1000);
    s.setAlpha(0.85);
    s.setBlock_nums(4);
    for (size_t i = 1; i < argc; i++) {
        if (string(argv[i]) == "-t") {
            s.setTrace(true);
        } else if (string(argv[i]) == "-m") {
            int tmpt = stoi(string(argv[i + 1]));
            if (tmpt <= 0) {
                usage();
                exit(1);
            }
            s.setMaxIterations(tmpt);
            i++;
        } else if (string(argv[i]) == "-a") {
            long double tmpa = atof(argv[i + 1]);
            if (tmpa > 1 || tmpa < 0) {
                usage();
                exit(1);
            }
            s.setAlpha(tmpa);
            i++;
        } else if (string(argv[i]) == "-b") {
            int tmpt = stoi(string(argv[i + 1]));
            if (tmpt <= 0 || tmpt >= 100) {
                usage();
                exit(1);
            }
            s.setBlock_nums(tmpt);
            i++;
        } else if (string(argv[i]) == "-c") {
            long double tmpc = atof(argv[i + 1]);
            s.setConvergence(tmpc);
            i++;
        } else if (i == argc - 1) {
            input = argv[i];
        } else {
            usage();
            exit(1);
        }
    }
    if (input == "") {
        usage();
        exit(1);
    }
    // input="WikiData.txt";
    s.readForRange(input);
    // s.setConvergence(1e-7);
    //s.setTrace(true);
    //cout<<"ss"<<endl;
    auto start = steady_clock::now();
    s.readFile(input);
    s.print();
    s.PageRank();
    auto end = steady_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    cout << "use "
         << double(duration.count()) * microseconds::period::num / microseconds::period::den
         << " s" << endl;

    return 0;
}