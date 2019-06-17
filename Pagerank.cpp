#include "Pagerank.h"

#define FormChecks(x) if(!FormCheck(x)) {cout<<"form error"<<endl;cerr<<"form error "<<x<<endl;exit(1);}
using namespace std;

void Pagerank::insert_item(int a, int b) {
    FormChecks(a);
    FormChecks(b);
    if (this->Matrix.find(a) != this->Matrix.end()) {
        Matrix[a].insert(b);
    } else {
        Matrix.insert({a, set<int>()});
        Matrix[a].insert({b});
    }
}


Pagerank::Pagerank(int _MinItem, int _MaxItem, int _MaxIterations, long double _alpha, bool _trace, int _block_nums,
                   long double _convergence) : MaxItem(_MaxItem), MinItem(_MinItem), MaxIterations(_MaxIterations),
                                          alpha(_alpha), trace(_trace), block_nums(_block_nums),
                                          convergence(_convergence) {
    this->Matrix.clear();
    //this->MatrixCount.clear();
    this->node_index.clear();
    this->index_node.clear();
    // cout << "init completed" << endl;
}

bool Pagerank::FormCheck(int tmp) {
    if (MinItem == -1 || MaxItem == -1) return false;
    return tmp >= MinItem && tmp <= MaxItem;
}

void Pagerank::PrintItem(int limit) {
    int cc = 0;
    for (auto i : Matrix) {
        for (auto j :i.second) {
            cout << i.first << " " << j << endl;
        }
        cc++;
        if (cc >= limit) break;
    }
}

void Pagerank::setMaxIterations(int MaxIterations) {
    Pagerank::MaxIterations = MaxIterations;
}

int Pagerank::getMinItem() const {
    return MinItem;
}

int Pagerank::getMaxItem() const {
    return MaxItem;
}

void Pagerank::setMaxItem(int MaxItem) {
    Pagerank::MaxItem = MaxItem;
}

void Pagerank::setMinItem(int MinItem) {
    Pagerank::MinItem = MinItem;
}

int Pagerank::getMaxIterations() const {
    return MaxIterations;
}

long double Pagerank::getAlpha() const {
    return alpha;
}

void Pagerank::setAlpha(long double alpha) {
    Pagerank::alpha = alpha;
}

void Pagerank::readFile(string input) {
    total_edges = 0;
    if (MinItem == -1) {
        cerr << "You should set MinItem first" << endl;
        exit(1);
    }
    if (MaxItem == -1) {
        cerr << "You should set MaxItem first" << endl;
        exit(1);
    }
    if (MinItem > MaxItem) {
        cerr << "MaxItem should be larger than MinItem" << endl;
        exit(1);
    }
    bool *used = new bool[MaxItem + 2];

    ifstream Input;
    Input.open(input.c_str());
    int a, b;
    while (Input >> a >> b) {
        insert_item(a, b);
        used[a] = true;
        used[b] = true;
        total_edges++;
    }
    Input.close();
    int len = 0;
    for (int i = MinItem; i <= MaxItem; i++) {
        if (used[i]) {
            node_index[i] = len;
            index_node[len] = i;
            len++;
        }
    }
    delete[] used;
    total = len;
    if (trace) {
        for (auto i:node_index) {
            cout << i.first << " -> " << i.second << endl;
        }
    }

    pr.resize((size_t) total);
    degree.resize((size_t) total);
    pr[0] = 1;
    // int maxdegree=0;
    vector<int> deadpoints;
    for (int i = 0; i < total; i++) {
        if (Matrix.find(index_node[i]) == Matrix.end()) {
            if (trace) deadpoints.push_back(i);
            degree[i] = 0;
        } else {
            degree[i] = (int) Matrix[index_node[i]].size();
            // maxdegree=max(maxdegree,degree[i]);
        }
    }
    // cout<<"------------"<<endl;
    // cout<<maxdegree<<" "<<deadpoints.size()<<endl;
    if (trace) {
        cout << "Dead Points:" << endl;
        for (auto i:deadpoints) {
            cout << index_node[i] << endl;
        }
        cout << deadpoints.size() << endl;
        cout << endl;
    }
    char file_name[100][300];
    ofstream Output[100];
    for (int i = 0; i < block_nums; i++) {
        sprintf(file_name[i], "block-%d.txt", i);
        Output[i].open(file_name[i]);
    }
    int tmpt = 0;
    for (auto i:Matrix) {
        for (auto j:i.second) {
            int jj = node_index[j];
            int ii = node_index[i.first];
            tmpt = jj / (total / (block_nums)+1);
            // cout<<"total"<<total/(total/block_nums+1)<<" "<<(total-1)/(total/block_nums+1)<<endl;
            // if((total/block_nums)*tmpt>jj||(total/block_nums)*(tmpt+1)<=jj) cout<<jj<<endl;
            Output[tmpt] << ii << " " << degree[ii] << " " << jj << endl;
        }
    }
    for (int i = 0; i < block_nums; i++) {
        Output[i].close();
    }
    Output[0].open("node_index.txt");
    for (auto i:node_index) {
        Output[0] << i.first << "->" << i.second << endl;
    }
    Output[0].close();
    cout << "Read completed" << endl;
}

bool Pagerank::isTrace() const {
    return trace;
}

void Pagerank::setTrace(bool trace) {
    Pagerank::trace = trace;
}

int Pagerank::getBlock_nums() const {
    return block_nums;
}

void Pagerank::setBlock_nums(int block_nums) {
    Pagerank::block_nums = block_nums;
}

void Pagerank::PageRank() {
    Prepare();
    long double beta = 1 - alpha;
    vector<long double> oldpr;
    oldpr.resize((size_t) total);
    for (auto &i:pr) {
        i = 1.0 / total;
    }

    long double diff = 1e9;
    int counts = 0;
    int DD = total/block_nums+1;
    while (diff > convergence && counts < MaxIterations) {
        long double sum_pr = 0;
        long double sum_dead = 0;
        for (int i = 0; i < total; i++) {
            if (degree[i] == 0) {
                sum_dead += pr[i];
            }
            sum_pr += pr[i];
        }

        for (int i = 0; i < total; i++) {
            oldpr[i] = pr[i] / sum_pr;
        }

        long double add_dead = alpha * sum_dead / total / sum_pr;
        long double add_t = beta * 1.0 / total;
        sum_pr = 0.0;
        diff = 0.0;
//#pragma omp parallel  for
        for (int B = 0; B < block_nums; B++) {
            char file_name[300];
            fstream Input;
            sprintf(file_name, "block-%d.txt", B);
            Input.open(file_name);
            vector<long double> tmp_pr;
            tmp_pr.resize((size_t) min(DD * (B + 1), total) - DD * (B) + 1);
            for (auto &i:tmp_pr) {
                i = 0.0;
            }
            int Min = DD * (B);
            int Max = min(DD * (B + 1), total);

            int a, b, c;
            while (Input >> a >> b >> c) {
                tmp_pr[c - Min] += alpha * 1.0 / degree[a] * oldpr[a];
            }
            for (int i = Min; i < Max; i++) {
                tmp_pr[i - Min] += add_dead + add_t;
            }
            for (int i = Min; i < Max; i++) {
                pr[i] = tmp_pr[i - Min];
                sum_pr += pr[i];
            }
            Input.close();
        }
        for (auto &i:pr) {
            i /= sum_pr;
        }
        for (int i = 0; i < total; i++) {
            diff += fabs(pr[i] - oldpr[i]);
        }
        cout<<diff<<endl;
        counts++;
    }
    vector<pair<long double, int>> oo;
    oo.resize(total);
    for (int i = 0; i < total; i++) {
        oo[i] = {-pr[i], i};
    }
    sort(oo.begin(), oo.end());
    ofstream OO;
    OO.open("out.txt");
    for (int i = 0; i < total; i++) {
        if (i >= 100) break;
        cout << index_node[oo[i].second] << " " << pr[oo[i].second] << endl;
        OO << index_node[oo[i].second] << " " << pr[oo[i].second] << endl;
    }
    OO.close();
    cout << counts << endl;
}

void Pagerank::Prepare() {

}

long double Pagerank::getConvergence() const {
    return convergence;
}

void Pagerank::setConvergence(long double convergence) {
    Pagerank::convergence = convergence;
}

void Pagerank::readForRange(string input) {
    ifstream Input;
    Input.open(input.c_str());
    int a=9999999, b=-1;
    int aa=999999,bb=-1;
    while (Input >> a >> b) {
       aa=min(aa,min(a,b));
       bb=max(bb,max(a,b));
    }
    MinItem=aa;
    MaxItem=bb;
    Input.close();
}

void Pagerank::print() {
    cout<<total<<" "<<total_edges<<" "<<MinItem<< " "<<MaxIterations<<" "<<MaxItem<<endl;

}
