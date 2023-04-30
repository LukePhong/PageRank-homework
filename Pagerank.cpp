#include "Pagerank.h"

#define FormChecks(x) if(!FormCheck(x)) {cout<<"form error"<<endl;cerr<<"form error "<<x<<endl;exit(1);}

void Pagerank::insert_item(int a, int b) {
    FormChecks(a);  // 输入越界检查
    FormChecks(b);
    if (this->Matrix.find(a) != this->Matrix.end()) {
        Matrix[a].insert(b);    // 已有直接插入
    } else {
        Matrix.insert({a, set<int>()}); // 没有则新建一个set
        Matrix[a].insert({b});
    }
}


Pagerank::Pagerank(int _MinItem, int _MaxItem, int _MaxIterations, long double _alpha, bool _trace, int _block_nums,
                   long double _convergence, bool _separate_pr) : MaxItem(_MaxItem), MinItem(_MinItem), MaxIterations(_MaxIterations),
                                          alpha(_alpha), trace(_trace), block_nums(_block_nums),
                                          convergence(_convergence), separate_pr(_separate_pr) {
    this->Matrix.clear();
    this->node_index.clear();
    this->index_node.clear();
}

bool Pagerank::FormCheck(int tmp) {
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
    bool *used = new bool[MaxItem + 2]; // TODO: 是否存在该节点（过大内存开销可能）

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
    for (int i = MinItem; i <= MaxItem; i++) {  // 使所有点处于从0开始的连续编号上
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

    // 初始化page rank值、出度值
    pr.resize((size_t) total);
    degree.resize((size_t) total);
    // pr[0] = 1;
    // int maxdegree=0;
    vector<int> deadpoints;
    for (int i = 0; i < total; i++) {
        if (Matrix.find(index_node[i]) == Matrix.end()) {   // 这个点没有出边，是死点
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

    // 文件分块存储
    char file_name[100][300];
    ofstream blockOutput[100];
    ofstream prOutput[100];
    for (int i = 0; i < block_nums; i++) {
        sprintf(file_name[i], "block-%d.txt", i);
        blockOutput[i].open(file_name[i]);
    }
    for (int i = 0; i < block_nums; i++) {
        sprintf(file_name[i], "pr-%d.txt", i);
        prOutput[i].open(file_name[i]);
    }
    int tmpt = 0;
    for (auto i:Matrix) {
        for (auto j:i.second) {         // 对于每一个文件，from信息是全的，按照to的编号划分，也就是把邻接矩阵横向划分
            int to_index = node_index[j];
            int from_index = node_index[i.first];
            tmpt = to_index / (total / (block_nums)+1);
            blockOutput[tmpt] << from_index << " " << degree[from_index] << " " << to_index << endl; // TODO: 输出了重复的信息
        }
    }
    for (int i = 0; i < block_nums; i++) {
        blockOutput[i].close();
        prOutput[i].close();
    }

    // 输出重定位信息
    if(trace){
        blockOutput[0].open("node_index.txt");
        for (auto i:node_index) {
            blockOutput[0] << i.first << "->" << i.second << endl;
        }
        blockOutput[0].close();
    }

    cout << "Read completed" << endl;
}

bool Pagerank::isTrace() const {
    return trace;
}

void Pagerank::setTrace(bool trace) {
    Pagerank::trace = trace;
}

bool Pagerank::prSeparate() const {
    return separate_pr;
}

void Pagerank::setSeparatePr(bool sep) {
    Pagerank::separate_pr = sep;
}

int Pagerank::getBlock_nums() const {
    return block_nums;
}

void Pagerank::setBlock_nums(int block_nums) {
    Pagerank::block_nums = block_nums;
}

void Pagerank::PageRank() {
    // Prepare();
    long double beta = 1 - alpha;
    vector<long double> oldpr;  // 上一轮的page rank向量
    // 读取之前写好的pr初始值
    // for (size_t i = 0; i < block_nums; i++)
    // {
    //     /* code */
    // }
    
    oldpr.resize((size_t) total);
    for (auto &i:pr) {          // 初始化每个节点的PageRank值为1/total
        i = 1.0 / total;
    }

    long double diff = 1e9;     // 两轮之间差异
    int counts = 0;             // 已经进行迭代次数
    int blockSize = total/block_nums + 1;

    long double sum_pr = 1;             // 全部page rank值的和
    long double sum_dead = 0;           // 端点page rank值的和
    for (int i = 0; i < total; i++) {   // 统计sum_pr和sum_dead
        if (degree[i] == 0) {
            sum_dead += pr[i];
        }
        // sum_pr += pr[i];
    }
    // 两个sum可以放到循环外面 就不用再次统计了
    while (diff > convergence && counts < MaxIterations) {

        // 第一遍读取，需要统计全部的pr分块

        // long double sum_pr = 0;             // 全部page rank值的和
        // long double sum_dead = 0;           // 端点page rank值的和
        // for (int i = 0; i < total; i++) {   // 统计sum_pr和sum_dead
        //     if (degree[i] == 0) {
        //         sum_dead += pr[i];
        //     }
        //     sum_pr += pr[i];
        // }

        for (int i = 0; i < total; i++) {   // 更新oldpr向量为当前pr向量除以sum_pr，保证oldpr的总和为1
            oldpr[i] = pr[i];// / sum_pr;   // 没有必要
        }

        long double add_dead = alpha * sum_dead / total / sum_pr;   // 用于防止端点“泄露”page rank值
        long double add_t = beta * 1.0 / total;         // (1-d)*(1/N)用于引入随机性
        sum_pr = 0.0;
        sum_dead = 0;
        diff = 0.0;
//#pragma omp parallel  for
        for (int currBlock = 0; currBlock < block_nums; currBlock++) {
            char file_name[300];
            fstream Input;
            sprintf(file_name, "block-%d.txt", currBlock);
            Input.open(file_name);

            vector<long double> tmp_pr;     // 存储当前块的PageRank值的临时向量
            tmp_pr.resize((size_t) min(blockSize * (currBlock + 1), total) - blockSize * (currBlock) + 1);
            for (auto &i:tmp_pr) {
                i = 0.0;
            }
            int Min = blockSize * currBlock;// 当前block的范围
            int Max = min(blockSize * (currBlock + 1), total);

            int from_index, degree_from, to_index;
            while (Input >> from_index >> degree_from >> to_index) {
                tmp_pr[to_index - Min] += alpha * 1.0 / degree[from_index] * oldpr[from_index]; // 算出公式后半部分PR(to)那一横行
            }
            for (int i = Min; i < Max; i++) {
                tmp_pr[i - Min] += add_dead + add_t;    // 光进行上面的运算不算完 需要添加：1 随机漫步项 2 作为端点没有degree，但是他们所在的列不能按照0来计算
            }
            for (int i = Min; i < Max; i++) {
                pr[i] = tmp_pr[i - Min];    //更新pr
                sum_pr += pr[i];        // 重新统计sum_pr
                if(degree[i] == 0){  // 重新统计sum_dead
                    sum_dead += pr[i];
                }
            }
            Input.close();
        }
        // for (auto &i:pr) {  // page rank值归一化（没必要，）
        //     i /= sum_pr;
        // }
        for (int i = 0; i < total; i++) {
            diff += fabs(pr[i] - oldpr[i]);
        }
        cout<<"iteration: "<<counts<<"\tdiff: "<<diff<<endl;
        counts++;
    }

}

void Pagerank::outputFile(string output){
    // 将计算结果按PageRank值从大到小排序，并输出到文件
    vector<pair<long double, int>> to_out;
    to_out.resize(total);
    for (int i = 0; i < total; i++) {
        to_out[i] = {-pr[i], i};
    }
    sort(to_out.begin(), to_out.end());
    ofstream outFile;
    outFile.open(output);
    for (int i = 0; i < total; i++) {
        if (i >= 100) break;
        // cout << index_node[to_out[i].second] << " " << pr[to_out[i].second] << endl;
        outFile << index_node[to_out[i].second] << " " << pr[to_out[i].second] << endl;
    }
    outFile.close();
}

// void Pagerank::Prepare() {

// }

long double Pagerank::getConvergence() const {
    return convergence;
}

void Pagerank::setConvergence(long double convergence) {
    Pagerank::convergence = convergence;
}

void Pagerank::readForRange(string input) {
    ifstream Input;
    Input.open(input.c_str());
    int a=INT32_MAX, b=-1;
    int aa=INT32_MAX,bb=-1;
    while (Input >> a >> b) {
       aa=min(aa,min(a,b));
       bb=max(bb,max(a,b));
    }
    MinItem=aa;
    MaxItem=bb;
    Input.close();
}

void Pagerank::printBasicGraphInfo() {
    cout<<total<<" "<<total_edges<<" "<<MinItem<< " "<<MaxIterations<<" "<<MaxItem<<endl;
}
