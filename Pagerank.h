#ifndef PAGERANK_PAGERANK_H
#define PAGERANK_PAGERANK_H

#include <bits/stdc++.h>

using namespace std;

class Pagerank {
private:
    map<int, set<int>> Matrix;           // 储存稀疏矩阵，从a到set中节点的有向边 自动去重
    int MaxItem, MinItem;                // 最大编号和最小编号
    int MaxIterations;                   // 最大迭代轮数
    long double alpha;                  // 随机游走系数
    map<int, int> node_index;             // 初始编号与连续化编号的映射
    map<int, int> index_node;             // 连续化后编号与初始编号的映射
    int total;                          // 连续化后的编号总数，离散化后编号为[0,total)
    bool trace;                         // 是否输出过程信息，调试用
    vector<long double> pr;             // PageRank数组
    vector<int> degree;                 // 连续化后每一个点的出度数
    int block_nums;                     // 分块块数
    int total_edges;                    // 输入数据一共多少边
    long double convergence;            // 收敛限界
    bool separate_pr;
public:
    Pagerank(int _MinItem = -1, int _MaxItem = -1, int _MaxIterations = -1,
        long double _alpha = 0, bool _trace = false, int block_nums = 1, long double convergence = 1e-6, bool separate_pr = false);

    ~Pagerank() = default;

    long double getConvergence() const;

    void setConvergence(long double convergence);

    bool isTrace() const;

    void setTrace(bool trace);

    bool prSeparate() const;

    void setSeparatePr(bool sep);

    int getBlock_nums() const;

    void setBlock_nums(int block_nums);

    long double getAlpha() const;

    void setAlpha(long double alpha);

    int getMaxItem() const;

    void setMaxItem(int MaxItem);

    void setMinItem(int MinItem);

    int getMaxIterations() const;

    int getMinItem() const;

    void setMaxIterations(int MaxIterations);

    void PrintItem(int limit = 10);

    void insert_item(int a, int b);

    bool FormCheck(int tmp);

    void printBasicGraphInfo();

    void readFile(string input);

    void PageRank();

    // void Prepare();

    void readForRange(string input);    // 读取文件一次，找出最大的和最小的节点号

    void outputFile(string output);
};


#endif //PAGERANK_PAGERANK_H
