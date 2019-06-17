# PageRank-homework
大数据课程最后一项大作业实现pagerank算法

* 语言：C/C++/JAVA/Python
* 考虑dead ends 和spider trap 节点
* 优化稀疏矩阵
* 实现分块计算
* 程序需要迭代至收敛
* 不可直接调用接口，例如实现pagerank时，调用Python的networkx包
* 结果格式(.txt文件):[NodeID] [Score]

---

main.cpp       // 程序入口

Pagerank.h     // 头文件

Pagerank.cpp   // 具体实现

CMakeLists.txt // CMake文件

---

```C++
class Pagerank {
    private:
        map<int, set<int>> Matrix;           // 储存稀疏矩阵
        int MaxItem, MinItem;                // 最大编号和最小编号
        int MaxIterations;                   // 最大迭代轮数
        long double alpha;                  // 随机游走系数
        map<int, int> node_index;             // 初始编号与离散化后编号的映射
        map<int, int> index_node;             // 离散化后编号与初始编号的映射
        int total;                          // 离散化后的编号总数，离散化后编号为[0,total)
        bool trace;                         // 是否输出过程信息，调试用
        vector<long double> pr;             // PageRank数组
        vector<int> degree;                 // 离散化后每一个点的出度数
        int block_nums;                     // 分块块数
        int total_edges;                    // 输入数据一共多少边
        long double convergence;            // 收敛限界
    public:
        ... // 一些私有变量的get和set接口
        ... // 一些调试和缺省的默认函数，PrintItem等
        
        Pagerank(int _MinItem = -1, int _MaxItem = -1, int _MaxIterations = -1,
            long double _alpha = 0, bool _trace = false, int block_nums = 1, long double convergence = 1e-5); // 构造函数
        
        void insert_item(int a, int b); // 向Matrix插入边

        bool FormCheck(int tmp); // 简单的格式检查

        void readFile(string input); // 读入函数

        void PageRank(); // 核心算法

        void readForRange(string input); // 读入预处理
}
```
