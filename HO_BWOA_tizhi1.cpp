#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include <sstream>
#include<random>
#include <limits>
#include<algorithm>
#include<numeric>
#include<cmath>
#include "Variant.hpp"
using namespace std;

/*生成随机排列的函数（与matlab中的randperm是同一个功能）*/
std::vector<int> randperm(int n, int k) {
    std::vector<int> indices(n);
    std::iota(indices.begin(), indices.end(), 0);  // 填充 0,1,2,...,n-1
    std::shuffle(indices.begin(), indices.end(), std::default_random_engine(std::random_device{}()));
    if (k < n) {
        indices.resize(k);
    }
    return indices;
}

/*模板化求和函数：支持 int、double、float 等任意数值类型*/
template <typename T>  // T 为通用类型参数（可推导为 int、double 等）
T sum_vec1d(const std::vector<T>& vec_1d) {
    return accumulate(vec_1d.begin(), vec_1d.end(), T{});
}

/*处理二维double vector：n* m → 1 * m（第一个维度求和，保留第二个维度），与matlab中的sum类似。*/
template <typename T>   //支持int,double,float类型
vector<T> sum_vec2d(const std::vector<std::vector<T>>& vec_2d) {
    // 边界检查：输入不能为空
    if (vec_2d.empty()) {
        throw std::invalid_argument("2D vector is empty (cannot sum first dimension)");
    }

    // 确定结果的列数（第二个维度大小），并初始化结果为0
    size_t cols = vec_2d[0].size();
    std::vector<T> result(cols, T{});

    // 遍历每个行（第一个维度），对应列相加
    for (const auto& row : vec_2d) {
        // 检查所有行的列数一致
        if (row.size() != cols) {
            throw std::invalid_argument("2D vector has inconsistent column counts");
        }
        for (size_t j = 0; j < cols; ++j) {
            result[j] += row[j];
        }
    }

    return result;
}

/*处理三维double vector：n* m* p → m* p（第一个维度求和，保留第二、三维度，与matlab中的sum类似。）*/
template <typename T>   //支持int,double,float类型
std::vector<std::vector<T>> sum_vec3d(const std::vector<std::vector<std::vector<T>>>& vec_3d) {
    // 边界检查：输入不能为空
    if (vec_3d.empty()) {
        throw std::invalid_argument("3D vector is empty (cannot sum first dimension)");
    }

    // 确定结果的行数（第二个维度）和列数（第三个维度），初始化结果为0
    size_t rows = vec_3d[0].size();
    size_t cols = vec_3d[0][0].size();
    std::vector<std::vector<T>> result(rows, std::vector<T>(cols, T{}));

    // 遍历每个"切片"（第一个维度），对应位置相加
    for (const auto& slice : vec_3d) {
        // 检查所有切片的行数一致
        if (slice.size() != rows) {
            throw std::invalid_argument("3D vector has inconsistent row counts across slices");
        }
        for (size_t i = 0; i < rows; ++i) {
            const auto& row = slice[i];
            // 检查所有行的列数一致
            if (row.size() != cols) {
                throw std::invalid_argument("3D vector has inconsistent column counts across rows");
            }
            for (size_t j = 0; j < cols; ++j) {
                result[i][j] += row[j];
            }
        }
    }

    return result;
}

/*支持一维vector，数据类型支持int,double,float*/
template <typename T>
bool writeMatrixToFile(const string& file_path, const vector<T>& vec)
{
    ofstream file(file_path);
    if (!file.is_open()) {
        cerr << "错误：无法打开文件" << file_path << endl;
        return false;
    }
    
    //遍历vector的每个元素
    for (size_t i = 0; i < vec.size(); i++)
    {
        file << vec[i];
        if (i != vec.size() - 1) {
            file << " ";
        }
    }
    file << endl;
    
    file.close();
    cout << "一维vector已成功写入=>" << file_path << endl;
    return true;
}

/*用于评估搜索代理优劣的目标函数*/
fobjReturn fobj0811v3(vector<double> Pos, int dim, vector<int> para_GRemain, vector<int> para_Galloc, int user_num, int sate_num, int slot_num, int Th_G)
{
    //初始化返回结构体
    fobjReturn fobj_return(user_num, sate_num, slot_num, 0.0);
    vector<vector<vector<double>>> M_3D = make_vec3(user_num, sate_num, slot_num, 0.0);
    //将一维数组重建为三维数组（主意这里matlab的reshape是列主序！！！）
    for (int index = 0; index < Pos.size(); index++)
    {
        int i = index / (sate_num * slot_num);
        int j = (index % (sate_num * slot_num)) / slot_num;
        int k = (index % (sate_num * slot_num)) % slot_num;
        M_3D[i][j][k] = Pos[index];
    }

    int lambda1 = 1e4;  //成功率项权重
    int lambda2 = 1;    //负载均衡项权重
    int lambda3 = 1;    //剩余时间项权重

    //处理约束1： 群组接入1颗星
    //解释：惩罚项pen1，超出1颗星的群组按照超出的数量的平方累加惩罚值，得到所有群组的总惩罚值。
    //作用：惩罚违反约束的解，让优化过程倾向于满足约束。
    double pen1 = 0;    
    int muy1 = 1e5;
    //vector<double> ceshi1 = sum_vec2d(M_3D[0]);
    //double ceshi2 = sum_vec1d(sum_vec2d(M_3D[0]));
    for (int g = 0; g < user_num; g++)
    {
        vector<double> M_3D_sum1 = sum_vec2d(M_3D[g]);  //对第二维（卫星）求和
        double M_3D_sum2 = sum_vec1d(M_3D_sum1);    //对第二、三维求和
        bool isExceed = 0;
        for (int i = 0; i < M_3D_sum1.size(); i++)
        {
            if (M_3D_sum1[i] > 1)
            {
                isExceed = 1;
            }
        }
        if (isExceed)
        {
            pen1 = pen1 + pow((M_3D_sum2 - 1), 2);
        }
    }
    pen1 = pen1 * muy1;

    //处理约束2：单星负载限制，先简单用用户数量估计。
    double pen2 = 0;
    double muy2 = 1;
    vector<int> loadS(sate_num, 0);    //记录sat负载，超载则直接记当前负载为门限值
    vector<bool> isLoadOver(sate_num, 0);  //标记sat是否超负载
    for (int s = 0; s < sate_num; s++)
    {
        for (int g = 0; g < user_num; g++)
        {
            //确保两次切换之间仅分配一个时隙给该群组，或者无分配，无分配则无需加该群组用户。
            if (sum_vec1d(M_3D[g][s]) == 1)
            {
                loadS[s] += para_Galloc[g]; //加该群组用户数
            }
        }
        if (loadS[s] > Th_G) //如果负载超出设定门限
        {
            isLoadOver[s] = 1;
            loadS[s] = Th_G;
            pen2 += pow((loadS[s] - Th_G), 2);
        }
    }
    pen2 = pen2 * muy2;

    //处理约束3：单星各时隙接入群组数量限制
	double pen3 = 0;
	double muy3 = 1e5;
    for (int s = 0; s < sate_num; s++)
    {
        for (int t = 0; t < slot_num; t++)
        {
            //计算该星该时隙接入了多少群组
            double sum_group_num = 0;
            for (int g = 0; g < user_num; g++)
            {
                sum_group_num += M_3D[g][s][t];
            }
            if (sum_group_num > 1)
            {
                pen3 += pow((sum_group_num - 1), 2);
            }
        }
    }
    pen3 = pen3 * muy3;


    //step1：处理约束1:--若决策指示多星/多时隙接入，保留时机最靠前判决
    //过程：对于超出约束的群组，M_3D中仅保留时间最靠前的接入情况（置1），其他都置0。
    for (int gID = 0; gID < user_num; gID++)
    {
        if (sum_vec1d(sum_vec2d(M_3D[gID])) > 1)   //单群组存在多星/多时隙接入情况
        {
            //先遍历时间，则可以快速找到时机最靠前判决
            bool hasFind = 0;
            for (int t = 0; t < slot_num; t++)
            {
                for (int s = 0; s < sate_num; s++)
                {
                    if (M_3D[gID][s][t] == 1){  //找到了
                        //保留时机最靠前判决
                        M_3D[gID] = make_vec2<double>(sate_num, slot_num, 0.0);
                        M_3D[gID][s][t] = 1;
                        hasFind = 1;
                        break;
                    }
                }
                if (hasFind) {
                    break;
                }
            }
        }
    }

    //step2：处理约束3--若s星t时隙接入群组数超出，则保留用户数较多的群组。
    for (int s = 0; s < sate_num; s++)
    {
        for (int t = 0; t < slot_num; t++)
        {
            //计算该星该时隙接入了多少群组
            double sum_group_num = 0;
            int max_user = 0;
            double max_user_group_idx = 0;  //存储用户数最多的群组的ID
            for (int g = 0; g < user_num; g++)
            {
                double pos = M_3D[g][s][t];
                sum_group_num += pos;
                if (pos == 1) {
                    //找最大值及其位置
                    if (para_Galloc[g] > max_user) {
                        max_user = para_Galloc[g];
                        max_user_group_idx = g;
                    }
                }
                //顺便将原来M_3D中s星t时隙的数据清0，方便之后仅保留用户数较多的群组
                M_3D[g][s][t] = 0;
            }
            if (sum_group_num >= 1) //只要原来有接入群组，就：
            {
                //保留用户数较多的群组
                M_3D[max_user_group_idx][s][t] = 1;
            }
        }
    }

    //step3：处理约束2-单星容量，记录在loadS（实际和上面约束2计算过程类似）
    fill(loadS.begin(), loadS.end(), 0.0);
    fill(isLoadOver.begin(), isLoadOver.end(), 0.0);
    for (int s = 0; s < sate_num; s++)
    {
        for (int g = 0; g < user_num; g++)
        {
            //确保两次切换之间仅分配一个时隙给该群组，或者无分配，无分配则无需加该群组用户。
            if (sum_vec1d(M_3D[g][s]) == 1)
            {
                loadS[s] += para_Galloc[g]; //加该群组用户数
            }
        }
        if (loadS[s] > Th_G) //如果负载超出设定门限
        {
            isLoadOver[s] = 1;
            loadS[s] = Th_G;    //超载直接记为门限
            //pen2 += pow((loadS[s] - Th_G), 2);
        }
    }

    //计算成功率项【待修改】
    int HO_success = sum_vec1d(loadS);

    //计算负载均衡项 越均衡值越小（其实是算负载数据的方差）
    double HO_load_avg = HO_success / sate_num;
    double HO_load = 0; //负载均衡项
    for (int s = 0; s < sate_num; s++)
    {
        HO_load += pow((loadS[s] - HO_load_avg), 2);
    }

    //【新增】计算剩余时间项，剩余时间越短的越早切换 越小越好
    //过程：遍历群组，找到所有接入位置（ss星tt时隙），由于之前处理了约束3，使得
    //单星单时隙仅有一个群组，所以只会找到一个接入位置。累积（该时间*10-群组剩余时间）。
    double HO_T = 0;
    for (int g = 0; g < user_num; g++)
    {
        //寻找该群组的接入位置
        bool hasFind = 0;
        for (int ss = 0; ss < sate_num; ss++)
        {
            for (int tt = 0; tt < slot_num; tt++)
            {
                if (M_3D[g][ss][tt] == 1)   //按理来说只有一个
                {
                    HO_T += abs((tt + 1) * 10 - para_GRemain[g]);
                    break;
                }
            }
            if (hasFind == 1) {
                break;
            }
        }
    }

    //load越小越好
    double objf = -lambda1 * HO_success + lambda2 * HO_load + lambda3 * HO_T;

    fobj_return.M_3D = M_3D;    //新位置矩阵
    fobj_return.HO_Sucnow = HO_success; //成功率项
    fobj_return.fitness = objf; //综合各因素计算的适应度值
    return fobj_return;
}



/*二进制鲸鱼优化算法*/
void BWOA()
{
    /*HO_BWOA_tizhi1.m*/
    /*读取*/
    //存储读取数据
    int user_num = 0;
    vector<int> para_GueN;  //群组内用户数
    vector<int> para_GRemain;   //群首剩余时间
    /*备注：体制修改参数。需要修改的是这个读取的文件的名字
    体制1：group_input/group1_group_info.txt
    体制2：group_input/group2_group_info.txt
    体制3：group_input/group2_group_info.txt*/
    //读取群组信息
    string group_info_txt = "group_input/group1_group_info.txt";
    ifstream inFile(group_info_txt, ios::in);
    if(!inFile.is_open())
    {
        cout << "Failed to open group_input/group1_group_info.txt" << endl;
        return;
	}
    string line;
    while (getline(inFile, line))
    {
        user_num++;
        istringstream iss(line);
        int val1, val2;
        if (iss >> val1 >> val2)
        {  
            para_GueN.push_back(val1);  //存储群组内用户
            para_GRemain.push_back(val2);   //存储群首剩余时间
        }
    }

    /*备注：体制修改参数。需要修改以下四个参数：para_Gsize，sate_num，para_T，frequency_type
    体制1:16,4,250,1
    体制2:8,1,300,2
    体制3:8,1,300,3*/
    /*参数*/
    const int frequency_type = 1;    //体制类型，1为体制1，2为体制2，3为体制3
    int para_Gsize = 16;    //【调整】群组内用户数量
    int sate_num = 4;  //【调整】卫星数量，假设分布3/2/1
    int para_T = 250;	//【调整】判决间隔250ms
    int para_preamble = 8; //【调整】前导资源数量，码分复用（单时隙可分配资源数量）
    int slot_num = para_T / 10;  //时隙数量（每个时隙10ms），时分复用
    int Th_G = 400;     //单星负载限制（用户数量）
    int SAT_begin_index[3] = {0,4,5};   //各个体制的卫星编号起始索引
    //vector<int> para_Galloc(user_num, para_Gsize); //创建一个长度为群组数量的vector,值都是给定的群组内用户数量。
    //各群组用户数量，先用14-16随机数生成代替（现在是从文件中读取）
    // 创建随机数生成器
    /*std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis_GueN(15, 17);
    std::vector<int> para_GueN(user_num);
    for (int i = 0; i < user_num; i++) {
        para_GueN[i] = dis_GueN(gen);
    }*/
    //各群组剩余时间，随机数生成（现在是从文件中读取）
    /*std::uniform_int_distribution<> dis_GRemain(400, 500);
    std::vector<int> para_GRemain(user_num);
    for (int i = 0; i < user_num; i++) {
        para_GRemain[i] = dis_GRemain(gen);
    }*/
    int BWOA_num = 1;
    int dim = user_num * sate_num * slot_num;  //群组数量*卫星数量*时隙数
    int SearchAgents_no = 80;   //使用80个搜索代理（个体）来寻找最优解。
    int Max_iter = 30;  //迭代次数（30，测试时暂改1）

    /*初始化*/
    //创建一个dim列的全零向量，用于记录当前群组在当前时隙最优解（所有搜索个体的leader）的位置。
    vector<double> Leader_pos(dim, 1);
    //记录Leader的当前和上一轮的的得分。
    double Leader_score = std::numeric_limits<double>::infinity();  //初始化为正无穷
    double Leader_score_pre = Leader_score;
    //停止算法运行的收敛容忍度
    double delta = 1e-6;
    int todoTol = 0;
    int Flag = 0;
    //（随机）初始化各个搜索个体的初始搜索位置。（之后用于计算每个卫星每个时隙的每个群组的最优解）
    vector<vector<double>> Positions(SearchAgents_no, vector<double>(dim));
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis_pos(0.0, 1.0);
    for (int i = 0; i < SearchAgents_no; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            if (dis_pos(gen) <= 0.5) {
                Positions[i][j] = 0;
            }
            else {
                Positions[i][j] = 1;
            }
        }
    }
    //初始化时生成一个较可行解M_decision（遍历时隙和卫星，按顺序给所有群组先初始分配一份资源，每个都分配过一份或资源都分完了则结束）。
    //std::vector<std::vector<std::vector<double>>> RRR(user_num, std::vector<std::vector<double>>(sate_num, std::vector<double>(slot_num, 0.0)));
    auto RRR = make_vec3<double>(user_num, sate_num, slot_num, 0.0);   //三维
    int Gid = 0;
    for (int i = 0; i < slot_num; i++)
    {
        for (int j = 0; j < sate_num; j++)
        {
            RRR[Gid][j][i] = 1;
            Gid = Gid + 1;
            if (Gid >= user_num) { //直到给所有群组都分配过资源后，分配完毕。
                break;
            }
        }
        if (Gid >= user_num) {
            break;
        }
    }
    //新建一个和RRR等维度的RRR_new
    auto RRR_new = make_vec3<double>(user_num, sate_num, slot_num, 0.0);   //三维
    int K = std::min(user_num, sate_num * slot_num);
    for (int k = 0; k < int(SearchAgents_no / 2); k++)
    {
        vector<int> P = randperm(user_num, K); //随机打乱1~user_num的数，取前K个。
        //打乱是为了将前面按群组顺序分配的资源变为随机分配给各个群组，直到分配完资源，或者已经分配给所有群组一份。
        for (int i = 0; i < K; i++)
        {
            // 将 RRR[i] 的内容复制到 RRR_new[P[i]]
            RRR_new[P[i]] = RRR[i];
        }
        //将RRR_new抻平为一维向量（注意matlab的reshape是列主序，要反过来内部先遍历行）
        int index = 0;
        for (int m = 0; m < user_num; m++)
        {
            for (int n = 0; n < sate_num; n++)
            {
                for (int q = 0; q < slot_num; q++)
                {
                    Positions[k][index++] = RRR_new[m][n][q];
                }
            }
        }
    }
    vector<double> Convergence_curve(Max_iter, 0);
    auto D_X = make_vec2<double>(SearchAgents_no, dim, 0.0);    //二维
    
    //循环计数器
    int iter = 0;

    //主循环
    vector<double> old_pos = Positions[0];
    while (iter < Max_iter && Flag <= 3)
    {
        for (int i = 0; i < Positions.size(); i++)  //遍历每个SearchAgent
        {
            // 对所有参与搜索的个体（search agent），逐一计其在目标函数（object function），也就是适应度函数中的取值。（用于评估搜索个体的优劣）
            fobjReturn fobj_return = fobj0811v3(Positions[i], dim, para_GRemain, para_GueN, user_num, sate_num, slot_num, Th_G);
            vector<vector<vector<double>>> M_3D = fobj_return.M_3D;
            double HO_Sucnow = fobj_return.HO_Sucnow;
            double fitness = fobj_return.fitness;

            //更新Leader
            //cout << fitness << endl;
            if (fitness < Leader_score)   //fitness越小越好，故最小化得分以寻找最优解，解决优化问题。
            {
                Leader_score = fitness; //更新最佳得分
                //最佳位置也要更新成新的M_3D抻平后的结果
                int pos_index = 0;
                for (int g = 0; g < user_num; g++)
                {
                    for (int s = 0; s < sate_num; s++)
                    {
                        for (int t = 0; t < slot_num; t++)
                        {
                            Leader_pos[pos_index++] = M_3D[g][s][t];
                        }
                    }
                }
                int HO_Success = HO_Sucnow; //似乎只是存出来但是没有用到
            }
        }
        if (old_pos != Leader_pos)
        {
            printf("第%d次迭代Pos更新\n", iter);
            old_pos = Leader_pos;
        }

        //从2到0线性递减。
        double a = 2 - iter * (2 / Max_iter);  

        //从-1到-2线性递减，用来计算t
        double a2 = -1 + iter * ((-1) / Max_iter);

        //更新每个搜索代理的位置Positions
        for (int i = 0; i < SearchAgents_no; i++)   //Positions.size()
        {
            double r1 = dis_pos(gen);   //[0-1]随机数
            double r2 = dis_pos(gen);   //[0-1]随机数

            double A = 2 * a * r1 - a;
            double C = 2 * r2;

            double p = dis_pos(gen);

            for (int j = 0; j < dim; j++)   //Positions[0].size()
            {
                if (p < 0.5)
                {
                    //选择“随机猎物搜索+收缩包围”策略（迭代前期abs(A)多大于1，后期多小于1。前期多探索，后期多开发。）
                    if (abs(A) >= 1)
                    {
                        //搜索猎物（探索阶段），该鲸鱼随机选择跟随一只鲸鱼，全局范围寻找最优“猎物”，避免陷入局部最优解。
                        int rand_leader_index = floor(SearchAgents_no * dis_pos(gen));
                        vector<double> X_rand = Positions[rand_leader_index];
                        double D_X_rand = abs(C * X_rand[j] - Positions[i][j]); //C是随机扰动参数[0~2]
                        D_X[i][j] = D_X_rand;   //存储当前鲸鱼到到随机参考鲸鱼对应维度的距离。

                    }
                    else if (abs(A) < 1)
                    {
                        //收缩包围（开发阶段），该鲸鱼向头鲸靠近，在最优解附近精细化搜索。
                        double D_Leader = abs(C * Leader_pos[j] - Positions[i][j]); //C是随机扰动参数[0~2]
                        D_X[i][j] = D_Leader;   //存储当前鲸鱼到头鲸之间的距离。
                    }
                }
                else if (p >= 0.5)
                {
                    //选择“螺旋路径逼近策略”
                    //该鲸鱼绕头鲸做螺旋逼近，在最优解附近做弧形搜索，兼顾开发和探索。
                    double distance2Leader = abs(Leader_pos[j] - Positions[i][j]);
                    D_X[i][j] = distance2Leader;    //存储当前鲸鱼到头鲸之间的距离。
                }

                //在9种通过实际距离D_X和系数A转化成[0~1]区间概率的转移函数中选择一种。
                //s越大，原位置Pos置1或翻转的概率越大。
                double s = 0;
                if (BWOA_num == 1)
                {
                    s = 1 / (1 + exp(-2 * A * D_X[i][j]));
                }
                else if (BWOA_num == 2)
                {
                    s = 1 / (1 + exp(-A * D_X[i][j]));
                }
                else if (BWOA_num == 3)
                {
                    s = 1 / (1 + exp(-A * D_X[i][j]) / 2);
                }
                else if (BWOA_num == 4)
                {
                    s = 1 / (1 + exp(-A * D_X[i][j]) / 3);
                }
                //对于前四种，采用S形转移函数
                if (BWOA_num <= 4)
                {
                    if (dis_pos(gen) < s) {
                        Positions[i][j] = 1;
                    }
                    else {
                        Positions[i][j] = 0;
                    }
                }

                if (BWOA_num == 5)
                {
                    s = abs(erf((sqrt(pi) / 2) * A * D_X[i][j]));   //v1转移函数
                }
                else if (BWOA_num == 6)
                {
                    s = abs(tanh(A * D_X[i][j]));   //v2
                }
                else if (BWOA_num == 7)
                {
                    s = abs(A * D_X[i][j] / sqrt(1 + A * pow(D_X[i][j], 2)));   //v3
                }
                else if (BWOA_num == 8)
                {
                    s = abs((2 / pi) * atan((pi / 2) * A * D_X[i][j])); //v4
                }
                else if (BWOA_num == 9)
                {
                    s = 1 / (1 + exp(-10 * (A * D_X[i][j]) - 0.5));
                }
                //对于后5种，采用V形转移函数
                if (BWOA_num > 4 && BWOA_num <= 9)
                {
                    if (dis_pos(gen) < s) {
                        Positions[i][j] = (Positions[i][j] == 1) ? 0 : 1;   //0-1翻转
                    }
                }
            }

            //【尝试】若1的数量多于sate_num * slot_num，则随机删去多余的1。
            double pos1sum = sum_vec1d(Positions[i]);
            if (pos1sum > sate_num * slot_num)
            {
                // 先找到Positions[i]中所有值为1的位置（索引），随机取其中的sate_num * slot_num个置1。
                vector<int> index;
                int pos_size = Positions[i].size();
                for (int k = 0; k < pos_size; k++)
                {
                    // 存储所有值为1的位置
                    if (Positions[i][k] == 1) {
                        index.emplace_back(k);
                    }
                }
                // 在0~index.size()的数中，随机取其中sate_num * slot_num个。
                vector<int> indexnew = randperm(index.size(), sate_num * slot_num);
                // 先将原来的Positions[i]中的数置0。
                vector<double> pos_zero(Positions[i].size(), 0.0);
                Positions[i] = pos_zero;
                // 再遍历indexnew，仅给指定sate_num * slot_num个位置赋值1。
                int ceshi_size = indexnew.size();
                vector<double> PPP(Positions[i].size(), 0.0);
                for (int k = 0; k < sate_num * slot_num; k++)
                {
                    PPP[index[indexnew[k]]] = 1;
                }
                int ceshi_1_num = sum_vec1d(PPP);
                Positions[i] = PPP;
            }
        }
        //记录每次迭代后的头鲸的适应度函数的得分（越小越好）
        Convergence_curve[iter] = Leader_score;
        //迭代次数增加1
        iter = iter + 1;

        fobjReturn fobj_return = fobj0811v3(Leader_pos, dim, para_GRemain, para_GueN, user_num, sate_num, slot_num, Th_G);
        printf("迭代次数%d,Score=%f,HO_success=%d\n", iter, Leader_score, (int)fobj_return.HO_Sucnow);

        if (todoTol == 1 && abs(Leader_score - Leader_score_pre) < delta)
        {
            Flag += 1;
            Convergence_curve = vector<double>(Convergence_curve.begin(), Convergence_curve.begin() + iter);
        }
        Leader_score_pre = Leader_score;
    }

    //将最优的分配结果存起来
    vector<vector<vector<double>>> ansHoD1 = make_vec3(user_num, sate_num, slot_num, 0.0);
    for (int index = 0; index < Leader_pos.size(); index++)
    {
        int i = index / (sate_num * slot_num);
        int j = (index % (sate_num * slot_num)) / slot_num;
        int k = (index % (sate_num * slot_num)) % slot_num;
        ansHoD1[i][j][k] = Leader_pos[index];
    }

    /*HO_decision_tizhi1.m*/
    //当前体制卫星编号的起始索引
    int deltaSAT = SAT_begin_index[frequency_type - 1];
    //读取各UE分群结果与剩余时间
    string user_info_txt = "group_input/all_users_info.txt";
    //int all_UEnum = 1500;
    vector<vector<int>> GID;
    ifstream inFile2(user_info_txt, ios::in);
    if (!inFile2.is_open())
    {
        cout << "Failed to open all_users_info.txt" << endl;
        return;
    }
    string line2;
    int line_idx = 0;
    while (getline(inFile2, line2))
    {
        istringstream iss(line2);
        int val1, val2; //val1：体制1/2/3；val2：分群结果。
        if (iss >> val1 >> val2)
        {
            GID.push_back({ val1,val2 });
            line_idx++;
        }
    }
    int all_UEnum = line_idx;   //用户数量

    //提取群组判决结果
    vector<double> G_Tsat(user_num, 0.0);
    vector<double> G_Ttime(user_num, 0.0);
    //遍历每个g，提取各个群组有效分配（s星和t时隙）结果到两数组。
    for (int g = 0; g < user_num; g++)
    {
        bool isFind = 0;
        int first_t_s_idx[2];
        //遍历这个群组的分配情况，如果能找得到则只取第一个（优先遍历时隙）。
        for (int t = 0; t < slot_num; t++)
        {
            for (int s = 0; s < sate_num; s++)
            {
                if (ansHoD1[g][s][t] == 1) {
                    isFind = 1;
                    first_t_s_idx[0] = t;//存储下标
                    first_t_s_idx[1] = s;//存储下标
                }
            }
        }

        //如果能找到，仅取第一个；找不到则存-1
        if (isFind) {
            G_Ttime[g] = first_t_s_idx[0];
            G_Tsat[g] = first_t_s_idx[1];
        }
        else {
            G_Ttime[g] = -1;
            G_Tsat[g] = -1;
        }
    }

    //匹配各用户判决结果
    // 说明：基于群组分配结果来匹配，但是仅在用户所在群组有资源的情况下
    // 对体制为1的用户分配。
    vector<int> U_Tsat(all_UEnum, -1);  //目标卫星
    vector<int> U_Ttime(all_UEnum, -1); //接入发起时间
    vector<int> U_prb(all_UEnum, -1);   //接入前导
    vector<int> gtime(user_num, 0);//记录每个群组内已分配接入时间的用户数量（前导资源有限）
    vector<int> U_no;   //记录没有被分配资源的用户id
    for (int u = 0; u < all_UEnum; u++)
    {
        if (GID[u][0] == frequency_type) //仅对体制为1的用户分配
        {
            int nowg = GID[u][1];   //UE所属群组的index
            gtime[nowg] += 1;

            //判断该用户所在群组是否有被分配到切换资源
            if (G_Tsat[nowg] != -1)  //群组有资源
            {
                U_Tsat[u] = G_Tsat[nowg] + 1 + deltaSAT;   //存为从1开始
                //群组接入单星单时隙，其前导资源数量作为群组内可用资源的用户数量，资源有限
                if (gtime[nowg] <= para_preamble)//还有前导资源
                {
                    U_Ttime[u] = G_Ttime[nowg] * 10;
                    U_prb[u] = gtime[nowg]; //第几份前导资源
                }
                else    //无前导资源了，延后5ms发起接入，并复用前导资源
                {
                    U_Ttime[u] = G_Ttime[nowg] * 10 + 5;
                    U_prb[u] = gtime[nowg] - para_preamble;
                }
            }
            else    //群组无资源
            {
                U_no.emplace_back(u);
            }     
        }
    }
    cout << "所在群组没有卫星时隙资源的用户数量：" << U_no.size() << endl;

    //输出txt
    string file_path1 = "U_Tsat1.txt";
    writeMatrixToFile(file_path1, U_Tsat);
    string file_path2 = "U_Ttime1.txt";
    writeMatrixToFile(file_path2, U_Ttime);
    string file_path3 = "U_prb1.txt";
    writeMatrixToFile(file_path3, U_prb);
}


int main() {
    BWOA();
    return 0;
}