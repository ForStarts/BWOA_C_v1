#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <random>
#include <limits>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <omp.h>
#include "Variant.hpp"
#include "nlohmann/json.hpp"
using json = nlohmann::json;
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
fobjReturn fobj0811v3(vector<double> Pos, int dim, vector<int> para_GRemain, vector<int> para_Galloc, int para_Gnum, int para_Snum, int para_period, int Th_G)
{
    //初始化返回结构体
    fobjReturn fobj_return(para_Gnum, para_Snum, para_period, 0.0);
    vector<vector<vector<double>>> M_3D = make_vec3(para_Gnum, para_Snum, para_period, 0.0);
    //将一维数组重建为三维数组（主意这里matlab的reshape是列主序！！！）
    for (int index = 0; index < Pos.size(); index++)
    {
        int i = index / (para_Snum * para_period);
        int j = (index % (para_Snum * para_period)) / para_period;
        int k = (index % (para_Snum * para_period)) % para_period;
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
    for (int g = 0; g < para_Gnum; g++)
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
    vector<int> loadS(para_Snum, 0);    //记录sat负载，超载则直接记当前负载为门限值
    vector<bool> isLoadOver(para_Snum, 0);  //标记sat是否超负载
    for (int s = 0; s < para_Snum; s++)
    {
        for (int g = 0; g < para_Gnum; g++)
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
    for (int s = 0; s < para_Snum; s++)
    {
        for (int t = 0; t < para_period; t++)
        {
            //计算该星该时隙接入了多少群组
            double sum_group_num = 0;
            for (int g = 0; g < para_Gnum; g++)
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
    for (int gID = 0; gID < para_Gnum; gID++)
    {
        if (sum_vec1d(sum_vec2d(M_3D[gID])) > 1)   //单群组存在多星/多时隙接入情况
        {
            //先遍历时间，则可以快速找到时机最靠前判决
            bool hasFind = 0;
            for (int t = 0; t < para_period; t++)
            {
                for (int s = 0; s < para_Snum; s++)
                {
                    if (M_3D[gID][s][t] == 1){  //找到了
                        //保留时机最靠前判决
                        M_3D[gID] = make_vec2<double>(para_Snum, para_period, 0.0);
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
    for (int s = 0; s < para_Snum; s++)
    {
        for (int t = 0; t < para_period; t++)
        {
            //计算该星该时隙接入了多少群组
            double sum_group_num = 0;
            int max_user = 0;
            double max_user_group_idx = 0;  //存储用户数最多的群组的ID
            for (int g = 0; g < para_Gnum; g++)
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
    for (int s = 0; s < para_Snum; s++)
    {
        for (int g = 0; g < para_Gnum; g++)
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
    double HO_load_avg = HO_success / para_Snum;
    double HO_load = 0; //负载均衡项
    for (int s = 0; s < para_Snum; s++)
    {
        HO_load += pow((loadS[s] - HO_load_avg), 2);
    }

    //【新增】计算剩余时间项，剩余时间越短的越早切换 越小越好
    //过程：遍历群组，找到所有接入位置（ss星tt时隙），由于之前处理了约束3，使得
    //单星单时隙仅有一个群组，所以只会找到一个接入位置。累积（该时间*10-群组剩余时间）。
    double HO_T = 0;
    for (int g = 0; g < para_Gnum; g++)
    {
        //寻找该群组的接入位置
        bool hasFind = 0;
        for (int ss = 0; ss < para_Snum; ss++)
        {
            for (int tt = 0; tt < para_period; tt++)
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


/*其他模块的，切换需要调用的函数*/
double calculateSINR(User user, vector<Satellite> sateInfo, vector<baseStation> baseInfo)
{
    // 根据用户当前位置、连接的实体、干扰模型计算SINR
    // 考虑路径损耗、阴影衰落、多径效应等
    return 0.0;
}


void init_para(vector<User>& userInfo, vector<Satellite>& sateInfo, 
               vector<Beam>& beamInfo, vector<baseStation>& baseInfo)
{
    // 清空现有的数据
    userInfo.clear();
    sateInfo.clear();
    beamInfo.clear();
    baseInfo.clear();
    
    try {
        // 读取JSON文件
        std::ifstream json_file("entity_data.json");
        if (!json_file.is_open()) {
            std::cerr << "错误：无法打开 entity_data.json 文件" << std::endl;
            std::cerr << "当前工作目录: " << std::filesystem::current_path() << std::endl;
            return;
        }
        
        // 解析JSON
        json json_data;
        json_file >> json_data;
        json_file.close();
        
        // 初始化基站信息
        if (json_data.contains("baseStations") && json_data["baseStations"].is_array()) {
            for (const auto& base_json : json_data["baseStations"]) {
                baseStation base;
                
                // 读取ID
                if (base_json.contains("id")) {
                    // JSON中的id是数字，直接读取
                    base.id = base_json["id"];
                    // 不需要额外转换
                }
                
                // 读取位置信息
                if (base_json.contains("x") && base_json["x"].is_number()) {
                    base.x = base_json["x"];
                }
                if (base_json.contains("y") && base_json["y"].is_number()) {
                    base.y = base_json["y"];
                }
                
                // 读取覆盖半径
                if (base_json.contains("coverage_radius") && base_json["coverage_radius"].is_number()) {
                    base.coverage_radius = base_json["coverage_radius"];
                }
                
                // 读取容量信息
                if (base_json.contains("capacity") && base_json["capacity"].is_number_integer()) {
                    base.capacity = base_json["capacity"];
                }
                if (base_json.contains("current_load") && base_json["current_load"].is_number_integer()) {
                    base.current_load = base_json["current_load"];
                }
                
                // 读取前导码数量
                if (base_json.contains("preamble_num") && base_json["preamble_num"].is_number_integer()) {
                    base.preamble_num = base_json["preamble_num"];
                }
                
                baseInfo.push_back(base);
            }
            std::cout << "基站初始化完成: " << baseInfo.size() << " 个基站" << std::endl;
        }
        
        // 初始化波束信息
        if (json_data.contains("beams") && json_data["beams"].is_array()) {
            for (const auto& beam_json : json_data["beams"]) {
                Beam beam;
                
                // 读取ID
                if (beam_json.contains("id")) {
                    // JSON中的id是数字，直接读取
                    beam.id = beam_json["id"];
                }
                
                // 读取位置信息
                if (beam_json.contains("x") && beam_json["x"].is_number()) {
                    beam.x = beam_json["x"];
                }
                if (beam_json.contains("y") && beam_json["y"].is_number()) {
                    beam.y = beam_json["y"];
                }
                
                // 读取覆盖半径
                if (beam_json.contains("coverage_radius") && beam_json["coverage_radius"].is_number()) {
                    beam.coverage_radius = beam_json["coverage_radius"];
                }
                
                // 读取容量信息
                if (beam_json.contains("capacity") && beam_json["capacity"].is_number_integer()) {
                    beam.capacity = beam_json["capacity"];
                }
                if (beam_json.contains("current_load") && beam_json["current_load"].is_number_integer()) {
                    beam.current_load = beam_json["current_load"];
                }
                
                // 读取前导码数量
                if (beam_json.contains("preamble_num") && beam_json["preamble_num"].is_number_integer()) {
                    beam.preamble_num = beam_json["preamble_num"];
                }
                
                beamInfo.push_back(beam);
            }
            std::cout << "波束初始化完成: " << beamInfo.size() << " 个波束" << std::endl;
        }
        
        // 初始化用户信息
        if (json_data.contains("users") && json_data["users"].is_array()) {
            for (const auto& user_json : json_data["users"]) {
                User user;
                
                // 读取位置信息
                if (user_json.contains("x") && user_json["x"].is_number()) {
                    user.x = user_json["x"];
                }
                if (user_json.contains("y") && user_json["y"].is_number()) {
                    user.y = user_json["y"];
                }
                
                // 读取连接状态
                if (user_json.contains("is_connected") && user_json["is_connected"].is_boolean()) {
                    user.is_connected = user_json["is_connected"];
                }
                
                // 读取接入实体类型
                if (user_json.contains("access_entity") && user_json["access_entity"].is_string()) {
                    user.access_entity = user_json["access_entity"].get<std::string>();
                }
                
                // 注意：根据您的要求，User结构体需要添加access_entity_id字段
                // 您需要在Variant.hpp中的User结构体中添加这个字段
                // 如果您还没有添加，请先在Variant.hpp中添加：
                // struct User {
                //     double x;
                //     double y;
                //     bool is_connected;
                //     std::string access_entity;
                //     int access_entity_id;  // 添加这一行
                // };
                
                // 读取接入实体ID
                if (user_json.contains("access_entity_id") && user_json["access_entity_id"].is_number_integer()) {
                    user.access_entity_id = user_json["access_entity_id"];
                }
                
                userInfo.push_back(user);
            }
            std::cout << "用户初始化完成: " << userInfo.size() << " 个用户" << std::endl;
        }
        
        // 初始化卫星信息（暂时为空，但保留结构）
        // 根据JSON的metadata获取卫星数量，或者使用默认值
        int satellite_count = 6; // 默认6颗卫星
        
        if (json_data.contains("metadata")) {
            const auto& metadata = json_data["metadata"];
            // 如果有卫星数量的元数据，可以在这里读取
            // 例如：if (metadata.contains("satellite_count")) { satellite_count = metadata["satellite_count"]; }
        }
        
        for (int i = 0; i < satellite_count; i++) {
            Satellite sate;
            sateInfo.push_back(sate);
        }
        std::cout << "卫星初始化完成: " << sateInfo.size() << " 颗卫星" << std::endl;
        
        // 输出连接统计（从metadata中读取）
        if (json_data.contains("metadata")) {
            const auto& metadata = json_data["metadata"];
            if (metadata.contains("connection_stats")) {
                const auto& stats = metadata["connection_stats"];
                std::cout << "\n连接统计信息:" << std::endl;
                if (stats.contains("connected_to_base")) {
                    std::cout << "  连接到基站: " << stats["connected_to_base"] << std::endl;
                }
                if (stats.contains("connected_to_beam")) {
                    std::cout << "  连接到波束: " << stats["connected_to_beam"] << std::endl;
                }
                if (stats.contains("not_connected")) {
                    std::cout << "  未连接: " << stats["not_connected"] << std::endl;
                }
                if (stats.contains("total_users")) {
                    std::cout << "  总用户数: " << stats["total_users"] << std::endl;
                }
            }
        }
        
        std::cout << "\n数据初始化全部完成!" << std::endl;
        
    } catch (const json::exception& e) {
        std::cerr << "JSON解析错误: " << e.what() << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "初始化错误: " << e.what() << std::endl;
    }
}




/*二进制鲸鱼优化算法
输入：userInfo 待切换的用户信息
输入：sateInfo 卫星信息
输入：baseInfo 基站信息*/
void BWOA(double current_time, vector<User> userInfo, vector<Satellite> sateInfo, vector<Beam> beamInfo, vector<baseStation> baseInfo, double threshold)
{
    int usernum = userInfo.size();
    int satenum = sateInfo.size();
    int basenum = baseInfo.size();  //基站数

    //并行化计算所有用户到其当前连接对象的SINR
    vector<double> user_sinr_vec(usernum);
    #pragma omp parallel for
    for(int i = 0; i < usernum; i++)
    {
        user_sinr_vec[i] = calculateSINR(userInfo[i], sateInfo, baseInfo);
    }


    /*检测需要切换的用户*/
    vector<User> handoff_candidates;
    vector<int> handoff_user_indices; // 保存候选用户在原始userInfo中的索引
    for(int i = 0; i < usernum; i++)
    {
        //仅考虑已连接的用户
        if(userInfo[i].is_connected)
        {
            double current_sinr = calculateSINR(userInfo[i], sateInfo, baseInfo);

            if(current_sinr < threshold)
            { 
                //筛选后的候选用户加入列表
                handoff_candidates.emplace_back(userInfo[i]);
                handoff_user_indices.push_back(i);
            }
        }
    }
    int num_candidates = handoff_candidates.size();
    cout << "检测到 " << num_candidates << " 个需要切换的用户" << endl;
    
    if (num_candidates == 0) {
        cout << "没有需要切换的用户，程序结束" << endl;
        return;
    }


    /*计算每个待切换用户的可见波束和基站，存到2个表中*/
    // 第一个表：vector<vector<Beam>> visible_beams_per_user
    // 第二维存储每个候选用户的可见波束
    vector<vector<Beam>> visible_beams_per_user(num_candidates);
    
    // 第二个表：vector<vector<baseStation>> visible_bases_per_user
    // 第二维存储每个候选用户的可见基站
    vector<vector<baseStation>> visible_bases_per_user(num_candidates);

    // 并行计算每个候选用户的可见实体
    #pragma omp parallel for
    for (int i = 0; i < num_candidates; i++) {
        const User& user = handoff_candidates[i];
        
        // 计算用户与所有波束的距离，判断可见波束
        for (const auto beam : beamInfo) {
            // 计算用户与波束的欧氏距离
            double dx = user.x - beam.x;
            double dy = user.y - beam.y;
            double distance = sqrt(dx*dx + dy*dy);
            
            // 如果距离小于波束覆盖半径，则波束可见
            if (distance <= beam.coverage_radius) {
                visible_beams_per_user[i].push_back(beam);
            }
        }
        
        // 计算用户与所有基站的距离，判断可见基站
        for (const auto base : baseInfo) {
            // 计算用户与基站的欧氏距离
            double dx = user.x - base.x;
            double dy = user.y - base.y;
            double distance = sqrt(dx*dx + dy*dy);
            
            // 如果距离小于基站覆盖半径，则基站可见
            if (distance <= base.coverage_radius) {
                visible_bases_per_user[i].push_back(base);
            }
        }
    }


    /*输出每个候选用户的可见实体统计信息*/
    cout << "\n候选用户可见实体统计:" << endl;
    int total_visible_beams = 0;
    int total_visible_bases = 0;
    int users_with_no_visible = 0;
    
    for (int i = 0; i < num_candidates; i++) {
        int num_beams = visible_beams_per_user[i].size();
        int num_bases = visible_bases_per_user[i].size();
        
        total_visible_beams += num_beams;
        total_visible_bases += num_bases;
        
        if (num_beams == 0 && num_bases == 0) {
            users_with_no_visible++;
        }
        
        // 可选：输出每个用户的详细情况
        
        cout << "  候选用户 " << i << " (原始索引 " << handoff_user_indices[i] << "): "
             << "可见波束=" << num_beams 
             << ", 可见基站=" << num_bases << endl;
        
    }
    
    cout << "  平均每个用户可见波束数: " << (double)total_visible_beams / num_candidates << endl;
    cout << "  平均每个用户可见基站数: " << (double)total_visible_bases / num_candidates << endl;
    cout << "  没有任何可见实体的用户数: " << users_with_no_visible << endl;


    /*需要先将当前的待切换用户从接入的波束或基站中移除*/
    for (const auto& user : handoff_candidates) {
        if (user.access_entity == "beam") {
            for (auto& beam : beamInfo) {
                if (beam.id == user.access_entity_id) {
                    beam.current_load--;
                }
            }
        }
        else if (user.access_entity == "base") {
            for (auto& base : baseInfo) {
                if (base.id == user.access_entity_id) {
                    base.current_load--;
                }
            }
        }
    }
    

    /*执行切换操作，0-1优化求解切换问题*/
    
}


int main() {
    vector<User> userInfo;
    vector<Satellite> sateInfo;
    vector<Beam> beamInfo;
    vector<baseStation> baseInfo;

    //测试：初始化数据
    init_para(userInfo, sateInfo, beamInfo, baseInfo);

    // 验证数据是否正确加载
    if (!userInfo.empty()) {
        std::cout << "\n前5个用户信息:" << std::endl;
        for (int i = 0; i < std::min(5, (int)userInfo.size()); i++) {
            const auto& user = userInfo[i];
            std::cout << "  用户" << i << ": x=" << user.x << ", y=" << user.y 
                      << ", 连接状态=" << (user.is_connected ? "已连接" : "未连接")
                      << ", 接入实体=" << user.access_entity
                      << ", 实体ID=" << user.access_entity_id << std::endl;
        }
    }
    
    if (!baseInfo.empty()) {
        std::cout << "\n前3个基站信息:" << std::endl;
        for (int i = 0; i < std::min(3, (int)baseInfo.size()); i++) {
            const auto& base = baseInfo[i];
            std::cout << "  基站" << i << ": x=" << base.x << ", y=" << base.y 
                      << ", 容量=" << base.capacity << ", 当前负载=" << base.current_load 
                      << ", 负载率=" << (base.capacity > 0 ? (double)base.current_load/base.capacity*100 : 0) << "%" << std::endl;
        }
    }
    
    if (!beamInfo.empty()) {
        std::cout << "\n前3个波束信息:" << std::endl;
        for (int i = 0; i < std::min(3, (int)beamInfo.size()); i++) {
            const auto& beam = beamInfo[i];
            std::cout << "  波束" << i << ": x=" << beam.x << ", y=" << beam.y 
                      << ", 容量=" << beam.capacity << ", 当前负载=" << beam.current_load 
                      << ", 负载率=" << (beam.capacity > 0 ? (double)beam.current_load/beam.capacity*100 : 0) << "%" << std::endl;
        }
    }

    double threshold = 0.5; //SINR切换阈值
    double current_time = 0.0;

    //注意：这里仅传入需要切换的用户，不需要切换的可以不传
    //注意：切换仅可切换覆盖该用户的基站或卫星。
    BWOA(current_time, userInfo, sateInfo, beamInfo, baseInfo, threshold);
    return 0;
}