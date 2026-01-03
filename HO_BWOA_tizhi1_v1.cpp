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
fobjReturn fobj0811v3(vector<double> Pos, int dim, int handoff_user_num, int entity_num, vector<Beam>& beamInfo, vector<baseStation>& baseInfo)
{
    int beam_num = beamInfo.size();
    int base_num = baseInfo.size();
    //初始化返回结构体
    fobjReturn fobj_return(handoff_user_num, entity_num, 0.0);
    vector<vector<double>> M_2D = make_vec2(handoff_user_num, entity_num, 0.0);
    // 将一维输入Pos重建为【用户×实体】二维矩阵
    for (int index = 0; index < Pos.size(); index++)
    {
        int user_idx = index / entity_num;    // 行：用户索引
        int entity_idx = index % entity_num;  // 列：实体索引
        M_2D[user_idx][entity_idx] = Pos[index];
    }

    int lambda1 = 1e4;  //成功率项权重
    int lambda2 = 1e2;    //负载均衡项权重

    // ====================== 约束1：每个用户只能接入1个实体 ======================
    //解释：惩罚项pen1，超出1颗星的用户按照超出的数量的平方累加惩罚值，得到所有用户的总惩罚值。
    //作用：惩罚违反约束的解，让优化过程倾向于满足约束。
    double pen1 = 0;    
    int muy1 = 1e5;
    for (int u = 0; u < handoff_user_num; u++)
    {
        double entity_access_sum = sum_vec1d(M_2D[u]);  // 该用户接入的实体数总和
        if (entity_access_sum > 1.0)  // 接入超过1个实体，触发惩罚
        {
            pen1 += pow(entity_access_sum - 1.0, 2);
        }
    }
    pen1 *= muy1;

    // ====================== 约束2：单个实体接入用户数不超限 ======================
    double pen2 = 0;
    double muy2 = 1e3;
    std::vector<int> entity_load(entity_num, 0);  // 每个实体的接入用户数
    std::vector<bool> is_load_over(entity_num, false);  // 标记实体是否超载
    //初始化entity_load为实体原有负载
    for (int e = 0; e < entity_num; e++)
    {
        if(e < beam_num)
        {
            entity_load[e] = beamInfo[e].current_load;
        }
        else
        {
            entity_load[e] = baseInfo[e - beam_num].current_load;
        }
    }
    for (int e = 0; e < entity_num; e++)
    {
        for (int u = 0; u < handoff_user_num; u++)
        {
            // 二元关系：M_2D[u][e]=1表示用户u接入实体e
            if (M_2D[u][e] == 1.0)
            {
                entity_load[e] += 1;
            }
        }
        // 检查是否超载，计算惩罚
        int Th_G = (e < beam_num) ? beamInfo[e].capacity : baseInfo[e - beam_num].capacity;
        if (entity_load[e] > Th_G)
        {
            is_load_over[e] = true;
            pen2 += pow(entity_load[e] - Th_G, 2);
        }
    }
    pen2 = pen2 * muy2;


    //step1：处理约束1:--若接入结果为多实体接入，保留单实体接入判决
    //过程：对于超出约束的用户，M_3D中仅保留最靠前的接入情况（置1），其他都置0。
    for (int u = 0; u < handoff_user_num; u++)
    {
        double access_sum = sum_vec1d(M_2D[u]);
        if (access_sum > 1.0)
        {
            // 找到用户最靠前接入的实体并保留
            for(int e = 0; e < entity_num; e++)
            {
                if(M_2D[u][e] == 1.0) { //找到了
                    M_2D[u] = vector<double>(entity_num, 0.0);
                    M_2D[u][e] = 1.0;
                    break;
                }
            }
        }
    }

    //step2：处理约束2-单星容量，记录在loadS（实际和上面约束2计算过程类似）
    std::fill(entity_load.begin(), entity_load.end(), 0);
    fill(is_load_over.begin(), is_load_over.end(), false);
    //初始化entity_load为实体原有负载
    for (int e = 0; e < entity_num; e++)
    {
        if(e < beam_num)
        {
            entity_load[e] = beamInfo[e].current_load;
        }
        else
        {
            entity_load[e] = baseInfo[e - beam_num].current_load;
        }
    }
    for (int e = 0; e < entity_num; e++)
    {
        for (int u = 0; u < handoff_user_num; u++)
        {
            // 二元关系：M_2D[u][e]=1表示用户u接入实体e
            if (M_2D[u][e] == 1.0)
            {
                entity_load[e] += 1;
            }
        }
        // 检查是否超载，超载直接记为门限
        int Th_G = (e < beam_num) ? beamInfo[e].capacity : baseInfo[e - beam_num].capacity;
        if (entity_load[e] > Th_G)
        {
            is_load_over[e] = true;
            entity_load[e] = Th_G;
        }
    }

    // 计算成功率项，即所有实体成功接入的用户总数
    int HO_success = sum_vec1d(entity_load);

    // 计算负载均衡项 越均衡值越小（修改成负载使用率的方差）
    double HO_load = 0.0;
    std::vector<double> load_utilization;  // 存储各实体的负载使用率（实际/上限）
    for (int e = 0; e < entity_num; e++)
    {
        int Th_G = (e < beam_num) ? beamInfo[e].capacity : baseInfo[e - beam_num].capacity;
        load_utilization.push_back(entity_load[e] / Th_G);
    }
    // 计算负载使用率的平均值
    double avg_utilization = 0.0;
    if (!load_utilization.empty())
    {
        avg_utilization = sum_vec1d(load_utilization) / load_utilization.size();
    }
    // 计算负载使用率的方差（方差越小，均衡性越好）
    for (double util : load_utilization)
    {
        HO_load += pow(util - avg_utilization, 2);
    }

    //load越小越好
    double objf = -lambda1 * HO_success + lambda2 * HO_load;

    fobj_return.M_2D = M_2D;    //新位置矩阵
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
void BWOA(double current_time, vector<User>& userInfo, vector<Satellite>& sateInfo, vector<Beam>& beamInfo, vector<baseStation>& baseInfo, double threshold)
{
    int user_num = userInfo.size();
    int sate_num = sateInfo.size();
    int base_num = baseInfo.size();
    int beam_num = beamInfo.size();

    //并行化计算所有用户到其当前连接对象的SINR
    vector<double> user_sinr_vec(user_num);
    #pragma omp parallel for
    for(int i = 0; i < user_num; i++)
    {
        user_sinr_vec[i] = calculateSINR(userInfo[i], sateInfo, baseInfo);
    }


    /*检测需要切换的用户*/
    vector<User> handoff_candidates;
    vector<int> handoff_user_indices; // 保存候选用户在原始userInfo中的索引
    for(int i = 0; i < user_num; i++)
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
    int handoff_user_num = handoff_candidates.size();
    cout << "检测到 " << handoff_user_num << " 个需要切换的用户" << endl;
    
    if (handoff_user_num == 0) {
        cout << "没有需要切换的用户，程序结束" << endl;
        return;
    }


    /*计算每个待切换用户的可见波束和基站，存到2个表中*/
    // 第一个表：vector<vector<int>> visible_beams_per_user
    // 第二维存储每个候选用户的可见波束索引
    vector<vector<int>> visible_beams_per_user(handoff_user_num);
    
    // 第二个表：vector<vector<int>> visible_bases_per_user
    // 第二维存储每个候选用户的可见基站索引
    vector<vector<int>> visible_bases_per_user(handoff_user_num);

    // 并行计算每个候选用户的可见实体
    #pragma omp parallel for
    for (int i = 0; i < handoff_user_num; i++) {
        const User& user = handoff_candidates[i];
        
        // 计算用户与所有波束的距离，判断可见波束
        for (int j = 0; j < beam_num; j++) {
            const Beam& beam = beamInfo[j];
            // 计算用户与波束的欧氏距离
            double dx = user.x - beam.x;
            double dy = user.y - beam.y;
            double distance = sqrt(dx*dx + dy*dy);
            
            // 如果距离小于波束覆盖半径，则波束可见
            if (distance <= beam.coverage_radius) {
                visible_beams_per_user[i].push_back(j);
            }
        }
        
        // 计算用户与所有基站的距离，判断可见基站
        for (int j = 0; j < base_num; j++) {
            const baseStation& base = baseInfo[j];
            // 计算用户与基站的欧氏距离
            double dx = user.x - base.x;
            double dy = user.y - base.y;
            double distance = sqrt(dx*dx + dy*dy);
            
            // 如果距离小于基站覆盖半径，则基站可见
            if (distance <= base.coverage_radius) {
                visible_bases_per_user[i].push_back(j);
            }
        }
    }


    /*输出每个候选用户的可见实体统计信息*/
    cout << "\n候选用户可见实体统计:" << endl;
    int total_visible_beams = 0;
    int total_visible_bases = 0;
    int users_with_no_visible = 0;
    
    for (int i = 0; i < handoff_user_num; i++) {
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
    
    cout << "  平均每个用户可见波束数: " << (double)total_visible_beams / handoff_user_num << endl;
    cout << "  平均每个用户可见基站数: " << (double)total_visible_bases / handoff_user_num << endl;
    cout << "  没有任何可见实体的用户数: " << users_with_no_visible << endl;


    /*需要先将当前的待切换用户从接入的波束或基站中移除，
    后续分配结果后再加到对应的实体上*/
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
    int entity_num = beam_num + base_num;   //可切换的总实体数
    int dim = handoff_user_num * entity_num;   //用户数*(波束数+基站数)*时隙数
    int SearchAgents_no = 80;   //使用80个搜索代理（个体）来寻找最优解。
    int Max_iter = 1;  //迭代次数（30，测试时暂改1）

    //创建一个dim列的全零向量，用于记录所有用户的切换最优解（所有搜索个体的leader）的位置。
    vector<double> Leader_pos(dim, 1);
    //记录Leader的当前和上一轮的的得分。
    double Leader_score = std::numeric_limits<double>::infinity();  //初始化为正无穷
    double Leader_score_pre = Leader_score;
    //停止算法运行的收敛容忍度
    double delta = 1e-6;
    int todoTol = 0;
    int Flag = 0;
    
    //（随机）初始化各个搜索个体的初始搜索位置。（之后用于计算每个卫星每个时隙的每个用户的最优解）
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


    /*初始化生成一个较可行解，考虑了用户的可见实体和实体的剩余容量来分配切换资源。
        思路：先获取所有切换实体的剩余容量。
        接着遍历所有待切换用户（打乱顺序遍历），在其可见实体中找到
        仍有剩余容量的，选择其中容量最大的一个，分配切换资源。*/
    
    vector<int> entity_remaining_capacity(entity_num, 0);  // 每个实体的剩余容量

    // 初始化波束的剩余容量（前beam_num个是波束）
    for (int i = 0; i < beam_num; i++) {
        entity_remaining_capacity[i] = beamInfo[i].capacity - beamInfo[i].current_load;
        // 剩余容量不能为负数
        if (entity_remaining_capacity[i] < 0) {
            throw std::runtime_error("波束的剩余容量不能为负数");
        }
    }
    
    // 初始化基站的剩余容量（后base_num个是基站）
    for (int i = 0; i < base_num; i++) {
        entity_remaining_capacity[beam_num + i] = baseInfo[i].capacity - baseInfo[i].current_load;
        // 剩余容量不能为负数
        if (entity_remaining_capacity[beam_num + i] < 0) {
            throw std::runtime_error("基站的剩余容量不能为负数");
        }
    }

    // 创建二维分配矩阵
    auto RRR = make_vec2(handoff_user_num, entity_num, 0.0);
    // 遍历所有的待切换用户
    // 先给出待切换用户的索引(0~handoff_user_num-1)的打乱顺序
    vector<int> shuffle_user_index = randperm(handoff_user_num, handoff_user_num);
    for (int k = 0; k < handoff_user_num; k++) { 
        int i = shuffle_user_index[k];  // 待切换用户的打乱索引
        // 遍历可见的波束，获取容量最大的一个的索引及其容量
        int max_beam_index = -1;
        int max_beam_capacity = -1;
        for (int j = 0; j < visible_beams_per_user[i].size(); j++) { 
            int beam_index = visible_beams_per_user[i][j];
            if (entity_remaining_capacity[beam_index] > max_beam_capacity) { 
                max_beam_index = beam_index;
                max_beam_capacity = entity_remaining_capacity[beam_index];
            }
        }

        //遍历可见的基站，获取容量最大的一个的索引（注意转成绝对索引）及其容量
        int max_base_index = -1;
        int max_base_capacity = -1;
        for (int j = 0; j < visible_bases_per_user[i].size(); j++) { 
            int base_index = visible_bases_per_user[i][j];
            if (entity_remaining_capacity[beam_num + base_index] > max_base_capacity) { 
                max_base_index = beam_num + base_index;
                max_base_capacity = entity_remaining_capacity[beam_num + base_index];
            }
        }

        //比较两者可见实体的容量，在两者都存在的情况下选择其中容量较大的实体分配资源，相同则随机选择一个
        //并更新剩余容量
        if (max_beam_index != -1 && max_base_index != -1) { 
            // 先判断两个实体的剩余容量是否都大于0（容量为0则不分配）
            if (entity_remaining_capacity[max_beam_index] > 0 && entity_remaining_capacity[max_base_index] > 0) {
                if (max_beam_capacity > max_base_capacity) { 
                    RRR[i][max_beam_index] = 1;
                    entity_remaining_capacity[max_beam_index] -= 1; // 分配后更新剩余容量
                }
                else if (max_beam_capacity < max_base_capacity) { 
                    RRR[i][max_base_index] = 1;
                    entity_remaining_capacity[max_base_index] -= 1; // 分配后更新剩余容量
                }
                else { // 容量相同则随机选择一个，且更新对应剩余容量
                    if (dis_pos(gen) <= 0.5) {
                        RRR[i][max_beam_index] = 1;
                        entity_remaining_capacity[max_beam_index] -= 1;
                    } else {
                        RRR[i][max_base_index] = 1;
                        entity_remaining_capacity[max_base_index] -= 1;
                    }
                }
            }
            // 仅beam容量>0，base容量=0：分配给beam
            else if (entity_remaining_capacity[max_beam_index] > 0) {
                RRR[i][max_beam_index] = 1;
                entity_remaining_capacity[max_beam_index] -= 1;
            }
            // 仅base容量>0，beam容量=0：分配给base
            else if (entity_remaining_capacity[max_base_index] > 0) {
                RRR[i][max_base_index] = 1;
                entity_remaining_capacity[max_base_index] -= 1;
            }
            // 两者容量都为0：不分配
        }
        else if (max_beam_index != -1) { 
            // 仅beam实体存在，且剩余容量>0时才分配
            if (entity_remaining_capacity[max_beam_index] > 0) {
                RRR[i][max_beam_index] = 1;
                entity_remaining_capacity[max_beam_index] -= 1;
            }
        }
        else if (max_base_index != -1) { 
            // 仅base实体存在，且剩余容量>0时才分配
            if (entity_remaining_capacity[max_base_index] > 0) {
                RRR[i][max_base_index] = 1;
                entity_remaining_capacity[max_base_index] -= 1;
            }
        }
    }

    // 只生成一半的搜索个体作为较可行解，另一半保持随机，将RRR抻平放入Positions中。
    int half_agents = SearchAgents_no / 2;
    for (int i = 0; i < half_agents; i++) 
    { 
        int index = 0;
        for(int j = 0; j < handoff_user_num; j++) 
        {
            for(int k = 0; k < entity_num; k++)
            {
                Positions[i][index++] = RRR[j][k];
            }
        }
    }
    vector<double> Convergence_curve(Max_iter, 0);
    auto D_X = make_vec2<double>(SearchAgents_no, dim, 0.0);    //二维
    
    //创建一个二维vector存储切换结果。
    vector<vector<double>> Switching_result = make_vec2<double>(handoff_user_num, entity_num, 0.0);

    //循环计数器
    int iter = 0;

    //主循环
    vector<double> old_pos = Positions[0];
    while (iter < Max_iter && Flag <= 3)
    {
        for (int i = 0; i < Positions.size(); i++)  //遍历每个SearchAgent
        {
            // 对所有参与搜索的个体（search agent），逐一计其在目标函数（object function），也就是适应度函数中的取值。（用于评估搜索个体的优劣）
            fobjReturn fobj_return = fobj0811v3(Positions[i], dim, handoff_user_num, entity_num, beamInfo, baseInfo);
            vector<vector<double>> M_2D = fobj_return.M_2D;
            double HO_Sucnow = fobj_return.HO_Sucnow;
            double fitness = fobj_return.fitness;

            //更新Leader
            //cout << fitness << endl;
            if (fitness < Leader_score)   //fitness越小越好，故最小化得分以寻找最优解，解决优化问题。
            {
                Leader_score = fitness; //更新最佳得分
                //最佳位置也要更新成新的M_3D抻平后的结果
                int pos_index = 0;
                for (int u = 0; u < handoff_user_num; u++) 
                { 
                    for(int e = 0; e < entity_num; e++)
                    {
                        Leader_pos[pos_index++] = M_2D[u][e];
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

            //若Positions[i]中1的数量多于handoff_user_num，则随机删去多余的1。
            double pos1sum = sum_vec1d(Positions[i]);
            if(pos1sum > handoff_user_num)
            {
                // 先找到Positions[i]中所有值为1的位置（索引），随机取其中的handoff_user_num个置1。
                vector<int> pos1_index;
                for (int k = 0; k < dim; k++)   //Positions[0].size()
                {
                    if (Positions[i][k] == 1)
                    {
                        pos1_index.push_back(k);
                    }
                }
                // 在0~index.size()的数中，随机取其中handoff_user_num个。
                vector<int> pos1_index_rand = randperm(pos1_index.size(), handoff_user_num);
                // 先将原来的Positions[i]中的数置0。
                vector<double> pos_zero(dim, 0.0);
                Positions[i] = pos_zero;
                // 再遍历pos1_index_rand，将pos1_index_rand中的索引对应的值置1。
                vector<double> PPP(dim, 0.0);
                for (int k = 0; k < pos1_index_rand.size(); k++) {
                    PPP[pos1_index[pos1_index_rand[k]]] = 1;
                }
                Positions[i] = PPP;
            }
        }
        //记录每次迭代后的头鲸的适应度函数的得分（越小越好）
        Convergence_curve[iter] = Leader_score;
        //迭代次数增加1
        iter = iter + 1;
        
        fobjReturn fobj_return = fobj0811v3(Leader_pos, dim, handoff_user_num, entity_num, beamInfo, baseInfo);
        Switching_result = fobj_return.M_2D;    //存储每次迭代后的切换结果
        printf("迭代次数%d,Score=%f,HO_success=%d\n", iter, Leader_score, (int)fobj_return.HO_Sucnow);

        if (todoTol == 1 && abs(Leader_score - Leader_score_pre) < delta)
        {
            Flag += 1;
            Convergence_curve = vector<double>(Convergence_curve.begin(), Convergence_curve.begin() + iter);
        }
        Leader_score_pre = Leader_score;
    }

    /*将最优的分配结果存起来。
        之前存储过所有待切换用户的原有索引handoff_user_indices，
        现在需要按照切换结果，赋值到原有的userInfo中作为切换输出，
        并且增加对应实体的current_load。*/
    int ceshi = 1;
    for (int i = 0; i < handoff_user_num; i++)
    {
        //计算切换用户的原索引
        int orginal_index = handoff_user_indices[i];
        vector<double> temp_vec = Switching_result[i];
        auto it = find(temp_vec.begin(), temp_vec.end(), 1);
        if (it == temp_vec.end()) {
            //没有找到1，则说明该用户没有被分配。
            //该用户变为未接入状态
            userInfo[orginal_index].is_connected = false;
            userInfo[orginal_index].access_entity = "";
            userInfo[orginal_index].access_entity_id = -1;
            continue;
        }

        //计算对应切换实体的索引
        int entity_index = it - temp_vec.begin();
        userInfo[orginal_index].is_connected = true;
        if (entity_index < beam_num)
        {
            userInfo[orginal_index].access_entity = "beam";
            userInfo[orginal_index].access_entity_id = beamInfo[entity_index].id;
            beamInfo[entity_index].current_load = beamInfo[entity_index].current_load + 1;
        }
        else
        {
            userInfo[orginal_index].access_entity = "baseStation";
            userInfo[orginal_index].access_entity_id = baseInfo[entity_index - beam_num].id;
            baseInfo[entity_index - beam_num].current_load = baseInfo[entity_index - beam_num].current_load + 1;
        }
    }
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