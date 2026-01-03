#include<iostream>
#include<vector>
using namespace std;

/*定义常数项*/
const double pi = 3.14159;

/*定义二维和三维vector的模板*/
//模板别名（支持任意类型 T）
template <typename T>
using Vec2 = vector<vector<T>>;

template <typename T>
using Vec3 = vector<Vec2<T>>;

// 模板函数：创建二维 vector（任意类型 T）
template <typename T>
Vec2<T> make_vec2(int rows, int cols, const T& val = T()) {
    return Vec2<T>(rows, vector<T>(cols, val));
}

// 模板函数：创建三维 vector（任意类型 T）
template <typename T>
Vec3<T> make_vec3(int dim1, int dim2, int dim3, const T& val = T()) {
    return Vec3<T>(dim1, make_vec2<T>(dim2, dim3, val));
}

/*fobj函数（适应度函数fitness objection function）的返回结构体*/
struct fobjReturn {
    vector<vector<vector<double>>> M_3D;
    double HO_Sucnow;
    double fitness;

    //默认构造函数
    fobjReturn() : HO_Sucnow(0.0), fitness(0.0) {}

    //带维度的构造函数
    fobjReturn(int dim1, int dim2, int dim3, double init_val = 0.0)
        :M_3D(dim1, vector<vector<double>>(dim2, vector<double>(dim3, init_val))),
        HO_Sucnow(0.0), fitness(0.0) {}
};

// 定义各个（卫星）体制的枚举
enum TIZHI {
    NR,
    DVB,
    GMR
};
#define NUM_TIZHI 3

// 用户信息结构体
struct User {
    // 其他信息比如位置信息
    double x;
    double y;
    bool is_connected;      //用户是否接入
    string access_entity;  //当前接入对象，包含"baseStation"或"beam"
    int access_entity_id;   // 当前接入对象的ID，如果未连接则为-1
};

// 卫星信息结构体
struct Satellite {
    // 其他信息比如位置信息
    int capacity;           //卫星的总容量
    int current_load;       //卫星的当前使用容量
    int preamble_num;       //卫星的前导码数量（公共？专用？）
};

// 波束信息结构体
struct Beam {
    int id;
    // 其他信息比如位置信息
    double x;
    double y;
    double coverage_radius; //波束的覆盖半径
    
    int capacity;           //波束的总容量
    int current_load;       //波束的当前使用容量
    int preamble_num;       //波束的前导码数量（公共？专用？）
};

struct baseStation {
    int id;
    // 其他信息比如位置信息
    double x;
    double y;
    double coverage_radius;    //基站的覆盖半径
    int capacity;           //基站的总容量
    int current_load;       //基站的当前使用容量
    int preamble_num;       //基站的前导码数量（公共？专用？）
};