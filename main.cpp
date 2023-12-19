// 包含头文件
#include <iostream>
#include <chrono>

// 自定义头文件
#include <fiberbundle.h>
#include <computedispersion.h>
#include <measure_utils.h>

using namespace std;
using namespace magic_sheep;

// 主函数
int main() {
    string input_filename;
    cout << "输入文件名：" << endl;
    cin >> input_filename;
    //filename = "1T_fiber.vtk";
    // 读取 纤维束 数据  .vtk
    fiberbundle myBundle;
    auto readTime = measureTimeMillis([&]() { myBundle.ReadFibers(input_filename); });
    std::cout << "读取文件耗时: " << readTime << " ms" << std::endl;
//    myBundle.ReadFibers(input_filename);  // 括号里是vtk文件名
    //filename = "fiber_1T_output.vtk";
    // 数据处理
    readTime = measureTimeMillis([&]() { computeDispersion(myBundle, 3, 3, "", 1, 1); });
    std::cout << "计算耗时: " << readTime << " ms" << std::endl;

//    computeDispersion(myBundle, 3, 3, "", 1, 1);
    /*int computeDispersion(fiberbundle& bundle,
    double scale,  半径
    unsigned int numberOfSamplingDirections, 采样方向数
        const std::string& outputFilename,
    unsigned int tractSubSampling = 1,  纤维束子采样
    unsigned int fiberPointSubSampling = 1);   点子采样  */

    // 输出处理结果  .vtk  文件
    string output_filename;
    cout << "输出文件名 vtk：" << endl;
    cin >> output_filename;
    myBundle.WriteFibers(output_filename, false, true);  //writeAscii,writeUnCompressed: true or false
    //myBundle.Print();
    myBundle.getTd();     // 输出 TD 值 ，保存为 TXT 文档
    myBundle.getFiberTd();   // 输出每根纤维对应的 TD 值，保存为 TXT 文档

    return 0;
}