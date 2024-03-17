﻿// 包含头文件
#include <iostream>
#include <chrono>

#include <mex.hpp>
#include <mexAdapter.hpp>

// 自定义头文件
#include <fiberbundle.h>
#include <computedispersion.h>
#include <filesystem>
#include <measure_utils.h>

using namespace std;

auto measureTimeMillis(const function<void()>& func) -> decltype(std::chrono::milliseconds().count());

class MexFunction : public matlab::mex::Function {
    matlab::data::ArrayFactory factory;
    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();

    vector<float> td;
    map<int, vector<float>> fiberTd;
public:
    /**
     * @example
     * @param outputs
     * @param inputs input_filename, output_filename, scale, numberOfSamplingDirections, tractSubSampling, fiberPointSubSampling
     */
    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) override {
        using namespace matlab::data;
        if (inputs.size() == 0 && outputs.size() == 2) {
            // 保存输出参数
            TypedArray<float> output_1{ factory.createArray<float>({td.size(), 1}, td.begin().base(), td.end().base()) };
            outputs[0] = std::move(output_1);
            outputs[1] = getMexArray(std::move(this->fiberTd));
            return ;
        }
        checkInput(inputs);
        if (outputs.size() != 2) {
            // 在matlab上打印错误日志
            matlabPtr->feval(u"error", 0, std::vector<Array>({ factory.createScalar(
                    "Error: the number of output should be less to 2. \n"
                    "rerun result this function with no input and two output. \n"
                    ) }));
        }
        string inputFilename;
        // 判断inputFilename的类型
        if (inputs[0].getType() == ArrayType::CHAR) {
            // 如果是char类型,则转换为string
            CharArray input_filename = inputs[0];
            inputFilename = input_filename.toAscii();
        } else if (inputs[0].getType() == ArrayType::MATLAB_STRING) {
            StringArray input_filename = inputs[0];
            inputFilename = matlab::engine::convertUTF16StringToUTF8String(input_filename[0]);
        } else {
            // 在matlab上打印错误日志
            matlabPtr->feval(u"error", 0, std::vector<Array>({ factory.createScalar(
                    "Error: the input file name should be char or string."
            ) }));
        }
        cout << "inputFilename: " << inputFilename << endl;
        // 判断对应路径的文件是否存在
        ifstream f(inputFilename);
        if (!f.good()) {
            // 在matlab上打印错误日志
            matlabPtr->feval(u"error", 0, std::vector<Array>({ factory.createScalar(
                    "Error: the input file is not exist. \n"
            ) }));
        }
        string outputFilename;
        if (inputs[1].getType() == ArrayType::CHAR)
        {
            // 如果是char类型,则转换为string
            CharArray output_filename = inputs[1];
            outputFilename = output_filename.toAscii();
        }
        else if (inputs[1].getType() == ArrayType::MATLAB_STRING) {
            StringArray output_filename = inputs[1];
            outputFilename = matlab::engine::convertUTF16StringToUTF8String(output_filename[0]);
        } else {
            // 在matlab上打印错误日志
            matlabPtr->feval(u"error", 0, std::vector<Array>({ factory.createScalar(
                    "Error: the output file name should be char or string."
            ) }));
        }
        // 判断父目录是否存在,不存在则创建
        string parentPath = outputFilename.substr(0, outputFilename.find_last_of("/\\"));
        if (!filesystem::exists(parentPath)) {
            filesystem::create_directories(parentPath);
        }
        cout << "outputFilename: " << outputFilename << endl;
        double scale = inputs[2][0];
        cout << "scale: " << scale << endl;
        uint numberOfSamplingDirections = inputs[3][0];
        cout << "numberOfSamplingDirections: " << numberOfSamplingDirections << endl;
        uint tractSubSampling = inputs[4][0];
        cout << "tractSubSampling: " << tractSubSampling << endl;
        uint fiberPointSubSampling = inputs[5][0];
        cout << "fiberPointSubSampling: " << fiberPointSubSampling << endl;

        fiberbundle myBundle;
        myBundle.ReadFibers(inputFilename);
        cout << "main compute cost: " << magic_sheep::measureTimeMillis([&]() {
            computeDispersion(myBundle, scale, numberOfSamplingDirections, outputFilename, tractSubSampling, fiberPointSubSampling);
        }) << "ms\n";
        // 输出计算完成

        myBundle.WriteFibers(outputFilename, false, true);

        if (outputs.size() != 2) {
            // 保存输出参数
            this->td = myBundle.getTd();
            this->fiberTd = myBundle.getFiberTd();

            matlabPtr->feval(u"error", 0, std::vector<Array>({ factory.createScalar(
                    "Error: the number of output should be equal to 2. \n"
            ) }));
        }
        // 输出参数第一个保存一个列向量
        auto output_1 = myBundle.getTd();
        matlab::data::TypedArray<float> output_1_array = factory.createArray<float>(
                {output_1.size(), 1}, output_1.begin().base(), output_1.end().base());
        outputs[0] = std::move(output_1_array);
        outputs[1] = getMexArray(myBundle.getFiberTd());
    }

    void
    checkInput(matlab::mex::ArgumentList &inputs)
    {
        using namespace matlab::data;
        if (inputs.size() != 6) {
            // 在matlab上打印错误日志
            matlabPtr->feval(u"error", 0,
                std::vector<Array>({ factory.createScalar(
                        "Use this function like this: \n"
                        "td, fiberTd = this_func(input_filename, output_filename, scale, numberOfSamplingDirections, tractSubSampling, fiberPointSubSampling) \n"
                        "example: this_func(\"xxx.vtk\", \"outxxx.vtk\", 3, int32(3), int32(1), int32(1))\n"
                ) }));
        }
    }

    matlab::data::Array getMexArray (const std::map<int, std::vector<float>>&& m) { // 定义一个将map转为matlab::data::Array的函数
        matlab::data::StructArray mx = factory.createStructArray ( { 1,1 }, { "keys", "values" }); // 创建一个1行1列的matlab::data::StructArray对象，包含两个字段：keys和values
        std::vector<int> keys; // 定义一个存储键的vector
        std::vector<matlab::data::Array> values; // 定义一个存储值的vector
        for (auto& p : m) { // 遍历map中的每一对键值
            keys.emplace_back(p.first); // 将键添加到vector中
            values.push_back (factory.createArray( { 1,p.second.size () }, p.second.begin (), p.second.end ())); // 将值的vector转为matlab::data::Array对象，并添加到vector中
        }
        mx [0] ["keys"] = factory.createArray ( { 1,m.size () }, keys.begin (), keys.end ()); // 将键的vector转为matlab::data::Array对象，并赋值给mx [0] ["keys"]
        mx [0] ["values"] = factory.createCellArray ( { 1,m.size () });
        matlab::data::TypedArray<matlab::data::Array> valuesArray = mx[0]["values"];
        for (size_t i{}; i < values.size(); ++i) {
            valuesArray[i] = std::move(values[i]);
        }// 将值的vector转为matlab::data::CellArray对象，并赋值给mx [0] ["values"]
        mx[0]["values"] = std::move(valuesArray);
        return mx; // 返回一个matlab::data::Array对象
    }
};

auto measureTimeMillis(const function<void()>& func) -> decltype(std::chrono::milliseconds().count())
{
    // 获取开始时间点
    auto start = std::chrono::high_resolution_clock::now();
    func();
    auto stop = std::chrono::high_resolution_clock::now();
    // 计算持续时间
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
    return duration;
}