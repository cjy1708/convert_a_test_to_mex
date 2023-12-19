// 包含头文件
#include <iostream>
#include <chrono>

#include <mex.hpp>
#include <mexAdapter.hpp>

// 自定义头文件
#include <fiberbundle.h>
#include <computedispersion.h>

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
        {
            TypedArray<char16_t> input_filename = inputs[0];
            u16string u16str{input_filename.begin(), input_filename.end()};
            inputFilename = matlab::engine::convertUTF16StringToUTF8String(u16str);
        }
        cout << "inputFilename: " << inputFilename << endl;
        string outputFilename;
        {
            TypedArray<char16_t> output_filename = inputs[1];
            u16string u16str{output_filename.begin(), output_filename.end()};
            outputFilename = matlab::engine::convertUTF16StringToUTF8String(u16str);
        }
        cout << "outputFilename: " << outputFilename << endl;
        uint scale = inputs[2][0];
        cout << "scale: " << scale << endl;
        uint numberOfSamplingDirections = inputs[3][0];
        cout << "numberOfSamplingDirections: " << numberOfSamplingDirections << endl;
        uint tractSubSampling = inputs[4][0];
        cout << "tractSubSampling: " << tractSubSampling << endl;
        uint fiberPointSubSampling = inputs[5][0];
        cout << "fiberPointSubSampling: " << fiberPointSubSampling << endl;

        fiberbundle myBundle;
        myBundle.ReadFibers(inputFilename);
        computeDispersion(myBundle, scale, numberOfSamplingDirections, outputFilename, tractSubSampling, fiberPointSubSampling);
        // 输出计算完成

        myBundle.WriteFibers(outputFilename, false, true);

        if (outputs.size() != 2) {
            // 保存输出参数
            this->td = myBundle.getTd();
            this->fiberTd = myBundle.getFiberTd();
        }
        // 输出参数第一个保存一个列向量
        auto output_1 = myBundle.getTd();
        matlab::data::TypedArray<double> output_1_array = factory.createArray<float>(
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
        for (size_t i{}; i < values.size(); ++i) {
            mx[0]["values"][i] = std::move(values[i]);
        }// 将值的vector转为matlab::data::CellArray对象，并赋值给mx [0] ["values"]
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