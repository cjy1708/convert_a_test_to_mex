#ifndef FIBERBUNDLE_H
#define FIBERBUNDLE_H
//#include<vtk>
//#include "vtkSmartPointer.h"
//#include "vtkPolyData.h"
//#include "vtkDataArray.h"
#include<vector>
#include<map>

#include"fiber.h"
#include <vtkType.h>
#include <vtkDataArray.h>
#include <vtkPolyData.h>

class fiberbundle {
public:
    typedef std::vector<fiber> FiberVector;

    void ReadFibers(const std::string &inputFibersFileName);

    void WriteFibers(const std::string &outputFibersFileName, bool writeAscii, bool writeUnCompressed);

    FiberVector &GetFibers() { return m_FiberBundle; }

    void Print();

    [[nodiscard]]
    std::vector<float> getTd();
    [[nodiscard]]
    std::map<int, std::vector<float>> getFiberTd();

private:
    template<typename TArray>
    bool AddArray(fiber &curFiber,
                  vtkIdType nPoints,
                  const vtkIdType *pts,
                  const std::string &curname,
                  vtkDataArray *curArray) {
        TArray *array = TArray::SafeDownCast(curArray);
        if (array == 0) {
            return false;
        }
        if (array->GetNumberOfComponents() == 1) {
            std::vector<float> curVec;
            for (vtkIdType curPoint = 0; curPoint < nPoints; ++curPoint) {
                curVec.push_back(array->GetValue(pts[curPoint]));
            }
            curFiber.Fields[curname] = curVec;
        } else if (array->GetNumberOfComponents() == 9) {
            stdMat_t curVec;
            for (vtkIdType curPoint = 0; curPoint < nPoints; ++curPoint) {
                const double *curval = array->GetTuple(pts[curPoint]);
                mat33_t mat;
                for (unsigned int i = 0, v = 0; i < 3; ++i) {
                    for (unsigned int j = 0; j < 3; ++j, ++v) {
                        mat(i, j) = curval[v];
                    }
                }
                curVec.push_back(mat);
            }
            curFiber.Tensors[curname] = curVec;
        }
        return true;
    }

    FiberVector m_FiberBundle;
    vtkSmartPointer<vtkPolyData> m_PolyData;
    std::string m_InputFibersFileName;
};

#endif