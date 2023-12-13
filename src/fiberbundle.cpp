#include "vtkPointData.h"
#include "vtkDataArray.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkPolyData.h"
#include "vtkCellArray.h"
#include "vtkLine.h"
#include "vtkSmartPointer.h"
#include <fstream>

#include <fiberbundle.h>
#include <vtkPolyDataIO.h>

typedef fiberbundle fiberbundle1;

void
fiberbundle1
::ReadFibers(const std::string& inputFibersFileName)
{
    this->m_InputFibersFileName = inputFibersFileName;
    //Read a vtk polydata file and convert to a internal usable structure for the algorithm.
    // and fill in the m_FiberBundle private elment.
    try
    {
        m_PolyData = vtkPolyDataIO::Read(inputFibersFileName.c_str());
    }
    catch (...)
    {
        std::cerr << "Can't read " << inputFibersFileName << std::endl;
        throw;
    }
    // one line per fiber
    int numLines = m_PolyData->GetNumberOfLines();

    vtkCellArray* curLines = m_PolyData->GetLines();
    curLines->InitTraversal();

    vtkPointData* pd = m_PolyData->GetPointData();

    int numComponents = pd->GetNumberOfComponents();

    for (vtkIdType curLine = 0; curLine < numLines; ++curLine)
    {
        fiber curFiber;
        vtkIdType nPoints;
        const vtkIdType *pts;
        if (curLines->GetNextCell(nPoints, pts) == 0)
        {
            break;
        }
        // one point array per polydata, and each line == a fiber
        for (vtkIdType curPoint = 0; curPoint < nPoints; ++curPoint)
        {
            vec3_t cur;
            const double* vtkCur = m_PolyData->GetPoint(pts[curPoint]);
            cur[0] = vtkCur[0];
            cur[1] = vtkCur[1];
            cur[2] = vtkCur[2];
            curFiber.Points.push_back(cur);
        }
        for (int curComp = 0; curComp < numComponents; ++curComp)
        {
            std::string curname = pd->GetArrayName(curComp);
            vtkDataArray* curArray = pd->GetArray(curname.c_str());
            if (!this->AddArray<vtkFloatArray>(curFiber, nPoints, pts, curname, curArray))
            {
                this->AddArray<vtkDoubleArray>(curFiber, nPoints, pts, curname, curArray);
            }
        }
        // TODO tensors
        m_FiberBundle.push_back(curFiber);
    }
}

/**
 *
 * @param outputFibersFileName
 * @param writeAscii
 * @param writeUnCompressed
 */
void
fiberbundle1
::WriteFibers(const std::string& outputFibersFileName, bool writeAscii, bool writeUnCompressed)
{
    //Write a vtk polydata.  Have to populate from the Fiber
    // allocate empty polyData
    this->m_PolyData = vtkSmartPointer<vtkPolyData>::New();

    vtkIdType num_points = 0;
    vtkIdType counter = 0;
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();

    for (auto & it : this->m_FiberBundle)
    {
        num_points += it.Points.size();
    }

    for (auto & it : this->m_FiberBundle)
    {
        int fiber_size = static_cast<int>(it.Points.size());
        auto* ids = new vtkIdType[fiber_size];
        for (vtkIdType curPoint = 0; curPoint < fiber_size; ++curPoint, ++counter)
        {
            vec3_t curPt = it.Points[curPoint];
            // add current point
            points->InsertNextPoint(curPt[0], curPt[1], curPt[2]);
            // add id of current point to id array
            ids[curPoint] = counter;
        }
        vtkSmartPointer<vtkLine> curLine = vtkSmartPointer<vtkLine>::New();
        curLine->Initialize(fiber_size, ids, points);
        lines->InsertNextCell(curLine);
        delete[] ids;
    }
    m_PolyData->SetPoints(points);
    m_PolyData->SetLines(lines);

    // do the per-point attributes
    fiber::FieldMapType AllDecorators;
    fiber::TensorMapType AllTensors;

    for (auto & it : this->m_FiberBundle)
    {
        for (fiber::FieldMapType::const_iterator it2 = it.Fields.begin();
            it2 != it.Fields.end(); ++it2)
        {
            const std::string& attName = it2->first;
            const std::vector<float>& array = it2->second;
            for (float it3 : array)
            {
                AllDecorators[attName].push_back(it3);
            }
        }
        for (fiber::TensorMapType::const_iterator it2 = it.Tensors.begin(); it2 != it.Tensors.end(); ++it2)
        {
            const std::string& attName = it2->first;
            const stdMat_t& tensors = it2->second;
            for (const auto & tensor : tensors)
            {
                AllTensors[attName].push_back(tensor);
            }
        }
    }

    vtkPointData* pd = this->m_PolyData->GetPointData();

    for (std::map<std::string, std::vector<float> >::const_iterator it = AllDecorators.begin();
        it != AllDecorators.end(); ++it)
    {
        vtkSmartPointer<vtkFloatArray> curAtt = vtkSmartPointer<vtkFloatArray>::New();
        curAtt->SetNumberOfComponents(1);
        curAtt->Allocate(it->second.size());
        curAtt->SetName(it->first.c_str());
        for (vtkIdType j = 0; j < static_cast<vtkIdType>(it->second.size()); ++j)
        {
            curAtt->InsertNextValue(it->second[j]);
        }
        if (it->first == "DDF")
        {
            pd->SetScalars(curAtt);
        }
        pd->AddArray(curAtt);
    }
    // TODO: do the tensors.
    for (std::map<std::string, stdMat_t>::const_iterator it = AllTensors.begin(); it != AllTensors.end(); ++it)
    {
        vtkSmartPointer<vtkFloatArray> curAtt = vtkSmartPointer<vtkFloatArray>::New();
        curAtt->SetNumberOfComponents(9);
        curAtt->Allocate(it->second.size() * 9);
        curAtt->SetName(it->first.c_str());
        float tmp[9];
        for (const auto & it2 : it->second)
        {
            for (unsigned int i = 0, v = 0; i < 3; ++i)
            {
                for (unsigned int j = 0; j < 3; ++j, ++v)
                {
                    tmp[v] = it2(i, j);
                }
            }
            curAtt->InsertNextTuple(tmp);
        }
        pd->SetTensors(curAtt);
    }

    try
    {
        vtkPolyDataIO::Write(m_PolyData,
            outputFibersFileName.c_str(),
            !writeAscii, // binary
            !writeUnCompressed  // compressed
        );
    }
    catch (...)
    {
        std::cerr << "Can't read " << outputFibersFileName << std::endl;
        throw;
    }
}

void
fiberbundle1
::Print()
{
    std::cerr << this->m_InputFibersFileName << ":" << std::endl;
    for (FiberVector::const_iterator it = this->m_FiberBundle.begin();
        it != this->m_FiberBundle.end(); ++it)
    {
        const fiber& curFiber = *it;
        std::cerr << "Fiber " << it - this->m_FiberBundle.begin() << ":" << std::endl;
        std::cerr << "  Points:" << std::endl;
        for (unsigned int i = 0; i < curFiber.Points.size(); ++i)
        {
            std::cerr << "    ["
                << curFiber.Points[i][0] << " "
                << curFiber.Points[i][1] << " "
                << curFiber.Points[i][2] << " "
                << "]" << std::endl;
        }
        for (fiber::FieldMapType::const_iterator curFiberIt = curFiber.Fields.begin();
            curFiberIt != curFiber.Fields.end(); ++curFiberIt)
        {
            std::cerr << "  " << curFiberIt->first << ":" << std::endl << "    ";
            for (unsigned int i = 0; i < curFiberIt->second.size(); ++i)
            {
                std::cerr << curFiberIt->second[i];
                if (i == curFiberIt->second.size() - 1)
                {
                    std::cerr << endl;
                }
                else
                    std::cerr << ' ';
            }
        }
    }
}

using namespace std;
vector<float>
fiberbundle1
::getTd()
{
    ofstream outfile("TD.txt", ios::trunc);
    
    for (const auto & curFiber : this->m_FiberBundle)
    {
        for (const auto & Field : curFiber.Fields)
        {
            if (Field.first == "DDF") {
                for (float i : Field.second)
                {
                    outfile << i << "\n";
                }
            }
        }
    }
    outfile.close();
}

vector<vector<float>>
fiberbundle1
::getFiberTd()
{
    ofstream outfile("fiber_TD.txt", ios::trunc);
    //outfile << this->m_InputFibersFileName << ":" << std::endl;
    for (FiberVector::const_iterator it = this->m_FiberBundle.begin();
        it != this->m_FiberBundle.end(); ++it)
    {
        const fiber& curFiber = *it;
        outfile << "Fiber " << it - this->m_FiberBundle.begin() << ":" << std::endl;
        for (fiber::FieldMapType::const_iterator curFiberIt = curFiber.Fields.begin();
            curFiberIt != curFiber.Fields.end(); ++curFiberIt)
        {
            if (curFiberIt->first == "DDF") {
                for (unsigned int i = 0; i < curFiberIt->second.size(); ++i)
                {
                    outfile << curFiberIt->second[i] << "  ";
                }
            }
        }
        outfile << "\n";
    }
    outfile.close();
}