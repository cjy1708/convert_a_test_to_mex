#include <computedispersion.h>
#include <Eigen/Dense>

#include <cmath>
#include <algorithm>
#include <utility>
#include <vector>
/* This function requires three compulsory input arguments.此函数需要三个强制输入参数。

      DDF＝compute_dispersion（TRACT，SCALE，NUMDIRS）将单元阵列TRACT作为输入，
        其中每个单元是3乘N矩阵。该矩阵的行对应于一根光纤的点坐标的x、y和z分量。N是光纤中的点数。     
      SCALE是一个数值，用于指定计算色散的比例。   
      NUMDIRS是在与局部切线向量正交的平面中对色散分布函数（DDF）进行采样所沿的方向的数目。
    
   DDF是一个（NUMDIRS+1）乘-P矩阵。P是TRACT中的总点数。
    DDF的列对应于纤维束的点，它们的顺序与调用cell2mat（tract）产生的顺序相同。
    DDF的行对应于在给定光纤点计算的DDF。
    第一行NUMDIRS给出了每个采样方向的DDF值，第（NUMDIRS+1）行给出了DDF的中值，
    这对应于Savadjiev等人的TD测量。
    此外，该函数最多允许三个可选的输入参数。
    DDF=compute_dispersion（TRACT，SCALE，NUMDIRS，SUB_TR）
    对TRACT中给定的一组光纤上的计算进行子采样，因此每个SUB_TR光纤中只有一根光纤执行计算。
    SUB_TR必须是正整数，默认值=1（无子采样）。
    DDF=compute_dispersion（TRACT、SCALE、NUMDIRS、SUB_TR、SUB_FB）
    对沿每条光纤的计算进行子采样，从而仅对沿光纤的每个SUB_FB点中的一个点执行计算。
    SUB_FB必须是正整数，默认值=1（无子采样）。
    DDF=compute_dispersion（TRACT、SCALE、NUMDIRS、SUB_TR、SUB_FB、OUTPUTFILENAME）
    将输出DDF矩阵保存到.mat文件中，该文件的文件名由参数OUTPUTFILENAME提供。
    默认行为：输出不写入文件。
*/

/*  PrintMat
    注意：在调试过程中，您可以调用此函数来显示特征矩阵的内容：（gdb）调用PrintMat（matrixName）
 */
void PrintMat(Eigen::Matrix<ukfPrecisionType, Eigen::Dynamic, Eigen::Dynamic>& mat, std::ostream& outfile = std::cerr)
{
    Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
    outfile << mat.format(CleanFmt);
}

void PrintMat(Eigen::Matrix<ukfPrecisionType, Eigen::Dynamic, 1>& mat, std::ostream& outfile = std::cerr)
{
    Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
    outfile << mat.format(CleanFmt);
}

namespace {
    typedef Eigen::Matrix<ukfPrecisionType, Eigen::Dynamic, Eigen::Dynamic> MatrixType;

    typedef std::vector<MatrixType> MatrixVector;

    typedef Eigen::Matrix<ukfPrecisionType, Eigen::Dynamic, 1> VectorType;

    /* 将矩阵的向量转换为单个矩阵。
      MATLAB cell2mat更通用，我认为——这确实是 computedispersion所需要的。
     光纤被存储为NRows、numPoints矩阵；输出必须是NRows、numPoints*numFibers
     */
    void cell2mat(const MatrixVector& fibers, MatrixType& tractMatrix)
    {
        const size_t fibersSize(fibers.size());
        unsigned int numCols(0);
        for (unsigned i = 0; i < fibersSize; i++) {
            numCols += fibers[i].cols();
        }
        tractMatrix.conservativeResize(fibers[0].rows(), numCols);
        for (unsigned int i = 0, offset = 0; i < fibersSize; ++i)
        {
            const MatrixType& curFiber = fibers[i];
            tractMatrix.block(0, offset, curFiber.rows(), curFiber.cols()) = curFiber;
            offset += curFiber.cols();
        }
    }

    MatrixType PCA(const MatrixType& mat)
    {
        MatrixType centered = mat.rowwise() - mat.colwise().mean();
        MatrixType cov = centered.adjoint() * centered;

        Eigen::SelfAdjointEigenSolver<MatrixType> eig(cov);

        VectorType eigenValues = eig.eigenvalues();
        MatrixType eigenVectors = eig.eigenvectors();

        MatrixType rVal(mat.cols(), mat.cols());
        for (unsigned i = 0; i < unsigned(mat.cols()); ++i)
        {
            rVal.row(i) = eigenVectors.row((mat.cols() - 1) - i).normalized();
        }
        return rVal;
    }

    /* Go from std::vector<std::vector<vec3_t> to std::vector<MatrixType>
     输入光纤类型在光纤束中定义为点向量的向量。对于局部计算，我们需要一个矩阵向量，而不是维数* 3 的大小。
     */
    void FiberVector2EigenFiberVector(const fiberbundle::FiberVector& fv /* in */, MatrixVector& lv)
    {
        const size_t fiberSize(fv.size());

        lv.clear();                   // 可能重复使用现有矢量
        lv.resize(fiberSize);         // 与输入纤维束尺寸相同

        for (unsigned long i = 0; i < fiberSize; ++i)
        {
            const stdVec_t& curFiber = fv[i].Points;
            const size_t curFiberSize(curFiber.size());

            MatrixType& curLocalFiber = lv[i];
            curLocalFiber.conservativeResize(3, curFiberSize);
            for (unsigned int j = 0; j < curFiberSize; ++j)
            {
                curLocalFiber.col(j) = curFiber[j];
            }
        }
    }

    /* FRAME  Compute a frame field along a curve
       FRAME沿曲线计算帧场
        [T, N, B] = frame(X, Y, Z, V)
      X、 Y和Z各自是 1 * M 的矢量。它们必须分别包含曲线点的x、y和z点坐标。M是沿着曲线的点数。
        V是一个 3 * 1 的矢量。
        T、 N和B是 3 * M 矩阵。
     此函数分别返回输出向量T、N和B中沿曲线的切线、法线和二法线向量场。
     注意，法线向量和二法线向量不一定是Frenet方程定义的向量。
     由于Frenet框架在实践中是非常不稳定的——Frenet的二重向量场在曲线的拐点处是不连续的，
     并且在直线段上消失——我们使用输入向量V来约束二重向量始终位于由切线向量和约束向量V定义的平面中。
     这导致曲线上的二法线（以及因此的法线）矢量场的行为更加平滑。
     */
    class Frame {
    public:
        typedef Eigen::Matrix<ukfPrecisionType, 3, 1> FixedVectorType;
        Frame(VectorType xCoords, VectorType yCoords, VectorType  zCoords,
            const VectorType& constraintVector) : m_Xcoords(std::move(xCoords)),
            m_Ycoords(std::move(yCoords)),
            m_Zcoords(std::move(zCoords)),
            m_ConstraintVector(constraintVector)
        {
        }
        void compute();
        const MatrixType& GetTangent()  const { return m_Tangent; }
        const MatrixType& GetNormal()   const { return m_Normal; }
        const MatrixType& GetBinormal() const { return m_Binormal; }
    private:
        VectorType       m_Xcoords;
        VectorType       m_Ycoords;
        VectorType       m_Zcoords;
        FixedVectorType  m_ConstraintVector;
        MatrixType       m_Tangent;
        MatrixType       m_Normal;
        MatrixType       m_Binormal;
    };

    void
        Frame
        ::compute()
    {
        unsigned int curveLength = this->m_Xcoords.rows();
        if (curveLength == 1)
        {
            this->m_Xcoords = this->m_Xcoords.transpose();
            this->m_Ycoords = this->m_Ycoords.transpose();
            this->m_Zcoords = this->m_Zcoords.transpose();
            curveLength = this->m_Xcoords.rows();
        }

        if (curveLength == 0)
        {
            this->m_Tangent(0) = 1.0;
            this->m_Tangent(1) = 0.0;
            this->m_Tangent(2) = 0.0;

            this->m_Normal(0) = 0.0;
            this->m_Normal(1) = 1.0;
            this->m_Normal(2) = 0.0;

            this->m_Binormal(0) = 0.0;
            this->m_Binormal(1) = 0.0;
            this->m_Binormal(2) = 1.0;
            return;
        }
        this->m_Tangent.resize(curveLength, 3);
        this->m_Normal.resize(curveLength, 3);
        this->m_Binormal.resize(curveLength, 3);
        this->m_Tangent.setZero();
        this->m_Normal.setZero();
        this->m_Binormal.setZero();

        MatrixType points(curveLength, 3);
        points.col(0) = this->m_Xcoords;
        points.col(1) = this->m_Ycoords;
        points.col(2) = this->m_Zcoords;

        for (unsigned i = 1; i < curveLength - 1; ++i)
        {
            this->m_Tangent.row(i) = (points.row(i + 1) - points.row(i - 1)) / 2;
            double tangentNorm = this->m_Tangent.row(i).norm();
            if (tangentNorm > 0.00001)
            {
                this->m_Tangent.row(i).normalize();
            }
            else
            {
                this->m_Tangent.row(i) = this->m_Tangent.row(i - 1);
            }
        }

        this->m_Tangent.row(0) = points.row(1) - points.row(0);
        double tangentNorm = this->m_Tangent.row(0).norm();
        if (tangentNorm > 0.00001)
        {
            this->m_Tangent.row(0).normalize();
        }
        else
        {
            this->m_Tangent.row(0) = this->m_Tangent.row(1);
        }

        this->m_Tangent.row(curveLength - 1) = points.row(curveLength - 1) -
            points.row(curveLength - 2);
        tangentNorm = this->m_Tangent.row(curveLength - 1).norm();
        if (tangentNorm > 0.00001)
        {
            this->m_Tangent.row(curveLength - 1).normalize();
        }
        else
        {
            this->m_Tangent.row(curveLength - 1) = this->m_Tangent.row(curveLength - 2);
        }

        for (unsigned i = 0; i < curveLength; ++i)
        {
            FixedVectorType curTanRow = this->m_Tangent.row(i);
            FixedVectorType orthogonalVector1 = this->m_ConstraintVector.cross(curTanRow);
            if (orthogonalVector1.norm() < 0.00001)
            {
                this->m_ConstraintVector *= 0.00002;
                this->m_ConstraintVector.normalize();
                orthogonalVector1 = this->m_ConstraintVector.cross(curTanRow);
            }
            FixedVectorType orthogonalVector2 = curTanRow.cross(orthogonalVector1).normalized();
            if (orthogonalVector2.dot(this->m_ConstraintVector) < 0.00001)
            {
                this->m_Binormal.row(i) = orthogonalVector2 * -1.0;
            }
            else
            {
                this->m_Binormal.row(i) = orthogonalVector2;
            }
        }

        for (unsigned int i = 0; i < curveLength; ++i)
        {
            FixedVectorType curTanRow = this->m_Tangent.row(i);
            FixedVectorType curBinorm = this->m_Binormal.row(i);
            this->m_Normal.row(i) = curBinorm.cross(curTanRow);
        }
    }

    /** RotateField
    通过旋转矩阵R旋转束。因此，整个纤维束在点currentPosition的局部坐标系中表示，使得局部纤维切线向量为x轴。
    返回rotatedPointCoordinates中的旋转点坐标，以及rotatedTangentVectorField中相应的旋转切线向量。
     */
    class RotateField
    {
    public:
        RotateField(const MatrixType& tractMatrix,
            const MatrixType& rotationMatrix,
            const MatrixType& currentPosition) :
            m_TractMatrix(tractMatrix),
            m_RotationMatrix(rotationMatrix),
            m_CurrentPosition(currentPosition)
        {
        }
        void compute();
        const MatrixType& GetRotatedPointCoordinates() const
        {
            return m_RotatedPointCoordinates;
        }
        const MatrixType& GetRotatedTangentVectorField() const
        {
            return m_RotatedTangentVectorField;
        }

    private:
        MatrixType m_RotatedPointCoordinates;
        MatrixType m_RotatedTangentVectorField;
        const MatrixType& m_TractMatrix;
        const MatrixType& m_RotationMatrix;
        const MatrixType& m_CurrentPosition;
    };

    void
        RotateField
        ::compute()
    {
        const unsigned long tractCols(this->m_TractMatrix.cols());
        MatrixType currentPositionMatrix(3, tractCols);
        for (unsigned i = 0; i < unsigned(this->m_TractMatrix.cols()); ++i)
        {
            currentPositionMatrix.col(i) = this->m_CurrentPosition;
        }
        MatrixType pointCoordinates = m_TractMatrix.block(0, 0, 3, m_TractMatrix.cols()) - currentPositionMatrix;
        this->m_RotatedPointCoordinates = this->m_RotationMatrix * pointCoordinates;
        this->m_RotatedPointCoordinates += currentPositionMatrix;
        this->m_RotatedTangentVectorField =
            this->m_RotationMatrix * this->m_TractMatrix.block(3, 0, 3, tractCols);
    }

    /* computeMeanVector
      计算以currentPosition位置为中心的具有半径刻度的圆盘形邻域内的平均矢量。
      参见Savadjiev等人MICCAI 2012，方程式2，以及第3.1节中的相关讨论。
     */
    bool computeMeanVector(const MatrixType& xCoordinates,
        const MatrixType& yCoordinates,
        const MatrixType& zCoordinates,
        const MatrixType& tangentVectorField,
        const MatrixType& currentPosition,
        double scale,
        Eigen::Vector3d& meanVector /* out */)
    {
        const double eps(2.2204e-16); // 这就是matlab所认为的eps。

        std::vector<bool> indexPointsInPlane(xCoordinates.cols(), false);
        unsigned int inPlanePointCount = 0;
        for (unsigned i = 0; i < unsigned(xCoordinates.cols()); ++i)
        {
            if (fabs(xCoordinates(i) - currentPosition(0)) < 0.5)
            {
                indexPointsInPlane[i] = true;
                ++inPlanePointCount;
            }
        }

        //这实现了matlab的疯狂：其中b是逻辑向量的A（b）是所有点Z的向量，其中b[Z]为真
        // b是一个逻辑向量。 A（b）是所有点Z的向量，其中b[Z]为真
        MatrixType pointCoordinates(3, inPlanePointCount);
        for (unsigned i = 0, curPoint = 0; i < unsigned(xCoordinates.cols()); ++i)
        {
            if (indexPointsInPlane[i])
            {
                pointCoordinates(0, curPoint) = xCoordinates(i) - currentPosition(0);
                pointCoordinates(1, curPoint) = yCoordinates(i) - currentPosition(1);
                pointCoordinates(2, curPoint) = zCoordinates(i) - currentPosition(2);
                ++curPoint;
            }
        }

        MatrixType pointDistance(1, inPlanePointCount);

        // indexPointsInDisk = pointDistance < scale*scale;
        std::vector<bool> indexPointsInDisk(inPlanePointCount);

        unsigned int inDiskPointCount = 0;
        for (unsigned int i = 0; i < inPlanePointCount; ++i)
        {
            double dotProd = pointCoordinates.col(i).dot(pointCoordinates.col(i));
            pointDistance(i) = dotProd;

            if (pointDistance(i) < scale * scale)
            {
                indexPointsInDisk[i] = true;
                ++inDiskPointCount;
            }
            else
            {
                indexPointsInDisk[i] = false;
            }
        }

        double n = 0.0;
        //这是对matlab代码的重新排列；如果磁盘中没有点，那么计算平均值就没有意义。
        if (inDiskPointCount > 10)
        {
            // vectorFieldInDisk = vectorFieldInDisk(:,indexPointsInDisk);
            // use second variable rather than overwrite in place
            // vectorFieldInDisk = tangentVectorField(:,indexPointsInPlane);
            MatrixType vectorFieldInDisk(3, inPlanePointCount);
            for (unsigned i = 0, curPoint = 0; i < unsigned(tangentVectorField.cols()); ++i)
            {
                if (indexPointsInPlane[i])
                {
                    vectorFieldInDisk.col(curPoint) = tangentVectorField.col(i);
                    curPoint++;
                }
            }
            MatrixType vectorFieldInDisk2(3, inDiskPointCount);
            for (unsigned int i = 0, curPoint = 0; i < inPlanePointCount; ++i)
            {
                if (indexPointsInDisk[i])
                {
                    vectorFieldInDisk2.col(curPoint) = vectorFieldInDisk.col(i);
                    ++curPoint;
                }
            }
            // indexNegativeOrientations = vectorFieldInDisk(1,:) < -eps;
            // vectorFieldInDisk(:,indexNegativeOrientations) = -vectorFieldInDisk(:,indexNegativeOrientations);
            // in other words, any  negative points multiply by scalar -1
            for (unsigned int i = 0; i < inDiskPointCount; ++i)
            {
                if (vectorFieldInDisk2(0, i) < -eps)
                {
                    vectorFieldInDisk2.col(i) *= -1.0;
                }
            }
            // meanVector = mean(vectorFieldInDisk,2);
            meanVector(0) = vectorFieldInDisk2.row(0).mean();
            meanVector(1) = vectorFieldInDisk2.row(1).mean();
            meanVector(2) = vectorFieldInDisk2.row(2).mean();
            n = meanVector.norm();
        }
        //
        // MATLAB代码将取空矩阵的平均值，这导致了充满NaN的平均值；
        // 我的解决方案是根本不尝试计算均值，因为Eigen不喜欢对空矩阵进行运算。
        // if(~any(isnan(meanVector)) && (length(indexPointsInDisk) > 10) )
        //     n = norm(meanVector);
        //     if(n>0.1)
        if (inDiskPointCount > 10 && n > 0.1)
        {
            meanVector /= n;
            return true;
        }
        // 如果平面中没有足够的点，或者磁盘中没有点，那么只返回状态0来指示这一点。
        meanVector = Eigen::Vector3d::Zero();
        return false;
    }

    double median(std::vector<double>& vec)
    {
        double rval;
        size_t size = vec.size();
        std::sort(vec.begin(), vec.end());
        if (size % 2 == 0)
        {
            rval = (vec[size / 2 - 1] + vec[size / 2]) / 2;
        }
        else
        {
            rval = vec[size / 2];
        }
        return rval;
    }
}


int computeDispersion(fiberbundle& bundle, double scale,
    unsigned int numberOfSamplingDirections,
    const std::string& outputFilename,
    unsigned int  tractSubSampling,
    unsigned int fiberPointSubSampling)
{
    // bundle.Print();
    const double Pi = 3.14;

    fiberbundle::FiberVector& fibers = bundle.GetFibers();

    // 本征矩阵的向量
    MatrixVector lFibers;
    FiberVector2EigenFiberVector(fibers, lFibers);

    MatrixType tractMatrix; // 3 x x matrix accumulating all fibers tracts
    cell2mat(lFibers, tractMatrix);

    MatrixType princomp = PCA(tractMatrix.transpose());
    VectorType constraintVector = princomp.col(2);
    
    /*计算每条光纤上每个点的色散
    注：如果tractSubSampling或fiberPointSubSamplication>1，则它打破了每个点一个离散度的假设
    注意：这遵循matlab将tan/nom/binorm存储在curPoints矩阵中的点之后。它们可以很容易地分开，并且可能是更清晰的代码。
    */
    for (unsigned int i = 0; i < fibers.size(); ++i)
    {
        MatrixType& curPoints = lFibers[i];

        const unsigned int curPointsSize(curPoints.cols());

        curPoints.conservativeResize(12, curPointsSize);

        VectorType x = curPoints.row(0);
        VectorType y = curPoints.row(1);
        VectorType z = curPoints.row(2);

        Frame frame(x, y, z, constraintVector);
        frame.compute();

        curPoints.block(3, 0, 3, curPointsSize) = frame.GetTangent().transpose();
        curPoints.block(6, 0, 3, curPointsSize) = frame.GetNormal().transpose();
        curPoints.block(9, 0, 3, curPointsSize) = frame.GetBinormal().transpose();
    }

    cell2mat(lFibers, tractMatrix);

    //保留每根纤维束对纤维束进行二次采样
    MatrixVector subSampledTract;
    for (unsigned int i = 0; i < lFibers.size(); i += tractSubSampling)
    {
        subSampledTract.push_back(lFibers[i]);
    }
    MatrixType subSampledTractMatrix;
    cell2mat(subSampledTract, subSampledTractMatrix);
    //  初始化DDF矩阵并在标准坐标系中计算采样方向，使得采样方向位于YZ平面中。
    //用-1s初始化DDF矩阵——当平面/磁盘中没有足够的点来给出平均结果时的重要值。
    MatrixType DistributionValues = MatrixType::Ones(numberOfSamplingDirections + 1, subSampledTractMatrix.cols()) * -1;
    MatrixType samplingDirections = MatrixType::Zero(3, numberOfSamplingDirections);
    // 计算采样方向
    double theta = (2.0 * Pi) / static_cast<double>(numberOfSamplingDirections);
    for (unsigned int j = 1; j <= numberOfSamplingDirections; ++j)
    {
        samplingDirections(1, j - 1) = cos(j * theta);
        samplingDirections(2, j - 1) = sin(j * theta);
    }
    // 计算每个纤维点的DDF,沿着（可能是二次采样的）束的每个纤维进行第*次采样。
    for (unsigned i = 0; i < unsigned(subSampledTractMatrix.cols()); i += fiberPointSubSampling)
    {
        MatrixType rotationMatrix(3, 3);
        rotationMatrix.row(0) = subSampledTractMatrix.block(3, i, 3, 1).transpose();
        rotationMatrix.row(1) = subSampledTractMatrix.block(6, i, 3, 1).transpose();
        rotationMatrix.row(2) = subSampledTractMatrix.block(9, i, 3, 1).transpose();
        MatrixType currentPosition = subSampledTractMatrix.block(0, i, 3, 1);
        RotateField rotateField(tractMatrix, rotationMatrix, currentPosition);
        rotateField.compute();
        const MatrixType& rotatedPointCoordinates = rotateField.GetRotatedPointCoordinates();
        const MatrixType& rotatedTangentVectorField = rotateField.GetRotatedTangentVectorField();
        MatrixType xCoordinates = rotatedPointCoordinates.row(0);
        MatrixType yCoordinates = rotatedPointCoordinates.row(1);
        MatrixType zCoordinates = rotatedPointCoordinates.row(2);
        Eigen::Vector3d referenceMeanVector;
        if (computeMeanVector(xCoordinates,
            yCoordinates,
            zCoordinates,
            rotatedTangentVectorField,
            currentPosition,
            scale,
            referenceMeanVector))
        {
            MatrixType samplingPosition(3, samplingDirections.cols());
            for (unsigned j = 0; j < unsigned(samplingDirections.cols()); j++)
            {
                samplingPosition.col(j) = currentPosition + (scale * samplingDirections.col(j));
            }
            for (unsigned int j = 0; j < unsigned(samplingDirections.cols()); j++)
            {
                Eigen::Vector3d meanVector;
                if (computeMeanVector(xCoordinates,
                    yCoordinates,
                    zCoordinates,
                    rotatedTangentVectorField,
                    samplingPosition.col(j),
                    scale,
                    meanVector))
                {
                    double dot = meanVector.dot(referenceMeanVector);
                    double acosDot = acos(dot);
                    if (acosDot < 0.0)
                    {
                        acosDot *= -1.0;
                    }
                    DistributionValues(j, i) = acosDot;
                }
            }
            /* 取计算平均值的中值   */
            MatrixType pointDDF = DistributionValues.block(0, i, numberOfSamplingDirections, 1);
            std::vector<double> nonNegDDF;
            for (unsigned int j = 0; j < unsigned(pointDDF.rows()); ++j)
            {
                if (pointDDF(j, 0) != -1)
                {
                    nonNegDDF.push_back(pointDDF(j, 0));
                }
            }
            if (!nonNegDDF.empty())
            {
                DistributionValues(numberOfSamplingDirections, i) = median(nonNegDDF);
            }
        }
    }

    MatrixType DDFOutput = DistributionValues.row(numberOfSamplingDirections);
    // 因此，处理子采样的“双关语”是跳过任何光纤
    for (unsigned int i = 0, curPoint = 0; i < fibers.size(); ++i)
    {
        fiber& curFiber = fibers[i];
        std::vector<float> curDDF;
        // 如果对该光纤进行了采样，则在该光纤的点处打印出DDF
        if (i % tractSubSampling == 0)
        {
            for (unsigned int j = 0; j < curFiber.Points.size(); ++j)
            {
                curDDF.push_back(DDFOutput(curPoint));
                ++curPoint;
            }
        }
        // 如果跳过了此光纤，则为每个点输出-1s。
        else
        {
            for (unsigned j = 0; j < curFiber.Points.size(); ++j)
            {
                curDDF.push_back(-1.0);
            }
        }
        curFiber.Fields["DDF"] = curDDF;
    }

    if (!outputFilename.empty())
    {
        std::ofstream outfile((outputFilename + ".mat").c_str());
        PrintMat(DDFOutput, outfile);
        outfile.close();
    }

    return 0; // success
}

