#ifndef _FIBER_H
#define _FIBER_H

#include <iostream>
#include <map>
#include <vector>
#include "linalg.h"

class fiber {
public:
	//访问张量的映射。键入张量名称；对3D矩阵的向量进行赋值。
	typedef std::map < std::string,
		stdMat_t,
		std::less<std::string>,
		Eigen::aligned_allocator<std::pair<const std::string, stdMat_t> > > TensorMapType;
	//以字段名为键的映射和字段内容的2D ukfPrecisionType矢量
	typedef std::map<std::string, std::vector<float> > FieldMapType;

	/** The points in RAS space of the fiber */
	stdVec_t Points;

	/** The fields of the fiber */
	FieldMapType Fields;

	/** The tensors of the fiber */
	TensorMapType Tensors;
};

#endif