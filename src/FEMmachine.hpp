/*
此模块负责各类FEM的实现。在每一类FEM中，包含：
    1. 将其他接口读入的网格信息转化为FEM所需的网格数据结构。 
    2. 相关组矩阵操作
*/

#pragma once

#include <petsc.h>
#include <string>
#include <vector>
#include <array>

#include "mesh.hpp"

namespace femm
{
    class LagrangeBase
    {
        public: 

    }
}