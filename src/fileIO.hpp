#include <iostream>

#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

namespace fileIO
{
    /*
    文件读写类，仅原生c++数据格式实现，不依赖PETSc
    */
    class FileIO
    {
    public:
        FileIO();

        void VecWriteCSV(double *vec, int size, std::string filename);
        void VecWriteCSV(double *vec, int size, double *coords, std::string filename);

        // 行优先
        void MatWriteCSV(double *mat, int m, int n, std::string filename);
        void MatWriteCSV(double *mat, int m, int n, double *coords, std::string filename);

        void HistWriteCSV(std::vector<double> &hist, std::string filename);
    };
}