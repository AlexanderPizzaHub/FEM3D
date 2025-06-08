#include <iostream>

#include "fileIO.hpp"

#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

namespace fileIO
{
    FileIO::FileIO() {};

    void FileIO::VecWriteCSV(double *vec, int size, std::string filename)
    {
        std::ofstream out(filename.c_str());
        out << "value" << std::endl;
        for (int i = 0; i < size; i++)
        {
            out << vec[i] << std::endl;
        }
        out.close();
    }

    void FileIO::VecWriteCSV(double *vec, int size, double *coords, std::string filename)
    {
        /*
        assert vec aligns with coords
        */
        std::ofstream out(filename.c_str());
        out << "x y value" << std::endl;
        for (int i = 0; i < size; i++)
        {
            out << coords[3 * i] << " " << coords[3 * i + 1] << " " << coords[3 * i + 2] << vec[i] << std::endl;
        }
        out.close();
    }

    void FileIO::MatWriteCSV(double *mat, int m, int n, std::string filename)
    {
        std::ofstream out(filename.c_str());
        out << "value" << std::endl;
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                out << mat[i * m + j] << " ";
            }
            out << std::endl;
        }
        out.close();
    }

    void FileIO::MatWriteCSV(double *mat, int m, int n, double *coords, std::string filename)
    {
        std::ofstream out(filename.c_str());
        out << "x y value" << std::endl;
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                out << coords[3 * (i * m + j)] << " " << coords[3 * (i * m + j) + 1] << " " << mat[i * m + j] << std::endl;
            }
        }
        out.close();
    }

    void FileIO::HistWriteCSV(std::vector<double> &hist, std::string filename)
    {
        std::ofstream out(filename.c_str());
        out << "value" << std::endl;
        for (size_t i = 0; i < hist.size(); i++)
        {
            out << hist[i] << std::endl;
        }
        out.close();
    }

}