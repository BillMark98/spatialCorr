#ifndef FILE_PATH_H
#define FILE_PATH_H
// class to safely add directory
#include <string>
using std::string;
class FilePath {
    private:
        string _pathName;
        // bool indicating whether the pathName ends with /
        bool _endSlash;
        // bool indicating whether the pathName begins with /
        bool _beginSlash;
    public:
        FilePath(): _pathName("./"), _endSlash(true), _beginSlash(false){};
        FilePath(const string & pathName);
        FilePath operator+(const FilePath & fp);
        void operator+=(const FilePath & fp);
        bool endSlash() const;

        // explicit typecast
        explicit operator std::string() const;
};
#endif