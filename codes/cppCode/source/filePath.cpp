#include "filePath.h"


FilePath::FilePath(const string & pathName): _pathName(pathName){
    if (_pathName.back() == '/') {
        _endSlash = true;
    }
    else {
        _endSlash = false;
    }

    if (_pathName[0] == '/') {
        _beginSlash = true;
    }
    else {
        _beginSlash = false;
    }
}

FilePath FilePath::operator+(const FilePath & fp) {
    string combinedPath;
    if (_endSlash) {
        if (!fp._beginSlash) {
            combinedPath = _pathName + fp._pathName;
        }
        else {
            if (_pathName.length() >= 2)
                combinedPath = _pathName.substr(0, _pathName.length() - 1) + fp._pathName;
            else {
                // _pathName = "/"
                combinedPath = fp._pathName;
            }
        }
    }
    else {
        if (!fp._beginSlash) {
            combinedPath = _pathName + "/" + fp._pathName;
        }
        else {
            combinedPath = _pathName + fp._pathName;
        }
    }
    return FilePath(combinedPath);
}

void FilePath::operator+=(const FilePath & fp) {
    if (_endSlash) {
        _pathName += fp._pathName;
    }
    else {
        _pathName += "/" + fp._pathName;
    }
    _endSlash = fp._endSlash;
}

bool FilePath::endSlash() const {
    return _endSlash;
}

FilePath::operator std::string() const {
    return _pathName;
}
