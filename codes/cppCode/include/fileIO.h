#ifndef FILE_IO_H
#define FILE_IO_H
#include <iostream>
#include <fstream>
#include <cstdint>
#include <string>

// #include "healpixCorrelation.h"

#include <bits/stdc++.h> 
#include <sys/stat.h> 
#include <sys/types.h> 
// Deprecated
// #include <experimental/filesystem>
// namespace fs = std::filesystem;
 
// void demo_exists(const fs::path& p, fs::file_status s = fs::file_status{})
// {
//     std::cout << p;
//     if(fs::status_known(s) ? fs::exists(s) : fs::exists(p))
//         std::cout << " exists\n";
//     else
//         std::cout << " does not exist\n";
// }

bool IsPathExist(const std::string &s);
template<typename Stream>
int createPath( mode_t mode, const std::string& rootPath, Stream&& path)
{
    struct stat st;
    const std::string pathName = std::forward<Stream>(path);
    for( auto iter = pathName.begin() ; iter != pathName.end(); )
    {
         auto newIter = std::find( iter, pathName.end(), '/' );
         auto newPath = rootPath + "/" + std::string( pathName.begin(), newIter);

         if( stat( newPath.c_str(), &st) != 0)
         {           
             if( mkdir( newPath.c_str(), mode) != 0 && errno != EEXIST )
             {
                std::cout << "cannot create folder [" << newPath << "] : " << strerror(errno) << std::endl;
                return -1;
             }
         }
         else
            if( !S_ISDIR(st.st_mode) )
             {
                 errno = ENOTDIR;
                 std:: cout << "path [" << newPath << "] not a dir " << std::endl;
                 return -1;
             }
             else
                 std::cout << "path [" << newPath << "] already exists " << std::endl;


         iter = newIter;
         if( newIter != pathName.end() )
             ++ iter;
    }
    return 0;
}


// std::string joinPath()

#endif