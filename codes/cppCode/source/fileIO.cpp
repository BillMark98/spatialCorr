#include "fileIO.h"

using namespace std; 
bool IsPathExist(const std::string &s)
{
  struct stat buffer;
  return (stat (s.c_str(), &buffer) == 0);
}

// template<typename Stream>
// int createPath( mode_t mode, const std::string& rootPath, Stream&& path)
// {
//     struct stat st;
//     const string pathName = std::forward<Stream>(path);
//     for( auto iter = pathName.begin() ; iter != pathName.end(); )
//     {
//          auto newIter = std::find( iter, pathName.end(), '/' );
//          auto newPath = rootPath + "/" + std::string( pathName.begin(), newIter);

//          if( stat( newPath.c_str(), &st) != 0)
//          {           
//              if( mkdir( newPath.c_str(), mode) != 0 && errno != EEXIST )
//              {
//                 std::cout << "cannot create folder [" << newPath << "] : " << strerror(errno) << std::endl;
//                 return -1;
//              }
//          }
//          else
//             if( !S_ISDIR(st.st_mode) )
//              {
//                  errno = ENOTDIR;
//                  std:: cout << "path [" << newPath << "] not a dir " << std::endl;
//                  return -1;
//              }
//              else
//                  std::cout << "path [" << newPath << "] already exists " << std::endl;


//          iter = newIter;
//          if( newIter != path.end() )
//              ++ iter;
//     }
//     return 0;
// }
// int createPath( mode_t mode, const std::string& rootPath, const std::string&& path )
// {
//     struct stat st;

//     for( std::string::const_iterator iter = path.begin() ; iter != path.end(); )
//     {
//          std::string::const_iterator newIter = std::find( iter, path.end(), '/' );
//          std::string newPath = rootPath + "/" + std::string( path.begin(), newIter);

//          if( stat( newPath.c_str(), &st) != 0)
//          {           
//              if( mkdir( newPath.c_str(), mode) != 0 && errno != EEXIST )
//              {
//                 std::cout << "cannot create folder [" << newPath << "] : " << strerror(errno) << std::endl;
//                 return -1;
//              }
//          }
//          else
//             if( !S_ISDIR(st.st_mode) )
//              {
//                  errno = ENOTDIR;
//                  std:: cout << "path [" << newPath << "] not a dir " << std::endl;
//                  return -1;
//              }
//              else
//                  std::cout << "path [" << newPath << "] already exists " << std::endl;


//          iter = newIter;
//          if( newIter != path.end() )
//              ++ iter;
//     }
//     return 0;
// }

// int main() 
  
// { 
  
//     // Creating a directory 
//     if (mkdir("geeksforgeeks", 0777) == -1) 
//         cerr << "Error :  " << strerror(errno) << endl; 
  
//     else
//         cout << "Directory created"; 
    
//     if (IsPathExist("geeksforgeeks")) {
//         cout << "path exists\n";
//     }
// } 