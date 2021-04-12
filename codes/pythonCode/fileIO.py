# import pandas as pd
import os
import re
import json
sizeInf = 1000000

# concentration , samplingrate, trajectory, molecule
def containsFiles(currentPath, extension = ".mat"):
    """
    check if the currentPath contains file with extension ".mat"
    ------
    Parameters:
    ----
    
    currentPath: str
    extension:str
    default ".mat"
    """
    
    for file in os.listdir(currentPath):
        if file.endswith(extension):
            return True
        
    return False

def getSubDir(currentPath, extension = ".mat", beginWith = '', mode = 'absolute', blacklist = []):
    """
    Get all subdirectories of the current path, that contains the ".mat" file
    ------
    Parameters
    -----
    currentPath: str
    string that indicate the current working directory
    beginWith: str
    only select those that begins with the suffix beginWith

    mode: str
        indicating the format of the returned path
        either absolute or relative (w.r.t currentPath)
    
    blacklist: list of str
        list of strings for directories not to be considered, meaning that only consider those possible directories that are not listed in the blacklist
        default empty        
    """

    dirLists = []
    print(currentPath)
    temp = os.walk(currentPath)
    print(temp)
    if (not os.path.exists(currentPath)):
        print("path not existing")
    for root, dirs,_ in os.walk(currentPath):

        for d in dirs:
            if containsFiles(os.path.join(root,d), extension):
                if(len(beginWith) == 0 or d.startswith(beginWith)):
                    if os.path.join(root,d) not in blacklist:
                        if (mode == 'absolute'):
                            dirLists.append(os.path.join(root,d))
                        elif (mode == 'relative'):
                            dirLists.append(d)
                    
    return dirLists    

def getFilesFromDirs(fileDirs, extension = ".mat", fileNumEachDir = sizeInf, getSubDirsBySelf = False, beginWith = ''):
    """
    given list of dirs, from each dir retrieve file that ends with extension, up to fileNumEachDir

    --------
    Param
    --------
    fileDirs: list of str

    getSubDirsBySelf: bool
        default is False, so will not generate subdirs in which file with extension exists,
        so already assume the fileDirs are the one we need

        if set to True, will use fileDirs to first find out the subDirs and then find files

    beginWith: str
        will be useful if getSubDirsBySelf is true (argument for the function getSubDir)
    """
    fileNames = []
    if (getSubDirsBySelf == True):
        fileDest = []
        for possibleFileDir in fileDirs:
            tempDir = getSubDir(possibleFileDir, extension, beginWith)
            fileDest += tempDir
        fileDirs = fileDest
    for fileDir in fileDirs:
        tempNames = getFiles(fileDir, extension, fileNum = fileNumEachDir)
        fileNames += tempNames
    return fileNames

def getFiles(fileDir, extension= ".mat", fileNum = sizeInf):
    """
    pick a random number of files (.xvg) from the current fileDir
    if no number is specified, all files will be returned
    ------
    Parameters
    -----
    fileDir: str
    the current fileDir

    extension: str
    extension of a file, default ".mat"

    fileNum: int
    the number of files extract

    return:
    a list of file names
    """
    
    # to do make it random
    fileLen = len(os.listdir(fileDir))
    fileTotalNum = min(fileLen, fileNum)
    count = 1
    fileNames = []
    for file in os.listdir(fileDir):
        if (count > fileTotalNum):
            break
        if file.endswith(extension):
            fileNames.append(os.path.join(fileDir,file))
            count += 1
    return fileNames

if __name__ == '__main__':
    # filePathPref = "/home/panwei/Desktop/bachelorThese/data/rss_superC/rss_measurement_evaluation_toolbox_Panwei/SuperC"
    # filePathPref = "/home/panwei/Desktop/bachelorThese/data/rss_superC/rss_measurement_evaluation_toolbox_Panwei/"

    # txfolder = "TX1_R27"
    # my_path = os.path.join(filePathPref,txfolder)

    # Example using rss_superC data
    # filePathPref = "/home/panwei/Desktop/bachelorThese/data/rss_superC/rss_measurement_evaluation_toolbox_Panwei/"
    # my_path = filePathPref
    # subDir = getSubDir(my_path, beginWith='TX1')

    # # Example using rss_raytracing_framework
    filePathPref = "/home/panwei/Desktop/bachelorThese/data/rss_raytracing_evaluation_framework/InputData/Frankfurt/"
    # filePathPref = "/home/panwei/Desktop/bachelorThese/data/rss_raytracing_evaluation_framework/InputData/Seoul/"
    # my_path = filePathPref
    # # subDir = getSubDir(my_path, mode = "relative", blacklist=["/home/panwei/Desktop/bachelorThese/data/rss_raytracing_evaluation_framework/InputData/Frankfurt/1564986", "/home/panwei/Desktop/bachelorThese/data/rss_raytracing_evaluation_framework/InputData/Frankfurt/4465066", "/home/panwei/Desktop/bachelorThese/data/rss_raytracing_evaluation_framework/InputData/Frankfurt/5955386", "/home/panwei/Desktop/bachelorThese/data/rss_raytracing_evaluation_framework/InputData/Frankfurt/1774496", "/home/panwei/Desktop/bachelorThese/data/rss_raytracing_evaluation_framework/InputData/Frankfurt/3163526", "/home/panwei/Desktop/bachelorThese/data/rss_raytracing_evaluation_framework/InputData/Frankfurt/3594326", "/home/panwei/Desktop/bachelorThese/data/rss_raytracing_evaluation_framework/InputData/Frankfurt/5751526", "/home/panwei/Desktop/bachelorThese/data/rss_raytracing_evaluation_framework/InputData/Frankfurt/2376016", "/home/panwei/Desktop/bachelorThese/data/rss_raytracing_evaluation_framework/InputData/Frankfurt/1625746", "/home/panwei/Desktop/bachelorThese/data/rss_raytracing_evaluation_framework/InputData/Frankfurt/1193486", "/home/panwei/Desktop/bachelorThese/data/rss_raytracing_evaluation_framework/InputData/Frankfurt/5203856", "/home/panwei/Desktop/bachelorThese/data/rss_raytracing_evaluation_framework/InputData/Frankfurt/5101656", "/home/panwei/Desktop/bachelorThese/data/rss_raytracing_evaluation_framework/InputData/Frankfurt/3552206", "/home/panwei/Desktop/bachelorThese/data/rss_raytracing_evaluation_framework/InputData/Frankfurt/5392746", "/home/panwei/Desktop/bachelorThese/data/rss_raytracing_evaluation_framework/InputData/Frankfurt/4581536", "/home/panwei/Desktop/bachelorThese/data/rss_raytracing_evaluation_framework/InputData/Frankfurt/182996", "/home/panwei/Desktop/bachelorThese/data/rss_raytracing_evaluation_framework/InputData/Frankfurt/2805226", "/home/panwei/Desktop/bachelorThese/data/rss_raytracing_evaluation_framework/InputData/Frankfurt/5285056", "/home/panwei/Desktop/bachelorThese/data/rss_raytracing_evaluation_framework/InputData/Frankfurt/6064536", "/home/panwei/Desktop/bachelorThese/data/rss_raytracing_evaluation_framework/InputData/Frankfurt/3453006", "/home/panwei/Desktop/bachelorThese/data/rss_raytracing_evaluation_framework/InputData/Frankfurt/4362816", "/home/panwei/Desktop/bachelorThese/data/rss_raytracing_evaluation_framework/InputData/Frankfurt/5836116", "/home/panwei/Desktop/bachelorThese/data/rss_raytracing_evaluation_framework/InputData/Frankfurt/5853826", "/home/panwei/Desktop/bachelorThese/data/rss_raytracing_evaluation_framework/InputData/Frankfurt/4475786", "/home/panwei/Desktop/bachelorThese/data/rss_raytracing_evaluation_framework/InputData/Frankfurt/2261176", "/home/panwei/Desktop/bachelorThese/data/rss_raytracing_evaluation_framework/InputData/Frankfurt/5355566", "/home/panwei/Desktop/bachelorThese/data/rss_raytracing_evaluation_framework/InputData/Frankfurt/3625286", "/home/panwei/Desktop/bachelorThese/data/rss_raytracing_evaluation_framework/InputData/Frankfurt/5194466", "/home/panwei/Desktop/bachelorThese/data/rss_raytracing_evaluation_framework/InputData/Frankfurt/3085836", "/home/panwei/Desktop/bachelorThese/data/rss_raytracing_evaluation_framework/InputData/Frankfurt/3725856", "/home/panwei/Desktop/bachelorThese/data/rss_raytracing_evaluation_framework/InputData/Frankfurt/5152246", "/home/panwei/Desktop/bachelorThese/data/rss_raytracing_evaluation_framework/InputData/Frankfurt/2245166", "/home/panwei/Desktop/bachelorThese/data/rss_raytracing_evaluation_framework/InputData/Frankfurt/4434416", "/home/panwei/Desktop/bachelorThese/data/rss_raytracing_evaluation_framework/InputData/Frankfurt/3251696", "/home/panwei/Desktop/bachelorThese/data/rss_raytracing_evaluation_framework/InputData/Frankfurt/5942206", "/home/panwei/Desktop/bachelorThese/data/rss_raytracing_evaluation_framework/InputData/Frankfurt/3821656", "/home/panwei/Desktop/bachelorThese/data/rss_raytracing_evaluation_framework/InputData/Frankfurt/3054566", "/home/panwei/Desktop/bachelorThese/data/rss_raytracing_evaluation_framework/InputData/Frankfurt/4063746", "/home/panwei/Desktop/bachelorThese/data/rss_raytracing_evaluation_framework/InputData/Frankfurt/5853136"])
    # subDir = getSubDir(my_path)
    # # to print string in double quotes (essential for the running of matlab code)
    # print(json.dumps(subDir))

    # #get folder name
    dataFolderPath = "/home/panwei/Desktop/bachelorThese/myCode/data/matlabData/frankfurt/hornAntennaAggregate/"
    # dataFolderPath = "/home/panwei/Desktop/bachelorThese/myCode/data/matlabData/frankfurt/hornAntenna/"
    # dataFolderPath = "/home/panwei/Desktop/bachelorThese/myCode/data/matlabData/seoul/phasedAntennaArrayAggregate/"
    # dataFolderPath = "/home/panwei/Desktop/bachelorThese/myCode/data/matlabData/seoul/phasedAntennaArray/"
    # dataFolderPath = "/home/panwei/Desktop/bachelorThese/myCode/data/matlabData/seoul/hornAntenna/"



    dataDir = getSubDir(dataFolderPath, ".csv", mode = 'relative')
    for dataDirName in dataDir :
        print "\"%s\"" % (dataDirName),
    # outStringRaw = json.dumps(dataDir)
    # outString = ' '.join(n for n in outStringRaw)
    # print(outString)
    # print(json.dumps(dataDir))
    # print("files:")
    # files = getFilesFromDirs(dataDir)
    # print(files)
