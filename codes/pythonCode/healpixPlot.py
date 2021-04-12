from mpl_toolkits import mplot3d
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
# matplotlib.use('AGG')
import numpy as np
import healpy as hp
import pandas as pd
import re
import json
from almCoeff import plotDiffCl
from almCoeff import plotCl
from almCoeff import plotMultiCls
from almCoeff import writeAngularPower

from dataPreProc import *
from plotHelpFunc import *
import os
import datetime
 
currentDT = datetime.datetime.now()
currentTime = currentDT.strftime("%H:%M:%S")

# dataDictSavedPath = "/home/panwei/Desktop/bachelorThese/myCode/pythonCode/dataDict.txt"
dataDictSavedPath = "/home/panwei/Desktop/bachelorThese/myCode/pythonCode/dataDict.json"


LOGLOGPLOT = False
SEMILOGYPLOT = False
SEMILOGXPLOT = False
HEALPIX_CORR = False
HEALPIX_AVERAGE = False
HEALPIX_ALL_ONE = False
COMPACT_TEST = False
SIGNAL_DATA = False
SPHERICAL = False
VEC_MAP = False

MY_DOUBLE = False

MATERN = False
PRINT_POISSON = False

NO_TITLE = True

PLOT_RECENT = False
CREATE_DIR = True
MULTI_RUN = True

RANDOM_MAP_INFO = True
if (HEALPIX_CORR == True):
    SIGNAL_DATA = False
if (SIGNAL_DATA == True):
    MULTI_RUN = False
OFFSET = False
ADD_TIME = False
pointType = 'b*'
#number of major ticks in a plot in the xaxis
SPHERE_MAJOR_TICKS_NUM = 1
MAJOR_TICKS_NUM = 4

datum = "7/5/"

NO_PLOT = True


LARGE_Y_RANGE = True
Y_MAX = 4
Y_MIN = -4
POS_TITLE = 0.98
plotNumThreshold = 20
plotted = 0
figureI = 0
# within the COMPACT_TEST
ALL_ONE_PLOT = False
HEAL_AVG_PLOT = True


# if set will output the saved plotted info
VERBOSE = True

# if HEALPIX_AVERAGE == True:
#     datum = "5/12/vecMod/healpixAverage/"
# elif HEALPIX_ALL_ONE == True:
#     datum += "AllOneMap/"
# # elif MATERN == True:
# #     if OFFSET == True:
# #         # datum = "5/12/vecMod/maternYesOffset/"
# #     else:
# #         # datum = "5/12/vecMod/maternNoOffset/"
# elif SIGNAL_DATA == True:
#     datum = "5/14/vecMod/dataTest/"
# # datum = "5/12/vecMod/healpixAverage/"
# datum = "4/13/healpixCorr/vecMod/"

rotateAngle = 90  # the angle that ylabel rotates
plotTogether = True

class FileNotFoundError(OSError):
    pass

def getItems(dataDict, config, threshold, \
    model = 'rayTrace', queries = ['correlation']):
    """
        given the dataDict structure, get the info wanted

    Parameters:

    config: str
        e.g. TX1_RX2
    
    threshold: str
        e.t. '-75'
    
    model: str
        e.g. 'rayTrace'
        
    queries: str
        e.g. 'correlation', 'healpixMap',
    """
    thresholdDict = dataDict[config]
    calcuPath = thresholdDict[threshold]
    fileOrigName = calcuPath[model]
    responses = {}
    dirName = os.path.dirname(fileOrigName)
    if 'correlation' in queries:
        responses['correlation'] = os.path.join(dirName, "SignalAverageHealpixMode4.txt")
    if 'healpixMap' in queries:
        responses['healpixMap'] = os.path.join(dirName, "SignalAverageHealpixMap.txt")
    if 'sphericalPoint' in queries:
        responses['sphericalPoint'] = os.path.join(dirName, "spherePoint.txt")
    
    infoName = model + "_info"
    responses['calcuInfo'] = calcuPath[infoName]
    return responses
    # get the parent dir name


def plot2Func(filePath1, filePath2, legend1, legend2,savePath):
    global figureI,plotted
    X1,Y1,X2,Y2 = [],[],[],[]
    for line in open(filePath1, 'r'):
    #   print(line)
        values = [float(s) for s in line.split()]
    #   print(values)
        if(len(values) >= 1):
            X1.append(values[0])
            Y1.append(values[1])
            
    # fig, ax = plt.subplots(num = figNum)
    X1 = np.array(X1,dtype= float)
    Y1 = np.array(Y1,dtype= float)
    
    if (Convert2Deg == True):
        X1 = X1 * 180 / np.pi
    # X1, Y1 = cleanData(X,Y)
    # if (CUT_RANGE == True):
    #     X1cut,Y1cut = cutData(X1,Y1)    
    for line in open(filePath2, 'r'):
    #   print(line)
        values = [float(s) for s in line.split()]
    #   print(values)
        if(len(values) >= 1):
            X2.append(values[0])
            Y2.append(values[1])
            
    # # fig, ax = plt.subplots(num = figNum)
    X2 = np.array(X2,dtype= float)
    Y2 = np.array(Y2,dtype= float)    

    if (Convert2Deg == True):
        X2 = X2 * 180 / np.pi
    fig,ax = plt.subplots(num = figureI)
    ax.plot(X1,Y1,label = legend1)
    ax.plot(X2,Y2,label = legend2)
    plt.title("LS estimator comparison with 500 pts")
    chartBox = ax.get_position()
    ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.8, chartBox.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), shadow=True, ncol=1)
    if (Convert2Deg == True):
        plt.xlabel(r'$\theta\:(degree)$',fontsize = 10)
    else:
        plt.xlabel(r'$\theta$',fontsize = 10)    
    plt.ylabel(r'$\omega(\theta)$',fontsize = 10,rotation = rotateAngle)    

    plt.savefig(savePath)


def plotMultiFunc(filePaths,legends,title,savePath):
    global figureI,plotted,cut_Val
    if (Convert2Deg == True):
        # 1.59 or 0.7
        cut_Val = 1.59 * 180/np.pi
    else:
        cut_Val = 0.7
    X = [None] * len(filePaths)
    Y = [None] * len(filePaths)
    count = 0
    for filePath in filePaths:
        v1Temp = []
        v2Temp = []
        for line in open(filePath, 'r'):
        #   print(line)
            values = [float(s) for s in line.split()]
        #   print(values)
            if(len(values) >= 1):
                v1Temp.append(values[0])
                v2Temp.append(values[1])
        X[count] = np.array(v1Temp)
        if (Convert2Deg == True):
            X[count] = X[count] * 180 / np.pi        
        Y[count] = np.array(v2Temp)
        count += 1
    fig, ax = plt.subplots(num = figureI)
    

    for i in range(count):
        X[i],Y[i] = cleanData(X[i],Y[i])
        X[i],Y[i] = cutData(X[i],Y[i])

    xMax = -100
    yMax = -100
    yMin = 100
    for i in range(count):
        if(xMax < max(X[i])):
            xMax = max(X[i])
        if(yMax < max(Y[i])):
            yMax = max(Y[i])
        if(yMin > min(Y[i])):
            yMin = min(Y[i])
    yMax = round(yMax) + 1
    yMin = round(yMin) - 0.1
    # X1 = np.array(X1,dtype= float)
    # Y1 = np.array(Y1,dtype= float)

    # X1, Y1 = cleanData(X,Y)
    # if (CUT_RANGE == True):
    #     X1cut,Y1cut = cutData(X1,Y1)    
    
            
    # # fig, ax = plt.subplots(num = figNum)
    # X2 = np.array(X2,dtype= float)
    # Y2 = np.array(Y1,dtype= float)    
    fig,ax = plt.subplots(num = figureI+10)
    # num = 210
    for i in range(count):
        # ax.plot(X[i][:num],Y[i][:num],label = legends[i])
        ax.plot(X[i],Y[i],label = legends[i])
        plt.savefig(savePath)

    ax.xaxis.set_major_locator(plt.MultipleLocator(xMax/MAJOR_TICKS_NUM))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(xMax/(MAJOR_TICKS_NUM * 2)))
    # ax.yaxis.set_major_locator(plt.MultipleLocator((yMax - yMin)/MAJOR_TICKS_NUM))
    # plt.xticks(np.arange(0,max(X),sxtep = 0.1))
    if (Convert2Deg == True):
        plt.xlabel(r'$\theta\:(degree)$',fontsize = 10)
    else:
        plt.xlabel(r'$\theta$',fontsize = 10)    
    plt.ylabel(r'$\omega(\theta)$',fontsize = 10,rotation = rotateAngle)
    plt.title(title)
    chartBox = ax.get_position()
    ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.8, chartBox.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), shadow=True, ncol=1)
    plt.savefig(savePath)

# def cutData(X,Y,Y2 = []):
#     global cut_Val
#     Y = Y[(X < cut_Val)]
#     if (len(Y2) >= 1):
#         Y2 = Y2[(X < cut_Val)]    
#     X = X[(X < cut_Val)]
#     return X,Y,Y2
# plot the data calculated by 4 estimators run by healpix Average
def plotEstimatorHealpixAverageAlgo(fileNamePrefix, namePrefix, 
    figureSaveNamePre, pointCount, 
    randomCount, deltaTheta, average, thetaBegin,thetaEnd, 
    figNums, compactPlot = True, **kwargs):
    '''
    fileNamePrefix: path to that directory
    namePrefix: prefix for the data
    figureSaveNamePre: prefix for the to be saved plotted figure
    '''
    global plotted,figureI
    peebles1 = fileNamePrefix + namePrefix + "1.txt"
    peebles2 = fileNamePrefix + namePrefix + "2.txt"
    hamilton = fileNamePrefix + namePrefix + "3.txt"
    landy = fileNamePrefix + namePrefix + "4.txt"


    additionName = ""
    titleAddition = ""
    count = 0
    for key,value in kwargs.items():
        additionName += key + ":" + str(value)
        titleAddition += key + ":" + str(value) + ", "
        count += 1
        if (count % 4 == 0) :
            titleAddition+= "\n"

    adjustTitlePos = False
    if (count >= 4):
        adjustTitlePos = True

    if (Convert2Deg == True):
        thetaBegin *= 180 / np.pi
        thetaEnd *= 180/np.pi
        deltaTheta *= 180/np.pi
    corrName = "thetaBegin:%.2f"%(thetaBegin)+"thetaEnd:%.2f"%(thetaEnd) + "points:" + str(pointCount) + "Random:" + \
            str(randomCount) + "Average:" + str(average) + " dT: %.4f" %(deltaTheta) + \
                additionName
    if LOGLOGPLOT == True:
        corrName += "loglogPlot"
    elif SEMILOGYPLOT == True:
        corrName += "semilogYPlot"
    elif SEMILOGXPLOT == True:
        corrName += "semilogXPlot"
    # figEstPre = figureSaveNamePre + namePrefix + corrName
    figEstPre = figureSaveNamePre + namePrefix
    figPeebles1 = figEstPre + "Peebles1.png"
    figPeebles2 = figEstPre+ "Peebles2.png"
    figHamil = figEstPre + "Hamilton.png"
    figLandy = figEstPre + "Landy.png"
    if (MATERN == False):
        title =  "thetaBegin:%.2f"%(thetaBegin) + " ," + "thetaEnd:%.2f, "%(thetaEnd)+ str(pointCount) +" data points, random " + \
            str(randomCount) + "points\n" + str(average) + " times average " + " dT: %.4f" %(deltaTheta) + \
                "\n" + titleAddition
    else :
        title =  "thetaBegin:%.2f"%(thetaBegin) + " ," + "thetaEnd:%.2f, "%(thetaEnd) + \
            str(randomCount) + "points\n" + str(average) + " times average " + " dT: %.4f" %(deltaTheta) + \
                "\n" + titleAddition        
    if LOGLOGPLOT == True:
        title += " loglogPlot"
    elif SEMILOGYPLOT == True:
        title += " semilogYPlot"
    elif SEMILOGXPLOT == True:
        title += " semilogXPlot"

    # title += '\n' +namePrefix
    title += "\n"
    titlePeebles1 = title + " PH"
    titlePeebles2 = title + " DP"
    titleHamilton = title + " Hamilton"
    titleLandy = title + " LS"

    fileNames = [peebles1, peebles2, hamilton,landy]
    figNames = [figPeebles1, figPeebles2, figHamil, figLandy]
    titleNames = [titlePeebles1, titlePeebles2, titleHamilton, titleLandy]
    counts = range(4)
    yvalues = [None] * 8
    xvalues = [None] * 8
    # xFinal = []
    for fileName, saveName, titleName ,figNum, count in zip(fileNames, figNames, titleNames,figNums,counts):
        fig, ax = plt.subplots(num = figNum)
        if (CUT_RANGE == True):
            xvalues[count],yvalues[count],xvalues[count + 4],yvalues[count+4] = printEstimator(fileName,titleName,saveName,figNum,ax,False,"",adjustTitlePos)
        else:
            xvalues[count], yvalues[count] = printEstimator(fileName,titleName,saveName,figNum,ax,False,"",adjustTitlePos)
        figureI += 1
        plotted += 1
        # if (count >= 1):
        #     xFinal = np.intersect1d(xvalues[count],xvalues[count - 1])
    if (CUT_RANGE == False):
        xvalues = xvalues[0:4]
        yvalues = yvalues[0:4]
    # i = 1
    # j = 2
    # fig, ax = plt.subplots(num = j)
    # printEstimator(fileNames[i],figNames[i],titleNames[i],j,ax,False,'')
    if(compactPlot == True):
        legendName = ["PH","DP", "Hamilton", "LS"]
        # plt.figure()
        # ax = plt.subplot(111)
            # plt.figure(figNum)
        fignum = figNums[-1] + 1
        fig, ax = plt.subplots(num = fignum)
        # fig, ax = plt.subplots()
        try:
            xvalues
        except NameError:
            for legend in legendName:
                printEstimator(fileName,title, figEstPre + "compactPlot.png",fignum,ax,True,legend,adjustTitlePos)
        else:
            sweepCompactPlot(xvalues,yvalues,title,figEstPre + "compactPlot.png", fignum, ax,legendName,adjustTitlePos)
        plotted += 1
        figureI += 1
        

def plotEstimatorMaternAlgo(fileNamePrefix, namePrefix, figureSaveNamePre, pointCount, randomCount, deltaTheta, average, spannedAngle, figNums, clusterCenterCount, 
        clusterTheta, clusterPointCount, totalPoint, thetaBegin,thetaEnd,dBegin, dEnd, compactPlot = False, useMathDouble = False, useOffsetPoisson = False):
    """
    Parameters:
    --------
    fileNamePrefix: str
        pre prefix for all data file ususally the path like ../data/6/04/matern/
    
    namePrefix: str
        prefix for the data file like "poissonHealpixMode"

    figureSaveNamePre: str
        prefix fot the saved figure, like ../../photodata/6/04/matern/

    """
    global plotted,figureI        
    peebles1 = fileNamePrefix + namePrefix + "1.txt"
    peebles2 = fileNamePrefix + namePrefix + "2.txt"
    hamilton = fileNamePrefix + namePrefix + "3.txt"
    landy = fileNamePrefix + namePrefix + "4.txt"
    if (Convert2Deg == True):
        thetaBegin *= 180 / np.pi
        thetaEnd *= 180/np.pi
        deltaTheta *= 180/np.pi
        clusterTheta *= 180/np.pi

    corrName = "spanned:%.2f"%(spannedAngle) + "points:" + str(pointCount) + "Random:" + str(randomCount) + "Average:" + str(average) + " dT: %.4f" %(deltaTheta)+ \
         "clusterCenter:" + str(clusterCenterCount) + "clusterTheta:" + str(clusterTheta) + "clusterPointCount:" + str(clusterCenterCount) + \
        "totalPoint:" + str(totalPoint) + "thetaBegin: " + str(thetaBegin) + "thetaEnd: "+ \
        str(thetaEnd) + "daughterThetaBegin:" + str(dBegin) + "daughterThetaEnd:" + str(dEnd)+ \
        "double:" + str(useMathDouble) + "offsetPoisson:" + str(useOffsetPoisson)

    corrNameTitle = "spanned:%.2f"%(spannedAngle) + "points:" + str(pointCount) + "Random:" + str(randomCount) + "Average:" + str(average) + " dT: %.4f" %(deltaTheta) \
        + "clusterCenter:" + str(clusterCenterCount) + "clusterTheta:" + str(clusterTheta) + "clusterPointCount:" + str(clusterCenterCount) + \
        "totalPoint:" + str(totalPoint) + "\nthetaBegin: " + str(thetaBegin) + "thetaEnd: " + \
        str(thetaEnd) + "daughterThetaBegin:" + str(dBegin) + "daughterThetaEnd:" + str(dEnd)
        # "double:" + str(useMathDouble) + "offsetPoisson:" + str(useOffsetPoisson)
    if LOGLOGPLOT == True:
        corrName += "loglogPlot"
    elif SEMILOGYPLOT == True:
        corrName += "semilogYPlot"
    elif SEMILOGXPLOT == True:
        corrName += "semilogXPlot"
    # figEstPre = figureSaveNamePre + namePrefix + corrName
    figEstPre = figureSaveNamePre + "useMath:" + str(useMathDouble) + "offsetPoisson:" + str(useOffsetPoisson)
    figPeebles1 = figEstPre + "Peebles1.png"
    figPeebles2 = figEstPre+ "Peebles2.png"
    figHamil = figEstPre + "Hamilton.png"
    figLandy = figEstPre + "Landy.png"
    title =  "spanned:%.2f"%(spannedAngle) + " ," + str(pointCount) +" poisson points, random " + str(randomCount) + "points\n" + str(average) + " times average " + " dT: %.4f" %(deltaTheta) \
        +"\nclusterCenter:" + str(clusterCenterCount) + ",thetaDist:" + str(clusterTheta) + ",points in cluster:" + str(clusterPointCount) + ",totalPoint:" + str(totalPoint) + "\nthetaBegin: " + str(thetaBegin) + "thetaEnd: " + str(thetaEnd) +\
            "clusterTheta:" + str(dEnd - dBegin)
            # str(useMathDouble) + "offsetPoisson:" + str(useOffsetPoisson)
    if LOGLOGPLOT == True:
        title += " loglogPlot"
    elif SEMILOGYPLOT == True:
        title += " semilogYPlot"
    elif SEMILOGXPLOT == True:
        title += " semilogXPlot"
    # title += '\n' +namePrefix
    titlePeebles1 = title + " PH"
    titlePeebles2 = title + " DP"
    titleHamilton = title + " Hamilton"
    titleLandy = title + " LS"

    fileNames = [peebles1, peebles2, hamilton,landy]
    figNames = [figPeebles1, figPeebles2, figHamil, figLandy]
    titleNames = [titlePeebles1, titlePeebles2, titleHamilton, titleLandy]
    counts = range(4)
    yvalues = [None] * 4
    
    for fileName, saveName, titleName ,figNum, count in zip(fileNames, figNames, titleNames,figNums,counts):
        fig, ax = plt.subplots(num = figNum)
        xvalues, yvalues[count] = printEstimator(fileName,titleName,saveName,figNum,ax)
        plotted += 1
        figureI+= 1
    if(compactPlot == True):
        legendName = ["Peebles","DP", "Hamilton", "LS"]
        plt.figure()
        fignum = figNums[-1] + 1
        fig, ax = plt.subplots(num = fignum)
        try:
            xvalues
        except NameError:
            for legend in legendName:
                printEstimator(fileName,title, figEstPre + "compactPlot.png",fignum,ax,True,legend)
        else:
            sweepCompactPlot(xvalues,yvalues,title,figEstPre + "compactPlot.png", fignum, ax,legendName)
        plotted += 1
        figureI += 1

def printEstimator(fileName, figTitle, saveName,figNum, ax,
    plotLegend = False, legend = "",adjustTitlePos = False):
    '''
    fileName: name for the data file
    figTitle: title of the figure
    saveName: name for the saved figure
    '''
    global figureI, plotted, cut_Val    
    X, Y , Xcut, Ycut,Y2,Y2cut= [], [],[],[],[],[]
    # try:

    for line in open(fileName, 'r'):
    #   print(line)
        values = [float(s) for s in line.split()]
    #   print(values)
        if(len(values) >= 1):
            X.append(values[0])
            Y.append(values[1])
            if(len(values) >= 3):
                Y2.append(values[2])
            
    # fig, ax = plt.subplots(num = figNum)
    X = np.array(X,dtype= float)
    if (Convert2Deg == True):
        X = X * 180 / np.pi
    Y = np.array(Y,dtype= float)
    Y2 = np.array(Y2,dtype=float)

    X, Y = cleanData(X,Y)
    if (CUT_RANGE == True):
        Xcut,Ycut = cutData(X,Y)
    # X , Y, Y2= cleanData(X,Y,Y2)
    # if (CUT_RANGE == True):
    #     Xcut,Ycut,Y2cut = cutData(X,Y,Y2)
    if (len(Y2) <= 1):
        if (CUT_RANGE == True):
            figureI += 2
            figTemp, axTemp = plt.subplots(num = figureI)
            if (len(Xcut)>= 1):
                xcutMax = round(max(Xcut),1)  
                if(plotLegend == False):

                    if LOGLOGPLOT == True:
                        axTemp.loglog(Xcut,Ycut,pointType)
                    elif SEMILOGYPLOT == True:
                        axTemp.semilogy(Xcut,Ycut,pointType)
                    elif SEMILOGXPLOT == True:
                        axTemp.semilogx(Xcut,Ycut,pointType)
                        insertSemiX = saveName.find("png")
                        saveName = saveName[:insertSemiX-1]+ "SEMIX"+saveName[insertSemiX-1:]                        
                    else:
                        axTemp.plot(Xcut, Ycut)
                else :
                    if LOGLOGPLOT == True:
                        axTemp.loglog(Xcut,Ycut,pointType)
                    elif SEMILOGYPLOT == True:
                        axTemp.semilogy(Xcut,Ycut,pointType)
                        insertSemiX = saveName.find("png")
                        saveName = saveName[:insertSemiX-1]+ "SEMIY"+saveName[insertSemiX-1:]
                    elif SEMILOGXPLOT == True:
                        axTemp.semilogx(Xcut,Ycut,pointType)
                        insertSemiX = saveName.find("png")
                        saveName = saveName[:insertSemiX-1]+ "SEMIX"+saveName[insertSemiX-1:]                        
                    else:
                        axTemp.plot(Xcut, Ycut,label = legend)       
                plotted += 1 
                figureI += 1
                axTemp.xaxis.set_major_locator(plt.MultipleLocator(xcutMax/MAJOR_TICKS_NUM))
                axTemp.xaxis.set_minor_locator(plt.MultipleLocator(xcutMax/(MAJOR_TICKS_NUM * 2))) 
                insertPos = saveName.find("png")
                saveTemp = saveName[:insertPos-1] + "CutRange" + saveName[insertPos-1:]
                if (Convert2Deg == True):
                    plt.xlabel(r'$\theta\:(degree)$',fontsize = 10)
                else:
                    plt.xlabel(r'$\theta$',fontsize = 10)
                plt.ylabel(r'$\omega(\theta)$',fontsize = 10,rotation = rotateAngle)
                if (adjustTitlePos == True):
                    plt.title(figTitle, fontsize = 7, y = POS_TITLE)
                else:
                    plt.title(figTitle,fontsize = 8)        
                plt.savefig(saveTemp)    
                YcutMax = max(max(Ycut),Y_MAX) + 1
                YcutMin = min(min(Ycut),Y_MIN) - 1     
                if (LARGE_Y_RANGE == True):
                    bottom,top = plt.ylim()
                    plt.ylim(YcutMin, YcutMax)
                    insertPos = saveTemp.find("png")
                    saveTempLarge = saveTemp[:insertPos-1] + "largeScale" + saveTemp[insertPos-1:]
                    plt.savefig(saveTempLarge)
                    plotted += 1
                    plt.ylim(bottom,top)   
            else:
                return [],[],[],[]    
        # plot the original data     
        if (len(X) >= 1):
            xMax = round(max(X),1)      
            figNew, axNew = plt.subplots(num = figureI + 1)
            figureI += 1
            if(plotLegend == False):
                if LOGLOGPLOT == True:
                    axNew.loglog(X,Y,pointType)
                elif SEMILOGYPLOT == True:
                    axNew.semilogy(X,Y,pointType)
                elif SEMILOGXPLOT == True:
                    axNew.semilogx(X,Y,pointType)
                else:
                    axNew.plot(X, Y)
            else :
                if LOGLOGPLOT == True:
                    axNew.loglog(X,Y,pointType)
                elif SEMILOGYPLOT == True:
                    axNew.semilogy(X,Y,pointType)
                elif SEMILOGXPLOT == True:
                    axNew.semilogx(X,Y,pointType)
                else:
                    axNew.plot(X, Y,label = legend)
            plotted += 1
            axNew.xaxis.set_major_locator(plt.MultipleLocator(xMax/MAJOR_TICKS_NUM))
            axNew.xaxis.set_minor_locator(plt.MultipleLocator(xMax/(MAJOR_TICKS_NUM * 2)))
            # plt.xticks(np.arange(0,max(X),sxtep = 0.1))
            if (Convert2Deg == True):
                plt.xlabel(r'$\theta\:(degree)$',fontsize = 10)
            else:
                plt.xlabel(r'$\theta$',fontsize = 10)

            plt.ylabel(r'$\omega(\theta)$',fontsize = 10,rotation = rotateAngle)

            if (adjustTitlePos == True):
                plt.title(figTitle, fontsize = 7, y = POS_TITLE)
            else:
                plt.title(figTitle,fontsize = 8)      
            # plt.title(figTitle,fontsize = 8)


            YMax = max(max(Y),Y_MAX) + 1
            YMin = min(min(Y),Y_MIN) - 1
            if (LARGE_Y_RANGE == True):
                bottom,top = plt.ylim()
                plt.ylim(YMin, YMax)
                insertPos = saveName.find("png")
                saveTemp = saveName[:insertPos-1] + "largeScale" + saveName[insertPos-1:]
                plt.savefig(saveTemp)
                plt.ylim(bottom,top)
                plotted += 1
            # plt.yticks(np.arange(min(Y)-0.1, max(Y)+0.1, 0.05))
            # plt.yticks(np.arange(round(min(Y) - 0.1,1), round(max(Y) + 0.1,1), 0.05))
            # title = "Pi healpix Peebles with data points: " + str(hpPointCount)
            
            # if(plotLegend == False):
            #     plt.tight_layout()
            plt.savefig(saveName)   #  xMax = round(max(X),1)      
            # except IOError:
            #     print("fileName starts with: " + fileName + " not Found")
            if(NO_PLOT == True):
                plt.close('all')
            if (CUT_RANGE == True):
                return X,Y,Xcut,Ycut
            else:
                return X,Y
        else:
            return [],[]
    else:
        if (CUT_RANGE == True):
            figureI += 2
            figTemp, axTemp = plt.subplots(num = figureI)
            figureI += 1
            xcutMax = round(max(Xcut),1)        
            if(plotLegend == False):

                if LOGLOGPLOT == True:
                    axTemp.loglog(Xcut,Ycut,pointType)
                    axTemp.loglog(Xcut,Y2cut,pointType)
                elif SEMILOGYPLOT == True:
                    axTemp.semilogy(Xcut,Ycut,pointType)
                    axTemp.semilogy(Xcut,Y2cut,pointType)
                elif SEMILOGXPLOT == True:
                    axTemp.semilogx(Xcut,Ycut,pointType)
                    axTemp.semilogx(Xcut,Y2cut,pointType)
                    insertSemiX = saveName.find("png")
                    saveName = saveName[:insertSemiX-1]+ "SEMIX"+saveName[insertSemiX-1:]                        
                else:
                    axTemp.plot(Xcut, Ycut, label = "original data")
                    axTemp.plot(Xcut, Y2cut, label = "poisson")
            else :
                if LOGLOGPLOT == True:
                    axTemp.loglog(Xcut,Ycut,pointType)
                    axTemp.loglog(Xcut,Y2cut,pointType)
                elif SEMILOGYPLOT == True:
                    axTemp.semilogy(Xcut,Ycut,pointType)
                    axTemp.semilogy(Xcut,Y2cut,pointType)
                elif SEMILOGXPLOT == True:
                    axTemp.semilogx(Xcut,Ycut,pointType)
                    axTemp.semilogx(Xcut,Y2cut,pointType)
                    insertSemiX = saveName.find("png")
                    saveName = saveName[:insertSemiX-1]+ "SEMIX"+saveName[insertSemiX-1:]                        
                else:
                    axTemp.plot(Xcut, Ycut,label = legend + "original data")   
                    axTemp.plot(Xcut, Y2cut,'--',label = legend + "poisson")    
            
            axTemp.xaxis.set_major_locator(plt.MultipleLocator(xcutMax/MAJOR_TICKS_NUM))
            axTemp.xaxis.set_minor_locator(plt.MultipleLocator(xcutMax/(MAJOR_TICKS_NUM * 2))) 
            insertPos = saveName.find("png")
            saveTemp = saveName[:insertPos-1] + "CutRange" + saveName[insertPos-1:]
            if (Convert2Deg == True):
                plt.xlabel(r'$\theta\:(degree)$',fontsize = 10)
            else:
                plt.xlabel(r'$\theta$',fontsize = 10)

            plt.ylabel(r'$\omega(\theta)$',fontsize = 10,rotation = rotateAngle)

            if (figTitle[-1] == "\n"):
                figTitle = figTitle[:-1]
            if (figTitle[-2] == ","):
                figTitle = figTitle[:-2]            

            if (adjustTitlePos == True):
                plt.title(figTitle, fontsize = 7, y = POS_TITLE)
            else:
                plt.title(figTitle,fontsize = 8)        
            plt.savefig(saveTemp)    
            plotted += 1
            YcutMax = max(maxVal(Ycut,Y2cut),Y_MAX)
            YcutMin = min(minVal(Ycut,Y2cut),Y_MIN)         
            if (LARGE_Y_RANGE == True):
                bottom,top = plt.ylim()
                plt.ylim(YcutMin, YcutMax)
                insertPos = saveTemp.find("png")
                saveTempLarge = saveTemp[:insertPos-1] + "largeScale" + saveTemp[insertPos-1:]
                plt.savefig(saveTempLarge)
                plotted += 1
                plt.ylim(bottom,top)   
                
        # plot the original data      
        xMax = round(max(X),1)      
        figNew, axNew = plt.subplots(num = figureI + 1)
        figureI += 1
        if(plotLegend == False):
            if LOGLOGPLOT == True:
                axNew.loglog(X,Y,pointType)
                axNew.loglog(X,Y2,'--',pointType)
            elif SEMILOGYPLOT == True:
                axNew.semilogy(X,Y,pointType)
                axNew.semilogy(X,Y2,'--',pointType)
            elif SEMILOGXPLOT == True:
                axNew.semilogx(X,Y,pointType)
                axNew.semilogx(X,Y2,'--',pointType)
            else:
                axNew.plot(X, Y)
                axNew.plot(X, Y2)
        else :
            if LOGLOGPLOT == True:
                axNew.loglog(X,Y,pointType)
                axNew.loglog(X,Y2,pointType)
            elif SEMILOGYPLOT == True:
                axNew.semilogy(X,Y,pointType)
                axNew.semilogy(X,Y2,pointType)
            elif SEMILOGXPLOT == True:
                axNew.semilogx(X,Y,pointType)
                axNew.semilogx(X,Y2,pointType)
            else:
                axNew.plot(X, Y,label = legend)
                axNew.plot(X, Y2,'--',label = legend)
        axNew.xaxis.set_major_locator(plt.MultipleLocator(xMax/MAJOR_TICKS_NUM))
        axNew.xaxis.set_minor_locator(plt.MultipleLocator(xMax/(MAJOR_TICKS_NUM * 2)))
        # plt.xticks(np.arange(0,max(X),sxtep = 0.1))
        if (Convert2Deg == True):
            plt.xlabel(r'$\theta\:(degree)$',fontsize = 10)
        else:
            plt.xlabel(r'$\theta$',fontsize = 10)

        plt.ylabel(r'$\omega(\theta)$',fontsize = 10,rotation = rotateAngle)

        if (adjustTitlePos == True):
            plt.title(figTitle, fontsize = 7, y = POS_TITLE)
        else:
            plt.title(figTitle,fontsize = 8)      
        # plt.title(figTitle,fontsize = 8)


        YMax = max(maxVal(Y,Y2),Y_MAX)
        YMin = min(minVal(Y,Y2),Y_MIN)
        if (LARGE_Y_RANGE == True):
            bottom,top = plt.ylim()
            plt.ylim(YMin, YMax)
            insertPos = saveName.find("png")
            saveTemp = saveName[:insertPos-1] + "largeScale" + saveName[insertPos-1:]
            plt.savefig(saveTemp)
            plotted += 1
            plt.ylim(bottom,top)

        # plt.yticks(np.arange(min(Y)-0.1, max(Y)+0.1, 0.05))
        # plt.yticks(np.arange(round(min(Y) - 0.1,1), round(max(Y) + 0.1,1), 0.05))
        # title = "Pi healpix Peebles with data points: " + str(hpPointCount)
        
        # if(plotLegend == False):
        #     plt.tight_layout()
        plt.savefig(saveName)
        plotted += 1
        # except IOError:
        #     print("fileName starts with: " + fileName + " not Found")
        if(NO_PLOT == True):
            plt.close('all')
        if (CUT_RANGE == True):
            return X,Y,Xcut,Ycut
        else:
            return X,Y


    # if(plotLegend == False):

    #     if LOGLOGPLOT == True:
    #         ax.loglog(X,Y,pointType)
    #     elif SEMILOGYPLOT == True:
    #         ax.semilogy(X,Y,pointType)
    #     elif SEMILOGXPLOT == True:
    #         ax.semilogx(X,Y,pointType)
    #     else:
    #         ax.plot(X, Y)
    # else :
    #     if LOGLOGPLOT == True:
    #         ax.loglog(X,Y,pointType)
    #     elif SEMILOGYPLOT == True:
    #         ax.semilogy(X,Y,pointType)
    #     elif SEMILOGXPLOT == True:
    #         ax.semilogx(X,Y,pointType)
    #     else:
    #         ax.plot(X, Y,label = legend)
    # ax.xaxis.set_major_locator(plt.MultipleLocator(xMax/MAJOR_TICKS_NUM))
    # ax.xaxis.set_minor_locator(plt.MultipleLocator(xMax/(MAJOR_TICKS_NUM * 2)))
    # # plt.xticks(np.arange(0,max(X),sxtep = 0.1))

    # plt.xlabel(r'$\theta$',fontsize = 10)
    # plt.ylabel(r'$\omega(\theta)$',fontsize = 10,rotation = rotateAngle)
    # plt.title(figTitle,fontsize = 8)


    # YMax = max(max(Y),Y_MAX)
    # YMin = min(min(Y),Y_MIN)
    # if (LARGE_Y_RANGE == True):
    #     bottom,top = plt.ylim()
    #     plt.ylim(YMin, YMax)
    #     insertPos = saveName.find("png")
    #     saveTemp = saveName[:insertPos-1] + "largeScale" + saveName[insertPos-1:]
    #     plt.savefig(saveTemp)
    #     plt.ylim(bottom,top)

    # # plt.yticks(np.arange(min(Y)-0.1, max(Y)+0.1, 0.05))
    # # plt.yticks(np.arange(round(min(Y) - 0.1,1), round(max(Y) + 0.1,1), 0.05))
    # # title = "Pi healpix Peebles with data points: " + str(hpPointCount)
    
    # # if(plotLegend == False):
    # #     plt.tight_layout()
    # plt.savefig(saveName)
    # # except IOError:
    # #     print("fileName starts with: " + fileName + " not Found")
    # if (CUT_RANGE == True):
    #     return X,Y,Xcut,Ycut
    # else:
    #     return X,Y
def multiTimeEstimatorPlot(fileNamePathPre, dataFileSuffix,figTitlePre,saveNamePre,ax, plotLegend = False,
    legend = "",adjustTitlePos = False) :
    """ plot the estimator function given multiple data file
    Paramters:
    ----------
    fileNamePathPre: list
        prefix for the string of the data file, e.g ["../data/result0","../data/result1"]
    
    numRange: list
        range of the sweeping number for the data file, e.g range(0,5), so 
        like result0.txt, result1.txt ...

    dataFileSuffix: str
        suffix of the data file, like "Pi_HealpixRandomMode4.txt"

    figTitlePre: str
        prefix for the title of the plot
    
    saveNamePre: str
        prefix for the dir to the saved figure, like "../photodata/landy"

    ax: obj
        axis of the plot
    
    """
    global figureI, plotted, cut_Val    
    setAxisRange = False
    xMax = 0
    yMax = -10000
    yMin = 10000
    # try:
    for filePrefix in fileNamePathPre:
        fileName = filePrefix + dataFileSuffix
        X, Y = [], []
        for line in open(fileName, 'r'):
        #   print(line)
            values = [float(s) for s in line.split()]
        #   print(values)
            X.append(values[0])
            Y.append(values[1])
        # fig, ax = plt.subplots(num = figNum)
        X = np.array(X,dtype= float)
        if (Convert2Deg == True):
            X = X * 180 / np.pi
        Y = np.array(Y,dtype= float)
        X,Y = cleanData(X,Y)
        if(plotLegend == False):

            if LOGLOGPLOT == True:
                ax.loglog(X,Y,pointType)
            elif SEMILOGYPLOT == True:
                ax.semilogy(X,Y,pointType)
            elif SEMILOGXPLOT == True:
                ax.semilogx(X,Y,pointType)
            else:
                ax.plot(X, Y)
        else :
            if LOGLOGPLOT == True:
                ax.loglog(X,Y,pointType)
            elif SEMILOGYPLOT == True:
                ax.semilogy(X,Y,pointType)
            elif SEMILOGXPLOT == True:
                ax.semilogx(X,Y,pointType)
            else:
                ax.plot(X, Y,label = legend)
        if (Convert2Deg == True):
            plt.xlabel(r'$\theta\:(degree)$',fontsize = 10)
        else:
            plt.xlabel(r'$\theta$',fontsize = 10)

        if (setAxisRange == False):
            xMax = round(max(X),1)
            setAxisRange = True
        
        if (max(Y) > yMax):
            yMax = max(Y)
        if (min(Y) < yMin):
            yMin = min(Y)

    ax.xaxis.set_major_locator(plt.MultipleLocator(xMax/MAJOR_TICKS_NUM))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(xMax/(MAJOR_TICKS_NUM * 2)))
    # plt.xticks(np.arange(0,max(X),sxtep = 0.1))
    plt.ylabel(r'$\omega(\theta)$',fontsize = 10,rotation = rotateAngle)
    # plt.yticks(np.arange(min(Y)-0.1, max(Y)+0.1, 0.05))
    # plt.yticks(np.arange(round(min(Y) - 0.1,1), round(max(Y) + 0.1,1), 0.05))
    # title = "Pi healpix Peebles with data points: " + str(hpPointCount)
    figTitle = figTitlePre
    if (adjustTitlePos == True):
        plt.title(figTitle, fontsize = 7, y = POS_TITLE)
    else:
        plt.title(figTitle,fontsize = 8)      
    # plt.title(figTitle,fontsize = 8)
    # if(plotLegend == False):
    #     plt.tight_layout()
    saveName = saveNamePre + "average:" + str(len(fileNamePathPre)) + ".png"
    if (LARGE_Y_RANGE == True):
        YMax = max(yMax,Y_MAX)
        YMin = min(yMin,Y_MIN)
        bottom,top = plt.ylim()
        plt.ylim(YMin, YMax)
        insertPos = saveName.find("png")
        saveTemp = saveName[:insertPos-1] + "largeScale" + saveName[insertPos-1:]
        plt.savefig(saveTemp)
        plotted += 1
        plt.ylim(bottom,top)

    plt.savefig(saveName)
    figureI += 1
    plotted += 1

def sweepMultiEstimator(fileNamePathPre, dataFileSuffix,figureSaveNamePre, pointCount, 
    randomCount, deltaTheta, average, spannedAngle, 
    figNums, compactPlot = True, adjustTitlePos = False,**kwargs):
    """ sweep the multi time estimator
    Paramters:
    ---------

    fileNamePathPre: list of strings
        prefix for the string of the data file, e.g ["../data/result0","../data/result1"]

    dataFileSuffix: str
        suffix of the data file, e.g "Pi_RandomMode"

    """ 
    global plotted, figureI
    peebles1 = dataFileSuffix + "1.txt"
    peebles2 = dataFileSuffix + "2.txt"
    hamilton = dataFileSuffix + "3.txt"
    landy = dataFileSuffix + "4.txt"
    additionName = ""
    for key,value in kwargs.items():
        additionName += key + ":" + str(value)
    
    corrName = "spanned:%.2f"%(spannedAngle) + "points:" + str(pointCount) + "Random:" + \
            str(randomCount) + "Average:" + str(average) + " dT: %.4f" %(deltaTheta) + \
                additionName
    strName = ["Peebles","DP","hamilton","LS"]
    for var in strName:
        var += corrName

    saveNameList = [figureSaveNamePre + var for var in strName]
    fileSuffix = [peebles1,peebles2,hamilton,landy]

    for suffix, figTitle,saveNamePre,fig in zip(fileSuffix,strName,saveNameList,figNums):
        fig,ax = plt.subplots(num = fig)
        multiTimeEstimatorPlot(fileNamePathPre,suffix,figTitle,saveNamePre,ax,False,"",adjustTitlePos)

def sweepCompactPlot(xvalues,yvalues, figTitle, saveName,figNum, ax,legends,adjustTitlePos = False):
    global figureI, plotted, cut_Val
    print("sweep Compact called: saveName: ", saveName)
    ## assume len(yvalues) = 4 / 8
    # for i, legend in zip(range(len(yvalues)), legends):
    if (len(xvalues[0]) == 0):
        return
    for i, legend in zip(range(4), legends):
        X = xvalues[i]
        Y = yvalues[i]
        if LOGLOGPLOT == True:
            ax.loglog(X,Y,pointType)
        elif SEMILOGYPLOT == True:
            ax.semilogy(X,Y,pointType)
        elif SEMILOGXPLOT == True:
            ax.semilogx(X,Y,pointType)
        else:
            ax.plot(X, Y,label = legend)
    if (Convert2Deg == True):
        plt.xlabel(r'$\theta\:(degree)$',fontsize = 10)
    else:
        plt.xlabel(r'$\theta$',fontsize = 10)
    xMax = round(max(X),1)
    
    ax.xaxis.set_major_locator(plt.MultipleLocator(xMax/MAJOR_TICKS_NUM))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(xMax/(MAJOR_TICKS_NUM * 2)))
    # plt.xticks(np.arange(0,max(X),step = 0.1))
    plt.ylabel(r'$\omega(\theta)$',fontsize = 10,rotation = rotateAngle)
    # plt.yticks(np.arange(min(Y)-0.1, max(Y)+0.1, 0.05))
    # plt.yticks(np.arange(round(min(Y) - 0.1,1), round(max(Y) + 0.1,1), 0.05))
    # title = "Pi healpix Peebles with data points: " + str(hpPointCount)
    if (figTitle[-1] == "\n"):
        figTitle = figTitle[:-1]
    if (figTitle[-2] == ","):
        figTitle = figTitle[:-2]

    if (adjustTitlePos == True):
        plt.title(figTitle, fontsize = 7, y = POS_TITLE)
    else:
        plt.title(figTitle,fontsize = 8)      
    # plt.title(figTitle,fontsize = 8)



    YMax = max(max(Y),Y_MAX)
    YMin = min(min(Y),Y_MIN)
    chartBox = ax.get_position()
    ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.85, chartBox.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), shadow=True, ncol=1)
    if (LARGE_Y_RANGE == True):
        bottom,top = plt.ylim()
        plt.ylim(YMin, YMax)
        insertPos = saveName.find("png")
        saveTemp = saveName[:insertPos-1] + "largeScale" + saveName[insertPos-1:]
        plt.savefig(saveTemp)
        plotted += 1
        plt.ylim(bottom,top)
    
    # if(plotLegend == False):
    #     plt.tight_layout()

    plt.savefig(saveName)
    plotted += 1

    # cut range
    if (CUT_RANGE == True):
        figureI += 2
        figTemp, axTemp = plt.subplots(num = figureI)
        figureI += 1
        for i, legend in zip(range(4,8), legends):
            X = xvalues[i]
            Y = yvalues[i]
            if LOGLOGPLOT == True:
                axTemp.loglog(X,Y,pointType)
            elif SEMILOGYPLOT == True:
                axTemp.semilogy(X,Y,pointType)
            elif SEMILOGXPLOT == True:
                axTemp.semilogx(X,Y,pointType)
            else:
                axTemp.plot(X, Y,label = legend)
        if (Convert2Deg == True):
            plt.xlabel(r'$\theta\:(degree)$',fontsize = 10)
        else:
            plt.xlabel(r'$\theta$',fontsize = 10)
        xMax = round(max(X),1)
        
        axTemp.xaxis.set_major_locator(plt.MultipleLocator(xMax/MAJOR_TICKS_NUM))
        axTemp.xaxis.set_minor_locator(plt.MultipleLocator(xMax/(MAJOR_TICKS_NUM * 2)))
        # plt.xticks(np.arange(0,max(X),step = 0.1))
        plt.ylabel(r'$\omega(\theta)$',fontsize = 10,rotation = rotateAngle)
        # plt.yticks(np.arange(min(Y)-0.1, max(Y)+0.1, 0.05))
        # plt.yticks(np.arange(round(min(Y) - 0.1,1), round(max(Y) + 0.1,1), 0.05))
        # title = "Pi healpix Peebles with data points: " + str(hpPointCount)

        if (adjustTitlePos == True):
            plt.title(figTitle, fontsize = 7, y = POS_TITLE)
        else:
            plt.title(figTitle,fontsize = 8)              
        # plt.title(figTitle,fontsize = 8)

        YMax = max(max(Y),Y_MAX)
        YMin = min(min(Y),Y_MIN)
        chartBox = axTemp.get_position()
        axTemp.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.85, chartBox.height])
        axTemp.legend(loc='center left', bbox_to_anchor=(1, 0.5), shadow=True, ncol=1)

        insertPos = saveName.find("png")
        saveName = saveName[:insertPos-1] + "CutRange" + saveName[insertPos-1:]
        if (LARGE_Y_RANGE == True):
            bottom,top = plt.ylim()
            plt.ylim(YMin, YMax)
            insertPos = saveName.find("png")
            saveTemp = saveName[:insertPos-1] + "largeScale" + saveName[insertPos-1:]
            plt.savefig(saveTemp)
            plotted += 1
            plt.ylim(bottom,top)        
    plt.savefig(saveName)
    plotted += 1

def plotMultiLines(fileLists, legendLists,figureI, 
    saveFilePath, title, suffixes = ['.pdf','.png'], 
    saveName = '', thresholdPlotted = 2,subdirectoryPrefix = '', method = "correlation"):
    """
    plot multi lines in a figure

    ---
    Parameters:
    ---

    subdirectoryPrefix: str
        added for angular power spectrum, so it is equal to "angularPower_",
        so the final directory where the plot is saved is then saveFilePath/angularPower_pdf/...
    """
    fig,ax = plt.subplots(num = figureI)
    # count how many curves actually plotted, because there may be null array
    realPlotted = 0
    for file,legend in zip(fileLists,legendLists):
        x = []
        y = []
        for line in open(file,'r'):
            values = [float(s) for s in line.split()]
            if (len(values) >= 1):
                x.append(values[0])
                y.append(values[1])
        x = np.array(x, dtype = float)
        y = np.array(y, dtype = float)
        if (method == "correlation"):
            if (Convert2Deg == True):
                x = x * 180 / np.pi
        x, y = cleanData(x, y)
        if(len(x) >= 1):
            ax.plot(x,y,label = legend)
            realPlotted += 1
        else:
            legendLists.remove(legend)
    # chartBox = ax.get_position()
    # ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.8, chartBox.height])
    # ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), shadow=True, ncol=1)    
    adjustLegendPadding(plt,ax)
    # plt.legend()
    plt.title(title)
    if (method == "correlation"):
        if (Convert2Deg == True):
            plt.xlabel(r'$\theta\:(degree)$',fontsize = 10)
        else:
            plt.xlabel(r'$\theta$',fontsize = 10)    
        plt.ylabel(r'$\omega(\theta)$',fontsize = 10,rotation = rotateAngle)    
    elif (method == "angularPower"):
        angularPowerAxisLabel(plt,ax)
    elif (method == "angularPowerScaled"):
        angularPowerAxisLabel(plt,ax, scaled=True)
    else:
        print("Unknown method in plotMultiLines\n")
        return

    # only save figure if there are at least 2 curves
    if (realPlotted >= thresholdPlotted):
        for suffix in suffixes:
            ending = suffix.lstrip('.')
            saveFilePathFinal = os.path.join(saveFilePath,subdirectoryPrefix + ending)
            if not os.path.exists(saveFilePathFinal):
                os.makedirs(saveFilePathFinal)
            if (saveName == ''):
                saveName = title
            if (VERBOSE):
                savedPath = os.path.join(saveFilePathFinal, saveName + suffix)
                print("savedfig path: {}".format(savedPath))            

            plt.savefig(os.path.join(saveFilePathFinal, saveName + suffix))
    plt.close('all')

def plotRayTraceMeas(rayTraceFile, measFile, figureI, 
    savePath, title,
    legends = ["ray_tracing","measured"], 
    suffixes = ['.pdf','.png'], 
    method = "correlation"):
    """
    plot the estimators for ray tracing and measured

    ---
    Parameters::
    ---
    method: str
        either "correlation" or "angularPower"

    """
    # plotMultiLines([rayTraceFile, measFile],legends , figureI,savePath, title, thresholdPlotted = 1)
    DoNotPlot = False
    if (method == "correlation"):
        x_rt = []
        y_rt = []
        x_meas = []
        y_meas = []
        for line in open(rayTraceFile,'r'):
            values = [float(s) for s in line.split()]
            if (len(values) >= 1):
                x_rt.append(values[0])
                y_rt.append(values[1])
        
        # read in measured
        for line in open(measFile, 'r'):
            values = [float(s) for s in line.split()]
            if (len(values) >= 1):
                x_meas.append(values[0])
                y_meas.append(values[1])
        fig,ax = plt.subplots(num = figureI)

        x_rt = np.array(x_rt, dtype = float)
        x_meas = np.array(x_meas, dtype = float)
        if (Convert2Deg == True):
            x_rt = x_rt * 180 / np.pi
            x_meas = x_meas * 180 / np.pi

        x_rt, y_rt = cleanData(np.array(x_rt), np.array(y_rt))
        x_meas, y_meas = cleanData(np.array(x_meas),np.array(y_meas))
        ax.plot(x_rt, y_rt,'r',label = legends[0])
        ax.plot(x_meas,y_meas,'b',label = legends[1])
        chartBox = ax.get_position()
        ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.8, chartBox.height])
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), shadow=True, ncol=1)    

        if (Convert2Deg == True):
            plt.xlabel(r'$\theta\:(degree)$',fontsize = 10)
        else:
            plt.xlabel(r'$\theta$',fontsize = 10)    
        plt.ylabel(r'$\omega(\theta)$',fontsize = 10,rotation = rotateAngle)    

    else:
        DoNotPlot = plotMultiCls([rayTraceFile,measFile],figureI,legendLists=legends,colors = ["r","b"],figureTitleName=title,saveName=savePath,fignum=figureI,plt = plt, writeOutFile = True)

    if (DoNotPlot == False):
        plt.title(title)
        for suffix in suffixes:
            savePathFinal = os.path.join(savePath, method + suffix.lstrip('.'))
            if not os.path.exists(savePathFinal):
                os.makedirs(savePathFinal)
            saveFileName = os.path.join(savePathFinal, title + suffix)
            plt.savefig(saveFileName)

        


def getCalcuInfo(calculationInfoPath):
    """
    given the calculationinfo.txt, return a dict
    """

                #     outFile << pointsCount << endl;
                # outFile << randomCount << endl;
                # outFile << deltaTheta << endl;
                # outFile << AVERAGE_TIME << endl;
                # outFile << spannedAngle << endl;
                # outFile << endTheta << endl;
                # outFile << thetaBegin << endl;
                # outFile << thetaEnd << endl;
                # outFile << threshold << endl;
                # outFile << FACTOR_SIGVAL << endl;
                # outFile << outFilePre << endl;
    paramList = ["pointsCount","randomCount","deltaTheta","AVERAGE_TIME","spannedAngle",\
        "endTheta","thetaBegin","thetaEnd","threshold","FACTOR_SIGVAL","outFilePre"]
    # print("paramList length")
    # print(len(paramList))
    count = 0
    infoDict = {}
    try:
        for line in open(calculationInfoPath,'r'):
            values = [s for s in line.split()]
            if (len(values) >= 1):
                infoDict[paramList[count]] = values[0]
            count += 1
    except IOError:
        print("calculationInfo.txt does not exist")
        return None
    # print("count:")
    # print(count)
    return infoDict

def simplifyName(name):
    """
    simplify the name, like HM_RXA_horn_MAXRSS to 1ah
    """
    segments = name.split('_')
    hTr = segments[2]
    if (hTr == "horn"):
        antenna = 'h'
    else:
        antenna = 'a'

    rx = segments[1]
    if (len(rx) == 3):
        rxName = '1' + rx[-1]
    else :
        rxName = '2' + rx[-1]
    finalName = rxName + antenna
    return finalName

def getDataDictCompact(dataFolderPath, estimator = '4', method = "correlation", directReadFromFile = False, \
    saveDict = False, saveDictPath = "./dataDict.txt",\
    readDictPath = "./dataDict.txt",\
    generateAngularPower = False, angularPowerFileName = "angularPower.txt"):
    """ 
    return the compact data dictionary
    {
        'HM_RXA_horn_MAXRSS' : {
            '75': {
                'meas' : path
                'meas_info': dictionary
                'rayTrace' : pathname
                'rayTrace_info': dictionary
            }    
        }
    }

    ----
    Parameter:
    ---

    dataFolderPath: str
        indicate the path to the data folder

    method: str
        either "correlation" or "angularPower"

    directReadFromFile: bool
        used to specify whether or not the dictionary will be generated from parsing the calulation data folder,
        or directly load from a file where the dictionary is already saved, note this value
        will also change the meaning of dataFolderPath
        default: False
        if set to true and readDictPath and saveDictPath is the same, then automatically saveDict will be set to false

    saveDict: bool
        used to specify whether or not the dictionary will be saved,
        format json
        default false,

    saveDictPath: str
        used to specify the path to save the dictionary,
        only useful if saveDict is set to be True 
        default "./dataDict.txt"
    
    readDictPath: str
        used to read dictionary data from.
        Only useful if directReadFromFile is set to be true
        default: "./dataDict.txt"

    generateAngularPower: boolean
        if set, will generate angularPower spectrum
    
    angularPowerFileName: str
        specify the fillename which stores the angular power spectrum
        default "angularPower.txt"
        Note that this is the name for the pure angular power spectrum, without scaling, even 
        when the method is set to angularPowerScaled, the name of the pure angular power spectrum file is till angularPowerFileName
    ----
    Convention
    ----

    The data folder should have naming structure like
    TX1_RXA_raytracing or TX1_RXA_meas
    """
    dataDict = {}

    if (directReadFromFile == False):
        for root, direct, files in os.walk(dataFolderPath):
            # if(len(direct) >= 1 and direct[0][0] == 'm'):
            #     # in the mode directory
            if (len(files) >= 1):
                # in the file folder
                calcuInfoDict = getCalcuInfo(os.path.join(root,"calculationInfo.txt"))
                if (calcuInfoDict is None ):
                    continue
                # calcuInfo = pd.read_csv(os.path.join(root,"calculationInfo.txt"), header = None)
                # simName = calcuInfo.iloc[-1][0]

                # get the threshold folder 
                thresholdVal = calcuInfoDict['threshold']
                # if not thresholdVal in dataDict:
                #     dataDict[thresholdVal] = {}
                simName = calcuInfoDict['outFilePre']
                simNamePrune = "_".join(simName.split("_")[:-1])
                if (simName.split("_")[-1] == "meas"):
                    raise Exception("Need to copy correctly from raytracing")
                    cutIndex = simName.rfind("_")
                    simNamePrune = simName[0:cutIndex]
                    # print("measured, simNamePrune")
                    # print(simNamePrune)
                    if not simNamePrune in dataDict:
                        dataDict[simNamePrune] = {}
                    if not thresholdVal in dataDict[simNamePrune]:
                        dataDict[simNamePrune][thresholdVal] = {}
                #     if (method == "correlation"):
                #         dataDict[simNamePrune][thresholdVal]['meas'] = os.path.join(root,"SignalAverageHealpixMode"+estimator+".txt")
                #     else:
                #         dataDict[simNamePrune][thresholdVal]['meas'] = os.path.join(root,"SignalAverageHealpixMap.txt")

                #     dataDict[simNamePrune][thresholdVal]['meas_info'] = calcuInfoDict                
                    if (method == "correlation"):
                        filePath = os.path.join(root,"SignalAverageHealpixMode"+estimator+".txt")
                        dataDict[simNamePrune][thresholdVal]['meas'] = filePath
                        if (generateAngularPower == True):
                            writeAngularPower([filePath])
                    elif (method == "angularPower"):
                        raise Exception("Need to consider angularPowerScaled")
                        filePath = os.path.join(root,"SignalAverageHealpixMap.txt")
                        angularPowerPath = os.path.join(root, angularPowerFileName)
                        # generate angularPower txt
                        writtenNum = -1
                        if (generateAngularPower == True or not os.path.exists(angularPowerPath)):
                            writtenNum = writeAngularPower([filePath], outName=angularPowerFileName)
                        if (writtenNum == 1):
                            dataDict[simNamePrune][thresholdVal]['meas'] = os.path.join(angularPowerPath)

                    dataDict[simNamePrune][thresholdVal]['meas_info'] = calcuInfoDict
                    
                elif (simName.split("_")[-1] == "raytracing"):
                    if not simNamePrune in dataDict:
                        dataDict[simNamePrune] = {}
                    if not thresholdVal in dataDict[simNamePrune]:
                        dataDict[simNamePrune][thresholdVal] = {}
                    if (method == "correlation"):
                        filePath = os.path.join(root,"SignalAverageHealpixMode"+estimator+".txt")
                        dataDict[simNamePrune][thresholdVal]['rayTrace'] = filePath
                        if (generateAngularPower == True):
                            writeAngularPower([filePath])
                    elif (method == "angularPower" or method == "angularPowerScaled"):
                        filePath = os.path.join(root,"SignalAverageHealpixMap.txt")
                        angularPowerPath = os.path.join(root, angularPowerFileName)
                        # print("treat angular power file as the original one, not scaled")
                        # generate angularPower txt
                        writtenNum = -1
                        angularPowerScaledName = angularPowerFileName.split(".")[0] + "_scaled.txt"
                        angularPowerScaledPath = os.path.join(root,angularPowerScaledName)
                        # assume coexistence of angularPowerFile and angularPowerScaledPath
                        if (generateAngularPower == True or not os.path.exists(angularPowerPath) or not os.path.exists(angularPowerScaledPath)):
                            writtenNum = writeAngularPower([filePath], outName=angularPowerFileName)
                        if (writtenNum == 1 or os.path.exists(angularPowerPath)):
                            if (method == "angularPower"):
                                dataDict[simNamePrune][thresholdVal]['rayTrace'] = os.path.join(angularPowerPath)
                            else:
                                # method == "angularPowerScaled"
                                dataDict[simNamePrune][thresholdVal]['rayTrace'] = os.path.join(angularPowerScaledPath)


                    dataDict[simNamePrune][thresholdVal]['rayTrace_info'] = calcuInfoDict
    else:
        dataDict = json.load(open(readDictPath,"r"))
        if (readDictPath == saveDictPath):
            saveDict = False
        # need to modify the corresponding entries if necessary
        NEED_MODIFY = True
        corrlelationFileName = "SignalAverageHealpixMode"+estimator+".txt"
        angularPowerScaledName = angularPowerFileName.split(".")[0] + "_scaled.txt"
        for config, thresholdDict in dataDict.items():
            for threshold, calcuPath in thresholdDict.items():
                models = []
                if 'rayTrace' in calcuPath.keys():
                    models.append('rayTrace')

                if 'meas' in calcuPath.keys():
                    models.append('meas')
                
                if len(models) < 1:
                    raise Exception("Unknown calcuPath keys!")
                for model in models:
                    fileName = calcuPath[model].split("/")[-1]
                    rootPath = os.path.dirname(fileName)
                    if (method == "angularPower" and fileName == angularPowerFileName 
                    or method == "angularPowerScaled" and fileName == angularPowerScaledName or method == "correlation" and fileName == corrlelationFileName):
                        NEED_MODIFY = False
                        break
                    if (method == "correlation"):
                        calcuPath[model] = os.path.join(rootPath,corrlelationFileName)
                    elif (method == "angularPower"):
                        calcuPath[model] = os.path.join(rootPath,angularPowerFileName)
                    elif (method == "angularPower_scaled"):
                        calcuPath[model] = os.path.join(rootPath,angularPowerScaledName)
                    else:
                        raise Exception("Unknown method!")
                if (NEED_MODIFY == False):
                    break
            if (NEED_MODIFY == False):
                break

    if (saveDict == True):
        with open(saveDictPath, 'w') as file:
            file.write(json.dumps(dataDict))
    return dataDict    

def plotFixConfigDiffThreshold(dataDict,
    savePathPre, estimator = '4', modelList = ['meas','rayTrace'], 
    method = "correlation"):
    """n
    plot the figure for fixed configuration with different threshold
    """

    for config, thresholdDict in dataDict.items():
        configName = simplifyName(config)
        for model in modelList:
            thresholdPathDict = {}
            legendLists = []
            fileLists = []
            for threshold, calcuPath in thresholdDict.items():
                thresholdPathDict[threshold] = calcuPath[model]
                fileLists.append(calcuPath[model])
                legendLists.append(threshold + "dBm")
            if (method == "correlation"):
                plotMultiLines(fileLists, legendLists, figureI, \
                    os.path.join(savePathPre,"sameConfigDiffDBM"), \
                        configName +"_"+ model)
            else:
                # change all files in the fileLists to "angularPower.txt"
                ToPlot = True
                angularPowerFiles = []
                angularPowerLegends = []
                for file,legend in zip(fileLists,legendLists):
                    slicedDir = file.split("/")
                    subdirect = "/".join(slicedDir[:-1])
                    fileNew = os.path.join(subdirect,"angularPower.txt")
                    if (not os.path.exists(fileNew)):
                        if (VERBOSE):
                            print("-----------------")
                            print("-----------------")
                            print("{} does not exist\n the file will be dropped".format(file))
                            print("the config: {}, model: {}".format(config,model))
                            print("removed file: {}\n removed legend{}".format(file,legend))
                        else:
                            1+1
                        # ToPlot = False
                        # break
                        # fileLists.remove(file)
                        # legendLists.remove(legend)
                    else:
                        if (VERBOSE):
                            print("++++++++++++++++++++++")
                            print("++++++++++++++++++++++")
                            print("{} exists\n the file will be added".format(file))
                            print("the config: {}, model: {}".format(config,model))
                            print("added file: {}\n added legend{}".format(file,legend))
                        
                        # change file to fileNew
                        angularPowerFiles.append(fileNew)
                        angularPowerLegends.append(legend)
                if (ToPlot):
                    if (VERBOSE):
                        print("config: {} will be plotted\n".format(config))
                    plotMultiLines(angularPowerFiles,angularPowerLegends,figureI,\
                    os.path.join(savePathPre, "sameConfigDiffDBM"),\
                        configName + "_" + model,subdirectoryPrefix="angularPower",method = "angularPower")
            
                # plotMultiCls(fileLists,figureI,legendLists=legendLists,\
                #     colors = ['r','b','g'],configName + "_" + model,\
                #         os.path.join(savePathPre,"sameConfigDiffDBM"),fignum = figureI,plt)

def plotFixModelDBMDiffLocation(dataDict, 
    savePathPre, estimator = '4', models = ['meas','rayTrace'], 
    thresholds = ['-100','-95','-90','-85','-80','-75','-70'],
    method = "correlation"):
    """
    plot the figure with fixed model (rayTracing or measurement), fixed dBm value,
    but different location,
    possible combination:
    I) 1Ah, 1Bh, 1Ch
    II) 1Aa,1Ba, 1Ca
    III) 2Ba, 2Da
    IV) 1Aa,1Ba,1Ca,2Ba,2Da
    V) 1Ba,2Ba
    """

    # configuration 
    configurationLists = dataDict.keys()
    config2NameDict = {}
    Name2configDict = {}
    for key in configurationLists:
        config2NameDict[key] = simplifyName(key)
        Name2configDict[simplifyName(key)] = key
    # print("config2NameDict")
    # print(config2NameDict)
    # print("Name2configDict")
    # print(Name2configDict)

    combinationLists = [
        ['1Ah','1Bh','1Ch'],
        ['1Aa','1Ba','1Ca'],
        ['2Ba','2Da'],
        ['1Aa','1Ba','1Ca','2Ba','2Da'],
        ['1Ba','2Ba']
    ]

    for combination in combinationLists:
        for model in models:
            for threshold in thresholds:
                filePaths = []
                legends = []                
                for name in combination:
                    config = Name2configDict[name]
                    filePaths.append(dataDict[config][threshold][model])
                    legends.append(name)
                # all names collected, ready to plot
                title = model + '_' + threshold + "dBm"
                if (method == "correlation"):
                    plotMultiLines(filePaths, legends,figureI,\
                        os.path.join(savePathPre, "sameModeldBmDiffLoc"),\
                            title, suffixes= ['.pdf','.png'], saveName= title + ''.join(combination))
                elif (method == "angularPower"):
                    ToPlot = True
                    angularPowerFiles = []
                    angularPowerLegends = []
                    for file,legend in zip(filePaths,legends):
                        slicedDir = file.split("/")
                        subdirect = "/".join(slicedDir[:-1])
                        fileNew = os.path.join(subdirect,"angularPower.txt")
                        if (not os.path.exists(fileNew)):
                            if (VERBOSE):
                                print("-----------------")
                                print("-----------------")
                                print("{} does not exist\n the file will be dropped".format(file))
                                print("the combination: {}, model: {}".format(combination,model))
                                print("removed file: {}\n removed legend{}".format(file,legend))
                            else:
                                1+1
                            # ToPlot = False
                            # break
                            # fileLists.remove(file)
                            # legendLists.remove(legend)
                        else:
                            if (VERBOSE):
                                print("++++++++++++++++++++++")
                                print("++++++++++++++++++++++")
                                print("{} exists\n the file will be added".format(file))
                                print("the config: {}, model: {}".format(config,model))
                                print("added file: {}\n added legend{}".format(file,legend))
                            
                            # change file to fileNew
                            angularPowerFiles.append(fileNew)
                            angularPowerLegends.append(legend)
                    if (ToPlot):
                        if (VERBOSE):
                            print("config: {} will be plotted\n".format(config))
                        plotMultiLines(angularPowerFiles,angularPowerLegends,figureI,\
                            os.path.join(savePathPre, "sameModeldBmDiffLoc"),\
                            title,suffixes= ['.pdf','.png'],saveName= title + ''.join(combination), subdirectoryPrefix="angularPower",method = "angularPower")

def plotFixedModelTypeDiffPos(dataDict, savePathPre, estimator = '4',\
    model = 'rayTrace',\
    thresholdValues = ['-75','-85','-90','-100'],\
    matchPatterns = ["TX.*RX[0-9]{2}"],\
    method = "correlation",
    antennaPattern = "hornAntenna",
    singleThresholdPlot = True,
    plot_HealpixMap = False,
    plot_SphericalPoints = False):
    """
    plot the estimator for one type, e.g. either raytracing or experiment,
    use different threshold values, different configuration
    e.g.
    TX1_RX1, TX2_RX1, thresholdValues = ['100'] will compare the estimator
    between the two configurations, with fixed rx and varying tx.

    ----
    Parameters:
    ----
    model: str
        used to specify the type of the model, e.g., experiment or raytracing
        possible value: 'meas', 'rayTrace'

    matchPatterns: list of strs
        each str used to match certain configuration

    method: str
        "correlation"(default) or "angularPower" or "angularPowerScaled"

    antennaPattern: str
        indicate the pattern of the antenna, default "hornantenna"

    singleThresholdPlot: boolean
        if set true (default), will also plot with single threshold value
    
    plot_HealpixMap: boolean
        if set true, then will also plot the healpix map for each single configuration
        default False
    
    plot_SphericalPoint: boolean
        if set true, will also plot the spherical point for each single configuration,
        default False
    """

    global figureI

    figureI = 1
    for matchPattern in matchPatterns:
        fileLists = []
        legendLists = []
        fileListsSingleThd = {}
        legendListsSingleThd = {}        
        for config, thresholdDict in dataDict.items():
            x = re.search(matchPattern, config)

            if x:
                if VERBOSE:
                    print("{0:s} matches with the pattern {1:s}\n".format(config, matchPattern))
                for threshold, calcuPath in thresholdDict.items():
                    # check if threshold in thresholdValues, and if calcuPath has the key model and if it does, is the fileName not empty
                    if (threshold in thresholdValues and model in calcuPath.keys() and len(calcuPath[model]) >= 1):
                        fileLists.append(calcuPath[model])
                        legendLists.append(config + threshold + "dBm")
                        if (singleThresholdPlot == True):
                            for singleThd in thresholdValues:
                                if threshold == singleThd:
                                    if singleThd not in fileListsSingleThd.keys():
                                        fileListsSingleThd[singleThd] = []
                                        legendListsSingleThd[singleThd] = []
                                    fileListsSingleThd[singleThd].append(calcuPath[model])
                                    legendListsSingleThd[singleThd].append(config + threshold + "dBm")
                        if(plot_HealpixMap == True):
                            print("print the healpixMap for {0:s}, threshold: {1:s}".format(config, threshold))
                            printHealpixMapWrapper(dataDict, savePathPre, config, threshold, model)
                                        
        if (len(fileLists) >= 1):
            if (figureI >= 18):
                figureI = 1
            else:
                figureI += 1
            print("matchPattern: {}".format(matchPattern))
            print("legendLists:")
            print(legendLists)
            plotMultiLines(fileLists, legendLists, figureI,\
                os.path.join(savePathPre,method),\
                    model + " " + antennaPattern, saveName = matchPattern + str(thresholdValues), method = method)
            # since fileLists has cardinality >= each singleThd lists, only in this case is it possible that the singleThd is not empty
            if (singleThresholdPlot == True):
                for singleThd in thresholdValues:
                    # need to test first if singleThd is a key in fileListsSingleThd
                    if (singleThd in  fileListsSingleThd.keys() and len(fileListsSingleThd[singleThd]) >= 1):
                        plotMultiLines(fileListsSingleThd[singleThd], legendListsSingleThd[singleThd], figureI,\
                            os.path.join(savePathPre,method),\
                            model+ " " + antennaPattern, saveName = matchPattern + str(singleThd), method = method)
                        figureI += 1
                        
        else:
            print("Pattern: {0:s} does not match any file".format(matchPattern))
        

        # clean the file and legend lists
        fileLists = []
        legendLists = []
        fileListsSingleThd = {}
        legendListsSingleThd = {}        

def plotWholeRayMeas(dataFolderPath,estimator = '4', method = "correlation"):
    """
    plot the estimators for ray tracing and measurements

    --- 
    Parameters:
    ---
    method: str
        either "correlation" or "angularPower"

    """
    # dataDict has key the threshold (as string)
    # each value is a dictionary, with two keys, the one with rayTracing, 
    # and the other with meas
    # like {-75: { HM_RXA_horn_MAXRSS : a dictionary,
#                 HM_RXB_TRX_MAXRSS : a dictionry,...},
        #   -80: { ...},... }
    # dataDict = {}

    # for root, direct, files in os.walk(dataFolderPath):
    #     # if(len(direct) >= 1 and direct[0][0] == 'm'):
    #     #     # in the mode directory
    #     if (len(files) >= 1):
    #         # in the file folder
    #         calcuInfoDict = getCalcuInfo(os.path.join(root,"calculationInfo.txt"))
    #         # calcuInfo = pd.read_csv(os.path.join(root,"calculationInfo.txt"), header = None)
    #         # simName = calcuInfo.iloc[-1][0]

    #         # get the threshold folder 
    #         thresholdVal = calcuInfoDict['threshold']
    #         if not thresholdVal in dataDict:
    #             dataDict[thresholdVal] = {}
    #         simName = calcuInfoDict['outFilePre']
    #         if (simName.split("_")[-1] == "meas"):
    #             cutIndex = simName.rfind("_")
    #             simNamePrune = simName[0:cutIndex]
    #             print("measured, simNamePrune")
    #             print(simNamePrune)
    #             if not simNamePrune in dataDict:
    #                 dataDict[thresholdVal][simNamePrune] = {}

    #             dataDict[thresholdVal][simNamePrune]['meas'] = os.path.join(root,"SignalAverageHealpixMode"+estimator+".txt")
    #             dataDict[thresholdVal][simNamePrune]['meas_info'] = calcuInfoDict                
    #         else:
    #             if not simName in dataDict[thresholdVal]:
    #                 dataDict[thresholdVal][simName] = {}
    #             dataDict[thresholdVal][simName]['rayTrace'] = os.path.join(root,"SignalAverageHealpixMode"+estimator+".txt")
    #             dataDict[thresholdVal][simName]['rayTrace_info'] = calcuInfoDict
    
    # figureI = 1
    # for key, value in dataDict.items():
    #     dbValue = key
    #     for configuration, infoDict in value.items():
    #         titleName = dbValue + "_" + simplifyName(configuration)
    #         plotRayTraceMeas(infoDict['rayTrace'], infoDict['meas'], figureI, "/home/panwei/Desktop/bachelorThese/photoData/Langenfeld", titleName)
    #     figureI += 1

    # dataDict = {}

    # for root, direct, files in os.walk(dataFolderPath):
    #     # if(len(direct) >= 1 and direct[0][0] == 'm'):
    #     #     # in the mode directory
    #     if (len(files) >= 1):
    #         # in the file folder
    #         calcuInfoDict = getCalcuInfo(os.path.join(root,"calculationInfo.txt"))
    #         # calcuInfo = pd.read_csv(os.path.join(root,"calculationInfo.txt"), header = None)
    #         # simName = calcuInfo.iloc[-1][0]

    #         # get the threshold folder 
    #         thresholdVal = calcuInfoDict['threshold']
    #         # if not thresholdVal in dataDict:
    #         #     dataDict[thresholdVal] = {}
    #         simName = calcuInfoDict['outFilePre']
    #         if (simName.split("_")[-1] == "meas"):
    #             cutIndex = simName.rfind("_")
    #             simNamePrune = simName[0:cutIndex]
    #             print("measured, simNamePrune")
    #             print(simNamePrune)
    #             if not simNamePrune in dataDict:
    #                 dataDict[simNamePrune] = {}
    #             if not thresholdVal in dataDict[simNamePrune]:
    #                 dataDict[simNamePrune][thresholdVal] = {}

    #             dataDict[simNamePrune][thresholdVal]['meas'] = os.path.join(root,"SignalAverageHealpixMode"+estimator+".txt")
    #             dataDict[simNamePrune][thresholdVal]['meas_info'] = calcuInfoDict                
    #         else:
    #             if not simName in dataDict:
    #                 dataDict[simName] = {}
    #             if not thresholdVal in dataDict[simName]:
    #                 dataDict[simName][thresholdVal] = {}
    #             dataDict[simName][thresholdVal]['rayTrace'] = os.path.join(root,"SignalAverageHealpixMode"+estimator+".txt")
    #             dataDict[simName][thresholdVal]['rayTrace_info'] = calcuInfoDict
    
    dataDict = getDataDictCompact(dataFolderPath = dataFolderPath, estimator = estimator, method = method)

    figureI = 1
    for configuration, value in dataDict.items():
        for dbValue, infoDict in value.items():
            titleName = dbValue + "dBm_" + simplifyName(configuration)
            plotRayTraceMeas(infoDict['rayTrace'], infoDict['meas'], figureI, "/home/panwei/Desktop/bachelorThese/photoData/Langenfeld/rayTraceExperiplot", titleName, method = method)
            plt.close('all')
    

    
    
# healpixMapFileName the file which stores the data
# figureSaveNamePre: prefix for the title
# saveNamePre: prefix for the name of the saved figure

# plot the envelope
# def plotEnvelope(fileNamePathPre, dataFileSuffix,ax,lineColor, lineType,plotLegend = True,legend = ""):
#     """ plot the envelope of several curves
#     assume all the xvalues are equal
#     for filePrefix in file

#     """
#     X,Ymin,Ymax = [],[]
#     for filePrefix in fileNamePathPre:
#         fileName = filePrefix + dataFileSuffix
#         xT,yT = [],[]
#         count = 0
#         for line in open(fileName, 'r'):
#         #   print(line)
#             values = [float(s) for s in line.split()]
#         #   print(values)
#             xT.append(values[0])
#             yT.append(values[1])
#             if(Ymin[count] > values[1]):
#                 Ymin[count] = values[1]
#             if (Ymax[count] < values[1]):
#                 Ymax[count] = values[1]

#         # fig, ax = plt.subplots(num = figNum)
#         if (len(X) != len(xT)):
#             print("len not equal")
#             return None
    
#         if(plotLegend == False):

#             if LOGLOGPLOT == True:
#                 ax.loglog(X,Ymax,color = lineColor, marker = lineType)
#                 ax.loglog(X,Ymin, color = lineColor, marker = lineType)
#             elif SEMILOGYPLOT == True:
#                 ax.semilogy(X,Ymax,color = lineColor, marker = lineType)
#                 ax.semilogy(X,Ymin, color = lineColor,marker = lineType)
#             elif SEMILOGXPLOT == True:
#                 ax.semilogx(X,Ymax,color = lineColor,marker = lineType)
#                 ax.semilogx(X,Ymin,color = lineColor, marker = lineType)
#             else:
#                 ax.plot(X, Ymax, color = lineColor, marker = lineType)
#                 ax.plot(X,Ymin, color = lineColor, marker = lineType)
#         else :
#             if LOGLOGPLOT == True:
#                 ax.loglog(X,Y,pointType)
#             elif SEMILOGYPLOT == True:
#                 ax.semilogy(X,Y,pointType)
#             elif SEMILOGXPLOT == True:
#                 ax.semilogx(X,Y,pointType)
#             else:
#                 ax.plot(X, Y,label = legend)
#         plt.xlabel(r'$\theta$',fontsize = 10)
#         if (setAxisRange == False):
#             xMax = round(max(X),1)
#             setAxisRange = True
#         if (max(Y) > yMax):
#             yMax = max(Y)
#         if (min(Y) < yMin):
#             yMin = min(Y)

#     ax.xaxis.set_major_locator(plt.MultipleLocator(xMax/MAJOR_TICKS_NUM))
#     ax.xaxis.set_minor_locator(plt.MultipleLocator(xMax/(MAJOR_TICKS_NUM * 2)))
#     # plt.xticks(np.arange(0,max(X),sxtep = 0.1))
#     plt.ylabel(r'$\omega(\theta)$',fontsize = 10,rotation = rotateAngle)
#     # plt.yticks(np.arange(min(Y)-0.1, max(Y)+0.1, 0.05))
#     # plt.yticks(np.arange(round(min(Y) - 0.1,1), round(max(Y) + 0.1,1), 0.05))
#     # title = "Pi healpix Peebles with data points: " + str(hpPointCount)
#     figTitle = figTitlePre

#     if (adjustTitlePos == True):
#         plt.title(figTitle, fontsize = 7, y = POS_TITLE)
#     else:
#         plt.title(figTitle,fontsize = 8)      

#     plt.title(figTitle,fontsize = 8)
#     # if(plotLegend == False):
#     #     plt.tight_layout()
#     saveName = saveNamePre + "average:" + str(len(fileNamePathPre)) + ".png"
#     if (LARGE_Y_RANGE == True):
#         YMax = max(yMax,Y_MAX)
#         YMin = min(yMin,Y_MIN)
#         bottom,top = plt.ylim()
#         plt.ylim(YMin, YMax)
#         insertPos = saveName.find("png")
#         saveTemp = saveName[:insertPos-1] + "largeScale" + saveName[insertPos-1:]
#         plt.savefig(saveTemp)
#         plt.ylim(bottom,top)

#     plt.savefig(saveName)
#     figureI += 1
#     plotted += 1        

def printHealpixMap(healpixMapFileName,figureTitleNamePre,\
    saveNamePre,pointCount, spannedAngle,saveFotoPath = ".", **kwargs):
    '''
    healpixMapFileName: data file name 
    figureTitleNamePre: the prefix for the title of the graph
    saveNamePre: the prefix for the to be saved figure

    pointCount : int

    spannedAngle : double
    saveFotoPath: str
        the path where the foto is saved, default "."
    '''
    global figureI, plotted, cut_Val

    corrName = "_spanned:%.4f"%(spannedAngle) +\
         " points:" + str(pointCount)
    additionName = ""
    for key,value in kwargs.items():
        additionName += key + ":" + str(value)
    
    figureName = figureTitleNamePre + ' ' + corrName + "\n" + additionName
    saveNameSuffix =   saveNamePre + corrName + additionName +  ".png"
    # in case the folder does not exist create one
    if not os.path.exists(saveFotoPath):
        os.makedirs(saveFotoPath)
    saveName = os.path.join(saveFotoPath, saveNameSuffix)
    hp1 = []
    # try:
    for line in open(healpixMapFileName,'r'):   
        if(SIGNAL_DATA == False or MY_DOUBLE == False):
            values = [float(s) for s in line.split()]
        else:
            values = [float(s) for s in line.split()]
        if (len(values) >= 1):
            hp1.append(values[0])

    hpOrder = hp1[0]
    hpPointCount = hp1[1]
    pixVal = hp1[2:]
    pixValLen = len(pixVal)
    Nside = 2**hpOrder
    Npixel = hp.nside2npix(Nside)
    pixVal.extend([0] * (int(Npixel) - pixValLen))
    pixVal = np.array(pixVal)
    hp.mollview(pixVal, title = '')
    hp.graticule()
    # plt.title()
    plt.savefig(saveName)
    figureI += 1
    plotted += 1

    plt.close('all')
     
    # except IOError:
    #     print("fileName starts with: " + healpixMapFileName + " not Found")

def printHealpixMapWrapper(dataDict, saveFotoPath, \
    config, threshold,\
        model = 'rayTrace'):
    """
    wrapper for the printHealpixMap

    """

    # first get the dictionary item
    if (config not in dataDict.keys()):
        raise Exception("The key {} not in dataDict fail to call plotSphericalPoint".format(config))
    
    responses = getItems(dataDict, config, threshold,model,
    queries=['healpixMap'])

    if (len(responses) < 1):
        raise Exception("No entry for the responses dictionary, cannot proceed!")
    
    titleNamePre = config + "_" + threshold + "_dBm"
    saveFotoPath = os.path.join(saveFotoPath, "HealpixMap")
    if not os.path.exists(saveFotoPath):
        os.makedirs(saveFotoPath)
    
    calcuInfo = responses['calcuInfo']
    pointsCount = calcuInfo['pointsCount']
    spannedAngle = calcuInfo['spannedAngle']
    if (pointsCount == '0'):
        print("zero points, no need to plot! config: {0:}  threshold: {1:}".format(config, threshold))
        return
    printHealpixMap(responses['healpixMap'], titleNamePre, \
        titleNamePre, pointCount=int(pointsCount),spannedAngle=float(spannedAngle),\
            saveFotoPath=saveFotoPath)


def plotSphericalPoint(sphericalPointDataFileName,figureTitleNamePre, \
    saveNamePre,pointCount, spannedAngle,**kwargs):
    """ Plot the 3d spherical point

    Parameters
    --------
    sphericalPointDataFileName: str
        name of the file that stores the data points

    figureTitleNamePre: str
        prefix of the title of the plot

    saveNamePre: str
        prefix of the saved figure
    """
    global figureI, plotted, cut_Val
    corrName = "spanned:%.4f"%(spannedAngle) + " points:" + str(pointCount)
    additionName = ""
    titleAddition = "spanned:%.4f"%(spannedAngle) + ", points:" + str(pointCount) + ", "
    count = 0
    for key,value in kwargs.items():
        if (count == 0):
            titleAddition += "\n"
        additionName += key + ":" + str(value)
        titleAddition += key + ":" + str(value) + ", "
        count += 1
        if (count % 3 == 0) :
            titleAddition+= "\n"
    X,Y,Z = [],[],[]
    titleName = figureTitleNamePre + "\n" + titleAddition
    if (titleName[-1] == "\n"):
        titleName = titleName[:-1]
    if (titleName[-2] == ","):
        titleName = titleName[:-2]
    with open(sphericalPointDataFileName,'r') as file_obj:
        for line in file_obj:
            values = [float(s) for s in line.split()]
            X.append(values[0])
            Y.append(values[1])
            Z.append(values[2])
    plt.figure(figureI + 1)
    figureI += 1
    ax = plt.axes(projection = '3d')
    ax.scatter3D(X,Y,Z,c = 'r',s = 2)
    ax.set_xlabel("x",labelpad = 15)
    ax.set_ylabel("y",labelpad = 15)
    ax.set_zlabel("z",labelpad = 10)
    if(len(X) <= 1):
        return
    saveName = saveNamePre + "sphericalPoint" + corrName + additionName

    zlimMin = min(Z) - 0.1
    zlimMax = max(Z)+ 0.1
    ax.set_xlim3d(-1.1,1.1)
    ax.set_ylim3d(-1.1,1.1)
    ax.set_zlim3d(zlimMin, zlimMax)
    # xMax = round(max(X),1)
    # yMax = round(max(Y),1)
    xMax = 1.0
    yMax = 1.0
    ax.xaxis.set_major_locator(plt.MultipleLocator(xMax/SPHERE_MAJOR_TICKS_NUM ))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(xMax/(SPHERE_MAJOR_TICKS_NUM  * 2)))
    ax.yaxis.set_major_locator(plt.MultipleLocator(yMax/SPHERE_MAJOR_TICKS_NUM ))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(yMax/(SPHERE_MAJOR_TICKS_NUM  * 2)))

    ax.set_aspect('equal', adjustable='box')
    ax.tick_params(axis='both', which='major', pad= 2)
    ax.tick_params(axis = 'z', which = 'major', pad = 5)
    # plt.title(titleName, fontsize = 8)
    if (NO_TITLE == False):
        ax.set_title(titleName, fontsize = 10, pad = 15)
    else :
        saveNameTemp = saveName + "noTitle"
    changeTick = False
    changeTitle = False
    for i in xrange(0,120,30):
        for ii in xrange(0,120,30):
            saveNameTemp = saveName + "elev:" + str(i) + "azim" + str(ii) + "*.png"
            ax.view_init(elev=i, azim=ii)
            if (i < 30) :
                ax.xaxis.set_major_locator(plt.MultipleLocator(xMax/1))
                ax.xaxis.set_minor_locator(plt.MultipleLocator(xMax/(1 * 2)))
                ax.yaxis.set_major_locator(plt.MultipleLocator(yMax/1))
                ax.yaxis.set_minor_locator(plt.MultipleLocator(yMax/(1* 2)))
                # ax.set_title(titleName, fontsize = 10, pad = 10)
                if (ii < 60):
                    # ax.set_title(titleName, fontsize = 10, pad = 5)
                    changeTitle = True
                else :
                    # ax.set_title(titleName, fontsize = 10, pad = 30)
                    changeTitle = False
                changeTick = True
                if (ii < 30):
                    ax.get_xaxis().set_ticks([])
                    ax.set_xlabel('x',labelpad = 0)
                elif (ii > 80):
                    ax.get_yaxis().set_ticks([])
                    ax.set_ylabel('y',labelpad = 0)
            # else :
                # ax.set_title(titleName, fontsize = 10, pad = 20)
            elif (i > 80):
                # ax.get_zaxis().set_ticks([])
                # ax.autoscale()
                ax.w_zaxis.line.set_lw(0.)
                ax.set_zticks([])
                ax.set_zlabel('z',labelpad = 0)
                if (ii > 80):
                    ax.set_xlabel('x',labelpad = 25)
                    # ax.set_title(titleName, fontsize = 10, pad = 10)
                    changeTitle = True
                if (ii < 30):
                    ax.set_ylabel('y',labelpad = 25)
                    # ax.set_title(titleName, fontsize = 10, pad = 10)
                    changeTitle = True
            else:
                if (ii < 30):
                    ax.tick_params(axis='y', which='major', pad= 1)
                    ax.set_ylabel('y',labelpad = 5)
                    changeTick = True
                if (ii > 80):
                    ax.tick_params(axis='x', which='major', pad= 1)
                    ax.set_xlabel('x',labelpad = 5)
                    changeTick = True
            plt.savefig(saveNameTemp,bbox_inches = 'tight')
            plotted += 1
            
            if (changeTick == True):
                ax.xaxis.set_major_locator(plt.MultipleLocator(xMax/SPHERE_MAJOR_TICKS_NUM ))
                ax.xaxis.set_minor_locator(plt.MultipleLocator(xMax/(SPHERE_MAJOR_TICKS_NUM  * 2)))
                ax.yaxis.set_major_locator(plt.MultipleLocator(yMax/SPHERE_MAJOR_TICKS_NUM ))
                ax.yaxis.set_minor_locator(plt.MultipleLocator(yMax/(SPHERE_MAJOR_TICKS_NUM  * 2)))
                ax.set_xlabel("x",labelpad = 8)
                ax.set_ylabel("y",labelpad = 8)
                ax.set_zlabel("z",labelpad = 10)
                changeTick = False
            # if (changeTitle == True):
            #     ax.set_title(titleName, fontsize = 10, pad = 20)
            #     changeTitle = False
    # plt.savefig(saveName)
    # figureI += 1
    # plotted += 1

def plotSphericalPointWrapper(dataDict, saveFotoPath, \
    config, threshold,\
        model = 'rayTrace'):
    """
    wrapper for the plotSphericalPoint

    """

    # first get the dictionary item
    if (config not in dataDict.keys()):
        raise Exception("The key {} not in dataDict fail to call plotSphericalPoint".format(config))
    
    responses = getItems(dataDict, config, threshold,model,
    queries=['sphericalPoint'])

def readConfig(configFilePath):
    """
    read in config files

    """
    #
    print ("to be programmed...")


def mainFunc(dataFolderPath,saveFotoPath, directReadFromFile = False, \
    saveDict = False, saveDictPath = "./dataDict.txt",\
    readDictPath = "./dataDict.txt",\
    generateAngularPower = False, \
        angularPowerFileName = "angularPower.txt", \
    estimator = '4',\
    models = ['rayTrace'],\
    thresholdValues = ['-75','-85','-90','-100'],\
    matchPatterns = ["TX.*RX[0-9]{2}"],\
    methods = ["correlation"],\
    antennaPattern = "hornAntenna",\
    singleThresholdPlot = True,\
    plot_HealpixMap = False,\
    plot_SphericalPoints = False):
    """ main function for plotting

    Parameters:
    models: list of strings
        either 'rayTrace', or 'meas'
    
    methods: list of strings
        either 'correlation' or 'angularPower' or 'angularPowerScaled'

    antennaPattern: str
        indicate the antenna pattern

    singleThresholdPlot: boolean
        if set true (default), will also plot with single threshold value
    
    plot_HealpixMap: boolean
        if set true, then will also plot the healpix map for each single configuration
        default False
    
    plot_SphericalPoint: boolean
        if set true, will also plot the spherical point for each single configuration,
        default False
    """

    for method in methods:
        for model in models:
            dataDict = getDataDictCompact(dataFolderPath,estimator=estimator, method = method,\
                directReadFromFile=directReadFromFile, saveDict=saveDict,\
                    saveDictPath=saveDictPath, readDictPath=readDictPath,\
                        generateAngularPower=generateAngularPower,\
                    angularPowerFileName=angularPowerFileName)
            
            plotFixedModelTypeDiffPos(dataDict, saveFotoPath, estimator=estimator,\
                model=model, thresholdValues=thresholdValues, matchPatterns=matchPatterns,\
                    method=method, singleThresholdPlot=singleThresholdPlot,\
                        plot_HealpixMap=plot_HealpixMap, plot_SphericalPoints=plot_SphericalPoints,
                        antennaPattern = antennaPattern)



def AllPlot(fileNamePrefix, stringFile, dataFilePrefix, healpixMapSuffix, saveNamePrefix,compactPlot = True,**kwargs):
    """ plot the four estimators, compact plot, the angular power spectrum and the healpix map

    Parameters:
    -------------
    fileNamePrefix: str
        prefix for all the data file in that directory, ususally the path to that dir

    stringFile: str
        prefix for the file that stores the name of the datafile, e.g stringName.txt

    dataFilePrefix: str
        the prefix for the data file which stores the correlation data like Pi_HealpixRandomMode

    healpixMapSuffix: str
        the suffix for the data that stores the map info, like Pi_HealpixAverageMap.txt

    saveNamePrefix: str
        the prefix for the saved plotted figure, e.g ../../photodata/6/18/
    """
    global figureI, plotted, cut_Val
    stringHealAvg = fileNamePrefix + stringFile
    nameHealAvg = []
    for line in open(stringHealAvg, 'r'):
        values = [s for s in line.split()]
        if (len(values) >= 1):
            nameHealAvg.append(values[0])
        
    if(ADD_TIME == True):
        saveNamePrefix += currentTime
    # figureSaveNamePre = saveNamePreix + 
    if (PLOT_RECENT == True):
        nameHealAvg = [nameHealAvg[-1]]

    multiRunPlot = False
    for strVar in nameHealAvg:
        if (SIGNAL_DATA == False):
            calculationInfo = fileNamePrefix + strVar + "calculationInfo.txt"
        else :
            calculationInfo = fileNamePrefix + strVar + "/calculationInfo.txt"
        info = []
    # try:
        for line in open(calculationInfo,'r'):
            values = [float(s) for s in line.split()]
            info.append(values[0])

        pointCount = int(info[0])
        randomCount = int(info[1])
        deltaTheta = info[2]
        average = int(info[3])
        spannedAngle = info[4]
        endTheta = info[5]
        if (Convert2Deg == True):
            cut_Val = spannedAngle * 180/np.pi
        else:
            cut_Val = spannedAngle
        kwd = {}
        saveNamePath = saveNamePrefix
        if (CREATE_DIR == True):
            saveNamePath += strVar + "/"
            if not os.path.exists(saveNamePath):
                os.makedirs(saveNamePath)
        dataPre = fileNamePrefix + strVar
        healpixMapFileName = dataPre + healpixMapSuffix
        healpixTitlePre = "healpixMap:"
        healpixSavePre = saveNamePath + "HealpixMap"
        # if (CUT_RANGE == True):
        #     saveNamePath += "cutRange"
        if (RANDOM_MAP_INFO == True):
            randomMapFile = dataPre + "randomMap.txt"
            randomSavePre = saveNamePath + "RandomMap"

        if (MATERN == True):
            indexBegin = 6
            _clusterCenterCount = int(info[indexBegin])
            _clusterTheta = info[indexBegin + 1]
            _clusterPointCount = int(info[indexBegin + 2])
            _totalPoint = int(info[indexBegin + 3])
            _thetaBegin = float(info[indexBegin + 4])
            _thetaEnd = float(info[indexBegin + 5])
            _dthetaBegin = float(info[indexBegin + 6])
            _dthetaEnd = float(info[indexBegin + 7])
            _useMathDouble = bool(int(info[indexBegin + 8]))
            _useOffsetPoisson = bool(int(info[indexBegin + 9]))
            _compactPlot = True
            kwd['clusterTheta'] = _clusterTheta
            kwd['clusterPointCount'] = _clusterPointCount
            kwd['totalPoint'] = _totalPoint
            mapTitlePre = "maternProcess"
            sphericalFile = dataPre + "spherePoint.txt"
            figureTitlePre = "spherical points"

            # plotEstimatorMaternAlgo(dataPre,namePrefix,saveNamePath,pointCount,randomCount,deltaTheta,average,spannedAngle,[figureI+1,figureI+2,figureI+3,figureI+4],_clusterCenterCount,
            #     _clusterTheta,_clusterPointCount,_totalPoint,_thetaBegin,_thetaEnd,_dthetaBegin, _dthetaEnd, _compactPlot, _useMathDouble, _useOffsetPoisson)
            
            # no need for total point, it's the same as pointCount
            if (Convert2Deg == True):
                _clusterTheta *= 180 / np.pi
                _dthetaEnd *= 180/np.pi
                _clusterTheta = float("{:.2f}".format(_clusterTheta))
                _dthetaEnd = float("{:.2f}".format(_dthetaEnd))
            plotEstimatorHealpixAverageAlgo(dataPre,namePrefix,saveNamePath,pointCount,randomCount,deltaTheta,average,_thetaBegin,_thetaEnd,[figureI+1,figureI+2,figureI+3,figureI+4],\
                compactPlot=_compactPlot,\
                clusterCenter=_clusterCenterCount,\
                clusterTheta = _clusterTheta, clusterPointCount = _clusterPointCount,\
                dTEnd = _dthetaEnd, data = "matern process") #dTBeg=_dthetaBegin, , #_useMathDouble, poissonOffset = _useOffsetPoisson, 
                        
            # printHealpixMap(healpixMapFileName,mapTitlePre,healpixSavePre,pointCount,spannedAngle, 
            #     clusterCenter = _clusterCenterCount, Theta = _clusterTheta,
            #     clusterPoints = _clusterPointCount, totalPoint = _totalPoint, thetaBegin = _thetaBegin, thetaEnd = _thetaEnd, daughterThetaBegin = _dthetaBegin,
            #     daughterThetaEnd = _dthetaEnd)
            # #Theta = _clusterTheta,
            # plotCl(healpixMapFileName,mapTitlePre,saveNamePath, _totalPoint, SIGNAL_DATA, MY_DOUBLE,spannedAngle,figureI, clusterCenter = _clusterCenterCount,  clusterPoints = _clusterPointCount, thetaBegin = _thetaBegin, thetaEnd = _thetaEnd)
            figureI += 5
            # the originial data has used  pointsCount to denote the parent process counts, not wanted
            # plotSphericalPoint(sphericalFile, figureTitlePre, saveNamePath,_totalPoint,spannedAngle,clusterCenter = _clusterCenterCount, dist = _clusterTheta, clusterPoints = _clusterPointCount, thetaEnd = _thetaEnd, clusterTheta= _dthetaEnd)

            # if (RANDOM_MAP_INFO == True):
            #     figureI += 5
            #     # printHealpixMap(randomMapFile,"randomly generated map", randomSavePre, pointCount,spannedAngle)
            #     plotDiffCl(healpixMapFileName, randomMapFile,"maternProcess", randomSavePre, _totalPoint,SIGNAL_DATA,MY_DOUBLE,spannedAngle, figureI)
            #     figureI += 5
            #     plotted += 5
            # daughterThetaBegin = _dthetaBegin,
        elif (HEALPIX_CORR == True):
            indexBegin = 6
            thetaBegin = info[indexBegin]
            thetaEnd = info[indexBegin + 1]
            plotEstimatorHealpixAverageAlgo(dataPre,dataFilePrefix,saveNamePath, pointCount,randomCount,deltaTheta,
                average,thetaBegin,thetaEnd,[figureI+1,figureI+2,figureI+3,figureI+4],compactPlot, data="generated poisson points")
            printHealpixMap(healpixMapFileName, healpixTitlePre,healpixSavePre,pointCount,spannedAngle)
            # plotCl(healpixMapFileName, "poissonPointProcess",healpixSavePre, pointCount, SIGNAL_DATA, MY_DOUBLE,spannedAngle, figureI)
            # figureI += 5
            # plotted += 4
            # plotSphericalPoint(dataPre + "spherePoint.txt", "spherical points ", saveNamePath,pointCount,
            # spannedAngle)

            # if (RANDOM_MAP_INFO == True):
            #     printHealpixMap(randomMapFile,"randomly generated map", randomSavePre,pointCount,spannedAngle)
            #     plotDiffCl(healpixMapFileName, randomMapFile, "poissonPointProcess", randomSavePre,pointCount,SIGNAL_DATA,MY_DOUBLE,spannedAngle, figureI)
            #     figureI += 5
            #     plotted += 5
        elif (HEALPIX_AVERAGE == True):
            indexBegin = 6
            thetaBegin = info[indexBegin]
            thetaEnd = info[indexBegin + 1]
            # phiBegin = info[indexBegin + 2]
            # phiEnd = info[indexBegin + 3]
            # plotEstimatorHealpixAverageAlgo(dataPre,dataFilePrefix,saveNamePath, pointCount,randomCount,deltaTheta,
            #     average,0.0,endTheta,[figureI+1,figureI+2,figureI+3,figureI+4],compactPlot, data="randomly generated healpix map")
            plotEstimatorHealpixAverageAlgo(dataPre,dataFilePrefix,saveNamePath, pointCount,randomCount,deltaTheta,
                average,thetaBegin,thetaEnd,[figureI+1,figureI+2,figureI+3,figureI+4],compactPlot, data="randomly generated healpix map")
            # printHealpixMap(healpixMapFileName, healpixTitlePre,healpixSavePre,pointCount,spannedAngle)
            # plotCl(healpixMapFileName, "healpixRandom",healpixSavePre, pointCount, SIGNAL_DATA, MY_DOUBLE,spannedAngle, figureI)
            # if (RANDOM_MAP_INFO == True):
            #     figureI += 5
            #     printHealpixMap(randomMapFile,"randomly generated map", randomSavePre,pointCount,spannedAngle)
            #     plotDiffCl(healpixMapFileName, randomMapFile, "poissonPointProcess", randomSavePre,pointCount,SIGNAL_DATA,MY_DOUBLE,spannedAngle, figureI)
            #     figureI += 5
            #     plotted += 5
            
        elif (SIGNAL_DATA == True):
            indexBegin = 6
            thetaBegin = info[indexBegin]
            thetaEnd = info[indexBegin + 1]
            threshold = info[indexBegin + 2]
            # factor = int(info[indexBegin + 3])
            dataPre += "/"
            healpixMapDataFile = dataPre + healpixMapSuffix
            randomMapDataFile = dataPre + "randomMap.txt"
            plotEstimatorHealpixAverageAlgo(dataPre,dataFilePrefix,saveNamePath, pointCount,randomCount,deltaTheta,
                average,thetaBegin,thetaEnd,[figureI+1,figureI+2,figureI+3,figureI+4],compactPlot, data="experiment data",
                threshold = str(threshold)+"dB")  #factor = factor
            # printHealpixMap(healpixMapDataFile, healpixTitlePre,healpixSavePre,pointCount,spannedAngle)
            # plotCl(healpixMapDataFile, "experiment data:",healpixSavePre, pointCount, SIGNAL_DATA, MY_DOUBLE,spannedAngle,figureI,threshold = threshold)
            # plotSphericalPoint(dataPre + "spherePoint.txt", "spherical points ", saveNamePath,pointCount,
            # spannedAngle)
            # if (RANDOM_MAP_INFO == True):
            #     figureI += 5
            #     # printHealpixMap(randomMapDataFile,"randomly generated map", randomSavePre,pointCount,spannedAngle)
            #     plotDiffCl(healpixMapDataFile, randomMapDataFile, "experiment data", randomSavePre,pointCount,SIGNAL_DATA,MY_DOUBLE,spannedAngle, figureI)
            #     # figureI += 5
            #     # plotted += 5


        else:            
            plotEstimatorHealpixAverageAlgo(dataPre,dataFilePrefix,saveNamePath, pointCount,randomCount,deltaTheta,
                average,spannedAngle,[figureI+1,figureI+2,figureI+3,figureI+4],compactPlot)
            printHealpixMap(healpixMapFileName, healpixTitlePre,healpixSavePre,pointCount,spannedAngle)
            plotCl(healpixMapFileName, healpixTitlePre,healpixSavePre, pointCount, SIGNAL_DATA, MY_DOUBLE,spannedAngle)
            plotSphericalPoint(dataPre + "spherePoint.txt", "spherical points ", saveNamePath,pointCount,
            spannedAngle)            
        figureI = figureI + 6
        plotted += 5
        if (MULTI_RUN == True and not multiRunPlot):
            saveMultiPath = saveNamePrefix + "/averageMulti/"
            if not os.path.exists(saveMultiPath):
                os.mkdir(saveMultiPath)
            fileMultiName = [fileNamePrefix + var for var in nameHealAvg]
            sweepMultiEstimator(fileMultiName,dataFilePrefix,saveMultiPath,pointCount,randomCount,deltaTheta,
            average,spannedAngle,[figureI+1,figureI+2,figureI+3,figureI+4])
            multiRunPlot = True

        if (NO_PLOT == True):
            plt.close('all')
            figureI = 0
            plotted = 0            


fileNamePrefix = "../data/" + datum
figureSaveNamePre = "../../photoData/" + datum
calculationInfo = fileNamePrefix + "calculationInfo.txt"
info = []
# for line in open(calculationInfo,'r'):
#     values = [float(s) for s in line.split()]
#     info.append(values[0])


# pointCount = int(info[0])
# randomCount = int(info[1])
# deltaTheta = info[2]
# average = int(info[3])
# spannedAngle = info[4]


# if(HEALPIX_CORR == True):
#     # fileNamePrefix += "greatAverage/poissonTest/"

#     fileNamePrefix += "mapVersion/rightAverage/poissonTest/" #"rightAverage/poissonTest/"
#     namePrefix = "poissonHealpixRandomMode"
#     comparePrefix = "compareHealpixRandomMode"
#     healpixMapFileName = "poissonHealpixRandomMap.txt"
#     healpixFigPre = "poissonHealpix"
#     healpixSavePre =  figureSaveNamePre  + "mapVersion/rightAverage/poissonTest/"#"rightAverage/poissonTest/"
#     #"greatAverage/poissonTest/"

#     healpixCompFileName = fileNamePrefix + "poissonHealpixRandomPixelMap.txt"
#     healpixCompPre = "pixMap_Healpix"
#     healpixComSavePre = figureSaveNamePre + "pixMap_Healpix"

#     AllPlot(fileNamePrefix,"stringName.txt",namePrefix,healpixMapFileName,healpixSavePre)
#     # plotEstimatorHealpixAverageAlgo(fileNamePrefix,namePrefix,figureSaveNamePre,pointCount,randomCount,deltaTheta,
#     #     average,spannedAngle,[figureI+1,figureI+2,figureI+3,figureI+4],plotTogether)
#     # printHealpixMap(healpixMapFileName,healpixFigPre,healpixSavePre,pointCount,spannedAngle)
#     # the comparison
#     # plotEstimatorHealpixAverageAlgo(fileNamePrefix,comparePrefix,figureSaveNamePre,pointCount,randomCount,deltaTheta,average,spannedAngle,[6,7,8,9])
#     # printHealpixMap(healpixCompFileName, healpixCompFileName, healpixComSavePre, pointCount, spannedAngle)

# elif (VEC_MAP == True):
#     # plotSphericalPoint("../data/"+datum + "vecMap/spherePoint.txt","vector","../../photoData/"+datum+"vecMap/",
#     #     866,100)
#     printHealpixMap("../data/"+datum + "map/phiConstraint/hpMap.txt","title","../../photoData/"+datum+"map/phiConstraint/",866,100 )
# elif(HEALPIX_AVERAGE == True):
#     # subPath = "wholeSphere/averageRR/eigenPixMult/nonCum/HealpixAverage/"
#     subPath = "wholeSphere/rightAverageRR/HeapixRandom/"  # phiConstraint/"
#     fileNamePrefix += subPath     # wholeSphere/rightAverageRR/HeapixRandom/" #"wholeSphere/greatAverageRR/HeapixRandom/"
#     namePrefix = "Pi_HealpixAverageRandomMode"
#     healpixMapFileName = "Pi_HealpixAverageRandomMap.txt"
#     healpixTitlePre = "healpix_Random:"
#     healpixSavePre = "../../photoData/" + datum + subPath # "wholeSphere/greatAverageRR/HeapixRandom/"

#     AllPlot(fileNamePrefix, "stringName.txt",namePrefix,healpixMapFileName,healpixSavePre)
#     # plotEstimatorHealpixAverageAlgo(fileNamePrefix,healpixAverageNamePre,figureSaveNamePre, pointCount,randomCount,deltaTheta,
#     #     average,spannedAngle,[figureI+1,figureI+2,figureI+3,figureI+4])
#     # printHealpixMap(healpixMapFileName, healpixTitlePre,healpixSavePre,pointCount,spannedAngle)
# elif(HEALPIX_ALL_ONE == True):
#     healpixAverageNamePre = "Pi_AllOneHealpixRandomMode"
#     healpixMapFileName = fileNamePrefix +  "Pi_AllOneHealpixRandomMap.txt"
#     healpixTitlePre = "AllOneHealpix_Random:"
#     healpixSavePre = figureSaveNamePre + "Pi_AllOneHealpixMap"
#     plotEstimatorHealpixAverageAlgo(fileNamePrefix,healpixAverageNamePre,figureSaveNamePre, pointCount,randomCount,deltaTheta,
#         average,spannedAngle,[figureI+1,figureI+2,figureI+3,figureI+4],compactPlot=True)
#     printHealpixMap(healpixMapFileName, healpixTitlePre,healpixSavePre,pointCount,spannedAngle)
#     plotCl(healpixMapFileName, healpixAverageNamePre,figureSaveNamePre, pointCount, SIGNAL_DATA, MY_DOUBLE,spannedAngle)
#     plotted += 4
#     figureI += 4
# elif(COMPACT_TEST == True):
#     allOneMapPre = "nonCum/AllOneMap/"
#     # "wholeSphere/averageRR/eigenPixMult/nonCum/geomSpace/HealpixAverage/"
#     healAvgPre = "wholeSphere/rightAverageRR/HeapixRandom/"
#     stringFilePre = "stringName.txt"
#     # all one map
#     if (ALL_ONE_PLOT == True):
#         allOneStr = fileNamePrefix + allOneMapPre
#         stringAllName = allOneStr + stringFilePre
#         nameAllOne = []
#         for line in open(stringAllName,'r'):
#             values = [s for s in line.split()]
#             nameAllOne.append(values[0])
        
#         # count = 0
#         healpixAverageNamePre = "Pi_AllOneHealpixRandomMode"
#         # nameAllOne = nameAllOne[5:]
#         figureSaveNamePreAllOne = figureSaveNamePre + allOneMapPre
#         for strVar in nameAllOne:

#             calculationInfo = allOneStr + strVar + "calculationInfo.txt"
#             info = []
#             try:
#                 for line in open(calculationInfo,'r'):
#                     values = [float(s) for s in line.split()]
#                     info.append(values[0])

#                 pointCount = int(info[0])
#                 randomCount = int(info[1])
#                 deltaTheta = info[2]
#                 average = int(info[3])
#                 spannedAngle = info[4]
#                 endTheta = info[5]

#                 healpixMapDataPre = allOneStr + strVar
#                 healpixMapFileName = healpixMapDataPre+  "Pi_AllOneHealpixRandomMap.txt"
#                 healpixTitlePre = "AllOneHealpix_Random:"
#                 healpixSavePre = healpixMapDataPre + "Pi_AllOneHealpixMap"

#                 plotEstimatorHealpixAverageAlgo(healpixMapDataPre,healpixAverageNamePre,figureSaveNamePreAllOne, pointCount,randomCount,deltaTheta,
#                     average,spannedAngle,[figureI+1,figureI+2,figureI+3,figureI+4],compactPlot=True)
#                 printHealpixMap(healpixMapFileName, healpixTitlePre,figureSaveNamePreAllOne,pointCount,spannedAngle)
#                 plotCl(healpixMapFileName, healpixAverageNamePre,figureSaveNamePreAllOne, pointCount, SIGNAL_DATA, MY_DOUBLE,spannedAngle)
#                 figureI = figureI + 6
#                 plotted += 4
#                 if (NO_PLOT == True):
#                     plt.close('all')
#                     figureI = 0
#                     plotted = 0
#             except TclError:
#                 print("dunno")
#             except IOError:
#                 print("in Compact Plot fileName starts with: " + strVar + " not Found")
             
#     # healAvg

#     if (HEAL_AVG_PLOT == True):
#         # healAvgStr = fileNamePrefix + healAvgPre
#         # stringHealAvg = healAvgStr + stringFilePre
#         # nameHealAvg = []
#         # for line in open(stringHealAvg, 'r'):
#         #     values = [s for s in line.split()]
#         #     nameHealAvg.append(values[0])
            
#         # # nameHealAvg = nameHealAvg[5:10]
#         # healpixAverageNamePre = "Pi_HealpixAverageRandomMode"    
#         # figureSaveNamePreHealAvg = figureSaveNamePre + healAvgPre
#         # for strVar in nameHealAvg:
#         #     calculationInfo = healAvgStr + strVar + "calculationInfo.txt"
#         #     info = []
#         # # try:
#         #     for line in open(calculationInfo,'r'):
#         #         values = [float(s) for s in line.split()]
#         #         info.append(values[0])

#         #     pointCount = int(info[0])
#         #     randomCount = int(info[1])
#         #     deltaTheta = info[2]
#         #     average = int(info[3])
#         #     spannedAngle = info[4]
#         #     endTheta = info[5]

#         #     healpixMapDataPre = healAvgStr + strVar
#         #     healpixMapFileName = healpixMapDataPre + "Pi_HealpixAverageRandomMap.txt"
#         #     healpixTitlePre = "healpix_Random:"
#         #     healpixSavePre = healpixMapDataPre + "Pi_HealpixAverageRandomMap"

#         #     plotEstimatorHealpixAverageAlgo(healpixMapDataPre,healpixAverageNamePre,figureSaveNamePreHealAvg, pointCount,randomCount,deltaTheta,
#         #         average,spannedAngle,[figureI+1,figureI+2,figureI+3,figureI+4])
#         #     printHealpixMap(healpixMapFileName, healpixTitlePre,figureSaveNamePreHealAvg,pointCount,spannedAngle)
#         #     plotCl(healpixMapFileName, healpixAverageNamePre,figureSaveNamePreHealAvg, pointCount, SIGNAL_DATA, MY_DOUBLE,spannedAngle)
#         #     figureI = figureI + 6
#         #     plotted += 4
#         #     if (NO_PLOT == True):
#         #         plt.close('all')
#         #         figureI = 0
#         #         plotted = 0    

#         AllPlot(fileNamePrefix + healAvgPre, stringFilePre, "Pi_HealpixAverageRandomMode","Pi_HealpixAverageRandomMap.txt",
#         figureSaveNamePre + healAvgPre)

# # except IOError:
#         #     print("in Compact plot fileName starts with: " + strVar + " not Found")

# elif(MATERN == True):
#     fileNamePrefix += "matern/"

#     subPath = "rightAverage/matern/linSpace/noPoissonOffset/"
#     namePrefix = "poissonHealpixRandomMode"
#     comparePrefix = "compareHealpixRandomMode"
#     healpixMapFileName = "poissonHealpixRandomMap.txt"
#     healpixFigPre = "poissonHealpix"
#     figureSaveNamePre += "matern/"
#     healpixSavePre = figureSaveNamePre + "poissonHealpix"


#     # "../../photoData/" + datum   "greatAverage/matern/linSpace/noPoissonOffset/" "greatAverage/matern/linSpace/withPoissonOffset/"
#     AllPlot("../data/" + datum + subPath,"stringName.txt",namePrefix,healpixMapFileName,"../../photoData/" + datum + subPath)
#     # plotEstimatorMaternAlgo(fileNamePrefix,namePrefix,figureSaveNamePre,pointCount,randomCount,deltaTheta,average,spannedAngle,[figureI+1,figureI+2,figureI+3,figureI+4],_clusterCenterCount,
#     #     _clusterTheta,_clusterPointCount,_totalPoint,_compactPlot, _useMathDouble, _useOffsetPoisson)
#     # printHealpixMap(healpixMapFileName,healpixFigPre,healpixSavePre,pointCount,spannedAngle, 
#     #     clusterCenter = _clusterCenterCount, Theta = _clusterTheta,
#     #     points = _clusterPointCount)
#     # plotCl(healpixMapFileName, healpixAverageNamePre,figureSaveNamePre, pointCount, SIGNAL_DATA, MY_DOUBLE,spannedAngle)
# elif(SIGNAL_DATA == True):
#     # thres = float(info[5])
#     subPath = "dataTest/linSpace/91Threshold/MAXthetaFinal/"#"dataTest/linSpace/91Threshold/MAXthetaFinal/"
#     namePrefix = "SignalAverageHealpixMode"
#     healpixMapFileName = "SignalAverageHealpixMap.txt"
#     healpixTitlePre = "signalData:"
#     healpixSavePre = figureSaveNamePre + "SignalAverageHealpixMap"

#     AllPlot("../data/" + datum + subPath, "stringName.txt", namePrefix, healpixMapFileName, "../../photoData/" + datum + subPath)
#     # plotEstimatorHealpixAverageAlgo(fileNamePrefix,healpixAverageNamePre,figureSaveNamePre, pointCount,randomCount,
#     #     deltaTheta,average,spannedAngle,[figureI+1,figureI+2,figureI+3,figureI+4],compactPlot = True, threshold = thres)
#     # printHealpixMap(healpixMapFileName, healpixTitlePre,healpixSavePre,pointCount,spannedAngle,threshold = thres)
#     # plotCl(healpixMapFileName, healpixAverageNamePre,figureSaveNamePre, pointCount, SIGNAL_DATA, MY_DOUBLE,spannedAngle)
# elif (SPHERICAL == True):
#     # thres = float(info[5])
#     healpixAverageNamePre = "Pi_AverageHealpixMode"
#     healpixMapFileName = fileNamePrefix +  "Pi_AverageHealpixRandomMap.txt"
#     healpixTitlePre = "signalData:"
#     healpixSavePre = figureSaveNamePre + "Pi_AverageHealpixRandomMap"
#     # plotEstimatorHealpixAverageAlgo(fileNamePrefix,healpixAverageNamePre,figureSaveNamePre, pointCount,randomCount,
#     #     deltaTheta,average,spannedAngle,[figureI+1,figureI+2,figureI+3,figureI+4],compactPlot = True, threshold = thres)
#     printHealpixMap(healpixMapFileName, healpixTitlePre,healpixSavePre,pointCount,spannedAngle)
#     plotCl(healpixMapFileName, healpixAverageNamePre,figureSaveNamePre, pointCount, SIGNAL_DATA, MY_DOUBLE,spannedAngle)

# if (not NO_PLOT and plotted <= plotNumThreshold) :
#     plt.show()


# plot2Func("../data/"+datum + "rawTest/yes/landy.txt", "../data/"+datum + "rawTest/no/landy.txt",\
# "averaged map", "non-averaged","../../photoData/" + datum + "rawTest/compare.png")

# tempPath = "/home/panwei/Desktop/bachelorThese/photoData/7/3/wholeSphere/averageRR/eigenPixMult/nonCum/HealpixAverage/points1000avgPoints:2000avgTime:20dT0.01000_sp3.14159_endT:3.14159/runNum:4"

# fP = "../data/" + datum + "dataTest/linSpace/"
# f104Path = "104Threshold/MAXthetaFinal/"
# f91Path = "91Threshold/MAXthetaFinal/"
# fileN1 = "mode:1/points306avgTime:10nr:612thd:-104dTtB:0.523599tE:1.57080.01000_endT:2.74889ft:1"
# fileN2 = "mode:1/points219avgTime:10nr:438thd:-95dTtB:0.523599tE:1.57080.01000_endT:2.74889ft:1"

# fileM0 = "mode:1/points151avgTime:10nr:302thd:-89dTtB:0.523599tE:1.57080.01000_endT:2.74889"
# fileM1 = "mode:1/points101avgTime:10nr:202thd:-85dTtB:0.523599tE:1.57080.01000_endT:2.74889"
# fileM2 = "mode:1/points62avgTime:10nr:124thd:-80dTtB:0.523599tE:1.57080.01000_endT:2.74889"
# fileM3 = "mode:1/points27avgTime:10nr:54thd:-75dTtB:0.523599tE:1.57080.01000_endT:2.74889ft:1"
# file1 = fP + f104Path + fileN1 + "/SignalAverageHealpixMode4.txt"
# file2 = fP + f104Path + fileN2 + "/SignalAverageHealpixMode4.txt"
# file4 = fP + f91Path + fileM1 + "/SignalAverageHealpixMode4.txt"
# file5 = fP + f91Path + fileM2 + "/SignalAverageHealpixMode4.txt"
# file6 = fP + f91Path + fileM3 + "/SignalAverageHealpixMode4.txt"
# file3 = fP + f91Path + fileM0 + "/SignalAverageHealpixMode4.txt"
# savePathMulti = "../../photoData/" + datum + "dataTest/linSpace/" + "combi0.png"
# plotMultiFunc([file1,file2,file3,file4,file5,file6],["-104dB","-95dB","-89dB","-85dB","-80dB","-75dB"], "comparison between different\nthresholds",savePathMulti)


if __name__ == "__main__":
    # plotWholeRayMeas("/home/panwei/Desktop/bachelorThese/myCode/data/10/27/dataTest/LangenfeldData",method = "angularPower")
    # plotWholeRayMeas("/home/panwei/Desktop/bachelorThese/myCode/data/10/27/dataTest/testFolder")
    # dataDict = getDataDictCompact("/home/panwei/Desktop/bachelorThese/myCode/data/10/27/dataTest/LangenfeldData")

    # rawDataDir = "/home/panwei/Desktop/bachelorThese/myCode/data/3/1/dataTest/superCData"
    # rawDataDir = "/home/panwei/Desktop/bachelorThese/myCode/data/3/31/dataTest/frankfurt"
    # rawDataDir = "/home/panwei/Desktop/bachelorThese/myCode/data/3/31/dataTest/seoul/hornAntenna"
    # savedFotoPath = "/home/panwei/Desktop/bachelorThese/photoData/seoul/hornAntenna"

    # rawDataDir = "/home/panwei/Desktop/bachelorThese/myCode/data/3/31/dataTest/seoul/phasedAntennaArray"
    # savedFotoPath = "/home/panwei/Desktop/bachelorThese/photoData/seoul/phasedAntennaArray"
    # rawDataDir = "/home/panwei/Desktop/bachelorThese/myCode/data/3/31/dataTest/seoul/hornAntennaAggregate"
    # savedFotoPath = "/home/panwei/Desktop/bachelorThese/photoData/seoul/hornAntennaAggregate"
    
    # rawDataDir = "/home/panwei/Desktop/bachelorThese/myCode/data/3/31/dataTest/seoul/phasedAntennaArrayAggregate"
    # savedFotoPath = "/home/panwei/Desktop/bachelorThese/photoData/seoul/phasedAntennaArrayAggregate"
    

    
    # rawDataDir = "/home/panwei/Desktop/bachelorThese/myCode/data/3/31/dataTest/frankfurt/phasedAntennaArray"
    # savedFotoPath = "/home/panwei/Desktop/bachelorThese/photoData/frankfurt/phasedAntennaArray"

    # rawDataDir = "/home/panwei/Desktop/bachelorThese/myCode/data/3/31/dataTest/frankfurt/hornAntennaAggregate"
    # savedFotoPath = "/home/panwei/Desktop/bachelorThese/photoData/frankfurt/hornAntennaAggregate"

    # rawDataDir = "/home/panwei/Desktop/bachelorThese/myCode/data/3/31/dataTest/frankfurt/phasedAntennaArrayNewMode"
    # savedFotoPath = "/home/panwei/Desktop/bachelorThese/photoData/frankfurt/phasedAntennaArrayNewMode"

    # rawDataDir = "/home/panwei/Desktop/bachelorThese/myCode/data/3/31/dataTest/frankfurt/phasedAntennaArrayAggregateNewMode"
    # savedFotoPath = "/home/panwei/Desktop/bachelorThese/photoData/frankfurt/phasedAntennaArrayAggregateNewMode"

    rawDataDir = "/home/panwei/Desktop/bachelorThese/myCode/data/3/31/dataTest/seoul/phasedAntennaArrayNewMode"
    savedFotoPath = "/home/panwei/Desktop/bachelorThese/photoData/seoul/phasedAntennaArrayNewMode"

    # rawDataDir = "/home/panwei/Desktop/bachelorThese/myCode/data/3/31/dataTest/seoul/phasedAntennaArrayAggregateNewMode"
    # savedFotoPath = "/home/panwei/Desktop/bachelorThese/photoData/seoul/phasedAntennaArrayAggregateNewMode"



    rawDataDir = "/home/panwei/Desktop/bachelorThese/myCode/data/3/31/dataTest/frankfurt/hornAntennaNewMode"
    savedFotoPath = "/home/panwei/Desktop/bachelorThese/photoData/frankfurt/hornAntennaNewMode"

    # rawDataDir = "/home/panwei/Desktop/bachelorThese/myCode/data/3/31/dataTest/frankfurt/hornAntennaAggregateNewMode"
    # savedFotoPath = "/home/panwei/Desktop/bachelorThese/photoData/frankfurt/hornAntennaAggregateNewMode"

    # rawDataDir = "/home/panwei/Desktop/bachelorThese/myCode/data/3/31/dataTest/seoul/hornAntennaNewMode"
    # savedFotoPath = "/home/panwei/Desktop/bachelorThese/photoData/seoul/hornAntennaNewMode"

    # rawDataDir = "/home/panwei/Desktop/bachelorThese/myCode/data/3/31/dataTest/seoul/hornAntennaAggregateNewMode"
    # savedFotoPath = "/home/panwei/Desktop/bachelorThese/photoData/seoul/hornAntennaAggregateNewMode"




    antennaPattern = "hornAntenna"
    # antennaPattern = "phasedAntennaArray"

    # dataDict = getDataDictCompact(rawDataDir,  \
    #     method = "correlation", \
    #         saveDict=True)
    # dataDict = getDataDictCompact(dataDictSavedPath,directReadFromFile=True)

    # print(dataDict)

    # matchPatterns = ["TX[1-2]_RX1(1|2)$", "TX1_RX2[0-9]$", "TX1_RX[A-Z]$","TX1_RX(I|V|X){1,2,3}$","TX2_RX[A-Z]$","TX2_RX[a-g].*$"]

    # matchPatterns = ["TX_4652916_RX_.*$", "TX_[0-9]*_RX_355_355", "TX_4[0-9]*_RX_.*", "TX_5[0-9]*_RX_.*", "TX_1[0-9]*_RX_.*", "TX_1852836_RX.*","TX_[0-9]*_RX_380_400", "TX_32[0-9]*_RX_.*", "TX_6[0-9]*_RX_.*"]
    # frankfurt
    # matchPatterns = ["TX_32[0-9]*_RX_355_355", 
    #     "TX_4(3|4)[0-9]*_RX_300_.*",
    #     "TX_43[0-9]*_RX_(300|380)_.*", 
    #     "TX_43[0-9]*_RX_(380|400)_.*", 
    #     "TX_5(6|7)[0-9]*_RX_600.*",
    #     "TX_5(8|9)[0-9]*_RX_600.*", 
    #     "TX_1[0-9]*_RX_(120|300).*", 
    #     "TX_15[0-9]*_RX_(300|380).*", 
    #     "TX_15[0-9]*_RX_(380|400).*", 
    #     "TX_1[0-9]*_RX_380_400",
    #     "TX_2[0-9]*_RX_380_400",
    #     "TX_3[0-9]*_RX_380_400",
    #     "TX_4[0-9]*_RX_380_400",
    #     "TX_4[0-9]*_RX_380_400"]

    # # aggregate
    # matchPatterns = ["RX_(120|300).*",
    # "RX_120.*",
    # "RX_300.*",
    # "RX_(300|380).*",
    # "RX_(300|355).*",
    # "RX_(355|380).*",
    # "RX_(380|400).*",
    # "RX_(400|600).*"]


    # seoul
    matchPatterns = ["TX_32[0-9]*_RX_355_355", 
        "TX_4(3|4)[0-9]*_RX_.*", 
        "TX_46[0-9]*_RX_.*", 
        "TX_5(6|7)[0-9]*_RX_600.*",
        "TX_5(8|9)[0-9]*_RX_600.*", 
        "TX_1[0-9]*_RX_.*", 
        "TX_1852836_RX.*",
        "TX_1[0-9]*_RX_380_400",
        "TX_2[0-9]*_RX_380_400",
        "TX_32[0-9]*_RX_380_400",
        "TX_33[0-9]*_RX_380_400",
        "TX_35[0-9]*_RX_380_400",
        "TX_4[0-9]*_RX_380_400",
        "TX_4[0-9]*_RX_380_400", 
        "TX_3271906_RX_(120|355).*", 
        "TX_3271906_RX_(120|380).*", 
        "TX_3271906_RX_(380|400).*", 
        "TX_6[0-9]*_RX_.*"]
    # matchPatterns = ["TX_4652916_RX_(380_400|400_380)", "TX_4652916_RX_380_400", "TX_4652916_RX_400_380"]
    # matchPatterns = ["TX_4652916_RX_380_400", "TX_4652916_RX_400_380"]

    # thresholdValues = ['-85','-75','-95']
    thresholdValues = ['-90','-95']
    # thresholdValues = ['-85','-80']
    # thresholdValues = ["-75","-80"]


    # plotFixedModelTypeDiffPos(dataDict, "/home/panwei/Desktop/bachelorThese/photoData/superCData",method = "correlation", thresholdValues= ["-100"] ,matchPattern="TX[1-2]_RX1(1|2)$")
    # plotFixedModelTypeDiffPos(dataDict, "/home/panwei/Desktop/bachelorThese/photoData/superCData",method = "correlation", thresholdValues= thresholdValues ,matchPatterns = matchPatterns)
    # savedFotoPath = "/home/panwei/Desktop/bachelorThese/photoData/frankfurt"

    # plotFixedModelTypeDiffPos(dataDict, savedFotoPath ,\
    #     model = "rayTrace", \
    #         thresholdValues= thresholdValues ,\
    #     matchPatterns = matchPatterns,\
    #          method = "correlation")
    # plotFixedModelTypeDiffPos(dataDict, savedFotoPath ,method = "angularPower", thresholdValues= thresholdValues ,matchPatterns = matchPatterns)

    mainFunc(rawDataDir,savedFotoPath,directReadFromFile=False, models = ['rayTrace'], \
        thresholdValues=thresholdValues, \
            saveDict=True,\
            matchPatterns=matchPatterns,\
        methods=['correlation','angularPower','angularPowerScaled'], plot_HealpixMap=False, generateAngularPower=False,singleThresholdPlot=True,
        antennaPattern = antennaPattern)

    # plotFixConfigDiffThreshold(dataDict,"/home/panwei/Desktop/bachelorThese/photoData/Langenfeld",method = "angularPower")
    # plotFixModelDBMDiffLocation(dataDict, "/home/panwei/Desktop/bachelorThese/photoData/Langenfeld", method = "angularPower")




    # hornAntenna/angularPower/pdf hornAntenna/angularPowerScaled/pdf hornAntenna/correlation/pdf hornAntennaAggregate/angularPower/pdf hornAntennaAggregate/angularPowerScaled/pdf hornAntennaAggregate/correlation/pdf phasedAntennaArray/angularPower/pdf phasedAntennaArray/angularPowerScaled/pdf phasedAntennaArray/correlation/pdf 
