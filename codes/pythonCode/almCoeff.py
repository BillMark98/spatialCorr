import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
import os
from plotHelpFunc import *
# from healpixPlot import plotted,figureI
CL_MAJOR_TICKS_NUM = 8

figX = 8
figY = 5


XThres = 75
def setYlim(plt, data):
    minVal = min(data)
    maxVal = max(data)
    plt.ylim(minVal, maxVal + 0.5 * maxVal)

def annot_max(x,y, ax=None):
    xmax = x[np.argmax(y)]
    ymax = y.max()
    # text= "x={:.3f}, y={:.3f}".format(xmax, ymax)
    # if not ax:
    #     ax=plt.gca()
    # bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
    # arrowprops=dict(arrowstyle="->",connectionstyle="angle,angleA=0,angleB=60")
    # kw = dict(xycoords='data',textcoords="axes fraction",
    #           arrowprops=arrowprops, bbox=bbox_props, ha="right", va="top")
    # ax.annotate(text, xy=(xmax, ymax), xytext=(0.94,0.96), **kw)

def annotPoint(x,y,ax):
    [bottom,top] = ax.get_ylim()
    diff = (top-bottom)/5
    for xp,yp in zip(x,y):
        text = "l={:d}".format(xp)
        xT = xp + 15
        yT = yp + diff
        if (yT > top):
            yT = yp + top/4
            xT = xT + 30
        bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
        arrowprops = dict(arrowstyle = "->", connectionstyle="arc3") #facecolor = 'black', shrink = 0.05, 
        kw = dict(xycoords = 'data',textcoords = 'data',arrowprops = arrowprops,bbox = bbox_props,ha="center",va="center")
        ax.annotate(text,xy=(xp,yp), xytext=(xT, yT), **kw)
    

def localMax(x,y,ax = None):
    global XThres
    indexX = []
    yThres = 0
    upBound = min(len(y),XThres)
    [bottom,top] = ax.get_ylim()
    ycut = (top - bottom)/(100)
    zeroHit = False
    indexMax = 0
    for i in range(1, upBound):
        if (y[i] >= y[i-1] and y[i] >= y[i+1] and (y[i] >= yThres or (zeroHit == True and (y[i] >= yThres/3 or indexMax <= 10)))):
            if(i >= 2):
                if(y[i] >= y[i-2] and y[i] >= y[i+2]):
                    indexX.append(i)
                    yThres = y[i]
                    indexMax = i
                    zeroHit = False
            else:
                indexX.append(i)
                yThres = y[i]
                zeroHit = False                
        if (y[i] < ycut):
            zeroHit = True
    annotPoint(x[indexX],y[indexX],ax)
    return indexX

def plotCl(healpixMapFileName,figureTitleNamePre,
    saveNamePre,pointCount, 
    SIGNAL_DATA, MY_DOUBLE, 
    spannedAngle,fignum,plotting = True, **kwargs):
    """ plot the angular power spectrum coeff

    SIGNAL_DATA, MY_DOUBLE no usage

    Parameters:
    --------
    healpixMapFileName: str
        name of the data file

    figureTitleNamePre: str
        prefix of the title
        
    saveNamePre: str
        prefix of the to be saved figures
    
    return:
        dictionary, {'Nside', 'cl','pointsCount'}
    """
    global figX, figY
    if (plotting == True):
        corrName = "spanned:%.4f"%(spannedAngle) + " points:" + str(pointCount)

        additionName = ""
        titleAddition = ""
        count = 0
        for key,value in kwargs.items():
            additionName += key + ":" + str(value)
            titleAddition += key + ":" + str(value) + ", "
            count += 1
            if (count % 3 == 0) :
                titleAddition+= "\n"
        
        figureName = figureTitleNamePre + ' ' + corrName + "\n" + titleAddition
        saveName =   saveNamePre + corrName + additionName
    hp1 = []
    for line in open(healpixMapFileName,'r'):   
        # if(SIGNAL_DATA == False or MY_DOUBLE == False):
        #     values = [float(s) for s in line.split()]
        # else:
        values = [float(s) for s in line.split()]
        hp1.append(values[0])

    hpOrder = int(hp1[0])
    hpPointCount = int(hp1[1])
    pixVal = hp1[2:]
    pixValLen = len(pixVal)
    Nside = 2**hpOrder
    Npixel = hp.nside2npix(Nside)
    pixVal.extend([0] * (int(Npixel) - pixValLen))
    pixVal = np.array(pixVal)
    if(hpPointCount <= 1):
        return None
    LMAX = Nside * 2
    cl = hp.anafast(pixVal,lmax = LMAX)
    alm = hp.map2alm(pixVal,lmax = LMAX)
    ell = np.arange(len(cl))
    coeff = (cl * (2 * ell + 1))/((hpPointCount)**2) * 4 * np.pi - 1/( hpPointCount)

    coeff2 = (cl * (2 * ell + 1)) / ((hpPointCount)**2)

    if (plotting == True):
        xMax = LMAX

        plt.figure(figsize = (figX, figY))
        fig1,ax1 = plt.subplots(num = fignum + 1)
        plt.xlim(0,LMAX)

        ax1.xaxis.set_major_locator(plt.MultipleLocator(xMax/CL_MAJOR_TICKS_NUM))
        ax1.xaxis.set_minor_locator(plt.MultipleLocator(xMax/(CL_MAJOR_TICKS_NUM * 2)))
        ax1.plot(ell, ell * (ell + 1) * cl)
        plt.xlabel("$\ell$")
        plt.ylabel("$\ell(\ell+1)C_{\ell}$")
        plt.grid()
        saveName1 = saveName + "LMAX:" + str(LMAX) + "CLcoeff_llp1.png"
        figureName1 = figureName + "LMAX:" + str(LMAX) +  " l(l+1)_normalized"
        plt.title(figureName1,fontsize = 8, y = 0.99)
        plt.savefig(saveName1)
        # print(alm[0])

        plt.figure(figsize= (figX, figY))
        fig2,ax2 = plt.subplots(num = fignum + 2)
        ax2.xaxis.set_major_locator(plt.MultipleLocator(xMax/CL_MAJOR_TICKS_NUM))
        ax2.xaxis.set_minor_locator(plt.MultipleLocator(xMax/(CL_MAJOR_TICKS_NUM * 2)))

        normedCl = (2 * ell + 1) * cl
        ax2.plot(ell, normedCl)
        plt.xlabel("$\ell$")
        plt.ylabel("$(2\ell+1)C_{\ell}$")
        plt.grid()
        saveName2= saveName + "LMAX:" + str(LMAX) + "CLcoeff_2lp1.png"
        figureName2 = figureName + "LMAX:" + str(LMAX) +  " (2l+1)_normalized"
        plt.title(figureName2,fontsize = 8, y = 0.99)

        #find local max
        indexX = localMax(ell,normedCl,ax2)
        maxXVal = ell[indexX]
        maxYVal = normedCl[indexX]
        plt.stem(maxXVal, maxYVal, linefmt = 'b--', markerfmt = 'ro', basefmt = ' ')
        plt.savefig(saveName2)



        plt.figure(figsize= (figX, figY))
        fig3,ax3 = plt.subplots(num = fignum + 3)
        ax3.xaxis.set_major_locator(plt.MultipleLocator(xMax/CL_MAJOR_TICKS_NUM))
        ax3.xaxis.set_minor_locator(plt.MultipleLocator(xMax/(CL_MAJOR_TICKS_NUM * 2)))
        ax3.plot(ell, coeff)
        plt.xlabel("$\ell$")
        plt.ylabel("$\widehat{c_{\ell}}$")
        plt.grid()
        saveName3= saveName + "LMAX:" + str(LMAX) + "CLcoeff_coeff.png"
        figureName3 = figureName + "LMAX:" + str(LMAX) +  " legendre coeff"
        plt.title(figureName3,fontsize = 8, y = 0.99)
        plt.savefig(saveName3)



        plt.figure(figsize= (figX, figY))
        fig4,ax4 = plt.subplots(num = fignum + 4)
        ax4.xaxis.set_major_locator(plt.MultipleLocator(xMax/CL_MAJOR_TICKS_NUM))
        ax4.xaxis.set_minor_locator(plt.MultipleLocator(xMax/(CL_MAJOR_TICKS_NUM * 2)))
        ax4.plot(ell, cl)
        plt.xlabel("$\ell$")
        plt.ylabel("$C_{\ell}$")
        plt.grid()
        saveName4= saveName + "LMAX:" + str(LMAX) + "CLcoeff_rawCL.png"
        figureName4 = figureName + "LMAX:" + str(LMAX) +  " Cl"
        plt.title(figureName4,fontsize = 8, y = 0.99)
        plt.savefig(saveName4)
    else:
        return {'Nside': Nside, 'cl':cl, 'pointsCount':hpPointCount}

def plotDiffCl(healpixMapFileName, randomMapFileName,figureTitleNamePre,saveNamePre,pointCount, SIGNAL_DATA, MY_DOUBLE, spannedAngle,fignum,**kwargs):
    """ plot the angular power spectrum coeff subtracted from a randomly generated map

    Parameters:
    --------
    healpixMapFileName: str
        name of the data file

    randomMapFileName: str
        name of the random Map file

    figureTitleNamePre: str
        prefix of the title
        
    saveNamePre: str
        prefix of the to be saved figures
    """
    global figX, figY
    corrName = "spanned:%.4f"%(spannedAngle) + " points:" + str(pointCount)

    additionName = ""
    titleAddition = ""
    count = 0
    for key,value in kwargs.items():
        additionName += key + ":" + str(value)
        titleAddition += key + ":" + str(value) + ", "
        count += 1
        if (count % 3 == 0) :
            titleAddition+= "\n"
    
    figureName = figureTitleNamePre + ' ' + corrName + "\n" + titleAddition
    saveName =   saveNamePre + corrName + additionName
    hp1 = []
    for line in open(healpixMapFileName,'r'):   
        if(SIGNAL_DATA == False or MY_DOUBLE == False):
            values = [float(s) for s in line.split()]
        else:
            values = [float(s) for s in line.split()]
        if(len(values) >= 1):
            hp1.append(values[0])

    hpOrder = int(hp1[0])
    hpPointCount = int(hp1[1])
    if(hpPointCount <= 1):
        print("healpix map empty, fileName: " + healpixMapFileName)
        return
    pixVal = hp1[2:]
    pixValLen = len(pixVal)
    Nside = 2**hpOrder
    Npixel = hp.nside2npix(Nside)
    pixVal.extend([0] * (int(Npixel) - pixValLen))
    pixVal = np.array(pixVal)

    LMAX = Nside * 2

    # read in random Map
    hp2 = []
    for line in open(randomMapFileName, 'r'):
        values = [float(s) for s in line.split()]
        if(len(values) >= 1):
            hp2.append(values[0])
    hp2Order = int(hp2[0])
    if (not hpOrder == hp2Order):
        print("two maps not of the same order")
        return
    hp2PointCount = int(hp2[1])
    if(hp2PointCount <= 1):
        print("random map empty, fileName: " + randomMapFileName)
        return
    pix2Val = hp2[2:]
    pix2ValLen = len(pix2Val)
    pix2Val.extend([0] * (int(Npixel) - pix2ValLen))
    pix2Val = np.array(pix2Val)


    cl = hp.anafast(pixVal,lmax = LMAX)
    alm = hp.map2alm(pixVal,lmax = LMAX)
    ell = np.arange(len(cl))
    # coeff = (cl * (2 * ell + 1))/((hpPointCount)**2) * 4 * np.pi - 1/( hpPointCount)

    # coeff2 = (cl * (2 * ell + 1)) / ((hpPointCount)**2)

    cl2 = hp.anafast(pix2Val,lmax = LMAX)
    # ell2 = np.arange(len(cl2))
    # coeffRandom = (cl * (2*ell + 1))/((hp2PointCount)**2)*
    if(len(cl) != len(cl2)):
        print("len(cl) !=  len(cl2)")
        return

    clFinal = cl - cl2
    xMax = LMAX


    plt.figure(figsize = (figX, figY))
    fig1,ax1 = plt.subplots(num = fignum + 1)
    ax1.xaxis.set_major_locator(plt.MultipleLocator(xMax/CL_MAJOR_TICKS_NUM))
    ax1.xaxis.set_minor_locator(plt.MultipleLocator(xMax/(CL_MAJOR_TICKS_NUM * 2)))
    ax1.plot(ell, ell * (ell + 1) * clFinal)
    plt.xlabel("$\ell$")
    plt.ylabel("$\ell(\ell+1)C_{\ell}$")
    plt.grid()

    annot_max(ell, ell * (ell + 1) * clFinal, ax1)
    saveName1 = saveName + "Diff_LMAX:" + str(LMAX) + "CLcoeff_llp1.png"
    figureName1 = figureName + "LMAX:" + str(LMAX) +  " l(l+1)_normalized\nrandom Map offset"
    plt.title(figureName1,fontsize = 8, y = 0.99)
    plt.savefig(saveName1)
    # print(alm[0])





    plt.figure(figsize = (figX, figY))
    fig2,ax2 = plt.subplots(num = fignum + 2)
    ax2.xaxis.set_major_locator(plt.MultipleLocator(xMax/CL_MAJOR_TICKS_NUM))
    ax2.xaxis.set_minor_locator(plt.MultipleLocator(xMax/(CL_MAJOR_TICKS_NUM * 2)))
    ax2.plot(ell, (2 * ell + 1) * clFinal)
    normedCl = (2*ell+1)*clFinal

    plt.xlabel("$\ell$")
    plt.ylabel("$(2\ell+1)C_{\ell}$")
    annot_max(ell, (2 *ell + 1) * clFinal, ax2)
    plt.grid()
    saveName2= saveName + "Diff_LMAX:" + str(LMAX) + "CLcoeff_2lp1.png"
    figureName2 = figureName + "LMAX:" + str(LMAX) +  " (2l+1)_normalized\nrandom Map offset"
    plt.title(figureName2,fontsize = 8, y = 0.99)

        #find local max
    indexX = localMax(ell,normedCl,ax2)
    maxXVal = ell[indexX]
    maxYVal = normedCl[indexX]
    plt.stem(maxXVal, maxYVal, linefmt = 'b--', markerfmt = 'ro', basefmt = ' ')
    plt.savefig(saveName2)


    plt.savefig(saveName2)

    # plt.figure(figsize= (figX, figY))
    # plt.plot(ell, coeff)
    # plt.xlabel("$\ell$")
    # plt.ylabel("$\widehat{c_{\ell}}$")
    # plt.grid()
    # saveName3= saveName + "LMAX:" + str(LMAX) + "CLcoeff_coeff.png"
    # figureName3 = figureName + "LMAX:" + str(LMAX) +  " legendre coeff"
    # plt.title(figureName3,fontsize = 8, y = 0.99)
    # plt.savefig(saveName3)


    plt.figure(figsize = (figX, figY))
    fig3,ax3 = plt.subplots(num = fignum + 3)
    ax3.xaxis.set_major_locator(plt.MultipleLocator(xMax/CL_MAJOR_TICKS_NUM))
    ax3.xaxis.set_minor_locator(plt.MultipleLocator(xMax/(CL_MAJOR_TICKS_NUM * 2)))
    ax3.plot(ell, clFinal)
    plt.xlabel("$\ell$")
    plt.ylabel("$C_{\ell}$")

    annot_max(ell, clFinal, ax3)
    plt.grid()
    saveName4= saveName + "Diff_LMAX:" + str(LMAX) + "CLcoeff_rawCL.png"
    figureName4 = figureName + "LMAX:" + str(LMAX) +  " Cl\nrandom Map offset"
    plt.title(figureName4,fontsize = 8, y = 0.99)
    plt.savefig(saveName4)    

def writeAngularPower(healpixMapFiles, \
    outName = "angularPower.txt", rescaled = True):
    """
    given lists of healpixMaps, write out the angular power spectrum
    within each healpixMap folder, with the name outName

    ----
    Parameters:
    ----

    healpixMapFiles: list of strs

    outName: str
        name of the file names
    
    rescaled: boolean
        boolean indicating whether the factor will be divided by 1/(4 pi zeta)

    return:
        number of angular power files written
    """

    count = 0
    for healpixFile in healpixMapFiles:
        
        returnedDict = plotCl(healpixFile,\
            '',\
            saveNamePre='',\
                pointCount=1,\
            SIGNAL_DATA = True, MY_DOUBLE = True,\
                spannedAngle = 1, fignum = 1, plotting = False)
        


        if returnedDict is not None:
            Nside = returnedDict['Nside']
            cl = returnedDict['cl']
            pointsCount = returnedDict['pointsCount']            
            ell = np.arange(len(cl))



            # ax.xaxis.set_major_locator(plt.MultipleLocator(xMax/CL_MAJOR_TICKS_NUM))
            # ax.xaxis.set_minor_locator(plt.MultipleLocator(xMax/(CL_MAJOR_TICKS_NUM * 2)))
            # ax.plot(ell, ell * (ell + 1) * cl)
            # plt.xlabel("$\ell$")
            # plt.ylabel("$\ell(\ell+1)C_{\ell}$")
            # plt.grid()
            # saveName1 = saveNamePre + "LMAX:" + str(LMAX) + "CLcoeff_llp1.png"
            # figureName1 = figureName + "LMAX:" + str(LMAX) +  " l(l+1)_normalized"        

            normedCl = (2 * ell + 1) * cl
            # find the subdirectory
            subdirs = healpixFile.split("/")
            subdirect = "/".join(subdirs[:-1])
            print("write angular power to {}".format(os.path.join(subdirect,outName)))
            with open(os.path.join(subdirect, outName),"w") as fileWrite:
                for l,coeff in zip(ell,normedCl):
                    fileWrite.write("{} {}\n".format(l,coeff))
            if (rescaled == True):
                normedCl_scaled = normedCl * 4 * np.pi / (pointsCount * pointsCount)
                fileName = outName.split('.')[0]
                outName_scaled = fileName + "_scaled.txt"
                with open(os.path.join(subdirect, outName_scaled), "w") as fileWrite:
                    for l,coeff in zip(ell,normedCl_scaled):
                        fileWrite.write("{} {}\n".format(l,coeff))

            count += 1

    return count
        
def plotMultiCls(healpixMapFiles, figureI,legendLists,colors,
    figureTitleName,
    saveName,fignum,plt, writeOutFile):
    """
    ---
    Parameters:
    ----

    writeOutFile: bool
        if True, then will write out the normed Cls to the file "angularPower.txt" in the directory where the 
        input healpixMapFile is
    """

    plt.figure(figsize = (figX, figY))

    fig,ax = plt.subplots(num = figureI)

    DoNotPlot = False

    for healpixFile,legend,color in zip(healpixMapFiles,legendLists,colors):
        Nside, cl = plotCl(healpixFile,\
            figureTitleName,\
            saveNamePre=saveName,\
                pointCount=1,\
            SIGNAL_DATA = True, MY_DOUBLE = True,\
                spannedAngle = 1, fignum = fignum, plotting = False)
        if (Nside > -1):
            ell = np.arange(len(cl))
            LMAX = Nside * 2

            xMax = LMAX

            plt.xlim(0,LMAX)

            # ax.xaxis.set_major_locator(plt.MultipleLocator(xMax/CL_MAJOR_TICKS_NUM))
            # ax.xaxis.set_minor_locator(plt.MultipleLocator(xMax/(CL_MAJOR_TICKS_NUM * 2)))
            # ax.plot(ell, ell * (ell + 1) * cl)
            # plt.xlabel("$\ell$")
            # plt.ylabel("$\ell(\ell+1)C_{\ell}$")
            # plt.grid()
            # saveName1 = saveNamePre + "LMAX:" + str(LMAX) + "CLcoeff_llp1.png"
            # figureName1 = figureName + "LMAX:" + str(LMAX) +  " l(l+1)_normalized"        


            normedCl = (2 * ell + 1) * cl
            ax.plot(ell, normedCl, label = legend,color = color)
            if (writeOutFile):
                # find the subdirectory
                subdirs = healpixFile.split("/")
                subdirect = "/".join(subdirs[:-1])
                print("write angular power to {}".format(os.path.join(subdirect,outNam)))
                with open(os.path.join(subdirect, outNam),"w") as fileWrite:
                    for l,coeff in zip(ell,normedCl):
                        fileWrite.write("{} {}\n".format(l,coeff))
        else:
            DoNotPlot = True
        
    if (DoNotPlot == False):
        # plt.xlabel("$\ell$")
        # plt.ylabel("$(2\ell+1)C_{\ell}$")

        # plt.grid()
        angularPowerAxisLabel(plt,ax)
        # x0, x1, y0, y1 = plt.axis()
        # plt.axis((x0 - plot_margin,
        #   x1 + plot_margin,
        #   y0 - plot_margin,
        #   y1 + plot_margin))
        # chartBox = ax.get_position()
        # ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.8, chartBox.height])
        # ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), shadow=True, ncol=1)        
        adjustLegendPadding(plt,ax)
    return DoNotPlot
    # saveName = saveNamePre + "LMAX:" + str(LMAX) + "CLcoeff_2lp1.png"
    # figureName = figureTitleNamePre + "LMAX:" + str(LMAX) +  " (2l+1)_normalized"

    
    # plt.title(figureTitleNamePre,fontsize = 8, y = 0.99)        
    # plt.savefig(saveNamePre)

if __name__ == "__main__":
    healpixMaps = ["/home/panwei/Desktop/bachelorThese/myCode/data/10/27/dataTest/LangenfeldData/HM_RXC_horn_MAXRSS_meas/100/points250avgTime_10nr_500thd_-100dTtB_0.523599tE_2.09440.01000_endT_1.88496ft_1/randomMap.txt", "/home/panwei/Desktop/bachelorThese/myCode/data/10/27/dataTest/LangenfeldData/HM_RXC_horn_MAXRSS/100/points122avgTime_10nr_244thd_-100dTtB_0.523599tE_2.09440.01000_endT_1.88496ft_1/SignalAverageHealpixMap.txt"]
    # healpixMaps = ["/home/panwei/Desktop/bachelorThese/myCode/data/10/27/dataTest/LangenfeldData/HM_RXA_horn_MAXRSS_meas/70/points3avgTime_10nr_6thd_-70dTtB_0.523599tE_1.57080.01000_endT_1.88496ft_1/SignalAverageHealpixMap.txt","/home/panwei/Desktop/bachelorThese/myCode/data/10/27/dataTest/LangenfeldData/HM_RXA_TRX_MAXRSS/70/points0avgTime_10nr_0thd_-70dTtB_0.523599tE_2.09440.01000_endT_1.88496ft_1/SignalAverageHealpixMap.txt"]
    plotMultiCls(healpixMaps,1,["experiment","rayTracing"],["r","b"],"test","/home/panwei/Desktop/bachelorThese/myCode/pythonCode",1,plt,writeOutFile=True)
    plt.savefig("/home/panwei/Desktop/bachelorThese/myCode/pythonCode/testMultiCl2.png")