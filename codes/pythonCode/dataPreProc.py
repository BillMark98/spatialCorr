import numpy as np


Convert2Deg = True



CUT_RANGE = True
if (Convert2Deg == True):
    cut_Val = 360
else:
    cut_Val = 3.14


def maxVal(x,y):
    return max(max(x),max(y))
def minVal(x,y):
    return min(min(x),min(y))

def cleanData(X,Y):
    # yValidIndex= np.isfinite(Y)
    # print(len(yValidIndex))

    X = X[~np.isnan(Y)]
    Y = Y[~np.isnan(Y)]
    X = X[np.isfinite(Y)]
    Y = Y[np.isfinite(Y)]
    return X,Y

# def cleanData(X,Y,Y2 = []):
#     # yValidIndex= np.isfinite(Y)
#     # print(len(yValidIndex))
#     X = X[np.isfinite(Y)]
#     if (len(Y2) >= 1):
#         Y2 = Y2[np.isfinite(Y)]
#     Y = Y[np.isfinite(Y)]
#     X = X[~np.isnan(Y)]
#     if (len(Y2) >= 1):
#         Y2 = Y2[~np.isnan(Y)]
#     Y = Y[~np.isnan(Y)]

#     if (len(Y2) >= 1):
#         X = X[np.isfinite(Y2)]
#         Y = Y[np.isfinite(Y2)]
#         Y2 = Y2[np.isfinite(Y2)]

#         X = X[~np.isnan(Y2)]
#         Y = Y[~np.isnan(Y2)]
#         Y2 = Y2[~np.isnan(Y2)]
#     return X,Y,Y2

def cutData(X,Y):
    global cut_Val
    Y = Y[(X < cut_Val)]
    X = X[(X < cut_Val)]
    return X,Y