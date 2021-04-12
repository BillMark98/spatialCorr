import scipy.io
import numpy as np
import os
import fileIO
namePos = [str(i) for i in range(1,5)]
print(namePos)
# dirPath = "../data/matlabData/LangenfeldData/"
dirPath = "../data/matlabData/superCData/measurement"

# filePathPref = "../../data/Langenfeld TX1_RXA/SimulationData_15 dBi_horn_thresh_91/HM_RXA_horn_MAXRSS"#"../../data/Langenfeld TX1_RXA/SimulationData_15 dBi_horn_thresh_104/HM_RXA_horn_TXPOS"
# filePathPref = "../../data/LangenfeldData/"
filePathPref = "~/Desktop/bachelorThese/data/rss_superC/rss_measurement_evaluation_toolbox_Panwei/SuperC"


def mat2csv(filePaths = [""]):
	traversePath = []
	if (len(filePaths) == 1 and len(filePaths[0]) == 0):
		traversePath = os.listdir(filePathPref)
	else:
		traversePath = fileIO.getFilesFromDirs(filePaths)
	for file in traversePath:
		if file.endswith(".mat"):
			print("file name: " + file)
		data = scipy.io.loadmat(os.path.join(filePathPref,file))
		fileName = file.split(".")[0]
		saveName = os.path.join(dirPath,fileName)
		if not os.path.exists(saveName):
			os.mkdir(saveName)
			
		print("saved at:")
		print(saveName)
		for i in data:
			if '__' not in i and 'readme' not in i:
				print(i)
				print("data:")
				print(data[i])
				print("dim: ")
				print(data[i].ndim)
				if (data[i].ndim < 3):
					iName = i.split("_")[0]
					if (not iName.endswith("Final")):
						iName += "Final"
					np.savetxt((os.path.join(saveName, (iName + ".csv"))), data[i], delimiter=' ')
				else:
					print("dim: 3, not processed")
					print("name: {0:}".format(i))

	# filepath = filePathPref + ".mat"
	# data = scipy.io.loadmat(filepath)
	# savePath = dirPath
	# for i in data:
	# 	if '__' not in i and 'readme' not in i:
	# 		np.savetxt((savePath+i+".csv"),data[i],delimiter=' ')
	# for index in range(4):
	# 	filepath = filePathPref + namePos[index] + ".mat"
	# 	data = scipy.io.loadmat(filepath)
	# 	savePath = dirPath + namePos[index]
	# 	for i in data:
	# 		if '__' not in i and 'readme' not in i:
	# 			np.savetxt((savePath+i+".csv"),data[i],delimiter=' ')



if __name__ == "__main__":
	# mat2csv()
	# os.chdir(dirPath)
	# for root, dirname, filename in os.walk('.'):
	# 	print("root:")
	# 	print(root)
	# 	print("dir")
	# 	print(dirname)
	# 	print("filename")
	# 	print(filename)

	print("test")
	# mat2csv(["/home/panwei/Desktop/bachelorThese/data/rss_superC/rss_measurement_evaluation_toolbox_Panwei/SuperC/TX1_RXV"])