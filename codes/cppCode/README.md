This cpp code is used to generate 2pCF analysis. The code
is based on the software [Healpix](https://healpix.sourceforge.io/).
To use this software, `Healpix` has to be preinstalled, the installation guide can be found at [here](https://healpix.sourceforge.io/downloads.php#SECTION_360_12000000000000000)

# Installation

***!!! Need to preinstall healpix successfully to proceed !!!***

Go to the `CMakeLists.txt` to change the corresponding path:

```cmake
set(HEALPIX_SOURCE_DIR "/home/panwei/Documents/Healpix_3.50/src/cxx/basic_gcc/")
set(CFITSIO_SOURCE_DIR "/home/panwei/Documents/cfitsio-3.47/")
```

Set the corresponding path to the `Healpix` source and the `cfitsio`.

```bash
mkdir build
cd build
cmake ..
make exampleIOtest
```

The corresponding configuration file can be found in `testConfig/` folder. You could specify the content for your own usage.
The exact meaning of each entry is the following:

```txt
"/home/panwei/Desktop/bachelorThese/myCode/data/matlabData/seoul/phasedAntennaArrayAggregate/"   --->  the parent path to the data folder
"RX_450_600_raytracing" "RX_400_380_raytracing" ---> the folder names where the input data is stored
"/home/panwei/Desktop/bachelorThese/myCode/data/3/31/dataTest/seoul/phasedAntennaArrayAggregateNewMode/"  --> the folder path to output the calculation result
"RX_450_600_raytracing" "RX_400_380_raytracing" --->  the folder path where the calculation result is stored
-95 -92 ---> the threshold values
4 ---> the mode for the correlation (see `Healpix_Correlation::SetPixValue` for more details)
```

## Results

### Log files

The explicit log files can be found under the folder `./build/logs/`
The `calcuInfo.txt` gives out the information for the calculation of the correlation estimators.
The `DR.txt` gives the calculed `DD`, `DR` and `RR` vectors, the first column indicates the bin width (in radians)

### result folders

The result folders consist of the following files
- `calculationInfo.txt` - the parameters for calculation correlation
    - the variables has the following meanings:
    1. `pointsCount` - the total point of the `Healpix_Map`
    2. `randomCount` - the total point of the randomly generated map
    3. `deltaTheta` - the bin width
    4. `AVERAGE_TIME` - the times doing an average for the calculation
    5. `spannedAngle` - the angle spanned of the considered calculation region
    6. `endTheta` - the end Theta of the considered region (set by the user)
    7. `thetaBegin` - the begin Theta of the considered region
    8. `thetaEnd` - the actual (note the difference to `endTheta`) maximal theta angle achieved by the data set
    9. `threshold` - the threshold value set
    10. `FACTOR_SIGVAL` - the extra factor multiplied when converting signal level to data point (default is 1, can be set in `include/myDefineConst.h`)
    11. `outFilePre` - the name of the folder where the input data is coming from
- `SignalAverageHealplixMap.txt` - the corresponding `Healpix` map points
- `SignalAverageHealpixMode*.txt` where `*` is a number from 1 - 4. This is the resultant estimator, for the estimators 1 - 4. See my thesis for detail about the naming. Typicall the estimator 4 which is the Landy-Szalay estimator is used.
- `spherePoint.txt` - file indicating the 3D point location (cartesian x,y,z coordinates)
- `randomMap.txt` - the file for the randomly generated map for the calculation

