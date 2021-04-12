#!/bin/bash
# matlab -nodisplay -nosplash -nodesktop -r "GpTheta(pi/4,pi/2);exit;"
matlab -nodisplay -nosplash -nodesktop -r "run('./mainFunc.m');exit;" # | tail -n +11