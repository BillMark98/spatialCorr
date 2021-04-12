
# Introduction

This is the software used for the spation correlation analysis.
There two main parts of the repository.

* code  - which is for the programming usage
* data  - which contains some data both for the generation of analysis data as well as post calculation, so can be used to visualize, plotting etc.


***!!! Important !!!***

The `data` folder is quite large (~ 2.9GB). If you only want to pull down the code part. It is advisable to use `sparse-checkout`.
Example:

```bash
git init
git remote add -f origin https://github.com/BillMark98/spatialCorr.git

git config core.sparseCheckout true
echo "codes/" >> .git/info/sparse-checkout
git pull origin master
```

For some detailed information. It might be helpful to have a look on [my thesis](Thesis.pdf).

# Code part

The code part divides into three stages. Preprocessing, processing and postprocessing (what else can it be :-) )

## preprocessing

The code in the `matlabCode` is mainly used for the preprocessing. See the subfolder for more details.
The main usage is to generate the data that can be further processed. The incoming data is the raytracing experiment data.
This repository already contains certain data that you can start with, so you dont have to do everything from the scratch. The
needed data for further processing can be found in `data/matlabData`. For further info see that [subfolder](codes/matlabCode/README.md).

## Processing

The `cpp` code is mainly used for the processing. It takes data input for the raytracing data, and then generate corresponding 
estimator correlation. For further details, see the corresponding [subfolder](codes/cppCode/README.md).

The processing step can also plays with itself, in particular, you can generate random `Poisson` process, `Matern` processes etc to get a
grasp of the behaviror of the estimators.

## Post processing

The `python` code is mainly used for the post processing. In this step. The `angular Power Spectrum` can also be generated. See the [subfolder](codes/pythonCode/README.md) for more details.
The idea is that the `python` code takes the generated data and plot the figures in the desired form.

# Data part

This folder contains some example data that you can play with, ranging from the data I used to generate analysis for my thesis, to raytracing data from the 
`superC`, `Langenfeld` , `Frankfurt` and `Seoul`. For more details, please have a look at the [subfolder](data/README.md).

