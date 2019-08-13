# Name explanation
Most of the functions and scripts starts with an 'A_' because I want to distinguish mine functions and yours. There shouldn't be any conflict.

# Folders expmalnation

## simulation
This folder is the simulation scripts for _old model_ and _perfused model_

2 '__parameter sensitivity__' folders are to test the influence of parameters on the perfused model or shuyu's changed model. To run for all three parameters, just need to run the '__\_main.m__'script.

'__shuyu changed model__' is the source code that associated with his model.

## real data fitting
This folder is to fit on the real data.
It relys deeply on the folder 'data'.

### old model

_A\_old\_fitting_ is using old model to fit for real data.

### perfused model

_A\_perfused\_fitting_ is using perfused model, able to compare contorl, tumor, noise.

_A\_fitting\_real\_perfused_ is using perfused model, able to compare intra, extra vascular.

## data
This is the folder where all the data is in.
_Don't_ change the reletive relationship between '_data_' and '_Andy testing_'. 'real data fitting' is accessing this folder.

Also, about the data. To make the script runs faster, I first transform the xlsx into '.mat' file. So the simulation scripts only need to load '.mat' file.

The '.mat' data is in 'Andy_fitted'.

### Andy_data
This folder contains original xlsx. Remember to put different tumors in 3 different folders: '_786O_','_A498_','_uok262_'.

Now that few of them contains the VIF information for the perfused model, I picked them up and put in the 

## helper_funcs
This folder contains helper functions. They are used to get the name of the data, or to load the data.

## presentation
This folder contains all the figures I generated throughout the summer. Take a look if you need to refresh about the progress in the past.
