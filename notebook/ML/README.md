## Overview
Currently, the only signal being used is ggHH. 

## Overview of the different notebooks in this directory
Sequential_DNN.ipynb: The original DNN script in case I (Stephanie) mess something up when trying things out on the other ones. 

Sequential_DNN_2018.ipynb: The DNN with 2018 samples that include the additional training variables of the DeepCSV b-tagging score but not the leading object pt to diobject pt ratios.

Sequential_DNN_2018.ipynb: The DNN with 2018 samples that includes the leading objet pt to diobject mass ratios in the training variables. This is also the notebook in which I have been trying to convert the DNN results and variables of the test set to a Root tree. 

Sequential_DNN_combine.ipynb: The DNN with 2018 samples where I tried to combine some of the backgrounds together into one large background. 

(Not added yet) Sequential_DNN_optimize_trial.ipynb: The DNN (with full recon==1 samples) where I am trying to optimize the hyperparameters. 

Note: at the end, I will likely try to combine these into a single notebook.

## Overview of plots and other saved files
BackgroundVsSignal_Plts: This folder contains the plots comparing the training variable distributions of the background to the signal. 

DNN_Score_Nums: Text files containing the remaining signal and background yields at various DNN score cuts.

DNN_Score_Plts: The distribution of DNN score for the background and signal. This folder also contains the plots of the significance against the DNN score cutoff.

Efficiency: The text files showing the efficiency at background acceptance of 1% and 0.1%.

Epoch_Plts: The plots of the loss, validation loss, accuracy, and validation accuracy over the training epochs. 

Mass_Scuplt_Plts: The mass sculplting plots.

Models: Saved DNN model files.

ROC_Plts: The ROC curve plots. 
