## Overview
Currently, the only signal being used is ggHH. 

## Overview of the different notebooks in this directory

Multiclass_DNN.ipynb: The multiclass DNN that uses integer values for the classes and includes down-sizing for the large diphoton+jets background. 

Multiclass_DNN_int.ipynb: The multiclass DNN that uses integer values for the classes and does not include down-sizing. 

Multiclass_DNN_one_hot.ipynb: The multiclass DNN that uses one-hot encording for the classes and does not include down-sizing.

Sequential_DNN_ggHH.ipynb: The DNN used for distinguishing the ggHH signal from the backgrounds (with VBFHH included as a background). The hyperparameters of the model has been optimized by Bayesian optimization. The model is trained and tested on the combined MC background and signal samples, and then applied to the data. 

## Overview of plots and other saved files
BackgroundVsSignal_Plts: This folder contains the plots comparing the training variable distributions of the background to the signal. 

DNN_Score_Nums: Text files containing the remaining signal and background yields at various DNN score cuts.

DNN_Score_Plts: The distribution of DNN score for the background and signal. This folder also contains the plots of the significance against the DNN score cutoff.

Efficiency: The text files showing the efficiency at background acceptance of 1% and 0.1%.

Epoch_Plts: The plots of the loss, validation loss, accuracy, and validation accuracy over the training epochs. 

Mass_Scuplt_Plts: The mass sculplting plots.

Models: Saved DNN model files.

ROC_Plts: The ROC curve plots. 
