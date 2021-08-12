## Overview
Currently, the only signal being used is ggHH. 

## Overview of the different notebooks in this directory
Sequential_DNN.ipynb: The original DNN script in case I (Stephanie) mess something up when trying things out on the other ones. 

Sequential_DNN_2018.ipynb: The DNN with 2018 samples that include the additional training variables of the DeepCSV b-tagging score but not the leading object pt to diobject pt ratios.

Sequential_DNN_2018.ipynb: The DNN with 2018 samples that includes the leading objet pt to diobject mass ratios in the training variables. This is also the notebook in which I have been trying to convert the DNN results and variables of the test set to a Root tree. 

Sequential_DNN_combine.ipynb: The DNN with 2018 samples where I tried to combine some of the backgrounds together into one large background. 

Sequential_DNN_optimize_talos.ipynb: The DNN (with full recon==1 samples) where I am trying to optimize the hyperparameters using an imported talos library (https://github.com/autonomio/talos). However, it takes quite a long time to run, and the accuracy seems worse than when I was using the DNN without optimization. I'm also trying to extract the parameters for the best model, though having some difficulty interpreting the output for that.

Sequence_DNN_optimize_RandomCV.ipynb: The DNN (with full recon==1 samples) where I am trying to optimize the hyperparameters using RandomCV from sklearn. This also takes quite some time to run. I have been running into two bugs though. The test scores of nan seems to occure when the model is not a classifier/does not predict anything, though that is not the case here. There may also be something with how the predict_classes method is no longer available in sklearn for Sequential models, but I haven't figured out how to address that yet. The other bug is the inability to clone due to n_neurons. I have tried turning the arrays into tuples as suggested on Stack Exchange, but that did not change anything. I also cannot tell why it is n_neurons specifically that is throwing the area. A post also said it might have been a bug with version 0.22.2, but the version I updated to was 0.24.2. 

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
