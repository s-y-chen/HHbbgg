This is a directory for the retired or archived notebooks that are no longer in use. 

## Overview of the different notebooks in this directory

Sequential_DNN.ipynb: The original DNN script in case I (Stephanie) mess something up when trying things out on the other ones. 

Sequential_DNN_2018.ipynb: The DNN with 2018 samples that include the additional training variables of the DeepCSV b-tagging score but not the leading object pt to diobject pt ratios.

Sequential_DNN_combine_nodownsampling.ipynb: The DNN with all the backgrounds combined together into one large background without downsampling of the background.

Sequential_DNN_combine.ipynb: The DNN with all the backgrounds combined together into one large background. 

Sequential_DNN_data_separate.ipynb: The DNN using the data outside the 115-135 diphoton mass window as the background, separated by year. This was due to the inability to successfully run 2016 data through the selection code.

Sequential_DNN_data.ipynb: The DNN using the data outside the 115-135 diphoton mass window as the background.

Sequential_DNN_mass.ipynb: The DNN that includes the leading objet pt to diobject mass ratios in the training variables. This is also the notebook in which I have been trying to convert the DNN results and variables of the test set to a Root tree, which is now working.

Sequential_DNN_optimize_bayesian.ipynb: The DNN (with full recon==1 samples) optimized using Bayesian methods. As this optimization process seems to work (though does not significantly improve performance), this is likely the optimization method that will be kept moving forward. 

Sequential_DNN_optimize_talos.ipynb: The DNN (with full recon==1 samples) where I am trying to optimize the hyperparameters using an imported talos library (https://github.com/autonomio/talos). However, it takes quite a long time to run, and the accuracy seems worse than when I was using the DNN without optimization. I'm also trying to extract the parameters for the best model, though having some difficulty interpreting the output for that.

Sequence_DNN_optimize_RandomCV.ipynb: The DNN (with full recon==1 samples) where I am trying to optimize the hyperparameters using RandomCV from sklearn. This also takes quite some time to run. I have been running into two bugs though. The test scores of nan seems to occure when the model is not a classifier/does not predict anything, though that is not the case here. There may also be something with how the predict_classes method is no longer available in sklearn for Sequential models, but I haven't figured out how to address that yet. The other bug is the inability to clone due to n_neurons. I have tried turning the arrays into tuples as suggested on Stack Exchange, but that did not change anything. I also cannot tell why it is n_neurons specifically that is throwing the area. A post also said it might have been a bug with version 0.22.2, but the version I updated to was 0.24.2. 

ttH_DNN-SYC_cp_combine.ipynb: A DNN using the ttHScore model structure but applied on all MC backgrounds combined with the same training variables as used in all other models for this part of the analysis. 

ttH_DNN-SYC_cp.ipynb: A DNN using the ttHScore model structure but applied on an individual, selected MC background with the same training variables as used in all other models for this part of the analysis. 
