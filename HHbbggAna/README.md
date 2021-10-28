cmsrel CMSSW_9_4_2

cd CMSSW_9_4_2/src

cmsenv

make clean; make

To run the code

./analyzeHHbbgg testfile/lists_2018/GluGluToHHTo2B2G_node_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8.txt  out.root GluGluToHHTo2B2G_node_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8 F T 2018

OR 

./analyzeHHbbgg testfile/lists_2016/data_ex.txt data T F 2016

nanoAOD:

https://cms-nanoaod-integration.web.cern.ch/integration/cms-swCMSSW_10_6_19/mc102X_doc.html#Photon
https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD#Photons

cms datasets:

https://twiki.cern.ch/twiki/bin/view/CMS/PdmVDataReprocessing
https://twiki.cern.ch/twiki/bin/view/CMS/PdmVDataReprocessingNanoAODv7

data

/EGamma/Run2018A-02Apr2020-v1/NANOAOD
/EGamma/Run2018B-02Apr2020-v1/NANOAOD
/EGamma/Run2018C-02Apr2020-v1/NANOAOD
/EGamma/Run2018D-02Apr2020-v1/NANOAOD

https://cmsweb.cern.ch/das/request?view=plain&limit=50&instance=prod%2Fglobal&input=dataset%3D%2F*DoubleEG*%2F*Run201*02Apr2020*%2FNANOAOD

/DoubleEG/Run2016B-02Apr2020_ver1-v1/NANOAOD
/DoubleEG/Run2016B-02Apr2020_ver2-v1/NANOAOD
/DoubleEG/Run2016C-02Apr2020-v1/NANOAOD
/DoubleEG/Run2016D-02Apr2020-v1/NANOAOD
/DoubleEG/Run2016E-02Apr2020-v1/NANOAOD
/DoubleEG/Run2016F-02Apr2020-v1/NANOAOD
/DoubleEG/Run2016G-02Apr2020-v1/NANOAOD
/DoubleEG/Run2016H-02Apr2020-v1/NANOAOD
/DoubleEG/Run2017B-02Apr2020-v1/NANOAOD
/DoubleEG/Run2017C-02Apr2020-v1/NANOAOD
/DoubleEG/Run2017D-02Apr2020-v1/NANOAOD
/DoubleEG/Run2017E-02Apr2020-v1/NANOAOD
/DoubleEG/Run2017F-02Apr2020-v1/NANOAOD

CMS DAS

 The DAS command line tool dasgoclient is available in any CMSSW releases.

    How can I use DAS CLI?

See help section of the DAS CLI tool:

dasgoclient --help # Go-based DAS CLI tool (recommended)
dasgoclient --examples # provides list of queries it supports

Here are few examples of DAS cli usage:

dasgoclient --query="dataset=/EG/Run2010A*/AOD"
dasgoclient --query="dataset=/EG/Run2010A*/AOD" --verbose=1
dasgoclient --query="dataset=/EG/Run2010A*/AOD | grep dataset.name"
dasgoclient --query="dataset=/EG/Run2010A*/AOD | grep dataset.name" --format=json

dasgoclient --query="dataset=/TTGG_0Jets_TuneCP5_13TeV_amcatnlo_madspin_pythia8/RunIIAutumn18NanoAODv7*/NANOAODSIM"

dasgoclient --query="file dataset=/TTGG_0Jets_TuneCP5_13TeV_amcatnlo_madspin_pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM"
