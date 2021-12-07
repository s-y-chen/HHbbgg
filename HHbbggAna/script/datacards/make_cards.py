# import ROOT as rt
# import uproot
# import numpy as np
# import pandas as pd

# def samp_to_df(samp_name):
#     file_name = '/storage/af/user/schen7/CMSSW_9_4_2/src/Higgs/HHbbgg/notebook/ML/DNN_Trees/combine_sequential_DNN/{0}_combine_seqDNN.root'.format(sampe_name)
#     samp_file = uproot.open(file_name)
#     samp_array = samp_file['tree'].arrays()
#     samp_df = pd.DataFrame(samp_array)
#     return samp_df

def PrintDatacard(categories, signals_proc, backgrounds_proc, signals_pdf, backgrounds_pdf, rate_lst, ofname):
    dcof = open(ofname, "w")

    number_of_categories = len(categories)
    number_of_backgrounds = 0
    
    number_of_backgrounds = len(backgrounds_proc)
    number_of_signals = len(signals_proc)
    
    dcof.write("imax {0}\n".format(number_of_categories))
    dcof.write("jmax {0}\n".format(number_of_backgrounds + number_of_signals - 1))
    dcof.write("kmax *\n")
    dcof.write("---------------\n")
        
    backgrounds = set(backgrounds_proc)
    signals = set(signals_proc)
    number_of_backgrounds = len(backgrounds_proc)
    number_of_signals = len(signals_proc)
    analysis_categories = list(set([c for c in categories]))
        
    print("analysis_categories: ",analysis_categories)

    #for cat in categories:
    #    dcof.write("shapes * {0} {1} $PROCESS $PROCESS__$SYSTEMATIC\n".format(
    #        cat,
    #        os.path.basename(filenames[cat])
    #    )

    dcof.write("---------------\n")

    dcof.write("bin\t" +  "\t".join(analysis_categories) + "\n")
    dcof.write("observation\t" + "\t".join("-1" for _ in analysis_categories) + "\n")
    dcof.write("---------------\n")

    bins        = []
    processes_0 = []
    processes_1 = []
    rates       = []
    icat = 0;
    for cat in categories:
                   
        for isig in range(len(signals)):
            bins.append(cat)
            processes_0.append(signals_proc[isig])
            i_sample = -isig
            processes_1.append(str(i_sample))         
            rates.append("{0}".format(rate_lst[icat][isig]))
            dcof.write("shapes\t"+signals_proc[isig]+"\t"+cat+"\t"+signals_pdf[icat*len(signals_proc)+isig]+ "\t"+"pdf:model"+signals_proc[isig]+cat+"\n")
               
        for ibkg in range(len(backgrounds)):
            bins.append(cat)
            processes_0.append(backgrounds_proc[ibkg])
            processes_1.append(str(ibkg+1))         
            rates.append("{0}".format(rate_lst[icat][ibkg+len(signals)]))  
            dcof.write("shapes\t"+backgrounds_proc[ibkg]+"\t"+cat+"\t"+backgrounds_pdf[icat*len(backgrounds_proc)+ibkg]+"\t"+"pdf:model"+backgrounds_proc[ibkg]+cat+"\n")
 
        #data
        dcof.write("shapes\tdata_obs \t"+cat+"\t"+backgrounds_pdf[icat*len(backgrounds_proc)+len(backgrounds)-1]+"\t" + "pdf:Data_13TeV_nonresonant_"+cat+"\n")

        icat+=1
        
    #Write process lines (names and IDs)
    dcof.write("bin\t"+"\t".join(bins)+"\n")
    dcof.write("process\t"+"\t".join(processes_0)+"\n")
    dcof.write("process\t"+"\t".join(processes_1)+"\n")
    dcof.write("rate\t"+"\t".join(rates)+"\n")
    dcof.write("---------------\n")
    for cat in categories:
        dcof.write("modelnonresonant"+cat+"_norm flatParam \n")

def main():
    print("making datacards")
    
    backgrounds = ["ggh","tth","vh","nonresonant"]
    signals = ["gghh"] # Todo add "VBFHH"
    
    backgrounds_pdf = ["pdfs_dnn_bjet/wsinput.GaussiangghggHHcat1.root",
                       "pdfs_dnn_bjet/wsinput.GaussiantthggHHcat1.root",
                       "pdfs_dnn_bjet/wsinput.GaussianvhggHHcat1.root",
                       "pdfs_dnn_bjet/wsinput.BernnonresonantggHHcat1.root",
                       "pdfs_dnn_bjet/wsinput.GaussiangghggHHcat2.root",
                       "pdfs_dnn_bjet/wsinput.GaussiantthggHHcat2.root",
                       "pdfs_dnn_bjet/wsinput.GaussianvhggHHcat2.root",
                       "pdfs_dnn_bjet/wsinput.BernnonresonantggHHcat2.root",
                        "pdfs_dnn_bjet/wsinput.GaussiangghggHHcat3.root",
                       "pdfs_dnn_bjet/wsinput.GaussiantthggHHcat3.root",
                       "pdfs_dnn_bjet/wsinput.GaussianvhggHHcat3.root",
                       "pdfs_dnn_bjet/wsinput.BernnonresonantggHHcat3.root",
                         "pdfs_dnn_bjet/wsinput.GaussiangghggHHcat4.root",
                       "pdfs_dnn_bjet/wsinput.GaussiantthggHHcat4.root",
                       "pdfs_dnn_bjet/wsinput.GaussianvhggHHcat4.root",
                       "pdfs_dnn_bjet/wsinput.BernnonresonantggHHcat4.root"
                      ]
    
    signals_pdf = ["pdfs_dnn_bjet/wsinput.GaussiangghhggHHcat1.root",
                   "pdfs_dnn_bjet/wsinput.GaussiangghhggHHcat2.root",
                   "pdfs_dnn_bjet/wsinput.GaussiangghhggHHcat3.root",
                   "pdfs_dnn_bjet/wsinput.GaussiangghhggHHcat4.root"
                  ]
    
    categories = ["ggHHcat1","ggHHcat2", "ggHHcat3",  "ggHHcat4" ]
                   
    ofname = "HHbbgg_datacard.txt"
    
    cat1_rates = [0.6206, 363.8, 1.972, 58.67, 1.0]
    cat2_rates = [0.5091, 80.09, 5.798, 17.57, 1.0]
    cat3_rates = [0.7004, 1.694, 0.06017, 1.036, 1.0]
    cat4_rates = [0.3302, 0.7709, 0.1766, 0.3322, 1.0]
    
    rate_lst = [cat1_rates, cat2_rates, cat3_rates, cat4_rates]
    
    PrintDatacard(categories, signals, backgrounds, signals_pdf, backgrounds_pdf, rate_lst, ofname)    

if __name__=="__main__":
    main()
