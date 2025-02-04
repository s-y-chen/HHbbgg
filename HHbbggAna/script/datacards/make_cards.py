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
    
    backgrounds = ["resonant", "nonresonant"]
    signals = ["gghh"] # Todo add "VBFHH"
    
    backgrounds_pdf = [
                       "pdfs_dnn0507_bjet/wsinput.GaussianresonantggHHcat1.root",
                       "pdfs_dnn0507_bjet/wsinput.BernnonresonantggHHcat1.root",
                       #
                       "pdfs_dnn0507_bjet/wsinput.GaussianresonantggHHcat2.root",
                       "pdfs_dnn0507_bjet/wsinput.BernnonresonantggHHcat2.root",
                       #
                       "pdfs_dnn0507_bjet/wsinput.GaussianresonantggHHcat3.root",
                       "pdfs_dnn0507_bjet/wsinput.BernnonresonantggHHcat3.root",
                       #
                       "pdfs_dnn0507_bjet/wsinput.GaussianresonantggHHcat4.root",
                       "pdfs_dnn0507_bjet/wsinput.BernnonresonantggHHcat4.root",
                       #
                       "pdfs_dnn0507_bjet/wsinput.Gaussian_2resonantggHHcat5.root",
                       "pdfs_dnn0507_bjet/wsinput.BernnonresonantggHHcat5.root",
                       #
                       "pdfs_dnn0507_bjet/wsinput.GaussianresonantggHHcat6.root",
                       "pdfs_dnn0507_bjet/wsinput.BernnonresonantggHHcat6.root",
                      ]
    
    signals_pdf = ["pdfs_dnn0507_bjet/wsinput.GaussiangghhggHHcat1.root",
                   "pdfs_dnn0507_bjet/wsinput.GaussiangghhggHHcat2.root",
                   "pdfs_dnn0507_bjet/wsinput.GaussiangghhggHHcat3.root",
                   "pdfs_dnn0507_bjet/wsinput.GaussiangghhggHHcat4.root",
                    "pdfs_dnn0507_bjet/wsinput.GaussiangghhggHHcat5.root",
                   "pdfs_dnn0507_bjet/wsinput.GaussiangghhggHHcat6.root"
                  ]
    
    categories = ["ggHHcat1","ggHHcat2", "ggHHcat3", "ggHHcat4", "ggHHcat5", "ggHHcat6"]
                   
    ofname = "HHbbgg_datacard.txt"
    
    cat1_rates = [0.6206, 451.0, 1.0]
    cat2_rates = [0.5091, 109.8, 1.0]
    cat3_rates = [0.1829, 1.5405, 1.0]
    cat4_rates = [0.1284, 0.6112, 1.0]
    cat5_rates = [0.5175, 1.374, 1.0]
    cat6_rates = [0.2019, 0.6855, 1.0]
    
    rate_lst = [cat1_rates, cat2_rates, cat3_rates, cat4_rates, cat5_rates, cat6_rates]
    
    PrintDatacard(categories, signals, backgrounds, signals_pdf, backgrounds_pdf, rate_lst, ofname)    

if __name__=="__main__":
    main()
