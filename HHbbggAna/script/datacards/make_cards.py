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
            dcof.write("shapes\t"+"\t"+signals_proc[isig]+"\t"+cat+"\t"+signals_pdf[icat*len(signals_proc)+isig]+ "\t"+ "pdf:model" + "\n")
               
        for ibkg in range(len(backgrounds)):
            bins.append(cat)
            processes_0.append(backgrounds_proc[ibkg])
            processes_1.append(str(ibkg+1))         
            rates.append("{0}".format(rate_lst[icat][ibkg+len(signals)]))  
            dcof.write("shapes\t"+"\t"+backgrounds_proc[ibkg]+"\t"+cat+"\t"+backgrounds_pdf[icat*len(backgrounds_proc)+ibkg]+"\t" + "pdf:model"+"\n")
        icat+=1
        
    #Write process lines (names and IDs)
    dcof.write("bin\t"+"\t".join(bins)+"\n")
    dcof.write("process\t"+"\t".join(processes_0)+"\n")
    dcof.write("process\t"+"\t".join(processes_1)+"\n")
    dcof.write("rate\t"+"\t".join(rates)+"\n")
    dcof.write("---------------\n")

def main():
    print("making datacards")
    
    backgrounds = ["ggh","vbfh","tth","VH","nonresonant"]
    signals = ["ggHH"] # Todo add "VBFHH"
    
#     backgrounds_pdf = ["pdfs/wsinput.CBgghggHHcat1.root",
#                        "pdfs/wsinput.CBvbfhggHHcat1.root",
#                        "pdfs/wsinput.CBtthggHHcat1.root",
#                        "pdfs/wsinput.CBvhggHHcat1.root",
#                        "pdfs/wsinput.BernnonresonantggHHcat1.root",
#                        "pdfs/wsinput.CBgghggHHcat2.root",
#                        "pdfs/wsinput.CBvbfhggHHcat2.root",
#                        "pdfs/wsinput.CBtthggHHcat2.root",
#                        "pdfs/wsinput.CBvhggHHcat2.root",
#                        "pdfs/wsinput.BernnonresonantggHHcat2.root"
#                       ]
    
#     signals_pdf = ["pdfs/wsinput.CBgghhggHHcat1.root",
#                    "pdfs/wsinput.CBgghhggHHcat2.root"                  
#                   ]
    backgrounds_pdf = ["pdfs/wsinput.GaussiangghggHHcat1.root",
                       "pdfs/wsinput.GaussianvbfhggHHcat1.root",
                       "pdfs/wsinput.GaussiantthggHHcat1.root",
                       "pdfs/wsinput.GaussianvhggHHcat1.root",
                       "pdfs/wsinput.BernnonresonantggHHcat1.root",
                       "pdfs/wsinput.GaussiangghggHHcat2.root",
                       "pdfs/wsinput.GaussianvbfhggHHcat2.root",
                       "pdfs/wsinput.GaussiantthggHHcat2.root",
                       "pdfs/wsinput.GaussianvhggHHcat2.root",
                       "pdfs/wsinput.BernnonresonantggHHcat2.root"
                      ]
    
    signals_pdf = ["pdfs/wsinput.GaussiangghhggHHcat1.root",
                   "pdfs/wsinput.GaussiangghhggHHcat2.root"                  
                  ]
    
    categories = ["hbbhgg_gghhbin1_13TeV","hbbhgg_gghhbin2_13TeV"]
                   
    ofname = "HHbbgg_datacard.txt"
    
    # get rates
#     samples = ['sig', 'GluGluHtoGG', 'VBFHToGG', 'ttHToGG', 'VHToGG', 'DiPhotonJetsBox']
#     dfs = []
#     for samp in samples:
#         dfs.append(samp_to_df(samp))
#     cat1_dfs = []
#     cat2_dfs = []
#     for df in dfs:
#         c1_df = df[df[b'DNN_score'] >= 0.5]
#         c2_df = df[df[b'DNN_score'] < 0.5]
#         cat1_dfs.append(c1_df)
#         cat2_dfs.append(c2_df)
#     cat1_rates = []
#     cat2_rates = []
#     for i in range(len(samples)):
#         c1_yield = np.sum(cat1_dfs[i].loc[:,b'genweight_scale'].values)
#         c2_yield = np.sum(cat2_dfs[i].loc[:, b'genweight_scale'].values)
#         cat1_rates.append(c1_yield)
#         cat2_rates.append(c2_yield)
    cat1_rates = [1.909, 3.181, 2.227, 10.35, 0.03114, 5486]
    cat2_rates = [0.2509, 30.18, 5.860, 66.83, 0.01578, 112100]
    rate_lst = [cat1_rates, cat2_rates]
    
    PrintDatacard(categories, signals, backgrounds, signals_pdf, backgrounds_pdf, rate_lst, ofname)    

if __name__=="__main__":
    main()
