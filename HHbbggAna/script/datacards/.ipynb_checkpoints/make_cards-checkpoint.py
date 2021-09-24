def PrintDatacard(categories, signals_proc, backgrounds_proc, signals_pdf, backgrounds_pdf, ofname):
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

    tmp = 1.0
    icat = 0;
    for cat in categories:
                   
        for isig in range(len(signals)):
            bins.append(cat)
            processes_0.append(signals_proc[isig])
            i_sample = -isig
            processes_1.append(str(i_sample))         
            rates.append("{0}".format(tmp))
            dcof.write("shapes\t"+"\t"+signals_proc[isig]+"\t"+cat+"\t"+signals_pdf[icat*len(signals_proc)+isig]+ "pdf:model" + "\n")
               
        for ibkg in range(len(backgrounds)):
            bins.append(cat)
            processes_0.append(backgrounds_proc[ibkg])
            processes_1.append(str(ibkg))         
            rates.append("{0}".format(tmp))  
            dcof.write("shapes\t"+"\t"+backgrounds_proc[ibkg]+"\t"+cat+"\t"+backgrounds_pdf[icat*len(backgrounds_proc)+ibkg]+"pdf:model"+"\n")
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
    
    backgrounds_pdf = ["pdfs/wsinput.CBgghggHHcat1.root",
                       "pdfs/wsinput.CBvbfhggHHcat1.root",
                       "pdfs/wsinput.CBtthggHHcat1.root",
                       "pdfs/wsinput.CBvhggHHcat1.root",
                       "pdfs/wsinput.BernnonresonantggHHcat1.root",
                       "pdfs/wsinput.CBgghggHHcat2.root",
                       "pdfs/wsinput.CBvbfhggHHcat2.root",
                       "pdfs/wsinput.CBtthggHHcat2.root",
                       "pdfs/wsinput.CBvhggHHcat2.root",
                       "pdfs/wsinput.BernnonresonantggHHcat2.root"
                      ]
    
    signals_pdf = ["pdfs/wsinput.CBgghhggHHcat1.root",
                   "pdfs/wsinput.CBgghhggHHcat2.root"                  
                  ]
    
    categories = ["hbbhgg_gghhbin1_13TeV","hbbhgg_gghhbin2_13TeV"]
                   
    ofname = "HHbbgg_datacard.txt"
       
    PrintDatacard(categories, signals, backgrounds, signals_pdf, backgrounds_pdf, ofname)    

if __name__=="__main__":
    main()
