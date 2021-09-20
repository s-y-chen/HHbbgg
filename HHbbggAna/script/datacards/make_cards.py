def PrintDatacard(categories, signals, backgrounds, filenames, ofname):
    dcof = open(ofname, "w")

    number_of_bins = len(categories)
    number_of_backgrounds = 0
    
    for cat in categories:
        for proc in cat.processes:
            if (proc in cat.signal_processes):
                signals += [proc]
            elif (proc in remove_proc):
                continue
            else:
                backgrounds += [proc]

    backgrounds = set(backgrounds)
    signals = set(signals)
    number_of_backgrounds = len(backgrounds)
    number_of_signals = len(signals)
    analysis_categories = list(set([c.full_name for c in categories]))
    dcof.write("imax {0}\n".format(number_of_bins))
    dcof.write("jmax {0}\n".format(number_of_backgrounds + number_of_signals - 1))
    dcof.write("kmax *\n")
    dcof.write("---------------\n")

    for cat in categories:
#old format
#        dcof.write("shapes * {0} {1} $PROCESS__$CHANNEL $PROCESS__$CHANNEL__$SYSTEMATIC\n".format(
        dcof.write("shapes * {0} {1} $PROCESS $PROCESS__$SYSTEMATIC\n".format(
            cat.full_name,
            os.path.basename(filenames[cat.full_name])
        ))

    dcof.write("---------------\n")

    dcof.write("bin\t" +  "\t".join(analysis_categories) + "\n")
    dcof.write("observation\t" + "\t".join("-1" for _ in analysis_categories) + "\n")
    dcof.write("---------------\n")

    bins        = []
    processes_0 = []
    processes_1 = []
    rates       = []

    for cat in categories:
        for i_sample, sample in enumerate(cat.processes):
            if (sample in remove_proc):
                continue
            else:
                bins.append(cat.full_name)
                processes_0.append(sample)
                if sample in cat.signal_processes:
                    i_sample = -i_sample
                processes_1.append(str(i_sample))
                rates.append("{0}".format(event_counts[sample]))

    #Write process lines (names and IDs)
    dcof.write("bin\t"+"\t".join(bins)+"\n")
    dcof.write("process\t"+"\t".join(processes_0)+"\n")
    dcof.write("process\t"+"\t".join(processes_1)+"\n")
    dcof.write("rate\t"+"\t".join(rates)+"\n")
    dcof.write("---------------\n")

def main():
    print("making datacards")
    
    backgrounds_proc = ["ggh","vbfh","tth","VH","nonresonant"]
    signals_proc = ["ggHH","VBFHH"]
    backgrounds_pdf = ["ggh","vbfh","tth","VH","nonresonant"]
    signals_pdf = ["ggHH","VBFHH"]
    categories = ["passDNN","failDNN"]
       
    PrintDatacard()    


if __name__=="__main__":
    main()
