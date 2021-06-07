import sys
import matplotlib.pyplot as plt
import numpy as np





def main():
    g2c = {}

    annot_file = sys.argv[1]
    sample_name = sys.argv[3]
    fusion_info = {}
    with open(annot_file,"r") as hand:
        for line in hand:
            fields = line.rstrip().split("\t")
            genes = sorted(fields[7].split("::"))
            fin = float(fields[8])
            annot = fields[9]
            gene_name = "::".join(genes)
            ffigf = float(fields[18])

            fusion_info[gene_name] = (genes,fin,annot,ffigf,fields[9])
     
    separated_fusions = {}
    for x,y in fusion_info.items():
        t = y[4]
        if "segdup" in y[4]:
            t = "FAIL:segdup"
            continue
        elif "lowsup" in y[4]:
            t = "FAIL:lowsup"
        elif "lowfin" in y[4]:
            t = "FAIL:lowfin"
        elif "noncoding" in y[4]:
            t = "FAIL:noncoding"
            continue
        elif "overlaps" in y[4]:
            t = "FAIL:overlaps"
            continue
    #    else:
        
        if t not in separated_fusions:
            separated_fusions[t] = ([],[])
        separated_fusions[t][0].append(y[1])
        separated_fusions[t][1].append(y[3]) 


        
    plt.style.use('ggplot')
    fig = plt.figure(figsize=(12,12))
    #ax1,ax2 = fig.subplots(1,2,sharey='row')


    ax1 = fig.add_subplot()

    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_ylim(bottom=1,top=10**4)
    colors = range(0,len(list(separated_fusions.items())))
    for c,kv in zip(colors,separated_fusions.items()):
        k,v =kv
        X = v[0]
        Y = v[1]
        ax1.scatter(X,Y,label=k,alpha=0.75)
        #ax2.scatter(X,Y,label=k)
        
        
        ######
    plt.legend(fontsize=16)
    plt.title(sample_name,fontsize=20)
        
    fig.text(0.5, 0.075, 'FiN', ha='center',fontsize=16)
    #plt.xlabel("FiN",fontsize=18)
    ax1.set_ylabel("ff-igf",fontsize=16)

    plt.savefig(sys.argv[2],dpi=300)
    return 0



if __name__ == "__main__":
    exit(main())
