import xlrd
import itertools
import math

aa={}
aa["F"]=["TTT","TTC"]
aa["L"]=["TTA","TTG","CTT","CTC","CTA","CTG"]
aa["I"]=["ATT","ATC","ATA"]
aa["M"]=["ATG"]
aa["V"]=["GTT","GTC","GTA","GTG"]
aa["S"]=["TCT","TCC","TCA","TCG","AGT","AGC"]
aa["P"]=["CCT","CCC","CCA","CCG"]
aa["T"]=["ACT","ACC","ACA","ACG"]
aa["A"]=["GCT","GCC","GCA","GCG"]
aa["Y"]=["TAT","TAC"]
aa["H"]=["CAT","CAC"]
aa["Q"]=["CAA","CAG"]
aa["N"]=["AAT","AAC"]
aa["K"]=["AAA","AAG"]
aa["D"]=["GAT","GAC"]
aa["E"]=["GAA","GAG"]
aa["C"]=["TGT","TGC"]
aa["W"]=["TGG"]
aa["R"]=["CGT","CGC","CGA","CGG","AGA","AGG"]
aa["G"]=["GGT","GGC","GGA","GGG"]
aa["*"]=["TAA","TAG","TGA"]
codons={}
codons["TTY"]="F"
codons["TTR"]="L"
codons["CTN"]="L"
codons["ATH"]="I"
codons["ATG"]="M"
codons["GTN"]="V"
codons["TCN"]="S"
codons["AGY"]="S"
codons["CCN"]="P"
codons["ACN"]="T"
codons["GCN"]="A"
codons["TAY"]="Y"
codons["CAY"]="H"
codons["CAR"]="Q"
codons["AAY"]="N"
codons["AAR"]="K"
codons["GAY"]="D"
codons["GAR"]="E"
codons["TGY"]="C"
codons["TGG"]="W"
codons["CGN"]="R"
codons["AGR"]="R"
codons["GGN"]="G"
codons["TAR"]="*"
codons["TGA"]="*"
iupac={}
iupac["A"]={"A"}
iupac["G"]={"G"}
iupac["C"]={"C"}
iupac["T"]={"T"}
iupac["R"]={"A","G"}
iupac["Y"]={"C","T"}
iupac["S"]={"G","C"}
iupac["W"]={"A","T"}
iupac["K"]={"G","T"}
iupac["M"]={"A","C"}
iupac["B"]={"C","G","T"}
iupac["D"]={"A","G","T"}
iupac["H"]={"A","C","T"}
iupac["V"]={"A","C","G"}
iupac["N"]={"N"}

wb = xlrd.open_workbook("MutaGenie.xlsx")
sheet = wb.sheet_by_index(0)

def translate(plasmid,base):
    seqa=""
    iupac={}
    iupac["A"]=["A","R","H","N"]
    iupac["G"]=["G","R","N"]
    iupac["T"]=["T","Y","H","N"]
    iupac["C"]=["C","Y","H","N"]
    while base+2<len(plasmid):
        codon=plasmid[base:base+3]
        n=0
        while codon[:2]+iupac[codon[2]][n] not in codons.keys():
            n+=1
        if codons[codon[:2]+iupac[codon[2]][n]] =="*":
            seqa+=codons[codon[:2]+iupac[codon[2]][n]]
            break
        else:
            seqa+=codons[codon[:2]+iupac[codon[2]][n]]
        base+=3
    return(seqa)

def find_all(template,find):
    n=0
    pos=[]
    while True:
        n=template.find(find,n)
        if n==-1: break
        pos+=n,
        n+=len(find)
    return(pos)

def enzyme(site):
    b=find_all(site,"(")
    front=[]
    back=[]
    for i in range(len(b)):
        if b[i]==0:
            front+=site[1:site.find(")")],
        else:
            back+=site[b[i]+1:-1],
    for i in range(len(front)):
        front[i]=math.floor(sum([int(i) for i in front[i].split("/")])/len(front[i].split("/")))
    for i in range(len(back)):
        back[i]=math.floor(sum([int(i) for i in back[i].split("/")])/len(back[i].split("/")))
    site="".join([i for i in site if i.isalpha()])
    while site[0]=="N":
        site=site[1:]
        front=list([i+1 for i in front])
    while site[-1]=="N":
        site=site[:-1]
        back=list([i+1 for i in back])
    front=list([-i for i in front])
    back=list([i+len(site) for i in back])
    if len(front+back)==0:
        front=[math.floor(len(site)/2)]
    site=list([iupac[i] for i in site])
    if {"N"} in site:
        site.insert(site.index({"N"}),[str(site.count({"N"}))])
    site=list([i for i in site if not i=={"N"}])
    site=list(itertools.product(*site))
    site=list(["".join(i) for i in site])
    return([site,front+back])

def cut(enz,template):
    cuts=[]
    if enz[0][0].isalpha():
        for i in enz[0]:
            cuts+=find_all(template,i)
    else:
        num="".join([i for i in enz[0][0] if i.isnumeric()])
        for i in enz[0]:
            f,b=i.split(num)
            temp=find_all(template,f)
            cuts+=[c for c in temp if template[c+len(f)+int(num):c+len(f+b)+int(num)]==b]
    cuts.sort()
    cuts=[c+f for c in enz[1] for f in cuts]
    return(cuts)

def different(enz,wt,mut): #[[str],[int]] , str , [str]
    cuts=[]     #[[int]] list of cut patterns
    slices={}   #dict[str([int])]=[str]
    wtcut=cut(enz,wt) #list of cut sites
    for m in mut: #str in [str]
        temp=cut(enz,m) #list of cut sites
        if not temp==wtcut:
            if temp in cuts:
                slices[str(temp)]+=m,
            else:
                cuts+=temp,
                slices[str(temp)]=[m] #dict[str([int])]=[str]
    if len(cuts)==0:
        cuts=wtcut
    for key in slices.keys(): #each different cut pattern
        temp=list(range(len(wt))) #positions in string
        tempa=[] #[[int]] list of dif positions for strings that cut differently
        for string in slices[key]: #strings that yield different cut
            tempa+=[pos for pos in range(len(wt)) if not wt[pos]==string[pos]],
        mutmin=min([len(i) for i in tempa]) #min required mutations
        slices[key]=[slices[key][i] for i in range(len(tempa)) if len(tempa[i])==mutmin]
        tempa=[i for i in tempa if len(i)==mutmin]
        for string in slices[key]:
            temp=[pos for pos in temp if not wt[pos]==string[pos]]
        bases={}
        alt=[]
        for pos in temp:
            bases[pos]=set([string[pos] for string in slices[key]])
            bases[pos]=[letter for letter,letters in iupac.items() if letters==bases[pos]][0]
        if not len(temp)==mutmin:
            alt=[[int("".join([letter for letter in split if letter.isnumeric()])) for split in i.split(",")] for i in set([str([a for a in i if a not in temp]) for i in tempa])]
            altbases=[[set([string[pos] for string in slices[key] if pos in tempa[slices[key].index(string)]]) for pos in positions] for positions in alt]
            altbases=[[[letter for letter,letters in iupac.items() if letters==pos][0] for pos in positions] for positions in altbases]
        slices[key]=[wt[pos]+str(plasmid.find(wt)+pos+1)+bases[pos] for pos in bases.keys()]
        for altpos in range(len(alt)):
            slices[key]+=[wt[alt[altpos][i]]+str(plasmid.find(wt)+alt[altpos][i]+1)+altbases[altpos][i] for i in range(len(alt[altpos]))],
    return[wtcut,cuts,slices]

def run_gel(plasmid,cutsites):
    cutsites=list(cutsites)
    cutsites=cutsites+[0,len(plasmid)]
    cutsites.sort()
    bands=[]
    for i in range(len(cutsites)-1):
        bands+=cutsites[i+1]-cutsites[i],
    try:
        bands[0]=bands[0]+bands.pop()
        bands.sort()
    except:
        bands=["uncut"]
    return(bands)

def compare(wt,mut): #lists of integers
    old=list(wt)
    new=list(mut)
    for i in range(len(wt)):
        temp=old.pop(0)
        try:
            new.remove(temp)
        except:
            old+=temp,
    if not "uncut" in mut and not "uncut" in wt and max(old+new)>=minband:
        score=max([min([abs(math.log10(o)-math.log10(m)) for m in mut if type(m)==int]) for o in old if o>=minband and type(o)==int]+[min([abs(math.log10(n)-math.log10(w)) for w in wt if type(w)==int]) for n in new if n>=minband and type(n)==int])
    try:
        score=max([min([abs(math.log10(o)-math.log10(m)) for m in mut if type(m)==int]) for o in old if o>=minband and type(o)==int]+[min([abs(math.log10(n)-math.log10(w)) for w in wt if type(w)==int]) for n in new if n>=minband and type(n)==int])
    except:
        score=max([o for o in old if type(o)==int]+[n for n in new if type(n)==int])
        if score<minband:
            score=0
        else:
            score=10
    return(score,old,new)

def summarise(enzymes):
    results={}
    for enz in enzymes.keys():
        wtcut,mutcuts,mutations=different(enzymes[enz],wtslice,mutslice)
        if len(mutations)>0:
            results[enz]={}
            wtsites=set(cut(enzymes[enz],plasmid)+[i-100 for i in cut(enzymes[enz],plasmid[-100:]+plasmid[:100])])
            temp=set([i for i in wtsites if i<1])
            wtsites=wtsites.difference(temp)
            temp=set([len(plasmid)+i for i in temp])
            wtsites.update(temp)
            wtbands=run_gel(plasmid,wtsites)
            for mutcut in mutcuts:
                mutsites=set(wtsites)
                for i in wtcut:
                    mutsites.remove(i+plasmid.find(wtslice))
                for i in mutcut:
                    mutsites.add(i+plasmid.find(wtslice))
                mutbands=run_gel(plasmid,mutsites)
                results[enz][str(mutcut)]=list(compare(wtbands,mutbands))+[mutations[str(mutcut)],wtbands,mutbands]
    return(results)

def make_table(results):
    lines=[]
    temp=[]
    for enz in results.keys():
        for pattern in results[enz].keys():
            lines+=[enz]+results[enz][pattern],
    ranking=[line[1] for line in lines if line[1]>minsep]
    ranking.sort()
    ranking.reverse()
    for score in ranking:
        temp+=[[len(temp)+1]+line for line in lines if line[1]==score]
        lines=[line for line in lines if not line[1]==score]
    lines=list(temp)
    b=max([len(str(i[3])+str(i[4])) for i in lines]+[12])
    e=max([len(i[1]) for i in lines]+[6])
    output=["Rank|"+"Enzyme".center(e," ")+"|Log10 separation|"+"Band changes".center(b," ")+"|Mutations|Mutation range"]
    for line in lines:
        if line[2]==10:
            line[2]="N/A"
        else:
            line[2]=round(line[2],2)
        band=str(line[3])[1:-1]+" -> "+str(line[4])[1:-1]
        band=band.replace("'","")
        mut=str(len([i for i in line[5] if type(i)==str])+min([len(o) for o in line[5] if type(o)==list], default=0))
        mutrange=str(1+min([max(e)-min(e) for e in [[int(i[1:-1]) for i in line[5] if type(i)==str]+[[int(o[1:-1]) for o in a][0] for a in line[5] if type(a)==list]]]))
        output+=[str(line[0]).center(4," ")+"|"+line[1].center(e," ")+"|"+str(line[2]).center(16," ")+"|"+band.center(b," ")+"|"+mut.center(9," ")+"|"+mutrange.center(14," ")]
    for line in output:
        print("|"+line+"|")
    print("")
    return(lines)

yes={}
no={}
for row in range(len(sheet.col_values(0,1))):
    if sheet.cell_value(1+row, 2).upper()=="Y":
        yes[sheet.cell_value(1+row,0)]=enzyme(sheet.cell_value(1+row,1))
    elif sheet.cell_value(1+row, 2).upper()=="N":
        no[sheet.cell_value(1+row,0)]=enzyme(sheet.cell_value(1+row,1))
plasmid=sheet.cell_value(1,6).upper()
cds=sheet.cell_value(4,6).upper()
mut=sheet.cell_value(7,6).upper()
mutrange=int(sheet.cell_value(1,4))
minband=int(sheet.cell_value(4,4))
minsep=float(sheet.cell_value(7,4))

base=atg=plasmid.find(cds)
seqa=translate(plasmid,base)
mutslice=[]

if seqa[int(mut[1:-1])-1]==mut[0]:
    muta=int(mut[1:-1])-1
    mutb=atg+3*(muta)
    seqa=seqa[:muta]+mut[-1]+seqa[muta+1:]
    wtslice=plasmid[mutb-30-mutrange*3:mutb+30+mutrange*3]
for pos in range(mutrange):
    inserta=seqa[muta-pos:muta+mutrange-pos]
    lst=[]
    for i in inserta:
        lst+=aa[i],
    insert=list(itertools.product(*lst))
    for i in range(len(insert)):
        insert[i]="".join(insert[i])
        mutslice+=wtslice[:30+3*(mutrange-pos)]+"".join(insert[i])+wtslice[30+6*mutrange-3*pos:],
results_yes=summarise(yes)
results_no=summarise(no)
results_both=results_yes.copy()
lines=make_table(results_yes)

print("Include unowned enzymes? (Y/N)")
while True:
    x=input(">>> ")
    if x.upper()=="Q":
        break
    elif x.upper()=="Y":
        results_both.update(results_no)
        lines=make_table(results_both)
        print("Type the name of an enzyme for more details, 'IUPAC' for a list of ambiguity codes or 'Q' to quit")
    elif x.upper()=="N":
        print("Type the name of an enzyme for more details, 'IUPAC' for a list of ambiguity codes or 'Q' to quit")
    elif True in [x.upper() in i.upper() for i in results_both.keys()]:
        temp=[x.upper()==i.upper() for i in results_both.keys()]
        if not True in temp:
            temp=[x.upper() in i.upper() for i in results_both.keys()]
        if temp.count(True)==1:
            found=False
            enz=list(results_both.keys())[temp.index(True)]
            for line in lines:
                if line[1]==enz:
                    found=True
                    temp=[type(i) for i in line[5]]
                    if str in temp:
                        muts=" and ".join([i for i in line[5] if type(i)==str])
                        if list in temp:
                            muts+=" and "+" or ".join([i[0] for i in line[5] if type(i)==list])
                    else:
                        muts=" or ".join([i[0] for i in line[5] if type(i)==list])
                    print("Rank:",line[0],"/// Mutations:",muts,"/// Old -> new bands:",str(line[6])[1:-1].replace("'",""),"->",str(line[7])[1:-1].replace("'",""))
            if not found:
                print("This enzyme has been filtered out because its Log10 separation score is "+str(round(max([results_both[enz][i][0] for i in results_both[enz].keys()]),3)))
        else:
            print("Did you mean "+" or ".join([list(results_both.keys())[i] for i in range(len(temp)) if temp[i]==True])+"?")
    elif x.upper()=="IUPAC":
        print("|IUPAC code|"+"Meaning".center(16," ")+"|")
        for code in iupac.keys():
            if code=="N":
                print("|"+"N".center(10," ")+"|A or G or C or T|")
            else:
                print("|"+code.center(10," ")+"|"+" or ".join(iupac[code]).center(16," ")+"|")
