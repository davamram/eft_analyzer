#########################################################################
##              Script to write reweighting card for                   ##
##                      genproductions MG5                             ##
#########################################################################



import argparse

parser = argparse.ArgumentParser(description='Script to generate reweighting datacard for genproduction with MG5')
args = parser.parse_args()

def BranchWritter(name,cptb,cptbi,ctw,ctwi,cbw,cbwi):
    reweightCard.write("launch --rwgt_name="+ name +"\n")
    reweightCard.write("set DIM6 8 " + str(cptb) + "\n")
    reweightCard.write("set DIM6 9 " + str(cptbi) + "\n")
    reweightCard.write("set DIM6 10 " + str(ctw) + "\n")
    reweightCard.write("set DIM6 12 " + str(ctwi) + "\n")
    reweightCard.write("set DIM6 14 " + str(cbw) + "\n")
    reweightCard.write("set DIM6 15 " + str(cbwi) + "\n")
    reweightCard.write("set decay 6 auto \n \n")

def BranchNamer(Wname,Wvalue):
    if Wvalue == 0: return ""
    minusplus = "_p"
    if Wvalue<0: minusplus = "_m"
    return "_" + Wname + minusplus + str(Wvalue).removeprefix("-")


reweightCard = open("reweight_crad.dat","w+")

reweightCard.write("#****************************************************************** \n")
reweightCard.write("#                       Reweight Module                           * \n")
reweightCard.write("#****************************************************************** \n \n")
reweightCard.write("change rwgt_dir rwgt \n \n")

###################### Only in Grid Points ######################

# for i in range(-10,11,5):
#     if i==0: continue
#     minusplus = "p"
#     value = str(i)
#     if i<0: minusplus = "m"
#     if len(str(i)) > 1: value = str(i).removeprefix("-")
#     name = "cptb_" + minusplus + value
#     BranchWritter(name,str(i),"0","0","0","0","0")

# for i in range(-10,11,5):
#     if i==0: continue
#     minusplus = "p"
#     value = str(i)
#     if i<0: minusplus = "m"
#     if len(str(i)) > 1: value = str(i).removeprefix("-")
#     name = "cptbi_" + minusplus + value
#     BranchWritter(name,"0",str(i),"0","0","0","0")

# for i in range(-5,6,1):
#     if (i==0 or i==-4 or i==-3 or i==3 or i==4): continue
#     minusplus = "p"
#     value = str(i)
#     if i<0: minusplus = "m"
#     if len(str(i)) > 1: value = str(i).removeprefix("-")
#     name = "ctwi_" + minusplus + value
#     BranchWritter(name,"0","0","0",str(i),"0","0")

# for i in range(-2,3,1):
#     if i==0: continue
#     minusplus = "p"
#     value = str(i)
#     if i<0: minusplus = "m"
#     if len(str(i)) > 1: value = str(i).removeprefix("-")
#     name = "ctw_" + minusplus + value
#     BranchWritter(name,"0","0",str(i),"0","0","0")

# for i in range(-5,6,1):
#     if (i==0 or i==-4 or i==-3 or i==3 or i==4): continue
#     minusplus = "p"
#     value = str(i)
#     if i<0: minusplus = "m"
#     if len(str(i)) > 1: value = str(i).removeprefix("-")
#     name = "cbwi_" + minusplus + value
#     BranchWritter(name,"0","0","0","0","0",str(i))

# for i in range(-2,3,1):
#     if i==0: continue
#     minusplus = "p"
#     value = str(i)
#     if i<0: minusplus = "m"
#     if len(str(i)) > 1: value = str(i).removeprefix("-")
#     name = "cbw_" + minusplus + value
#     BranchWritter(name,"0","0","0","0",str(i),"0")

###################### Whole Grid Points ######################
# In this part, we create the points in the parameter space   #
#  for ctw, ctwi, cbw, cbwi = {-5,-4,-3,-2,-1,0,1,2,3,4,5}    #
#         for cptb and cptbi = {-10,-5,0,5,10}                #
#      and taking acount of the cuts defined below            #
###############################################################




count = 0
cuts1 = [-4,-3,-1,1,3,4] #Cuts on values of ctw, ctwi, cbw, cbwi
cuts2 = [] #Cuts on values of cptb and cptbi
listOfProcesses = []

for cptb in range(-10,11,5):
    if cptb in cuts2:
        name_cptb = ""
        cptb = 0
        continue
    name_cptb = BranchNamer("cptb",cptb)

    for cptbi in range(-10,11,5):
        if cptbi in cuts2:
            name_cptbi = ""
            cptbi = 0
            continue
        name_cptbi = BranchNamer("cptbi",cptbi)

        for ctw in range(-5,6,1):
            if ctw in cuts1:
                name_ctw = ""
                ctw = 0
                continue
            name_ctw = BranchNamer("ctw",ctw)

            for ctwi in range(-5,6,1):
                if ctwi in cuts1:
                    name_ctwi = ""
                    ctwi = 0
                    continue
                name_ctwi = BranchNamer("ctwi",ctwi)

                for cbw in range(-5,6,1):
                    if cbw in cuts1:
                        name_cbw = ""
                        cbw = 0
                        continue
                    name_cbw = BranchNamer("cbw",cbw)

                    for cbwi in range(-5,6,1):
                        if cbwi in cuts1:
                            name_cbwi = ""
                            cbwi = 0
                            continue
                        name_cbwi = BranchNamer("cbwi",cbwi)

                        process = (name_cptb + name_cptbi + name_ctw + name_ctwi + name_cbw + name_cbwi).removeprefix("_")

                        if process == "":
                            process = "SM"
                            ctwi = "1e-10"
                        
                        reweightCard.write("#Process nb = " + str(count+1) + "\n")

                        BranchWritter(process,cptb,cptbi,ctw,ctwi,cbw,cbwi)

                        count += 1

print("Number of Processes = ", count)


######################Make your own process manually######################

# switch = input("Do you want to add a new branch? [yes/no]: ")
# WilsonCoeff = ["cptb","cptbi","ctw","ctwi","cbw","cbwi"]
# ParamCardCorresp = ["8","9","10","12","14","15"]
# while switch == "yes":
#     switch = ""
#     branchName = input("Enter branch name:")
#     reweightCard.write("launch --rwgt_name="+branchName+"\n")

#     for i in range(len(WilsonCoeff)):
#         CoeffValue = input(WilsonCoeff[i] + " = ")
#         reweightCard.write("set DIM6 " + ParamCardCorresp[i] + " " + CoeffValue + "    #" + WilsonCoeff[i] + "\n")

#     reweightCard.write("set decay 6 auto \n \n")
#     switch = input("Do you want to add a new branch? [yes/no]")

reweightCard.close()

            

