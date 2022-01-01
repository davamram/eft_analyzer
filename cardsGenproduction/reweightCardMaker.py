import argparse

parser = argparse.ArgumentParser(description='Script to generate reweighting datacard for genproduction with MG5')
args = parser.parse_args()

reweightCard = open("_reweight_crad.dat","w+")

reweightCard.write("#****************************************************************** \n")
reweightCard.write("#                       Reweight Module                           * \n")
reweightCard.write("#****************************************************************** \n \n")
reweightCard.write("change rwgt_dir rwgt \n \n")

switch = "yes"

WilsonCoeff = ["cptb","cptbi","ctw","ctwi","cbw","cbwi"]
ParamCardCorresp = ["8","9","10","12","14","15"]



while switch == "yes":
    switch = ""
    branchName = input("Enter branch name:")
    reweightCard.write("launch --rwgt_name="+branchName+"\n")

    for i in range(len(WilsonCoeff)):
        CoeffValue = input(WilsonCoeff[i] + " = ")
        reweightCard.write("set DIM6 " + ParamCardCorresp[i] + " " + CoeffValue + "    #" + WilsonCoeff[i] + "\n")

    reweightCard.write("set decay 6 auto \n \n")
    switch = input("Do you want to add a new branch? [yes/no]")

reweightCard.close()

            

