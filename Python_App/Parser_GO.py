import csv
import os
import re

os.chdir("/home/izem/Documents")
file = open("/home/izem/Downloads/go.obo", "r")
file.readline()

list_BP_myo=[]
list_BP_angio=[]


for term in file.read().split("[Term]"):
    if term.__contains__("biological_process") and term.lower().__contains__("muscle") and term.__contains__("obsolete")==False : #i use lower() to avoid missing Muscle with M
        list_BP_myo.append(term[5:16].strip("\n")) #i remove le saut de ligne and i take only the go id locate at 5 to 16
    if term.__contains__("biological_process") and term.lower().__contains__("angiogenesis") and term.__contains__("obsolete")==False :
        list_BP_angio.append(term[5:16].strip("\n"))



with open('muscle.csv', 'w') as f:
    write = csv.writer(f)
    write.writerow(list_BP_myo)
    write.writerow(list_BP_myo)

with open('angio.csv', 'w') as f:
    write = csv.writer(f)
    write.writerow(list_BP_angio)
    write.writerow(list_BP_angio)



