import numpy as np
import matplotlib.pyplot as plt
import os
import pandas

instances="20_USA-road-d.BAY.gr 20_USA-road-d.COL.gr 20_USA-road-d.NY.gr 40_USA-road-d.BAY.gr 40_USA-road-d.COL.gr 40_USA-road-d.NY.gr 60_USA-road-d.BAY.gr 60_USA-road-d.COL.gr 60_USA-road-d.NY.gr 80_USA-road-d.BAY.gr 80_USA-road-d.COL.gr 80_USA-road-d.NY.gr 100_USA-road-d.BAY.gr 100_USA-road-d.COL.gr 100_USA-road-d.NY.gr 120_USA-road-d.BAY.gr 120_USA-road-d.COL.gr 120_USA-road-d.NY.gr 140_USA-road-d.BAY.gr 140_USA-road-d.COL.gr 140_USA-road-d.NY.gr 160_USA-road-d.BAY.gr 160_USA-road-d.COL.gr 160_USA-road-d.NY.gr 180_USA-road-d.BAY.gr 180_USA-road-d.COL.gr 180_USA-road-d.NY.gr 200_USA-road-d.BAY.gr 200_USA-road-d.COL.gr 200_USA-road-d.NY.gr 250_USA-road-d.BAY.gr 250_USA-road-d.COL.gr 250_USA-road-d.NY.gr 300_USA-road-d.BAY.gr 300_USA-road-d.COL.gr 300_USA-road-d.NY.gr 350_USA-road-d.BAY.gr 350_USA-road-d.COL.gr 350_USA-road-d.NY.gr"

tab_files= instances.split(" ")

tab_methods=['bc','s','d','pc','pck']

dico_methods=dict()

for s in tab_methods:
	value=dict()
	time=dict()
	list_files= os.listdir('./ECMA/results/'+s)
	for file in list_files:
		title=file.split('_')[0]+'_'+file.split('_')[1]
		if title in value.keys():
			value[title].append(float(file.split('_')[3]))
		else:
			value[title]=[float(file.split('_')[3])]

		t=float((file.split('_')[4])[:-4])
		if t>1200:
			t=1200
		if title in time.keys():
			time[title].append(t)
		else:
			time[title]=[t]

	dico_methods[s+"_val"]=np.zeros(len(tab_files))
	dico_methods[s+"_time"]=np.zeros(len(tab_files))

	i=0
	for file in tab_files:
		if file in value.keys():
			dico_methods[s+"_val"][i]=sum(value[file])*1./len(value[file])
			dico_methods[s+"_time"][i]=sum(time[file])*1./len(time[file])
		i+=1

	id_val=np.argsort(dico_methods[s+"_val"])
	id_t=np.argsort(dico_methods[s+"_time"])
	
	plt.plot( (dico_methods[s+"_time"])[id_t], np.arange(len(tab_files)), label=s)

	
#resolution exacte
file=open("ECMA/results/re/res_exact.txt",'r')
lines=file.readlines()
file.close()

tab_val=[]
tab_t=[]
for file in tab_files:
	for line in lines:
		file_line=line.split(' ')[0]
		if file_line == file:
			t=float(line.split(' ')[1])
			val=float(line.split(' ')[2])
			tab_val.append(val)
			tab_t.append(t)
dico_methods['re_val']=np.array(tab_val)
dico_methods['re_time']=np.array(tab_t)

id_val=np.argsort(dico_methods['re_val'])
id_t=np.argsort(dico_methods['re_time'])

plt.plot( (dico_methods['re_time'])[id_t], np.arange(len(lines)),label='resolution exacte' )

#resolution exacte
file=open("ECMA/results/h/res_h.txt",'r')
lines=file.readlines()
file.close()

tab_val=[]
tab_t=[]
for file in tab_files:
	for line in lines:
		file_line=line.split(' ')[0]
		if file_line == file:
			t=float(line.split(' ')[1])
			val=float(line.split(' ')[2])
			tab_val.append(val)
			tab_t.append(t)
dico_methods['h_val']=np.array(tab_val)
dico_methods['h_time']=np.array(tab_t)

id_val=np.argsort(dico_methods['h_val'])
id_t=np.argsort(dico_methods['h_time'])

plt.plot( (dico_methods['h_time'])[id_t], np.arange(len(lines)) , label='heuristique')

plt.xlabel('Temps (s)')
plt.ylabel("Nombre d'instances résolues")
plt.legend(bbox_to_anchor=(1,1), loc="upper left")
plt.show()

"""
tab_methods+=['re','h']

PR=np.zeros((len(tab_methods),len(tab_files)))
for i in range(len(tab_methods)):
	PR[i]=100*np.abs(dico_methods[tab_methods[i]+'_val']-dico_methods['d_val'])*1./dico_methods['d_val']
			
			


dataframe = pandas.DataFrame(PR[1]) 
dataframe.to_csv("PR.csv")

dataframe = pandas.DataFrame(PR) 
dataframe.to_csv("gap.csv")
"""
	
