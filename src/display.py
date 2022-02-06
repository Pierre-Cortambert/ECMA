import numpy as np
import matplotlib.pyplot as plt
import os

tab_files= os.listdir('../instances')

tab_methods=['bc','s','d','pc','pck']

dico_methods=dict()

for s in tab_methods:
	value=dict()
	time=dict()
	list_files= os.listdir('./ECMA/results/'+s)
	for file in list_files:
		if file != '.directory':

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

	dico_methods[s+"_val"]=np.zeros(len(tab_files)-1)
	dico_methods[s+"_time"]=np.zeros(len(tab_files)-1)

	i=0
	for file in tab_files:
		if file != '.directory':
			if file in value.keys():
				dico_methods[s+"_val"][i]=sum(value[file])*1./len(value[file])
				dico_methods[s+"_time"][i]=sum(time[file])*1./len(time[file])
			i+=1

	id_val=np.argsort(dico_methods[s+"_val"])
	id_t=np.argsort(dico_methods[s+"_time"])
	
	plt.plot( (dico_methods[s+"_time"])[id_t], np.arange(len(tab_files)-1), label=s)

	
#resolution exacte
file=open("../results/re/res_exact.txt",'r')
lines=file.readlines()
file.close()

tab_val=[]
tab_t=[]
for line in lines:
	t=float(line.split(' ')[1])
	val=float(line.split(' ')[2])
	tab_val.append(val)
	tab_t.append(t)
dico_methods['re_val']=np.array(tab_val)
dico_methods['re_time']=np.array(tab_t)

id_val=np.argsort(dico_methods['re_val'])
id_t=np.argsort(dico_methods['re_time'])

plt.plot( (dico_methods['re_time'])[id_t], np.arange(len(lines)),label='re' )

#resolution exacte
file=open("../results/h/res_h.txt",'r')
lines=file.readlines()
file.close()

tab_val=[]
tab_t=[]
for line in lines:
	t=float(line.split(' ')[1])
	val=float(line.split(' ')[2])
	tab_val.append(val)
	tab_t.append(t)
dico_methods['re_val']=np.array(tab_val)
dico_methods['re_time']=np.array(tab_t)

id_val=np.argsort(dico_methods['re_val'])
id_t=np.argsort(dico_methods['re_time'])

plt.plot( (dico_methods['re_time'])[id_t], np.arange(len(lines)) , label='h')

plt.xlabel('Temps (s)')
plt.ylabel("Nombre d'instances résolues")
plt.legend(bbox_to_anchor=(1,1), loc="upper left")
#plt.legend(tab_methods+['re','h'])
plt.show()
	
	
