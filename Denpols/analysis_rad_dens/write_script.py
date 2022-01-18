import os	

gen=[3,4,5]
bb=[10,16,25,36,63,100]
#box_g3=[30,48,75,108,189,500]
#box_g4=[36.45,64.8,75,72,189,500] 
#box_g5=box_g3
lp=10
path=os.getcwd()
print("Current working directory is %s"%path)

g=open("radial_exec.sh","w")
for i in gen:
	count=0
	for j in bb:
		foldername=path+"/G{}{}".format(i,j)
		if i==3:
			g.write("python3 analysis_rho_r.py %s %d %d %d 1\n\n"%(foldername,i,j,lp))
			count+=1	
		elif i==4:
			g.write("python3 analysis_rho_r.py %s %d %d %d 1\n\n"%(foldername,i,j,lp))
			count+=1
		elif i==5:
			g.write("python3 analysis_rho_r.py %s %d %d %d 1\n\n"%(foldername,i,j,lp))
			count+=1
		
		
g.close()
