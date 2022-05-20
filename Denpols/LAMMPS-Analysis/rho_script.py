import os 

nbb=[10,16,25,36,63,100]
#Ls=[30,48,75,108,189,500] #for G3 and G5
#Ls=[36.45,64.8,75,72,189,500] #For G4
lp=10
gen=int(input("Generation?"))
file="G{}_sq.sh".format(gen)
g=open(file,"w")
frames=1000
#framestep=[1e5,1e5,1e5,1e5,1e5,1e6] #For G3
#framestep=[1e4,1e4,1e5,1e5,2e5,1e6] #For G4
framestep=[2e4,2e4,1e6,1.4e5,1e7,1e6] #For G5

skip=2**(gen+1)+1
foldername='./radial_density'
if not os.path.exists(foldername):
    os.makedirs(foldername)
for i in range(len(nbb)):
	rdf_name=foldername +'/G{}{}_rdf'.format(gen,nbb[i])
	sq_name=foldername +'/G{}{}'.format(gen,nbb[i])
	rdf_name_lin=foldername +'/G{}{}_rdf_lin'.format(gen,nbb[i])
	sq_name_lin=foldername +'/G{}{}_lin'.format(gen,nbb[i])



	myL=Ls[i]
	filename="./sq_trjs/G{}{}.trj".format(gen,nbb[i])
	line="python3 quick_rdf.py {} {} {} {} {} {} {}".format(filename,myL,1,int(framestep[i]),frames,rdf_name,sq_name)
	line1="python3 quick_rdf.py {} {} {} {} {} {} {}".format(filename,myL,skip,int(framestep[i]),frames,rdf_name_lin,sq_name_lin)


	g.write(line)
	g.write("\n")
	g.write(line1)
	g.write("\n")