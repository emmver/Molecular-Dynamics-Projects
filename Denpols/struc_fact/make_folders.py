import os 



path=os.getcwd()
print("Current working directory is %s"%path)
gen=[3,4,5]
bb=[10,16,25,36,63,100]

for i in gen:
	for j in bb:
		foldername=path+"/G{}{}".format(i,j)


		try: 
			os.mkdir(foldername)
		except OSError: 
			print("Creation of directory failed %s"%(foldername))
		else:
			print("Success %s"%(foldername))
			
			
