
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import os
import re
os.chdir('/home/yangfang/mywork/bacteria3/')
#open a file
csvfile = open('abb_600_names_taxonomy.txt')
#initialize a list get the class
ret = []
# read this file by tab
for line in csvfile:
	x = line.split("\t")
	for i in x[2:]:
		if i != 'NA':
			name = i 
			break
	#get the class names
	ret.append(name)
#get the unique class names

ret = list(set(ret))

csvfile.close()


#open the color file
color = open('color_hcl.txt')
color_list = []
for line in color:
	col = line.split("\t")
	color_list.append(col[0].strip())

color_list = color_list[:len(ret)]
color.close()
# create a dict 
color_dict = dict(zip(ret,color_list))

# a function to get the color code
def ch_color(phylum):
	for (k,v) in color_dict.items():
		if phylum == k:
			return v

#open files			
csvfile2 = open('271euk_taxonomy.txt')
f = open("color_hcl_abb.txt","w")
#write the file
for line in csvfile2:
	x = line.split("\t")
	for i in x[1:]:
		if i != 'NA':
			name = i 
			break
	color_code = ch_color(name)
	# write by format: Homo_sapiens   range   #00ff00   important1
	f.write("%s\t%s\t%s\t%s\n"%(x[0],"range",color_code,name))

csvfile.close()
f.close()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


