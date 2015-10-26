import urllib2
import re
url = "http://www.kegg.jp/kegg/catalog/org_list.html"
content = urllib2.urlopen(url).read()

#Extract the organisms abbreviation 
abb = re.findall("(?<=\\?org=)(.*?)(?='>)",content)

#Extract the organisms fullname
fullname = re.findall("(?<=www_bfind\\?T[\\d]{5}'>).*(?=</a>)",content)
list1 = []
fp = open("/home/yangfang/organisms_abb_fulllnames.txt","w")
for i in range(len(abb)):
	list1.append(str(abb[i])+"\t"+str(fullname[i]))
	
fp.write("\n".join(list1))
fp.close()

		

