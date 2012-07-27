import os

files = os.listdir("input")

total = str(len(files))

for i,filename in enumerate(files):
	print filename
	print str(i)+" of "+total
	os.system("./Kavosh 6 "+filename)

