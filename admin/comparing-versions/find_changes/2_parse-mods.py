import csv

print("Preparing file comparison lists (basicfiles.txt AND modfiles.txt)")

f = open('changedfiles.txt', "r")
files = f.read()
files = files.split('\n')

basicf = open('basicfiles.txt', 'a')
modf = open('modfiles.txt', 'a')
for line in range(0,len(files)):
	data = files[line].split(' ')
	if len(data) == 5:
		fileA = data[1]
		fileB = data[3]
		basicf.write(fileA+'\n')
		modf.write(fileB+'\n')
print("Files ready to compare")
