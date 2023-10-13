
print("Parsing list of differing files into \"uniquefiles.txt\" and \"changedfiles.txt\" ")

f = open('file_diffs.log', "r")
files = f.read()
files = files.split('\n')

unqs = open('uniquefiles.txt', 'a')
chng = open('changedfiles.txt', 'a')
for line in range(0,len(files)):
	if files[line].startswith("Only"):
		unqs.write(files[line]+'\n')
	else:
		chng.write(files[line]+'\n')
unqs.write('\nAll files checked')
chng.write('\nAll files checked')
