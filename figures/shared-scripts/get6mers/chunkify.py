import re, random
from sys import exit
from sys import argv

n = 6
k = 2

f = open(argv[1],"r")
#f = open("lang/samples/english.txt","r")
#f = open("lang/samples/xhosa.txt","r")


ls = ""

for l in f:
	l = re.sub(r'[^a-z]+', '_', l.lower())
	l = re.sub(r'^[^a-z]+', '', l)
	l = re.sub(r'[^a-z]+$', '', l)
	l = l.strip()
	if l != "" :
		ls = ls + l + "_"

chk = []

i = 0
while i+n < len(ls) :
#	print(ls[i:i+n])
	i = i+n
	chk.append(ls[i:i+n]) 

	
# print( len(set(chk)) )

x = list( set( chk ) )
random.shuffle( x )

chk = x[0:1000]	

for c in chk:
	print(c)

exit()

def sim(a,b):
	r = 0
	for i in range( min( [len(a),len(b)] ) ):
		if a[i] != b[i] :
			r = 0
		else :
			r += 1
	return r

print( "graph { " )
for i in range( len( chk ) ):
	for j in range( i+1, len(chk) ):
		if sim( chk[i], chk[j] ) >= 3 :
			print( chk[i] + " -- "+ chk[j] )
print( "}" )