import os, sys, time, collections
from optparse import OptionParser
from Bio import AlignIO


# Get counts for characters in columns of MSA
################################################################################################################################
# Parameters

usage = "usage: %prog [option]"
parser = OptionParser(usage=usage)

parser.add_option("-i","--inputfile",
                  dest = "inputfile",
                  default = "",
                  metavar = "path",
                  help = "Path to input alignment [required]")

parser.add_option("-o","--outputpath",
                  dest = "outputpath",
                  default = "",
                  metavar = "path",
                  help = "Path to save output file in [required]")

parser.add_option("-n","--normalise",
                  dest = "normalise",
                  action = "store_true",
                  help = "Normalise each column in the alignment [default: %default] [optional]")

(options,args) = parser.parse_args()

inputfile  = os.path.abspath(options.inputfile)
outputpath = os.path.abspath(options.outputpath)+"/"
normalise  = options.normalise
alphabet   = ["A", "C", "G", "T", "R", "Y", "W", "S", "K", "M", "D", "V", "H", "B", "N", "-", "?"]


################################################################################################################################
start = time.time()

if (not os.path.exists(outputpath)):
    os.mkdir(outputpath)

outfile = open(outputpath+inputfile[inputfile.rfind("/")+1:inputfile.rfind(".")]+".hist.csv","w")
outfile.write(",".join(alphabet)+"\n")    

msa = AlignIO.read(inputfile, "fasta")
for i in range(0, msa.get_alignment_length()):
	col = msa[:,i]

	# Old slow way to get column histogram
 	#coldict = dict((c, col.count(c)) for c in col)

 	# Faster way
	coldict = collections.Counter(col.upper())
 	
 	# Mising characters not in alphabet
	for c in coldict.keys():
		if (c not in alphabet):
			sys.stdout.write(c+" not in alphabet!\n")

	total = 0
	line  = []
	for c in alphabet:
		if (c in coldict.keys()):
			line.append(coldict[c])
			total += coldict[c]
		else:
			line.append(0)

	if (normalise):
		for i in range(0,len(line)):
			line[i] /= total

	outfile.write(",".join(map(str, line))+"\n")
	#sys.stdout.write("\t".join(map(str, line))+"\t"+str(total)+"\n")
#
outfile.close()

end = time.time()
sys.stdout.write("Total time taken: "+str(end-start)+" seconds\n")



