import os, sys, time
import yaml
import numpy as np

from fnmatch       import fnmatch
from datetime      import datetime
from optparse      import OptionParser
from Bio           import SeqIO, Seq, AlignIO
from Bio.SeqRecord import SeqRecord

from dateutils import getDoubleDate
from fileutils import parseList, getMapping


###########################################################################################################
# Extract pol sequences from csv and fasta file (if provided) according to metadata criteria
#   - Extract ALL matching sequences
#   - Assumes missing values are indicated by "NA" in the csv file
#   - Use "*" to select all within a group (except missing values)
#   - Wildcards do not apply for included, extra and excluded ids (-I, -X and -E options). For these use exact ids.
#   - Selection criteria can be passed as either input arguments or a yaml configuration file
#   - Input arguments supersede the yaml configuration file
#
# THIS IS TO SELECT SEQUENCES THAT FULFILL CERTAIN CRITERIA, NOT TO SUBSAMPLE THEM!
#
# Select genomes by:
#
#   -I, --include    : Included ids                      (file with ids in the first column) OR (comma-separated list)
#   -E, --exclude    : Excluded ids                      (file with ids in the first column) OR (comma-separated list)
#   -D, --daterange  : Date range                        [yyyy/mm/dd-yyyy/mm/dd) OR [float-float)
#   -C, --country    : Country                           (comma-separated list)
#   -P, --province   : Province                          (comma-separated list)
#   -L, --length     : Minimum ungapped sequence length  (integer)
#
#
# Ignored metadata (perhaps added in the future)
#
#
#
################################################################################################################################
# Parameters
################################################################################################################################


usage = "usage: %prog [option]"
parser = OptionParser(usage=usage)


################################################################################################################################
# General options

parser.add_option("-i","--inputfile",
                  dest = "inputfile",
                  default = "",
                  metavar = "path",
                  help = "Path to input csv file [required]")

parser.add_option("-a","--alignfile",
                  dest = "alignfile",
                  default = "",
                  metavar = "path",
                  help = "Path to input alignment file [optional]")

parser.add_option("-o","--outputpath",
                  dest = "outputpath",
                  default = "",
                  metavar = "path",
                  help = "Path to save output file in [required]")

parser.add_option("-p","--prefix",
                  dest = "prefix",
                  default = "",
                  metavar = "path",
                  help = "Prefix for the output file(s) [required]")

parser.add_option("-c","--configfile",
                  dest = "configfile",
                  default = "",
                  metavar = "yaml file",
                  help = "Configure file to use for options [default = %default]")

parser.add_option("-d","--dbsep",
                  dest = "dbsep",
                  default = ",",
                  metavar = "string",
                  help = "Separator between columns in output database [default: %default] [optional]")

parser.add_option("-f","--fastasep",
                  dest = "fasep",
                  default = "|",
                  metavar = "string",
                  help = "Separator between fields in fasta sequence id's [default: %default] [optional]")

parser.add_option("-m","--msaonly",
                  dest = "msaonly",
                  default = False,
                  action = "store_true",
                  help = "Write only the output alignment, not the csv files [default: %default] [optional]")

parser.add_option("-u","--updateids",
                  dest = "updateids",
                  default = False,
                  action = "store_true",
                  help = "Update sequence ids [default: %default] [optional]")

parser.add_option("-s","--seqidformat",
                  dest = "seqidformat",
                  default = "EBOV|{id}|{accession}|{country}|{province}|{date}",
                  metavar = "string",
                  help = "Format of the output sequence ids if updateids is True [default: %default] [optional]")

################################################################################################################################
# Selection criteria

parser.add_option("-I","--include",
                  dest = "include",
                  default = "",
                  metavar = "comma-separated list of strings/path",
                  help = "Comma-separated list of sequence ids to select or csv file with sequence ids in first column [optional]")

parser.add_option("-E","--exclude",
                  dest = "exclude",
                  default = "",
                  metavar = "comma-separated list of strings/path",
                  help = "Comma-separated list of sequence ids to exclude or csv file with sequence ids in first column [optional]")

parser.add_option("-D","--daterange",
                  dest = "daterange",
                  default = "",
                  metavar = "range",
                  help = "Range of 2 dates to define lower (inclusive) and upper (exclusive) bounds. Either floating point or yyyy/mm/dd (e.g. 2011.5-2012.5 or 2011/01/01-2012/02/10) [optional]")

parser.add_option("-C","--country",
                  dest = "country",
                  default = "",
                  metavar = "comma-separated list of strings",
                  help = "Comma-separated list of countries to select [optional]")

parser.add_option("-P","--province",
                  dest = "province",
                  default = "",
                  metavar = "comma-separated list of strings",
                  help = "Comma-separated list of provinces to select [optional]")

parser.add_option("-L","--length",
                  dest = "length",
                  default = "",
                  metavar = "integer",
                  help = "Minimum length of sequences to select (ungapped, inclusive, use alignment if possible) [optional]")

(options,args) = parser.parse_args()

inputfile  = os.path.abspath(options.inputfile)
alignfile  = os.path.abspath(options.alignfile) if options.alignfile != "" else ""
outputpath = os.path.abspath(options.outputpath)+"/"
prefix     = options.prefix
dbsep      = options.dbsep
fasep      = options.fasep
msaonly    = options.msaonly
updateids  = options.updateids
seqidformat= options.seqidformat

###########################

configfile = os.path.abspath(options.configfile) if options.configfile != "" else ""

inputpars  = dict()
inputpars["include"]   = options.include
inputpars["exclude"]   = options.exclude
inputpars["daterange"] = options.daterange
inputpars["country"]    = options.country
inputpars["province"]  = options.province
inputpars["length"]    = options.length




# Column mapping
fieldmapping = {"id"            : "sample_id",  
                "accession"     : "accession",              
                "date"          : "date",                                
                "country"       : "country",
                "province"      : "location",
                "virus"         : "virus",
                "sample_lab"    : "sample_lab",
                "sequencing_lab": "sequencing_lab",
                "platform"      : "platform",
                "GP82"          : "GP82",
                "reference"     : "reference"} 

datefmt    = "%Y/%m/%d"
datedigits = 5



################################################################################################################################


def readSeqIds(inputfile, idx=0, sep=','):
      """Read ids from file

      """

      idfile = open(inputfile,"r")

      ids = []
      for line in idfile:
            parts = line.strip().split(sep)
            ids.append(parts[idx])

      idfile.close()
      return ids
#


def readAlignment(inputfile, fmt, sep="|"):
      """Read alignment from file

      """

      alignfile = open(inputfile,"r")
      alignment = AlignIO.read(alignfile,fmt)
      alignfile.close()

      seqdict = dict()
      for seqrecord in alignment:
          idx   = seqrecord.id.find(sep)          
          seqid = seqrecord.id[idx+1:seqrecord.id.find(sep, idx+1)]
          #seq   = seqrecord.seq
          seqdict[seqid] = seqrecord

      return seqdict
#


def fuzzyStringMatch(candidate, matchlist):
    """Equivalent to (candidate in matchlist) but with fnmatch for wildcards in the list

    """

    return ( sum(list(map( lambda x : fnmatch(candidate,x) , matchlist ))) > 0 )
#


def getSequenceFields(datalist, mapping):
    """Return dictionary with fields of the datalist mapped with the mapping

    """

    seqdict = dict()
    for k in mapping.keys():
        seqdict[k] = datalist[mapping[k]]

    return seqdict
#


def printCriteria(selection):


    sys.stdout.write("\nSELECTION CRITERIA\n")
    for par in selection:

        if (isinstance(selection[par],tuple)):
            valuestr = "[%.5f,%.5f)" % selection[par]
        elif (isinstance(selection[par], list)):
            valuestr = ",".join(selection[par])
        elif (isinstance(selection[par], int)):
            valuestr = int(selection[par])
        else:
            valuestr = selection[par]

        sys.stdout.write("%15s : %s\n" % (par, valuestr))
#




################################################################################################################################
start = time.time()

##################################################
# Load config file
if (configfile != ""):
    config     = open(configfile,'r').read().replace("\t"," ")
    configpars = yaml.load(config)  

    #print(configpars)

    for par in inputpars.keys():
        if (inputpars[par] != ""):
            configpars[par] = inputpars[par]
else:
    configpars = inputpars


##################################################
# Process selection criteria
selection = dict()
for par in configpars:

    if (configpars[par] != ""):

        if (par == "include" or par == "exclude"):
            idfile = os.path.abspath(configpars[par])
            if (os.path.exists(idfile) and os.path.isfile(idfile)):
                #print("loading file")
                selection[par] = readSeqIds(idfile, idx=0, sep=dbsep)

            else:
                selection[par] = parseList(configpars[par])

        elif (par == "daterange"):

            try:
                floatdates = parseList(configpars["daterange"], float, sep="-")
            except Exception as A: 
                print("something")
                floatdates = parseList(configpars["daterange"], lambda f: getDoubleDate(datetime.strptime(f,datefmt)), sep="-")
            except Exception as B:
                sys.stdout.write("ERROR! Cannot parse date range\n")

            selection[par] = (min(floatdates), max(floatdates))

        elif (par == "length"):
            selection[par] = int(configpars[par])

        else:
            selection[par] = parseList(configpars[par])

printCriteria(selection)


# Read alignment
if (alignfile != ""):
     msaseqs = readAlignment(alignfile, "fasta")

# Make output folder
if (not os.path.exists(outputpath)):
    os.mkdir(outputpath)

evaluated = 0
matched   = 0
skipped   = 0
sequences = []

db        = open(inputfile,'r')
header    = db.readline().strip().split(dbsep)
mapping   = getMapping(header, fieldmapping)


if (prefix == ""):
    prefix = inputfile[inputfile.rfind('/')+1:inputfile.rfind('.')]

if (msaonly != True):
    outputdb  = open(outputpath+prefix+".csv", "w")
    skipfile  = open(outputpath+prefix+".skipped.csv","w")
    outputdb.write(dbsep.join(header)+"\n")
    skipfile.write(dbsep.join(header)+"\n")

for line in db:
    evaluated += 1
    row = line.strip().split(dbsep)
    seqid = row[mapping["id"]]
    
    ##################################################
    # Check selection criteria
    skip  = False
    extra = False
    for par in selection.keys():
        #sys.stdout.write("%s\t%s\t%s\n" % (par, row[mapping[par]], str(selection[par])))
    
        # ID in excludes
        if (par == "exclude"):
            #if (fuzzyStringMatch(row[mapping["id"]], selection[par]) == True):
            if (row[mapping["id"]] in selection[par]):
                skip = True
                continue

        # ID in includes
        elif (par == "include"):
            #if (fuzzyStringMatch(row[mapping["id"]], selection[par]) == False):
            if (row[mapping["id"]] not in selection[par]):
                skip = True
                continue

        # Date range
        elif (par == "daterange"):
            try: 
                date = getDoubleDate(datetime.strptime(row[mapping["date"]], "%Y-%m-%d"))

                if (date < selection["daterange"][0] or date >= selection["daterange"][1]):
                    raise ValueError()
            except: 
                skip = True
                continue
        
        # Length
        elif (par == 'length'):

            if (alignfile != ""):
                seq = msaseqs[seqid].seq.ungap("-").ungap("N") 
            else:
                sys.stdout.write("ERROR! No alignment file specified, cannot select by sequence length\n")
                raise RuntimeError
           

            if (len(seq) < selection["length"]):
                skip = True
                continue

        # Other comma-separated fields
        else:

            if (row[mapping[par]] == 'NA' or fuzzyStringMatch(row[mapping[par]], selection[par]) == False):
                skip = True
                continue

            #if (selection[par] != [] and row[mapping[par]] not in selection[par]):
            #    skip = True
            #    continue

    # Sequence is in extra list and is always included
    if (extra):
        skip = False

    ##################################################
    # Skip or save
    if (skip):
        if (msaonly == False):
            skipfile.write(dbsep.join(row)+"\n")
        skipped += 1
    else:
        matched += 1

        # Save metadata
        if (msaonly != True):
            outputdb.write(dbsep.join(row)+"\n")

        # Save alignment
        if (alignfile != ""):            
            seq = msaseqs[seqid]

            # Update sequence ids
            if (updateids):
                seqdict = getSequenceFields(row, mapping)
                newid   = seqidformat.format(**seqdict)
                seq.id  = newid
                seq.description = ""
            sequences.append(seq)


db.close()

if (msaonly != True):
    outputdb.close()
    skipfile.close()

# Write sequence file
if (len(sequences) > 0):
    SeqIO.write(sequences,outputpath+prefix+".fas","fasta")


# Print some statistics
sys.stdout.write("\nRESULTS SUMMARY\n")
sys.stdout.write("%10d sequence(s) evaluated\n" % evaluated)
sys.stdout.write("%10d matching sequences found\n" % matched)
sys.stdout.write("%10d sequences skipped\n" % skipped)
sys.stdout.write("%10d matching sequences extracted from alignment\n\n" % len(sequences))


end = time.time()
sys.stdout.write("Total time taken: "+str(end-start)+" seconds\n")
