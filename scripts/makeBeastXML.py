import os, sys, io, yaml, json, datetime
import numpy as np
from fnmatch import fnmatch
from optparse import OptionParser
from Bio import SeqIO
from dateutils import getDoubleDate

usage = "usage: %prog [option]"
parser = OptionParser(usage=usage)

parser.add_option("-i","--inputpath",
                  dest = "inputpath",
                  default = "",
                  metavar = "path",
                  help = "Path to input files [required]")

parser.add_option("-c","--config",
                  dest = "config",
                  default = "*.cfg",
                  metavar = "",
                  help = "Pattern to match for config files [default = %default]")

parser.add_option("-x","--template",
                  dest = "template",
                  default = "",
                  metavar = "path",
                  help = "Path to template XML file [required]")

parser.add_option("-o","--outputpath",
                  dest = "outputpath",
                  default = "",
                  metavar = "path",
                  help = "Path to save output file in [required]")

parser.add_option("-n","--name",
                  dest = "name",
                  default = "",
                  metavar = "path",
                  help = "Name of the runs [required]")

parser.add_option("-s","--seed",
                  dest = "seed",
                  default = "127",
                  metavar = "integer",
                  help = "Seed or comma separated list of seeds to use [required]")

parser.add_option("-t","--threads",
                  dest = "threads",
                  default = "2",
                  metavar = "integer",
                  help = "Threads to use [required]")

parser.add_option("-q","--queue",
                  dest = "queue",
                  default = "24",
                  metavar = "integer",
                  help = "Queue to submit to [required]")

(options,args) = parser.parse_args()

if (options.inputpath != ""):
	config         = options.config
	inputpath      = os.path.abspath(options.inputpath)+"/"
else:
	config         = options.config[options.config.rfind("/")+1:]
	inputpath      = os.path.abspath(options.config[:options.config.rfind("/")])+"/"

################################################################################################################################  


def writeAlignment(filename, output, msa_id="msa"):
      """
      	Writes alignment in BEAST format to output buffer
      	Return dictionary with sequence ids and dates
      """

      datedict = dict()
      tipdates = dict()
      for seq in SeqIO.parse(filename, "fasta"):		

            # Process date
            datestring = seq.id.split('|')[-1]
            (seqdate, sequpper, seqlower)  = getYearDate(datestring, datefmt = "%Y-%m-%d")

            #print("%s\t%f\t%f\t%f\n" % (datestring, seqdate, sequpper, seqlower))

            if (sequpper - seqlower != 0): 
                sys.stdout.write("Estimating sampling date for: %s\n" % seq.id)
                tipdates[seq.id] = [sequpper, seqlower]
            #

            datedict[seq.id] = seqdate

            # Write sequence
            output.write('\t\t\t<sequence id="%s:%s" taxon="%s" totalcount="4" value="%s"/>\n' % (seq.id, msa_id, seq.id, seq.seq))

      #

      return({'datetrait' : datedict, 'tipdates' : tipdates})
#

def checkDateTraitsEqual(datetrait1, datetrait2):

      if (len(datetrait1.keys()) != len(datetrait2.keys())):
            sys.stdout.write("Warning: Alignments do not have the same number of sequences\n")
            return False 

      for seq in datetrait1.keys():
            if (datetrait1[seq] != datetrait2[seq]):
                  sys.stdout.write("Dates not equal for %s (%.5f vs %.5f)\n")
                  return False

      return True
#

def writeDateTrait(dates, output):
      """
      	Writes dateTrait in BEAST format for sequence ids in date dictionary
      """

      maxdate = max(dates.values())
      mindate = min(dates.values())

      traits = []
      for seqid in dates:
            traits.append('\n\t\t\t\t\t%s=%.13f' %(seqid, dates[seqid]))
            if (dates[seqid] == maxdate): 
                  sys.stdout.write("Most recent sample: %s\n" % seqid)
            if (dates[seqid] == mindate): 
                  sys.stdout.write("Oldest sample: %s\n" % seqid)
      output.write(','.join(traits)+'\n')

#

def writeTipDatesSampling(tipdates, distOutput, opOutput, logOutput, maxDate=None):
      """
            Writes distribution and operators for estimating tip dates
      """

      for seqid in tipdates.keys():

            # Max date cannot be in the future
            if (maxDate != None and maxDate < tipdates[seqid][1]):
                  tipdates[seqid][1] = maxDate

            distOutput.write('\n\t\t\t<distribution id="tipDates:%s" monophyletic="false" spec="beast.math.distributions.MRCAPrior" tipsonly="true" tree="@Tree.t:tree">\n' % seqid)
            distOutput.write('\t\t\t\t<taxonset id="TaxonSet:%s" spec="TaxonSet">\n' % seqid)
            distOutput.write('\t\t\t\t\t<taxon id="%s" spec="Taxon"/>\n' % seqid)
            distOutput.write('\t\t\t\t</taxonset>\n')
            distOutput.write('\t\t\t\t<distr lower="%f" offset="0.0" spec="beast.math.distributions.Uniform" upper="%f"/>\n' % (tipdates[seqid][0], tipdates[seqid][1]))
            distOutput.write('\t\t\t</distribution>\n')

            #opOutput.write('\n\t<operator windowSize="1" spec="TipDatesRandomWalker" taxonset="@TaxonSet:%s" tree="@Tree.t:tree" weight="1.0"/>\n' % seqid)
            #opOutput.write('\n\t<operator windowSize="1" spec="SampledNodeDateRandomWalkerForZeroBranchSATrees" taxonset="@TaxonSet:%s" tree="@Tree.t:tree" weight="1.0"/>\n' % seqid)
            opOutput.write('\n\t\t<operator windowSize="1" spec="SampledNodeDateRandomWalker" taxonset="@TaxonSet:%s" tree="@Tree.t:tree" weight="1.0"/>\n' % seqid)

            logOutput.write('\t\t<log idref="tipDates:%s"/>\n' % seqid)
      #
#

def writeDatesCSV(dates, tipdates, output):

      maxdate = max(dates.values())

      output.write("id\tdate\tlower\tupper\n")
      for seqid in dates:

            if (seqid in tipdates.keys()):
                  upper = min(tipdates[seqid][1], maxdate)
                  lower = tipdates[seqid][0]
            else:
                  upper = dates[seqid]
                  lower = dates[seqid]

            output.write("%s\t%f\t%f\t%f\n" % (seqid, dates[seqid], lower, upper))
      #
#


# def getDoubleDate(date): 
#
#      dec31   = datetime.datetime.strptime(str(date.year)+"-12-31", "%Y-%m-%d")
#
#      datett  = date.timetuple()
#      dec31tt = dec31.timetuple()
#
#      return (date.year + float(datett.tm_yday-1)/dec31tt.tm_yday)
##



def getYearDate(datestring, datefmt = "%Y-%m-%d"):
      """
      	Takes in date in format 2014-01-01 and returns as year followed by fraction, eg. 2014.0

   		Also accounts for unknown month or day to return upper and lower bounds
      """

      if (datestring.count("-") == 2):
            date  = datetime.datetime.strptime(datestring, datefmt)
            upper = date 
            lower = date 
      # Only month
      elif (datestring.count("-") == 1):            
            date  = datetime.datetime.strptime(datestring+"-15", datefmt)
            upper = datetime.datetime.strptime("%d-%d-01" % (date.year, date.month), datefmt)
            if (date.month == 12):
                  day = 31
            else:
                  day = (datetime.date(date.year, date.month+1, 1)-datetime.date(date.year, date.month, 1)).days
            lower = datetime.datetime.strptime("%d-%d-%d" % (date.year, date.month, day), datefmt)
      # Only year
      elif (datestring.count("-") == 0):            
            date  = datetime.datetime.strptime(datestring+"-07-01", datefmt)
            upper = datetime.datetime.strptime(datestring+"-01-01", datefmt)
            lower = datetime.datetime.strptime(datestring+"-12-31", datefmt)
      else:
            sys.stdout.write("Unknown date format - %s\n" % datestring)

      #return ({'date' : getDoubleDate(date), 'upper' : getDoubleDate(upper), 'lower' : getDoubleDate(lower)})
      return (getDoubleDate(date), getDoubleDate(upper), getDoubleDate(lower))
#


def makeXMLFile(pars, template, outputpath=""):

	sys.stdout.write(pars["name"]+"...\n")
	formatpars = dict()
	for par in pars:
		formatpars['$'+par] = pars[par]
	output = template.format(**formatpars)

	if (outputpath == ""):
		outputpath = pars['outputpath']

	if (not os.path.exists(outputpath)):
		os.mkdir(outputpath)

	outfile = open(outputpath+"/"+pars["name"]+".xml", 'w')
	outfile.write(output)
	outfile.close()
#


def formatPars(pars):

      output_align_ig  = io.StringIO()
      output_align_cds = io.StringIO()
 
      dates_cds = writeAlignment(pars["alignment_cds"], output_align_cds, msa_id="cds")
      pars["alignment_cds"] = output_align_cds.getvalue()      

      # Check that dates and sequence id's are equal
      if ("alignment_ig" in pars.keys()): 
            dates_ig  = writeAlignment(pars["alignment_ig"],  output_align_ig, msa_id="ig")
            checkDateTraitsEqual(dates_ig['datetrait'], dates_cds['datetrait'])
            pars["alignment_ig"]  = output_align_ig.getvalue()      
          
      
      # Tip dates
      output_dates = io.StringIO()

      output_tipDates    = io.StringIO()
      output_opTipDates  = io.StringIO()
      output_logTipDates = io.StringIO()
      if ("datetrait" in pars.keys()):
            pars["datetrait"] = open(pars["datetrait"]).read().replace("\n",",\n")

            pars["tipdatesPriors"]    = ""
            pars["tipdatesOperators"] = ""
            pars["tipdatesLoggers"]   = ""
      else:
            writeDateTrait(dates_cds['datetrait'], output_dates)
            writeTipDatesSampling(dates_cds['tipdates'], output_tipDates, output_opTipDates, output_logTipDates, maxDate=max(dates_cds['datetrait'].values()))

            pars["datetrait"]  = output_dates.getvalue()

            pars["tipdatesPriors"]    = output_tipDates.getvalue()
            pars["tipdatesOperators"] = output_opTipDates.getvalue()
            pars["tipdatesLoggers"]   = output_logTipDates.getvalue()

      # Typetraits
      # pars["typetrait"] = open(pars["typetrait"]).read().replace("\t","=").replace("\n",",\n\t\t\t\t")

      if ("tmrca_mean" in pars.keys()):
            maxdate = max(dates['datetrait'].values())
            pars['tmrca_mean'] = maxdate - pars['tmrca_mean']

      output_align_ig.close()
      output_align_cds.close()
      output_dates.close()
      output_tipDates.close()
      output_opTipDates.close()
      output_logTipDates.close()


      # BDSKY parameters
      pars["samplingSliceDim"]   = pars["samplingProportionDimension"]-1
      if ("samplingProportion" not in pars.keys()):
            pars["samplingProportion"] = 0.001
      pars["samplingProportion"] = "0.0"+(" %.5f" % pars["samplingProportion"])*(pars["samplingProportionDimension"]-1)

      
      # Set origin limits (if origin present)
      if ("origin_max" in pars.keys()):
            pars["origin_min"]  = max(dates_cds['datetrait'].values())-min(dates_cds['datetrait'].values())
            pars["origin_max"]  = max(dates_cds['datetrait'].values())-getYearDate(str(pars['origin_max']))[0] 
            pars["origin_init"] = max(dates_cds['datetrait'].values())-getYearDate(str(pars['origin_init']))[0]
#


################################################################################################################################  

for filename in sorted(os.listdir(inputpath)):
      if (fnmatch(filename,config)):

            # Load config file
            configfile = open(inputpath+filename, 'r').read().replace("\t"," ")
            pars 	     = yaml.load(configfile)	

            # Set BEAST specific parameters
            queue      = int(options.queue)
            threads    = int(options.threads)
            seeds      = list(map(int, options.seed.split(',')))
            basename   = pars["name"] if options.name == '' else options.name
            outputpath = os.path.abspath(pars["outputpath"]  if options.outputpath == '' else options.outputpath)
            template   = open(os.path.abspath(pars["template"] if options.template == '' else options.template), 'r').read()

            # Make output directory
            if (not os.path.exists(outputpath)):
                  os.makedirs(outputpath)
            
            # Set parameters not in the config file
            formatPars(pars)

            # Replace config file and save
            #pars["name"] = basename+"_"+str(i)
            makeXMLFile(pars, template, outputpath=outputpath)

            # Output scripts
            scriptfile = open(outputpath+"/"+basename+".sh",'a')
            eulerfile  = open(outputpath+"/"+basename+".euler.sh",'a')

            scriptfile.write("# %s \n" % filename)
            eulerfile.write("# %s \n" % filename)

            for seed in seeds:
                  cmd ="java -jar -Xms2G -Xmx4G $1 -overwrite -seed %d -threads %d %s" % (seed, threads, pars["name"]+".xml") 
                  scriptfile.write("nohup %s > %s_%d.out &\n" % (cmd, pars["name"], seed))
                  eulerfile.write("bsub -W%s:0 -n %d -o %s_%d.euler.out -R 'rusage[mem=4096]' %s\n" % (queue, threads, pars["name"], seed, cmd))

            scriptfile.write("\n")
            eulerfile.write("\n")
            scriptfile.close()
            eulerfile.close()
      #
#
