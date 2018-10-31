import sys, datetime

########################################################
# Some functions for working with and converting dates #
########################################################


def str2int(s):
    try: 
      i = int(s)  
    except:
      i = -1
    return(i)
#


def getDoubleDate(date): 
    """Given datetime object return the decimal date for it

    """

    dec31   = datetime.datetime.strptime(str(date.year)+"-12-31", "%Y-%m-%d")

    datett  = date.timetuple()
    dec31tt = dec31.timetuple()

    return (date.year + float(datett.tm_yday-1)/dec31tt.tm_yday)
##


def getDateString(year, month=None, day=None, calendar=True, digits=None):
    """Return the date as a string

    Examples:

        Calendar dates:
          getDate(year="2015", month="1", day="15", calendar=True) -> "2015-01-15"
          getDate(year="2015", month="1", calendar=True) -> "2015-01"
          getDate(year="2015", month="NA", calendar=True) -> "2015"
          getDate(year="NA", month="NA", calendar=True) -> "NA"

        Decimal dates:
          getDate(year="2015", month="1", day="15", calendar=False) -> "2015.0383561643835"
          getDate(year="2015", month="1", calendar=False) -> "2015.0410958904108"
          getDate(year="2015", month="NA", calendar=False) -> "2015.5"
          getDate(year="NA", month="NA", calendar=False) -> "NA"

    """

    year  = str2int(year)
    month = str2int(month)
    day   = str2int(day)

    if (year < 0):
        # Nothing defined
        date = "NA"
    elif (month < 0):
        # Only year defined
        if (calendar):
            date = "%04d" % year
        else:
            date = "%04d.5" % year

    else:
        if (day < 0):
          # Year and month defined
    
          if (calendar):
              date = "%04d-%02d" % (year, month)
          else:
              monthstart = datetime.date(year, month, 1)
              
              if (month == 12):
                  day = 31
              else:
                  day = (datetime.date(year, month+1, 1)-monthstart).days
              monthend   = datetime.date(year,month,day)

              ddate = (getDoubleDate(monthstart)+getDoubleDate(monthend))/2.0
              if (digits == None): 
                  date = "%f" % ddate
              else:
                  date = ("%." + "%df" % digits) % ddate

        else: 
          # All defined

          if (calendar):
              date = "%04d-%02d-%02d" % (year, month, day)          
          else:
              ddate = getDoubleDate(datetime.date(year,month,day))
              if (digits == None): 
                  date = "%f" % ddate
              else:
                  date = ("%." + "%df" % digits) % ddate              

    return(date)
#




################################################################################################################################  
# Some tests

if __name__ == "__main__":

  years  = ["2015", "2015", "2000", "2015", "2016", "2016"]
  months = ["01",   "12",   "05",   "02",   "02",   "NA"]
  days   = ["01",   "31",   "15",   "NA",   "NA",   "NA"]

  for i in range(0,len(years)):

      sys.stdout.write("%s-%s-%s = " % (years[i], months[i], days[i]))
      sys.stdout.write("%-10s =" % getDateString(years[i],months[i],days[i],calendar=True))
      sys.stdout.write("%19s =" % getDateString(years[i],months[i],days[i],calendar=False))
      sys.stdout.write("%9s\n" % getDateString(years[i],months[i],days[i],calendar=False, digits=3))










