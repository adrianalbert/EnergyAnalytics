# Python script to download and extract data from the CalISO website.
# Developed using Python 2.6
#
# It templates a get request that emulates the submission of their
# wacky form data.
# See "Prices" tab of http://oasis.caiso.com/
# You will probably need to use IE to see the correct interface...
#
# Use and modify as you wish, but please include this full header 
# in your version
#
# written by Sam Borgeson (sborgeson@berkeley.edu)
#
# 8/23/2011
#   Initial working version with zipFile extraction info 
#   and procedural/automated file naming.
#
import os, zipfile
import urllib, urllib2 
from StringIO import StringIO

def retrieveISOData(start,end):
  fn = getFileName(start,end)
  if(os.path.exists(fn)): 
    print '%s alred on disk. Skipping download.' % (fn)
    return
  unzipData(retrieveDailyZip(start,end),fn)

def getLoad(start,end):
  tacArea = 'CA ISO-TAC'
  label = 'Total Actual Hourly Integrated Load'
  day = datetime.timedelta(1)
  curDate = start
  import csv
  fn = getFileName(start,end)
  dataRows = []
  print 'Loading data from %s' % fn
  with open(fn, 'r') as f:
    reader = csv.reader(f)
    try:
      head = reader.next() + ['date']
      if(head[0].startswith('<?xml')):
        print('XML error message. Skipping.')
        return ([],[])
      #print head
      labelIdx   = head.index('LABEL')
      tacAreaIdx = head.index('TAC_AREA_NAME')
      for row in reader:
        if(row[tacAreaIdx] == tacArea and row[labelIdx] == label):
          dataRows.append(row + [curDate.strftime('%Y-%m-%d')])
          curDate = curDate + day
    except csv.Error, e:
        sys.exit('file %s, line %d: %s' % (fn, reader.line_num, e))
  return (head,dataRows)

def retrieveDailyZip(startDate="20110802",endDate='20110802'):
	# dates are in the format YYYYmmdd, as in 20110803
  urlTemplate = 'http://oasis.caiso.com/mrtu-oasis/SingleZip?resultformat=6&queryname=SLD_FCST&startdate=%s&enddate=%s'
  url = urlTemplate % (startDate.strftime('%Y%m%d'), endDate.strftime('%Y%m%d'))
  print url
  urllib.urlretrieve(url,"temp.zip")
  return open("temp.zip", "rb")

def getFileName(start,end):
  return "%s_%s_ISO.csv" % (start.strftime('%Y%m%d'),end.strftime('%Y%m%d'))
  
def unzipData(inFile, outFileName, fileName=None):
  outFile = None
  try:
    #print zipStream.read()
    zipData = zipfile.ZipFile(inFile) # StringIO(zipStream))
    if fileName is None:
      fileName = zipData.namelist()[0]
    #for name in myzipfile.namelist():
    print "Writing %s to %s" % (fileName, outFileName)
    # write out the file with the right name...
    outFile = open(outFileName,"wb")
    outFile.write(zipData.read(fileName))
    #zipData.extract(fileName,outFileName)
  finally:
    try:
      outFile.close()
    except: pass
    try:
      zipData.close()
    except: pass
    try:
      nm = inFile.name
      inFile.close()
      os.remove(nm)
    except: pass

if __name__ == "__main__":
  import csv
  # usage: retrieveISOData will download and write out data from the 
  # ISO using the data and hour supplied (1-25 inclusive) to name the 
  # resulting CSV file
  # ex: for Aug. 2, 2011 hr 2, use:
  #retrieveISOData("20110802",2)

  # OR for a scripted batch over a range of datas and hours:
  #dates = ["20110702","20110703","20110704","20110705"]
  #hours = [1,2,3,4,5,6,7]
  #for d in dates:
    #for h in hours:
    #retrieveISOData(d,h)
  
  # or run through a full year of dates...
  import datetime
  import math
  start = datetime.datetime.strptime('2008-12-22','%Y-%m-%d') # first reliable date in CAISO db
  end   = datetime.datetime.strptime('2013-06-05','%Y-%m-%d')
  dt    = datetime.timedelta(days=20)
  dateRange = end - start
  dlDates = [start] + [start + dt * x for x in range(1,int(math.ceil(dateRange.days / float(dt.days))))] + [end]
  for i in range(len(dlDates)-1): # download the data
    s = dlDates[i]
    e = dlDates[i+1]-datetime.timedelta(1)
    retrieveISOData(s, e) # writes the data file to disk
  with open('caiso.csv', 'wb') as f: # parse the data
    writer = csv.writer(f)
    writeHead = True
    for i in range(len(dlDates)-1):
      s = dlDates[i]
      e = dlDates[i+1]-datetime.timedelta(1)
      (head,readings) = getLoad(s,e)
      if writeHead: 
        writer.writerow(head)
        writeHead = False
      writer.writerows(readings)
  
    
  '''td = datetime.date.today()
  
  start = td - dt * n
  day = start
  data = []
  for i in range(0,n):
    dayStr = day.strftime("%Y%m%d")
    retrieveISOData(dayStr)     # writes the data file to disk
    readings = getLoad(dayStr)  # opens the data file to extract its data
    for (i,val) in enumerate(readings):
      hr = datetime.datetime.combine( day, datetime.time(i, 0) )
      data.append([hr.strftime("%m/%d/%Y %H:%M"),val])
      #print "%s, %s" % (hr.strftime("%m/%d/%Y %H:%M"),val)
    day = day + dt
  with open('CAISO.csv', 'wb') as f:
    writer = csv.writer(f)
    writer.writerows(data)'''