import urllib
import sys
import time
import re
from itertools import groupby
from datetime import date, timedelta as td

#Global Constants
ATTEMPTS     = 5
WAIT_ERR_ZIP = 0.5
WAIT_ERR_DAY = 0.5
WAIT_OK_DAY  = 0.5

# Open input file
if(len(sys.argv) > 1):
    input = open(sys.argv[1], 'r').readlines()
else:
    input = open("input", 'r').readlines()
input = [x.strip() for x in input]
input = input[1:len(input)]

output_dir = ''
if (len(sys.argv)>2):
	output_dir = sys.argv[2]

print('input:  ' + sys.argv[1])
print('output: ' + output_dir)

# function to get nearby station names given a zipcode
def latlong(station):
	# Read webpage to find station name for a given zipcode
	station_url = "http://www.wunderground.com/weatherstation/WXDailyHistory.asp?ID=%s" % (station)
	nAttempts = 0
	zip_ok = False
	while nAttempts < ATTEMPTS:
		try:
			webfile = urllib.urlopen(station_url)
			data = webfile.read()	
			zip_ok = True					
			break
		except IOError:
			print('Error - retrying...' + str(nAttempts) + '/' + str(ATTEMPTS))
			time.sleep(WAIT_ERR_ZIP)
			continue            
			
	# Hack to find lat/long
	startIndex = data.find("Lat:")
	startIndex = data.find('(', startIndex, startIndex + 50) + 2
	endIndex   = data.find('&', startIndex, startIndex + 100)-1
	lat = float(data[startIndex:endIndex])
	startIndex = data.find("Lon:")
	startIndex = data.find('(', startIndex, startIndex + 50) + 2
	endIndex   = data.find('&', startIndex, startIndex + 100)-1
	lon = float(data[startIndex:endIndex])
	
	info = [lat,lon]
	return info
		
# test data request function		
# tmp = request_data(my_stations, dates_all)

# Main loop[
k = 0
N = len(input)
zipcode_old = ''
data_tot = ''
print("Station,Lat,Long")
for line in input:
	
	k = k + 1
	station = line
	info = latlong(station)
	
	print( "%s,%f,%f"%(station,info[0],info[1]) )
	
	
# Script END
