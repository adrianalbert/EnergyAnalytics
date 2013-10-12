import urllib
import sys
import time
import re
from itertools import groupby

#Global Constants
month_lengths = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
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
def stations(zipcode):
	# Read webpage to find station name for a given zipcode
	zipcode_url = "http://www.wunderground.com/cgi-bin/findweather/getForecast?query=%s" % (zipcode)
	nAttempts = 0
	zip_ok = False
	while nAttempts < ATTEMPTS:
		try:
			webfile = urllib.urlopen(zipcode_url)
			data = webfile.read()	
			zip_ok = True					
			break
		except IOError:
			print('Error - retrying...' + str(nAttempts) + '/' + str(ATTEMPTS))
			time.sleep(WAIT_ERR_ZIP)
			continue            
			
	# Hack to find current station
	startIndex = data.find("pwsid") + 7
	endIndex   = data.find('"', startIndex, startIndex + 20)
	station_id = data[startIndex:endIndex]
	
	# Hack to find nearby stations
	pattern1 = re.compile(r'class="pwsrt"\spwsid="^[A-Za-z0-9_-]*"') 
	pattern2 = re.compile(r'class="pwsrt"\spwsid="([^"]*)"')	
	idxRm    = data.find("Station Location")
	dataTrnc = data[idxRm:len(data)]
	nearby   = [station_id]
	for match in pattern2.finditer(dataTrnc):
		nearby.append(match.group(1))
		
	# return stations
	nearby = [ key for key,_ in groupby(nearby)]
	return nearby
	

#my_stations = stations('98720')

#http://www.wunderground.com/weatherstation/WXDailyHistory.asp ID=KCAPORTO6&day=23&month=07&year=2012
#http://www.wunderground.com/weatherstation/WXDailyHistory.asp?ID=KCAPORTO6&day=23&month=10&year=2010&format=1

# generates a list of dates of interest
def generate_dates(day):
	if (len(day) == 0):
		dates = []
		for month in range(3,11):
			for day in range(1,(month_lengths[month-1]+1)):
				dates.append([month,day])
	else:
		year,month,day = day.split('-')
		dates = [[month,day]]
						
	return dates

# dates = generate_dates('2010-04-05')
# dates_all = generate_dates('')

# for given list of stations and dates, send requests for data iteratively
def request_data(nearby_stations, dates):
	cur_data = ''
	for date in dates:			# loop over dates
		month,day = date
		date_url = "&day=" + str(day) + "&month=" + str(month) + "&year=2010&format=1"  
		station_ok = ''
		print(str(month) + ' - ' + str(day))
		for station_id in nearby_stations:
			station_url = "http://www.wunderground.com/weatherstation/WXDailyHistory.asp?ID=%s&" % (station_id)		
			url = station_url + date_url			          
			print('  -->'+url)
			nAttempts = 0	
			data_ok = False		
			while nAttempts < ATTEMPTS:
				while True:
					try:
						webfile = urllib.urlopen(url)
						data = webfile.read()						
						break
					except IOError:
						print('Error - retrying...')
						time.sleep(WAIT_ERR_DAY)
						continue    
				        
				if (len(data)>220):
					data_ok = True
					break
				else:
					print("     ~~ ! Null reading, trying again " + str(nAttempts+1) + "/" + str(ATTEMPTS) )
					nAttempts = nAttempts + 1
					time.sleep(0.1)				
					
			if (data_ok):
				station_ok = station_id
				break
			else:
				print('!!! Error at station  ' + station_id)
				print('    Trying nearby station...')
			
		# process data for dumping to file
		print(' OK: successfully retrieved data from ' + station_ok )
		data = data.replace("\nTime", "Time")
		data = data.replace("<br>", "")
		data = data.replace(",\n", "")
		data = data[0:len(data)-1]
		if (len(cur_data)>0):
			start= data.find('Time')
			stop = data.find('DateUTC') + 7
			data = data.replace(data[start:stop],'')
		
		# append to data	
		cur_data = cur_data + data
		time.sleep(WAIT_OK_DAY)
		
	return cur_data
	

# test data request function		
# tmp = request_data(my_stations, dates_all)

# Main loop[
k = 0
N = len(input)
zipcode_old = ''
data_tot = ''
for line in input:
	
	# treat input
	k = k + 1
	elems = line.split(",")	
	if (len(elems) == 2):
		zipcode,day = elems
		print(str(k) + '/' + str(N) + ': ' + zipcode + '; ' + day)
	else:
		zipcode = elems[0]
		day     = ''
		print(str(k) + '/' + str(N) + ': ' + zipcode)
	
	# generate dates of interest
	dates_list = generate_dates(day)
	
	# get nearby stations
	nearby_stations = stations(zipcode)
	#	print(nearby_stations)
	
	# send data requests 
	data_cur = request_data(nearby_stations, dates_list)
	if (len(data_cur) == 0):
		print(str(zipcode) + ': No data available!!!')
		continue
	
	# store data to file
	if (day == ''):  	# if pulling all data for a zipcode store to own file
		myfile = open(output_dir + "%s.csv" % (zipcode), "w")
		myfile.write(data_cur)
		myfile.close()	
		time.sleep(WAIT_OK_DAY)
	else:				# if working on a list of days store to a single file per zipcode
		if (zipcode == zipcode_old): # if it's same zipcode, append to existing data
			# first remove first line of header
			if (len(data_tot)>0):
				idxFirstLine = data_cur.find('DateUTC\n') + 7
				data_cur     = data_cur[idxFirstLine:len(data_cur)]
			
			data_tot = data_tot + data_cur
		else:						 # if different zipcode, dump previous data
			if (len(data_tot)>0):
				cur_file = "%s%s.csv" %(output_dir,zipcode_old)
				myfile = open(cur_file, "w")
				myfile.write(data_tot)
				myfile.close()	
				print('*** Dumped zipcode %s data to file %s'%(zipcode_old, cur_file))

			data_tot = data_cur
			zipcode_old = zipcode
			print('   New zipcode: ' + str(zipcode))

# Script END
