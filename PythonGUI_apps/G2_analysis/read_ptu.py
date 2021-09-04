#MODIFIED FROM 
#https://github.com/PicoQuant/PicoQuant-Time-Tagged-File-Format-Demos/blob/master/PTU/Python/Read_PTU.py


# Read_PTU.py    Read PicoQuant Unified Histogram Files
# This is demo code. Use at your own risk. No warranties.
# Keno Goertz, PicoQUant GmbH, February 2018

# Note that marker events have a lower time resolution and may therefore appear 
# in the file slightly out of order with respect to regular (photon) event records.
# This is by design. Markers are designed only for relatively coarse 
# synchronization requirements such as image scanning. 

# T Mode data are written to an output file [filename]
# We do not keep it in memory because of the huge amout of memory
# this would take in case of large files. Of course you can change this, 
# e.g. if your files are not too big. 
# Otherwise it is best process the data on the fly and keep only the results.

import time
import sys
import struct
import io
import pandas as pd
import numpy as np
import csv
import matplotlib.pyplot as plt
from helper_funcs import refresh_text_box

# Tag Types
tyEmpty8      = struct.unpack(">i", bytes.fromhex("FFFF0008"))[0]
tyBool8       = struct.unpack(">i", bytes.fromhex("00000008"))[0]
tyInt8        = struct.unpack(">i", bytes.fromhex("10000008"))[0]
tyBitSet64    = struct.unpack(">i", bytes.fromhex("11000008"))[0]
tyColor8      = struct.unpack(">i", bytes.fromhex("12000008"))[0]
tyFloat8      = struct.unpack(">i", bytes.fromhex("20000008"))[0]
tyTDateTime   = struct.unpack(">i", bytes.fromhex("21000008"))[0]
tyFloat8Array = struct.unpack(">i", bytes.fromhex("2001FFFF"))[0]
tyAnsiString  = struct.unpack(">i", bytes.fromhex("4001FFFF"))[0]
tyWideString  = struct.unpack(">i", bytes.fromhex("4002FFFF"))[0]
tyBinaryBlob  = struct.unpack(">i", bytes.fromhex("FFFFFFFF"))[0]

# Record types
rtPicoHarpT3     = struct.unpack(">i", bytes.fromhex('00010303'))[0]
rtPicoHarpT2     = struct.unpack(">i", bytes.fromhex('00010203'))[0]
rtHydraHarpT3    = struct.unpack(">i", bytes.fromhex('00010304'))[0]
rtHydraHarpT2    = struct.unpack(">i", bytes.fromhex('00010204'))[0]
rtHydraHarp2T3   = struct.unpack(">i", bytes.fromhex('01010304'))[0]
rtHydraHarp2T2   = struct.unpack(">i", bytes.fromhex('01010204'))[0]
rtTimeHarp260NT3 = struct.unpack(">i", bytes.fromhex('00010305'))[0]
rtTimeHarp260NT2 = struct.unpack(">i", bytes.fromhex('00010205'))[0]
rtTimeHarp260PT3 = struct.unpack(">i", bytes.fromhex('00010306'))[0]
rtTimeHarp260PT2 = struct.unpack(">i", bytes.fromhex('00010206'))[0]
rtMultiHarpT3    = struct.unpack(">i", bytes.fromhex('00010307'))[0]
rtMultiHarpT2    = struct.unpack(">i", bytes.fromhex('00010207'))[0]

# global variables
global ptufile
global inputfile
global outputfile
global recNum
global oflcorrection
global truensync
global dlen
global isT2
global globRes
global numRecords
global tagDataList
global tagNames
global tagValues

def open_ptu():
	"""
	Open ptu file and determine if valid file.
	"""
	global ptufile, inputfile, outputfile, tagNames, tagValues, numRecords, globRes
	inputfile = open(ptufile, "rb")
	# The following is needed for support of wide strings
	outputfile = io.open(ptufile + "_output222.txt", "w+", encoding="utf-16le")

	# Check if inputfile is a valid PTU file
	# Python strings don't have terminating NULL characters, so they're stripped
	magic = inputfile.read(8).decode("utf-8").strip('\0')
	if magic != "PQTTTR":
		print("ERROR: Magic invalid, this is not a PTU file.")
		inputfile.close()
		outputfile.close()
		exit(0)

def read_settings():
	"""
	Read settings in header.
	"""
	global inputfile, outputfile, tagDataList, tagNames, tagValues, numRecords, globRes
	version = inputfile.read(8).decode("utf-8").strip('\0')
	outputfile.write("Tag version: %s\n" % version)

	# Write the header data to outputfile and also save it in memory.
	# There's no do ... while in Python, so an if statement inside the while loop
	# breaks out of it
	tagDataList = []    # Contains tuples of (tagName, tagValue)
	while True:
		tagIdent = inputfile.read(32).decode("utf-8").strip('\0')
		tagIdx = struct.unpack("<i", inputfile.read(4))[0]
		tagTyp = struct.unpack("<i", inputfile.read(4))[0]
		if tagIdx > -1:
			evalName = tagIdent + '(' + str(tagIdx) + ')'
		else:
			evalName = tagIdent
		outputfile.write("\n%-40s" % evalName)
		if tagTyp == tyEmpty8:
			inputfile.read(8)
			outputfile.write("<empty Tag>")
			tagDataList.append((evalName, "<empty Tag>"))
		elif tagTyp == tyBool8:
			tagInt = struct.unpack("<q", inputfile.read(8))[0]
			if tagInt == 0:
				outputfile.write("False")
				tagDataList.append((evalName, "False"))
			else:
				outputfile.write("True")
				tagDataList.append((evalName, "True"))
		elif tagTyp == tyInt8:
			tagInt = struct.unpack("<q", inputfile.read(8))[0]
			outputfile.write("%d" % tagInt)
			tagDataList.append((evalName, tagInt))
		elif tagTyp == tyBitSet64:
			tagInt = struct.unpack("<q", inputfile.read(8))[0]
			outputfile.write("{0:#0{1}x}".format(tagInt,18))
			tagDataList.append((evalName, tagInt))
		elif tagTyp == tyColor8:
			tagInt = struct.unpack("<q", inputfile.read(8))[0]
			outputfile.write("{0:#0{1}x}".format(tagInt,18))
			tagDataList.append((evalName, tagInt))
		elif tagTyp == tyFloat8:
			tagFloat = struct.unpack("<d", inputfile.read(8))[0]
			outputfile.write("%-3E" % tagFloat)
			tagDataList.append((evalName, tagFloat))
		elif tagTyp == tyFloat8Array:
			tagInt = struct.unpack("<q", inputfile.read(8))[0]
			outputfile.write("<Float array with %d entries>" % tagInt/8)
			tagDataList.append((evalName, tagInt))
		elif tagTyp == tyTDateTime:
			tagFloat = struct.unpack("<d", inputfile.read(8))[0]
			tagTime = int((tagFloat - 25569) * 86400)
			tagTime = time.gmtime(tagTime)
			outputfile.write(time.strftime("%a %b %d %H:%M:%S %Y", tagTime))
			tagDataList.append((evalName, tagTime))
		elif tagTyp == tyAnsiString:
			tagInt = struct.unpack("<q", inputfile.read(8))[0]
			tagString = inputfile.read(tagInt).decode("utf-8").strip("\0")
			outputfile.write("%s" % tagString)
			tagDataList.append((evalName, tagString))
		elif tagTyp == tyWideString:
			tagInt = struct.unpack("<q", inputfile.read(8))[0]
			tagString = inputfile.read(tagInt).decode("utf-16le", errors="ignore").strip("\0")
			outputfile.write(tagString)
			tagDataList.append((evalName, tagString))
		elif tagTyp == tyBinaryBlob:
			tagInt = struct.unpack("<q", inputfile.read(8))[0]
			outputfile.write("<Binary blob with %d bytes>" % tagInt)
			tagDataList.append((evalName, tagInt))
		else:
			print("ERROR: Unknown tag type")
			exit(0)
		if tagIdent == "Header_End":
			break

	# Reformat the saved data for easier access
	tagNames = [tagDataList[i][0] for i in range(0, len(tagDataList))]
	tagValues = [tagDataList[i][1] for i in range(0, len(tagDataList))]

	# get important variables from headers
	numRecords = tagValues[tagNames.index("TTResult_NumberOfRecords")]
	globRes = tagValues[tagNames.index("MeasDesc_GlobalResolution")]
	
def gotOverflow(count):
	global outputfile, recNum
	outputfile.write("%u OFL -1 0 %2x\n" % (recNum, count))

def gotMarker(timeTag, markers):
	global outputfile, recNum
	outputfile.write("%u MAR %2x %u\n" % (recNum, markers, timeTag))

def gotPhoton(timeTag, channel, dtime):
	global outputfile, isT2, recNum, globRes
	if isT2:
		outputfile.write("%u CHN %1x %u %8.0lf\n" % (recNum, channel, timeTag,\
						 (timeTag * globRes * 1e12)))
	else:
		outputfile.write("%u CHN %1x %u %8.0lf %10u\n" % (recNum, channel,\
						 timeTag, (timeTag * globRes * 1e9), dtime))

def readPT3(textBrowser = None):
	"""
	Read ptu file for T3 measurement.
	"""
	global inputfile, outputfile, recNum, oflcorrection, dlen, numRecords
	T3WRAPAROUND = 65536
	for recNum in range(0, numRecords):
		# The data is stored in 32 bits that need to be divided into smaller
		# groups of bits, with each group of bits representing a different
		# variable. In this case, channel, dtime and nsync. This can easily be
		# achieved by converting the 32 bits to a string, dividing the groups
		# with simple array slicing, and then converting back into the integers.
		try:
			recordData = "{0:0{1}b}".format(struct.unpack("<I", inputfile.read(4))[0], 32)
		except:
			warning_msg = "The file ended earlier than expected, at record %d/%d."\
				  % (recNum, numRecords)
			if textBrowser is not None:
				refresh_text_box(textBrowser, warning_msg)
			exit(0)

		channel = int(recordData[0:4], base=2)
		dtime = int(recordData[4:16], base=2)
		nsync = int(recordData[16:32], base=2)
		if channel == 0xF: # Special record
			if dtime == 0: # Not a marker, so overflow
				gotOverflow(1)
				oflcorrection += T3WRAPAROUND
			else:
				truensync = oflcorrection + nsync
				gotMarker(truensync, dtime)
		else:
			if channel == 0 or channel > 4: # Should not occur
				print("Illegal Channel: #%1d %1u" % (dlen, channel))
				outputfile.write("\nIllegal channel ")
			truensync = oflcorrection + nsync
			gotPhoton(truensync, channel, dtime)
			dlen += 1
		if recNum % 100000 == 0:
			progress_msg = "\rProgress: %.1f%%" % (float(recNum)*100/float(numRecords))
			if textBrowser is not None:
				refresh_text_box(textBrowser, progress_msg)
			sys.stdout.write(progress_msg)
			sys.stdout.flush()

def readPT2(textBrowser = None):
	"""
	Read ptu file for PT2 measurement.
	"""
	global inputfile, outputfile, recNum, oflcorrection, numRecords
	T2WRAPAROUND = 210698240
	for recNum in range(0, numRecords):
		try:
			recordData = "{0:0{1}b}".format(struct.unpack("<I", inputfile.read(4))[0], 32)
		except:
			warning_msg ="The file ended earlier than expected, at record %d/%d.\n"\
				  % (recNum, numRecords)
			if textBrowser is not None:
				refresh_text_box(textBrowser, warning_msg)
			exit(0)

		channel = int(recordData[0:4], base=2)
		time = int(recordData[4:32], base=2)
		if channel == 0xF: # Special record
			# lower 4 bits of time are marker bits
			markers = int(recordData[28:32], base=2)
			if markers == 0: # Not a marker, so overflow
				gotOverflow(1)
				oflcorrection += T2WRAPAROUND
			else:
				# Actually, the lower 4 bits for the time aren't valid because
				# they belong to the marker. But the error caused by them is
				# so small that we can just ignore it.
				truetime = oflcorrection + time
				gotMarker(truetime, markers)
		else:
			if channel > 4: # Should not occur
				print("Illegal Channel: #%1d %1u" % (recNum, channel))
				outputfile.write("\nIllegal channel ")
			truetime = oflcorrection + time
			gotPhoton(truetime, channel, time)
		if recNum % 100000 == 0:
			if textBrowser is not None:
					refresh_text_box(textBrowser, "Progress: %.1f%%" % (float(recNum)*100/float(numRecords)))
				
			progress_msg = "\rProgress: %.1f%%\n" % (float(recNum)*100/float(numRecords))
			sys.stdout.write(progress_msg)
			sys.stdout.flush()
			
def setup_conversion():
	"""
	Run pre-read functions.
	"""
	open_ptu()
	read_settings()
	
def generate_txt(textBrowser = None):
	"""
	Convert ptu to txt.
	"""
	global oflcorrection, dlen, inputfile, outputfile, tagValues, tagNames, numRecords, isT2
	oflcorrection = 0
	dlen = 0
	outputfile.write("\n-----------------------\n")
	recordType = tagValues[tagNames.index("TTResultFormat_TTTRRecType")]
	start_msg = "Writing %d records, this may take a while..." % (numRecords)
	if textBrowser is not None:
		refresh_text_box(textBrowser, start_msg + "\n")
	if recordType == rtPicoHarpT2:
		isT2 = True
		print("PicoHarp T2 data")
		outputfile.write("PicoHarp T2 data\n")
		outputfile.write("\nrecord# type chan   nsync truetime/ps\n")
		readPT2(textBrowser)
		if textBrowser is not None:
			refresh_text_box(textBrowser, start_msg + "\n")
			refresh_text_box(textBrowser, "Read complete.\n")
	elif recordType == rtPicoHarpT3:
		isT2 = False
		print("PicoHarp T3 data")
		outputfile.write("PicoHarp T3 data\n")
		outputfile.write("\nrecord# type chan   nsync truetime/ns dtime\n")
		readPT3(textBrowser)
		if textBrowser is not None:
			refresh_text_box(textBrowser, start_msg + "\n")
			refresh_text_box(textBrowser, "Read complete.\n")
	else:
		print("ERROR: Unknown record type")
		exit(0)
	inputfile.close()
	outputfile.close()
	return outputfile.name