import sys, math
from numpy import *

def ExtractWaveOS_output(infilename):
	print ('Open Gussian output')
	ifile = open(infilename,'r') #open file for reading
	lines = ifile.readlines()
	ifile.close()
	print ('Close Gussian output')

	WaveLength = []
	V_OS = []

	for line in lines:
		if line.find("Excited State  ") >=0:
			line_StateInfo = line.split()
	#		print (line_StateInfo[6])
			WaveLength.append(float(line_StateInfo[6]))
			OS_info = line_StateInfo[8].split('=')
	#		print(OS_info[1])
			V_OS.append(float(OS_info[1]))

	print (WaveLength)

	return WaveLength, V_OS

def GauFilter(WaveLength, V_OS, x):

	N_data = len(WaveLength)

	g = 0.0
	sigma = 10.0

	for i in range(N_data):
		g = g + V_OS[i]*math.exp(-(x-WaveLength[i])**2/(2*sigma**2))/(math.sqrt(2*math.pi)*sigma)

	return g

def GetSpectrum(WaveLength, V_OS):

        interval = 1.00
        Init_value = 200.0
        Final_value = 1000.0

        Integral = 0.

        Index_wave = []
        Intensity = []
        Integral_value = []

        Num_Bin = int((Final_value-Init_value)/interval)

        for i in range(Num_Bin):
                Index_wave.append(Init_value+i*interval)
                Intensity.append(GauFilter(WaveLength, V_OS, Init_value+i*interval))
                Integral = Integral+Intensity[i]*interval

        for i in range(Num_Bin):
                Intensity[i] = Intensity[i]/Integral

        Integral = 0.
        for i in range(Num_Bin):
                Integral = Integral+Intensity[i]
                Integral_value.append(Integral)

        return Index_wave, Intensity, Integral_value



#usage ='Usage; %s infile' % sys.argv[0]
#
#try:
#   infilename = sys.argv[1]
#except:
#   print (usage); sys.exit()

#WaveLength, V_OS = ExtractWaveOS_output(infilename)

#Index_Wave, Intensity = GetSpectrum(WaveLength, V_OS)

#outfilename = 'WaveLength_f.txt'
#ofile = open(outfilename,'w')
#
#for j in range(len(Index_Wave)):
#	ofile.write('% 10.5f  % 10.5f \n' % (Index_Wave[j], Intensity[j]))
#
#ofile.close()



	
