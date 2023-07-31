import tkinter as tk
import os
import sys
from threading import Thread
import serial
import matplotlib.pyplot as plt
from drawnow import *
import math
import cmath
import numpy as np
import Adafruit_GPIO.SPI as SPI
import Adafruit_MCP3008
import RPi.GPIO as GPIO
import time
from time import sleep
import skrf as rf
from skrf.calibration import SOLT
from skrf.calibration import OnePort
import pylab
import scipy
import scipy.signal as sps
import numpy as np
from math import factorial
from pylab import *

port = serial.Serial("/dev/ttyACM0", baudrate = 115200, timeout=3.0)
noOfPoints = 52

master = tk.Tk()
master.geometry("900x525")
master.resizable(False,False)

def AveragingThreadS11():
    Thread(target=AveragingS11, daemon=True).start()

def CorrectedDataThreadS11():
    Thread(target=CorrectDataS11, daemon=True).start()
   
def PlotCalibratedThreadS11():
    Thread(target=PlotCalibratedS11, daemon=True).start()
    
def ResultSmoothingThreadS11():
    Thread(target=ResultSmoothingS11, daemon=True).start()

def AveragingThreadS21():
    Thread(target=AveragingS21, daemon=True).start()

def CorrectedDataThreadS21():
    Thread(target=CorrectDataS21, daemon=True).start()
   
def PlotCalibratedThreadS21():
    Thread(target=PlotCalibratedS21, daemon=True).start()
    
def ResultSmoothingThreadS21():
    Thread(target=ResultSmoothingS21, daemon=True).start()
    
def AveragingThreadS12():
    Thread(target=AveragingS12, daemon=True).start()

def CorrectedDataThreadS12():
    Thread(target=CorrectDataS12, daemon=True).start()
   
def PlotCalibratedThreadS12():
    Thread(target=PlotCalibratedS12, daemon=True).start()
    
def ResultSmoothingThreadS12():
    Thread(target=ResultSmoothingS12, daemon=True).start()
    
def AveragingThreadS22():
    Thread(target=AveragingS22, daemon=True).start()

def CorrectedDataThreadS22():
    Thread(target=CorrectDataS22, daemon=True).start()
   
def PlotCalibratedThreadS22():
    Thread(target=PlotCalibratedS22, daemon=True).start()
    
def ResultSmoothingThreadS22():
    Thread(target=ResultSmoothingS22, daemon=True).start()
    
def close_plot_thread():
    Thread(target=close_plot, daemon=True).start()

##################################S11#################################
###############################Averaging##############################
def AveragingS11():
    S11Display.config(text = "Collecting data...please wait")
    #-------------------------------------------#
    gainC = []
    ##plt.ion()
    cnt=0
    file = 0
    magSum = 0
    #-------------------------------------------#
    ##Hardware SPI configuration:
    SPI_PORT   = 0
    SPI_DEVICE = 0
    mcp = Adafruit_MCP3008.MCP3008(spi=SPI.SpiDev(SPI_PORT, SPI_DEVICE))
    #-------------------------------------------#
    def plotgainC():
        plt.ylim(-40,10)
        plt.xlim(0,3)
        plt.title('S11 (Input Return Loss)')
        plt.grid(True)
        plt.ylabel('Magnitude (dB)')
        plt.xlabel('Frequency (GHz)')
        plt.plot(np.linspace(0.1,2.7,noOfPoints),gainC, 'rx-', label='dB')
        plt.legend(loc='upper right')
        plt.show()
        

    for turn in range(15):
        
        a = str(turn)
        b = ('result')
        c = ('.s1p')
        nameOfTxt = (b+a+c)
    ##    print("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\nOpening", b+a+c)
        textFile = open('S11 52 points/Noise averaging/'+nameOfTxt, "w+")
        initialFreq = 50000000
        textFile = open('S11 52 points/Noise averaging/'+nameOfTxt, "a")
        
        textFile.write("!Output from S11\n")
        textFile.write("!!S1P File: Measurement: S11:\n")
        textFile.write("# Hz S  dB   R 50 \n")

        for i in range(0,noOfPoints): #pre-load dummy data
            gainC.append(0)

        port.write(b"C0E1r1l100.0u2700.0s50.0t50.000[8.00]8.00^1")
        port.write(b"g")     

        for x in range(noOfPoints):
            
            initialFreq = initialFreq + 50000000
            textFile.write(str(initialFreq))
            textFile.write('\t')        
            
            values = [0]*8
            for i in range(8):
                values[i] = mcp.read_adc(i)
                    
            ##values[0] is the Vmag, values[1] is the Vphase, check the MCP3008 pin config
            
            
            #Magnitude----------------------------------#
            #-------------------------------------------#
            vMag = values[0]
            magRatio = ((vMag-558.82)/(17.503)) ##formula for line equation from readings
            r = round(magRatio, 2) ##round up value for ease of reading
            
    ##        print('Magnitude Ratio = ', r, 'dB')
            
            volRatio = round((10**(r/20)),2)
    ##        print('Voltage Ratio = ', volRatio)
            
            #-------------------------------------------# 
            #Phase--------------------------------------#
            #-------------------------------------------#
            vPhase = values[1]
    ##      phaseDiff = 0.1*(810-((vPhase/1024)*1.8))
            phaseDiff = ((vPhase/1024.0)*180)-90 #found this equation via pages.cs.wisc.edu/~timc/e/spna/index.html 
            phaseDiffInRad = (phaseDiff/180)*3.14159265359
            phi = round(phaseDiffInRad, 2)##round value up for ease of reading
            
    ##        print('Phase           = ', phi, 'Radians')
            
            #-------------------------------------------#
            #Converting to complex form-----------------#
            #-------------------------------------------#       
            realForm = volRatio * math.cos(phi)
            imaginaryForm = volRatio * math.sin(phi)
    ##        print('In complex form -> ',round(realForm,2), '+ ', round(imaginaryForm,2), 'j\n')
            
            textFile.write(str(magRatio))
            textFile.write('\t')
            textFile.write(str(phaseDiff)) 
            textFile.write('\n')

            
            gainC.append(r)
            gainC.pop(0)
                    
            time.sleep(0.05)
        
    ##  print(gainC)
    ##    drawnow(plotgainC)
        textFile.close()
    ##    print(str(round(((turn/15)*100),0))+'% complete')
        S11Display.config(text = str(round(((turn/15)*100),0))+'% complete')
        time.sleep(0.5)
        gainC = []

    S11Display.config(text = 'Data collection complete')
    ##print('Data collection complete\n\n\n\n')        
    ##switch off outputs of SynthHD
    port.write(b"C0E0r0C1E0r0")

    avgFile = open('S11 52 points/Noise averaging/averaged.s1p', "w+")
    initialFreq = 50000000
    newFreq = 50000000
    avgFile = open('S11 52 points/Noise averaging/averaged.s1p', "a")
        
    avgFile.write("!Output from S11\n")
    avgFile.write("!!S1P File: Measurement: S11:\n")
    avgFile.write("# Hz S  dB   R 50 \n")

    ##print("Initiating averaging")
        
    for row in range(noOfPoints):
        magSum = 0
        phsSum = 0
        for file in range(15):
            
            a = str(file)
            b = ('result')
            c = ('.s1p')
            nameOfTxt = (b+a+c)
    ##        print(b+a+c)

            currLines = []
            currFile = open('S11 52 points/Noise averaging/' + nameOfTxt, "r")
            currLines = [_.split() for _ in currFile.readlines()]
            
    ##        print(currLines[row+3][1])
            
            mag = float(currLines[row+3][1])
            phs = float(currLines[row+3][2])
            
            magSum = mag + magSum
            phsSum = phs + phsSum
    ##        print(magSum)
            
    ##        time.sleep(0.05)
            
        newFreq = newFreq + initialFreq
        avgFile.write(str(newFreq))
        avgFile.write("\t")    
        avg = 0
        avg = magSum/15
    ##    print(avg)
    ##    print('\n\n')
        avgFile.write(str(avg)+'\t')
        avg = 0
        avg = phsSum/15
    ##    print(avg)
    ##    print('\n\n')
        avgFile.write(str(avg)+'\n')
      
    S11Display.config(text = 'Averaging complete')  
    avgFile.close()

	
###################################S11################################
###############################CorrectData############################	
def CorrectDataS11():
    
    S11Display.config(text = "Processing. Please wait.")
    gainC = []
    plt.ion()
    cnt=0
    
    nameOfTxt = "attenuatorSampleS11.s1p"

    #-------------------------------------------#

    textFile = open("S11 52 points/"+nameOfTxt, "w+")
    textFile.write("!Output from S11\n")
    textFile.write("!!S1P File: Measurement: S11:\n")
    textFile.write("# Hz S  dB   R 50 \n")
    textFile.close()

    #-------------------------------------------#

    ##Hardware SPI configuration:
    SPI_PORT   = 0
    SPI_DEVICE = 0
    mcp = Adafruit_MCP3008.MCP3008(spi=SPI.SpiDev(SPI_PORT, SPI_DEVICE))

    #-------------------------------------------#
    def plotgainC():
        plt.ylim(-40,10)
        plt.xlim(0,3)
        plt.title('S11 (Input Return Loss)')
        plt.grid(True)
        plt.ylabel('Magnitude (dB)')
        plt.xlabel('Frequency (GHz)')
        plt.plot(np.linspace(0.1,2.7,noOfPoints),gainC, 'rx-', label='dB')
        plt.legend(loc='upper right')
        plt.show()


    initialFreq = 0
    # open results from previous data, and create a matrix with them prevLines[x][y]
    prevLines = []
    textFile = open("S11 52 points/Noise averaging/averaged.s1p", "r")
    prevLines = [_.split() for _ in textFile.readlines()]
    

        
    # read the data file in as a list
    fin = open("S11 52 points/"+nameOfTxt, "r" )
    data_list = fin.readlines()
    fin.close()
         
    # remove list items from index 3 to 5 (inclusive)
    del data_list[3:noOfPoints+10]
        
    # write the changed data (list) to a file
    fout = open("S11 52 points/"+nameOfTxt, "w")
    fout.writelines(data_list)
    fout.close()

    textFile = open("S11 52 points/"+nameOfTxt, "a")

    for i in range(0,noOfPoints): #pre-load dummy data
        gainC.append(0)

    port.write(b"C0E1r1l100.0u2700.0s50.0t50.000[8.00]8.00^1")
    port.write(b"g")     

    for x in range(noOfPoints):
        initialFreq = initialFreq + 50000000
        print(initialFreq)
        textFile.write(str(initialFreq))
        textFile.write('\t')        
            
        values = [0]*8
        for i in range(8):
            values[i] = mcp.read_adc(i)
                    
    ##values[0] is the Vmag, values[1] is the Vphase, check the MCP3008 pin config
            
            
        #Magnitude----------------------------------#
        #-------------------------------------------#
        vMag = values[0]
        magRatio = ((vMag-558.82)/(17.503)) ##formula for line equation from readings
            
        tempR = prevLines[x+3][1] #temporary real value from result.s1p
        tempI = prevLines[x+3][2] #temporary imaginary value from result.s1p
        ftempR = float(tempR) #converting to float
        ftempI = float(tempI) #converting to float
    ##        prevMag = math.sqrt((ftempR*ftempR)+(ftempI*ftempI))
            
        print('\n\n\n\n\nMag = ', magRatio, 'dB')
        print('Prev = ', ftempR, 'dB')
            
        if magRatio < 0: #if value is negative
            if ftempR > 0:
                corr = magRatio + ftempR
                print('Magnitude Ratio = ', corr, 'dB')
            if ftempR < 0:
                corr = magRatio - ftempR
                print('Magnitude Ratio = ', corr, 'dB')                
        if magRatio > 0: #if value is positive
            if ftempR > 0:
                corr = magRatio - ftempR
                print('Magnitude Ratio = ', corr, 'dB')
            if ftempR < 0:
                corr = magRatio + ftempR
                print('Magnitude Ratio = ', corr, 'dB')   


        r = round(corr, 2) ##round up value for ease of reading
                   
        #-------------------------------------------# 
        #Phase--------------------------------------#
        #-------------------------------------------#
        vPhase = values[1]
    ##      phaseDiff = 0.1*(810-((vPhase/1024)*1.8))
        phaseDiff = ((vPhase/1024.0)*180)-90 #found this equation via pages.cs.wisc.edu/~timc/e/spna/index.html 
        phaseDiffInRad = (phaseDiff/180)*3.14159265359
        phi = round(phaseDiffInRad, 2)##round value up for ease of reading
            
        print('Phase           = ', phi, 'Radians')
            
        #-------------------------------------------#
        #Converting to complex form-----------------#
        #-------------------------------------------#       
        realForm = r * math.cos(phi)
        imaginaryForm = r * math.sin(phi) 
            
        textFile.write(str(magRatio))
        textFile.write('\t')
        textFile.write(str(phaseDiff)) 
        textFile.write('\n')
           
        print('Magnitude RatioII = ', r, 'dB')
        gainC.append(r)
        gainC.pop(0)
                    
        time.sleep(0.05)
    initialFreq = 0
    drawnow(plotgainC)
        
    textFile.close()
        
    S11Display.config(text = "Plot Displayed")
         
    gainC = []
        
port.write(b"C0E0r0C1E0r0")	

##################################S11#################################
############################PlotCalibrated############################
def PlotCalibratedS11():
    my_ideals = [\
            rf.Network('S11 52 points/IdealCalibVal/IdealShort.s1p'),
            rf.Network('S11 52 points/IdealCalibVal/IdealOpen.s1p'),
            rf.Network('S11 52 points/IdealCalibVal/IdealBBLoad.s1p'),
            ]

    my_measured = [\
            rf.Network('S11 52 points/MeasCalibVal/smoothMeasShort.s1p'),
            rf.Network('S11 52 points/MeasCalibVal/smoothMeasOpen.s1p'),
            rf.Network('S11 52 points/MeasCalibVal/smoothMeasBBLoad.s1p'),
            ]

    cal = rf.OnePort(\
            ideals = my_ideals,
            measured = my_measured,
            )

    cal.run()
    dut = rf.Network('S11 52 points/attenuatorSampleS11.s1p')

    dut_caled = cal.apply_cal(dut)
    dut_caled.name =  dut.name + ' Calibrated'
    dut_caled.write_touchstone()

    nameOfCalSignal =  dut.name + ' Calibrated.s1p' 

    ro_ori = rf.Network('S11 52 points/attenuatorSampleS11.s1p')
    ro_cal = rf.Network(nameOfCalSignal)

    figure()

    ro_ori.plot_s_db(label='Original')
    ro_cal.plot_s_db(label='Calibrated')

    draw();show();
    
    nameOfCalTxt = "attenuatorSampleS11 Calibrated.s1p"
    nameOfCaldBTxt = "Calibrated&IndBFormatS11.s1p"
    noOfPoints = 52

    realList = []
    imagList = []
    freqList = []

    #--------------------------------------------------#
    textFile = open(nameOfCalTxt, "r")
    original = [_.split() for _ in textFile.readlines()]

    for a in range(noOfPoints):
        real = original[a+3][1]
        realList.append(real) #this list holds all noise filled original values for Magnitude
        
    for b in range(noOfPoints):
        imag = original[b+3][2]
        imagList.append(imag) #this list holds all noise filled original values for Phase
        
    for c in range(noOfPoints):
        freq = original[c+3][0]
        freqList.append(freq)#this list holds all frequenciesfor a in (noOfPoints):

    rA = np.array(realList)
    iA = np.array(imagList)
    fA = np.array(freqList)

    realArray = rA.astype(np.float)
    imagArray = iA.astype(np.float)
    freqArray = fA.astype(np.float)
    #--------------------------------------------------#

    textFile = open(nameOfCaldBTxt, "w+")
    textFile.write("!Output from S11\n")
    textFile.write("!!S1P File: Measurement: S11:\n")
    textFile.write("# Hz S  dB   R 50 \n")

    for d in range(noOfPoints):
    ##    print(freqArray[d])
        rl = realArray[d]
        im = imagArray[d]
        mag = 20 * math.log10(math.sqrt((rl*rl)+(im*im)))
        phs = math.degrees(math.atan(im/rl))
        
        magS = str(mag)
        phsS = str(phs)
        
        textFile.write(str(freqArray[d]))
        textFile.write("\t")    
        textFile.write(magS)
        textFile.write("\t")
        textFile.write(phsS)
        textFile.write("\n")
        
    textFile.close()

##################################S11#################################
############################ResultSmoothing###########################
def ResultSmoothingS11():
    nameOfTxt = "smoothCalibrated&IndBFormatS11.s1p"
    noOfPoints = 52 #there is actually 52 points read

    values = []
    magList = []
    phsList = []
    freqList = []
    valm = []
    valp = []
    freq = []
    smoothed = []
    textFile = open("Calibrated&IndBFormatS11.s1p", "r")
    values = [_.split() for _ in textFile.readlines()] # open results from previous data, and create a matrix with them prevLines[x][y]
    a = 0
    b = 0
    c = 0
    d = 0

    for a in range(noOfPoints):
        valM = values[a+3][1]
        print(values[a+3][1])
        magList.append(valM) #this list holds all noise filled original values for Magnitude
        
    for b in range(noOfPoints):
        valP = values[b+3][2]
        phsList.append(valP) #this list holds all noise filled original values for Phase
        
    for c in range(noOfPoints):
        freq = values[c+3][0]
        freqList.append(freq)#this list holds all frequencies
        
    #converting all three into numpy arrays   
    magArray = np.asarray(magList)
    phsArray = np.asarray(phsList)
    freqArray = np.asarray(freqList)

    ##print(magArray)
    ##print(phsArray)
    ##print(freqArray)

    smoothedMag = sps.savgol_filter(magArray,21,8)
    ##smoothedPhs = sps.savgol_filter(phsArray,3,1)

    textFile = open(nameOfTxt, "w+")
    textFile.write("!Output from S11\n")
    textFile.write("!!S1P File: Measurement: S11:\n")
    textFile.write("# Hz S  dB   R 50 \n")

    ##print(freqArray[5])

    for d in range(noOfPoints):
        textFile.write(freqArray[d])
        textFile.write("\t")
        textFile.write(str(smoothedMag[d]))
        textFile.write("\t")
        textFile.write(str(phsArray[d]))
        textFile.write("\n")

    textFile.close()

    plt.ylim(-40,10)
    plt.xlim(0,3)
    plt.title('S11 (Input Return Loss)')
    plt.grid(True)
    plt.ylabel('Magnitude (dB)')
    plt.xlabel('Frequency (GHz)')
    plt.plot(np.linspace(0.1,2.7,noOfPoints),gainC, 'rx-', label='dB')
    plt.legend(loc='upper right')
    plt.show()
    
    

##################################S21#################################
###############################Averaging##############################
def AveragingS21():

    S21Display.config(text = "Collecting data...please wait")
    #-------------------------------------------#
    gainC = []
    plt.ion()
    cnt=0
    file = 0
    magSum = 0
    #-------------------------------------------#
    ##Hardware SPI configuration:
    SPI_PORT   = 0
    SPI_DEVICE = 0
    mcp = Adafruit_MCP3008.MCP3008(spi=SPI.SpiDev(SPI_PORT, SPI_DEVICE))
    #-------------------------------------------#
    def plotgainC():
        plt.ylim(-40,10)
        plt.xlim(0,3)
        plt.title('S21 (Forward Gain/Isolation)')
        plt.grid(True)
        plt.ylabel('Magnitude (dB)')
        plt.xlabel('Frequency (GHz)')
        pltplot(np.linspace(0.1,2.7,noOfPoints),gainC, 'rx-', label='dB')
        plt.legend(loc='upper right')
        plt.show()

    for turn in range(15):
        
        a = str(turn)
        b = ('result')
        c = ('.s2p')
        nameOfTxt = (b+a+c)
        textFile = open('S21 52 points/Noise averaging/'+nameOfTxt, "w+")
        initialFreq = 50000000
        textFile = open('S21 52 points/Noise averaging/'+nameOfTxt, "a")
        
        textFile.write("!Output from S21\n")
        textFile.write("!!S1P File: Measurement: S21:\n")
        textFile.write("# Hz S  dB   R 50 \n")

        for i in range(0,noOfPoints): #pre-load dummy data
            gainC.append(0)

        port.write(b"C0E1r1l100.0u2700.0s50.0t50.000[8.00]8.00^1")
        port.write(b"g")     

        for x in range(noOfPoints):
            
            initialFreq = initialFreq + 50000000
            textFile.write(str(initialFreq))
            textFile.write('\t')        
            
            values = [0]*8
            for i in range(8):
                values[i] = mcp.read_adc(i)
            #Magnitude----------------------------------#
            #-------------------------------------------#
            vMag = values[0]
            magRatio = ((vMag-558.82)/(17.503)) ##formula for line equation from readings
            r = round(magRatio, 2) ##round up value for ease of reading
            volRatio = round((10**(r/20)),2)
            #-------------------------------------------# 
            #Phase--------------------------------------#
            #-------------------------------------------#
            vPhase = values[1]
            phaseDiff = ((vPhase/1024.0)*180)-90 #found this equation via pages.cs.wisc.edu/~timc/e/spna/index.html 
            phaseDiffInRad = (phaseDiff/180)*3.14159265359
            phi = round(phaseDiffInRad, 2)##round value up for ease of reading
            #-------------------------------------------#
            #Converting to complex form-----------------#
            #-------------------------------------------#       
            realForm = volRatio * math.cos(phi)
            imaginaryForm = volRatio * math.sin(phi)
            
            textFile.write('0.0')
            textFile.write('\t')
            textFile.write('0.0')
            textFile.write('\t')
            
            textFile.write(str(magRatio))
            textFile.write('\t')
            textFile.write(str(phaseDiff)) 
            textFile.write('\t')
            
            textFile.write('0.0')
            textFile.write('\t')
            textFile.write('0.0')
            textFile.write('\t')
            
            textFile.write('0.0')
            textFile.write('\t')
            textFile.write('0.0')
            textFile.write('\n')

            gainC.append(r)
            gainC.pop(0)
                    
            time.sleep(0.05)

        textFile.close()
        S21Display.config(text = str(round(((turn/15)*100),2))+'% complete')
        time.sleep(0.5)
        gainC = []

    port.write(b"C0E0r0C1E0r0")

    newFile = open('S21 52 points/Noise averaging/averaged.s2p', "w+")
    initialFreq = 50000000
    newFreq = 50000000
    newFile = open('S21 52 points/Noise averaging/averaged.s2p', "a")
        
    newFile.write("!Output from S21\n")
    newFile.write("!!S1P File: Measurement: S21:\n")
    newFile.write("# Hz S  dB   R 50 \n")
        
    for row in range(noOfPoints):
        magSum = 0
        phsSum = 0
        for file in range(15):
            
            a = str(file)
            b = ('result')
            c = ('.s2p')
            nameOfTxt = (b+a+c)

            prevLines = []
            textFile = open('S21 52 points/Noise averaging/'+nameOfTxt, "r")
            prevLines = [_.split() for _ in textFile.readlines()]
            
            mag = float(prevLines[row+3][3])
            phs = float(prevLines[row+3][4])
            
            magSum = mag + magSum
            phsSum = phs + phsSum
            
        newFreq = newFreq + initialFreq
        newFile.write(str(newFreq))
        newFile.write("\t")
        
        newFile.write('0.0')
        newFile.write('\t')
        newFile.write('0.0')
        newFile.write('\t')
        
        avg = 0
        avg = magSum/15
        newFile.write(str(avg)+'\t')
        avg = 0
        avg = phsSum/15
        newFile.write(str(avg)+'\t')
      
        newFile.write('0.0')
        newFile.write('\t')
        newFile.write('0.0')
        newFile.write('\t')
            
        newFile.write('0.0')
        newFile.write('\t')
        newFile.write('0.0')
        newFile.write('\n')  
      
    newFile.close()
    S21Display.config(text = "Averaging complete")


###################################S21################################
###############################CorrectData############################
def CorrectDataS21():

    S21Display.config(text = "Processing. Please wait.")

    gainC = []
    plt.ion()
    cnt=0
    corr = 0

    nameOfTxt = "attenuatorSampleS21.s2p"
    #-------------------------------------------#

    textFile = open("S21 52 points/"+nameOfTxt, "w+")
    textFile.write("!Output from S21\n")
    textFile.write("!!S1P File: Measurement: S21:\n")
    textFile.write("# Hz S  dB   R 50 \n")
    textFile.close()

    #-------------------------------------------#

    ##Hardware SPI configuration:
    SPI_PORT   = 0
    SPI_DEVICE = 0
    mcp = Adafruit_MCP3008.MCP3008(spi=SPI.SpiDev(SPI_PORT, SPI_DEVICE))

    #-------------------------------------------#
    def plotgainC():
        plt.ylim(-40,10)
        plt.xlim(0,3)
        plt.title('S21 (Forward Gain/Isolation)')
        plt.grid(True)
        plt.ylabel('Magnitude (dB)')
        plt.xlabel('Frequency (GHz)')
        plt.plot(np.linspace(0.1,2.7,noOfPoints),gainC, 'rx-', label='dB')
        plt.legend(loc='upper right')
        plt.show()


    initialFreq = 0
    # open results from previous data, and create a matrix with them prevLines[x][y]
    prevLines = []
    textFile = open("S21 52 points/Noise averaging/averaged.s2p", "r")
    prevLines = [_.split() for _ in textFile.readlines()]
            
    # read the data file in as a list
    fin = open("S21 52 points/"+nameOfTxt, "r" )
    data_list = fin.readlines()
    fin.close()
         
    # remove list items from index 3 to 5 (inclusive)
    del data_list[3:noOfPoints+10]
        
    # write the changed data (list) to a file
    fout = open("S21 52 points/"+nameOfTxt, "w")
    fout.writelines(data_list)
    fout.close()

    textFile = open("S21 52 points/"+nameOfTxt, "a")

    for i in range(0,noOfPoints): #pre-load dummy data
        gainC.append(0)

    port.write(b"C0E1r1l100.0u2700.0s50.0t50.000[8.00]8.00^1")
    port.write(b"g")     

    for x in range(noOfPoints):
        initialFreq = initialFreq + 50000000
        textFile.write(str(initialFreq))
        textFile.write('\t')        
            
        values = [0]*8
        for i in range(8):
            values[i] = mcp.read_adc(i)
                    
    ##values[0] is the Vmag, values[1] is the Vphase, check the MCP3008 pin config
            
            
        #Magnitude----------------------------------#
        #-------------------------------------------#
        vMag = values[0]
        magRatio = ((vMag-558.82)/(17.503)) ##formula for line equation from readings
             
        tempR = prevLines[x+3][3] #temporary real value from result.s1p
        tempI = prevLines[x+3][4] #temporary imaginary value from result.s1p
        ftempR = float(tempR) #converting to float
        ftempI = float(tempI) #converting to float
    ##        prevMag = math.sqrt((ftempR*ftempR)+(ftempI*ftempI))

        if magRatio < 0: #if value is negative
            if ftempR > 0:
                corr = magRatio + ftempR
            if ftempR < 0:
                corr = magRatio - ftempR              
        if magRatio > 0: #if value is positive
            if ftempR > 0:
                corr = magRatio - ftempR
            if ftempR < 0:
                corr = magRatio + ftempR
        if magRatio == 0:
            corr = magRatio
        
        print(magRatio)
        print(corr)
        r = round(magRatio, 2) ##round up value for ease of reading
        print(r)
                   
        #-------------------------------------------# 
        #Phase--------------------------------------#
        #-------------------------------------------#
        vPhase = values[1]
    ##  phaseDiff = 0.1*(810-((vPhase/1024)*1.8))
        phaseDiff = ((vPhase/1024.0)*180)-90 #found this equation via pages.cs.wisc.edu/~timc/e/spna/index.html 
        phaseDiffInRad = (phaseDiff/180)*3.14159265359
        phi = round(phaseDiffInRad, 2)##round value up for ease of reading
            
        #-------------------------------------------#
        #Converting to complex form-----------------#
        #-------------------------------------------#       
        realForm = r * math.cos(phi)
        imaginaryForm = r * math.sin(phi)
            
        textFile.write('0.0')
        textFile.write('\t')
        textFile.write('0.0')
        textFile.write('\t')
        
        textFile.write(str(magRatio))
        textFile.write('\t')
        textFile.write(str(phaseDiff)) 
        textFile.write('\t')
        
        textFile.write('0.0')
        textFile.write('\t')
        textFile.write('0.0')
        textFile.write('\t')
            
        textFile.write('0.0')
        textFile.write('\t')
        textFile.write('0.0')
        textFile.write('\n')
            
        gainC.append(r)
        gainC.pop(0)
                    
        time.sleep(0.05)
    initialFreq = 0
    ##  print(gainC)
    drawnow(plotgainC)
        
    textFile.close()
        
    S21Display.config(text = "Plot Displayed")
            
    port.write(b"C0E0r0C1E0r0")
        
    gainC = []


##################################S21#################################
############################PlotCalibrated############################
def PlotCalibratedS21():
    # a list of Network types, holding 'ideal' responses
    my_ideals = [
        rf.Network('S21 52 points/IdealCalibVal/IdealThru.s2p'),
        rf.Network('S21 52 points/IdealCalibVal/IdealShort.s2p'),
        rf.Network('S21 52 points/IdealCalibVal/IdealOpen.s2p'),
        rf.Network('S21 52 points/IdealCalibVal/IdealBBLoad.s2p'),
        ]

    # a list of Network types, holding 'measured' responses
    my_measured = [
        rf.Network('S21 52 points/MeasCalibVal/MeasThru.s2p'),
        rf.Network('S21 52 points/MeasCalibVal/MeasShort.s2p'),
        rf.Network('S21 52 points/MeasCalibVal/MeasOpen.s2p'),
        rf.Network('S21 52 points/MeasCalibVal/MeasBBLoad.s2p'),
        ]

    ## create a SOLT instance
    cal = SOLT(
        ideals = my_ideals,
        measured = my_measured,
        )

    ## run, and apply calibration to a DUT
    # run calibration algorithm
    cal.run()

    # apply it to a dut
    dut = rf.Network('S21 52 points/attenuatorSampleS21.s2p')
    dut_caled = cal.apply_cal(dut)
    dut_caled.name =  'S21 52 points/' + dut.name + ' Calibrated'

    # save results
    dut_caled.write_touchstone()

    nameOfCalTxt = "attenuatorSampleS21 Calibrated.s2p"
    nameOfCaldBTxt = "S21 52 points/Calibrated&IndBFormatS21.s2p"
    noOfPoints = 52

    freqList = []
    lineList = []

    #--------------------------------------------------#
    textFile = open("S21 52 points/"+nameOfCalTxt, "r")
    original = [_.split() for _ in textFile.readlines()]

    for c in range(noOfPoints):
        freq = original[c+3][0]
        freqList.append(freq)#this list holds all frequenciesfor a in (noOfPoints):
        
        freqArray = np.array(freqList)
        
    textFile = open(nameOfCaldBTxt, "w+")
    textFile.write("!Output from S21\n")
    textFile.write("!!S1P File: Measurement: S21:\n")
    textFile.write("# Hz S  dB   R 50 \n")

    for a in range(noOfPoints):
        
        textFile.write(str(freqArray[a]))
        textFile.write("\t")
        
        lineList = []
        lA = []
        lArray = []

        for column in range(8):
            line = original[a+3][column+1]
            lineList.append(line)
            
        lA = np.array(lineList)
            
        lArray = lA.astype(np.float)
            
        rl1 = lArray[0]
        im1 = lArray[1]
        rl2 = lArray[2]
        im2 = lArray[3]
        rl3 = lArray[4]
        im3 = lArray[5]
        rl4 = lArray[6]
        im4 = lArray[7]
        
        mag1 = 20 * math.log10(math.sqrt((rl1*rl1)+(im1*im1)))
        phs1 = math.degrees(math.atan(im1/rl1))
            
        mag2 = 20 * math.log10(math.sqrt((rl2*rl2)+(im2*im2)))
        phs2 = math.degrees(math.atan(im2/rl2))
            
        mag3 = 20 * math.log10(math.sqrt((rl3*rl3)+(im3*im3)))
        phs3 = math.degrees(math.atan(im3/rl3))
        
        mag4 = 20 * math.log10(math.sqrt((rl4*rl4)+(im4*im4)))
        phs4 = math.degrees(math.atan(im4/rl4))
        
        magS1 = str(mag1)
        phsS1 = str(phs1)
        magS2 = str(mag2)
        phsS2 = str(phs2)
        magS3 = str(mag3)
        phsS3 = str(phs3)
        magS4 = str(mag4)
        phsS4 = str(phs4)          
        
        textFile.write(magS1)
        textFile.write("\t")
        textFile.write(phsS1)
        textFile.write("\t")
            
        textFile.write(magS2)
        textFile.write("\t")
        textFile.write(phsS2)
        textFile.write("\t")
            
        textFile.write(magS3)
        textFile.write("\t")
        textFile.write(phsS3)
        textFile.write("\t")
            
        textFile.write(magS4)
        textFile.write("\t")
        textFile.write(phsS4)
        textFile.write("\n")
        
        mag1 = 0
        phs1 = 0
        mag2 = 0
        phs2 = 0
        mag3 = 0
        phs3 = 0
        mag4 = 0
        phs4 = 0
        
    textFile.close()

    indBFormat = rf.Network('S21 52 points/Calibrated&IndBFormatS21.s2p')

    pylab.figure(2)
    pylab.title('indB, Mag')
    ##dut.plot_s_db(m=1,n=0)
    indBFormat.plot_s_db(m=1,n=0) # m,n are S-Matrix indices
    # show the plots
    pylab.ylim(-40,10)
    pylab.show()


##################################S21#################################
############################ResultSmoothing###########################
def ResultSmoothingS21():
    nameOfTxt = "Calibrated&SmoothS21.s2p"
    noOfPoints = 52 #there are actually 52 points to read

    values = []
    magList1 = []
    phsList1 = []
    magList2 = []
    phsList2 = []
    magList3 = []
    phsList3 = []
    magList4 = []
    phsList4 = []
    freqList = []
    valM1 = []
    valP1 = []
    valM2 = []
    valP2 = []
    valM3 = []
    valP3 = []
    valM4 = []
    valP4 = []
    freq = []
    smoothed = []
    ##textFile = open("/home/pi/Desktop/S21 52 points/Calibrated&IndBFormat.s1p", "r")
    textFile = open("S21 52 points/Calibrated&IndBFormatS21.s2p", "r")
    values = [_.split() for _ in textFile.readlines()] # open results from previous data, and create a matrix with them prevLines[x][y]
    a = 0

    for a in range(noOfPoints):
        valM1 = values[a+3][1]
        magList1.append(valM1)
        valP1 = values[a+3][2]
        phsList1.append(valP1)
        
        valM2 = values[a+3][3]
        magList2.append(valM2)
        valP2 = values[a+3][4]
        phsList2.append(valP2)
        
        valM3 = values[a+3][5]
        magList3.append(valM3)
        valP3 = values[a+3][6]
        phsList3.append(valP3)
        
        valM4 = values[a+3][7]
        magList4.append(valM4)
        valP4 = values[a+3][8]
        phsList4.append(valP4)
        
        freq = values[a+3][0]
        freqList.append(freq)#this list holds all frequencies
        
    #converting all three into numpy arrays   
    magArray1 = np.asarray(magList1)
    phsArray1 = np.asarray(phsList1)
    magArray2 = np.asarray(magList2)
    phsArray2 = np.asarray(phsList2)
    magArray3 = np.asarray(magList3)
    phsArray3 = np.asarray(phsList3)
    magArray4 = np.asarray(magList4)
    phsArray4 = np.asarray(phsList4)
    freqArray = np.asarray(freqList)

    smoothedMag1 = sps.savgol_filter(magArray1,21,8)
    smoothedMag2 = sps.savgol_filter(magArray2,21,7)
    smoothedMag3 = sps.savgol_filter(magArray3,21,8)
    smoothedMag4 = sps.savgol_filter(magArray4,21,8)
    ##smoothedPhs = sps.savgol_filter(phsArray,3,1)

    textFile = open("S21 52 points/"+nameOfTxt, "w+")
    textFile.write("!Output from S21\n")
    textFile.write("!!S1P File: Measurement: S21:\n")
    textFile.write("# Hz S  dB   R 50 \n")

    ##print(freqArray[5])

    for a in range(noOfPoints):
        textFile.write(freqArray[a])
        textFile.write("\t")
        textFile.write(str(smoothedMag1[a]))
        textFile.write("\t")
        textFile.write(str(phsArray1[a]))
        textFile.write("\t")
        textFile.write(str(smoothedMag2[a]))
        textFile.write("\t")
        textFile.write(str(phsArray2[a]))
        textFile.write("\t")
        textFile.write(str(smoothedMag3[a]))
        textFile.write("\t")
        textFile.write(str(phsArray3[a]))
        textFile.write("\t")
        textFile.write(str(smoothedMag4[a]))
        textFile.write("\t")
        textFile.write(str(phsArray4[a]))
        textFile.write("\n")

    textFile.close()

    plt.ylim(-40,10)
    plt.xlim(0,3)
    plt.title('S21 (Forward Gain/Isolation)')
    plt.grid(True)
    plt.ylabel('Magnitude (dB)')
    plt.xlabel('Frequency (GHz)')
    plt.legend(loc='upper right')
    ##plt.plot(np.linspace(0.1,2.7,noOfPoints),smoothedMag1, label='dBm', color = 'blue')
    plt.plot(np.linspace(0.1,2.7,noOfPoints),smoothedMag2, label='dBm', color = 'orange')
    ##plt.plot(np.linspace(0.1,2.7,noOfPoints),smoothedMag3, label='dBm', color = 'green')
    ##plt.plot(np.linspace(0.1,2.7,noOfPoints),smoothedMag4, label='dBm', color = 'red')

    plt.show()
    S21Display.config(text = "Smoothing complete")


##################################S12#################################
###############################Averaging##############################
def AveragingS12():

	S12Display.config(text = "Collecting data...please wait")
    #-------------------------------------------#
	gainC = []
	plt.ion()
	cnt=0
	file = 0
	magSum = 0
    #-------------------------------------------#
    ##Hardware SPI configuration:
	SPI_PORT   = 0
	SPI_DEVICE = 0
	mcp = Adafruit_MCP3008.MCP3008(spi=SPI.SpiDev(SPI_PORT, SPI_DEVICE))
    #-------------------------------------------#
	def plotgainC():
            plt.ylim(-40,10)
            plt.xlim(0,3)
            plt.title('S12 (Reverse Gain/Isolation)')
            plt.grid(True)
            plt.ylabel('Magnitude (dB)')
            plt.xlabel('Frequency (GHz)')
            plt.plot(np.linspace(0.1,2.7,noOfPoints),gainC, 'rx-', label='dB')
            plt.legend(loc='upper right')
            plt.show()

	for turn in range(15):
        
		a = str(turn)
		b = ('result')
		c = ('.s2p')
		nameOfTxt = (b+a+c)
		textFile = open('S12 52 points/Noise averaging/'+nameOfTxt, "w+")
		initialFreq = 50000000
		textFile = open('S12 52 points/Noise averaging/'+nameOfTxt, "a")
        
		textFile.write("!Output from S21\n")
		textFile.write("!!S1P File: Measurement: S12:\n")
		textFile.write("# Hz S  dB   R 50 \n")

		for i in range(0,noOfPoints): #pre-load dummy data
			gainC.append(0)

		port.write(b"C1E1r1l100.0u2700.0s50.0t50.000[0.00]0.00^1")
		port.write(b"g")     

		for x in range(noOfPoints):
            
			initialFreq = initialFreq + 50000000
			textFile.write(str(initialFreq))
			textFile.write('\t')        
            
			values = [0]*8
			for i in range(8):
				values[i] = mcp.read_adc(i)
            #Magnitude----------------------------------#
            #-------------------------------------------#
			vMag = values[2]
			magRatio = ((vMag-559.03)/(17.915)) ##formula for line equation from readings
			r = round(magRatio, 2) ##round up value for ease of reading
			volRatio = round((10**(r/20)),2)
            #-------------------------------------------# 
            #Phase--------------------------------------#
            #-------------------------------------------#
			vPhase = values[3]
			phaseDiff = ((vPhase/1024.0)*180)-90 #found this equation via pages.cs.wisc.edu/~timc/e/spna/index.html 
			phaseDiffInRad = (phaseDiff/180)*3.14159265359
			phi = round(phaseDiffInRad, 2)##round value up for ease of reading
            #-------------------------------------------#
            #Converting to complex form-----------------#
            #-------------------------------------------#       
			realForm = volRatio * math.cos(phi)
			imaginaryForm = volRatio * math.sin(phi)

			textFile.write('0.0')
			textFile.write('\t')
			textFile.write('0.0')
			textFile.write('\t')
       
			textFile.write('0.0')
			textFile.write('\t')
			textFile.write('0.0')
			textFile.write('\t')
        
			textFile.write(str(magRatio))
			textFile.write('\t')
			textFile.write(str(phaseDiff)) 
			textFile.write('\t')
        
			textFile.write('0.0')
			textFile.write('\t')
			textFile.write('0.0')
			textFile.write('\n')

			gainC.append(r)
			gainC.pop(0)
                    
			time.sleep(0.05)

		textFile.close()
		S21Display.config(text = str(round(((turn/15)*100),2))+'% complete')
		time.sleep(0.5)
		gainC = []

	port.write(b"C0E0r0C1E0r0")

	newFile = open('S12 52 points/Noise averaging/averaged.s2p', "w+")
	initialFreq = 50000000
	newFreq = 50000000
	newFile = open('S12 52 points/Noise averaging/averaged.s2p', "a")
        
	newFile.write("!Output from S12\n")
	newFile.write("!!S1P File: Measurement: S12:\n")
	newFile.write("# Hz S  dB   R 50 \n")
        
	for row in range(noOfPoints):
		magSum = 0
		phsSum = 0
		for file in range(15):
            
			a = str(file)
			b = ('result')
			c = ('.s2p')
			nameOfTxt = (b+a+c)

			prevLines = []
			textFile = open('S12 52 points/Noise averaging/'+nameOfTxt, "r")
			prevLines = [_.split() for _ in textFile.readlines()]
            
			mag = float(prevLines[row+3][3])
			phs = float(prevLines[row+3][4])
            
			magSum = mag + magSum
			phsSum = phs + phsSum
            
		newFreq = newFreq + initialFreq
		newFile.write(str(newFreq))
		newFile.write("\t")
        
		newFile.write('0.0')
		newFile.write('\t')
		newFile.write('0.0')
		newFile.write('\t')
          
		newFile.write('0.0')
		newFile.write('\t')
		newFile.write('0.0')
		newFile.write('\t')
            
		avg = 0
		avg = magSum/15
		print(avg)
		print('\n\n')
		newFile.write(str(avg)+'\t')
		avg = 0
		avg = phsSum/15
		print(avg)
		print('\n\n')
		newFile.write(str(avg)+'\t')

		newFile.write('0.0')
		newFile.write('\t')
		newFile.write('0.0')
		newFile.write('\n')    
      
	newFile.close()
	S12Display.config(text = "Averaging complete")	
	
###################################S12################################
###############################CorrectData############################
def CorrectDataS12():

	S12Display.config(text = "Processing. Please wait.")

	gainC = []
	plt.ion()
	cnt=0
	corr = 0

	nameOfTxt = "attenuatorSampleS12.s2p"
    #-------------------------------------------#

	textFile = open("S12 52 points/"+nameOfTxt, "w+")
	textFile.write("!Output from S12\n")
	textFile.write("!!S1P File: Measurement: S12:\n")
	textFile.write("# Hz S  dB   R 50 \n")
	textFile.close()

    #-------------------------------------------#

    ##Hardware SPI configuration:
	SPI_PORT   = 0
	SPI_DEVICE = 0
	mcp = Adafruit_MCP3008.MCP3008(spi=SPI.SpiDev(SPI_PORT, SPI_DEVICE))

    #-------------------------------------------#
	def plotgainC():
            plt.ylim(-40,10)
            plt.xlim(0,3)
            plt.title('S12 (Reverse Gain/Isolation)')
            plt.grid(True)
            plt.ylabel('Magnitude (dB)')
            plt.xlabel('Frequency (GHz)')
            plt.plot(np.linspace(0.1,2.7,noOfPoints),gainC, 'rx-', label='dB')
            plt.legend(loc='upper right')
            plt.show()


	initialFreq = 0
    # open results from previous data, and create a matrix with them prevLines[x][y]
	prevLines = []
	textFile = open("S12 52 points/Noise averaging/averaged.s2p", "r")
	prevLines = [_.split() for _ in textFile.readlines()]
            
    # read the data file in as a list
	fin = open("S12 52 points/"+nameOfTxt, "r" )
	data_list = fin.readlines()
	fin.close()
         
    # remove list items from index 3 to 5 (inclusive)
	del data_list[3:noOfPoints+10]
        
    # write the changed data (list) to a file
	fout = open("S12 52 points/"+nameOfTxt, "w")
	fout.writelines(data_list)
	fout.close()

	textFile = open("S12 52 points/"+nameOfTxt, "a")

	for i in range(0,noOfPoints): #pre-load dummy data
		gainC.append(0)

	port.write(b"C1E1r1l100.0u2700.0s50.0t50.000[0.00]0.00^1")
	port.write(b"g")     

	for x in range(noOfPoints):
		initialFreq = initialFreq + 50000000
		textFile.write(str(initialFreq))
		textFile.write('\t')        
            
		values = [0]*8
		for i in range(8):
			values[i] = mcp.read_adc(i)
                    
    ##values[0] is the Vmag, values[1] is the Vphase, check the MCP3008 pin config
            
            
        #Magnitude----------------------------------#
        #-------------------------------------------#
		vMag = values[2]
		magRatio = ((vMag-559.03)/(17.915)) ##formula for line equation from readings
             
		tempR = prevLines[x+3][5] #temporary real value from result.s1p
		tempI = prevLines[x+3][6] #temporary imaginary value from result.s1p
		ftempR = float(tempR) #converting to float
		ftempI = float(tempI) #converting to float
    ##        prevMag = math.sqrt((ftempR*ftempR)+(ftempI*ftempI))

		if magRatio < 0: #if value is negative
			if ftempR > 0:
				corr = magRatio + ftempR
			if ftempR < 0:
				corr = magRatio - ftempR              
		if magRatio > 0: #if value is positive
			if ftempR > 0:
				corr = magRatio - ftempR
			if ftempR < 0:
				corr = magRatio + ftempR
		if magRatio == 0:
			corr = magRatio
        
		print(magRatio)
		print(corr)
		r = round(magRatio, 2) ##round up value for ease of reading
		print(r)
                   
        #-------------------------------------------# 
        #Phase--------------------------------------#
        #-------------------------------------------#
		vPhase = values[3]
    ##  phaseDiff = 0.1*(810-((vPhase/1024)*1.8))
		phaseDiff = ((vPhase/1024.0)*180)-90 #found this equation via pages.cs.wisc.edu/~timc/e/spna/index.html 
		phaseDiffInRad = (phaseDiff/180)*3.14159265359
		phi = round(phaseDiffInRad, 2)##round value up for ease of reading
            
        #-------------------------------------------#
        #Converting to complex form-----------------#
        #-------------------------------------------#       
		realForm = r * math.cos(phi)
		imaginaryForm = r * math.sin(phi)
            
		textFile.write('0.0')
		textFile.write('\t')
		textFile.write('0.0')
		textFile.write('\t')
       
		textFile.write('0.0')
		textFile.write('\t')
		textFile.write('0.0')
		textFile.write('\t')

		textFile.write(str(magRatio))
		textFile.write('\t')
		textFile.write(str(phaseDiff)) 
		textFile.write('\t')
            
		textFile.write('0.0')
		textFile.write('\t')
		textFile.write('0.0')
		textFile.write('\n')
            
		gainC.append(r)
		gainC.pop(0)
                    
		time.sleep(0.05)
	initialFreq = 0
    ##  print(gainC)
	drawnow(plotgainC)
        
	textFile.close()
        
	S12Display.config(text = "Plot Displayed")
            
	port.write(b"C0E0r0C1E0r0")
        
	gainC = []	

##################################S12#################################
############################PlotCalibrated############################
def PlotCalibratedS12():
    # a list of Network types, holding 'ideal' responses
    my_ideals = [
        rf.Network('S12 52 points/IdealCalibVal/IdealThru.s2p'),
        rf.Network('S12 52 points/IdealCalibVal/IdealShort.s2p'),
        rf.Network('S12 52 points/IdealCalibVal/IdealOpen.s2p'),
        rf.Network('S12 52 points/IdealCalibVal/IdealBBLoad.s2p'),
        ]

    # a list of Network types, holding 'measured' responses
    my_measured = [
        rf.Network('S12 52 points/MeasCalibVal/MeasThru.s2p'),
        rf.Network('S12 52 points/MeasCalibVal/MeasShort.s2p'),
        rf.Network('S12 52 points/MeasCalibVal/MeasOpen.s2p'),
        rf.Network('S12 52 points/MeasCalibVal/MeasBBLoad.s2p'),
        ]

    ## create a SOLT instance
    cal = SOLT(
        ideals = my_ideals,
        measured = my_measured,
        )

    ## run, and apply calibration to a DUT
    # run calibration algorithm
    cal.run()

    # apply it to a dut
    dut = rf.Network('S12 52 points/attenuatorSampleS12.s2p')
    dut_caled = cal.apply_cal(dut)
    dut_caled.name =  'S12 52 points/' + dut.name + ' calibrated'

    # save results
    dut_caled.write_touchstone()

    nameOfCalTxt = "attenuatorSampleS12 calibrated.s2p"
    nameOfCaldBTxt = "S12 52 points/Calibrated&IndBFormatS12.s2p"
    noOfPoints = 52

    freqList = []
    lineList = []

    #--------------------------------------------------#
    textFile = open("S12 52 points/"+nameOfCalTxt, "r")
    original = [_.split() for _ in textFile.readlines()]

    for c in range(noOfPoints):
        freq = original[c+3][0]
        freqList.append(freq)#this list holds all frequenciesfor a in (noOfPoints):
        
        freqArray = np.array(freqList)
        
    textFile = open(nameOfCaldBTxt, "w+")
    textFile.write("!Output from S12\n")
    textFile.write("!!S1P File: Measurement: S12:\n")
    textFile.write("# Hz S  dB   R 50 \n")

    for a in range(noOfPoints):
        
        textFile.write(str(freqArray[a]))
        textFile.write("\t")
        
        lineList = []
        lA = []
        lArray = []

        for column in range(8):
            line = original[a+3][column+1]
            lineList.append(line)
            
        lA = np.array(lineList)
            
        lArray = lA.astype(np.float)
            
        rl1 = lArray[0]
        im1 = lArray[1]
        rl2 = lArray[2]
        im2 = lArray[3]
        rl3 = lArray[4]
        im3 = lArray[5]
        rl4 = lArray[6]
        im4 = lArray[7]
        
        mag1 = 20 * math.log10(math.sqrt((rl1*rl1)+(im1*im1)))
        phs1 = math.degrees(math.atan(im1/rl1))
            
        mag2 = 20 * math.log10(math.sqrt((rl2*rl2)+(im2*im2)))
        phs2 = math.degrees(math.atan(im2/rl2))
            
        mag3 = 20 * math.log10(math.sqrt((rl3*rl3)+(im3*im3)))
        phs3 = math.degrees(math.atan(im3/rl3))
        
        mag4 = 20 * math.log10(math.sqrt((rl4*rl4)+(im4*im4)))
        phs4 = math.degrees(math.atan(im4/rl4))
        
        magS1 = str(mag1)
        phsS1 = str(phs1)
        magS2 = str(mag2)
        phsS2 = str(phs2)
        magS3 = str(mag3)
        phsS3 = str(phs3)
        magS4 = str(mag4)
        phsS4 = str(phs4)          
        
        textFile.write(magS1)
        textFile.write("\t")
        textFile.write(phsS1)
        textFile.write("\t")
            
        textFile.write(magS2)
        textFile.write("\t")
        textFile.write(phsS2)
        textFile.write("\t")
            
        textFile.write(magS3)
        textFile.write("\t")
        textFile.write(phsS3)
        textFile.write("\t")
            
        textFile.write(magS4)
        textFile.write("\t")
        textFile.write(phsS4)
        textFile.write("\n")
        
        mag1 = 0
        phs1 = 0
        mag2 = 0
        phs2 = 0
        mag3 = 0
        phs3 = 0
        mag4 = 0
        phs4 = 0
        
    textFile.close()

    indBFormat = rf.Network('S12 52 points/Calibrated&IndBFormatS12.s2p')

    pylab.figure(2)
    pylab.title('indB, Mag')
    ##dut.plot_s_db(m=1,n=0)
    indBFormat.plot_s_db(m=1,n=0) # m,n are S-Matrix indices
    # show the plots
    pylab.ylim(-40,10)
    pylab.show()	
	
##################################S12#################################
############################ResultSmoothing###########################
def ResultSmoothingS12():
    nameOfTxt = "MeasBBLoadSmoothS12.s2p"
    noOfPoints = 52 #there are actually 52 points read

    values = []
    magList1 = []
    phsList1 = []
    magList2 = []
    phsList2 = []
    magList3 = []
    phsList3 = []
    magList4 = []
    phsList4 = []
    freqList = []
    valM1 = []
    valP1 = []
    valM2 = []
    valP2 = []
    valM3 = []
    valP3 = []
    valM4 = []
    valP4 = []
    freq = []
    smoothed = []
    ##textFile = open("/home/pi/Desktop/S12 52 points/Calibrated&IndBFormat.s1p", "r")
    textFile = open("S12 52 points/Calibrated&IndBFormatS12.s2p", "r")
    values = [_.split() for _ in textFile.readlines()] # open results from previous data, and create a matrix with them prevLines[x][y]
    a = 0

    for a in range(noOfPoints):
        valM1 = values[a+3][1]
        magList1.append(valM1)
        valP1 = values[a+3][2]
        phsList1.append(valP1)
        
        valM2 = values[a+3][3]
        magList2.append(valM2)
        valP2 = values[a+3][4]
        phsList2.append(valP2)
        
        valM3 = values[a+3][5]
        magList3.append(valM3)
        valP3 = values[a+3][6]
        phsList3.append(valP3)
        
        valM4 = values[a+3][7]
        magList4.append(valM4)
        valP4 = values[a+3][8]
        phsList4.append(valP4)
        
        freq = values[a+3][0]
        freqList.append(freq)#this list holds all frequencies
        
    #converting all three into numpy arrays   
    magArray1 = np.asarray(magList1)
    phsArray1 = np.asarray(phsList1)
    magArray2 = np.asarray(magList2)
    phsArray2 = np.asarray(phsList2)
    magArray3 = np.asarray(magList3)
    phsArray3 = np.asarray(phsList3)
    magArray4 = np.asarray(magList4)
    phsArray4 = np.asarray(phsList4)
    freqArray = np.asarray(freqList)

    smoothedMag1 = sps.savgol_filter(magArray1,21,8)
    smoothedMag2 = sps.savgol_filter(magArray2,9,3)
    smoothedMag3 = sps.savgol_filter(magArray3,21,8)
    smoothedMag4 = sps.savgol_filter(magArray4,21,8)
    ##smoothedPhs = sps.savgol_filter(phsArray,3,1)

    textFile = open("S12 52 points/"+nameOfTxt, "w+")
    textFile.write("!Output from S12\n")
    textFile.write("!!S1P File: Measurement: S12:\n")
    textFile.write("# Hz S  dB   R 50 \n")

    ##print(freqArray[5])

    for a in range(noOfPoints):
        textFile.write(freqArray[a])
        textFile.write("\t")
        textFile.write(str(smoothedMag1[a]))
        textFile.write("\t")
        textFile.write(str(phsArray1[a]))
        textFile.write("\t")
        textFile.write(str(smoothedMag2[a]))
        textFile.write("\t")
        textFile.write(str(phsArray2[a]))
        textFile.write("\t")
        textFile.write(str(smoothedMag3[a]))
        textFile.write("\t")
        textFile.write(str(phsArray3[a]))
        textFile.write("\t")
        textFile.write(str(smoothedMag4[a]))
        textFile.write("\t")
        textFile.write(str(phsArray4[a]))
        textFile.write("\n")

    textFile.close()

    plt.ylim(-40,10)
    plt.xlim(0,3)
    plt.title('S12 (Reverse Gain/Isolation)')
    plt.grid(True)
    plt.ylabel('Magnitude (dB)')
    plt.xlabel('Frequency (GHz)')
    plt.legend(loc='upper right')
    ##plt.plot(np.linspace(0.1,2.7,noOfPoints),smoothedMag1, label='dBm', color = 'blue')
    plt.plot(np.linspace(0.1,2.7,noOfPoints),smoothedMag2, label='dB', color = 'orange')
    ##plt.plot(np.linspace(0.1,2.7,noOfPoints),smoothedMag3, label='dBm', color = 'green')
    ##plt.plot(np.linspace(0.1,2.7,noOfPoints),smoothedMag4, label='dBm', color = 'red')

    plt.show()
    S12Display.config(text = "Smoothing complete")


##################################S22#################################
###############################Averaging##############################
def AveragingS11():
    S22Display.config(text = "Collecting data...please wait")
    #-------------------------------------------#
    gainC = []
    ##plt.ion()
    cnt=0
    file = 0
    magSum = 0
    #-------------------------------------------#
    ##Hardware SPI configuration:
    SPI_PORT   = 0
    SPI_DEVICE = 0
    mcp = Adafruit_MCP3008.MCP3008(spi=SPI.SpiDev(SPI_PORT, SPI_DEVICE))
    #-------------------------------------------#
    def plotgainC():
        plt.ylim(-40,10)
        plt.xlim(0,3)
        plt.title('S22 (Output Return Loss)')
        plt.grid(True)
        plt.ylabel('Magnitude (dB)')
        plt.xlabel('Frequency (GHz)')
        pltplot(np.linspace(0.1,2.7,noOfPoints),gainC, 'rx-', label='dB')
        plt.legend(loc='upper right')
        plt.show()

    for turn in range(15):
        
        a = str(turn)
        b = ('result')
        c = ('.s1p')
        nameOfTxt = (b+a+c)
    ##    print("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\nOpening", b+a+c)
        textFile = open('S22 52 points/Noise averaging/'+nameOfTxt, "w+")
        initialFreq = 50000000
        textFile = open('S22 52 points/Noise averaging/'+nameOfTxt, "a")
        
        textFile.write("!Output from S22\n")
        textFile.write("!!S1P File: Measurement: S22:\n")
        textFile.write("# Hz S  dB   R 50 \n")

        for i in range(0,noOfPoints): #pre-load dummy data
            gainC.append(0)

        port.write(b"C1E1r1l100.0u2700.0s50.0t50.000[0.00]0.00^1")
        port.write(b"g")     

        for x in range(noOfPoints):
            
            initialFreq = initialFreq + 50000000
            textFile.write(str(initialFreq))
            textFile.write('\t')        
            
            values = [0]*8
            for i in range(8):
                values[i] = mcp.read_adc(i)
                    
            ##values[0] is the Vmag, values[1] is the Vphase, check the MCP3008 pin config
            
            
            #Magnitude----------------------------------#
            #-------------------------------------------#
            vMag = values[2]
            magRatio = ((vMag-559.03)/(17.915)) ##formula for line equation from readings
            r = round(magRatio, 2) ##round up value for ease of reading
            
    ##        print('Magnitude Ratio = ', r, 'dB')
            
            volRatio = round((10**(r/20)),2)
    ##        print('Voltage Ratio = ', volRatio)
            
            #-------------------------------------------# 
            #Phase--------------------------------------#
            #-------------------------------------------#
            vPhase = values[3]
    ##      phaseDiff = 0.1*(810-((vPhase/1024)*1.8))
            phaseDiff = ((vPhase/1024.0)*180)-90 #found this equation via pages.cs.wisc.edu/~timc/e/spna/index.html 
            phaseDiffInRad = (phaseDiff/180)*3.14159265359
            phi = round(phaseDiffInRad, 2)##round value up for ease of reading
            
    ##        print('Phase           = ', phi, 'Radians')
            
            #-------------------------------------------#
            #Converting to complex form-----------------#
            #-------------------------------------------#       
            realForm = volRatio * math.cos(phi)
            imaginaryForm = volRatio * math.sin(phi)
    ##        print('In complex form -> ',round(realForm,2), '+ ', round(imaginaryForm,2), 'j\n')
            
            textFile.write(str(magRatio))
            textFile.write('\t')
            textFile.write(str(phaseDiff)) 
            textFile.write('\n')

            
            gainC.append(r)
            gainC.pop(0)
                    
            time.sleep(0.05)
        
    ##  print(gainC)
    ##    drawnow(plotgainC)
        textFile.close()
    ##    print(str(round(((turn/15)*100),0))+'% complete')
        S22Display.config(text = str(round(((turn/15)*100),0))+'% complete')
        time.sleep(0.5)
        gainC = []

    S11Display.config(text = 'Data collection complete')
    ##print('Data collection complete\n\n\n\n')        
    ##switch off outputs of SynthHD
    port.write(b"C0E0r0C1E0r0")

    avgFile = open('S22 52 points/Noise averaging/averaged.s1p', "w+")
    initialFreq = 50000000
    newFreq = 50000000
    avgFile = open('S22 52 points/Noise averaging/averaged.s1p', "a")
        
    avgFile.write("!Output from S22\n")
    avgFile.write("!!S1P File: Measurement: S22:\n")
    avgFile.write("# Hz S  dB   R 50 \n")

    ##print("Initiating averaging")
        
    for row in range(noOfPoints):
        magSum = 0
        phsSum = 0
        for file in range(15):
            
            a = str(file)
            b = ('result')
            c = ('.s1p')
            nameOfTxt = (b+a+c)
    ##        print(b+a+c)

            currLines = []
            currFile = open('S22 52 points/Noise averaging/' + nameOfTxt, "r")
            currLines = [_.split() for _ in currFile.readlines()]
            
    ##        print(currLines[row+3][1])
            
            mag = float(currLines[row+3][1])
            phs = float(currLines[row+3][2])
            
            magSum = mag + magSum
            phsSum = phs + phsSum
    ##        print(magSum)
            
    ##        time.sleep(0.05)
            
        newFreq = newFreq + initialFreq
        avgFile.write(str(newFreq))
        avgFile.write("\t")    
        avg = 0
        avg = magSum/15
    ##    print(avg)
    ##    print('\n\n')
        avgFile.write(str(avg)+'\t')
        avg = 0
        avg = phsSum/15
    ##    print(avg)
    ##    print('\n\n')
        avgFile.write(str(avg)+'\n')
      
    S22Display.config(text = 'Averaging complete')  
    avgFile.close()

	
###################################S22################################
###############################CorrectData############################	
def CorrectDataS22():
    
    S22Display.config(text = "Processing. Please wait.")
    gainC = []
    plt.ion()
    cnt=0
    
    nameOfTxt = "attenuatorSampleS22.s1p"

    #-------------------------------------------#

    textFile = open("S22 52 points/"+nameOfTxt, "w+")
    textFile.write("!Output from S11\n")
    textFile.write("!!S1P File: Measurement: S22:\n")
    textFile.write("# Hz S  dB   R 50 \n")
    textFile.close()

    #-------------------------------------------#

    ##Hardware SPI configuration:
    SPI_PORT   = 0
    SPI_DEVICE = 0
    mcp = Adafruit_MCP3008.MCP3008(spi=SPI.SpiDev(SPI_PORT, SPI_DEVICE))

    #-------------------------------------------#
    def plotgainC():
        plt.ylim(-40,10)
        plt.xlim(0,3)
        plt.title('S22 (Output Return Loss)')
        plt.grid(True)
        plt.ylabel('Magnitude (dB)')
        plt.xlabel('Frequency (GHz)')
        pltplot(np.linspace(0.1,2.7,noOfPoints),gainC, 'rx-', label='dB')
        plt.legend(loc='upper right')
        plt.show()


    initialFreq = 0
    # open results from previous data, and create a matrix with them prevLines[x][y]
    prevLines = []
    textFile = open("S22 52 points/Noise averaging/averaged.s1p", "r")
    prevLines = [_.split() for _ in textFile.readlines()]
    

        
    # read the data file in as a list
    fin = open("S22 52 points/"+nameOfTxt, "r" )
    data_list = fin.readlines()
    fin.close()
         
    # remove list items from index 3 to 5 (inclusive)
    del data_list[3:noOfPoints+10]
        
    # write the changed data (list) to a file
    fout = open("S22 52 points/"+nameOfTxt, "w")
    fout.writelines(data_list)
    fout.close()

    textFile = open("S22 52 points/"+nameOfTxt, "a")

    for i in range(0,noOfPoints): #pre-load dummy data
        gainC.append(0)

    port.write(b"C1E1r1l100.0u2700.0s50.0t50.000[0.00]0.00^1")
    port.write(b"g")     

    for x in range(noOfPoints):
        initialFreq = initialFreq + 50000000
        print(initialFreq)
        textFile.write(str(initialFreq))
        textFile.write('\t')        
            
        values = [0]*8
        for i in range(8):
            values[i] = mcp.read_adc(i)
                    
    ##values[0] is the Vmag, values[1] is the Vphase, check the MCP3008 pin config
            
            
        #Magnitude----------------------------------#
        #-------------------------------------------#
        vMag = values[2]
        magRatio = ((vMag-559.03)/(17.915)) ##formula for line equation from readings
            
        tempR = prevLines[x+3][1] #temporary real value from result.s1p
        tempI = prevLines[x+3][2] #temporary imaginary value from result.s1p
        ftempR = float(tempR) #converting to float
        ftempI = float(tempI) #converting to float
    ##        prevMag = math.sqrt((ftempR*ftempR)+(ftempI*ftempI))
            
        print('\n\n\n\n\nMag = ', magRatio, 'dB')
        print('Prev = ', ftempR, 'dB')
            
        if magRatio < 0: #if value is negative
            if ftempR > 0:
                corr = magRatio + ftempR
                print('Magnitude Ratio = ', corr, 'dB')
            if ftempR < 0:
                corr = magRatio - ftempR
                print('Magnitude Ratio = ', corr, 'dB')                
        if magRatio > 0: #if value is positive
            if ftempR > 0:
                corr = magRatio - ftempR
                print('Magnitude Ratio = ', corr, 'dB')
            if ftempR < 0:
                corr = magRatio + ftempR
                print('Magnitude Ratio = ', corr, 'dB')   


        r = round(corr, 2) ##round up value for ease of reading
                   
        #-------------------------------------------# 
        #Phase--------------------------------------#
        #-------------------------------------------#
        vPhase = values[3]
    ##      phaseDiff = 0.1*(810-((vPhase/1024)*1.8))
        phaseDiff = ((vPhase/1024.0)*180)-90 #found this equation via pages.cs.wisc.edu/~timc/e/spna/index.html 
        phaseDiffInRad = (phaseDiff/180)*3.14159265359
        phi = round(phaseDiffInRad, 2)##round value up for ease of reading
            
        print('Phase           = ', phi, 'Radians')
            
        #-------------------------------------------#
        #Converting to complex form-----------------#
        #-------------------------------------------#       
        realForm = r * math.cos(phi)
        imaginaryForm = r * math.sin(phi) 
            
        textFile.write(str(magRatio))
        textFile.write('\t')
        textFile.write(str(phaseDiff)) 
        textFile.write('\n')
           
        print('Magnitude RatioII = ', r, 'dB')
        gainC.append(r)
        gainC.pop(0)
                    
        time.sleep(0.05)
    initialFreq = 0
    drawnow(plotgainC)
        
    textFile.close()
        
    S22Display.config(text = "Plot Displayed")
         
    gainC = []
        
port.write(b"C0E0r0C1E0r0")	

##################################S22#################################
############################PlotCalibrated############################
def PlotCalibratedS22():
    my_ideals = [\
            rf.Network('S22 52 points/IdealCalibVal/IdealShort.s1p'),
            rf.Network('S22 52 points/IdealCalibVal/IdealOpen.s1p'),
            rf.Network('S22 52 points/IdealCalibVal/IdealBBLoad.s1p'),
            ]

    my_measured = [\
            rf.Network('S22 52 points/MeasCalibVal/smoothMeasShort.s1p'),
            rf.Network('S22 52 points/MeasCalibVal/smoothMeasOpen.s1p'),
            rf.Network('S22 52 points/MeasCalibVal/smoothMeasBBLoad.s1p'),
            ]

    cal = rf.OnePort(\
            ideals = my_ideals,
            measured = my_measured,
            )

    cal.run()
    dut = rf.Network('S22 52 points/attenuatorSampleS22.s1p')

    dut_caled = cal.apply_cal(dut)
    dut_caled.name =  dut.name + ' Calibrated'
    dut_caled.write_touchstone()

    nameOfCalSignal =  dut.name + ' Calibrated.s1p' 

    ro_ori = rf.Network('S22 52 points/attenuatorSampleS22.s1p')
    ro_cal = rf.Network(nameOfCalSignal)

    figure()

    ro_ori.plot_s_db(label='Original')
    ro_cal.plot_s_db(label='Calibrated')

    draw();show();
    
    nameOfCalTxt = "attenuatorSampleS22 Calibrated.s1p"
    nameOfCaldBTxt = "Calibrated&IndBFormatS22.s1p"
    noOfPoints = 52

    realList = []
    imagList = []
    freqList = []

    #--------------------------------------------------#
    textFile = open(nameOfCalTxt, "r")
    original = [_.split() for _ in textFile.readlines()]

    for a in range(noOfPoints):
        real = original[a+3][1]
        realList.append(real) #this list holds all noise filled original values for Magnitude
        
    for b in range(noOfPoints):
        imag = original[b+3][2]
        imagList.append(imag) #this list holds all noise filled original values for Phase
        
    for c in range(noOfPoints):
        freq = original[c+3][0]
        freqList.append(freq)#this list holds all frequenciesfor a in (noOfPoints):

    rA = np.array(realList)
    iA = np.array(imagList)
    fA = np.array(freqList)

    realArray = rA.astype(np.float)
    imagArray = iA.astype(np.float)
    freqArray = fA.astype(np.float)
    #--------------------------------------------------#

    textFile = open(nameOfCaldBTxt, "w+")
    textFile.write("!Output from S22\n")
    textFile.write("!!S1P File: Measurement: S22:\n")
    textFile.write("# Hz S  dB   R 50 \n")

    for d in range(noOfPoints):
    ##    print(freqArray[d])
        rl = realArray[d]
        im = imagArray[d]
        mag = 20 * math.log10(math.sqrt((rl*rl)+(im*im)))
        phs = math.degrees(math.atan(im/rl))
        
        magS = str(mag)
        phsS = str(phs)
        
        textFile.write(str(freqArray[d]))
        textFile.write("\t")    
        textFile.write(magS)
        textFile.write("\t")
        textFile.write(phsS)
        textFile.write("\n")
        
    textFile.close()

##################################S22#################################
############################ResultSmoothing###########################
def ResultSmoothingS22():
    nameOfTxt = "smoothCalibrated&IndBFormatS22.s1p"
    noOfPoints = 52 #there is actually 52 points read

    values = []
    magList = []
    phsList = []
    freqList = []
    valm = []
    valp = []
    freq = []
    smoothed = []
    textFile = open("Calibrated&IndBFormatS22.s1p", "r")
    values = [_.split() for _ in textFile.readlines()] # open results from previous data, and create a matrix with them prevLines[x][y]
    a = 0
    b = 0
    c = 0
    d = 0

    for a in range(noOfPoints):
        valM = values[a+3][1]
        print(values[a+3][1])
        magList.append(valM) #this list holds all noise filled original values for Magnitude
        
    for b in range(noOfPoints):
        valP = values[b+3][2]
        phsList.append(valP) #this list holds all noise filled original values for Phase
        
    for c in range(noOfPoints):
        freq = values[c+3][0]
        freqList.append(freq)#this list holds all frequencies
        
    #converting all three into numpy arrays   
    magArray = np.asarray(magList)
    phsArray = np.asarray(phsList)
    freqArray = np.asarray(freqList)

    ##print(magArray)
    ##print(phsArray)
    ##print(freqArray)

    smoothedMag = sps.savgol_filter(magArray,21,8)
    ##smoothedPhs = sps.savgol_filter(phsArray,3,1)

    textFile = open(nameOfTxt, "w+")
    textFile.write("!Output from S22\n")
    textFile.write("!!S1P File: Measurement: S22:\n")
    textFile.write("# Hz S  dB   R 50 \n")

    ##print(freqArray[5])

    for d in range(noOfPoints):
        textFile.write(freqArray[d])
        textFile.write("\t")
        textFile.write(str(smoothedMag[d]))
        textFile.write("\t")
        textFile.write(str(phsArray[d]))
        textFile.write("\n")

    textFile.close()

    plt.ylim(-40,10)
    plt.xlim(0,3)
    plt.title('S22 (Output Return Loss)')
    plt.grid(True)
    plt.ylabel('Magnitude (dB)')
    plt.xlabel('Frequency (GHz)')
    pltplot(np.linspace(0.1,2.7,noOfPoints),gainC, 'rx-', label='dB')
    plt.legend(loc='upper right')
    plt.show()




def close_plot():
    plt.clf()
    plt.close('all')
    pylab.close('all')
    
def restart_program():
    python = sys.executable
    os.execl(python, python, * sys.argv)

def close_window(): 
    master.destroy()
    sys.exit()


######################################################################
###############################Interface##############################
######################################################################

tk.Label(master, text = "VNA", font = 'Gentona 20 bold').grid(pady=10, padx=10, row=0, column=1,sticky="E")
tk.Label(master, text = "INTERFACE", font = 'Gentona 20 bold').grid(pady=10, padx=10, row=0, column=2,sticky="W")



##################################S11#################################
S11Sign = tk.Label(master, text = "S11")
S11Sign.grid(pady=10, padx=10, ipadx=10, row=1, column=0)
S11Display = tk.Label(master, text = "Choose option to continue")
S11Display.grid(pady=10, padx=10, ipadx=10, row=2, column=0)

tk.Button(master, text='Averaging', command=AveragingThreadS11, width=15).grid(pady=10, padx=20, row=3, column=0)
tk.Button(master, text='Plot once', command=CorrectedDataThreadS11, width=15).grid(pady=10, padx=20, row=4, column=0)
tk.Button(master, text='Plot calibrated', command=PlotCalibratedThreadS11, width=15).grid(pady=10, padx=20, row=5, column=0)
tk.Button(master, text='Result smoothing', command=ResultSmoothingThreadS11, width=15).grid(pady=10, padx=20, row=6, column=0)
###################################################################### 

##################################S21#################################
S21Sign = tk.Label(master, text = "S21")
S21Sign.grid(pady=10, padx=10, ipadx=10, row=1, column=1)
S21Display = tk.Label(master, text = "Choose option to continue")
S21Display.grid(pady=10, padx=10, ipadx=10, row=2, column=1)

tk.Button(master, text='Averaging', command=AveragingThreadS21, width=15).grid(pady=10, padx=20, row=3, column=1)
tk.Button(master, text='Plot once', command=CorrectedDataThreadS21, width=15).grid(pady=10, padx=20, row=4, column=1)
tk.Button(master, text='Plot calibrated', command=PlotCalibratedThreadS21, width=15).grid(pady=10, padx=20, row=5, column=1)
tk.Button(master, text='Result smoothing', command=ResultSmoothingThreadS21, width=15).grid(pady=10, padx=20, row=6, column=1)
tk.Button(master, text='Close plot', command=close_plot, width=15).grid(pady=10, padx=20, row=7, column=1)
######################################################################

##################################S12#################################
S12Sign = tk.Label(master, text = "S12")
S12Sign.grid(pady=10, padx=10, ipadx=10, row=1, column=2)
S12Display = tk.Label(master, text = "Choose option to continue")
S12Display.grid(pady=10, padx=10, ipadx=10, row=2, column=2)

tk.Button(master, text='Averaging', command=AveragingThreadS12, width=15).grid(pady=10, padx=20, row=3, column=2)
tk.Button(master, text='Plot once', command=CorrectedDataThreadS12, width=15).grid(pady=10, padx=20, row=4, column=2)
tk.Button(master, text='Plot calibrated', command=PlotCalibratedThreadS12, width=15).grid(pady=10, padx=20, row=5, column=2)
tk.Button(master, text='Result smoothing', command=ResultSmoothingThreadS12, width=15).grid(pady=10, padx=20, row=6, column=2)
tk.Button(master, text='Close plot', command=close_plot_thread, width=15).grid(pady=10, padx=20, row=7, column=2)
######################################################################

##################################S22#################################
S22Sign = tk.Label(master, text = "S22")
S22Sign.grid(pady=10, padx=10, ipadx=10, row=1, column=3)
S22Display = tk.Label(master, text = "Choose option to continue")
S22Display.grid(pady=10, padx=10, ipadx=10, row=2, column=3)

tk.Button(master, text='Averaging', command=AveragingThreadS22, width=15).grid(pady=10, padx=20, row=3, column=3)
tk.Button(master, text='Plot once', command=CorrectedDataThreadS22, width=15).grid(pady=10, padx=20, row=4, column=3)
tk.Button(master, text='Plot calibrated', command=PlotCalibratedThreadS22, width=15).grid(pady=10, padx=20, row=5, column=3)
tk.Button(master, text='Result smoothing', command=ResultSmoothingThreadS22, width=15).grid(pady=10, padx=20, row=6, column=3)
######################################################################


ExitButton = tk.Button(master, text='Exit', command=close_window, width=15).grid(pady=10, padx=20, row=8, column=1)
ExitButton = tk.Button(master, text='Restart', command=close_window, width=15).grid(pady=10, padx=20, row=8, column=2)
tk.Label(master, text='NOTE:\nFor S11 ans S22 measurements, cable A goes to yellow module and cable B goes to blue module.\nFor S12 ans S21 measurements, cable B goes to yellow module and cable A goes to blue module.\nIf [Plot Once] does not work on the first try, press again.\nSame concept applies to the [Close plot] button.', justify = 'left').place(x=10, y=430)

master.mainloop( )

