import sys
import numpy as np
import os
import glob

# Building the Weight Matrix Ground Truth
LTP2 = 1.0e-3
LTD2 = -1.0e-2
Init_MFDCN = 0.07

Weight_Matrix = []
MF_Activity=np.loadtxt("MF_12.dat")
PC_Activity=np.loadtxt("PC_40.dat")

for t in range(12):
    t_array=MF_Activity[t]
    count=0
    LTD_amount = []
    LTD_time = []
    LTP_amount = []
    LTP_time = []
    CurrentWeight = []
    LTD_index = 0
    for t2 in range(40):
        t2_array=PC_Activity[t2]
        for n,SpikesMF in enumerate(t_array):
            LTP_time.append(SpikesMF)
            LTP_amount.append(LTP2)
            for m,SpikesPC in enumerate(t2_array):
                sd = SpikesMF - SpikesPC
                if sd > 0 and sd <= 10:
                    LTD_amount.append(LTD2*np.exp(np.fabs(400.0*sd/1000.0))*pow(np.cos(sd/1000.0),2))
                    LTD_time.append(SpikesMF)
                elif sd <=0 and sd >=-10:
                    LTD_amount.append(LTD2*np.exp(np.fabs(400.0*sd/1000.0))*pow(np.cos(sd/1000.0),2))
                    if n+1<len(t_array):
                        next_MF_spike=t_array[n+1]
                        posto=n+1
                        while(next_MF_spike < SpikesPC and posto<len(t_array)):
                            next_MF_spike=t_array[posto]
                            posto+=1
                        LTD_time.append(next_MF_spike)
                        if posto>=len(t_array):
                            LTD_amount.pop()
                            LTD_time.pop()
                    else:
                        LTD_amount.pop()   
      
        LTD_amount=np.array([x for _,x in sorted(zip(LTD_time,LTD_amount))])
        LTD_time = sorted(LTD_time)
        LTD_time_u = np.unique(LTD_time)
        amou = []
        for LTD_times_u in LTD_time_u:
             ind = np.where(LTD_time == LTD_times_u)[0]
             amou.append(np.sum(LTD_amount[ind]))
        CurrentWeight = Init_MFDCN
        for x,LTPtimes in enumerate(LTP_time):
             CurrentWeight += LTP_amount[x]
             print(LTPtimes, LTD_times_u)
             CurrentLTD = amou[np.where(LTPtimes == LTD_times_u)[0]]
             if not len(CurrentLTD):
                 CurrentWeight += CurrentLTD
             print(t+1, 53, t_array[n]+1, round(CurrentWeight,3))

'''
    count = 0
    for n,LTD_times_u in enumerate(LTD_time_u):
        ind =  np.where(t_array>LTD_times_u)[0]
        for index in ind:
            CurrentWeight[index]=CurrentWeight[index]+amou[n]
    for n,w in enumerate(CurrentWeight): 
        Weight_Matrix.append([t+1, 65601, t_array[n]+1, round(w,3)])

Weight_Matrix=np.array(Weight_Matrix)

SimResults1 = np.loadtxt("PFPC1.csv")
SimResults4 = np.loadtxt("PFPC4.csv")

print( ["Length Ground Truth: " +  str(len(Weight_Matrix)) + " Length Simulation1 Results: " + str(len(SimResults1)) + " Length Simulation4 Results: " + str(len(SimResults4)) ])

SimResults1=np.matrix(SimResults1)
SimResults4=np.matrix(SimResults4)
Weight_Matrix=np.matrix(Weight_Matrix)

Weight_Matrix.sort(axis=0)
SimResults1.sort(axis=0)
SimResults4.sort(axis=0)

Difference1 = SimResults1-Weight_Matrix
# Get rid of approximation errors
Difference1[np.where(np.bitwise_and(Difference1 < 1.01e-3, Difference1 > -1.01e-3 ))]=0.0
Error1 = np.sum(Difference1)

Difference4 = SimResults4-Weight_Matrix
# Get rid of approximation errors
Difference4[np.where(np.bitwise_and(Difference4 < 1.01e-3, Difference4 > -1.01e-3 ))]=0.0
Error4 = np.sum(Difference4)

Error = Error1 + Error4

if Error == 0.0:
    sys.exit(0)
else:
   sys.exit(-1)
'''
