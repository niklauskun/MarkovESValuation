import os
import pandas as pd
import numpy as np
from scipy.io import loadmat
import timeit

initial = timeit.default_timer()

def transition_matrix(df, state_num):
    n = state_num
    M = [[0]*n for _ in range(n)]
    for x in range(df.shape[1]):
        transitions = df[x]
        for (i,j) in zip(transitions,transitions[1:]):
            M[i][j] += 1

    #now convert to probabilities:
    for row in M:
        s = sum(row)
        if s > 0:
            row[:] = [f/s for f in row]
    return M


#RTP data
BASE_DIRECTORY = 'C:\\Users\\wenmi\\Desktop\\MarkovESValuation'
RESULTS_DIRECTORY = os.path.join(BASE_DIRECTORY, 'transition_matrix')
mat = loadmat(os.path.join(BASE_DIRECTORY, 'RTP_NYC_2010_2019.mat'))
RTP = pd.DataFrame(mat['RTP'])
resolution = 12 #slices in one hour
price_bar = 200 #price greater than price bar is price spike
state_gap = 10 #state gap
state_num = price_bar//state_gap + 2 #state number]
transition_matrix_num = 24

for i in range(transition_matrix_num):
    start = i*resolution
    end = (i+1)*resolution
    hour_data = (RTP[start : end]//state_gap+1).astype(int)
    hour_data[hour_data < 0] = 0
    hour_data[hour_data > state_num-2] = state_num-1
    m = transition_matrix(hour_data, state_num)
    np.savetxt(os.path.join(RESULTS_DIRECTORY,'matrix'+str(i)+".csv"), m, delimiter = ',')

stop = timeit.default_timer()

print('Time: ', stop - initial)

#import datetime

#initial_date = datetime.datetime.strptime('2010-01-01 00:00:00', '%Y-%m-%d %H:%M:%S')
#time_resolution = 300 #data resolution in sec, defalut by 5 mins#
#states_resolution = 1000

#for i,j in [(i,j) for i in range(df.shape[0]) for j in range(df.shape[1])]:
#    delta = datetime.timedelta(days = j,seconds = i * time_resolution)
#    t = initial_date + delta
#    d = df.iloc[i,j]
#    l = d//states_resolution
#    demand_rt = demand_rt.append({'time': t, 'demand(MW)': d, 'demand_level': l}, ignore_index=True)
#demand_rt.sort_values('time', inplace=True)
#demand_rt.to_csv(os.path.join(BASE_DIRECTORY,"real_time_demand.csv"), index=False)

