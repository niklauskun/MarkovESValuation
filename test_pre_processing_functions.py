import timeit
import os
import pandas as pd
import numpy as np
from scipy.io import loadmat

start = timeit.default_timer()

def transition_matrix(hour_data):
    n = states_num
    M = [[0]*n for _ in range(n)]
    for x in range(df.shape[1]+1):
        transitions = df.iloc[x]
        for (i,j) in zip(transitions,transitions[1:]):
            M[i][j] += 1

    #now convert to probabilities:
    for row in M:
        s = sum(row)
        if s > 0:
            row[:] = [f/s for f in row]
    return M

#test data
test_data = {'1':[0,0,2,2],
       '2':[0,2,1,0],
       '3':[0,1,1,1]}
states_num = 3 #number of states
df = pd.DataFrame(test_data)

m = transition_matrix(df)

#RTP data
BASE_DIRECTORY = 'C:\\Users\\wenmi\\Desktop\\MarkovESValuation'
mat = loadmat(os.path.join(BASE_DIRECTORY, 'RTP_NYC_2010_2019.mat'))
mat_df = pd.DataFrame(mat['RTP'])
#time_resolution = 60 #time resolution in mins
hour0 = mat_df[0:12]

for i in range(24):
    hour_data = mat_df[i*12 : (i+1)*12]
    #add code for ground to 21 nodes



    m = transition_matrix(hour_data)
    np.savetxt('new.csv', m, delimiter = ',')
    

stop = timeit.default_timer()

print('Time: ', stop - start)

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

