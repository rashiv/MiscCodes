import numpy as np
import pylab as plt
import matplotlib as mpl
from matplotlib import rc
from scipy.optimize import curve_fit

def func(x, a, b,c):
    return a*(1.0/(1+np.exp(-b*(x-c))) )

def plotWithFill(dataAll,nruns):
    clasifications = {}
    mean_nonmono = 0.0;
    for currseed in dataAll:
        data = dataAll[currseed][1]
        data['Success %'] = (data['Success'])/nruns
        data['Time'] = data['TimeSinceNucleation']
        data['ErrD'] = (data['Success']-data['SuccessMin'])/nruns; data['ErrU'] = (data['SuccessMax']-data['Success'])/nruns
        max_time = np.min(data['Time'][ data['Success %'] > 0.95]); min_time = np.min(data['Time'][ data['Success %'] > 0.05])

        mean_err = np.mean( (data['ErrU']+data['ErrD'])/2.0 )
        transition_time = (max_time-min_time);

        data['Time'] = data['Time']/transition_time
        data['TimeAlt'] = data['Time']*6;
        data['TimeAlt'] = data['TimeAlt'].astype(np.int)

        data_binned = data.groupby(['TimeAlt']).mean()
        data_binned.index = data_binned.index/6.0

        succ_arr = np.array(data_binned['Success %'])
        min_change = np.min( succ_arr[1:]-succ_arr[:-1] )

        parameter, covariance_matrix = curve_fit(func, data['Time'], data['Success %'], maxfev=64500,p0=[1.0,-0.5,1.0])
        err_curr = ( np.mean( (data['Success %'] - func(data['Time'], *parameter))**2.0 ) )
        err_curr >
        if min_change < -0.02 and :
            clasifications[currseed] = True
            mean_nonmono += 1.0
        else:
            clasifications[currseed] = False
    mean_nonmono = mean_nonmono/float(len(dataAll.keys()))
    return clasifications

