from __future__ import division
import numpy as np
import pandas
from datetime import datetime
import statsmodels.api as sm
import scipy.stats as ss
from matplotlib import pyplot as plt

def generateWeather(n):
    '''Generates n-years of mean temperatures and precipitation amounts in 3 consecutive 30-day periods post planting date'''

    # range of historical data
    dates = pandas.date_range(start=datetime(1980,1,1),end=datetime(2010,12,31))
    
    # load database with lat/long coordinates and associated plant dates
    database = np.loadtxt('AfricaDatabase2.csv',delimiter=',',skiprows=1)
        
    for j in range(np.shape(database)[0]):
        #print(j)
        data = np.loadtxt('Weather/agmerra/txt/lat' + str(np.round(database[j,1],3)) + '_long' + str(np.round(database[j,2],3)) + '.txt', \
            skiprows=4, usecols=[4,5])
        
        # find historical monthly mean temperatures and total precipitation amounts in 3 consecutive 30-day periods post planting
        # first find how many of the 3 30-day periods post plant date are complete and size monthly data matrices accordingly
        # if plant date is Feb 29th, change to Feb 28th since there is no Feb 29th in 2010
        if database[j,3] == 2 and database[j,4] == 29:
            database[j,4] = 28

        lastPlantDate = np.where(dates==datetime(2010,int(database[j,3]),int(database[j,4])))[0][0]
        remainingDays = len(dates) - lastPlantDate + 1
        if remainingDays < 30:
            month1 = np.zeros([30,2])
            month2 = np.zeros([30,2])
            month3 = np.zeros([30,2])
        elif remainingDays >= 30 and remainingDays < 60:
            month1 = np.zeros([31,2])
            month2 = np.zeros([30,2])
            month3 = np.zeros([30,2])
        elif remainingDays >= 60 and remainingDays < 90:
            month1 = np.zeros([31,2])
            month2 = np.zeros([31,2])
            month3 = np.zeros([30,2])
        else:
            month1 = np.zeros([31,2])
            month2 = np.zeros([31,2])
            month3 = np.zeros([31,2])
            
        # populate monthly data matrices
        for i in range(np.shape(month1)[0]):
            day1 = np.where(dates==datetime(1980+i,int(database[j,3]),int(database[j,4])))[0][0]
            month1[i,0] = np.mean(data[day1:(day1+30),0])
            month1[i,1] = np.sum(data[day1:(day1+30),1])
            
        for i in range(np.shape(month2)[0]):
            day1 = np.where(dates==datetime(1980+i,int(database[j,3]),int(database[j,4])))[0][0]
            month2[i,0] = np.mean(data[(day1+30):(day1+60),0])
            month2[i,1] = np.sum(data[(day1+30):(day1+60),1])
            
        for i in range(np.shape(month3)[0]):
            day1 = np.where(dates==datetime(1980+i,int(database[j,3]),int(database[j,4])))[0][0]
            month3[i,0] = np.mean(data[(day1+60):(day1+90),0])
            month3[i,1] = np.sum(data[(day1+60):(day1+90),1])
            
        # fit linear trend to temperature data and calculate temperature anomalies from trend
        months = [month1, month2, month3]
        pars = []
        for i in range(len(months)):
            months[i][:,0], tempPars = removeTrend(months[i][:,0])
            pars.append(tempPars)
            
        # fit normal distribution to temperature anomalies, gamma distribution to total precipitation
        # and use a Gaussian copula to capture correlation across the 2 variables each month
        mus = np.zeros(3) # should be zero because calculating anomalies
        sigmas = np.zeros(3)
        for i in range(len(months)):
            mus[i], sigmas[i] = ss.norm.fit(months[i][:,0], floc=0) # freeze location parameter at 0 because fitting anomalies
        
        alphas = np.zeros(3)
        betas = np.zeros(3)
        for i in range(len(months)):
            # change any monthly total precipitation amounts of 0 to 0.1 so there is no error in estimation
            months[i][:,1] = np.where(months[i][:,1]>0.0, months[i][:,1], 0.1)
            alphas[i], loc, betas[i] = ss.gamma.fit(months[i][:,1], floc=0) # freeze location parameter at 0
            
        # fit bivariate Gaussian copula to monthly temperature anomalies and total monthly precipitation in each 30-day period post plant date
        percentiles = []
        for i in range(len(months)):
            tempPcts = ss.norm.cdf(months[i][:,0], loc=0, scale=sigmas[i]) # percentiles of normal marginal
            pptPcts = ss.gamma.cdf(months[i][:,1], a=alphas[i], loc=0, scale=betas[i]) # percentiles of gamma marginal     
            percentiles.append(np.array([tempPcts,pptPcts]))
            
        zscores = []
        rhos = np.zeros(3)
        for i in range(len(percentiles)):
            temp_zs = ss.norm.ppf(percentiles[i][0], loc=0, scale=1) # z-scores of standard normal
            ppt_zs = ss.norm.ppf(percentiles[i][1], loc=0, scale=1) # z-scores of standard normal
            zscores.append(np.array([temp_zs,ppt_zs]))
            rhos[i] = np.corrcoef(percentiles[i][0],percentiles[i][1])[0][1] # Spearman rank correlation
        
        # generate bivariate 3-period time series of z-scores from bivariate Gaussian copula
        # generate z1 (1st period) unconditionally (so each year assumed independent), then calculate z2 and z3 from
        # z2 = A1*z1 + B1*eps ; z3 = A2*z2 + B2*eps (so within-year, period-to-period correlation is preserved)
        A1, B1 = calcA_and_B(zscores[0], zscores[1])
        A2, B2 = calcA_and_B(zscores[1], zscores[2])
        
        z1 = np.random.multivariate_normal(np.zeros(2),np.array([[1,rhos[0]],[rhos[0],1]]),n)
        z2 = np.zeros(np.shape(z1))
        for i in range(n):
            z2[i] = np.dot(A1,z1[i,:]) + np.dot(B1,np.random.normal(0,1,2))
            
        z3 = np.zeros(np.shape(z2))
        for i in range(n):
            z3[i] = np.dot(A2,z2[i,:]) + np.dot(B2,np.random.normal(0,1,2))
        
        # transform z-scores back to temperature anomalies and precipitation amounts        
        zs = [z1, z2, z3]
        tempData = []
        precipData = []
        for i in range(len(zs)):
            tempPcts = ss.norm.cdf(zs[i][:,0], loc=0, scale=1) # convert z-scores to percentiles
            temps = ss.norm.ppf(tempPcts, loc=0, scale = sigmas[i]) # convert percentiles to temp anomalies
            tempData.append(temps)
            
            pptPcts = ss.norm.cdf(zs[i][:,1], loc=0, scale=1) # convert z-scores to percentiles
            precip = ss.gamma.ppf(pptPcts, a=alphas[i], loc=0, scale=betas[i])
            precipData.append(precip)
                
            #compareWeather(months, tempData, precipData, sigmas, alphas, betas, j)
            
        # add mean temperature predicted by trend for 2017 to simulated temperature anomalies
        for i in range(len(tempData)):
            tempData[i] = tempData[i] + pars[i][0]*37 + pars[i][1] # 37 years since 1980
            
        # write weather data to file - one file for mean monthly temperatures and one for total monthly precipitation amounts
        # each row is a different simulated year with 3 columns for the 3 30-day periods post plant date
        f = open('Weather/synthetic/lat' + str(np.round(database[j,1],3)) + '_long' + str(np.round(database[j,2],3)) + '_precip.txt','w')
        for i in range(n):
            f.write(str(np.round(precipData[0][i],1)) + ' ' + str(np.round(precipData[1][i],1)) + ' ' + str(np.round(precipData[2][i],1)) + '\n')
            
        f.close()
        
        f = open('Weather/synthetic/lat' + str(np.round(database[j,1],3)) + '_long' + str(np.round(database[j,2],3)) + '_temp.txt','w')
        for i in range(n):
            f.write(str(np.round(tempData[0][i],1)) + ' ' + str(np.round(tempData[1][i],1)) + ' ' + str(np.round(tempData[2][i],1)) + '\n')
            
        f.close()        
        
    return None

def compareWeather(months, tempData, precipData, sigmas, alphas, betas, rowNo):
    fig = plt.figure()
    ax = fig.add_subplot(3,2,1) # temp on left, precip on right; rows = months
    ax.hist(months[0][:,0], color='b', alpha=0.5, density=True)
    ax.hist(tempData[0], color='g', alpha=0.5, density=True)
    maxT = max(np.max(tempData[0]),np.max(months[0][:,0]))
    minT = min(np.min(tempData[0]),np.min(months[0][:,0]))
    x = np.arange(minT, maxT, (maxT-minT)/100.0)
    ax.plot(x,ss.norm.pdf(x, loc=0, scale=sigmas[0]),c='k',linewidth=2)
    ax.set_title('Temperature Month 1')
    
    ax = fig.add_subplot(3,2,2) # temp on left, precip on right; rows = months
    ax.hist(months[0][:,1], color='b', alpha=0.5, density=True)
    ax.hist(precipData[0], color='g', alpha=0.5, density=True)
    maxP = max(np.max(precipData[0]),np.max(months[0][:,1]))
    minP = min(np.min(precipData[0]),np.min(months[0][:,1]))
    x = np.arange(minP, maxP, (maxP-minP)/100.0)
    ax.plot(x,ss.gamma.pdf(x, a=alphas[0], loc=0, scale=betas[0]),c='k',linewidth=2)
    ax.set_title('Precipitation Month 1')
    
    ax = fig.add_subplot(3,2,3) # temp on left, precip on right; rows = months
    ax.hist(months[1][:,0], color='b', alpha=0.5, density=True)
    ax.hist(tempData[1], color='g', alpha=0.5, density=True)
    maxT = max(np.max(tempData[1]),np.max(months[1][:,0]))
    minT = min(np.min(tempData[1]),np.min(months[1][:,0]))
    x = np.arange(minT, maxT, (maxT-minT)/100.0)
    ax.plot(x,ss.norm.pdf(x, loc=0, scale=sigmas[1]),c='k',linewidth=2)
    ax.set_title('Temperature Month 2')
    
    ax = fig.add_subplot(3,2,4) # temp on left, precip on right; rows = months
    ax.hist(months[1][:,1], color='b', alpha=0.5, density=True)
    ax.hist(precipData[1], color='g', alpha=0.5, density=True)
    maxP = max(np.max(precipData[1]),np.max(months[1][:,1]))
    minP = min(np.min(precipData[1]),np.min(months[1][:,1]))
    x = np.arange(minP, maxP, (maxP-minP)/100.0)
    ax.plot(x,ss.gamma.pdf(x, a=alphas[1], loc=0, scale=betas[1]),c='k',linewidth=2)
    ax.set_title('Precipitation Month 2')
    
    ax = fig.add_subplot(3,2,5) # temp on left, precip on right; rows = months
    ax.hist(months[2][:,0], color='b', alpha=0.5, density=True)
    ax.hist(tempData[2], color='g', alpha=0.5, density=True)
    maxT = max(np.max(tempData[2]),np.max(months[2][:,0]))
    minT = min(np.min(tempData[2]),np.min(months[2][:,0]))
    x = np.arange(minT, maxT, (maxT-minT)/100.0)
    ax.plot(x,ss.norm.pdf(x, loc=0, scale=sigmas[2]),c='k',linewidth=2)
    ax.set_title('Temperature Month 3')
    
    ax = fig.add_subplot(3,2,6) # temp on left, precip on right; rows = months
    ax.hist(months[2][:,1], color='b', alpha=0.5, density=True, label='Historical')
    ax.hist(precipData[2], color='g', alpha=0.5, density=True, label='Synthetic')
    maxP = max(np.max(precipData[2]),np.max(months[2][:,1]))
    minP = min(np.min(precipData[2]),np.min(months[2][:,1]))
    x = np.arange(minP, maxP, (maxP-minP)/100.0)
    ax.plot(x,ss.gamma.pdf(x, a=alphas[2], loc=0, scale=betas[2]),c='k',linewidth=2)
    ax.set_title('Precipitation Month 3')
    
    handles, labels = plt.gca().get_legend_handles_labels()
    fig.subplots_adjust(bottom=0.15)
    fig.legend(handles, labels, loc='lower center',ncol=2,fontsize=16,frameon=True)
    fig.set_size_inches([10.5,11.5])
    fig.savefig('Weather/validation/Row' + str(rowNo) + '.png')
    fig.clf()
        
    return None
    
def removeTrend(y):
	X = np.arange(1,len(y)+1,1)
	X = sm.add_constant(X,prepend=False)
	trend = sm.OLS(y,X).fit()
	pars = trend.params
	z = np.zeros(len(y))
	for i in range(len(z)):
		z[i] = y[i] - pars[1] - pars[0]*(i+1)

	return z, pars
 
def calcA_and_B(X1, X2):
    X = np.transpose(np.concatenate((X1,X2),1)) # combine both months to calculate covariance matrix
    S = np.cov(np.transpose(X)) # temp and precip covariance matrix across both months
    
    S1 = np.zeros(np.shape(S)) # lag=one-month temp and precip covariance matrix
    # if less columns in second month because it didn't complete before end of 2010, remove last row of X1
    if np.shape(X2)[1] < np.shape(X1)[1]:
        X1 = X1[:,0:-1]

    S1[0,0] = np.cov(X2[0,:],X1[0,:])[0][1] # temp w/ previous temp
    S1[0,1] = np.cov(X2[0,:],X1[1,:])[0][1] # temp w/ previous precip
    S1[1,0] = np.cov(X2[1,:],X1[0,:])[0][1] # precip w/ previous temp
    S1[1,1] = np.cov(X2[1,:],X1[1,:])[0][1] # precip w/ previous precip
    
    A = np.dot(S1,np.linalg.inv(S))
    # check if positive definite
    if isPD(S - np.dot(A, np.transpose(S1))) == True:
        B = np.linalg.cholesky(S - np.dot(A, np.transpose(S1)))
    else:
        B = nearestPD(S - np.dot(A, np.transpose(S1)))
    
    return A, B

# fix non-positive definite matrices (from https://stackoverflow.com/questions/43238173/python-convert-matrix-to-positive-semi-definite)
def nearestPD(A):
    '''Find the nearest positive-definite matrix to input

    A Python/Numpy port of John D'Errico's `nearestSPD` MATLAB code [1], which
    credits [2].

    [1] https://www.mathworks.com/matlabcentral/fileexchange/42885-nearestspd

    [2] N.J. Higham, "Computing a nearest symmetric positive semidefinite
    matrix" (1988): https://doi.org/10.1016/0024-3795(88)90223-6'''

    B = (A + A.T) / 2
    _, s, V = np.linalg.svd(B)

    H = np.dot(V.T, np.dot(np.diag(s), V))

    A2 = (B + H) / 2

    A3 = (A2 + A2.T) / 2

    if isPD(A3):
        return A3

    spacing = np.spacing(np.linalg.norm(A))
    # The above is different from [1]. It appears that MATLAB's `chol` Cholesky
    # decomposition will accept matrixes with exactly 0-eigenvalue, whereas
    # Numpy's will not. So where [1] uses `eps(mineig)` (where `eps` is Matlab
    # for `np.spacing`), we use the above definition. CAVEAT: our `spacing`
    # will be much larger than [1]'s `eps(mineig)`, since `mineig` is usually on
    # the order of 1e-16, and `eps(1e-16)` is on the order of 1e-34, whereas
    # `spacing` will, for Gaussian random matrixes of small dimension, be on
    # othe order of 1e-16. In practice, both ways converge, as the unit test
    # below suggests.
    I = np.eye(A.shape[0])
    k = 1
    while not isPD(A3):
        mineig = np.min(np.real(np.linalg.eigvals(A3)))
        A3 += I * (-mineig * k**2 + spacing)
        k += 1

    return A3

def isPD(B):
    '''Returns true when input is positive-definite, via Cholesky decomposition'''
    try:
        _ = np.linalg.cholesky(B)
        return True
    except np.linalg.LinAlgError:
        return False
    
generateWeather(1000)