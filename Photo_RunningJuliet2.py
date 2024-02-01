import juliet
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.append('../')
import utils

############INPUT############
target = 'WASP6'
LightKurves = ['TESS_binnedOOT', 'ASASSNg_binnedClipped', 'ASASSNv_binnedClipped'] #save the data in 'GitHub/twin-planets/PhotometricMonitoring/PhotoMonReduc.ipynb'
Dists = [['uniform', 'fixed', 'uniform', 'uniform', 'loguniform'], None, None]
Hyperps = [[[-0.5,0.5], 1., [-200, 200], [-.1,.1], [1e-5,1e5]], None, None]
runName = 'binnedTESSnASASSNpass1'
PlotCorn, PlotGPfit, PlotPhaseFold = True, True, True
phase_num = 1
FontSize, TickSize = 20, 17
windowFuncTest = True
############

times, fluxes, fluxes_error = {},{},{}  #Create dictionaries for LC data
params, dists, hyperps = [], [], []
CombName =''
for LK in range(len(LightKurves)):
    insturment = LightKurves[LK].split('_')[0]
    CombName += insturment+'_'
    Data = np.load(target+'/'+target+'_'+LightKurves[LK]+'.npy')

    if windowFuncTest:
        FlUX = np.ones(len(Data[1]))
    else:
        FlUX = Data[1]
    times[insturment], fluxes[insturment], fluxes_error[insturment] = Data[0], FlUX, Data[2] # Save data into LC dictionaries

    # Name of the parameters to be fit:
    params += ['mflux_'+insturment, 'mdilution_'+insturment, 'sigma_w_'+insturment, 'GP_B_'+insturment, 'GP_C_'+insturment]#,'GP_L_TESS', 'GP_Prot_TESS']
    
    # Distribution for each of the parameters:
    #IMPORTANT!!!! ['GP_B','GP_L','GP_Prot','GP_C'] = ["amplitude", "timescale", "period", "factor"] Quasi-periodic kernels hyperparameters
    if Dists[LK]:
        dists += Dists[LK]
    else: #default values
        dists += ['uniform', 'fixed', 'uniform', 'uniform', 'uniform']

    # Hyperparameters of the distributions (mean and standard-deviation for normal
    # distributions, lower and upper limits for uniform and loguniform distributions, and
    # fixed values for fixed "distributions", which assume the parameter is fixed)
    if Hyperps[LK]:
        hyperps += Hyperps[LK]
    else:
        hyperps += [[-0.1,0.1], 1., [0, 8e4], [-0.1, 0.1], [0, 9e3]]

## Now to add the combined parameters of all photometric data
params +=  ['GP_L_'+CombName[:-1], 'GP_Prot_'+CombName[:-1]]  #['GP_L_TESS', 'GP_Prot_TESS']
dists +=  ['loguniform', 'uniform'] 
hyperps +=  [[1e-6,1e6], [5, 60]] 


OutFolder = 'Juliet/'+target+'/'+runName
if not os.path.isdir(OutFolder):
    os.makedirs(OutFolder)

priors = {}
for param, dist, hyperp in zip(params, dists, hyperps):
    priors[param] = {}
    priors[param]['distribution'], priors[param]['hyperparameters'] = dist, hyperp
dataset = juliet.load(priors=priors, t_lc=times, y_lc=fluxes,yerr_lc=fluxes_error, GP_regressors_lc=times,out_folder=OutFolder)

# Fit:
results = dataset.fit(n_live_points=1000)

############# To plot corner plots #############
if PlotCorn:
    import corner
    first_time = True
    posterior_names = []
    for i in range(len(params)):
        if dists[i] != 'fixed':
            if first_time:
                posterior_data = results.posteriors['posterior_samples'][params[i]]
                posterior_names.append(params[i])
                first_time = False
            else:
                posterior_data  = np.vstack((posterior_data, results.posteriors['posterior_samples'][params[i]]))
                posterior_names.append(params[i])
    posterior_data = posterior_data.T    
    print ("posterior_data.shape", posterior_data.shape)
    figure = corner.corner(posterior_data, labels=posterior_names)
    plt.savefig(OutFolder+'/Corner.png')
    plt.close()

############# To plot the data with the GP fit overplotted for each light curve #############
if PlotGPfit:
    def PlotLCnGPfit(subpltNum, title, mod_inf, LCs, TimeBnds): #function that will be called to make all the subplots of data
        #mod_inf =  [modelTimes, modelFlux, [1sig_upperBnd, 1sig_LowerBnd], [2sig_upperBnd, 2sig_LowerBnd], [3sig_upperBnd, 3sig_LowerBnd]]
        #To write the JD date in reduced form
        LenDiff = len(str(int(np.floor(LCs[0][0]-LCs[0][-1])))) #to get the order of mag number difference from the earliest and lastest date
        startTime = str(int(np.floor(LCs[0][0])))
        JDsubtract = int(startTime[:len(startTime)-LenDiff] + ('0'*LenDiff)) #use that order of mag difference as to subtract from the JD dates so tiimes don't have such high numbers        
        plt.subplot(subpltNum)
        Ylim = (np.min(LCs[1]-LCs[2]), np.max(LCs[1]+LCs[2]))
        Xlim = (TimeBnds[0]-5-JDsubtract, TimeBnds[1]+5-JDsubtract)
        # plt.title(title, fontweight='bold', fontsize=FontSize)
        plt.text(Xlim[0]+Xlim[0]*.01,Ylim[1]-Ylim[1]*.01, title, fontweight='bold', fontsize=FontSize, zorder=20) #if want to print the title on the figure, rather than an actual title (for the paper)
        plt.errorbar(LCs[0]-JDsubtract, LCs[1], LCs[2], ls='', marker='o', mfc='green', ms=7, mew=1, zorder=8) 
        plt.plot(mod_inf[0]-JDsubtract, mod_inf[1], 'r--', linewidth = 5, zorder=10)
        plt.fill_between(mod_inf[0]-JDsubtract,  mod_inf[2][0],  mod_inf[2][1], color='grey', alpha=0.5)
        plt.fill_between(mod_inf[0]-JDsubtract,  mod_inf[3][0],  mod_inf[3][1], color='grey', alpha=0.3)
        plt.fill_between(mod_inf[0]-JDsubtract,  mod_inf[4][0],  mod_inf[4][1], color='grey', alpha=0.1)
        plt.xlim(Xlim)
        # plt.ylim((np.mean(flux2)-3*np.std(flux2), np.mean(flux2)+3*np.std(flux2)))
        plt.ylim(Ylim)
        plt.ylabel('relative Flux', fontsize=FontSize,fontweight='bold')
        plt.xlabel('JD - '+str(JDsubtract), fontsize=FontSize,fontweight='bold')
        plt.xticks(fontsize=TickSize, fontweight='bold')
        plt.yticks(fontsize=TickSize, fontweight='bold')
        plt.tight_layout()
        return None

    #To first figure out how many light curves there will be:
    FigBounds = [] #to keep track of the upper and lower bound of a given Figure
    subLKs = [1]*len(LightKurves) #to keep track of how many sub light curves there are for each LightKurves dataset
    for curve in range(len(LightKurves)): 
        Data = np.load(target+'/'+target+'_'+LightKurves[curve]+'.npy')
        Difs = np.diff(Data[0])
        diffThreshold = 100*np.mean(Difs)
        BigDiffs = Difs[np.where(Difs > diffThreshold)[0]] #if the difference between sequential observations is more than 100 times the average time sampling, then split this data into two parts
        subLKs[curve] += len(BigDiffs)
        if len(BigDiffs) == 0: #if there are no big gaps, then the data range goes from the first observation to the last
            FigBounds.append([Data[0][0], Data[0][-1]])
        else:
            startT = Data[0][0]
            for BD in BigDiffs: 
                Gap = np.where(Difs == BD)[0][0] #since there will be only one of this, get only the int index 
                endT = Data[0][Gap]
                FigBounds.append([startT, endT])
                startT = Data[0][Gap+1]
            #need an extra iteration outside of the loop, because there will be one more number of time bounds per gaps. i.e. two gaps yields 3 time ranges, etc.
            FigBounds.append([startT, Data[0][-1]])
    n = np.sum(subLKs)
    plt.figure(1, figsize=(18,3+(n*3))) #add an additional 3 in height for every new row
    subIdxBase = str(np.sum(n))+'1' #always have one column of data and number of rows equal to the total number of subplots
    totCnt = 0
    #Now to make a subplot of the data and the model for EACH light curve regions:
    for F in range(len(LightKurves)): 
        Data = np.load(target+'/'+target+'_'+LightKurves[F]+'.npy')
        InstName = LightKurves[F].split('_')[0]
        modelTimes = dataset.times_lc[InstName]
        model, error68_up, error68_down = results.lc.evaluate(InstName, return_err=True) #1-sigma
        model, error95_up, error95_down = results.lc.evaluate(InstName, return_err=True, alpha=0.95) #2-sigma
        model, error995_up, error995_down = results.lc.evaluate(InstName, return_err=True, alpha=0.995) #3-sigma (a little less. 0.997 is 3sig, but that doesn't work)
        ModelInfo = [modelTimes, model, [error68_up, error68_down], [error95_up, error95_down], [error995_up, error995_down]]
        if subLKs[F] == 1:
            Titl = InstName
            PlotLCnGPfit(int(subIdxBase+str(totCnt+1)), Titl, ModelInfo, Data, FigBounds[totCnt])
            totCnt +=1
        else:
            for subs in range(subLKs[F]):
                Titl = InstName+'_group'+str(subs)
                PlotLCnGPfit(int(subIdxBase+str(totCnt+1)), Titl, ModelInfo, Data, FigBounds[totCnt])
                totCnt +=1
    plt.savefig(OutFolder+'/PhotometricMonitoringGPfit.png')
    print ("Saved 'PhotometricMonitoringGPfit.png' in path '"+OutFolder+"'")
    plt.close()

# ############# To plot the phase folded data, folded on the best period and over plot a fitted sine curve with that period #############
# if PlotGPfit:
#     #stuff over here;

