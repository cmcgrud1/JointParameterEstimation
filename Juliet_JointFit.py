import juliet
import numpy as np
import matplotlib.pyplot as plt
import lightkurve as lk
import os
import sys
import corner
from uncertainties import ufloat #used for error propergation. It uses linear error propagation theory
from uncertainties import umath #documentation: https://pythonhosted.org/uncertainties/_downloads/uncertaintiesPythonPackage.pdf
sys.path.append('/Users/chimamcgruder/Research/GitHub/twin-planets/')
import utils

############### Target input info ###############
#General info
ID, Survey = '6', 'WASP' #number ID of star, #name of survey that discovered the planet
Period, Perd_unc = 3.3610060,0.001
t0, t0_unc = 2454596.43267,0.0011
fit_ecc, Nlive = False, 300

#Plotting Flags
plot_corner, plot_single_trans, plot_phase_trans, plot_RV = True, True, True, True

#Stellar Parameters
Rs, Ms = [0.856, 0.028], [0.846, 0.032] #[value,uncertainity] #units R_sun, M_sun

#RV info
RV_insturments = ['CORALIE', 'HARPS']
sys_vel_range = {'CORALIE':[-10,10], 'HARPS':[-10,10]}
K_rang = [50.,90.] #assuming unifrom priors in [m/s], assuming systemic velocity is already subtracted

#### PATH INFO
target = Survey+'-'+ID
RV_dat_file = str(Survey)+str(ID)+'/'+str(Survey[0])+str(ID)+'RV_sysCorr_noRM.dat'
Fi = 'fix'
if fit_ecc:
    Fi = 'fit'
final_path = str(Survey)+str(ID)+'/'+str(Survey[0])+str(ID)+'_jointfit_'+Fi+'E'

# #TESS info. For now just using TESS data, but will change that later
# TESS_data = [utils.get_tess_lc(target, obs=1), utils.get_tess_lc(target, obs=7)] #This doesn't work in Hydra. Have to save TESS data as .npy file
TESS_file = str(Survey)+str(ID)+'/'+str(Survey[0])+str(ID)+'_TESS.npy'

#############################################



####### First define the TRANSIT priors (for now just using TESS data) #######
priors = {}

params = ['P_p1','t0_p1','r1_p1','r2_p1','q1_TESS','q2_TESS','ecc_p1','omega_p1',\
              'rho', 'mdilution_TESS', 'mflux_TESS', 'sigma_w_TESS'] # Name of the parameters to be fit

if fit_ecc:
    ecc_prior, omgea_prior = 'uniform', 'uniform'
    ecc_bnds, omgea_bnds = [0,.1], [0,360] #assuming if not 0, eccentricity, close to 0. omega can go from 0 to 360, so allow that to float
else:
    ecc_prior, omgea_prior = 'fixed', 'fixed'
    ecc_bnds, omgea_bnds = 0, 90    
dists = ['normal','normal','uniform','uniform','uniform','uniform',ecc_prior, omgea_prior,\
                 'loguniform', 'fixed', 'normal', 'loguniform'] # Distribution for each of the parameters

# Hyperparameters of the distributions (mean and standard-deviation for normal
# distributions, lower and upper limits for uniform and loguniform distributions, and
# fixed values for fixed "distributions", which assume the parameter is fixed). Note prior
# on t0 has an added 2457000 to convert from TESS JD to JD:
hyperps = [[Period, Perd_unc], [t0, t0_unc], [0.,1], [0.,1.], [0.,1.], [0.,1.], ecc_bnds, omgea_bnds,\
                   [100., 10000.], 1.0, [0.,0.1], [0.1, 1000.]]

# Populate the priors dictionary:
for param, dist, hyperp in zip(params, dists, hyperps):
    priors[param] = {}
    priors[param]['distribution'], priors[param]['hyperparameters'] = dist, hyperp

    
####### Now define the RV priors #######
params, dists, hyperps = ['K_p1'], ['uniform'], [K_rang,]
for inst in RV_insturments:
    params.append('mu_'+inst), params.append('sigma_w_'+inst)
    dists.append('uniform'), dists.append('loguniform') #Distributions
    hyperps.append(sys_vel_range[inst]), hyperps.append([1e-3, 100.]) # Hyperparameters


# Populate the priors dictionary:
for param, dist, hyperp in zip(params, dists, hyperps):
    priors[param] = {}
    priors[param]['distribution'], priors[param]['hyperparameters'] = dist, hyperp

# First put TESS photometric data in dictionaries, add 2457000 to the times to convert from TESS JD to JD:
times, fluxes, fluxes_error = {},{},{}
# for tess in TESS_data:
#   times['TESS'] = np.concatenate([TESS_data[0].time.value, TESS_data[1].time.value]) + 2457000
#   fluxes['TESS'] = np.concatenate([TESS_data[0].flux.value, TESS_data[1].flux.value])
#   fluxes_error['TESS'] = np.concatenate([TESS_data[0].flux_err.value, TESS_data[1].flux_err.value])
#Have to save the TESS data in .npy form, using that data here
Data = np.load(TESS_file)
times['TESS'], fluxes['TESS'], fluxes_error['TESS'] = Data[0,:], Data[1,:], Data[2,:]

# RV data is given in a file, so let's just pass the filename to juliet and load the dataset:
dataset = juliet.load(priors=priors, t_lc=times, y_lc=fluxes, yerr_lc=fluxes_error,\
                      rvfilename=RV_dat_file, out_folder=final_path) #RV data

# And now let's fit it!
results = dataset.fit(n_live_points=Nlive)

####### To plot the joint fit corner plot #######
if plot_corner:
    first_time = True
    posterior_names = []
    for par in list(priors.keys()):
        if priors[par]['distribution'] != 'fixed':
            if first_time:
                posterior_data = results.posteriors['posterior_samples'][par]
                posterior_names.append(par)
                first_time = False
            else:
                posterior_data  = np.vstack((posterior_data, results.posteriors['posterior_samples'][par]))
                posterior_names.append(par)
    posterior_data = posterior_data.T    
    figure = corner.corner(posterior_data, labels=posterior_names, quantiles=[0.15865,0.50, 0.84135])
    plt.savefig(final_path+'/Corner.png')
    plt.close()


####### To plot the Transit results #######
if plot_single_trans: #first plot each individual transit
    # Extract median model and the ones that cover the 68% credibility band around it:
    transit_model, transit_up68, transit_low68  = results.lc.evaluate('TESS', return_err=True)

    #Here to plot the TESS data
    PstSum =np.loadtxt(final_path+'/posteriors.dat', 'str') #.dat file containing the summary of the posterior results

    #To get joint parameters:
    P, t0 = np.array(PstSum[np.where(PstSum[:,0] == 'P_p1')[0][0],1:],dtype=np.float), np.array(PstSum[np.where(PstSum[:,0] == 't0_p1')[0][0],1:],dtype=np.float)
    t0_cent = t0[0]+(P[0]/2)#mid_transit+(P[0]/2). So centered the phased data on the tranist (when RV == 0)

    # To get transit parameters:
    RpRs, aRs =np.array(PstSum[np.where(PstSum[:,0] == 'p_p1')[0][0],1:],dtype=np.float), np.array(PstSum[np.where(PstSum[:,0] == 'a_p1')[0][0],1:],dtype=np.float)
    inc, rho_s = np.array(PstSum[np.where(PstSum[:,0] == 'inc_p1')[0][0],1:],dtype=np.float), np.array(PstSum[np.where(PstSum[:,0] == 'rho')[0][0],1:],dtype=np.float)
    b = np.array(PstSum[np.where(PstSum[:,0] == 'b_p1')[0][0],1:],dtype=np.float)

    T_dur_min = utils.TransDur(ufloat(P[0],np.max(P[1:])), ufloat(Rs[0],Rs[1]), ufloat(b[0],np.max(b[1:])), ufloat(RpRs[0],np.max(RpRs[1:])), ufloat(aRs[0],np.max(aRs[1:])))/60 #to convert to minutes
    print ("Estimated transit duration:", str(np.round(T_dur_min.n/(60),3))+"[hr] \u00B1", str(np.round(T_dur_min.s,2))+"[min]")

    #To 1st determine how many transit are observed (obs_trans)
    thresh = 15 #the limiting number of datapoints within a transit time to be considered that the transit was observed 
    mid_trans, mid_t =[], t0[0]#NOTE: make sure using first transit as t0 in fit!!!
    transits = (dataset.times_lc['TESS'][-1]-t0[0])/P[0] # to get the number of transits in our observations
    transits_int, T_dur_day= int(np.ceil(transits)), T_dur_min/(60*24)
    for t in range(transits_int): 
        transitRange = [mid_t-(T_dur_day.n/2)+(T_dur_day.s*2), mid_t+(T_dur_day.n/2)-(T_dur_day.s*2)] #to find the time frame where the transit is occuring. subtracting 2*std to make sure we are well within the transit
        inTtimes = np.where( (transitRange[0] < dataset.times_lc['TESS']) & (dataset.times_lc['TESS'] < transitRange[1]) )[0] # to get all the times that are in that transit
        if len(inTtimes) > thresh:
            mid_trans.append(mid_t)
        mid_t+=P[0]
    print ("observed transits:", len(mid_trans))

    #Now to plot those transits
    # total_datapoints = 70 #how many datapoints in the LC. if more than this, bin the data
    rows = int(np.ceil(len(mid_trans)/2)) #want only 2 columns
    fig = plt.figure(figsize=(18,5+(rows*5))) #make the figure large enough to fit all plots
    TwoPrtTitle = [target,'photometric transits'] # to split the title in 2 parts, since it will be printed on both headers
    for T in range(len(mid_trans)):
    #     if T < 2: #print the title on the 1st row (1st 2 columns) #not working!!!
    #         plt.title(TwoPrtTitle[T], fontsize=20, fontweight='bold')
        halfTransDur = (T_dur_day.n/2)+(T_dur_day.s) #half of the transit duration, plus the uncertainity, to ensure encompases full transit
        transit_Range = [mid_trans[T]-(2*halfTransDur), mid_trans[T]+(2*halfTransDur)] #need same amount of time in transit as out, hence the 2*
        Ttimes = np.where( (transit_Range[0] < dataset.times_lc['TESS']) & (dataset.times_lc['TESS'] < transit_Range[1]) )[0]
        tim, flx, err = dataset.times_lc['TESS'][Ttimes], dataset.data_lc['TESS'][Ttimes],  dataset.errors_lc['TESS'][Ttimes]
        position = ((T//2)*4, T%2) # position of main lc figure. (row, column): (row poisiton is index*4 because each lc takes 4 rows, even indeces are on the right side of the table)
        #To plot the main LC:
        plt.subplot2grid((4*rows, 2), position, rowspan=3, colspan=1)
        if position[0]%2 == 0: #only plot y-axis on left columns of LC plots
            plt.ylabel("relative flux")
    ### Actually don't need to bin the non-phased data
    #     if len(tim) > total_datapoints: #this will never happen with TESS data, just perparing for when I use other data
    #         plt.plot(tim, flx,'k.', alpha=0.1) #1st plot unbinned data
    #         bin_s = int(np.round(len(tim)/total_datapoints)) # to get the size of a bin
    #         tim, flx, err = utils.BinDat(tim, flx, err, Bin=bin_s)
        mod = transit_model[Ttimes]
        plt.errorbar(tim, flx, yerr=err, fmt='o', mec='cornflowerblue', ecolor='cornflowerblue', elinewidth=3, mfc='white', ms=7, label='t0 = '+str(round(mid_trans[T],5))+'[days]')
        plt.plot(tim, mod, color='black',zorder=10) # Plot the median model
        #To plot the residuals:
        pos_resd = (position[0]+3, position[1])
        plt.subplot2grid((4*rows, 2), pos_resd, rowspan=1, colspan=1)
    #     mod = np.interp(tim, dataset.times_lc['TESS'][Ttimes], mod)
        resid = flx-mod
        FitStats = '$\sigma$ = %d [ppm]'%(np.std(resid)*1e6)
        plt.errorbar(tim, resid, yerr=err, fmt='o', mec='red', ecolor='red', elinewidth=3, mfc='white', ms=4, label=FitStats)
        plt.plot(tim, np.zeros(len(tim)), color='black')
        plt.legend()
        if position[0]%2 == 0: #only plot y-axis on left columns of residuals plots
            plt.ylabel('Residuals')
        if T > len(mid_trans)-3: #only plot x-axis on last two residals
            plt.xlabel("JD date [days]")
    # plt.show()
    plt.savefig(final_path+'/'+Survey[0]+str(ID)+'_transits.png')
    plt.close()


if plot_phase_trans: # Get phases
    plt.figure(figsize=(18,7))
    #To plot data:
    plt.subplot2grid((4, 1), (0,0), rowspan=3, colspan=1)
    X_phased, LC_phased, LC_phasedErr, Pha_ord = utils.PhaseFold(dataset.times_lc['TESS'],dataset.data_lc['TESS'], dataset.errors_lc['TESS'], P[0], PhaseNums=1, t0=t0_cent)
    plt.plot(X_phased, LC_phased, '.', color ='cornflowerblue', alpha = 0.3, zorder=2)
    plt.plot(X_phased,transit_model[Pha_ord], color='black',zorder=10)
    plt.fill_between(X_phased,transit_up68[Pha_ord],transit_low68[Pha_ord],color='red',alpha=0.95,zorder=5)
    paz, flx, err = utils.BinDat(X_phased, LC_phased, LC_phasedErr,Bin=30)
    plt.errorbar(paz, flx, yerr=err, fmt='o', mec='blueviolet', ecolor='blueviolet', elinewidth=3, mfc='white', ms=4,zorder=10)
    sorted_flx = np.sort(dataset.data_lc['TESS']) #to pick the 3rd largest and 3rd smallest flx counts for plot y limit
    plt.ylim((sorted_flx[3],sorted_flx[-3]))
    plt.xlim((-.05,.05))
    #To plot residuals:
    plt.subplot2grid((4, 1), (3,0), rowspan=1, colspan=1)
    paz, bin_modl, zero_err = utils.BinDat(X_phased, transit_model[Pha_ord], np.ones(len(Pha_ord)),Bin=30)
    paz, flx, err = utils.BinDat(X_phased, dataset.data_lc['TESS'][Pha_ord], dataset.errors_lc['TESS'][Pha_ord],Bin=30)
    resid = flx-bin_modl
    FitStats = '$\sigma$ = %d [ppm]'%(np.std(resid)*1e6)
    plt.plot(paz, np.zeros(len(paz)), 'k-')
    plt.errorbar(paz, resid, yerr=err, fmt='o', mec='blueviolet', ecolor='blueviolet', elinewidth=3, mfc='white', ms=4, label=FitStats)
    plt.legend()
    plt.xlim((-.05,.05))
    # plt.show()
    plt.savefig(final_path+'/'+Survey[0]+str(ID)+'_phasefoldedTrans.png')
    plt.close()

# To get RV parameters:
K = np.array(PstSum[np.where(PstSum[:,0] == 'K_p1')[0][0],1:],dtype=np.float)
insturs = list(dataset.times_rv.keys())

#to count how many different plots there should be
Mx_timeGap = 30 #max distance in time one observation can be from another in order to make a new figure
TotalTimes, Mus = np.array([]), []
for inst in insturs:
    TotalTimes = np.concatenate((TotalTimes, dataset.times_rv[inst]))
    Mus.append(np.median(results.posteriors['posterior_samples']['mu_'+inst]))
    print ("For inst: mu="+str(Mus[-1])+", length="+str(len(dataset.times_rv[inst])))

colors = ['cornflowerblue', 'orangered', 'chartreuse', 'slategrey', 'blueviolet', 'darkgoldenrod1'] #color charts: https://www.webucator.com/article/python-color-constants-module/

# if plot_single_RV: #Don't think I need this. If I decide I do, i'll have to edit the code
#     #Now to plot the data with the gaps at the proper locations
#     Diff = np.diff(TotalTimes)
#     Gaps = np.where(Diff > Mx_timeGap)[0]
#     for G in Gaps:
#         if G == Gaps[-1]: #if at last Gap.....?
#         idx_cnt = 0 #to keep track of the global index (index after combining all data)
#         ins_cnt = 0 #to keep track of which insturment we on
#         for ins in insturs:
#             if idx_cnt <= G+1 and G+1 < idx_cnt+len(dataset.times_rv[ins]):
#                 timeRV, RV, RVerr = dataset.times_rv[ins][G+1:],dataset.data_rv[ins][G+1:], dataset.errors_rv[ins][G+1:]
#                 plt.errorbar(timeRV, RV-Mus[ins_cnt], yerr=RVerr, fmt='o', mec=colors[ins_cnt], ecolor=colors[ins_cnt],\
#                  elinewidth=3, mfc='white', ms=7, label=ins, zorder=10)
#             elif G == Gaps[-1]: #if at last Gap.....?
#             idx_cnt += len(dataset.times_rv[ins])
#             ins_cnt += 1
#         plt.legend()
#         plt.show()
#         plt.close()

if plot_RV:
    #just make a model from the last used insturment (assuming inputted in chronilogical order). Shouldn't matter much, since a sine function
    model_times = np.linspace(np.min(dataset.times_rv[insturs[-1]])-1,np.max(dataset.times_rv[insturs[-1]])+1,4000) 
    keplerian = results.rv.evaluate(insturs[-1], t=model_times) - np.mean(Mus)
    idx0 = utils.find_nearest(keplerian,0)
    print ("model_times[idx0]:", model_times[idx0])
    t0_cent = t0[0]+(P[0]/2)#mid_transit+(P[0]/2). So centered on the tranist 

    #To plot the phase folded data
    fig = plt.figure(figsize=(18,6))
    plt.title("Radial Velocity vs. Phase", fontsize=20, fontweight='bold')
    ins_cnt, Max_rv_amp = 0, 0  #to keep track of which insturment we on, and the max RV signal
    for InSt in insturs: #to plot each insturment seperately 
        X_phased, RV_phased, RV_phasedErr, Pha_ord = utils.PhaseFold(dataset.times_rv[InSt],dataset.data_rv[InSt], dataset.errors_rv[InSt], P[0], PhaseNums=1, t0=t0_cent)
        plt.errorbar(X_phased, RV_phased-Mus[ins_cnt], yerr=RV_phasedErr, fmt='o', mec=colors[ins_cnt], ecolor=colors[ins_cnt],\
                     elinewidth=3, mfc='white', ms=7, label=InSt, zorder=10)    
        max_i = np.max(np.abs(RV_phased-Mus[ins_cnt]))
        if Max_rv_amp < max_i:
            Max_rv_amp = max_i
        ins_cnt +=1

    #To plot the best fit lc
    mX_phased = np.linspace(-.5, .5, num=len(model_times))
    B, phase_shift = 2*np.pi, 0 #Period = 2pi/B. Remember Period = 1, since phased
    plt.plot(mX_phased, -K[0]*np.sin(B*(mX_phased+phase_shift)), 'k') #in case you forget how sine functions work:  https://www.mathsisfun.com/algebra/amplitude-period-frequency-phase-shift.html
    plt.plot(np.linspace(-.5, .5, 10), np.ones(10), 'k--', alpha =.5)
    print_results = "K="+str(np.round(K[0],2))+"\u00B1"+str(np.round(np.max(K[1:]),2))+"[m/s] \nP="+str(np.round(P[0],7))+"\u00B1"+str(np.round(np.max(P[1:]),7))+"[days]"
    plt.text(-.45, -Max_rv_amp+10, print_results, fontsize=20)

    plt.legend(fontsize=20)
    plt.ylabel('RV (m/s)', fontsize=20)
    plt.xlabel('Phase', fontsize=20)
    plt.xlim([-.5, .5])
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    # plt.show()
    plt.savefig(final_path+'/'+Survey[0]+str(ID)+'_RV.png')
    plt.close()

####### To estimate the planetary parameters #######
#TODO: add asymmetric uncertainties. To do this, basically do the calculation twice, one for the upper and lower bound
def Mass_P(P, Ms, K, i, aRs, Rs):#to return the mass of the planet in jupiter masses
    #assuming units of: P[days], Ms[M_sun], K[m/s], i[deg], aRs[a/Rs]
    # to calculate mass of planet from RV: https://socratic.org/questions/5929fc3711ef6b28473ff602#:~:text=The%20radial%20velocity%20of%20the,p%3D2%CF%80ap%20.&text=Substituting%20values%20gives%20vp,the%20planet%20mass%20and%20velocity.&text=A%20solar%20mass%20is%201.99%E2%8B%851030%20kg.
    #m_p = (P*Ms*K*sin(i)/2*pi*a)
    days2sec, deg2rad = 86400, np.pi/180 #sec in a day, degrees to rads
    M_jup, M_sun, R_sun = 1.898e27, 1.989e30, 6.957e8 #kg, kg, meter
    aRs_2_a = lambda aRs, Rs: aRs*Rs*R_sun #meter
    num = (P*days2sec)*(Ms*M_sun)*(K*umath.sin(i*deg2rad))
    denom = 2*np.pi*aRs_2_a(aRs,Rs)
    return (num/denom)/M_jup #to convert from SI to jupiter masses

Per, Ms, K1 = ufloat(P[0],np.max(P[1:])), ufloat(Ms[0], Ms[1]), ufloat(K[0],np.max(K[1:]))
i, Rs, AU = ufloat(inc[0],np.max(inc[1:])), ufloat(Rs[0], Rs[1]), 1.496e11 #meters
a_Rs = ufloat(aRs[0],np.max(aRs[1:]))
print ("Per", Per, "\nMs", Ms, "\nK1", K1, "\ni", i, "\nRpRs",ufloat(RpRs[0],np.max(RpRs[1:])))
print ("Rs", Rs, "\na_Rs", a_Rs, "\nb", ufloat(b[0],np.max(b[1:])), "\nrho_*",ufloat(rho_s[0],np.max(rho_s[1:])))
Massp = Mass_P(Per, Ms, K1, i, a_Rs, Rs)
print ("mass of planet:", Massp.n, '\u00B1', Massp.s)
