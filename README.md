# JointParameterEstimation
Routines used to jointly estimate orbital and planetary parameters for an exoplanet with radial velocity and transit observations. As well as jointly using multiple photometric observations to estimate the rotational period of the host star. These observations include out-of-transit TESS photometric observations - which is relatively high cadence but short duration monitoring, and ASAS-SN photometric observations - which have lower cadence and precision but are a conglomeration of data spanning years. The joint probably distributions are sampled using Juliet [https://github.com/nespinoza/juliet](url)

This routine was implemented and utilized in McGruder et. al. 2022 (ADS link: [https://ui.adsabs.harvard.edu/abs/2023ApJ...944L..56M/abstract](url))

Required packages (excluding anaconda standard packages):
1) lightkurve - [https://docs.lightkurve.org/](url)
2) uncertainties - [https://pypi.org/project/uncertainties/](url)
3) juliet - [https://github.com/nespinoza/juliet](url)
4) corner - [https://github.com/dfm/corner.py](url)


Scripts: 
1) utils.py - Internal module that stores routines that are needed for both 'Juliet_JointFit.py' and 'Photo_RunningJuliet2.py.' This is also the same 
2) Juliet_JointFit.py - script used to estimate the system parameters which are the planetary orbital period, stellar mass, RV semi-major amplitude, inclination, planetary radius (relative to stellar radius), Radius of star, semi-major axis of planetary orbit, transit impact parameter, stellar density, and mass of planet.
3) Photo_RunningJuliet2.py - script used to estimate the host star's rotation period utilizing the multiple different photometric data.

This routine was used for many transiting systems, but for an example, the data and outputs are for the star WASP-6 and its known transiting hot Jupiter. 

Input data (included)- Juliet_JointFit.py:
1) "W6_fullRV_sysCorr_noRM.dat" RV data, with the exclusion of the RV data when the planet transits since the code cannot yet properly model the Rossiter–McLaughlin effect.
2) "W6_TESS.npy" is the TESS data coming from the TESS reduction pipeline.
3) Additional inputs must be written in the code directly, in the section titled "Target input info," for specific settings desired.

Output (included)- Juliet_JointFit.py: 
1) "W6_jointfit_fixE" (the output folder name is specified in the script "Target input info") is a folder containing the posterior samples and results, Corner plots, and best fit RV and transit models overplotted on their corresponding datasets. 
2) The system parameters, which are printed out by the script.

Input data (included)- Photo_RunningJuliet2.py:
1) "WASP6_ASASSNg_binnedClipped.npy" is the binned ASAS-SN data with the g-band observations.
2) "WASP6_ASASSNv_binnedClipped.npy" is the binned ASAS-SN data with the v-band observations. The g and v band observations are analyzed separately because of different systematics and sensitivities, which need to be modeled differently. 
3) "WASP6_TESS_binnedOOT.npy" is reduced TESS data, which is also binned, and all data during the WASP-6 transit removed.
4) Additional inputs must be written in the code directly, in the section titled "INPUT," for specific settings desired.

Output (included)- Photo_RunningJuliet2.py:
"binnedTESSnASASSNpass1" (the output folder name is specified in the script "INPUT") is a folder containing the posterior samples and results, Corner plots, and best fit models overplotted on their corresponding datasets.
   Note: The folder 'Juliet' is created to store all photometric monitoring data, and a subfolder with the name of the target (which is an "INPUT" parameter) is created as the code runs.
