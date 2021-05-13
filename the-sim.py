#!/Applications/anaconda/bin/python

# 2016-03-03 refining model, energy dependence in Thomas-Imel param alpha
# 2021-05-11 simplifying code to just provide LUX / LUX-like sim, and
#            adding crystalline case
#             

import sys
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from scipy.optimize import minimize
from scipy.interpolate import interp1d
from matplotlib import colors
from matplotlib import rc # import matplotlib.colors as colors

rc('text', usetex=True)
rc('font',**{'family':'serif','serif':['Palatino']})
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

############################################# my imports
#import myfunctions as mf
#mf.np = np
def gaussian(x,xx):
    y = x[0]*np.exp(-np.power((xx - x[1])/x[2], 2.)/2) /(np.sqrt(2.*np.pi)*x[2])
    return y
def gaussian_min(x,xx,counts):
	return np.sum((gaussian(x,xx)-counts)**2)

############################################# define stuff

def xxfunction(particle,Ed,er):
	if particle=='electron':
		if Ed==180: # tailored to LUX
			xx = 0.032+ 0.004/(1+np.exp((er-5)/0.5))
		elif Ed==530: # tailored to XENON100
			xx = 0.036*(np.exp(-(er-3.5)/12))
			xx[er<3.5] = xx[er==3.5]			
			xx[er>=17] = xx[er==17]			
		elif Ed==730: # tailored to XENON10
			xx = 0.034*np.exp(-(er-2)/13)
			xx[er<2]= xx[er==2]
			xx[er>=17] = xx[er==17]			
		elif Ed==3900: # tailored to ZEPLINIII
			xx = 0.025*(np.exp(-(er-0)/14))
			xx[er>=17] = xx[er==17]
		else:
			xx = x[2]*np.ones(len(er))			
	elif particle=='neutron':
		x2 = 0.044
		if Ed==180: # tailored to LUX
			xx = x2+ 0.002*(1/(1+np.exp((er-15)/2))-1)
		elif Ed==530: # tailored to XENON100
			xx = x2+ 0.005*(1/(1+np.exp((er-15)/2))-1)
		elif Ed==730: # tailored to XENON10
			xx = x2+ 0.005*(1/(1+np.exp((er-15)/2))-1)
		elif Ed==3900: # tailored to ZEPLINIII
			xx = x2+ 0.005*(1/(1+np.exp((er-15)/2))-1)
			#xx = (x2-0.015) # ad-hoc to address suspected low-E multiple scatter contamination
		else:
			xx = x2*np.ones(len(er))			
	return xx

def signalYields(particle,Ed,er,detector):
	wq = 0.0138 # average energy required to create a single quanta, keV
	if particle == 'electron':
		x0 = np.array([1,0.06])
		fn = 1 # quenching is defined relative to electron recoils
		NexOverNi = np.abs(x0[1]) 
		alphaOveraSquaredv = xxfunction(particle,Ed,er)
		Nt = er * fn / wq # total quanta created
		Ni = Nt/(1+NexOverNi) # ions created
		if (detector=='crystal'): # more e- escape recombination due to e- mobility x2 increase
			Ni = Ni*1.15
		xi = Ni/4 * alphaOveraSquaredv # Thomas-Imel parameter
		#r = 1 - np.log(1+xi)/xi # everybody else uses this so define it :)
		f0 = np.log(1+xi)/xi # Thomas-Imel fraction of ions that escape recombination
	elif particle == 'neutron':
		x0 = np.array([0.166,1.1])
		Z = 54;A = 131.3
		epsilon = 11.5 * er / Z**(7/3) # Lindhard reduced energy
		#k = 0.133 * Z**(2/3)*A**(-1/2) # Lindhard approx electronic stopping
		g = 3*epsilon**0.15 + 0.7*epsilon**0.6 + epsilon # Lindhard empirical spline function
		k = np.abs(x0[0])
		NexOverNi = np.abs(x0[1])
		alphaOveraSquaredv = xxfunction(particle,Ed,er)
		fn = k * g / (1 + k*g)
		Nt = er * fn / wq # total quanta created
		Ni = Nt/(1+NexOverNi)
		xi = Ni/4 * alphaOveraSquaredv 
		f0 = np.log(1+xi)/xi
	return (Nt,Ni,f0,x0[1])

def pmtresponse(sigma,samples,softwareThreshold_phe,double_phe_fraction):
	# note that since mu == 1, sigma and resolution (sigma/mu) are interchangeable here
	n_single_phe = np.sum( np.random.uniform(0,1,samples)>double_phe_fraction ) 
	n_double_phe = samples - n_single_phe
	a_single = np.random.normal(1,sigma,n_single_phe)
	a_double = np.random.normal(2,np.sqrt(2)*sigma,n_double_phe)
	a = np.concatenate((a_single,a_double),axis=0)
	mask = a>softwareThreshold_phe #0
	a = a*mask 
	return a

def bandsim(particle,Ed,er,dp,pmtps1,pmtps2,jj,detector='None'):
	ii = len(er)
	Ne = np.zeros((ii,jj)); S1 = np.zeros((ii,jj)); S2 = np.zeros((ii,jj)); Ng = np.zeros((ii,jj)); Nee = np.zeros((ii,jj));
	Nt = np.zeros(ii); Ni = np.zeros(ii); f0 = np.zeros(ii);
	(Nt,Ni,f0,x0) = signalYields(particle,Ed,er,detector)
	for i in range(1,ii): # energies
		for j in range(0,jj): # events
			if particle=="electron":
				tauS1_ns = 27 # ns
				if f0[i]>0 and f0[i]<=1:					
					p=-1
					while p<=0 or p>1:
						if er[i]<30:							
							p = np.random.normal(f0[i],(x0*(1-np.exp(-(er[i])/1))) )
#  							p = np.random.normal(f0[i],x0)
					ne0 = np.random.binomial((Ni[i]),p) 
				else:
					ne0 = 0
			elif particle=="neutron":
				tauS1_ns = 4 # ns
				ne0 = np.random.binomial(Ni[i],f0[i])
				ne1 = 0
			Ne[i][j] = ne0 #+ ne1
			
			if 0: # DO NOT USE. This is not quite correct, has been used in the past, and is copied here for reference only	
				S2[i][j] = eta_extraction * dp[1] * (np.sqrt(Ne[i][j]) * np.random.normal(0,dp[4],1)+Ne[i][j] )
			else: # 2 Feb 2016 +
				Nee[i][j] = np.random.binomial(Ne[i][j],dp[11]) # Number of electrons extracted
				hitpattern2 = np.random.binomial(dp[3]*Nee[i][j], pmtps2*dp[2])
				S2[i][j] = np.sum( pmtresponse(dp[10],np.sum(hitpattern2),dp[9],dp[7]) )
							
			Ng[i][j] = (Nt[i] - Ne[i][j])
			if Ng[i][j]<=0:
				S1[i][j]=0
			else:
				hitpattern = np.random.binomial(Ng[i][j], pmtps1*dp[0]) # distribute Ng photons onto the PMTs, with specified probability to create a phe
				candidateS1 = pmtresponse(dp[10],np.sum(hitpattern),softwareThreshold_phe,dp[7]) #
				if np.sum(candidateS1>0)>1: # apply decay time logic
					phe_times = np.sort(np.random.exponential(tauS1_ns,np.sum(candidateS1>0) ))
					d_times = np.diff(phe_times)
					min_dt = np.min(d_times)
					in_window = (min_dt<dp[8])
				else:
					in_window = True # this can be left as True bc subsequent logic will fail if no candidate phe
				
				if (np.sum(candidateS1>0)>=dp[6]) and in_window:
					S1[i][j] = np.sum( candidateS1 )
				else:
					S1[i][j] = 0
	S1[S1==0]=1e-6; S2[S2==0]=1e-6; #S1[S1==0] = np.finfo(float).eps
	S2[S2<dp[5]] = 0
	if 1:#(detector=='LUX'):
		S2 = S2/2 # to compare with LUX data assume only using bottom array so half the S2
	
	yy=np.log10(S2/S1); yy[np.isnan(yy)]=0; yy[np.isinf(yy)]=0; yy[yy>5]=0;
	return (Ne,Ng,S1,S2,yy)


def getBandMuSigma(xe_s1,xe_be,S1,yy,detector):
	sim_mu = np.zeros((len(xe_s1),1))
	sim_1s = np.zeros((len(xe_s1),1))
	sim_bc = np.zeros((len(xe_s1),1))
	showplot=0
	if particle=="electron" and showplot:
		plt.rcParams["figure.figsize"] = (12,8)
		fig=plt.figure(16);plt.clf()
	for i in range(0,len(xe_be)-1):
		mask0 = np.logical_and(S1>xe_be[i],S1<xe_be[i+1])
		mask1 = np.logical_and(mask0,yy<5)
		mask = np.logical_and(mask1,yy>0.5)
		sim_bc[i] = np.sum(mask)
		if (sim_bc[i]>1):
			if 0: #detector=="XENON100":
				sim_mu[i] = np.mean(yy[mask])
				sim_1s[i] = np.std(yy[mask])
			else:			
				counts,bine=np.histogram(yy[mask],np.arange(1.,3.5,0.03));binc = (bine[0:-1] + bine[1:])/2
				b0 = np.array([np.sum(counts)/10,np.mean(yy[mask]),1*np.std(yy[mask])]) #b0 = np.array([15,2.37,0.142])
				#b0 = [100,2,0.1]
				res = ( minimize(gaussian_min,b0,args=(binc,counts),method='nelder-mead',
								options={'xtol': 1e-8, 'disp': False}) )
				sim_mu[i] = res.x[1]
				sim_1s[i] = res.x[2]
			if particle=="electron" and showplot:# and detector!="XENON100":
				fig=plt.figure(16)
				ax = plt.subplot(3,5,i+1); plt.cla()
				plt.hist(yy[mask],bine,edgecolor="None",facecolor='darkslateblue',alpha=0.5)
				x=np.arange(0.8,3.5,0.01)
				#plt.plot(x,gaussian(b0,x),'g-')
				if (detector=='LUX'):
					plt.plot(x,gaussian([res.x[0],xermu[i],xer1s[i]],x),'k-',color='black',linewidth=1)
				plt.plot(x,gaussian(res.x,x),'--',color='darkslateblue',linewidth=0.75)

				plt.yscale('log');ax.set_ylim(1,300)
				plt.text(1.0,80,"$S1: [%.1f , %.1f]$" % (xe_be[i],xe_be[i+1]) ,fontsize=8)
				plt.text(1.0,150,"$\sigma/\mu = %.3f/%.2f$" % (res.x[2],res.x[1]) ,fontsize=10)
				plt.show(); plt.draw()
				plt.savefig("figs/fits0.png")

	return (sim_mu,sim_1s,sim_bc)

def getLeakage(xe_s1,sim_er_mu,sim_er_1s,sim_nr_mu):
	leakage_frac = np.zeros(len(xe_s1))
	uncert_frac = np.zeros(len(xe_s1))
	y = np.arange(0.5,4.5,1e-4)
	for i in range(0,len(xe_s1)):
		g = gaussian([1,sim_er_mu[i],sim_er_1s[i]],y)
		cut = y<sim_nr_mu[i]
		integral = np.sum(g)
		leakage = np.sum(g[cut])
		leakage_frac[i] = leakage/integral
		uncert_frac[i] = np.sqrt(leakage)/integral
	return (leakage_frac, uncert_frac)


############################################# do stuff


if 1: # plot bands
	if 1:
		if 1:
			detector = 'LUX'
			g1 = 0.115
			eta_extraction = 0.50
		if 0:
			detector = 'LUX-like'
			g1 = 0.115
			eta_extraction = 0.95
		if 1:
			detector = 'crystal'
			g1 = 0.115
			eta_extraction = 0.50
		
		Ed = 180
		double_phe_fraction = 0.2
		npmttop = 61
		npmtbot = 61
		resolution_pmt = 0.37
		g2prime = 0.15 # actual geometric * QE light collection
		gasgain = 160
		g2 = g2prime*gasgain
		resolution_e = np.sqrt(g2)/g2
		nco = 3
		softwareThreshold_phe = 1/3.
		S1window_ns = 100 # ns
		S2Threshold_phe = g2*8.1 # based on 200/24.6

		pmtptop = np.ones(npmttop); pmtptop.fill(0.20/npmttop)
		pmtpbot = np.ones(npmtbot); pmtpbot.fill(0.80/npmtbot)
		pmtps1 = np.concatenate((pmtptop,pmtpbot),axis=0)

		# simple 50:50 assumption for S2
		pmtps2 = np.ones(npmttop+npmtbot); pmtps2.fill(1.0/(npmttop+npmtbot))

		# assign params to array for convenience
		dp = np.zeros(12)
		dp[0] = g1
		dp[1] = g2 # not used, = g2prime*gasgain
		dp[2] = g2prime
		dp[3] = gasgain
		dp[4] = resolution_e # not used, is defined via g2prime and gasgain
		dp[5] = S2Threshold_phe 
		dp[6] = nco
		dp[7] = double_phe_fraction
		dp[8] = S1window_ns
		dp[9] = softwareThreshold_phe
		dp[10] = resolution_pmt
		dp[11] = eta_extraction

		# LUX-specific params digitized from papers
		xe_be = np.append([np.arange(2,4.1,1)],[np.arange(5,55,5)]) # 80
		xe_s1 = (xe_be[1:]+xe_be[0:-1])/2 #np.append([np.arange(2.5,5,1)],[np.arange(7.5,50,5)]) # 75
		xermu = np.array([2.284,2.175,2.097,1.951,1.814,1.731,1.670,1.621,1.585,1.552,1.524,1.500])
		xer128s = np.array([2.109,2.009,1.936,1.796,1.667,1.586,1.529,1.482,1.448,1.417,1.393,1.371])
		xer1s = (xermu-xer128s)/1.28
		xer3s = xermu-xer1s*3
		xnrmu = np.array([1.830,1.720,1.650,1.549,1.481,1.429,1.383,1.349,1.311,1.287,1.262,1.236])
		xnr128s = np.array([1.637,1.525,1.442,1.357,1.330,1.292,1.266,1.237,1.206,1.188,1.163,1.138])
		xnr1s = (xnrmu-xnr128s)/1.28


## plot S2/S1 vs S1
	plt.rcParams["figure.figsize"] = (12,5)
	fig=plt.figure(9); plt.clf();#  axes=fig.add_axes([0.15, 0.15, 0.8, 0.8])
	if 1:
		jj = 500 # number of simulated events per energy
		particle='electron'
		er = np.arange(0.75,12,0.25) # keV
		(Ne,Ng,S1,S2,yy) = bandsim(particle,Ed,er,dp,pmtps1,pmtps2,jj,detector)
		(sim_er_mu,sim_er_1s,sim_bc) = getBandMuSigma(xe_s1,xe_be,S1,yy,detector)
		plt.figure(9); plt.subplot(1,2,1); plt.cla()
		plt.plot(S1,yy,'k.',color='darkslateblue',markeredgecolor="None",alpha=0.13) 
		plt.plot(xe_s1,sim_er_mu,'ko',ms=4,markeredgewidth=0.75,markerFaceColor='None')
		plt.plot(xe_s1,sim_er_mu-sim_er_1s,'k+',ms=6,markeredgewidth=0.75)
		plt.plot(xe_s1,sim_er_mu+sim_er_1s,'k+',ms=6,markeredgewidth=0.75)

		if detector=='LUX':		# digitized LUX data to compare
			plt.plot(xe_s1,xermu,'k.-',ms=5,markerfacecolor="None",alpha=1.,lineWidth=0.5)
			plt.plot(xe_s1,xermu-xer1s,'k.--',ms=5,markerfacecolor="None",alpha=1.,lineWidth=0.5)
			plt.plot(xe_s1,xermu+xer1s,'k.--',ms=5,markerfacecolor="None",alpha=1.,lineWidth=0.5)

		plt.xlabel('S1',fontsize=16)
		plt.ylabel('$log_{10}$(S2/S1)',fontsize=16)
		#plt.grid(True)
		plt.xticks(np.arange(0,51,10),fontsize=16)
		plt.yticks(np.arange(0.5,3.5,0.5),fontsize=16)
		plt.axis([0,50,0.4,3.3])
		plt.minorticks_on()
		plt.show(); plt.draw()

	if 1:
		jj = 100 # number of simulated events per energy
		particle='neutron'
		er = np.arange(0.75,35,0.25) # keV
		(Ne,Ng,S1,S2,yy) = bandsim(particle,Ed,er,dp,pmtps1,pmtps2,jj,detector)
		(sim_nr_mu,sim_nr_1s,sim_bc) = getBandMuSigma(xe_s1,xe_be,S1,yy,detector)
		plt.figure(9); plt.subplot(1,2,2); plt.cla()
		plt.plot(S1,yy,'k.',color='chocolate',markeredgecolor="None",alpha=0.13) 
		plt.plot(xe_s1,sim_nr_mu,'ko',ms=4,markeredgewidth=0.75,markerFaceColor='None')
		plt.plot(xe_s1,sim_nr_mu-sim_nr_1s,'k+',ms=6,markeredgewidth=0.75)
		plt.plot(xe_s1,sim_nr_mu+sim_nr_1s,'k+',ms=6,markeredgewidth=0.75)

		if detector=='LUX':		# digitized LUX data to compare
			plt.plot(xe_s1,xnrmu,'k.-',ms=4,markerfacecolor="None",alpha=1.,lineWidth=0.5)
			plt.plot(xe_s1,xnrmu-xnr1s,'k.--',ms=6,markerfacecolor="None",alpha=1.,lineWidth=0.5)
			plt.plot(xe_s1,xnrmu+xnr1s,'k.--',ms=6,markerfacecolor="None",alpha=1.,lineWidth=0.5)

		if 0:
			energies = [1,2,5,10,15,20]
			s1s = np.arange(0.01,50,0.01)
			for i in range(len(energies)):
				s2s = (energies[i]/0.0138 - s1s/g1)*g2
				yy = np.log10((s2s)/s1s)
				plt.plot(s1s,yy,'k:')
	
		plt.xlabel('S1',fontsize=16)
		plt.ylabel('$log_{10}$(S2/S1)',fontsize=16)
		#plt.grid(True)
		plt.xticks(np.arange(0,51,10),fontsize=16)
		plt.yticks(np.arange(0.5,3.5,0.5),fontsize=16)
		plt.axis([0,50,0.4,3.3])
		plt.minorticks_on()
		plt.show(); plt.draw()


		plt.savefig("figs/bands.png")

	if 1:
		plt.rcParams["figure.figsize"] = (7,5)
		fig=plt.figure(10); plt.clf();  axes=fig.add_axes([0.15, 0.15, 0.75, 0.8])
		(leakage_frac,uncert_frac) = getLeakage(xe_s1,sim_er_mu,sim_er_1s,sim_nr_mu)
		plt.errorbar(xe_s1,leakage_frac,yerr=uncert_frac,xerr=(xe_s1-xe_be[0:-1]), fmt='o',color='darkslateblue',markeredgecolor='darkslateblue',markerfacecolor="white",label='simulation')
		plt.yscale('log')
		if 1: # plot LUX data for comparison
			discx = np.arange(1.5,49.5+1,1); discy = 1e-4*np.array([48.4,44.5,11.7,2.5,0,2.84,2.83,5.52,8.28,2.83,8.26,0,8.91,5.78,17.8,14.9,23.9,21.3,12.5,50.7,9.50,3.22,20.1,45.1,34.9,29.1,31.7,36.5,14.6,15.9,22.7,19.8,15.3,20.5,12.7,21.6,31.2,22.4,32.9,9.17,4.97,10.1,39.8,37.1,10.8,33.9,0,5.90,12.4])
			plt.plot(discx,discy,'ks',color='gray',markerfacecolor="None",label='LUX data')
		
		ytik = [1e-6,1e-5,1e-4,1e-3,1e-2]
		xx = 51
		for ii in range(2,len(ytik)):
			plt.text(xx,ytik[ii]*0.90,100*(1-ytik[ii]),fontsize=16)
		plt.text(52,1.0e-6,'discrimination (\%)',rotation=90,fontsize=14,verticalalignment='bottom',horizontalalignment='left')
				
		plt.plot(np.array([0,50]),np.array([1,1])*5e-3,'k--')
		plt.axis([0,50,1e-6,3e-2])
		plt.xlabel('detected scintillation photons',fontsize=16)
		plt.ylabel('leakage fraction',fontsize=16)
		plt.xticks(np.arange(0,51,10),fontsize=16)
		plt.yticks(ytik,fontsize=16)
		plt.grid(True)
		plt.legend(loc='lower right')
		plt.show(); plt.draw()
		plt.savefig("figs/discrim0.png")
	
		# to compare
		try:
			tmp = np.load('tmp.npy')
		except:
			print('no file to load')
			
		comp_leakage = tmp/leakage_frac
		print(comp_leakage)

		tmp = leakage_frac; np.save('tmp',tmp)

