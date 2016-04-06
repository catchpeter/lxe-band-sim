#!/Applications/anaconda/bin/python

# 2016-03-03 refining model, energy dependence in Thomas-Imel param alpha

import sys
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from scipy.optimize import minimize
from scipy.interpolate import interp1d
from matplotlib import colors
from matplotlib import rc # import matplotlib.colors as colors

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})

#plt.ioff()

############################################# my imports
import myfunctions as mf
mf.np = np

############################################# define stuff

def xxfunction(particle,Ed,er):
	if particle=='electron':
		if Ed==180:
			xx = 0.032+ 0.004/(1+np.exp((er-5)/0.5))
		elif Ed==530:
			xx = 0.036*(np.exp(-(er-3.5)/12))
			xx[er<3.5] = xx[er==3.5]			
			xx[er>=17] = xx[er==17]			
		elif Ed==730:
			xx = 0.034*np.exp(-(er-2)/13)
			xx[er<2]= xx[er==2]
  			xx[er>=17] = xx[er==17]			
		elif Ed==3900:
			xx = 0.025*(np.exp(-(er-0)/14))
			xx[er>=17] = xx[er==17]
		else:
			xx = x[2]*np.ones(len(er))			
	elif particle=='neutron':
		x2 = 0.044
		if Ed==180:
			xx = x2+ 0.002*(1/(1+np.exp((er-15)/2))-1)
		elif Ed==530:
			xx = x2+ 0.005*(1/(1+np.exp((er-15)/2))-1)
		elif Ed==730:
			xx = x2+ 0.005*(1/(1+np.exp((er-15)/2))-1)
		elif Ed==3900:
			xx = x2+ 0.005*(1/(1+np.exp((er-15)/2))-1)
			#xx = (x2-0.015) # ad-hoc to address suspected low-E multiple scatter contamination
		else:
			xx = x2*np.ones(len(er))			
	return xx

def signalYields(particle,Ed,er):
	wq = 0.0138 # average energy required to create a single quanta, keV
	if particle == 'electron':
		x0 = np.array([1,0.06])
		fn = 1 # quenching is defined relative to electron recoils
		NexOverNi = np.abs(x0[1]) 
		alphaOveraSquaredv = xxfunction(particle,Ed,er)
		Nt = er * fn / wq # total quanta created
		Ni = Nt/(1+NexOverNi) # ions created
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
	return (Nt,Ni,f0)

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

# dp[0] = g1
# dp[1] = g2 # not used, = g2prime*gasgain
# dp[2] = g2prime
# dp[3] = gasgain
# dp[4] = resolution_e
# dp[5] = S2Threshold_phe
# dp[6] = nco
# dp[7] = double_phe_fraction
# dp[8] = S1window_ns
# dp[9] = resolution_pmt


def bandsim(particle,Ed,er,dp,pmtps1,pmtps2,jj):
	ii = len(er)
	if Ed==180: # LUX electric field
		x0 = 0.06
	elif Ed==530: # XENON100 electric field
		x0 = 0.07
	elif Ed==730: # XENON10 electric field
		x0 = 0.07		
	Ne = np.zeros((ii,jj)); S1 = np.zeros((ii,jj)); S2 = np.zeros((ii,jj)); Ng = np.zeros((ii,jj)); Nee = np.zeros((ii,jj));
	Nt = np.zeros(ii); Ni = np.zeros(ii); f0 = np.zeros(ii);
	(Nt,Ni,f0) = signalYields(particle,Ed,er)
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
	yy=np.log10(S2/S1); yy[np.isnan(yy)]=0; yy[np.isinf(yy)]=0; yy[yy>5]=0;
	return (Ne,Ng,S1,S2,yy)


def getBandMuSigma(xe_s1,xe_be,S1,yy,detector):
	sim_mu = np.zeros((len(xe_s1),1))
	sim_1s = np.zeros((len(xe_s1),1))
	sim_bc = np.zeros((len(xe_s1),1))
	showplot=0
	if particle=="electron" and showplot:
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
				res = ( minimize(mf.gaussian_min,b0,args=(binc,counts),method='nelder-mead',
								options={'xtol': 1e-8, 'disp': False}) )
				sim_mu[i] = res.x[1]
				sim_1s[i] = res.x[2]
			if particle=="electron" and showplot:# and detector!="XENON100":
				fig=plt.figure(16)
				ax = plt.subplot(3,5,i+1); plt.cla()
				plt.hist(yy[mask],bine,edgecolor="0.2",facecolor='None',alpha=0.1)
				x=np.arange(0.8,3.5,0.01)
				#plt.plot(x,mf.gaussian(b0,x),'g-')
				plt.plot(x,mf.gaussian(res.x,x),'k-.')
				plt.plot(x,mf.gaussian([res.x[0],xermu[i],xer1s[i]],x),'b-')
				plt.yscale('log');ax.set_ylim(1,300)
				plt.text(0.6,50,"$S1: [%.1f , %.1f]$" % (xe_be[i],xe_be[i+1]) ,fontsize=8)
				plt.text(0.6,100,"$\mu = %.2f$ \n $\sigma = %.3f$" % (res.x[1],res.x[2]) ,fontsize=10)
				plt.show(0); plt.draw()
	return (sim_mu,sim_1s,sim_bc)

def getLeakage(xe_s1,sim_er_mu,sim_er_1s,sim_nr_mu):
	leakage_frac = np.zeros(len(xe_s1))
	y = np.arange(0.5,4.5,1e-4)
	for i in range(0,len(xe_s1)):
		g = mf.gaussian([1,sim_er_mu[i],sim_er_1s[i]],y)
		cut = y<sim_nr_mu[i]
		integral = np.sum(g)
		leakage = np.sum(g[cut])
		leakage_frac[i] = leakage/integral
	return leakage_frac
		
############################################# do stuff


if 1: # plot bands
	if 1:
		detector = 'LUX'
		Ed = 180
		#Ed = 1
		er = np.arange(0.75,12,0.25) # keV
	elif 0:
		detector = 'XENON100'
		Ed = 530
		er = np.arange(0.75,25,0.25) # keV
	elif 0:
		detector = 'XENON10'
		Ed = 730
		er = np.arange(0.75,20,0.25) # keV
	elif 0:
		detector = 'ZEPLIN3'
		Ed = 3900
		er = np.arange(0.75,30,0.25) # keV

	execfile("define-detector.py")

## plot S2/S1 vs S1
	if 0:
		jj = 400 # number of simulated events per energy
		particle='electron'
		#x0 = define_x0(particle,Ed)
		(Ne,Ng,S1,S2,yy) = bandsim(particle,Ed,er,dp,pmtps1,pmtps2,jj)
		(sim_er_mu,sim_er_1s,sim_bc) = getBandMuSigma(xe_s1,xe_be,S1,yy,detector)
		plt.figure(9); ax = plt.subplot(1,2,1); plt.cla()
		plt.plot(S1,yy,'k.',color='darkslateblue',markeredgecolor="None",alpha=0.13) 
		execfile("add-to-plot.py")
		plt.plot(xe_s1,sim_er_mu,'k+',ms=6,markerfacecolor="None",markeredgewidth=0.75)
		plt.plot(xe_s1,sim_er_mu-sim_er_1s,'k+',ms=6,markerfacecolor="None",markeredgewidth=0.75)
# 		for bb in range(0,len(xe_s1)):
# 			plt.plot([xe_be[bb], xe_be[bb+1]],[sim_er_mu[bb], sim_er_mu[bb]],'k-')
# 			#plt.plot([xe_s1[bb], xe_s1[bb]],[sim_er_1s[bb], sim_er_1s[bb]]/np.sqrt(sim_bc[bb]),'k-')
# 			plt.plot([xe_be[bb], xe_be[bb+1]],[sim_er_mu[bb]-sim_er_1s[bb], sim_er_mu[bb]-sim_er_1s[bb]],'k-')

	if 0:
		jj = 100 # number of simulated events per energy
		particle='neutron'
		er = np.arange(0.75,35,0.25) # keV
		#x0 = define_x0(particle,Ed)
		(Ne,Ng,S1,S2,yy) = bandsim(particle,Ed,er,dp,pmtps1,pmtps2,jj)
		(sim_nr_mu,sim_nr_1s,sim_bc) = getBandMuSigma(xe_s1,xe_be,S1,yy,detector)
		plt.figure(9);ax = plt.subplot(1,2,2); plt.cla()
		plt.plot(S1,yy,'k.',color='chocolate',markeredgecolor="None",alpha=0.13) 
		execfile("add-to-plot.py")
		plt.plot(xe_s1,sim_nr_mu,'k+',ms=6,markerfacecolor="None",markeredgewidth=0.75)
		plt.plot(xe_s1,sim_nr_mu-sim_nr_1s,'k+',ms=6,markerfacecolor="None",markeredgewidth=0.75)
# 		for bb in range(0,len(xe_s1)):
# 			plt.plot([xe_be[bb], xe_be[bb+1]],[sim_nr_mu[bb], sim_nr_mu[bb]],'k-')
# 			#plt.plot([xe_s1[bb], xe_s1[bb]],[sim_nr_1s[bb], sim_nr_1s[bb]]/np.sqrt(sim_bc[bb]),'k-')
# 			plt.plot([xe_be[bb], xe_be[bb+1]],[sim_nr_mu[bb]-sim_nr_1s[bb], sim_nr_mu[bb]-sim_nr_1s[bb]],'k-')

		plt.savefig("figs/bands.pdf")

	if 1:
		jj=400
		execfile("get-leakage.py")	
		f5=plt.figure(10); plt.clf();
		plt.semilogy(xe_s1,leakage_frac,'r+',markersize=10) #,color='firebrick'
		if detector == 'LUX':
			discx = np.arange(1.5,49.5+1,1); discy = 1e-4*np.array([48.4,44.5,11.7,2.5,0,2.84,2.83,5.52,8.28,2.83,8.26,0,8.91,5.78,17.8,14.9,23.9,21.3,12.5,50.7,9.50,3.22,20.1,45.1,34.9,29.1,31.7,36.5,14.6,15.9,22.7,19.8,15.3,20.5,12.7,21.6,31.2,22.4,32.9,9.17,4.97,10.1,39.8,37.1,10.8,33.9,0,5.90,12.4])
		elif detector == 'XENON10':
			discx = np.array([5.5,7.7,9.9,12.1,15.4,19.8,24.2]); discy = 1e-4*np.array([7.9,17.4,10.9,41.3,43.1,43.6,72.1])
		try:
			plt.plot(discx,discy,'ks',color='gray')
		except:
			print "* oops, no data to plot"

		if 1: # for proposal
			execfile("define-detector.py")
			g1 = 0.116*2; dp[0]=g1
			eta_extraction = 0.95; dp[11] = eta_extraction
			execfile("get-leakage.py")
			plt.semilogy(xe_s1,leakage_frac,'go',markersize=10)

		plt.plot(45,4e-6,'b*',markersize=10,markeredgecolor='b')
		plt.plot([45,45],[2e-6 , 8e-6],'b-')
		plt.text(33,4e-6,'typical uncertainty:',verticalalignment='center',color='blue')
		
		ytik = [1e-6,1e-5,1e-4,1e-3,1e-2]
		xx = 51
		for ii in range(2,len(ytik)):
			plt.text(xx,ytik[ii]*0.90,100*(1-ytik[ii]),fontsize=16)
		plt.text(52,1.5e-6,'discrimination (\%):',rotation=90,fontsize=14,verticalalignment='bottom',horizontalalignment='left')
				
		plt.plot(np.array([0,50]),np.array([1,1])*5e-3,'k--')
		plt.axis([0,50,1e-6,3e-2])
		plt.xlabel('detected scintillation photons',fontsize=18)
		plt.ylabel('leakage fraction',fontsize=18)
		plt.xticks(np.arange(0,51,10),fontsize=18)
		plt.yticks(ytik,fontsize=18)
		plt.grid(True)
		plt.show(0); plt.draw()
		plt.savefig("figs/discrim0.pdf")

if 0: # plot recombination energy dependence
	f10=plt.figure(11); plt.clf();
# 	axt = plt.subplot(2,1,1); plt.cla()

	particle='electron'; Ed=180
	xx=xxfunction(particle,Ed,er)
	e180, = plt.plot(er[er<10],xx[er<10],'k-',color='darkslateblue',label='electron 180 V/cm')

	particle='electron'; Ed=730
	xx=xxfunction(particle,Ed,er)
	e730, = plt.plot(er[er<15],xx[er<15],'k--',color='darkslateblue',label='electron 730 V/cm')

	particle='electron'; Ed=530
	xx=xxfunction(particle,Ed,er)
	e530, = plt.plot(er[er<17],xx[er<17],'k-.',color='darkslateblue',label='electron 530 V/cm')

# 		particle='electron'; Ed=3900
# 		x0 = define_x0(particle,Ed)
# 		xx=xxfunction(particle,Ed,er)
# 		e3900, = plt.plot(er,xx,'k:',color='darkslateblue',label='electron 3900 V/cm')

# 	plt.ylabel(r'$\alpha / a^2 v$',fontsize=18)
# 	plt.xticks(np.arange(0,51,5),fontsize=18)
# 	plt.yticks(np.arange(0,0.07,.01),fontsize=18)
# 	plt.minorticks_on()
# 	plt.grid(True)
# 	plt.axis([0,20,0,0.05])
# 
# 	axb = plt.subplot(2,1,2); plt.cla()

	particle='neutron'; Ed=180
	xx=xxfunction(particle,Ed,er)
	n180, = plt.plot(er,xx,'k-',color='chocolate',label='neutron 180 V/cm')

	particle='neutron'; Ed=730
	xx=xxfunction(particle,Ed,er)
	n730, = plt.plot(er,xx,'k--',color='chocolate',label='neutron 730 V/cm')

	plt.xlabel('Energy~/~keV',fontsize=18)
	plt.ylabel(r'$\alpha / a^2 v$',fontsize=18)
	plt.xticks(np.arange(0,51,5),fontsize=18)
	plt.yticks(np.arange(0,0.07,.01),fontsize=18)
	plt.minorticks_on()
	#plt.grid(True)
	plt.axis([0.5,50,0.01,0.05])
	plt.xscale('log');
	#plt.yscale('log');

	plt.legend(handles=[n180,n730,e180,e530,e730],loc='lower left')
	plt.show(0); plt.draw()
	plt.savefig("figs/yields.pdf")

if 0: # plot discrimination versus x
	f4=plt.figure(4); plt.clf(); #ax7=f7.add_axes([0.1, 0.1, 0.8, 0.8]); plt.clf();# plt.hold(True)
	detector = 'LUX'

	Ed = 180
	execfile("define-detector.py")
	g1 = 0.08; dp[0]=g1 #0.116*1.5
	execfile("get-leakage.py")
	plt.semilogy(xe_s1,leakage_frac,'rx')

	Ed = 180	
	execfile("define-detector.py")
	execfile("get-leakage.py")
	lux180g11, = plt.semilogy(xe_s1,leakage_frac,'ko')

	Ed = 180
	execfile("define-detector.py")
	eta_extraction = 0.95; dp[11] = eta_extraction
	execfile("get-leakage.py")
	lux180g11eta95, = plt.semilogy(xe_s1,leakage_frac,'k+')

	Ed = 730
	execfile("define-detector.py")
	execfile("get-leakage.py")
	lux730g11, = plt.semilogy(xe_s1,leakage_frac,'b^')


	plt.plot(np.array([0,50]),np.array([1,1])*5e-3,'k--')
	plt.axis([0,50,1e-6,3e-2])
	plt.xlabel('$S1$',fontsize=18)
	plt.ylabel('$leakage~fraction$',fontsize=18)
	plt.xticks(np.arange(0,51,10),fontsize=18)
	plt.yticks([1e-6,1e-5,1e-4,1e-3,1e-2],fontsize=18)
	plt.grid(True)
	plt.show(0); plt.draw()
	plt.savefig("figs/disc7.pdf")


if 0:
	f5=plt.figure(5); plt.clf();
	detector = 'LUX'
	jj = 400
	Ed = 180	
	execfile("define-detector.py")
	eta_extraction = 0.95; dp[11] = eta_extraction
	g1 = 0.057; dp[0]=g1
	execfile("get-leakage.py")
	plt.semilogy(xe_s1,leakage_frac,'kx',color='orange')
	
	Ed = 180	
	execfile("define-detector.py")
	execfile("get-leakage.py")
	plt.semilogy(xe_s1,leakage_frac,'ko')

	Ed = 180	
	execfile("define-detector.py")
	eta_extraction = 0.95; dp[11] = eta_extraction
	execfile("get-leakage.py")
	plt.semilogy(xe_s1,leakage_frac,'k+')

	Ed = 730
	execfile("define-detector.py")
	eta_extraction = 0.95; dp[11] = eta_extraction
	execfile("get-leakage.py")
	plt.semilogy(xe_s1,leakage_frac,'k^',color='aqua')

	Ed = 180
	execfile("define-detector.py")
	eta_extraction = 0.95; dp[11] = eta_extraction
	g1 = 0.117*2; dp[0]=g1
	execfile("get-leakage.py")
	plt.semilogy(xe_s1,leakage_frac,'k*',color='indigo')


	plt.plot(np.array([0,50]),np.array([1,1])*5e-3,'k--')
	plt.axis([0,50,1e-6,3e-2])
	plt.xlabel('$S1$',fontsize=18)
	plt.ylabel('$leakage~fraction$',fontsize=18)
	plt.xticks(np.arange(0,51,10),fontsize=18)
	plt.yticks([1e-6,1e-5,1e-4,1e-3,1e-2],fontsize=18)
	plt.grid(True)
	plt.show(0); plt.draw()
	plt.savefig("figs/disc5.pdf")

