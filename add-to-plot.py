#
if detector=='XENON10': # xenon10 from my thesis
	x10S1 = np.array([3,5,7,9,11.5,14.5,17.5,20.5,23.5,27.5,32.5,37.5])
	x10yy = np.array([2.89,2.74,2.66,2.60,2.55,2.50,2.46,2.44,2.43,2.40,2.39,2.37])
	plt.plot(x10S1,x10yy,'b.-',ms=7)#,markeredgecolor="None")
if 1: # band mu from PR article
	plt.plot(xe_s1,xermu,'b.-',ms=7,markerfacecolor="None",alpha=1.,markeredgewidth=1)
	plt.plot(xe_s1,xermu-xer1s,'b.--',ms=7,markerfacecolor="None",alpha=1.,markeredgewidth=1)
	plt.plot(xe_s1,xermu+xer1s,'b.--',ms=7,markerfacecolor="None",alpha=1.,markeredgewidth=1)
	plt.plot(xe_s1,xnrmu,'r.-',ms=7,markerfacecolor="None",alpha=1.,markeredgewidth=1)
	plt.plot(xe_s1,xnrmu-xnr1s,'r.--',ms=7,markerfacecolor="None",alpha=1.,markeredgewidth=1)
	plt.plot(xe_s1,xnrmu+xnr1s,'r.--',ms=7,markerfacecolor="None",alpha=1.,markeredgewidth=1)


energies = [1,2,5,10,15,20]
s1s = np.arange(0.01,50,0.01)
for i in range(len(energies)):
	s2s = (energies[i]/0.0138 - s1s/g1)*g2
	yy = np.log10((s2s)/s1s)
	plt.plot(s1s,yy,'k:')
	
plt.xlabel('$S1$',fontsize=18)
plt.ylabel('$log_{10}(S2/S1)$',fontsize=18)
#plt.grid(True)
plt.xticks(np.arange(0,51,10),fontsize=18)
plt.yticks(np.arange(0.5,3.5,0.5),fontsize=18)
plt.axis([0,50,0.4,3.3])
plt.minorticks_on()
plt.show(0); plt.draw()
