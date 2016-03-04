#
if detector == 'XENON10':
	xe_be = np.append([2],[np.arange(5,55,5)]) # 80
	xe_s1 = np.append([3.5],[np.arange(7.5,50,5)]) # 75
	xermu = np.array([2.77,2.644,2.535,2.464,2.416,2.385,2.367,2.356,2.349,2.347])#,2.346,2.347,2.351,2.352,2.357])
	xer3s = np.array([2.383,2.258,2.143,2.068,2.020,1.986,1.963,1.949,1.937,1.933])#,1.929,1.925,1.924,1.924,1.924])
	xer1s = (xermu-xer3s)/3
	xnrmu = np.array([2.398,2.241,2.160,2.108,2.069,2.045,2.010,1.969,1.939,1.935])#,1.904,1.892,1.877,1.846,1.844])
	xnr3s = np.array([1.866,1.805,1.846,1.833,1.788,1.819,1.788,1.760,1.739,1.741])
	xnr1s = (xnrmu-xnr3s)/3
	
	g1 = 0.08
	double_phe_fraction = 0.2
	npmttop = 48
	npmtbot = 41
	resolution_pmt = 0.59#*np.ones(npmttop+npmtbot)
	g2prime = 0.15 # actual geometric * QE light collection
	gasgain = 160
	g2 = gasgain*g2prime # = 24
	resolution_e = np.sqrt(g2)/g2
		
	# simple 50:50 assumption for S2
	pmtps2 = np.ones(npmttop+npmtbot); pmtps2.fill(1.0/(npmttop+npmtbot))
		
	nco = 2
	softwareThreshold_phe = 1/3.
	S2Threshold_phe = 300
	S1window_ns = 300 # ns
	eta_extraction = 0.95

elif detector == 'LUX':
	xe_be = np.append([np.arange(2,4.1,1)],[np.arange(5,55,5)]) # 80
	xe_s1 = np.append([np.arange(2.5,5,1)],[np.arange(7.5,50,5)]) # 75
	xermu = np.array([2.284,2.175,2.097,1.951,1.814,1.731,1.670,1.621,1.585,1.552,1.524,1.500])
	xer128s = np.array([2.109,2.009,1.936,1.796,1.667,1.586,1.529,1.482,1.448,1.417,1.393,1.371])
	xer1s = (xermu-xer128s)/1.28
	xer3s = xermu-xer1s*3
	xnrmu = np.array([1.830,1.720,1.650,1.549,1.481,1.429,1.383,1.349,1.311,1.287,1.262,1.236])
	xnr128s = np.array([1.637,1.525,1.442,1.357,1.330,1.292,1.266,1.237,1.206,1.188,1.163,1.138])
	xnr1s = (xnrmu-xnr128s)/1.28
	
	if 1:
		g1 = 0.115
		eta_extraction = 0.50
	elif 1:
		g1 = 0.23
		eta_extraction = 0.950
	elif 0:
		g1 = 0.5
		eta_extraction = 0.950

	double_phe_fraction = 0.2

	npmttop = 61
	npmtbot = 61
	resolution_pmt = 0.37
	g2prime = 0.15 # actual geometric * QE light collection
	gasgain = 160
	g2 = gasgain*g2prime/2 # = 12
	resolution_e = np.sqrt(g2)/g2

	# simple 50:50 assumption for S2, but using only bottom array for S2
	pmtps2 = np.ones(npmtbot); pmtps2.fill(0.5/(npmtbot))

	nco = 3
	softwareThreshold_phe = 1/3.
	S2Threshold_phe = g2*8.1 # based on 200/24.6
	S1window_ns = 100 # ns
	
elif detector == 'XENON100':
	xe_be = np.append([2],[np.arange(5,55,5)]) # bands from Aprile talk at WONDER 2010, LNGS
	xe_s1 = np.append([3.5],[np.arange(7.5,50,5)]) # 75
	xermu = np.array([2.714,2.550,2.438,2.378,2.349,2.329,2.318,2.316,2.312,2.310])
	xer3s = np.array([2.274,2.115,2.012,1.959,1.929,1.913,1.912,1.905,1.899,1.899])
	xer1s = (xermu-xer3s)/3
	xnrmu = np.array([2.357,2.159,2.086,2.045,2.009,1.975,1.941,1.908,1.878,1.854])
	xnr1s = np.zeros(len(xnrmu))
	
	g1 = 0.057 # 2.2/0.58*0.015
	double_phe_fraction = 0.2
	npmttop = 98
	npmtbot = 80
	resolution_pmt = 0.59
	g2prime = 0.15 # actual geometric * QE light collection
	gasgain = 119
	g2 = gasgain*g2prime # = 17.8 in early operation 2010, see Lang talk at UC Davis
	resolution_e = np.sqrt(g2)/g2
		
	# simple 50:50 assumption for S2
	pmtps2 = np.ones(npmttop+npmtbot); pmtps2.fill(1.0/(npmttop+npmtbot))
		
	nco = 2
	softwareThreshold_phe = 1/3.
	S2Threshold_phe = 300
	S1window_ns = 20 # ns
	eta_extraction = 0.95

elif detector == 'ZEPLIN3':
	secretFactor = 1.8/738 # this factor obtained via email from T Sumner, it is phe/keVee at 122 keVee for S1 and S2
 	#value on ZEPLIN3 plot translates to usual coords as np.log10(10**(value)/secretFactor)
 	xe_be = 1.8*np.arange(2,16.5,1)
 	xe_s1 = 1.8*np.arange(2.5,16,1) 
 	xermu = np.log10(10**np.array([0.224,0.121,0.058,0.001,-0.026,-0.053,-0.077,-0.100,-0.111,-0.125,-0.131,-0.131,-0.146,-0.139]) / secretFactor)
 	xer1s = np.array([0.154,0.148,0.155,0.160,0.159,0.167,0.167,0.170,0.171,0.176,0.176,0.181,0.182,0.181])
 	xnrmu = np.log10(10**np.array([-0.249,-0.330,-0.384,-0.431,-0.469,-0.507,-0.534,-0.561,-0.583,-0.607,-0.626,-0.642,-0.661,-0.676]) / secretFactor)
 	xnr1s = (xnrmu - np.log10(10**np.array([-0.384,-0.451,-0.499,-0.542,-0.575,-0.606,-0.630,-0.654,-0.674,-0.692,-0.710,-0.725,-0.740,-0.754])/ secretFactor) )
	
	g1 = 0.075 # 5*0.015
	double_phe_fraction = 0.2
	npmttop = 0
	npmtbot = 31
	resolution_pmt = 0.63 # arxiv:0905.2523
	g2prime = 0.15 # actual geometric * QE light collection
	gasgain = 204
	g2 = gasgain*g2prime # = 30.6
	resolution_e = np.sqrt(g2)/g2
		
	# simple 50:50 assumption for S2
	pmtps2 = np.ones(npmttop+npmtbot); pmtps2.fill(1.0/(npmttop+npmtbot))
		
	nco = 3
	softwareThreshold_phe = 1/3.
	# smallest observable event has ~15 e-
	S2Threshold_phe = g2*11 # stated 11 e- trigger, so assume higher software threshold (a few electrons to turn on)
	S1window_ns = 300 # ns
	#eta_extraction = 0.83 # stated value, but does not jive?
	#eta_extraction = 0.75 # ER agreement, but at 730 V/cm !?
	eta_extraction = 0.52# NR ~agreement, assuming stated single e- size
	
elif detector == 'LZ':
	g1 = 0.075
	double_phe_fraction = 0.2
	npmttop = 247
	npmtbot = 241
	resolution_pmt = 0.37#*np.ones(npmttop+npmtbot)
	g2 = 25
	resolution_e = np.sqrt(g2)/g2
	nco = 3
	softwareThreshold_phe = 1/3.
	S1window_ns = 100 # ns
	eta_extraction = 0.95


# probability for photon to impinge on each PMT (average assumption)
try:
	pmtptop = np.ones(npmttop); pmtptop.fill(0.20/npmttop)
except: #ZeroDivisionError:
	erro = sys.exc_info()[0]
	print erro
	pmtptop = []
pmtpbot = np.ones(npmtbot); pmtpbot.fill(0.80/npmtbot)
pmtps1 = np.concatenate((pmtptop,pmtpbot),axis=0)


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