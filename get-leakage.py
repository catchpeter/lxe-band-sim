
particle='electron'
er = np.arange(0.75,20,0.25) # keV
#x0 = define_x0(particle,Ed)
(Ne,Ng,S1,S2,yy) = bandsim(particle,Ed,er,dp,pmtps1,pmtps2,jj)
S2[S2<S2Threshold_phe] = 0
(sim_er_mu,sim_er_1s,sim_bc) = getBandMuSigma(xe_s1,xe_be,S1,yy,detector)

particle='neutron'
er = np.arange(1,70,0.25) # keV
#x0 = define_x0(particle,Ed)
(Ne,Ng,S1,S2,yy) = bandsim(particle,Ed,er,dp,pmtps1,pmtps2,jj)
S2[S2<S2Threshold_phe] = 0
(sim_nr_mu,sim_nr_1s,sim_bc) = getBandMuSigma(xe_s1,xe_be,S1,yy,detector)

leakage_frac = getLeakage(xe_s1,sim_er_mu,sim_er_1s,sim_nr_mu)
