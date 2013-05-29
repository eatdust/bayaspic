!write down the prior associated with the aspic models. Syntax is of
!the form: (p1min,p1max,maptype1,p2min,p2max,maptype2...). If the end
!of inflation is an additional parameter, it is assumed to be the last
!one. Maptype denotes the relation between the aspic parameters and
!the one we are sampling on, "flat" means they are identical, "log"
!means we are sampling on the log10, "ln" means we are sampling on the
!log Neper, "mlog" means -log() as "mln" for -ln(). The prior bounds are
!assumed to be on the sampling parameters.

ZEROPRIOR(hi)

ONEPRIOR(rchi,5,20,flat)
ONEPRIOR(lfi,4,4,flat)
ONEPRIOR(rcmi,-6,-3,log)
ONEPRIOR(rcqi,-3,-0.1,log)
ONEPRIOR(ni,0,3,log)
ONEPRIOR(esi,-3,1,log)
ONEPRIOR(pli,0.,sqrt(2./3.),flat)
ONEPRIOR(kmii,0.382,4,log)
ONEPRIOR(hf1i,-3,3,log)
ONEPRIOR(li,-0.1,0.1,flat)
ONEPRIOR(rpi1,1,1.5,flat)
ONEPRIOR(dwi,log(2.*sqrt(2.)),3,log)
ONEPRIOR(mhi,-2,2,log)
ONEPRIOR(rgi,-4,4,log)
ONEPRIOR(mssmi,-3,-2,log)
ONEPRIOR(ripi,-5,-4,log)
ONEPRIOR(ai,-3,-0.3,log)
ONEPRIOR(cnai,-4,-0.2,log)
ONEPRIOR(cnbi,-5,-1.4,log)
ONEPRIOR(osti,1,3,log)
ONEPRIOR(wrhi,-3,3,log)

!with xend

TWOPRIORS(rpi2,1,1.5,flat,1,2,log)
TWOPRIORS(ii,-1,1,log,0,2,log)
TWOPRIORS(twi,-3,1,log,2,20,flat)
TWOPRIORS(bsusybi,-1.,-1.,log,-200,-0,flat)
TWOPRIORS(csi,-5,0,log,-10,10,flat)
TWOPRIORS(cnci,-5,0,log,0,2,log)
TWOPRIORS(imi,1,6,flat,0,2,log)

!witout xend

TWOPRIORS(cwi,4._kp*exp(1._kp),4._kp*exp(1._kp),flat,-5,-3,log)
TWOPRIORS(sfi,2,10,flat,-0.2,3,log)
TWOPRIORS(kmiii,4,6,log,10,14,log)
TWOPRIORS(lmi1,0.01,1,flat,-2,2,log)
TWOPRIORS(gmssmi,0.95,1.05,flat,5e-5,5e-4,flat)
TWOPRIORS(gripi,0.95,1.05,flat,-2,0,log)
TWOPRIORS(ti,0.4,0.6,flat,-6,-2,log)
TWOPRIORS(bei,-3,3,log,-1.5,3.,log)
TWOPRIORS(psni,-6,0,log,-1,2,log)
TWOPRIORS(ncki,-3,3,log,-4,0,log)
TWOPRIORS(oi,-3,-1,log,-4,0,log)
TWOPRIORS(sbi,0,1,flat,1e-5,0.1,flat)
TWOPRIORS(ssbi1,0,10,flat,-3,1,log)
TWOPRIORS(ssbi2,-5,-1,log,-3,3,log)
TWOPRIORS(ssbi3,0,10,flat,-3,3,mlog)
TWOPRIORS(ssbi4,-4,1,log,-5,-3,mlog)
TWOPRIORS(ssbi5,-10,0,flat,-6,-4,log)
TWOPRIORS(ssbi6,-10,0,flat,-5,0,log)

!with xend

THREEPRIORS(lmi2,0.5,0.99,flat,0.1,10,flat,10,100,flat)
THREEPRIORS(rmi1,1e-3,0.2,flat,-3,0,log,0.36788,1,flat)
THREEPRIORS(rmi2,1e-3,0.2,flat,-3,0,log,1,2.7183,flat)
THREEPRIORS(rmi3,-1e3,-0.2,flat,-3,0,log,0.36788,1,flat)
THREEPRIORS(rmi4,-1e3,-0.2,flat,-3,0,log,1,2.7183,flat)
THREEPRIORS(vhi,1.01,6,flat,-0.2,3,log,0,1,flat)
THREEPRIORS(dsi,1.01,6,flat,-10,0,log,0,3,log)
THREEPRIORS(cndi,-2,-1,log,-1,1,log,0,3,log)

!witout xend

THREEPRIORS(gmlfi,1,5,flat,1,5,flat,-3,3,log)
THREEPRIORS(lpi1,4,4,flat,3,3,flat,-3,3,log)
THREEPRIORS(lpi2,4,4,flat,3,3,flat,-3,-3,log)
THREEPRIORS(lpi3,4,4,flat,3,3,flat,-3,-3,log)
