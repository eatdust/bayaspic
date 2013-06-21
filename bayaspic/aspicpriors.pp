!write down the prior associated with the aspic models. Syntax is of
!the form: (p1min,p1max,maptype1,p2min,p2max,maptype2...). If the end
!of inflation is an additional parameter, it is assumed to be the last
!one. Maptype denotes the relation between the aspic parameters and
!the one we are sampling on, "flat" means they are identical, "log"
!means we are sampling on the log10, "ln" means we are sampling on the
!log Neper, "mlog" means -log() as "mln" for -ln(), "inv" means
!inverse 1/(), "invsqrt" means the square root of the inverse 1/sqrt()
!etc...  The prior bounds are assumed to be on the sampling
!parameters.

ZEROPRIOR(hi)

ONEPRIOR(rchi,-48,-20,flat)
ONEPRIOR(rchiB,-65,100,flat)

ONEPRIOR(lfi,0.2,5,flat)
ONEPRIOR(lfi23,2._kp/3._kp,2._kp/3._kp,flat)
ONEPRIOR(lfi1,1.,1.,flat)
ONEPRIOR(lfi2,2.,2.,flat)
ONEPRIOR(lfi3,3.,3.,flat)
ONEPRIOR(lfi4,4.,4.,flat)

ONEPRIOR(rcmi,-7,-3,log)

ONEPRIOR(rcqi,-3,-0.1,log)

ONEPRIOR(ni,0,2.5,log)

ONEPRIOR(esi,0.1,6.,flat)
!ONEPRIOR(esi,0.5,2.,invsqrt)
!ONEPRIOR(esi,-3.,1.,log)
!ONEPRIOR(esi,sqrt(2._kp),sqrt(2._kp),flat)
!ONEPRIOR(esi,sqrt(2._kp/3._kp),sqrt(2._kp/3._kp),flat)


ONEPRIOR(pli,0.,1.1,flat)
!ONEPRIOR(pli,-4.,-0.1,log)

ONEPRIOR(kmii,0.382,4,log)
!ONEPRIOR(kmii,2.4095,2.7183,flat)

ONEPRIOR(hf1i,-3,3,log)

ONEPRIOR(li,-0.1_kp,0.1_kp,flat)
ONEPRIOR(li+,-3._kp,-1._kp,log)
ONEPRIOR(li-,-2.5_kp,-1._kp,mlog)


ONEPRIOR(rpi1,1.,1.5,flat)

ONEPRIOR(dwi,log(2._kp*sqrt(2._kp)),3,log)

ONEPRIOR(mhi,-2,2,log)
!ONEPRIOR(mhi,-2,0,log)
!ONEPRIOR(mhi,0,2,log)

ONEPRIOR(rgi,-4,4,log)
!ONEPRIOR(rgi,-4,0,log)
!ONEPRIOR(rgi,0,4,log)

ONEPRIOR(mssmi,-3,3,log)
!ONEPRIOR(mssmi,0.00002,0.0002,flat)

ONEPRIOR(ripi,-3,3,log)
!ONEPRIOR(ripi,0.00002,0.0002,flat)

ONEPRIOR(ai,-3,-0.29,log)

ONEPRIOR(cnai,-4,-0.18,log)

ONEPRIOR(cnbi,-5,-1.4,log)

ONEPRIOR(osti,1,4,log)

ONEPRIOR(wri,-3,3,log)
!ONEPRIOR(wri,1.,1.,flat)

!with xend

TWOPRIORS(rpi2,1.,1.5,flat,1,2,log)

TWOPRIORS(ii,0,10,flat,-1,4,log)
!TWOPRIORS(ii,0,1,iif,-1,4,log)
TWOPRIORS(iif,0,1,inv,-1,4,log)
!TWOPRIORS(ii,0.1,4,iiloglambda,-1,4,log)
TWOPRIORS(iilambda,0.1,4,log,-1,4,log)

! Warning: the second param here is phiend/phi0 (changed to phiend/Mp in map_aspic_param in wraspic.f90)
!TWOPRIORS(twi,-4,-1,twi,log(2._kp),log(20._kp),twi)
TWOPRIORS(twi,-4,-1,log,log(2._kp),log(20._kp),log)

TWOPRIORS(bsusybi,0.,0.3,flat,-200,-0,flat)
!TWOPRIORS(bsusybi,-3.,-1.,log,-200,-0,flat)

!TWOPRIORS(csi,-5,-1,log,0.,1.,csi)
TWOPRIORS(csi,-5,-1,log,0.,1.,flat)

!TWOPRIORS(cnci,-5,-1,log,1,100,cnci)
TWOPRIORS(cnci,-5,-1,log,1,100,flat)

!TWOPRIORS(imi,1,6,flat,1,100,imi)
TWOPRIORS(imi,1,6,flat,1,100,flat)

!witout xend

TWOPRIORS(cwi,4._kp*exp(1._kp),4._kp*exp(1._kp),flat,-5,-3,log)
!TWOPRIORS(cwi,4._kp*exp(1._kp),4._kp*exp(1._kp),flat,5*10.**(-5.),5*10.**(-4.),flat)

!TWOPRIORS(sfi,2,10,flat,-1.,2.,log)
!TWOPRIORS(sfi,2,10,flat,-1.,0.,log)
!TWOPRIORS(sfi,2,10,flat,0,2.,log)
!TWOPRIORS(sfi,1,1,flat,1.1,3.,log)
!TWOPRIORS(sfi,1,1,flat,-1.,0.,log)
!TWOPRIORS(sfi,1,1,flat,0.,2.,log)
!TWOPRIORS(sfi,2,2,flat,-1.,2.,log)
!TWOPRIORS(sfi,2,2,flat,-1.,0.,log)
!TWOPRIORS(sfi,2,2,flat,0.,2.,log)
!TWOPRIORS(sfi,3,3,flat,-1.,2.,log)
!TWOPRIORS(sfi,3,3,flat,-1.,0.,log)
!TWOPRIORS(sfi,3,3,flat,0.,2.,log)
!TWOPRIORS(sfi,4,4,flat,-1.,2.,log)
!TWOPRIORS(sfi,4,4,flat,-1.,0.,log)
TWOPRIORS(sfi,4,4,flat,0.,2.,log)


TWOPRIORS(kmiii,2,16,log,9,15,log)

TWOPRIORS(lmi1,0.01,0.99,flat,-2,2,log)
TWOPRIORS(lmio,1.1,6,flat,-2,2,flat)
!TWOPRIORS(lmi1,1.1,6,lmio,-2,2,lmio)


! Warning: the fine tuning condition on alpha must be activated
! when the first param lies in [-16,-1] (cf check_aspic_hardprior in aspicpriors.F08),
! and desactivated otherwise
!TWOPRIORS(gmssmi,-4.,-1.,logOfxMinusOne,-3.,3.,log)
!TWOPRIORS(gmssmi,-16.,-1.,logOfxMinusOne,-3.,3.,log)
!TWOPRIORS(gmssmi,-4.,-1.,logOfxMinusOne,0.00002,0.0002,flat)
!TWOPRIORS(gmssmi,-16.,-1.,logOfxMinusOne,0.00002,0.0002,flat)
!TWOPRIORS(gmssmi,-4.,-1.,logOfOneMinusX,-3.,3.,log)
!TWOPRIORS(gmssmi,-16.,-1.,logOfOneMinusX,-3.,3.,log)
!TWOPRIORS(gmssmi,-4.,-1.,logOfOneMinusX,0.00002,0.0002,flat)
!TWOPRIORS(gmssmi,-16.,-1.,logOfOneMinusX,0.00002,0.0002,flat)
!TWOPRIORS(gmssmi,0.9,1.1,flat,-3.,3.,log)
!TWOPRIORS(gmssmi,0.9,1.1,flat,0.00002,0.0002,flat)

TWOPRIORS(gmssmiA,-4.,-1.,log,-3.,3.,log)
TWOPRIORS(gmssmiB,-16.,-1.,log,-3.,3.,log)
TWOPRIORS(gmssmiC,-4.,-1.,log,0.00002,0.0002,flat)
TWOPRIORS(gmssmiD,-16.,-1.,log,0.00002,0.0002,flat)
TWOPRIORS(gmssmiE,-4.,-1.,log,-3.,3.,log)
TWOPRIORS(gmssmiF,-16.,-1.,log,-3.,3.,log)
TWOPRIORS(gmssmiG,-4.,-1.,log,0.00002,0.0002,flat)
TWOPRIORS(gmssmiH,-16.,-1.,log,0.00002,0.0002,flat)
TWOPRIORS(gmssmiI,0.9,1.1,flat,-3.,3.,log)
TWOPRIORS(gmssmiJ,0.9,1.1,flat,0.00002,0.0002,flat)


!!! Warning: the fine tuning condition on alpha must be activated
!!! when the first param lies in [-16,-1] (cf check_aspic_hardprior in aspicpriors.F08),
!!! and desactivated otherwise
!TWOPRIORS(gripi,-4.,-1.,logOfxMinusOne,-3.,3.,log)
!TWOPRIORS(gripi,-16.,-1.,logOfxMinusOne,-3.,3.,log)
!TWOPRIORS(gripi,-4.,-1.,logOfxMinusOne,0.00002,0.0002,flat)
!TWOPRIORS(gripi,-16.,-1.,logOfxMinusOne,0.00002,0.0002,flat)
!TWOPRIORS(gripi,-4.,-1.,logOfOneMinusX,-3.,3.,log)
!TWOPRIORS(gripi,-16.,-1.,logOfOneMinusX,-3.,3.,log)
!TWOPRIORS(gripi,-4.,-1.,logOfOneMinusX,0.00002,0.0002,flat)
!TWOPRIORS(gripi,-16.,-1.,logOfOneMinusX,0.00002,0.0002,flat)
!TWOPRIORS(gripi,0.9,1.1,flat,-3.,3.,log)
!TWOPRIORS(gripi,0.9,1.1,flat,0.00002,0.0002,flat)

TWOPRIORS(gripiA,-4.,-1.,log,-3.,3.,log)
TWOPRIORS(gripiB,-16.,-1.,log,-3.,3.,log)
TWOPRIORS(gripiC,-4.,-1.,log,0.00002,0.0002,flat)
TWOPRIORS(gripiD,-16.,-1.,log,0.00002,0.0002,flat)
TWOPRIORS(gripiE,-4.,-1.,log,-3.,3.,log)
TWOPRIORS(gripiF,-16.,-1.,log,-3.,3.,log)
TWOPRIORS(gripiG,-4.,-1.,log,0.00002,0.0002,flat)
TWOPRIORS(gripiH,-16.,-1.,log,0.00002,0.0002,flat)
TWOPRIORS(gripiI,0.9,1.1,flat,-3.,3.,log)
TWOPRIORS(gripiJ,0.9,1.1,flat,0.00002,0.0002,flat)

!!! Warning: the fine tuning condition on alpha must be activated
!!! when the first param lies in [0.5,0.5±10.**(-7.)] (cf check_aspic_hardprior in aspicpriors.F08),
!!! and desactivated otherwise
!TWOPRIORS(ti,0.5,0.6,flat,-5,-3,log)
!TWOPRIORS(ti,0.5,0.5+1d-7,flat,-5,-3,log)
!TWOPRIORS(ti,0.4,0.5,flat,-5,-3,log)
!TWOPRIORS(ti,0.5-1d-7,0.5,flat,-5,-3,log)
!TWOPRIORS(ti,0.4,0.6,flat,-5,-3,log)
!TWOPRIORS(ti,0.5-1d-7,0.5+1d-7,flat,-5,-3,log)

TWOPRIORS(tiA,0.5_kp,0.6_kp,flat,-5,-3,log)
TWOPRIORS(tiB,0.5_kp,0.5_kp+1d-7,flat,-5,-3,log)
TWOPRIORS(tiC,0.4_kp,0.5_kp,flat,-5,-3,log)
TWOPRIORS(tiD,0.5_kp-1d-7,0.5_kp,flat,-5,-3,log)
TWOPRIORS(tiE,0.4_kp,0.6_kp,flat,-5,-3,log)
TWOPRIORS(tiF,0.5_kp-1d-7,0.5_kp+1d-7,flat,-5,-3,log)


TWOPRIORS(bei,-3,3,log,-1.5,3.,log)

TWOPRIORS(psni,-5,-1,log,-2,1,log)
!!! Warning: activate fine tune by the reject condition,
!!! cf check_aspic_hardprior in aspicpriors.F08
TWOPRIORS(psniB,-9,1,log,-2,1,log)


TWOPRIORS(ncki,-5,-1,log,0.1,5.,flat)
!TWOPRIORS(ncki,-5,-1,log,-5.,-0.1,flat)

TWOPRIORS(oi,-3,-1,log,-3,-1,log)

TWOPRIORS(sbi,-5.,-1.,log,-5.,-1.,log)

TWOPRIORS(ssbi1,-5,1,log,-5,1,log)
TWOPRIORS(ssbi2,-5,-1,mlog,-5,-1,mlog)
TWOPRIORS(ssbi3,-5,1,log,-5,1,mlog)
TWOPRIORS(ssbi4,-5,1,log,-5,1,mlog)
TWOPRIORS(ssbi5,-5,-1,mlog,-5,-1,log)
TWOPRIORS(ssbi6,-5,1,mlog,-5,1,log)

!with xend

THREEPRIORS(lmi2,0.01,0.99,flat,-2,2,log,-5,5,log)
!THREEPRIORS(lmi2,1.1,6,lmio,-2,2,lmio,-5,5,log)

THREEPRIORS(rmi1,0.05,0.2,flat,-2,0,log,1./exp(1.),1.,flat)
THREEPRIORS(rmi2,0.05,0.2,flat,-2,0,log,1,exp(1.),flat)
THREEPRIORS(rmi3,-0.2,-0.05,flat,-2,0,log,1./exp(1.),1,flat)
THREEPRIORS(rmi4,-0.2,-0.05,flat,-2,0,log,1,exp(1.),flat)

THREEPRIORS(vhi,1.01,6,flat,0.,3,log,0,1,flat)

THREEPRIORS(dsi,1.01,6,flat,-9,-7,log,-1,4,log)

!last params is xend/xendmax
THREEPRIORS(cndi,-2,-1,log,1,6,flat,0,1,flat)

!witout xend

THREEPRIORS(gmlfi,2,2,flat,2,2,flat,-5,1,log)
!THREEPRIORS(gmlfi,1,6,flat,1,6,flat,-5,1,log)
!THREEPRIORS(gmlfi,2,2,flat,1,1,flat,-5,1,log)
!THREEPRIORS(gmlfi,2,2,flat,3,3,flat,-5,1,log)
!THREEPRIORS(gmlfi,3,3,flat,1,1,flat,-5,1,log)
!THREEPRIORS(gmlfi,3,3,flat,2,2,flat,-5,1,log)
!THREEPRIORS(gmlfi,3,3,flat,3,3,flat,-5,1,log)
!THREEPRIORS(gmlfi,4,4,flat,1,1,flat,-5,1,log)
!THREEPRIORS(gmlfi,4,4,flat,2,2,flat,-5,1,log)
!THREEPRIORS(gmlfi,4,4,flat,3,3,flat,-5,1,log)

THREEPRIORS(lpi1,1,6,flat,1,6,flat,-3,3,log)
!THREEPRIORS(lpi1,4,4,flat,1,1,flat,-3,3,log)
!THREEPRIORS(lpi1,4,4,flat,3,3,flat,-3,3,log)

THREEPRIORS(lpi2,1,6,flat,1,6,flat,-3,3,log)
!THREEPRIORS(lpi2,4,4,flat,1,1,flat,-3,3,log)
!THREEPRIORS(lpi2,4,4,flat,3,3,flat,-3,3,log)

THREEPRIORS(lpi3,1,6,flat,1,6,flat,-3,3,log)
!THREEPRIORS(lpi3,4,4,flat,1,1,flat,-3,3,log)
!THREEPRIORS(lpi3,4,4,flat,3,3,flat,-3,3,log)
