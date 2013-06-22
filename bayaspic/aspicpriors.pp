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

ONEPRIOR(rchio,-48._kp,-20._kp,flat)
ONEPRIOR(rchi,-65._kp,100._kp,flat)

ONEPRIOR(lfi,0.2_kp,5._kp,flat)
ONEPRIOR(lfi23,2._kp/3._kp,2._kp/3._kp,flat)
ONEPRIOR(lfi1,1._kp,1._kp,flat)
ONEPRIOR(lfi2,2._kp,2._kp,flat)
ONEPRIOR(lfi3,3._kp,3._kp,flat)
ONEPRIOR(lfi4,4._kp,4._kp,flat)

ONEPRIOR(rcmi,-7._kp,-3._kp,log)

ONEPRIOR(rcqi,-3._kp,-0.1_kp,log)

ONEPRIOR(ni,0._kp,2.5_kp,log)

ONEPRIOR(esi,0.1_kp,6._kp,flat)
ONEPRIOR(esio,0.5_kp,2._kp,invsqrt)
ONEPRIOR(esil,-3._kp,1._kp,log)
ONEPRIOR(esisqrt2,sqrt(2._kp),sqrt(2._kp),flat)
ONEPRIOR(esisqrt23,sqrt(2._kp/3._kp),sqrt(2._kp/3._kp),flat)


ONEPRIOR(plip,0._kp,1.1_kp,flat)
ONEPRIOR(pli,-4._kp,0._kp,log)

ONEPRIOR(kmii,0.382_kp,4._kp,log)
ONEPRIOR(kmiiV>0,2.4095_kp,2.7183_kp,flat)

ONEPRIOR(hf1i,-3._kp,3._kp,log)

ONEPRIOR(li,-0.1_kp,0.1_kp,flat)
ONEPRIOR(li+,-3._kp,-1._kp,log)
ONEPRIOR(li-,-2.8_kp,-1._kp,mlog)

ONEPRIOR(rpi1,1._kp,1.5_kp,flat)

ONEPRIOR(dwi,log(2._kp*sqrt(2._kp)),5._kp,log)

ONEPRIOR(mhi,-2._kp,2._kp,log)
ONEPRIOR(mhis,-2._kp,0._kp,log)
ONEPRIOR(mhil,0._kp,2._kp,log)

ONEPRIOR(rgi,-4._kp,4._kp,log)
ONEPRIOR(rgi1Over16,0.0625_kp,0.0625_kp,log)
ONEPRIOR(rgis,-4._kp,0._kp,log)
ONEPRIOR(rgil,0._kp,4._kp,log)

ONEPRIOR(mssmip,-3._kp,3._kp,log)
ONEPRIOR(mssmio,0.00002_kp,0.0002_kp,flat)

ONEPRIOR(ripip,-3._kp,3._kp,log)
ONEPRIOR(ripio,0.00002_kp,0.0002_kp,flat)

ONEPRIOR(ai,-3._kp,log(0.512378_kp)/log(10._kp),log)

ONEPRIOR(cnai,-4._kp,-0.27_kp,log)

ONEPRIOR(cnbi,-5._kp,-1.4_kp,log)

ONEPRIOR(osti,1._kp,4._kp,log)

ONEPRIOR(wrig,-3._kp,3._kp,log)
ONEPRIOR(wrio,1._kp,1._kp,flat)

!with xend

TWOPRIORS(rpi2,1.,1.5,flat,1,2,log)

TWOPRIORS(ii,0._kp,10._kp,flat,-1._kp,4._kp,log)
TWOPRIORS(iif,0.1_kp,1._kp,inv,-1._kp,4._kp,log)
TWOPRIORS(iilambda,0.1_kp,4._kp,log,-1._kp,4._kp,log)

! Warning: the second param here is phiend/phi0 (changed to phiend/Mp in redefine_aspic_params located in aspicpriors.F8)
TWOPRIORS(twi,-4._kp,-1._kp,log,log(2._kp)/log(10._kp),log(20._kp)/log(10._kp),log)

TWOPRIORS(bsusybif,0._kp,0.3_kp,flat,-200._kp,-0._kp,flat)
TWOPRIORS(bsusybil,-3._kp,-1._kp,log,-200._kp,-0._kp,flat)

TWOPRIORS(csi,-5,-1,log,0.,1.,flat)

TWOPRIORS(cnci,-5,-1,log,1,100,flat)

!TWOPRIORS(imi,1,6,flat,1,100,imi)
TWOPRIORS(imi,1,6,flat,1,100,flat)

!witout xend

TWOPRIORS(cwil,4._kp*exp(1._kp),4._kp*exp(1._kp),flat,-5._kp,-3._kp,log)
TWOPRIORS(cwif,4._kp*exp(1._kp),4._kp*exp(1._kp),flat,5._kp*10.**(-5._kp),5._kp*10.**(-4._kp),flat)

TWOPRIORS(sfi,2._kp,10._kp,flat,-1._kp,2._kp,log)
TWOPRIORS(sfis,2._kp,10._kp,flat,-1._kp,0._kp,log)
TWOPRIORS(sfil,2._kp,10._kp,flat,0._kp,2._kp,log)
TWOPRIORS(sfi1,1._kp,1._kp,flat,1.1_kp,3._kp,log)
TWOPRIORS(sfi2,2._kp,2._kp,flat,-1._kp,2._kp,log)
TWOPRIORS(sfi2s,2._kp,2._kp,flat,-1._kp,0._kp,log)
TWOPRIORS(sfi2l,2._kp,2._kp,flat,0._kp,2._kp,log)
TWOPRIORS(sfi3,3._kp,3._kp,flat,-1._kp,2._kp,log)
TWOPRIORS(sfi3s,3._kp,3._kp,flat,-1._kp,0._kp,log)
TWOPRIORS(sfi3l,3._kp,3._kp,flat,0._kp,2._kp,log)
TWOPRIORS(sfi4,4._kp,4._kp,flat,-1._kp,2._kp,log)
TWOPRIORS(sfi4s,4._kp,4._kp,flat,-1._kp,0._kp,log)
TWOPRIORS(sfi4l,4._kp,4._kp,flat,0._kp,2._kp,log)


TWOPRIORS(kmiii,2._kp,16._kp,log,9._kp,15._kp,log)

TWOPRIORS(lmi1p,0.1_kp,1._kp,flat,-2._kp,2._kp,log)
TWOPRIORS(lmi1o,1._kp,6._kp,flat,-2._kp,2._kp,flat)


TWOPRIORS(gmssmiA,-4._kp,-1._kp,log,-3._kp,3._kp,log)
TWOPRIORS(gmssmiB,-20._kp,-1._kp,log,-3._kp,3._kp,log)
TWOPRIORS(gmssmiC,-4._kp,-1._kp,log,0.00002_kp,0.0002_kp,flat)
TWOPRIORS(gmssmiD,-28._kp,-16._kp,log,0.00002_kp,0.0002_kp,flat)
TWOPRIORS(gmssmiE,-4._kp,-1._kp,log,-3._kp,3._kp,log)
TWOPRIORS(gmssmiF,-20._kp,-1._kp,log,-3._kp,3._kp,log)
TWOPRIORS(gmssmiG,-4._kp,-1._kp,log,0.00002_kp,0.0002_kp,flat)
TWOPRIORS(gmssmiH,-28._kp,-1._kp,log,0.00002_kp,0.0002_kp,flat)
TWOPRIORS(gmssmiI,0.9_kp,1.1_kp,flat,-3._kp,3._kp,log)
TWOPRIORS(gmssmiJ,0.9_kp,1.1_kp,flat,0.00002_kp,0.0002_kp,flat)


TWOPRIORS(gripiA,-4._kp,-1._kp,log,-3._kp,3._kp,log)
TWOPRIORS(gripiB,-20._kp,-1._kp,log,-3._kp,3._kp,log)
TWOPRIORS(gripiC,-4._kp,-1._kp,log,0.00002_kp,0.0002_kp,flat)
TWOPRIORS(gripiD,-28._kp,-16._kp,log,0.00002_kp,0.0002_kp,flat)
TWOPRIORS(gripiE,-4._kp,-1._kp,log,-3._kp,3._kp,log)
TWOPRIORS(gripiF,-20._kp,-1._kp,log,-3._kp,3._kp,log)
TWOPRIORS(gripiG,-4._kp,-1._kp,log,0.00002_kp,0.0002_kp,flat)
TWOPRIORS(gripiH,-28._kp,-1._kp,log,0.00002_kp,0.0002_kp,flat)
TWOPRIORS(gripiI,0.9_kp,1.1_kp,flat,-3._kp,3._kp,log)
TWOPRIORS(gripiJ,0.9_kp,1.1_kp,flat,0.00002_kp,0.0002_kp,flat)


TWOPRIORS(tiA,0.5_kp,0.6_kp,flat,-5._kp,-3._kp,log)
TWOPRIORS(tiB,0.5_kp,0.5_kp+10._kp**(-7._kp),flat,-5._kp,-3._kp,log)
TWOPRIORS(tiC,0.4_kp,0.5_kp,flat,-5._kp,-3._kp,log)
TWOPRIORS(tiD,0.5_kp-10._kp**(-7._kp),0.5_kp,flat,-5._kp,-3._kp,log)
TWOPRIORS(tiE,0.4_kp,0.6_kp,flat,-5._kp,-3._kp,log)
TWOPRIORS(tiF,0.5_kp-10._kp**(-7._kp),0.5_kp+10._kp**(-7._kp),flat,-5._kp,-3._kp,log)
TWOPRIORS(tiOneHalf,0.5_kp,0.5_kp,flat,-5._kp,-3._kp,log)


TWOPRIORS(bei,-3._kp,3._kp,log,-1.5_kp,3._kp,log)

TWOPRIORS(psni,-5._kp,-1._kp,log,-2._kp,1._kp,log)
!! Warning: the first parameter here is alpha/f**2 (changed to alpha in redefine_aspic_params located in aspicpriors.F8)
TWOPRIORS(psnift,-5._kp,-1._kp,log,-2._kp,1._kp,log)


TWOPRIORS(ncki,-5,-1,log,0.1,5.,flat)
!TWOPRIORS(ncki,-5,-1,log,-5.,-0.1,flat)

TWOPRIORS(oi,-3._kp,-1._kp,log,-3._kp,-1._kp,log)

TWOPRIORS(sbi,-5._kp,-1._kp,log,-5._kp,-1._kp,log)

TWOPRIORS(ssbi1,-5,1,log,-5,1,log)
TWOPRIORS(ssbi2,-5,-1,mlog,-5,-1,mlog)
TWOPRIORS(ssbi3,-5,1,log,-5,1,mlog)
TWOPRIORS(ssbi4,-5,1,log,-5,1,mlog)
TWOPRIORS(ssbi5,-5,-1,mlog,-5,-1,log)
TWOPRIORS(ssbi6,-5,1,mlog,-5,1,log)

!with xend

THREEPRIORS(lmi2p,0.1_kp,0.99_kp,flat,-2._kp,2._kp,log,-5._kp,10._kp,log)
THREEPRIORS(lmi2o,1.1_kp,6._kp,flat,-2._kp,2._kp,flat,-5._kp,10._kp,log)


THREEPRIORS(rmi1,0.05,0.2,flat,-2,0,log,1./exp(1.),1.,flat)
THREEPRIORS(rmi2,0.05,0.2,flat,-2,0,log,1,exp(1.),flat)
THREEPRIORS(rmi3,-0.2,-0.05,flat,-2,0,log,1./exp(1.),1,flat)
THREEPRIORS(rmi4,-0.2,-0.05,flat,-2,0,log,1,exp(1.),flat)

THREEPRIORS(vhi,1.01,6,flat,0.,3,log,0,1,flat)

THREEPRIORS(dsi,1.01,6,flat,-9,-7,log,-1,4,log)

THREEPRIORS(cndi,-2,-1,log,1.1,6,flat,0,100.,flat)

!witout xend

THREEPRIORS(gmlfi22,2,2,flat,2,2,flat,-5,1,log)
!THREEPRIORS(gmlfi,1,6,flat,1,6,flat,-5,1,log)
!THREEPRIORS(gmlfi21,2,2,flat,1,1,flat,-5,1,log)
!THREEPRIORS(gmlfi23,2,2,flat,3,3,flat,-5,1,log)
!THREEPRIORS(gmlfi31,3,3,flat,1,1,flat,-5,1,log)
!THREEPRIORS(gmlfi32,3,3,flat,2,2,flat,-5,1,log)
!THREEPRIORS(gmlfi33,3,3,flat,3,3,flat,-5,1,log)
!THREEPRIORS(gmlfi41,4,4,flat,1,1,flat,-5,1,log)
!THREEPRIORS(gmlfi42,4,4,flat,2,2,flat,-5,1,log)
!THREEPRIORS(gmlfi43,4,4,flat,3,3,flat,-5,1,log)

THREEPRIORS(lpi1,1,6,flat,1,6,flat,-3,3,log)
!THREEPRIORS(lpi1,4,4,flat,1,1,flat,-3,3,log)
!THREEPRIORS(lpi1,4,4,flat,3,3,flat,-3,3,log)

THREEPRIORS(lpi2,1,6,flat,1,6,flat,-3,3,log)
!THREEPRIORS(lpi2,4,4,flat,1,1,flat,-3,3,log)
!THREEPRIORS(lpi2,4,4,flat,3,3,flat,-3,3,log)

THREEPRIORS(lpi3,1,6,flat,1,6,flat,-3,3,log)
!THREEPRIORS(lpi3,4,4,flat,1,1,flat,-3,3,log)
!THREEPRIORS(lpi3,4,4,flat,3,3,flat,-3,3,log)
