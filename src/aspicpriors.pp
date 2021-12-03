!write down the prior associated with the aspic models. Syntax is of
!the form:

!   NPRIOR(modelname,p1min,p1max,maptype1,p2min,p2max,maptype2...)

!If the end of inflation is an additional parameter, it is assumed to
!come *just after* the potential parameters. Any other additional and
!non-potential parameters comes next.

!The "maptype" entry denotes the relation between the aspic parameters
!and the one we are sampling on, "flat" means they are identical,
!"log" means we are sampling on the log10(aspicparam), "ln" means we
!are sampling on the log Neper, "mlog" means -log() as "mln" for
!-ln(), "inv" means inverse 1/(), "invsqrt" means the square root of
!the inverse 1/sqrt() etc...

!The prior bounds are assumed to be on the actual *sampled
!parameters*.  Therefore, if maptype is "ln", and the prior bounds are
![0,1], it means that the aspic parameters will vary in the range
![1,e].

!Redefinition of aspic parameters can made in the function
!redefine_aspic_params() in aspicpriors.F08. This can be used if the
!observable parameter is a combination of the fundamental aspic
!parameters. This function is called in map_aspic_params() (see
!wraspic.f90) and assumes that the parameters are the aspic ones, not
!the ones on which we sample. Therefore, it is called after the
!sampled parameter inverse rescaling.

ZEROPRIOR(si)

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
ONEPRIOR(esi8,8._kp,8._kp,flat)

ONEPRIOR(plip,0._kp,1.1_kp,flat)
ONEPRIOR(pli,-4._kp,0._kp,log)

ONEPRIOR(kmii,0.382_kp,4._kp,log)
ONEPRIOR(kmiivp,2.4095_kp,2.7183_kp,flat)

ONEPRIOR(hf1i,-3._kp,3._kp,log)

ONEPRIOR(li,-0.3_kp,0.3_kp,flat)
ONEPRIOR(lip,log10(0.003_kp),log10(0.3_kp),log)
ONEPRIOR(lin,-0.4_kp,-0.1_kp,flat)

ONEPRIOR(rpi1,1._kp+epsilon(1._kp),1.5_kp,flat)
ONEPRIOR(rpi3,0.8_kp,1._kp-epsilon(1._kp),flat)

ONEPRIOR(dwi,log(2._kp*sqrt(2._kp)),5._kp,log)
ONEPRIOR(dwi25,25._kp,25._kp,flat)

ONEPRIOR(mhi,-2._kp,2._kp,log)
ONEPRIOR(mhis,-2._kp,0._kp,log)
ONEPRIOR(mhil,0._kp,2._kp,log)
ONEPRIOR(mhixs,0.01_kp,0.01_kp,flat)

ONEPRIOR(rgi,-4._kp,4._kp,log)
ONEPRIOR(rgi116,0.0625_kp,0.0625_kp,log)
ONEPRIOR(rgis,-4._kp,0._kp,log)
ONEPRIOR(rgil,0._kp,4._kp,log)

!ripi requires quad-precision
ONEPRIOR(mssmip,-3._kp,3._kp,log)
ONEPRIOR(mssmio,0.00002_kp,0.0002_kp,flat)

!ripi requires quad-precision
ONEPRIOR(ripip,-3._kp,3._kp,log)
ONEPRIOR(ripio,0.00002_kp,0.0002_kp,flat)
ONEPRIOR(ripiS,10._kp,50._kp,flat)

!ONEPRIOR(ai,-3._kp,log10(0.512378_kp),log)
ONEPRIOR(ai,-3._kp,log10(0.51_kp),log)

ONEPRIOR(cnai,-4._kp,-0.27_kp,log)

ONEPRIOR(cnbi,-5._kp,-1.4_kp,log)

ONEPRIOR(osti,1._kp,4._kp,log)

ONEPRIOR(wrig,-3._kp,3._kp,log)
ONEPRIOR(wrio,1._kp,1._kp,flat)

!with xend

TWOPRIORS(rpi2,1._kp+epsilon(1._kp),1.2_kp,flat,0.5_kp,2._kp,log)

TWOPRIORS(ii,0._kp,10._kp,flat,-1._kp,4._kp,log)
TWOPRIORS(iif,0.1_kp,1._kp,inv,-1._kp,4._kp,log)
TWOPRIORS(iilambda,0.1_kp,4._kp,log,-1._kp,4._kp,log)

!second param phiend/phi0
TWOPRIORS(twiA1,-4._kp,-1._kp,log,log10(2._kp),log10(20._kp),log)
TWOPRIORS(twiA2,-4._kp,-1._kp,log,log10(2._kp),log10(40._kp),log)
!first param is N=10^(-10) (phi0/Mp)^(-2)
TWOPRIORS(twiB1,1._kp,100._kp,flat,log10(2._kp),log10(20._kp),log)
TWOPRIORS(twiB2,1._kp,100._kp,flat,log10(2._kp),log10(40._kp),log)

TWOPRIORS(bsusybif,0._kp,0.3_kp,flat,-200._kp,-0._kp,flat)
TWOPRIORS(bsusybil,-3._kp,-1._kp,log,-200._kp,-0._kp,flat)

TWOPRIORS(csi,-5,-1,log,0,1,flat)

TWOPRIORS(cnci,-5,-1,log,1,10,flat)

TWOPRIORS(imi,1._kp,6._kp,flat,1._kp,100._kp,flat)
TWOPRIORS(imi1,1._kp,1._kp,flat,1._kp,100._kp,flat)
TWOPRIORS(imi2,2._kp,2._kp,flat,1._kp,100._kp,flat)
TWOPRIORS(imi3,3._kp,3._kp,flat,1._kp,100._kp,flat)
TWOPRIORS(imi4,4._kp,4._kp,flat,1._kp,100._kp,flat)
TWOPRIORS(imi5,5._kp,5._kp,flat,1._kp,100._kp,flat)
TWOPRIORS(imi6,6._kp,6._kp,flat,1._kp,100._kp,flat)




!witout xend

TWOPRIORS(cwil,4._kp*exp(1._kp),4._kp*exp(1._kp),flat,-5._kp,-3._kp,log)
TWOPRIORS(cwif,4._kp*exp(1._kp),4._kp*exp(1._kp),flat,5._kp*10._kp**(-5._kp),5._kp*10._kp**(-4._kp),flat)

TWOPRIORS(sfi,2._kp,10._kp,flat,-1._kp,2._kp,log)
TWOPRIORS(sfis,2._kp,10._kp,flat,-1._kp,0._kp,log)
TWOPRIORS(sfil,2._kp,10._kp,flat,0._kp,2._kp,log)
TWOPRIORS(sfi1,1._kp,1._kp,flat,1.1_kp,3._kp,log)
TWOPRIORS(sfi2,2._kp,2._kp,flat,-1._kp,2._kp,log)
TWOPRIORS(sfi2l,2._kp,2._kp,flat,0._kp,2._kp,log)
TWOPRIORS(sfi3,3._kp,3._kp,flat,-1._kp,2._kp,log)
TWOPRIORS(sfi3s,3._kp,3._kp,flat,-1._kp,0._kp,log)
TWOPRIORS(sfi3l,3._kp,3._kp,flat,0._kp,2._kp,log)
TWOPRIORS(sfi4,4._kp,4._kp,flat,-1._kp,2._kp,log)
TWOPRIORS(sfi4s,4._kp,4._kp,flat,-1._kp,0._kp,log)
TWOPRIORS(sfi4l,4._kp,4._kp,flat,0._kp,2._kp,log)

!first param is alpha/(beta V) and 2nd is V
TWOPRIORS(kmiii,0.2_kp,5._kp,flat,5._kp,7._kp,log)

TWOPRIORS(lmi1p,0.1_kp,1._kp,flat,0.01_kp,10._kp,flat)
TWOPRIORS(lmi1o,1._kp,6._kp,flat,-2._kp,2._kp,log)

!gmssmi requires aspic quad-precision, and lnR>0
TWOPRIORS(gmssmi,0.9_kp,1.1_kp,flat,-5._kp,5._kp,log)
TWOPRIORS(gmssmiopA,-28._kp,-23._kp,log,0.00002_kp,0.0002_kp,flat)
TWOPRIORS(gmssmiopB,-28._kp,-21.5_kp,log,0.00002_kp,0.0002_kp,flat)
TWOPRIORS(gmssmiomA,-28._kp,-23._kp,log,0.00002_kp,0.0002_kp,flat)
TWOPRIORS(gmssmiomB,-28._kp,-20._kp,log,0.00002_kp,0.0002_kp,flat)

!grippi requires aspic quad-precision lnR>0
TWOPRIORS(gripi,0.9_kp,1.1_kp,flat,-5._kp,5._kp,log)
TWOPRIORS(gripiopA,-15._kp,-10._kp,log,0.02_kp,0.2_kp,flat)
TWOPRIORS(gripiopB,-15._kp,-8._kp,log,0.02_kp,0.2_kp,flat)
TWOPRIORS(gripiomA,-15._kp,-10._kp,log,0.02_kp,0.2_kp,flat)
TWOPRIORS(gripiomB,-15._kp,-8._kp,log,0.02_kp,0.2_kp,flat)
TWOPRIORS(gripiS,0.98_kp,1.02_kp,flat,10._kp,50._kp,flat)


!TWOPRIORS(tiA,0.5_kp,0.6_kp,flat,-5._kp,-3._kp,log)
TWOPRIORS(tiB1,0.5_kp,0.5_kp+10._kp**(-7._kp),flat,-5._kp,-3._kp,log)
TWOPRIORS(tiB2,0.5_kp,0.5_kp+2.*10._kp**(-7._kp),flat,-5._kp,-3._kp,log)
TWOPRIORS(tiB3,0.5_kp,0.5_kp+10._kp**(-6._kp),flat,-5._kp,-3._kp,log)
!TWOPRIORS(tiC,0.4_kp,0.5_kp,flat,-5._kp,-3._kp,log)
TWOPRIORS(tiD1,0.5_kp-10._kp**(-7._kp),0.5_kp,flat,-5._kp,-3._kp,log)
TWOPRIORS(tiD2,0.5_kp-10._kp**(-6._kp),0.5_kp,flat,-5._kp,-3._kp,log)
TWOPRIORS(tiD3,0.5_kp-10._kp**(-5._kp),0.5_kp,flat,-5._kp,-3._kp,log)
!TWOPRIORS(tiE,0.4_kp,0.6_kp,flat,-5._kp,-3._kp,log)
TWOPRIORS(tiF1,0.5_kp-10._kp**(-7._kp),0.5_kp+10._kp**(-7._kp),flat,-5._kp,-3._kp,log)
TWOPRIORS(tiF2,0.5_kp-10._kp**(-6._kp),0.5_kp+10._kp**(-6._kp),flat,-5._kp,-3._kp,log)
TWOPRIORS(tiF3,0.5_kp-10._kp**(-5._kp),0.5_kp+10._kp**(-5._kp),flat,-5._kp,-3._kp,log)
TWOPRIORS(ti12,0.5_kp,0.5_kp,flat,-5._kp,-3._kp,log)

TWOPRIORS(bei,-3._kp,3._kp,log,-1.5_kp,3._kp,log)

!first param alpha/f**2
TWOPRIORS(psnift1,-5._kp,-1._kp,log,-2._kp,1._kp,log)
TWOPRIORS(psnift2,-3._kp,-1._kp,log,-2._kp,1._kp,log)
TWOPRIORS(psnift3,-2._kp,-1._kp,log,-2._kp,1._kp,log)
TWOPRIORS(psnift4,-1.5_kp,-1._kp,log,-2._kp,1._kp,log)

!first param is alpha + hard prior such that alpha/f^2 < 0.1
TWOPRIORS(psnioA,-7._kp,-1._kp,log,-2._kp,1._kp,log)
TWOPRIORS(psnioB,-5._kp,-1._kp,log,-2._kp,1._kp,log)
TWOPRIORS(psnioC,-3._kp,-1._kp,log,-2._kp,1._kp,log)


TWOPRIORS(nckip,-4._kp,-1._kp,log,0.02_kp,0.2_kp,flat)
TWOPRIORS(nckim,-4._kp,-1._kp,log,-0.1_kp,-0.02_kp,flat)

TWOPRIORS(oi,-3._kp,-1._kp,log,-3._kp,-1._kp,log)

TWOPRIORS(sbi,-5._kp,-2._kp,log,-4._kp,-1._kp,log)
!alpha is set to alphamin(beta)
TWOPRIORS(sbialphamin,1._kp,1._kp,flat,-4._kp,-1._kp,log)

TWOPRIORS(ssbi1,-3._kp,1._kp,log,-5._kp,1._kp,log)
TWOPRIORS(ssbi2,-5._kp,-1._kp,mlog,-5._kp,1._kp,mlog)
TWOPRIORS(ssbi3,-5._kp,1._kp,log,-5._kp,1._kp,mlog)
TWOPRIORS(ssbi4,-5._kp,1._kp,log,-5._kp,1._kp,mlog)
TWOPRIORS(ssbi5,-5._kp,-1._kp,mlog,-5._kp,-1._kp,log)
TWOPRIORS(ssbi6,-5._kp,1._kp,mlog,-5._kp,1._kp,log)

TWOPRIORS(ssbi1f,0.001,10._kp,flat,0.00001_kp,10._kp,flat)
TWOPRIORS(ssbi2f,-0.1_kp,-0.00001_kp,flat,-10._kp,-0.00001_kp,flat)
TWOPRIORS(ssbi3f,0.00001_kp,10._kp,flat,-10._kp,-0.00001_kp,flat)
TWOPRIORS(ssbi4f,0.00001,10._kp,flat,-10._kp,-0.00001_kp,flat)
TWOPRIORS(ssbi5f,-0.1_kp,-0.00001_kp,flat,0.00001,0.1_kp,flat)
TWOPRIORS(ssbi6f,-10._kp,-0.00001_kp,flat,0.00001,10._kp,flat)

TWOPRIORS(nfi1,0._kp,10._kp,flat,1._kp+epsilon(1._kp),10._kp,flat)
TWOPRIORS(nfi3n,-10._kp,-epsilon(1._kp),flat,epsilon(1._kp),1._kp-epsilon(1._kp),flat)
TWOPRIORS(nfi3p,epsilon(1._kp),10._kp,flat,-10._kp,-epsilon(1._kp),flat)

TWOPRIORS(saai,-3._kp,3._kp,log,1._kp/3._kp,2._kp,flat)
TWOPRIORS(sabi13,-3._kp,3._kp,log,1._kp/3._kp, 1._kp/3._kp,flat)
TWOPRIORS(sabi12,-3._kp,3._kp,log,0.5_kp,0.5_kp,flat)
TWOPRIORS(sabi1,-3._kp,3._kp,log,1._kp,1._kp,flat)
TWOPRIORS(sabi32,-3._kp,3._kp,log,1.5_kp,1.5_kp,flat)
TWOPRIORS(sabi2,-3._kp,3._kp,log,2._kp,2._kp,flat)


!with xend
THREEPRIORS(lmi2o,1.1_kp,6._kp,flat,-2._kp,2._kp,log,0._kp,2._kp,log)
THREEPRIORS(lmi2p,0.1_kp,0.99_kp,flat,0.01_kp,10._kp,flat,0._kp,2._kp,log)

!rmi1* requires quadruple precision for aspic for some of these
!two (four). Inflation can occur at the top of the potential with hyper-fine
!tuned conditions (and ns=1...)
THREEPRIORS(rmi1,0.01,0.2,flat,-2,0,log,exp(-1._kp),1._kp,flat)
THREEPRIORS(rmi1l,-3.,-1.,log,-2,0,log,exp(-1._kp),1._kp,flat)

THREEPRIORS(rmi2,0.01,0.2,flat,-2,0,log,1._kp,exp(1._kp),flat)
THREEPRIORS(rmi2l,-3,-1,log,-2,0,log,1._kp,exp(1._kp),flat)

THREEPRIORS(rmi3,-0.2,-0.01,flat,-2,0,log,exp(-1._kp),1._kp,flat)
THREEPRIORS(rmi3l,-3,-1,mlog,-2,0,log,exp(-1._kp),1._kp,flat)

!the last params is y such that xend = 1 + y(xendmax-1)
THREEPRIORS(rmi4,-0.2,-0.01,flat,-2,0,log,0._kp,1._kp,flat)
THREEPRIORS(rmi4l,-3,-1,mlog,-2,0,log,0._kp,1._kp,flat)

!last params is y with xend=xendmin+y*(xendmax-xendmin)
THREEPRIORS(vhi,1+epsilon(1._kp),6,flat,0,3,log,0,1,flat)

!vhi l1 requires quadruple precision for aspic
THREEPRIORS(vhil1,0,0.9,flat,0,3,log,0,1,flat)
THREEPRIORS(vhi12,0.5_kp,0.5_kp,flat,0,3,log,0,1,flat)
THREEPRIORS(vhi1,1,1,flat,0,3,log,0,1,flat)
THREEPRIORS(vhi2,2,2,flat,0,3,log,0,1,flat)
THREEPRIORS(vhi3,3,3,flat,0,3,log,0,1,flat)
THREEPRIORS(vhi4,4,4,flat,0,3,log,0,1,flat)


!second param is mu/mumax, last is y with xend=xendmin+y*(xendmax-xendmin)
THREEPRIORS(dsi,1+epsilon(1._kp),6,flat,-5,0,log,-5,-0.7,log)
!second param is mu/Mpl + hard prior mu<mumax, last param is y, as for
! dsi
THREEPRIORS(dsio,1+epsilon(1._kp),6,flat,1d-9,1d-6,flat,-5,-0.7,log)
THREEPRIORS(dsi2,2,2,flat,1d-9,1d-6,flat,-5,-0.7,log)

!last params is xend/xendmax
THREEPRIORS(cndi,-2._kp,-1._kp,log,1._kp,6._kp,flat,0._kp,1._kp,flat)

!last aprams is y such that xend = xendmin+y(xendmax-xendmin)
THREEPRIORS(nfi2,-10._kp,0._kp,flat,1._kp+epsilon(1._kp),10._kp,flat,0._kp,1._kp,flat)

!last aprams is y such that xend = xendmin+y(xendmax-xendmin)
THREEPRIORS(nfi4p,epsilon(1._kp),10._kp,flat,epsilon(1._kp),1._kp-epsilon(1._kp),flat,0._kp,1._kp,flat)
THREEPRIORS(nfi4n,-10._kp,-epsilon(1._kp),flat,-10._kp,0._kp,flat,0._kp,1._kp,flat)

!witout xend

THREEPRIORS(gmlfi,1,6,flat,1,6,flat,-5,1,log)
THREEPRIORS(gmlfi2313,2._kp/3._kp,2._kp/3._kp,flat,1._kp/3._kp,1._kp/3._kp,flat,-5,1,log)
THREEPRIORS(gmlfi2343,2._kp/3._kp,2._kp/3._kp,flat,4._kp/3._kp,4._kp/3._kp,flat,-5,1,log)
THREEPRIORS(gmlfi11,1,1,flat,1,1,flat,-5,1,log)
THREEPRIORS(gmlfi12,1,1,flat,2,2,flat,-5,1,log)
THREEPRIORS(gmlfi13,1,1,flat,3,3,flat,-5,1,log)
THREEPRIORS(gmlfi21,2,2,flat,1,1,flat,-5,1,log)
THREEPRIORS(gmlfi22,2,2,flat,2,2,flat,-5,1,log)
THREEPRIORS(gmlfi23,2,2,flat,3,3,flat,-5,1,log)
THREEPRIORS(gmlfi31,3,3,flat,1,1,flat,-5,1,log)
THREEPRIORS(gmlfi32,3,3,flat,2,2,flat,-5,1,log)
THREEPRIORS(gmlfi33,3,3,flat,3,3,flat,-5,1,log)

THREEPRIORS(lpi1,1,6,flat,1,6,flat,-3,3,log)
THREEPRIORS(lpi141,4,4,flat,1,1,flat,-3,3,log)
THREEPRIORS(lpi142,4,4,flat,2,2,flat,-3,3,log)
THREEPRIORS(lpi143,4,4,flat,3,3,flat,-3,3,log)


THREEPRIORS(lpi22,1,6,flat,2,2,flat,2,5,log)
THREEPRIORS(lpi24,1,6,flat,4,4,flat,2,5,log)
THREEPRIORS(lpi26,1,6,flat,6,6,flat,2,5,log)

THREEPRIORS(lpi32,1,6,flat,2,2,flat,2,5,log)
THREEPRIORS(lpi34,1,6,flat,4,4,flat,2,5,log)
THREEPRIORS(lpi36,1,6,flat,6,6,flat,2,5,log)


!with xend and xini

!bi/kklti within string scenarios. The parameters are p=4, mu<<Mpl, xstg and
!the last is phiuv=xuv*mu
FOURPRIORS(bistg,4._kp,4._kp,flat,-6._kp,log10(2._kp),log,0._kp,3._kp,log,-2._kp,log10(2._kp),log)
FOURPRIORS(kkltistg,4._kp,4._kp,flat,-6._kp,log10(2._kp),log,0._kp,3._kp,log,-2._kp,log10(2._kp),log)

!bi/kklti pheno, xstg and xuv have no role and they are fixed to
!unphysical values (xstg=-1,xuv=-1). Inflation ends at xeps1

FOURPRIORS(kklti,2._kp,10._kp,flat,-3._kp,3._kp,log,-1._kp,-1._kp,flat,-1._kp,-1._kp,flat)
FOURPRIORS(kkltis,2._kp,10._kp,flat,-3._kp,0._kp,log,-1._kp,-1._kp,flat,-1._kp,-1._kp,flat)

FOURPRIORS(bi,2._kp,10._kp,flat,-3._kp,3._kp,log,-1._kp,-1._kp,flat,-1._kp,-1._kp,flat)
FOURPRIORS(bis,2._kp,10._kp,flat,-3._kp,0._kp,log,-1._kp,-1._kp,flat,-1._kp,-1._kp,flat)
FOURPRIORS(bi1s,1._kp,1._kp,flat,-3._kp,0._kp,log,-1._kp,-1._kp,flat,-1._kp,-1._kp,flat)
FOURPRIORS(bi2s,2._kp,2._kp,flat,-3._kp,0._kp,log,-1._kp,-1._kp,flat,-1._kp,-1._kp,flat)
FOURPRIORS(bi3s,3._kp,3._kp,flat,-3._kp,0._kp,log,-1._kp,-1._kp,flat,-1._kp,-1._kp,flat)
FOURPRIORS(bi4s,4._kp,4._kp,flat,-3._kp,0._kp,log,-1._kp,-1._kp,flat,-1._kp,-1._kp,flat)
FOURPRIORS(bi5s,5._kp,5._kp,flat,-3._kp,0._kp,log,-1._kp,-1._kp,flat,-1._kp,-1._kp,flat)
FOURPRIORS(bi6s,6._kp,6._kp,flat,-3._kp,0._kp,log,-1._kp,-1._kp,flat,-1._kp,-1._kp,flat)


