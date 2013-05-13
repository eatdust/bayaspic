! write down here the prior associated with the aspic models. Syntax
! is of the form: (p1min,p1max,maptype1,p2min,p2max,maptype2...). If
! the end of inflation is an additional parameter, it is assumed to be
! the last one. Maptype denotes the relation between the aspic
! parameters and the one we are sampling on, "flat" means they are
! identical, "log" means we are sampling on the log10, "ln" means we are
! sampling on the log neper. The prior bounds are assumed to be on
! the sampling parameters.

ONEPRIOR(rchi,-1,1,flat)
ONEPRIOR(lfi,0.1,5,flat)
ONEPRIOR(rcmi,-6,-3,flat)
ONEPRIOR(rcqi,-3,0,flat)
ONEPRIOR(ni,0,3,flat)
ONEPRIOR(esi,-3,1,flat)
ONEPRIOR(pli,-4,0,flat)

TWOPRIORS(sfi,1.5,10,flat,-0.2,3,log)
TWOPRIORS(ii,-1,1,flat,0,2,log)

THREEPRIORS(vhi,1.01,5,flat,-0.2,3,log,0,1,flat)
