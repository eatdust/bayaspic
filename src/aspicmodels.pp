!   This file is part of bayaspic
!
!   Copyright (C) 2013-2021 C. Ringeval
!
!   bayaspic is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   bayaspic is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with bayaspic.  If not, see <https://www.gnu.org/licenses/>.


!write down the list of aspic models having a *standard reheating
!history* sorted according to the total number of ASPIC parameters,
!append _XEND if the field value at the end of inflation is one of
!them. Add the directive NON_MINIMAL if the model has a non-minimal
!coupling to gravity (required to compute rhoEndInf in the Jordan
!Frame)
  
ZERO(si)


ONE(rchi)
ONE(lfi)
ONE(rcmi)
ONE(rcqi)
ONE(ni)
ONE(esi)
ONE(pli)
ONE(kmii)
ONE(hf1i)
ONE(li)
ONE(rpi1)
ONE(rpi3)
ONE(dwi)
ONE(mhi)
ONE(rgi)
ONE(mssmi)
ONE(ripi)
ONE(ai)
ONE(cnai)
ONE(cnbi)
ONE(osti)
ONE(wri)
ONE(saai)
ONE(ccsi1)
ONE(ccsi3)
ONE(pai)
ONE(sbki)
ONE(ahi)


TWO(cwi)
TWO(sfi)
TWO(kmiii)
TWO(lmi1)
TWO(gmssmi)
TWO(gripi)
TWO(ti)
TWO(bei)
TWO(psni)
TWO(ncki)
TWO(oi)
TWO(sbi)
TWO(ssbi1)
TWO(ssbi2)
TWO(ssbi3)
TWO(ssbi4)
TWO(ssbi5)
TWO(ssbi6)
TWO(nfi1)
TWO(nfi3)
TWO(vfmi)
TWO(hbi)
TWO(shi)
TWO(sabi)
TWO(sati)
TWO(fi)
TWO(hni1)
TWO(saii1)
TWO(saii2)
TWO(gdwi)

TWO_XEND(rpi2)
TWO_XEND(ii)
TWO_XEND(twi)
TWO_XEND(bsusybi)
TWO_XEND(csi)
TWO_XEND(cnci)
TWO_XEND(imi)
TWO_XEND(ccsi2)
TWO_XEND(sdi)



THREE(gmlfi)
THREE(lpi1)
THREE(lpi2)
THREE(lpi3)
THREE(ncli)
THREE(saiii1)
THREE(saiii2)
THREE(saiii3)


THREE_XEND(lmi2)
THREE_XEND(rmi1)
THREE_XEND(rmi2)
THREE_XEND(rmi3)
THREE_XEND(rmi4)
THREE_XEND(vhi)
THREE_XEND(dsi)
THREE_XEND(cndi)
THREE_XEND(nfi2)
THREE_XEND(nfi4)

THREE_XEND(bi)
THREE_XEND(kklti)
THREE_XEND(hni2)
