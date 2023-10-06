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
!them. The logical "true" or "false" says if the model is
!non-minimally coupled to gravity, or not.

ZERO(si,.true.)

ONE(rchi,.false.)
ONE(lfi,.false.)
ONE(rcmi,.false.)
ONE(rcqi,.false.)
ONE(ni,.false.)
ONE(esi,.false.)
ONE(pli,.false.)
ONE(kmii,.false.)
ONE(hf1i,.false.)
ONE(li,.false.)
ONE(rpi1,.true.)
ONE(rpi3,.true.)
ONE(dwi,.false.)
ONE(mhi,.false.)
ONE(rgi,.false.)
ONE(mssmi,.false.)
ONE(ripi,.false.)
ONE(ai,.false.)
ONE(cnai,.false.)
ONE(cnbi,.false.)
ONE(osti,.false.)
ONE(wri,.false.)
ONE(saai,.false.)
ONE(ccsi1,.true.)
ONE(ccsi3,.true.)
ONE(pai,.false.)
ONE(sbki,.false.)
ONE(ahi,.false.)


TWO(cwi,.false.)
TWO(sfi,.false.)
TWO(kmiii,.false.)
TWO(lmi1,.false.)
TWO(gmssmi,.false.)
TWO(gripi,.false.)
TWO(ti,.false.)
TWO(bei,.false.)
TWO(psni,.false.)
TWO(ncki,.false.)
TWO(oi,.false.)
TWO(sbi,.false.)
TWO(ssbi1,.false.)
TWO(ssbi2,.false.)
TWO(ssbi3,.false.)
TWO(ssbi4,.false.)
TWO(ssbi5,.false.)
TWO(ssbi6,.false.)
TWO(nfi1,.false.)
TWO(nfi3,.false.)
TWO(vfmi,.false.)
TWO(hbi,.false.)
TWO(shi,.false.)
TWO(sabi,.false.)
TWO(sati,.false.)
TWO(fi,.false.)
TWO(hni1,.false.)
TWO(saii1,.false.)
TWO(saii2,.false.)
TWO(gdwi,.false.)

TWO_XEND(rpi2,.true.)
TWO_XEND(ii,.false.)
TWO_XEND(twi,.false.)
TWO_XEND(bsusybi,.false.)
TWO_XEND(csi,.false.)
TWO_XEND(cnci,.false.)
TWO_XEND(imi,.false.)
TWO_XEND(ccsi2,.true.)
TWO_XEND(sdi,.false.)



THREE(gmlfi,.false.)
THREE(lpi1,.false.)
THREE(lpi2,.false.)
THREE(lpi3,.false.)
THREE(ncli,.false.)

THREE_XEND(lmi2,.false.)
THREE_XEND(rmi1,.false.)
THREE_XEND(rmi2,.false.)
THREE_XEND(rmi3,.false.)
THREE_XEND(rmi4,.false.)
THREE_XEND(vhi,.false.)
THREE_XEND(dsi,.false.)
THREE_XEND(cndi,.false.)
THREE_XEND(nfi2,.false.)
THREE_XEND(nfi4,.false.)

THREE_XEND(bi,.false.)
THREE_XEND(kklti,.false.)
THREE_XEND(hni2,.false.)
