C     =================================================
      Subroutine DefEv (iunit,lmax)     !### changed 14/04/08
c     Valence energies in P subspace (Eps1) are set equal to:
c     Kval= 0...................... HF energies from DAT file (iunit);
c           1...................... HF energies of the first valence
c                                   orbital with the same l and j;
c           2...................... Energies are taken from the INPUT
c                                   file and do not depend on the QN n;
c     Valence energies in Q subspace (Eps2) are set equal to HF
c     energies from DAT file (unit=iunit).
C     - - - - - - - - - - - - - - - - - - - - - - - - -
      include "hfd.par"
      INCLUDE "conf.par"
       implicit real*8 (a-h,j,o-z)
       common /Ns/Ns/Nso/Nso/Nsv/Nsv/Nmax/Nmax/Ii/Ii
     >        /Kval/Kval/Kout/Kout
       common /Lj/Lj(IPs)/Ll/Ll(IPs)/Eps1/Eps1(IPs)/Eps2/Eps2(IPs)
       character*4 eval
       dimension ee(20),p(IP6),q(IP6)
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        rewind(10)                                !### search for Eval
 200    read(10,5,end=700) eval                   !#### part of the input
 5      format(A4)
        if (eval.NE.'Eval') goto 200
        read(10,15) Kval                          !### defines valence
15      format(5X,I3)
        if (Kval.EQ.2) then
          do l1=1,lmax+1
            if (l1.EQ.1) then
              read(10,*,err=710,end=710) l,ee(2*l1)
            else
              read(10,*,err=710,end=710) l,ee(2*l1-1),ee(2*l1)
            end if
            if (l.NE.l1-1) goto 710
          end do
        end if
        call ReadF (iunit,1,p,q,2)
        ns1=p(2)+0.1d0
        Nmax=min(Nmax,ns1)
        do ni=1,Nmax
          li=Ll(ni)
          lji=Lj(ni)
          call ReadF (iunit,ni+4,p,q,1)
          Eps2(ni)=-p(Ii+1)
          if (Kval.EQ.0) then
            Eps1(ni)=-p(Ii+1)
          end if
          if (Kval.EQ.1) then
            do 10 nk=Nso+1,Ns
              lk=Ll(nk)
              ljk=Lj(nk)
              if(li.NE.lk.OR.lji.NE.ljk) goto 10
              call ReadF (iunit,nk+4,p,q,1)
              Eps1(ni)=-p(Ii+1)
              goto 210
 10         end do
 210        continue
          end if
          if (Kval.EQ.2) then                                          	
            ie=2*li+1
            if (lji.GT.2*li) ie=ie+1
            Eps1(ni)=ee(ie)
          end if
        end do
        nmin=min(Nso+1,Nsv)
        if (Kout.GT.1) write(11,25) (i,Eps1(i),Eps2(i),i=Nmin,Nmax)
 25     format(5X,'Energies of valence orbitals in P and Q spaces:',
     >       /3(I4,2F8.3))
       Return
 700    write (*,*) 'DefEv: can not find "Eval" in INPUT'
       Stop
 710    write (*,*) 'DefEv: can not find energies for l=',l1-1
       Stop
      end