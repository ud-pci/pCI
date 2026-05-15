      function cre(kap1,k,kap2)
      implicit doubleprecision(a-h,o-z)
c
c    computes the relativistic reduced matrix element
c           (j1 // c(k) // j2)
c    defined by eq.(5.15),i.p.grant, advances in physics(1970),
c    vol.19,p.762.
c
c         kap1,kap2 are the kappa values corresponding
c         to j1,j2
c
c    the triangle conditions are tested by the 3j-coefficient
c    routine, clrx
c
c  subroutines called: clrx
c
c
      cre=0d0
      k1=iabs(kap1)
      k2=iabs(kap2)
      cre=2d0*dsqrt(dble(k1*k2))*clrx(kap1,k,kap2)
      if (mod(k1,2).eq.1) cre =-cre
      return
      end
