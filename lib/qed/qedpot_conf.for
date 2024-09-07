c       =================================================
        program qedpot
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       This code compute self-energy model, vacuum
c       polarization Uehling and Wichmann-Kroll potentials
c       and store them to the file qedmod.dat
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       This code use following subroutines
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c  grid        - Subroutine initializes Semi-logarithmic
c                radial grid, using change of variables
c                rho(r)=al*r+bt*ln(r). The variable rho(r)
c                is given in the uninform grid
c                rho_i=rho_1+h*(i-1), where h is the grid
c                spacing, rho_1=al*r_1+bt*ln(r_1) and
c                r_1 > 0 is the first grid point closest
c                to the point of origin.
c  init_se -     This subroutine initialize some data needed to
c                constract Self-Energy (SE) potential. It also
c                compute the diagonal and non-diagonal SE matrix
c                elements using subroutine FSE_dat
c
c  FSE_dat  -    This subroutine computes diagonal and non-diagonal
c                self-energy (SE) matrix elements, using
c                interpolation procedure and tabulated data  taken
c                from the paper:
c                V.M. Shabaev, I.I. Tupitsyn and V. A. Yerokhin,
c                PRA, 88, 012513 (2013).
c
c   local_se_pot - This subroutine compute the local part
c                of the Self-Energy (SE) potential.
c                Local potential multiplied on r is store in the
c                array 'vsemi_loc(i,kappa)' for each value of the
c                relativistic quantum number kappa, where 'i' is
c                a number of grid point. 
c
c   nonlocal_se_pot -This subroutine compute the nonlocal
c                (separable) part of the Self-Energy (SE) potential.
c                This potential has the form:
c                V_nonloc = \sum_{a,b} |a> D_ab <b|, where a,b
c                are projected wave functions.
c                The projected wave functions and matrix D_ab are
c                stored in the file se_pot.dat
c
c    uehling  -  This subroutine computes Uehling vacuum-polarization
c                potential and store it to the file 'uehling.dat'.
C                The subroutine use accurate approximation method,
c                described in the paper:
c                L.W.Fullerton and G.A. Ringer, PRA 13, 1283 (1976)
c
c     wk  -      This subroutine computes Wichmann-Kroll point nuclear
c                vacuum-polarization potential and store it to the
c                file 'wk.dat'.
C                The subroutine use analytical approximation,
c                obtained in the paper:
c                A G Fainshtein, N L Manakov and A A Nekipelov,
c                J.Phys.E: At.Mol.Opt.Phys. v.23 (1990) 559-569. 
c
c       nucl  -  This subroutine computes nonpoint part of the
c                nuclear potential.
c
c       Dirac -  This subroutine computes radial Dirac
c                wavefunctions
c
c       tint  -  This subroutine computes the matrix elements
c                of r^k tint= <p,q|r^k|a,b>
c
c       sint  -  This subroutine computes the matrix elements of
c                the potential v(r). ds=<p,q|v|a,b>.
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Some variables:
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       z     ..... Nuclear charge
c       knucl ..... Nuclear model (0 - point, 1 - sphere,
c                   2 - Fermi model, 3 - Gauss)
c       cl    ..... Speed of light (in atomic units)
c       ii    ..... Number of points (ii=maxii-20), where
c                   maxii is given in the include file
c       h     ..... grid step
c       r(i)  ..... non-unform radial grid (i=1,ii)
c       p(i)  ..... Large component of the wavefunction
c       q(i)  ..... Small component of the wavefunction
c       y(i)  ..... Coulomb (or external) potential multiplied
c       unuc(i) ... Nonpoint part of the Nuclear potential
c                   multiplied by r(i)
c       uehl(i) ... Uehling potential multiplied by r(i).
c       wk(i)   ... Wichmann-Kroll potential multiplied by r(i).
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        parameter (maxii=4*512,maxns=5*128)
        parameter(maxff=2*(2*maxns+1))
        include "hfd.par"
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /cl/cl/z/z/ns/ns/kvar/kvar
        common /ii/ii/h/h/r1/r1/r2/r2/al/al/bt/bt
        common /p/p(maxii)/q/q(maxii)/a/a(maxii)/b/b(maxii)
        common /y/y(maxii)/cp/cp(maxii)/cq/cq(maxii)
        common /ca/ca(maxii)/cb/cb(maxii)
        common /c/c(maxii)/ww/ww(maxii)
        common /unuc/unuc(maxii)/dnuc/dnuc(maxii)
        common /r/r(maxii)/v/v(maxii)
        common /vloc/vloc(maxii)
        common /vsemi_loc/vsemi_loc(maxii,6)
     &  /wsemi_loc/wsemi_loc(maxii,6)
        common /uehl/uehl(maxii)/v_wk/v_wk(maxii)
        common /nmax/nmax/knucl/knucl/inucl/inucl
     1  /rnucl/rnucl
        common /niter/niter/mi/m1,m2,m3/nit/nit
        common /ns_proj/ns_proj
        common /nn/nn(maxns)/ll/ll(maxns)/jj/jj(maxns)
     &  /kk/kk(maxns)
        common /nn_proj/nn_proj(maxns)/ll_proj/ll_proj(maxns)
     &  /jj_proj/jj_proj(maxns)/kk_proj/kk_proj(maxns)
        common /uq/uq(maxii,maxns)/up/up(maxii,maxns)
        common /de_qed0/de_qed0(maxns)/ee_loc/ee_loc(maxns)
        common /qed_matrix/qed_matrix(maxns,maxns)
        common /fname/fname
        common /Dab/Dab(maxns,maxns)
        common /v_magn/v_magn(maxii)
        real*8 Gab(maxns,maxns),sab(maxns,maxns),vab(maxns,maxns)
        real*8 p_c(maxii),q_c(maxii),r_c(maxii),v_c(maxii),
     &  pq(4*maxii)
        real*8 v_lf(maxii),v_electr(maxii)
        real*8 ee_ue(maxns),ee_wk(maxns)
        real*8 ee_lf(maxns),ee_mn(maxns),ee_el(maxns)
        dimension IQN(4*maxii)
        equivalence (IQN(1),PQ(21))
        common /conf_dat/r1_c,r2_c,h_c,al_c,bt_c,rnucl_c,
     &  ns_c,ii_c,ispl,nbytes
        common /ff/ff(maxii*maxff)
        character*1 let(11),fname*12
        data
     1  let /'s','p','d','f','g','h','i','k','l','m','n'/
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        nmax=9
        knucl=2
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        cl=137.03599911d0
        pi=3.1415926535897932385d0
        alpha=1.d0/cl
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       file xx:qedpot.res
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        open(unit=11,file='qedpot.res',status='unknown')
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        kvar=1
        fname='CONF.DAT'
        open(unit=10,file='qedpot.inp',status='old',err=200)
c        read(10,*) z
c        read(10,*) knucl
        read(10,*,err=200,end=200) kvar
        read(10,*,err=200,end=200) fname
 200    close(unit=10)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        call conf_init(Ip6,maxii,p_c)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (Ns_c.gt.maxns) then
          write( *,5) ns_c,maxns
          write(11,5) ns_c,maxns
 5        format(/2x,'Error: Ns > maxns',/2x,'Ns=',i7,/2x,'Maxns=',i4,
     &    /2x,'Encrease maxns parmeter in the top of main program.') 
          stop
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        write( *,'(/2x,a,f10.4,/2x,a,i2,/2x,a,i2)') 
     &  'Z=',z,'Knucl=',knucl,'Kvar=',kvar
        write(11,'(/2x,a,f10.4,/2x,a,i2,/2x,a,i2)') 
     &  'Z=',z,'Knucl=',knucl,'Kvar=',kvar
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        call atom_data
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        alz=z/cl
        ns=ns_c
        rnucl=rnucl_c
        r2=r2_c
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        call grid(maxii,r,v)
        call nucl(maxii,r,unuc,dnuc)
        call init_se(maxns,nn_proj,ll_proj,jj_proj,
     &  kk_proj,qed_matrix)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Hydrogen like wave functions
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        write( *,15)
        write(11,15)
15      format(/2x,'Hydrogen-like wavefunctions')
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        write( *,25)
        write(11,25)
25      format(2x,47('=')/10x,'jj',8x,'ee',12x,'d(0)',
     &  6x,'Niter')
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do ni=1,ns_proj
          n=nn_proj(ni)
          l=ll_proj(ni)+1
          kappa=kk_proj(ni)
          call dirac(n,kappa,maxii,p,q,a,b,c,y,cp,cq,
     &    r,v,unuc,ww)
c
          e=0.5d0*p(ii+1)
          d=p(ii+5)**2+q(ii+5)**2
          write( *,35) ni,n,let(l),jj_proj(ni),e,d,niter
          write(11,35) ni,n,let(l),jj_proj(ni),e,d,niter
35        format (i3,i4,a1,i2,'/2',f16.8,e15.6,i5)
          call write_func(ni,p,q,2,maxii,ff)
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        write( *,45)
        write(11,45)
45      format(2x,47('='))
c       - - - - - - - - - - - - - - - - - - - - - - - - -
	write( *,'(/2x,a)') "Uehling potential..." 
        call uehling(maxii,r,unuc,uehl)
	write( *,'(2x,a)') "Wichmann-Kroll potential..." 
        call wk(maxii,r,v_wk)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        write( *,55)
        write(11,55)
55      format(/23x,'Ue[a.u.]',8x,'WK[a.u.]')
        do ni=1,ns_proj
          n=nn_proj(ni)
          kappa=kk_proj(ni)
          j=2*iabs(kappa)-1
          l=(iabs(2*kappa+1)-1)/2
          coef=alpha/pi*alz**4/n**3*cl**2
          call read_func(ni,p,q,2,maxii,ff)
          call sint(de_ue,uehl,p,q,p,q,r,v)
          call sint(de_wk,v_wk,p,q,p,q,r,v)
          write( *,65) ni,n,let(l+1),j,de_ue,de_wk
          write(11,65) ni,n,let(l+1),j,de_ue,de_wk
65        format(2x,i3,i5,a1,i2,'/2',2f16.8)
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do i=1,maxii
          v_magn(i)=0.d0
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (kvar.eq.5) then
          goto 210
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
	write( *,'(/2x,a)') 
     &  "Local part of Self-Energy (SE) potential..."     
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (kvar.eq.1) then
          call local_se_pot(maxns,nn_proj,ll_proj,
     &    jj_proj,kk_proj,ee_loc,maxii,p,q,y,r,v,vloc,ww,
     &    vsemi_loc,wsemi_loc,ff)
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (kvar.ne.1) then
          call local_pot_flambaum(maxns,ns_proj,nn_proj,
     &    ll_proj,jj_proj,kk_proj,maxii,p,q,cp,cq,y,unuc,
     &    dnuc,r,v,vloc,ww,vsemi_loc,wsemi_loc,v_lf,v_magn,
     &    v_electr,ee_lf,ee_mn,ee_el,ff)
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (kvar.ne.3) then
          write( *,'(/2x,a)') 
     &    "Nonlocal part of Self-Energy (SE) potential..." 
          call nonlocal_se_pot(maxii,maxns,nn_proj,
     &    ll_proj,jj_proj,kk_proj,de_qed0,dab,sab,vab,
     &    qed_matrix,p,q,a,b,cp,cq,r,v,vloc,v_magn,ww,
     &    vsemi_loc,wsemi_loc,up,uq,ff)
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
210     write( *,75)
        write(11,75)
75     format(/2x,'Falz=F(alpha*z)',1x,
     &  '[Johnson, Soff, At.Data.Nucl.Data Tables, 33, 405 (1985)]')
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        call print_var(kvar)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        write( *,125)
        write(11,125)
125     format(23x,'Ue_Falz',8x,'WK_Falz',8x,'SE_Falz',
     &  4x,'SE_Falz_Table')
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do 10 ni=1,ns_proj
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        n=nn_proj(ni)
        kappa=kk_proj(ni)
        if (kvar.eq.5.and.kappa.ne.-1) goto 10
        j=2*iabs(kappa)-1
        l=(iabs(2*kappa+1)-1)/2
        coef=alpha/pi*alz**4/n**3*cl**2
        call read_func(ni,p,q,2,maxii,ff)
        if (kvar.ne.5) then
          call qed_int(kappa,maxns,kk_proj,dab,maxii,
     &    p,q,p,q,r,v,cp,cq,ca,cb,uq,up,uehl,v_wk,vloc,
     &    vsemi_loc,v_magn,de_se,de_ue,de_wk)
        else
          ei=0.5d0*p(ii+1)
          call semi_empirical(ei,ei,de_se)
          call sint(de_ue,uehl,p,q,p,q,r,v)
          call sint(de_wk,v_wk,p,q,p,q,r,v)
        endif
        falz_ue=de_ue/coef
        falz_wk=de_wk/coef
        falz_se=de_se/coef
        fse_int=qed_matrix(ni,ni)/coef
        write( *,135) ni,n,let(l+1),j,falz_ue,falz_wk,
     &  falz_se,fse_int
        write(11,135) ni,n,let(l+1),j,falz_ue,falz_wk,
     &  falz_se,fse_int
135     format(2x,i3,i5,a1,i2,'/2',4f15.6)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
10      continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        call qed_test(kappa,maxns,kk_proj,dab,maxii,
     &  p,q,a,b,c,r,v,cp,cq,ca,cb,uq,up,uehl,v_wk,vloc,
     &  vsemi_loc,v_magn,de_se,de_ue,de_wk,y,unuc,ww)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        call conf_read(maxns,kk_proj,nn,ll,jj,kk,dab,
     &  maxii,p,q,a,b,cp,cq,ca,cb,uq,up,uehl,v_wk,
     &  vloc,vsemi_loc,v_magn,r,v,Ip6,p_c,q_c,r_c,v_c,
     &  pq,IQN,ff)
c
        call form_sgc
        call print_var(kvar)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        stop
        end
C       =================================================
        Subroutine conf_init(Ip6,maxii,p_c)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /z/z/ns/ns/fname/fname
        common /conf_dat/r1_c,r2_c,h_c,al_c,bt_c,rnucl_c,
     &  ns_c,ii_c,ispl,nbytes
        dimension p_c(maxii)
        character fname*12
C       - - - - - - - - - - - - - - - - - - - - - - - - -
        c1=0.01d0
C       - - - - - - - - - - - - - - - - - - - - - - - - -
        call recunit(nbytes)
        write(  *,5) nbytes
5       format(/4x,'Record unit =',i2,' bytes')
        lrec=IP6*(8/nbytes)
        lrec=20*(8/nbytes)
C       - - - - - - - - - - - - - - - - - - - - - - - - -
        open (12,file=fname,status='OLD',
     >        access='DIRECT',recl=lrec,err=700)
        read (12,rec=1) (p_c(i),i=1,20)
C       - - - - - - - - - - - - - - - - - - - - - - - - -
        z   = p_c(1)
        Ns_c= p_c(2)+c1
        ii_c= p_c(3)+c1
        r1_c= p_c(4)
        r2_c= p_c(5)
        h_c = p_c(6)
        bt_c= p_c(7)
        al_c= p_c(8)
        Rnucl_c=p_c(13)
        write( 6,15) Z,Ns_c,Ip6
        write(11,15) Z,Ns_c,Ip6
15      format (/2X,'Z  =',F6.2,/2X,'Ns =',I4,/2x,'Ip6=',i4)
c
        if (dabs(bt_c-1.d0).lt.1.d-6) then
          ispl=0
          write( 6,25) ii_c,r1_c,r2_c,h_c,al_c,bt_c
          write(11,25) ii_c,r1_c,r2_c,h_c,al_c,bt_c
25        format (/2x,'Conf Semi-logariphmic grid:',/2x,
     &    'II =',I5,/2x,'R1 =',e13.5,/2x,'R2 =',e13.5,
     &    /2x,'H  =',e13.5,/2x,'Al =',e13.5,/2x,'Bt =',e13.5)
        else
          ispl=1
          al_c=0.d0
          bt_c=1.d0
          write( 6,35) ii_c,r1_c,r2_c,h_c
          write(11,35) ii_c,r1_c,r2_c,h_c
35        format (/2x,'Conf Johnson grid:',/2x,
     &    'II =',I5,/2x,'R1 =',e13.5,/2x,'R2 =',e13.5,
     &    /2x,'H  =',e13.5)
        endif
C       - - - - - - - - - - - - - - - - - - - - - - - - -
        close(unit=12)
        return
C       - - - - - - - - - - - - - - - - - - - - - - - - -
 700    write( 6,105) fname
        write(11,105) fname
 105    format(/2X,'file ',A12,' is absent'/)
        Stop
C       - - - - - - - - - - - - - - - - - - - - - - - - -
        end
C       =================================================
        Subroutine conf_read(maxns,kk_proj,nn,ll,jj,kk,
     &  dab,maxii,p,q,a,b,cp,cq,ca,cb,uq,up,uehl,v_wk,
     &  vloc,vsemi_loc,v_magn,r,v,Ip6,p_c,q_c,r_c,v_c,pq,
     &  IQN,ff)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        parameter(IPs=128)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /ii/ii/r2/r2/cl/cl/al/al/bt/bt
        common /kvar/kvar/fname/fname
        integer kk_proj(maxns)
        integer nn(maxns),ll(maxns),jj(maxns),kk(maxns)
        real*8 Dab(maxns,maxns)
        real*8 p(*),q(*),a(*),b(*),cp(*),cq(*),ca(*),cb(*),
     &  r(*),v(*)
        real*8 uehl(maxii),v_wk(maxii)
        real*8 uq(maxii,maxns),up(maxii,maxns)
        real*8 vsemi_loc(maxii,6),vloc(maxii)
        common /Z/Z/Rnucl/Rnucl
        common /conf_dat/r1_c,r2_c,h_c,al_c,bt_c,rnucl_c,
     &  ns_c,ii_c,ispl,nbytes
c        dimension p_c(IP6),q_c(IP6),p1(IP6),q1(IP6),pq(4*IP6)
c        dimension r_c(IP6),v_c(IP6)
        dimension p_c(maxii),q_c(maxii),pq(4*maxii)
        dimension r_c(maxii),v_c(maxii)
        logical longbasis
        real*8 ff(*)
c        dimension IQN(4*1024)
        dimension IQN(4*1024)
c        equivalence (IQN(1),PQ(21))
c        equivalence (p_c(1),pq(1)), (q_c(1),pq(IP6+1)),
c     >      (p1(1),pq(2*IP6+1)), (q1(1),pq(3*IP6+1))
        character*1 let(11),fname*12
        data let/'S','P','D','F','G','H','I','K','L','M','N'/
C       - - - - - - - - - - - - - - - - - - - - - - - - -
        c1=0.01d0
        nmax=9
C       - - - - - - - - - - - - - - - - - - - - - - - - -
        lrec=IP6*(8/nbytes)
        open (12,file=fname,status='OLD',
     >        access='DIRECT',recl=lrec)
        read (12,rec=1) (p_c(i),i=1,Ip6)
        read (12,rec=2) (q_c(i),i=1,Ip6)
        read (12,rec=3) (r_c(i),i=1,Ip6)
        read (12,rec=4) (v_c(i),i=1,Ip6)
        do i=1,Ip6
          pq(i)=p_c(i)
          pq(i+Ip6)=q_c(i)
        enddo
        longbasis=dabs(p_c(20)-0.98765d0).LT.1.d-6
        read (12,rec=5) (p_c(i),i=1,Ip6)
        read (12,rec=6) (q_c(i),i=1,Ip6)
        do i=1,Ip6
          pq(i+2*Ip6)=p_c(i)
          pq(i+3*Ip6)=q_c(i)
        enddo
C       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (longbasis) then
          write( *,'(/2x,a)') 'Using variant for long basis '
          write(11,'(/2x,a)') 'Using variant for long basis '
c
          do ni=1,Ns_c
            Nn(ni)=IQN(4*ni-3)
            Ll(ni)=IQN(4*ni-2)
            Kk(ni)=IQN(4*ni-1)
            Jj(ni)=IQN(4*ni)
          end do          
        else
          if=20
          do ni=1,Ns_c
            if=if+1
            Nn(ni)=pq(if)+c1
            if=if+1
            Ll(ni)=pq(if)+c1
            if=if+3
            c2=dsign(c1,pq(if))
            Kk(ni)=pq(if)+c2
            if=if+1
            c2=dsign(c1,pq(if))
            Jj(ni)=pq(if)+c2
          end do
        endif
C       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (ispl.eq.1 ) then
          r0=v_c(1)-r_c(1)
          ro_0=dlog(v_c(1)/r0)
          d0=ro_0/h_c
          i0=d0+0.5d0
        endif
        write( *,'(/23x,a,9x,a)') 'ee[a.u.]','Norm'
        write(11,'(/23x,a,9x,a)') 'ee[a.u.]','Norm'
C       - - - - - - - - - - - - - - - - - - - - - - - - -
        do 10 ni=1,Ns_c
C       - - - - - - - - - - - - - - - - - - - - - - - - -
        read(12,rec=2*ni+7) (p_c(i),i=1,Ip6)
        read(12,rec=2*ni+8) (q_c(i),i=1,Ip6)
        if (ispl.eq.1 ) then
          p_c(ii_c+3)=ii_c
          q_c(ii_c+3)=ii_c
        endif
        imax_c=p_c(ii_c+3)+c1
        if (imax_c.le.0.or.imax_c.gt.ii_c) p_c(ii_c+3)=ii_c
        q_c(ii_c+3)=p_c(ii_c+3)
        q_c(ii_c+4)=p_c(ii_c+4)
        n=Nn(ni)
        l=Ll(ni)
        j=Jj(ni)
        e=p_c(ii_c+1)
        imax_c=p_c(ii_c+3)+c1
        do i=1,ii
          ri=r(i)
          if (ri.le.r_c(imax_c)) then
            if (ispl.eq.0) then
              p(i)=flagr(ri,p_c,r_c,ii_c,h_c,al_c,bt_c)
              q(i)=flagr(ri,q_c,r_c,ii_c,h_c,al_c,bt_c)
            else
              p(i)=flagr_spl(ri,p_c,r_c,ii_c,h_c,r0,i0)
              q(i)=flagr_spl(ri,q_c,r_c,ii_c,h_c,r0,i0)
            endif
            imax=i
          else
            p(i)=0.d0
            q(i)=0.d0
          endif
        enddo
        p(ii+1)=p_c(ii_c+1)
        p(ii+3)=imax
        p(ii+4)=p_c(ii_c+4)
        q(ii+1)=q_c(ii_c+1)
        q(ii+3)=imax
        q(ii+4)=q_c(ii_c+4)
        do m=0,nmax
          p(ii+5+m)=p_c(ii_c+5+m)
          q(ii+5+m)=q_c(ii_c+5+m)
        enddo
        dn=tint(0,p,q,p,q,r,v)
        dn=dsqrt(dn)
        write( *,'(i5,i4,a1,i3,a,f18.8,f14.8)') ni,n,
     &  let(l+1),j,'/2',e,dn
        write(11,'(i5,i4,a1,i3,a,f18.8,f14.8)') ni,n,
     &  let(l+1),j,'/2',e,dn
        do i=1,ii
          p(i)=p(i)/dn
          q(i)=q(i)/dn
        enddo
        do m=0,nmax
          i=ii+5+m
          p(i)=p(i)/dn
          q(i)=q(i)/dn
        enddo
C       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Test
C       - - - - - - - - - - - - - - - - - - - - - - - - -
        gam=p_c(ii_c+4)
        p0=p_c(ii_c+5)
        q0=q_c(ii_c+5)
        i1=1
        i2=3
        i3=5
        do i=1,10
          x1=r_c(i1)
          x2=r_c(i2)
          x3=r_c(i3)
          ri=r_c(i)
          f1=p_c(i1)/x1**gam-p0
          f2=p_c(i2)/x2**gam-p0
          f3=p_c(i3)/x3**gam-p0
          g1=q_c(i1)/x1**gam-q0
          g2=q_c(i2)/x2**gam-q0
          g3=q_c(i3)/x3**gam-q0
          ci1=ri*(ri-x2)*(ri-x3)/(x1*(x1-x2)*(x1-x3))
          ci2=ri*(ri-x1)*(ri-x3)/(x2*(x2-x1)*(x2-x3))
          ci3=ri*(ri-x1)*(ri-x2)/(x3*(x3-x1)*(x3-x2))
          fi=(p0+ci1*f1+ci2*f2+ci3*f3)*ri**gam
          gi=(q0+ci1*g1+ci2*g2+ci3*g3)*ri**gam
          ri=r_c(i)
          pi=0.d0
          qi=0.d0
          do m=0,nmax
            j=ii_c+5+m
            pi=pi+p_c(j)*(ri/r_c(1))**m
            qi=qi+q_c(j)*(ri/r_c(1))**m
          enddo
          pi=ri**gam*pi
          qi=ri**gam*qi
c         write( *,'(i5,3e17.8,2x,2e17.8)') i,ri,p_c(i),pi,q_c(i),qi
c         write( *,'(i5,3e17.8,2x,2e17.8)') i,ri,p_c(i),fi,q_c(i),gi
        enddo
        call write_func(ni,p,q,2,maxii,ff)
C       - - - - - - - - - - - - - - - - - - - - - - - - -
 10     continue
C       - - - - - - - - - - - - - - - - - - - - - - - - -
        close(unit=12)
C       - - - - - - - - - - - - - - - - - - - - - - - - -
        pi=3.1415926535897932385d0
        alpha=1.d0/cl
        alz=z/cl
        open(unit=13,file='qedpot.dat',status='unknown')
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        write( *,5)
        write(11,5)
 5      format(/15x,'SE [a.u.]',2x,'Uehling[a.u]',3x,
     &  'WK[a.u.]',3x,'F_SE(alZ)',3x,'F_U(alZ)',
     &  4x,'F_WK(alZ)')
C       - - - - - - - - - - - - - - - - - - - - - - - - -
        do 20 ni=1,Ns_c
        do 30 nj=ni,Ns_c
C       - - - - - - - - - - - - - - - - - - - - - - - - -
        qed_ij=0.d0
        if (kk(ni).ne.kk(nj)) goto 30
        if (ll(ni).gt.2) goto 200
        n=nn(ni)
        kappa=kk(ni)
        j=jj(ni)
        l=ll(ni)
        call read_func(ni,p,q,2,maxii,ff)
        call read_func(nj,a,b,2,maxii,ff)
        if (kvar.ne.5) then
          call qed_int(kappa,maxns,kk_proj,dab,maxii,
     &    p,q,a,b,r,v,cp,cq,ca,cb,uq,up,uehl,v_wk,vloc,
     &    vsemi_loc,v_magn,de_se,de_ue,de_wk)
        else
          de_se=0.d0
          if (kappa.eq.-1) then
            ei=0.5d0*p(ii+1)
            ej=0.5d0*a(ii+1)
            if (ei.gt.0.and.ej.gt.0)
     &      call semi_empirical(ei,ej,de_se)
            call sint(de_ue,uehl,p,q,a,b,r,v)
            call sint(de_wk,v_wk,p,q,a,b,r,v)
          endif
        endif
        qed_ij=de_se+de_ue+de_wk
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (nj.ne.ni) goto 200
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        coef=pi/alpha/alz**4*nn(ni)**3/cl**2
        fse=de_se*coef
        fue=de_ue*coef
        fwk=de_wk*coef
        l=ll(ni)
        write( *,15) ni,nn(ni),let(l+1),jj(ni),de_se,
     &  de_ue,de_wk,fse,fue,fwk
        write( 11,15) ni,nn(ni),let(l+1),jj(ni),de_se,
     &  de_ue,de_wk,fse,fue,fwk
15      format (i3,i3,a1,i2,'/2',2e13.4,e12.4,3e12.4)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        ds=tint(0,p,q,a,b,r,v)
c        call sint(de,uehl,p,q,a,b,r,v)
c        call sint(dw,v_wk,p,q,a,b,r,v)
c        dse=tint(0,p,q,cp,cq,r,v)
c
c        falz_ue=de/coef
c        falz_wk=dw/coef*cl**2
c        write( *,65) ni,n,let(l+1),j,falz_ue,falz_wk,
c     &  falz_se,fse_int
c        write(11,65) ni,n,let(l+1),j,falz_ue,falz_wk,
c     &  falz_se,fse_int
c65      format(2x,i3,i5,a1,i2,'/2',4f15.6)
C       - - - - - - - - - - - - - - - - - - - - - - - - -
 200    write(13,'(2i5,d16.6)') ni,nj,qed_ij
C       - - - - - - - - - - - - - - - - - - - - - - - - -
 30     continue
 20     continue
C       - - - - - - - - - - - - - - - - - - - - - - - - -
        close(unit=13)
        Return
C       - - - - - - - - - - - - - - - - - - - - - - - - -
        End
c       =================================================
        subroutine form_sgc
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
        parameter (maxii=4*512,maxns=5*128)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /kvar/kvar
        common /nn/nn(maxns)/ll/ll(maxns)/jj/jj(maxns)
     &  /kk/kk(maxns)
        character*3 str0,str1,str(9)*1,lbl*4
        character*128 string
        equivalence (str(1),str1)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        write(*,'(/2x,a)') 'Forming file SGCqed.CON...'
        open(unit=13,file='qedpot.dat',status='old',err=710)
        open(unit=14,file='SGC.CON',status='old',err=500)
        open(unit=15,file='SGCqed.CON',status='unknown')
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        read(14,'(a)') string
        read(string,15) Nmax,lmax,Kt,Kval,Khot,C_SMS,Klow,Kqed
 15     format(1X,' Nmax=',I3,'lmax=',I1,'  kt=',I2,'  kval=',I2,
     >         '  Khot=',I2,' C_SMS=',F9.6,' Klow=',I1,' Kqed=',I1)
        if (kqed.NE.0) then
          write(*,*) ' In SGC.CON kqed =',kqed
          write(*,*) ' File SGC.CON already includes QED corrections!'
          stop
        end if
        Kqed=kvar
        write(15,25) Nmax,lmax,Kt,Kval,Khot,C_SMS,Klow,Kqed
 25     format(1X,' Nmax=',I3,'lmax=',I1,'  kt=',I2,'  kval=',I2,
     >         '  Khot=',I2,' C_SMS=',F9.6,' Klow=',I1,' Kqed=',I1)

 100    read(14,35,err=300) irpt,ka,kb,sg,dsg,esg
 35     format(1X,I4,1X,I3,1X,I3,1X,F12.8,2(1X,E12.5))
c        write(*,35) irpt,ka,kb,sg,dsg,esg

 200    read(13,*) ma,mb,sgq
c        write(*,*) ma,mb,sgq
        if (ma.GT.ka) then
          write(*,45) ka,kb,ma,mb
 45       format(2x,'form_sgc error: looking for ME '2I4,
     &    ' got ',2I4)
          Stop
        end if
        if (ma.LT.ka) goto 200
        if (mb.LT.kb) goto 200
        if (ma.NE.ka.OR.mb.NE.kb) then
          write(*,55) ka,kb,ma,mb
 55       format(2x,'form_sgc error: expecting indeces '2I4,
     &    ' got ',2I3)
          Stop
        else
        if (sg.NE.0.d0) sg=sg+sgq     ! no QED corrections for zero MEs
          write(15,65) irpt,ka,kb,sg,dsg,esg
c          write(*, 65) irpt,ka,kb,sg,dsg,esg
 65       format(1X,I4,1X,I3,1X,I3,1X,F12.8,2(1X,E12.5))
          goto 100
        end if

 300    backspace(14)
        str0='scr'
        read(14,'(9a1)') (str(i),i=1,9)
        do i=1,9
          k=i
          if (str(i).NE.' ') goto 400
        end do
 400    do i=k,9
          str(i-k+1)=str(i)
        end do
        if (str1.NE.str0) then
          write(*,75) str0,str1
 75       format(2x,'form_sgc error: expecting string: 'a3,
     &    ' got: ',a3)
          Stop
        end if
        write(15,*)' screening coefficients:'
        do k=0,9
          read (14,85) k1,Chot
          write(15,85) k1,Chot
 85       format(5X,I2,5X,F8.5)
        end do
        write(15,*)' >>>>>>>>> END <<<<<<<<<'
        close(unit=13)
        close(unit=14)
        close(unit=15)
        write(*,'(2x,a)') 'File SGCqed.CON formed using SGC.CON'
        Return
c       - - - - - - - - - - - - - - - - - - - - - - - - -
 500    open(unit=14,file='MBPT.INP',status='old',err=720)
        open(unit=15,file='SGCqed.CON',status='unknown')
        write(*,*) ' SGC.CON is abscent, using MBPT.INP'
C     - - - - - - - - - - - - - - - - - - - - - - - - -
 600    read (14,95) lbl                         !### search for MBPT
 95     format(A4)                               !#### part of the input
        if (lbl.NE.'MBPT') goto 600
        call reads(14,string)
        read (string,*) Nso                      !### CI core
        call reads(14,string)
        read (string,*) Nsh                      !### defines SCF field
        call reads(14,string)
        read (string,*) nss                      !### last virtual shell
        call reads(14,string)
        read (string,*) nsv                      !### first vacant shell
        call reads(14,string)
        read (string,*) Nmax
        call reads(14,string)
        read (string,*) lmax
        read (14,*)
        call reads(14,string)
        read (string,*) kt1
        close(unit=14)

        write(6, 25) Nmax,lmax,Kt1,0,0,0.d0,0,kvar
        write(15,25) Nmax,lmax,Kt1,0,0,0.d0,0,kvar

        irpt=0
        do ka=Nso+1,Nmax
          na =nn(ka)
          la =ll(ka)
          ja =jj(ka)
          if (la.LE.lmax) then
            do kb=ka,Nmax
              nb =nn(kb)
              lb =ll(kb)
              jb =jj(kb)
              if (lb.EQ.la.AND.jb.EQ.ja) then
                irpt=irpt+1
 700            read(13,*) ma,mb,sgq
                if (ma.GT.ka) then
                  write(*,45) ka,kb,ma,mb
                  Stop
                end if
                if (ma.LT.ka) goto 700
                if (mb.LT.kb) goto 700
                if (ma.NE.ka.OR.mb.NE.kb) then
                  write(*,55) ka,kb,ma,mb
                  Stop
                else
                  write(15,65) irpt,ka,kb,sgq,0.d0,0.d0
                end if
              end if
            end do
          end if
        end do
        write(15,*)' screening coefficients:'
        do k=1,10
          write(15,85) (k-1),1.d0
        end do
        write(15,*)' >>>>>>>>> END <<<<<<<<<'
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        close(unit=13)
        close(unit=15)
        write(*,'(2x,a)') 'File SGCqed.CON formed using MBPT.INP'
        return
c       - - - - - - - - - - - - - - - - - - - - - - - - -
 710    write(*,'(2x,a)') 'form_sgc error: No file qedpot.dat'
        Stop
 720    write(*,'(2x,a)') 'form_sgc error: No files SGC.CON or MBPT.INP'
        Stop
        End
C       =================================================
        Subroutine reads(kan,string)
        character*128 string
        integer kan
C       - - - - - - - - - - - - - - - - - - - - - - - - -
        read(kan,'(a)') string
        string(1:5)='     '
        return
        end
c       =================================================
        subroutine print_var(kvar)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        write( *,5) kvar
        write(11,5) kvar
5       format(/2x,'Variant:',i2)
c
        if (kvar.eq.1) then
        write( *,15)
        write(11,15)
15      format(2x,'QEDMOD [V.M. Shabaev, I.I. Tupitsyn V.A. Yerokhin',
     &  /10x,'Phys.rev.A, 88, 012513 (2013) (2013)]')
        endif
c
        if (kvar.eq.2) then
        write( *,25)
        write(11,25)
25      format(2x,
     &  'Flambaum Local potential+QEDMOD non-local correction')
        endif
c
        if (kvar.eq.3) then
        write( *,35)
        write(11,35)
35      format(2x,'Local potential: [V.V. Flambaum and J.S.M. Ginges,',
     &  /20x,'Phys.rev.A, 72, 052115 (2005)]')
        endif
c
        if (kvar.eq.4) then
        write( *,45)
        write(11,45)
45      format(2x,'QEDPOT: [I.I. Tupitsyn and E.V. Berseneva',
     &  /11x,'Optics and Spectroscopy, 114, 682 (2013)]')
        endif
c
        if (kvar.eq.5) then
        write( *,55)
        write(11,55)
55      format(2x,'Semi-empirical method:')
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
1000    return
        end
c       =================================================
        subroutine local_se_pot(maxns,nn_proj,ll_proj,
     &  jj_proj,kk_proj,ee_loc,maxii,p,q,y,r,v,vloc,ww,
     &  vsemi_loc,wsemi_loc,ff)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       This subroutine compute the local part of the
c       Self-Energy (SE) potential.
c       Local potential multiplied on r is store in the
c       array 'vsemi_loc(i,kappa)' for each value of the
c       relativistic quantum number kappa, where 'i' is
c       a number of grid point. 
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /cl/cl/z/z/r1/r1/ii/ii/h/h
        common /ns_proj/ns_proj/num_kappa/num_kappa        
        common /kk_se/kk_se(6)/ll_se/ll_se(6)/nn_se/nn_se(6)
        common /se_loc/se_loc(6)
        common /nshell/nshell(-5:5,10)        
        integer nn_proj(maxns),ll_proj(maxns),jj_proj(maxns),
     &  kk_proj(maxns)
        real*8 ee_loc(maxns)
        real*8 p(maxii),q(maxii),y(maxii)
        real*8 r(maxii),v(maxii)
        real*8 vloc(maxii),ww(maxii)
        real*8 vsemi_loc(maxii,6),wsemi_loc(maxii,6)
        real*8 ff(*)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        parameter (pi=3.1415926535897932385d0)
        alpha=1.d0/cl
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do il=1,num_kappa
          do i=1,maxii
            vsemi_loc(i,il)=0.d0
            wsemi_loc(i,il)=0.d0
          enddo
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do 10 il=1,num_kappa
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        kappa=kk_se(il)
        l=ll_se(il)
        n=nn_se(il)
        ni=nshell(kappa,n)
        if (ni.eq.0) goto 10
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do i=1,ii
          ri=r(i)
          vloc(i)=dexp(-ri*cl)*ri
          ww(i)=  dexp(-ri*cl)*ri
        enddo
        do i=ii,1,-1
          if (dabs(vloc(i)).gt.1.d-40) goto 200
          imax=i-1
          vloc(i)=0.d0
        enddo
c
 200    vloc(ii+3)=imax
        vloc(ii+4)=1.d0
        ww(ii+3)=ii
        ww(ii+4)=1.d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        call read_func(ni,p,q,2,maxii,ff)
        call sint(de,vloc,p,q,p,q,r,v)
        vloc0=se_loc(il)/de
        call sint(de,ww,p,q,p,q,r,v)
        w0=1.d0/de
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do i=1,ii
          vloc(i)=vloc0*vloc(i)
          ww(i)=w0*ww(i)
        enddo
        do i=1,maxii
          vsemi_loc(i,il)=vloc(i)
          wsemi_loc(i,il)=ww(i)
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
10      continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
c       =================================================
        subroutine local_pot_flambaum(maxns,ns_proj,nn_proj,
     &  ll_proj,jj_proj,kk_proj,maxii,p,q,cp,cq,y,unuc,dnuc,
     &  r,v,vloc,ww,vsemi_loc,wsemi_loc,v_lf,v_magn,
     &  v_electr,ee_lf,ee_mn,ee_el,ff)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /cl/cl/z/z/r1/r1/ii/ii/h/h
        common /nmax/nmax/knucl/knucl
        common /gauss/wg2(128),tg2(128)
        real*8 p(maxii),q(maxii),cp(maxii),cq(maxii),
     &  y(maxii),r(maxii),v(maxii)
        integer nn(maxns),ll(maxns),jj(maxns),kk(maxns)
        integer nn_proj(maxns),ll_proj(maxns),jj_proj(maxns),
     &  kk_proj(maxns)
        real*8 unuc(maxii),dnuc(maxii)
        real*8 v_lf(maxii),v_magn(maxii),v_electr(maxii)
        real*8 vloc(maxii),ww(maxii)
        real*8 vsemi_loc(maxii,6),wsemi_loc(maxii,6)
        real*8 ff(*)
        real*8 ee_lf(maxns),ee_mn(maxns),ee_el(maxns)
        character*1 let(11)
        data
     1  let /'s','p','d','f','g','h','i','k','l','m','n'/
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        parameter (pi=3.1415926535897932385d0)
        alpha=1.d0/cl
        nodes=64
c       aB = 5.2917720859**(-11)
	par = 4*dlog(1/(z*alpha)+0.5d0)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
	write( *,'(2x,a)') "Low frequency term..."     
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (z.le.83.d0) then
          Bz = 0.074d0 + 0.35d0*z*alpha
        end if
        if (z.gt.83.d0.and.z.lt.92.d0) then
          Bz=0.095d0+0.35d0*z*alpha
        end if
        if(z.ge.92.d0)then
          Bz=0.411d0*(1.d0-0.010d0*(110.d0-z))
        end if
	do i=1,ii
          ri=r(i)
          v_lf(i) = Bz*z**4*alpha**5*cl**2*dexp(-Z*ri)*ri
        enddo
        v_lf(ii+3)=ii
        v_lf(ii+4)=1.d0
        call origin_coef(ii,maxii,v_lf,r)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do ni=1,ns_proj
          de=0.d0
          call read_func(ni,p,q,2,maxii,ff)
          call sint(de,v_lf,p,q,p,q,r,v)
          ee_lf(ni)=de
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - - 
        write( *,'(2x,a)') "Magnetic form factor..." 
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do i=1,ii
          ri=r(i)
          x=2*ri*cl
          if (dabs(x).lt.7000.d0) then
c           call ITIKB(x,ti1,tk1)
            call TKB(x,tk1)
            bk0=0.d0
            bk1=0.d0
	    if (x.lt.700) call IK01B(x,bi0,di0,bi1,di1,bk0,dk0,bk1,dk1)
c
c           call K1_my(X,dbk1)
c
c           tk2=x*(bk1-tk1)
c           v_m = alpha*Z/(4*pi*ri)*(alpha*(tk2 - 1.d0)/ri + 2*tk1)
c
            vn=z-unuc(i)
            if (knucl.eq.0) then
              v_m = alpha/(2*pi*ri)*z*(bk1-alpha/(2*ri))
            else
              v_m = alpha/(2*pi*ri)*((vn+dnuc(i))*(bk1-alpha/(2*ri))-
     &              dnuc(i)*tk1)
            endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            cp(i)= alpha*ri/(2*pi)*((-z+unuc(i))/ri)*
     &             (bk1-tk1-alpha/(2*ri))
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            v_magn(i) = v_m*ri
          else
            v_magn(i)=0.d0
          endif
        enddo
c
        v_magn(ii+3)=ii
        v_magn(ii+4)=1.d0
        call origin_coef(ii,maxii,v_magn,r)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c           Test
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        do i=2,665
c          ri=r(i)
c          vm1=(cp(i+1)-cp(i-1))/v(i)/(2*h)
c          v_m1=v_magn(i)/ri
c          v_magn(i)=vm1*ri
c          write(*,'(i5,3e14.6)') i,ri,v_m1,vm1
c        enddo
c        do m=0,nmax
c           i=ii+5+m
c           v_magn(i)=0.d0
c         enddo
c         pause
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do ni=1,ns_proj
          call read_func(ni,p,q,2,maxii,ff)
          call sint(de,v_magn,p,q,q,p,r,v)
c          ee_mn(ni)=kk(ni)*de
          ee_mn(ni)=-de
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        write( *,'(2x,a)') "Electric form factor..." 
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do i=1,ii
          ri=r(i)
       	  rm = ri*cl
          x=2*ri*cl
c
          ro_max= 50.d0
          ro_min=-100.d0
       	  Np3 = 800
	  h3 = (ro_max-ro_min)/Np3
	  v_ln3 = 0.d0
	  ai = ro_min
	  bi = ai + h3
	  do j = 1, Np3
	    call D01BAZ(ai, bi, nodes, wg2, tg2, IFAIL)
	    do in = 1, nodes
	       xi = tg2(in)
               v_ln3 = v_ln3 + wg2(in)*f_ln3(xi, 2*rm)
	    enddo
	    ai = bi
	    bi = ai + h3
          enddo	
c
c         call ITIKB(x,ti1,tk1)
          call TKB(x,tk1)
            bk0=0.d0
            bk1=0.d0
	    if (x.lt.700) call IK01B(x,bi0,di0,bi1,di1,bk0,dk0,bk1,dk1)
c
          tk2=x*(bk1-tk1)
          v_el=v_ln3 +par*(bk0-0.5d0*tk2)- 1.5d0 * bk0 + tk2
          xi=(z-80.d0)*alpha
          az=(1.071-1.976*xi**2-2.128d0*xi**3+0.169d0*xi**4)*rm
          az=az/(rm+0.07*(alpha*z)**2)
          v_electr(i) =az*v_el*alpha*(z-unuc(i))/(ri*pi)*ri
        enddo
c
         v_electr(ii+3)=ii
         v_electr(ii+4)=1.d0
         call origin_coef(ii,maxii,v_electr,r)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do ni=1,ns_proj
          call read_func(ni,p,q,2,maxii,ff)
          call sint(de,v_electr,p,q,p,q,r,v)
          ee_el(ni)=de
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Total potential
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do i=1,ii
          vloc(i)=v_electr(i)+v_lf(i)
        enddo
        vloc(ii+3)=ii
        vloc(ii+4)=0.d0
        call origin_coef(ii,maxii,vloc,r)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        write( *,15)
        write(11,15)
15      format(/2x,'Loc_pot:',3x,'Electr.[a.u.]',2x,
     &  'Low.fr[a.u.]',2x,'Magn.[a.u.]',3x,'Total[a.u.]')
        do ni=1,ns_proj
          n=nn_proj(ni)
          l=ll_proj(ni)
          j=jj_proj(ni)
          e_lf=ee_lf(ni)
          e_mn=ee_mn(ni)
          e_el=ee_el(ni)
          total=e_lf+e_mn+e_el
          write( *,25) n,let(l+1),j,e_el,e_lf,e_mn,total
          write(11,25) n,let(l+1),j,e_el,e_lf,e_mn,total
 25       format(i4,a1,i3,'/2',4f14.7)
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        write( *,35)
        write(11,35)
35      format(/2x,'Loc_pot:',4x,'FSE_Electr.',3x,'FSE_Low.fr',
     &  5x,'FSE_Magn.',6x,'FSE_Tot')
        do ni=1,ns_proj
          n=nn_proj(ni)
          l=ll_proj(ni)
          j=jj_proj(ni)
          coef=alpha/pi*(alpha*z)**4/n**3*cl**2
          f_lf=ee_lf(ni)/coef
          f_mn=ee_mn(ni)/coef
          f_el=ee_el(ni)/coef
          total=f_lf+f_mn+f_el
          write( *,45) n,let(l+1),j,f_el,f_lf,f_mn,total
          write(11,45) n,let(l+1),j,f_el,f_lf,f_mn,total
 45       format(i4,a1,i3,'/2',4f14.7)
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do i=1,maxii
          ww(i)=vloc(i)
        enddo
        do il=1,6
          do i=1,maxii
            vsemi_loc(i,il)=vloc(i)
            wsemi_loc(i,il)=ww(i)
          enddo
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
c       =================================================
	real*8 function f_ln3(ro,x)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
	implicit  real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       f_ln3 function with ln and exp for electric form
c       factor with substitution t^2 - 1 = dexp(ro) 
c       as we use DO1BAZ exp is to be added in equation
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        p=dexp(0.5d0*ro)
	t=dsqrt(p*p+1.d0)
	t2=t*t
c        d2=t2-1.d0
c        if (d2.le.1.d-40) d2=1.d-100
c	f_ln3 = dsqrt(d2)*(1-1.d0/(2*t2))/(2*t)*
c     &  ro*dexp(-t*x)
c	f_ln3 = dsqrt(d2)*(1-1.d0/(2*t2))/(2*t)*
c     &  dlog(d2)*dexp(-t*x)
	f_ln3 = ro*p*(1-1.d0/(2*t2))/(2*t)*dexp(-t*x)
c       - - - - - - - - - - - - - - - - - - - - - - - - -	
        return
	end
C       =======================================================
        SUBROUTINE TKB(X,TK)
C       =======================================================
C       Purpose: Integrate Bessel functions I0(t) and K0(t)
C                with respect to t from 0 to x
C       Input :  x  --- Upper limit of the integral ( x Ô 0 )
C       Output:  TI --- Integration of I0(t) from 0 to x
C                TK --- Integration of K0(t) from x to infty
C       =======================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PI=3.141592653589793D0
        IF (X.EQ.0.0D0) THEN
           TK=PI/2
        ELSE IF (X.LE.2.0D0) THEN
           T1=X/2.0D0
           T=T1*T1
           TK=((((((.116D-5*T+.2069D-4)*T+.62664D-3)*T
     &        +.01110118D0)*T+.11227902D0)*T+.50407836D0)*T
     &        +.84556868D0)*T1
cc
           T1=X/5.0D0
           T=T1*T1
           TI=((((((((.59434D-3*T+.4500642D-2)*T
     &        +.044686921D0)*T+.300704878D0)*T+1.471860153D0)
     &        *T+4.844024624D0)*T+9.765629849D0)*T
     &        +10.416666367D0)*T+5.0D0)*T1
cc
              TK=TK-DLOG(X/2.0D0)*TI

           TK=PI/2-TK
        ELSE IF (X.GT.2.0.AND.X.LE.4.0D0) THEN
           T=2.0D0/X
           TK=(((.0160395D0*T-.0781715D0)*T+.185984D0)*T
     &        -.3584641D0)*T+1.2494934D0
           TK=TK*DEXP(-X)/DSQRT(X)
        ELSE IF (X.GT.4.0.AND.X.LE.7.0D0) THEN
           T=4.0D0/X
           TK=(((((.37128D-2*T-.0158449D0)*T+.0320504D0)*T
     &        -.0481455D0)*T+.0787284D0)*T-.1958273D0)*T
     &        +1.2533141D0
           TK=TK*DEXP(-X)/DSQRT(X)
        ELSE
           T=7.0D0/X
           TK=(((((.33934D-3*T-.163271D-2)*T+.417454D-2)*T
     &        -.933944D-2)*T+.02576646D0)*T-.11190289D0)*T
     &        +1.25331414D0
           TK=TK*DEXP(-X)/DSQRT(X)
        ENDIF
        RETURN
        END
C       =========================================================
        SUBROUTINE K1_my(X,DBK1)
C       =========================================================
C       Purpose: Compute modified Bessel functions I0(x), I1(1),
C                K0(x) and K1(x), and their derivatives
C       Input :  x   --- Argument ( x Ô 0 )
C       Output:  BI1 --- I1(x)
C                BK1 --- K1(x)-1/x
C       =========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        IF (X.EQ.0.0D0) THEN
           BI1=0.0D0
           DBK1=0.d0
           RETURN
        ELSE IF (X.LE.3.75D0) THEN
           T=X/3.75D0
           T2=T*T
           BI1=X*((((((.00032411D0*T2+.00301532D0)*T2
     &         +.02658733D0)*T2+.15084934D0)*T2+.51498869D0)
     &         *T2+.87890594D0)*T2+.5D0)
        ELSE
           T=3.75D0/X
           BI1=((((((((-.420059D-2*T+.01787654D0)*T
     &         -.02895312D0)*T+.02282967D0)*T-.01031555D0)*T
     &         +.163801D-2)*T-.00362018D0)*T-.03988024D0)*T
     &         +.39894228D0)*DEXP(X)/DSQRT(X)
        ENDIF
        IF (X.LE.2.0D0) THEN
           T=X/2.0D0
           T2=T*T
           DBK1=((((((-.00004686D0*T2-.00110404D0)*T2
     &         -.01919402D0)*T2-.18156897D0)*T2-.67278579D0)
     &         *T2+.15443144D0)*T2)/X+BI1*DLOG(T)
        ELSE
           T=2.0D0/X
           T2=T*T
           DBK1=((((((-.00068245D0*T+.00325614D0)*T
     &         -.00780353D0)*T+.01504268D0)*T-.0365562D0)*T+
     &         .23498619D0)*T+1.25331414D0)*DEXP(-X)/DSQRT(X)-1.d0/X
        ENDIF
        RETURN
        END
c       =================================================
        SUBROUTINE IK01B(X,BI0,DI0,BI1,DI1,BK0,DK0,BK1,DK1)
C       =========================================================
C       Purpose: Compute modified Bessel functions I0(x), I1(1),
C                K0(x) and K1(x), and their derivatives
C       Input :  x   --- Argument ( x Ô 0 )
C       Output:  BI0 --- I0(x)
C                DI0 --- I0'(x)
C                BI1 --- I1(x)
C                DI1 --- I1'(x)
C                BK0 --- K0(x)
C                DK0 --- K0'(x)
C                BK1 --- K1(x)
C                DK1 --- K1'(x)
C       =========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        IF (X.EQ.0.0D0) THEN
           BI0=1.0D0
           BI1=0.0D0
           BK0=1.0D+300
           BK1=1.0D+300
           DI0=0.0D0
           DI1=0.5D0
           DK0=-1.0D+300
           DK1=-1.0D+300
           RETURN
        ELSE IF (X.LE.3.75D0) THEN
           T=X/3.75D0
           T2=T*T
           BI0=(((((.0045813D0*T2+.0360768D0)*T2+.2659732D0)
     &         *T2+1.2067492D0)*T2+3.0899424D0)*T2
     &         +3.5156229D0)*T2+1.0D0
           BI1=X*((((((.00032411D0*T2+.00301532D0)*T2
     &         +.02658733D0)*T2+.15084934D0)*T2+.51498869D0)
     &         *T2+.87890594D0)*T2+.5D0)
        ELSE
           T=3.75D0/X
           BI0=((((((((.00392377D0*T-.01647633D0)*T
     &         +.02635537D0)*T-.02057706D0)*T+.916281D-2)*T
     &         -.157565D-2)*T+.225319D-2)*T+.01328592D0)*T
     &         +.39894228D0)*DEXP(X)/DSQRT(X)
           BI1=((((((((-.420059D-2*T+.01787654D0)*T
     &         -.02895312D0)*T+.02282967D0)*T-.01031555D0)*T
     &         +.163801D-2)*T-.00362018D0)*T-.03988024D0)*T
     &         +.39894228D0)*DEXP(X)/DSQRT(X)
        ENDIF
        IF (X.LE.2.0D0) THEN
           T=X/2.0D0
           T2=T*T
           BK0=(((((.0000074D0*T2+.0001075D0)*T2+.00262698D0)
     &         *T2+.0348859D0)*T2+.23069756D0)*T2+.4227842D0)
     &         *T2-.57721566D0-BI0*DLOG(T)
           BK1=((((((-.00004686D0*T2-.00110404D0)*T2
     &         -.01919402D0)*T2-.18156897D0)*T2-.67278579D0)
     &         *T2+.15443144D0)*T2+1.0D0)/X+BI1*DLOG(T)
        ELSE
           T=2.0D0/X
           T2=T*T
           BK0=((((((.00053208D0*T-.0025154D0)*T+.00587872D0)
     &         *T-.01062446D0)*T+.02189568D0)*T-.07832358D0)
     &         *T+1.25331414D0)*DEXP(-X)/DSQRT(X)
           BK1=((((((-.00068245D0*T+.00325614D0)*T
     &         -.00780353D0)*T+.01504268D0)*T-.0365562D0)*T+
     &         .23498619D0)*T+1.25331414D0)*DEXP(-X)/DSQRT(X)
        ENDIF
        DI0=BI1
        DI1=BI0-BI1/X
        DK0=-BK1
        DK1=-BK0-BK1/X
        RETURN
        END
c       =================================================
        SUBROUTINE D01BAZ(A, B, NPTS, WEIGHT, ABSCIS, IFAIL)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       podprogramma iz NAG, vychislyaet uzly polinomov Legeandr'a
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        DOUBLE PRECISION A, B
        INTEGER IFAIL, NPTS
        DOUBLE PRECISION ABSCIS(NPTS), WEIGHT(NPTS)
        DOUBLE PRECISION HFRNGE, PNTMID
        INTEGER I, IIJJ, N, NL, NN, NPTSA
        DOUBLE PRECISION ABST(136), WTST(136)
        INTEGER NSTOR(16)
        DATA WTST( 1),WTST( 2),WTST( 3),WTST( 4),WTST( 5)/+
     *0.200000000000000000000000000000D   1,+
     *0.100000000000000000000000000000D   1,+
     *0.555555555555555555555555555555D   0,+
     *0.888888888888888888888888888888D   0,+
     *0.347854845137453857373063949221D   0/
      DATA WTST( 6),WTST( 7),WTST(08),WTST(09),WTST(10)/+
     *0.652145154862546142626936050778D   0,+
     *0.236926885056189087514264040719D   0,+
     *0.478628670499366468041291514835D   0,+
     *0.568888888888888888888888888888D   0,+
     *0.171324492379170345040296142172D   0/
      DATA WTST(11),WTST(12),WTST(13),WTST(14),WTST(15)/+
     *0.360761573048138607569833513837D   0,+
     *0.467913934572691047389870343989D   0,+
     *0.101228536290376259152531354309D   0,+
     *0.222381034453374470544355994426D   0,+
     *0.313706645877887287337962201986D   0/
      DATA WTST(16),WTST(17),WTST(18),WTST(19),WTST(20)/+
     *0.362683783378361982965150449277D   0,+
     *0.666713443086881375935688098933D  -1,+
     *0.149451349150580593145776339657D   0,+
     *0.219086362515982043995534934228D   0,+
     *0.269266719309996355091226921569D   0/
      DATA WTST(21),WTST(22),WTST(23),WTST(24),WTST(25)/+
     *0.295524224714752870173892994651D   0,+
     *0.471753363865118271946159614850D  -1,+
     *0.106939325995318430960254718193D   0,+
     *0.160078328543346226334652529543D   0,+
     *0.203167426723065921749064455809D   0/
      DATA WTST(26),WTST(27),WTST(28),WTST(29),WTST(30)/+
     *0.233492536538354808760849898924D   0,+
     *0.249147045813402785000562436042D   0,+
     *0.351194603317518630318328761381D  -1,+
     *0.801580871597602098056332770628D  -1,+
     *0.121518570687903184689414809072D   0/
      DATA WTST(31),WTST(32),WTST(33),WTST(34),WTST(35)/+
     *0.157203167158193534569601938623D   0,+
     *0.185538397477937813741716590125D   0,+
     *0.205198463721295603965924065661D   0,+
     *0.215263853463157790195876443316D   0,+
     *0.271524594117540948517805724560D  -1/
      DATA WTST(36),WTST(37),WTST(38),WTST(39),WTST(40)/+
     *0.622535239386478928628438369943D  -1,+
     *0.951585116824927848099251076022D  -1,+
     *0.124628971255533872052476282192D   0,+
     *0.149595988816576732081501730547D   0,+
     *0.169156519395002538189312079030D   0/
      DATA WTST(41),WTST(42),WTST(43),WTST(44),WTST(45)/+
     *0.182603415044923588866763667969D   0,+
     *0.189450610455068496285396723208D   0,+
     *0.176140071391521183118619623518D  -1,+
     *0.406014298003869413310399522749D  -1,+
     *0.626720483341090635695065351870D  -1/
      DATA WTST(46),WTST(47),WTST(48),WTST(49),WTST(50)/+
     *0.832767415767047487247581432220D  -1,+
     *0.101930119817240435036750135480D   0,+
     *0.118194531961518417312377377711D   0,+
     *0.131688638449176626898494499748D   0,+
     *0.142096109318382051329298325067D   0/
      DATA WTST(51),WTST(52),WTST(53),WTST(54),WTST(55)/+
     *0.149172986472603746787828737001D   0,+
     *0.152753387130725850698084331955D   0,+
     *0.123412297999871995468056670700D  -1,+
     *0.285313886289336631813078159518D  -1,+
     *0.442774388174198061686027482113D  -1/
      DATA WTST(56),WTST(57),WTST(58),WTST(59),WTST(60)/+
     *0.592985849154367807463677585001D  -1,+
     *0.733464814110803057340336152531D  -1,+
     *0.861901615319532759171852029837D  -1,+
     *0.976186521041138882698806644642D  -1,+
     *0.107444270115965634782577342446D   0/
      DATA WTST(61),WTST(62),WTST(63),WTST(64),WTST(65)/+
     *0.115505668053725601353344483906D   0,+
     *0.121670472927803391204463153476D   0,+
     *0.125837456346828296121375382511D   0,+
     *0.127938195346752156974056165224D   0,+
     *0.701861000947009660040706373885D  -2/
      DATA WTST(66),WTST(67),WTST(68),WTST(69),WTST(70)/+
     *0.162743947309056706051705622063D  -1,+
     *0.253920653092620594557525897892D  -1,+
     *0.342738629130214331026877322523D  -1,+
     *0.428358980222266806568786466061D  -1,+
     *0.509980592623761761961632446895D  -1/
      DATA WTST( 71),WTST( 72),WTST( 73),WTST( 74),WTST( 75)/+
     *0.586840934785355471452836373001D  -1,+
     *0.658222227763618468376500637069D  -1,+
     *0.723457941088485062253993564784D  -1,+
     *0.781938957870703064717409188283D  -1,+
     *0.833119242269467552221990746043D  -1/
      DATA WTST( 76),WTST( 77),WTST( 78),WTST( 79),WTST( 80)/+
     *0.876520930044038111427714627518D  -1,+
     *0.911738786957638847128685771116D  -1,+
     *0.938443990808045656391802376681D  -1,+
     *0.956387200792748594190820022041D  -1,+
     *0.965400885147278005667648300635D  -1/
      DATA WTST( 81),WTST( 82),WTST( 83),WTST( 84),WTST( 85)/+
     *0.315334605230583863267731154389D  -2,+
     *0.732755390127626210238397962178D  -2,+
     *0.114772345792345394895926676090D  -1,+
     *0.155793157229438487281769558344D  -1,+
     *0.196161604573555278144607196522D  -1/
      DATA WTST( 86),WTST( 87),WTST( 88),WTST( 89),WTST( 90)/+
     *0.235707608393243791405193013784D  -1,+
     *0.274265097083569482000738362625D  -1,+
     *0.311672278327980889020657568463D  -1,+
     *0.347772225647704388925485859638D  -1,+
     *0.382413510658307063172172565237D  -1/
      DATA WTST( 91),WTST( 92),WTST( 93),WTST( 94),WTST( 95)/+
     *0.415450829434647492140588223610D  -1,+
     *0.446745608566942804194485871258D  -1,+
     *0.476166584924904748259066234789D  -1,+
     *0.503590355538544749578076190878D  -1,+
     *0.528901894851936670955050562646D  -1/
      DATA WTST( 96),WTST( 97),WTST( 98),WTST( 99),WTST(100)/+
     *0.551995036999841628682034951916D  -1,+
     *0.572772921004032157051502346847D  -1,+
     *0.591148396983956357464748174335D  -1,+
     *0.607044391658938800529692320278D  -1,+
     *0.620394231598926639041977841375D  -1/
      DATA WTST(101),WTST(102),WTST(103),WTST(104),WTST(105)/+
     *0.631141922862540256571260227502D  -1,+
     *0.639242385846481866239062018255D  -1,+
     *0.644661644359500822065041936577D  -1,+
     *0.647376968126839225030249387365D  -1,+
     *0.178328072169643294729607914497D  -2/
      DATA WTST(106),WTST(107),WTST(108),WTST(109),WTST(110)/+
     *0.414703326056246763528753572855D  -2,+
     *0.650445796897836285611736039998D  -2,+
     *0.884675982636394772303091465973D  -2,+
     *0.111681394601311288185904930192D  -1,+
     *0.134630478967186425980607666859D  -1/
      DATA WTST(111),WTST(112),WTST(113),WTST(114),WTST(115)/+
     *0.157260304760247193219659952975D  -1,+
     *0.179517157756973430850453020011D  -1,+
     *0.201348231535302093723403167285D  -1,+
     *0.222701738083832541592983303841D  -1,+
     *0.243527025687108733381775504090D  -1/
      DATA WTST(116),WTST(117),WTST(118),WTST(119),WTST(120)/+
     *0.263774697150546586716917926252D  -1,+
     *0.283396726142594832275113052002D  -1,+
     *0.302346570724024788679740598195D  -1,+
     *0.320579283548515535854675043478D  -1,+
     *0.338051618371416093915654821107D  -1/
      DATA WTST(121),WTST(122),WTST(123),WTST(124),WTST(125)/+
     *0.354722132568823838106931467152D  -1,+
     *0.370551285402400460404151018095D  -1,+
     *0.385501531786156291289624969468D  -1,+
     *0.399537411327203413866569261283D  -1,+
     *0.412625632426235286101562974736D  -1/
      DATA WTST(126),WTST(127),WTST(128),WTST(129),WTST(130)/+
     *0.424735151236535890073397679088D  -1,+
     *0.435837245293234533768278609737D  -1,+
     *0.445905581637565630601347100309D  -1,+
     *0.454916279274181444797709969712D  -1,+
     *0.462847965813144172959532492322D  -1/
      DATA WTST(131),WTST(132),WTST(133),WTST(134),WTST(135)/+
     *0.469681828162100173253262857545D  -1,+
     *0.475401657148303086622822069442D  -1,+
     *0.479993885964583077281261798713D  -1,+
     *0.483447622348029571697695271580D  -1,+
     *0.485754674415034269347990667839D  -1/
      DATA WTST(136)/+0.486909570091397203833653907347D  -1/
      DATA ABST( 1),ABST( 2),ABST( 3),ABST( 4),ABST( 5)/+
     *0.000000000000000000000000000000D   0,+
     *0.577350269189625764509148780501D   0,+
     *0.774596669241483377035853079956D   0,+
     *0.000000000000000000000000000000D   0,+
     *0.861136311594052575223946488892D   0/
      DATA ABST( 6),ABST( 7),ABST( 8),ABST( 9),ABST(10)/+
     *0.339981043584856264802665759103D   0,+
     *0.906179845938663992797626878299D   0,+
     *0.538469310105683091036314420700D   0,+
     *0.000000000000000000000000000000D   0,+
     *0.932469514203152027812301554493D   0/
      DATA ABST(11),ABST(12),ABST(13),ABST(14),ABST(15)/+
     *0.661209386466264513661399595019D   0,+
     *0.238619186083196908630501721680D   0,+
     *0.960289856497536231683560868569D   0,+
     *0.796666477413626739591553936475D   0,+
     *0.525532409916328985817739049189D   0/
      DATA ABST(16),ABST(17),ABST(18),ABST(19),ABST(20)/+
     *0.183434642495649804939476142360D   0,+
     *0.973906528517171720077964012084D   0,+
     *0.865063366688984510732096688423D   0,+
     *0.679409568299024406234327365114D   0,+
     *0.433395394129247190799265943165D   0/
      DATA ABST(21),ABST(22),ABST(23),ABST(24),ABST(25)/+
     *0.148874338981631210884826001129D   0,+
     *0.981560634246719250690549090149D   0,+
     *0.904117256370474856678465866119D   0,+
     *0.769902674194304687036893833212D   0,+
     *0.587317954286617447296702418940D   0/
      DATA ABST(26),ABST(27),ABST(28),ABST(29),ABST(30)/+
     *0.367831498998180193752691536643D   0,+
     *0.125233408511468915472441369463D   0,+
     *0.986283808696812338841597266704D   0,+
     *0.928434883663573517336391139377D   0,+
     *0.827201315069764993189794742650D   0/
      DATA ABST(31),ABST(32),ABST(33),ABST(34),ABST(35)/+
     *0.687292904811685470148019803019D   0,+
     *0.515248636358154091965290718551D   0,+
     *0.319112368927889760435671824168D   0,+
     *0.108054948707343662066244650219D   0,+
     *0.989400934991649932596154173450D   0/
      DATA ABST(36),ABST(37),ABST(38),ABST(39),ABST(40)/+
     *0.944575023073232576077988415534D   0,+
     *0.865631202387831743880467897712D   0,+
     *0.755404408355003033895101194847D   0,+
     *0.617876244402643748446671764048D   0,+
     *0.458016777657227386342419442983D   0/
      DATA ABST(41),ABST(42),ABST(43),ABST(44),ABST(45)/+
     *0.281603550779258913230460501460D   0,+
     *0.950125098376374401853193354249D  -1,+
     *0.993128599185094924786122388471D   0,+
     *0.963971927277913791267666131197D   0,+
     *0.912234428251325905867752441203D   0/
      DATA ABST(46),ABST(47),ABST(48),ABST(49),ABST(50)/+
     *0.839116971822218823394529061701D   0,+
     *0.746331906460150792614305070355D   0,+
     *0.636053680726515025452836696226D   0,+
     *0.510867001950827098004364050955D   0,+
     *0.373706088715419560672548177024D   0/
      DATA ABST(51),ABST(52),ABST(53),ABST(54),ABST(55)/+
     *0.227785851141645078080496195368D   0,+
     *0.765265211334973337546404093988D  -1,+
     *0.995187219997021360179997409700D   0,+
     *0.974728555971309498198391993008D   0,+
     *0.938274552002732758523649001708D   0/
      DATA ABST(56),ABST(57),ABST(58),ABST(59),ABST(60)/+
     *0.886415527004401034213154341982D   0,+
     *0.820001985973902921953949872669D   0,+
     *0.740124191578554364243828103099D   0,+
     *0.648093651936975569252495786910D   0,+
     *0.545421471388839535658375617218D   0/
      DATA ABST(61),ABST(62),ABST(63),ABST(64),ABST(65)/+
     *0.433793507626045138487084231913D   0,+
     *0.315042679696163374386793291319D   0,+
     *0.191118867473616309158639820757D   0,+
     *0.640568928626056260850430826247D  -1,+
     *0.997263861849481563544981128665D   0/
      DATA ABST(66),ABST(67),ABST(68),ABST(69),ABST(70)/+
     *0.985611511545268335400175044630D   0,+
     *0.964762255587506430773811928118D   0,+
     *0.934906075937739689170919134835D   0,+
     *0.896321155766052123965307243719D   0,+
     *0.849367613732569970133693004967D   0/
      DATA ABST( 71),ABST( 72),ABST( 73),ABST( 74),ABST( 75)/+
     *0.794483795967942406963097298970D   0,+
     *0.732182118740289680387426665091D   0,+
     *0.663044266930215200975115168663D   0,+
     *0.587715757240762329040745476401D   0,+
     *0.506899908932229390023747474377D   0/
      DATA ABST( 76),ABST( 77),ABST( 78),ABST( 79),ABST( 80)/+
     *0.421351276130635345364119436172D   0,+
     *0.331868602282127649779916805730D   0,+
     *0.239287362252137074544603209165D   0,+
     *0.144471961582796493485186373598D   0,+
     *0.483076656877383162348125704405D  -1/
      DATA ABST( 81),ABST( 82),ABST( 83),ABST( 84),ABST( 85)/+
     *0.998771007252426118600541491563D   0,+
     *0.993530172266350757547928750849D   0,+
     *0.984124583722826857744583600026D   0,+
     *0.970591592546247250461411983800D   0,+
     *0.952987703160430860722960666025D   0/
      DATA ABST( 86),ABST( 87),ABST( 88),ABST( 89),ABST( 90)/+
     *0.931386690706554333114174380101D   0,+
     *0.905879136715569672822074835671D   0,+
     *0.876572020274247885905693554805D   0,+
     *0.843588261624393530711089844519D   0,+
     *0.807066204029442627082553043024D   0/
      DATA ABST( 91),ABST( 92),ABST( 93),ABST( 94),ABST( 95)/+
     *0.767159032515740339253855437522D   0,+
     *0.724034130923814654674482233493D   0,+
     *0.677872379632663905211851280675D   0,+
     *0.628867396776513623995164933069D   0,+
     *0.577224726083972703817809238540D   0/
      DATA ABST( 96),ABST( 97),ABST( 98),ABST( 99),ABST(100)/+
     *0.523160974722233033678225869137D   0,+
     *0.466902904750958404544928861650D   0,+
     *0.408686481990716729916225495814D   0,+
     *0.348755886292160738159817937270D   0,+
     *0.287362487355455576735886461316D   0/
      DATA ABST(101),ABST(102),ABST(103),ABST(104),ABST(105)/+
     *0.224763790394689061224865440174D   0,+
     *0.161222356068891718056437390783D   0,+
     *0.970046992094626989300539558536D  -1,+
     *0.323801709628693620333222431521D  -1,+
     *0.999305041735772139456905624345D   0/
      DATA ABST(106),ABST(107),ABST(108),ABST(109),ABST(110)/+
     *0.996340116771955279346924500676D   0,+
     *0.991013371476744320739382383443D   0,+
     *0.983336253884625956931299302156D   0,+
     *0.973326827789910963741853507352D   0,+
     *0.961008799652053718918614121897D   0/
      DATA ABST(111),ABST(112),ABST(113),ABST(114),ABST(115)/+
     *0.946411374858402816062481491347D   0,+
     *0.929569172131939575821490154559D   0,+
     *0.910522137078502805756380668008D   0,+
     *0.889315445995114105853404038272D   0,+
     *0.865999398154092819760783385070D   0/
      DATA ABST(116),ABST(117),ABST(118),ABST(119),ABST(120)/+
     *0.840629296252580362751691544695D   0,+
     *0.813265315122797559741923338086D   0,+
     *0.783972358943341407610220525213D   0,+
     *0.752819907260531896611863774885D   0,+
     *0.719881850171610826848940217831D   0/
      DATA ABST(121),ABST(122),ABST(123),ABST(124),ABST(125)/+
     *0.685236313054233242563558371031D   0,+
     *0.648965471254657339857761231993D   0,+
     *0.611155355172393250248852971018D   0,+
     *0.571895646202634034283878116659D   0,+
     *0.531279464019894545658013903544D   0/
      DATA ABST(126),ABST(127),ABST(128),ABST(129),ABST(130)/+
     *0.489403145707052957478526307021D   0,+
     *0.446366017253464087984947714758D   0,+
     *0.402270157963991603695766771260D   0,+
     *0.357220158337668115950442615046D   0,+
     *0.311322871990210956157512698560D   0/
      DATA ABST(131),ABST(132),ABST(133),ABST(134),ABST(135)/+
     *0.264687162208767416373964172510D   0,+
     *0.217423643740007084149648748988D   0,+
     *0.169644420423992818037313629748D   0,+
     *0.121462819296120554470376463492D   0,+
     *0.729931217877990394495429419403D  -1/
      DATA ABST(136)/+0.243502926634244325089558428537D  -1/
      DATA NSTOR(1), NSTOR(2), NSTOR(3), NSTOR(4) /1,2,3,4/
      DATA NSTOR(5), NSTOR(6), NSTOR(7), NSTOR(8) /5,6,8,10/
      DATA NSTOR(9), NSTOR(10), NSTOR(11), NSTOR(12) /12,14,16,20/
      DATA NSTOR(13), NSTOR(14), NSTOR(15), NSTOR(16) /24,32,48,64/
      DO 20 I=1,NPTS
         WEIGHT(I) = 0.0D0
         ABSCIS(I) = 0.0D0
   20 CONTINUE
      N = 0
      NPTSA = 0
      IFAIL = 0
      DO 60 I=1,16
         IF (NPTS.LT.NSTOR(I)) GO TO 80
         N = N + (NPTSA+1)/2
         NPTSA = NSTOR(I)
         IF (NPTS.EQ.NSTOR(I)) GO TO 100
   60 CONTINUE
   80 IFAIL = 1
  100 HFRNGE = 0.5D0*(B-A)
      PNTMID = 0.5D0*(A+B)
      NL = NPTSA/2
      IF (NL.LT.1) GO TO 140
      DO 120 NN=1,NL
         N = N + 1
         IIJJ = NPTSA + 1 - NN
         ABSCIS(NN) = HFRNGE*ABST(N) + PNTMID
         WEIGHT(NN) = HFRNGE*WTST(N)
         ABSCIS(IIJJ) = -HFRNGE*ABST(N) + PNTMID
         WEIGHT(IIJJ) = HFRNGE*WTST(N)
  120 CONTINUE
  140 IF (NPTSA.LE.(NL+NL)) GO TO 160
      N = N + 1
      ABSCIS(NL+1) = HFRNGE*ABST(N) + PNTMID
      WEIGHT(NL+1) = HFRNGE*WTST(N)
  160 RETURN
      END      
c       =================================================
        subroutine init_se(maxns,nn_proj,ll_proj,jj_proj,
     &  kk_proj,qed_matrix)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       This subroutine initialize some data needed to
c       constract Self-Energy (SE) potential. It also
c       compute the diagonal and non-diagonal SE matrix
c       elements using subroutine FSE_dat and store them
c       in 'qed_matrix' array.
c
c       The quantum numbers listed below define the quantum
c       numbers of the projected wave functions used by the
c       nonlocal part of SE model potential.
c
c       cl          - Speed of light (in atomic units)
c       ns_proj     - number of atomic shells
c       maxns       - maximal number of atomic shells
c                     (see file 'qedmod.inc')
c       ni          - shell number (ni=1,,,ns)
c       nn_proj(ni) - principal quantum numbers
c       ll_proj(ni) - orbital quantum numbers
c       jj_proj(ni) - angular quantum numbers
c       kk_proj(ni) - relativistic quantum numbers
c       num_kappa   - number of different relativistic
c                     quantum numbers
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /cl/cl/z/z/ns_proj/ns_proj/knucl/knucl
        common /num_kappa/num_kappa/kvar/kvar
        common /kk_se/kk_se(6)/ll_se/ll_se(6)/nn_se/nn_se(6)
        common /se_loc/se_loc(6)/num_loc/num_loc(-3:3)
        common /nshell/nshell(-5:5,10)
        common /nn_max/nn_max(0:3)
        integer nn_proj(maxns),ll_proj(maxns),jj_proj(maxns),
     &  kk_proj(maxns)
        real*8 qed_matrix(maxns,maxns)
        parameter (np=5)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        parameter (pi=3.1415926535897932385d0)
        alpha=1.d0/cl
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Default value of the max principal quantum number
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        nma=3
        if (kvar.eq.4) nma=4
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do l=0,3
          if (nma.gt.0) then
            nn_max(l)=nma
          else
            nn_max(l)=iabs(nma)+l
          endif
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        kk_se(1)=-1
        kk_se(2)=1
        kk_se(3)=-2
        kk_se(4)=2
        kk_se(5)=-3
        kk_se(6)=3
        num_kappa=5
        ni=0
        do il=1,num_kappa
          kappa=kk_se(il)
          l=(iabs(2*kappa+1)-1)/2
          ll_se(il)=l
          nn_se(il)=l+1
          do n=l+1,nn_max(l)
            ni=ni+1
            kk_proj(ni)=kappa
            ll_proj(ni)=l
            jj_proj(ni)=2*iabs(kappa)-1
            nn_proj(ni)=n
          enddo
        enddo
        ns_proj=ni
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        num_loc(-3)=5
        num_loc(-2)=3
        num_loc(-1)=1
        num_loc( 0)=0
        num_loc( 1)=2
        num_loc( 2)=4
        num_loc( 3)=6
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do kappa=-5,5
        do n=1,10
          nshell(kappa,n)=0
        enddo
        enddo
        do ni=1,ns_proj
          kappa=kk_proj(ni)
          n=nn_proj(ni)
          nshell(kappa,n)=ni
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Calculation of the qed_matrix
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        call se_pnt_stor()
        call se_fn_stor()
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do ni=1,ns_proj
        do nj=ni,ns_proj
          qed_matrix(ni,nj)=0.d0
          kappa=kk_proj(ni)
          nn1=nn_proj(ni)
          nn2=nn_proj(nj)
          if (kk_proj(nj).eq.kappa) then
            call qed_se_interpol(nn1,nn2,kappa,se_ij)
c            iz=z+0.5d0
c            call FSE_dat(kappa,nn1,nn2,iz,FSE_pnt,FSE_ext)
c            dn1=nn1
c            dn2=nn2
c            dn=dsqrt(dn1*dn2)
c            coef=alpha/pi*(z*alpha)**4/dn**3*cl**2
c            if (knucl.eq.0) then
c              qed_matrix(ni,nj)=fse_pnt*coef
c              qed_matrix(nj,ni)=fse_pnt*coef
c            else
c              qed_matrix(ni,nj)=fse_ext*coef
c              qed_matrix(nj,ni)=fse_ext*coef
               qed_matrix(ni,nj)=se_ij
               qed_matrix(nj,ni)=se_ij
            do il=1,num_kappa
              n=nn_se(il)
              if (kappa.eq.kk_se(il).and.nn1.eq.n.and.nn2.eq.n) then
                se_loc(il)=qed_matrix(ni,nj)
              endif
            enddo
          endif
        enddo
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
c       =================================================
        subroutine qed_se_interpol(nn1,nn2,kappa,se_ij)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /cl/cl/z/z/knucl/knucl
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        parameter (pi=3.1415926535897932385d0)
        alpha=1.d0/cl
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Calculation of the qed_matrix
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        call se_pnt_stor()
        call se_fn_stor()
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        se_ij=0.d0
        iz=z+0.5d0
        call FSE_dat(kappa,nn1,nn2,iz,FSE_pnt,FSE_ext)
        dn1=nn1
        dn2=nn2
        dn=dsqrt(dn1*dn2)
        coef=alpha/pi*(z*alpha)**4/dn**3*cl**2
        if (knucl.eq.0) then
          se_ij=fse_pnt*coef
        else
          se_ij=fse_ext*coef
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
c       =================================================
        subroutine nonlocal_se_pot(maxii,maxns,nn_proj,
     &  ll_proj,jj_proj,kk_proj,de_qed0,dab,sab,vab,
     &  qed_matrix,p,q,a,b,cp,cq,r,v,vloc,v_magn,ww,
     &  vsemi_loc,wsemi_loc,up,uq,ff)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       This subroutine compute the nonlocal (separable)
c       part of the Self-Energy (SE) potential. This
c       potential has the form:
c       V_nonloc = \sum_{a,b} |a> D_ab <b|, where a,b
c       are projected wave functions.
c       The projected wave functions and matrix D_ab are
c       stored in the file se_pot.dat
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /cl/cl/z/z/r1/r1/ns_proj/ns_proj
        common /ii/ii/h/h/al/al/bt/bt
        common /kvar/kvar
        common /num_loc/num_loc(-3:3)
        common /nshell/nshell(-5:5,10)
        integer nn_proj(maxns),ll_proj(maxns),jj_proj(maxns),
     &  kk_proj(maxns)
        real*8 de_qed0(maxns)
        real*8 Dab(maxns,maxns),Gab(maxns,maxns),
     &  sab(maxns,maxns),vab(maxns,maxns)
        real*8 qed_matrix(maxns,maxns)
        real*8 vloc(maxii),ww(maxii),v_magn(maxii)
        real*8 vsemi_loc(maxii,6),wsemi_loc(maxii,6)
        real*8 p(maxii),q(maxii),a(maxii),b(maxii),
     &  cp(maxii),cq(maxii),r(maxii),v(maxii)
        real*8 uq(maxii,maxns),up(maxii,maxns)
        integer iwrk1(maxns),iwrk2(maxns)
        real*8 ff(*)
        parameter (np=5)
        character*1 let(11)
        data
     1  let /'s','p','d','f','g','h','i','k','l','m','n'/
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        parameter (pi=3.1415926535897932385d0)
        alpha=1.d0/cl
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Functions up, uq
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do ni=1,ns_proj
          kappa=kk_proj(ni)
          n=nn_proj(ni)
          l=ll_proj(ni)
          call read_func(ni,p,q,2,maxii,ff)
          il=num_loc(kappa)
          do i=1,maxii
            vloc(i)=vsemi_loc(i,il)       
            ww(i)=wsemi_loc(i,il)
          enddo
          do i=ii+1,maxii
            a(i)=p(i)
            b(i)=q(i)
          enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
          if (kvar.eq.1.or.kvar.eq.2) then
            k=nn_proj(ni)-ll_proj(ni)
              do i=1,ii
                dw=dexp(-r(i)*z*2)*r(i)
                factor_p=(1.d0-(-1.d0)**k)/2.d0
                factor_q=(1.d0+(-1.d0)**k)/2.d0
                a(i)= factor_p*dw*p(i)
                b(i)= factor_q*dw*q(i)
            enddo
          endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Weight function does not include the magnetic part
c       - - - - - - - - - - - - - - - - - - - - - - - - -
          if (kvar.eq.4) then
            do i=1,ii
c             a(i)=vloc(i)*p(i)-v_magn(i)*q(i)
c             b(i)=vloc(i)*q(i)-v_magn(i)*p(i)
              a(i)=vloc(i)*p(i)
              b(i)=vloc(i)*q(i)
            enddo
          endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
          gam=p(ii+4)+ww(ii+4)
          a(ii+4)=gam
          b(ii+4)=gam
          a(ii+3)=ii
          b(ii+3)=ii
c
c
          do i=1,ii
            up(i,ni)=a(i)
            uq(i,ni)=b(i)
          enddo
          imax=ii
          do i=ii,1,-1
            if (dabs(up(i,ni))+dabs(uq(i,ni)).gt.1.d-40) goto 200
            imax=i-1
            up(i,ni)=0.d0
            uq(i,ni)=0.d0
          enddo
 200      up(ii+3,ni)=imax
          uq(ii+3,ni)=imax
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Matrix Sab and Vab
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do 10 ni=1,ns_proj
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        kappa=kk_proj(ni)
        il=num_loc(kappa)
        do i=1,maxii
          vloc(i)=vsemi_loc(i,il)
          ww(i)=wsemi_loc(i,il)
        enddo
        call read_func(ni,p,q,2,maxii,ff)
        do nj=1,ns_proj
          sab(ni,nj)=0.d0
          vab(ni,nj)=0.d0
          if (kk_proj(ni).eq.kk_proj(nj)) then
            call read_func(nj,cp,cq,2,maxii,ff)
            call sint(de,vloc,p,q,cp,cq,r,v)
            vab(ni,nj)=de
            if (kvar.eq.2.or.kvar.eq.3.or.kvar.eq.4) then
              call sint(dm,v_magn,p,q,cq,cp,r,v)
              vab(ni,nj)=vab(ni,nj)-dm
            endif
            do i=1,maxii
              a(i)=up(i,nj)
              b(i)=uq(i,nj)
            enddo
            sab(ni,nj)=tint(-1,p,q,a,b,r,v)
          endif
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
10      continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Matrix Inversion
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do ni=1,ns_proj
        do nj=1,ns_proj
          Gab(ni,nj)=Sab(ni,nj)
        enddo
        enddo
        call dminv(Gab,ns_proj,maxns,det,iwrk1,iwrk2)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Inversion test
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        del=0.d0
        do ni=1,ns_proj
        do nj=1,ns_proj
          aij=0.d0
          do nk=1,ns_proj
            aij=aij+Gab(ni,nk)*Sab(nk,nj)
          enddo
          if (ni.eq.nj) aij=aij-1.d0
          del=del+dabs(aij)/ns_proj
        enddo
        enddo
        write( *,'(/2x,a,e12.2)') 'Inversion test:',del
        write(11,'(/2x,a,e12.2)') 'Inversion test:',del
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       D-matrix
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        call dab_matrix(maxns,kk_proj,dab,gab,vab,sab,
     &  qed_matrix)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do ni=1, ns_proj
          delt=0.d0
          do nj=1, ns_proj
          do nk=1, ns_proj
            delt=delt+Sab(ni,nj)*Dab(nj,nk)*Sab(ni,nk)
          enddo
          enddo
          de_qed0(ni)=Vab(ni,ni)+delt
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
 1000   return
        end
c       =================================================
        subroutine dab_matrix(maxns,kk_proj,dab,gab,vab,
     &  sab,qed_matrix)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /kvar/kvar
        common /ns_proj/ns_proj
        integer kk_proj(maxns)
        real*8 Dab(maxns,maxns),Gab(maxns,maxns),
     &  vab(maxns,maxns),sab(maxns,maxns)
        real*8 qed_matrix(maxns,maxns)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       D-matrix
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do ni=1,ns_proj
        do nj=1,ns_proj
          dab(ni,nj)=0.d0
        enddo
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do 10 ni=1,ns_proj
        do 20 nj=1,ns_proj
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        Dab(ni,nj)=0.d0
        if (kk_proj(nj).ne.kk_proj(ni)) goto 20
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (kvar.eq.4) then
          di=qed_matrix(ni,ni)-vab(ni,ni)
          dj=qed_matrix(nj,nj)-vab(nj,nj)
          Dab(ni,nj)=0.5d0*(di/Sab(ni,ni)*Gab(ni,nj)+
     &    dj/Sab(nj,nj)*Gab(ni,nj))
          goto 20
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do 30 nk=1,ns_proj
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (kk_proj(nk).ne.kk_proj(ni)) goto 30
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do 40 nl=1,ns_proj
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (kk_proj(nl).ne.kk_proj(ni)) goto 40
        qkl=qed_matrix(nk,nl)
        qlk=qed_matrix(nl,nk)
        dkl=0.5d0*(qkl+qlk)-vab(nk,nl)
        Dab(ni,nj)=Dab(ni,nj)+Gab(ni,nk)*dkl*Gab(nj,nl)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
40      continue
30      continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
20      continue
10      continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
c       =================================================
        subroutine semi_empirical(ei,ej,de_se)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /ii/ii/z/z/cl/cl
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        pi=3.1415926535897932385d0
        alpha=1.d0/cl
        alz=z/cl
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        zi=z
        kappa=-1
        nn_max=5
        coef=alpha/pi*alz**4/nn_max**3*cl**2
        call qed_se_interpol(nn_max,nn_max,kappa,se_ij)
        falz_se=se_ij/coef
        de_se=(4.d0*ei*ej)**0.75d0*alpha**3*z**2*falz_se/pi/zi
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
c       =================================================
        subroutine qed_test(kappa,maxns,kk_proj,
     &  dab,maxii,p,q,a,b,c,r,v,cp,cq,ca,cb,uq,up,uehl,v_wk,
     &  vloc,vsemi_loc,v_magn,de_se,de_ue,de_wk,y,unuc,ww)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /ii/ii/z/z/cl/cl
        common /nmax/nmax/kvar/kvar
        common /ns_proj/ns_proj/num_kappa/num_kappa
        common /kk_se/kk_se(6)
        integer kk_proj(maxns)
        real*8 Dab(maxns,maxns)
        real*8 p(maxii),q(maxii),a(maxii),b(maxii),c(maxii)
        real*8 y(maxii),unuc(maxii),ww(maxii)
        real*8 r(maxii),v(maxii)
        real*8 cp(maxii),cq(maxii),ca(maxii),cb(maxii)
        real*8 uehl(maxii),v_wk(maxii)
        real*8 uq(maxii,maxns),up(maxii,maxns)
        real*8 vsemi_loc(maxii,6),vloc(maxii),v_magn(maxii)
        real*8 dproj(maxns)
        integer nn(20),kk(20)
        character*1 let(11)
        data
     1  let /'s','p','d','f','g','h','i','k','l','m','n'/
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        ns=10
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       4s-state
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        nn(1)=4
        kk(1)=-1
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       5s-state
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        nn(2)=5
        kk(2)=-1
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       4p_1/2-state
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        nn(3)=4
        kk(3)=1
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       5p_1/2-state
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        nn(4)=5
        kk(4)=1
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       4p_3/2-state
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        nn(5)=4
        kk(5)=-2
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       5p_3/2-state
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        nn(6)=5
        kk(6)=-2
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       4d_3/2-state
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        nn(7)=4
        kk(7)=2
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       5d_3/2-state
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        nn(8)=5
        kk(8)=2
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       4d_5/2-state
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        nn(9)=4
        kk(9)=-3
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       5d_5/2-state
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        nn(10)=5
        kk(10)=-3
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        pi=3.1415926535897932385d0
        alpha=1.d0/cl
        alz=z/cl
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Hydrogen like wave functions
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        write( *,15)
        write(11,15)
15      format(/4x,'QED Test:',/23x,'Se_Falz',4x,'SE_Falz_Table',
     &  3x,'Relative Error')
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do 10 ni=1,ns
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        n=nn(ni)
        kappa=kk(ni)
        if (kappa.ne.-1.and.kvar.eq.5) goto 10
        j=2*iabs(kappa)-1
          l=(iabs(2*kappa+1)-1)/2
          coef=alpha/pi*alz**4/n**3*cl**2
          call dirac(n,kappa,maxii,p,q,a,b,c,y,cp,cq,
     &    r,v,unuc,ww)
          e=0.5d0*p(ii+1)
          if (kvar.ne.5) then
            call qed_int(kappa,maxns,kk_proj,dab,maxii,
     &      p,q,p,q,r,v,cp,cq,ca,cb,uq,up,uehl,v_wk,vloc,
     &      vsemi_loc,v_magn,de_se,de_ue,de_wk)
          else
            ei=0.5d0*p(ii+1)
            call semi_empirical(ei,ei,de_se)
          endif
          call qed_se_interpol(n,n,kappa,se_ij)
          falz_se=de_se/coef
          fse_int=se_ij/coef
          del=(falz_se-fse_int)/fse_int
          write( *,135) ni,n,let(l+1),j,falz_se,fse_int,del
          write(11,135) ni,n,let(l+1),j,falz_se,fse_int,del
135       format(2x,i3,i5,a1,i2,'/2',2f15.6,e16.2)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
 10     continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
c       =================================================
        subroutine qed_int(kappa,maxns,kk_proj,dab,maxii,
     &  p,q,a,b,r,v,cp,cq,ca,cb,uq,up,uehl,v_wk,vloc,
     &  vsemi_loc,v_magn,de_se,de_ue,de_wk)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /ns/ns/ii/ii/z/z/cl/cl
        common /nmax/nmax/kvar/kvar
        common /ns_proj/ns_proj/num_kappa/num_kappa
        common /kk_se/kk_se(6)
        integer kk_proj(maxns)
        real*8 Dab(maxns,maxns)
        real*8 p(maxii),q(maxii),a(maxii),b(maxii)
        real*8 r(maxii),v(maxii)
        real*8 cp(maxii),cq(maxii),ca(maxii),cb(maxii)
        real*8 uehl(maxii),v_wk(maxii)
        real*8 uq(maxii,maxns),up(maxii,maxns)
        real*8 vsemi_loc(maxii,6),vloc(maxii),v_magn(maxii)
        real*8 dproj(maxns)
        character*1 let(11)
        data
     1  let /'s','p','d','f','g','h','i','k','l','m','n'/
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        call sint(de_ue,uehl,p,q,a,b,r,v)
        call sint(de_wk,v_wk,p,q,a,b,r,v)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        il=0
        do i=1,num_kappa
          if (kappa.eq.kk_se(i)) then
            il=i
          endif
        enddo
        if (il.gt.0) then
          do i=1,maxii
            vloc(i)=vsemi_loc(i,il)
          enddo
        else
          do i=1,maxii
            vloc(i)=0.d0
          enddo
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Local
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do i=1,ii
          ca(i)=vloc(i)*a(i)
          cb(i)=vloc(i)*b(i)
        enddo
        ca(ii+3)=ii
        cb(ii+3)=ii
        ca(ii+4)=vloc(ii+4)+a(ii+4)
        cb(ii+4)=vloc(ii+4)+b(ii+4)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        call origin_coef(ii,maxii,ca,r)
        call origin_coef(ii,maxii,cb,r)
        de_loc=tint(-1,p,q,ca,cb,r,v)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Magnetic term
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Small component has opposite sign.
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (kvar.ne.1.and.kvar.ne.5) then
          call sint(de_magn,v_magn,p,q,b,a,r,v)
          de_loc=de_loc-de_magn
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        de_unloc=0.d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (kvar.eq.3) goto 200
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Non-local potential
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do i=1,ii
          ca(i)=0.d0
          cb(i)=0.d0
        enddo
        do nj=1,ns_proj
          dproj(nj)=0.d0
          if (kk_proj(nj).eq.kappa) then
            do i=1,maxii
              cp(i)=up(i,nj)
              cq(i)=uq(i,nj)
            enddo
            dproj(nj)=tint(-1,a,b,cp,cq,r,v)
          endif
        enddo
c
        do nj=1,ns_proj
          if (kk_proj(nj).eq.kappa) then
            dnonloc=0.d0
            do nk=1,ns_proj
              if (kk_proj(nk).eq.kappa) then
                dnonloc=dnonloc+Dab(nj,nk)*dproj(nk)
              endif
            enddo
            do i=1,ii
              ca(i)=ca(i)+dnonloc*up(i,nj)
              cb(i)=cb(i)+dnonloc*uq(i,nj)
            enddo
            ca(ii+4)=vloc(ii+4)+up(ii+4,nj)
            cb(ii+4)=vloc(ii+4)+uq(ii+4,nj)
          endif
        enddo
c
        ca(ii+3)=ii
        cb(ii+3)=ii
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        call origin_coef(ii,maxii,ca,r)
        call origin_coef(ii,maxii,cb,r)
        de_unloc=tint(-1,p,q,ca,cb,r,v)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
 200    de_se=de_loc+de_unloc
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
c       =================================================
        subroutine uehling(maxii,r,unuc,uehl)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       This subroutine computes Uehling vacuum-polarization
c       potential and store it to the file 'uehling.dat'.
C       The subroutine use accurate approximation method,
c       described in the paper:
c       L.W.Fullerton and G.A. Ringer, PRA 13, 1283 (1976)
c       This approximation may be applied to point-charge
c       or to finite charge distribution nuclear models.
c
c       r(i)    .... Radial grid
c       unuc(i) .... Nonpoint part of the nuclear potential
c                    multiplied by r(i).
c       uehl(i) .... Uehling potential multiplied by r(i).
c       cl     ..... Speed of light (in atomic units)
c       z      ..... Nuclear charge
c       knucl  ..... Nuclear model (0 - point, 1 - uniform sphere,
c                    2 - Fermi model, 3 - Gauss)
c       ii     ..... Number of grid points
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /cl/cl/z/z/ii/ii/h/h/al/al/bt/bt
        common /inucl/inucl/knucl/knucl
        real*8 r(maxii),unuc(maxii),uehl(maxii)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        parameter (pi=3.1415926535897932385d0)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Calculation Uehling potential times on r
c       - - - - - - - - - - - - - - - - - - - - - - - - -
	Do i=1,maxii
          uehl(i)=0.D0
	End do
c       - - - - - - - - - - - - - - - - - - - - - - - - -
	coef1=-2.d0/3.d0/cl**2/(4*pi)
	coef2=-z*8.d0/3.d0/cl/(4*pi)
	ir=0
 500    ir=ir+1
        rr=r(ir)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (knucl.ne.0) then
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Finite charge nuclear distribution
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        npoints=maxii
        r0=r(1)
        rmax=r(inucl)
        tmax=dlog(rmax)
        tmin=dlog(r0)
        hi=(tmax-tmin)/(npoints-1)
        fint=0.d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do i=1,npoints
          t=hi*(i-1)+tmin
          ri=dexp(t)
          x1=2*cl*dabs(rr-ri)
          x2=2*cl*dabs(rr+ri)
c
          d1=dk0(x1)
          d2=dk0(x2)
c
          fu=d1-d2
          fr=fu*ro_nucl(ri)*ri
          fr=fr*ri
          fint=fint+fr*hi
          if (i.eq.1) then
            f0=fr
          endif
          if (i.eq.2) then
            f1=fr
          endif
          if (i.eq.2) then
            f2=fr
          endif
          if (i.eq.npoints-2) then
            fr1=fr
          endif
          if (i.eq.npoints-1) then
            fr2=fr
          endif
          if (i.eq.npoints) then
            fr3=fr
          endif
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        dfr=(3*fr3-4*fr2+fr1)/(2*hi)
        df0=(3*f2-4*f1+f0)/(2*hi)
        fint=fint-0.5d0*hi*(fr3+f0)-hi**2/12.d0*(dfr-df0)
c        write(*,'(5e16.5)') fint,fint1,fint0,
c     &  -hi**2/12.d0*dfr,hi**2/12.d0*df0
c
        x1=2*cl*dabs(rr-r0)
        x2=2*cl*dabs(rr+r0)
        d1=dk0(x1)
        d2=dk0(x2)
        fu=d1-d2
        fr=fu*ro_nucl(r0)*r0
        fint0=fr*r0/3.d0
        fint=fint+fint0
c
        Uehl(ir)=fint*coef1
c
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        else
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Point nuclear distribution
c       - - - - - - - - - - - - - - - - - - - - - - - - -
          x=2*cl*rr
          Uehl(ir)=dk1(x)*coef2
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
	If (dabs(uehl(ir)/r(ir)).gt.1.D-40.and.ir.lt.ii) go to 500
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        imax=ir
        uehl(ii+3)=imax
        uehl(ii+4)=0.d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        open (unit=15,file='uehling.dat',status='unknown')
        write(15,'(f22.14,16x,a)') z,'z'
        write(15,'(i22,16x,a)') ii,'ii'
        write(15,'(d22.14,16x,a)')  h,'h'
        write(15,'(d22.14,16x,a)') al,'al'
        write(15,'(d22.14,16x,a)') bt,'bt'
        write(15,*)
        write(15,'(17x,a,20x,a,17x,a,20x,a)') 'r','v_uehling',
     &  'unuc'
        do i=1,ii
          write (15,'(i6,3d24.14)') i,r(i),uehl(i),unuc(i)
        enddo
        write(15,*)
        write(15,'(d22.14,16x,a)') uehl(ii+3),'uehl(ii+3)'
        write(15,'(d22.14,16x,a)') uehl(ii+4),'uehl(ii+4)'
        close(unit=15)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
c       =================================================
        subroutine uehling1(maxii,r,unuc,uehl)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       This subroutine computes Uehling vacuum-polarization
c       potential and store it to the file 'uehling.dat'.
C       The subroutine use accurate approximations
c       of the Modified Bessel functions
c       This approximation may be applied to point-charge
c       or to finite charge distribution nuclear models.
c
c       r(i)    .... Radial grid
c       unuc(i) .... Nonpoint part of the nuclear potential
c                    multiplied by r(i).
c       uehl(i) .... Uehling potential multiplied by r(i).
c       cl     ..... Speed of light (in atomic units)
c       z      ..... Nuclear charge
c       knucl  ..... Nuclear model (0 - point, 1 - uniform sphere,
c                    2 - Fermi model)
c       ii     ..... Number of grid points
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /cl/cl/z/z/ii/ii/h/h/al/al/bt/bt
        common /inucl/inucl/knucl/knucl
        real*8 r(maxii),unuc(maxii),uehl(maxii)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        parameter (pi=3.1415926535897932385d0)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Calculation Uehling potential times on r
c       - - - - - - - - - - - - - - - - - - - - - - - - -
	Do i=1,maxii
          uehl(i)=0.D0
	End do
c       - - - - - - - - - - - - - - - - - - - - - - - - - 
        alpha=1.d0/cl
	ir=0
 500    ir=ir+1
        rr=r(ir)
        coef1=-2*alpha/(3*pi)
        coef2=-2*alpha/(3*pi)/cl/4.d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (knucl.eq.0) then
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Point nuclei
c       - - - - - - - - - - - - - - - - - - - - - - - - -
          x=2*rr*cl
          uehl(ir)=z*coef1*dk0_bessel(x)
        else
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Finite charge nuclear distribution
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        npoints=maxii
        r0=r(1)
        rmax=r(inucl)
        tmax=dlog(rmax)
        tmin=dlog(r0)
        hi=(tmax-tmin)/(npoints-1)
        fint=0.d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do i=1,npoints
c       - - - - - - - - - - - - - - - - - - - - - - - - -
          t=hi*(i-1)+tmin
          ri=dexp(t)
          x1=2*cl*dabs(rr-ri)
          x2=2*cl*dabs(rr+ri)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
          d1=dk1_bessel(x1)
          d2=dk1_bessel(x2)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
          fu=d1-d2
          fr=fu*ro_nucl(ri)*ri
          fr=fr*ri
          fint=fint+fr*hi
          if (i.eq.1) then
            f0=fr
          endif
          if (i.eq.2) then
            f1=fr
          endif
          if (i.eq.2) then
            f2=fr
          endif
          if (i.eq.npoints-2) then
            fr1=fr
          endif
          if (i.eq.npoints-1) then
            fr2=fr
          endif
          if (i.eq.npoints) then
            fr3=fr
          endif
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        dfr=(3*fr3-4*fr2+fr1)/(2*hi)
        df0=(3*f2-4*f1+f0)/(2*hi)
        fint=fint-0.5d0*hi*(fr3+f0)-hi**2/12.d0*(dfr-df0)
c
        x1=2*cl*dabs(rr-r0)
        x2=2*cl*dabs(rr+r0)
        d1=dk1_bessel(x1)
        d2=dk1_bessel(x2)
        fu=d1-d2
        fr=fu*ro_nucl(r0)*r0
        fint0=fr*r0/3.d0
        fint=fint+fint0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        Uehl(ir)=fint*coef2
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
	If (dabs(uehl(ir)/r(ir)).gt.1.D-40.and.ir.lt.ii) go to 500
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        imax=ir
        uehl(ii+3)=imax
        uehl(ii+4)=0.d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
c       =================================================
        real*8 function dk0(x)
        implicit real*8(a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Fullerton'a and Ringer'a, PRA 13, 1283 (1976)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        dk0=0.d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        a0= 0.88357293375d0
        a1=-0.28259817381d0
        a2=-0.58904879578d0
        a3= 0.12500133434d0
        a4=-0.032729913852d0
        a5= 0.0082888574511d0
        a6=-0.0000103277658d0
        a7= 0.0000636436689d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        b0=-319.999594323d0
        b1= 2.53900995981d0
        b2= 1.d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        c0=-319.999594333d0
        c1= 2.53901020662d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        d0=5.018065179d0
        d1=71.51891262d0
        d2=211.6209929d0
        d3=31.40327478d0
        d4=-1.d0 
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        e0=2.669207401d0
        e1=51.72549669d0
        e2=296.9809720d0
        e3=536.4324164d0
        e4=153.5335924d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if(x.ge.0.d0.and.x.le.1.d0) then       
          dp7=a0+a1*x+a2*x**2+a3*x**3+a4*x**4+a5*x**5
     #        +a6*x**6+a7*x**7
          dp21=(b0+b1*x**2+b2*x**4)/(c0+c1*x**2)
          if (x.gt.0) then
            dk0=dp7+x*dp21*dlog(x)
          else
            dk0=dp7
          endif
        else
          dk0=dexp(-x)*(d0+d1/x+d2/x**2+d3/x**3+d4/x**4)
     #        /(dsqrt(x))**3/(e0+e1/x+e2/x**2+e3/x**3+e4/x**4)
        end if
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
c       =================================================
        real*8 function dk1(x)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8(a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Fullerton'a and Ringer'a, PRA 13, 1283 (1976)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        dk1=0.d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        a0=-0.71740181754d0
        a1= 1.1780972274d0
        a2=-0.37499963087d0
        a3= 0.13089675530d0
        a4=-0.038258286439d0
        a5=-0.0000242972873d0
        a6=-0.0003592014867d0
        a7=-0.0000171700907d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        b0=-64.0514843293d0
        b1= 0.711722714285d0
        b2= 1.d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        c0= 64.0514843293d0
        c1= -0.711722686403d0
        c2=0.0008042207748
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        d0=217.2386409d0
        d1=1643.364528d0
        d2=2122.244512d0
        d3=-45.12004044d0
        d4=1.d0 
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        e0=115.5589983d0
        e1=1292.191441d0
        e2=3831.198012d0
        e3=2904.410075d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if(x.ge.0.d0.and.x.le.1.d0) then       
          dp7=a0+a1*x+a2*x**2+a3*x**3+a4*x**4+a5*x**5
     #        +a6*x**6+a7*x**7
          dp21=(b0+b1*x**2+b2*x**4)/(c0+c1*x**2+c2*x**4)
          dk1=dp7+dp21*dlog(x)
        else
          dk1=dexp(-x)*(d0+d1/x+d2/x**2+d3/x**3+d4/x**4)
     #        /(dsqrt(x))**3/(e0+e1/x+e2/x**2+e3/x**3)
        end if
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
c       =================================================
        real*8 function dk0_bessel(x)
        implicit real*8(a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (dabs(x).lt.700.d0) then
c	  call ITIKB(x,ti1,tk1)
          call TKB(x,tk1)
          bk0=0.d0
          bk1=0.d0
          if (x.lt.700) call IK01B(x,bi0,di0,bi1,di1,bk0,dk0,bk1,dk1)
c
          tk2=x*(bk1-tk1)
          tk3=(-x*tk2+tk1+x*bk0)/2.d0
          tk4=(-x*tk3+2*tk2+x*tk1)/3.d0
          dk0_bessel=bk0-0.5d0*tk2-0.5d0*tk4
        else
          dk0_bessel=0.d0
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
c       =================================================
        real*8 function dk1_bessel(x)
        implicit real*8(a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (dabs(x).lt.700.d0) then
c	  call ITIKB(x,ti1,tk1)
          call TKB(x,tk1)
          bk0=0.d0
          bk1=0.d0
          if (x.lt.700) call IK01B(x,bi0,di0,bi1,di1,bk0,dk0,bk1,dk1)
c
          tk2=x*(bk1-tk1)
          tk3=(-x*tk2+tk1+x*bk0)/2.d0
          tk4=(-x*tk3+2*tk2+x*tk1)/3.d0
          tk5=(-x*tk4+3*tk3+x*tk2)/4.d0
          dk1_bessel=tk1-0.5d0*tk3-0.5d0*tk5
        else
          dk1_bessel=0.d0
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
c       =================================================
        subroutine wk(maxii,r,v_wk)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       This subroutine computes Wichmann-Kroll point nuclear
c       vacuum-polarization potential and store it to the
c       file 'wk.dat'.
C       The subroutine use analytical approximation,
c       obtained in the paper:
c       A G Fainshtein, N L Manakov and A A Nekipelov,
c       J.Phys.E: At.Mol.Opt.Phys. v.23 (1990) 559-569. 
c
c       r(i)    .... Radial grid
c       wk(i)   .... Wichmann-Kroll potential multiplied by r(i).
c       cl     ..... Speed of light (in atomic units)
c       z      ..... Nuclear charge
c       ii     ..... Number of grid points
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /z/z/cl/cl/ii/ii/h/h/al/al/bt/bt
        real*8 r(maxii),v_wk(maxii)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Calculation Wichmann-Kroll potential times on r
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        call wk_init
c
	Do i=1,ii
          x=r(i)*cl
          v_wk(i)=u_wk(x)*r(i)*cl**2
	End do
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        imax=ii
        do i=ii,1,-1
	  If (dabs(v_wk(i)/r(i)).gt.1.D-40) go to 200
          imax=i
          v_wk(i)=0.d0
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
 200    v_wk(ii+3)=imax
        v_wk(ii+4)=0.d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        open (unit=15,file='wk.dat',status='unknown')
        write(15,'(f22.14,16x,a)') z,'z'
        write(15,'(i22,16x,a)') ii,'ii'
        write(15,'(d22.14,16x,a)')  h,'h'
        write(15,'(d22.14,16x,a)') al,'al'
        write(15,'(d22.14,16x,a)') bt,'bt'
        write(15,*)
        write(15,'(17x,a,23x,a,19x)') 'r','wk'
        do i=1,ii
          write (15,'(i6,3d24.14)') i,r(i),v_wk(i)
        enddo
        write(15,*)
        write(15,'(d22.14,16x,a)') v_wk(ii+3),'wk(ii+3)'
        write(15,'(d22.14,16x,a)') v_wk(ii+4),'wk(ii+4)'
        close(unit=15)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
c       =================================================
        subroutine wk_init
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       A G Fainshtein, N L Manakov and A A Nekipelov,
c       J.Phys.E: At.Mol.Opt.Phys. v.23 (1990) 559-569. 
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Initialization
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        real*8 f,p,q,a
        integer s
        common /wk_coeff/f(15,5),p(0:35),q(0:15),
     &  a(4,0:7,4),s(4,4)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do m=1,15        
        do n=1,5
          f(m,n) = 0.d0
        enddo
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        p(0)=  2.094 001 968 517 539 d-02
        p(1)= -6.328 682 056 509 706 d-02
        p(2)=  1.656 675 047 469 226 d-01
        p(3)=  6.218 254 052 794 603 d-02
        p(4)=  8.561 450 880 149 488 d-01
        p(5)= -1.324 099 530 525 646 d+00
        p(6)=  1.022 894 990 055 742 d-01
        p(7)=  1.818 708 880 809 046 d-01
        p(8)=  5.002 561 255 243 687 d-04
        p(9)=  2.167 778 174 657 628 d-01
        p(10)= 5.453 835 857 080 859 d-01
        p(11)= 3.137 663 571 113 079 d-01
        p(12)=-1.619 774 880 547 481 d-01
        p(13)=-6.749 903 602 516 955 d-02
        p(14)=-2.935 144 138 820 853 d-04
        p(15)= 9.755 800 734 714 044 d-02
        p(16)=-1.994 016 822 313 688 d-01
        p(17)=-4.777 623 414 220 549 d-02
        p(18)= 7.503 472 520 308 676 d-03
        p(19)= 3.440 742 164 594 281 d-05

        p(20)= 2.829 420 9807 d-03
        p(21)= 1.419 551 8237 d-04
        p(22)= 1.348 983 3978 d-05
        p(23)= 2.013 400 8556 d-07
        p(24)= 8.602 918 4536 d-06
        p(25)=-1.787 016 5724 d-05
        p(26)= 2.326 960 7856 d-05
        p(27)=-1.779 846 3760 d-05
        p(28)= 8.727 474 4996 d-06
        p(29)=-2.923 560 6110 d-06
        p(30)= 6.892 936 0802 d-07
        p(31)=-1.149 169 9126 d-07
        p(32)= 1.330 283 5807 d-08
        p(33)=-1.019 098 8753 d-09
        p(34)= 4.650 804 1566 d-11
        p(35)=-9.578 592 4366 d-13
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        q(0)  = 2.09400223235d-02
        q(1)  = 1.42425723499d-02
        q(2)  = 0.76419630941d-02
        q(3)  = 0.36124746165d-02
        q(4)  = 0.15118743043d-02
        q(5)  = 0.05247168246d-02
        q(6)  = 0.01797758276d-02
        q(7)  =-0.04761316284d-02
        q(8)  = 0.19713924639d-02
        q(9)  =-0.62421312791d-02
        q(10) = 1.29259217859d-02
        q(11) =-1.89507086023d-02
        q(12) = 1.89089237003d-02
        q(13) =-1.24611625668d-02
        q(14) = 0.48793792307d-02
        q(15) =-0.08858438257d-02
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        f(5,1) = 0.3812769d+00
        f(5,2) = 0.0266795d+00
        f(5,3) = 0.0000860d+00

        f(6,1) = 0.4946161d+00
        f(6,2) = 0.0458093d+00
        f(6,3) = 0.0002918d+00
        f(6,4) = 6.d-08

        f(7,1) = 0.6150374d+00
        f(7,2) = 0.0707397d+00
        f(7,3) = 0.0007216d+00
        f(7,4) = 5.d-07

        f(8,1) = 0.7412290d+00
        f(8,2) = 0.1014082d+00
        f(8,3) = 0.0014705d+00
        f(8,4) = 23.0d-07
        f(8,5) = 1.0d-10

        f(9,1) =0.8721898d+00
        f(9,2) =0.1376368d+00
        f(9,3) =0.0026288d+00
        f(9,4) =69.d-07
        f(9,5) =2.0d-09

        f(10,1) =1.0071432d+00
        f(10,2) =0.1791838d+00
        f(10,3) =0.0042776d+00
        f(10,4) =167.0d-07
        f(10,5) =8.0d-09

        f(11,1) = 1.1454776d+00
        f(11,2) = 0.2557770d+00
        f(11,3) = 0.0064866d+00
        f(11,4) = 346.0d-07
        f(11,5) = 3.0d-08


        f(12,1) = 1.2867042d+00
        f(12,2) = 0.2771342d+00
        f(12,3) = 0.0093134d+00
        f(12,4) = 643.0d-07
        f(12,5) = 9.0d-08

        f(13,1) = 1.4304276d+00
        f(13,1) = 0.3329753d+00
        f(13,1) = 0.0128043d+00
        f(13,1) = 1097.0d-07
        f(13,1) = 2.0d-07

        f(14,1) = 1.5763244d+00
        f(14,1) = 0.3930296d+00
        f(14,1) = 0.0169953d+00
        f(14,1) = 1750.0d-07
        f(14,1) = 5.0d-07

        f(15,1) = 1.7241274d+00
        f(15,1) = 0.4570401d+00
        f(15,1) = 0.0219132d+00
        f(15,1) = 2644.0d-07
        f(15,1) = 9.0d-07
c       - - - - - - - - - - - - - - - - - - - - - - - - -      
        do i=1,4
        do m=0,7
        do n=1,4
          a(i,m,n)=0.d0
        enddo
        enddo
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -     
        a(1,0,1) = 2.9675d+00
        a(1,0,2) = 7.8894d-01
        a(1,0,3) =-3.0201d-04
        a(1,0,4) = 4.9725d-04

        a(1,1,1) =-4.4288d+00
        a(1,1,2) =-2.9946d+01
        a(1,1,3) =-9.0441d-01
        a(1,1,4) = 2.4826d+00

        a(1,2,1) = 5.9339d01
        a(1,2,2) = 8.4881d+02
        a(1,2,3) = 5.8388d+01
        a(1,2,4) =-1.0745d+02

        a(1,3,1) =-5.6820d+02
        a(1,3,2) =-1.4888d+04
        a(1,3,3) =-1.5899d+03
        a(1,3,4) = 2.7128d+03

        a(1,4,1) = 4.6131d+03
        a(1,4,2) = 1.6003d+05
        a(1,4,3) = 2.2880d+04
        a(1,4,4) =-3.8249d+04

        a(1,5,1) =-2.6282d+04
        a(1,5,2) =-1.0259d+06
        a(1,5,3) =-1.8036d+05
        a(1,5,4) = 3.0100d+05

        a(1,6,1) = 8.9448d+04
        a(1,6,2) = 3.5980d+06
        a(1,6,3) = 7.3716d+05
        a(1,6,4) =-1.2368d+06

        a(1,7,1) =-1.3412d+05
        a(1,7,2) =-5.3067d+06
        a(1,7,3) =-1.2228d+06
        a(1,7,4) = 2.0678d+06
c       - - - - - - - - - - - - - - - - - - - - - - - - -     
        a(2,0,1) =-6.64673d-01
        a(2,0,2) = 4.68647d-05
        a(2,0,3) = 2.56839d-02

        a(2,1,1) = 5.70540d-01
        a(2,1,2) =-9.83836d-04
        a(2,1,3) =-3.21021d-03

        a(2,2,1) =-9.43955d-01
        a(2,2,2) = 8.28746d-03
        a(2,2,3) = 1.01163d-02

        a(2,3,1) = 1.18405d+00
        a(2,3,2) = 5.94372d-02
        a(2,3,3) =-1.62093d-02

        a(2,4,1) =-8.91689d-01
        a(2,4,2) =-5.31431d-02
        a(2,4,3) = 1.17762d-02

        a(2,5,1) = 3.64707d-01
        a(2,5,2) = 7.52653d-03
        a(2,5,3) =-4.69814d-03

        a(2,6,1) =-6.17902d-02
        a(2,6,2) =-4.54081d-04
        a(2,6,3) = 7.76469d-04
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        a(3,0,1) = 1.93743d-02
        a(3,0,2) = 1.64381d-01
        a(3,0,3) = 2.53250d-02

        a(3,1,1) =-6.43399d-01
        a(3,1,2) =-3.05149d-01
        a(3,1,3) = 3.56651d-03

        a(3,2,1) = 2.56565d-01
        a(3,2,2) = 2.31560d-01
        a(3,2,3) =-8.14980d-03

        a(3,3,1) =-1.09235d-01
        a(3,3,2) =-9.49948d-02
        a(3,3,3) = 6.21168d-03

        a(3,4,1) = 4.14226d-02
        a(3,4,2) = 2.22450d-02
        a(3,4,3) =-3.29527d-03

        a(3,5,1) =-7.28744d-03
        a(3,5,2) =-2.80168d-03
        a(3,5,3) = 7.60885d-04

        a(3,6,1) = 4.54914d-04
        a(3,6,2) = 1.47382d-04
        a(3,6,3) =-6.03993d-05
c       - - - - - - - - - - - - - - - - - - - - - - - - -     
        a(4,0,1) =-2.4734758d+09
        a(4,0,2) =-1.5115262d+04

        a(4,1,1) = 3.2270784d+09
        a(4,1,2) = 1.4945586d+04

        a(4,2,1) =-1.7207915d+09
        a(4,2,2) =-5.8361389d+03

        a(4,3,1) = 4.7635740d+08
        a(4,3,2) = 1.1476935d+03

        a(4,4,1) =-7.1414759d+07
        a(4,4,2) =-1.2098290d+02

        a(4,5,1) = 5.4305859d+06
        a(4,5,2) = 6.5439613d+00

        a(4,6,1) =-1.6479928d+05
        a(4,6,2) =-1.4296421d-01
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        s(1,1) = 1
        s(1,2) = 1
        s(1,3) = 0
        s(1,4) = 0

        s(2,1) = 0
        s(2,2) =-3
        s(2,3) = 0
        s(2,4) = 0

        s(3,1) =-1
        s(3,2) = 2
        s(3,3) = 0
        s(3,4) = 0

        s(4,1) =-13
        s(4,2) =-6
        s(4,3) = 0
        s(4,4) = 0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
c       =================================================
        real*8 function u_wk(x)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       A G Fainshtein, N L Manakov and A A Nekipelov,
c       J.Phys.E: At.Mol.Opt.Phys. v.23 (1990) 559-569. 
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /z/z/cl/cl
        real*8 break (5)
        parameter (pi=3.1415926535 8979323846d0)
        real*8 f,p,q,a
        integer s
        common /wk_coeff/f(15,5),p(0:35),q(0:15),
     &  a(4,0:7,4),s(4,4)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        alz=z/cl
        dlam=1.d0-dsqrt(1.d0-alz)*dsqrt(1.d0+alz)
c       - - - - - - - - - - - - - - - - - - - - - - - - -      
        if (x.lt.10.d0) then
          i=1
          if (x.gt.4.d0) then
            i=4
          else
            if (x.gt.1.5d0) then
              i=3
            else
              if (x.gt.0.15d0) i=2
            endif
          endif
c
          v3=0.d0
          if (x.le.4.25d0) then   
            do n=0,8
              v3=v3+p(n)*x**(n-1)
            enddo
            do n=0,5
              v3=v3+p(n+9)*x**(n+2)*dlog(x)
            enddo
            do n=0,4
              v3=v3+p(n+15)*x**(n+3)*dlog(x)**2
            enddo
          else
            do n=0,15
              v3=v3+p(n+20)*(1.d1/x)**(2*n)
            enddo
            v3=v3/x**5
          endif
          v3=v3/cl*alz**3
c
          u_wk=v3
c
          if (x.lt.0.15d0) then
            q5=0.d0
            do n=1,15                              
              q5=q5+q(n)*dlam**n
            enddo
            q5=q5*alz**3/cl
            u_wk=u_wk+q5/x
          endif
c
          sum=1.d0
          do n=1,4
          do m=0,7
            sum=sum+a(i,m,n)*(dlam**n)*x**(m+s(i,n))
          enddo
          enddo
          u_wk=u_wk/sum
c
        else
c
          u_wk=2.d0/225.d0/x**5+59.d0/1323.d0/x**7+(1977.d0+
     &    20.d0*alz*alz)/4725.d0/x**9+(1258680.d0+34144.d0*alz*alz)/
     &    190575/x**11+(1960420032.d0+93618070.d0*alz*alz+
     &    96096.d0*alz**4)/12297285.d0/x**13
          u_wk=u_wk*alz**3    
c      
          prod=576.d0/x**13
          do m=5,15
            prod=prod*m*m/(x*x)
            sum=0.d0
            do n=1,5
              sum=sum+f(m,n)*alz**(2*n+1)
            enddo
            u_wk=u_wk+sum*prod
          enddo
          u_wk=u_wk/(pi*cl)
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -      
        return
        end
c       =================================================
        subroutine dirac(n,kappa,maxii,p,q,a,b,c,y,cp,cq,
     &  r,v,unuc,w)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       This subroutine computes radial Dirac wavefunctions
c       n ......... Principal quantum number
c       kappa ..... Relativistic quantum number
c       cl    ..... Speed of light (in atomic units)
c       e     ..... one-electron energy
c       ii    ..... Number of grid points
c       h     ..... grid step
c       i=1,ii
c       p(i)  ..... Large component of the wavefunction
c       q(i)  ..... Small component of the wavefunction
c       r(i)  ..... radial distance
c       y(i)  ..... Coulomb (or external) potential multiplied
c                   by r(i)
c       unuc(i) ... Nonpoint part of the Nuclear potential
c                   multiplied by r(i)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /cl/cl/z/z
        common /ii/ii/h/h/r1/r1/al/al/bt/bt
        common /nmax/nmax
        common /knucl/knucl/inucl/inucl/rnucl/rnucl
     &  /vnuc/vnuc(20)
        common /niter/niter/nit/nit/mi/m1,m2,m3
        real*8 p(maxii),q(maxii),a(maxii),b(maxii),c(maxii)
        real*8 cp(maxii),cq(maxii),y(maxii),unuc(maxii)
        real*8 r(maxii),v(maxii),w(maxii)
        real*8 fac(30)
        character*1 let(11)
        data
     1  let /'s','p','d','f','g','h','i','k','l','m','n'/
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        nmax=9
        l=(iabs(2*kappa+1)-1)/2
        k=kappa
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do i=1,ii
          y(i)=-z
        enddo
        do m=0,nmax
          i=ii+5+m
          y(i)=0.d0
        enddo
        y(ii+5)=0.d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        eps0=1.d-10*(256.d0/maxii)**2
        if (eps0.lt.1.d-14) eps0=1.d-14
        tol=1.d-50
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        fac(1)=0.d0
        do i=2,30
          d=i-1.d0
          fac(i)=fac(i-1)+dlog(d)
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do i=1,maxii
          p(i)=0.d0
          q(i)=0.d0
        enddo
        e=(z/n)**2
        p(ii+1)=e
        p(ii+2)=1.d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       m1 - First turning point
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        dm=z*z-e*(l+0.5d0)**2
        if (dm.gt.0.1) then
          rm1=(l+0.5d0)**2/(z+dsqrt(dm))
          m1=(al*(rm1-r(1))+bt*dlog(rm1/r(1)))/h+1.d0
        else
          m1=11
        endif
        if (m1.lt.3) m1=3
        m1=(m1/2)*2+1
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        p(ii+15)=m1
        p(ii+3)=ii
        q(ii+3)=ii
        gam=dsqrt(1.d0-(vnuc(1)/cl)**2)
        p(ii+4)=gam
        q(ii+4)=gam
        p(ii+5)=0.d0
        q(ii+5)=0.d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        dz=z
        if (z.lt.0.1d0) dz=1.d0
        a0=dexp(dlog(2.d0*dz/n)+(l+0.5d0)*
     1  dlog(2.d0*dz/n)-
     1  fac(l+l+2)+0.5d0*(fac(n+l+1)-fac(n-l)))/dsqrt(2.d0*n)
        if (kappa.gt.0)
     1  q(ii+5)=-a0*(l+kappa+1)/(2*cl)
        if (kappa.lt.0) p(ii+5)=a0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        c1=0.01d0
        eps=eps0*0.1d0
        hh=0.5d0*h
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        ee1=0.d0
        dd1=0.d0
        niter=0
800     niter=niter+1
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (niter.gt.200) then
          write( *,'(/2x,a)') 'Niter > 200 in dirac'
          write(11,'(/2x,a)') 'Niter > 200 in dirac'
          call exit1
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do i=1,ii
          cp(i)=0.d0
          cq(i)=0.d0
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Boundary conditions at the origin
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        call origin_loc(kappa,p,q,y,unuc,r)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        gam=p(ii+4)
        e=p(ii+1)
        e0=e
        m1=p(ii+15)+c1
        imax=ii
        iemax=0
        emax=0.d0
        iemin=0
        emin=0.d0
        iflag=0
        ifail=0
        istep=0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Calculation of the homogeneous part
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        tn=dsqrt(1.d0/p(ii+2))
        t=hh/cl
        t1=t*tn
        s=0.5d0*e
        hk=hh*kappa
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do 20 i=1,ii
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        c(i)=hk
        yi=y(i)+unuc(i)
        a(i)=t*(s*r(i)+yi)
        d=v(i)/r(i)
        c(i)=c(i)*d
        a(i)=a(i)*d
c       - - - - - - - - - - - - - - - - - - - - - - - - -
20      continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Calculation of the point m3
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        m3=ii
        s=1.d0/(h*h)
        do 30 i=1,ii
        j=ii+1-i
        if (e*v(j)*v(j).ge.s) goto 30
        m3=j
        goto 320
30      continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
320     if (niter.gt.1) then
          if (knucl.ne.1.or.inucl.le.1) then
            call correct(1,ii-2,tn,p,q,cp,cq)
          else
            call correct(1,inucl,tn,p,q,cp,cq)
            call correct(inucl,ii-2,tn,p,q,cp,cq)
          endif
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        n1=0
        n2=0
        n3=0
        nt1=0
        nit=0
        nt3=0
        nj=0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       The begining of the iteration procedure
c       - - - - - - - - - - - - - - - - - - - - - - - - -
600     if (nit.le.380) goto 500
        if(nt1.eq.0.and.nt3.eq.0)
     1  write( *,5) iflag,m1,m2,m3,nj,e0,e
        if (nt1.eq.0.and.nt3.eq.0)
     1  write(11,5) iflag,m1,m2,m3,nj,e0,e
5       format(' Flag=',i1,2x,'m1=',i4,2x,'m2=',i4,2x,
     1  'm3=',i4,2x,'nj =',i2,2x,'e =',2e15.8)
        if (nt1.ne.0.or.nt3.ne.0)
     1  write( *,15) iflag,m1,m2,m3,nj,e0,dp
        if (nt1.ne.0.or.nt3.ne.0)
     1  write(11,15) iflag,m1,m2,m3,nj,e0,dp
15      format(' Flag=',i1,2x,'m1=',i4,2x,'m2=',i4,2x,
     1  'm3=',i4,2x,'nj=',i2,2x,'e =',e15.8,2x,'df=',e12.4)
        if (nit.le.450) goto 500
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do i=1,m2
          write(11,'(6e14.6)') r(i),p(i),q(i),y(i),cp(i),cq(i)
        enddo
        j=2*iabs(kappa)-1
        write( *,25) n,let(l+1),j
        write(11,25) n,let(l+1),j
25      format(/2x,'Divergancy for the shell:',i3,a1,i2,'/2')
        write( *,35) niter
        write(11,35) niter
35      format(2x,'nit>250',2x,'Niter=',i5)
        close (unit=11)
        call exit1
c       - - - - - - - - - - - - - - - - - - - - - - - - -
500     nit=nit+1
        nj=0
        d=0.5d0*(e-e0)*hh/cl
        s=cl*h
        st=0.d0
        do i=1,ii
          t=a(i)+d*v(i)
          a(i)=t
          b(i)=s*v(i)-a(i)
          w(i)=c(i)*c(i)+a(i)*b(i)
          if (st.gt.t) st=t
        enddo
        if (st.lt.0.d0) goto 350
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       iflag=1
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        iflag=1
        n1=n1+1
        if (iemax.eq.0) emax=e
        if (iemax.eq.1.and.e.lt.emax) emax=e
        iemax=1
        e0=e
        e=e*(1.d0-0.5d0*n1/(n1+1.d0))
        if (iemin.eq.1.and.e.le.emin) e=0.5d0*(e0+emin)
        goto 600
350     n1=0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       m2 - Second turning point
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (nt1.gt.0.and.nt3.gt.0) goto 380
        im=m1+m3
        s=0.5d0*kappa*(kappa+1)*hh/cl
        do 70 i=m1,m3
        j=im-i
        rj=r(j)
        df=a(j)+s*v(j)/(rj*rj)
        if (df.gt.0.d0.and.w(j).gt.0.d0) goto 70
        m2=j
        goto 370
70      continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       iflag=2
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        iflag=2
        if (iemax.eq.0) emax=e
        if (iemax.eq.1.and.e.lt.emax) emax=e
        iemax=1
        e0=e
        e=e*0.8d0
        if (iemin.eq.1.and.e.le.emin) e=0.5d0*(e0+emin)
        goto 600
370     if (m2.lt.m3-2) goto 380
        if (iemin.eq.0) emin=e
        if (iemin.eq.1.and.e.gt.emin) emin=e
        iemin=1
        e0=e
        e=e*1.2d0
        if (iemax.eq.1.and.e.ge.emax) e=0.5d0*(e0+emax)
        goto 600
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Boundary conditions at the origin
c       - - - - - - - - - - - - - - - - - - - - - - - - -
380     i0=1
        pi=0.d0
        qi=0.d0
        do m=0,nmax
          j=ii+5+m
          pi=pi+p(j)
          qi=qi+q(j)
        enddo
        st=r(i0)**gam
        pi=pi*st
        qi=qi*st
        p(i0)=(1.d0+c(i0))*pi+b(i0)*qi
        q(i0)=(1.d0-c(i0))*qi+a(i0)*pi
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Calculation functions p(i),q(i) at points (i=1,m2)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        i1=i0+1
        pj=p(i0)
        qj=q(i0)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do 120 i=i1,m2
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        j=i-1
        dw=w(j)
        dc=c(j)
        da=a(j)
        db=b(j)
        dt=0.5d0*(1.d0-dw)
        dq=((dw+dc)*qj-da*pj)/dt+cq(j)
        dp=((dw-dc)*pj-db*qj)/dt+cp(j)
        pj=pj+dp
        qj=qj+dq
        p(i)=pj
        q(i)=qj
c       - - - - - - - - - - - - - - - - - - - - - - - - -
120     continue
        pm2=pj
        qm2=qj
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Testing number of nodes.
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        nj=n-l-1
        do i=2,m2
          if (p(i)*p(i-1).lt.0) nj=nj-1
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (nj.eq.0) goto 390
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       iflag=3
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        iflag=3
        i=1
        if (nj.lt.0) i=-1
        n2=n2+1
        e0=e
        j=iabs(nj)
        if (j.lt.n2) j=n2
        e=e*(1.d0-0.04d0*(i*j)/(j+2.d0))
        if (e.gt.2*cl**2) e=0.5d0*(e0+2*cl**2)
c        if (e.gt.2*cl**2) e=2*cl**2
        if (iemin.eq.1.and.e.le.emin) e=0.5d0*(emin+e0)
        if (iemax.eq.1.and.e.ge.emax) e=0.5d0*(emax+e0)
        if (nj.gt.0.and.iemax.eq.0) emax=e0
        if (nj.gt.0.and.iemax.eq.1.and.e.lt.emax) emax=e0
        if (nj.gt.0) iemax=1
        if (nj.lt.0.and.iemin.eq.0) emin=e0
        if (nj.lt.0.and.iemin.eq.1.and.e.gt.emin) emin=e0
        if (nj.lt.0) iemin=1
        if (iemin.eq.0.or.iemax.eq.0) goto 600
        if (emax-emin.gt.eps*emax*0.00001) goto 600
        ifail=1
        goto 390
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Boundary conditions at the infinity
c       - - - - - - - - - - - - - - - - - - - - - - - - -
390     im=ii-2
        i1=im+imax
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do 130 i=im,imax
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        j=i1-i
        dw=w(j)
        if (dw.gt.0.d0) then
          dw=dsqrt(dw)
        else
          write( *,45)
          write(11,45)
45        format(2x,'Warning!!  R2 too small or',
     1    ' Bound state does not exist.')
          call exit1
        endif
        db=b(j)
        dc=c(j)
        pj=(dw-dc)/db
        qj=cq(j)/(dw+cp(j))
        p(j)=pj
        q(j)=qj
c       - - - - - - - - - - - - - - - - - - - - - - - - -
130     continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Calculation the ratio q(i)/p(i) at points (i=m2,ii)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        qj=qj*(1.d0-dc-pj*db)
        p(im)=pj
        q(im)=qj
        j=im-1
400     dc=c(j)
        da=a(j)
        db=b(j)
        dw=w(j)
        dq=dc+dc
        dt=1.d0+dw+dq+2*db*pj
        if (j.le.m3) dp=(da+da+(1.d0+dw-dq)*pj)/dt
        if (j.gt.m3) dp=(dsqrt(dw)-dc)/db
        qj=(qj+pj*cp(j)-cq(j))/dt*(1.d0-dw)
        pj=dp
        if (j.eq.m2) goto 410
        p(j)=pj
        q(j)=qj
        j=j-1
        goto 400
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Sewing the wave functions at the point m2
c       - - - - - - - - - - - - - - - - - - - - - - - - -
410     dq=qm2
        qm2=pj*pm2+qj
        dp=qm2-dq
        dq=dq+kappa/r(m2)*pm2/(2*cl)
        t=dabs(dp)/(dabs(dq)+eps/cl)
c       if (t.lt.5*eps.and.nit.gt.2) goto 450
        if (t.lt.5*eps) goto 450
        if (dp*b(m2).le.0.d0) goto 420
        j=1
        e3=e
        dt3=dp
        nt3=1
        goto 430
420     e1=e
        dt1=dp
        nt1=1
        j=-1
430     if (nt1.eq.nt3) goto 440
        istep=0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       iflag=4
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        iflag=4
        i=n-l-1
        n3=n3+1
        d=i-2*(i/2)-0.5d0
        t=0.01d0
        e0=e
        d=t*(j*n3)*d/(n3+3.d0)
        e=e*(1.d0+d)
c        if (e.gt.2*cl**2) e=0.5d0*(e0+2*cl**2)
         goto 600
        if (iemin.eq.1.and.e.le.emin) e=0.5d0*(e0+emin)
        if (iemax.eq.1.and.e.ge.emax) e=0.5d0*(e0+emax)
        if (d.le.0.d0.and.iemax.eq.0) emax=e0
        if (d.le.0.d0.and.iemax.eq.1.and.e.le.emax) emax=e0
        if (d.le.0.d0) iemax=1
        if (d.ge.0.d0.and.iemin.eq.0) emin=e0
        if (d.ge.0.d0.and.iemin.eq.1.and.e.ge.emin) emin=e0
        if (d.ge.0.d0) iemin=1
        if (iemin.eq.0.or.iemax.eq.0) goto 600
        if (emax-emin.gt.eps*emax*0.00001) goto 600
        ifail=1
        goto 450
440     n3=0
        e0=e
        e=(dt3*e1-dt1*e3)/(dt3-dt1)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        istep=istep+1
        if (istep.gt.3) then
           e=0.5d0*(e1+e3)
           istep=0
        endif
        if (nt3.gt.5.or.nt1.gt.5) then
          e=0.5d0*(e1+e3)
          nt3=1
          nt1=1
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       iflag=5
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        iflag=5
        goto 600
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Calculation functions p(i),q(i) at points (i=m2,ii)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
450     i2=m2
        im=ii-2
        do 140 i=m2+1,im
        j=i-1
        d=c(j)+c(j)
        s=b(j)+b(j)
        t=((1.d0+w(j)-d)*p(j)-s*q(j))/(1.d0-w(j))
        pj=t+cp(j)
        qj=p(i)*pj+q(i)
        if (dabs(pj).lt.1.d-20) goto 460
        i2=i
        p(i)=pj
        q(i)=qj
140     continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Reconstruction functions p(i),q(i) at points (i=1,m2)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
460     do j=1,i2
          pj=p(j)
          qj=q(j)
          t=1.d0-w(j)
          p(j)=((1.d0-c(j))*pj-b(j)*qj)/t
          q(j)=((1.d0+c(j))*qj-a(j)*pj)/t
        enddo
        if (i2.lt.im) goto 470
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        i2=im
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do 160 i=im+1,ii
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        d=(r(i)-r(i-1))/(hh*v(i))
        s=dsqrt(w(i))
        t=cp(i)-s
        t1=t*d
        if (dabs(t1).lt.0.01d0) t2=-d
        if (dabs(t1).ge.0.01d0) t2=(1.d0-dexp(t1))/t
        t1=p(i-1)*dexp(-s*d)
        t2=t2*b(i)*q(i)
        pj=t1+t2
        if (dabs(pj).lt.tol) goto 470
        i2=i
        t2=t1
        q(i)=pj*p(i)+q(i)
        p(i)=pj
c       - - - - - - - - - - - - - - - - - - - - - - - - -
160     continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        goto 1000
c       - - - - - - - - - - - - - - - - - - - - - - - - -
470     imax=i2
        if (imax.eq.ii) goto 1000
        do i=imax+1,ii
          p(i)=0.d0
          q(i)=0.d0
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
1000    p(ii+1)=e
        p(ii+3)=imax
        p(ii+18)=m2
        q(ii+1)=e
        q(ii+3)=imax
        if (ifail.ne.0) nit=-nit
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (dabs(p(imax)).gt.1.d-3) then
        write( *,65) kappa,p(imax)
        write(11,65) kappa,p(imax)
65      format(2x,'Warning!!  R2 too small. Kappa =',i2,
     1  '  p(imax) =',e9.2)
        if (ifail.ne.0) then
        write( *,75)
        write(11,75)
75      format(2x,'Warning!!  Failed to find eigenvalue.',
     1  '  Kappa =',i2)
        endif
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Normalization
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        d=tint(0,p,q,p,q,r,v)
        d=1.d0/dsqrt(d)
        do i=1,imax
          p(i)=p(i)*d
          q(i)=q(i)*d
        enddo
        do m=0,nmax
          i=m+ii+5
          p(i)=p(i)*d
          q(i)=q(i)*d
        enddo
        p(ii+2)=1.d0
        q(ii+2)=1.d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        ee0=0.5d0*p(ii+1)
        dd0=(p(ii+5)**2+q(ii+5)**2)/p(ii+2)
        de=dabs((ee0-ee1)/(ee0+ee1+eps))
        dn=dabs((dd0-dd1)/(dd0+dd1+eps))
        del=de
        if (dn.gt.de) del=dn
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        ee1=ee0
        dd1=dd0
        if (niter.ge.3) then
          if (del.gt.1.d-10.or.niter.le.2) goto 800
        else
          goto 800
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
c       =================================================
        subroutine correct(imin,imax,tn,p,q,cp,cq)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /h/h
        real*8 p(*),q(*),cp(*),cq(*)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        s=tn/3.d0
        t=tn/8.d0
        d=tn/120.d0
        i0=imin+2
        im=imax-3
        do 10 i=i0,im
        j =i+1
        j1=j+1
        j2=j1+1
        i1=i-1
        i2=i1-1
        tp=s*(p(j)-p(i))-t*(p(j1)-p(i1))+d*(p(j2)-p(i2))
        tq=s*(q(j)-q(i))-t*(q(j1)-q(i1))+d*(q(j2)-q(i2))
c        tp=tn/12.d0*(-p(j1)+3.d0*p(j)-3.d0*p(i)+p(i1))
c        tq=tn/12.d0*(-q(j1)+3.d0*q(j)-3.d0*q(i)+q(i1))
        cp(i)=cp(i)+tp
        cq(i)=cq(i)+tq
10      continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        i=imin+1
        j =i+1
        j1=j+1
        j2=j1+1
        j3=j2+1
        i1=i-1
        tp=d*(9.d0*p(i1)-25.d0*p(i)+20.d0*p(j)-
     1  5.d0*p(j2)+p(j3))
        tq=d*(9.d0*q(i1)-25.d0*q(i)+20.d0*q(j)-
     1  5.d0*q(j2)+q(j3))
        cp(i)=cp(i)+tp
        cq(i)=cq(i)+tq
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        i=imin
        j =i+1
        j1=j+1
        j2=j1+1
        j3=j2+1
        j4=j3+1
        tp=d*(29.d0*p(i)-115.d0*p(j)+180.d0*p(j1)-
     1  140.d0*p(j2)+55.d0*p(j3)-9.d0*p(j4))
        tq=d*(29.d0*q(i)-115.d0*q(j)+180.d0*q(j1)-
     1  140.d0*q(j2)+55.d0*q(j3)-9.d0*q(j4))
        cp(i)=cp(i)+tp
        cq(i)=cq(i)+tq
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        i=imax-2
        i1=i+1
        i2=i1+1
        j =i-1
        j1=j-1
        j2=j1-1
        tp=d*(9.d0*p(i2)-25.d0*p(i1)+20.d0*p(i)-
     1  5.d0*p(j1)+p(j2))
        tq=d*(9.d0*q(i2)-25.d0*q(i1)+20.d0*q(i)-
     1  5.d0*q(j1)+q(j2))
        cp(i)=cp(i)-tp
        cq(i)=cq(i)-tq
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        i=imax-1
        i1=i+1
        j =i-1
        j1=j-1
        j2=j1-1
        j3=j2-1
        tp=d*(29.d0*p(i1)-115.d0*p(i)+180.d0*p(j)-
     1  140.d0*p(j1)+55.d0*p(j2)-9.d0*p(j3))
        tq=d*(29.d0*q(i1)-115.d0*q(i)+180.d0*q(j)-
     1  140.d0*q(j1)+55.d0*q(j2)-9.d0*q(j3))
        cp(i)=cp(i)-tp
        cq(i)=cq(i)-tq
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
c       =================================================
        subroutine origin_loc(kappa,p,q,y,unuc,r)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /ii/ii/cl/cl/z/z
        common /nmax/nmax
        real*8 p(*),q(*),y(*),unuc(*),r(*)
        dimension yy(20),p1(20),q1(20),vp(20),vq(20)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Taylor series expansion of the wavefunctions around r=0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        c1=0.01d0
        e=p(ii+1)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        yy(1)=0.d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        ny=y(ii+4)+c1
        do n=0,nmax
          i=n+1
          yy(i+ny)= y(ii+n+5)/r(1)**n
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        vp(1)=unuc(ii+5)-z
        vq(1)=unuc(ii+5)-z
        do m=1,nmax
          i=ii+5+m
          ui=unuc(i)/r(1)**m
          vp(m+1)=yy(m+1)+ui
          vq(m+1)=yy(m+1)+ui
        enddo
c
        vp(2)=vp(2)+(0.5d0*e-2*cl*cl)
        vq(2)=vq(2)+0.5d0*e
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        gam=dsqrt(kappa**2-(vp(1)/cl)**2)
        t=iabs(kappa)
        s=gam+t
        d=dsqrt((p(ii+5)**2+q(ii+5)**2)*0.5d0*s/
     1  (t*p(ii+2)))
        t=d*dabs(vp(1))/(s*cl)
        if (kappa.gt.0) goto 350
        p1(1)=d
        q1(1)=t
        goto 360
350     p1(1)=t
        q1(1)=-d
360     do n=1,nmax
          i=n+1
          d=n*(n+2*gam)*cl*cl
          t1=(n+gam-kappa)*cl
          t2=(n+gam+kappa)*cl
          s1=0.d0
          s2=0.d0
          do m=1,n
            j=m+1
            s1=s1+p1(i-m)*vq(j)
            s2=s2+q1(i-m)*vp(j)
          enddo
          p1(i)=(-vp(1)*s1+t1*s2)/d
          q1(i)=(-vq(1)*s2-t2*s1)/d
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        p(ii+4)=gam
        q(ii+4)=gam
        do n=0,nmax
          i=n+1
          p(ii+5+n)=p1(i)*r(1)**n
          q(ii+5+n)=q1(i)*r(1)**n
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
1000    return
        end
c       =================================================
        subroutine atom_data
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /z/z/r2/r2/cl/cl
        common /name/name
        common /knucl/knucl/rnucl/rnucl
        character*1 let(11)
        common /Rrms/Rrms(120)/name_at/name_at(120)
        character*2 name,name_at
        data
     1  let /'s','p','d','f','g','h','i','k','l','m','n'/
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        call init_atom_data
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        iz=z+0.5d0
        name=name_at(iz)
        write( *,5) name
        write(11,5) name
5       format (/2x,a)
        write( *,15) z
        write(11,15) z
15      format (2x,'Z=',f10.4)
        rms=Rrms(iz)
        write( *,25) rms
        write(11,25) rms
25      format (2x,'Rms=',f8.4,' fm')
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Defult value of the lihgt velocity.
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c        cl=137.03599911d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Default value of Atomic box.
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       r2=150.0/z
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        write( *,'(2x,a,f16.9)') 'Cl=',cl
        write(11,'(2x,a,f16.9)') 'Cl=',cl
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Nuclear model:
c       Knucl=1  Uniform sphere model
c       Knucl=2  Fermi model
c       Knucl=3  Gauss model
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (iz.le.15) knucl=1
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (dabs(rms).lt.1.d-10) then
          knucl=0
          rnucl=0.d0
          write( *,'(/2x,a)') 'Point nuuclei.'
          write(11,'(/2x,a)') 'Point nuuclei.'
        else
          fermi2at=1.d-5/0.529177249d0
          rfermi=rms*dsqrt(5.d0/3.d0)
          rnucl=rfermi*fermi2at
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (knucl.lt.0.or.knucl.gt.3) then
          write( *,'(/2x,a,3x,a,i4)') 'Wrong Nuclear model !!!',
     &    'Nucl=',knucl
          write(11,'(/2x,a,3x,a,i4)') 'Wrong Nuclear model !!!',
     &    'Nucl=',knucl
         call exit1
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
 1000   return
        call exit1
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        end
c       =================================================
        subroutine init_atom_data
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /Rrms/Rrms(120)/name_at/name_at(120)
        common /Am/Am(120)
        character*2 name_at
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        name_at( 1)=' H'
        name_at( 2)='He'
        name_at( 3)='Li'
        name_at( 4)='Be'
        name_at( 5)=' B'
        name_at( 6)=' C'
        name_at( 7)=' N'
        name_at( 8)=' O'
        name_at( 9)=' F'
        name_at(10)='Ne'
        name_at(11)='Na'
        name_at(12)='Mg'
        name_at(13)='Al'
        name_at(14)='Al'
        name_at(15)=' P'
        name_at(16)=' S'
        name_at(17)='Cl'
        name_at(18)='Ar'
        name_at(19)='K '
        name_at(20)='Ca'
        name_at(21)='Sc'
        name_at(22)='Ti'
        name_at(23)=' V'
        name_at(24)='Cr'
        name_at(25)='Mn'
        name_at(26)='Fe'
        name_at(27)='Co'
        name_at(28)='Ni'
        name_at(29)='Cu'
        name_at(30)='Zn'
        name_at(31)='Ga'
        name_at(32)='Ge'
        name_at(33)='As'
        name_at(34)='Se'
        name_at(35)='Br'
        name_at(36)='Kr'
        name_at(37)='Rb'
        name_at(38)='Sr'
        name_at(39)=' Y'
        name_at(40)='Zr'
        name_at(41)='Nb'
        name_at(42)='Mo'
        name_at(43)='Tc'
        name_at(44)='Ru'
        name_at(45)='Rh'
        name_at(46)='Pd'
        name_at(47)='Ag'
        name_at(48)='Cd'
        name_at(49)='In'
        name_at(50)='Sn'
        name_at(51)='Sb'
        name_at(52)='Te'
        name_at(53)=' I'
        name_at(54)='Xe'
        name_at(55)='Cs'
        name_at(56)='Ba'
        name_at(57)='La'
        name_at(58)='Ce'
        name_at(59)='Pr'
        name_at(60)='Nd'
        name_at(61)='Pm'
        name_at(62)='Sm'
        name_at(63)='Eu'
        name_at(64)='Gd'
        name_at(65)='Tb'
        name_at(66)='Dy'
        name_at(67)='Ho'
        name_at(68)='Er'
        name_at(69)='Tm'
        name_at(70)='Yb'
        name_at(71)='Lu'
        name_at(72)='Hf'
        name_at(73)='Ta'
        name_at(74)=' W'
        name_at(75)='Re'
        name_at(76)='Os'
        name_at(77)='Ir'
        name_at(78)='Pt'
        name_at(79)='Au'
        name_at(80)='Hg'
        name_at(81)='Tl'
        name_at(82)='Pb'
        name_at(83)='Bi'
        name_at(84)='Po'
        name_at(85)='At'
        name_at(86)='Rn'
        name_at(87)='Fr'
        name_at(88)='Ra'
        name_at(89)='Ac'
        name_at(90)='Th'
        name_at(91)='Pa'
        name_at(92)=' U'
        name_at(93)='Np'
        name_at(94)='Pu'
        name_at(95)='Am'
        name_at(96)='Cm'
        name_at(97)='Bk'
        name_at(98)='Cf'
        name_at(99)='Es'
        name_at(100)='Fm'
        name_at(101)='Md'
        name_at(102)='No'
        name_at(103)='Lr'
c       Rutherfordium
        name_at(104)='Rf'
        name_at(105)='Db'
        name_at(106)='Sg'
        name_at(107)='Bh'
        name_at(108)='Hs'
        name_at(109)='Mt'
        name_at(110)='Ds'
        name_at(111)='Rg'
        name_at(112)='Cn' 
        name_at(113)='Ut'
        name_at(114)='Fl'
        name_at(115)='Up'
        name_at(116)='Lv'
        name_at(117)='Us'
        name_at(118)='Uo'
        name_at(119)='Ue'
        name_at(120)='Ub'
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       From NIST: http://www.nist.gov/pml/data/comp.cfm
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        Am(  1) =  1.0
        Am(  2) =  4.0
        Am(  3) =  6.0
        Am(  4) =  9.0
        Am(  5) =  10.0
        Am(  6) =  12.0
        Am(  7) =  14.0
        Am(  8) =  16.0
        Am(  9) =  19.0
        Am( 10) =  20.0
        Am( 11) =  23.0
        Am( 12) =  24.0
        Am( 13) =  27.0
        Am( 14) =  28.0
        Am( 15) =  31.0
        Am( 16) =  32.0
        Am( 17) =  35.0
        Am( 18) =  38.0
        Am( 19) =  39.0
        Am( 20) =  40.0
        Am( 21) =  45.0
        Am( 22) =  48.0
        Am( 23) =  51.0
        Am( 24) =  52.0
        Am( 25) =  55.0
        Am( 26) =  56.0
        Am( 27) =  59.0
        Am( 28) =  60.0
        Am( 29) =  63.0
        Am( 30) =  66.0
        Am( 31) =  69.0
        Am( 32) =  74.0
        Am( 33) =  75.0
        Am( 34) =  80.0
        Am( 35) =  79.0
        Am( 36) =  86.0
        Am( 37) =  87.0
        Am( 38) =  88.0
        Am( 39) =  89.0
        Am( 40) =  90.0
        Am( 41) =  93.0
        Am( 42) =  92.0
        Am( 43) =  97.0
        Am( 44) = 104.0
        Am( 45) = 103.0
        Am( 46) = 108.0
        Am( 47) = 109.0
        Am( 48) = 114.0
        Am( 49) = 115.0
        Am( 50) = 120.0
        Am( 51) = 121.0
        Am( 52) = 130.0
        Am( 53) = 127.0
        Am( 54) = 136.0
        Am( 55) = 133.0
        Am( 56) = 138.0
        Am( 57) = 139.0
        Am( 58) = 140.0
        Am( 59) = 141.0
        Am( 60) = 142.0
        Am( 61) = 145.0
        Am( 62) = 144.0
        Am( 63) = 145.0
        Am( 64) = 160.0
        Am( 65) = 159.0
        Am( 66) = 148.0
        Am( 67) = 165.0
        Am( 68) = 170.0
        Am( 69) = 169.0
        Am( 70) = 176.0
        Am( 71) = 175.0
        Am( 72) = 178.0
        Am( 73) = 181.0
        Am( 74) = 184.0
        Am( 75) = 185.0
        Am( 76) = 192.0
        Am( 77) = 191.0
        Am( 78) = 194.0
        Am( 79) = 197.0
        Am( 80) = 198.0
        Am( 81) = 205.0
        Am( 82) = 208.0
        Am( 83) = 209.0
        Am( 84) = 208.0
        Am( 85) = 210.0
        Am( 86) = 212.0
        Am( 87) = 212.0
        Am( 88) = 214.0
        Am( 89) = 227.0
        Am( 90) = 232.0
        Am( 91) = 232.0
        Am( 92) = 238.0
        Am( 93) = 237.0
        Am( 94) = 239.0
        Am( 95) = 243.0
        Am( 96) = 244.0
        Am( 97) = 247.0
        Am( 98) = 251.0
        Am( 99) = 252.0
        Am(100) = 257.0
        Am(101) = 258.0
        Am(102) = 259.0
        Am(103) = 262.0
        Am(104) = 265.0
        Am(105) = 268.0
        Am(106) = 271.0
        Am(107) = 272.0
        Am(108) = 270.0
        Am(109) = 276.0
        Am(110) = 281.0
        Am(111) = 280.0
        Am(112) = 285.0
        Am(113) = 284.0
        Am(114) = 289.0
        Am(115) = 288.0
        Am(116) = 293.0
        Am(117) = 292.0
        Am(118) = 294.0
        Am(119) = 295.0
        Am(120) = 295.0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do iz=1,120
          if (iz.le.90) then
            rms=(0.836d0*Am(iz)**(1.d0/3.d0)+0.570d0)
          else
            rms=(0.77d0*Am(iz)**(1.d0/3.d0)+0.980d0)
          endif
          Rrms(iz)=rms
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Nuclear R_rms in fm.
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        Rrms  (  1) =         0.8783 d0 !#       A=1
        Rrms  (  2) =         1.6755 d0 !#       A=4        
        Rrms  (  3) =         2.5890 d0 !#       A=6
        Rrms  (  4) =         2.5190 d0 !#       A=9
        Rrms  (  5) =         2.4277 d0 !#       A=10
        Rrms  (  6) =         2.4702 d0 !#       A=12
        Rrms  (  7) =         2.5582 d0 !#       A=14
        Rrms  (  8) =         2.6991 d0 !#       A=16
        Rrms  (  9) =         2.8976 d0 !#       A=19
        Rrms  ( 10) =         3.0055 d0 !#       A=20
        Rrms  ( 11) =         2.9936 d0 !#       A=23
        Rrms  ( 12) =         3.0570 d0 !#       A=24
        Rrms  ( 13) =         3.0610 d0 !#       A=27
        Rrms  ( 14) =         3.1224 d0 !#       A=28
        Rrms  ( 15) =         3.1889 d0 !#       A=31
        Rrms  ( 16) =         3.3468 d0 !#       A=32
        Rrms  ( 17) =         3.3654 d0 !#       A=35
        Rrms  ( 18) =         3.4028 d0 !#       A=38
        Rrms  ( 19) =         3.4349 d0 !#       A=39
        Rrms  ( 20) =         3.4776 d0 !#       A=40
        Rrms  ( 21) =         3.5459 d0 !#       A=45
        Rrms  ( 22) =         3.5921 d0 !#       A=48
        Rrms  ( 23) =         3.6002 d0 !#       A=51
        Rrms  ( 24) =         3.6452 d0 !#       A=52
        Rrms  ( 25) =         3.7057 d0 !#       A=55
        Rrms  ( 26) =         3.7377 d0 !#       A=56
        Rrms  ( 27) =         3.7875 d0 !#       A=59
        Rrms  ( 28) =         3.8118 d0 !#       A=60
        Rrms  ( 29) =         3.8823 d0 !#       A=63
        Rrms  ( 30) =         3.9491 d0 !#       A=66
        Rrms  ( 31) =         3.9973 d0 !#       A=69
        Rrms  ( 32) =         4.0742 d0 !#       A=74
        Rrms  ( 33) =         4.0968 d0 !#       A=75
        Rrms  ( 34) =         4.1400 d0 !#       A=80
        Rrms  ( 35) =         4.1629 d0 !#       A=79
        Rrms  ( 36) =         4.1835 d0 !#       A=86
        Rrms  ( 37) =         4.1989 d0 !#       A=87
        Rrms  ( 38) =         4.2240 d0 !#       A=88
        Rrms  ( 39) =         4.2430 d0 !#       A=89
        Rrms  ( 40) =         4.2694 d0 !#       A=90
        Rrms  ( 41) =         4.3240 d0 !#       A=93
        Rrms  ( 42) =         4.3151 d0 !#       A=92
c       Rrms  ( 43) =         0.0000 d0 !#       A=97
        Rrms  ( 44) =         4.5098 d0 !#       A=104
        Rrms  ( 45) =         4.4945 d0 !#       A=103
        Rrms  ( 46) =         4.5563 d0 !#       A=108
        Rrms  ( 47) =         4.5638 d0 !#       A=109
        Rrms  ( 48) =         4.6087 d0 !#       A=114
        Rrms  ( 49) =         4.6156 d0 !#       A=115
        Rrms  ( 50) =         4.6519 d0 !#       A=120
        Rrms  ( 51) =         4.6802 d0 !#       A=121
        Rrms  ( 52) =         4.7423 d0 !#       A=130
        Rrms  ( 53) =         4.7964 d0 !#       A=127
        Rrms  ( 54) =         4.5098 d0 !#       A=136
        Rrms  ( 55) =         4.8041 d0 !#       A=133
        Rrms  ( 56) =         4.8378 d0 !#       A=138
        Rrms  ( 57) =         4.8550 d0 !#       A=139
        Rrms  ( 58) =         4.8771 d0 !#       A=140
        Rrms  ( 59) =         4.8919 d0 !#       A=141
        Rrms  ( 60) =         4.9123 d0 !#       A=142
c       Rrms  ( 61) =         0.0000 d0 !#       A=145
        Rrms  ( 62) =         4.9524 d0 !#       A=144
        Rrms  ( 63) =         4.9663 d0 !#       A=145
        Rrms  ( 64) =         5.1734 d0 !#       A=160
        Rrms  ( 65) =         5.0600 d0 !#       A=159
        Rrms  ( 66) =         5.0455 d0 !#       A=148
        Rrms  ( 67) =         5.2022 d0 !#       A=165
        Rrms  ( 68) =         5.2789 d0 !#       A=170
        Rrms  ( 69) =         5.2256 d0 !#       A=169
        Rrms  ( 70) =         5.3215 d0 !#       A=176
        Rrms  ( 71) =         5.3770 d0 !#       A=175
        Rrms  ( 72) =         5.3371 d0 !#       A=178
        Rrms  ( 73) =         5.3507 d0 !#       A=181
        Rrms  ( 74) =         5.3658 d0 !#       A=184
        Rrms  ( 75) =         5.3596 d0 !#       A=185
        Rrms  ( 76) =         5.4126 d0 !#       A=192
        Rrms  ( 77) =         5.3968 d0 !#       A=191
        Rrms  ( 78) =         5.4236 d0 !#       A=194
        Rrms  ( 79) =         5.4371 d0 !#       A=197
        Rrms  ( 80) =         5.4463 d0 !#       A=198
        Rrms  ( 81) =         5.4759 d0 !#       A=205
        Rrms  ( 82) =         5.5012 d0 !#       A=208
        Rrms  ( 83) =         5.5211 d0 !#       A=209
        Rrms  ( 84) =         5.5584 d0 !#       A=208
c       Rrms  ( 85) =         0.0000 d0 !#       A=210
        Rrms  ( 86) =         5.5915 d0 !#       A=212
        Rrms  ( 87) =         5.5915 d0 !#       A=212
        Rrms  ( 88) =         5.6079 d0 !#       A=214
c       Rrms  ( 89) =         0.0000 d0 !#       A=227
        Rrms  ( 90) =         5.7848 d0 !#       A=232
c       Rrms  ( 91) =         0.0000 d0 !#
        Rrms  ( 92) =         5.8571 d0 !#       A=238
c       Rrms  ( 93) =         0.00000d0 !#
        Rrms  ( 94) =         5.8601 d0 !#       A=239
        Rrms  ( 95) =         5.9048 d0 !#       A=243
        Rrms  ( 96) =         5.8420 d0 !#       A=244
        Rrms  (100) =         5.8570 d0 !#  fm 
        Rrms  (105) =         5.9190 d0 !#  fm 
        Rrms  (110) =         5.9930 d0 !#  fm 
        Rrms  (115) =         6.0880 d0 !#  fm 
        Rrms  (120) =         6.1750 d0 !#  fm 
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
c       =================================================
        subroutine FSE_dat(ka,n1,n2,iz,SE_pnt,SE_ext)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       This subroutine computes diagonal and non-diagonal
c       self-energy (SE) matrix elements, using interpolation
c       procedure and tabulated data taken from the paper:
c       V.M. Shabaev, I.I. Tupitsyn and V. A. Yerokhin,
c       PRA, 88, 012513 (2013).
c
c       ka    .... Relativistic quantum number
c       n1,n2 .... principal quntum numbers
c       SE_pnt ... SE matrix element value, calculated
c                  with point nuclear model
c       SE_ext ... SE matrix element value, calculated
c                  with extended nuclear model
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        real*8 Func(120)
        parameter (TINY = 1.d-16)
        parameter (acl = 137.0359895d0)
c
        real*8 FSE(120,-3:2,5,5),errtot(120,-3:2,5,5)
        common /re_stor/FSE,errtot

        real*8 FSEfns(120,-3:2,5,5),errfns(120,-3:2,5,5)
        common /re_fns_stor/FSEfns,errfns
c
        Func = 0.d0
        do i = 1,120
          Func(i) = FSE(i,KA,N1,N2)
          if (KA.eq.-1.and.abs(Func(i)).gt.TINY) then
            Func(i) = Func(i) - 4.d0/3.d0*(-2)*log(i/acl)
          end if
        enddo
        SE_pnt = DInterpolation(iz,Func)
        if (KA.eq.-1) SE_pnt = SE_pnt + 4.d0/3.d0*(-2)*log(iz/acl)
c
        Func = 0.d0
        do i = 1,120
          Func(i) = FSE(i,KA,N1,N2) + FSEfns(i,KA,N1,N2)
          if (KA.eq.-1.and.abs(Func(i)).gt.TINY) then
            Func(i) = Func(i) - 4.d0/3.d0*(-2)*log(i/acl)
          end if
       enddo
      SE_ext = DInterpolation(iz,Func)
      if (KA.eq.-1) SE_ext = SE_ext + 4.d0/3.d0*(-2)*log(iz/acl)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
      return
      end
c===========================================================================
c
c===========================================================================
      real*8 function DInterpolation(iz,Func)
c
c 4- or 3-point interpolation of an array Func(Z) with "holes"
c
      implicit real*8 (a-h,o-z)
      real*8 Func(120),F(0:10)
      integer indx(120)
      parameter (TINY = 1.d-16)
c---
      if (dabs(Func(iz)).gt.TINY) then
        DInterpolation = Func(iz)
        return
      end if
c
c-
c     indexing
c--
      l = 1
      do i = 1,120
        if (dabs(Func(i)).gt.TINY) then
          indx(l) = i
          l = l+1
        end if
      enddo
      lmax = l-1
      if (lmax.eq.0) stop 'DInterpolation: nothing to interpolate'

      if (iz.lt.indx(1)) stop 'DInterpolation: iz too small'
      if (iz.gt.indx(lmax)) stop 'DInterpolation: iz too large'
c
c-
c     determining points of interpolation
c--
      do l = 1,lmax
        if (indx(l).gt.iz) then
          imid = l
          go to 111
        end if
      enddo
111   continue

      if (imid.le.2) then
        imin = 1
        imax = 3
      else if (imid.ge.lmax) then
        imin = lmax-2
        imax = lmax
      else
        imin = imid-2
        imax = imid+1
      end if

      do j = imin,imax
        F(j-imin) = Func(indx(j))
      enddo
c
c-
c     Interpolation
c---
      do j = imin+1,imax
         a = F(j-imin-1)
         b = indx(j-1)
         do i = j,imax
           F(i-imin) = (a- F(i-imin))/ (b- indx(i))
         enddo
      enddo

      FFF = F(imax-imin)
      do i = imax-1,imin,-1
         delt = iz-indx(i)
         FFF = F(i-imin)+ delt* FFF
      enddo

      DInterpolation = FFF
c--
      return
      end
c=======================================================================
c
c=======================================================================
      subroutine se_fn_stor()
      implicit real*8 (a-h,o-z)

      real*8 FSEfns(120,-3:2,5,5),errtot(120,-3:2,5,5)
      common /re_fns_stor/FSEfns,errtot
*
* ns
*
      FSEfns( 10, -1, 1, 1) =        -0.0000332422 d0 !# Fermi    3.005 fm
      errtot( 10, -1, 1, 1) =         0.0000000005 d0 !# 
      FSEfns( 10, -1, 1, 2) =        -0.0000329499 d0 !# Fermi    3.005 fm
      errtot( 10, -1, 1, 2) =         0.0000000007 d0 !# 
      FSEfns( 10, -1, 1, 3) =        -0.0000327626 d0 !# Fermi    3.005 fm
      errtot( 10, -1, 1, 3) =         0.0000000113 d0 !# 
      FSEfns( 10, -1, 1, 4) =        -0.0000326765 d0 !# Fermi    3.005 fm
      errtot( 10, -1, 1, 4) =         0.0000000719 d0 !# 
      FSEfns( 10, -1, 1, 5) =        -0.0000325876 d0 !# Fermi    3.005 fm
      errtot( 10, -1, 1, 5) =         0.0000000222 d0 !# 
      FSEfns( 10, -1, 2, 2) =        -0.0000325388 d0 !# Fermi    3.005 fm
      errtot( 10, -1, 2, 2) =         0.0000000010 d0 !# 
      FSEfns( 10, -1, 2, 3) =        -0.0000323342 d0 !# Fermi    3.005 fm
      errtot( 10, -1, 2, 3) =         0.0000000047 d0 !# 
      FSEfns( 10, -1, 2, 4) =        -0.0000322250 d0 !# Fermi    3.005 fm
      errtot( 10, -1, 2, 4) =         0.0000000115 d0 !# 
      FSEfns( 10, -1, 2, 5) =        -0.0000321482 d0 !# Fermi    3.005 fm
      errtot( 10, -1, 2, 5) =         0.0000000068 d0 !# 
      FSEfns( 10, -1, 3, 3) =        -0.0000321143 d0 !# Fermi    3.005 fm
      errtot( 10, -1, 3, 3) =         0.0000000379 d0 !# 
      FSEfns( 10, -1, 3, 4) =        -0.0000319885 d0 !# Fermi    3.005 fm
      errtot( 10, -1, 3, 4) =         0.0000000149 d0 !# 
      FSEfns( 10, -1, 3, 5) =        -0.0000319218 d0 !# Fermi    3.005 fm
      errtot( 10, -1, 3, 5) =         0.0000000277 d0 !# 
      FSEfns( 10, -1, 4, 4) =        -0.0000318745 d0 !# Fermi    3.005 fm
      errtot( 10, -1, 4, 4) =         0.0000000315 d0 !# 
      FSEfns( 10, -1, 4, 5) =        -0.0000318093 d0 !# Fermi    3.005 fm
      errtot( 10, -1, 4, 5) =         0.0000003396 d0 !# 
      FSEfns( 10, -1, 5, 5) =        -0.0000320003 d0 !# Fermi    3.005 fm
      errtot( 10, -1, 5, 5) =         0.0000023940 d0 !# 
      FSEfns( 15, -1, 1, 1) =        -0.0000601284 d0 !# Fermi    3.189 fm
      errtot( 15, -1, 1, 1) =         0.0000000002 d0 !# 
      FSEfns( 15, -1, 1, 2) =        -0.0000597668 d0 !# Fermi    3.189 fm
      errtot( 15, -1, 1, 2) =         0.0000000028 d0 !# 
      FSEfns( 15, -1, 1, 3) =        -0.0000593837 d0 !# Fermi    3.189 fm
      errtot( 15, -1, 1, 3) =         0.0000000059 d0 !# 
      FSEfns( 15, -1, 1, 4) =        -0.0000591536 d0 !# Fermi    3.189 fm
      errtot( 15, -1, 1, 4) =         0.0000000310 d0 !# 
      FSEfns( 15, -1, 1, 5) =        -0.0000589855 d0 !# Fermi    3.189 fm
      errtot( 15, -1, 1, 5) =         0.0000000104 d0 !# 
      FSEfns( 15, -1, 2, 2) =        -0.0000591026 d0 !# Fermi    3.189 fm
      errtot( 15, -1, 2, 2) =         0.0000000015 d0 !# 
      FSEfns( 15, -1, 2, 3) =        -0.0000586694 d0 !# Fermi    3.189 fm
      errtot( 15, -1, 2, 3) =         0.0000000083 d0 !# 
      FSEfns( 15, -1, 2, 4) =        -0.0000584297 d0 !# Fermi    3.189 fm
      errtot( 15, -1, 2, 4) =         0.0000000012 d0 !# 
      FSEfns( 15, -1, 2, 5) =        -0.0000582710 d0 !# Fermi    3.189 fm
      errtot( 15, -1, 2, 5) =         0.0000000132 d0 !# 
      FSEfns( 15, -1, 3, 3) =        -0.0000582003 d0 !# Fermi    3.189 fm
      errtot( 15, -1, 3, 3) =         0.0000000204 d0 !# 
      FSEfns( 15, -1, 3, 4) =        -0.0000579280 d0 !# Fermi    3.189 fm
      errtot( 15, -1, 3, 4) =         0.0000000099 d0 !# 
      FSEfns( 15, -1, 3, 5) =        -0.0000577684 d0 !# Fermi    3.189 fm
      errtot( 15, -1, 3, 5) =         0.0000000420 d0 !# 
      FSEfns( 15, -1, 4, 4) =        -0.0000576602 d0 !# Fermi    3.189 fm
      errtot( 15, -1, 4, 4) =         0.0000000246 d0 !# 
      FSEfns( 15, -1, 4, 5) =        -0.0000574819 d0 !# Fermi    3.189 fm
      errtot( 15, -1, 4, 5) =         0.0000000875 d0 !# 
      FSEfns( 15, -1, 5, 5) =        -0.0000571881 d0 !# Fermi    3.189 fm
      errtot( 15, -1, 5, 5) =         0.0000003894 d0 !# 
      FSEfns( 20, -1, 1, 1) =        -0.0001029398 d0 !# Fermi    3.476 fm
      errtot( 20, -1, 1, 1) =         0.0000000007 d0 !# 
      FSEfns( 20, -1, 1, 2) =        -0.0001028204 d0 !# Fermi    3.476 fm
      errtot( 20, -1, 1, 2) =         0.0000000060 d0 !# 
      FSEfns( 20, -1, 1, 3) =        -0.0001021307 d0 !# Fermi    3.476 fm
      errtot( 20, -1, 1, 3) =         0.0000000068 d0 !# 
      FSEfns( 20, -1, 1, 4) =        -0.0001016689 d0 !# Fermi    3.476 fm
      errtot( 20, -1, 1, 4) =         0.0000000131 d0 !# 
      FSEfns( 20, -1, 1, 5) =        -0.0001013497 d0 !# Fermi    3.476 fm
      errtot( 20, -1, 1, 5) =         0.0000000204 d0 !# 
      FSEfns( 20, -1, 2, 2) =        -0.0001020626 d0 !# Fermi    3.476 fm
      errtot( 20, -1, 2, 2) =         0.0000000040 d0 !# 
      FSEfns( 20, -1, 2, 3) =        -0.0001012607 d0 !# Fermi    3.476 fm
      errtot( 20, -1, 2, 3) =         0.0000000097 d0 !# 
      FSEfns( 20, -1, 2, 4) =        -0.0001007934 d0 !# Fermi    3.476 fm
      errtot( 20, -1, 2, 4) =         0.0000000031 d0 !# 
      FSEfns( 20, -1, 2, 5) =        -0.0001004774 d0 !# Fermi    3.476 fm
      errtot( 20, -1, 2, 5) =         0.0000000135 d0 !# 
      FSEfns( 20, -1, 3, 3) =        -0.0001003779 d0 !# Fermi    3.476 fm
      errtot( 20, -1, 3, 3) =         0.0000000193 d0 !# 
      FSEfns( 20, -1, 3, 4) =        -0.0000998355 d0 !# Fermi    3.476 fm
      errtot( 20, -1, 3, 4) =         0.0000000084 d0 !# 
      FSEfns( 20, -1, 3, 5) =        -0.0000995160 d0 !# Fermi    3.476 fm
      errtot( 20, -1, 3, 5) =         0.0000000316 d0 !# 
      FSEfns( 20, -1, 4, 4) =        -0.0000992965 d0 !# Fermi    3.476 fm
      errtot( 20, -1, 4, 4) =         0.0000000395 d0 !# 
      FSEfns( 20, -1, 4, 5) =        -0.0000989307 d0 !# Fermi    3.476 fm
      errtot( 20, -1, 4, 5) =         0.0000000352 d0 !# 
      FSEfns( 20, -1, 5, 5) =        -0.0000985979 d0 !# Fermi    3.476 fm
      errtot( 20, -1, 5, 5) =         0.0000005790 d0 !# 
      FSEfns( 25, -1, 1, 1) =        -0.0001599385 d0 !# Fermi    3.706 fm
      errtot( 25, -1, 1, 1) =         0.0000000005 d0 !# 
      FSEfns( 25, -1, 1, 2) =        -0.0001608116 d0 !# Fermi    3.706 fm
      errtot( 25, -1, 1, 2) =         0.0000000055 d0 !# 
      FSEfns( 25, -1, 1, 3) =        -0.0001597261 d0 !# Fermi    3.706 fm
      errtot( 25, -1, 1, 3) =         0.0000000049 d0 !# 
      FSEfns( 25, -1, 1, 4) =        -0.0001589237 d0 !# Fermi    3.706 fm
      errtot( 25, -1, 1, 4) =         0.0000000083 d0 !# 
      FSEfns( 25, -1, 1, 5) =        -0.0001583700 d0 !# Fermi    3.706 fm
      errtot( 25, -1, 1, 5) =         0.0000000191 d0 !# 
      FSEfns( 25, -1, 2, 2) =        -0.0001605573 d0 !# Fermi    3.706 fm
      errtot( 25, -1, 2, 2) =         0.0000000069 d0 !# 
      FSEfns( 25, -1, 2, 3) =        -0.0001592569 d0 !# Fermi    3.706 fm
      errtot( 25, -1, 2, 3) =         0.0000000041 d0 !# 
      FSEfns( 25, -1, 2, 4) =        -0.0001584533 d0 !# Fermi    3.706 fm
      errtot( 25, -1, 2, 4) =         0.0000000059 d0 !# 
      FSEfns( 25, -1, 2, 5) =        -0.0001578948 d0 !# Fermi    3.706 fm
      errtot( 25, -1, 2, 5) =         0.0000000212 d0 !# 
      FSEfns( 25, -1, 3, 3) =        -0.0001578029 d0 !# Fermi    3.706 fm
      errtot( 25, -1, 3, 3) =         0.0000000188 d0 !# 
      FSEfns( 25, -1, 3, 4) =        -0.0001568532 d0 !# Fermi    3.706 fm
      errtot( 25, -1, 3, 4) =         0.0000000196 d0 !# 
      FSEfns( 25, -1, 3, 5) =        -0.0001562851 d0 !# Fermi    3.706 fm
      errtot( 25, -1, 3, 5) =         0.0000000178 d0 !# 
      FSEfns( 25, -1, 4, 4) =        -0.0001559109 d0 !# Fermi    3.706 fm
      errtot( 25, -1, 4, 4) =         0.0000000683 d0 !# 
      FSEfns( 25, -1, 4, 5) =        -0.0001552462 d0 !# Fermi    3.706 fm
      errtot( 25, -1, 4, 5) =         0.0000000402 d0 !# 
      FSEfns( 25, -1, 5, 5) =        -0.0001547911 d0 !# Fermi    3.706 fm
      errtot( 25, -1, 5, 5) =         0.0000006040 d0 !# 
      FSEfns( 30, -1, 1, 1) =        -0.0002389360 d0 !# Fermi    3.929 fm
      errtot( 30, -1, 1, 1) =         0.0000000006 d0 !# 
      FSEfns( 30, -1, 1, 2) =        -0.0002422039 d0 !# Fermi    3.929 fm
      errtot( 30, -1, 1, 2) =         0.0000000057 d0 !# 
      FSEfns( 30, -1, 1, 3) =        -0.0002406004 d0 !# Fermi    3.929 fm
      errtot( 30, -1, 1, 3) =         0.0000000033 d0 !# 
      FSEfns( 30, -1, 1, 4) =        -0.0002392910 d0 !# Fermi    3.929 fm
      errtot( 30, -1, 1, 4) =         0.0000000180 d0 !# 
      FSEfns( 30, -1, 1, 5) =        -0.0002383593 d0 !# Fermi    3.929 fm
      errtot( 30, -1, 1, 5) =         0.0000000015 d0 !# 
      FSEfns( 30, -1, 2, 2) =        -0.0002436753 d0 !# Fermi    3.929 fm
      errtot( 30, -1, 2, 2) =         0.0000000061 d0 !# 
      FSEfns( 30, -1, 2, 3) =        -0.0002416886 d0 !# Fermi    3.929 fm
      errtot( 30, -1, 2, 3) =         0.0000000020 d0 !# 
      FSEfns( 30, -1, 2, 4) =        -0.0002403755 d0 !# Fermi    3.929 fm
      errtot( 30, -1, 2, 4) =         0.0000000039 d0 !# 
      FSEfns( 30, -1, 2, 5) =        -0.0002394333 d0 !# Fermi    3.929 fm
      errtot( 30, -1, 2, 5) =         0.0000000298 d0 !# 
      FSEfns( 30, -1, 3, 3) =        -0.0002394353 d0 !# Fermi    3.929 fm
      errtot( 30, -1, 3, 3) =         0.0000000122 d0 !# 
      FSEfns( 30, -1, 3, 4) =        -0.0002378587 d0 !# Fermi    3.929 fm
      errtot( 30, -1, 3, 4) =         0.0000000265 d0 !# 
      FSEfns( 30, -1, 3, 5) =        -0.0002368978 d0 !# Fermi    3.929 fm
      errtot( 30, -1, 3, 5) =         0.0000000131 d0 !# 
      FSEfns( 30, -1, 4, 4) =        -0.0002362964 d0 !# Fermi    3.929 fm
      errtot( 30, -1, 4, 4) =         0.0000000936 d0 !# 
      FSEfns( 30, -1, 4, 5) =        -0.0002351690 d0 !# Fermi    3.929 fm
      errtot( 30, -1, 4, 5) =         0.0000000892 d0 !# 
      FSEfns( 30, -1, 5, 5) =        -0.0002342814 d0 !# Fermi    3.929 fm
      errtot( 30, -1, 5, 5) =         0.0000008817 d0 !# 
      FSEfns( 35, -1, 1, 1) =        -0.0003509700 d0 !# Fermi    4.163 fm
      errtot( 35, -1, 1, 1) =         0.0000000012 d0 !# 
      FSEfns( 35, -1, 1, 2) =        -0.0003591894 d0 !# Fermi    4.163 fm
      errtot( 35, -1, 1, 2) =         0.0000000025 d0 !# 
      FSEfns( 35, -1, 1, 3) =        -0.0003568998 d0 !# Fermi    4.163 fm
      errtot( 35, -1, 1, 3) =         0.0000000030 d0 !# 
      FSEfns( 35, -1, 1, 4) =        -0.0003548049 d0 !# Fermi    4.163 fm
      errtot( 35, -1, 1, 4) =         0.0000000247 d0 !# 
      FSEfns( 35, -1, 1, 5) =        -0.0003532721 d0 !# Fermi    4.163 fm
      errtot( 35, -1, 1, 5) =         0.0000000490 d0 !# 
      FSEfns( 35, -1, 2, 2) =        -0.0003647582 d0 !# Fermi    4.163 fm
      errtot( 35, -1, 2, 2) =         0.0000000046 d0 !# 
      FSEfns( 35, -1, 2, 3) =        -0.0003618143 d0 !# Fermi    4.163 fm
      errtot( 35, -1, 2, 3) =         0.0000000163 d0 !# 
      FSEfns( 35, -1, 2, 4) =        -0.0003597109 d0 !# Fermi    4.163 fm
      errtot( 35, -1, 2, 4) =         0.0000000217 d0 !# 
      FSEfns( 35, -1, 2, 5) =        -0.0003581513 d0 !# Fermi    4.163 fm
      errtot( 35, -1, 2, 5) =         0.0000000009 d0 !# 
      FSEfns( 35, -1, 3, 3) =        -0.0003584288 d0 !# Fermi    4.163 fm
      errtot( 35, -1, 3, 3) =         0.0000000069 d0 !# 
      FSEfns( 35, -1, 3, 4) =        -0.0003558694 d0 !# Fermi    4.163 fm
      errtot( 35, -1, 3, 4) =         0.0000000387 d0 !# 
      FSEfns( 35, -1, 3, 5) =        -0.0003542773 d0 !# Fermi    4.163 fm
      errtot( 35, -1, 3, 5) =         0.0000000639 d0 !# 
      FSEfns( 35, -1, 4, 4) =        -0.0003533391 d0 !# Fermi    4.163 fm
      errtot( 35, -1, 4, 4) =         0.0000001260 d0 !# 
      FSEfns( 35, -1, 4, 5) =        -0.0003514290 d0 !# Fermi    4.163 fm
      errtot( 35, -1, 4, 5) =         0.0000003211 d0 !# 
      FSEfns( 35, -1, 5, 5) =        -0.0003500550 d0 !# Fermi    4.163 fm
      errtot( 35, -1, 5, 5) =         0.0000014354 d0 !# 
      FSEfns( 40, -1, 1, 1) =        -0.0004810963 d0 !# Fermi    4.270 fm
      errtot( 40, -1, 1, 1) =         0.0000000016 d0 !# 
      FSEfns( 40, -1, 1, 2) =        -0.0004977758 d0 !# Fermi    4.270 fm
      errtot( 40, -1, 1, 2) =         0.0000000013 d0 !# 
      FSEfns( 40, -1, 1, 3) =        -0.0004947724 d0 !# Fermi    4.270 fm
      errtot( 40, -1, 1, 3) =         0.0000000099 d0 !# 
      FSEfns( 40, -1, 1, 4) =        -0.0004916555 d0 !# Fermi    4.270 fm
      errtot( 40, -1, 1, 4) =         0.0000000179 d0 !# 
      FSEfns( 40, -1, 1, 5) =        -0.0004893123 d0 !# Fermi    4.270 fm
      errtot( 40, -1, 1, 5) =         0.0000000582 d0 !# 
      FSEfns( 40, -1, 2, 2) =        -0.0005110657 d0 !# Fermi    4.270 fm
      errtot( 40, -1, 2, 2) =         0.0000000043 d0 !# 
      FSEfns( 40, -1, 2, 3) =        -0.0005070380 d0 !# Fermi    4.270 fm
      errtot( 40, -1, 2, 3) =         0.0000000303 d0 !# 
      FSEfns( 40, -1, 2, 4) =        -0.0005038920 d0 !# Fermi    4.270 fm
      errtot( 40, -1, 2, 4) =         0.0000000199 d0 !# 
      FSEfns( 40, -1, 2, 5) =        -0.0005014962 d0 !# Fermi    4.270 fm
      errtot( 40, -1, 2, 5) =         0.0000000139 d0 !# 
      FSEfns( 40, -1, 3, 3) =        -0.0005023451 d0 !# Fermi    4.270 fm
      errtot( 40, -1, 3, 3) =         0.0000000150 d0 !# 
      FSEfns( 40, -1, 3, 4) =        -0.0004984775 d0 !# Fermi    4.270 fm
      errtot( 40, -1, 3, 4) =         0.0000000414 d0 !# 
      FSEfns( 40, -1, 3, 5) =        -0.0004960325 d0 !# Fermi    4.270 fm
      errtot( 40, -1, 3, 5) =         0.0000000778 d0 !# 
      FSEfns( 40, -1, 4, 4) =        -0.0004946894 d0 !# Fermi    4.270 fm
      errtot( 40, -1, 4, 4) =         0.0000003456 d0 !# 
      FSEfns( 40, -1, 4, 5) =        -0.0004916343 d0 !# Fermi    4.270 fm
      errtot( 40, -1, 4, 5) =         0.0000006611 d0 !# 
      FSEfns( 40, -1, 5, 5) =        -0.0004897984 d0 !# Fermi    4.270 fm
      errtot( 40, -1, 5, 5) =         0.0000017338 d0 !# 
      FSEfns( 45, -1, 1, 1) =        -0.0006899295 d0 !# Fermi    4.494 fm
      errtot( 45, -1, 1, 1) =         0.0000000027 d0 !# 
      FSEfns( 45, -1, 1, 2) =        -0.0007226410 d0 !# Fermi    4.494 fm
      errtot( 45, -1, 1, 2) =         0.0000000020 d0 !# 
      FSEfns( 45, -1, 1, 3) =        -0.0007185546 d0 !# Fermi    4.494 fm
      errtot( 45, -1, 1, 3) =         0.0000000104 d0 !# 
      FSEfns( 45, -1, 1, 4) =        -0.0007136817 d0 !# Fermi    4.494 fm
      errtot( 45, -1, 1, 4) =         0.0000000073 d0 !# 
      FSEfns( 45, -1, 1, 5) =        -0.0007099282 d0 !# Fermi    4.494 fm
      errtot( 45, -1, 1, 5) =         0.0000000132 d0 !# 
      FSEfns( 45, -1, 2, 2) =        -0.0007512565 d0 !# Fermi    4.494 fm
      errtot( 45, -1, 2, 2) =         0.0000000032 d0 !# 
      FSEfns( 45, -1, 2, 3) =        -0.0007455098 d0 !# Fermi    4.494 fm
      errtot( 45, -1, 2, 3) =         0.0000000429 d0 !# 
      FSEfns( 45, -1, 2, 4) =        -0.0007405573 d0 !# Fermi    4.494 fm
      errtot( 45, -1, 2, 4) =         0.0000000341 d0 !# 
      FSEfns( 45, -1, 2, 5) =        -0.0007366872 d0 !# Fermi    4.494 fm
      errtot( 45, -1, 2, 5) =         0.0000000101 d0 !# 
      FSEfns( 45, -1, 3, 3) =        -0.0007387361 d0 !# Fermi    4.494 fm
      errtot( 45, -1, 3, 3) =         0.0000000307 d0 !# 
      FSEfns( 45, -1, 3, 4) =        -0.0007325972 d0 !# Fermi    4.494 fm
      errtot( 45, -1, 3, 4) =         0.0000000907 d0 !# 
      FSEfns( 45, -1, 3, 5) =        -0.0007286547 d0 !# Fermi    4.494 fm
      errtot( 45, -1, 3, 5) =         0.0000000422 d0 !# 
      FSEfns( 45, -1, 4, 4) =        -0.0007266961 d0 !# Fermi    4.494 fm
      errtot( 45, -1, 4, 4) =         0.0000009636 d0 !# 
      FSEfns( 45, -1, 4, 5) =        -0.0007215272 d0 !# Fermi    4.494 fm
      errtot( 45, -1, 4, 5) =         0.0000009620 d0 !# 
      FSEfns( 45, -1, 5, 5) =        -0.0007190190 d0 !# Fermi    4.494 fm
      errtot( 45, -1, 5, 5) =         0.0000017977 d0 !# 
      FSEfns( 50, -1, 1, 1) =        -0.0009606070 d0 !# Fermi    4.654 fm
      errtot( 50, -1, 1, 1) =         0.0000000016 d0 !# 
      FSEfns( 50, -1, 1, 2) =        -0.0010198793 d0 !# Fermi    4.654 fm
      errtot( 50, -1, 1, 2) =         0.0000000054 d0 !# 
      FSEfns( 50, -1, 1, 3) =        -0.0010145370 d0 !# Fermi    4.654 fm
      errtot( 50, -1, 1, 3) =         0.0000000220 d0 !# 
      FSEfns( 50, -1, 1, 4) =        -0.0010071015 d0 !# Fermi    4.654 fm
      errtot( 50, -1, 1, 4) =         0.0000000444 d0 !# 
      FSEfns( 50, -1, 1, 5) =        -0.0010012478 d0 !# Fermi    4.654 fm
      errtot( 50, -1, 1, 5) =         0.0000000359 d0 !# 
      FSEfns( 50, -1, 2, 2) =        -0.0010752132 d0 !# Fermi    4.654 fm
      errtot( 50, -1, 2, 2) =         0.0000000148 d0 !# 
      FSEfns( 50, -1, 2, 3) =        -0.0010672745 d0 !# Fermi    4.654 fm
      errtot( 50, -1, 2, 3) =         0.0000000240 d0 !# 
      FSEfns( 50, -1, 2, 4) =        -0.0010596473 d0 !# Fermi    4.654 fm
      errtot( 50, -1, 2, 4) =         0.0000000433 d0 !# 
      FSEfns( 50, -1, 2, 5) =        -0.0010535496 d0 !# Fermi    4.654 fm
      errtot( 50, -1, 2, 5) =         0.0000000342 d0 !# 
      FSEfns( 50, -1, 3, 3) =        -0.0010578180 d0 !# Fermi    4.654 fm
      errtot( 50, -1, 3, 3) =         0.0000000500 d0 !# 
      FSEfns( 50, -1, 3, 4) =        -0.0010483094 d0 !# Fermi    4.654 fm
      errtot( 50, -1, 3, 4) =         0.0000002631 d0 !# 
      FSEfns( 50, -1, 3, 5) =        -0.0010421234 d0 !# Fermi    4.654 fm
      errtot( 50, -1, 3, 5) =         0.0000000779 d0 !# 
      FSEfns( 50, -1, 4, 4) =        -0.0010394333 d0 !# Fermi    4.654 fm
      errtot( 50, -1, 4, 4) =         0.0000020847 d0 !# 
      FSEfns( 50, -1, 4, 5) =        -0.0010309845 d0 !# Fermi    4.654 fm
      errtot( 50, -1, 4, 5) =         0.0000007844 d0 !# 
      FSEfns( 50, -1, 5, 5) =        -0.0010270474 d0 !# Fermi    4.654 fm
      errtot( 50, -1, 5, 5) =         0.0000015529 d0 !# 
      FSEfns( 55, -1, 1, 1) =        -0.0013335911 d0 !# Fermi    4.804 fm
      errtot( 55, -1, 1, 1) =         0.0000000021 d0 !# 
      FSEfns( 55, -1, 1, 2) =        -0.0014370597 d0 !# Fermi    4.804 fm
      errtot( 55, -1, 1, 2) =         0.0000000259 d0 !# 
      FSEfns( 55, -1, 1, 3) =        -0.0014301303 d0 !# Fermi    4.804 fm
      errtot( 55, -1, 1, 3) =         0.0000000224 d0 !# 
      FSEfns( 55, -1, 1, 4) =        -0.0014187550 d0 !# Fermi    4.804 fm
      errtot( 55, -1, 1, 4) =         0.0000000784 d0 !# 
      FSEfns( 55, -1, 1, 5) =        -0.0014096232 d0 !# Fermi    4.804 fm
      errtot( 55, -1, 1, 5) =         0.0000000613 d0 !# 
      FSEfns( 55, -1, 2, 2) =        -0.0015386464 d0 !# Fermi    4.804 fm
      errtot( 55, -1, 2, 2) =         0.0000000503 d0 !# 
      FSEfns( 55, -1, 2, 3) =        -0.0015276926 d0 !# Fermi    4.804 fm
      errtot( 55, -1, 2, 3) =         0.0000000074 d0 !# 
      FSEfns( 55, -1, 2, 4) =        -0.0015158787 d0 !# Fermi    4.804 fm
      errtot( 55, -1, 2, 4) =         0.0000000522 d0 !# 
      FSEfns( 55, -1, 2, 5) =        -0.0015062469 d0 !# Fermi    4.804 fm
      errtot( 55, -1, 2, 5) =         0.0000000488 d0 !# 
      FSEfns( 55, -1, 3, 3) =        -0.0015145509 d0 !# Fermi    4.804 fm
      errtot( 55, -1, 3, 3) =         0.0000000853 d0 !# 
      FSEfns( 55, -1, 3, 4) =        -0.0014997633 d0 !# Fermi    4.804 fm
      errtot( 55, -1, 3, 4) =         0.0000006116 d0 !# 
      FSEfns( 55, -1, 3, 5) =        -0.0014900295 d0 !# Fermi    4.804 fm
      errtot( 55, -1, 3, 5) =         0.0000001935 d0 !# 
      FSEfns( 55, -1, 4, 4) =        -0.0014864435 d0 !# Fermi    4.804 fm
      errtot( 55, -1, 4, 4) =         0.0000037224 d0 !# 
      FSEfns( 55, -1, 4, 5) =        -0.0014728644 d0 !# Fermi    4.804 fm
      errtot( 55, -1, 4, 5) =         0.0000002960 d0 !# 
      FSEfns( 55, -1, 5, 5) =        -0.0014659606 d0 !# Fermi    4.804 fm
      errtot( 55, -1, 5, 5) =         0.0000021184 d0 !# 
      FSEfns( 60, -1, 1, 1) =        -0.0018307500 d0 !# Fermi    4.912 fm
      errtot( 60, -1, 1, 1) =         0.0000000054 d0 !# 
      FSEfns( 60, -1, 1, 2) =        -0.0020049152 d0 !# Fermi    4.912 fm
      errtot( 60, -1, 1, 2) =         0.0000000497 d0 !# 
      FSEfns( 60, -1, 1, 3) =        -0.0019960492 d0 !# Fermi    4.912 fm
      errtot( 60, -1, 1, 3) =         0.0000000150 d0 !# 
      FSEfns( 60, -1, 1, 4) =        -0.0019787368 d0 !# Fermi    4.912 fm
      errtot( 60, -1, 1, 4) =         0.0000000937 d0 !# 
      FSEfns( 60, -1, 1, 5) =        -0.0019646143 d0 !# Fermi    4.912 fm
      errtot( 60, -1, 1, 5) =         0.0000000692 d0 !# 
      FSEfns( 60, -1, 2, 2) =        -0.0021832739 d0 !# Fermi    4.912 fm
      errtot( 60, -1, 2, 2) =         0.0000001052 d0 !# 
      FSEfns( 60, -1, 2, 3) =        -0.0021682585 d0 !# Fermi    4.912 fm
      errtot( 60, -1, 2, 3) =         0.0000000260 d0 !# 
      FSEfns( 60, -1, 2, 4) =        -0.0021500008 d0 !# Fermi    4.912 fm
      errtot( 60, -1, 2, 4) =         0.0000000625 d0 !# 
      FSEfns( 60, -1, 2, 5) =        -0.0021348830 d0 !# Fermi    4.912 fm
      errtot( 60, -1, 2, 5) =         0.0000000531 d0 !# 
      FSEfns( 60, -1, 3, 3) =        -0.0021501875 d0 !# Fermi    4.912 fm
      errtot( 60, -1, 3, 3) =         0.0000002127 d0 !# 
      FSEfns( 60, -1, 3, 4) =        -0.0021272933 d0 !# Fermi    4.912 fm
      errtot( 60, -1, 3, 4) =         0.0000012268 d0 !# 
      FSEfns( 60, -1, 3, 5) =        -0.0021120615 d0 !# Fermi    4.912 fm
      errtot( 60, -1, 3, 5) =         0.0000006369 d0 !# 
      FSEfns( 60, -1, 4, 4) =        -0.0021073997 d0 !# Fermi    4.912 fm
      errtot( 60, -1, 4, 4) =         0.0000056250 d0 !# 
      FSEfns( 60, -1, 4, 5) =        -0.0020861511 d0 !# Fermi    4.912 fm
      errtot( 60, -1, 4, 5) =         0.0000022123 d0 !# 
      FSEfns( 60, -1, 5, 5) =        -0.0020738903 d0 !# Fermi    4.912 fm
      errtot( 60, -1, 5, 5) =         0.0000025490 d0 !# 
      FSEfns( 65, -1, 1, 1) =        -0.0025637408 d0 !# Fermi    5.060 fm
      errtot( 65, -1, 1, 1) =         0.0000000075 d0 !# 
      FSEfns( 65, -1, 1, 2) =        -0.0028570956 d0 !# Fermi    5.060 fm
      errtot( 65, -1, 1, 2) =         0.0000000639 d0 !# 
      FSEfns( 65, -1, 1, 3) =        -0.0028454345 d0 !# Fermi    5.060 fm
      errtot( 65, -1, 1, 3) =         0.0000000139 d0 !# 
      FSEfns( 65, -1, 1, 4) =        -0.0028183548 d0 !# Fermi    5.060 fm
      errtot( 65, -1, 1, 4) =         0.0000001067 d0 !# 
      FSEfns( 65, -1, 1, 5) =        -0.0027959862 d0 !# Fermi    5.060 fm
      errtot( 65, -1, 1, 5) =         0.0000001159 d0 !# 
      FSEfns( 65, -1, 2, 2) =        -0.0031688653 d0 !# Fermi    5.060 fm
      errtot( 65, -1, 2, 2) =         0.0000001646 d0 !# 
      FSEfns( 65, -1, 2, 3) =        -0.0031476209 d0 !# Fermi    5.060 fm
      errtot( 65, -1, 2, 3) =         0.0000000686 d0 !# 
      FSEfns( 65, -1, 2, 4) =        -0.0031185340 d0 !# Fermi    5.060 fm
      errtot( 65, -1, 2, 4) =         0.0000000766 d0 !# 
      FSEfns( 65, -1, 2, 5) =        -0.0030941812 d0 !# Fermi    5.060 fm
      errtot( 65, -1, 2, 5) =         0.0000000908 d0 !# 
      FSEfns( 65, -1, 3, 3) =        -0.0031221001 d0 !# Fermi    5.060 fm
      errtot( 65, -1, 3, 3) =         0.0000004314 d0 !# 
      FSEfns( 65, -1, 3, 4) =        -0.0030858122 d0 !# Fermi    5.060 fm
      errtot( 65, -1, 3, 4) =         0.0000018549 d0 !# 
      FSEfns( 65, -1, 3, 5) =        -0.0030612573 d0 !# Fermi    5.060 fm
      errtot( 65, -1, 3, 5) =         0.0000009851 d0 !# 
      FSEfns( 65, -1, 4, 4) =        -0.0030543298 d0 !# Fermi    5.060 fm
      errtot( 65, -1, 4, 4) =         0.0000078246 d0 !# 
      FSEfns( 65, -1, 4, 5) =        -0.0030208782 d0 !# Fermi    5.060 fm
      errtot( 65, -1, 4, 5) =         0.0000034049 d0 !# 
      FSEfns( 65, -1, 5, 5) =        -0.0030001571 d0 !# Fermi    5.060 fm
      errtot( 65, -1, 5, 5) =         0.0000032217 d0 !# 
      FSEfns( 70, -1, 1, 1) =        -0.0037301886 d0 !# Fermi    5.311 fm
      errtot( 70, -1, 1, 1) =         0.0000000158 d0 !# 
      FSEfns( 70, -1, 1, 2) =        -0.0042358562 d0 !# Fermi    5.311 fm
      errtot( 70, -1, 1, 2) =         0.0000000830 d0 !# 
      FSEfns( 70, -1, 1, 3) =        -0.0042195309 d0 !# Fermi    5.311 fm
      errtot( 70, -1, 1, 3) =         0.0000000456 d0 !# 
      FSEfns( 70, -1, 1, 4) =        -0.0041751454 d0 !# Fermi    5.311 fm
      errtot( 70, -1, 1, 4) =         0.0000001267 d0 !# 
      FSEfns( 70, -1, 1, 5) =        -0.0041381687 d0 !# Fermi    5.311 fm
      errtot( 70, -1, 1, 5) =         0.0000001898 d0 !# 
      FSEfns( 70, -1, 2, 2) =        -0.0047918360 d0 !# Fermi    5.311 fm
      errtot( 70, -1, 2, 2) =         0.0000002132 d0 !# 
      FSEfns( 70, -1, 2, 3) =        -0.0047599370 d0 !# Fermi    5.311 fm
      errtot( 70, -1, 2, 3) =         0.0000001449 d0 !# 
      FSEfns( 70, -1, 2, 4) =        -0.0047112483 d0 !# Fermi    5.311 fm
      errtot( 70, -1, 2, 4) =         0.0000000567 d0 !# 
      FSEfns( 70, -1, 2, 5) =        -0.0046702143 d0 !# Fermi    5.311 fm
      errtot( 70, -1, 2, 5) =         0.0000001487 d0 !# 
      FSEfns( 70, -1, 3, 3) =        -0.0047219516 d0 !# Fermi    5.311 fm
      errtot( 70, -1, 3, 3) =         0.0000006215 d0 !# 
      FSEfns( 70, -1, 3, 4) =        -0.0046617976 d0 !# Fermi    5.311 fm
      errtot( 70, -1, 3, 4) =         0.0000026881 d0 !# 
      FSEfns( 70, -1, 3, 5) =        -0.0046203283 d0 !# Fermi    5.311 fm
      errtot( 70, -1, 3, 5) =         0.0000014642 d0 !# 
      FSEfns( 70, -1, 4, 4) =        -0.0046091153 d0 !# Fermi    5.311 fm
      errtot( 70, -1, 4, 4) =         0.0000112682 d0 !# 
      FSEfns( 70, -1, 4, 5) =        -0.0045543444 d0 !# Fermi    5.311 fm
      errtot( 70, -1, 4, 5) =         0.0000050217 d0 !# 
      FSEfns( 70, -1, 5, 5) =        -0.0045183050 d0 !# Fermi    5.311 fm
      errtot( 70, -1, 5, 5) =         0.0000042695 d0 !# 
      FSEfns( 75, -1, 1, 1) =        -0.0051106650 d0 !# Fermi    5.339 fm
      errtot( 75, -1, 1, 1) =         0.0000000327 d0 !# 
      FSEfns( 75, -1, 1, 2) =        -0.0059214635 d0 !# Fermi    5.339 fm
      errtot( 75, -1, 1, 2) =         0.0000000734 d0 !# 
      FSEfns( 75, -1, 1, 3) =        -0.0058994600 d0 !# Fermi    5.339 fm
      errtot( 75, -1, 1, 3) =         0.0000000464 d0 !# 
      FSEfns( 75, -1, 1, 4) =        -0.0058306355 d0 !# Fermi    5.339 fm
      errtot( 75, -1, 1, 4) =         0.0000001408 d0 !# 
      FSEfns( 75, -1, 1, 5) =        -0.0057729713 d0 !# Fermi    5.339 fm
      errtot( 75, -1, 1, 5) =         0.0000002576 d0 !# 
      FSEfns( 75, -1, 2, 2) =        -0.0068422119 d0 !# Fermi    5.339 fm
      errtot( 75, -1, 2, 2) =         0.0000002356 d0 !# 
      FSEfns( 75, -1, 2, 3) =        -0.0067964399 d0 !# Fermi    5.339 fm
      errtot( 75, -1, 2, 3) =         0.0000001655 d0 !# 
      FSEfns( 75, -1, 2, 4) =        -0.0067191270 d0 !# Fermi    5.339 fm
      errtot( 75, -1, 2, 4) =         0.0000000425 d0 !# 
      FSEfns( 75, -1, 2, 5) =        -0.0066537731 d0 !# Fermi    5.339 fm
      errtot( 75, -1, 2, 5) =         0.0000002201 d0 !# 
      FSEfns( 75, -1, 3, 3) =        -0.0067428296 d0 !# Fermi    5.339 fm
      errtot( 75, -1, 3, 3) =         0.0000008649 d0 !# 
      FSEfns( 75, -1, 3, 4) =        -0.0066484879 d0 !# Fermi    5.339 fm
      errtot( 75, -1, 3, 4) =         0.0000037224 d0 !# 
      FSEfns( 75, -1, 3, 5) =        -0.0065823318 d0 !# Fermi    5.339 fm
      errtot( 75, -1, 3, 5) =         0.0000021117 d0 !# 
      FSEfns( 75, -1, 4, 4) =        -0.0065653209 d0 !# Fermi    5.339 fm
      errtot( 75, -1, 4, 4) =         0.0000152098 d0 !# 
      FSEfns( 75, -1, 4, 5) =        -0.0064805313 d0 !# Fermi    5.339 fm
      errtot( 75, -1, 4, 5) =         0.0000071704 d0 !# 
      FSEfns( 75, -1, 5, 5) =        -0.0064215581 d0 !# Fermi    5.339 fm
      errtot( 75, -1, 5, 5) =         0.0000050562 d0 !# 
      FSEfns( 80, -1, 1, 1) =        -0.0072903626 d0 !# Fermi    5.463 fm
      errtot( 80, -1, 1, 1) =         0.0000000618 d0 !# 
      FSEfns( 80, -1, 1, 2) =        -0.0086303557 d0 !# Fermi    5.463 fm
      errtot( 80, -1, 1, 2) =         0.0000000557 d0 !# 
      FSEfns( 80, -1, 1, 3) =        -0.0085976786 d0 !# Fermi    5.463 fm
      errtot( 80, -1, 1, 3) =         0.0000000440 d0 !# 
      FSEfns( 80, -1, 1, 4) =        -0.0084855894 d0 !# Fermi    5.463 fm
      errtot( 80, -1, 1, 4) =         0.0000001473 d0 !# 
      FSEfns( 80, -1, 1, 5) =        -0.0083915131 d0 !# Fermi    5.463 fm
      errtot( 80, -1, 1, 5) =         0.0000003934 d0 !# 
      FSEfns( 80, -1, 2, 2) =        -0.0102004289 d0 !# Fermi    5.463 fm
      errtot( 80, -1, 2, 2) =         0.0000002872 d0 !# 
      FSEfns( 80, -1, 2, 3) =        -0.0101297301 d0 !# Fermi    5.463 fm
      errtot( 80, -1, 2, 3) =         0.0000002314 d0 !# 
      FSEfns( 80, -1, 2, 4) =        -0.0100005107 d0 !# Fermi    5.463 fm
      errtot( 80, -1, 2, 4) =         0.0000000262 d0 !# 
      FSEfns( 80, -1, 2, 5) =        -0.0098914184 d0 !# Fermi    5.463 fm
      errtot( 80, -1, 2, 5) =         0.0000003354 d0 !# 
      FSEfns( 80, -1, 3, 3) =        -0.0100490494 d0 !# Fermi    5.463 fm
      errtot( 80, -1, 3, 3) =         0.0000012824 d0 !# 
      FSEfns( 80, -1, 3, 4) =        -0.0098936307 d0 !# Fermi    5.463 fm
      errtot( 80, -1, 3, 4) =         0.0000054465 d0 !# 
      FSEfns( 80, -1, 3, 5) =        -0.0097830863 d0 !# Fermi    5.463 fm
      errtot( 80, -1, 3, 5) =         0.0000032543 d0 !# 
      FSEfns( 80, -1, 4, 4) =        -0.0097559817 d0 !# Fermi    5.463 fm
      errtot( 80, -1, 4, 4) =         0.0000213624 d0 !# 
      FSEfns( 80, -1, 4, 5) =        -0.0096183237 d0 !# Fermi    5.463 fm
      errtot( 80, -1, 4, 5) =         0.0000108680 d0 !# 
      FSEfns( 80, -1, 5, 5) =        -0.0095178598 d0 !# Fermi    5.463 fm
      errtot( 80, -1, 5, 5) =         0.0000062689 d0 !# 
      FSEfns( 85, -1, 1, 1) =        -0.0103954397 d0 !# Fermi    5.539 fm
      errtot( 85, -1, 1, 1) =         0.0000001288 d0 !# 
      FSEfns( 85, -1, 1, 2) =        -0.0125903933 d0 !# Fermi    5.539 fm
      errtot( 85, -1, 1, 2) =         0.0000000239 d0 !# 
      FSEfns( 85, -1, 1, 3) =        -0.0125387635 d0 !# Fermi    5.539 fm
      errtot( 85, -1, 1, 3) =         0.0000000322 d0 !# 
      FSEfns( 85, -1, 1, 4) =        -0.0123550234 d0 !# Fermi    5.539 fm
      errtot( 85, -1, 1, 4) =         0.0000001226 d0 !# 
      FSEfns( 85, -1, 1, 5) =        -0.0122010608 d0 !# Fermi    5.539 fm
      errtot( 85, -1, 1, 5) =         0.0000005494 d0 !# 
      FSEfns( 85, -1, 2, 2) =        -0.0152431342 d0 !# Fermi    5.539 fm
      errtot( 85, -1, 2, 2) =         0.0000003471 d0 !# 
      FSEfns( 85, -1, 2, 3) =        -0.0151304015 d0 !# Fermi    5.539 fm
      errtot( 85, -1, 2, 3) =         0.0000002650 d0 !# 
      FSEfns( 85, -1, 2, 4) =        -0.0149126387 d0 !# Fermi    5.539 fm
      errtot( 85, -1, 2, 4) =         0.0000000553 d0 !# 
      FSEfns( 85, -1, 2, 5) =        -0.0147296409 d0 !# Fermi    5.539 fm
      errtot( 85, -1, 2, 5) =         0.0000005324 d0 !# 
      FSEfns( 85, -1, 3, 3) =        -0.0150061273 d0 !# Fermi    5.539 fm
      errtot( 85, -1, 3, 3) =         0.0000019347 d0 !# 
      FSEfns( 85, -1, 3, 4) =        -0.0147484747 d0 !# Fermi    5.539 fm
      errtot( 85, -1, 3, 4) =         0.0000080209 d0 !# 
      FSEfns( 85, -1, 3, 5) =        -0.0145629340 d0 !# Fermi    5.539 fm
      errtot( 85, -1, 3, 5) =         0.0000050855 d0 !# 
      FSEfns( 85, -1, 4, 4) =        -0.0145194264 d0 !# Fermi    5.539 fm
      errtot( 85, -1, 4, 4) =         0.0000298119 d0 !# 
      FSEfns( 85, -1, 4, 5) =        -0.0142948651 d0 !# Fermi    5.539 fm
      errtot( 85, -1, 4, 5) =         0.0000165308 d0 !# 
      FSEfns( 85, -1, 5, 5) =        -0.0141237981 d0 !# Fermi    5.539 fm
      errtot( 85, -1, 5, 5) =         0.0000073506 d0 !# 
      FSEfns( 90, -1, 1, 1) =        -0.0154050611 d0 !# Fermi    5.710 fm
      errtot( 90, -1, 1, 1) =         0.0000002578 d0 !# 
      FSEfns( 90, -1, 1, 2) =        -0.0191143509 d0 !# Fermi    5.710 fm
      errtot( 90, -1, 1, 2) =         0.0000000427 d0 !# 
      FSEfns( 90, -1, 1, 3) =        -0.0190228068 d0 !# Fermi    5.710 fm
      errtot( 90, -1, 1, 3) =         0.0000000423 d0 !# 
      FSEfns( 90, -1, 1, 4) =        -0.0187070286 d0 !# Fermi    5.710 fm
      errtot( 90, -1, 1, 4) =         0.0000000483 d0 !# 
      FSEfns( 90, -1, 1, 5) =        -0.0184438929 d0 !# Fermi    5.710 fm
      errtot( 90, -1, 1, 5) =         0.0000007472 d0 !# 
      FSEfns( 90, -1, 2, 2) =        -0.0237389688 d0 !# Fermi    5.710 fm
      errtot( 90, -1, 2, 2) =         0.0000003655 d0 !# 
      FSEfns( 90, -1, 2, 3) =        -0.0235442622 d0 !# Fermi    5.710 fm
      errtot( 90, -1, 2, 3) =         0.0000003182 d0 !# 
      FSEfns( 90, -1, 2, 4) =        -0.0231588853 d0 !# Fermi    5.710 fm
      errtot( 90, -1, 2, 4) =         0.0000000979 d0 !# 
      FSEfns( 90, -1, 2, 5) =        -0.0228378112 d0 !# Fermi    5.710 fm
      errtot( 90, -1, 2, 5) =         0.0000007908 d0 !# 
      FSEfns( 90, -1, 3, 3) =        -0.0233382052 d0 !# Fermi    5.710 fm
      errtot( 90, -1, 3, 3) =         0.0000031005 d0 !# 
      FSEfns( 90, -1, 3, 4) =        -0.0228904573 d0 !# Fermi    5.710 fm
      errtot( 90, -1, 3, 4) =         0.0000123461 d0 !# 
      FSEfns( 90, -1, 3, 5) =        -0.0225649649 d0 !# Fermi    5.710 fm
      errtot( 90, -1, 3, 5) =         0.0000083350 d0 !# 
      FSEfns( 90, -1, 4, 4) =        -0.0224912801 d0 !# Fermi    5.710 fm
      errtot( 90, -1, 4, 4) =         0.0000429935 d0 !# 
      FSEfns( 90, -1, 4, 5) =        -0.0221081017 d0 !# Fermi    5.710 fm
      errtot( 90, -1, 4, 5) =         0.0000260841 d0 !# 
      FSEfns( 90, -1, 5, 5) =        -0.0218051428 d0 !# Fermi    5.710 fm
      errtot( 90, -1, 5, 5) =         0.0000085540 d0 !# 
      FSEfns( 95, -1, 1, 1) =        -0.0232784867 d0 !# Fermi    5.905 fm
      errtot( 95, -1, 1, 1) =         0.0000005364 d0 !# 
      FSEfns( 95, -1, 1, 2) =        -0.0296287584 d0 !# Fermi    5.905 fm
      errtot( 95, -1, 1, 2) =         0.0000004502 d0 !# 
      FSEfns( 95, -1, 1, 3) =        -0.0294523500 d0 !# Fermi    5.905 fm
      errtot( 95, -1, 1, 3) =         0.0000003652 d0 !# 
      FSEfns( 95, -1, 1, 4) =        -0.0288946627 d0 !# Fermi    5.905 fm
      errtot( 95, -1, 1, 4) =         0.0000000723 d0 !# 
      FSEfns( 95, -1, 1, 5) =        -0.0284341521 d0 !# Fermi    5.905 fm
      errtot( 95, -1, 1, 5) =         0.0000009889 d0 !# 
      FSEfns( 95, -1, 2, 2) =        -0.0378010798 d0 !# Fermi    5.905 fm
      errtot( 95, -1, 2, 2) =         0.0000003947 d0 !# 
      FSEfns( 95, -1, 2, 3) =        -0.0374449535 d0 !# Fermi    5.905 fm
      errtot( 95, -1, 2, 3) =         0.0000005737 d0 !# 
      FSEfns( 95, -1, 2, 4) =        -0.0367434569 d0 !# Fermi    5.905 fm
      errtot( 95, -1, 2, 4) =         0.0000002227 d0 !# 
      FSEfns( 95, -1, 2, 5) =        -0.0361657521 d0 !# Fermi    5.905 fm
      errtot( 95, -1, 2, 5) =         0.0000012097 d0 !# 
      FSEfns( 95, -1, 3, 3) =        -0.0370841663 d0 !# Fermi    5.905 fm
      errtot( 95, -1, 3, 3) =         0.0000051727 d0 !# 
      FSEfns( 95, -1, 3, 4) =        -0.0362851180 d0 !# Fermi    5.905 fm
      errtot( 95, -1, 3, 4) =         0.0000193916 d0 !# 
      FSEfns( 95, -1, 3, 5) =        -0.0356999279 d0 !# Fermi    5.905 fm
      errtot( 95, -1, 3, 5) =         0.0000139890 d0 !# 
      FSEfns( 95, -1, 4, 4) =        -0.0355706659 d0 !# Fermi    5.905 fm
      errtot( 95, -1, 4, 4) =         0.0000625930 d0 !# 
      FSEfns( 95, -1, 4, 5) =        -0.0349004028 d0 !# Fermi    5.905 fm
      errtot( 95, -1, 4, 5) =         0.0000415532 d0 !# 
      FSEfns( 95, -1, 5, 5) =        -0.0343531076 d0 !# Fermi    5.905 fm
      errtot( 95, -1, 5, 5) =         0.0000093011 d0 !# 
      FSEfns(100, -1, 1, 1) =        -0.0340600903 d0 !# Fermi    5.857 fm
      errtot(100, -1, 1, 1) =         0.0000007681 d0 !# 
      FSEfns(100, -1, 1, 2) =        -0.0445273923 d0 !# Fermi    5.857 fm
      errtot(100, -1, 1, 2) =         0.0000011675 d0 !# 
      FSEfns(100, -1, 1, 3) =        -0.0441872273 d0 !# Fermi    5.857 fm
      errtot(100, -1, 1, 3) =         0.0000009604 d0 !# 
      FSEfns(100, -1, 1, 4) =        -0.0432288627 d0 !# Fermi    5.857 fm
      errtot(100, -1, 1, 4) =         0.0000005734 d0 !# 
      FSEfns(100, -1, 1, 5) =        -0.0424468314 d0 !# Fermi    5.857 fm
      errtot(100, -1, 1, 5) =         0.0000014824 d0 !# 
      FSEfns(100, -1, 2, 2) =        -0.0584441368 d0 !# Fermi    5.857 fm
      errtot(100, -1, 2, 2) =         0.0000009660 d0 !# 
      FSEfns(100, -1, 2, 3) =        -0.0577991551 d0 !# Fermi    5.857 fm
      errtot(100, -1, 2, 3) =         0.0000009351 d0 !# 
      FSEfns(100, -1, 2, 4) =        -0.0565553052 d0 !# Fermi    5.857 fm
      errtot(100, -1, 2, 4) =         0.0000006623 d0 !# 
      FSEfns(100, -1, 2, 5) =        -0.0555452021 d0 !# Fermi    5.857 fm
      errtot(100, -1, 2, 5) =         0.0000016883 d0 !# 
      FSEfns(100, -1, 3, 3) =        -0.0571727383 d0 !# Fermi    5.857 fm
      errtot(100, -1, 3, 3) =         0.0000083692 d0 !# 
      FSEfns(100, -1, 3, 4) =        -0.0557863562 d0 !# Fermi    5.857 fm
      errtot(100, -1, 3, 4) =         0.0000292268 d0 !# 
      FSEfns(100, -1, 3, 5) =        -0.0547648029 d0 !# Fermi    5.857 fm
      errtot(100, -1, 3, 5) =         0.0000224319 d0 !# 
      FSEfns(100, -1, 4, 4) =        -0.0545433578 d0 !# Fermi    5.857 fm
      errtot(100, -1, 4, 4) =         0.0000864909 d0 !# 
      FSEfns(100, -1, 4, 5) =        -0.0534055503 d0 !# Fermi    5.857 fm
      errtot(100, -1, 4, 5) =         0.0000628191 d0 !# 
      FSEfns(100, -1, 5, 5) =        -0.0524495098 d0 !# Fermi    5.857 fm
      errtot(100, -1, 5, 5) =         0.0000099057 d0 !# 
      FSEfns(105, -1, 1, 1) =        -0.0523252115 d0 !# Fermi    5.919 fm
      errtot(105, -1, 1, 1) =         0.0000040693 d0 !# 
      FSEfns(105, -1, 1, 2) =        -0.0703378615 d0 !# Fermi    5.919 fm
      errtot(105, -1, 1, 2) =         0.0000033229 d0 !# 
      FSEfns(105, -1, 1, 3) =        -0.0696221657 d0 !# Fermi    5.919 fm
      errtot(105, -1, 1, 3) =         0.0000024747 d0 !# 
      FSEfns(105, -1, 1, 4) =        -0.0678777394 d0 !# Fermi    5.919 fm
      errtot(105, -1, 1, 4) =         0.0000007360 d0 !# 
      FSEfns(105, -1, 1, 5) =        -0.0664763188 d0 !# Fermi    5.919 fm
      errtot(105, -1, 1, 5) =         0.0000028243 d0 !# 
      FSEfns(105, -1, 2, 2) =        -0.0951103852 d0 !# Fermi    5.919 fm
      errtot(105, -1, 2, 2) =         0.0000027795 d0 !# 
      FSEfns(105, -1, 2, 3) =        -0.0938403600 d0 !# Fermi    5.919 fm
      errtot(105, -1, 2, 3) =         0.0000080720 d0 !# 
      FSEfns(105, -1, 2, 4) =        -0.0915035519 d0 !# Fermi    5.919 fm
      errtot(105, -1, 2, 4) =         0.0000005171 d0 !# 
      FSEfns(105, -1, 2, 5) =        -0.0896373586 d0 !# Fermi    5.919 fm
      errtot(105, -1, 2, 5) =         0.0000041405 d0 !# 
      FSEfns(105, -1, 3, 3) =        -0.0926516227 d0 !# Fermi    5.919 fm
      errtot(105, -1, 3, 3) =         0.0000030002 d0 !# 
      FSEfns(105, -1, 3, 4) =        -0.0901193048 d0 !# Fermi    5.919 fm
      errtot(105, -1, 3, 4) =         0.0000055463 d0 !# 
      FSEfns(105, -1, 3, 5) =        -0.0882294727 d0 !# Fermi    5.919 fm
      errtot(105, -1, 3, 5) =         0.0000054205 d0 !# 
      FSEfns(105, -1, 4, 4) =        -0.0877491601 d0 !# Fermi    5.919 fm
      errtot(105, -1, 4, 4) =         0.0000287486 d0 !# 
      FSEfns(105, -1, 4, 5) =        -0.0857790545 d0 !# Fermi    5.919 fm
      errtot(105, -1, 4, 5) =         0.0000980557 d0 !# 
      FSEfns(105, -1, 5, 5) =        -0.0840228749 d0 !# Fermi    5.919 fm
      errtot(105, -1, 5, 5) =         0.0000269353 d0 !# 
      FSEfns(110, -1, 1, 1) =        -0.0829910966 d0 !# Fermi    5.993 fm
      errtot(110, -1, 1, 1) =         0.0000054704 d0 !# 
      FSEfns(110, -1, 1, 2) =        -0.1148127810 d0 !# Fermi    5.993 fm
      errtot(110, -1, 1, 2) =         0.0000096329 d0 !# 
      FSEfns(110, -1, 1, 3) =        -0.1132291640 d0 !# Fermi    5.993 fm
      errtot(110, -1, 1, 3) =         0.0000067639 d0 !# 
      FSEfns(110, -1, 1, 4) =        -0.1099257440 d0 !# Fermi    5.993 fm
      errtot(110, -1, 1, 4) =         0.0000074484 d0 !# 
      FSEfns(110, -1, 1, 5) =        -0.1073223624 d0 !# Fermi    5.993 fm
      errtot(110, -1, 1, 5) =         0.0000051531 d0 !# 
      FSEfns(110, -1, 2, 2) =        -0.1601463597 d0 !# Fermi    5.993 fm
      errtot(110, -1, 2, 2) =         0.0000078291 d0 !# 
      FSEfns(110, -1, 2, 3) =        -0.1574996664 d0 !# Fermi    5.993 fm
      errtot(110, -1, 2, 3) =         0.0000030038 d0 !# 
      FSEfns(110, -1, 2, 4) =        -0.1529315378 d0 !# Fermi    5.993 fm
      errtot(110, -1, 2, 4) =         0.0000100212 d0 !# 
      FSEfns(110, -1, 2, 5) =        -0.1493536549 d0 !# Fermi    5.993 fm
      errtot(110, -1, 2, 5) =         0.0000063673 d0 !# 
      FSEfns(110, -1, 3, 3) =        -0.1551018295 d0 !# Fermi    5.993 fm
      errtot(110, -1, 3, 3) =         0.0000064415 d0 !# 
      FSEfns(110, -1, 3, 4) =        -0.1502657803 d0 !# Fermi    5.993 fm
      errtot(110, -1, 3, 4) =         0.0000092625 d0 !# 
      FSEfns(110, -1, 3, 5) =        -0.1466615841 d0 !# Fermi    5.993 fm
      errtot(110, -1, 3, 5) =         0.0000058706 d0 !# 
      FSEfns(110, -1, 4, 4) =        -0.1457543813 d0 !# Fermi    5.993 fm
      errtot(110, -1, 4, 4) =         0.0000538675 d0 !# 
      FSEfns(110, -1, 4, 5) =        -0.1420770982 d0 !# Fermi    5.993 fm
      errtot(110, -1, 4, 5) =         0.0001548883 d0 !# 
      FSEfns(110, -1, 5, 5) =        -0.1387489851 d0 !# Fermi    5.993 fm
      errtot(110, -1, 5, 5) =         0.0000294270 d0 !# 
      FSEfns(115, -1, 1, 1) =        -0.1370180896 d0 !# Fermi    6.088 fm
      errtot(115, -1, 1, 1) =         0.0000230103 d0 !# 
      FSEfns(115, -1, 1, 2) =        -0.1951714261 d0 !# Fermi    6.088 fm
      errtot(115, -1, 1, 2) =         0.0000298042 d0 !# 
      FSEfns(115, -1, 1, 3) =        -0.1914893533 d0 !# Fermi    6.088 fm
      errtot(115, -1, 1, 3) =         0.0000150155 d0 !# 
      FSEfns(115, -1, 1, 4) =        -0.1849295101 d0 !# Fermi    6.088 fm
      errtot(115, -1, 1, 4) =         0.0000179736 d0 !# 
      FSEfns(115, -1, 1, 5) =        -0.1798765340 d0 !# Fermi    6.088 fm
      errtot(115, -1, 1, 5) =         0.0000129208 d0 !# 
      FSEfns(115, -1, 2, 2) =        -0.2811342014 d0 !# Fermi    6.088 fm
      errtot(115, -1, 2, 2) =         0.0000257442 d0 !# 
      FSEfns(115, -1, 2, 3) =        -0.2752797658 d0 !# Fermi    6.088 fm
      errtot(115, -1, 2, 3) =         0.0000412525 d0 !# 
      FSEfns(115, -1, 2, 4) =        -0.2659139858 d0 !# Fermi    6.088 fm
      errtot(115, -1, 2, 4) =         0.0000205667 d0 !# 
      FSEfns(115, -1, 2, 5) =        -0.2587389687 d0 !# Fermi    6.088 fm
      errtot(115, -1, 2, 5) =         0.0000185591 d0 !# 
      FSEfns(115, -1, 3, 3) =        -0.2701051898 d0 !# Fermi    6.088 fm
      errtot(115, -1, 3, 3) =         0.0000256563 d0 !# 
      FSEfns(115, -1, 3, 4) =        -0.2604362976 d0 !# Fermi    6.088 fm
      errtot(115, -1, 3, 4) =         0.0000097992 d0 !# 
      FSEfns(115, -1, 3, 5) =        -0.2532548152 d0 !# Fermi    6.088 fm
      errtot(115, -1, 3, 5) =         0.0000172053 d0 !# 
      FSEfns(115, -1, 4, 4) =        -0.2514528005 d0 !# Fermi    6.088 fm
      errtot(115, -1, 4, 4) =         0.0000913055 d0 !# 
      FSEfns(115, -1, 4, 5) =        -0.2442302399 d0 !# Fermi    6.088 fm
      errtot(115, -1, 4, 5) =         0.0000492932 d0 !# 
      FSEfns(115, -1, 5, 5) =        -0.2377012762 d0 !# Fermi    6.088 fm
      errtot(115, -1, 5, 5) =         0.0000227588 d0 !# 
      FSEfns(120, -1, 1, 1) =        -0.2373357540 d0 !# Fermi    6.175 fm
      errtot(120, -1, 1, 1) =         0.0000721569 d0 !# 
      FSEfns(120, -1, 1, 2) =        -0.3480506083 d0 !# Fermi    6.175 fm
      errtot(120, -1, 1, 2) =         0.0001033663 d0 !# 
      FSEfns(120, -1, 1, 3) =        -0.3390523486 d0 !# Fermi    6.175 fm
      errtot(120, -1, 1, 3) =         0.0000524825 d0 !# 
      FSEfns(120, -1, 1, 4) =        -0.3252957280 d0 !# Fermi    6.175 fm
      errtot(120, -1, 1, 4) =         0.0000713717 d0 !# 
      FSEfns(120, -1, 1, 5) =        -0.3149743336 d0 !# Fermi    6.175 fm
      errtot(120, -1, 1, 5) =         0.0000387187 d0 !# 
      FSEfns(120, -1, 2, 2) =        -0.5181665881 d0 !# Fermi    6.175 fm
      errtot(120, -1, 2, 2) =         0.0000871056 d0 !# 
      FSEfns(120, -1, 2, 3) =        -0.5043772409 d0 !# Fermi    6.175 fm
      errtot(120, -1, 2, 3) =         0.0000951408 d0 !# 
      FSEfns(120, -1, 2, 4) =        -0.4841123727 d0 !# Fermi    6.175 fm
      errtot(120, -1, 2, 4) =         0.0000793829 d0 !# 
      FSEfns(120, -1, 2, 5) =        -0.4689516532 d0 !# Fermi    6.175 fm
      errtot(120, -1, 2, 5) =         0.0000590038 d0 !# 
      FSEfns(120, -1, 3, 3) =        -0.4924167100 d0 !# Fermi    6.175 fm
      errtot(120, -1, 3, 3) =         0.0000744635 d0 !# 
      FSEfns(120, -1, 3, 4) =        -0.4720298101 d0 !# Fermi    6.175 fm
      errtot(120, -1, 3, 4) =         0.0000458188 d0 !# 
      FSEfns(120, -1, 3, 5) =        -0.4569886322 d0 !# Fermi    6.175 fm
      errtot(120, -1, 3, 5) =         0.0000525072 d0 !# 
      FSEfns(120, -1, 4, 4) =        -0.4531494130 d0 !# Fermi    6.175 fm
      errtot(120, -1, 4, 4) =         0.0001857706 d0 !# 
      FSEfns(120, -1, 4, 5) =        -0.4383367202 d0 !# Fermi    6.175 fm
      errtot(120, -1, 4, 5) =         0.0000702882 d0 !# 
      FSEfns(120, -1, 5, 5) =        -0.4248184784 d0 !# Fermi    6.175 fm
      errtot(120, -1, 5, 5) =         0.0000251700 d0 !# 
*
* np_1/2
*
      FSEfns( 60,  1, 2, 2) =        -0.0000460267 d0 !# Fermi    4.912 fm
      errtot( 60,  1, 2, 2) =         0.0000000040 d0 !# 
      FSEfns( 60,  1, 2, 3) =        -0.0000510722 d0 !# Fermi    4.912 fm
      errtot( 60,  1, 2, 3) =         0.0000000023 d0 !# 
      FSEfns( 60,  1, 2, 4) =        -0.0000517128 d0 !# Fermi    4.912 fm
      errtot( 60,  1, 2, 4) =         0.0000000031 d0 !# 
      FSEfns( 60,  1, 2, 5) =        -0.0000518387 d0 !# Fermi    4.912 fm
      errtot( 60,  1, 2, 5) =         0.0000000030 d0 !# 
      FSEfns( 60,  1, 3, 3) =        -0.0000560101 d0 !# Fermi    4.912 fm
      errtot( 60,  1, 3, 3) =         0.0000000032 d0 !# 
      FSEfns( 60,  1, 3, 4) =        -0.0000578867 d0 !# Fermi    4.912 fm
      errtot( 60,  1, 3, 4) =         0.0000000536 d0 !# 
      FSEfns( 60,  1, 3, 5) =        -0.0000581169 d0 !# Fermi    4.912 fm
      errtot( 60,  1, 3, 5) =         0.0000000235 d0 !# 
      FSEfns( 60,  1, 4, 4) =        -0.0000593483 d0 !# Fermi    4.912 fm
      errtot( 60,  1, 4, 4) =         0.0000002868 d0 !# 
      FSEfns( 60,  1, 4, 5) =        -0.0000600701 d0 !# Fermi    4.912 fm
      errtot( 60,  1, 4, 5) =         0.0000000703 d0 !# 
      FSEfns( 60,  1, 5, 5) =        -0.0000608095 d0 !# Fermi    4.912 fm
      errtot( 60,  1, 5, 5) =         0.0000000940 d0 !# 
      FSEfns( 65,  1, 2, 2) =        -0.0000883856 d0 !# Fermi    5.060 fm
      errtot( 65,  1, 2, 2) =         0.0000000091 d0 !# 
      FSEfns( 65,  1, 2, 3) =        -0.0000980020 d0 !# Fermi    5.060 fm
      errtot( 65,  1, 2, 3) =         0.0000000067 d0 !# 
      FSEfns( 65,  1, 2, 4) =        -0.0000994580 d0 !# Fermi    5.060 fm
      errtot( 65,  1, 2, 4) =         0.0000000034 d0 !# 
      FSEfns( 65,  1, 2, 5) =        -0.0000997816 d0 !# Fermi    5.060 fm
      errtot( 65,  1, 2, 5) =         0.0000000051 d0 !# 
      FSEfns( 65,  1, 3, 3) =        -0.0001074436 d0 !# Fermi    5.060 fm
      errtot( 65,  1, 3, 3) =         0.0000000091 d0 !# 
      FSEfns( 65,  1, 3, 4) =        -0.0001107595 d0 !# Fermi    5.060 fm
      errtot( 65,  1, 3, 4) =         0.0000000973 d0 !# 
      FSEfns( 65,  1, 3, 5) =        -0.0001112452 d0 !# Fermi    5.060 fm
      errtot( 65,  1, 3, 5) =         0.0000000433 d0 !# 
      FSEfns( 65,  1, 4, 4) =        -0.0001134524 d0 !# Fermi    5.060 fm
      errtot( 65,  1, 4, 4) =         0.0000004890 d0 !# 
      FSEfns( 65,  1, 4, 5) =        -0.0001146174 d0 !# Fermi    5.060 fm
      errtot( 65,  1, 4, 5) =         0.0000001412 d0 !# 
      FSEfns( 65,  1, 5, 5) =        -0.0001158314 d0 !# Fermi    5.060 fm
      errtot( 65,  1, 5, 5) =         0.0000002716 d0 !# 
      FSEfns( 70,  1, 2, 2) =        -0.0001704541 d0 !# Fermi    5.311 fm
      errtot( 70,  1, 2, 2) =         0.0000000164 d0 !# 
      FSEfns( 70,  1, 2, 3) =        -0.0001889211 d0 !# Fermi    5.311 fm
      errtot( 70,  1, 2, 3) =         0.0000000146 d0 !# 
      FSEfns( 70,  1, 2, 4) =        -0.0001919592 d0 !# Fermi    5.311 fm
      errtot( 70,  1, 2, 4) =         0.0000000027 d0 !# 
      FSEfns( 70,  1, 2, 5) =        -0.0001926249 d0 !# Fermi    5.311 fm
      errtot( 70,  1, 2, 5) =         0.0000000080 d0 !# 
      FSEfns( 70,  1, 3, 3) =        -0.0002070891 d0 !# Fermi    5.311 fm
      errtot( 70,  1, 3, 3) =         0.0000000192 d0 !# 
      FSEfns( 70,  1, 3, 4) =        -0.0002130078 d0 !# Fermi    5.311 fm
      errtot( 70,  1, 3, 4) =         0.0000001669 d0 !# 
      FSEfns( 70,  1, 3, 5) =        -0.0002139224 d0 !# Fermi    5.311 fm
      errtot( 70,  1, 3, 5) =         0.0000000756 d0 !# 
      FSEfns( 70,  1, 4, 4) =        -0.0002179293 d0 !# Fermi    5.311 fm
      errtot( 70,  1, 4, 4) =         0.0000008389 d0 !# 
      FSEfns( 70,  1, 4, 5) =        -0.0002198055 d0 !# Fermi    5.311 fm
      errtot( 70,  1, 4, 5) =         0.0000002486 d0 !# 
      FSEfns( 70,  1, 5, 5) =        -0.0002219483 d0 !# Fermi    5.311 fm
      errtot( 70,  1, 5, 5) =         0.0000002352 d0 !# 
      FSEfns( 75,  1, 2, 2) =        -0.0003035259 d0 !# Fermi    5.339 fm
      errtot( 75,  1, 2, 2) =         0.0000000227 d0 !# 
      FSEfns( 75,  1, 2, 3) =        -0.0003362634 d0 !# Fermi    5.339 fm
      errtot( 75,  1, 2, 3) =         0.0000000212 d0 !# 
      FSEfns( 75,  1, 2, 4) =        -0.0003418416 d0 !# Fermi    5.339 fm
      errtot( 75,  1, 2, 4) =         0.0000000057 d0 !# 
      FSEfns( 75,  1, 2, 5) =        -0.0003429585 d0 !# Fermi    5.339 fm
      errtot( 75,  1, 2, 5) =         0.0000000135 d0 !# 
      FSEfns( 75,  1, 3, 3) =        -0.0003685455 d0 !# Fermi    5.339 fm
      errtot( 75,  1, 3, 3) =         0.0000000310 d0 !# 
      FSEfns( 75,  1, 3, 4) =        -0.0003782520 d0 !# Fermi    5.339 fm
      errtot( 75,  1, 3, 4) =         0.0000002721 d0 !# 
      FSEfns( 75,  1, 3, 5) =        -0.0003797105 d0 !# Fermi    5.339 fm
      errtot( 75,  1, 3, 5) =         0.0000001279 d0 !# 
      FSEfns( 75,  1, 4, 4) =        -0.0003864830 d0 !# Fermi    5.339 fm
      errtot( 75,  1, 4, 4) =         0.0000013449 d0 !# 
      FSEfns( 75,  1, 4, 5) =        -0.0003891701 d0 !# Fermi    5.339 fm
      errtot( 75,  1, 4, 5) =         0.0000004327 d0 !# 
      FSEfns( 75,  1, 5, 5) =        -0.0003924099 d0 !# Fermi    5.339 fm
      errtot( 75,  1, 5, 5) =         0.0000003472 d0 !# 
      FSEfns( 80,  1, 2, 2) =        -0.0005545012 d0 !# Fermi    5.463 fm
      errtot( 80,  1, 2, 2) =         0.0000000354 d0 !# 
      FSEfns( 80,  1, 2, 3) =        -0.0006140147 d0 !# Fermi    5.463 fm
      errtot( 80,  1, 2, 3) =         0.0000000284 d0 !# 
      FSEfns( 80,  1, 2, 4) =        -0.0006241286 d0 !# Fermi    5.463 fm
      errtot( 80,  1, 2, 4) =         0.0000000147 d0 !# 
      FSEfns( 80,  1, 2, 5) =        -0.0006258120 d0 !# Fermi    5.463 fm
      errtot( 80,  1, 2, 5) =         0.0000000232 d0 !# 
      FSEfns( 80,  1, 3, 3) =        -0.0006728207 d0 !# Fermi    5.463 fm
      errtot( 80,  1, 3, 3) =         0.0000000420 d0 !# 
      FSEfns( 80,  1, 3, 4) =        -0.0006889994 d0 !# Fermi    5.463 fm
      errtot( 80,  1, 3, 4) =         0.0000004674 d0 !# 
      FSEfns( 80,  1, 3, 5) =        -0.0006911319 d0 !# Fermi    5.463 fm
      errtot( 80,  1, 3, 5) =         0.0000002309 d0 !# 
      FSEfns( 80,  1, 4, 4) =        -0.0007029403 d0 !# Fermi    5.463 fm
      errtot( 80,  1, 4, 4) =         0.0000022408 d0 !# 
      FSEfns( 80,  1, 4, 5) =        -0.0007066137 d0 !# Fermi    5.463 fm
      errtot( 80,  1, 4, 5) =         0.0000008063 d0 !# 
      FSEfns( 80,  1, 5, 5) =        -0.0007113658 d0 !# Fermi    5.463 fm
      errtot( 80,  1, 5, 5) =         0.0000005630 d0 !# 
      FSEfns( 85,  1, 2, 2) =        -0.0010044909 d0 !# Fermi    5.539 fm
      errtot( 85,  1, 2, 2) =         0.0000000504 d0 !# 
      FSEfns( 85,  1, 2, 3) =        -0.0011115153 d0 !# Fermi    5.539 fm
      errtot( 85,  1, 2, 3) =         0.0000000398 d0 !# 
      FSEfns( 85,  1, 2, 4) =        -0.0011290666 d0 !# Fermi    5.539 fm
      errtot( 85,  1, 2, 4) =         0.0000000229 d0 !# 
      FSEfns( 85,  1, 2, 5) =        -0.0011310872 d0 !# Fermi    5.539 fm
      errtot( 85,  1, 2, 5) =         0.0000000400 d0 !# 
      FSEfns( 85,  1, 3, 3) =        -0.0012174922 d0 !# Fermi    5.539 fm
      errtot( 85,  1, 3, 3) =         0.0000000821 d0 !# 
      FSEfns( 85,  1, 3, 4) =        -0.0012437587 d0 !# Fermi    5.539 fm
      errtot( 85,  1, 3, 4) =         0.0000008063 d0 !# 
      FSEfns( 85,  1, 3, 5) =        -0.0012462736 d0 !# Fermi    5.539 fm
      errtot( 85,  1, 3, 5) =         0.0000004253 d0 !# 
      FSEfns( 85,  1, 4, 4) =        -0.0012667057 d0 !# Fermi    5.539 fm
      errtot( 85,  1, 4, 4) =         0.0000037085 d0 !# 
      FSEfns( 85,  1, 4, 5) =        -0.0012709666 d0 !# Fermi    5.539 fm
      errtot( 85,  1, 4, 5) =         0.0000015088 d0 !# 
      FSEfns( 85,  1, 5, 5) =        -0.0012773232 d0 !# Fermi    5.539 fm
      errtot( 85,  1, 5, 5) =         0.0000008361 d0 !# 
      FSEfns( 90,  1, 2, 2) =        -0.0018832560 d0 !# Fermi    5.710 fm
      errtot( 90,  1, 2, 2) =         0.0000000721 d0 !# 
      FSEfns( 90,  1, 2, 3) =        -0.0020817064 d0 !# Fermi    5.710 fm
      errtot( 90,  1, 2, 3) =         0.0000000929 d0 !# 
      FSEfns( 90,  1, 2, 4) =        -0.0021119769 d0 !# Fermi    5.710 fm
      errtot( 90,  1, 2, 4) =         0.0000000431 d0 !# 
      FSEfns( 90,  1, 2, 5) =        -0.0021130957 d0 !# Fermi    5.710 fm
      errtot( 90,  1, 2, 5) =         0.0000000779 d0 !# 
      FSEfns( 90,  1, 3, 3) =        -0.0022786544 d0 !# Fermi    5.710 fm
      errtot( 90,  1, 3, 3) =         0.0000003917 d0 !# 
      FSEfns( 90,  1, 3, 4) =        -0.0023215471 d0 !# Fermi    5.710 fm
      errtot( 90,  1, 3, 4) =         0.0000014542 d0 !# 
      FSEfns( 90,  1, 3, 5) =        -0.0023229715 d0 !# Fermi    5.710 fm
      errtot( 90,  1, 3, 5) =         0.0000008227 d0 !# 
      FSEfns( 90,  1, 4, 4) =        -0.0023594727 d0 !# Fermi    5.710 fm
      errtot( 90,  1, 4, 4) =         0.0000063522 d0 !# 
      FSEfns( 90,  1, 4, 5) =        -0.0023625669 d0 !# Fermi    5.710 fm
      errtot( 90,  1, 4, 5) =         0.0000029236 d0 !# 
      FSEfns( 90,  1, 5, 5) =        -0.0023698002 d0 !# Fermi    5.710 fm
      errtot( 90,  1, 5, 5) =         0.0000012469 d0 !# 
      FSEfns( 95,  1, 2, 2) =        -0.0036011916 d0 !# Fermi    5.905 fm
      errtot( 95,  1, 2, 2) =         0.0000000471 d0 !# 
      FSEfns( 95,  1, 2, 3) =        -0.0039742375 d0 !# Fermi    5.905 fm
      errtot( 95,  1, 2, 3) =         0.0000000574 d0 !# 
      FSEfns( 95,  1, 2, 4) =        -0.0040245664 d0 !# Fermi    5.905 fm
      errtot( 95,  1, 2, 4) =         0.0000000571 d0 !# 
      FSEfns( 95,  1, 2, 5) =        -0.0040200777 d0 !# Fermi    5.905 fm
      errtot( 95,  1, 2, 5) =         0.0000001286 d0 !# 
      FSEfns( 95,  1, 3, 3) =        -0.0043452971 d0 !# Fermi    5.905 fm
      errtot( 95,  1, 3, 3) =         0.0000004560 d0 !# 
      FSEfns( 95,  1, 3, 4) =        -0.0044134236 d0 !# Fermi    5.905 fm
      errtot( 95,  1, 3, 4) =         0.0000026813 d0 !# 
      FSEfns( 95,  1, 3, 5) =        -0.0044082665 d0 !# Fermi    5.905 fm
      errtot( 95,  1, 3, 5) =         0.0000016307 d0 !# 
      FSEfns( 95,  1, 4, 4) =        -0.0044743118 d0 !# Fermi    5.905 fm
      errtot( 95,  1, 4, 4) =         0.0000110163 d0 !# 
      FSEfns( 95,  1, 4, 5) =        -0.0044698694 d0 !# Fermi    5.905 fm
      errtot( 95,  1, 4, 5) =         0.0000057085 d0 !# 
      FSEfns( 95,  1, 5, 5) =        -0.0044735859 d0 !# Fermi    5.905 fm
      errtot( 95,  1, 5, 5) =         0.0000018198 d0 !# 
      FSEfns(100,  1, 2, 2) =        -0.0066963332 d0 !# Fermi    5.857 fm
      errtot(100,  1, 2, 2) =         0.0000000355 d0 !# 
      FSEfns(100,  1, 2, 3) =        -0.0073717497 d0 !# Fermi    5.857 fm
      errtot(100,  1, 2, 3) =         0.0000005209 d0 !# 
      FSEfns(100,  1, 2, 4) =        -0.0074459787 d0 !# Fermi    5.857 fm
      errtot(100,  1, 2, 4) =         0.0000000773 d0 !# 
      FSEfns(100,  1, 2, 5) =        -0.0074222279 d0 !# Fermi    5.857 fm
      errtot(100,  1, 2, 5) =         0.0000002499 d0 !# 
      FSEfns(100,  1, 3, 3) =        -0.0080454070 d0 !# Fermi    5.857 fm
      errtot(100,  1, 3, 3) =         0.0000027178 d0 !# 
      FSEfns(100,  1, 3, 4) =        -0.0081420420 d0 !# Fermi    5.857 fm
      errtot(100,  1, 3, 4) =         0.0000047517 d0 !# 
      FSEfns(100,  1, 3, 5) =        -0.0081145808 d0 !# Fermi    5.857 fm
      errtot(100,  1, 3, 5) =         0.0000031056 d0 !# 
      FSEfns(100,  1, 4, 4) =        -0.0082291928 d0 !# Fermi    5.857 fm
      errtot(100,  1, 4, 4) =         0.0000182078 d0 !# 
      FSEfns(100,  1, 4, 5) =        -0.0081996691 d0 !# Fermi    5.857 fm
      errtot(100,  1, 4, 5) =         0.0000105535 d0 !# 
      FSEfns(100,  1, 5, 5) =        -0.0081854818 d0 !# Fermi    5.857 fm
      errtot(100,  1, 5, 5) =         0.0000023247 d0 !# 
      FSEfns(105,  1, 2, 2) =        -0.0131673991 d0 !# Fermi    5.919 fm
      errtot(105,  1, 2, 2) =         0.0000003382 d0 !# 
      FSEfns(105,  1, 2, 3) =        -0.0144435510 d0 !# Fermi    5.919 fm
      errtot(105,  1, 2, 3) =         0.0000005980 d0 !# 
      FSEfns(105,  1, 2, 4) =        -0.0145386068 d0 !# Fermi    5.919 fm
      errtot(105,  1, 2, 4) =         0.0000002982 d0 !# 
      FSEfns(105,  1, 2, 5) =        -0.0144542828 d0 !# Fermi    5.919 fm
      errtot(105,  1, 2, 5) =         0.0000006908 d0 !# 
      FSEfns(105,  1, 3, 3) =        -0.0157193791 d0 !# Fermi    5.919 fm
      errtot(105,  1, 3, 3) =         0.0000092416 d0 !# 
      FSEfns(105,  1, 3, 4) =        -0.0158395088 d0 !# Fermi    5.919 fm
      errtot(105,  1, 3, 4) =         0.0000088679 d0 !# 
      FSEfns(105,  1, 3, 5) =        -0.0157428982 d0 !# Fermi    5.919 fm
      errtot(105,  1, 3, 5) =         0.0000062708 d0 !# 
      FSEfns(105,  1, 4, 4) =        -0.0159493280 d0 !# Fermi    5.919 fm
      errtot(105,  1, 4, 4) =         0.0000313142 d0 !# 
      FSEfns(105,  1, 4, 5) =        -0.0158438865 d0 !# Fermi    5.919 fm
      errtot(105,  1, 4, 5) =         0.0000202036 d0 !# 
      FSEfns(105,  1, 5, 5) =        -0.0157688074 d0 !# Fermi    5.919 fm
      errtot(105,  1, 5, 5) =         0.0000027837 d0 !# 
      FSEfns(110,  1, 2, 2) =        -0.0270370673 d0 !# Fermi    5.993 fm
      errtot(110,  1, 2, 2) =         0.0000013861 d0 !# 
      FSEfns(110,  1, 2, 3) =        -0.0295024189 d0 !# Fermi    5.993 fm
      errtot(110,  1, 2, 3) =         0.0000021870 d0 !# 
      FSEfns(110,  1, 2, 4) =        -0.0295607652 d0 !# Fermi    5.993 fm
      errtot(110,  1, 2, 4) =         0.0000026391 d0 !# 
      FSEfns(110,  1, 2, 5) =        -0.0292927021 d0 !# Fermi    5.993 fm
      errtot(110,  1, 2, 5) =         0.0000017169 d0 !# 
      FSEfns(110,  1, 3, 3) =        -0.0319711919 d0 !# Fermi    5.993 fm
      errtot(110,  1, 3, 3) =         0.0000169058 d0 !# 
      FSEfns(110,  1, 3, 4) =        -0.0320471952 d0 !# Fermi    5.993 fm
      errtot(110,  1, 3, 4) =         0.0000170474 d0 !# 
      FSEfns(110,  1, 3, 5) =        -0.0317437148 d0 !# Fermi    5.993 fm
      errtot(110,  1, 3, 5) =         0.0000129010 d0 !# 
      FSEfns(110,  1, 4, 4) =        -0.0321152395 d0 !# Fermi    5.993 fm
      errtot(110,  1, 4, 4) =         0.0000552878 d0 !# 
      FSEfns(110,  1, 4, 5) =        -0.0317903845 d0 !# Fermi    5.993 fm
      errtot(110,  1, 4, 5) =         0.0000393514 d0 !# 
      FSEfns(110,  1, 5, 5) =        -0.0315280209 d0 !# Fermi    5.993 fm
      errtot(110,  1, 5, 5) =         0.0000024861 d0 !# 
      FSEfns(115,  1, 2, 2) =        -0.0587211601 d0 !# Fermi    6.088 fm
      errtot(115,  1, 2, 2) =         0.0000058611 d0 !# 
      FSEfns(115,  1, 2, 3) =        -0.0635871710 d0 !# Fermi    6.088 fm
      errtot(115,  1, 2, 3) =         0.0000134671 d0 !# 
      FSEfns(115,  1, 2, 4) =        -0.0633263601 d0 !# Fermi    6.088 fm
      errtot(115,  1, 2, 4) =         0.0000105344 d0 !# 
      FSEfns(115,  1, 2, 5) =        -0.0624915145 d0 !# Fermi    6.088 fm
      errtot(115,  1, 2, 5) =         0.0000048769 d0 !# 
      FSEfns(115,  1, 3, 3) =        -0.0684636708 d0 !# Fermi    6.088 fm
      errtot(115,  1, 3, 3) =         0.0000063875 d0 !# 
      FSEfns(115,  1, 3, 4) =        -0.0681787386 d0 !# Fermi    6.088 fm
      errtot(115,  1, 3, 4) =         0.0000342403 d0 !# 
      FSEfns(115,  1, 3, 5) =        -0.0672470060 d0 !# Fermi    6.088 fm
      errtot(115,  1, 3, 5) =         0.0000277651 d0 !# 
      FSEfns(115,  1, 4, 4) =        -0.0679209152 d0 !# Fermi    6.088 fm
      errtot(115,  1, 4, 4) =         0.0000980558 d0 !# 
      FSEfns(115,  1, 4, 5) =        -0.0669423946 d0 !# Fermi    6.088 fm
      errtot(115,  1, 4, 5) =         0.0000783455 d0 !# 
      FSEfns(115,  1, 5, 5) =        -0.0660991092 d0 !# Fermi    6.088 fm
      errtot(115,  1, 5, 5) =         0.0000025238 d0 !# 
      FSEfns(120,  1, 2, 2) =        -0.1367388667 d0 !# Fermi    6.175 fm
      errtot(120,  1, 2, 2) =         0.0000240498 d0 !# 
      FSEfns(120,  1, 2, 3) =        -0.1464018402 d0 !# Fermi    6.175 fm
      errtot(120,  1, 2, 3) =         0.0000023019 d0 !# 
      FSEfns(120,  1, 2, 4) =        -0.1446238134 d0 !# Fermi    6.175 fm
      errtot(120,  1, 2, 4) =         0.0000278291 d0 !# 
      FSEfns(120,  1, 2, 5) =        -0.1419595534 d0 !# Fermi    6.175 fm
      errtot(120,  1, 2, 5) =         0.0000166647 d0 !# 
      FSEfns(120,  1, 3, 3) =        -0.1560815728 d0 !# Fermi    6.175 fm
      errtot(120,  1, 3, 3) =         0.0000190062 d0 !# 
      FSEfns(120,  1, 3, 4) =        -0.1541428160 d0 !# Fermi    6.175 fm
      errtot(120,  1, 3, 4) =         0.0000713092 d0 !# 
      FSEfns(120,  1, 3, 5) =        -0.1512230030 d0 !# Fermi    6.175 fm
      errtot(120,  1, 3, 5) =         0.0000630258 d0 !# 
      FSEfns(120,  1, 4, 4) =        -0.1523627479 d0 !# Fermi    6.175 fm
      errtot(120,  1, 4, 4) =         0.0001773380 d0 !# 
      FSEfns(120,  1, 4, 5) =        -0.1493759579 d0 !# Fermi    6.175 fm
      errtot(120,  1, 4, 5) =         0.0001616935 d0 !# 
      FSEfns(120,  1, 5, 5) =        -0.1467001003 d0 !# Fermi    6.175 fm
      errtot(120,  1, 5, 5) =         0.0000197602 d0 !# 
*
* np_3/2
*
      FSEfns( 75, -2, 2, 2) =        -0.0000743395 d0 !# Fermi    5.339 fm
      errtot( 75, -2, 2, 2) =         0.0000000285 d0 !# 
      FSEfns( 75, -2, 2, 3) =        -0.0000892171 d0 !# Fermi    5.339 fm
      errtot( 75, -2, 2, 3) =         0.0000001036 d0 !# 
      FSEfns( 75, -2, 2, 4) =        -0.0000906607 d0 !# Fermi    5.339 fm
      errtot( 75, -2, 2, 4) =         0.0000000159 d0 !# 
      FSEfns( 75, -2, 2, 5) =        -0.0000908090 d0 !# Fermi    5.339 fm
      errtot( 75, -2, 2, 5) =         0.0000001513 d0 !# 
      FSEfns( 75, -2, 3, 3) =        -0.0000999990 d0 !# Fermi    5.339 fm
      errtot( 75, -2, 3, 3) =         0.0000000049 d0 !# 
      FSEfns( 75, -2, 3, 4) =        -0.0001074154 d0 !# Fermi    5.339 fm
      errtot( 75, -2, 3, 4) =         0.0000000504 d0 !# 
      FSEfns( 75, -2, 3, 5) =        -0.0001086764 d0 !# Fermi    5.339 fm
      errtot( 75, -2, 3, 5) =         0.0000003975 d0 !# 
      FSEfns( 75, -2, 4, 4) =        -0.0001111104 d0 !# Fermi    5.339 fm
      errtot( 75, -2, 4, 4) =         0.0000000385 d0 !# 
      FSEfns( 75, -2, 4, 5) =        -0.0001153890 d0 !# Fermi    5.339 fm
      errtot( 75, -2, 4, 5) =         0.0000011434 d0 !# 
      FSEfns( 75, -2, 5, 5) =        -0.0001169400 d0 !# Fermi    5.339 fm
      errtot( 75, -2, 5, 5) =         0.0000006128 d0 !# 
      FSEfns( 80, -2, 2, 2) =        -0.0001179414 d0 !# Fermi    5.463 fm
      errtot( 80, -2, 2, 2) =         0.0000013282 d0 !# 
      FSEfns( 80, -2, 2, 3) =        -0.0001427576 d0 !# Fermi    5.463 fm
      errtot( 80, -2, 2, 3) =         0.0000000984 d0 !# 
      FSEfns( 80, -2, 2, 4) =        -0.0001459271 d0 !# Fermi    5.463 fm
      errtot( 80, -2, 2, 4) =         0.0000000344 d0 !# 
      FSEfns( 80, -2, 2, 5) =        -0.0001465424 d0 !# Fermi    5.463 fm
      errtot( 80, -2, 2, 5) =         0.0000000944 d0 !# 
      FSEfns( 80, -2, 3, 3) =        -0.0001609586 d0 !# Fermi    5.463 fm
      errtot( 80, -2, 3, 3) =         0.0000000208 d0 !# 
      FSEfns( 80, -2, 3, 4) =        -0.0001731653 d0 !# Fermi    5.463 fm
      errtot( 80, -2, 3, 4) =         0.0000000439 d0 !# 
      FSEfns( 80, -2, 3, 5) =        -0.0001756446 d0 !# Fermi    5.463 fm
      errtot( 80, -2, 3, 5) =         0.0000002480 d0 !# 
      FSEfns( 80, -2, 4, 4) =        -0.0001794966 d0 !# Fermi    5.463 fm
      errtot( 80, -2, 4, 4) =         0.0000000387 d0 !# 
      FSEfns( 80, -2, 4, 5) =        -0.0001864410 d0 !# Fermi    5.463 fm
      errtot( 80, -2, 4, 5) =         0.0000006686 d0 !# 
      FSEfns( 80, -2, 5, 5) =        -0.0001891511 d0 !# Fermi    5.463 fm
      errtot( 80, -2, 5, 5) =         0.0000004453 d0 !# 
      FSEfns( 85, -2, 2, 2) =        -0.0001802039 d0 !# Fermi    5.539 fm
      errtot( 85, -2, 2, 2) =         0.0000032375 d0 !# 
      FSEfns( 85, -2, 2, 3) =        -0.0002210303 d0 !# Fermi    5.539 fm
      errtot( 85, -2, 2, 3) =         0.0000000395 d0 !# 
      FSEfns( 85, -2, 2, 4) =        -0.0002272783 d0 !# Fermi    5.539 fm
      errtot( 85, -2, 2, 4) =         0.0000000547 d0 !# 
      FSEfns( 85, -2, 2, 5) =        -0.0002287741 d0 !# Fermi    5.539 fm
      errtot( 85, -2, 2, 5) =         0.0000004566 d0 !# 
      FSEfns( 85, -2, 3, 3) =        -0.0002507821 d0 !# Fermi    5.539 fm
      errtot( 85, -2, 3, 3) =         0.0000000649 d0 !# 
      FSEfns( 85, -2, 3, 4) =        -0.0002704517 d0 !# Fermi    5.539 fm
      errtot( 85, -2, 3, 4) =         0.0000000183 d0 !# 
      FSEfns( 85, -2, 3, 5) =        -0.0002749658 d0 !# Fermi    5.539 fm
      errtot( 85, -2, 3, 5) =         0.0000012054 d0 !# 
      FSEfns( 85, -2, 4, 4) =        -0.0002809513 d0 !# Fermi    5.539 fm
      errtot( 85, -2, 4, 4) =         0.0000000313 d0 !# 
      FSEfns( 85, -2, 4, 5) =        -0.0002918884 d0 !# Fermi    5.539 fm
      errtot( 85, -2, 4, 5) =         0.0000030370 d0 !# 
      FSEfns( 85, -2, 5, 5) =        -0.0002964508 d0 !# Fermi    5.539 fm
      errtot( 85, -2, 5, 5) =         0.0000024249 d0 !# 
      FSEfns( 90, -2, 2, 2) =        -0.0002802014 d0 !# Fermi    5.710 fm
      errtot( 90, -2, 2, 2) =         0.0000123356 d0 !# 
      FSEfns( 90, -2, 2, 3) =        -0.0003457960 d0 !# Fermi    5.710 fm
      errtot( 90, -2, 2, 3) =         0.0000000133 d0 !# 
      FSEfns( 90, -2, 2, 4) =        -0.0003577965 d0 !# Fermi    5.710 fm
      errtot( 90, -2, 2, 4) =         0.0000001737 d0 !# 
      FSEfns( 90, -2, 2, 5) =        -0.0003609026 d0 !# Fermi    5.710 fm
      errtot( 90, -2, 2, 5) =         0.0000005205 d0 !# 
      FSEfns( 90, -2, 3, 3) =        -0.0003949395 d0 !# Fermi    5.710 fm
      errtot( 90, -2, 3, 3) =         0.0000001456 d0 !# 
      FSEfns( 90, -2, 3, 4) =        -0.0004272502 d0 !# Fermi    5.710 fm
      errtot( 90, -2, 3, 4) =         0.0000000123 d0 !# 
      FSEfns( 90, -2, 3, 5) =        -0.0004352042 d0 !# Fermi    5.710 fm
      errtot( 90, -2, 3, 5) =         0.0000013810 d0 !# 
      FSEfns( 90, -2, 4, 4) =        -0.0004448596 d0 !# Fermi    5.710 fm
      errtot( 90, -2, 4, 4) =         0.0000000021 d0 !# 
      FSEfns( 90, -2, 4, 5) =        -0.0004617223 d0 !# Fermi    5.710 fm
      errtot( 90, -2, 4, 5) =         0.0000034602 d0 !# 
      FSEfns( 90, -2, 5, 5) =        -0.0004695545 d0 !# Fermi    5.710 fm
      errtot( 90, -2, 5, 5) =         0.0000030992 d0 !# 
      FSEfns( 95, -2, 2, 2) =        -0.0004260804 d0 !# Fermi    5.905 fm
      errtot( 95, -2, 2, 2) =         0.0000232596 d0 !# 
      FSEfns( 95, -2, 2, 3) =        -0.0005349687 d0 !# Fermi    5.905 fm
      errtot( 95, -2, 2, 3) =         0.0000012317 d0 !# 
      FSEfns( 95, -2, 2, 4) =        -0.0005574350 d0 !# Fermi    5.905 fm
      errtot( 95, -2, 2, 4) =         0.0000000540 d0 !# 
      FSEfns( 95, -2, 2, 5) =        -0.0005639367 d0 !# Fermi    5.905 fm
      errtot( 95, -2, 2, 5) =         0.0000009001 d0 !# 
      FSEfns( 95, -2, 3, 3) =        -0.0006152096 d0 !# Fermi    5.905 fm
      errtot( 95, -2, 3, 3) =         0.0000024184 d0 !# 
      FSEfns( 95, -2, 3, 4) =        -0.0006682180 d0 !# Fermi    5.905 fm
      errtot( 95, -2, 3, 4) =         0.0000004456 d0 !# 
      FSEfns( 95, -2, 3, 5) =        -0.0006830948 d0 !# Fermi    5.905 fm
      errtot( 95, -2, 3, 5) =         0.0000024517 d0 !# 
      FSEfns( 95, -2, 4, 4) =        -0.0006974712 d0 !# Fermi    5.905 fm
      errtot( 95, -2, 4, 4) =         0.0000000032 d0 !# 
      FSEfns( 95, -2, 4, 5) =        -0.0007256855 d0 !# Fermi    5.905 fm
      errtot( 95, -2, 4, 5) =         0.0000062577 d0 !# 
      FSEfns( 95, -2, 5, 5) =        -0.0007386798 d0 !# Fermi    5.905 fm
      errtot( 95, -2, 5, 5) =         0.0000060652 d0 !# 
      FSEfns(100, -2, 2, 2) =        -0.0005967179 d0 !# Fermi    5.857 fm
      errtot(100, -2, 2, 2) =         0.0000388523 d0 !# 
      FSEfns(100, -2, 2, 3) =        -0.0007700984 d0 !# Fermi    5.857 fm
      errtot(100, -2, 2, 3) =         0.0000069047 d0 !# 
      FSEfns(100, -2, 2, 4) =        -0.0008083335 d0 !# Fermi    5.857 fm
      errtot(100, -2, 2, 4) =         0.0000000934 d0 !# 
      FSEfns(100, -2, 2, 5) =        -0.0008201205 d0 !# Fermi    5.857 fm
      errtot(100, -2, 2, 5) =         0.0000015131 d0 !# 
      FSEfns(100, -2, 3, 3) =        -0.0008922198 d0 !# Fermi    5.857 fm
      errtot(100, -2, 3, 3) =         0.0000137370 d0 !# 
      FSEfns(100, -2, 3, 4) =        -0.0009725804 d0 !# Fermi    5.857 fm
      errtot(100, -2, 3, 4) =         0.0000002547 d0 !# 
      FSEfns(100, -2, 3, 5) =        -0.0009972737 d0 !# Fermi    5.857 fm
      errtot(100, -2, 3, 5) =         0.0000043257 d0 !# 
      FSEfns(100, -2, 4, 4) =        -0.0010178243 d0 !# Fermi    5.857 fm
      errtot(100, -2, 4, 4) =         0.0000000108 d0 !# 
      FSEfns(100, -2, 4, 5) =        -0.0010591682 d0 !# Fermi    5.857 fm
      errtot(100, -2, 4, 5) =         0.0000112299 d0 !# 
      FSEfns(100, -2, 5, 5) =        -0.0010793846 d0 !# Fermi    5.857 fm
      errtot(100, -2, 5, 5) =         0.0000117771 d0 !# 
      FSEfns(105, -2, 2, 2) =        -0.0008529545 d0 !# Fermi    5.919 fm
      errtot(105, -2, 2, 2) =         0.0000645388 d0 !# 
      FSEfns(105, -2, 2, 3) =        -0.0011186768 d0 !# Fermi    5.919 fm
      errtot(105, -2, 2, 3) =         0.0000234403 d0 !# 
      FSEfns(105, -2, 2, 4) =        -0.0011847694 d0 !# Fermi    5.919 fm
      errtot(105, -2, 2, 4) =         0.0000000076 d0 !# 
      FSEfns(105, -2, 2, 5) =        -0.0012067664 d0 !# Fermi    5.919 fm
      errtot(105, -2, 2, 5) =         0.0000000726 d0 !# 
      FSEfns(105, -2, 3, 3) =        -0.0013057455 d0 !# Fermi    5.919 fm
      errtot(105, -2, 3, 3) =         0.0000405841 d0 !# 
      FSEfns(105, -2, 3, 4) =        -0.0014295692 d0 !# Fermi    5.919 fm
      errtot(105, -2, 3, 4) =         0.0000000454 d0 !# 
      FSEfns(105, -2, 3, 5) =        -0.0014726024 d0 !# Fermi    5.919 fm
      errtot(105, -2, 3, 5) =         0.0000002647 d0 !# 
      FSEfns(105, -2, 4, 4) =        -0.0015002175 d0 !# Fermi    5.919 fm
      errtot(105, -2, 4, 4) =         0.0000000298 d0 !# 
      FSEfns(105, -2, 4, 5) =        -0.0015644808 d0 !# Fermi    5.919 fm
      errtot(105, -2, 4, 5) =         0.0000005624 d0 !# 
      FSEfns(105, -2, 5, 5) =        -0.0015963036 d0 !# Fermi    5.919 fm
      errtot(105, -2, 5, 5) =         0.0000006757 d0 !# 
      FSEfns(110, -2, 2, 2) =        -0.0011791200 d0 !# Fermi    5.993 fm
      errtot(110, -2, 2, 2) =         0.0001109911 d0 !# 
      FSEfns(110, -2, 2, 3) =        -0.0015822778 d0 !# Fermi    5.993 fm
      errtot(110, -2, 2, 3) =         0.0000159903 d0 !# 
      FSEfns(110, -2, 2, 4) =        -0.0017008694 d0 !# Fermi    5.993 fm
      errtot(110, -2, 2, 4) =         0.0000000131 d0 !# 
      FSEfns(110, -2, 2, 5) =        -0.0017414424 d0 !# Fermi    5.993 fm
      errtot(110, -2, 2, 5) =         0.0000005503 d0 !# 
      FSEfns(110, -2, 3, 3) =        -0.0018577248 d0 !# Fermi    5.993 fm
      errtot(110, -2, 3, 3) =         0.0000482506 d0 !# 
      FSEfns(110, -2, 3, 4) =        -0.0020537203 d0 !# Fermi    5.993 fm
      errtot(110, -2, 3, 4) =         0.0000001206 d0 !# 
      FSEfns(110, -2, 3, 5) =        -0.0021282392 d0 !# Fermi    5.993 fm
      errtot(110, -2, 3, 5) =         0.0000012714 d0 !# 
      FSEfns(110, -2, 4, 4) =        -0.0021613984 d0 !# Fermi    5.993 fm
      errtot(110, -2, 4, 4) =         0.0000000690 d0 !# 
      FSEfns(110, -2, 4, 5) =        -0.0022597332 d0 !# Fermi    5.993 fm
      errtot(110, -2, 4, 5) =         0.0000039836 d0 !# 
      FSEfns(110, -2, 5, 5) =        -0.0023090774 d0 !# Fermi    5.993 fm
      errtot(110, -2, 5, 5) =         0.0000049095 d0 !# 
*
      FSEfns(115, -2, 2, 2) =        -0.0016 d0 ! to be replaced!
*
      FSEfns(115, -2, 2, 3) =        -0.0021541076 d0 !# Fermi    6.088 fm
      errtot(115, -2, 2, 3) =         0.0000213650 d0 !# 
*
      FSEfns(115, -2, 2, 4) =         -0.0023562978d0 ! to be replaced!
*
      FSEfns(115, -2, 2, 5) =        -0.0024303191 d0 !# Fermi    6.088 fm
      errtot(115, -2, 2, 5) =         0.0000001429 d0 !# 
      FSEfns(115, -2, 3, 3) =        -0.0025463145 d0 !# Fermi    6.088 fm
      errtot(115, -2, 3, 3) =         0.0000281289 d0 !# 
      FSEfns(115, -2, 3, 4) =        -0.0028329163 d0 !# Fermi    6.088 fm
      errtot(115, -2, 3, 4) =         0.0000002174 d0 !# 
      FSEfns(115, -2, 3, 5) =        -0.0029598288 d0 !# Fermi    6.088 fm
      errtot(115, -2, 3, 5) =         0.0000000341 d0 !# 
      FSEfns(115, -2, 4, 4) =        -0.0029901482 d0 !# Fermi    6.088 fm
      errtot(115, -2, 4, 4) =         0.0000000330 d0 !# 
      FSEfns(115, -2, 4, 5) =        -0.0031327820 d0 !# Fermi    6.088 fm
      errtot(115, -2, 4, 5) =         0.0000010758 d0 !# 
      FSEfns(115, -2, 5, 5) =        -0.0032061888 d0 !# Fermi    6.088 fm
      errtot(115, -2, 5, 5) =         0.0000012377 d0 !# 
      FSEfns(120, -2, 2, 2) =        -0.0018358470 d0 !# Fermi    6.175 fm
      errtot(120, -2, 2, 2) =         0.0000391651 d0 !# 
      FSEfns(120, -2, 2, 3) =        -0.0026749798 d0 !# Fermi    6.175 fm
      errtot(120, -2, 2, 3) =         0.0000444820 d0 !# 
      FSEfns(120, -2, 2, 4) =        -0.0030169655 d0 !# Fermi    6.175 fm
      errtot(120, -2, 2, 4) =         0.0000042415 d0 !# 
      FSEfns(120, -2, 2, 5) =        -0.0031509181 d0 !# Fermi    6.175 fm
      errtot(120, -2, 2, 5) =         0.0000000235 d0 !# 
      FSEfns(120, -2, 3, 3) =        -0.0031669163 d0 !# Fermi    6.175 fm
      errtot(120, -2, 3, 3) =         0.0001086365 d0 !# 
      FSEfns(120, -2, 3, 4) =        -0.0035663027 d0 !# Fermi    6.175 fm
      errtot(120, -2, 3, 4) =         0.0000084947 d0 !# 
      FSEfns(120, -2, 3, 5) =        -0.0037780591 d0 !# Fermi    6.175 fm
      errtot(120, -2, 3, 5) =         0.0000007637 d0 !# 
      FSEfns(120, -2, 4, 4) =        -0.0037749587 d0 !# Fermi    6.175 fm
      errtot(120, -2, 4, 4) =         0.0000135474 d0 !# 
      FSEfns(120, -2, 4, 5) =        -0.0039644485 d0 !# Fermi    6.175 fm
      errtot(120, -2, 4, 5) =         0.0000002098 d0 !# 
      FSEfns(120, -2, 5, 5) =        -0.0040652475 d0 !# Fermi    6.175 fm
      errtot(120, -2, 5, 5) =         0.0000000484 d0 !# 
*
* nd_3/2
*
      FSEfns(100,  2, 3, 3) =        -0.0000302133 d0 !# Fermi    5.857 fm
      errtot(100,  2, 3, 3) =         0.0000023792 d0 !# 
      FSEfns(100,  2, 3, 4) =        -0.0000402735 d0 !# Fermi    5.857 fm
      errtot(100,  2, 3, 4) =         0.0000001013 d0 !# 
      FSEfns(100,  2, 3, 5) =        -0.0000419748 d0 !# Fermi    5.857 fm
      errtot(100,  2, 3, 5) =         0.0000004043 d0 !# 
      FSEfns(100,  2, 4, 4) =        -0.0000475623 d0 !# Fermi    5.857 fm
      errtot(100,  2, 4, 4) =         0.0000000020 d0 !# 
      FSEfns(100,  2, 4, 5) =        -0.0000542426 d0 !# Fermi    5.857 fm
      errtot(100,  2, 4, 5) =         0.0000013179 d0 !# 
      FSEfns(100,  2, 5, 5) =        -0.0000574322 d0 !# Fermi    5.857 fm
      errtot(100,  2, 5, 5) =         0.0000018192 d0 !# 
      FSEfns(105,  2, 3, 3) =        -0.0000562106 d0 !# Fermi    5.919 fm
      errtot(105,  2, 3, 3) =         0.0000059643 d0 !# 
      FSEfns(105,  2, 3, 4) =        -0.0000761065 d0 !# Fermi    5.919 fm
      errtot(105,  2, 3, 4) =         0.0000000212 d0 !# 
      FSEfns(105,  2, 3, 5) =        -0.0000802645 d0 !# Fermi    5.919 fm
      errtot(105,  2, 3, 5) =         0.0000004761 d0 !# 
      FSEfns(105,  2, 4, 4) =        -0.0000898778 d0 !# Fermi    5.919 fm
      errtot(105,  2, 4, 4) =         0.0000000039 d0 !# 
      FSEfns(105,  2, 4, 5) =        -0.0001027565 d0 !# Fermi    5.919 fm
      errtot(105,  2, 4, 5) =         0.0000015436 d0 !# 
      FSEfns(105,  2, 5, 5) =        -0.0001088988 d0 !# Fermi    5.919 fm
      errtot(105,  2, 5, 5) =         0.0000022943 d0 !# 
      FSEfns(110,  2, 3, 3) =        -0.0001046033 d0 !# Fermi    5.993 fm
      errtot(110,  2, 3, 3) =         0.0000092689 d0 !# 
      FSEfns(110,  2, 3, 4) =        -0.0001439051 d0 !# Fermi    5.993 fm
      errtot(110,  2, 3, 4) =         0.0000002855 d0 !# 
      FSEfns(110,  2, 3, 5) =        -0.0001536611 d0 !# Fermi    5.993 fm
      errtot(110,  2, 3, 5) =         0.0000010164 d0 !# 
      FSEfns(110,  2, 4, 4) =        -0.0001697844 d0 !# Fermi    5.993 fm
      errtot(110,  2, 4, 4) =         0.0000000048 d0 !# 
      FSEfns(110,  2, 4, 5) =        -0.0001943270 d0 !# Fermi    5.993 fm
      errtot(110,  2, 4, 5) =         0.0000016602 d0 !# 
      FSEfns(110,  2, 5, 5) =        -0.0002059997 d0 !# Fermi    5.993 fm
      errtot(110,  2, 5, 5) =         0.0000026678 d0 !# 
      FSEfns(115,  2, 3, 3) =        -0.0002035607 d0 !# Fermi    6.088 fm
      errtot(115,  2, 3, 3) =         0.0000426949 d0 !# 
      FSEfns(115,  2, 3, 4) =        -0.0002719449 d0 !# Fermi    6.088 fm
      errtot(115,  2, 3, 4) =         0.0000003874 d0 !# 
      FSEfns(115,  2, 3, 5) =        -0.0002944958 d0 !# Fermi    6.088 fm
      errtot(115,  2, 3, 5) =         0.0000005159 d0 !# 
      FSEfns(115,  2, 4, 4) =        -0.0003200959 d0 !# Fermi    6.088 fm
      errtot(115,  2, 4, 4) =         0.0000001607 d0 !# 
      FSEfns(115,  2, 4, 5) =        -0.0003665518 d0 !# Fermi    6.088 fm
      errtot(115,  2, 4, 5) =         0.0000013203 d0 !# 
      FSEfns(115,  2, 5, 5) =        -0.0003885120 d0 !# Fermi    6.088 fm
      errtot(115,  2, 5, 5) =         0.0000020054 d0 !# 
      FSEfns(120,  2, 3, 3) =        -0.0003717512 d0 !# Fermi    6.175 fm
      errtot(120,  2, 3, 3) =         0.0001088767 d0 !# 
      FSEfns(120,  2, 3, 4) =        -0.0005058913 d0 !# Fermi    6.175 fm
      errtot(120,  2, 3, 4) =         0.0000093326 d0 !# 
      FSEfns(120,  2, 3, 5) =        -0.0005580444 d0 !# Fermi    6.175 fm
      errtot(120,  2, 3, 5) =         0.0000001531 d0 !# 
      FSEfns(120,  2, 4, 4) =        -0.0005925089 d0 !# Fermi    6.175 fm
      errtot(120,  2, 4, 4) =         0.0000176351 d0 !# 
      FSEfns(120,  2, 4, 5) =        -0.0006778614 d0 !# Fermi    6.175 fm
      errtot(120,  2, 4, 5) =         0.0000000825 d0 !# 
      FSEfns(120,  2, 5, 5) =        -0.0007176534 d0 !# Fermi    6.175 fm
      errtot(120,  2, 5, 5) =         0.0000000848 d0 !# 
*
*
      return
      end
c=======================================================================
c
c=======================================================================
      subroutine se_pnt_stor()
      implicit real*8 (a-h,o-z)

      real*8 FSE(120,-3:2,5,5),errtot(120,-3:2,5,5)
      common /re_stor/FSE,errtot

*
* ns
*
      FSE   ( 10, -1, 1, 1) =         4.6541616501 d0 ! #
      errtot( 10, -1, 1, 1) =         0.0000001912 d0 
      FSE   ( 10, -1, 1, 2) =         4.7960625271 d0 ! #
      errtot( 10, -1, 1, 2) =         0.0000005780 d0 
      FSE   ( 10, -1, 1, 3) =         4.8144711901 d0 ! #
      errtot( 10, -1, 1, 3) =         0.0000002123 d0 
      FSE   ( 10, -1, 1, 4) =         4.8192789034 d0 ! #
      errtot( 10, -1, 1, 4) =         0.0000015563 d0 
      FSE   ( 10, -1, 1, 5) =         4.8210462454 d0 ! #
      errtot( 10, -1, 1, 5) =         0.0000682876 d0 
      FSE   ( 10, -1, 2, 2) =         4.8944153180 d0 ! #
      errtot( 10, -1, 2, 2) =         0.0000059602 d0 
      FSE   ( 10, -1, 2, 3) =         4.9324970060 d0 ! #
      errtot( 10, -1, 2, 3) =         0.0000006169 d0 
      FSE   ( 10, -1, 2, 4) =         4.9437144559 d0 ! #
      errtot( 10, -1, 2, 4) =         0.0000010781 d0 
      FSE   ( 10, -1, 2, 5) =         4.9480354607 d0 ! #
      errtot( 10, -1, 2, 5) =         0.0000358082 d0 
      FSE   ( 10, -1, 3, 3) =         4.9524137201 d0 ! #
      errtot( 10, -1, 3, 3) =         0.0000005633 d0 
      FSE   ( 10, -1, 3, 4) =         4.9676882753 d0 ! #
      errtot( 10, -1, 3, 4) =         0.0000058940 d0 
      FSE   ( 10, -1, 3, 5) =         4.9739563386 d0 ! #
      errtot( 10, -1, 3, 5) =         0.0000058935 d0 
      FSE   ( 10, -1, 4, 4) =         4.9749198812 d0 ! #
      errtot( 10, -1, 4, 4) =         0.0000017730 d0 
      FSE   ( 10, -1, 4, 5) =         4.9824412050 d0 ! #
      errtot( 10, -1, 4, 5) =         0.0000132481 d0 
      FSE   ( 10, -1, 5, 5) =         4.9858002595 d0 ! #
      errtot( 10, -1, 5, 5) =         0.0000217368 d0 
      FSE   ( 15, -1, 1, 1) =         3.8014099586 d0 ! #
      errtot( 15, -1, 1, 1) =         0.0000005909 d0 
      FSE   ( 15, -1, 1, 2) =         3.9462555448 d0 ! #
      errtot( 15, -1, 1, 2) =         0.0000001240 d0 
      FSE   ( 15, -1, 1, 3) =         3.9639442216 d0 ! #
      errtot( 15, -1, 1, 3) =         0.0000001592 d0 
      FSE   ( 15, -1, 1, 4) =         3.9680725706 d0 ! #
      errtot( 15, -1, 1, 4) =         0.0000006148 d0 
      FSE   ( 15, -1, 1, 5) =         3.9693322058 d0 ! #
      errtot( 15, -1, 1, 5) =         0.0000293885 d0 
      FSE   ( 15, -1, 2, 2) =         4.0508754709 d0 ! #
      errtot( 15, -1, 2, 2) =         0.0000006640 d0 
      FSE   ( 15, -1, 2, 3) =         4.0885402456 d0 ! #
      errtot( 15, -1, 2, 3) =         0.0000000632 d0 
      FSE   ( 15, -1, 2, 4) =         4.0991852096 d0 ! #
      errtot( 15, -1, 2, 4) =         0.0000003758 d0 
      FSE   ( 15, -1, 2, 5) =         4.1030385036 d0 ! #
      errtot( 15, -1, 2, 5) =         0.0000168553 d0 
      FSE   ( 15, -1, 3, 3) =         4.1082320300 d0 ! #
      errtot( 15, -1, 3, 3) =         0.0000012845 d0 
      FSE   ( 15, -1, 3, 4) =         4.1229389963 d0 ! #
      errtot( 15, -1, 3, 4) =         0.0000023742 d0 
      FSE   ( 15, -1, 3, 5) =         4.1287529130 d0 ! #
      errtot( 15, -1, 3, 5) =         0.0000034940 d0 
      FSE   ( 15, -1, 4, 4) =         4.1296289809 d0 ! #
      errtot( 15, -1, 4, 4) =         0.0000007631 d0 
      FSE   ( 15, -1, 4, 5) =         4.1366932336 d0 ! #
      errtot( 15, -1, 4, 5) =         0.0000037949 d0 
      FSE   ( 15, -1, 5, 5) =         4.1396217404 d0 ! #
      errtot( 15, -1, 5, 5) =         0.0000026819 d0 
      FSE   ( 20, -1, 1, 1) =         3.2462544375 d0 ! #
      errtot( 20, -1, 1, 1) =         0.0000007592 d0 
      FSE   ( 20, -1, 1, 2) =         3.3945569178 d0 ! #
      errtot( 20, -1, 1, 2) =         0.0000003353 d0 
      FSE   ( 20, -1, 1, 3) =         3.4114354162 d0 ! #
      errtot( 20, -1, 1, 3) =         0.0000004557 d0 
      FSE   ( 20, -1, 1, 4) =         3.4147822141 d0 ! #
      errtot( 20, -1, 1, 4) =         0.0000008577 d0 
      FSE   ( 20, -1, 1, 5) =         3.4154581897 d0 ! #
      errtot( 20, -1, 1, 5) =         0.0000252995 d0 
      FSE   ( 20, -1, 2, 2) =         3.5066469903 d0 ! #
      errtot( 20, -1, 2, 2) =         0.0000005291 d0 
      FSE   ( 20, -1, 2, 3) =         3.5438492743 d0 ! #
      errtot( 20, -1, 2, 3) =         0.0000010351 d0 
      FSE   ( 20, -1, 2, 4) =         3.5538369584 d0 ! #
      errtot( 20, -1, 2, 4) =         0.0000003272 d0 
      FSE   ( 20, -1, 2, 5) =         3.5571508258 d0 ! #
      errtot( 20, -1, 2, 5) =         0.0000088432 d0 
      FSE   ( 20, -1, 3, 3) =         3.5633041687 d0 ! #
      errtot( 20, -1, 3, 3) =         0.0000014114 d0 
      FSE   ( 20, -1, 3, 4) =         3.5773380655 d0 ! #
      errtot( 20, -1, 3, 4) =         0.0000011616 d0 
      FSE   ( 20, -1, 3, 5) =         3.5826271291 d0 ! #
      errtot( 20, -1, 3, 5) =         0.0000013220 d0 
      FSE   ( 20, -1, 4, 4) =         3.5834053942 d0 ! #
      errtot( 20, -1, 4, 4) =         0.0000009558 d0 
      FSE   ( 20, -1, 4, 5) =         3.5899114212 d0 ! #
      errtot( 20, -1, 4, 5) =         0.0000063208 d0 
      FSE   ( 20, -1, 5, 5) =         3.5923101823 d0 ! #
      errtot( 20, -1, 5, 5) =         0.0000140388 d0 
      FSE   ( 25, -1, 1, 1) =         2.8501040428 d0 ! #
      errtot( 25, -1, 1, 1) =         0.0000006390 d0 
      FSE   ( 25, -1, 1, 2) =         3.0023383931 d0 ! #
      errtot( 25, -1, 1, 2) =         0.0000002023 d0 
      FSE   ( 25, -1, 1, 3) =         3.0183482345 d0 ! #
      errtot( 25, -1, 1, 3) =         0.0000006501 d0 
      FSE   ( 25, -1, 1, 4) =         3.0208336985 d0 ! #
      errtot( 25, -1, 1, 4) =         0.0000049806 d0 
      FSE   ( 25, -1, 1, 5) =         3.0208606704 d0 ! #
      errtot( 25, -1, 1, 5) =         0.0000248337 d0 
      FSE   ( 25, -1, 2, 2) =         3.1229584751 d0 ! #
      errtot( 25, -1, 2, 2) =         0.0000001361 d0 
      FSE   ( 25, -1, 2, 3) =         3.1596648335 d0 ! #
      errtot( 25, -1, 2, 3) =         0.0000003609 d0 
      FSE   ( 25, -1, 2, 4) =         3.1689283771 d0 ! #
      errtot( 25, -1, 2, 4) =         0.0000013050 d0 
      FSE   ( 25, -1, 2, 5) =         3.1716302411 d0 ! #
      errtot( 25, -1, 2, 5) =         0.0000059229 d0 
      FSE   ( 25, -1, 3, 3) =         3.1788759407 d0 ! #
      errtot( 25, -1, 3, 3) =         0.0000009879 d0 
      FSE   ( 25, -1, 3, 4) =         3.1921521035 d0 ! #
      errtot( 25, -1, 3, 4) =         0.0000007211 d0 
      FSE   ( 25, -1, 3, 5) =         3.1968544647 d0 ! #
      errtot( 25, -1, 3, 5) =         0.0000016859 d0 
      FSE   ( 25, -1, 4, 4) =         3.1975190641 d0 ! #
      errtot( 25, -1, 4, 4) =         0.0000014393 d0 
      FSE   ( 25, -1, 4, 5) =         3.2033918589 d0 ! #
      errtot( 25, -1, 4, 5) =         0.0000044279 d0 
      FSE   ( 25, -1, 5, 5) =         3.2051933149 d0 ! #
      errtot( 25, -1, 5, 5) =         0.0000217181 d0 
      FSE   ( 30, -1, 1, 1) =         2.5520150977 d0 ! #
      errtot( 30, -1, 1, 1) =         0.0000005460 d0 
      FSE   ( 30, -1, 1, 2) =         2.7086705915 d0 ! #
      errtot( 30, -1, 1, 2) =         0.0000007167 d0 
      FSE   ( 30, -1, 1, 3) =         2.7237721173 d0 ! #
      errtot( 30, -1, 1, 3) =         0.0000014791 d0 
      FSE   ( 30, -1, 1, 4) =         2.7253271034 d0 ! #
      errtot( 30, -1, 1, 4) =         0.0000136348 d0 
      FSE   ( 30, -1, 1, 5) =         2.7246493453 d0 ! #
      errtot( 30, -1, 1, 5) =         0.0000079518 d0 
      FSE   ( 30, -1, 2, 2) =         2.8388374559 d0 ! #
      errtot( 30, -1, 2, 2) =         0.0000007938 d0 
      FSE   ( 30, -1, 2, 3) =         2.8750230361 d0 ! #
      errtot( 30, -1, 2, 3) =         0.0000013067 d0 
      FSE   ( 30, -1, 2, 4) =         2.8834827856 d0 ! #
      errtot( 30, -1, 2, 4) =         0.0000012466 d0 
      FSE   ( 30, -1, 2, 5) =         2.8855154299 d0 ! #
      errtot( 30, -1, 2, 5) =         0.0000037311 d0 
      FSE   ( 30, -1, 3, 3) =         2.8939800779 d0 ! #
      errtot( 30, -1, 3, 3) =         0.0000006942 d0 
      FSE   ( 30, -1, 3, 4) =         2.9064132272 d0 ! #
      errtot( 30, -1, 3, 4) =         0.0000004561 d0 
      FSE   ( 30, -1, 3, 5) =         2.9104687717 d0 ! #
      errtot( 30, -1, 3, 5) =         0.0000016400 d0 
      FSE   ( 30, -1, 4, 4) =         2.9110012014 d0 ! #
      errtot( 30, -1, 4, 4) =         0.0000014455 d0 
      FSE   ( 30, -1, 4, 5) =         2.9161712036 d0 ! #
      errtot( 30, -1, 4, 5) =         0.0000015248 d0 
      FSE   ( 30, -1, 5, 5) =         2.9173068110 d0 ! #
      errtot( 30, -1, 5, 5) =         0.0000066504 d0 
      FSE   ( 35, -1, 1, 1) =         2.3199760963 d0 ! #
      errtot( 35, -1, 1, 1) =         0.0000002987 d0 
      FSE   ( 35, -1, 1, 2) =         2.4815928312 d0 ! #
      errtot( 35, -1, 1, 2) =         0.0000008490 d0 
      FSE   ( 35, -1, 1, 3) =         2.4957592654 d0 ! #
      errtot( 35, -1, 1, 3) =         0.0000007701 d0 
      FSE   ( 35, -1, 1, 4) =         2.4963125948 d0 ! #
      errtot( 35, -1, 1, 4) =         0.0000059436 d0 
      FSE   ( 35, -1, 1, 5) =         2.4948774588 d0 ! #
      errtot( 35, -1, 1, 5) =         0.0000008123 d0 
      FSE   ( 35, -1, 2, 2) =         2.6223351290 d0 ! #
      errtot( 35, -1, 2, 2) =         0.0000011966 d0 
      FSE   ( 35, -1, 2, 3) =         2.6579635799 d0 ! #
      errtot( 35, -1, 2, 3) =         0.0000014913 d0 
      FSE   ( 35, -1, 2, 4) =         2.6655590526 d0 ! #
      errtot( 35, -1, 2, 4) =         0.0000006147 d0 
      FSE   ( 35, -1, 2, 5) =         2.6668560399 d0 ! #
      errtot( 35, -1, 2, 5) =         0.0000029094 d0 
      FSE   ( 35, -1, 3, 3) =         2.6766590495 d0 ! #
      errtot( 35, -1, 3, 3) =         0.0000007264 d0 
      FSE   ( 35, -1, 3, 4) =         2.6881585310 d0 ! #
      errtot( 35, -1, 3, 4) =         0.0000002657 d0 
      FSE   ( 35, -1, 3, 5) =         2.6915039451 d0 ! #
      errtot( 35, -1, 3, 5) =         0.0000017337 d0 
      FSE   ( 35, -1, 4, 4) =         2.6918855015 d0 ! #
      errtot( 35, -1, 4, 4) =         0.0000017549 d0 
      FSE   ( 35, -1, 4, 5) =         2.6962727289 d0 ! #
      errtot( 35, -1, 4, 5) =         0.0000008460 d0 
      FSE   ( 35, -1, 5, 5) =         2.6966717704 d0 ! #
      errtot( 35, -1, 5, 5) =         0.0000046493 d0 
      FSE   ( 40, -1, 1, 1) =         2.1352284872 d0 ! #
      errtot( 40, -1, 1, 1) =         0.0000003270 d0 
      FSE   ( 40, -1, 1, 2) =         2.3024092805 d0 ! #
      errtot( 40, -1, 1, 2) =         0.0000005039 d0 
      FSE   ( 40, -1, 1, 3) =         2.3156170744 d0 ! #
      errtot( 40, -1, 1, 3) =         0.0000056186 d0 
      FSE   ( 40, -1, 1, 4) =         2.3151045289 d0 ! #
      errtot( 40, -1, 1, 4) =         0.0000065072 d0 
      FSE   ( 40, -1, 1, 5) =         2.3128525683 d0 ! #
      errtot( 40, -1, 1, 5) =         0.0000023201 d0 
      FSE   ( 40, -1, 2, 2) =         2.4548280564 d0 ! #
      errtot( 40, -1, 2, 2) =         0.0000012656 d0 
      FSE   ( 40, -1, 2, 3) =         2.4898570353 d0 ! #
      errtot( 40, -1, 2, 3) =         0.0000005512 d0 
      FSE   ( 40, -1, 2, 4) =         2.4965109103 d0 ! #
      errtot( 40, -1, 2, 4) =         0.0000005937 d0 
      FSE   ( 40, -1, 2, 5) =         2.4970012422 d0 ! #
      errtot( 40, -1, 2, 5) =         0.0000022906 d0 
      FSE   ( 40, -1, 3, 3) =         2.5082725945 d0 ! #
      errtot( 40, -1, 3, 3) =         0.0000009079 d0 
      FSE   ( 40, -1, 3, 4) =         2.5187348954 d0 ! #
      errtot( 40, -1, 3, 4) =         0.0000004838 d0 
      FSE   ( 40, -1, 3, 5) =         2.5212997491 d0 ! #
      errtot( 40, -1, 3, 5) =         0.0000011928 d0 
      FSE   ( 40, -1, 4, 4) =         2.5215067733 d0 ! #
      errtot( 40, -1, 4, 4) =         0.0000008800 d0 
      FSE   ( 40, -1, 4, 5) =         2.5250239025 d0 ! #
      errtot( 40, -1, 4, 5) =         0.0000008025 d0 
      FSE   ( 40, -1, 5, 5) =         2.5246055804 d0 ! #
      errtot( 40, -1, 5, 5) =         0.0000031349 d0 
      FSE   ( 45, -1, 1, 1) =         1.9859439071 d0 ! #
      errtot( 45, -1, 1, 1) =         0.0000002739 d0 
      FSE   ( 45, -1, 1, 2) =         2.1593837478 d0 ! #
      errtot( 45, -1, 1, 2) =         0.0000010759 d0 
      FSE   ( 45, -1, 1, 3) =         2.1716101043 d0 ! #
      errtot( 45, -1, 1, 3) =         0.0000060308 d0 
      FSE   ( 45, -1, 1, 4) =         2.1699505427 d0 ! #
      errtot( 45, -1, 1, 4) =         0.0000006984 d0 
      FSE   ( 45, -1, 1, 5) =         2.1668201449 d0 ! #
      errtot( 45, -1, 1, 5) =         0.0000011817 d0 
      FSE   ( 45, -1, 2, 2) =         2.3246889530 d0 ! #
      errtot( 45, -1, 2, 2) =         0.0000008584 d0 
      FSE   ( 45, -1, 2, 3) =         2.3590691278 d0 ! #
      errtot( 45, -1, 2, 3) =         0.0000013967 d0 
      FSE   ( 45, -1, 2, 4) =         2.3646839476 d0 ! #
      errtot( 45, -1, 2, 4) =         0.0000011536 d0 
      FSE   ( 45, -1, 2, 5) =         2.3642837636 d0 ! #
      errtot( 45, -1, 2, 5) =         0.0000020935 d0 
      FSE   ( 45, -1, 3, 3) =         2.3771661860 d0 ! #
      errtot( 45, -1, 3, 3) =         0.0000006923 d0 
      FSE   ( 45, -1, 3, 4) =         2.3864706180 d0 ! #
      errtot( 45, -1, 3, 4) =         0.0000004269 d0 
      FSE   ( 45, -1, 3, 5) =         2.3881672292 d0 ! #
      errtot( 45, -1, 3, 5) =         0.0000009391 d0 
      FSE   ( 45, -1, 4, 4) =         2.3881732709 d0 ! #
      errtot( 45, -1, 4, 4) =         0.0000008491 d0 
      FSE   ( 45, -1, 4, 5) =         2.3907183864 d0 ! #
      errtot( 45, -1, 4, 5) =         0.0000008335 d0 
      FSE   ( 45, -1, 5, 5) =         2.3893908800 d0 ! #
      errtot( 45, -1, 5, 5) =         0.0000022818 d0 
      FSE   ( 50, -1, 1, 1) =         1.8642745205 d0 ! #
      errtot( 50, -1, 1, 1) =         0.0000002887 d0 
      FSE   ( 50, -1, 1, 2) =         2.0447802539 d0 ! #
      errtot( 50, -1, 1, 2) =         0.0000014294 d0 
      FSE   ( 50, -1, 1, 3) =         2.0559984853 d0 ! #
      errtot( 50, -1, 1, 3) =         0.0000040291 d0 
      FSE   ( 50, -1, 1, 4) =         2.0531038085 d0 ! #
      errtot( 50, -1, 1, 4) =         0.0000023455 d0 
      FSE   ( 50, -1, 1, 5) =         2.0490186966 d0 ! #
      errtot( 50, -1, 1, 5) =         0.0000025161 d0 
      FSE   ( 50, -1, 2, 2) =         2.2243367444 d0 ! #
      errtot( 50, -1, 2, 2) =         0.0000003526 d0 
      FSE   ( 50, -1, 2, 3) =         2.2579941618 d0 ! #
      errtot( 50, -1, 2, 3) =         0.0000006797 d0 
      FSE   ( 50, -1, 2, 4) =         2.2624558314 d0 ! #
      errtot( 50, -1, 2, 4) =         0.0000004518 d0 
      FSE   ( 50, -1, 2, 5) =         2.2610663353 d0 ! #
      errtot( 50, -1, 2, 5) =         0.0000016283 d0 
      FSE   ( 50, -1, 3, 3) =         2.2757219816 d0 ! #
      errtot( 50, -1, 3, 3) =         0.0000005245 d0 
      FSE   ( 50, -1, 3, 4) =         2.2837187383 d0 ! #
      errtot( 50, -1, 3, 4) =         0.0000003948 d0 
      FSE   ( 50, -1, 3, 5) =         2.2844479785 d0 ! #
      errtot( 50, -1, 3, 5) =         0.0000007600 d0 
      FSE   ( 50, -1, 4, 4) =         2.2842199096 d0 ! #
      errtot( 50, -1, 4, 4) =         0.0000010261 d0 
      FSE   ( 50, -1, 4, 5) =         2.2856711833 d0 ! #
      errtot( 50, -1, 4, 5) =         0.0000009585 d0 
      FSE   ( 50, -1, 5, 5) =         2.2833188520 d0 ! #
      errtot( 50, -1, 5, 5) =         0.0000038029 d0 
      FSE   ( 55, -1, 1, 1) =         1.7648305921 d0 ! #
      errtot( 55, -1, 1, 1) =         0.0000003199 d0 
      FSE   ( 55, -1, 1, 2) =         1.9533491099 d0 ! #
      errtot( 55, -1, 1, 2) =         0.0000025293 d0 
      FSE   ( 55, -1, 1, 3) =         1.9635223677 d0 ! #
      errtot( 55, -1, 1, 3) =         0.0000008434 d0 
      FSE   ( 55, -1, 1, 4) =         1.9592850285 d0 ! #
      errtot( 55, -1, 1, 4) =         0.0000011246 d0 
      FSE   ( 55, -1, 1, 5) =         1.9541567767 d0 ! #
      errtot( 55, -1, 1, 5) =         0.0000020503 d0 
      FSE   ( 55, -1, 2, 2) =         2.1487259469 d0 ! #
      errtot( 55, -1, 2, 2) =         0.0000004373 d0 
      FSE   ( 55, -1, 2, 3) =         2.1815623455 d0 ! #
      errtot( 55, -1, 2, 3) =         0.0000005570 d0 
      FSE   ( 55, -1, 2, 4) =         2.1847292144 d0 ! #
      errtot( 55, -1, 2, 4) =         0.0000010059 d0 
      FSE   ( 55, -1, 2, 5) =         2.1822303260 d0 ! #
      errtot( 55, -1, 2, 5) =         0.0000016775 d0 
      FSE   ( 55, -1, 3, 3) =         2.1988437247 d0 ! #
      errtot( 55, -1, 3, 3) =         0.0000004031 d0 
      FSE   ( 55, -1, 3, 4) =         2.2053513981 d0 ! #
      errtot( 55, -1, 3, 4) =         0.0000008728 d0 
      FSE   ( 55, -1, 3, 5) =         2.2049901972 d0 ! #
      errtot( 55, -1, 3, 5) =         0.0000005746 d0 
      FSE   ( 55, -1, 4, 4) =         2.2044844007 d0 ! #
      errtot( 55, -1, 4, 4) =         0.0000009325 d0 
      FSE   ( 55, -1, 4, 5) =         2.2046992491 d0 ! #
      errtot( 55, -1, 4, 5) =         0.0000002705 d0 
      FSE   ( 55, -1, 5, 5) =         2.2011982143 d0 ! #
      errtot( 55, -1, 5, 5) =         0.0000031490 d0 
      FSE   ( 60, -1, 1, 1) =         1.6838360339 d0 ! #
      errtot( 60, -1, 1, 1) =         0.0000003833 d0 
      FSE   ( 60, -1, 1, 2) =         1.8814853270 d0 ! #
      errtot( 60, -1, 1, 2) =         0.0000013520 d0 
      FSE   ( 60, -1, 1, 3) =         1.8905603113 d0 ! #
      errtot( 60, -1, 1, 3) =         0.0000017218 d0 
      FSE   ( 60, -1, 1, 4) =         1.8848478833 d0 ! #
      errtot( 60, -1, 1, 4) =         0.0000002619 d0 
      FSE   ( 60, -1, 1, 5) =         1.8785713346 d0 ! #
      errtot( 60, -1, 1, 5) =         0.0000015871 d0 
      FSE   ( 60, -1, 2, 2) =         2.0945167604 d0 ! #
      errtot( 60, -1, 2, 2) =         0.0000016366 d0 
      FSE   ( 60, -1, 2, 3) =         2.1264006925 d0 ! #
      errtot( 60, -1, 2, 3) =         0.0000007606 d0 
      FSE   ( 60, -1, 2, 4) =         2.1280938050 d0 ! #
      errtot( 60, -1, 2, 4) =         0.0000005165 d0 
      FSE   ( 60, -1, 2, 5) =         2.1243373003 d0 ! #
      errtot( 60, -1, 2, 5) =         0.0000012352 d0 
      FSE   ( 60, -1, 3, 3) =         2.1431258045 d0 ! #
      errtot( 60, -1, 3, 3) =         0.0000008024 d0 
      FSE   ( 60, -1, 3, 4) =         2.1479196839 d0 ! #
      errtot( 60, -1, 3, 4) =         0.0000004274 d0 
      FSE   ( 60, -1, 3, 5) =         2.1463158562 d0 ! #
      errtot( 60, -1, 3, 5) =         0.0000002324 d0 
      FSE   ( 60, -1, 4, 4) =         2.1454791319 d0 ! #
      errtot( 60, -1, 4, 4) =         0.0000004771 d0 
      FSE   ( 60, -1, 4, 5) =         2.1442817241 d0 ! #
      errtot( 60, -1, 4, 5) =         0.0000002423 d0 
      FSE   ( 60, -1, 5, 5) =         2.1394733541 d0 ! #
      errtot( 60, -1, 5, 5) =         0.0000026837 d0 
      FSE   ( 65, -1, 1, 1) =         1.6186365349 d0 ! #
      errtot( 65, -1, 1, 1) =         0.0000015369 d0 
      FSE   ( 65, -1, 1, 2) =         1.8267466265 d0 ! #
      errtot( 65, -1, 1, 2) =         0.0000002503 d0 
      FSE   ( 65, -1, 1, 3) =         1.8346468951 d0 ! #
      errtot( 65, -1, 1, 3) =         0.0000007981 d0 
      FSE   ( 65, -1, 1, 4) =         1.8272958330 d0 ! #
      errtot( 65, -1, 1, 4) =         0.0000006643 d0 
      FSE   ( 65, -1, 1, 5) =         1.8197404100 d0 ! #
      errtot( 65, -1, 1, 5) =         0.0000005427 d0 
      FSE   ( 65, -1, 2, 2) =         2.0596099487 d0 ! #
      errtot( 65, -1, 2, 2) =         0.0000010862 d0 
      FSE   ( 65, -1, 2, 3) =         2.0903599411 d0 ! #
      errtot( 65, -1, 2, 3) =         0.0000006889 d0 
      FSE   ( 65, -1, 2, 4) =         2.0903508289 d0 ! #
      errtot( 65, -1, 2, 4) =         0.0000005178 d0 
      FSE   ( 65, -1, 2, 5) =         2.0851529020 d0 ! #
      errtot( 65, -1, 2, 5) =         0.0000010282 d0 
      FSE   ( 65, -1, 3, 3) =         2.1063791790 d0 ! #
      errtot( 65, -1, 3, 3) =         0.0000010888 d0 
      FSE   ( 65, -1, 3, 4) =         2.1091797153 d0 ! #
      errtot( 65, -1, 3, 4) =         0.0000003729 d0 
      FSE   ( 65, -1, 3, 5) =         2.1061438054 d0 ! #
      errtot( 65, -1, 3, 5) =         0.0000002559 d0 
      FSE   ( 65, -1, 4, 4) =         2.1049102412 d0 ! #
      errtot( 65, -1, 4, 4) =         0.0000002776 d0 
      FSE   ( 65, -1, 4, 5) =         2.1020854956 d0 ! #
      errtot( 65, -1, 4, 5) =         0.0000001382 d0 
      FSE   ( 65, -1, 5, 5) =         2.0957739778 d0 ! #
      errtot( 65, -1, 5, 5) =         0.0000031428 d0 
      FSE   ( 70, -1, 1, 1) =         1.5674077683 d0 ! #
      errtot( 70, -1, 1, 1) =         0.0000007665 d0 
      FSE   ( 70, -1, 1, 2) =         1.7875768251 d0 ! #
      errtot( 70, -1, 1, 2) =         0.0000003134 d0 
      FSE   ( 70, -1, 1, 3) =         1.7941919818 d0 ! #
      errtot( 70, -1, 1, 3) =         0.0000001228 d0 
      FSE   ( 70, -1, 1, 4) =         1.7849949906 d0 ! #
      errtot( 70, -1, 1, 4) =         0.0000002494 d0 
      FSE   ( 70, -1, 1, 5) =         1.7759984338 d0 ! #
      errtot( 70, -1, 1, 5) =         0.0000000904 d0 
      FSE   ( 70, -1, 2, 2) =         2.0428907653 d0 ! #
      errtot( 70, -1, 2, 2) =         0.0000006310 d0 
      FSE   ( 70, -1, 2, 3) =         2.0722632075 d0 ! #
      errtot( 70, -1, 2, 3) =         0.0000005256 d0 
      FSE   ( 70, -1, 2, 4) =         2.0702592323 d0 ! #
      errtot( 70, -1, 2, 4) =         0.0000005207 d0 
      FSE   ( 70, -1, 2, 5) =         2.0633889312 d0 ! #
      errtot( 70, -1, 2, 5) =         0.0000008580 d0 
      FSE   ( 70, -1, 3, 3) =         2.0873719121 d0 ! #
      errtot( 70, -1, 3, 3) =         0.0000007525 d0 
      FSE   ( 70, -1, 3, 4) =         2.0878272222 d0 ! #
      errtot( 70, -1, 3, 4) =         0.0000001979 d0 
      FSE   ( 70, -1, 3, 5) =         2.0831201705 d0 ! #
      errtot( 70, -1, 3, 5) =         0.0000008325 d0 
      FSE   ( 70, -1, 4, 4) =         2.0814083739 d0 ! #
      errtot( 70, -1, 4, 4) =         0.0000011112 d0 
      FSE   ( 70, -1, 4, 5) =         2.0766895951 d0 ! #
      errtot( 70, -1, 4, 5) =         0.0000007409 d0 
      FSE   ( 70, -1, 5, 5) =         2.0686333460 d0 ! #
      errtot( 70, -1, 5, 5) =         0.0000025486 d0 
      FSE   ( 75, -1, 1, 1) =         1.5289842235 d0 ! #
      errtot( 75, -1, 1, 1) =         0.0000011676 d0 
      FSE   ( 75, -1, 1, 2) =         1.7631528976 d0 ! #
      errtot( 75, -1, 1, 2) =         0.0000012896 d0 
      FSE   ( 75, -1, 1, 3) =         1.7683247857 d0 ! #
      errtot( 75, -1, 1, 3) =         0.0000001014 d0 
      FSE   ( 75, -1, 1, 4) =         1.7570171145 d0 ! #
      errtot( 75, -1, 1, 4) =         0.0000007612 d0 
      FSE   ( 75, -1, 1, 5) =         1.7463747794 d0 ! #
      errtot( 75, -1, 1, 5) =         0.0000011053 d0 
      FSE   ( 75, -1, 2, 2) =         2.0441145434 d0 ! #
      errtot( 75, -1, 2, 2) =         0.0000034286 d0 
      FSE   ( 75, -1, 2, 3) =         2.0717811241 d0 ! #
      errtot( 75, -1, 2, 3) =         0.0000002732 d0 
      FSE   ( 75, -1, 2, 4) =         2.0674026270 d0 ! #
      errtot( 75, -1, 2, 4) =         0.0000005958 d0 
      FSE   ( 75, -1, 2, 5) =         2.0585661712 d0 ! #
      errtot( 75, -1, 2, 5) =         0.0000007161 d0 
      FSE   ( 75, -1, 3, 3) =         2.0857024829 d0 ! #
      errtot( 75, -1, 3, 3) =         0.0000004452 d0 
      FSE   ( 75, -1, 3, 4) =         2.0833650859 d0 ! #
      errtot( 75, -1, 3, 4) =         0.0000001006 d0 
      FSE   ( 75, -1, 3, 5) =         2.0766826162 d0 ! #
      errtot( 75, -1, 3, 5) =         0.0000012075 d0 
      FSE   ( 75, -1, 4, 4) =         2.0743905006 d0 ! #
      errtot( 75, -1, 4, 4) =         0.0000008461 d0 
      FSE   ( 75, -1, 4, 5) =         2.0674449638 d0 ! #
      errtot( 75, -1, 4, 5) =         0.0000021889 d0 
      FSE   ( 75, -1, 5, 5) =         2.0573411038 d0 ! #
      errtot( 75, -1, 5, 5) =         0.0000027030 d0 
      FSE   ( 80, -1, 1, 1) =         1.5027776215 d0 ! #
      errtot( 80, -1, 1, 1) =         0.0000012280 d0 
      FSE   ( 80, -1, 1, 2) =         1.7533339380 d0 ! #
      errtot( 80, -1, 1, 2) =         0.0000004719 d0 
      FSE   ( 80, -1, 1, 3) =         1.7568347793 d0 ! #
      errtot( 80, -1, 1, 3) =         0.0000001814 d0 
      FSE   ( 80, -1, 1, 4) =         1.7430737646 d0 ! #
      errtot( 80, -1, 1, 4) =         0.0000001615 d0 
      FSE   ( 80, -1, 1, 5) =         1.7305240572 d0 ! #
      errtot( 80, -1, 1, 5) =         0.0000009220 d0 
      FSE   ( 80, -1, 2, 2) =         2.0639050024 d0 ! #
      errtot( 80, -1, 2, 2) =         0.0000054140 d0 
      FSE   ( 80, -1, 2, 3) =         2.0894223504 d0 ! #
      errtot( 80, -1, 2, 3) =         0.0000003084 d0 
      FSE   ( 80, -1, 2, 4) =         2.0821707929 d0 ! #
      errtot( 80, -1, 2, 4) =         0.0000003243 d0 
      FSE   ( 80, -1, 2, 5) =         2.0709895837 d0 ! #
      errtot( 80, -1, 2, 5) =         0.0000004391 d0 
      FSE   ( 80, -1, 3, 3) =         2.1017815683 d0 ! #
      errtot( 80, -1, 3, 3) =         0.0000001831 d0 
      FSE   ( 80, -1, 3, 4) =         2.0960767499 d0 ! #
      errtot( 80, -1, 3, 4) =         0.0000005114 d0 
      FSE   ( 80, -1, 3, 5) =         2.0870270072 d0 ! #
      errtot( 80, -1, 3, 5) =         0.0000023117 d0 
      FSE   ( 80, -1, 4, 4) =         2.0840250929 d0 ! #
      errtot( 80, -1, 4, 4) =         0.0000017351 d0 
      FSE   ( 80, -1, 4, 5) =         2.0744331033 d0 ! #
      errtot( 80, -1, 4, 5) =         0.0000038732 d0 
      FSE   ( 80, -1, 5, 5) =         2.0618992117 d0 ! #
      errtot( 80, -1, 5, 5) =         0.0000060514 d0 
      FSE   ( 85, -1, 1, 1) =         1.4887627951 d0 ! #
      errtot( 85, -1, 1, 1) =         0.0000004567 d0 
      FSE   ( 85, -1, 1, 2) =         1.7586920635 d0 ! #
      errtot( 85, -1, 1, 2) =         0.0000008825 d0 
      FSE   ( 85, -1, 1, 3) =         1.7601950857 d0 ! #
      errtot( 85, -1, 1, 3) =         0.0000002954 d0 
      FSE   ( 85, -1, 1, 4) =         1.7435294945 d0 ! #
      errtot( 85, -1, 1, 4) =         0.0000002409 d0 
      FSE   ( 85, -1, 1, 5) =         1.7287339420 d0 ! #
      errtot( 85, -1, 1, 5) =         0.0000010241 d0 
      FSE   ( 85, -1, 2, 2) =         2.1038766530 d0 ! #
      errtot( 85, -1, 2, 2) =         0.0000023680 d0 
      FSE   ( 85, -1, 2, 3) =         2.1266345768 d0 ! #
      errtot( 85, -1, 2, 3) =         0.0000012495 d0 
      FSE   ( 85, -1, 2, 4) =         2.1158453643 d0 ! #
      errtot( 85, -1, 2, 4) =         0.0000003997 d0 
      FSE   ( 85, -1, 2, 5) =         2.1018237365 d0 ! #
      errtot( 85, -1, 2, 5) =         0.0000004452 d0 
      FSE   ( 85, -1, 3, 3) =         2.1369215511 d0 ! #
      errtot( 85, -1, 3, 3) =         0.0000002779 d0 
      FSE   ( 85, -1, 3, 4) =         2.1271013981 d0 ! #
      errtot( 85, -1, 3, 4) =         0.0000010817 d0 
      FSE   ( 85, -1, 3, 5) =         2.1151723078 d0 ! #
      errtot( 85, -1, 3, 5) =         0.0000042337 d0 
      FSE   ( 85, -1, 4, 4) =         2.1112951549 d0 ! #
      errtot( 85, -1, 4, 4) =         0.0000031085 d0 
      FSE   ( 85, -1, 4, 5) =         2.0985205920 d0 ! #
      errtot( 85, -1, 4, 5) =         0.0000054367 d0 
      FSE   ( 85, -1, 5, 5) =         2.0830662886 d0 ! #
      errtot( 85, -1, 5, 5) =         0.0000080577 d0 
      FSE   ( 90, -1, 1, 1) =         1.4875420586 d0 ! #
      errtot( 90, -1, 1, 1) =         0.0000004300 d0 
      FSE   ( 90, -1, 1, 2) =         1.7806483878 d0 ! #
      errtot( 90, -1, 1, 2) =         0.0000010236 d0 
      FSE   ( 90, -1, 1, 3) =         1.7796814226 d0 ! #
      errtot( 90, -1, 1, 3) =         0.0000004753 d0 
      FSE   ( 90, -1, 1, 4) =         1.7595080436 d0 ! #
      errtot( 90, -1, 1, 4) =         0.0000001648 d0 
      FSE   ( 90, -1, 1, 5) =         1.7420199017 d0 ! #
      errtot( 90, -1, 1, 5) =         0.0000012466 d0 
      FSE   ( 90, -1, 2, 2) =         2.1668847406 d0 ! #
      errtot( 90, -1, 2, 2) =         0.0000105096 d0 
      FSE   ( 90, -1, 2, 3) =         2.1860456587 d0 ! #
      errtot( 90, -1, 2, 3) =         0.0000005600 d0 
      FSE   ( 90, -1, 2, 4) =         2.1708190214 d0 ! #
      errtot( 90, -1, 2, 4) =         0.0000011344 d0 
      FSE   ( 90, -1, 2, 5) =         2.1532965114 d0 ! #
      errtot( 90, -1, 2, 5) =         0.0000011021 d0 
      FSE   ( 90, -1, 3, 3) =         2.1935610669 d0 ! #
      errtot( 90, -1, 3, 3) =         0.0000007152 d0 
      FSE   ( 90, -1, 3, 4) =         2.1786364496 d0 ! #
      errtot( 90, -1, 3, 4) =         0.0000019933 d0 
      FSE   ( 90, -1, 3, 5) =         2.1631482345 d0 ! #
      errtot( 90, -1, 3, 5) =         0.0000044863 d0 
      FSE   ( 90, -1, 4, 4) =         2.1581849937 d0 ! #
      errtot( 90, -1, 4, 4) =         0.0000047549 d0 
      FSE   ( 90, -1, 4, 5) =         2.1415268947 d0 ! #
      errtot( 90, -1, 4, 5) =         0.0000101996 d0 
      FSE   ( 90, -1, 5, 5) =         2.1225129529 d0 ! #
      errtot( 90, -1, 5, 5) =         0.0000131343 d0 
      FSE   ( 95, -1, 1, 1) =         1.5005126055 d0 ! #
      errtot( 95, -1, 1, 1) =         0.0000003994 d0 
      FSE   ( 95, -1, 1, 2) =         1.8217585036 d0 ! #
      errtot( 95, -1, 1, 2) =         0.0000004259 d0 
      FSE   ( 95, -1, 1, 3) =         1.8176315255 d0 ! #
      errtot( 95, -1, 1, 3) =         0.0000004096 d0 
      FSE   ( 95, -1, 1, 4) =         1.7931254929 d0 ! #
      errtot( 95, -1, 1, 4) =         0.0000000844 d0 
      FSE   ( 95, -1, 1, 5) =         1.7723426843 d0 ! #
      errtot( 95, -1, 1, 5) =         0.0000019557 d0 
      FSE   ( 95, -1, 2, 2) =         2.2575389745 d0 ! #
      errtot( 95, -1, 2, 2) =         0.0000024094 d0 
      FSE   ( 95, -1, 2, 3) =         2.2719274874 d0 ! #
      errtot( 95, -1, 2, 3) =         0.0000020979 d0 
      FSE   ( 95, -1, 2, 4) =         2.2510170928 d0 ! #
      errtot( 95, -1, 2, 4) =         0.0000006596 d0 
      FSE   ( 95, -1, 2, 5) =         2.2290932886 d0 ! #
      errtot( 95, -1, 2, 5) =         0.0000015286 d0 
      FSE   ( 95, -1, 3, 3) =         2.2756972006 d0 ! #
      errtot( 95, -1, 3, 3) =         0.0000013323 d0 
      FSE   ( 95, -1, 3, 4) =         2.2543335077 d0 ! #
      errtot( 95, -1, 3, 4) =         0.0000033403 d0 
      FSE   ( 95, -1, 3, 5) =         2.2343634831 d0 ! #
      errtot( 95, -1, 3, 5) =         0.0000099693 d0 
      FSE   ( 95, -1, 4, 4) =         2.2280361304 d0 ! #
      errtot( 95, -1, 4, 4) =         0.0000070436 d0 
      FSE   ( 95, -1, 4, 5) =         2.2065668056 d0 ! #
      errtot( 95, -1, 4, 5) =         0.0000108643 d0 
      FSE   ( 95, -1, 5, 5) =         2.1831431153 d0 ! #
      errtot( 95, -1, 5, 5) =         0.0000121584 d0 
      FSE   (100, -1, 1, 1) =         1.5302003919 d0 ! #
      errtot(100, -1, 1, 1) =         0.0000008280 d0 
      FSE   (100, -1, 1, 2) =         1.8862532902 d0 ! #
      errtot(100, -1, 1, 2) =         0.0000007429 d0 
      FSE   (100, -1, 1, 3) =         1.8779384620 d0 ! #
      errtot(100, -1, 1, 3) =         0.0000003045 d0 
      FSE   (100, -1, 1, 4) =         1.8479417809 d0 ! #
      errtot(100, -1, 1, 4) =         0.0000003909 d0 
      FSE   (100, -1, 1, 5) =         1.8230307446 d0 ! #
      errtot(100, -1, 1, 5) =         0.0000031527 d0 
      FSE   (100, -1, 2, 2) =         2.3831206493 d0 ! #
      errtot(100, -1, 2, 2) =         0.0000072678 d0 
      FSE   (100, -1, 2, 3) =         2.3910450962 d0 ! #
      errtot(100, -1, 2, 3) =         0.0000021206 d0 
      FSE   (100, -1, 2, 4) =         2.3626782831 d0 ! #
      errtot(100, -1, 2, 4) =         0.0000003177 d0 
      FSE   (100, -1, 2, 5) =         2.3350910500 d0 ! #
      errtot(100, -1, 2, 5) =         0.0000020766 d0 
      FSE   (100, -1, 3, 3) =         2.3896821669 d0 ! #
      errtot(100, -1, 3, 3) =         0.0000021186 d0 
      FSE   (100, -1, 3, 4) =         2.3600318193 d0 ! #
      errtot(100, -1, 3, 4) =         0.0000049740 d0 
      FSE   (100, -1, 3, 5) =         2.3342966151 d0 ! #
      errtot(100, -1, 3, 5) =         0.0000098409 d0 
      FSE   (100, -1, 4, 4) =         2.3262306382 d0 ! #
      errtot(100, -1, 4, 4) =         0.0000102251 d0 
      FSE   (100, -1, 4, 5) =         2.2986847808 d0 ! #
      errtot(100, -1, 4, 5) =         0.0000212702 d0 
      FSE   (100, -1, 5, 5) =         2.2696895509 d0 ! #
      errtot(100, -1, 5, 5) =         0.0000219507 d0 
      FSE   (105, -1, 1, 1) =         1.5809115478 d0 ! #
      errtot(105, -1, 1, 1) =         0.0000001781 d0 
      FSE   (105, -1, 1, 2) =         1.9810700368 d0 ! #
      errtot(105, -1, 1, 2) =         0.0000009347 d0 
      FSE   (105, -1, 1, 3) =         1.9669958379 d0 ! #
      errtot(105, -1, 1, 3) =         0.0000003824 d0 
      FSE   (105, -1, 1, 4) =         1.9298263377 d0 ! #
      errtot(105, -1, 1, 4) =         0.0000017423 d0 
      FSE   (105, -1, 1, 5) =         1.8995927341 d0 ! #
      errtot(105, -1, 1, 5) =         0.0000062058 d0 
      FSE   (105, -1, 2, 2) =         2.5552975921 d0 ! #
      errtot(105, -1, 2, 2) =         0.0000108931 d0 
      FSE   (105, -1, 2, 3) =         2.5542549748 d0 ! #
      errtot(105, -1, 2, 3) =         0.0000009450 d0 
      FSE   (105, -1, 2, 4) =         2.5158264066 d0 ! #
      errtot(105, -1, 2, 4) =         0.0000006259 d0 
      FSE   (105, -1, 2, 5) =         2.4807443599 d0 ! #
      errtot(105, -1, 2, 5) =         0.0000039201 d0 
      FSE   (105, -1, 3, 3) =         2.5457268792 d0 ! #
      errtot(105, -1, 3, 3) =         0.0000026036 d0 
      FSE   (105, -1, 3, 4) =         2.5051462485 d0 ! #
      errtot(105, -1, 3, 4) =         0.0000056716 d0 
      FSE   (105, -1, 3, 5) =         2.4717994700 d0 ! #
      errtot(105, -1, 3, 5) =         0.0000164459 d0 
      FSE   (105, -1, 4, 4) =         2.4614719391 d0 ! #
      errtot(105, -1, 4, 4) =         0.0000139592 d0 
      FSE   (105, -1, 4, 5) =         2.4260678678 d0 ! #
      errtot(105, -1, 4, 5) =         0.0000252154 d0 
      FSE   (105, -1, 5, 5) =         2.3898587241 d0 ! #
      errtot(105, -1, 5, 5) =         0.0000223418 d0 
      FSE   (110, -1, 1, 1) =         1.6600624972 d0 ! #
      errtot(110, -1, 1, 1) =         0.0000024132 d0 
      FSE   (110, -1, 1, 2) =         2.1179340107 d0 ! #
      errtot(110, -1, 1, 2) =         0.0000006996 d0 
      FSE   (110, -1, 1, 3) =         2.0956039877 d0 ! #
      errtot(110, -1, 1, 3) =         0.0000022235 d0 
      FSE   (110, -1, 1, 4) =         2.0487054891 d0 ! #
      errtot(110, -1, 1, 4) =         0.0000057516 d0 
      FSE   (110, -1, 1, 5) =         2.0113588498 d0 ! #
      errtot(110, -1, 1, 5) =         0.0000136278 d0 
      FSE   (110, -1, 2, 2) =         2.7935929168 d0 ! #
      errtot(110, -1, 2, 2) =         0.0000068428 d0 
      FSE   (110, -1, 2, 3) =         2.7797252362 d0 ! #
      errtot(110, -1, 2, 3) =         0.0000026453 d0 
      FSE   (110, -1, 2, 4) =         2.7272304999 d0 ! #
      errtot(110, -1, 2, 4) =         0.0000029971 d0 
      FSE   (110, -1, 2, 5) =         2.6818710510 d0 ! #
      errtot(110, -1, 2, 5) =         0.0000025775 d0 
      FSE   (110, -1, 3, 3) =         2.7609230876 d0 ! #
      errtot(110, -1, 3, 3) =         0.0000011041 d0 
      FSE   (110, -1, 3, 4) =         2.7054581563 d0 ! #
      errtot(110, -1, 3, 4) =         0.0000045068 d0 
      FSE   (110, -1, 3, 5) =         2.6617214814 d0 ! #
      errtot(110, -1, 3, 5) =         0.0000209815 d0 
      FSE   (110, -1, 4, 4) =         2.6483654846 d0 ! #
      errtot(110, -1, 4, 4) =         0.0000231317 d0 
      FSE   (110, -1, 4, 5) =         2.6024810801 d0 ! #
      errtot(110, -1, 4, 5) =         0.0000364987 d0 
      FSE   (110, -1, 5, 5) =         2.5566292634 d0 ! #
      errtot(110, -1, 5, 5) =         0.0000269712 d0 
      FSE   (115, -1, 1, 1) =         1.7810909606 d0 ! #
      errtot(115, -1, 1, 1) =         0.0000007916 d0 
      FSE   (115, -1, 1, 2) =         2.3179621550 d0 ! #
      errtot(115, -1, 1, 2) =         0.0000052814 d0 
      FSE   (115, -1, 1, 3) =         2.2831867646 d0 ! #
      errtot(115, -1, 1, 3) =         0.0000092259 d0 
      FSE   (115, -1, 1, 4) =         2.2224270395 d0 ! #
      errtot(115, -1, 1, 4) =         0.0000174230 d0 
      FSE   (115, -1, 1, 5) =         2.1751111450 d0 ! #
      errtot(115, -1, 1, 5) =         0.0000321070 d0 
      FSE   (115, -1, 2, 2) =         3.1330543617 d0 ! #
      errtot(115, -1, 2, 2) =         0.0000036116 d0 
      FSE   (115, -1, 2, 3) =         3.1000552795 d0 ! #
      errtot(115, -1, 2, 3) =         0.0000138342 d0 
      FSE   (115, -1, 2, 4) =         3.0269519179 d0 ! #
      errtot(115, -1, 2, 4) =         0.0000151203 d0 
      FSE   (115, -1, 2, 5) =         2.9668150676 d0 ! #
      errtot(115, -1, 2, 5) =         0.0000095315 d0 
      FSE   (115, -1, 3, 3) =         3.0659228841 d0 ! #
      errtot(115, -1, 3, 3) =         0.0000091061 d0 
      FSE   (115, -1, 3, 4) =         2.9892812337 d0 ! #
      errtot(115, -1, 3, 4) =         0.0000050363 d0 
      FSE   (115, -1, 3, 5) =         2.9307159873 d0 ! #
      errtot(115, -1, 3, 5) =         0.0000141026 d0 
      FSE   (115, -1, 4, 4) =         2.9131317323 d0 ! #
      errtot(115, -1, 4, 4) =         0.0000103875 d0 
      FSE   (115, -1, 4, 5) =         2.8526514077 d0 ! #
      errtot(115, -1, 4, 5) =         0.0000446237 d0 
      FSE   (115, -1, 5, 5) =         2.7933222292 d0 ! #
      errtot(115, -1, 5, 5) =         0.0000316181 d0 
      FSE   (120, -1, 1, 1) =         1.9708509708 d0 ! #
      errtot(120, -1, 1, 1) =         0.0000105678 d0 
      FSE   (120, -1, 1, 2) =         2.6233456086 d0 ! #
      errtot(120, -1, 1, 2) =         0.0000239671 d0 
      FSE   (120, -1, 1, 3) =         2.5684888224 d0 ! #
      errtot(120, -1, 1, 3) =         0.0000360244 d0 
      FSE   (120, -1, 1, 4) =         2.4865606278 d0 ! #
      errtot(120, -1, 1, 4) =         0.0000536514 d0 
      FSE   (120, -1, 1, 5) =         2.4242886646 d0 ! #
      errtot(120, -1, 1, 5) =         0.0000816646 d0 
      FSE   (120, -1, 2, 2) =         3.6437480112 d0 ! #
      errtot(120, -1, 2, 2) =         0.0000388883 d0 
      FSE   (120, -1, 2, 3) =         3.5803676174 d0 ! #
      errtot(120, -1, 2, 3) =         0.0000567411 d0 
      FSE   (120, -1, 2, 4) =         3.4749711888 d0 ! #
      errtot(120, -1, 2, 4) =         0.0000624803 d0 
      FSE   (120, -1, 2, 5) =         3.3920911941 d0 ! #
      errtot(120, -1, 2, 5) =         0.0000543726 d0 
      FSE   (120, -1, 3, 3) =         3.5218836171 d0 ! #
      errtot(120, -1, 3, 3) =         0.0000551109 d0 
      FSE   (120, -1, 3, 4) =         3.4130980091 d0 ! #
      errtot(120, -1, 3, 4) =         0.0000484251 d0 
      FSE   (120, -1, 3, 5) =         3.3319491069 d0 ! #
      errtot(120, -1, 3, 5) =         0.0000192737 d0 
      FSE   (120, -1, 4, 4) =         3.3080182144 d0 ! #
      errtot(120, -1, 4, 4) =         0.0000380108 d0 
      FSE   (120, -1, 4, 5) =         3.2258749095 d0 ! #
      errtot(120, -1, 4, 5) =         0.0000047700 d0 
      FSE   (120, -1, 5, 5) =         3.1464394270 d0 ! #
      errtot(120, -1, 5, 5) =         0.0000415234 d0 
*
* np_1/2
*
      FSE   ( 10,  1, 2, 2) =        -0.1148458245 d0 ! #
      errtot( 10,  1, 2, 2) =         0.0000000885 d0 
      FSE   ( 10,  1, 2, 3) =        -0.0962280593 d0 ! #
      errtot( 10,  1, 2, 3) =         0.0000003697 d0 
      FSE   ( 10,  1, 2, 4) =        -0.0946095046 d0 ! #
      errtot( 10,  1, 2, 4) =         0.0000497708 d0 
      FSE   ( 10,  1, 2, 5) =        -0.0945936418 d0 ! #
      errtot( 10,  1, 2, 5) =         0.0000090288 d0 
      FSE   ( 10,  1, 3, 3) =        -0.1020426709 d0 ! #
      errtot( 10,  1, 3, 3) =         0.0000458951 d0 
      FSE   ( 10,  1, 3, 4) =        -0.0940795276 d0 ! #
      errtot( 10,  1, 3, 4) =         0.0000025885 d0 
      FSE   ( 10,  1, 3, 5) =        -0.0923252750 d0 ! #
      errtot( 10,  1, 3, 5) =         0.0000026149 d0 
      FSE   ( 10,  1, 4, 4) =        -0.0963648332 d0 ! #
      errtot( 10,  1, 4, 4) =         0.0001304792 d0 
      FSE   ( 10,  1, 4, 5) =        -0.0921072660 d0 ! #
      errtot( 10,  1, 4, 5) =         0.0000004403 d0 
      FSE   ( 10,  1, 5, 5) =        -0.0932360890 d0 ! #
      errtot( 10,  1, 5, 5) =         0.0000076741 d0 
      FSE   ( 15,  1, 2, 2) =        -0.1045475395 d0 ! #
      errtot( 15,  1, 2, 2) =         0.0000022414 d0 
      FSE   ( 15,  1, 2, 3) =        -0.0851346702 d0 ! #
      errtot( 15,  1, 2, 3) =         0.0000005565 d0 
      FSE   ( 15,  1, 2, 4) =        -0.0832540226 d0 ! #
      errtot( 15,  1, 2, 4) =         0.0000012201 d0 
      FSE   ( 15,  1, 2, 5) =        -0.0831089106 d0 ! #
      errtot( 15,  1, 2, 5) =         0.0000036383 d0 
      FSE   ( 15,  1, 3, 3) =        -0.0900510646 d0 ! #
      errtot( 15,  1, 3, 3) =         0.0000006544 d0 
      FSE   ( 15,  1, 3, 4) =        -0.0817698594 d0 ! #
      errtot( 15,  1, 3, 4) =         0.0000215948 d0 
      FSE   ( 15,  1, 3, 5) =        -0.0798659093 d0 ! #
      errtot( 15,  1, 3, 5) =         0.0000087056 d0 
      FSE   ( 15,  1, 4, 4) =        -0.0837328151 d0 ! #
      errtot( 15,  1, 4, 4) =         0.0000097039 d0 
      FSE   ( 15,  1, 4, 5) =        -0.0793524935 d0 ! #
      errtot( 15,  1, 4, 5) =         0.0000002639 d0 
      FSE   ( 15,  1, 5, 5) =        -0.0803491028 d0 ! #
      errtot( 15,  1, 5, 5) =         0.0000061152 d0 
      FSE   ( 20,  1, 2, 2) =        -0.0925166642 d0 ! #
      errtot( 20,  1, 2, 2) =         0.0000011219 d0 
      FSE   ( 20,  1, 2, 3) =        -0.0721457319 d0 ! #
      errtot( 20,  1, 2, 3) =         0.0000008742 d0 
      FSE   ( 20,  1, 2, 4) =        -0.0699330937 d0 ! #
      errtot( 20,  1, 2, 4) =         0.0000007903 d0 
      FSE   ( 20,  1, 2, 5) =        -0.0696403771 d0 ! #
      errtot( 20,  1, 2, 5) =         0.0000016452 d0 
      FSE   ( 20,  1, 3, 3) =        -0.0760140998 d0 ! #
      errtot( 20,  1, 3, 3) =         0.0000004413 d0 
      FSE   ( 20,  1, 3, 4) =        -0.0673682830 d0 ! #
      errtot( 20,  1, 3, 4) =         0.0000078853 d0 
      FSE   ( 20,  1, 3, 5) =        -0.0652656115 d0 ! #
      errtot( 20,  1, 3, 5) =         0.0000037151 d0 
      FSE   ( 20,  1, 4, 4) =        -0.0689762728 d0 ! #
      errtot( 20,  1, 4, 4) =         0.0000009135 d0 
      FSE   ( 20,  1, 4, 5) =        -0.0644259140 d0 ! #
      errtot( 20,  1, 4, 5) =         0.0000033944 d0 
      FSE   ( 20,  1, 5, 5) =        -0.0652620568 d0 ! #
      errtot( 20,  1, 5, 5) =         0.0000098380 d0 
      FSE   ( 25,  1, 2, 2) =        -0.0790647217 d0 ! #
      errtot( 25,  1, 2, 2) =         0.0000009127 d0 
      FSE   ( 25,  1, 2, 3) =        -0.0575867968 d0 ! #
      errtot( 25,  1, 2, 3) =         0.0000006797 d0 
      FSE   ( 25,  1, 2, 4) =        -0.0549905748 d0 ! #
      errtot( 25,  1, 2, 4) =         0.0000009247 d0 
      FSE   ( 25,  1, 2, 5) =        -0.0545296984 d0 ! #
      errtot( 25,  1, 2, 5) =         0.0000010689 d0 
      FSE   ( 25,  1, 3, 3) =        -0.0603006086 d0 ! #
      errtot( 25,  1, 3, 3) =         0.0000009016 d0 
      FSE   ( 25,  1, 3, 4) =        -0.0512411402 d0 ! #
      errtot( 25,  1, 3, 4) =         0.0000009762 d0 
      FSE   ( 25,  1, 3, 5) =        -0.0489176050 d0 ! #
      errtot( 25,  1, 3, 5) =         0.0000010711 d0 
      FSE   ( 25,  1, 4, 4) =        -0.0524648533 d0 ! #
      errtot( 25,  1, 4, 4) =         0.0000011303 d0 
      FSE   ( 25,  1, 4, 5) =        -0.0477295312 d0 ! #
      errtot( 25,  1, 4, 5) =         0.0000027484 d0 
      FSE   ( 25,  1, 5, 5) =        -0.0483818558 d0 ! #
      errtot( 25,  1, 5, 5) =         0.0000149363 d0 
      FSE   ( 30,  1, 2, 2) =        -0.0643291535 d0 ! #
      errtot( 30,  1, 2, 2) =         0.0000010023 d0 
      FSE   ( 30,  1, 2, 3) =        -0.0416066936 d0 ! #
      errtot( 30,  1, 2, 3) =         0.0000008246 d0 
      FSE   ( 30,  1, 2, 4) =        -0.0385788456 d0 ! #
      errtot( 30,  1, 2, 4) =         0.0000006399 d0 
      FSE   ( 30,  1, 2, 5) =        -0.0379331722 d0 ! #
      errtot( 30,  1, 2, 5) =         0.0000005742 d0 
      FSE   ( 30,  1, 3, 3) =        -0.0430774157 d0 ! #
      errtot( 30,  1, 3, 3) =         0.0000003144 d0 
      FSE   ( 30,  1, 3, 4) =        -0.0335639412 d0 ! #
      errtot( 30,  1, 3, 4) =         0.0000006195 d0 
      FSE   ( 30,  1, 3, 5) =        -0.0309982269 d0 ! #
      errtot( 30,  1, 3, 5) =         0.0000013065 d0 
      FSE   ( 30,  1, 4, 4) =        -0.0343843944 d0 ! #
      errtot( 30,  1, 4, 4) =         0.0000007662 d0 
      FSE   ( 30,  1, 4, 5) =        -0.0294518301 d0 ! #
      errtot( 30,  1, 4, 5) =         0.0000017486 d0 
      FSE   ( 30,  1, 5, 5) =        -0.0299209914 d0 ! #
      errtot( 30,  1, 5, 5) =         0.0000092775 d0 
      FSE   ( 35,  1, 2, 2) =        -0.0483413223 d0 ! #
      errtot( 35,  1, 2, 2) =         0.0000008678 d0 
      FSE   ( 35,  1, 2, 3) =        -0.0242369752 d0 ! #
      errtot( 35,  1, 2, 3) =         0.0000004412 d0 
      FSE   ( 35,  1, 2, 4) =        -0.0207321613 d0 ! #
      errtot( 35,  1, 2, 4) =         0.0000003800 d0 
      FSE   ( 35,  1, 2, 5) =        -0.0198880505 d0 ! #
      errtot( 35,  1, 2, 5) =         0.0000005398 d0 
      FSE   ( 35,  1, 3, 3) =        -0.0243881601 d0 ! #
      errtot( 35,  1, 3, 3) =         0.0000007961 d0 
      FSE   ( 35,  1, 3, 4) =        -0.0143896285 d0 ! #
      errtot( 35,  1, 3, 4) =         0.0000011645 d0 
      FSE   ( 35,  1, 3, 5) =        -0.0115622247 d0 ! #
      errtot( 35,  1, 3, 5) =         0.0000013799 d0 
      FSE   ( 35,  1, 4, 4) =        -0.0147948202 d0 ! #
      errtot( 35,  1, 4, 4) =         0.0000003426 d0 
      FSE   ( 35,  1, 4, 5) =        -0.0096607996 d0 ! #
      errtot( 35,  1, 4, 5) =         0.0000012083 d0 
      FSE   ( 35,  1, 5, 5) =        -0.0099482400 d0 ! #
      errtot( 35,  1, 5, 5) =         0.0000061966 d0 
      FSE   ( 40,  1, 2, 2) =        -0.0310493856 d0 ! #
      errtot( 40,  1, 2, 2) =         0.0000007508 d0 
      FSE   ( 40,  1, 2, 3) =        -0.0054223348 d0 ! #
      errtot( 40,  1, 2, 3) =         0.0000004984 d0 
      FSE   ( 40,  1, 2, 4) =        -0.0013967404 d0 ! #
      errtot( 40,  1, 2, 4) =         0.0000004230 d0 
      FSE   ( 40,  1, 2, 5) =        -0.0003438548 d0 ! #
      errtot( 40,  1, 2, 5) =         0.0000004626 d0 
      FSE   ( 40,  1, 3, 3) =        -0.0041845172 d0 ! #
      errtot( 40,  1, 3, 3) =         0.0000013216 d0 
      FSE   ( 40,  1, 3, 4) =         0.0063234510 d0 ! #
      errtot( 40,  1, 3, 4) =         0.0000009802 d0 
      FSE   ( 40,  1, 3, 5) =         0.0094281879 d0 ! #
      errtot( 40,  1, 3, 5) =         0.0000013422 d0 
      FSE   ( 40,  1, 4, 4) =         0.0063376529 d0 ! #
      errtot( 40,  1, 4, 4) =         0.0000004898 d0 
      FSE   ( 40,  1, 4, 5) =         0.0116722370 d0 ! #
      errtot( 40,  1, 4, 5) =         0.0000016814 d0 
      FSE   ( 40,  1, 5, 5) =         0.0115596485 d0 ! #
      errtot( 40,  1, 5, 5) =         0.0000047554 d0 
      FSE   ( 45,  1, 2, 2) =        -0.0123298655 d0 ! #
      errtot( 45,  1, 2, 2) =         0.0000007725 d0 
      FSE   ( 45,  1, 2, 3) =         0.0149683084 d0 ! #
      errtot( 45,  1, 2, 3) =         0.0000003140 d0 
      FSE   ( 45,  1, 2, 4) =         0.0195573147 d0 ! #
      errtot( 45,  1, 2, 4) =         0.0000002376 d0 
      FSE   ( 45,  1, 2, 5) =         0.0208258759 d0 ! #
      errtot( 45,  1, 2, 5) =         0.0000001660 d0 
      FSE   ( 45,  1, 3, 3) =         0.0176588768 d0 ! #
      errtot( 45,  1, 3, 3) =         0.0000012547 d0 
      FSE   ( 45,  1, 3, 4) =         0.0286981790 d0 ! #
      errtot( 45,  1, 3, 4) =         0.0000010804 d0 
      FSE   ( 45,  1, 3, 5) =         0.0320907224 d0 ! #
      errtot( 45,  1, 3, 5) =         0.0000008178 d0 
      FSE   ( 45,  1, 4, 4) =         0.0291278478 d0 ! #
      errtot( 45,  1, 4, 4) =         0.0000005430 d0 
      FSE   ( 45,  1, 4, 5) =         0.0346522785 d0 ! #
      errtot( 45,  1, 4, 5) =         0.0000013698 d0 
      FSE   ( 45,  1, 5, 5) =         0.0347031247 d0 ! #
      errtot( 45,  1, 5, 5) =         0.0000032706 d0 
      FSE   ( 50,  1, 2, 2) =         0.0080125804 d0 ! #
      errtot( 50,  1, 2, 2) =         0.0000007058 d0 
      FSE   ( 50,  1, 2, 3) =         0.0371420374 d0 ! #
      errtot( 50,  1, 2, 3) =         0.0000002151 d0 
      FSE   ( 50,  1, 2, 4) =         0.0423358740 d0 ! #
      errtot( 50,  1, 2, 4) =         0.0000003022 d0 
      FSE   ( 50,  1, 2, 5) =         0.0438228258 d0 ! #
      errtot( 50,  1, 2, 5) =         0.0000002475 d0 
      FSE   ( 50,  1, 3, 3) =         0.0413566916 d0 ! #
      errtot( 50,  1, 3, 3) =         0.0000010378 d0 
      FSE   ( 50,  1, 3, 4) =         0.0529356196 d0 ! #
      errtot( 50,  1, 3, 4) =         0.0000009711 d0 
      FSE   ( 50,  1, 3, 5) =         0.0566225034 d0 ! #
      errtot( 50,  1, 3, 5) =         0.0000007801 d0 
      FSE   ( 50,  1, 4, 4) =         0.0537692089 d0 ! #
      errtot( 50,  1, 4, 4) =         0.0000005915 d0 
      FSE   ( 50,  1, 4, 5) =         0.0594655897 d0 ! #
      errtot( 50,  1, 4, 5) =         0.0000013150 d0 
      FSE   ( 50,  1, 5, 5) =         0.0596595276 d0 ! #
      errtot( 50,  1, 5, 5) =         0.0000027320 d0 
      FSE   ( 55,  1, 2, 2) =         0.0302532197 d0 ! #
      errtot( 55,  1, 2, 2) =         0.0000007093 d0 
      FSE   ( 55,  1, 2, 3) =         0.0613899047 d0 ! #
      errtot( 55,  1, 2, 3) =         0.0000001652 d0 
      FSE   ( 55,  1, 2, 4) =         0.0672281575 d0 ! #
      errtot( 55,  1, 2, 4) =         0.0000004476 d0 
      FSE   ( 55,  1, 2, 5) =         0.0689312325 d0 ! #
      errtot( 55,  1, 2, 5) =         0.0000005831 d0 
      FSE   ( 55,  1, 3, 3) =         0.0672027494 d0 ! #
      errtot( 55,  1, 3, 3) =         0.0000010417 d0 
      FSE   ( 55,  1, 3, 4) =         0.0793233396 d0 ! #
      errtot( 55,  1, 3, 4) =         0.0000006702 d0 
      FSE   ( 55,  1, 3, 5) =         0.0833060563 d0 ! #
      errtot( 55,  1, 3, 5) =         0.0000014498 d0 
      FSE   ( 55,  1, 4, 4) =         0.0805394320 d0 ! #
      errtot( 55,  1, 4, 4) =         0.0000007779 d0 
      FSE   ( 55,  1, 4, 5) =         0.0863800646 d0 ! #
      errtot( 55,  1, 4, 5) =         0.0000011244 d0 
      FSE   ( 55,  1, 5, 5) =         0.0866596520 d0 ! #
      errtot( 55,  1, 5, 5) =         0.0000121012 d0 
      FSE   ( 60,  1, 2, 2) =         0.0547635950 d0 ! #
      errtot( 60,  1, 2, 2) =         0.0000007458 d0 
      FSE   ( 60,  1, 2, 3) =         0.0881036496 d0 ! #
      errtot( 60,  1, 2, 3) =         0.0000008912 d0 
      FSE   ( 60,  1, 2, 4) =         0.0946232882 d0 ! #
      errtot( 60,  1, 2, 4) =         0.0000005283 d0 
      FSE   ( 60,  1, 2, 5) =         0.0965325981 d0 ! #
      errtot( 60,  1, 2, 5) =         0.0000005652 d0 
      FSE   ( 60,  1, 3, 3) =         0.0955973180 d0 ! #
      errtot( 60,  1, 3, 3) =         0.0000010111 d0 
      FSE   ( 60,  1, 3, 4) =         0.1082523692 d0 ! #
      errtot( 60,  1, 3, 4) =         0.0000010916 d0 
      FSE   ( 60,  1, 3, 5) =         0.1125198275 d0 ! #
      errtot( 60,  1, 3, 5) =         0.0000009113 d0 
      FSE   ( 60,  1, 4, 4) =         0.1098166237 d0 ! #
      errtot( 60,  1, 4, 4) =         0.0000008624 d0 
      FSE   ( 60,  1, 4, 5) =         0.1157613657 d0 ! #
      errtot( 60,  1, 4, 5) =         0.0000010451 d0 
      FSE   ( 60,  1, 5, 5) =         0.1161158095 d0 ! #
      errtot( 60,  1, 5, 5) =         0.0000112917 d0 
      FSE   ( 65,  1, 2, 2) =         0.0820377333 d0 ! #
      errtot( 65,  1, 2, 2) =         0.0000005604 d0 
      FSE   ( 65,  1, 2, 3) =         0.1178030189 d0 ! #
      errtot( 65,  1, 2, 3) =         0.0000004673 d0 
      FSE   ( 65,  1, 2, 4) =         0.1250350411 d0 ! #
      errtot( 65,  1, 2, 4) =         0.0000006132 d0 
      FSE   ( 65,  1, 2, 5) =         0.1271314329 d0 ! #
      errtot( 65,  1, 2, 5) =         0.0000003069 d0 
      FSE   ( 65,  1, 3, 3) =         0.1270724372 d0 ! #
      errtot( 65,  1, 3, 3) =         0.0000009536 d0 
      FSE   ( 65,  1, 3, 4) =         0.1402384333 d0 ! #
      errtot( 65,  1, 3, 4) =         0.0000010000 d0 
      FSE   ( 65,  1, 3, 5) =         0.1447728420 d0 ! #
      errtot( 65,  1, 3, 5) =         0.0000007959 d0 
      FSE   ( 65,  1, 4, 4) =         0.1421039688 d0 ! #
      errtot( 65,  1, 4, 4) =         0.0000010204 d0 
      FSE   ( 65,  1, 4, 5) =         0.1480902126 d0 ! #
      errtot( 65,  1, 4, 5) =         0.0000036257 d0 
      FSE   ( 65,  1, 5, 5) =         0.1484711433 d0 ! #
      errtot( 65,  1, 5, 5) =         0.0000109822 d0 
      FSE   ( 70,  1, 2, 2) =         0.1127326348 d0 ! #
      errtot( 70,  1, 2, 2) =         0.0000006249 d0 
      FSE   ( 70,  1, 2, 3) =         0.1511762488 d0 ! #
      errtot( 70,  1, 2, 3) =         0.0000006463 d0 
      FSE   ( 70,  1, 2, 4) =         0.1591436156 d0 ! #
      errtot( 70,  1, 2, 4) =         0.0000010007 d0 
      FSE   ( 70,  1, 2, 5) =         0.1613933834 d0 ! #
      errtot( 70,  1, 2, 5) =         0.0000002455 d0 
      FSE   ( 70,  1, 3, 3) =         0.1623335497 d0 ! #
      errtot( 70,  1, 3, 3) =         0.0000010647 d0 
      FSE   ( 70,  1, 3, 4) =         0.1759682431 d0 ! #
      errtot( 70,  1, 3, 4) =         0.0000010570 d0 
      FSE   ( 70,  1, 3, 5) =         0.1807344969 d0 ! #
      errtot( 70,  1, 3, 5) =         0.0000005999 d0 
      FSE   ( 70,  1, 4, 4) =         0.1780673094 d0 ! #
      errtot( 70,  1, 4, 4) =         0.0000011563 d0 
      FSE   ( 70,  1, 4, 5) =         0.1840242688 d0 ! #
      errtot( 70,  1, 4, 5) =         0.0000032523 d0 
      FSE   ( 70,  1, 5, 5) =         0.1843508630 d0 ! #
      errtot( 70,  1, 5, 5) =         0.0000102767 d0 
      FSE   ( 75,  1, 2, 2) =         0.1477297082 d0 ! #
      errtot( 75,  1, 2, 2) =         0.0000008940 d0 
      FSE   ( 75,  1, 2, 3) =         0.1891427513 d0 ! #
      errtot( 75,  1, 2, 3) =         0.0000007266 d0 
      FSE   ( 75,  1, 2, 4) =         0.1978526252 d0 ! #
      errtot( 75,  1, 2, 4) =         0.0000004447 d0 
      FSE   ( 75,  1, 2, 5) =         0.2002036549 d0 ! #
      errtot( 75,  1, 2, 5) =         0.0000012477 d0 
      FSE   ( 75,  1, 3, 3) =         0.2023214811 d0 ! #
      errtot( 75,  1, 3, 3) =         0.0000009779 d0 
      FSE   ( 75,  1, 3, 4) =         0.2163548612 d0 ! #
      errtot( 75,  1, 3, 4) =         0.0000009484 d0 
      FSE   ( 75,  1, 3, 5) =         0.2212957182 d0 ! #
      errtot( 75,  1, 3, 5) =         0.0000006792 d0 
      FSE   ( 75,  1, 4, 4) =         0.2185927622 d0 ! #
      errtot( 75,  1, 4, 4) =         0.0000011407 d0 
      FSE   ( 75,  1, 4, 5) =         0.2244133283 d0 ! #
      errtot( 75,  1, 4, 5) =         0.0000028048 d0 
      FSE   ( 75,  1, 5, 5) =         0.2245869371 d0 ! #
      errtot( 75,  1, 5, 5) =         0.0000094193 d0 
      FSE   ( 80,  1, 2, 2) =         0.1882259703 d0 ! #
      errtot( 80,  1, 2, 2) =         0.0000028946 d0 
      FSE   ( 80,  1, 2, 3) =         0.2329475983 d0 ! #
      errtot( 80,  1, 2, 3) =         0.0000007514 d0 
      FSE   ( 80,  1, 2, 4) =         0.2423831586 d0 ! #
      errtot( 80,  1, 2, 4) =         0.0000000447 d0 
      FSE   ( 80,  1, 2, 5) =         0.2447529662 d0 ! #
      errtot( 80,  1, 2, 5) =         0.0000009297 d0 
      FSE   ( 80,  1, 3, 3) =         0.2483089413 d0 ! #
      errtot( 80,  1, 3, 3) =         0.0000008835 d0 
      FSE   ( 80,  1, 3, 4) =         0.2626304613 d0 ! #
      errtot( 80,  1, 3, 4) =         0.0000007453 d0 
      FSE   ( 80,  1, 3, 5) =         0.2676557055 d0 ! #
      errtot( 80,  1, 3, 5) =         0.0000005625 d0 
      FSE   ( 80,  1, 4, 4) =         0.2648737415 d0 ! #
      errtot( 80,  1, 4, 4) =         0.0000010902 d0 
      FSE   ( 80,  1, 4, 5) =         0.2704103810 d0 ! #
      errtot( 80,  1, 4, 5) =         0.0000024437 d0 
      FSE   ( 80,  1, 5, 5) =         0.2702975127 d0 ! #
      errtot( 80,  1, 5, 5) =         0.0000088725 d0 
      FSE   ( 85,  1, 2, 2) =         0.2358861272 d0 ! #
      errtot( 85,  1, 2, 2) =         0.0000006329 d0 
      FSE   ( 85,  1, 2, 3) =         0.2843101933 d0 ! #
      errtot( 85,  1, 2, 3) =         0.0000009836 d0 
      FSE   ( 85,  1, 2, 4) =         0.2944134952 d0 ! #
      errtot( 85,  1, 2, 4) =         0.0000004175 d0 
      FSE   ( 85,  1, 2, 5) =         0.2966765211 d0 ! #
      errtot( 85,  1, 2, 5) =         0.0000006754 d0 
      FSE   ( 85,  1, 3, 3) =         0.3020485026 d0 ! #
      errtot( 85,  1, 3, 3) =         0.0000009619 d0 
      FSE   ( 85,  1, 3, 4) =         0.3164875542 d0 ! #
      errtot( 85,  1, 3, 4) =         0.0000006738 d0 
      FSE   ( 85,  1, 3, 5) =         0.3214572218 d0 ! #
      errtot( 85,  1, 3, 5) =         0.0000001037 d0 
      FSE   ( 85,  1, 4, 4) =         0.3185457530 d0 ! #
      errtot( 85,  1, 4, 4) =         0.0000010332 d0 
      FSE   ( 85,  1, 4, 5) =         0.3235937953 d0 ! #
      errtot( 85,  1, 4, 5) =         0.0000021941 d0 
      FSE   ( 85,  1, 5, 5) =         0.3230102379 d0 ! #
      errtot( 85,  1, 5, 5) =         0.0000094843 d0 
      FSE   ( 90,  1, 2, 2) =         0.2930727872 d0 ! #
      errtot( 90,  1, 2, 2) =         0.0000029047 d0 
      FSE   ( 90,  1, 2, 3) =         0.3456610655 d0 ! #
      errtot( 90,  1, 2, 3) =         0.0000002823 d0 
      FSE   ( 90,  1, 2, 4) =         0.3563069170 d0 ! #
      errtot( 90,  1, 2, 4) =         0.0000007688 d0 
      FSE   ( 90,  1, 2, 5) =         0.3582706392 d0 ! #
      errtot( 90,  1, 2, 5) =         0.0000004091 d0 
      FSE   ( 90,  1, 3, 3) =         0.3660106922 d0 ! #
      errtot( 90,  1, 3, 3) =         0.0000008552 d0 
      FSE   ( 90,  1, 3, 4) =         0.3803048337 d0 ! #
      errtot( 90,  1, 3, 4) =         0.0000005425 d0 
      FSE   ( 90,  1, 3, 5) =         0.3850041057 d0 ! #
      errtot( 90,  1, 3, 5) =         0.0000008809 d0 
      FSE   ( 90,  1, 4, 4) =         0.3819021082 d0 ! #
      errtot( 90,  1, 4, 4) =         0.0000010027 d0 
      FSE   ( 90,  1, 4, 5) =         0.3861726454 d0 ! #
      errtot( 90,  1, 4, 5) =         0.0000020297 d0 
      FSE   ( 90,  1, 5, 5) =         0.3848610924 d0 ! #
      errtot( 90,  1, 5, 5) =         0.0000095417 d0 
      FSE   ( 95,  1, 2, 2) =         0.3632535536 d0 ! #
      errtot( 95,  1, 2, 2) =         0.0000006534 d0 
      FSE   ( 95,  1, 2, 3) =         0.4205449474 d0 ! #
      errtot( 95,  1, 2, 3) =         0.0000005027 d0 
      FSE   ( 95,  1, 2, 4) =         0.4314927703 d0 ! #
      errtot( 95,  1, 2, 4) =         0.0000003331 d0 
      FSE   ( 95,  1, 2, 5) =         0.4328582529 d0 ! #
      errtot( 95,  1, 2, 5) =         0.0000007838 d0 
      FSE   ( 95,  1, 3, 3) =         0.4437830839 d0 ! #
      errtot( 95,  1, 3, 3) =         0.0000007326 d0 
      FSE   ( 95,  1, 3, 4) =         0.4575249396 d0 ! #
      errtot( 95,  1, 3, 4) =         0.0000003650 d0 
      FSE   ( 95,  1, 3, 5) =         0.4616201283 d0 ! #
      errtot( 95,  1, 3, 5) =         0.0000012794 d0 
      FSE   ( 95,  1, 4, 4) =         0.4582515669 d0 ! #
      errtot( 95,  1, 4, 4) =         0.0000010512 d0 
      FSE   ( 95,  1, 4, 5) =         0.4613299507 d0 ! #
      errtot( 95,  1, 4, 5) =         0.0000017494 d0 
      FSE   ( 95,  1, 5, 5) =         0.4589180147 d0 ! #
      errtot( 95,  1, 5, 5) =         0.0000100868 d0 
      FSE   (100,  1, 2, 2) =         0.4517110710 d0 ! #
      errtot(100,  1, 2, 2) =         0.0000005386 d0 
      FSE   (100,  1, 2, 3) =         0.5143245938 d0 ! #
      errtot(100,  1, 2, 3) =         0.0000006610 d0 
      FSE   (100,  1, 2, 4) =         0.5251348746 d0 ! #
      errtot(100,  1, 2, 4) =         0.0000010737 d0 
      FSE   (100,  1, 2, 5) =         0.5254289934 d0 ! #
      errtot(100,  1, 2, 5) =         0.0000008064 d0 
      FSE   (100,  1, 3, 3) =         0.5407683544 d0 ! #
      errtot(100,  1, 3, 3) =         0.0000005327 d0 
      FSE   (100,  1, 3, 4) =         0.5533142661 d0 ! #
      errtot(100,  1, 3, 4) =         0.0000001757 d0 
      FSE   (100,  1, 3, 5) =         0.5562805569 d0 ! #
      errtot(100,  1, 3, 5) =         0.0000009773 d0 
      FSE   (100,  1, 4, 4) =         0.5525435008 d0 ! #
      errtot(100,  1, 4, 4) =         0.0000010258 d0 
      FSE   (100,  1, 4, 5) =         0.5538170049 d0 ! #
      errtot(100,  1, 4, 5) =         0.0000033976 d0 
      FSE   (100,  1, 5, 5) =         0.5497578146 d0 ! #
      errtot(100,  1, 5, 5) =         0.0000107705 d0 
      FSE   (105,  1, 2, 2) =         0.5668829959 d0 ! #
      errtot(105,  1, 2, 2) =         0.0000029886 d0 
      FSE   (105,  1, 2, 3) =         0.6355044937 d0 ! #
      errtot(105,  1, 2, 3) =         0.0000004238 d0 
      FSE   (105,  1, 2, 4) =         0.6453794957 d0 ! #
      errtot(105,  1, 2, 4) =         0.0000010314 d0 
      FSE   (105,  1, 2, 5) =         0.6438321133 d0 ! #
      errtot(105,  1, 2, 5) =         0.0000007102 d0 
      FSE   (105,  1, 3, 3) =         0.6654876076 d0 ! #
      errtot(105,  1, 3, 3) =         0.0000005543 d0 
      FSE   (105,  1, 3, 4) =         0.6757905573 d0 ! #
      errtot(105,  1, 3, 4) =         0.0000006775 d0 
      FSE   (105,  1, 3, 5) =         0.6767781235 d0 ! #
      errtot(105,  1, 3, 5) =         0.0000013688 d0 
      FSE   (105,  1, 4, 4) =         0.6725256372 d0 ! #
      errtot(105,  1, 4, 4) =         0.0000005710 d0 
      FSE   (105,  1, 4, 5) =         0.6710564572 d0 ! #
      errtot(105,  1, 4, 5) =         0.0000039398 d0 
      FSE   (105,  1, 5, 5) =         0.6645105542 d0 ! #
      errtot(105,  1, 5, 5) =         0.0000130066 d0 
      FSE   (110,  1, 2, 2) =         0.7231010809 d0 ! #
      errtot(110,  1, 2, 2) =         0.0000024477 d0 
      FSE   (110,  1, 2, 3) =         0.7984133614 d0 ! #
      errtot(110,  1, 2, 3) =         0.0000003074 d0 
      FSE   (110,  1, 2, 4) =         0.8058805341 d0 ! #
      errtot(110,  1, 2, 4) =         0.0000008478 d0 
      FSE   (110,  1, 2, 5) =         0.8011809556 d0 ! #
      errtot(110,  1, 2, 5) =         0.0000006865 d0 
      FSE   (110,  1, 3, 3) =         0.8322064792 d0 ! #
      errtot(110,  1, 3, 3) =         0.0000010801 d0 
      FSE   (110,  1, 3, 4) =         0.8384901011 d0 ! #
      errtot(110,  1, 3, 4) =         0.0000023469 d0 
      FSE   (110,  1, 3, 5) =         0.8360715404 d0 ! #
      errtot(110,  1, 3, 5) =         0.0000015930 d0 
      FSE   (110,  1, 4, 4) =         0.8310625186 d0 ! #
      errtot(110,  1, 4, 4) =         0.0000070889 d0 
      FSE   (110,  1, 4, 5) =         0.8253565828 d0 ! #
      errtot(110,  1, 4, 5) =         0.0000092692 d0 
      FSE   (110,  1, 5, 5) =         0.8149667269 d0 ! #
      errtot(110,  1, 5, 5) =         0.0000165014 d0 
      FSE   (115,  1, 2, 2) =         0.9467890585 d0 ! #
      errtot(115,  1, 2, 2) =         0.0000004168 d0 
      FSE   (115,  1, 2, 3) =         1.0292351054 d0 ! #
      errtot(115,  1, 2, 3) =         0.0000028649 d0 
      FSE   (115,  1, 2, 4) =         1.0314488443 d0 ! #
      errtot(115,  1, 2, 4) =         0.0000041366 d0 
      FSE   (115,  1, 2, 5) =         1.0212266008 d0 ! #
      errtot(115,  1, 2, 5) =         0.0000048884 d0 
      FSE   (115,  1, 3, 3) =         1.0668034094 d0 ! #
      errtot(115,  1, 3, 3) =         0.0000042372 d0 
      FSE   (115,  1, 3, 4) =         1.0658634277 d0 ! #
      errtot(115,  1, 3, 4) =         0.0000053553 d0 
      FSE   (115,  1, 3, 5) =         1.0574963459 d0 ! #
      errtot(115,  1, 3, 5) =         0.0000052884 d0 
      FSE   (115,  1, 4, 4) =         1.0513011787 d0 ! #
      errtot(115,  1, 4, 4) =         0.0000036270 d0 
      FSE   (115,  1, 4, 5) =         1.0387997103 d0 ! #
      errtot(115,  1, 4, 5) =         0.0000090940 d0 
      FSE   (115,  1, 5, 5) =         1.0222265036 d0 ! #
      errtot(115,  1, 5, 5) =         0.0000223526 d0 
      FSE   (120,  1, 2, 2) =         1.2926467686 d0 ! #
      errtot(120,  1, 2, 2) =         0.0000098632 d0 
      FSE   (120,  1, 2, 3) =         1.3816381666 d0 ! #
      errtot(120,  1, 2, 3) =         0.0000152017 d0 
      FSE   (120,  1, 2, 4) =         1.3726430057 d0 ! #
      errtot(120,  1, 2, 4) =         0.0000165318 d0 
      FSE   (120,  1, 2, 5) =         1.3522046884 d0 ! #
      errtot(120,  1, 2, 5) =         0.0000233789 d0 
      FSE   (120,  1, 3, 3) =         1.4218772035 d0 ! #
      errtot(120,  1, 3, 3) =         0.0000195628 d0 
      FSE   (120,  1, 3, 4) =         1.4073774769 d0 ! #
      errtot(120,  1, 3, 4) =         0.0000211269 d0 
      FSE   (120,  1, 3, 5) =         1.3881257162 d0 ! #
      errtot(120,  1, 3, 5) =         0.0000205520 d0 
      FSE   (120,  1, 4, 4) =         1.3798169806 d0 ! #
      errtot(120,  1, 4, 4) =         0.0000313866 d0 
      FSE   (120,  1, 4, 5) =         1.3557326497 d0 ! #
      errtot(120,  1, 4, 5) =         0.0000090408 d0 
      FSE   (120,  1, 5, 5) =         1.3285516178 d0 ! #
      errtot(120,  1, 5, 5) =         0.0000423980 d0 
*
* np_{3/2}
*
      FSE   ( 10, -2, 2, 2) =         0.1303555156 d0 ! #
      errtot( 10, -2, 2, 2) =         0.0000093092 d0 
      FSE   ( 10, -2, 2, 3) =         0.1334816968 d0 ! #
      errtot( 10, -2, 2, 3) =         0.0000008618 d0 
      FSE   ( 10, -2, 2, 4) =         0.1315747788 d0 ! #
      errtot( 10, -2, 2, 4) =         0.0000005350 d0 
      FSE   ( 10, -2, 2, 5) =         0.1301409832 d0 ! #
      errtot( 10, -2, 2, 5) =         0.0000088665 d0 
      FSE   ( 10, -2, 3, 3) =         0.1420915638 d0 ! #
      errtot( 10, -2, 3, 3) =         0.0000024251 d0 
      FSE   ( 10, -2, 3, 4) =         0.1454332627 d0 ! #
      errtot( 10, -2, 3, 4) =         0.0000064408 d0 
      FSE   ( 10, -2, 3, 5) =         0.1456678050 d0 ! #
      errtot( 10, -2, 3, 5) =         0.0000087682 d0 
      FSE   ( 10, -2, 4, 4) =         0.1473972946 d0 ! #
      errtot( 10, -2, 4, 4) =         0.0000013900 d0 
      FSE   ( 10, -2, 4, 5) =         0.1496278285 d0 ! #
      errtot( 10, -2, 4, 5) =         0.0000025735 d0 
      FSE   ( 10, -2, 5, 5) =         0.1502972841 d0 ! #
      errtot( 10, -2, 5, 5) =         0.0000050186 d0 
      FSE   ( 15, -2, 2, 2) =         0.1365632182 d0 ! #
      errtot( 15, -2, 2, 2) =         0.0000011674 d0 
      FSE   ( 15, -2, 2, 3) =         0.1398708292 d0 ! #
      errtot( 15, -2, 2, 3) =         0.0000014896 d0 
      FSE   ( 15, -2, 2, 4) =         0.1379811945 d0 ! #
      errtot( 15, -2, 2, 4) =         0.0000027376 d0 
      FSE   ( 15, -2, 2, 5) =         0.1365453973 d0 ! #
      errtot( 15, -2, 2, 5) =         0.0000044155 d0 
      FSE   ( 15, -2, 3, 3) =         0.1490479918 d0 ! #
      errtot( 15, -2, 3, 3) =         0.0000003106 d0 
      FSE   ( 15, -2, 3, 4) =         0.1524640027 d0 ! #
      errtot( 15, -2, 3, 4) =         0.0000012122 d0 
      FSE   ( 15, -2, 3, 5) =         0.1527213536 d0 ! #
      errtot( 15, -2, 3, 5) =         0.0000051851 d0 
      FSE   ( 15, -2, 4, 4) =         0.1546012986 d0 ! #
      errtot( 15, -2, 4, 4) =         0.0000008067 d0 
      FSE   ( 15, -2, 4, 5) =         0.1568681995 d0 ! #
      errtot( 15, -2, 4, 5) =         0.0000059429 d0 
      FSE   ( 15, -2, 5, 5) =         0.1576144678 d0 ! #
      errtot( 15, -2, 5, 5) =         0.0000046057 d0 
      FSE   ( 20, -2, 2, 2) =         0.1438355353 d0 ! #
      errtot( 20, -2, 2, 2) =         0.0000005601 d0 
      FSE   ( 20, -2, 2, 3) =         0.1473472725 d0 ! #
      errtot( 20, -2, 2, 3) =         0.0000012930 d0 
      FSE   ( 20, -2, 2, 4) =         0.1454558714 d0 ! #
      errtot( 20, -2, 2, 4) =         0.0000007110 d0 
      FSE   ( 20, -2, 2, 5) =         0.1440081367 d0 ! #
      errtot( 20, -2, 2, 5) =         0.0000019370 d0 
      FSE   ( 20, -2, 3, 3) =         0.1571968950 d0 ! #
      errtot( 20, -2, 3, 3) =         0.0000029437 d0 
      FSE   ( 20, -2, 3, 4) =         0.1607270902 d0 ! #
      errtot( 20, -2, 3, 4) =         0.0000005165 d0 
      FSE   ( 20, -2, 3, 5) =         0.1610081639 d0 ! #
      errtot( 20, -2, 3, 5) =         0.0000026873 d0 
      FSE   ( 20, -2, 4, 4) =         0.1630496200 d0 ! #
      errtot( 20, -2, 4, 4) =         0.0000001933 d0 
      FSE   ( 20, -2, 4, 5) =         0.1653755804 d0 ! #
      errtot( 20, -2, 4, 5) =         0.0000042754 d0 
      FSE   ( 20, -2, 5, 5) =         0.1661945229 d0 ! #
      errtot( 20, -2, 5, 5) =         0.0000062230 d0 
      FSE   ( 25, -2, 2, 2) =         0.1519158671 d0 ! #
      errtot( 25, -2, 2, 2) =         0.0000012243 d0 
      FSE   ( 25, -2, 2, 3) =         0.1556473684 d0 ! #
      errtot( 25, -2, 2, 3) =         0.0000006404 d0 
      FSE   ( 25, -2, 2, 4) =         0.1537460826 d0 ! #
      errtot( 25, -2, 2, 4) =         0.0000009068 d0 
      FSE   ( 25, -2, 2, 5) =         0.1522774224 d0 ! #
      errtot( 25, -2, 2, 5) =         0.0000012831 d0 
      FSE   ( 25, -2, 3, 3) =         0.1662722968 d0 ! #
      errtot( 25, -2, 3, 3) =         0.0000014675 d0 
      FSE   ( 25, -2, 3, 4) =         0.1699410275 d0 ! #
      errtot( 25, -2, 3, 4) =         0.0000000869 d0 
      FSE   ( 25, -2, 3, 5) =         0.1702433646 d0 ! #
      errtot( 25, -2, 3, 5) =         0.0000012730 d0 
      FSE   ( 25, -2, 4, 4) =         0.1724766858 d0 ! #
      errtot( 25, -2, 4, 4) =         0.0000008696 d0 
      FSE   ( 25, -2, 4, 5) =         0.1748645221 d0 ! #
      errtot( 25, -2, 4, 5) =         0.0000029269 d0 
      FSE   ( 25, -2, 5, 5) =         0.1757615095 d0 ! #
      errtot( 25, -2, 5, 5) =         0.0000158693 d0 
      FSE   ( 30, -2, 2, 2) =         0.1606465252 d0 ! #
      errtot( 30, -2, 2, 2) =         0.0000009740 d0 
      FSE   ( 30, -2, 2, 3) =         0.1646129571 d0 ! #
      errtot( 30, -2, 2, 3) =         0.0000005751 d0 
      FSE   ( 30, -2, 2, 4) =         0.1626896565 d0 ! #
      errtot( 30, -2, 2, 4) =         0.0000006555 d0 
      FSE   ( 30, -2, 2, 5) =         0.1611896809 d0 ! #
      errtot( 30, -2, 2, 5) =         0.0000004396 d0 
      FSE   ( 30, -2, 3, 3) =         0.1761268854 d0 ! #
      errtot( 30, -2, 3, 3) =         0.0000007241 d0 
      FSE   ( 30, -2, 3, 4) =         0.1799380603 d0 ! #
      errtot( 30, -2, 3, 4) =         0.0000007380 d0 
      FSE   ( 30, -2, 3, 5) =         0.1802581848 d0 ! #
      errtot( 30, -2, 3, 5) =         0.0000006997 d0 
      FSE   ( 30, -2, 4, 4) =         0.1827097121 d0 ! #
      errtot( 30, -2, 4, 4) =         0.0000003704 d0 
      FSE   ( 30, -2, 4, 5) =         0.1851675581 d0 ! #
      errtot( 30, -2, 4, 5) =         0.0000006789 d0 
      FSE   ( 30, -2, 5, 5) =         0.1861484465 d0 ! #
      errtot( 30, -2, 5, 5) =         0.0000072246 d0 
      FSE   ( 35, -2, 2, 2) =         0.1699002568 d0 ! #
      errtot( 35, -2, 2, 2) =         0.0000008928 d0 
      FSE   ( 35, -2, 2, 3) =         0.1741323791 d0 ! #
      errtot( 35, -2, 2, 3) =         0.0000004047 d0 
      FSE   ( 35, -2, 2, 4) =         0.1721747463 d0 ! #
      errtot( 35, -2, 2, 4) =         0.0000002875 d0 
      FSE   ( 35, -2, 2, 5) =         0.1706323519 d0 ! #
      errtot( 35, -2, 2, 5) =         0.0000003379 d0 
      FSE   ( 35, -2, 3, 3) =         0.1866373215 d0 ! #
      errtot( 35, -2, 3, 3) =         0.0000004922 d0 
      FSE   ( 35, -2, 3, 4) =         0.1906063345 d0 ! #
      errtot( 35, -2, 3, 4) =         0.0000002672 d0 
      FSE   ( 35, -2, 3, 5) =         0.1909421252 d0 ! #
      errtot( 35, -2, 3, 5) =         0.0000015518 d0 
      FSE   ( 35, -2, 4, 4) =         0.1936364304 d0 ! #
      errtot( 35, -2, 4, 4) =         0.0000011093 d0 
      FSE   ( 35, -2, 4, 5) =         0.1961714708 d0 ! #
      errtot( 35, -2, 4, 5) =         0.0000009194 d0 
      FSE   ( 35, -2, 5, 5) =         0.1972414670 d0 ! #
      errtot( 35, -2, 5, 5) =         0.0000045323 d0 
      FSE   ( 40, -2, 2, 2) =         0.1795943967 d0 ! #
      errtot( 40, -2, 2, 2) =         0.0000006774 d0 
      FSE   ( 40, -2, 2, 3) =         0.1841258050 d0 ! #
      errtot( 40, -2, 2, 3) =         0.0000003504 d0 
      FSE   ( 40, -2, 2, 4) =         0.1821221829 d0 ! #
      errtot( 40, -2, 2, 4) =         0.0000002379 d0 
      FSE   ( 40, -2, 2, 5) =         0.1805263057 d0 ! #
      errtot( 40, -2, 2, 5) =         0.0000001013 d0 
      FSE   ( 40, -2, 3, 3) =         0.1977264637 d0 ! #
      errtot( 40, -2, 3, 3) =         0.0000001830 d0 
      FSE   ( 40, -2, 3, 4) =         0.2018706511 d0 ! #
      errtot( 40, -2, 3, 4) =         0.0000010736 d0 
      FSE   ( 40, -2, 3, 5) =         0.2022212587 d0 ! #
      errtot( 40, -2, 3, 5) =         0.0000010659 d0 
      FSE   ( 40, -2, 4, 4) =         0.2051683362 d0 ! #
      errtot( 40, -2, 4, 4) =         0.0000011927 d0 
      FSE   ( 40, -2, 4, 5) =         0.2078007247 d0 ! #
      errtot( 40, -2, 4, 5) =         0.0000010403 d0 
      FSE   ( 40, -2, 5, 5) =         0.2089664850 d0 ! #
      errtot( 40, -2, 5, 5) =         0.0000029729 d0 
      FSE   ( 45, -2, 2, 2) =         0.1896627924 d0 ! #
      errtot( 45, -2, 2, 2) =         0.0000007180 d0 
      FSE   ( 45, -2, 2, 3) =         0.1945358292 d0 ! #
      errtot( 45, -2, 2, 3) =         0.0000002032 d0 
      FSE   ( 45, -2, 2, 4) =         0.1924757530 d0 ! #
      errtot( 45, -2, 2, 4) =         0.0000002291 d0 
      FSE   ( 45, -2, 2, 5) =         0.1908149175 d0 ! #
      errtot( 45, -2, 2, 5) =         0.0000004161 d0 
      FSE   ( 45, -2, 3, 3) =         0.2093414681 d0 ! #
      errtot( 45, -2, 3, 3) =         0.0000001379 d0 
      FSE   ( 45, -2, 3, 4) =         0.2136840691 d0 ! #
      errtot( 45, -2, 3, 4) =         0.0000009603 d0 
      FSE   ( 45, -2, 3, 5) =         0.2140414445 d0 ! #
      errtot( 45, -2, 3, 5) =         0.0000007253 d0 
      FSE   ( 45, -2, 4, 4) =         0.2172805518 d0 ! #
      errtot( 45, -2, 4, 4) =         0.0000004553 d0 
      FSE   ( 45, -2, 4, 5) =         0.2200081812 d0 ! #
      errtot( 45, -2, 4, 5) =         0.0000001949 d0 
      FSE   ( 45, -2, 5, 5) =         0.2212733911 d0 ! #
      errtot( 45, -2, 5, 5) =         0.0000014397 d0 
      FSE   ( 50, -2, 2, 2) =         0.2000534182 d0 ! #
      errtot( 50, -2, 2, 2) =         0.0000006843 d0 
      FSE   ( 50, -2, 2, 3) =         0.2053206429 d0 ! #
      errtot( 50, -2, 2, 3) =         0.0000001992 d0 
      FSE   ( 50, -2, 2, 4) =         0.2031954359 d0 ! #
      errtot( 50, -2, 2, 4) =         0.0000004572 d0 
      FSE   ( 50, -2, 2, 5) =         0.2014585055 d0 ! #
      errtot( 50, -2, 2, 5) =         0.0000006829 d0 
      FSE   ( 50, -2, 3, 3) =         0.2214475603 d0 ! #
      errtot( 50, -2, 3, 3) =         0.0000006580 d0 
      FSE   ( 50, -2, 3, 4) =         0.2260098330 d0 ! #
      errtot( 50, -2, 3, 4) =         0.0000008001 d0 
      FSE   ( 50, -2, 3, 5) =         0.2263729412 d0 ! #
      errtot( 50, -2, 3, 5) =         0.0000010924 d0 
      FSE   ( 50, -2, 4, 4) =         0.2299329174 d0 ! #
      errtot( 50, -2, 4, 4) =         0.0000003969 d0 
      FSE   ( 50, -2, 4, 5) =         0.2327659934 d0 ! #
      errtot( 50, -2, 4, 5) =         0.0000004965 d0 
      FSE   ( 50, -2, 5, 5) =         0.2341348903 d0 ! #
      errtot( 50, -2, 5, 5) =         0.0000013189 d0 
      FSE   ( 55, -2, 2, 2) =         0.2107243878 d0 ! #
      errtot( 55, -2, 2, 2) =         0.0000006771 d0 
      FSE   ( 55, -2, 2, 3) =         0.2164500994 d0 ! #
      errtot( 55, -2, 2, 3) =         0.0000004224 d0 
      FSE   ( 55, -2, 2, 4) =         0.2142533245 d0 ! #
      errtot( 55, -2, 2, 4) =         0.0000005341 d0 
      FSE   ( 55, -2, 2, 5) =         0.2124296766 d0 ! #
      errtot( 55, -2, 2, 5) =         0.0000006648 d0 
      FSE   ( 55, -2, 3, 3) =         0.2340237976 d0 ! #
      errtot( 55, -2, 3, 3) =         0.0000003382 d0 
      FSE   ( 55, -2, 3, 4) =         0.2388320337 d0 ! #
      errtot( 55, -2, 3, 4) =         0.0000011135 d0 
      FSE   ( 55, -2, 3, 5) =         0.2392019390 d0 ! #
      errtot( 55, -2, 3, 5) =         0.0000008547 d0 
      FSE   ( 55, -2, 4, 4) =         0.2431117748 d0 ! #
      errtot( 55, -2, 4, 4) =         0.0000007830 d0 
      FSE   ( 55, -2, 4, 5) =         0.2460603758 d0 ! #
      errtot( 55, -2, 4, 5) =         0.0000006776 d0 
      FSE   ( 55, -2, 5, 5) =         0.2475374374 d0 ! #
      errtot( 55, -2, 5, 5) =         0.0000012777 d0 
      FSE   ( 60, -2, 2, 2) =         0.2216408285 d0 ! #
      errtot( 60, -2, 2, 2) =         0.0000008811 d0 
      FSE   ( 60, -2, 2, 3) =         0.2279029565 d0 ! #
      errtot( 60, -2, 2, 3) =         0.0000005394 d0 
      FSE   ( 60, -2, 2, 4) =         0.2256320457 d0 ! #
      errtot( 60, -2, 2, 4) =         0.0000002632 d0 
      FSE   ( 60, -2, 2, 5) =         0.2237115311 d0 ! #
      errtot( 60, -2, 2, 5) =         0.0000005942 d0 
      FSE   ( 60, -2, 3, 3) =         0.2470573740 d0 ! #
      errtot( 60, -2, 3, 3) =         0.0000005501 d0 
      FSE   ( 60, -2, 3, 4) =         0.2521459374 d0 ! #
      errtot( 60, -2, 3, 4) =         0.0000014052 d0 
      FSE   ( 60, -2, 3, 5) =         0.2525219684 d0 ! #
      errtot( 60, -2, 3, 5) =         0.0000009271 d0 
      FSE   ( 60, -2, 4, 4) =         0.2568115232 d0 ! #
      errtot( 60, -2, 4, 4) =         0.0000002141 d0 
      FSE   ( 60, -2, 4, 5) =         0.2598939255 d0 ! #
      errtot( 60, -2, 4, 5) =         0.0000023561 d0 
      FSE   ( 60, -2, 5, 5) =         0.2614808830 d0 ! #
      errtot( 60, -2, 5, 5) =         0.0000015444 d0 
      FSE   ( 65, -2, 2, 2) =         0.2327727938 d0 ! #
      errtot( 65, -2, 2, 2) =         0.0000005068 d0 
      FSE   ( 65, -2, 2, 3) =         0.2396639955 d0 ! #
      errtot( 65, -2, 2, 3) =         0.0000006036 d0 
      FSE   ( 65, -2, 2, 4) =         0.2373208381 d0 ! #
      errtot( 65, -2, 2, 4) =         0.0000001857 d0 
      FSE   ( 65, -2, 2, 5) =         0.2352944385 d0 ! #
      errtot( 65, -2, 2, 5) =         0.0000002083 d0 
      FSE   ( 65, -2, 3, 3) =         0.2605459255 d0 ! #
      errtot( 65, -2, 3, 3) =         0.0000006188 d0 
      FSE   ( 65, -2, 3, 4) =         0.2659545223 d0 ! #
      errtot( 65, -2, 3, 4) =         0.0000008747 d0 
      FSE   ( 65, -2, 3, 5) =         0.2663381674 d0 ! #
      errtot( 65, -2, 3, 5) =         0.0000013144 d0 
      FSE   ( 65, -2, 4, 4) =         0.2710402295 d0 ! #
      errtot( 65, -2, 4, 4) =         0.0000003290 d0 
      FSE   ( 65, -2, 4, 5) =         0.2742669703 d0 ! #
      errtot( 65, -2, 4, 5) =         0.0000043468 d0 
      FSE   ( 65, -2, 5, 5) =         0.2759766025 d0 ! #
      errtot( 65, -2, 5, 5) =         0.0000039268 d0 
      FSE   ( 70, -2, 2, 2) =         0.2440924375 d0 ! #
      errtot( 70, -2, 2, 2) =         0.0000006353 d0 
      FSE   ( 70, -2, 2, 3) =         0.2517228169 d0 ! #
      errtot( 70, -2, 2, 3) =         0.0000005347 d0 
      FSE   ( 70, -2, 2, 4) =         0.2493142998 d0 ! #
      errtot( 70, -2, 2, 4) =         0.0000002296 d0 
      FSE   ( 70, -2, 2, 5) =         0.2471743201 d0 ! #
      errtot( 70, -2, 2, 5) =         0.0000002565 d0 
      FSE   ( 70, -2, 3, 3) =         0.2744918445 d0 ! #
      errtot( 70, -2, 3, 3) =         0.0000008976 d0 
      FSE   ( 70, -2, 3, 4) =         0.2802655595 d0 ! #
      errtot( 70, -2, 3, 4) =         0.0000009240 d0 
      FSE   ( 70, -2, 3, 5) =         0.2806617171 d0 ! #
      errtot( 70, -2, 3, 5) =         0.0000003972 d0 
      FSE   ( 70, -2, 4, 4) =         0.2858108925 d0 ! #
      errtot( 70, -2, 4, 4) =         0.0000004710 d0 
      FSE   ( 70, -2, 4, 5) =         0.2892001983 d0 ! #
      errtot( 70, -2, 4, 5) =         0.0000020362 d0 
      FSE   ( 70, -2, 5, 5) =         0.2910405501 d0 ! #
      errtot( 70, -2, 5, 5) =         0.0000034446 d0 
      FSE   ( 75, -2, 2, 2) =         0.2555731103 d0 ! #
      errtot( 75, -2, 2, 2) =         0.0000004438 d0 
      FSE   ( 75, -2, 2, 3) =         0.2640717150 d0 ! #
      errtot( 75, -2, 2, 3) =         0.0000008928 d0 
      FSE   ( 75, -2, 2, 4) =         0.2616109532 d0 ! #
      errtot( 75, -2, 2, 4) =         0.0000004245 d0 
      FSE   ( 75, -2, 2, 5) =         0.2593515636 d0 ! #
      errtot( 75, -2, 2, 5) =         0.0000006784 d0 
      FSE   ( 75, -2, 3, 3) =         0.2889014876 d0 ! #
      errtot( 75, -2, 3, 3) =         0.0000010291 d0 
      FSE   ( 75, -2, 3, 4) =         0.2950937186 d0 ! #
      errtot( 75, -2, 3, 4) =         0.0000007755 d0 
      FSE   ( 75, -2, 3, 5) =         0.2955102787 d0 ! #
      errtot( 75, -2, 3, 5) =         0.0000012166 d0 
      FSE   ( 75, -2, 4, 4) =         0.3011424426 d0 ! #
      errtot( 75, -2, 4, 4) =         0.0000002069 d0 
      FSE   ( 75, -2, 4, 5) =         0.3047150182 d0 ! #
      errtot( 75, -2, 4, 5) =         0.0000020062 d0 
      FSE   ( 75, -2, 5, 5) =         0.3066963797 d0 ! #
      errtot( 75, -2, 5, 5) =         0.0000006256 d0 
      FSE   ( 80, -2, 2, 2) =         0.2671860633 d0 ! #
      errtot( 80, -2, 2, 2) =         0.0000004998 d0 
      FSE   ( 80, -2, 2, 3) =         0.2767026253 d0 ! #
      errtot( 80, -2, 2, 3) =         0.0000003754 d0 
      FSE   ( 80, -2, 2, 4) =         0.2742103327 d0 ! #
      errtot( 80, -2, 2, 4) =         0.0000004793 d0 
      FSE   ( 80, -2, 2, 5) =         0.2718283867 d0 ! #
      errtot( 80, -2, 2, 5) =         0.0000005024 d0 
      FSE   ( 80, -2, 3, 3) =         0.3037815869 d0 ! #
      errtot( 80, -2, 3, 3) =         0.0000010970 d0 
      FSE   ( 80, -2, 3, 4) =         0.3104547855 d0 ! #
      errtot( 80, -2, 3, 4) =         0.0000003789 d0 
      FSE   ( 80, -2, 3, 5) =         0.3109023172 d0 ! #
      errtot( 80, -2, 3, 5) =         0.0000016346 d0 
      FSE   ( 80, -2, 4, 4) =         0.3170562776 d0 ! #
      errtot( 80, -2, 4, 4) =         0.0000002585 d0 
      FSE   ( 80, -2, 4, 5) =         0.3208368880 d0 ! #
      errtot( 80, -2, 4, 5) =         0.0000040729 d0 
      FSE   ( 80, -2, 5, 5) =         0.3229705444 d0 ! #
      errtot( 80, -2, 5, 5) =         0.0000034505 d0 
      FSE   ( 85, -2, 2, 2) =         0.2788981621 d0 ! #
      errtot( 85, -2, 2, 2) =         0.0000006659 d0 
      FSE   ( 85, -2, 2, 3) =         0.2896041461 d0 ! #
      errtot( 85, -2, 2, 3) =         0.0000002717 d0 
      FSE   ( 85, -2, 2, 4) =         0.2871122544 d0 ! #
      errtot( 85, -2, 2, 4) =         0.0000002836 d0 
      FSE   ( 85, -2, 2, 5) =         0.2846082319 d0 ! #
      errtot( 85, -2, 2, 5) =         0.0000003253 d0 
      FSE   ( 85, -2, 3, 3) =         0.3191370091 d0 ! #
      errtot( 85, -2, 3, 3) =         0.0000011196 d0 
      FSE   ( 85, -2, 3, 4) =         0.3263624743 d0 ! #
      errtot( 85, -2, 3, 4) =         0.0000003693 d0 
      FSE   ( 85, -2, 3, 5) =         0.3268580210 d0 ! #
      errtot( 85, -2, 3, 5) =         0.0000009994 d0 
      FSE   ( 85, -2, 4, 4) =         0.3335724677 d0 ! #
      errtot( 85, -2, 4, 4) =         0.0000003846 d0 
      FSE   ( 85, -2, 4, 5) =         0.3375899257 d0 ! #
      errtot( 85, -2, 4, 5) =         0.0000016892 d0 
      FSE   ( 85, -2, 5, 5) =         0.3398894420 d0 ! #
      errtot( 85, -2, 5, 5) =         0.0000010019 d0 
      FSE   ( 90, -2, 2, 2) =         0.2906677659 d0 ! #
      errtot( 90, -2, 2, 2) =         0.0000006102 d0 
      FSE   ( 90, -2, 2, 3) =         0.3027577290 d0 ! #
      errtot( 90, -2, 2, 3) =         0.0000004236 d0 
      FSE   ( 90, -2, 2, 4) =         0.3003097995 d0 ! #
      errtot( 90, -2, 2, 4) =         0.0000002855 d0 
      FSE   ( 90, -2, 2, 5) =         0.2976882429 d0 ! #
      errtot( 90, -2, 2, 5) =         0.0000011723 d0 
      FSE   ( 90, -2, 3, 3) =         0.3349646230 d0 ! #
      errtot( 90, -2, 3, 3) =         0.0000011120 d0 
      FSE   ( 90, -2, 3, 4) =         0.3428240782 d0 ! #
      errtot( 90, -2, 3, 4) =         0.0000001368 d0 
      FSE   ( 90, -2, 3, 5) =         0.3433900746 d0 ! #
      errtot( 90, -2, 3, 5) =         0.0000019959 d0 
      FSE   ( 90, -2, 4, 4) =         0.3507043847 d0 ! #
      errtot( 90, -2, 4, 4) =         0.0000002932 d0 
      FSE   ( 90, -2, 4, 5) =         0.3549898461 d0 ! #
      errtot( 90, -2, 4, 5) =         0.0000069259 d0 
      FSE   ( 90, -2, 5, 5) =         0.3574710657 d0 ! #
      errtot( 90, -2, 5, 5) =         0.0000069556 d0 
      FSE   ( 95, -2, 2, 2) =         0.3024393229 d0 ! #
      errtot( 95, -2, 2, 2) =         0.0000004550 d0 
      FSE   ( 95, -2, 2, 3) =         0.3161296602 d0 ! #
      errtot( 95, -2, 2, 3) =         0.0000003401 d0 
      FSE   ( 95, -2, 2, 4) =         0.3137846012 d0 ! #
      errtot( 95, -2, 2, 4) =         0.0000004346 d0 
      FSE   ( 95, -2, 2, 5) =         0.3110556871 d0 ! #
      errtot( 95, -2, 2, 5) =         0.0000008298 d0 
      FSE   ( 95, -2, 3, 3) =         0.3512454342 d0 ! #
      errtot( 95, -2, 3, 3) =         0.0000009196 d0 
      FSE   ( 95, -2, 3, 4) =         0.3598308336 d0 ! #
      errtot( 95, -2, 3, 4) =         0.0000008915 d0 
      FSE   ( 95, -2, 3, 5) =         0.3604962072 d0 ! #
      errtot( 95, -2, 3, 5) =         0.0000020918 d0 
      FSE   ( 95, -2, 4, 4) =         0.3684540184 d0 ! #
      errtot( 95, -2, 4, 4) =         0.0000021704 d0 
      FSE   ( 95, -2, 4, 5) =         0.3730370141 d0 ! #
      errtot( 95, -2, 4, 5) =         0.0000021435 d0 
      FSE   ( 95, -2, 5, 5) =         0.3757177429 d0 ! #
      errtot( 95, -2, 5, 5) =         0.0000058347 d0 
      FSE   (100, -2, 2, 2) =         0.3141340050 d0 ! #
      errtot(100, -2, 2, 2) =         0.0000013199 d0 
      FSE   (100, -2, 2, 3) =         0.3296634093 d0 ! #
      errtot(100, -2, 2, 3) =         0.0000003163 d0 
      FSE   (100, -2, 2, 4) =         0.3274971896 d0 ! #
      errtot(100, -2, 2, 4) =         0.0000009798 d0 
      FSE   (100, -2, 2, 5) =         0.3246787141 d0 ! #
      errtot(100, -2, 2, 5) =         0.0000016223 d0 
      FSE   (100, -2, 3, 3) =         0.3679302963 d0 ! #
      errtot(100, -2, 3, 3) =         0.0000008706 d0 
      FSE   (100, -2, 3, 4) =         0.3773413766 d0 ! #
      errtot(100, -2, 3, 4) =         0.0000008818 d0 
      FSE   (100, -2, 3, 5) =         0.3781464646 d0 ! #
      errtot(100, -2, 3, 5) =         0.0000024070 d0 
      FSE   (100, -2, 4, 4) =         0.3867759314 d0 ! #
      errtot(100, -2, 4, 4) =         0.0000020525 d0 
      FSE   (100, -2, 4, 5) =         0.3917005015 d0 ! #
      errtot(100, -2, 4, 5) =         0.0000058757 d0 
      FSE   (100, -2, 5, 5) =         0.3946000176 d0 ! #
      errtot(100, -2, 5, 5) =         0.0000111605 d0 
      FSE   (105, -2, 2, 2) =         0.3256377411 d0 ! #
      errtot(105, -2, 2, 2) =         0.0000005442 d0 
      FSE   (105, -2, 2, 3) =         0.3432572878 d0 ! #
      errtot(105, -2, 2, 3) =         0.0000005157 d0 
      FSE   (105, -2, 2, 4) =         0.3413670121 d0 ! #
      errtot(105, -2, 2, 4) =         0.0000008815 d0 
      FSE   (105, -2, 2, 5) =         0.3384860615 d0 ! #
      errtot(105, -2, 2, 5) =         0.0000014820 d0 
      FSE   (105, -2, 3, 3) =         0.3849185070 d0 ! #
      errtot(105, -2, 3, 3) =         0.0000006792 d0 
      FSE   (105, -2, 3, 4) =         0.3952605807 d0 ! #
      errtot(105, -2, 3, 4) =         0.0000006836 d0 
      FSE   (105, -2, 3, 5) =         0.3962572037 d0 ! #
      errtot(105, -2, 3, 5) =         0.0000002067 d0 
      FSE   (105, -2, 4, 4) =         0.4055825646 d0 ! #
      errtot(105, -2, 4, 4) =         0.0000019612 d0 
      FSE   (105, -2, 4, 5) =         0.4108856038 d0 ! #
      errtot(105, -2, 4, 5) =         0.0000011261 d0 
      FSE   (105, -2, 5, 5) =         0.4140235138 d0 ! #
      errtot(105, -2, 5, 5) =         0.0000176800 d0 
      FSE   (110, -2, 2, 2) =         0.3367752137 d0 ! #
      errtot(110, -2, 2, 2) =         0.0000004056 d0 
      FSE   (110, -2, 2, 3) =         0.3567356085 d0 ! #
      errtot(110, -2, 2, 3) =         0.0000009897 d0 
      FSE   (110, -2, 2, 4) =         0.3552428521 d0 ! #
      errtot(110, -2, 2, 4) =         0.0000008439 d0 
      FSE   (110, -2, 2, 5) =         0.3523385474 d0 ! #
      errtot(110, -2, 2, 5) =         0.0000009292 d0 
      FSE   (110, -2, 3, 3) =         0.4020161009 d0 ! #
      errtot(110, -2, 3, 3) =         0.0000006346 d0 
      FSE   (110, -2, 3, 4) =         0.4133922399 d0 ! #
      errtot(110, -2, 3, 4) =         0.0000005466 d0 
      FSE   (110, -2, 3, 5) =         0.4146434491 d0 ! #
      errtot(110, -2, 3, 5) =         0.0000052556 d0 
      FSE   (110, -2, 4, 4) =         0.4246742488 d0 ! #
      errtot(110, -2, 4, 4) =         0.0000020074 d0 
      FSE   (110, -2, 4, 5) =         0.4303883246 d0 ! #
      errtot(110, -2, 4, 5) =         0.0000019065 d0 
      FSE   (110, -2, 5, 5) =         0.4337819291 d0 ! #
      errtot(110, -2, 5, 5) =         0.0000312746 d0 
      FSE   (115, -2, 2, 2) =         0.3472694455 d0 ! #
      errtot(115, -2, 2, 2) =         0.0000002748 d0 
      FSE   (115, -2, 2, 3) =         0.3697908655 d0 ! #
      errtot(115, -2, 2, 3) =         0.0000018570 d0 
      FSE   (115, -2, 2, 4) =         0.3688439095 d0 ! #
      errtot(115, -2, 2, 4) =         0.0000010530 d0 
      FSE   (115, -2, 2, 5) =         0.3659709743 d0 ! #
      errtot(115, -2, 2, 5) =         0.0000017340 d0 
      FSE   (115, -2, 3, 3) =         0.4188617233 d0 ! #
      errtot(115, -2, 3, 3) =         0.0000029467 d0 
      FSE   (115, -2, 3, 4) =         0.4313548628 d0 ! #
      errtot(115, -2, 3, 4) =         0.0000008340 d0 
      FSE   (115, -2, 3, 5) =         0.4329400815 d0 ! #
      errtot(115, -2, 3, 5) =         0.0000032238 d0 
      FSE   (115, -2, 4, 4) =         0.4436543598 d0 ! #
      errtot(115, -2, 4, 4) =         0.0000017806 d0 
      FSE   (115, -2, 4, 5) =         0.4497937916 d0 ! #
      errtot(115, -2, 4, 5) =         0.0000040763 d0 
      FSE   (115, -2, 5, 5) =         0.4534514048 d0 ! #
      errtot(115, -2, 5, 5) =         0.0000565282 d0 
      FSE   (120, -2, 2, 2) =         0.3566565647 d0 ! #
      errtot(120, -2, 2, 2) =         0.0000003431 d0 
      FSE   (120, -2, 2, 3) =         0.3818623479 d0 ! #
      errtot(120, -2, 2, 3) =         0.0000013928 d0 
      FSE   (120, -2, 2, 4) =         0.3816390071 d0 ! #
      errtot(120, -2, 2, 4) =         0.0000008193 d0 
      FSE   (120, -2, 2, 5) =         0.3788735523 d0 ! #
      errtot(120, -2, 2, 5) =         0.0000016952 d0 
      FSE   (120, -2, 3, 3) =         0.4347728689 d0 ! #
      errtot(120, -2, 3, 3) =         0.0000028434 d0 
      FSE   (120, -2, 3, 4) =         0.4484090429 d0 ! #
      errtot(120, -2, 3, 4) =         0.0000031094 d0 
      FSE   (120, -2, 3, 5) =         0.4504190434 d0 ! #
      errtot(120, -2, 3, 5) =         0.0000004959 d0 
      FSE   (120, -2, 4, 4) =         0.4617438408 d0 ! #
      errtot(120, -2, 4, 4) =         0.0000030405 d0 
      FSE   (120, -2, 4, 5) =         0.4682866352 d0 ! #
      errtot(120, -2, 4, 5) =         0.0000068442 d0 
      FSE   (120, -2, 5, 5) =         0.4721931021 d0 ! #
      errtot(120, -2, 5, 5) =         0.0001051288 d0 
*
* nd3/2
*
      FSE   ( 10,  2, 3, 3) =        -0.0426945670 d0 ! #
      errtot( 10,  2, 3, 3) =         0.0000738791 d0 
      FSE   ( 10,  2, 3, 4) =        -0.0351668516 d0 ! #
      errtot( 10,  2, 3, 4) =         0.0000206238 d0 
      FSE   ( 10,  2, 3, 5) =        -0.0343847138 d0 ! #
      errtot( 10,  2, 3, 5) =         0.0000836646 d0 
      FSE   ( 10,  2, 4, 4) =        -0.0406850693 d0 ! #
      errtot( 10,  2, 4, 4) =         0.0000328404 d0 
      FSE   ( 10,  2, 4, 5) =        -0.0371244931 d0 ! #
      errtot( 10,  2, 4, 5) =         0.0000368428 d0 
      FSE   ( 10,  2, 5, 5) =        -0.0395046782 d0 ! #
      errtot( 10,  2, 5, 5) =         0.0000546671 d0 
      FSE   ( 15,  2, 3, 3) =        -0.0424044600 d0 ! #
      errtot( 15,  2, 3, 3) =         0.0000096153 d0 
      FSE   ( 15,  2, 3, 4) =        -0.0348591086 d0 ! #
      errtot( 15,  2, 3, 4) =         0.0000004291 d0 
      FSE   ( 15,  2, 3, 5) =        -0.0340660500 d0 ! #
      errtot( 15,  2, 3, 5) =         0.0000054587 d0 
      FSE   ( 15,  2, 4, 4) =        -0.0403331380 d0 ! #
      errtot( 15,  2, 4, 4) =         0.0000111853 d0 
      FSE   ( 15,  2, 4, 5) =        -0.0367680444 d0 ! #
      errtot( 15,  2, 4, 5) =         0.0000124107 d0 
      FSE   ( 15,  2, 5, 5) =        -0.0391287090 d0 ! #
      errtot( 15,  2, 5, 5) =         0.0000399985 d0 
      FSE   ( 20,  2, 3, 3) =        -0.0420147408 d0 ! #
      errtot( 20,  2, 3, 3) =         0.0000006731 d0 
      FSE   ( 20,  2, 3, 4) =        -0.0344694485 d0 ! #
      errtot( 20,  2, 3, 4) =         0.0000007400 d0 
      FSE   ( 20,  2, 3, 5) =        -0.0336923562 d0 ! #
      errtot( 20,  2, 3, 5) =         0.0000013108 d0 
      FSE   ( 20,  2, 4, 4) =        -0.0398833227 d0 ! #
      errtot( 20,  2, 4, 4) =         0.0000064670 d0 
      FSE   ( 20,  2, 4, 5) =        -0.0363142666 d0 ! #
      errtot( 20,  2, 4, 5) =         0.0000038228 d0 
      FSE   ( 20,  2, 5, 5) =        -0.0386606406 d0 ! #
      errtot( 20,  2, 5, 5) =         0.0000160294 d0 
      FSE   ( 25,  2, 3, 3) =        -0.0415464226 d0 ! #
      errtot( 25,  2, 3, 3) =         0.0000014006 d0 
      FSE   ( 25,  2, 3, 4) =        -0.0339913029 d0 ! #
      errtot( 25,  2, 3, 4) =         0.0000014660 d0 
      FSE   ( 25,  2, 3, 5) =        -0.0332401766 d0 ! #
      errtot( 25,  2, 3, 5) =         0.0000013836 d0 
      FSE   ( 25,  2, 4, 4) =        -0.0393420962 d0 ! #
      errtot( 25,  2, 4, 4) =         0.0000146051 d0 
      FSE   ( 25,  2, 4, 5) =        -0.0357599028 d0 ! #
      errtot( 25,  2, 4, 5) =         0.0000008170 d0 
      FSE   ( 25,  2, 5, 5) =        -0.0380901635 d0 ! #
      errtot( 25,  2, 5, 5) =         0.0000065511 d0 
      FSE   ( 30,  2, 3, 3) =        -0.0409967365 d0 ! #
      errtot( 30,  2, 3, 3) =         0.0000002949 d0 
      FSE   ( 30,  2, 3, 4) =        -0.0334302225 d0 ! #
      errtot( 30,  2, 3, 4) =         0.0000029740 d0 
      FSE   ( 30,  2, 3, 5) =        -0.0326926269 d0 ! #
      errtot( 30,  2, 3, 5) =         0.0000009788 d0 
      FSE   ( 30,  2, 4, 4) =        -0.0387022608 d0 ! #
      errtot( 30,  2, 4, 4) =         0.0000003582 d0 
      FSE   ( 30,  2, 4, 5) =        -0.0350922902 d0 ! #
      errtot( 30,  2, 4, 5) =         0.0000019277 d0 
      FSE   ( 30,  2, 5, 5) =        -0.0373972556 d0 ! #
      errtot( 30,  2, 5, 5) =         0.0000077586 d0 
      FSE   ( 35,  2, 3, 3) =        -0.0403533441 d0 ! #
      errtot( 35,  2, 3, 3) =         0.0000013058 d0 
      FSE   ( 35,  2, 3, 4) =        -0.0327588032 d0 ! #
      errtot( 35,  2, 3, 4) =         0.0000005887 d0 
      FSE   ( 35,  2, 3, 5) =        -0.0320416628 d0 ! #
      errtot( 35,  2, 3, 5) =         0.0000014802 d0 
      FSE   ( 35,  2, 4, 4) =        -0.0379443238 d0 ! #
      errtot( 35,  2, 4, 4) =         0.0000002981 d0 
      FSE   ( 35,  2, 4, 5) =        -0.0343035786 d0 ! #
      errtot( 35,  2, 4, 5) =         0.0000016878 d0 
      FSE   ( 35,  2, 5, 5) =        -0.0365751334 d0 ! #
      errtot( 35,  2, 5, 5) =         0.0000067819 d0 
      FSE   ( 40,  2, 3, 3) =        -0.0396145187 d0 ! #
      errtot( 40,  2, 3, 3) =         0.0000017729 d0 
      FSE   ( 40,  2, 3, 4) =        -0.0319676660 d0 ! #
      errtot( 40,  2, 3, 4) =         0.0000012071 d0 
      FSE   ( 40,  2, 3, 5) =        -0.0312697740 d0 ! #
      errtot( 40,  2, 3, 5) =         0.0000012369 d0 
      FSE   ( 40,  2, 4, 4) =        -0.0370543185 d0 ! #
      errtot( 40,  2, 4, 4) =         0.0000007537 d0 
      FSE   ( 40,  2, 4, 5) =        -0.0333683971 d0 ! #
      errtot( 40,  2, 4, 5) =         0.0000017361 d0 
      FSE   ( 40,  2, 5, 5) =        -0.0356055707 d0 ! #
      errtot( 40,  2, 5, 5) =         0.0000036520 d0 
      FSE   ( 45,  2, 3, 3) =        -0.0387562875 d0 ! #
      errtot( 45,  2, 3, 3) =         0.0000012114 d0 
      FSE   ( 45,  2, 3, 4) =        -0.0310405888 d0 ! #
      errtot( 45,  2, 3, 4) =         0.0000007828 d0 
      FSE   ( 45,  2, 3, 5) =        -0.0303602290 d0 ! #
      errtot( 45,  2, 3, 5) =         0.0000007623 d0 
      FSE   ( 45,  2, 4, 4) =        -0.0360123024 d0 ! #
      errtot( 45,  2, 4, 4) =         0.0000003299 d0 
      FSE   ( 45,  2, 4, 5) =        -0.0322679342 d0 ! #
      errtot( 45,  2, 4, 5) =         0.0000006742 d0 
      FSE   ( 45,  2, 5, 5) =        -0.0344629915 d0 ! #
      errtot( 45,  2, 5, 5) =         0.0000024868 d0 
      FSE   ( 50,  2, 3, 3) =        -0.0377689999 d0 ! #
      errtot( 50,  2, 3, 3) =         0.0000010336 d0 
      FSE   ( 50,  2, 3, 4) =        -0.0299580533 d0 ! #
      errtot( 50,  2, 3, 4) =         0.0000008338 d0 
      FSE   ( 50,  2, 3, 5) =        -0.0292925527 d0 ! #
      errtot( 50,  2, 3, 5) =         0.0000003229 d0 
      FSE   ( 50,  2, 4, 4) =        -0.0347959482 d0 ! #
      errtot( 50,  2, 4, 4) =         0.0000003516 d0 
      FSE   ( 50,  2, 4, 5) =        -0.0309778852 d0 ! #
      errtot( 50,  2, 4, 5) =         0.0000011392 d0 
      FSE   ( 50,  2, 5, 5) =        -0.0331240919 d0 ! #
      errtot( 50,  2, 5, 5) =         0.0000018402 d0 
      FSE   ( 55,  2, 3, 3) =        -0.0366327635 d0 ! #
      errtot( 55,  2, 3, 3) =         0.0000010044 d0 
      FSE   ( 55,  2, 3, 4) =        -0.0286977579 d0 ! #
      errtot( 55,  2, 3, 4) =         0.0000007833 d0 
      FSE   ( 55,  2, 3, 5) =        -0.0280395326 d0 ! #
      errtot( 55,  2, 3, 5) =         0.0000002017 d0 
      FSE   ( 55,  2, 4, 4) =        -0.0333804078 d0 ! #
      errtot( 55,  2, 4, 4) =         0.0000006756 d0 
      FSE   ( 55,  2, 4, 5) =        -0.0294700043 d0 ! #
      errtot( 55,  2, 4, 5) =         0.0000011295 d0 
      FSE   ( 55,  2, 5, 5) =        -0.0316005062 d0 ! #
      errtot( 55,  2, 5, 5) =         0.0000175736 d0 
      FSE   ( 60,  2, 3, 3) =        -0.0353288622 d0 ! #
      errtot( 60,  2, 3, 3) =         0.0000009720 d0 
      FSE   ( 60,  2, 3, 4) =        -0.0272342273 d0 ! #
      errtot( 60,  2, 3, 4) =         0.0000002768 d0 
      FSE   ( 60,  2, 3, 5) =        -0.0265822080 d0 ! #
      errtot( 60,  2, 3, 5) =         0.0000009346 d0 
      FSE   ( 60,  2, 4, 4) =        -0.0317367205 d0 ! #
      errtot( 60,  2, 4, 4) =         0.0000008137 d0 
      FSE   ( 60,  2, 4, 5) =        -0.0277135459 d0 ! #
      errtot( 60,  2, 4, 5) =         0.0000007179 d0 
      FSE   ( 60,  2, 5, 5) =        -0.0297738809 d0 ! #
      errtot( 60,  2, 5, 5) =         0.0000164018 d0 
      FSE   ( 65,  2, 3, 3) =        -0.0338348325 d0 ! #
      errtot( 65,  2, 3, 3) =         0.0000008020 d0 
      FSE   ( 65,  2, 3, 4) =        -0.0255381831 d0 ! #
      errtot( 65,  2, 3, 4) =         0.0000001063 d0 
      FSE   ( 65,  2, 3, 5) =        -0.0248826581 d0 ! #
      errtot( 65,  2, 3, 5) =         0.0000008586 d0 
      FSE   ( 65,  2, 4, 4) =        -0.0298317190 d0 ! #
      errtot( 65,  2, 4, 4) =         0.0000009208 d0 
      FSE   ( 65,  2, 4, 5) =        -0.0256703405 d0 ! #
      errtot( 65,  2, 4, 5) =         0.0000031204 d0 
      FSE   ( 65,  2, 5, 5) =        -0.0276495831 d0 ! #
      errtot( 65,  2, 5, 5) =         0.0000156335 d0 
      FSE   ( 70,  2, 3, 3) =        -0.0321246775 d0 ! #
      errtot( 70,  2, 3, 3) =         0.0000008372 d0 
      FSE   ( 70,  2, 3, 4) =        -0.0235769044 d0 ! #
      errtot( 70,  2, 3, 4) =         0.0000001392 d0 
      FSE   ( 70,  2, 3, 5) =        -0.0229079881 d0 ! #
      errtot( 70,  2, 3, 5) =         0.0000004550 d0 
      FSE   ( 70,  2, 4, 4) =        -0.0276280230 d0 ! #
      errtot( 70,  2, 4, 4) =         0.0000010320 d0 
      FSE   ( 70,  2, 4, 5) =        -0.0233085222 d0 ! #
      errtot( 70,  2, 4, 5) =         0.0000048648 d0 
      FSE   ( 70,  2, 5, 5) =        -0.0251853851 d0 ! #
      errtot( 70,  2, 5, 5) =         0.0000145573 d0 
      FSE   ( 75,  2, 3, 3) =        -0.0301695234 d0 ! #
      errtot( 75,  2, 3, 3) =         0.0000008518 d0 
      FSE   ( 75,  2, 3, 4) =        -0.0213124518 d0 ! #
      errtot( 75,  2, 3, 4) =         0.0000007632 d0 
      FSE   ( 75,  2, 3, 5) =        -0.0206182825 d0 ! #
      errtot( 75,  2, 3, 5) =         0.0000005156 d0 
      FSE   ( 75,  2, 4, 4) =        -0.0250830736 d0 ! #
      errtot( 75,  2, 4, 4) =         0.0000010494 d0 
      FSE   ( 75,  2, 4, 5) =        -0.0205551020 d0 ! #
      errtot( 75,  2, 4, 5) =         0.0000011846 d0 
      FSE   ( 75,  2, 5, 5) =        -0.0223316096 d0 ! #
      errtot( 75,  2, 5, 5) =         0.0000136805 d0 
      FSE   ( 80,  2, 3, 3) =        -0.0279365352 d0 ! #
      errtot( 80,  2, 3, 3) =         0.0000007053 d0 
      FSE   ( 80,  2, 3, 4) =        -0.0187037808 d0 ! #
      errtot( 80,  2, 3, 4) =         0.0000007641 d0 
      FSE   ( 80,  2, 3, 5) =        -0.0179675727 d0 ! #
      errtot( 80,  2, 3, 5) =         0.0000008292 d0 
      FSE   ( 80,  2, 4, 4) =        -0.0221482632 d0 ! #
      errtot( 80,  2, 4, 4) =         0.0000010371 d0 
      FSE   ( 80,  2, 4, 5) =        -0.0173881170 d0 ! #
      errtot( 80,  2, 4, 5) =         0.0000044876 d0 
      FSE   ( 80,  2, 5, 5) =        -0.0190316354 d0 ! #
      errtot( 80,  2, 5, 5) =         0.0000143240 d0 
      FSE   ( 85,  2, 3, 3) =        -0.0253880383 d0 ! #
      errtot( 85,  2, 3, 3) =         0.0000007702 d0 
      FSE   ( 85,  2, 3, 4) =        -0.0157018433 d0 ! #
      errtot( 85,  2, 3, 4) =         0.0000004181 d0 
      FSE   ( 85,  2, 3, 5) =        -0.0149030353 d0 ! #
      errtot( 85,  2, 3, 5) =         0.0000000778 d0 
      FSE   ( 85,  2, 4, 4) =        -0.0187675851 d0 ! #
      errtot( 85,  2, 4, 4) =         0.0000011024 d0 
      FSE   ( 85,  2, 4, 5) =        -0.0137229322 d0 ! #
      errtot( 85,  2, 4, 5) =         0.0000034912 d0 
      FSE   ( 85,  2, 5, 5) =        -0.0152204560 d0 ! #
      errtot( 85,  2, 5, 5) =         0.0000138834 d0 
      FSE   ( 90,  2, 3, 3) =        -0.0224815986 d0 ! #
      errtot( 90,  2, 3, 3) =         0.0000008025 d0 
      FSE   ( 90,  2, 3, 4) =        -0.0122508737 d0 ! #
      errtot( 90,  2, 3, 4) =         0.0000002935 d0 
      FSE   ( 90,  2, 3, 5) =        -0.0113640737 d0 ! #
      errtot( 90,  2, 3, 5) =         0.0000010548 d0 
      FSE   ( 90,  2, 4, 4) =        -0.0148773719 d0 ! #
      errtot( 90,  2, 4, 4) =         0.0000009303 d0 
      FSE   ( 90,  2, 4, 5) =        -0.0094967721 d0 ! #
      errtot( 90,  2, 4, 5) =         0.0000054307 d0 
      FSE   ( 90,  2, 5, 5) =        -0.0108243848 d0 ! #
      errtot( 90,  2, 5, 5) =         0.0000151062 d0 
      FSE   ( 95,  2, 3, 3) =        -0.0191691018 d0 ! #
      errtot( 95,  2, 3, 3) =         0.0000006726 d0 
      FSE   ( 95,  2, 3, 4) =        -0.0082878399 d0 ! #
      errtot( 95,  2, 3, 4) =         0.0000001263 d0 
      FSE   ( 95,  2, 3, 5) =        -0.0072804008 d0 ! #
      errtot( 95,  2, 3, 5) =         0.0000010458 d0 
      FSE   ( 95,  2, 4, 4) =        -0.0104045765 d0 ! #
      errtot( 95,  2, 4, 4) =         0.0000010190 d0 
      FSE   ( 95,  2, 4, 5) =        -0.0046284796 d0 ! #
      errtot( 95,  2, 4, 5) =         0.0000029491 d0 
      FSE   ( 95,  2, 5, 5) =        -0.0057584517 d0 ! #
      errtot( 95,  2, 5, 5) =         0.0000151748 d0 
      FSE   (100,  2, 3, 3) =        -0.0153961211 d0 ! #
      errtot(100,  2, 3, 3) =         0.0000005677 d0 
      FSE   (100,  2, 3, 4) =        -0.0037411059 d0 ! #
      errtot(100,  2, 3, 4) =         0.0000000553 d0 
      FSE   (100,  2, 3, 5) =        -0.0025707845 d0 ! #
      errtot(100,  2, 3, 5) =         0.0000008620 d0 
      FSE   (100,  2, 4, 4) =        -0.0052664750 d0 ! #
      errtot(100,  2, 4, 4) =         0.0000008523 d0 
      FSE   (100,  2, 4, 5) =         0.0009747135 d0 ! #
      errtot(100,  2, 4, 5) =         0.0000049931 d0 
      FSE   (100,  2, 5, 5) =         0.0000746515 d0 ! #
      errtot(100,  2, 5, 5) =         0.0000181975 d0 
      FSE   (105,  2, 3, 3) =        -0.0111025995 d0 ! #
      errtot(105,  2, 3, 3) =         0.0000003583 d0 
      FSE   (105,  2, 3, 4) =         0.0014689646 d0 ! #
      errtot(105,  2, 3, 4) =         0.0000000230 d0 
      FSE   (105,  2, 3, 5) =         0.0028560792 d0 ! #
      errtot(105,  2, 3, 5) =         0.0000016140 d0 
      FSE   (105,  2, 4, 4) =         0.0006284954 d0 ! #
      errtot(105,  2, 4, 4) =         0.0000008030 d0 
      FSE   (105,  2, 4, 5) =         0.0074139479 d0 ! #
      errtot(105,  2, 4, 5) =         0.0000021225 d0 
      FSE   (105,  2, 5, 5) =         0.0067804282 d0 ! #
      errtot(105,  2, 5, 5) =         0.0000214135 d0 
      FSE   (110,  2, 3, 3) =        -0.0062243224 d0 ! #
      errtot(110,  2, 3, 3) =         0.0000002229 d0 
      FSE   (110,  2, 3, 4) =         0.0074277815 d0 ! #
      errtot(110,  2, 3, 4) =         0.0000002150 d0 
      FSE   (110,  2, 3, 5) =         0.0091039796 d0 ! #
      errtot(110,  2, 3, 5) =         0.0000012846 d0 
      FSE   (110,  2, 4, 4) =         0.0073779115 d0 ! #
      errtot(110,  2, 4, 4) =         0.0000022577 d0 
      FSE   (110,  2, 4, 5) =         0.0147961458 d0 ! #
      errtot(110,  2, 4, 5) =         0.0000021945 d0 
      FSE   (110,  2, 5, 5) =         0.0144709514 d0 ! #
      errtot(110,  2, 5, 5) =         0.0000295358 d0 
      FSE   (115,  2, 3, 3) =        -0.0006991989 d0 ! #
      errtot(115,  2, 3, 3) =         0.0000022060 d0 
      FSE   (115,  2, 3, 4) =         0.0142184868 d0 ! #
      errtot(115,  2, 3, 4) =         0.0000010402 d0 
      FSE   (115,  2, 3, 5) =         0.0162787413 d0 ! #
      errtot(115,  2, 3, 5) =         0.0000008118 d0 
      FSE   (115,  2, 4, 4) =         0.0150745744 d0 ! #
      errtot(115,  2, 4, 4) =         0.0000007413 d0 
      FSE   (115,  2, 4, 5) =         0.0232225816 d0 ! #
      errtot(115,  2, 4, 5) =         0.0000007847 d0 
      FSE   (115,  2, 5, 5) =         0.0232472091 d0 ! #
      errtot(115,  2, 5, 5) =         0.0000450997 d0 
      FSE   (120,  2, 3, 3) =         0.0055184937 d0 ! #
      errtot(120,  2, 3, 3) =         0.0000024593 d0 
      FSE   (120,  2, 3, 4) =         0.0218988224 d0 ! #
      errtot(120,  2, 3, 4) =         0.0000054275 d0 
      FSE   (120,  2, 3, 5) =         0.0244719476 d0 ! #
      errtot(120,  2, 3, 5) =         0.0000019589 d0 
      FSE   (120,  2, 4, 4) =         0.0237641674 d0 ! #
      errtot(120,  2, 4, 4) =         0.0000067330 d0 
      FSE   (120,  2, 4, 5) =         0.0327510062 d0 ! #
      errtot(120,  2, 4, 5) =         0.0000040739 d0 
      FSE   (120,  2, 5, 5) =         0.0331637266 d0 ! #
      errtot(120,  2, 5, 5) =         0.0000767216 d0 
*
* nd5/2
*
      FSE   ( 10, -3, 3, 3) =         0.0407667086 d0 ! #
      errtot( 10, -3, 3, 3) =         0.0000013252 d0 
      FSE   ( 10, -3, 3, 4) =         0.0376710564 d0 ! #
      errtot( 10, -3, 3, 4) =         0.0000062916 d0 
      FSE   ( 10, -3, 3, 5) =         0.0354113675 d0 ! #
      errtot( 10, -3, 3, 5) =         0.0000247526 d0 
      FSE   ( 10, -3, 4, 4) =         0.0428180287 d0 ! #
      errtot( 10, -3, 4, 4) =         0.0000024583 d0 
      FSE   ( 10, -3, 4, 5) =         0.0424007731 d0 ! #
      errtot( 10, -3, 4, 5) =         0.0000074780 d0
      FSE   ( 10, -3, 5, 5) =         0.0440 d0 ! #
      errtot( 10, -3, 5, 5) =         0.d0
      FSE   ( 15, -3, 3, 3) =         0.0411846556 d0 ! #
      errtot( 15, -3, 3, 3) =         0.0000129797 d0 
      FSE   ( 15, -3, 3, 4) =         0.0380676083 d0 ! #
      errtot( 15, -3, 3, 4) =         0.0000020409 d0 
      FSE   ( 15, -3, 3, 5) =         0.0358258648 d0 ! #
      errtot( 15, -3, 3, 5) =         0.0000132237 d0 
      FSE   ( 15, -3, 4, 4) =         0.0433335786 d0 ! #
      errtot( 15, -3, 4, 4) =         0.0000053484 d0 
      FSE   ( 15, -3, 4, 5) =         0.0428853427 d0 ! #
      errtot( 15, -3, 4, 5) =         0.0000191085 d0 
      FSE   ( 15, -3, 5, 5) =         0.0445266060 d0 ! #
      errtot( 15, -3, 5, 5) =         0.0000013919 d0 
      FSE   ( 20, -3, 3, 3) =         0.0417489358 d0 ! #
      errtot( 20, -3, 3, 3) =         0.0000014089 d0 
      FSE   ( 20, -3, 3, 4) =         0.0386191102 d0 ! #
      errtot( 20, -3, 3, 4) =         0.0000005283 d0 
      FSE   ( 20, -3, 3, 5) =         0.0363304824 d0 ! #
      errtot( 20, -3, 3, 5) =         0.0000022518 d0 
      FSE   ( 20, -3, 4, 4) =         0.0439660809 d0 ! #
      errtot( 20, -3, 4, 4) =         0.0000041632 d0 
      FSE   ( 20, -3, 4, 5) =         0.0435316301 d0 ! #
      errtot( 20, -3, 4, 5) =         0.0000038674 d0 
      FSE   ( 20, -3, 5, 5) =         0.0452150438 d0 ! #
      errtot( 20, -3, 5, 5) =         0.0000050023 d0 
      FSE   ( 25, -3, 3, 3) =         0.0424344728 d0 ! #
      errtot( 25, -3, 3, 3) =         0.0000001033 d0 
      FSE   ( 25, -3, 3, 4) =         0.0392915945 d0 ! #
      errtot( 25, -3, 3, 4) =         0.0000007590 d0 
      FSE   ( 25, -3, 3, 5) =         0.0369818699 d0 ! #
      errtot( 25, -3, 3, 5) =         0.0000011427 d0 
      FSE   ( 25, -3, 4, 4) =         0.0447514104 d0 ! #
      errtot( 25, -3, 4, 4) =         0.0000088477 d0 
      FSE   ( 25, -3, 4, 5) =         0.0443255611 d0 ! #
      errtot( 25, -3, 4, 5) =         0.0000009382 d0 
      FSE   ( 25, -3, 5, 5) =         0.0460362107 d0 ! #
      errtot( 25, -3, 5, 5) =         0.0000151693 d0 
      FSE   ( 30, -3, 3, 3) =         0.0432349890 d0 ! #
      errtot( 30, -3, 3, 3) =         0.0000007535 d0 
      FSE   ( 30, -3, 3, 4) =         0.0400805154 d0 ! #
      errtot( 30, -3, 3, 4) =         0.0000000685 d0 
      FSE   ( 30, -3, 3, 5) =         0.0377484863 d0 ! #
      errtot( 30, -3, 3, 5) =         0.0000010437 d0 
      FSE   ( 30, -3, 4, 4) =         0.0456776298 d0 ! #
      errtot( 30, -3, 4, 4) =         0.0000007852 d0 
      FSE   ( 30, -3, 4, 5) =         0.0452591015 d0 ! #
      errtot( 30, -3, 4, 5) =         0.0000014207 d0 
      FSE   ( 30, -3, 5, 5) =         0.0470172909 d0 ! #
      errtot( 30, -3, 5, 5) =         0.0000102828 d0 
      FSE   ( 35, -3, 3, 3) =         0.0441489645 d0 ! #
      errtot( 35, -3, 3, 3) =         0.0000002210 d0 
      FSE   ( 35, -3, 3, 4) =         0.0409792394 d0 ! #
      errtot( 35, -3, 3, 4) =         0.0000009255 d0 
      FSE   ( 35, -3, 3, 5) =         0.0386210780 d0 ! #
      errtot( 35, -3, 3, 5) =         0.0000011405 d0 
      FSE   ( 35, -3, 4, 4) =         0.0467369786 d0 ! #
      errtot( 35, -3, 4, 4) =         0.0000023860 d0 
      FSE   ( 35, -3, 4, 5) =         0.0463286301 d0 ! #
      errtot( 35, -3, 4, 5) =         0.0000004824 d0 
      FSE   ( 35, -3, 5, 5) =         0.0481433930 d0 ! #
      errtot( 35, -3, 5, 5) =         0.0000045925 d0 
      FSE   ( 40, -3, 3, 3) =         0.0451704883 d0 ! #
      errtot( 40, -3, 3, 3) =         0.0000008231 d0 
      FSE   ( 40, -3, 3, 4) =         0.0419901151 d0 ! #
      errtot( 40, -3, 3, 4) =         0.0000012273 d0 
      FSE   ( 40, -3, 3, 5) =         0.0396040659 d0 ! #
      errtot( 40, -3, 3, 5) =         0.0000008424 d0 
      FSE   ( 40, -3, 4, 4) =         0.0479314290 d0 ! #
      errtot( 40, -3, 4, 4) =         0.0000008815 d0 
      FSE   ( 40, -3, 4, 5) =         0.0475368385 d0 ! #
      errtot( 40, -3, 4, 5) =         0.0000009332 d0 
      FSE   ( 45, -3, 3, 3) =         0.0463008560 d0 ! #
      errtot( 45, -3, 3, 3) =         0.0000014141 d0 
      FSE   ( 45, -3, 3, 4) =         0.0431138407 d0 ! #
      errtot( 45, -3, 3, 4) =         0.0000006100 d0 
      FSE   ( 45, -3, 3, 5) =         0.0406967365 d0 ! #
      errtot( 45, -3, 3, 5) =         0.0000007355 d0 
      FSE   ( 45, -3, 4, 4) =         0.0492622535 d0 ! #
      errtot( 45, -3, 4, 4) =         0.0000009478 d0 
      FSE   ( 45, -3, 4, 5) =         0.0488848438 d0 ! #
      errtot( 45, -3, 4, 5) =         0.0000013266 d0 
      FSE   ( 50, -3, 3, 3) =         0.0475390690 d0 ! #
      errtot( 50, -3, 3, 3) =         0.0000011340 d0 
      FSE   ( 50, -3, 3, 4) =         0.0443495017 d0 ! #
      errtot( 50, -3, 3, 4) =         0.0000006748 d0 
      FSE   ( 50, -3, 3, 5) =         0.0419006388 d0 ! #
      errtot( 50, -3, 3, 5) =         0.0000002351 d0 
      FSE   ( 50, -3, 4, 4) =         0.0507304254 d0 ! #
      errtot( 50, -3, 4, 4) =         0.0000007729 d0 
      FSE   ( 50, -3, 4, 5) =         0.0503739566 d0 ! #
      errtot( 50, -3, 4, 5) =         0.0000014937 d0 
      FSE   ( 50, -3, 5, 5) =         0.0524278775 d0 ! #
      errtot( 50, -3, 5, 5) =         0.0000065710 d0 
      FSE   ( 55, -3, 3, 3) =         0.0488932398 d0 ! #
      errtot( 55, -3, 3, 3) =         0.0000008571 d0 
      FSE   ( 55, -3, 3, 4) =         0.0456982186 d0 ! #
      errtot( 55, -3, 3, 4) =         0.0000007071 d0 
      FSE   ( 55, -3, 3, 5) =         0.0432160433 d0 ! #
      errtot( 55, -3, 3, 5) =         0.0000006516 d0 
      FSE   ( 55, -3, 4, 4) =         0.0523371381 d0 ! #
      errtot( 55, -3, 4, 4) =         0.0000008211 d0 
      FSE   ( 55, -3, 4, 5) =         0.0520052381 d0 ! #
      errtot( 55, -3, 4, 5) =         0.0000011553 d0 
      FSE   ( 55, -3, 5, 5) =         0.0541504469 d0 ! #
      errtot( 55, -3, 5, 5) =         0.0000064686 d0 
      FSE   ( 60, -3, 3, 3) =         0.0503478493 d0 ! #
      errtot( 60, -3, 3, 3) =         0.0000009038 d0 
      FSE   ( 60, -3, 3, 4) =         0.0471565118 d0 ! #
      errtot( 60, -3, 3, 4) =         0.0000010405 d0 
      FSE   ( 60, -3, 3, 5) =         0.0446423389 d0 ! #
      errtot( 60, -3, 3, 5) =         0.0000009164 d0 
      FSE   ( 60, -3, 4, 4) =         0.0540809961 d0 ! #
      errtot( 60, -3, 4, 4) =         0.0000005016 d0 
      FSE   ( 60, -3, 4, 5) =         0.0537888842 d0 ! #
      errtot( 60, -3, 4, 5) =         0.0000048588 d0 
      FSE   ( 60, -3, 5, 5) =         0.0560254080 d0 ! #
      errtot( 60, -3, 5, 5) =         0.0000063489 d0 
      FSE   ( 65, -3, 3, 3) =         0.0519077234 d0 ! #
      errtot( 65, -3, 3, 3) =         0.0000008953 d0 
      FSE   ( 65, -3, 3, 4) =         0.0487329643 d0 ! #
      errtot( 65, -3, 3, 4) =         0.0000017032 d0 
      FSE   ( 65, -3, 3, 5) =         0.0461805464 d0 ! #
      errtot( 65, -3, 3, 5) =         0.0000005649 d0 
      FSE   ( 65, -3, 4, 4) =         0.0559629481 d0 ! #
      errtot( 65, -3, 4, 4) =         0.0000004511 d0 
      FSE   ( 65, -3, 4, 5) =         0.0557048871 d0 ! #
      errtot( 65, -3, 4, 5) =         0.0000058711 d0 
      FSE   ( 65, -3, 5, 5) =         0.0580532481 d0 ! #
      errtot( 65, -3, 5, 5) =         0.0000063221 d0 
      FSE   ( 70, -3, 3, 3) =         0.0535705997 d0 ! #
      errtot( 70, -3, 3, 3) =         0.0000009324 d0 
      FSE   ( 70, -3, 3, 4) =         0.0504147510 d0 ! #
      errtot( 70, -3, 3, 4) =         0.0000003920 d0 
      FSE   ( 70, -3, 3, 5) =         0.0478274516 d0 ! #
      errtot( 70, -3, 3, 5) =         0.0000004471 d0 
      FSE   ( 70, -3, 4, 4) =         0.0579811623 d0 ! #
      errtot( 70, -3, 4, 4) =         0.0000004718 d0 
      FSE   ( 70, -3, 4, 5) =         0.0577623984 d0 ! #
      errtot( 70, -3, 4, 5) =         0.0000040286 d0 
      FSE   ( 70, -3, 5, 5) =         0.0602322030 d0 ! #
      errtot( 70, -3, 5, 5) =         0.0000062225 d0 
      FSE   ( 75, -3, 3, 3) =         0.0553332229 d0 ! #
      errtot( 75, -3, 3, 3) =         0.0000009645 d0 
      FSE   ( 75, -3, 3, 4) =         0.0522054522 d0 ! #
      errtot( 75, -3, 3, 4) =         0.0000011986 d0 
      FSE   ( 75, -3, 3, 5) =         0.0495798771 d0 ! #
      errtot( 75, -3, 3, 5) =         0.0000007848 d0 
      FSE   ( 75, -3, 4, 4) =         0.0601328064 d0 ! #
      errtot( 75, -3, 4, 4) =         0.0000003916 d0 
      FSE   ( 75, -3, 4, 5) =         0.0599575762 d0 ! #
      errtot( 75, -3, 4, 5) =         0.0000072673 d0 
      FSE   ( 75, -3, 5, 5) =         0.0625599630 d0 ! #
      errtot( 75, -3, 5, 5) =         0.0000045668 d0 
      FSE   ( 80, -3, 3, 3) =         0.0571913711 d0 ! #
      errtot( 80, -3, 3, 3) =         0.0000007936 d0 
      FSE   ( 80, -3, 3, 4) =         0.0540987831 d0 ! #
      errtot( 80, -3, 3, 4) =         0.0000022712 d0 
      FSE   ( 80, -3, 3, 5) =         0.0514332405 d0 ! #
      errtot( 80, -3, 3, 5) =         0.0000001812 d0 
      FSE   ( 80, -3, 4, 4) =         0.0624136856 d0 ! #
      errtot( 80, -3, 4, 4) =         0.0000005945 d0 
      FSE   ( 80, -3, 4, 5) =         0.0622862108 d0 ! #
      errtot( 80, -3, 4, 5) =         0.0000087302 d0 
      FSE   ( 80, -3, 5, 5) =         0.0650314634 d0 ! #
      errtot( 80, -3, 5, 5) =         0.0000038898 d0 
      FSE   ( 85, -3, 3, 3) =         0.0591392161 d0 ! #
      errtot( 85, -3, 3, 3) =         0.0000008889 d0 
      FSE   ( 85, -3, 3, 4) =         0.0560884432 d0 ! #
      errtot( 85, -3, 3, 4) =         0.0000008988 d0 
      FSE   ( 85, -3, 3, 5) =         0.0533820187 d0 ! #
      errtot( 85, -3, 3, 5) =         0.0000025557 d0 
      FSE   ( 85, -3, 4, 4) =         0.0648172620 d0 ! #
      errtot( 85, -3, 4, 4) =         0.0000007386 d0 
      FSE   ( 85, -3, 4, 5) =         0.0647447117 d0 ! #
      errtot( 85, -3, 4, 5) =         0.0000113573 d0 
      FSE   ( 85, -3, 5, 5) =         0.0676412281 d0 ! #
      errtot( 85, -3, 5, 5) =         0.0000037392 d0 
      FSE   ( 90, -3, 3, 3) =         0.0611702183 d0 ! #
      errtot( 90, -3, 3, 3) =         0.0000006708 d0 
      FSE   ( 90, -3, 3, 4) =         0.0581673791 d0 ! #
      errtot( 90, -3, 3, 4) =         0.0000001621 d0 
      FSE   ( 90, -3, 3, 5) =         0.0554191686 d0 ! #
      errtot( 90, -3, 3, 5) =         0.0000023569 d0 
      FSE   ( 90, -3, 4, 4) =         0.0673359234 d0 ! #
      errtot( 90, -3, 4, 4) =         0.0000006335 d0 
      FSE   ( 90, -3, 4, 5) =         0.0673258575 d0 ! #
      errtot( 90, -3, 4, 5) =         0.0000089432 d0 
      FSE   ( 90, -3, 5, 5) =         0.0703808374 d0 ! #
      errtot( 90, -3, 5, 5) =         0.0000031232 d0 
      FSE   ( 95, -3, 3, 3) =         0.0632756882 d0 ! #
      errtot( 95, -3, 3, 3) =         0.0000005925 d0 
      FSE   ( 95, -3, 3, 4) =         0.0603261584 d0 ! #
      errtot( 95, -3, 3, 4) =         0.0000009832 d0 
      FSE   ( 95, -3, 3, 5) =         0.0575330868 d0 ! #
      errtot( 95, -3, 3, 5) =         0.0000012000 d0 
      FSE   ( 95, -3, 4, 4) =         0.0699593651 d0 ! #
      errtot( 95, -3, 4, 4) =         0.0000015628 d0 
      FSE   ( 95, -3, 4, 5) =         0.0700169138 d0 ! #
      errtot( 95, -3, 4, 5) =         0.0000147459 d0 
      FSE   ( 95, -3, 5, 5) =         0.0732383471 d0 ! #
      errtot( 95, -3, 5, 5) =         0.0000030050 d0 
      FSE   (100, -3, 3, 3) =         0.0654458661 d0 ! #
      errtot(100, -3, 3, 3) =         0.0000005306 d0 
      FSE   (100, -3, 3, 4) =         0.0625542461 d0 ! #
      errtot(100, -3, 3, 4) =         0.0000038384 d0 
      FSE   (100, -3, 3, 5) =         0.0597122507 d0 ! #
      errtot(100, -3, 3, 5) =         0.0000026797 d0 
      FSE   (100, -3, 4, 4) =         0.0726751556 d0 ! #
      errtot(100, -3, 4, 4) =         0.0000029972 d0 
      FSE   (100, -3, 4, 5) =         0.0728009040 d0 ! #
      errtot(100, -3, 4, 5) =         0.0000103797 d0 
      FSE   (100, -3, 5, 5) =         0.0761976114 d0 ! #
      errtot(100, -3, 5, 5) =         0.0000065757 d0 
      FSE   (105, -3, 3, 3) =         0.0676693445 d0 ! #
      errtot(105, -3, 3, 3) =         0.0000005064 d0 
      FSE   (105, -3, 3, 4) =         0.0648382477 d0 ! #
      errtot(105, -3, 3, 4) =         0.0000030810 d0 
      FSE   (105, -3, 3, 5) =         0.0619425181 d0 ! #
      errtot(105, -3, 3, 5) =         0.0000006331 d0 
      FSE   (105, -3, 4, 4) =         0.0754780746 d0 ! #
      errtot(105, -3, 4, 4) =         0.0000050689 d0 
      FSE   (105, -3, 4, 5) =         0.0756683991 d0 ! #
      errtot(105, -3, 4, 5) =         0.0000071640 d0 
      FSE   (105, -3, 5, 5) =         0.0792453310 d0 ! #
      errtot(105, -3, 5, 5) =         0.0000151773 d0 
      FSE   (110, -3, 3, 3) =         0.0699324264 d0 ! #
      errtot(110, -3, 3, 3) =         0.0000009595 d0 
      FSE   (110, -3, 3, 4) =         0.0671634257 d0 ! #
      errtot(110, -3, 3, 4) =         0.0000029982 d0 
      FSE   (110, -3, 3, 5) =         0.0642093856 d0 ! #
      errtot(110, -3, 3, 5) =         0.0000022965 d0 
      FSE   (110, -3, 4, 4) =         0.0783302409 d0 ! #
      errtot(110, -3, 4, 4) =         0.0000050259 d0 
      FSE   (110, -3, 4, 5) =         0.0786011659 d0 ! #
      errtot(110, -3, 4, 5) =         0.0000053496 d0 
      FSE   (110, -3, 5, 5) =         0.0823608883 d0 ! #
      errtot(110, -3, 5, 5) =         0.0000301164 d0 
      FSE   (115, -3, 3, 3) =         0.0722233420 d0 ! #
      errtot(115, -3, 3, 3) =         0.0000017132 d0 
      FSE   (115, -3, 3, 4) =         0.0695117018 d0 ! #
      errtot(115, -3, 3, 4) =         0.0000023485 d0 
      FSE   (115, -3, 3, 5) =         0.0664897140 d0 ! #
      errtot(115, -3, 3, 5) =         0.0000027397 d0 
      FSE   (115, -3, 4, 4) =         0.0812226726 d0 ! #
      errtot(115, -3, 4, 4) =         0.0000047085 d0 
      FSE   (115, -3, 4, 5) =         0.0815643827 d0 ! #
      errtot(115, -3, 4, 5) =         0.0000185802 d0 
      FSE   (115, -3, 5, 5) =         0.0855092813 d0 ! #
      errtot(115, -3, 5, 5) =         0.0000551350 d0 
      FSE   (120, -3, 3, 3) =         0.0745275901 d0 ! #
      errtot(120, -3, 3, 3) =         0.0000019079 d0 
      FSE   (120, -3, 3, 4) =         0.0718664331 d0 ! #
      errtot(120, -3, 3, 4) =         0.0000072836 d0 
      FSE   (120, -3, 3, 5) =         0.0687710355 d0 ! #
      errtot(120, -3, 3, 5) =         0.0000079537 d0 
      FSE   (120, -3, 4, 4) =         0.0841346914 d0 ! #
      errtot(120, -3, 4, 4) =         0.0000095044 d0 
      FSE   (120, -3, 4, 5) =         0.0845558147 d0 ! #
      errtot(120, -3, 4, 5) =         0.0000361896 d0 
      FSE   (120, -3, 5, 5) =         0.0886748746 d0 ! #
      errtot(120, -3, 5, 5) =         0.0001055026 d0 
*
*
      return
      end
c       =================================================
        subroutine grid(maxii,r,v)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       This subroutine initializes Semi-logarithmic
c       radial grid, using change of variables
c       rho(r)=al*r+bt*ln(r). The variable rho(r) is given
c       in the uninform grid rho_i=rho_1+h*(i-1), where h
c       is the grid spacing, rho_1=al*r_1+bt*ln(r_1) and
c       r_1 > 0 is the first grid point closest to the
c       point of origin.
c       ii   .... number of points (ii=maxii-20), where
c                 maxii is given in the include file
c                 'qedmod.inc'
c       r(i) .... non-unform radial grid (i=1,ii)
c       z    .... Nuclear charge
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /z/z
        common /ii/ii/r1/r1/r2/r2/h/h/al/al/bt/bt
        common /knucl/knucl/inucl/inucl/rnucl/rnucl
        real*8 r(maxii),v(maxii)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        r1=dexp(-6.d0)*(256.d0/maxii)**2/5.d0
        if (z.gt.0.99) then
          r1=r1/z
        else
          r1=r1/100.d0
        endif
        if (rnucl.gt.0.and.r1.gt.0.10d0*rnucl) r1=0.10d0*rnucl
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        rmax=dabs(r2)
        imax=maxii-20
        imax=((imax+1)/2)*2-1
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Semi-logarithmic grid: ro=al*r+bt*ln(r)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        ii=imax
        bt=1.d0
        h=bt*dlog(r2/r1)/((ii-1)*0.67)
        al=(h*(imax-1)-bt*dlog(r2/r1))/(r2-r1)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (al.ge.0.d0) goto 210
        ii=(bt*dlog(r2/r1))/h+1
        write( *,99) imax,ii,h,al,bt,r1,r2
99      format(/'  ii > Imax'/'  Imax=',i6/'  ii  =',i5/
     1  '  h   =',e13.5,/'  al  =',e13.5/,'  bt  =',e13.5/
     2  '  r1  =',e13.5,'  r2  =',e13.5)
        call exit1
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Nuclear correction
c       - - - - - - - - - - - - - - - - - - - - - - - - -
210     inucl=0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (knucl.eq.1) then
          ro1=al*r1+bt*dlog(r1)
          ron=al*rnucl+bt*dlog(rnucl)
          d=ron-ro1
          inucl=(d/h+0.0001d0)+1
          if (inucl.lt.11) then
            inucl=11
            dr=d-(inucl-1)*h
            ro1=ro1+dr
            r0=dexp(ro1)/bt
            r1=rr(ro1,r0)
          endif
          dr=(rnucl-r1)/(r2-r1)
          h=(bt*dlog(rnucl/r1)-bt*dlog(r2/r1)*dr)/
     &    ((inucl-1)-(ii-1)*dr)
          al=((inucl-1)*h-bt*dlog(rnucl/r1))/(rnucl-r1)
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do i=1,maxii
          r(i)=0.d0
          v(i)=0.d0
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        call tbr(ii,r1,h,al,bt,r,v)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        write( *,5) ii,r1,r2,h,al,bt
        write(11,5) ii,r1,r2,h,al,bt
5       format (/2x,'Semi-logarithmic radial grid:',
     1  /2x,'ii =',i5,/2x,'r1 =',e13.5,
     2  /2x,'r2 =',e13.5,/2x,'h  =',e13.5,
     3  /2x,'al =',e13.5,/2x,'bt =',e13.5)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
1000    return
        end
c       =================================================
        subroutine tbr(ii,r1,h,al,bt,r,v)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        dimension r(*),v(*)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        r(1)=r1
        d=r1
        t=al*d+bt
        v(1)=d/t
        p=al*d+bt*dlog(d)
        do 10 i=2,ii
        p=p+h
200     t=al*d+bt*dlog(d)
        t=(p-t)/(al*d+bt)
        d=d*(1.d0+t)
        if (dabs(t).gt.0.5d-11) goto 200
        t=al*d+bt
        v(i)=d/t
        r(i)=d
10      continue
        return
        end
c       =================================================
        function rr(ro,r0)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /al/al/bt/bt
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        r=r0
200     t=al*r+bt*dlog(r)
        t=(ro-t)/(al*r+bt)
        if (t.le.-1.d0) t=-0.5d0
        r=r*(1.d0+t)
        if (dabs(t).gt.0.1d-8) goto 200
        rr=r
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
c       =================================================
        subroutine nucl(maxii,r,unuc,dnuc)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       This subroutine computes nonpoint part of the
c       nuclear potential.
c
c       unuc(i) .... Nonpoint part of the nuclear potential
c                    multiplied by r(i).
c       cl     ..... Speed of light (in atomic units)
c       z      ..... Nuclear charge
c       knucl  ..... Nuclear model (0 - point, 1 - uniform sphere,
c                    2 - Fermi model, 3 - Gauss distribution)
c       ii     ..... Number of grid points
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /z/z/cl/cl/ii/ii
        real*8 r(maxii),unuc(maxii),dnuc(maxii)
        common /inucl/inucl/knucl/knucl/rnucl/rnucl
     1  /vnuc/vnuc(20)
        common /nmax/nmax
        common /fermi/aa,tt,cc,ro0,rnucl0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        parameter (pi=3.1415926535897932385d0)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        rnucl0=rnucl/dsqrt(5.d0/3.d0)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do n=0,nmax
          i=n+1
          vnuc(i)=0.d0
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Point nuclear
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (knucl.eq.0) then
c       - - - - - - - - - - - - - - - - - - - - - - - - -
          write( *,5)
          write(11,5)
5         format(/2x,'Point Nuclear Charge'/2x,20('*'))
c       - - - - - - - - - - - - - - - - - - - - - - - - -
          vnuc(1)=-z
          goto 1000
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        fermi2at=1.d-5/0.529177249d0
        rfermi=rnucl/fermi2at
        rms_au=rnucl/dsqrt(5.d0/3.d0)
        rms_fermi=rfermi/dsqrt(5.d0/3.d0)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Volume distribution
c       - - - - - - - - - - - - - - - - - - - - - - - - -
200     if (knucl.eq.1) then
c       - - - - - - - - - - - - - - - - - - - - - - - - -
          write( *,15)
          write(11,15)
15        format(/2x,'Volume Distribution of the Nuclear Charge'
     &    /2x,41('*'))
          write( *,25) rnucl,rfermi,rms_au,rms_fermi
          write(11,25) rnucl,rfermi,rms_au,rms_fermi
25        format 
     &    ( 2x,'Rnucl     =',e13.6,1x,'[au] =',f8.4,1x,'[Fermi]',
     &     /2x,'Sqrt<R^2> =',e13.6,1x,'[au] =',f8.4,1x,'[Fermi]')
c
          vnuc(2)=-1.5d0*z/rnucl
          vnuc(4)= 0.5d0*z/rnucl**3
          goto 1000
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (knucl.eq.2) then
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Fermi distribution
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       1 fermi = 10**(-13) cm
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        fermi2at=1.d-5/0.529177249d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        write( *,35)
        write(11,35)
35      format(/2x,'Fermi Distribution of the Nuclear Charge',
     &    /2x,40('*'))
c
        write( *,25) rnucl,rfermi,rms_au,rms_fermi
        write(11,25) rnucl,rfermi,rms_au,rms_fermi
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Setting nuclear parameters
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        tt=2.30d0
        aa=tt/(4.d0*dlog(3.d0))*fermi2at
        cc=5.d0/3.d0*rnucl0**2-7.d0/3.d0*pi**2*aa**2
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (cc.le.0.d0) then
          write( *,'(/2x,a,1x,a/2x,a)') '*** Warning!!!',
     &    'Z too small. Fermi Distribution does not work ***.',
     &    'The program use Sphere Ditribution.' 
          write(11,'(/2x,a,1x,a/2x,a)') '*** Warning!!!',
     &    'Z too small. Fermi Distribution does not work ***.',
     &    'The program use Sphere Ditribution.'
          knucl=1
          goto 200
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        cc=dsqrt(cc)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        write( *,45) cc,cc/fermi2at,tt*fermi2at,tt
        write(11,45) cc,cc/fermi2at,tt*fermi2at,tt
45      format (/2x,'Parametr c=',e13.6,1x,'[a.u.] =',f8.4,
     &  1x,'[Fermi]',/2x,'Parametr t=',e13.6,1x,
     &  '[a.u.] =',f8.4,1x,'[Fermi]')
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        dnorm=cc**3*(1+(pi*aa/cc)**2-6*(aa/cc)**3*sk(3,-cc/aa))/3
        ro0=z/dnorm
        write( *,55) ro0,ro0*fermi2at**3
        write(11,55) ro0,ro0*fermi2at**3
55      format (2x,'Ro_nucl(0)=',e13.6,1x,'[a.u.] =',f8.4,1x,
     1  '[Fermi]')
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Calculation of the Tailor coefficients
c       - - - - - - - - - - - - - - - - - - - - - - - - -
	B0=ro0/(1+dexp(-cc/aa))
	B1=-ro0/((1+dexp(-cc/aa))**2)*dexp(-cc/aa)/aa
        B2=B1/(1+dexp(-cc/aa))*(1.d0-dexp(-cc/aa))/aa
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       RM(n)=<ro(r)*r^n>/ro0;  ro0=2*ro(c)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
	RM1=RM(1,cc,aa)
        vnuc(2)=-RM1*ro0
c        vnuc(4)=-(B0/3.d0-B0/2.d0)
c        vnuc(5)=-(B1/4.d0-B1/3.d0)
        vnuc(4)=B0/6.d0
        vnuc(5)=B1/12.d0
        vnuc(6)=B2/40.d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (knucl.eq.3) then
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Gaussian distribution
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       V(r)=-Z*erf(dsqrt(1.5)*r/Rms)/r
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        write( *,65)
        write(11,65)
65      format(/2x,'Gaussian Distribution of the Nuclear Charge'
     &  /2x,43('*'))
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        dn0=-2*z/dsqrt(pi)
        a0=dsqrt(1.5d0)/rnucl0
        vnuc(2)= dn0*a0
        vnuc(4)=-dn0*a0**3/3.d0
        vnuc(6)= dn0*a0**5/10.d0
        vnuc(8)=-dn0*a0**7/42.d0
        vnuc(10)=dn0*a0**9/216.d0
c       vnuc(12)=-dn0*a0**11/1320.d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (dabs(vnuc(1)).gt.dabs(cl)) then
          write( *,'(/2x,a)') '*** Error: Vnuc(1) > |Cl| ***'
          write(11,'(/2x,a)') '*** Error: Vnuc(1) > |Cl| ***'
        call exit1
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
1000    do i=1,maxii
          unuc(i)=0.d0
          dnuc(i)=0.d0
        enddo
        if (knucl.eq.0) return
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (knucl.eq.1) then
c
          if (inucl.lt.1) return
c
          do i=1,inucl
            vi=0.d0
            do n=0,nmax
              vi=vi+vnuc(n+1)*r(i)**n
            enddo
            unuc(i)=vi+z
          enddo
c
          do i=1,inucl
            dvi=0.d0
            do n=0,nmax-1
              dvi=dvi+n*vnuc(n+1)*r(i)**n
            enddo
c           dvi1=-1.5d0*z/r(inucl)*(1.d0-r(i)**2/r(inucl)**2)
            dnuc(i)=dvi
          enddo
c
          goto 1100
c
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Fermi and Gaussian distributions        
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do i=1,ii
          if (knucl.eq.2) vi=vpot_fermi(z,cc,aa,r(i))
          if (knucl.eq.3) vi=vpot_gauss(z,rnucl0,r(i))
          unuc(i)=z+vi
        enddo
        do i=1,ii
          if (knucl.eq.2) dnuc(i)=dV_dR(z,cc,aa,r(i))*r(i)
          if (knucl.eq.3) dnuc(i)=dV_dR_gauss(z,rnucl0,r(i))*r(i)
        enddo
        do i=ii,1,-1
          if (dabs(unuc(i))+dabs(dnuc(i)).gt.1.d-13) goto 1100
          unuc(i)=0.d0
          dnuc(i)=0.d0
          inucl=i
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
 1100   unuc(ii+3)=inucl
        unuc(ii+4)=0.d0
        unuc(ii+5)=z+vnuc(1)
        do m=1,nmax-1
          i=ii+5+m
          unuc(i)=vnuc(m+1)*r(1)**m
        enddo
        unuc(ii+5+nmax)=0.d0
c
        dnuc(ii+3)=inucl
        dnuc(ii+4)=0.d0
        do m=0,nmax
          i=ii+5+m
          dnuc(i)=m*unuc(i)
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
c       =================================================
        real*8 function ro_nucl(r)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8(a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /knucl/knucl/rnucl/rnucl/z/z
        common /fermi/aa,tt,cc,ro0,rnucl0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        parameter (pi=3.1415926535897932385d0)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Volume distribution
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (knucl.eq.1) then
          if(r.lt.rnucl+1.d-10)then
            ro_nucl=3.d0*z/rnucl**3
          else
            ro_nucl=0.d0
          endif
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Fermi distribution
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (knucl.eq.2) then
          x=(r-cc)/aa
          if (x.le.200) then 
            ro_nucl=ro0/(1.d0+dexp(x))
          else
            ro_nucl=0.d0
          endif
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Gaussian distribution
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (knucl.eq.3) then
          ro_nucl=z*4*pi*(1.5d0/(pi*rnucl0**2))**1.5d0*
     &             dexp(-1.5d0*(r/rnucl0)**2)
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
c       =================================================              
        function vpot_fermi(Z,C,A,R)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (A-H,O-Z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        parameter (pi=3.1415926535897932385d0)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        dn=1+(pi*a/c)**2-6d0*((a/c)**3)*SK(3,-c/a)
c       - - - - - - - - - - - - - - - - - - - - - - - - -      
        IF(R.LT.C) THEN
        B1=-6d0*((a/c)**3)*SK(3,-c/a)+6d0*((a/c)**3)*SK(3,(r-c)/a)
        B2=r/c*(3d0/2d0+((pi*a/c)**2)/2d0-3d0*((a/c)**2)*SK(2,(r-c)/a))
        B3=-((r/c)**3)/2d0
        V=Z/dn*(B1+B2+B3)
c       - - - - - - - - - - - - - - - - - - - - - - - - -      
        ELSE
c       - - - - - - - - - - - - - - - - - - - - - - - - -      
        B1=dn
        B2=3d0*((a/c)**2)*(r/c*SK(2,(c-r)/a)+2d0*a/c*SK(3,(c-r)/a))
        V=Z/dn*(B1+B2)
        END IF
c       - - - - - - - - - - - - - - - - - - - - - - - - -      	
        vpot_fermi=-V
c       - - - - - - - - - - - - - - - - - - - - - - - - -      
	RETURN
        END
c       =================================================              
        function vpot_gauss(z,rms,r)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (A-H,O-Z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        external erf
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        vpot_gauss=-z*erf(dsqrt(1.5d0)*r/rms)
c       - - - - - - - - - - - - - - - - - - - - - - - - -      
	return
        end
c       =================================================
        FUNCTION dV_dR(Z,C,A,R)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (A-H,O-Z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        parameter (pi=3.1415926535897932385d0)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        dn=1+(pi*a/c)**2-6d0*((a/c)**3)*SK(3,-c/a)
c       - - - - - - - - - - - - - - - - - - - - - - - - -      
        IF(R.LT.C) THEN
          B1=6d0*((a/c)**3)*SK(2,(r-c)/a)/a
          B2=1d0/c*(3d0/2d0+((pi*a/c)**2)/2d0-3d0*((a/c)**2)*
     &    SK(2,(r-c)/a))
          B3=r/c*(-3d0*((a/c)**2)*SK(1,(r-c)/a))/a
          B4=-((r/c)**3)/2d0
	  B4=-3d0/2d0/(c**3)*(r**2)
          V=Z/dn*(B1+B2+B3+B4)
c       - - - - - - - - - - - - - - - - - - - - - - - - -      
        ELSE
c       - - - - - - - - - - - - - - - - - - - - - - - - -      
          B1=3d0*((a/c)**2)*(1d0/c*SK(2,(c-r)/a))
          B2=-3d0*((a/c)**2)*(r/c*SK(1,(c-r)/a)+2d0*a/c*
     &    SK(2,(c-r)/a))/a
          V=Z/dn*(B1+B2)
        END IF
c       - - - - - - - - - - - - - - - - - - - - - - - - -      	
        dV_dR=-V
c       - - - - - - - - - - - - - - - - - - - - - - - - -      
	RETURN
        END
c       =================================================              
        function dV_dR_gauss(z,rms,r)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (A-H,O-Z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        parameter (pi=3.1415926535897932385d0)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        dv=-z*2/dsqrt(pi)*dexp(-1.5d0*(r/rms)**2)
        dV_dR_gauss=dsqrt(1.5d0)/rms*dv
c       - - - - - - - - - - - - - - - - - - - - - - - - -      
	return
        end
c       =================================================
        function sk(k,x)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        parameter (error=1d-20)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        dk=K
        N=0
        sk=0d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
10      N=N+1
        dn=N
        add=DEXP(dN*X)/(dN**dK)*(-1)**N
	if(N.eq.1) add1=add
        sk=sk+add
        if (add.eq.0d0) goto 1000
        if (dabs(add/add1).lt.error) goto 1000
        goto 10
c       - - - - - - - - - - - - - - - - - - - - - - - - -      
1000    return
        end
c       =================================================
        FUNCTION RM(M,C,A)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (A-H,O-Z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        parameter (error=1.d-4)
        parameter (pi=3.1415926535897932385d0)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (m.eq.1) then
        rm=c*c/2+a**2*sk(2,-c/a)+2*a**2*(pi**2/12)
        return
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        RM=c**(M+1)/(M+1)
        RM1=c**(M+1)/(M+1)
        N=0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c        write(*,*)
c        write(*,5551)'N,ERR1,add=                                      '
c5551    format(1X,A,$)
c        ibacsp = 8
c       - - - - - - - - - - - - - - - - - - - - - - - - -
10      N=N+1
        dn=N
ccc      add=((-1)**N)*DEXP(-dN*C/A)*(fxneax(M,dN/A,C)-fxneax(M,dN/A,0d0))
       add=((-1)**N)*(fxneax(M,dN/A,C)-DEXP(-dN*C/A)*fxneax(M,dN/A,0d0))
        RM=RM+add
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c        write(*,5552) (ibacsp,isp=1,39) ,N,dabs(add/RM1),add
c5552    format(39A1,I11,F14.11,F14.10,$)
c       - - - - - - - - - - - - - - - - - - - - - - - - -	 
        if(dabs(add/RM1).lt.error) goto 20
        GOTO 10
c       - - - - - - - - - - - - - - - - - - - - - - - - -       
20      N=0
c        write(*,*)
c        write(*,5551)'N,ERR2,add=                                      '
        addd=add
c       - - - - - - - - - - - - - - - - - - - - - - - - -             
30      N=N+1
        dN=N
ccc      add=((-1)**(N+1))*DEXP(dN*C/A)*(-fxneax(M,-dN/A,C))
        add=((-1)**(N+1))*(-fxneax(M,-dN/A,C))
        IF(N.eq.1) RM1=add
        RM=RM+add
c       - - - - - - - - - - - - - - - - - - - - - - - - -             
c        write(*,5552) (ibacsp,isp=1,39),N,dabs(add/RM1),add
c       - - - - - - - - - - - - - - - - - - - - - - - - -             
        if(dabs(add).lt.dabs(addd)) goto 1000
        GOTO 30
c       - - - - - - - - - - - - - - - - - - - - - - - - -
1000    RETURN
        END
c       =================================================
        FUNCTION fxneax(N,A,X)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (A-H,O-Z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        SUM=0d0
        DO 10 l=0,N
c       SUM=SUM+dfac(n)/dfac(n-l)*((A*X)**(n-l))*(-1)**l
        SUM=SUM+dfac(n)/dfac(n-l)*(STPDI(A*X,n-l))*(-1)**l
10      CONTINUE
c       - - - - - - - - - - - - - - - - - - - - - - - - -
ccc      fxneax=SUM*DEXP(A*X)/(A**(N+1))
        fxneax=SUM/(A**(N+1))
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        RETURN
        END
c       =================================================
        function dfac(n)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit double precision (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        parameter (one=1.0D+00)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if(n.eq.0) then
        dfac=one
	return
        end if
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        df= one
        ns=n
        do 5, i = ns , 1, -1
        df=df * i
   5    continue
        dfac=df
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
c       =================================================       
        FUNCTION STPDI(X,I)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        IMPLICIT REAL*8 (A-H,O-Z)      
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        IF(I.EQ.0) THEN
        stpdi=1.d0
        goto 1000
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        ELSE
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        stpdi=x**i
        goto 1000
        END IF 
c       - - - - - - - - - - - - - - - - - - - - - - - - -
1000    return
        END
c       =================================================              
        DOUBLE PRECISION FUNCTION erf(x)
c       - - - - - - - - - - - - - - - - - - - - - - - - -      
C             EVALUATION OF THE REAL ERROR FUNCTION
C       .. Scalar Arguments ..
        DOUBLE PRECISION x
C       ..
C       .. Local Scalars ..
        DOUBLE PRECISION ax,bot,c,t,top,x2
C       ..
C       .. Local Arrays ..
        DOUBLE PRECISION a(5),b(3),p(8),q(8),r(5),s(4)
C       ..
C       .. Intrinsic Functions ..
        INTRINSIC abs,exp,sign
c       - - - - - - - - - - - - - - - - - - - - - - - - -      
C       Data statements ..
c       - - - - - - - - - - - - - - - - - - - - - - - - -      
        DATA c/.564189583547756D0/
        DATA a(1)/.771058495001320D-04/,a(2)/-.133733772997339D-02/,
     +       a(3)/.323076579225834D-01/,a(4)/.479137145607681D-01/,
     +       a(5)/.128379167095513D+00/
        DATA b(1)/.301048631703895D-02/,b(2)/.538971687740286D-01/,
     +       b(3)/.375795757275549D+00/
        DATA p(1)/-1.36864857382717D-07/,p(2)/5.64195517478974D-01/,
     +       p(3)/7.21175825088309D+00/,p(4)/4.31622272220567D+01/,
     +       p(5)/1.52989285046940D+02/,p(6)/3.39320816734344D+02/,
     +       p(7)/4.51918953711873D+02/,p(8)/3.00459261020162D+02/
        DATA q(1)/1.00000000000000D+00/,q(2)/1.27827273196294D+01/,
     +       q(3)/7.70001529352295D+01/,q(4)/2.77585444743988D+02/,
     +       q(5)/6.38980264465631D+02/,q(6)/9.31354094850610D+02/,
     +       q(7)/7.90950925327898D+02/,q(8)/3.00459260956983D+02/
        DATA r(1)/2.10144126479064D+00/,r(2)/2.62370141675169D+01/,
     +       r(3)/2.13688200555087D+01/,r(4)/4.65807828718470D+00/,
     +       r(5)/2.82094791773523D-01/
        DATA s(1)/9.41537750555460D+01/,s(2)/1.87114811799590D+02/,
     +       s(3)/9.90191814623914D+01/,s(4)/1.80124575948747D+01/
c       - - - - - - - - - - - - - - - - - - - - - - - - -      
C       Executable Statements ..
c       - - - - - - - - - - - - - - - - - - - - - - - - -      
        ax = abs(x)
        IF (ax.GT.0.5D0) GO TO 10
        t = x*x
        top = ((((a(1)*t+a(2))*t+a(3))*t+a(4))*t+a(5)) + 1.0D0
        bot = ((b(1)*t+b(2))*t+b(3))*t + 1.0D0
        erf = x* (top/bot)
        RETURN
c       - - - - - - - - - - - - - - - - - - - - - - - - -      
   10   IF (ax.GT.4.0D0) GO TO 20
        top = ((((((p(1)*ax+p(2))*ax+p(3))*ax+p(4))*ax+p(5))*ax+
     +        p(6))*ax+p(7))*ax + p(8)
        bot = ((((((q(1)*ax+q(2))*ax+q(3))*ax+q(4))*ax+q(5))*ax+
     +        q(6))*ax+q(7))*ax + q(8)
        erf = 0.5D0 + (0.5D0-exp(-x*x)*top/bot)
        IF (x.LT.0.0D0) erf = -erf
        RETURN
c       - - - - - - - - - - - - - - - - - - - - - - - - -      
   20   IF (ax.GE.5.8D0) GO TO 30
        x2 = x*x
        t = 1.0D0/x2
        top = (((r(1)*t+r(2))*t+r(3))*t+r(4))*t + r(5)
        bot = (((s(1)*t+s(2))*t+s(3))*t+s(4))*t + 1.0D0
        erf = (c-top/ (x2*bot))/ax
        erf = 0.5D0 + (0.5D0-exp(-x2)*erf)
        IF (x.LT.0.0D0) erf = -erf
        RETURN
c       - - - - - - - - - - - - - - - - - - - - - - - - -      
   30   erf = sign(1.0D0,x)
        RETURN
        END
c       =================================================
        subroutine origin_coef(ii,maxii,f,r)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /nmax/nmax
        real*8 f(maxii),r(maxii)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        ih=maxii/512
        if ((ih/2)*2.ne.ih) ih=ih+1
        do m=0,nmax
          i=ii+5+m
          f(i)=0.d0
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        gam=f(ii+4)
        i1=1
        i2=i1+ih
        i3=i2+ih
        x1=r(i1)
        x2=r(i2)
        x3=r(i3)
        f1=f(i1)/r(i1)**gam
        f2=f(i2)/r(i2)**gam
        f3=f(i3)/r(i3)**gam
        h1=x2-x1
        h2=x3-x2
        dh=1.d0/(h1*h2*(h1+h2))
        p0=dh*(h2*x2*x3*f1-(h1+h2)*x1*x3*f2+h1*x1*x2*f3)
        p1=dh*(-h2*(x2+x3)*f1+(h1+h2)*(x1+x3)*f2-h1*(x1+x2)*f3)
        p2=dh*2*(h2*f1-(h1+h2)*f2+h1*f3)
        f(ii+5)=p0
        f(ii+6)=p1*r(1)
        f(ii+7)=0.5d0*p2*r(1)**2
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        dint1=0.d0
        dint2=0.d0
        do i=1,21,10*ih
          fi=0.d0
          do m=0,nmax
            j=m+ii+5
            fi=fi+f(j)*(r(i)/r(1))**m
          enddo
          fi=fi*r(i)**gam
          dint1=dint1+f(i)*(h1+h2)/2
          dint2=dint2+fi*(h1+h2)/2
        enddo
c        write( *,'(a,2e16.8)') 'Test origin:',dint1,dint2
c        write(11,'(a,2e16.8)') 'Test origin:',dint1,dint2
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
c       =================================================
        subroutine read_func(ni,v1,v2,nrec,maxii,ff)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        real*8 ff(*)
        integer ni,nrec
        real*8 v1(maxii),v2(maxii)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        nr1=2*ni-1
        nr2=nr1+1
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do i=1,maxii
          iff=maxii*(nr1-1)+i
c         v1(i)=ff(i,nr1)
          v1(i)=ff(iff)
        enddo
        if (nrec.eq.1) goto 1000
        do i=1,maxii
          iff=maxii*(nr2-1)+i
c         v2(i)=ff(i,nr2)
          v2(i)=ff(iff)
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
1000    return
        end
c       =================================================
        function tint(k,p,q,a,b,r,v)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /ii/ii/h/h/al/al/bt/bt/nmax/nmax
        dimension r(*),v(*)
        dimension p(*),q(*),a(*),b(*)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        i0=1
        imax=p(ii+3)+0.01d0
        r0=r(i0)
        v0=v(i0)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        g=p(ii+4)+a(ii+4)+k
        t0=0.d0
        p0=0.d0
        do m=0,nmax
          i=ii+5+m
          d=0.d0
          do n=0,m
            j=ii+5+n
            d=d+p(j)*a(i-n)+q(j)*b(i-n)
          enddo
          t0=t0+d/(g+m+1)
          p0=p0+d*(g+m)
        enddo
        t0=t0*r0**(g+1)
        p0=p0*r0**(g-1)
        f0=(p(i0)*a(i0)+q(i0)*b(i0))*v0*r0**k
        p0=p0*v0*v0+f0*bt/(al*r0+bt)**2
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        dt=0.d0
        do i=i0,imax
          fi=(p(i)*a(i)+q(i)*b(i))*v(i)
          if (k.ne.0) fi=fi*r(i)**k
          dt=dt+fi
        enddo
        dt=dt*h
        dt=t0+dt-h*(0.5d0*f0-h/12.d0*p0)-0.5d0*h*fi
c       - - - - - - - - - - - - - - - - - - - - - - - - -
1000    tint=dt
        return
        end
c       =================================================
        function tint_sympson(k,p,q,a,b,ii,h,r,v)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       This subroutine computes the matrix elements of r^k
c       tint= <p,q|r^k|a,b>
c       Numerical Integration is performed using Simpson's rule
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        real*8 p(*),q(*),a(*),b(*),r(*),v(*)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        i0=1
        r0=r(i0)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        imax=ii
        do i=ii,1,-1
          d=dabs(p(i))+dabs(q(i))+dabs(a(i))+dabs(b(i))
          if (d.gt.1.d-20) goto 200
          imax=i
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
 200    g=p(ii+4)+a(ii+4)+k
        ro0=(p(i0)*a(i0)+q(i0)*b(i0))*r0**k
        t0=ro0/(g+1)*r0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        ds=0.d0
        i1=1
        i2=imax
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       If the number of points is even we use four-point
c       Simpson's rule for first 4 points. The number
c       of the rest points is odd
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (2*(imax/2).eq.imax) then
          i=i1
          f1=(p(i)*a(i)+q(i)*b(i))*v(i)*r(i)**k
          i=i1+1
          f2=(p(i)*a(i)+q(i)*b(i))*v(i)*r(i)**k
          i=i1+2
          f3=(p(i)*a(i)+q(i)*b(i))*v(i)*r(i)**k
          i=i1+3
          f4=(p(i)*a(i)+q(i)*b(i))*v(i)*r(i)**k
          ds=ds+3.d0/8.d0*h*(f1+3*f2+3*f3+f4)
          i1=i1+3
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Simpson's rule for the odd number of points.
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        dt1=0.d0
        do i=i1,i2,2
          fi=(p(i)*a(i)+q(i)*b(i))*v(i)*r(i)**k
          dt1=dt1+fi
        enddo
        f1=(p(i1)*a(i1)+q(i1)*b(i1))*v(i1)*r(i1)**k
        f2=(p(i2)*a(i2)+q(i2)*b(i2))*v(i2)*r(i2)**k
        dt1=2*dt1-f1-f2
c
        dt2=0.d0
        do i=i1+1,i2-1,2
          fi=(p(i)*a(i)+q(i)*b(i))*v(i)*r(i)**k
          dt2=dt2+fi
        enddo
        dt2=4*dt2
        dt=h*(dt1+dt2)/3.d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        tint_sympson=ds+dt+t0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
c       =================================================
        subroutine sint(ds,c,p,q,a,b,r,v)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       This subroutine computes the matrix elements of
c       the potential v(r). ds=<p,q|v|a,b>.
c       The array c(r) is the potential v(r) multiplied on r
c       c(r)=v(r)*r. Numerical Integration is performed
c       using Simpson's rule
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Numerical Integration. Simpson's rule.
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /h/h/ii/ii
        real*8 c(*),p(*),q(*),a(*),b(*),r(*),v(*)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        i0=1
        r0=r(i0)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        imax=ii
        do i=ii,1,-1
          d=dabs(p(i))+dabs(q(i))+dabs(a(i))+dabs(b(i))+dabs(c(i))
          if (d.gt.1.d-20) goto 200
          imax=i-1
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
 200    g=p(ii+4)+a(ii+4)+c(ii+4)-1.d0
        ro0=c(i0)*(p(i0)*a(i0)+q(i0)*b(i0))/r(i0)
        t0=ro0/(g+1)*r0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        ds=0.d0
        i1=1
        i2=imax
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       If the number of points is even we use four-point
c       Simpson's rule for first 4 points. The number
c       of the rest points is odd
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (2*(imax/2).eq.imax) then
          i=i1
          f1=(p(i)*a(i)+q(i)*b(i))*c(i)*v(i)/r(i)
          i=i1+1
          f2=(p(i)*a(i)+q(i)*b(i))*c(i)*v(i)/r(i)
          i=i1+2
          f3=(p(i)*a(i)+q(i)*b(i))*c(i)*v(i)/r(i)
          i=i1+3
          f4=(p(i)*a(i)+q(i)*b(i))*c(i)*v(i)/r(i)
          ds=ds+3.d0/8.d0*h*(f1+3*f2+3*f3+f4)
          i1=i1+3
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Simpson's rule for the odd number of points.
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        dt1=0.d0
        do i=i1,i2,2
          fi=(p(i)*a(i)+q(i)*b(i))*c(i)*v(i)/r(i)
          dt1=dt1+fi
        enddo
        f1=(p(i1)*a(i1)+q(i1)*b(i1))*c(i1)*v(i1)/r(i1)
        f2=(p(i2)*a(i2)+q(i2)*b(i2))*c(i2)*v(i2)/r(i2)
        dt1=2*dt1-f1-f2
c
        dt2=0.d0
        do i=i1+1,i2-1,2
          fi=(p(i)*a(i)+q(i)*b(i))*c(i)*v(i)/r(i)
          dt2=dt2+fi
        enddo
        dt2=4*dt2
        dt=h*(dt1+dt2)/3.d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        ds=ds+dt+t0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
c       =================================================
        subroutine write_func(ni,v1,v2,nrec,maxii,ff)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        real*8 ff(*)
        integer ni,nrec
        real*8 v1(maxii),v2(maxii)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        nr1=2*ni-1
        nr2=nr1+1
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do i=1,maxii
          iff=maxii*(nr1-1)+i
c         ff(i,nr1)=v1(i)
          ff(iff)=v1(i)
        enddo
        if (nrec.eq.1) goto 1000
        do i=1,maxii
          iff=maxii*(nr2-1)+i
c         ff(i,nr2)=v2(i)
          ff(iff)=v2(i)
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
1000    return
        end
c       =================================================
        function flagr(x,f,r,ii,h,al,bt)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        real*8 f(*),r(*)
        common /nmax/nmax
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        es  =f(ii+1)
        imax=f(ii+3)+0.01d0
        gs  =f(ii+4)
        if (x.gt.r(imax)) goto 210
        if (x.lt.r(   1)) goto 220
        ro=al*(x-r(1))+bt*dlog(x/r(1))
        i3=3
        if (x.le.r(i3)) goto 200
        i3=imax-3
        if (x.ge.r(i3)) goto 200
        i3=ro/h+1
200     p=(ro-(i3-1)*h)/h
        i2=i3-1
        i1=i2-1
        i4=i3+1
        i5=i4+1
        i6=i5+1
        pm1=p-1.d0
        pm2=p-2.d0
        pm3=p-3.d0
        pp1=p+1.d0
        pp2=p+2.d0
        fi=0.1d0*pm2*pm1*p*pp1*(pp2*f(i6)-pm3*f(i1))
        fi=fi+0.5d0*pm3*pm1*p*pp2*(pm2*f(i2)-pp1*f(i5))
        fi=fi+pm3*pm2*pp1*pp2*(p*f(i4)-pm1*f(i3))
        flagr=fi/12.d0
        return
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Asimptotics at the infinity x > r(imax)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
210     flagr=0.d0
        return
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Asimptotics at the origin pri x < r(1)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
220     fi=0.d0
        do m=0,nmax
          i=ii+5+m
          fi=fi+f(i)*(x/r(1))**m
        enddo
        flagr=x**gs*fi
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
c       =================================================
        function flagr_spl(x,f,r,ii,h,r0,i0)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        real*8 f(*),r(*)
        common /nmax/nmax
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        es  =f(ii+1)
        imax=f(ii+3)+0.01d0
        gs  =f(ii+4)
        if (x.gt.r(imax)) goto 210
        if (x.lt.r(   1)) goto 220
        ro=dlog(x/r0+1.d0)-i0*h
c       ro=al*(x-r(1))+bt*dlog(x/r(1))
        i3=3
        if (x.le.r(i3)) goto 200
        i3=imax-3
        if (x.ge.r(i3)) goto 200
        i3=ro/h+1
200     p=(ro-(i3-1)*h)/h
c        write(*,*) ro,x,r(i3),r(i3+1),i0
        i2=i3-1
        i1=i2-1
        i4=i3+1
        i5=i4+1
        i6=i5+1
        pm1=p-1.d0
        pm2=p-2.d0
        pm3=p-3.d0
        pp1=p+1.d0
        pp2=p+2.d0
        fi=0.1d0*pm2*pm1*p*pp1*(pp2*f(i6)-pm3*f(i1))
        fi=fi+0.5d0*pm3*pm1*p*pp2*(pm2*f(i2)-pp1*f(i5))
        fi=fi+pm3*pm2*pp1*pp2*(p*f(i4)-pm1*f(i3))
        flagr_spl=fi/12.d0
        return
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Asimptotics at the infinity x > r(imax)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
210     flagr_spl=0.d0
        return
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Asimptotics at the origin pri x < r(1)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
220     fi=0.d0
        do m=0,nmax
          i=ii+5+m
          fi=fi+f(i)*(x/r(1))**m
        enddo
        flagr_spl=x**gs*fi
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
c       =================================================
        subroutine dminv(a,n,nmax,d,l,m)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Invertion of arbitrary matrix:  C=A**(-1),
c       using Gauss-Gordan method
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       l,m - work arrays (dimension - n)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        dimension a(nmax,nmax),l(nmax),m(nmax)
c       double precision a,d,biga,hold
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Searching maximal element.
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        d=1.d0
        do 80 k=1,n
        l(k)=k
        m(k)=k
        biga=a(k,k)
        do 20 j=k,n
        do 20 i=k,n
10      if (dabs(biga)-dabs(a(i,j))) 15,20,20
15      biga=a(i,j)
        l(k)=i
        m(k)=j
20      continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Changing the strings
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        j=l(k)
        if (j-k) 35,35,25
25      do 30 i=1,n
        hold=-a(k,i)
        a(k,i)=a(j,i)
30      a(j,i)=hold
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Changing the columns
c       - - - - - - - - - - - - - - - - - - - - - - - - -
35      i=m(k)
        if (i-k) 45,45,38
38      jp=n*(i-1)
        do 40 j=1,n
        hold=-a(j,k)
        a(j,k)=a(j,i)
40      a(j,i)=hold
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Deviding column on the leader element.
c       The leader element is in 'biga'.
c       - - - - - - - - - - - - - - - - - - - - - - - - -
45      if (biga) 48,46,48
46      d=0.d0
        return
48      biga=-biga
        do 55 i=1,n
        if (i-k) 50,55,50
50      a(i,k)=a(i,k)/biga
55      continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        biga=-biga
        do 65 i=1,n
        hold=a(i,k)
        ij=i-n
        do 65 j=1,n
        ij=ij+n
        if (i-k) 60,65,60
60      if (j-k) 62,65,62
62      kj=ij-i+k
        a(i,j)=hold*a(k,j)+a(i,j)
65      continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Deviding the string on the leader element.
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        kj=k-n
        do 75 j=1,n
        kj=kj+n
        if (j-k) 70,75,70
70      a(k,j)=a(k,j)/biga
75      continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Composition of the leader elements.
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        d=d*biga
        a(k,k)=1.d0/biga
80      continue
        k=n
100     k=k-1
        if (k) 150,150,105
105     i=l(k)
        if (i-k) 120,120,108
108     do 110 j=1,n
        hold=a(j,k)
        a(j,k)=-a(j,i)
110     a(j,i)=hold
120     j=m(k)
        if (j-k) 100,100,125
125     do 130 i=1,n
        hold=a(k,i)
        a(k,i)=-a(j,i)
130     a(j,i)=hold
        goto 100
c       - - - - - - - - - - - - - - - - - - - - - - - - -
150     return
        end
c       =================================================
        subroutine recunit(nbytes)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        character*8 d1,t1,d2,t2
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Automatic calculating of the record unit length.
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        t1='abcdefgh'
        d1='        '
        t2='hgfedcba'
        d2='        '
        lrec=0
        iflag=1
200     lrec=lrec+1
        if (lrec.gt.8) then
          write(*,*)  'lrec > 8'
          call exit1
        endif
        open(unit=13,file='test.tmp',status='unknown',
     1  access='direct',recl=lrec)
        write(13,rec=1,err=210) t1
        write(13,rec=2,err=210) t2
        read(13,rec=1,err=210) d1
        read(13,rec=2,err=210) d2
        if (d1.ne.t1) goto 210
        if (d2.ne.t2) goto 210
        iflag=0
210     close(unit=13,status='delete')
        if (iflag.ne.0) goto 200
        nbytes=8/lrec
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
c       =================================================
        subroutine exit1
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        call exit(1)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
