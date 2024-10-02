c     ==================   21/04/96   =====================
      Program RPA_DTM             !### last update 15/10/17
c#    Uses files RPA_**.INT to change integrals in the file DTM.INT
      INCLUDE "rpa.par"
       implicit real*8 (a-h,o-z)
       common /Ns/Ns/Nso/Nso/Nsx/Nsx/Z/Z/Nint/Nint
     1        /Nr/Nr/Nv/Nv/Kl/Kl/Ita/Ita
       common /Nn/Nn(IPs)/Kk/Kk(IPs)/Ll/Ll(IPs)/Jj/Jj(IPs)
     >        /Ki/Ki(IPcs)
     >        /Ivi/Ivi(IPv)   /Ivk/Ivk(IPv)
     >        /Isi/Isi(IPv)   /Isk/Isk(IPv)
     >        /Iti/Iti(IPv)   /Itk/Itk(IPv)
     >        /Fp/Fp(IPv)     /Fm/Fm(IPv)     /F0v/F0v(IPv)
     >        /RPAsbt/RPAsbt(IPv)
     >        /SgRPA/SgRPA(IPv)               /Str_Rad/Str_Rad(IPv)
       common /Let/Let(9),Alet(IPcs),Blet(5)
       dimension ita1(IPcs)
       character*1 Let
       character*4 Alet,Blet
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        data ita1/+1,+1,+1,+1,-1,-1,-1,+1,+1,+1,+1,+1,+1/
        data Let/'s','p','d','f','g','h','i','k','l'/
        data Alet/'A_hf','B_hf','E1_L','EDM ','PNC ','E1_V','AM  ',
     >            'MQM ',' M1 ','E2  ','E3  ','M2  ','M3  '/
        data Blet/'Rint','RPA1','RPA2','RPA3','RPA4'/
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        write (*,5)
 5      format (/3X,' RPA_DTM  substitutes radial integrals',
     >         /3X,' of zero order operators from  DTM.INT',
     >         /3X,' by radial integrals calculated by',
     >         /3X,' RPA, RPA_sbt, Sg_RPA and Str_Rad',/)
        call DTM_in(IPcs,NPh)
        is=0
        num=0
        Open (11,file='RPA_DTM.RES',status='UNKNOWN')
        do Kl=1,IPcs
          if (Ki(Kl).EQ.1) then
            Ita= ita1(Kl)      !### - real(+1)/imaginary(-1)
            call RPA_in(kv)
            if (kv.NE.0) then
              is=is+1
              call Change(kv,nnew,nold,ndis,NPh)
              write (*,15) Alet(Kl),nnew,nold,ndis
 15           format (5X,'Operator ',A4,': ',I4,' integrals changed ',
     >                I4,' unchanged',I4,' skipped')
c              read (*,*)
              if (nnew.GT.0) Ki(Kl)=kv+1
              num=num+nnew
            end if
          end if
        end do
        if (num.EQ.0) then
          write (*,*) ' no new integrals found'
c          read(*,*)
        else
          write (*,*) ' saving changes to DTM.INT...'
          call DTM_out(IPcs)
        end if
        write (11,*) '     >>>>>>>>>> END <<<<<<<<<<'
       stop
      end
C     =================================================
      include "check_dim.inc"
C     =================================================
      subroutine DTM_in(IPcs,NPh)     !### Reads DTM.INT
      INCLUDE "conf.par"
C     --------------------
      Parameter(IPcs1= 13)
C     --------------------
       implicit real*8 (a-h,o-z)
       common /Ns/Ns/Nso/Nso/Z/Z/Nint/Nint/Nx/Nx
       common /Ll/Ll(IPs)/Ki/Ki(IPcs1)
       common /Let/Let(9),Alet(IPcs1),Blet(5)
       common /Rnt/Rnt(IPh)/Int/Int(IPh)
       character*1 Let
       character*4 Alet,Blet
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        IF (IPcs.NE.IPcs1) THEN
          write(*,*) 'Subr. DTM_in:'
          write(*,15) IPcs1,IPcs
 15       Format(' Parameter IPcs1=',I3,' must be equal to IPcs=',I3)
          write(*,*) 'Program is terminated'
          STOP
        END IF

        NPh= IPh

        Nx=IPx                  !### Nx is used for indexation of Rnt
        open (13,file='DTM.INT',status='OLD',form='UNFORMATTED',
     >        err=200)
        read (13,end=201,err=201) Ns,Nso,Z
        write (*,*) ' Ns =',Ns,' Nso=',Nso
        read (13) Nint,(Ll(i),i=1,Ns)
        write (*,*) ' Nint =',Nint
        write(*,*) 'IPcs =',IPcs
        read (13) (Ki(i),i=1,IPcs)
        write (*,*) (Ki(i),i=1,IPcs)
        write (*,*) ' reading integrals...'
        read (13) (Rnt(i),Int(i),i=1,Nint)
        write (*, 5) Nint,(Alet(i),Blet(ki(i)),i=1,IPcs)
 5      format (/1X,'### Radial integrals from DTM.INT (',
     >             ' Nint =',I6,' ) ###',
     >             /(4X,A4,' calculated by ',A4))
        close (13)
        is=0
        do i=1,IPcs
         if (ki(i).EQ.1) is=is+1
        end do
        if (is.GT.0) write (*,*) is,' zero order operators found'
       return
 200    write (*,*) ' Cannot open file DTM.INT'
       stop
 201   write (*,*) ' Error in file DTM.INT'
       stop
      end
C     =================================================
      subroutine DTM_out(IPcs)    !### Writes to DTM.INT
      INCLUDE "conf.par"
C     --------------------
      Parameter(IPcs1= 13)
C     --------------------
       implicit real*8 (a-h,o-z)
       common /Ns/Ns/Nso/Nso/Z/Z/Nint/Nint
       common /Ll/Ll(IPs)/Ki/Ki(IPcs1)
       common /Let/Let(9),Alet(IPcs1),Blet(5)
       common /Rnt/Rnt(IPh)/Int/Int(IPh)
       character*1 Let
       character*4 Alet,Blet
C     - - - - - - - - - - - - - - - - - - - - - - - - -
       IF (IPcs.NE.IPcs1) THEN
         write(*,*) 'Subr. DTM_in:'
         write(*,15) IPcs1,IPcs
 15      Format(' Parameter IPcs1=',I3,' must be equal to IPcs=',I3)
         write(*,*) 'Program is terminated'
         STOP
       END IF
       
        open (13,file='DTM.INT',status='OLD',form='UNFORMATTED',
     >        err=200)
        read (13)
        read (13)
        write (13) (Ki(i),i=1,IPcs)
        write (13) (Rnt(i),Int(i),i=1,Nint)
        write (*, 5) Nint,(Alet(i),Blet(ki(i)),i=1,IPcs)
 5      format (/1X,'### Radial integrals written to DTM.INT (',
     >             ' Nint =',I6,' ) ###',
     >             /(4X,A4,' calculated by ',A4))
        close (13)
       return
 200    write (*,*) ' Error in file DTM.INT'
       stop
      end
C     =================================================
      Subroutine RPA_in(kv)   !### Reads file RPA_n.INT
      INCLUDE "rpa.par"
       implicit real*8 (a-h,o-z)
       common /Ns1/Ns1/Nso1/Nso1/Nr/Nr/Nv/Nv/Kl/Kl/Nvs/Nvs/Nts/Nts
       common /Ivi/Ivi(IPv)   /Ivk/Ivk(IPv)
     >        /Isi/Isi(IPv)   /Isk/Isk(IPv)
     >        /Iti/Iti(IPv)   /Itk/Itk(IPv)
     >        /Fp/Fp(IPv)     /Fm/Fm(IPv)     /F0v/F0v(IPv)
     >        /RPAsbt/RPAsbt(IPv)
     >        /SgRPA/SgRPA(IPv)               /Str_Rad/Str_Rad(IPv)
       character*1 cn(0:9),cf(9)
       character*9 FNAME
       equivalence (FNAME,cf(1))
        data cn /'0','1','2','3','4','5','6','7','8','9'/
        data cf /'R','P','A','_','n','.','I','N','T'/
C     - - - - - - - - - - - - - - - - - - - - - - - - -
        itest=IPs
        call Check_Dim(itest) ! IPs must be the same as in conf.par
!old        cf(5)=cn(Kl)
        if (Kl.LT.10) then
          cf(4)='_'
          cf(5)=cn(Kl)
        else
          cf(4)='1'
          if (Kl.EQ.10) then
            cf(5)='0'
          else
            cf(5)=cn(Kl-10)
          end if
        end if
        kv=0
        open(13,file=FNAME,status='OLD',form='UNFORMATTED',err=200)
        read (13,err=300) Ns1,Nso1,nsx,Nr,Nv
        if (Nr.GT.IPv.OR.Nv.GT.IPv) then
          write (*,*) '>>> ERROR: dimension is too high in ',FNAME
          write (*,*) ' Nr =',Nr,' Nv =',Nv
          write (*,*) ' current dimension: ',IPv
c          read (*,*)
         return
        end if
        read (13,err=300)
        read (13,err=300)
        if (Nv.GE.1) then
          read (13,err=400) (Ivi(i),Ivk(i),i=1,Nv)
          read (13,err=400) (F0v(i),Fp(i),Fm(i),i=1,Nv)
          write(*,15) FNAME,ns1,nso1,nsx,Nr,Nv
 15       format(5X,'RPA MEs from ',A9,' read OK: Ns=',I3,
     >           ' Nso=',I2,' Nsx=',I2,' Nr=',I4,' Nv=',I4)
          kv=1
          read (13,err=100,end=100) nvt,(RPAsbt(i),i=1,nvt)
          if (nvt.GT.0) then
             kv=2
             write (*,25) nvt
 25          format(5X,I4,' RPA_sbt MEs read OK.')
             if (nvt.NE.Nv)
     >          write (*,*) ' Warning: number of RPA_sbt MEs:',nvt
          end if
          read (13,err=100,end=100) Nvs,(Isi(i),Isk(i),i=1,Nvs)
          read (13,err=100,end=100) (SgRPA(i),i=1,nvs)
          kv=3
          write (*,35) Nvs
 35       format(5X,I4,' Sg_RPA MEs read OK.')
          read (13,err=100,end=100) Nts,(Iti(i),Itk(i),i=1,Nts)
          read (13,err=100,end=100) (Str_Rad(i),i=1,Nts)
          kv=4
          write (*,45) Nts
 45       format(5X,I4,' Str_Rad MEs read OK.')
        end if
 100    close(13)
        if (kv.GT.1) then
          write (*,*)
     >      ' Choose: RPA(1) + RPA_sbt(2) + Sg_RPA(3) + Str_Rad(4): '
          read (*,*) kv1
          kv=min(kv,kv1)
          if (kv.GT.2.AND.nvt.GT.0) then
            write (*,*) ' Include RPA_sbt(1) or skip it (0): '
            read (*,*) ksbt
            if (ksbt.EQ.0) kv=-kv
          end if
        end if
       return
 200    write (*,*) FNAME,' is absent...'
       return
 300    write (*,*) FNAME,' is damaged...'
       return
 400    write (*,*) ' no valence MEs in ',FNAME
       return
      end
C     =================================================
      Subroutine Change(kv,nnew,nold,ndis,NPh) !### substitutes zero
                                               !### order integrals
      INCLUDE "rpa.par"
C     ----------------------
      Parameter(IPh1= 350000)
C     ----------------------
       implicit real*8 (a-h,o-z)
       common /Ns1/Ns1/Nso1/Nso1/Nso/Nso/Nr/Nr/Nv/Nv/Nx/Nx
     >        /Nvs/Nvs/Nts/Nts/Kl/Kl/Nint/Nint/Ita/Ita
       common /Ivi/Ivi(IPv)   /Ivk/Ivk(IPv)
     >        /Isi/Isi(IPv)   /Isk/Isk(IPv)
     >        /Iti/Iti(IPv)   /Itk/Itk(IPv)
     >        /Fp/Fp(IPv)     /Fm/Fm(IPv)     /F0v/F0v(IPv)
     >        /RPAsbt/RPAsbt(IPv)
     >        /SgRPA/SgRPA(IPv)               /Str_Rad/Str_Rad(IPv)
       common /Rnt/Rnt(IPh1)/Int/Int(IPh1)
       common /Let/Let(9),Alet(IPcs),Blet(5)
       character*1 Let
       character*4 Alet,Blet
C     - - - - - - - - - - - - - - - - - - - - - - - - -
       IF (IPh1.NE.NPh) THEN
         write(*,*) 'Subr. Change:'
         write(*,65) IPh1,NPh
 65      Format(' Parameter IPh1=',I6,' must be equal to IPh=',I6)
         write(*,*) 'Program is terminated'
         STOP
       END IF

        if (kv.LT.0) then
          kv=-kv
          ksbt=0
        else
          ksbt=1
        end if
        write (11,5) Alet(Kl)
 5      format (65('='),/1X,'Operator ',A4,/4X,
     >   'i   k    DTM      RPA0     RPA1    RPA_sbt  Sg_RPA  Str_Rad')
        nnew=0
        nold=0
        ndis=0
        do i=1,Nint
           ind=Int(i)
           k=ind/(Nx*Nx)
           if (k.EQ.Kl) then
              nold=nold+1
              na=(ind-Nx*Nx*k)/Nx+Nso
              nb= ind-Nx*Nx*k-Nx*(na-Nso)+Nso
              r=Rnt(i)
              is=0
              rt=0.d0
              rs=0.d0
              rst=0.d0
              do n=1,Nv
                 if (Ivi(n).EQ.na.AND.Ivk(n).EQ.nb) then
                    is=n
                    r0=Ita*F0v(n)
                    rm=Ita*(Fp(n)-Fm(n))
                 end if
              end do
              if (is.NE.0) then
                 err=dabs(r-r0)/(dabs(r)+dabs(r0)+1.d-30)
                 if (err.LT.1.d-5) then
                    nnew=nnew+1
                    if (ksbt*kv.GE.2) rt=Ita*RPAsbt(is)
                    if (kv.GE.3) then
                       do n=1,Nvs
                          if (Isi(n).EQ.na.AND.Isk(n).EQ.nb) then
                             rs=Ita*SgRPA(n)
                          end if
                       end do
                    end if
                    if (kv.GE.4) then
                       do n=1,Nts
                          if (Iti(n).EQ.na.AND.Itk(n).EQ.nb) then
                             rst=Ita*Str_Rad(n)
                          end if
                       end do
                    end if
                    Rnt(i)=rm+rt+rs+rst
                    if (ksbt.EQ.1) then
                      write (11,15) na,nb,r,r0,rm,rt,rs,rst
 15                   format (1X,2I4,6F9.4)
                    else
                      write (11,25) na,nb,r,r0,rm,   rs,rst
 25                   format (1X,2I4,3F9.4,' skipped ',2F9.4)
                    end if
                 else
                    write (*,35) Alet(k),na,nb,r,err
 35                 format (1X,'Discrepancy for ',
     >                   A4,'(',I3,',',I3,') =',E15.8,' err =',E9.2)
                    write (11,45) na,nb,r,r0,err
 45                 format (1X,2I4,2F9.4,' err= ',E10.3,' Skipped!')
                    ndis=ndis+1
                 end if
              else
                write (11,55) na,nb,r
 55             format (1X,2I4,F9.4,5X,'Absent in RPA!')
              end if
           end if
        end do
        nold=nold-nnew
       return
      end
