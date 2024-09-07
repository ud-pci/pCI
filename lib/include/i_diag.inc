c     =================================================
c     Iterative diagonalization of diagonally-dominated matrix
c     =================================================
      subroutine I_Diag(nd,  ! real dimension
     >                  nd1, ! declared dimension
     >                  H,   ! initial matrix (not changed)
     >                  Z,   ! matrix of eigenvectors in columns
     >                  D,   ! eigenvalues
     >                  X,   ! reserved space
     >                  res) ! max residue
       implicit real*8 (a-h,o-z)
       common /Sf/Sf         !### damping factor
       dimension H(nd1,nd1),Z(nd1,nd1),D(nd1),X(nd1)
c     - - - - - - - - - - - - - - - - - - - - - - - - -
        if (Sf.LT.1.d-3) Sf=1.d0                 ! - no damping

        call Test_Z(nd,nd1,H,Z,D,res1)

        do n=1,nd

          s=0.d0                                 !### new eigenvalue:
          e=0.d0                                 !#### e_n=<X|Z_n>/<Z_n|Z_n>
          zx=0.d0
          nzx=0
          do k=1,nd
            xk=0.d0
            a=Z(k,n)
            a2=a*a
            if (a2.GT.zx) then
              zx=a2
              nzx=k                              !nzx - dominant projection
            end if
            s=s+a2
            do m=1,nd
              if(m.NE.k) then
                xk=xk+H(k,m)*Z(m,n)
              else
                xd=H(k,k)*a
              end if
            end do
            X(k)=xk                              !### X=(H-Diag)*Z_n
            e=e+a*(xk+xd)
          end do
          e=e/s                                  !### new eigenvalue
c          write (11,5) n,s,e,D(n)
 5        format(4X,'n=',I3,' <n|n>=',F14.8,' E_n=',2E14.7)
          D(n)=e

          do k=1,nd                              !### new eigenvector
            if (k.EQ.nzx) then
              X(k)=Z(k,n)
            else
              X(k)=Sf*X(k)/(e-H(k,k)) + (1.d0-Sf)*Z(k,n)
            end if
          end do

          do m=1,n-1                             !### orthogonalization
            s=0.d0
            do k=1,nd
              s=s+Z(k,m)*X(k)
            end do
c            write(11,15) n,m,s
 15         format(4X,'<',I3,'|',I3,'>=',E14.7)
            do k=1,nd
              X(k)=X(k)-s*Z(k,m)
            end do
          end do

          s=0.d0                                 !### normalization
          do k=1,nd
            a=X(k)
            s=s+a*a
          end do
c          write(11,15) n,n,s
          s=1.d0/dsqrt(s)

          do k=1,nd
            Z(k,n)=s*X(k)
          end do

        end do

        call Test_Z(nd,nd1,H,Z,D,res)
        write( *,25) res1,res
        write(11,25) res1,res
 25     format(4X,'I_Diag: residue was',E10.3,
     >           ' and now is',E10.3)
        if (res/res1.GT.0.9d0) then
          read(*,*)
          stop
        end if
       return
      end
c     =================================================
      subroutine Test_Z(nd,   ! real dimension
     >                  nd1,  ! declared dimension
     >                  H,    ! initial matrix (not changed)
     >                  Z,    ! matrix of eigenvectors in columns
     >                  D,    ! eigenvalues
     >                  x_max)! max residue for eigenvalue eq-n
       implicit real*8 (a-h,o-z)
       dimension H(nd1,nd1),Z(nd1,nd1),D(nd1)
c     - - - - - - - - - - - - - - - - - - - - - - - - -
        x_max=0.d0
        do n=1,nd
          x=0.d0
          write(11,5) n,D(n)
 5        format(4X,'Eivenvector ',I3,' E=',E15.8)
          do k=1,nd
            xk=0.d0
            do m=1,nd
              xk=xk+H(k,m)*Z(m,n)
            end do
            dxk=(xk-D(n)*Z(k,n))
            x=x+dxk*dxk
c            write(11,15) k, xk, dxk
 15         format(4X,'X(',I3,')=',E12.5,' dX=',E12.5)
          end do
          x=dsqrt(x)
          if (x_max.LT.x) x_max=x
          write(11,25) n,x
 25       format(4X,'Residue for n=',I3,' is ',E12.5)
        end do
       return
      end
