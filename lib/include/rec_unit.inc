c       =================================================
        subroutine recunit
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        character*8 d1,t1,d2,t2
        common /ipmr/ipmr
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Determination of the record unit.
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
          stop
        end if
        open(unit=13,file='test.tmp',status='unknown',
     >       access='direct',recl=lrec)
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
        ipmr=4/nbytes
C        write( *,5) nbytes
C5       format(/2x,'Record unit =',i2,' byte(s)')
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
