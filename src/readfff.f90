Module readfff

    Use params

    Implicit None

    Private

    Integer, Public :: Ird
    Real(dp), Public, Dimension(IP6, IPs) :: ArrP, ArrQ
    Public :: ReadF, WriteF, ReadFF

  Contains

    Subroutine ReadF (Kan,record,V1,V2,nrec)
        !FORTRAN-77, MS-DOS VERSION
        Implicit None

        Integer :: record,nrec, i, ii, nr1, nr2, kan
        Real(dp), Dimension(IP6) :: V1, V2

        ii=IP6
        nr1=2*record-1
        nr2=nr1+1
        Read(Kan,rec=nr1) (V1(I),i=1,ii)
        If (nrec.EQ.2) Read(Kan,rec=nr2) (V2(I),i=1,ii)

       Return
    End Subroutine ReadF

    Subroutine WriteF (kan,record,v1,v2,nrec)
        Implicit None
        ! fortran-77, ms-dos version

        Integer  :: record,nrec, i, ii, nr1, nr2, kan
        Real(dp), Dimension(IP6) :: V1, V2

        ii=IP6
        nr1=2*record-1
        nr2=nr1+1
        Write(kan,rec=nr1) (v1(i),i=1,ii)
        If (nrec.eq.2) Write(kan,rec=nr2) (v2(i),i=1,ii)

       Return
    End Subroutine WriteF

    Subroutine ReadFF (Kan,record,V1,V2,nrec)
        Implicit None

        Integer :: record, nrec, kan, ni, i
        Real(dp), Dimension(IP6) :: V1, V2

        If (Ird.NE.1) Then
            If (Ns.GT.IPs) Then
                Write(*,*) ' ReadFF: Basis set is too long.'
                Write(*,*) ' Max length:',IPs,' Ns=',Ns
                Stop
            End If
            Do ni=1,Ns
                Call ReadF(Kan,ni+4,V1,V2,2)
                Do i=1,IP6
                   ArrP(i,ni)=V1(i)
                   ArrQ(i,ni)=V2(i)
                End Do
            End Do
            Ird=1
        End If
        
        If (record.LE.4.OR.record.GT.Ns+4) Then
            Call ReadF (Kan,record,V1,V2,nrec)
        Else
            ni=record-4
            Do i=1,IP6
               V1(i)=ArrP(i,ni)
               V2(i)=ArrQ(i,ni)
            End Do
        End If

        Return
    End Subroutine ReadFF

End Module readfff