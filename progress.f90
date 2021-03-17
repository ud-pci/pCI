Module progress

    Implicit None
    
    Private
    
    Public :: ProgressData, ProgressStart, ProgressTick, ProgressEnd
    
    Type ProgressData
        Integer(Kind=8)     :: currentCount, totalCount
        Double Precision    :: fracPerCount, currentFrac, nextMinorTick, nextMajorTick
        Real                :: projTimeStep
        Real                :: t0, tProj
        Integer             :: ioUnit
    End Type ProgressData
    
    Real, Parameter         :: minorFraction = 0.02D0
    Real, Parameter         :: majorFraction = 0.10D0

  Contains
  
    Subroutine ProgressStart(this, name, ioUnit, totalCount, projTimeStep)
        Implicit None
        
        Type(ProgressData), Intent(InOut)   :: this
        Character(Len=*), Intent(In)        :: name
        Integer, Intent(In)                 :: ioUnit
        Integer(Kind=8), Intent(In)         :: totalCount
        Real, Intent(In)                    :: projTimeStep
        
        this%ioUnit = ioUnit
        this%currentCount = 0
        this%totalCount = totalCount
        this%fracPerCount = 1.0D0 / Dble(totalCount)
        this%currentFrac = 0.0D0
        this%nextMinorTick = minorFraction
        this%nextMajorTick = majorFraction
        If (projTimeStep > 0.0 .and. projTimeStep < 60.0) then
            this%projTimeStep = 1500.0
        Else
            this%projTimeStep = projTimeStep
        End If
        Call cpu_time(this%t0)
        this%tProj = this%t0 + this%projTimeStep
        If (ioUnit > 0) then
            Write(ioUnit,fmt='(A,A,A)') '* Progress monitor "', Trim(AdjustL(name)), '" {'
            Flush(ioUnit)
        End If
    End Subroutine
    
    Subroutine ProgressTick(this)
        Implicit None
        
        Type(ProgressData), Intent(InOut)           :: this
        Real                                        :: now, dTime
        Character(Len=24)                           :: timeStr
        Logical                                     :: haveWritten
        
        If (this%currentCount < this%totalCount) then
            this%currentCount = this%currentCount + 1
            this%currentFrac = this%currentFrac + this%fracPerCount
            If (this%ioUnit > 0) then
                haveWritten = .false.
                If (ABS(this%currentFrac - this%nextMajorTick) < 1.0D-3) then
                    Write(this%ioUnit, fmt='(A,X,I3,A)', advance='no') '|-', NInt(this%currentFrac * 1.0D2), '%  '
                    haveWritten = .true.
                    this%nextMajorTick = this%nextMajorTick + majorFraction
                Else If (ABS(this%currentFrac - this%nextMinorTick) < 1.0D-3) then
                    Write(this%ioUnit, fmt='(A)', advance='no') '|-       '
                    haveWritten = .true.
                    this%nextMinorTick = this%nextMinorTick + minorFraction
                    ! Advance the minor tick again if it's coincident with the major:
                    If (ABS(this%nextMajorTick - this%nextMinorTick) < 1.0D-3) then
                        this%nextMinorTick = this%nextMinorTick + minorFraction
                    EndIf
                End If
                ! Check the time?
                If (this%projTimeStep > 0) then
                    Call cpu_time(now)
                    If (now > this%tProj) then
                        dTime = now - this%t0
                        ! Predict how much longer this will take:
                        dTime = (dTime / this%currentFrac) - dTime
                        Write(timeStr, fmt='(F18.2)') dTime
                        Write(this%ioUnit, fmt='(A,X,A,X,A)', advance='no') '(approx', Trim(AdjustL(timeStr)), 'seconds remaining)'
                        haveWritten = .true.
                        this%tProj = now + this%projTimeStep
                    End If
                End If
                
                If (haveWritten) then
                    Write(this%ioUnit,fmt='(A)') ''
                    Flush(this%ioUnit)
                EndIf
            End If
        End If
    End Subroutine
    
    Subroutine ProgressEnd(this)
        Implicit None
        
        Type(ProgressData), Intent(InOut)           :: this
        Real                                        :: now
        Character(Len=24)                           :: timeStr
        
        Call cpu_time(now)
        this%currentCount = this%totalCount
        this%currentFrac = 1.0
        If (this%ioUnit > 0) then
            Write(this%ioUnit, fmt='(A)') '}'
            now = 1.0 / (now - this%t0)
            now = now * this%totalCount
            Write(timeStr, fmt='(F20.1)') now
            Write(this%ioUnit, fmt='(A,X,A,X,A)') '* Overall performance was', Trim(AdjustL(timeStr)), ' s^-1'
            Flush(this%ioUnit)
        End If
    End Subroutine

End Module
        