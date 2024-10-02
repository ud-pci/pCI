Module str_fmt

    Use, Intrinsic :: iso_fortran_env, Only : dp => real64, int64

    Implicit None
    
    Private

    Public :: startTimer, stopTimer, FormattedMemSize, FormattedTime

  Contains

    Subroutine startTimer(start_time)
        Implicit None
        Integer(Kind=int64), Intent(Out) :: start_time

        Call system_clock(start_time)
        
    End Subroutine startTimer

    Subroutine stopTimer(start_time, timeStr)
        Implicit None
        Integer(Kind=int64), Intent(In) :: start_time
        Integer(Kind=int64)             :: end_time, clock_rate
        Real                            :: total_time
        Character(Len=16), Intent(Out) :: timeStr

        Call system_clock(count_rate=clock_rate)
        Call system_clock(end_time)
        total_time = (end_time - start_time) / clock_rate
        Call FormattedTime(total_time, timeStr)
        
    End Subroutine stopTimer

    Subroutine FormattedMemSize(bytes, outString)
        Implicit none

        Integer(Kind=8), Intent(In)     :: bytes
        Character(Len=*), Intent(InOut) :: outString
        Character(Len=24)               :: strBuffer
        Character*6                     :: units
        Data units /'BKMGTP'/

        Integer                         :: unitIdx
        Real(Kind=8)                    :: calcBytes

        If (bytes < 0) Then
            Write(strBuffer,'(I24)') bytes
            Write(0, '(A,A)') 'FormattedMemSize: negative byte size encountered: ', AdjustL(Trim(strBuffer))
            Error Stop
        End If

        unitIdx = 1
        calcBytes = Real(bytes,8)
        Do While ( unitIdx < 6 .and. calcBytes >= 1024.0 )
            unitIdx = unitIdx + 1
            calcBytes = calcBytes / 1024.0
        EndDo
        If ( unitIdx == 1 ) Then
            Write(strBuffer, '(I8,A)') bytes,' B'
        Else
            Write(strBuffer, '(F8.1,A)') calcBytes,' '//units(unitIdx:unitIdx)//'iB'
        EndIf
        outString = AdjustL(Trim(strBuffer))
        Return
    End Subroutine FormattedMemSize

    Subroutine FormattedTime(inpSec, outString)
        Implicit none

        Real, Intent(In)     :: inpSec
        Character(Len=*), Intent(InOut) :: outString
        Character(Len=24)               :: strBuffer
        Integer :: intSec, intMin, intHour, intDay, i, iskip
        Real(dp) :: days, hours, minutes, seconds, secMin, secHour, secDay, hourSec, minSec, remSec
        Integer, dimension(4) :: arrTime
        Character(Len=3), dimension(4) :: charArrTime
        Character(Len=4)               :: units

        Data units /'dhms'/

        intSec = int(inpSec)

        ! If time took 0 seconds, return immediately with output "0s"
        If (intSec == 0) Then
            outString = "00s"
            Return
        End If

        secMin = 60
        secHour = 60 * secMin
        secDay = 24 * secHour

        ! extract days
        days = floor(inpSec/secDay)

        ! extract hours
        hourSec = modulo(inpSec,secDay)
        hours = floor(hourSec/secHour)

        ! extract minutes
        minSec = modulo(hourSec,secHour)
        minutes = floor(minSec/secMin)

        ! extract remaining seconds
        remSec = modulo(minSec,secMin)
        seconds = ceiling(remSec)

        ! convert times to integer format
        intSec = int(seconds)
        intMin = int(minutes)
        intHour = int(hours)
        intDay = int(days)

        ! initialize time arrays
        arrTime = (/intDay,intHour,intMin,intSec/)
        charArrTime = (/"","","",""/)

        ! convert array of times into a string with units
        iskip = 0
        Do i=1,4
            If (arrTime(i) /= 0) Then
                If (arrTime(i) < 10) Then
                    Write(charArrTime(i-iskip),'(A,I1,A)') '0',arrTime(i), units(i:i)
                Else
                    Write(charArrTime(i-iskip),'(I2,A)') arrTime(i), units(i:i)
                End If
            Else
                iskip = iskip + 1
            End If
        End Do

        strBuffer = charArrTime(1) // " " // charArrTime(2) // " " //charArrTime(3) // " " // charArrTime(4)
        outString = AdjustL(Trim(strBuffer))
        Return
    End Subroutine FormattedTime
    
End module
