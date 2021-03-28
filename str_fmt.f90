Module str_fmt

    Implicit None
    
    Private

    Public :: FormattedMemSize, FormattedTime

  Contains

    Subroutine FormattedMemSize(bytes, outString)
        Implicit none

        Integer(Kind=8), Intent(In)     :: bytes
        Character(Len=*), Intent(InOut) :: outString
        Character(Len=24)               :: strBuffer
        Character*6                     :: units
        Data units /'BKMGTP'/

        Integer                         :: unitIdx
        Real(Kind=8)                    :: calcBytes

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
        use params, only : dp
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
            outString = "0s"
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