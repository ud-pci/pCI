Module utils
    ! This module contains general utility functions and subroutines for the pCI software package

    Use, Intrinsic :: iso_fortran_env, Only : dp => real64, int64

    Implicit None

    Private

    Public :: CountSubstr, DetermineRecordLength, CheckRecordLength, ToUpperString, WriteSkeleton, &
              CsvReal, CsvSci

Contains

    Function CountSubstr(str, substr) Result(count)
        ! This function counts the number of occurences of a substring in a string
        Implicit None

        Integer :: i, position, count
        Character(Len=*), Intent(In) :: str
        Character(Len=*), Intent(In) :: substr

        count = 0
        i = 1
        Do 
            position = index(str(i:), substr)
            If (position == 0) Return
            count = count + 1
            i = i + position + len(substr) - 1
        End Do

    End Function CountSubstr

    Subroutine DetermineRecordLength(lrec)
        ! Determines the minimum usable record length for direct access files
        Integer, intent(Out) :: lrec

        Real(dp) :: tmp

        Inquire(IOLENGTH=lrec) tmp

    End Subroutine DetermineRecordLength

    Subroutine CheckRecordLength(unit, filename, recl)
        ! Verifies that the file's size is consistent with the expected record length.
        ! Throws an error if the file size is not a multiple of the expected record size.
        Integer, Intent(In) :: unit, recl
        Character(*), Intent(In) :: filename

        Integer(int64) :: file_size, bytes_per_record
        Integer :: lrec
        Real(dp) :: tmp

        Inquire(IOLENGTH=lrec) tmp
        bytes_per_record = int(recl, int64) * 8_int64 / int(lrec, int64)

        Inquire(unit=unit, SIZE=file_size)
        If (file_size < 0) Then
            Write(*,'(2A)') 'ERROR: INQUIRE SIZE not supported, cannot verify record length for ', trim(filename)
            Stop
        End If

        If (mod(file_size, bytes_per_record) /= 0) Then
            Write(*,'(A)')    'ERROR: record length inconsistency'
            Write(*,'(2A)')   '  File:                  ', trim(filename)
            Write(*,'(A,I0)') '  File size (bytes):     ', file_size
            Write(*,'(A,I0)') '  Expected record (bytes):', bytes_per_record
            Stop
        End If

    End Subroutine CheckRecordLength

    Pure Function ToUpperString(str) Result(upper)
        ! Converts every lowercase letter (a-z) in a string to uppercase (A-Z)
        Character(*), Intent(In) :: str
        Character(len(str)) :: upper
        Integer :: i, c

        Do i = 1, Len(str)
            c = IChar(str(i:i))
            If (c >= IChar('a') .and. c <= IChar('z')) Then
                upper(i:i) = Char(c - 32)
            Else
                upper(i:i) = str(i:i)
            End If
        End Do
    End Function ToUpperString

    Function CsvReal(x, d) Result(s)
        Implicit None
        Real(dp), Intent(In) :: x
        Integer, Intent(In) :: d
        Character(Len=30) :: s, tmp
        Character(Len=12) :: fmt
        Write(fmt, '("(F28.",I0,")")') d
        Write(tmp, fmt) x
        s = Trim(AdjustL(tmp))
    End Function CsvReal

    Function CsvSci(x, d) Result(s)
        Implicit None
        Real(dp), Intent(In) :: x
        Integer, Intent(In) :: d
        Character(Len=30) :: s, tmp
        Character(Len=12) :: fmt
        Write(fmt, '("(ES28.",I0,")")') d
        Write(tmp, fmt) x
        s = Trim(AdjustL(tmp))
    End Function CsvSci

    Subroutine WriteSkeleton(filename, lines)
        ! Write a skeleton key=value input file when the file is not found,
        ! then stop so the user can fill it in and rerun.
        Character(Len=*), Intent(In) :: filename
        Character(Len=*), Intent(In) :: lines(:)
        Integer :: i, u
        Open(newunit=u, file=filename, status='REPLACE')
        Do i = 1, Size(lines)
            Write(u, '(A)') Trim(lines(i))
        End Do
        Close(u)
        Write(*,'(3A)') Trim(filename), ' not found — wrote a skeleton. Fill it in and rerun.'
        Stop
    End Subroutine WriteSkeleton

End Module utils
