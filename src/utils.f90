Module utils
    ! This module contains general utility functions and subroutines for the pCI software package

    Implicit None

    Private

    Public :: CountSubstr, DetermineRecordLength

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

    Subroutine DetermineRecordLength(lrec, success)
        ! Determines the minimum usable record length for direct access files 
        Integer, intent(Out) :: lrec
        Logical, intent(Out) :: success

        Integer :: irecl, err_stat
        Character(Len=8) :: test_data, read_data

        test_data = 'abcdefgh'
        irecl = 1
        success = .false.

        Do
            ! Open a file to test writing and reading test_data with increasing record length
            Open(unit=10, file='test.tmp', status='UNKNOWN', access='DIRECT', recl=irecl, iostat=err_stat)
            If (err_stat /= 0) then
                Exit
            End If

            ! Try writing test_data
            Write(10, rec=1, iostat=err_stat) test_data
            If (err_stat /= 0) Then
                Close(unit=10, status='DELETE')
                irecl = irecl + 1
                Cycle
            End If

            ! Try reading read_data
            Read(10, rec=1, iostat=err_stat) read_data
            Close(unit=10, status='DELETE')
            If (err_stat == 0 .and. read_data == test_data) Then
                success = .true.
                lrec = irecl
                Exit
            End If
        End Do
    
        If (.not. success) Then
            lrec = -1
        End If

    End Subroutine DetermineRecordLength

End Module utils