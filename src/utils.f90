module utils
    ! This module contains general utility functions and subroutines for the pCI software package

    implicit none

    private

    public :: CountSubstr

contains

    function CountSubstr(str, substr) result(count)
        ! This function counts the number of occurences of a substring in a string
        implicit none

        integer :: i, position, count
        character(len=*), intent(in) :: str
        character(len=*), intent(in) :: substr

        count = 0
        i = 1
        do 
            position = index(str(i:), substr)
            if (position == 0) return
            count = count + 1
            i = i + position + len(substr) - 1
        end do

    end function CountSubstr

end module utils