program check_matrices
    use, intrinsic :: iso_fortran_env, only : dp => real64, int64

    implicit none

    character(len=16) :: filename

    print*, 'Enter matrix file (CONFp.HIJ or CONFp.JJJ) to check: '
    read(*, *) filename
    call read_matrix(filename)

contains

    subroutine read_matrix(filename)
        ! read matrix element file written by pconf program using a single process
        implicit none

        character(len=*), intent(in) :: filename
        integer(kind=int64) :: num_total_elements, num_zeros

        integer(kind=int64) :: i
        integer :: n, err_stat, num_procs, num_elements
        character(len=16) :: err_msg, count_str
        logical :: found_zero
        
        ! Read num_cores
        open(unit=15,file='nprocs.conf',status='UNKNOWN',form='unformatted',access='stream',iostat=err_stat,iomsg=err_msg)
        read(15) num_procs
        close(15)
    
        Write(count_str, '(I4)') num_procs
        print*, 'Number of processes used to write ', trim(adjustl(filename)), ': ', trim(adjustl(count_str))

        open(unit=15,file=filename,status='OLD',form='unformatted',access='stream',iostat=err_stat,iomsg=err_msg)
        num_total_elements = 0_int64
        found_zero = .false.
        do n = 1, num_procs
            read(15) num_elements
            if (num_elements == 0) then
                print*, 'core', n-1, 'has 0 matrix elements!'
                found_zero = .true.
            end if
            num_total_elements = num_total_elements + num_elements
        end do
        print*, 'Number of matrix elements in ', filename, ':', num_total_elements
        if (found_zero) print*, 'WARNING: check results of pconf to see if number of matrix elements match'

    end subroutine read_matrix

end program check_matrices