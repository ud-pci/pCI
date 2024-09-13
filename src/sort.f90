program sort
    use, intrinsic :: iso_fortran_env, only : dp => real64, int64
    use str_fmt, only : startTimer, stopTimer, FormattedTime

    implicit none

    integer(kind=int64) :: start_time
    character(Len=16)   :: time_str

    ! Start program timer
    call startTimer(start_time)

    call sort_matrix('CONFp.JJJ', 'CONF.JJJ')
    print*, ''
    call sort_matrix('CONFp.HIJ', 'CONF.HIJ')

    ! End program timer
    call stopTimer(start_time, time_str)
    write(*,'(2X,A)'), 'TIMING >>> Total computation time of sort was '// trim(time_str)

contains

    subroutine sort_matrix(original_filename, filename)
        implicit none

        character(len=*) :: filename, original_filename

        integer :: num_procs, err_stat, j, cnt, cnt2, num_dets
        integer(kind=int64) :: i, start_time, start_time2, num_total_elements, num_zeros
        character(len=16)   :: time_str, time_str2, err_msg, count_str

        type matrix
            integer, dimension(:), allocatable :: indices1, indices2
            real(dp), dimension(:), allocatable :: values
        end type matrix

        type(matrix) :: mat
    
        ! Read unsorted matrix elements
        call startTimer(start_time)
        call read_matrix_serial(original_filename, num_total_elements, num_zeros, mat%indices1, mat%indices2, mat%values)
        call stopTimer(start_time, time_str2)
        print*, trim(adjustl(original_filename)) // ' read in ' // trim(adjustl(time_str2))
    
        ! Print number of zero-valued matrix elements to be removed
        if (num_zeros > 0) then
            write(count_str, '(I16)') num_zeros
            print*, 'Number of zero-valued matrix elements to be removed: ' // trim(adjustl(count_str))
            write(count_str, '(I16)') num_total_elements - num_zeros
            print*, 'Number of matrix elements after zeros are removed: ' // trim(adjustl(count_str))
        end if
        
        ! Print number of determinants
        num_dets = maxval(mat%indices1)
        write(count_str,'(I16)') num_dets
        print*, 'Total number of determinants: ', trim(adjustl(count_str))
    
        ! Sort matrix elements
        call startTimer(start_time2)
        call quicksort_matrix(mat%indices1, mat%indices2, mat%values, 1, num_total_elements)
        call stopTimer(start_time2, time_str2)
        print*, trim(adjustl(filename)) // ' sorted in ' // trim(adjustl(time_str2))
        
        ! Write sorted matrix element file
        call startTimer(start_time2)
        open(unit=15,file=filename,status='UNKNOWN',form='unformatted',access='stream',iostat=err_stat,iomsg=err_msg)
        write(15) num_total_elements-num_zeros
        do i = 1_int64, num_total_elements
            if (mat%values(i) /= 0) write(15) mat%indices1(i), mat%indices2(i), mat%values(i)
        end do
        close(15)

        deallocate(mat%indices1, mat%indices2, mat%values)

        call stopTimer(start_time2, time_str)
        print*, 'Writing ' // trim(adjustl(filename)) // ' took ' // trim(adjustl(time_str))
    end subroutine sort_matrix

    subroutine read_matrix_serial(filename, num_total_elements, num_zeros, indices1, indices2, values)
        ! read parallel matrix element file written by pconf program using a single process
        implicit none

        character(len=*), intent(in) :: filename
        integer(kind=int64), intent(out) :: num_total_elements, num_zeros
        integer, dimension(:), allocatable, intent(out) :: indices1, indices2
        real(dp), dimension(:), allocatable, intent(out) :: values

        integer(kind=int64) :: i
        integer :: n, err_stat, num_procs, num_elements
        character(len=16) :: err_msg, count_str
        
        ! Read num_cores
        open(unit=15,file='nprocs.conf',status='UNKNOWN',form='unformatted',access='stream',iostat=err_stat,iomsg=err_msg)
        read(15) num_procs
        close(15)
    
        Write(count_str, '(I4)') num_procs
        print*, 'Number of processes used to write ', trim(adjustl(filename)), ': ', trim(adjustl(count_str))

        open(unit=15,file=filename,status='UNKNOWN',form='unformatted',access='stream',iostat=err_stat,iomsg=err_msg)
        num_total_elements = 0_int64
        do n = 1, num_procs
            read(15) num_elements
            num_total_elements = num_total_elements + num_elements
        end do

        ! Print total number of matrix elements read from original file
        write(count_str,'(I16)') num_total_elements
        print*, 'Total number of matrix elements to be read from ', trim(adjustl(filename)), ': ', trim(adjustl(count_str))
        
        allocate(indices1(num_total_elements), indices2(num_total_elements), values(num_total_elements))
        do i = 1, num_total_elements
            read(15) indices1(i)
        end do
        if (num_total_elements > 214748647) print*, 'first array of indices read'

        do i = 1, num_total_elements
            read(15) indices2(i)
        end do
        if (num_total_elements > 214748647) print*, 'second array of indices read'

        num_zeros = 0_dp
        do i = 1, num_total_elements
            read(15) values(i)
            if (values(i) == 0) num_zeros = num_zeros + 1
        end do
        if (num_total_elements > 214748647) print*, 'array of matrix element values read'
        close(15)

    end subroutine read_matrix_serial

    recursive subroutine quicksort(arr, left, right)
        ! serial implementation of the QuickSort algorithm
        ! 
        ! inputs:
        ! arr - integer array to sort
        ! left, right - integer positions of first and last elements of sub-array to sort
        !
        implicit none

        integer, dimension(:), allocatable, intent(inout) :: arr
        integer(kind=int64), intent(in) :: left, right
        integer(kind=int64) :: pivot, i, j
        integer :: temp

        if (left < right) then
            pivot = arr((left + right) / 2) 
            i = left
            j = right

            do while (i <= j)
                do while (arr(i) < pivot)
                    i = i + 1
                end do
                do while (arr(j) > pivot)
                    j = j - 1
                end do
                if (i <= j) then
                    temp = arr(i)
                    arr(i) = arr(j)
                    arr(j) = temp
                    i = i + 1
                    j = j - 1
                end if
            end do

            call quicksort(arr, left, j)
            call quicksort(arr, i, right)
        end if
    end subroutine quicksort

    recursive subroutine quicksort_matrix(arr1, arr2, arr3, left, right)
        ! serial implementation of the QuickSort algorithm applied to the custom Matrix datatype
        implicit none

        integer, dimension(:), allocatable, intent(inout) :: arr1, arr2
        real(dp), dimension(:), allocatable, intent(inout) :: arr3
        integer(kind=int64), intent(in) :: left, right
        integer(kind=int64) :: pivot1, pivot2, i, j
        integer :: temp1, temp2
        real(dp) :: temp3

        if (left < right) then
            pivot1 = arr1((left + right) / 2) 
            pivot2 = arr2((left + right) / 2)
            i = left
            j = right

            do while (i <= j)
                do while ((arr1(i) < pivot1) .or. (arr1(i) == pivot1 .and. arr2(i) < pivot2))
                    i = i + 1
                end do

                do while (arr1(j) > pivot1 .or. (arr1(j) == pivot1 .and. arr2(j) > pivot2))
                    j = j - 1
                end do

                if (i <= j) then
                    call swapI(arr1(i), arr1(j))
                    call swapI(arr2(i), arr2(j))
                    call swapR(arr3(i), arr3(j))
                    i = i + 1
                    j = j - 1
                end if

            end do

            call quicksort_matrix(arr1, arr2, arr3, left, j)
            call quicksort_matrix(arr1, arr2, arr3, i, right)
        end if
    end subroutine quicksort_matrix

    subroutine swapI(a, b)
        implicit none

        integer, intent(inout) :: a, b
        integer :: temp

        temp = a
        a = b
        b = temp
    end subroutine swapI

    subroutine swapR(a, b)
        implicit none

        real(dp), intent(inout) :: a, b
        real(dp) :: temp

        temp = a
        a = b
        b = temp
    end subroutine swapR

end program sort