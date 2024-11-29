Module wigner

    Use params 

    Implicit None

    Private
    
    Public :: FJ3, FJ6, FJ9, DELFJ, NGFJ, swap, sortArray

  Contains

    Real(dp) Function FJ3(a1, a2, a3, a4, a5, a6)
        Implicit None

        Real(dp), Intent(In) :: a1, a2, a3, a4, a5, a6

        Integer :: k, l, ku, km
        Real(dp) :: eps, un, all
        Real(dp), Dimension(9) :: x
        Integer, Dimension(9) :: lu, li
        Real(dp), Dimension(3) :: y, z
        Integer, Dimension(3) :: ma, na

        eps = 0.0001_dp

        If (dabs(a4 + a5 + a6) >= eps) Then
            fj3 = 0_dp
            Return
        End If

        x(1) = a1 + a2 - a3
        x(2) = a1 - a2 + a3
        x(3) = -a1 + a2 + a3
        x(4) = a1 + a4
        x(5) = a1 - a4
        x(6) = a2 + a5
        x(7) = a2 - a5
        x(8) = a3 + a6
        x(9) = a3 - a6
        x = Anint(x * 100_dp) / 100_dp
        
        Do k = 1, 9
            If (x(k) < 0_dp) Then
                fj3 = 0_dp
                Return
            Else
                lu(k) = Int(x(k))
                If (dabs(lu(k) - x(k)) >= eps) Then
                    fj3 = 0.0_dp
                    Return
                End If
            End If
        End Do
        
        y(1) = a1 + a2 - a3
        y(2) = a1 - a4
        y(3) = a2 + a5
        y = Anint(y * 100_dp) / 100_dp
        
        z(1) = 0_dp
        z(2) = -a3 + a2 - a4
        z(3) = -a3 + a1 + a5
        z = Anint(z * 100_dp) / 100_dp
        
        DO k = 1, 3
            ma(k) = Int(y(k))
            na(k) = Int(z(k))
            If (dabs(ma(k) - y(k)) >= eps .or. dabs(na(k) - z(k)) >= eps) Then
                fj3 = 0.0_dp
                Return
            End If
        End Do
        
        Call sortArray(ma, 3)
        Call sortArray(na, 3)
        
        If (ma(1) < 0 .or. ma(1) < na(3)) Then
            fj3 = 0.0_dp
            Return
        End If
        
        li(1) = a1 + a2 + a3 + 1_dp
        Do k = 1, 2
            li(k + 1) = ma(k + 1) - ma(1)
            li(k + 3) = ma(k + 1) - ma(1)
            li(k + 5) = na(3) - na(k)
            li(k + 7) = na(3) - na(k)
        End Do
        
        Call sortArray(lu, 9)
        Call sortArray(li, 9)
        
        all = 0_dp
        Do k = 1, 9
            If (lu(k) /= li(k)) Then
                If (lu(k) > li(k)) Then
                    Call NGFJ(1, li(k), lu(k), all)
                Else
                    Call NGFJ(-1, lu(k), li(k), all)
                End If
            End If
        End Do
        
        un = all / 2_dp
        km = ma(1) - na(3) + 1

        fj3 = 0.0_dp
        Do ku = 1, km
            k = na(3) + ku - 1
            all = 0.0_dp
            Call NGFJ(-1, 0, ma(1) - k, all)
            Call NGFJ(-1, ma(2) - ma(1), ma(2) - k, all)
            Call NGFJ(-1, ma(3) - ma(1), ma(3) - k, all)
            Call NGFJ(-1, 0, k - na(3), all)
            Call NGFJ(-1, na(3) - na(1), k - na(1), all)
            Call NGFJ(-1, na(3) - na(2), k - na(2), all)
            l = k + a1 - a2 - a6
            fj3 = fj3 + (-1)**l * dexp(all)
        End Do
        
        fj3 = dexp(un) * fj3

    End Function FJ3

    Real(dp) Function FJ6(a1, a2, a3, a4, a5, a6)
        Implicit None

        Real(dp), Intent(In) :: a1, a2, a3, a4, a5, a6

        Integer :: k, l, km
        Real(dp) :: eps, x1, x2, x3, x4, y1, y2, y3, un, all
        Integer, Dimension(4) :: ma
        Integer, Dimension(3) :: na
        Real(dp), Dimension(12) :: x
        Integer, Dimension(13) :: lu, li

        fj6 = 0.0_dp

        x1 = a1 + a2 + a3
        x2 = a1 + a5 + a6
        x3 = a4 + a2 + a6
        x4 = a4 + a5 + a3
        y1 = a1 + a2 + a4 + a5
        y2 = a2 + a3 + a5 + a6
        y3 = a3 + a1 + a6 + a4

        eps = 0.0001_dp
        ma(1) = x1 + eps
        ma(2) = x2 + eps
        ma(3) = x3 + eps
        ma(4) = x4 + eps
        na(1) = y1 + eps
        na(2) = y2 + eps
        na(3) = y3 + eps

        If (Any(dabs(ma - [x1, x2, x3, x4]) >= eps) .or. &
            Any(dabs(na - [y1, y2, y3]) >= eps)) Return

        Call sortArray(ma, 4)
        Call sortArray(na, 3)

        If (ma(4) > na(1)) Return
        
        x(1) = a1 + a2 - a3
        x(2) = a1 - a2 + a3
        x(3) = -a1 + a2 + a3
        x(4) = a1 + a5 - a6
        x(5) = a1 - a5 + a6
        x(6) = -a1 + a5 + a6
        x(7) = a4 + a2 - a6
        x(8) = a4 - a2 + a6
        x(9) = -a4 + a2 + a6
        x(10) = a4 + a5 - a3
        x(11) = a4 - a5 + a3
        x(12) = -a4 + a5 + a3
        x = Anint(x * 100_dp) / 100_dp

        Do k = 1, 12
            If (x(k) >= 0) lu(k) = x(k) + eps
        End Do
        lu(13) = ma(4) + 1
        li(1:3) = ma(1:3) + 1
        li(4:6) = ma(4) - ma(1:3)
        li(7:9) = li(4:6)
        li(10:11) = na(2:3) - na(1)
        li(12:13) = li(10:11)

        Call sortArray(lu, 13)
        Call sortArray(li, 13)

        all = 0._dp
        Do k = 1, 13
            If (lu(k) /= li(k)) Then
                If (lu(k) < li(k)) Then
                    Call NGFJ(-1, lu(k), li(k), all)
                Else
                    Call NGFJ(1, li(k), lu(k), all)
                End If
            End If
        End Do
    
        un = all / 2_dp
        km = na(1) - ma(4) + 1
        Do l = 1, km
            k = ma(4) + l - 1
            all = 0_dp
            Call NGFJ(1, ma(4) + 1, k + 1, all)
            Call NGFJ(-1, ma(4) - ma(1), k - ma(1), all)
            Call NGFJ(-1, ma(4) - ma(2), k - ma(2), all)
            Call NGFJ(-1, ma(4) - ma(3), k - ma(3), all)
            Call NGFJ(-1, 0, k - ma(4), all)
            Call NGFJ(-1, 0, na(1) - k, all)
            Call NGFJ(-1, na(2) - na(1), na(2) - k, all)
            Call NGFJ(-1, na(3) - na(1), na(3) - k, all)
            fj6 = fj6 + (-1)**k * dexp(all)
        End Do
    
        fj6 = dexp(un) * fj6
    End Function FJ6

    Real(dp) Function FJ9(a11, a12, a13, a21, a22, a23, a31, a32, a33)
        Implicit None

        Real(dp), Intent(In) :: a11, a12, a13, a21, a22, a23, a31, a32, a33

        Integer :: n, m, k
        Real(dp) :: x, y, a, d
        Real(dp), Dimension(3) :: up, down

        fj9 = 0.0_dp

        Call DELFJ(a11, a12, a13, n)
        If (n /= 1) Return
        Call DELFJ(a21, a22, a23, n)
        If (n /= 1) Return
        Call DELFJ(a31, a32, a33, n)
        If (n /= 1) Return
        Call DELFJ(a11, a21, a31, n)
        If (n /= 1) Return
        Call DELFJ(a12, a22, a32, n)
        If (n /= 1) Return
        Call DELFJ(a13, a23, a33, n)
        If (n /= 1) Return

        up(1) = a11 + a33
        down(1) = dabs(a11 - a33)
        up(2) = a32 + a21
        down(2) = dabs(a32 - a21)
        up(3) = a12 + a23
        down(3) = dabs(a12 - a23)

        x = Maxval(up)
        y = Minval(down)
        IF (x <= y) Return

        d = x - y + 1.0_dp
        n = Int(d)
        If (dabs(d - Real(n, dp)) > 0.00001_dp) Return

        Do k = 1, n
            a = y + Real(k - 1, dp)
            m = Int(2 * a + 0.00001_dp)
            fj9 = fj9 + (-1)**m * (m + 1) * &
                   FJ6(a11, a21, a31, a32, a33, a) * &
                   FJ6(a12, a22, a32, a21, a, a23) * &
                   FJ6(a13, a23, a33, a, a11, a12)
        End Do
    End Function FJ9

    Subroutine DELFJ(a, b, c, n)
        Implicit None

        Real(dp), Intent(In) :: a, b, c
        Integer, Intent(Out) :: n
        Real(dp) :: x, y

        x = a + b
        y = dabs(a - b)

        If (x < c) Then
            n = 2
        Else If (y < c) Then
            n = 1
        Else
            n = 2
        End If

    End Subroutine DELFJ

    Subroutine NGFJ(j, m, n, all)
        Implicit None

        Integer, Intent(In) :: j, m, n
        Real(dp), Intent(InOut) :: all
        Integer  :: l
        Real(dp) :: al

        If (m < n) Then
            Do l = m + 1, n
                al = log(1.0_dp * l)
                all = all + j * al
            End Do
        End If
        
    End Subroutine NGFJ

    Subroutine swap(m1, m2)
        Implicit None

        Integer, Intent(InOut) :: m1, m2
        Integer :: n

        n = m1
        m1 = m2
        m2 = n
        
    End Subroutine swap

    Subroutine sortArray(arr, n)
        Implicit None
        
        Integer, Intent(In) :: n
        Integer, Dimension(n), Intent(InOut) :: arr

        Integer :: i, j
        
        Do i = 1, n-1
            Do j = i+1, n
                If (arr(i) > arr(j)) Then
                    Call swap(arr(i), arr(j))
                End If
              End Do
        End Do

    End Subroutine sortArray

End Module