Module env_var

    Use, Intrinsic :: iso_fortran_env, Only : sp => real32, dp => real64, Int64

    Implicit None
    
    Private

    Public :: GetEnvIsPresent, GetEnvLogical, GetEnvInteger, GetEnvInteger64

  Contains

    Logical Function GetEnvIsPresent(varName)
        Implicit none
        
        Character(Len=*), Intent(In)    :: varName
        
        Character(Len=24)               :: varValue
    
        Call GetEnv(varName, varValue)
        GetEnvIsPresent = (Len_Trim(varValue) > 0)
        Return
    End Function
    
    Logical Function GetEnvLogical(varName, defaultValue)
        Implicit none
        
        Character(Len=*), Intent(In)    :: varName
        Logical, Intent(In), Optional   :: defaultValue
        
        Character(Len=24)               :: varValue
        Integer                         :: iValue, convStat
        
        GetEnvLogical = defaultValue
        Call GetEnv(varName, varValue)
        If (Len_Trim(varValue) > 0) then
            If (varValue == 'y' .or. varValue == 'Y' .or. varValue == 'yes' .or. varValue == 'YES' .or. varValue == 't' .or. varValue == 'T' .or. varValue == 'true' .or. varValue == 'TRUE') then
                GetEnvLogical = .true.
            Else If (varValue == 'n' .or. varValue == 'N' .or. varValue == 'no' .or. varValue == 'NO' .or. varValue == 'f' .or. varValue == 'F' .or. varValue == 'false' .or. varValue == 'FALSE') then
                GetEnvLogical = .false.
            Else
                Read(varValue,*,iostat=convStat) iValue
                If (convStat == 0) then
                    If (iValue == 0) then
                        GetEnvLogical = .false.
                    Else
                        GetEnvLogical = .true.
                    EndIf
                End If
            End If
        End If
        Return
    End Function
    
    Integer Function GetEnvInteger(varName, defaultValue)
        Implicit none
        
        Character(Len=*), Intent(In)    :: varName
        Integer, Intent(In), Optional   :: defaultValue
        
        Character(Len=24)               :: varValue
        Integer                         :: iValue, convStat
        
        GetEnvInteger = defaultValue
        Call GetEnv(varName, varValue)
        If (Len_Trim(varValue) > 0) then
            Read(varValue,*,iostat=convStat) iValue
            If (convStat == 0) GetEnvInteger = iValue
        End If
        Return
    End Function

    Integer(Kind=Int64) Function GetEnvInteger64(varName, defaultValue)
        Implicit none

        Character(Len=*), Intent(In)    :: varName
        Integer(Kind=Int64), Intent(In), Optional   :: defaultValue

        Character(Len=48)               :: varValue
        Integer                         :: convStat
        Integer(Kind=Int64)             :: iValue

        GetEnvInteger64 = defaultValue
        Call GetEnv(varName, varValue)
        If (Len_Trim(varValue) > 0) then
            Read(varValue,*,iostat=convStat) iValue
            If (convStat == 0) GetEnvInteger64 = iValue
        End If
        Return
    End Function
    
End Module
