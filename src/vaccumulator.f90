Module vaccumulator

    Use conf_variables
    
    Implicit None
    
    Private
    
    Public :: IVAccumulator, RVAccumulator
    Public :: IVAccumulatorInit, IVAccumulatorAdd, IVAccumulatorCopy, IVAccumulatorReset
    Public :: RVAccumulatorInit, RVAccumulatorAdd, RVAccumulatorCopy, RVAccumulatorReset
    Public :: IVAccumulatorContinue, RVAccumulatorContinue
    
    Type IVAccumulator
        Integer, Allocatable, Dimension(:)  :: vAccum
        Integer(kind=int64)                 :: vLen, vSize, vGrowBy
    End Type IVAccumulator
    
    Type RVAccumulator
        Real(kind=Hamil%knd), Allocatable, Dimension(:) :: vAccum
        Integer(kind=int64)                 :: vLen, vSize, vGrowBy
    End Type RVAccumulator

  Contains
  
    Subroutine IVAccumulatorInit(this, growBy)
        Implicit None
        
        Type(IVAccumulator), Intent(InOut)  :: this
        Integer, Intent(In)                 :: growBy
        
        this%vLen = 0
        If (growBy <= 0) then
            this%vGrowBy = 1000
        Else
            this%vGrowBy = growBy
        End If
    End Subroutine

    Subroutine IVAccumulatorContinue(this, growBy)
        Implicit None
        
        Type(IVAccumulator), Intent(InOut)  :: this
        Integer, Intent(In)                 :: growBy
        
        this%vLen = size(this%vAccum)
        this%vSize = size(this%vAccum)
        If (growBy <= 0) then
            this%vGrowBy = 1000
        Else
            this%vGrowBy = growBy
        End If
    End Subroutine
    
    Subroutine IVAccumulatorAdd(this, iVal)
        Implicit None
        
        Type(IVAccumulator), Intent(InOut)  :: this
        Integer, Intent(In)                 :: iVal
        Integer, Dimension(:), Allocatable  :: vTmp
        
        If (.not. Allocated(this%vAccum)) then
            ! Initial allocation of vAccum:
            Allocate(this%vAccum(this%vGrowBy))
            this%vLen = 0
            this%vSize = this%vGrowBy
        Else If (this%vLen == this%vSize) then
            ! Expand vAccum to the next larger size:
            Allocate(vTmp(this%vSize + this%vGrowBy))
            vTmp(1:this%vSize) = this%vAccum(1:this%vSize)
            Call Move_Alloc(vTmp, this%vAccum)
            this%vSize = this%vSize + this%vGrowBy
        End If
    
        this%vLen = this%vLen + 1_int64
        this%vAccum(this%vLen) = iVal
    End Subroutine
    
    Subroutine IVAccumulatorCopy(this, outV, outVLen)
        Implicit None
        
        Type(IVAccumulator), Intent(InOut)              :: this
        Integer, Dimension(:), Allocatable, Intent(Out) :: outV
        Integer, Intent(Out)                            :: outVLen
        
        outVLen = this%vLen
        If (this%vLen > 0) then
            Allocate(outV(this%vLen))
            outV = this%vAccum(1:this%vLen)
        End If
    End Subroutine
    
    Subroutine IVAccumulatorReset(this)
        Implicit None
        
        Type(IVAccumulator), Intent(InOut)  :: this
        
        If (Allocated(this%vAccum)) Deallocate(this%vAccum)
        this%vSize = 0
        this%vLen = 0
    End Subroutine

    Subroutine RVAccumulatorInit(this, growBy)
        Implicit None
        
        Type(RVAccumulator), Intent(InOut)  :: this
        Integer, Intent(In)                 :: growBy
        
        this%vLen = 0
        If (growBy <= 0) then
            this%vGrowBy = 1000
        Else
            this%vGrowBy = growBy
        End If
    End Subroutine

    Subroutine RVAccumulatorContinue(this, growBy)
        Implicit None
        
        Type(RVAccumulator), Intent(InOut)  :: this
        Integer, Intent(In)                 :: growBy
        
        this%vLen = size(this%vAccum)
        this%vSize = size(this%vAccum)
        If (growBy <= 0) then
            this%vGrowBy = 1000
        Else
            this%vGrowBy = growBy
        End If
    End Subroutine
    
    Subroutine RVAccumulatorAdd(this, rVal)
        Implicit None
        
        Type(RVAccumulator), Intent(InOut)  :: this
        Real(kind=Hamil%knd), Intent(In)                :: rVal
        Real(kind=Hamil%knd), Dimension(:), Allocatable :: vTmp
        
        If (.not. Allocated(this%vAccum)) then
            ! Initial allocation of vAccum:
            Allocate(this%vAccum(this%vGrowBy))
            this%vLen = 0
            this%vSize = this%vGrowBy
        Else If (this%vLen == this%vSize) then
            ! Expand vAccum to the next larger size:
            Allocate(vTmp(this%vSize + this%vGrowBy))
            vTmp(1:this%vSize) = this%vAccum(1:this%vSize)
            Call Move_Alloc(vTmp, this%vAccum)
            this%vSize = this%vSize + this%vGrowBy
        End If
    
        this%vLen = this%vLen + 1_int64
        this%vAccum(this%vLen) = rVal
    End Subroutine
    
    Subroutine RVAccumulatorCopy(this, outV, outVLen)
        Implicit None
        
        Type(RVAccumulator), Intent(InOut)                  :: this
        Real(kind=Hamil%knd), Dimension(:), Allocatable, Intent(Out)    :: outV
        Integer, Intent(Out)                                :: outVLen
        
        outVLen = this%vLen
        If (this%vLen > 0) then
            Allocate(outV(this%vLen))
            outV = this%vAccum(1:this%vLen)
        End If
    End Subroutine
    
    Subroutine RVAccumulatorReset(this)
        Implicit None
        
        Type(RVAccumulator), Intent(InOut)  :: this
        
        If (Allocated(this%vAccum)) Deallocate(this%vAccum)
        this%vSize = 0
        this%vLen = 0
    End Subroutine

End Module
