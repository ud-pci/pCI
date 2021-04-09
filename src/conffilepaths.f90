!
! Module:   ConfFilePaths
! Usage:    use conffilepaths
!
! This Fortran module manages a set of single-access file types with which the conf program performs
! i/o.  Each type of file has a default file path (relative to the working directory) that was the
! legacy name for it.  For each file type, an environment variable may be set to override the
! default path.  Finally, the command argument list is also consulted for values (that override
! both the default and a value coming from an environment variable).
!
! The 'ConfHamiltonianFileIndex' file type is slightly unique in that it can be split between multiple
! files.  The default path contains the token sequence "###" to be overridden with a zero-padded index
! value that defaults to 0.  A path coming from the environment variable or a command argument MUST
! also contain one or more "#" characters to be acceptable.
!
Module ConfFilePaths

    Implicit None

    Private

    !
    ! We made this module private by default, so expose those items that should be accessible by external
    ! consumers of this module.
    !
    ! If you add new file types, be sure to expose their index variables here.
    !
    Public :: ConfFileInit
    Public :: ConfFileGetIsOpen, ConfFileGetPathLen, ConfFileGetPath, ConfFileSetPath, ConfFileOpen, ConfFileClose, ConfFileGetOpenUnit
    Public :: ConfInputFileIndex, ConfResultsFileIndex, ConfRadialFnsFileIndex, ConfRadialIntsFileIndex
    Public :: ConfHamiltonianFileIndex, ConfDeterminantsFileIndex, ConfEigenvectorsFileIndex, ConfFormJRestartFileIndex
    Public :: ConfGauntFileIndex, ConfProgressFileIndex, ConfSelfEnergyFileIndex, ConfScreeningRadialIntsFileIndex
    Public :: ConfControlFileIndex

    !
    ! For each type of file we handle, an integer index must be assigned.  This integer indexes that file
    ! type in the ConfFiles and ConfFilesInfo arrays.
    !
    ! When adding new file types, be sure to use the index immediatley following the last-defined index
    ! and don't forget to increment ConfMaxFileIndex.
    !
    Integer, Parameter  :: ConfInputFileIndex = 1
    Integer, Parameter  :: ConfResultsFileIndex = 2
    Integer, Parameter  :: ConfRadialFnsFileIndex = 3
    Integer, Parameter  :: ConfRadialIntsFileIndex = 4
    Integer, Parameter  :: ConfHamiltonianFileIndex = 5
    Integer, Parameter  :: ConfDeterminantsFileIndex = 6
    Integer, Parameter  :: ConfEigenvectorsFileIndex = 7
    Integer, Parameter  :: ConfFormJRestartFileIndex = 8
    Integer, Parameter  :: ConfGauntFileIndex = 9
    Integer, Parameter  :: ConfProgressFileIndex = 10
    Integer, Parameter  :: ConfSelfEnergyFileIndex = 11
    Integer, Parameter  :: ConfScreeningRadialIntsFileIndex = 12
    Integer, Parameter  :: ConfControlFileIndex = 13
    
    Integer, Parameter  :: ConfMaxFileIndex = 13
    
    !
    ! The ConfFile type encapsulates optional file names and a unit number (once the file is
    ! actually opened).
    !
    ! A single array of this type holds state information for each of the file types
    ! we manage; only a single file of a type can be open at once.
    !
    Type ConfFile
        Character(Len=:), Allocatable   :: filePath, hijFilePath
        Integer                         :: ioUnit
    End Type ConfFile
    Type(ConfFile)                      :: ConfFiles(ConfMaxFileIndex)
    
    !
    ! Logical array that determines which of the files we manage is actually opened at this
    ! time.
    !
    Logical :: ConfFileIsOpen(ConfMaxFileIndex) = .false.
    
    !
    ! For each type of file we handle, define a set of strings for:
    !
    !     - default file name (e.g. CONF.INP)
    !     - environment variable that overrides default (e.g. CONF_FILE_INP)
    !     - command-line flag that overrides environment variable (e.g. --file-inp)
    !
    ! When adding new file types be sure to add a record to the ConfFilesInfo array.
    !
    Type ConfFileInfo
        Character(Len=16)           :: defaultName
        Character(Len=24)           :: envVar
        Character(Len=16)           :: cliFlag
    End Type ConfFileInfo
    Type(ConfFileInfo), Parameter   :: ConfFilesInfo(ConfMaxFileIndex) = (/                                          &
                                            ConfFileInfo('CONF.INP',     'CONF_FILE_INP',         '--file-inp'),     &
                                            ConfFileInfo('CONF.RES',     'CONF_FILE_RES',         '--file-res'),     &
                                            ConfFileInfo('CONF.DAT',     'CONF_FILE_DAT',         '--file-dat'),     &
                                            ConfFileInfo('CONF.INT',     'CONF_FILE_INT',         '--file-int'),     &
                                            ConfFileInfo('CONF_###.HIJ', 'CONF_FILE_HIJ_PATTERN', '--file-hij'),     &
                                            ConfFileInfo('CONF.DET',     'CONF_FILE_DET',         '--file-det'),     &
                                            ConfFileInfo('CONF.XIJ',     'CONF_FILE_XIJ',         '--file-xij'),     &
                                            ConfFileInfo('CONF.JJJ',     'CONF_FILE_JJJ',         '--file-jjj'),     &
                                            ConfFileInfo('CONF.GNT',     'CONF_FILE_GNT',         '--file-gnt'),     &
                                            ConfFileInfo('CONF.PROGRESS','CONF_FILE_PROGRESS',    '--file-progress'),&
                                            ConfFileInfo('SGC.CON',      'CONF_FILE_SGC_CON',     '--file-sgc-con'), &
                                            ConfFileInfo('SCRC.CON',     'CONF_FILE_SCRC_CON',    '--file-scrc-con'),&
                                            ConfFileInfo('c.in',         'CONF_FILE_CONTROL',     '--file-control')  &
                                        /)
  Contains
  
    Subroutine ConfFileInit(shouldSummarize)
        Implicit None
        
        Logical, Intent(In), Optional   :: shouldSummarize
        
        Character(Len=4096)             :: thePath
        Integer                         :: i, iMax, j, thePathLen, status
        
        ! Check for environment variables first:
        Do i = 1, ConfMaxFileIndex
            Call GetEnv(ConfFilesInfo(i)%envVar, thePath)
            If (Len_Trim(thePath) > 0 ) Call ConfFileSetPath(i, Trim(thePath))
        End Do
        
        ! Now check for CLI arguments:
#ifdef HAVE_GET_COMMAND_ARGUMENT
        i = 1
        iMax = Command_Argument_Count()
        Do While (i <= iMax)
            Call Get_Command_Argument(i, thePath, thePathLen, status)
            If (status /= 0) then
                If (status == -1) then
                    Write(0,fmt='(A,X,I,X,A)') 'ERROR: Argument at index', i, 'exceeded size of copy buffer'
                    Error Stop
                End If
                Write(0,fmt='(A,X,I)') 'ERROR: Failed to copy argument at index', i
                Error Stop
            End If
            If (thePathLen > 7) then
                If (thePath(1:7) == '--file-') then
                    Do j = 1, ConfMaxFileIndex
                        If ((thePathLen >= Len_Trim(ConfFilesInfo(j)%cliFlag)) .and. (thePath(1:Len_Trim(ConfFilesInfo(j)%cliFlag)) == ConfFilesInfo(j)%cliFlag)) then
                            thePath = thePath(Len_Trim(ConfFilesInfo(j)%cliFlag) + 1:thePathLen)
                            thePathLen = Len_Trim(thePath)
                            If (thePathLen > 0) then
                                If (thePath(1:1) == '=') then
                                    thePath = thePath(2:thePathLen)
                                    If (Len_Trim(thePath) > 0) then
                                        Call ConfFileSetPath(j, Trim(thePath))
                                    Else
                                        Call Get_Command_Argument(i, thePath)
                                        Write(0,fmt='(A,X,A)') 'ERROR: Empty path provided:', Trim(thePath)
                                        Error Stop
                                    End If
                                End If
                            Else
                                If (i == iMax) then
                                    Call Get_Command_Argument(i, thePath)
                                    Write(0,fmt='(A,X,A,X,A)') 'ERROR: No path specified with', Trim(thePath), 'flag'
                                    Error Stop
                                End If
                                i = i + 1
                                Call Get_Command_Argument(i, thePath, thePathLen, status)
                                If (status /= 0) then
                                    If (status == -1) then
                                        Write(0,fmt='(A,X,I,X,A)') 'ERROR: File argument at index', i, 'exceeded size of copy buffer'
                                        Error Stop
                                    End If
                                    Write(0,fmt='(A,X,I)') 'ERROR: Failed to copy file argument at index', i
                                    Error Stop
                                End If
                                If (Len_Trim(thePath) > 0) then
                                    Call ConfFileSetPath(j, Trim(thePath))
                                Else
                                    Call Get_Command_Argument(i - 1, thePath)
                                    Write(0,fmt='(A,X,A)') 'ERROR: Empty path provided to', Trim(thePath)
                                    Error Stop
                                End If
                            End If
                        End If
                    End Do
                End If
            End If
            i = i + 1
        End Do
#else
        i = 1
        iMax = IArgc()
        Do While (i <= iMax)
            Call GetArg(i, thePath)
            thePathLen = Len_Trim(thePath)
            If (thePathLen > 7) then
                If (thePath(1:7) == '--file-') then
                    Do j = 1, ConfMaxFileIndex
                        If ((thePathLen >= Len_Trim(ConfFilesInfo(j)%cliFlag)) .and. (thePath(1:Len_Trim(ConfFilesInfo(j)%cliFlag)) == ConfFilesInfo(j)%cliFlag)) then
                            thePath = thePath(Len_Trim(ConfFilesInfo(j)%cliFlag) + 1:thePathLen)
                            thePathLen = Len_Trim(thePath)
                            If (thePathLen > 0) then
                                If (thePath(1:1) == '=') then
                                    thePath = thePath(2:thePathLen)
                                    If (Len_Trim(thePath) > 0) then
                                        Call ConfFileSetPath(j, Trim(thePath))
                                    Else
                                        Call GetArg(i, thePath)
                                        Write(0,fmt='(A,X,A)') 'ERROR: Empty path provided:', Trim(thePath)
                                        Error Stop
                                    End If
                                End If
                            Else
                                If (i == iMax) then
                                    Call GetArg(i, thePath)
                                    Write(0,fmt='(A,X,A,X,A)') 'ERROR: No path specified with', Trim(thePath), 'flag'
                                    Error Stop
                                End If
                                i = i + 1
                                Call GetArg(i, thePath)
                                If (Len_Trim(thePath) > 0) then
                                    Call ConfFileSetPath(j, Trim(thePath))
                                Else
                                    Call GetArg(i - 1, thePath)
                                    Write(0,fmt='(A,X,A)') 'ERROR: Empty path provided to', Trim(thePath)
                                    Error Stop
                                End If
                            End If
                        End If
                    End Do
                End If
            End If
            i = i + 1
        End Do
#endif
        If (Present(shouldSummarize)) then
            If (shouldSummarize) then
                Write(*,fmt='(A)') 'Configured file paths:'
                Do i = 1, ConfMaxFileIndex
                    Call ConfFileGetPath(i, thePath)
                    Write(*,fmt='("    ",A," = ",A)') ConfFilesInfo(i)%envVar, Trim(thePath)
                End Do
            End If
        End If
    End Subroutine
  
    Logical Function ConfFileGetIsOpen(filePathIndex)
        Implicit None
        
        Integer, Intent(In)             :: filePathIndex
        
        If (filePathIndex < 1 .or. filePathIndex > ConfMaxFileIndex) then
            Write(0, fmt='(A,X,I2)') 'ERROR: Invalid filePathIndex to ConfFileGetPath():', filePathIndex
            Error Stop
        End If
        ConfFileGetIsOpen = ConfFileIsOpen(filePathIndex)
        Return
    End Function

    Integer Function ConfFileGetPathLen(filePathIndex)
        Implicit None
        
        Integer, Intent(In)             :: filePathIndex
        
        If (filePathIndex < 1 .or. filePathIndex > ConfMaxFileIndex) then
            Write(0, fmt='(A,X,I2)') 'ERROR: Invalid filePathIndex to ConfFileGetPath():', filePathIndex
            Error Stop
        End If
        If (Allocated(ConfFiles(filePathIndex)%hijFilePath)) then
            ConfFileGetPathLen = Len(ConfFiles(filePathIndex)%hijFilePath)
        Else If (Allocated(ConfFiles(filePathIndex)%filePath)) then
            ConfFileGetPathLen = Len(ConfFiles(filePathIndex)%filePath)
        Else
            ConfFileGetPathLen = Len(ConfFilesInfo(filePathIndex)%defaultName)
        End If
    End Function

    Subroutine ConfFileGetPath(filePathIndex, outFilePath)
        Implicit None
        
        Integer, Intent(In)             :: filePathIndex
        Character(Len=*), Intent(InOut) :: outFilePath
        
        If (filePathIndex < 1 .or. filePathIndex > ConfMaxFileIndex) then
            Write(0, fmt='(A,X,I2)') 'ERROR: Invalid filePathIndex to ConfFileGetPath():', filePathIndex
            Error Stop
        End If
        If (Allocated(ConfFiles(filePathIndex)%hijFilePath)) then
            outFilePath = ConfFiles(filePathIndex)%hijFilePath
        Else If (Allocated(ConfFiles(filePathIndex)%filePath)) then
            outFilePath = ConfFiles(filePathIndex)%filePath
        Else
            outFilePath = ConfFilesInfo(filePathIndex)%defaultName
        End If
    End Subroutine

    Subroutine ConfFileGetConfiguredPath(filePathIndex, outFilePath)
        Implicit None
        
        Integer, Intent(In)             :: filePathIndex
        Character(Len=*), Intent(InOut) :: outFilePath
        
        If (Allocated(ConfFiles(filePathIndex)%filePath)) then
            outFilePath = ConfFiles(filePathIndex)%filePath
        Else
            outFilePath = ConfFilesInfo(filePathIndex)%defaultName
        End If
    End Subroutine

    Subroutine ConfFileSetPath(filePathIndex, inFilePath)
        Implicit None
        
        Integer, Intent(In)             :: filePathIndex
        Character(Len=*), Intent(In)    :: inFilePath
        
        Integer                         :: inFilePathLen
        
        If (filePathIndex < 1 .or. filePathIndex > ConfMaxFileIndex) then
            Write(0, fmt='(A,X,I2)') 'ERROR: Invalid filePathIndex to ConfFileSetPath():', filePathIndex
            Error Stop
        End If
        If (Allocated(ConfFiles(filePathIndex)%filePath)) then
            Write(0, fmt='(A,X,I2,X,A)') 'ERROR: File at index', filePathIndex, 'already configured or open'
            Error Stop
        End If
        
        ! Make sure the Hamiltonian inFilePath is a template:
        If (filePathIndex == ConfHamiltonianFileIndex) then
            inFilePathLen = Index(inFilePath, '#')
            If (inFilePathLen == 0) then
                Write(0,fmt='(A,X,A)') 'ERROR: Invalid Hamiltonian file name (lacks "#" tokens):', Trim(inFilePath)
                Error Stop
            End If
        End If        
        
        ! Allocate a string and copy inFilePath to it:
        inFilePathLen = Len(inFilePath)
        Allocate(Character(Len=inFilePathLen) :: ConfFiles(filePathIndex)%filePath)
        ConfFiles(filePathIndex)%filePath = inFilePath
    End Subroutine
    
    Integer Function ConfFileOpen(filePathIndex, hijIndex, unit, form, status, access, recl, stat)
        Implicit None
        
        Integer, Intent(In)                     :: filePathIndex
        Integer, Intent(In), Optional           :: hijIndex
        Integer, Intent(In)                     :: unit
        Character(Len=*), Intent(In), Optional  :: status, form, access
        Integer, Intent(InOut), Optional        :: stat
        Integer, Intent(In), Optional           :: recl
        
        Integer                                 :: ioStat, ioHijIndex, fnStat
        Integer                                 :: hijNameLen, cIndex
        Character(Len=:), Allocatable           :: hijName
        Character, Dimension(10), Parameter     :: base10Digits = (/ '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' /)
        
        If (filePathIndex < 1 .or. filePathIndex > ConfMaxFileIndex) then
            Write(0, fmt='(A,X,I2)') 'ERROR: Invalid filePathIndex to ConfFileSetPath():', filePathIndex
            Error Stop
        End If
        If (ConfFileIsOpen(filePathIndex)) then
            Write(0, fmt='(A,X,I2,X,A)') 'ERROR: Conf file is already open (index', filePathIndex, ')'
            Error Stop
        End If
        
        ! Default return value:
        ConfFileOpen = -1
        
        If (filePathIndex == ConfHamiltonianFileIndex) then
            If (Present(hijIndex)) then
                ioHijIndex = hijIndex
            Else
                ioHijIndex = 0
            End If
            
            ! Copy the template into place:
            hijNameLen = ConfFileGetPathLen(filePathIndex)
            Allocate(Character(Len=hijNameLen) :: ConfFiles(filePathIndex)%hijFilePath)
            Call ConfFileGetConfiguredPath(filePathIndex, ConfFiles(filePathIndex)%hijFilePath)
            
            ! Now replace any '#' characters with the hijIndex digits:
            cIndex = Index(ConfFiles(filePathIndex)%hijFilePath, '#', .true.)
            If (cIndex == 0) then
                Write(0,fmt='(A,X,A)') 'ERROR: Invalid Hamiltonian file name (lacks "#" tokens):', ConfFiles(filePathIndex)%hijFilePath
                Error Stop
            End If
            Do While (cIndex > 0)
                ConfFiles(filePathIndex)%hijFilePath = ConfFiles(filePathIndex)%hijFilePath(1:cIndex - 1) // base10Digits(1 + Mod(ioHijIndex, 10)) // ConfFiles(filePathIndex)%hijFilePath(cIndex + 1:hijNameLen)
                ioHijIndex = ioHijIndex / 10
                cIndex = Index(ConfFiles(filePathIndex)%hijFilePath, '#', .true.)
            End Do
            
            ! Filename is ready, get it opened:
            If (Present(status)) then
                If (Present(form)) then
                    If (Present(access)) then
                        If (Present(recl)) then
                            fnStat = ConfFileOpenFSAR(unit=unit, file=ConfFiles(filePathIndex)%hijFilePath, status=status, form=form, access=access, recl=recl, iostat=ioStat)
                        Else
                            fnStat = ConfFileOpenFSA(unit=unit, file=ConfFiles(filePathIndex)%hijFilePath, status=status, form=form, access=access, iostat=ioStat)
                        End If
                    Else
                        If (Present(recl)) then
                            fnStat = ConfFileOpenFSR(unit=unit, file=ConfFiles(filePathIndex)%hijFilePath, status=status, form=form, recl=recl, iostat=ioStat)
                        Else
                            fnStat = ConfFileOpenFS(unit=unit, file=ConfFiles(filePathIndex)%hijFilePath, status=status, form=form, iostat=ioStat)
                        End If
                    End If
                Else
                    If (Present(access)) then
                        If (Present(recl)) then
                            fnStat = ConfFileOpenSAR(unit=unit, file=ConfFiles(filePathIndex)%hijFilePath, status=status, access=access, recl=recl, iostat=ioStat)
                        Else
                            fnStat = ConfFileOpenSA(unit=unit, file=ConfFiles(filePathIndex)%hijFilePath, status=status, access=access, iostat=ioStat)
                        End If
                    Else
                        If (Present(recl)) then
                            fnStat = ConfFileOpenSR(unit=unit, file=ConfFiles(filePathIndex)%hijFilePath, status=status, recl=recl, iostat=ioStat)
                        Else
                            fnStat = ConfFileOpenS(unit=unit, file=ConfFiles(filePathIndex)%hijFilePath, status=status, iostat=ioStat)
                        End If
                    End If
                End If
            Else
                If (Present(form)) then
                    If (Present(access)) then
                        If (Present(recl)) then
                            fnStat = ConfFileOpenFAR(unit=unit, file=ConfFiles(filePathIndex)%hijFilePath, form=form, access=access, recl=recl, iostat=ioStat)
                        Else
                            fnStat = ConfFileOpenFA(unit=unit, file=ConfFiles(filePathIndex)%hijFilePath, form=form, access=access, iostat=ioStat)
                        End If
                    Else
                        If (Present(recl)) then
                            fnStat = ConfFileOpenFR(unit=unit, file=ConfFiles(filePathIndex)%hijFilePath, form=form, recl=recl, iostat=ioStat)
                        Else
                            fnStat = ConfFileOpenF(unit=unit, file=ConfFiles(filePathIndex)%hijFilePath, form=form, iostat=ioStat)
                        End If
                    End If
                Else
                    If (Present(access)) then
                        If (Present(recl)) then
                            fnStat = ConfFileOpenAR(unit=unit, file=ConfFiles(filePathIndex)%hijFilePath, access=access, recl=recl, iostat=ioStat)
                        Else
                            fnStat = ConfFileOpenA(unit=unit, file=ConfFiles(filePathIndex)%hijFilePath, access=access, iostat=ioStat)
                        End If
                    Else
                        If (Present(recl)) then
                            fnStat = ConfFileOpenR(unit=unit, file=ConfFiles(filePathIndex)%hijFilePath, recl=recl, iostat=ioStat)
                        Else
                            fnStat = 0
                            Open(unit=unit, file=ConfFiles(filePathIndex)%hijFilePath, iostat=ioStat)
                            If (ioStat /= 0) fnStat = -1
                        End If
                    End If
                End If
            End If
            If (ioStat /= 0) Deallocate(ConfFiles(filePathIndex)%hijFilePath)
        Else
            If (Allocated(ConfFiles(filePathIndex)%filePath)) then
                If (Present(status)) then
                    If (Present(form)) then
                        If (Present(access)) then
                            If (Present(recl)) then
                                fnStat = ConfFileOpenFSAR(unit=unit, file=ConfFiles(filePathIndex)%filePath, status=status, form=form, access=access, recl=recl, iostat=ioStat)
                            Else
                                fnStat = ConfFileOpenFSA(unit=unit, file=ConfFiles(filePathIndex)%filePath, status=status, form=form, access=access, iostat=ioStat)
                            End If
                        Else
                            If (Present(recl)) then
                                fnStat = ConfFileOpenFSR(unit=unit, file=ConfFiles(filePathIndex)%filePath, status=status, form=form, recl=recl, iostat=ioStat)
                            Else
                                fnStat = ConfFileOpenFS(unit=unit, file=ConfFiles(filePathIndex)%filePath, status=status, form=form, iostat=ioStat)
                            End If
                        End If
                    Else
                        If (Present(access)) then
                            If (Present(recl)) then
                                fnStat = ConfFileOpenSAR(unit=unit, file=ConfFiles(filePathIndex)%filePath, status=status, access=access, recl=recl, iostat=ioStat)
                            Else
                                fnStat = ConfFileOpenSA(unit=unit, file=ConfFiles(filePathIndex)%filePath, status=status, access=access, iostat=ioStat)
                            End If
                        Else
                            If (Present(recl)) then
                                fnStat = ConfFileOpenSR(unit=unit, file=ConfFiles(filePathIndex)%filePath, status=status, recl=recl, iostat=ioStat)
                            Else
                                fnStat = ConfFileOpenS(unit=unit, file=ConfFiles(filePathIndex)%filePath, status=status, iostat=ioStat)
                            End If
                        End If
                    End If
                Else
                    If (Present(form)) then
                        If (Present(access)) then
                            If (Present(recl)) then
                                fnStat = ConfFileOpenFAR(unit=unit, file=ConfFiles(filePathIndex)%filePath, form=form, access=access, recl=recl, iostat=ioStat)
                            Else
                                fnStat = ConfFileOpenFA(unit=unit, file=ConfFiles(filePathIndex)%filePath, form=form, access=access, iostat=ioStat)
                            End If
                        Else
                            If (Present(recl)) then
                                fnStat = ConfFileOpenFR(unit=unit, file=ConfFiles(filePathIndex)%filePath, form=form, recl=recl, iostat=ioStat)
                            Else
                                fnStat = ConfFileOpenF(unit=unit, file=ConfFiles(filePathIndex)%filePath, form=form, iostat=ioStat)
                            End If
                        End If
                    Else
                        If (Present(access)) then
                            If (Present(recl)) then
                                fnStat = ConfFileOpenAR(unit=unit, file=ConfFiles(filePathIndex)%filePath, access=access, recl=recl, iostat=ioStat)
                            Else
                                fnStat = ConfFileOpenA(unit=unit, file=ConfFiles(filePathIndex)%filePath, access=access, iostat=ioStat)
                            End If
                        Else
                            If (Present(recl)) then
                                fnStat = ConfFileOpenR(unit=unit, file=ConfFiles(filePathIndex)%filePath, recl=recl, iostat=ioStat)
                            Else
                                fnStat = 0
                                Open(unit=unit, file=ConfFiles(filePathIndex)%filePath, iostat=ioStat)
                                If (ioStat /= 0) fnStat = -1
                            End If
                        End If
                    End If
                End If
            Else
                If (Present(status)) then
                    If (Present(form)) then
                        If (Present(access)) then
                            If (Present(recl)) then
                                fnStat = ConfFileOpenFSAR(unit=unit, file=ConfFilesInfo(filePathIndex)%defaultName, status=status, form=form, access=access, recl=recl, iostat=ioStat)
                            Else
                                fnStat = ConfFileOpenFSA(unit=unit, file=ConfFilesInfo(filePathIndex)%defaultName, status=status, form=form, access=access, iostat=ioStat)
                            End If
                        Else
                            If (Present(recl)) then
                                fnStat = ConfFileOpenFSR(unit=unit, file=ConfFilesInfo(filePathIndex)%defaultName, status=status, form=form, recl=recl, iostat=ioStat)
                            Else
                                fnStat = ConfFileOpenFS(unit=unit, file=ConfFilesInfo(filePathIndex)%defaultName, status=status, form=form, iostat=ioStat)
                            End If
                        End If
                    Else
                        If (Present(access)) then
                            If (Present(recl)) then
                                fnStat = ConfFileOpenSAR(unit=unit, file=ConfFilesInfo(filePathIndex)%defaultName, status=status, access=access, recl=recl, iostat=ioStat)
                            Else
                                fnStat = ConfFileOpenSA(unit=unit, file=ConfFilesInfo(filePathIndex)%defaultName, status=status, access=access, iostat=ioStat)
                            End If
                        Else
                            If (Present(recl)) then
                                fnStat = ConfFileOpenSR(unit=unit, file=ConfFilesInfo(filePathIndex)%defaultName, status=status, recl=recl, iostat=ioStat)
                            Else
                                fnStat = ConfFileOpenS(unit=unit, file=ConfFilesInfo(filePathIndex)%defaultName, status=status, iostat=ioStat)
                            End If
                        End If
                    End If
                Else
                    If (Present(form)) then
                        If (Present(access)) then
                            If (Present(recl)) then
                                fnStat = ConfFileOpenFAR(unit=unit, file=ConfFilesInfo(filePathIndex)%defaultName, form=form, access=access, recl=recl, iostat=ioStat)
                            Else
                                fnStat = ConfFileOpenFA(unit=unit, file=ConfFilesInfo(filePathIndex)%defaultName, form=form, access=access, iostat=ioStat)
                            End If
                        Else
                            If (Present(recl)) then
                                fnStat = ConfFileOpenFR(unit=unit, file=ConfFilesInfo(filePathIndex)%defaultName, form=form, recl=recl, iostat=ioStat)
                            Else
                                fnStat = ConfFileOpenF(unit=unit, file=ConfFilesInfo(filePathIndex)%defaultName, form=form, iostat=ioStat)
                            End If
                        End If
                    Else
                        If (Present(access)) then
                            If (Present(recl)) then
                                fnStat = ConfFileOpenAR(unit=unit, file=ConfFilesInfo(filePathIndex)%defaultName, access=access, recl=recl, iostat=ioStat)
                            Else
                                fnStat = ConfFileOpenA(unit=unit, file=ConfFilesInfo(filePathIndex)%defaultName, access=access, iostat=ioStat)
                            End If
                        Else
                            If (Present(recl)) then
                                fnStat = ConfFileOpenR(unit=unit, file=ConfFilesInfo(filePathIndex)%defaultName, recl=recl, iostat=ioStat)
                            Else
                                fnStat = 0
                                Open(unit=unit, file=ConfFilesInfo(filePathIndex)%defaultName, iostat=ioStat)
                                If (ioStat /= 0) fnStat = -1
                            End If
                        End If
                    End If
                End If
            End If
        End If
        If (ioStat == 0) then
            ConfFileIsOpen(filePathIndex) = .true.
            ConfFiles(filePathIndex)%ioUnit = unit
            ConfFileOpen = unit
        End If
        If (Present(stat)) stat = ioStat
        Return
    End Function
    
    Integer Function ConfFileClose(filePathIndex, status, stat)
        Implicit None
        
        Integer, Intent(In)                     :: filePathIndex
        Integer, Intent(InOut), Optional        :: stat
        Character(Len=*), Intent(In), Optional  :: status
        
        Integer                                 :: ioStat
        
        If (filePathIndex < 1 .or. filePathIndex > ConfMaxFileIndex) then
            Write(0, fmt='(A,X,I2)') 'ERROR: Invalid filePathIndex to ConfFileSetPath():', filePathIndex
            Error Stop
        End If
        If (ConfFileIsOpen(filePathIndex)) then
            ! Default return value:
            ConfFileClose = -1
        
            If (Present(status)) then
                Close(ConfFiles(filePathIndex)%ioUnit, status=status, iostat=ioStat)
            Else
                Close(ConfFiles(filePathIndex)%ioUnit, iostat=ioStat)
            End If
            If (Allocated(ConfFiles(filePathIndex)%hijFilePath)) Deallocate(ConfFiles(filePathIndex)%hijFilePath)
            ConfFiles(filePathIndex)%ioUnit = -1
            ConfFileIsOpen(filePathIndex) = .false.
            If (Present(stat)) stat=ioStat
            If (ioStat == 0) ConfFileClose = 0
        Else
            ConfFileClose = 0
        End If
        Return
    End Function
    
    Integer Function ConfFileGetOpenUnit(filePathIndex)
        Implicit None
        
        Integer, Intent(In)                     :: filePathIndex

        If (filePathIndex < 1 .or. filePathIndex > ConfMaxFileIndex) then
            Write(0, fmt='(A,X,I2)') 'ERROR: Invalid filePathIndex to ConfFileSetPath():', filePathIndex
            Error Stop
        End If
        If (.not. ConfFileIsOpen(filePathIndex)) then
            Write(0, fmt='(A,X,I2,X,A)') 'ERROR: Conf file is not open (index', filePathIndex, ')'
            Error Stop
        End If
        ConfFileGetOpenUnit = ConfFiles(filePathIndex)%ioUnit
        Return
    End Function
    
!
! The following 16 functions Open() an i/o unit with a varying set of arguments; since the Open()
! intrinsic is usually NOT a Fortran function, it doesn't understand Fortran optional arguments
! and thus cannot be called with argument values that were not presented to the function doing the
! Open().  Ugh.
!
    Integer Function ConfFileOpenFSAR(unit, file, form, status, access, recl, iostat)
        Implicit None
        
        Integer, Intent(In)           :: unit
        Character(Len=*), Intent(In)  :: file, status, form, access
        Integer, Intent(In)           :: recl
        Integer, Intent(InOut)        :: iostat
        
        ConfFileOpenFSAR = 0
        Open(unit=unit, file=file, status=status, form=form, access=access, recl=recl, iostat=iostat)
        If (iostat /= 0) ConfFileOpenFSAR = -1
        Return
    End Function
    
    Integer Function ConfFileOpenFSA(unit, file, form, status, access, iostat)
        Implicit None
        
        Integer, Intent(In)           :: unit
        Character(Len=*), Intent(In)  :: file, status, form, access
        Integer, Intent(InOut)        :: iostat
        
        ConfFileOpenFSA = 0
        Open(unit=unit, file=file, status=status, form=form, access=access, iostat=iostat)
        If (iostat /= 0) ConfFileOpenFSA = -1
        Return
    End Function
    
    Integer Function ConfFileOpenFSR(unit, file, form, status, recl, iostat)
        Implicit None
        
        Integer, Intent(In)           :: unit
        Character(Len=*), Intent(In)  :: file, status, form
        Integer, Intent(In)           :: recl
        Integer, Intent(InOut)        :: iostat
        
        ConfFileOpenFSR = 0
        Open(unit=unit, file=file, status=status, form=form, recl=recl, iostat=iostat)
        If (iostat /= 0) ConfFileOpenFSR = -1
        Return
    End Function
    
    Integer Function ConfFileOpenFAR(unit, file, form, access, recl, iostat)
        Implicit None
        
        Integer, Intent(In)           :: unit
        Character(Len=*), Intent(In)  :: file, form, access
        Integer, Intent(In)           :: recl
        Integer, Intent(InOut)        :: iostat
        
        ConfFileOpenFAR = 0
        Open(unit=unit, file=file, form=form, access=access, recl=recl, iostat=iostat)
        If (iostat /= 0) ConfFileOpenFAR = -1
        Return
    End Function
    
    Integer Function ConfFileOpenSAR(unit, file, status, access, recl, iostat)
        Implicit None
        
        Integer, Intent(In)           :: unit
        Character(Len=*), Intent(In)  :: file, status, access
        Integer, Intent(In)           :: recl
        Integer, Intent(InOut)        :: iostat
        
        ConfFileOpenSAR = 0
        Open(unit=unit, file=file, status=status, access=access, recl=recl, iostat=iostat)
        If (iostat /= 0) ConfFileOpenSAR = -1
        Return
    End Function
    
    Integer Function ConfFileOpenFS(unit, file, form, status, iostat)
        Implicit None
        
        Integer, Intent(In)           :: unit
        Character(Len=*), Intent(In)  :: file, status, form
        Integer, Intent(InOut)        :: iostat
        
        ConfFileOpenFS = 0
        Open(unit=unit, file=file, status=status, form=form, iostat=iostat)
        If (iostat /= 0) ConfFileOpenFS = -1
        Return
    End Function
    
    Integer Function ConfFileOpenFA(unit, file, form, access, iostat)
        Implicit None
        
        Integer, Intent(In)           :: unit
        Character(Len=*), Intent(In)  :: file, form, access
        Integer, Intent(InOut)        :: iostat
        
        ConfFileOpenFA = 0
        Open(unit=unit, file=file, form=form, access=access, iostat=iostat)
        If (iostat /= 0) ConfFileOpenFA = -1
        Return
    End Function
    
    Integer Function ConfFileOpenFR(unit, file, form, recl, iostat)
        Implicit None
        
        Integer, Intent(In)           :: unit
        Character(Len=*), Intent(In)  :: file, form
        Integer, Intent(In)           :: recl
        Integer, Intent(InOut)        :: iostat
        
        ConfFileOpenFR = 0
        Open(unit=unit, file=file, form=form, recl=recl, iostat=iostat)
        If (iostat /= 0) ConfFileOpenFR = -1
        Return
    End Function
    
    Integer Function ConfFileOpenSA(unit, file, status, access, iostat)
        Implicit None
        
        Integer, Intent(In)           :: unit
        Character(Len=*), Intent(In)  :: file, status, access
        Integer, Intent(InOut)        :: iostat
        
        ConfFileOpenSA = 0
        Open(unit=unit, file=file, status=status, access=access, iostat=iostat)
        If (iostat /= 0) ConfFileOpenSA = -1
        Return
    End Function
    
    Integer Function ConfFileOpenSR(unit, file, status, recl, iostat)
        Implicit None
        
        Integer, Intent(In)           :: unit
        Character(Len=*), Intent(In)  :: file, status
        Integer, Intent(In)           :: recl
        Integer, Intent(InOut)        :: iostat
        
        ConfFileOpenSR = 0
        Open(unit=unit, file=file, status=status, recl=recl, iostat=iostat)
        If (iostat /= 0) ConfFileOpenSR = -1
        Return
    End Function
    
    Integer Function ConfFileOpenAR(unit, file, access, recl, iostat)
        Implicit None
        
        Integer, Intent(In)           :: unit
        Character(Len=*), Intent(In)  :: file, access
        Integer, Intent(In)           :: recl
        Integer, Intent(InOut)        :: iostat
        
        ConfFileOpenAR = 0
        Open(unit=unit, file=file, access=access, recl=recl, iostat=iostat)
        If (iostat /= 0) ConfFileOpenAR = -1
        Return
    End Function
    
    Integer Function ConfFileOpenF(unit, file, form, iostat)
        Implicit None
        
        Integer, Intent(In)           :: unit
        Character(Len=*), Intent(In)  :: file, form
        Integer, Intent(InOut)        :: iostat
        
        ConfFileOpenF = 0
        Open(unit=unit, file=file, form=form, iostat=iostat)
        If (iostat /= 0) ConfFileOpenF = -1
        Return
    End Function
    
    Integer Function ConfFileOpenS(unit, file, status, iostat)
        Implicit None
        
        Integer, Intent(In)           :: unit
        Character(Len=*), Intent(In)  :: file, status
        Integer, Intent(InOut)        :: iostat
        
        ConfFileOpenS = 0
        Open(unit=unit, file=file, status=status, iostat=iostat)
        If (iostat /= 0) ConfFileOpenS = -1
        Return
    End Function
    
    Integer Function ConfFileOpenA(unit, file, access, iostat)
        Implicit None
        
        Integer, Intent(In)           :: unit
        Character(Len=*), Intent(In)  :: file, access
        Integer, Intent(InOut)        :: iostat
        
        ConfFileOpenA = 0
        Open(unit=unit, file=file, access=access, iostat=iostat)
        If (iostat /= 0) ConfFileOpenA = -1
        Return
    End Function
    
    Integer Function ConfFileOpenR(unit, file, recl, iostat)
        Implicit None
        
        Integer, Intent(In)           :: unit
        Character(Len=*), Intent(In)  :: file
        Integer, Intent(In)           :: recl
        Integer, Intent(InOut)        :: iostat
        
        ConfFileOpenR = 0
        Open(unit=unit, file=file, recl=recl, iostat=iostat)
        If (iostat /= 0) ConfFileOpenR = -1
        Return
    End Function
    
End Module

#ifdef CONFFILEPATHS_UNIT_TEST

Program ConfFilePathsUnitTest

    Use conffilepaths
    
    Implicit None

    Call ConfFileInit(shouldSummarize=.true.)

End

#endif
