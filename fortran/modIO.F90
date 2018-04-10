module IO

  use modPar
  
  implicit none

  integer(kind=iKind) :: logginglevel = 1
  !
  ! logginglevel = 0 ! No messages printed
  !              = 1 ! Routine messages (startup, first entry, computational warnings)
  !              = 2 ! Per outer timestep messages
  !              = 3 ! Inner timestep messages
  !              = 4 ! full debug
  !
  ! Note: messages with level -1 will ALWAYS get printed
  
  !
  ! CONSTANTS
  !
  integer(kind=4), parameter,private    :: FileNameLen = 256
  integer(kind=4), parameter,private    :: errMsgLen   = 256  
  integer(kind=4), parameter,private    :: lKind       = 4
  !
  type tFile
     !
     character(len=FileNameLen)  :: FileName = '' 
     integer(kind=4)             :: FileUnit = 0
     integer(kind=4)             :: IOStat   = 0
     character(len=3)            :: FileType = 'TXT'   !File type [BIN/TXT]
     character(len=10)           :: action   = 'READ'
     character(len=10)           :: access   = 'SEQUENTIAL'
     logical                     :: eof      = .false.
     integer                     :: recL     = 0
     integer(kind=8)             :: CurPos   = 0
     integer(kind=8)             :: NtPos    = 1
     logical                     :: fileopen = .false.
     !
  end type tFile
  character(len=errMsgLen)   :: errMsg
  character(len=FileNameLen) :: scratchFile = 'scratch.ray'
  type(tFile) ::  logfile
  type(tFile) ::  debugfile
  
  logical     ::  log_open = .false.

  type(tFile) ::  amp_file
  type(tFile) ::  drift_file
  type(tFile) :: outputFile
  !
  
  !
  type tOutput
     !
     type(tFile)                  :: file
     logical(lKind)               :: outputrequested = .false.
     integer                      :: iPar(4) = 0
     real(kind=rKind)             :: rPar(4) = 0.
     logical(kind=lKind)          :: lPar(4) = .false.
     real(kind=rKind),allocatable :: data(:)
     !
  end type tOutput
  type(tOutput) :: OutputRay
  !
  
  !
contains
  !
  ! 
  !*****************************************************************************
  !
  logical(kind=lKind) function openFile(File)
  ! 
  !*****************************************************************************

  !=============================================================================
  !                               DESCRIPTION
    !=============================================================================

    ! 1.) PURPOSE:
    !
    !
    ! 2.) METHOD:
    !
    !
    ! 3.) AUTHOR:
    !
    !   P.B. Smit      june 2010      New function
    !

    !
    !=============================================================================
    !                                INTERFACE
    !=============================================================================

    !
    !-----------------------------------------------------------------------------
    !                              DEPENDENCIES
    !-----------------------------------------------------------------------------
    !
    implicit none
    !
    !-----------------------------------------------------------------------------
    !                               VARIABLES
    !-----------------------------------------------------------------------------
    !


    !-----------------------------  ARGUMENT(s) ----------------------------------

    type(TFile), intent(inout)  :: File

    !------------------------------  LOCAL   -------------------------------------

    logical(kind=lKind)         :: exists
    integer(kind=iKind)         :: ioStat

    !=============================================================================
    !                                   SOURCE
    !=============================================================================

    exists = .false._lKind
    File%FileUnit = 0_iKind




    !
    ! Get a free file unit, if some error occurs, raise error
    !
    File%FileUnit = getFreeUnit()
    if (File%FileUnit < 0)     goto 200

    !
    inquire(file=trim(File%FileName),EXIST=exists)
    if ( exists ) then
       !
       if (trim(File%action) == 'WRITE') then
          !
          open(unit=File%FileUnit,file=trim(File%Filename),ACTION='READ')
          close(unit=File%FileUnit,status='delete',iostat=iostat)
          !
       endif
       !
    else         
       !
       If (Trim(File%Action) == 'Read') Then
          goto 100
       endif
       !
    endif
    !

    if     (File%FileType == 'TXT') then
       !
       open(unit=File%FileUnit  ,file=trim(File%Filename)  ,form='FORMATTED'     ,          &
            access='sequential' ,action=trim(File%action)  ,iostat=iostat,Err=300   )
       !
    elseif (File%FileType == 'BIN') then
       !
       open(unit=File%FileUnit  ,file=trim(File%Filename)  ,form='UNFORMATTED'   ,          &
            access=trim(File%access),action=trim(File%action),Err=300   )
    endif
    !    
    OpenFile = .true.
    return
    !

    !=============================================================================
    !                              ERRORS/FORMAT
    !=============================================================================

    !--------------------------------  ERROR(s) ----------------------------------
    !   File does not exist
    call wlog ('READTEXTVAR (io.F90) || Reading failed:' // trim(file%Filename) // ' File does not exist',-1)
100 OpenFile = .false.
    return
    !
    !   Could not find a free unit number
    call wlog ('READTEXTVAR (io.F90) || Reading failed:' // trim(file%Filename) // ' Could not find a free unit number',-1)
200 OpenFile = .false.
    return
    !
    !   Error during opening of file    
    write(errMsg,'(i5)') file%iostat
300 OpenFile = .false.
    call wlog ('READTEXTVAR (io.F90) || Reading failed:' // trim(file%Filename) // 'code: '  &
          // trim(errMsg) // ' Could not open file',-1)
    !
    return
    !

    !--------------------------------  FORMAT(s)----------------------------------

  end function OpenFile

  ! 
  !*****************************************************************************
  !
  integer(kind=iKind) function getFreeUnit()
    ! 
    !*****************************************************************************
    !

    !-----------------------------------------------------------------------------
    !                                   INTERFACE
    !-----------------------------------------------------------------------------

    !   --------- DESCRIPTION  --------

    !   Function returning an available unit number for i/o operations

    !   --------- DEPENDENCIES --------

    implicit none

    !   --------- VARIABLES --------
    !
    !   ARGUMENT
    !
    !    
    !   LOCAL
    !
    integer(kind=iKind),parameter :: startUnit   = 10
    integer(kind=iKind),parameter :: maxUnit     = 10000
    integer(kind=iKind)           :: iUnit
    logical(kind=lKind)           :: connected
    !  
    !-----------------------------------------------------------------------------
    !                                    SOURCE
    !-----------------------------------------------------------------------------
    !
    !
    connected     = .true.
    getFreeUnit   = -1
    !
    FindUnitLoop: do iUnit=startUnit,maxUnit
       !
       inquire(unit=iUnit,opened=connected)
       if (.not. connected) then           
          getFreeUnit = iUnit
          exit FindUnitLoop
       endif
       !
    enddo FindUnitLoop
    !
  end function getFreeUnit
  !

 subroutine log_int( msg , ival, ilevel )
    !
    character( len = * ),intent(in) :: msg
    integer(kind=iKind)    ,intent(in)   :: ival
    integer(kind = iKind),intent(in)  :: ilevel

    character( len = 100)             :: msg2

    write( msg2 , '(i10)' ) ival
    msg2 = trim( msg ) // trim( msg2 )
    
    call wlog( msg2 , ilevel )
    !
  end subroutine log_int
  
  subroutine log_real( msg , rval, ilevel )
    !
    character( len = * ),intent(in) :: msg
    real(kind=rKind)    ,intent(in)   :: rval
    integer(kind = iKind),intent(in)  :: ilevel

    character( len = 100)             :: msg2

    write( msg2 , '(E15.7)' ) rval
    msg2 = trim( msg ) // trim( msg2 )
    
    call wlog( msg2 , ilevel )
    !
  end subroutine log_real
  

  
  subroutine wlog( message ,ilevel )
    !
    character( len = * ) ,intent(in) :: message
    integer(kind = iKind),intent(in) :: ilevel
    !
    if ( ilevel <= loggingLevel ) then
       !
       if ( init_log() ) then
          !
          write( logfile%fileunit , * ) trim(message)
          flush( logfile%fileunit )
          !
       endif
       !
    endif    
    !
  end subroutine wlog

  
  logical function init_log()
    !
    if ( log_open ) then
       !
       init_log = .true.
       return
       !
    endif

    logfile%FileName = 'logfile.txt'
    logfile%FileType = 'TXT'
    logfile%action   = 'WRITE'

    if ( openFile(logfile) ) then
       !
       log_open = .true.
       init_log = .true.
       !
    else
       !
       log_open = .false.
       init_log = .false.
       !
    endif
    write( logfile%fileunit , '("%--------------------------------------------------------------------------------------")' )
    write( logfile%fileunit , '("%                             LOG FILE DRIFTERS LIBRARY                                ")' )
    write( logfile%fileunit , '("%--------------------------------------------------------------------------------------")' )
    flush( logfile%fileunit )
    !
  end function init_log
  !


  
  !  
end module IO
