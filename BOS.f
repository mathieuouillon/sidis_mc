C***********************************************************************
      SUBROUTINE CLASBOSFILL()
C***********************************************************************
      include "bcs.inc"
      include "names.inc"
      include "common.f"
c
c
      integer ierr,lunout
      parameter (lunout=33)
      INTEGER ind,indx,mctk,mcvx,mbank,part,boswrite
      INTEGER j,jj
      real charge
C ********************************************************
C       Write to the output event file
C ********************************************************
      call bdrop(iw,'E')
      call bgarb(iw)
      ind=mbank(iw,'HEAD',0,8,1)
      if(ind.ne.0)then
c        write(*,*) 'test HEAD bank'
        iw(ind+1)=1
        iw(ind+2)=1
        iw(ind+3)=IEVENT
        iw(ind+4)= 100           ! event time
        iw(ind+5)=-4
        iw(ind+6)=0
        iw(ind+7)=15
        iw(ind+8)=0
      endif
      ind= nbank('TGBI',0,4,1)
      if(ind.ne.0)then
        if (cl_pol.gt.0.5) then
         iw(ind+1) = 32768
        else
         iw(ind+1) = 0 
        endif
      endif
c
c  get mcvx,mctk bank pointers
ccc        mcvx=nbank('MCVX',0,5,1)
ccc        mctk=nbank('MCTK',0,11,4)

       jj=0
         do 998 j=1,N       
       if(k(j,1).gt.10.and..not.(k(j,2).eq.111)) goto 998
       if(k(j,2).eq.22.and.  k(k(j,3),2).eq.111) goto 998
       jj=jj+1
 998   continue
       ind=mbank(iw,'MCTK',0,11,jj) 
       mcvx=mbank(iw,'MCVX',0,5,jj)
       part= nbank('PART',0,13,jj) 
        jj=0
       if(ind.ne.0)then
         do 999 j=1,N
       if(k(j,1).gt.10.and..not.(k(j,2).eq.111)) goto 999
       if(k(j,2).eq.22.and.  k(k(j,3),2).eq.111) goto 999
       charge = 0
       select case (k(j,2))
C ... Electron, Gamma, positon
           case (11)
             charge = -1.
           case (22)
             charge = 0.
           case (-11)
             charge = 1.
C ... Pi0 
           case (111)
             charge = 0.
C ... Pi+ 
           case (211)
             charge = 1.
C ... Pi- 
           case (-211)
             charge = -1.
C ... K+   
           case (321)
             charge = 1.
C ... K-   
           case (-321)
             charge = -1.
C ... p    
           case (2212)
             charge = 1.
C ... p bar
           case (-2212)
             charge = -1.
C ... n    
           case (2112)
             charge = 0.
C ... n bar
           case (-2112)
             charge = 0.
C ... K long
           case (130)
             charge = 0.
C ... mu +
           case (-13)
             charge = 1.
C ... mu -
           case (13)
             charge = -1.
C ... neutrino e
           case (12)
             charge = 0.
C ... neutrino e bar
           case (-12)
             charge = 0.
           case default
             write(*,*) 'hey une part non identifie!'
             write(*,*) k(j,1),k(j,2)

         end select
         jj=jj+1
           indx=ind+(jj-1)*11
          pmom=p(j,1)*p(j,1)+p(j,2)*p(j,2)+p(j,3)*p(j,3)
          if(pmom.gt.0.0)then
             pmom=sqrt(pmom)
             rw(indx+1)=p(j,1)/pmom             !cx /plu(j,8) 
             rw(indx+2)=p(j,2)/pmom             !cy /plu(j,8) 
             rw(indx+3)=p(j,3)/pmom             !cz /plu(j,8)
           else
             rw(indx+1)=0.0 
             rw(indx+2)=0.0 
             rw(indx+3)=0.0
             rw(indx+4)=pmom                     !plu(j,8)
            endif 
           rw(indx+4)=pmom                     !plu(j,8)
c           rw(indx+5)=p(j,5)                ! mass of the particle
c           rw(indx+6)=plu(j,6)              ! charge
           rw(indx+6)= charge              ! charge
           iw(indx+7)=k(j,2)    ! PID LUND code 
           iw(indx+8)=k(j,1)    ! LUND Status
           iw(indx+9)=1         ! Beg. Vertex. 
           if(k(j,1).gt.10)then
             iw(indx+10)=1      ! End Vertex.
           else
             iw(indx+10)=0         
           endif
           iw(indx+11)=k(k(j,3),2)   ! LUND code for parent track
       if(part.ne.0) then
         iw(part+1)     = lund2geantid(k(j,2))                          ! particle ID (GEANT)
         rw(part+2)     = v(j,1)                                        ! x vertex position
         rw(part+3)     = v(j,2)                                        ! y vertex position
         rw(part+4)     = v(j,3)-25                                    ! z vertex position
         rw(part+5)     = sqrt(p(j,5)*p(j,5)+pmom*pmom)                 ! energy 
         rw(part+6)     = p(j,1)                                        ! px
         rw(part+7)     = p(j,2)                                        ! py
         rw(part+8)     = p(j,3)                                        ! pz
c         rw(part+9)     = plu(j,6)                                      ! charge
         rw(part+9)     = charge                                      ! charge
         iw(part+10)    = j                                             ! Track pointer
         rw(part+11)    = 0.
         rw(part+12)    = 0.
         iw(part+13)    = k(k(j,3),2)
         part=part+13
        endif

       if(mcvx.ne.0)then
           rw(mcvx+1)=v(j,1) 
           rw(mcvx+2)=v(j,2) 
           rw(mcvx+3)=v(j,3)-25
c           rw(mcvx+3)=v(j,3)-25+(2.*ranf(0)-1)*.001
           rw(mcvx+4)=0.
           iw(mcvx+5)=0
           mcvx=mcvx+5
       endif

 999    CONTINUE
       endif
c 

c
      call fwbos(iw,33,'E',ierr)        ! write banks to file
      call bdrop(iw,'E')                ! drop the bank, reclaim the space
      call bgarb(iw)                    ! garbage collection
c        ierr=bosWrite(lunout,iw,'HEADMCTKMCVXTAGR')
c
c$$$      targs	= t_targ*ran(idum)
c$$$
c$$$c     Change into proper coordinate system
c$$$      v(j,3)	= v(j,3) + targs - t_targ / 2.0
c$$$               z_pos=z_centr+2.*random_num()*z_width-z_width

      RETURN
      END     


      INTEGER FUNCTION LUND2GEANTID(I)
      PARAMETER (NSEL=44)
      INTEGER IPLUND(NSEL),IGE
      DATA IPLUND/
     +    22,   -11,    11,    12,   -13,    13,   111,   211,
     +  -211,   130,   321,  -321,  2112,  2212, -2212,   310,
     +   221,  3122,  3222,  3212,  3112,  3322,  3312,  3334,
     + -2112, -3122, -3112, -3212, -3222, -3322, -3312, -3334,
     +   -15,    15,   411,  -411,   421,  -421,   431,  -431,
     +  4122,    24,   -24,    23/

       DO IGE=1,NSEL
        IF(I.EQ.IPLUND(IGE)) THEN
         LUND2GEANTID=IGE
         RETURN
        ENDIF
       ENDDO
         LUND2GEANTID=0 
       RETURN
       END
c
      INTEGER FUNCTION GEANT2LUNDID(I)
      PARAMETER (NSEL=44)
      INTEGER IPLUND(NSEL),I
      DATA IPLUND/
     +    22,   -11,    11,    12,   -13,    13,   111,   211,
     +  -211,   130,   321,  -321,  2112,  2212, -2212,   310,
     +   221,  3122,  3222,  3212,  3112,  3322,  3312,  3334,
     + -2112, -3122, -3112, -3212, -3222, -3322, -3312, -3334,
     +   -15,    15,   411,  -411,   421,  -421,   431,  -431,
     +  4122,    24,   -24,    23/

      IF(I.GT.0.AND.I.LE.NSEL) THEN
       GEANT2LUNDID=IPLUND(I)
      ELSE
       GEANT2LUNDID=0
      ENDIF 
       RETURN
      END

C***********************************************************************
      SUBROUTINE CLASBOSEND(recname)
C***********************************************************************
      implicit none
      include "bcs.inc"
      character*8 recname
      integer ierr
      call fwbos(iw,33,'0',ierr)
      call bosta
      call fparm('CLOSE ' //recname)
c      call fclos()

      print *,'END CLASBOSEND',ierr

      RETURN
      END     


C***********************************************************************
      SUBROUTINE CLASBOSINIT(recname)
C***********************************************************************
      implicit none
      
      character*8 recname
      character*132 ddl_file
      include "bcs.inc"
      include "names.inc"
c
c       Open BOS input & output files
c
      CALL bnames(1000)
      call bos(iw,nbcs)
c
      call revinm ('CLAS_PARMS','clasbanks.ddl',ddl_file)
      CALL txt2bos(ddl_file) !'clasbanks.ddl')

      close (33)

      call fparm ('OPEN UNIT=33 FILE="'//bosout//'" 
     &      RECL=32760 ACTION=WRITE STATUS=NEW FORM=BINARY')

      RETURN
      END     




      SUBROUTINE TXT2BOS(FILENAME)
c
c_begin_doc
c  RCS ID string
c  $Id: clasDIS.F,v 1.3 2001/05/16 00:04:42 avakian Exp $
c
c  Documentation for subroutine TXT2BOS
c
c  Purpose: reads in an ascii file (FILENAME) and converts that info into
c  -------- BOS bank definition.
c
c  Input Parameters:  FILENAME :C*(*): ascii file contain BOS bank definitions
c  ----------------   Format is as follows:
c anything following a `*' or a `!' is considered a comment
c   *   this is a comment, the following is an example bank definition
c   *   leading blanks are ignored
C   * the first line should contain a 4 character  NAME followed by
c   * the name of the bank and ACTION word (CREATE, WRITE, DELETE, MODIFY)
c   *  TABLE BANKname fmt  ACTION1 ACTION2
c   * The NAME line is followed by a row ordered list of format statements
c   * as defined in page 10 of the BOS manual, for example
c   *           * run and event number 
c   *           * 2 group of 3 floating points
c   * 
c   * ACTION = CREATE implies MKFMT will be called
c   * ACTION = WRITE implies BLIST(IW,'E+','bname') will be called,
c   *          absence of WRITE implies BLIST(IW,'E-','bname') will be called
c   * ACTION = DELETE implies BLIST(IW,'R+','bname') will be called
c   * ACTION = MODIFY implies THIS definition should OVERRIDE definition in
c   *                 datafile
c   * ACTION = DISPLAY implies the bank will be `shipped' to the event store
c              buffer so that CED can display it
c
c  Output Parameters:
c  -----------------
c
c  Other routines:
c  ---------------
c
c  Notes:
c  ------
c
c  Author:   Arne Freyberger      Created:  Fri Oct  6 15:21:23 EDT 1995
c  -------
c
c  Major revisions:
c  ----------------
c     
c
c_end_doc
c
      IMPLICIT NONE
      SAVE
c
c_begin_inc
c  include files :
c  ---------------------
c BOS common block  uncomment the next line for BOS include file
      include "bcs.inc"
      include "bnkfmt.inc"
c                           CLAS control module
      include "clasmdl.inc"
c_end_inc
c
c_begin_var
c  input/output variables:
c  -----------------------
c
c  Local pre-defined variables:
c  ---------------------------
c  RCS information: 
      CHARACTER*132  CFILE, CREVIS, CSTATE, CDATE, CAUTHO, CRCSID
      PARAMETER (CFILE = '$rcsRCSfile$')
      PARAMETER (CREVIS = '$rcsRevision$')
      PARAMETER (CSTATE = '$rcsState$')
      PARAMETER (CDATE = '$rcsDate$')
      PARAMETER (CAUTHO = '$rcsAuthor$')
      DATA CRCSID/   
     1'$Id: clasDIS.F,v 1.3 2001/05/16 00:04:42 avakian Exp $'   
     2/   
c  Module information:
      CHARACTER*(*)  CRNAME, CRAUTH
      CHARACTER*100  CRMESS
      PARAMETER (CRNAME='TXT2BOS')
      PARAMETER (CRAUTH='Arne Freyberger')
c
c  Local User defined variables:
c  -----------------------------
      INTEGER   INDEXA, INDEXN
      EXTERNAL  INDEXA, INDEXN
      INTEGER ICLOC,  ISCAN, lenocc
      EXTERNAL ICLOC, ISCAN, lenocc
      character*132 spaces, strip
      external spaces, strip
c
      CHARACTER*(*) FILENAME
      character*512 LINE
      CHARACTER*256 CNAME, CWORD, CTMP
      CHARACTER*400 BFORMAT
      CHARACTER*4 BNAME
      CHARACTER*20 CWORDS(40)
      character*128  ctest
      LOGICAL LCREATE, LMODIFY, LWRITE, LDELETE
      LOGICAL LDONE, LCLEAR, LDISPLAY, LINIT, LFORMAT_DONE
      INTEGER ILUN, IRET, NLINE
      INTEGER IN_BEG, IN_END, IFORMAT_BEG, IFORMAT_END, IBEG, IEND
      INTEGER IBEGA, IBEGN, IEND1, IEND2
      INTEGER WCNT, IWORDS, ICOLUMN, iretlen
c_end_var
      DATA LINIT/.FALSE./
      DATA BFORMAT/'                              '/
c
c  executable code for routine TXT2BOS:
c----6----------------------------------------------------------------72
c
      IF (.NOT. LINIT) THEN
        LINIT = .TRUE.
        IDELETE = 0
        IWRITE = 0
        IREMOVE = 0
        ICREATE = 0
      ENDIF

c open the ascii file
 
      CALL RESLUN(CRNAME,ILUN,IRET)
      CALL REOPEN(FILENAME,ILUN,'OLD',IRET)
      if (iret .lt. 0) return
c 
      ICOLUMN = 0

c      write(crmess,1002)
 1002 format('Bankname : Format',T45,'WRITE DELETE DISPLAY')
c      call recmes(crname,'i',crmess)

c  do stuff with the file, read in line by line and 
c  parse each line into words

      LCLEAR = .TRUE.
10    CONTINUE
      READ(ILUN,11,ERR=999,END=999)LINE
11    FORMAT(A512)
      NLINE = NLINE + 1
      LINE = SPACES(LINE,1)
      IF (LINE(1:1) .EQ. '*' .OR. 
     1    LINE(1:1) .EQ. '!') GOTO 10

c  convert to upper case

      CALL CLTOU(LINE(1:256))

c construct the BOS bank, first break the line into words
c separated by blanks

      IEND1 = ICLOC('*',1,LINE,1,256)
      IEND2 = ICLOC('!',1,LINE,1,256)
      IEND = 256
      IF (IEND1*IEND2 .NE. 0) THEN
       IEND = MIN(IEND1,IEND2)
      ELSE
       IEND = MAX(IEND1,IEND2)
      ENDIF
      IF (IEND .EQ. 0) GOTO 10
      IBEGA = INDEXA(LINE(1:IEND))
      IBEGN = INDEXN(LINE(1:IEND))
      IF (IBEGN*IBEGA .EQ. 0) THEN
        IBEG = MAX(IBEGA,IBEGN)
      ELSE
        IBEG = MIN(IBEGA,IBEGN)
      ENDIF
      if (ibeg .eq. 0) goto 10
      IF (IEND .EQ. 0) IEND = 256
c
c  separate the line into words then loop over them
c
      CALL REWORD(LINE,IWORDS,CWORDS)
      IF (IWORDS .LE. 1) GOTO 10
      CWORD = CWORDS(1)
      IF (CWORD(1:1) .EQ. '!' .OR. CWORD(1:1) .EQ. '*') GOTO 10
c
c  look for TABLE as the START of the bank definition
c
      IF (CWORD(1:2) .EQ. 'TA') THEN 
        IF (.NOT. LCLEAR .and. IFORMAT_END .gt. 1) THEN
           IN_END = IFORMAT_END - 1
           iretlen = 1
           call fmcre(bformat(1:IN_END+2), ctest, iretlen)
c           WRITE(CRMESS,1001)
c     1     BNAME,BFORMAT(1:min(23,IFORMAT_END-1)), 
c     1     BNAME,ctest(1:30), 
c     2                LWRITE, LDELETE, LDISPLAY
c           CALL RECMES(' ','N',CRMESS)
c           IF (LCREATE) THEN
               if (iretlen.LT.6) then
                  CALL MKFMT(IW,BNAME,BFORMAT(1:IN_END))
               else
                  CALL MKFMT(IW,BNAME,ctest(1:iretlen))
               endif
             ICREATE = ICREATE + 1
             CREATE_NAME(ICREATE) = BNAME
             CREATE_FORMAT(ICREATE) = BFORMAT(1:IN_END)
             NCREATE_COL(ICREATE) = ICOLUMN
c           ENDIF
           IF (LWRITE) THEN
              CALL BLIST(IW,'E+',BNAME)
              IWRITE = IWRITE + 1
              CBANK_WRITE(1+4*(IWRITE-1):4+4*(IWRITE-1)) = BNAME
           ELSE
              CALL BLIST(IW,'E-',BNAME)
              IREMOVE = IREMOVE + 1
              CBANK_REMOVE(1+4*(IREMOVE-1):4+4*(IREMOVE-1)) = BNAME
           ENDIF
           IF (LDELETE) THEN 
              CALL BLIST(IW,'R+',BNAME)
              IDELETE = IDELETE + 1
              CBANK_DELETE(1+4*(IDELETE-1):4+4*(IDELETE-1)) = BNAME
           ENDIF
           IF (LDISPLAY) THEN
            NBANK_DISPLAY = NBANK_DISPLAY + 1
            CBANK_DISPLAY(NBANK_DISPLAY) = BNAME
            CBANK_FORMAT(NBANK_DISPLAY) = BFORMAT(1:IN_END)
           ENDIF
           LDONE = .FALSE.
           LFORMAT_DONE = .FALSE.
           IFORMAT_END = 0
           ICOLUMN = 0
           BNAME = '    '
           BFORMAT = '                '
        ENDIF

c  Done writting out previous bank now start over and initialize logicals

        LCLEAR = .FALSE.
        LCREATE = .FALSE.
        LDELETE = .FALSE.
        LWRITE = .FALSE.
        LMODIFY = .FALSE.
        LDISPLAY = .FALSE.
        LFORMAT_DONE = .FALSE.

c BANK name is the SECOND word!!

        BNAME = CWORDS(2)
c
        IF (BNAME .NE. '    ') THEN
c
        DO 23 WCNT=3,IWORDS
         CNAME = cwords(wcnt)
         IF (CNAME(1:3) .EQ. 'B16' .OR.
     1       CNAME(1:3) .EQ. 'B08' .OR.
     2       CNAME(1:3) .EQ. 'B32') THEN 
          BFORMAT = CNAME
          IFORMAT_END = LENOCC(BFORMAT) + 1
          LFORMAT_DONE = .TRUE.
         ELSEIF (CNAME(1:6) .EQ. 'CREATE') THEN
          LCREATE = .TRUE.
         ELSEIF (CNAME(1:5) .EQ. 'WRITE') THEN
          LWRITE = .TRUE.
         ELSEIF (CNAME(1:6) .EQ. 'DELETE') THEN
          LDELETE = .TRUE.
         ELSEIF (CNAME(1:6) .EQ. 'MODIFY') THEN
          LMODIFY = .TRUE.
         ELSEIF (CNAME(1:7) .EQ. 'DISPLAY') THEN
          LDISPLAY = .TRUE.
         ENDIF
23      CONTINUE
        ENDIF

c   check for END or END TABLE presence

      ELSEIF (CWORD(1:3) .EQ. 'END' .AND. BNAME .NE. '    ') THEN
          LDONE = .TRUE.
c
c  the following lines parse the FORMAT lines (not a TAble or END line)
c
      ELSEIF (BNAME .NE. '    ') THEN
         IFORMAT_BEG = IFORMAT_END + 1
         CNAME = CWORDS(3)
         IN_END = ISCAN(CNAME,' ') - 1
         if (in_end .eq. 0) goto 10
c
c  reject relational formats
c
         IF (CNAME(1:1) .EQ. 'D' .OR.
     1       CNAME(1:1) .EQ. 'R' .OR.
     2       CNAME(1:1) .EQ. 'M') GOTO 10
         ICOLUMN = ICOLUMN + 1
         ctmp = cwords(2)
         create_elements(icreate+1,icolumn)= ctmp(1:8)
         if (.not. lformat_done) then
          IN_BEG = 1
          IFORMAT_END = IFORMAT_BEG + IN_END - 1
          BFORMAT(IFORMAT_BEG:IFORMAT_END) = CNAME(1:IN_END)
          IFORMAT_END = IFORMAT_END + 1
          BFORMAT(IFORMAT_END:IFORMAT_END) = ','
         endif
       ENDIF
c
       IF (LDONE .AND. BNAME .NE. '    ') THEN
           IN_END = IFORMAT_END - 1
           CALL MKFMT(IW,BNAME,BFORMAT(1:IN_END))
           ICREATE = ICREATE + 1
           CREATE_NAME(ICREATE) = BNAME
           CREATE_FORMAT(ICREATE) = BFORMAT(1:IN_END)
           NCREATE_COL(ICREATE) = ICOLUMN

           IF (LWRITE) THEN
              CALL BLIST(IW,'E+',BNAME)
              IWRITE = IWRITE + 1
              CBANK_WRITE(1+4*(IWRITE-1):4+4*(IWRITE-1)) = BNAME
           ELSE
              CALL BLIST(IW,'E-',BNAME)
              IREMOVE = IREMOVE + 1
              CBANK_REMOVE(1+4*(IREMOVE-1):4+4*(IREMOVE-1)) = BNAME
           ENDIF
           IF (LDELETE) THEN 
              CALL BLIST(IW,'R+',BNAME)
              IDELETE = IDELETE + 1
              CBANK_DELETE(1+4*(IDELETE-1):4+4*(IDELETE-1)) = BNAME
           ENDIF
           IF (LDISPLAY) THEN
            NBANK_DISPLAY = NBANK_DISPLAY + 1
            CBANK_DISPLAY(NBANK_DISPLAY) = BNAME
            CBANK_FORMAT(NBANK_DISPLAY) = BFORMAT(1:IN_END)
           ENDIF
           iretlen = 1
           call fmcre(bformat, ctest, iretlen)
c           WRITE(CRMESS,1001)
c     1     BNAME,BFORMAT(1:min(23,IFORMAT_END-1)), 
c     1     BNAME,ctest, 
c     2                LWRITE, LDELETE, LDISPLAY
 1001      FORMAT('|',13x,A4,' : ',A30,' : ',T60,3(L1,6X))
c           CALL RECMES(' ','N',CRMESS)
           LDONE = .FALSE.
           LFORMAT_DONE = .FALSE.
           LCLEAR = .TRUE.
           IFORMAT_END = 0
           ICOLUMN = 0
           BNAME = '    '
           BFORMAT = '            '
       ENDIF
      GOTO 10
c
c  but format last bank before exiting
c
999   CONTINUE
      if (.not. lclear) then
       iretlen = 1
       call fmcre(bformat, ctest, iretlen)
c       WRITE(CRMESS,1001)
c     1     BNAME,BFORMAT(1:min(23,IFORMAT_END-1)), 
c     1     BNAME,ctest, 
c     2                LWRITE, LDELETE, LDISPLAY
c       CALL RECMES(' ','N',CRMESS)
       if (iretlen.LT.6) then
          CALL MKFMT(IW,BNAME,BFORMAT(1:IFORMAT_END-1))
       else
          CALL MKFMT(IW,BNAME,ctest(1:iretlen))
       endif
       ICREATE = ICREATE + 1
       CREATE_NAME(ICREATE) = BNAME
       CREATE_FORMAT(ICREATE) = BFORMAT(1:IFORMAT_END-1)
       NCREATE_COL(ICREATE) = ICOLUMN

       IF (LWRITE) THEN
           CALL BLIST(IW,'E+',BNAME)
           IWRITE = IWRITE + 1
           CBANK_WRITE(1+4*(IWRITE-1):4+4*(IWRITE-1)) = BNAME
       ELSE
           CALL BLIST(IW,'E-',BNAME)
           IREMOVE = IREMOVE + 1
           CBANK_REMOVE(1+4*(IREMOVE-1):4+4*(IREMOVE-1)) = BNAME
       ENDIF
       IF (LDELETE) THEN 
           CALL BLIST(IW,'R+',BNAME)
           IDELETE = IDELETE + 1
           CBANK_DELETE(1+4*(IDELETE-1):4+4*(IDELETE-1)) = BNAME
       ENDIF
       IF (LDISPLAY) THEN
         NBANK_DISPLAY = NBANK_DISPLAY + 1
         CBANK_DISPLAY(NBANK_DISPLAY) = BNAME
         CBANK_FORMAT(NBANK_DISPLAY) = BFORMAT(1:IN_END)
       ENDIF
c       iretlen = 30
c       call fmcre(bformat, ctest, iretlen)
c       WRITE(CRMESS,1001)
c     1     BNAME,BFORMAT(1:min(23,IFORMAT_END-1)), 
c     1     BNAME,ctest, 
c     2                LWRITE, LDELETE, LDISPLAY
c       CALL RECMES(' ','I',CRMESS)
       IFORMAT_END = 0
      endif
c      write(crmess,1002)
c      call recmes(crname,'i',crmess)
      CALL RESLUN(CRNAME,-ILUN,IRET)
c
c$$$      call recmes(crname,'i',
c$$$     1  'The following banks will be DELETED from memory before reading in the next event')
c$$$      if (idelete .ne. 0) 
c$$$     1      call recmes(crname,'i',cbank_delete(1:4+4*(IDELETE-1)))
c$$$      call recmes(crname,'i',
c$$$     1      'The following banks will NOT be written out')
c$$$      if (iremove .ne. 0) 
c$$$     1      call recmes(crname,'i',cbank_remove(1:4+4*(IREMOVE-1)))
c$$$      call recmes(crname,'i',
c$$$     1      'The following banks WILL be written out')
c$$$      if (iwrite .ne. 0) 
c$$$     1      call recmes(crname,'i',cbank_write(1:4+4*(IWRITE-1)))
      CLOSE(ILUN)
      RETURN
      END
