      subroutine init_fastmc(parfile)
c==   Read parameters from data base
      implicit none
      CHARACTER*30 parfile 
      include "read_par_clas12g.inc"

      PRINT *,'|==========> Cofiguration file: ',parfile
      OPEN(11,FILE=parfile,STATUS='old') ! Control par.file

      CALL readint4(11,'Cnfg_num  ',Cnfg_num) ! Configuration No 

      CALL readrea4(11,'Tmax      ',Tmax) ! Max Torus Current
      CALL readrea4(11,'Thf_lim   ',Thf_lim) ! Theta forward max

c Acceptance:

      CALL readrea4(11,'Thacc1    ',Thacc1)!First  theta acceptance parameter 
      CALL readrea4(11,'Thacc2    ',Thacc2)!Second theta acceptance parameter 
      CALL readrea4(11,'Thacc3    ',Thacc3)!Third  theta acceptance parameter 
      CALL readrea4(11,'Thacc4    ',Thacc4)!Fourth theta acceptance parameter 
      CALL readrea4(11,'Thacc5    ',Thacc5)!Fifth  theta acceptance parameter 
      CALL readrea4(11,'Thacc6    ',Thacc6)!Sixth  theta acceptance parameter 

      CALL readrea4(11,'Thacc01   ',Thacc01)!First  theta acceptance parameter(Theta_acc01) 
      CALL readrea4(11,'Thmom1    ',Thmom1 )!First  theta-momentum acc. parameter
      CALL readrea4(11,'Thmom2    ',Thmom2 )!Second theta-momentum acc. parameter
      CALL readrea4(11,'Thmom3    ',Thmom3 )!Third  theta-momentum acc. parameter
      CALL readrea4(11,'Thmom4    ',Thmom4 )!Fourth theta-momentum acc. parameter
      CALL readrea4(11,'Phimom1   ',Phimom1)!First  Phi-momentum acc. parameter        
      CALL readrea4(11,'Phimom2   ',Phimom2)!Second Phi-momentum acc. parameter        
      CALL readrea4(11,'Phimom3   ',Phimom3)!Third  Phi-momentum acc. parameter        
      CALL readrea4(11,'Phimom4   ',Phimom4)!Fourth Phi-momentum acc. parameter        
      CALL readrea4(11,'Phimom5   ',Phimom5)!Fifth  Phi-momentum acc. parameter 

c    FORWARD:
ckm may 23, 2005. new parameters added
      CALL readrea4(11,'FoMT1S1    ',FoMT1S1)!Forward 1-st param for sigma_1
      CALL readrea4(11,'FoMT1S2    ',FoMT1S2)!Forward 2-d param for sigma_1 
      CALL readrea4(11,'FoMT1S3    ',FoMT1S3)!Forward 3-d param for sigma_1 
      CALL readrea4(11,'FoMT2S1    ',FoMT2S1)!Forward 1-st param for sigma_1
      CALL readrea4(11,'FoMT2S2    ',FoMT2S2)!Forward 2-d param for sigma_1 
      CALL readrea4(11,'FoMT2S3    ',FoMT2S3)!Forward 3-d param for sigma_1 

      CALL readrea4(11,'FoInMoS1   ',FoInMoS1)!Forward inbending mometum sigma_1
      CALL readrea4(11,'FoInMoS2   ',FoInMoS2)!Forward inbending mometum sigma_2
      CALL readrea4(11,'FoOuMoS1   ',FoOuMoS1)!Forward outbending mometum sigma_1
      CALL readrea4(11,'FoOuMoS2   ',FoOuMoS2)!Forward outbending mometum sigma_2
      CALL readrea4(11,'FoNeMoS1   ',FoNeMoS1)!Forward neutral mometum sigma_1
      CALL readrea4(11,'FoNeMoS2   ',FoNeMoS2)!Forward neutral mometum sigma_2
                      
      CALL readrea4(11,'FoInThS1   ',FoInThS1)!Forward inbending  Theta sigma_1
      CALL readrea4(11,'FoInThS2   ',FoInThS2)!Forward inbending  Theta sigma_2
      CALL readrea4(11,'FoOuThS1   ',FoOuThS1)!Forward outbending Theta sigma_1
      CALL readrea4(11,'FoOuThS2   ',FoOuThS2)!Forward outbending Theta sigma_2
      CALL readrea4(11,'FoNeThS1   ',FoNeThS1)!Forward neutral    Theta sigma_1
      CALL readrea4(11,'FoNeThS2   ',FoNeThS2)!Forward neutral    Theta sigma_2
ckm May 28, 2005: Gail's phi-resolution.                   
ccc      CALL readrea4(11,'FoIOFiS1   ',FoIOFiS1)!Forward in(out)bending  phi sigma_1
ckm May 31, 2005: New fit from Gail:

      CALL readrea4(11,'FoA0FiS1   ',FoA0FiS1)!Forward in(out)bending  A0 phi sigma_1
      CALL readrea4(11,'FoA1FiS1   ',FoA1FiS1)!Forward in(out)bending  A1 phi sigma_1
      CALL readrea4(11,'FoA2FiS1   ',FoA2FiS1)!Forward in(out)bending  A2 phi sigma_1

      CALL readrea4(11,'FoA0FiS2   ',FoA0FiS2)!Forward in(out)bending  A0 phi sigma_2
      CALL readrea4(11,'FoA1FiS2   ',FoA1FiS2)!Forward in(out)bending  A1 phi sigma_2
      CALL readrea4(11,'FoA2FiS2   ',FoA2FiS2)!Forward in(out)bending  A1 phi sigma_2
      CALL readrea4(11,'FoNeFiS1   ',FoNeFiS1)!Forward neutral    phi sigma_1
      CALL readrea4(11,'FoNeFiS2   ',FoNeFiS2)!Forward neutral    phi sigma_2

c    CENTRAL:

      CALL readrea4(11,'CeA0ThS2   ',CeA0ThS2)!Central in(out)bending  A0  Theta sigma_2
      CALL readrea4(11,'CeA1ThS2   ',CeA1ThS2)!Central in(out)bending  A1  Theta sigma_2
      CALL readrea4(11,'CeA2ThS2   ',CeA2ThS2)!Central in(out)bending  A2  Theta sigma_2
      CALL readrea4(11,'CeInMoS1   ',CeInMoS1)!Central inbending mometum sigma_1
      CALL readrea4(11,'CeInMoS2   ',CeInMoS2)!Central inbending mometum sigma_2
      CALL readrea4(11,'CeOuMoS1   ',CeOuMoS1)!Central outbending mometum sigma_1
      CALL readrea4(11,'CeOuMoS2   ',CeOuMoS2)!Central outbending mometum sigma_2
      CALL readrea4(11,'CeNeMoS1   ',CeNeMoS1)!Central neutral mometum sigma_1
      CALL readrea4(11,'CeNeMoS2   ',CeNeMoS2)!Central neutral mometum sigma_2
                                                                   
      CALL readrea4(11,'CeInThS1   ',CeInThS1)!Central inbending  Theta sigma_1
      CALL readrea4(11,'CeInThS2   ',CeInThS2)!Central inbending  Theta sigma_2
      CALL readrea4(11,'CeOuThS1   ',CeOuThS1)!Central outbending Theta sigma_1
      CALL readrea4(11,'CeOuThS2   ',CeOuThS2)!Central outbending Theta sigma_2
      CALL readrea4(11,'CeNeThS1   ',CeNeThS1)!Central neutral    Theta sigma_1
      CALL readrea4(11,'CeNeThS2   ',CeNeThS2)!Central neutral    Theta sigma_2
                                                                   
      CALL readrea4(11,'CeInFiS1   ',CeInFiS1)!Central inbending  phi sigma_1
      CALL readrea4(11,'CeInFiS2   ',CeInFiS2)!Central inbending  phi sigma_2
      CALL readrea4(11,'CeOuFiS1   ',CeOuFiS1)!Central outbending phi sigma_1
      CALL readrea4(11,'CeOuFiS2   ',CeOuFiS2)!Central outbending phi sigma_2
      CALL readrea4(11,'CeNeFiS1   ',CeNeFiS1)!Central neutral    phi sigma_1
      CALL readrea4(11,'CeNeFiS2   ',CeNeFiS2)!Central neutral    phi sigma_2


c close data parameters file 
      CLOSE(11)
      return
      end
C================================================================
      SUBROUTINE readrea4(lun,key,rea4)
C---> Read a value (rea4) of variable (key) from unit (lun).
C     ======================================================
      CHARACTER*10 key
      CHARACTER*4  type
      CHARACTER*50 text
      call readtest(1,lun,key,type,ierr)
      IF (type.EQ.'REA4') THEN
        READ(lun,FMT='(18X,E10.4,A50)',ERR=13,END=13) rea4,text
        PRINT 2,   key,rea4,text
        WRITE(9,2) key,rea4,text
  2     FORMAT('      REA4|==> ',a10,'= ',E10.4,2x,a50)
        RETURN
      END IF
 13   call readtest(2,lun,key,type,ierr)
      RETURN
      END
C====
      SUBROUTINE readint4(lun,key,int4)
C---> Read a value (int4) of variable (key) from unit (lun).
C     ======================================================
      CHARACTER*10 key
      CHARACTER*4  type
      CHARACTER*50 text
      call readtest(1,lun,key,type,ierr)
      IF (type.EQ.'INT4') THEN
        READ(lun,FMT='(18X,I8,A50)',ERR=13,END=13) int4,text
        PRINT 2,   key,int4,text
        WRITE(9,2) key,int4,text
  2     FORMAT('      INT4|==> ',a10,'= ',I9,2x,a50)
        RETURN
      END IF
 13   call readtest(2,lun,key,type,ierr)
      RETURN
      END
C====
      SUBROUTINE readch50(lun,key,ch50)
C---> Read a value (ch50) of variable (key) from unit (lun).
C     ======================================================
      CHARACTER*10 key
      CHARACTER*4  type
      CHARACTER*50 ch50
      CHARACTER*12 text
      call readtest(1,lun,key,type,ierr)
      IF (type.EQ.'CH50') THEN
C        READ(lun,FMT='(18X,A50,A12)',ERR=13,END=13) ch50,text
        READ(lun,FMT='(18X,A50,A12)',ERR=13) ch50,text
        PRINT 2,   key,ch50,text
        WRITE(9,2) key,ch50,text
  2     FORMAT('      CH50|==> ',a10,'= ',A50,2x,a12)
        RETURN
      END IF
 13   call readtest(2,lun,key,type,ierr)
      RETURN
      END
C====
      SUBROUTINE readch30(lun,key,ch30)
C---> Read a value (ch30) of variable (key) from unit (lun).
C     ======================================================
      CHARACTER*10 key
      CHARACTER*4  type
      CHARACTER*30 ch30
      CHARACTER*32 text
      call readtest(1,lun,key,type,ierr)
      IF (type.EQ.'CH30') THEN
C        READ(lun,FMT='(18X,A30,A32)',ERR=13,END=13) ch30,text
        READ(lun,FMT='(18X,A30,A32)',ERR=13) ch30,text
        PRINT 2,   key,ch30,text
        WRITE(9,2) key,ch30,text
  2     FORMAT('      CH30|==> ',a10,'= ',A30,2x,a32)
        RETURN
      END IF
 13   call readtest(2,lun,key,type,ierr)
      RETURN
      END
C====
      SUBROUTINE readtest(mode,lun,key,type,ierr)
C---> Test of reading a variable (key) from unit (lun).
C     =================================================
      CHARACTER*10 key,test
      CHARACTER*4  type
      data istart/1/
      IF(istart.eq.1) THEN
        istart=0
        PRINT 10,    lun
        WRITE (9,10) lun
  10    FORMAT(' READ-PAR.|==> Control values from unit (',i2,')',/,
     &         50('='))
      END IF
      IF(mode.eq.1) THEN
      ierr=0
      REWIND(lun)
  1   READ(lun,FMT='(A10,2X,A4)',ERR=13,END=13) test,type
        IF(test.EQ.key) THEN
          BACKSPACE(lun)
          RETURN
        ELSE
          GOTO 1
        END IF
 13   ierr=1
      RETURN
      ELSE
      PRINT 14, key,test,type
 14   FORMAT(' *** ReadErr: key,test,type=',
     &       2a10,1x,a4,i9,/' *** STOP')
      STOP
      END IF
      END

C===============================================================
