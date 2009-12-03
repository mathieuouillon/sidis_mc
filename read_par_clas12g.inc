c    include file to read parameters
c      CHARACTER*30 parfile 
      INTEGER Cnfg_num
      REAL 
     &Tmax,Thf_lim,
     &Thacc1,Thacc2,Thacc3,Thacc4,Thacc5,Thacc6,
     &Thacc01,Thmom1,Thmom2,Thmom3,Thmom4,
     &Phimom1,Phimom2,Phimom3,Phimom4,Phimom5,
     +FoMT1S1,FoMT1S2,FoMT1S3,FoMT2S1,FoMT2S2,FoMT2S3,
     &FoInMoS1,FoInMoS2,FoOuMoS1,FoOuMoS2,FoNeMoS1,FoNeMoS2,
     &FoInThS1,FoInThS2,FoOuThS1,FoOuThS2,FoNeThS1,FoNeThS2,
c     &FoInFiS1,FoInFiS2,FoOuFiS1,FoOuFiS2,FoNeFiS1,FoNeFiS2,
cc     &FoIOFiS1,FoA0FiS2,FoA1FiS2,FoA2FiS2,FoNeFiS1,FoNeFiS2,
     &FoA0FiS1,FoA1FiS1,FoA2FiS1,
     +FoA0FiS2,FoA1FiS2,FoA2FiS2,FoNeFiS1,FoNeFiS2,
     *CeA0ThS2,CeA1ThS2,CeA2ThS2,
     &CeInMoS1,CeInMoS2,CeOuMoS1,CeOuMoS2,CeNeMoS1,CeNeMoS2,
     &CeInThS1,CeInThS2,CeOuThS1,CeOuThS2,CeNeThS1,CeNeThS2,
     &CeInFiS1,CeInFiS2,CeOuFiS1,CeOuFiS2,CeNeFiS1,CeNeFiS2

      COMMON/clas12data/Cnfg_num,
     &Tmax,Thf_lim,
     &Thacc1,Thacc2,Thacc3,Thacc4,Thacc5,Thacc6,
     &Thacc01,Thmom1,Thmom2,Thmom3,Thmom4,
     &Phimom1,Phimom2,Phimom3,Phimom4,Phimom5,
     +FoMT1S1,FoMT1S2,FoMT1S3,FoMT2S1,FoMT2S2,FoMT2S3,
     &FoInMoS1,FoInMoS2,FoOuMoS1,FoOuMoS2,FoNeMoS1,FoNeMoS2,
     &FoInThS1,FoInThS2,FoOuThS1,FoOuThS2,FoNeThS1,FoNeThS2,
c     &FoInFiS1,FoInFiS2,FoOuFiS1,FoOuFiS2,FoNeFiS1,FoNeFiS2,
cc     &FoIOFiS1,FoA0FiS2,FoA1FiS2,FoA2FiS2,FoNeFiS1,FoNeFiS2,
     &FoA0FiS1,FoA1FiS1,FoA2FiS1,
     +FoA0FiS2,FoA1FiS2,FoA2FiS2,FoNeFiS1,FoNeFiS2,
     *CeA0ThS2,CeA1ThS2,CeA2ThS2,
     &CeInMoS1,CeInMoS2,CeOuMoS1,CeOuMoS2,CeNeMoS1,CeNeMoS2,
     &CeInThS1,CeInThS2,CeOuThS1,CeOuThS2,CeNeThS1,CeNeThS2,
     &CeInFiS1,CeInFiS2,CeOuFiS1,CeOuFiS2,CeNeFiS1,CeNeFiS2
