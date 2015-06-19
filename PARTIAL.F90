
MODULE VARIABLES
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:)::XR,YR,ZR,XI,YI,ZI   !Q,ATOM,FRE(3*N)
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::OMEGA,QQ,PDOS,CV
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:)::GGX,GGY,GGZ
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::GX,GY,GZ
double precision,allocatable,dimension(:)::GT,NORMFACTOR,CV_ATOM
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::PDOSNTYP,WPDOSNTYP
INTEGER, ALLOCATABLE,DIMENSION(:)::NNTYPE,NNWTYPE
DOUBLE PRECISION,ALLOCATABLE, DIMENSION(:)::MASS,AMASS,AMASSW,EMEV
double precision,allocatable,dimension(:,:)::norm
DOUBLE PRECISION, ALLOCATABLE,DIMENSION(:)::BB,WEIGHT,TOTAL_NEUTRON_DOS
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::UX2,UY2,UZ2,UU2,UUW2,UWT
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::ENERGY_meV,BROAD_DOS,ANORMGX,ANORMGY,ANORMGZ,UUTOT
DOUBLE PRECISION::PI,S1,S2,TU2,ESTART,TSTART,TEND
CHARACTER::ABC*25,INFILE*100
DOUBLE PRECISION,DIMENSION(700,100)::DFREQ
INTEGER ::NATOM,IWFLAG,IUFLAG,NTYPE,NDF,NFU
INTEGER::NQPOINTS,NSTEP,NWTYPE
INTEGER ::I,i1,N1,N2,II,J,JJ,III,NATM,NWTYPE1,K
INTEGER::NTYPE1,KKK,ISIG,ITEMP
DOUBLE PRECISION::E1,E2,DELE,SUMM,DELE1,ESTEP,DELTAT,FRE
DOUBLE PRECISION::ANORMFACTOR,FWHM1,FWHM,SIG
!ESTART=-1.0  !
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::SPH,TTD,THETAD,THETAT
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::CVD
INTEGER::NDEBYE,JJJ,NDBSTEP
DOUBLE PRECISION::DTHETA,X,DELTASP,DELTADB,TSPSTART,TSPEND,TDBSTART,TDBEND
END MODULE VARIABLES













PROGRAM  PARTIAL
USE VARIABLES
IMPLICIT NONE

CALL READ_OPEN_FILE

IF(ITEMP==1)THEN

        OPEN(1,FILE="TEMP2")
        DO I=1,NSTEP
        READ(1,*)EMEV(I),(PDOSNTYP(J,I),J=1,NTYPE)
        END DO
        CLOSE(1)


        OPEN(1,FILE="TEMP1")
        OPEN(2,FILE="TEMP")
        DO I=1,NSTEP
        READ(1,*)EMEV(I),(PDOS(J,I),J=1,NATOM)
        READ(2,*)EMEV(I),gx(i,:),gy(i,:),gz(i,:),GT(I)
        END DO
        CLOSE(1)

        DELE1=EMEV(2)-EMEV(1)
       
         IF(DELE1.NE.DELE)THEN

          PRINT*,'DIFFERENCE IN ENERGY STEP PLEASE CHECK YOUR TEMP* FILES AND DELE IN PDOS_INP'
        !PAUSE
        END IF
        IF(IWFLAG==1)THEN        

        OPEN(1,FILE="TEMP4")
        DO I=1,NSTEP
        READ(1,*)EMEV(I),(WPDOSNTYP(J,I),J=1,NWTYPE)
        END DO
        CLOSE(1)
        END IF

        IF(IUFLAG==1)THEN
        
        CALL U2(TU2)
        CALL UUT(TSTART,TEND,DELTAT)
        END IF
        
        CALL NEUTRON_WEIGHTED_DOS

        CALL BROAD_DOS_CALC
        
      

     
ELSE


        CALL READ_PHONONQ


        CALL CALCULATE_PARTIAL
        IF(IUFLAG==1)THEN
        CALL U2(TU2)
        CALL UUT(TSTART,TEND,DELTAT)
        END IF

      


        CALL NEUTRON_WEIGHTED_DOS

        CALL BROAD_DOS_CALC
       
END IF          
IF(ITEMP.EQ.0) THEN

        CALL SPECIFIC_HEATV(TSPSTART,TSPEND,DELTASP)

        CALL DEBYE1(TDBSTART,DELTADB,NDBSTEP)
END IF
                               
 END PROGRAM PARTIAL





















!!!!!!!!!!!!!!!!!!!!!!!!!!!!@@@@@@@@@@@@@@@@@@@@@@@SUBROUTINES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE READ_OPEN_FILE
USE VARIABLES
IMPLICIT NONE

!THESE ARE THE THREE BASIC INPUT FILES WHICH HAS TO BE MODIFIED ACCORDING TO SYSTEM AND YOUR CALCULATION REQUIRMENTS   
!OPEN(1,FILE="ZRWO_888.d16")
OPEN(2,FILE="PDOS_INP")
OPEN(3,FILE="WPDOS_INP")
!OPEN(13,FILE="USQUARE_FILE")


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
READ(2,22)INFILE
22 FORMAT(A100)
OPEN(1,FILE=INFILE)
READ(2,*)NATOM,IWFLAG,IUFLAG,ITEMP
PRINT*,NATOM
READ(2,*)NFU
!PAUSE
READ(2,*)NTYPE
READ(2,*)NQPOINTS
READ(2,*)NSTEP,DELE,ESTART    ! DELE In meV mili electron Volt
READ(2,*)TSPSTART,TSPEND,DELTASP
READ(2,*)TDBSTART,DELTADB,NDBSTEP
IF(IWFLAG==1)THEN
READ(3,*)NWTYPE
END IF
IF(IUFLAG==1)THEN
OPEN(4,FILE="T2U")
READ(4,*)TU2,TSTART,TEND,DELTAT
PRINT*,TSTART,TEND
!CALL U2(TU2)
END IF                 !!!!!!!!!!!!!!!!!!!!!!!!!!!TU2 IS THE TEMPERATURE AT WHICH U2 WILL BE CALCULATES.


NDEBYE= INT((TSPEND-TSPSTART)/DELTASP)
NDEBYE=NDEBYE+1
PRINT*,NDEBYE
PAUSE
!READ(2,*)NEND


PI=4.0*atan(1.0)


ALLOCATE(THETAD(NDBSTEP),THETAT(NDBSTEP),SPH(NDEBYE),TTD(NDBSTEP),CVD(NDEBYE,NDBSTEP))
ALLOCATE(CV(NQPOINTS,3*NATOM),CV_ATOM(3*NATOM))
ALLOCATE(NNWTYPE(NWTYPE),WPDOSNTYP(NWTYPE,NSTEP))
ALLOCATE(UUW2(NSTEP,NWTYPE))
ALLOCATE(NORMFACTOR(NTYPE))
ALLOCATE(PDOSNTYP(NTYPE,NSTEP))
ALLOCATE(NNTYPE(NTYPE))
ALLOCATE(GGX(NSTEP,NQPOINTS,NATOM),GGY(NSTEP,NQPOINTS,NATOM),GGZ(NSTEP,NQPOINTS,NATOM),GX(NSTEP,NATOM),GY(NSTEP,NATOM),GZ(NSTEP,NATOM),GT(NSTEP))

ALLOCATE(PDOS(NATOM,NSTEP))
ALLOCATE(XR(NQPOINTS,3*NATOM,NATOM),XI(NQPOINTS,3*NATOM,NATOM),YR(NQPOINTS,3*NATOM,NATOM),YI(NQPOINTS,3*NATOM,NATOM),ZR(NQPOINTS,3*NATOM,NATOM),ZI(NQPOINTS,3*NATOM,NATOM))
ALLOCATE(OMEGA(NQPOINTS,3*NATOM))
ALLOCATE(QQ(NQPOINTS,3))
ALLOCATE(MASS(NATOM))
ALLOCATE(NORM(NQPOINTS,3*NATOM))
ALLOCATE(WEIGHT(NTYPE),BB(NTYPE))
ALLOCATE(TOTAL_NEUTRON_DOS(NSTEP),AMASS(NTYPE),BROAD_DOS(NSTEP))
ALLOCATE(ENERGY_meV(NSTEP))
ALLOCATE(UX2(NSTEP,NATOM),UY2(NSTEP,NATOM),UZ2(NSTEP,NATOM),UU2(NSTEP,NATOM))
ALLOCATE(ANORMGX(NATOM),ANORMGY(NATOM),ANORMGZ(NATOM))
ALLOCATE(AMASSW(NWTYPE))
ALLOCATE(EMEV(NSTEP),UUTOT(NSTEP))
ALLOCATE(UWT(2000,NATOM))
READ(2,*)NNTYPE
READ(2,*)AMASS
READ(2,*)BB       !!!!!!!!!!!!THIS IS THE SCATTERING CROSS SECTION INSTEAD OF BB
BB=SQRT(BB/(4*pi))
PRINT*,BB
PRINT*,BB**2
PAUSE

IF(IWFLAG==1)THEN
READ(3,*)NNWTYPE
PRINT*,NNWTYPE
READ(3,*)AMASSW
PRINT*,AMASSW
END IF


NDF=3*NATOM
PRINT*,NDF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
II=1
!PRINT*,NTYPE


DO I=1,NTYPE
DO J=II,NNTYPE(I)+II-1
MASS(J)=AMASS(I)
END DO
II=II+NNTYPE(I)
END DO
!PRINT*,MASS,AMASS
!PAUSE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

OPEN(4,FILE="TEMP")
OPEN(5,FILE="TEMP1")
OPEN(6,FILE="TEMP2")
OPEN(7,FILE="TEMP3")
OPEN(8,FILE="TEMP4")
OPEN(9,FILE="BROAD")
OPEN(10,FILE="MULTI-IN")
END SUBROUTINE







SUBROUTINE READ_PHONONQ
USE VARIABLES
IMPLICIT NONE
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!!!!!!!!!!!!!!!!!!RAEDING SECTION FOR PHONON FORMAT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

I=0


DO            !!!!!!!!!!!!!!!THIS LOOP FOR Q POINTS
N1=0
N2=0
       READ(1,1,END=99)ABC
         IF (ABC=="  Wave")THEN
!PAUSE

              I=I+1
              BACKSPACE(1)
              READ(1,2)QQ(I,:)
!              PRINT*,QQ(i,:),I

                II=0
             DO J=1,NATOM                             !!!!!!!!!!!!!!!!!!!!!!THIS LOOP IS FOR OMEGA =3*NATOM 
             N1=N2+1
             N2=N1+2
             READ(1,3) (OMEGA(I,II),II=N1,N2)
!PRINT*,OMEGA(i,:)
!PAUSE
            DO JJ=1,NATOM                            !!!!!!!!!!!!!!!!!!!!!!!THIS IS FOR EIGEN VECTOR FOR EACH ATOM FOR A GIVEN Q (I) AND OMEGA(II)

             READ(1,4)(XR(I,II,JJ),XI(I,II,JJ),II=N1,N2)
             READ(1,4)(YR(I,II,JJ),YI(I,II,JJ),II=N1,N2)
             READ(1,4)(ZR(I,II,JJ),ZI(I,II,JJ),II=N1,N2)





             END DO
             READ(1,*)

             END DO
!PRINT*,omega(I,:)
!PAUSE
         END IF
 END DO
99  CONTINUE
1             FORMAT(A6)
2             FORMAT(17X,3(F8.5,3x))
3            FORMAT(8X,3F18.3)
4           FORMAT(14X,6(F7.4,2x))
NQPOINTS=I
PRINT*,NQPOINTS
!#########################@@@@@@@@@@@@@@@@@@@@@######################$$$$$$$$$$$$$$$$$$$$$$$$$$@@@@@@@@@@@@@@@@@@@@@@@@@@@$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!!

END SUBROUTINE

SUBROUTINE CALCULATE_PARTIAL
USE VARIABLES
IMPLICIT NONE
!######@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@CALCULATION OF PARTIAL DOS TEMP CONTAIN PARTIAL DOS OF EACH ATOM WITH X Y and Z COMPONENT@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
I=0


              !0.0
!               DELE=0.2418  ! STEP IN THz    1 Thz=4.136  READ FROM INPUT FILE

               DO I=1,NSTEP   !!!!!!!!!!!!!!FIRST LOOP START

                 GGX=0
                 GGY=0
                 GGZ=0

                 E1=0.2418*ESTART+(I-1)*DELE*0.2418
                 E2=0.2418*ESTART+I*DELE*0.2418  
                 !print*,E1,E2
                     
                      DO II=1,NQPOINTS           !!!!!!!!!SECOND LOOP START

                             DO III=1,3*NATOM     !!!!!!!!!!!THIRD LOOP START
 
                                  IF(OMEGA(II,III).LT.E2.AND.OMEGA(II,III).GE.E1)THEN
                                  !print*,omega(ii,iii)
                                  GGX(I,II,:)=(XR(II,III,:)**2)+(XI(II,III,:)**2)
                                  GGY(I,II,:)=(YR(II,III,:)**2)+(YI(II,III,:)**2)
                                  GGZ(I,II,:)=(ZR(II,III,:)**2)+(ZI(II,III,:)**2)
                                  !print*,ggx(i,ii,:)
                                    !pause
                                  GX(I,:)=GX(I,:)+GGX(I,II,:)
                                  GY(I,:)=GY(I,:)+GGY(I,II,:)
                                  GZ(I,:)=GZ(I,:)+GGZ(I,II,:)
                                  END IF 
                            END DO            !!!!!!!!!!! THIRD LOOP END
                     END DO                   !!!!!!!!!! SECOND LOOP END



END DO   !!!!!!!!!!!!!THIRD  LOOP END 
                                      
                                  DO I1=1,NATOM
                                  ANORMGX(I1)=SUM(GX(:,I1))*DELE
                                  ANORMGY(I1)=SUM(GY(:,I1))*DELE
                                  ANORMGZ(I1)=SUM(GZ(:,I1))*DELE
                                  !PRINT*,ANORMGX(I1),ANORMGY(I1),ANORMGZ(I1)
                                  GX(:,I1)=GX(:,I1)/ANORMGX(I1)
                                  GY(:,I1)=GY(:,I1)/ANORMGY(I1)
                                  GZ(:,I1)=GZ(:,I1)/ANORMGZ(I1)


                                   END DO
                                
                                  




DO I=1,NSTEP    !!!!!!!!!!!AGAIN FIRST  LOOP FOR ENERGY START 

                                  GT(I)=SUM(GX(I,:)+GY(I,:)+GZ(I,:))/(3.0*NATOM)
                                  write(4,5)ESTART+(2*I-1)*DELE/2.0,gx(i,:),gy(i,:),gz(i,:),GT(I)

                                  DO NATM=1,NATOM
                                        PDOS(NATM,I)=(GX(I,NATM)+GY(I,NATM)+GZ(I,NATM))/3.0
                                       
                                        DFREQ(I,NATM)=PDOS(NATM,I)
                                  END DO
                                 
                                  WRITE(5,6)ESTART+(2*I-1)*DELE/2.0,PDOS(:,I),GT(I)/(NATOM)
  

                                  N1=1
                                  N2=0
                             DO NTYPE1=1,NTYPE
                                  N2=N2+NNTYPE(NTYPE1)
                                  !print*,n1,n2
                                        DO II=N1,N2
                                            PDOSNTYP(NTYPE1,I)=PDOSNTYP(NTYPE1,I)+PDOS(II,I)
                                        END DO
                                            N1=N2+1
                         
                          PDOSNTYP(NTYPE1,I)=PDOSNTYP(NTYPE1,I)/NNTYPE(NTYPE1)
                            END DO







END DO
             






                                 

5 format(1000(f10.5,2x))                           
6 FORMAT(100(F10.5,2X))
7 FORMAT(100(F10.5,2X))
! ############################################################################################################################################################################################################


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!HERE WE CREATE INPUT FILE FOR MULTIPHONON PROGRAM  WRITTEN BY. DR. CHAPLOT FILE NAME IS MULTI-IN


                                                     
                                  WRITE(10,57)DFREQ(:,:)
                                   57 FORMAT(8F10.6)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                  DO I=1,NSTEP
                                      E1=0.2418*ESTART+(I-1)*DELE*0.2418
                                      E2=0.2418*ESTART+I*DELE*0.2418
                                         PDOSNTYP(:,I)=PDOSNTYP(:,I)   !/NNTYPE(:)  !/(NORMFACTOR(:))
                                       WRITE(6,7)ESTART+(2*I-1)*0.5*DELE,PDOSNTYP(:,I),SUM(NNTYPE(:)*PDOSNTYP(:,I))/SUM(NNTYPE(:))
                                  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   




!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!THIS IS TO CREATE NEW
            !PARTIAL DOS ACCORDING TO DESIRED wycoff classification of atoms
            !like O1 O2 W1 W2 and so on and the output will be written in TEMP4
            !while in TEMP3
            !you will get PARTIAL DOS according to atom type not by atom site
            !the input for this calculation is  WPDOS_INP !!YOU can change
            !it..its not a big deal...cheers!!!!!!!!
!!!!!!!!SET IWFLAG IN PDOS_INP LAST LINE TO 1 IF YOU WANT TO SEPERATELY CALCULATE WYCOFF SITE PROJECTED PARTIAL DOS



IF (IWFLAG==1)THEN

DO I=1,NSTEP

N1=1
N2=0
                             DO NWTYPE1=1,NWTYPE

                                  N2=N2+NNWTYPE(NWTYPE1)
                                  !print*,n1,n2
                                        DO II=N1,N2
                                            WPDOSNTYP(NWTYPE1,I)=WPDOSNTYP(NWTYPE1,I)+PDOS(II,I)
                                        END DO
                                            N1=N2+1
                             END DO

                                 WRITE(8,7)ESTART+(2*I-1)*0.5*DELE,WPDOSNTYP(:,I)/NNWTYPE(:)

END DO
END IF
!##################################################################################################################################







END SUBROUTINE


SUBROUTINE NEUTRON_WEIGHTED_DOS
USE VARIABLES
IMPLICIT NONE
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!NEUTRON WEIGHTED DOS@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
                                 WEIGHT(:)=4*pi*NNTYPE(:)*(BB(:)**2)/(AMASS(:))
                                 PRINT*,BB(:),AMASS(:)  
                                 PRINT*,NNTYPE,WEIGHT !,SUM(WEIGHT(:))

                                  DO I=1,NSTEP
 
                                 TOTAL_NEUTRON_DOS(I)=SUM(PDOSNTYP(:,I)*WEIGHT(:))

                                 END DO




 ANORMFACTOR=SUM(TOTAL_NEUTRON_DOS(:)*DELE)
TOTAL_NEUTRON_DOS(:)=TOTAL_NEUTRON_DOS(:)/ANORMFACTOR

Do I=1,nstep
write(7,8)ESTART+(2*I-1)*DELE*0.5,TOTAL_neutron_dos(i)
end do
                                

ANORMFACTOR=SUM(TOTAL_NEUTRON_DOS(:)*DELE)
8 FORMAT(F10.4,4X,F15.4)
PRINT*,ANORMFACTOR
!PAUSE
!#######################################################################################################################



END SUBROUTINE


SUBROUTINE BROAD_DOS_CALC
USE VARIABLES
IMPLICIT NONE
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@APPLY BROADENNING IN NEUTRON WEIGTED DOS@@@@@@@@@@@@@@@@@@@@@@@@@@@@
       
!CALL BROAD(TOTAL_NEUTRON_DOS,BROAD_DOS,FWHM,NSTEP)
ESTEP=0.0

         DO I=1,NSTEP

ESTEP=(((2*I-1))*0.5*DELE+ESTART)/(DELE)
IF(ESTEP.LT.0)ESTEP=ABS(ESTEP)

IF(ESTEP.LT.40)THEN
FWHM1=0.15
ELSE
FWHM1=0.15
END IF
         FWHM=FWHM1*ESTEP
!PRINT*,DELE,DELE*4.136,ENERGY_meV(I)
!PAUSE
                 SIG=FWHM/(2*1.1774)
                 ISIG=2*SIG+1
                 IF(FWHM.NE.0.0)S2= 0.5/SIG**2
!         SUM=SQRT(2.0D0*3.141592D0)*SIG
!        DO 12 I=1,N
                 BROAD_DOS(I)=TOTAL_NEUTRON_DOS(I)
                 IF(FWHM==0.0)GO TO 12
                 BROAD_DOS(I)=0.0
                 SUMM=0
                 DO  J=-ISIG,ISIG
                 K=I+J
                 IF(K.LE.0)K=1
                 IF(K.GT.NSTEP)K=NSTEP
                 S1=-S2*J**2
                 S1=DEXP(S1)
                 BROAD_DOS(I)=BROAD_DOS(I)+TOTAL_NEUTRON_DOS(K)*S1
                 SUMM=SUMM+S1
                 END DO
12              CONTINUE
                 BROAD_DOS(I)=BROAD_DOS(I)/SUMM
!                 end if
!       write(56,8)ENERGY_meV(I),BROAD_DOS(I)
         END DO
                                  !NORMALAZATION OF BROAD DATA

                                      ANORMFACTOR=SUM(BROAD_DOS(:)*DELE)
!PRINT*,ANORMFACTOR
                                                                       
                                  DO I=1,NSTEP
                                       ESTEP=(((2*I-1))*0.5*DELE+ESTART)/(DELE)                               
                                       BROAD_DOS(I)=BROAD_DOS(I)/ANORMFACTOR
                                       WRITE(9,8)ESTEP*DELE,BROAD_DOS(I)
                                  END DO



ANORMFACTOR=SUM(BROAD_DOS(:)*dele)
!PRINT*,ANORMFACTOR
8 FORMAT(F10.4,4X,F15.4)
END SUBROUTINE



SUBROUTINE U2(TT)
USE VARIABLES
IMPLICIT NONE
DOUBLE PRECISION,EXTERNAL::NBOSE1
DOUBLE PRECISION::TT,H
OPEN(11,FILE="U2")
OPEN(12,FILE='UU2TOT')
H=6.626068E-26/1.66053886E-27

Do J=1,NSTEP
!DO I=1,NATOM

FRE=ESTART*0.2418+(2*J-1)*DELE*0.2418/2.0
UX2(J,:)=(NBOSE1(FRE,TT)+0.5)*H*GX(J,:)/(3*4*PI*PI*MASS(:)*FRE)
UY2(J,:)=(NBOSE1(FRE,TT)+0.5)*H*GY(J,:)/(3*4*PI*PI*MASS(:)*FRE)
UZ2(J,:)=(NBOSE1(FRE,TT)+0.5)*H*GZ(J,:)/(3*4*PI*PI*MASS(:)*FRE)
UU2(J,:)=UX2(J,:)+UY2(J,:)+UZ2(J,:)
!END DO
UU2(J,:)=(NBOSE1(FRE,TT)+0.5)*H*PDOS(:,J)/(4*PI*PI*MASS(:)*FRE)  

!!!!!!!!!!!UUTOT IS TOTAL MASS WEIGHTED U2!!!!!!!!!!!!!!!!!!IT CAN BE USED AS A GUIDELINE OF MINIMUM AMPLITUDE INITIAL GUESS FOR POTENTIAL CALCULATION IN FROZEN PHONON CALCULATION

UUTOT(J)=(NBOSE1(FRE,TT)+0.5)*H*GT(J)/(4*PI*PI*FRE)
!PRINT*,J,UX2(J,:)
!PAUSE
WRITE(11,11)FRE*4.136,UU2(J,:)
WRITE(12,11)FRE*4.136,UUTOT(J)
11 FORMAT(1000F12.6)

END DO


PRINT*,'TOTAL MASS WEIGHTED U2 ',SUM(UUTOT(:))*DELE
PRINT*,'UU2 ATOM WISE IN ANGSTROM SQUARE',(SUM(UU2(:,I)),I=1,NATOM)
PAUSE
!WRITE(13,13)TT,(SUM(UU2(:,I)*DELE),I=1,NATOM)
!13 FORMAT(I6,1000F10.4)



IF (IWFLAG==1)THEN
OPEN(12,FILE="U2WYCOFF")
DO I=1,NSTEP
FRE=ESTART*0.2418+(2*I-1)*DELE*0.2418/2.0

UUW2(I,:)=(NBOSE1(FRE,TT)+0.5)*H*WPDOSNTYP(:,I)/(4*PI*PI*AMASSW(:)*FRE) 
WRITE(12,11)FRE*4.136,UUW2(I,:)
END DO
END IF
END SUBROUTINE
!#######################################################################################

SUBROUTINE UUT(T1,T2,DELT)

 USE VARIABLES
IMPLICIT NONE
DOUBLE PRECISION,EXTERNAL::NBOSE1
DOUBLE PRECISION::TT,H ,T1,T2,DELT
!PRINT*,T1,T2
!PAUSE
TT=T1-DELT
H=6.626068E-26/1.66053886E-27
OPEN(14,FILE="USQUARET")

DO  
!IF(TT>T2)EXIT
TT=TT+DELT
!PRINT*,TT,T1,T2
!PAUSE
        IF(TT>T2)EXIT
        DO I=1,NSTEP
        FRE=ESTART*0.2418+(2*I-1)*DELE*0.2418/2.0

        UUW2(I,:)=(NBOSE1(FRE,TT)+0.5)*H*WPDOSNTYP(:,I)/(4*PI*PI*AMASSW(:)*FRE)


        END DO

        DO J=1,NWTYPE
        UWT(TT,J)=SUM(UUW2(:,J)) !*DELE 

        END DO
        WRITE(14,14)TT,UWT(TT,:)    

14      FORMAT(F10.5,1000F10.4)


END DO






END SUBROUTINE







!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



!!!!!!!BOSE POPULATON FACTOR CALCULATION!!!!!!!BOSE POPULATON FACTOR
!CALCULATION!!!!!!!BOSE POPULATON FACTOR CALCULATION!!!!!!!BOSE POPULATON FACTOR
!CALCULATION!!!!!!!BOSE POPULATON FACTOR CALCULATION
DOUBLE PRECISION  FUNCTION  NBOSE1(FRE,TT)
DOUBLE PRECISION::FRE, TT
BETA1=47.99237243603648754
BETAT=BETA1/TT
!PRINT*,TT,BETAT,BETA1 
!PAUSE
NBOSE1=1/(DEXP(BETAT*FRE)-1.0)
!PRINT*,TT,BETAT,FRE,NBOSE1
!PAUSE
END FUNCTION NBOSE1
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


!!!SPECIFIC HEAT CALCULATION!!!SPECIFIC HEAT CALCULATION!!!SPECIFIC HEAT
!CALCULATION!!!SPECIFIC HEAT CALCULATION!!!SPECIFIC HEAT CALCULATION!!!SPECIFIC
!HEAT CALCULATION!!!SPECIFIC HEAT CALCULATION
SUBROUTINE SPECIFIC_HEATV(T1,T2,DELT)
USE VARIABLES
IMPLICIT NONE
DOUBLE PRECISION::TEMP,T1,T2,DELT,CVT,H,KB,BETA1
DOUBLE PRECISION,EXTERNAL::DNBOSE1
H=6.626069300000000159e-22                     ! ORDER OF 12 LESS BECAUSE WHEN WE MULTIPLY W IT IS IN THZ
KB=1.380650500000000147e-23
BETA1=47.99237243603648754                                    !HERE WE HAVE MULTIPLY IT BY E12 SINCE FREQUENCY IS IN THZ

!OPEN(33,FILE="SPECIFIC_HEAT")
JJJ=0


DO II=T1,T2,DELT
     JJJ=JJJ+1
        TEMP=T1+(JJJ-1)*DELT

        DO J=1,NDF
        DO I=1,NQPOINTS

        FRE=OMEGA(I,J)

        IF(FRE==0.0)THEN
        CV(I,J)=0.0
        ELSE
        CV(I,J)=H*FRE*(DNBOSE1(FRE,TEMP))*6.023E23/NFU           !/(V0*0.14824687E-30)

        END IF
        END DO
        CV_ATOM(J)=SUM(CV(:,J))/NQPOINTS
        END DO
        CVT=SUM(CV_ATOM(:))
        SPH(JJJ)=CVT
END DO
END SUBROUTINE

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!!!!!!!!!DERIVATIVE OF BOSE FUNCTION!!!!!!!!!!!!!!!DERIVATIVE OF BOSE
!FUNCTION!!!!!!!!!!!!!!!DERIVATIVE OF BOSE FUNCTION!!!!!!!!!!!!!!!DERIVATIVE OF
!BOSE FUNCTION!!!!!!!!!!!!!!!DERIVATIVE OF BOSE FUNCTION!!!!!
DOUBLE PRECISION FUNCTION DNBOSE1(FRE,TT)
IMPLICIT NONE
DOUBLE PRECISION::FRE,BETAT,BETA1,TT
BETA1=47.99237243603648754
BETAT=BETA1/TT


IF(FRE==0)THEN
DNBOSE1=0
ELSE
!PRINT*,TT,BETAT,BETA1,FRE
DNBOSE1=BETAT*FRE*DEXP(BETAT*FRE)/(TT*(DEXP(FRE*BETAT)-1.0)*(DEXP(FRE*BETAT)-1.0))
!                         (DEXP(BETAT*FRE)-1.0)**2.0)
END IF
END FUNCTION DNBOSE1


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

SUBROUTINE DEBYE1(TDSTART,TDSTEP,NDSTEP)

USE VARIABLES
IMPLICIT NONE
INTEGER::NDSTEP
DOUBLE PRECISION::TD,T,XM,TDSTART,TDSTEP,SS1

OPEN(36,FILE='DEBYE')

OPEN(37,FILE="DEBYET")

!PRINT*,SPH

PRINT*,TDSTEP,TDSTART,NDSTEP 
PAUSE

JJJ=0

DO I=TSPSTART,TSPEND,DELTASP

JJJ=JJJ+1

!TTD(JJJ)=I

!PAUSE


       DO J=1,NDSTEP

        THETAD(J)=TDSTART+J*TDSTEP

        TD=THETAD(J)

        

        XM=REAL(TD/I)

        TTD(J)=XM

        JJ=0.0

        X=0

       IF(XM.GT.70)GO TO 199 

               DO

                JJ=JJ+1

                X=JJ*XM/10000  !-0.0005

!                X=X+0.01

 !               PRINT*,X

                IF(X.GT.XM) GO TO 199

                CVD(JJJ,J)=CVD(JJJ,J)+(X**4)*EXP(X)/(EXP(X)-1.0)**2
                
                END DO
PRINT*,JJ
                199 CONTINUE

          

                CVD(JJJ,J)=NATOM*CVD(JJJ,J)*9*8.31*0.0001*XM/(NFU*(XM**3))

!          PRINT*,X,XM,CVD(JJJ,J)

!          PAUSE

                END DO

  !              WRITE(37,36)I,(TTD(J),J=1,50)

                WRITE(36,36)I,(CVD(JJJ,J),J=1,NDSTEP)

36              FORMAT(I5,2x,20000F8.2)

37              FORMAT(I5,100F15.2)                  

END DO



JJJ=0

DO I=TSPSTART,TSPEND,DELTASP



JJJ=JJJ+1

DO J=1,NDSTEP

IF(SPH(JJJ).LT.CVD(JJJ,J))THEN
SS1=(SPH(JJJ)-CVD(JJJ,J-1))/(CVD(JJJ,J)-CVD(JJJ,J-1))

!DTHETA=THETAD(J-1)+SS1*(THETAD(J)-THETAD(J-1))  !J*20
DTHETA=THETAD(J)
END IF



END DO

WRITE(37,37)I,SPH(JJJ),DTHETA,I/DTHETA  !THETAD(JJJ)

END DO

!PRINT*,SPH



END SUBROUTINE
!!!!!!!!!!!!!!!!THIS PROGRAM HAS TO COMPLETE
