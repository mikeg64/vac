*DECK COPYRIGHT
************************************************************************
* DEMO program: @(#)timesteps2.F	2.1                      
*                                                                      *
* Author: A.J.C. Belien                                                *
*         FOM-Instituut voor Plasmafysica `Rijnhuizen'                 *
*         P.O. Box 1207                                                *
*         3430 BE Nieuwegein                                           *
*                                                                      *
* (C) 1994-1996: FOM-Institute for Plasma Physics `Rijnhuizen',        *
*                Academic Computing Services Amsterdam (SARA)          *
* (C) 1997-2000: FOM-Institute for Plasma Physics `Rijnhuizen'         *
************************************************************************
*DECK TIMESTEPS2
      PROGRAM TIMESTEPS2                                                
C                                                                       
#include "lover.h"
C                                                                       
C     * DATA MATRIX M AND INITIAL RESOLUTION RI.                        
C     * M(x,t) CONTAINS A SCALAR FIELD.                                 
      REAL    M(4)                                                      
      INTEGER RI(4)                                                     
C                                                                       
C     * INITIAL VALUES                                                  
      DATA RI/4,1,1,4/                                                  
C                                                                       
C     * UNIFORM (PRODUCT ARRAY) GRID -> x(I,J,K)=X(I), y(I,J,K)=Y(J), et
      REAL X(4),                                                        
     A     Y(4),                                                        
     B     Z(4),                                                        
     C     T(4)                                                         
C                                                                       
C     * CO: ARRAY CONTAINING WHICH DIMENSIONS WILL BE WRITTEN OUT.      
C     * SL: ARRAY CONTAINING THE PLANE NUMBERS OF THE SLICED DIMENSIONS.
      INTEGER CO(4), SL(4)                                              
C                                                                       
C     * ST: ARRAY CONTAINING THE RANKS OF THE FIELDS (STRUCTURES).      
C     * TY: ARRAY CONTAINING THE TYPES OF THE FIELDS (STRUCTURES).      
C     * FN: ARRAY CONTAINING THE NAMES OF THE FIELDS (STRUCTURES).      
C                                                                       
C     * ALL THREE ARRAYS MUST HAVE A DIMENSION OF NUMBER OF FIELDS + 1  
      INTEGER   ST(2)                                                   
      CHARACTER FN(2)*22, TY(2)*8                                       
C                                                                       
C     * COUNTERS                                                        
      INTEGER I,J
C                                                                       
C     * LOCAL VARIABLES                                                 
      REAL         WORK(4*1*4*1*1*1)                                         
      INTEGER      IWORK(4*1)                                                
      CHARACTER*23 FILENAME
C                                                                       
C     * INITIAL VALUES                                                  
      DATA TIMESTEPS/4/                                                 
C                                                                       
C     * THE SECOND COMPUTATIONAL DEPENDENCE AND THE THIRD DEPENDENCE    
C       ARE SLICED.                                                     
      CO(1)=1                                                           
C                                                                       
C     ** THE SLICED PLANES ARE DEFINED HERE.                            
      SL(2)=1                                                           
      SL(3)=1                                                           
C                                                                       
C     ** WE DEAL WITH ONE SCALAR FIELDS.                                
      ST(1)=VISRANKSCA                                                  
C                                                                       
C     ** THE NAMES OF THE SCALAR FIELD IS DEFINED HERE.                 
      FN(1)='firstfield'                                                
C                                                                       
C     ** THE PRIMITIVE TYPE OF THE FIELD IS DEFINED HERE.               
      TY(1)=VISTYPEDOU                                                  
C                                                                       
C     * M(x,t)=sin(k*x-omega*t); k=1.0, omega=4.0                       
      DO J=1,RI(1)                                                      
         X(J)=2*PI*(J-1)/3.0                                            
         Y(J)=0.0                                                       
         Z(J)=0.0                                                       
      END DO                                                            
C                                                                       
      DO I=1,TIMESTEPS                                                  
         T(I)=0.5*PI*(I-1)/3.0                                          
      END DO                                                            
C                                                                       
C     * SETTING ARGUMENTS                                               
      CALL VISSETARGS(1,4,1,1,4,RI,CO,SL,1,.TRUE.,1,VISFORMASC,
     A                VISGEOMAUX)         
C                                                                       
C     * WE WILL OUTPUT ONE TIMESTEP A TIME                              
      CALL VISONESTEP                                                   
C                                                                       
C     * PRODUCING DX INTERNAL ASCII FILE, FIRST TIME STEP               
      CALL VISGOFORIT('D',VISINTERNA,M,X,Y,Z,T,'TIMESTEPS2.dx',         
     A                FN,ST,TY,WORK,IWORK)                              
C                                                                       
      DO I=1,TIMESTEPS                                                  
C                                                                       
C        ** CALCULATING M FOR THIS TIMESTEP                             
         DO J=1,RI(1)                                                   
            M(J)=SIN(1.0*X(J)-4.0*T(I))                                 
         END DO                                                         
C                                                                       
C        ** WRITING OUT THE TIMESTEP                                    
         FILENAME='TIMESTEPS2.dx'
         CALL VISGOFORIT('D',VISINTERNA,M,X,Y,Z,T,FILENAME,
     A                   FN,ST,TY,WORK,IWORK)                           
      END DO                                                            
      END                                                               
