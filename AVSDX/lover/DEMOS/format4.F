*DECK COPYRIGHT
************************************************************************
* DEMO program: @(#)format4.F	2.1                      
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
*DECK FORMAT4
      PROGRAM FORMAT4                                                   
C                                                                       
#include "lover.h"
C                                                                       
C     * DATA MATRIX M AND INITIAL RESOLUTION RI                         
      REAL    M(3,4,4,4)                                                
      INTEGER RI(4)                                                     
C                                                                       
C     * RI CONTAINS THE RESOLUTION IN THE FOUR DIMENSIONS               
      DATA RI/3,3,3,1/                                                  
C                                                                       
C     * UNIFORM (PRODUCT ARRAY) GRID -> x(I,J,K)=X(I), y(I,J,K)=Y(J), et
      REAL X(4,4,4),                                                    
     A     Y(4,4,4),                                                    
     B     Z(4,4,4),                                                    
     C     T(1)                                                         
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
      INTEGER   ST(3)                                                   
      CHARACTER FN(3)*22, TY(3)*8                                       
C                                                                       
C     * COUNTERS                                                        
      INTEGER I,J,K                                                     
C                                                                       
C     * LOCAL VARIABLES                                                 
      REAL         WORK(4*3*4*4*4*1)                                         
      INTEGER      IWORK(1*2)                                                
      CHARACTER*23 FILENAME
C                                                                       
C     * M(1,x,y,z)=(x^2+y^2+z^2), M(2,x,y,z)=sqrt(x^2+y^2)              
      DO I=1,RI(3)                                                      
         DO J=1,RI(2)                                                   
            DO K=1,RI(1)                                                
               X(K,J,I)=(K-1)*0.5                                       
               Y(K,J,I)=(J-1)*1.0                                       
               Z(K,J,I)=(I-1)*2.0                                       
               M(1,K,J,I)=X(K,J,I)**2+Y(K,J,I)**2+Z(K,J,I)**2           
               M(2,K,J,I)=SQRT(X(K,J,I)**2+Y(K,J,I)**2)                 
            END DO                                                      
         END DO                                                         
      END DO                                                            
C                                                                       
C     * THE THREE COORDINATE DIRECTIONS WILL BE EXPORTED,               
C       THE TIME DIMENSION NEED NOT BE SLICED SINCE ITS LENGTH          
C       EQUALS 1.                                                       
      CO(1)=1                                                           
      CO(2)=2                                                           
      CO(3)=3                                                           
C                                                                       
C     ** WE DEAL WITH TWO SCALAR FIELDS.                                
      ST(1)=VISRANKSCA                                                  
      ST(2)=VISRANKSCA                                                  
C                                                                       
C     ** THE NAMES OF THE TWO SCALAR FIELDS ARE DEFINED HERE.           
      FN(1)='firstfield'                                                
      FN(2)='secondfield'                                               
C                                                                       
C     ** THE PRIMITIVE TYPES OF THE FIELDS ARE DEFINED HERE.            
      TY(1)=VISTYPEDOU                                                  
      TY(2)=VISTYPEDOU                                                  
C                                                                       
C     * SETTING ARGUMENTS                                               
      CALL VISSETARGS(3,4,4,4,1,RI,CO,SL,3,.FALSE.,2,VISFORMBIN,
     A                VISGEOMAUX)        
C                                                                       
C     * PRODUCING AVS INTERNAL FILE                                     
      FILENAME='AVSINTBIN'
      CALL VISGOFORIT('A',VISINTERNA,M,X,Y,Z,T,FILENAME,             
     A                FN,ST,TY,WORK,IWORK)                              
      END                                                               
