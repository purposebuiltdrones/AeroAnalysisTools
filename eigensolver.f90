!***********************************************************************
      SUBROUTINE EIGENSOLVER(IV,N,Ar,Br,EVAL,V,IERR)                  !*
!     SUBROUTINE EIGENSOLVER(IV,N,A,B,EVAL,V,IERR) Uncomment for complex A+B matrices
!***********************************************************************   
!* AUTHOR                                                              * 
!*         This main body of this code was written in 1979 by the      *
!*         COMPUTER SCIENCES CORPORATION of HAMPTON, VA.  The driver   *
!*         was modified in 2006 by N.R. Alley, Utah State University.  *
!*         Modification include:                                       *
!*             1. Conversion to FORTRAN 90 format                      *
!*             2. User friendy i/o:                                    *
!*                  -Input A+B matrices are real (this can be easily   *
!*                     changed to input complex A+B matrices by modi-  *
!*                     fying call line, and removing line 116)         *
!*                  -Eigenvalues are computed internally and sorted    *
!*                     in desceding order according to the magnitude   *
!*                     of the real component                           *
!*                  -Eivenvectors are defined to a magnitude of unity  *
!*             3. Ability to compute only the performane index (IV=3)  *
!*                  was removed                                        *
!*             4. Code modified to Double Precision                    *
!*                                                                     *
!* PURPOSE                                                             *
!*         TO CALCULATE THE EIGENVALUES AND, OPTIONALLY, THE EIGEN-    *
!*         VECTORS AND/OR PERFORMANCE INDEX OF THE GENERALIZED EIGEN-  *
!*         PROBLEM A*X = LAMBDA*B*X, WHERE A AND B ARE MATRICES OF     *
!*         ORDER N.                                                    *
!*                                                                     *
!*         CODE MUST BE COMPILED IN DOUBLE PRECISION                   *
!*                                                                     *
!* INPUTS															   *
!*                                                                     *
!*         IV  AN INPUT INTEGER OPTION PARAMETER                       *
!*                IV = 0, COMPUTE EIGENVALUES ONLY.                    *
!*                IV = 1, COMPUTE EIGENVALUES AND EIGENVECTORS.        *
!*                IV = 2, COMPUTE EIGENVALUES, EIGENVECTORS, AND       * 
!*                        PERFORMANCE INDEX.                           * 
!*                                                                     * 
!*          N  AN INPUT INTEGER SPECIFYING THE ORDER OF MATRICES       * 
!*             A AND B.                                                * 
!*                                                                     * 
!*         Ar  AN INPUT TWO-DIMENSIONAL ARRAY WITH DIMENSION NxN.      * 
!*             A MUST BE TYPED REAL AND ON INPUT, MUST CONTAIN THE     *
!*             ELEMENTS OF THE COEFFICIENT MATRIX A.  IF IV < 2, A     *
!*             IS DESTROYED ON OUTPUT.                                 * 
!*                                                                     * 
!*         Br  AN INPUT TWO-DIMENSIONAL ARRAY WITH DIMENSION NxN.      * 
!*             B MUST BE TYPED REAL AND ON INPUT, MUST CONTAIN THE     *
!*             ELEMENTS OF THE COEFFICIENT MATRIX B.  IF IV < 2, B     *
!*             IS DESTROYED ON OUTPUT.                                 * 
!*                                                                     *
!* OUTPUTS                                                             *
!*                                                                     *
!*       EVAL  AN OUTPUT ONE-DIMENSIONAL ARRAY CONTAINING THE EIGEN-   *
!*             VALUES WHICH MUST BE TYPED COMPLEX AND DIMENSIONED AT N *
!*                                                                     *
!*          V  AN INPUT/OUTPUT TWO-DIMENSIONAL ARRAY WHICH MUST BE     * 
!*             TYPED COMPLEX AND DIMENSIONED NxN.  IF IV = 1 OR 2,     *
!*             V IS OUTPUT AND THE j-TH COLUMN OF V CONTAINS THE       *
!*             EIGENVECTOR CORRESPONDING TO THE j-TH EIGENVALUE        *
!*             WHERE V(i,j)                                            *
!*             IF IV = 0, V MAY BE A DUMMY PARAMETER.                  *                               
!*                                                                     *
!*       IERR  AN OUTPUT INTEGER ERROR PARAMETER.                      * 
!*                =0  NO ERROR DETECTED.                               * 
!*                =J  THE ROUTINE FAILED TO CONVERGE TO THE J-TH       * 
!*                    EIGENVALUE AFTER 30 ITERATIONS.  EIGENVALUES     * 
!*                    J+1,J+2,...,N ARE CORRECT BUT EIGENVALUES        * 
!*                    1,2,...,J MAY BE INACCURATE.  EIGENVECTORS ARE   * 
!*                    NOT CORRECT.  PERFORMANCE INDEX IS SET TO 1000   * 
!*                    IF IT IS REQUESTED.                              *
!*                =-1 IV IS LESS THAN 0 OR GREATER THAN 2.             * 
!*                                                                     *
!* PARAMETERS                                                          *
!*                                                                     *
!*      ALPHA  AN INPUT/OUTPUT ONE-DIMENSIONAL ARRAY WHICH MUST BE     * 
!*             TYPED COMPLEX AND DIMENSIONED AT LEAST N.  IF IV < 3,   * 
!*             ALPHA IS OUTPUT AND CONTAINS THE NUMERATORS OF THE      * 
!*             RATIOS REPRESENTING THE EIGENVALUES.  THE J-TH EIGEN-   * 
!*             VALUE IS EQUAL TO ALPHA(J)/BETA(J) FOR J=1,2,...,N.     * 
!*             IF IV = 3, ALPHA IS INPUT AND MUST CONTAIN THE NUMER-   * 
!*             ATORS OF THE RATIOS DESCRIBING THE EIGENVALUES.         * 
!*                                                                     * 
!*       BETA  AN INPUT/OUTPUT ONE-DIMENSIONAL ARRAY WHICH MUST BE     * 
!*             TYPED COMPLEX AND DIMENSIONED AT LEAST N.  IF IV < 3,   * 
!*             BETA IS OUTPUT AND CONTAINS THE DENOMINATORS OF THE     *
!*             RATIOS REPRESENTING THE EIGENVALUES.  THE J-TH EIGEN-   * 
!*             VALUE IS EQUAL TO ALPHA(J)/BETA(J) FOR J=1,2,...,N.     * 
!*             IF IV = 3,  BETA IS INPUT AND MUST CONTAIN THE DENOM-   * 
!*             INATORS OF THE RATIOS DESCRIBING THE EIGENVALUES.       * 
!*                                                                     * 
!*         WK  A COMPLEX ONE-DIMENSIONAL WORKING STORAGE ARRAY.  IF    * 
!*             IV < 2, WK MAY BE A DUMMY PARAMETER AND NEED NOT BE     * 
!*             DIMENSIONED.  IF IV = 2, WK MUST CONTAIN AT LEAST       * 
!*             2*N**2 LOCATIONS.  IF IV = 3, WK MUST CONTAIN AT LEAST  * 
!*             ONE  LCATION.  IF IV > 1, THE REAL PART OF THE FIRST    * 
!*             ELEMENT OF WK WILL CONTAIN THE PERFORMANCE INDEX ON     * 
!*             OUTPUT.  THE ROUTINE IS CONSIDERED TO HAVE PERFORMED    * 
!*             (WELL, SATISFACTORILY, POORLY) IF THE PERFORMANCE INDEX * 
!*             IS (LESS THAN 1, BETWEEN 1 AND 100, GREATER THAN 100),  * 
!*             RESPECTIVELY.                                           * 
!*                                                                     *
!***********************************************************************
IMPLICIT NONE
INTEGER::N,JJOB,IERR,JER,KER,iMAX,IV,IVP1,I,J,K,L
REAL::Ar(N,N),Br(N,N)
COMPLEX::A(N,N),B(N,N),EVAL(N),V(N,N),CZERO,SUMR,SUMS
COMPLEX,ALLOCATABLE,DIMENSION(:)::ALPHA,BETA,HOLD2  
COMPLEX,ALLOCATABLE,DIMENSION(:,:)::WK                                      
REAL::ANORM,ASUM,PI,S,SUMR2,ZNORM,BSUM,REPS,ZERO,ONE,THOUS,BNORM
DATA CZERO/(0.d0,0.d0)/                              
DATA REPS/1.0e30/,ZERO/0.d0/,ONE/1.d0/,THOUS/1000.d0/  

!-ALLOCATE WORK ARRAYS
	ALLOCATE(WK(2*N**2,2*N**2),ALPHA(N),BETA(N),HOLD2(N))
!-INITIALIZE A and B WORK ARRAYS
	A=Ar ; B=Br              !Delete this line if A and B become inputs
!-INITIALIZE PARAMETERS 
	JJOB = 0 ; IERR = 0 ; JER = 0 ; KER = 0 ; iMAX=N
!-CHECK - IV OUT OF RANGE 
	IF(IV .GE. 0 .AND. IV .LE. 2) GO TO 15                            
		KER = -1 ; GO TO 105                                                         
   15 IF(IV .EQ. 1 .OR. IV .EQ. 2) JJOB=1                                
      IF(IV .EQ. 3) JJOB=2                                               
      IF(IV-2)30,20,35                                                   
   20 CONTINUE                                                           
!-SAVE INPUT A AND B IF IV = 2          
    do J=1,N                                                        
		do I=1,N                                                        
			WK(I,J) = A(I,J)                                                
			WK(I,J+N) = B(I,J)                                             
		enddo
	enddo
   30 CONTINUE                                                           
      CALL QXZ165(N,A,iMAX,B,iMAX,JJOB,V,iMAX)                         
   35 CALL QXZ166(N,A,iMAX,B,iMAX,V,iMAX,JJOB,ALPHA,BETA,IERR)         
      IVP1 = IV+1                                                       
      GO TO (105,85,40,50),IVP1                                          

   40 CONTINUE                                                          
!-MOVE ORIGINAL MATRICES BACK TO A AND B                          
      DO 45 J=1,N                                                        
      DO 45 I=1,N                                                        
         A(I,J) = WK(I,J)                                               
         B(I,J) = WK(I,J+N)                                             
   45 CONTINUE                                                          
      WK(1,1) = CMPLX(THOUS,ZERO)                                        
      IF(IERR .NE. 0) GO TO 105                                          
 50   ANORM = ZERO                                                       
      BNORM = ZERO                                                       
      DO 60 J=1,N                                                        
          ASUM = ZERO                                                    
          BSUM = ZERO                                                    
          DO 55 I=1,N                                                    
             ASUM = ASUM+ABS(A(I,J))                                    
             BSUM = BSUM+ABS(B(I,J))                                    
 55       CONTINUE                                                      
          ANORM = MAX(ANORM,ASUM)                                      
          BNORM = MAX(BNORM,BSUM)                                     
 60       CONTINUE                                                      
          IF (ANORM.EQ.ZERO) ANORM = ONE                                 
          IF (BNORM.EQ.ZERO) BNORM = ONE                                 
!-COMPUTE PERFORMANCE INDEX
      PI = ZERO                                                          
      DO 80 J=1,N                                                        
         ZNORM = ZERO                                                    
         DO 65 I=1,N                                                    
            ZNORM = MAX(ABS(V(I,J)),ZNORM)                               
   65    CONTINUE                                                        
         S = ZERO                                                        
         DO 75 L=1,N                                                     
            SUMR = CZERO                                                 
            SUMS = CZERO                                                 
            DO 70 K=1,N                                                  
               SUMR = SUMR+A(L,K)*V(K,J)                                
               SUMS = SUMS+B(L,K)*V(K,J)                                 
   70       CONTINUE                                                     
            SUMR = BETA(J)*SUMR-ALPHA(J)*SUMS                            
            S = MAX(S,ABS(SUMR))                                         
   75    CONTINUE                                                        
         SUMR2 = ABS(ALPHA(J))*BNORM*ZNORM                               
         SUMR2 = ABS(BETA(J))*ANORM*ZNORM+SUMR2                         
         PI = MAX(PI,S/SUMR2)                                            
   80 CONTINUE                                                           
      PI = PI/REPS                                                       
      WK(1,1) = CMPLX(PI,ZERO)                                           
   85 CONTINUE                                                          
                                                                       
!-NORMALIZE VECTORS TO 2-NORM OF ONE  
      DO 100 J=1,N                                                       
         S = ZERO                                                       
         DO 90 I=1,N                                                     
            S = S + REAL(V(I,J))**2 + IMAG(V(I,J))**2                   
   90    CONTINUE                                                       
         S = SQRT(S)                                                   
         SUMS = CMPLX(S,ZERO)                                        
         DO 95 I=1,N                                                   
            V(I,J) = V(I,J)/SUMS                                        
   95    CONTINUE                                                       
  100 CONTINUE                                                          
  105 CONTINUE                                                          
      IF (IERR .EQ. 0) IERR = KER                                        
	  EVAL=ALPHA/BETA
	  Ar=real(A) ; Br=real(B)
!-Sort Eigenvalues in descending order-!
	  call sorteg(N,EVAL,V)
!-Define Eigenvectors so that each Eigenvector has a magnitude of unity-!
!-See Equation 8.1.37, pg 704, Mechanics of Flight, W.F. Phillips, 2004-!
      call normcv(N,V)
	  DEALLOCATE(WK,ALPHA,BETA,HOLD2)
 9005 RETURN                                                            
      END      
!**********************************************************************
SUBROUTINE normcv(n,v)
!**********************************************************************
! normalize the eigenvectors to a root-mean-square norm of one,
! with the largest component real.
COMPLEX v(n,n),vmax
do j=1,n
	rmag=0.
	rmax=0.
	max=1
	do i=1,n
		vmag=abs(v(i,j))
		rmag=rmag+vmag*vmag
		if(vmag.gt.rmax)then
			rmax=vmag
			max=i
		endif
	end do
	rmag=sqrt(rmag)
	vmax=rmag*v(max,j)/abs(v(max,j))
	do i=1,n
		v(i,j)=v(i,j)/vmax
	end do
end do
return
end
!**********************************************************************
SUBROUTINE sorteg(n,e,v)
!**********************************************************************
! Sort the complex eigenvalues and eigenvectors by eigenvalues magnitude
COMPLEX e(n),v(n,n),tmp
do i=1,n-1
	max=i
	do j=i+1,n
		if(abs(e(j)).gt.abs(e(max)))then
			max=j
		else if(abs(abs(e(max))-abs(e(j))).lt..000001)then
			if(imag(e(j)).gt.imag(e(max)))max=j
		endif
	end do
	if(abs(abs(e(max))-abs(e(i))).lt..000001)then
		if(imag(e(i)).gt.imag(e(max)))max=i
	endif
	if(max.ne.i)then
		tmp=e(i)
		e(i)=e(max)
		e(max)=tmp
		do k=1,n
			tmp=v(k,i)
			v(k,i)=v(k,max)
			v(k,max)=tmp
		end do
	endif
end do
return
end                                                          
!**********************************************************************
      SUBROUTINE QXZ165 (N,A,IA,B,IB,IJOB,Z,IZ)                         
!**********************************************************************
!                                                                       !ELZHCL 4
!                                                                       !ELZHCL 5
!   FUNCTION            - REDUCE TWO COMPLEX MATRICES, A AND B, SIMUL-  !ELZHCL 6
!                           TANEOUSLY, A TO UPPER HESSENBERG AND B TO   !ELZHCL 7
!                           UPPER TRIANGULAR FORM.                      !ELZHCL 8
!   USAGE               - CALL QXZ165(N,A,IA,B,IB,IJOB,Z,IZ)            !ELZHCL 9
!   PARAMETERS   N      - THE ORDER OF A,B, AND Z. (INPUT)              !ELZHCL10
!                A      - ON INPUT, A IS THE COMPLEX MATRIX TO BE       !ELZHCL11
!                           REDUCED TO UPPER HESSENBERG FORM.           !ELZHCL12
!                         ON OUTPUT, A CONTAINS THE COMPLEX UPPER       !ELZHCL13
!                           HESSENBERG MATRIX, LAZ.                     !ELZHCL14
!                IA     - THE ROW DIMENSION OF A IN THE CALLING         !ELZHCL15
!                           PROGRAM. IA MUST NOT BE LESS THAN N.(INPUT) !ELZHCL16
!                B      - ON INPUT, B IS THE COMPLEX MATRIX TO BE       !ELZHCL17
!                           REDUCED TO UPPER TRIANGULAR FORM.           !ELZHCL18
!                         ON OUTPUT, B CONTAINS THE UPPER TRIANGULAR    !ELZHCL19
!                           MATRIX, LBZ.                                !ELZHCL20
!                IB     - THE ROW DIMENSION OF B IN THE CALLING         !ELZHCL21
!                           PROGRAM. IB MUST NOT BE LESS THAN N.(INPUT) !ELZHCL22
!                IJOB   - INPUT OPTION PARAMETER. WHEN IJOB = 1, THE    !ELZHCL23
!                           TRANSFORMATIONS NEEDED TO COMPUTE THE       !ELZHCL24
!                           EIGENVECTORS OF THE PROBLEM ARE RETURNED    !ELZHCL25
!                           IN Z. OTHERWISE, Z IS NOT USED.             !ELZHCL26
!                Z      - ON OUTPUT, WHEN IJOB = 1, THE COMPLEX MATRIX  !ELZHCL27
!                           Z CONTAINS THE TRANSFORMATION FROM THE      !ELZHCL28
!                           HESSENBERG REDUCTION.  LAZ IS UPPER HESSEN- !ELZHCL29
!                           BERG AND LBZ IS UPPER TRIANGULAR.  L IS NOT !ELZHCL30
!                           EXPLICITLY FORMED. WHEN IJOB IS NOT EQUAL   !ELZHCL31
!                           TO 1, Z IS NOT USED.                        !ELZHCL32
!                IZ     - THE ROW DIMENSION OF Z IN THE CALLING         !ELZHCL33
!                           PROGRAM. IF IJOB IS NOT EQUAL TO 1, THE     !ELZHCL34
!                           EIGENVECTORS ARE NOT COMPUTED AND Z IS      !ELZHCL35
!                           NOT USED. (INPUT)                           !ELZHCL36
!   PRECISION           - SINGLE                                        !ELZHCL37
!   LANGUAGE            - FORTRAN                                       !ELZHCL38
!-----------------------------------------------------------------------!ELZHCL39
!   LATEST REVISION     - APRIL 24, 1978                                !ELZHCL40
!                                                                       !ELZHCL41
      COMPLEX            Y,W,ZC,A(IA,N),B(IB,N),Z(IZ,*),ONE,ZERO        !ELZHCL42
      REAL               C,D,ZRO                                        !ELZHCL43
      DATA               ZERO/(0.0,0.0)/,ONE/(1.0,0.0)/                 !ELZHCL44
      DATA               ZRO/0.0/                                       !ELZHCL45
      IF(N.EQ.1) GO TO 135                                              !ELZHCL46
      NM1 = N-1                                                         !ELZHCL47
      NP1 = N+1                                                         !ELZHCL48
!                                  REDUCE B TO TRIANGULAR FORM USING    !ELZHCL49
!                                  ELEMENTARY TRANSFORMATIONS           !ELZHCL50
      DO 40 I=1,NM1                                                     !ELZHCL51
         D = ZRO                                                        !ELZHCL52
         IP1 = I+1                                                      !ELZHCL53
         DO 5 K=IP1,N                                                   !ELZHCL54
            Y = B(K,I)                                                  !ELZHCL55
            C = ABS(REAL(Y))+ABS(IMAG(Y))                              !ELZHCL56
            IF (C.LE.D) GO TO 5                                         !ELZHCL57
            D = C                                                       !ELZHCL58
            II = K                                                      !ELZHCL59
    5    CONTINUE                                                       !ELZHCL60
         IF (D.EQ.ZRO) GO TO 40                                         !ELZHCL61
         Y = B(I,I)                                                     !ELZHCL62
         IF (D.LE.ABS(REAL(Y))+ABS(IMAG(Y))) GO TO 20                  !ELZHCL63
!                                  MUST INTERCHANGE                     !ELZHCL64
         DO 10 J=1,N                                                    !ELZHCL65
            Y = A(I,J)                                                  !ELZHCL66
            A(I,J) = A(II,J)                                            !ELZHCL67
            A(II,J) = Y                                                 !ELZHCL68
   10    CONTINUE                                                       !ELZHCL69
         DO 15 J=I,N                                                    !ELZHCL70
            Y = B(I,J)                                                  !ELZHCL71
            B(I,J) = B(II,J)                                            !ELZHCL72
            B(II,J) = Y                                                 !ELZHCL73
   15    CONTINUE                                                       !ELZHCL74
   20    DO 35 J=IP1,N                                                  !ELZHCL75
            Y = B(J,I)/B(I,I)                                           !ELZHCL76
            IF (REAL(Y).EQ.ZRO.AND.IMAG(Y).EQ.ZRO) GO TO 35            !ELZHCL77
            DO 25 K=1,N                                                 !ELZHCL78
               A(J,K) = A(J,K)-Y*A(I,K)                                 !ELZHCL79
   25       CONTINUE                                                    !ELZHCL80
            DO 30 K=IP1,N                                               !ELZHCL81
               B(J,K) = B(J,K)-Y*B(I,K)                                 !ELZHCL82
   30       CONTINUE                                                    !ELZHCL83
   35    CONTINUE                                                       !ELZHCL84
         B(IP1,I) = ZERO                                                !ELZHCL85
   40 CONTINUE                                                          !ELZHCL86
!                                  INITIALIZE Z                         !ELZHCL87
      IF (IJOB.NE.1) GO TO 55                                           !ELZHCL88
      DO 50 I=1,N                                                       !ELZHCL89
         DO 45 J=1,N                                                    !ELZHCL90
            Z(I,J) = ZERO                                               !ELZHCL91
   45    CONTINUE                                                       !ELZHCL92
         Z(I,I) = ONE                                                   !ELZHCL93
   50 CONTINUE                                                          !ELZHCL94
!                                  REDUCE A TO UPPER HESSENBERG FORM    !ELZHCL95
   55 NM2 = N-2                                                         !ELZHCL96
      IF (NM2.LT.1) GO TO 135                                           !ELZHCL97
      DO 130 J=1,NM2                                                    !ELZHCL98
         JM2 = NM1-J                                                    !ELZHCL99
         JP1 = J+1                                                      !ELZHC100
         DO 125 II=1,JM2                                                !ELZHC101
            I = NP1-II                                                  !ELZHC102
            IM1 = I-1                                                   !ELZHC103
            IMJ = I-J                                                   !ELZHC104
            W = A(I,J)                                                  !ELZHC105
            ZC = A(IM1,J)                                               !ELZHC106
            IF (ABS(REAL(W))+ABS(IMAG(W)).LE.ABS(REAL(ZC))+ABS(IMAG(ZC))) GO TO 70  
!                                  MUST INTERCHANGE ROWS                !ELZHC109
            DO 60 K=J,N                                                 !ELZHC110
               Y = A(I,K)                                               !ELZHC111
               A(I,K) = A(IM1,K)                                        !ELZHC112
               A(IM1,K) = Y                                             !ELZHC113
   60       CONTINUE                                                    !ELZHC114
            DO 65 K=IM1,N                                               !ELZHC115
               Y = B(I,K)                                               !ELZHC116
               B(I,K) = B(IM1,K)                                        !ELZHC117
               B(IM1,K) = Y                                             !ELZHC118
   65       CONTINUE                                                    !ELZHC119
   70       ZC = A(I,J)                                                 !ELZHC120
            IF (REAL(ZC).EQ.ZRO.AND.IMAG(ZC).EQ.ZRO) GO TO 85          !ELZHC121
            Y = ZC/A(IM1,J)                                             !ELZHC122
            DO 75 K=JP1,N                                               !ELZHC123
               A(I,K) = A(I,K)-Y*A(IM1,K)                               !ELZHC124
   75       CONTINUE                                                    !ELZHC125
            DO 80 K=IM1,N                                               !ELZHC126
               B(I,K) = B(I,K)-Y*B(IM1,K)                               !ELZHC127
   80       CONTINUE                                                    !ELZHC128
!                                  TRANSFORMATION FROM THE RIGHT        !ELZHC129
   85       W = B(I,IM1)                                                !ELZHC130
            ZC = B(I,I)                                                 !ELZHC131
            IF (ABS(REAL(W))+ABS(IMAG(W)).LE.ABS(REAL(ZC))+ABS(IMAG(ZC))) GO TO 105    
!                                  MUST INTERCHANGE COLUMNS             !ELZHC134
            DO 90 K=1,I                                                 !ELZHC135
               Y = B(K,I)                                               !ELZHC136
               B(K,I) = B(K,IM1)                                        !ELZHC137
               B(K,IM1) = Y                                             !ELZHC138
   90       CONTINUE                                                    !ELZHC139
            DO 95 K=1,N                                                 !ELZHC140
               Y = A(K,I)                                               !ELZHC141
               A(K,I) = A(K,IM1)                                        !ELZHC142
               A(K,IM1) = Y                                             !ELZHC143
   95       CONTINUE                                                    !ELZHC144
            IF (IJOB.NE.1) GO TO 105                                    !ELZHC145
            DO 100 K=IMJ,N                                              !ELZHC146
               Y = Z(K,I)                                               !ELZHC147
               Z(K,I) = Z(K,IM1)                                        !ELZHC148
               Z(K,IM1) = Y                                             !ELZHC149
  100       CONTINUE                                                    !ELZHC150
  105       ZC = B(I,IM1)                                               !ELZHC151
            IF (REAL(ZC).EQ.ZRO.AND.IMAG(ZC).EQ.ZRO) GO TO 125         !ELZHC152
            Y = ZC/B(I,I)                                               !ELZHC153
            DO 110 K=1,IM1                                              !ELZHC154
               B(K,IM1) = B(K,IM1)-Y*B(K,I)                             !ELZHC155
  110       CONTINUE                                                    !ELZHC156
            B(I,IM1) = ZERO                                             !ELZHC157
            DO 115 K=1,N                                                !ELZHC158
               A(K,IM1) = A(K,IM1)-Y*A(K,I)                             !ELZHC159
  115       CONTINUE                                                    !ELZHC160
            IF (IJOB.NE.1) GO TO 125                                    !ELZHC161
            DO 120 K=IMJ,N                                              !ELZHC162
               Z(K,IM1) = Z(K,IM1)-Y*Z(K,I)                             !ELZHC163
  120       CONTINUE                                                    !ELZHC164
  125    CONTINUE                                                       !ELZHC165
         A(JP1+1,J) = ZERO                                              !ELZHC166
  130 CONTINUE                                                          !ELZHC167
  135 RETURN                                                            !ELZHC168
      END                                                               !ELZHC169
!**********************************************************************
SUBROUTINE QXZ166(N,A,IA,B,IB,Z,IZ,IJOB,EIGA,EIGB,IER)
!**********************************************************************
!                                                                       !ELZVCL 3
!                                                                       !ELZVCL 4
!   FUNCTION            - CALCULATE THE EIGENVALUES AND, OPTIONALLY,    !ELZVCL 5
!                           EIGENVECTORS OF THE SYSTEM A*Z=LAMBDA*B*Z   !ELZVCL 6
!                           WHERE COMPLEX MATRIX A IS UPPER HESSENBERG  !ELZVCL 7
!                           AND COMPLEX MATRIX B IS UPPER TRIANGULAR.   !ELZVCL 8
!   USAGE               - CALL QXZ166(N,A,IA,B,IB,Z,IZ,IJOB,EIGA,       !ELZVCL 9
!                           EIGB,INFER,IER)                             !ELZVCL10
!   PARAMETERS   N      - ORDER OF A, B, AND Z. LENGTH OF EIGA          !ELZVCL11
!                           AND EIGB. (INPUT)                           !ELZVCL12
!                A      - ON INPUT, A IS AN N BY N UPPER HESSENBERG     !ELZVCL13
!                           COMPLEX MATRIX. ON OUTPUT, A IS DESTROYED.  !ELZVCL14
!                IA     - ROW DIMENSION OF A IN THE CALLING PROGRAM.    !ELZVCL15
!                           (INPUT). IA MUST NOT BE LESS THAN N.        !ELZVCL16
!                B      - ON INPUT, B IS AN N BY N UPPER TRIANGULAR     !ELZVCL17
!                           COMPLEX MATRIX. ON OUTPUT, B IS DESTROYED.  !ELZVCL18
!                IB     - ROW DIMENSION OF B IN THE CALLING PROGRAM     !ELZVCL19
!                           (INPUT). IB MUST NOT BE LESS THAN N.        !ELZVCL20
!                Z      - ON INPUT, COMPLEX MATRIX Z SHOULD CONTAIN     !ELZVCL21
!                           EITHER THE IDENTITY MATRIX OR THE TRANS-    !ELZVCL22
!                           FORMATION MATRIX PRODUCED BY QXZ165.        !ELZVCL23
!                         ON OUTPUT, Z WILL CONTAIN THE EIGENVECTORS    !ELZVCL24
!                           OF THE SYSTEM GIVEN BY INPUT A AND B OR     !ELZVCL25
!                           THE EIGENVECTORS OF THE SYSTEM WHICH        !ELZVCL26
!                           WAS REDUCED BY QXZ165.                      !ELZVCL27
!                IZ     - ROW DIMENSION OF Z IN THE CALLING PROGRAM.    !ELZVCL28
!                           IF IJOB IS NOT EQUAL TO 1, THE EIGENVECTORS !ELZVCL29
!                           ARE NOT COMPUTED AND Z IS NOT USED. (INPUT) !ELZVCL30
!                IJOB   - INPUT OPTION PARAMETER. WHEN IJOB = 1, THE    !ELZVCL31
!                           TRANSFORMATIONS NEEDED TO COMPUTE THE       !ELZVCL32
!                           EIGENVECTORS OF THE PROBLEM ARE RETURNED    !ELZVCL33
!                           IN Z. OTHERWISE, Z IS NOT USED.             !ELZVCL34
!                           WHEN IJOB = 2, QXZ166 MERELY RENORMALIZES   !ELZVCL35
!                           EIGENVECTORS TO HAVE LARGEST COMPONENT      !ELZVCL36
!                           = 1.0 TO GET CONSISTENT PERFORMANCE INDEX.  !ELZVCL37
!                EIGA   - OUTPUT COMPLEX VECTORS OF LENGTH N.  EIGA     !ELZVCL38
!                EIGB       CONTAINS THE DIAGONAL OF THE TRIANGULARIZED !ELZVCL39
!                           A, AND EIGB CONTAINS THE DIAGONAL OF THE    !ELZVCL40
!                           TRIANGULARIZED B. SEE PROGRAMMING NOTES.    !ELZVCL41
!                         THE J-TH EIGENVALUE OF THE SYSTEM A*Z=        !ELZVCL42
!                           LAMBDA*B*Z IS EIGA(J)/EIGB(J).              !ELZVCL43
!                IER    - OUTPUT INTEGER ERROR CODE.                    !ELZVCL44
!                           IER = J, QXZ166 FAILED TO CONVERGE ON EIGEN-!ELZVCL45
!                             VALUE J WITHIN 30 ITERATIONS.  EIGEN-     !ELZVCL46
!                             VALUES J+1,J+2,...,N HAVE BEEN COMPUTED   !ELZVCL47
!                             CORRECTLY.  EIGENVALUES 1,2,...,N MAY BE  !ELZVCL48
!                             INACCURATE.  EIGENVECTORS ARE NOT CORRECT.!ELZVCL49
!                                                                       !ELZVCL50
!   PRECISION           - SINGLE                                        !ELZVCL51
!   LANGUAGE            - FORTRAN                                       !ELZVCL52
!-----------------------------------------------------------------------!ELZVCL53
!   LATEST REVISION     - APRIL 24, 1978                                !ELZVCL54
COMPLEX A(IA,N),B(IB,N),EIGA(N),EIGB(N),Z(IZ,*),ANNM1,ALFM
COMPLEX BETM,D,SL,DEN,NUM,ANM1M1,S,W,Y,ZC,ZERO                                  
REAL EPSA,EPSB,SS,R,ANORM,BNORM,ANI,BNI,C,D0,D1,D2,E0,E1,ONE,PT8,TWO,ZRO                 
DATA PT8/.8/,TWO/2.0/,ZRO/0.0/,ONE/1.0/,ZERO/(0.0,0.0)/                                
      IER = 0                                                           !ELZVCL62
      IF(IJOB .GT. 1) GO TO 202                                         !ELZVCL63
      NN = N                                                            !ELZVCL64
!                                  COMPUTE THE MACHINE PRECISION TIMES  !ELZVCL65
!                                    THE NORM OF A AND B                !ELZVCL66
      ANORM = ZRO                                                       !ELZVCL67
      BNORM = ZRO                                                       !ELZVCL68
      IM1 = 1                                                           !ELZVCL69
      DO 15 I=1,N                                                       !ELZVCL70
         ANI = ZRO                                                      !ELZVCL71
         IF (I.EQ.1) GO TO 5                                            !ELZVCL72
         Y = A(I,IM1)                                                   !ELZVCL73
         ANI = ANI+ABS(REAL(Y))+ABS(IMAG(Y))                           !ELZVCL74
    5    BNI = ZRO                                                      !ELZVCL75
         DO 10 J=I,N                                                    !ELZVCL76
            ANI = ANI+ABS(REAL(A(I,J)))+ABS(IMAG(A(I,J)))              !ELZVCL77
            BNI = BNI+ABS(REAL(B(I,J)))+ABS(IMAG(B(I,J)))              !ELZVCL78
   10    CONTINUE                                                       !ELZVCL79
         IF (ANI.GT.ANORM) ANORM = ANI                                  !ELZVCL80
         IF (BNI.GT.BNORM) BNORM = BNI                                  !ELZVCL81
         IM1 = I                                                        !ELZVCL82
   15 CONTINUE                                                          !ELZVCL83
      IF (ANORM.EQ.ZRO) ANORM = ONE                                     !ELZVCL84
      IF (BNORM.EQ.ZRO) BNORM = ONE                                     !ELZVCL85
      EPSB = BNORM                                                      !ELZVCL86
      EPSA = ANORM                                                      !ELZVCL87
   20 EPSA = EPSA/TWO                                                   !ELZVCL88
      EPSB = EPSB/TWO                                                   !ELZVCL89
      C = ANORM+EPSA                                                    !ELZVCL90
      IF (C.GT.ANORM) GO TO 20                                          !ELZVCL91
      IF (N.LE.1) GO TO 160                                             !ELZVCL92
   25 ITS = 0                                                           !ELZVCL93
      NM1 = NN-1                                                        !ELZVCL94
!                                  CHECK FOR NEGLIGIBLE                 !ELZVCL95
!                                    SUBDIAGONAL ELEMENTS               !ELZVCL96
   30 D2 = ABS(REAL(A(NN,NN)))+ABS(IMAG(A(NN,NN)))                     !ELZVCL97
      NNP2 = NN+2                                                       !ELZVCL98
      DO 35 LB=2,NN                                                     !ELZVCL99
         L = NNP2-LB                                                    !ELZVC100
         LM1 = L-1                                                      !ELZVC101
         SS = D2                                                        !ELZVC102
         Y = A(LM1,LM1)                                                 !ELZVC103
         D2 = ABS(REAL(Y))+ABS(IMAG(Y))                                !ELZVC104
         SS = SS+D2                                                     !ELZVC105
         Y = A(L,LM1)                                                   !ELZVC106
         R = SS+ABS(REAL(Y))+ABS(IMAG(Y))                              !ELZVC107
         IF (R.EQ.SS) GO TO 40                                          !ELZVC108
   35 CONTINUE                                                          !ELZVC109
      L = 1                                                             !ELZVC110
   40 IF (L.EQ.NN) GO TO 160                                            !ELZVC111
      IF (ITS.LT.30) GO TO 45                                           !ELZVC112
      IF(IER .EQ. 0) IER = NN                                           !ELZVC113
      IF (ABS(REAL(A(NN,NM1)))+ABS(IMAG(A(NN,NM1))).GT.PT8*ABS(REAL(ANNM1))+ABS(IMAG(ANNM1))) GO TO 9000                    !ELZVC115
   45 IF (ITS.EQ.10.OR.ITS.EQ.20) GO TO 55                              !ELZVC116
!                                  COMPUTE SHIFT AS EIGENVALUE          !ELZVC117
!                                    OF LOWER 2 BY 2                    !ELZVC118
      ANNM1 = A(NN,NM1)                                                 !ELZVC119
      ANM1M1 = A(NM1,NM1)                                               !ELZVC120
      S = A(NN,NN)*B(NM1,NM1)-ANNM1*B(NM1,NN)                           !ELZVC121
      W = ANNM1*B(NN,NN)*(A(NM1,NN)*B(NM1,NM1)-B(NM1,NN)*ANM1M1)        !ELZVC122
      Y = (ANM1M1*B(NN,NN)-S)/TWO                                       !ELZVC123
      ZC = sqrt(Y*Y+W)                                                  !ELZVC124
      IF (REAL(ZC).EQ.ZRO.AND.IMAG(ZC).EQ.ZRO) GO TO 50                !ELZVC125
      D0 = REAL(Y/ZC)                                                   !ELZVC126
      IF (D0.LT.ZRO) ZC = -ZC                                           !ELZVC127
   50 DEN = (Y+ZC)*B(NM1,NM1)*B(NN,NN)                                  !ELZVC128
      IF (REAL(DEN).EQ.ZRO.AND.IMAG(DEN).EQ.ZRO) DEN = CMPLX(EPSA,ZRO)                                                   !ELZVC130
      NUM = (Y+ZC)*S-W                                                  !ELZVC131
      GO TO 60                                                          !ELZVC132
!                                  AD-HOC SHIFT                         !ELZVC133
   55 Y = A(NM1,NN-2)                                                   !ELZVC134
      NUM = CMPLX(ABS(REAL(ANNM1))+ABS(IMAG(ANNM1)),ABS(REAL(Y))+ABS(IMAG(Y)))                                                    !ELZVC136
      DEN = CMPLX(ONE,ZRO)                                              !ELZVC137
!                                  CHECK FOR 2 CONSECUTIVE SMALL        !ELZVC138
!                                    SUBDIAGONAL ELEMENTS               !ELZVC139
   60 IF (NN.EQ.L+1) GO TO 70                                           !ELZVC140
      D2 = ABS(REAL(A(NM1,NM1)))+ABS(IMAG(A(NM1,NM1)))                 !ELZVC141
      E1 = ABS(REAL(ANNM1))+ABS(IMAG(ANNM1))                           !ELZVC142
      D1 = ABS(REAL(A(NN,NN)))+ABS(IMAG(A(NN,NN)))                     !ELZVC143
      NL = NN-(L+1)                                                     !ELZVC144
      DO 65 MB=1,NL                                                     !ELZVC145
         M = NN-MB                                                      !ELZVC146
         MN1 = M-1                                                      !ELZVC147
         E0 = E1                                                        !ELZVC148
         Y = A(M,MN1)                                                   !ELZVC149
         E1 = ABS(REAL(Y))+ABS(IMAG(Y))                                !ELZVC150
         D0 = D1                                                        !ELZVC151
         D1 = D2                                                        !ELZVC152
         Y = A(MN1,MN1)                                                 !ELZVC153
         D2 = ABS(REAL(Y))+ABS(IMAG(Y))                                !ELZVC154
         Y = A(M,M)*DEN-B(M,M)*NUM                                      !ELZVC155
         D0 = (D0+D1+D2)*(ABS(REAL(Y))+ABS(IMAG(Y)))                   !ELZVC156
         E0 = E0*E1*(ABS(REAL(DEN))+ABS(IMAG(DEN)))+D0                 !ELZVC157
         IF (E0.EQ.D0) GO TO 75                                         !ELZVC158
   65 CONTINUE                                                          !ELZVC159
   70 M = L                                                             !ELZVC160
   75 CONTINUE                                                          !ELZVC161
      ITS = ITS+1                                                       !ELZVC162
      W = A(M,M)*DEN-B(M,M)*NUM                                         !ELZVC163
      ZC = A(M+1,M)*DEN                                                 !ELZVC164
      D1 = ABS(REAL(ZC))+ABS(IMAG(ZC))                                 !ELZVC165
      D2 = ABS(REAL(W))+ABS(IMAG(W))                                   !ELZVC166
!                                  FIND L AND M AND SET A=LAM AND B=LBM !ELZVC167
      LOR1 = L                                                          !ELZVC168
      NNORN = NN                                                        !ELZVC169
      IF (IJOB.NE.1) GO TO 80                                           !ELZVC170
      LOR1 = 1                                                          !ELZVC171
      NNORN = N                                                         !ELZVC172
   80 IM1 = M                                                           !ELZVC173
      DO 155 I=M,NM1                                                    !ELZVC174
         J = I+1                                                        !ELZVC175
         JP1 = J+1                                                      !ELZVC176
!                                  FIND ROW TRANSFORMATIONS TO RESTORE  !ELZVC177
!                                    A TO UPPER HESSENBERG FORM.        !ELZVC178
!                                    APPLY TRANSFORMATIONS TO A AND B   !ELZVC179
         IF (I.EQ.M) GO TO 85                                           !ELZVC180
         W = A(I,IM1)                                                   !ELZVC181
         ZC = A(J,IM1)                                                  !ELZVC182
         D1 = ABS(REAL(ZC))+ABS(IMAG(ZC))                              !ELZVC183
         D2 = ABS(REAL(W))+ABS(IMAG(W))                                !ELZVC184
         IF (D1.EQ.ZRO) GO TO 30                                        !ELZVC185
   85    IF (D2.GT.D1) GO TO 95                                         !ELZVC186
!                                  MUST INTERCHANGE ROWS                !ELZVC187
         DO 90 K=I,NNORN                                                !ELZVC188
            Y = A(I,K)                                                  !ELZVC189
            A(I,K) = A(J,K)                                             !ELZVC190
            A(J,K) = Y                                                  !ELZVC191
            Y = B(I,K)                                                  !ELZVC192
            B(I,K) = B(J,K)                                             !ELZVC193
            B(J,K) = Y                                                  !ELZVC194
   90    CONTINUE                                                       !ELZVC195
         IF (I.GT.M) A(I,IM1) = A(J,IM1)                                !ELZVC196
         IF (D2.EQ.ZRO) GO TO 110                                       !ELZVC197
!                                  THE SCALING OF W AND Z IS DESIGNED   !ELZVC198
!                                    TO AVOID A DIVISION BY ZERO        !ELZVC199
!                                    WHEN THE DENOMINATOR IS SMALL      !ELZVC200
         Y = CMPLX(REAL(W)/D1,IMAG(W)/D1)/CMPLX(REAL(ZC)/D1,IMAG(ZC)/D1)                                                  !ELZVC202
         GO TO 100                                                      !ELZVC203
   95    Y = CMPLX(REAL(ZC)/D2,IMAG(ZC)/D2)/CMPLX(REAL(W)/D2,IMAG(W)/D2)                                                   !ELZVC205
  100    DO 105 K=I,NNORN                                               !ELZVC206
            A(J,K) = A(J,K)-Y*A(I,K)                                    !ELZVC207
            B(J,K) = B(J,K)-Y*B(I,K)                                    !ELZVC208
  105    CONTINUE                                                       !ELZVC209
  110    IF (I.GT.M) A(J,IM1) = ZERO                                    !ELZVC210
!                                  PERFORM TRANSFORMATIONS FROM RIGHT   !ELZVC211
!                                    TO RESTORE B TO TRIANGULAR FORM    !ELZVC212
!                                  APPLY TRANSFORMATIONS TO A           !ELZVC213
         ZC = B(J,I)                                                    !ELZVC214
         IM1 = I                                                        !ELZVC215
         W = B(J,J)                                                     !ELZVC216
         D2 = ABS(REAL(W))+ABS(IMAG(W))                                !ELZVC217
         D1 = ABS(REAL(ZC))+ABS(IMAG(ZC))                              !ELZVC218
         IF (D1.EQ.ZRO) GO TO 30                                        !ELZVC219
         IF (D2.GT.D1) GO TO 135                                        !ELZVC220
!                                  MUST INTERCHANGE COLUMNS             !ELZVC221
         DO 115 K=LOR1,J                                                !ELZVC222
            Y = A(K,J)                                                  !ELZVC223
            A(K,J) = A(K,I)                                             !ELZVC224
            A(K,I) = Y                                                  !ELZVC225
            Y = B(K,J)                                                  !ELZVC226
            B(K,J) = B(K,I)                                             !ELZVC227
            B(K,I) = Y                                                  !ELZVC228
  115    CONTINUE                                                       !ELZVC229
         IF (I.EQ.NM1) GO TO 120                                        !ELZVC230
         Y = A(JP1,J)                                                   !ELZVC231
         A(JP1,J) = A(JP1,I)                                            !ELZVC232
         A(JP1,I) = Y                                                   !ELZVC233
  120    IF (IJOB.NE.1) GO TO 130                                       !ELZVC234
         DO 125 K=1,N                                                   !ELZVC235
            Y = Z(K,J)                                                  !ELZVC236
            Z(K,J) = Z(K,I)                                             !ELZVC237
            Z(K,I) = Y                                                  !ELZVC238
  125    CONTINUE                                                       !ELZVC239
  130    B(J,I) = ZERO                                                  !ELZVC240
         IF (D2.EQ.ZRO) GO TO 155                                       !ELZVC241
         ZC = CMPLX(REAL(W)/D1,IMAG(W)/D1)/CMPLX(REAL(ZC)/D1,IMAG(ZC)/D1)                                                  !ELZVC243
         GO TO 140                                                      !ELZVC244
  135    ZC = CMPLX(REAL(ZC)/D2,IMAG(ZC)/D2)/CMPLX(REAL(W)/D2,IMAG(W)/D2)                                                   !ELZVC246
  140    DO 145 K=LOR1,J                                                !ELZVC247
            A(K,I) = A(K,I)-ZC*A(K,J)                                   !ELZVC248
            B(K,I) = B(K,I)-ZC*B(K,J)                                   !ELZVC249
  145    CONTINUE                                                       !ELZVC250
         B(J,I) = ZERO                                                  !ELZVC251
         IF (I.LT.NM1) A(JP1,I) = A(JP1,I)-ZC*A(JP1,J)                  !ELZVC252
         IF (IJOB.NE.1) GO TO 155                                       !ELZVC253
         DO 150 K=1,N                                                   !ELZVC254
            Z(K,I) = Z(K,I)-ZC*Z(K,J)                                   !ELZVC255
  150    CONTINUE                                                       !ELZVC256
  155 CONTINUE                                                          !ELZVC257
      GO TO 30                                                          !ELZVC258
  160 CONTINUE                                                          !ELZVC259
      EIGA(NN) = A(NN,NN)                                               !ELZVC260
      EIGB(NN) = B(NN,NN)                                               !ELZVC261
      IF (NN.EQ.1) GO TO 165                                            !ELZVC262
      NN = NM1                                                          !ELZVC263
      IF (NN.GT.1) GO TO 25                                             !ELZVC264
      GO TO 160                                                         !ELZVC265
!                                  FIND EIGENVECTORS USING B FOR        !ELZVC266
!                                    INTERMEDIATE STORAGE               !ELZVC267
  165 IF (IJOB.NE.1) GO TO 9000                                         !ELZVC268
      M = N                                                             !ELZVC269
  170 CONTINUE                                                          !ELZVC270
      ALFM = A(M,M)                                                     !ELZVC271
      BETM = B(M,M)                                                     !ELZVC272
      B(M,M) = CMPLX(ONE,ZRO)                                           !ELZVC273
      L = M-1                                                           !ELZVC274
      IF (L.EQ.0) GO TO 185                                             !ELZVC275
  175 CONTINUE                                                          !ELZVC276
      L1 = L+1                                                          !ELZVC277
      SL = ZERO                                                         !ELZVC278
      DO 180 J=L1,M                                                     !ELZVC279
         SL = SL+(BETM*A(L,J)-ALFM*B(L,J))*B(J,M)                       !ELZVC280
  180 CONTINUE                                                          !ELZVC281
      Y = BETM*A(L,L)-ALFM*B(L,L)                                       !ELZVC282
      IF (REAL(Y).EQ.ZRO.AND.IMAG(Y).EQ.ZRO) Y = CMPLX((EPSA+EPSB)/TWO,ZRO)                                        !ELZVC284
      B(L,M) = -SL/Y                                                    !ELZVC285
      L = L-1                                                           !ELZVC286
  185 IF (L.GT.0) GO TO 175                                             !ELZVC287
      M = M-1                                                           !ELZVC288
      IF (M.GT.0) GO TO 170                                             !ELZVC289
!                                  TRANSFORM TO ORIGINAL COORDINATE     !ELZVC290
!                                    SYSTEM                             !ELZVC291
      M = N                                                             !ELZVC292
  190 CONTINUE                                                          !ELZVC293
      DO 200 I=1,N                                                      !ELZVC294
         S = ZERO                                                       !ELZVC295
         DO 195 J=1,M                                                   !ELZVC296
            S = S+Z(I,J)*B(J,M)                                         !ELZVC297
  195    CONTINUE                                                       !ELZVC298
         Z(I,M) = S                                                     !ELZVC299
  200 CONTINUE                                                          !ELZVC300
      M = M-1                                                           !ELZVC301
      IF (M.GT.0) GO TO 190                                             !ELZVC302
!                                  NORMALIZE SO THAT LARGEST            !ELZVC303
!                                    COMPONENT = 1.0                    !ELZVC304
  202 CONTINUE                                                          !ELZVC305
      M = N                                                             !ELZVC306
  205 CONTINUE                                                          !ELZVC307
      SS = ZRO                                                          !ELZVC308
      DO 210 I=1,N                                                      !ELZVC309
         R = ABS(REAL(Z(I,M)))+ABS(IMAG(Z(I,M)))                       !ELZVC310
         IF (R.LT.SS) GO TO 210                                         !ELZVC311
         SS = R                                                         !ELZVC312
         D = Z(I,M)                                                     !ELZVC313
  210 CONTINUE                                                          !ELZVC314
      IF (SS.EQ.ZRO) GO TO 220                                          !ELZVC315
      DO 215 I=1,N                                                      !ELZVC316
         Z(I,M) = Z(I,M)/D                                              !ELZVC317
  215 CONTINUE                                                          !ELZVC318
  220 M = M-1                                                           !ELZVC319
      IF (M.GT.0) GO TO 205                                             !ELZVC320
      GO TO 9005                                                        !ELZVC321
 9000 CONTINUE                                                          !ELZVC322
 9005 RETURN                                                            !ELZVC323
      END                                                               !ELZVC324

!**********************************************************************
SUBROUTINE CXGCOIT (NiMAX,N,A,NRHS,B,CONLIM,IPIVOT,IFAC,UL,CON,WK,IERR)     
!**********************************************************************
!                                                                       !CXGCOIT4
!               F4.13   IS ROUTINE ID.                                  !CXGCOIT5
!***********************************************************************!CXGCOIT6
!*                                                                     *!CXGCOIT7
!* LANGUAGE:                                                           *!CXGCOIT8
!*         FORTRAN                                                     *!CXGCOIT9
!*                                                                     *!CXGCOI10
!* PURPOSE:                                                            *!CXGCOI11
!*         TO PERFORM THE LU DECOMPOSITION OF A COMPLEX MATRIX A AND   *!CXGCOI12
!*         SIMULTANEOUSLY COMPUTE AN ESTIMATE OF THE CONDITION NUMBER  *!CXGCOI13
!*         OF THE MATRIX. OPTIONALLY,A SYSTEM OF EQUATIONS AX=B MAY BE *!CXGCOI14
!*         SOLVED FOR ANY NUMBER OF COMPLEX VALUED RIGHT-HAND SIDE     *!CXGCOI15
!*         VECTORS B. AS A SUB-OPTION TO THIS OPTION, ITERATIVE        *!CXGCOI16
!*         IMPROVEMENT OF THE SOLUTION WILL BE PERFORMED IF THE ESTI-  *!CXGCOI17
!*         MATED CONDITION NUMBER IS ABOVE A USER SPECIFIED VALUE.     *!CXGCOI18
!*                                                                     *!CXGCOI19
!* USE:                                                                *!CXGCOI20
!*    CALL !CXGCOIT (NiMAX,N,A,NRHS,B,CONLIM,IPIVOT,                     *!CXGCOI21
!*                  IFAC,UL,CON,WK,IERR)                               *!CXGCOI22
!*                                                                     *!CXGCOI23
!* PARAMETERS:                                                         *!CXGCOI24
!*                                                                     *!CXGCOI25
!*       NiMAX  AN INPUT INTEGER SPECIFYING THE FIRST DIMENSION OF THE  *!CXGCOI26
!*             A,B, AND UL ARRAYS IN CALLING PROGRAM. NiMAX >= N        *!CXGCOI27
!*                                                                     *!CXGCOI28
!*          N  AN INPUT INTEGER SPECIFYING THE ORDER OF THE A MATRIX.  *!CXGCOI29
!*                                                                     *!CXGCOI30
!*          A  AN INPUT TWO-DIMENSIONAL COMPLEX ARRAY WITH FIRST       *!CXGCOI31
!*             DIMENSION EQUAL TO NiMAX AND SECOND DIMENSION AT LEAST   *!CXGCOI32
!*             N. A MUST BE TYPED COMPLEX IN THE CALLING PROGRAM. ON   *!CXGCOI33
!*             ENTRY, A MUST CONTAIN THE COEFFICIENT MATRIX TO BE      *!CXGCOI34
!*             DECOMPOSED. A IS NOT ALTERED BY GLMCOIT.                *!CXGCOI35
!*                                                                     *!CXGCOI36
!*       NRHS  AN INPUT INTEGER SPECIFYING THE NUMBER OF RIGHT HAND    *!CXGCOI37
!*             SIDE VECTORS OF THE SYSTEM OF EQUATIONS TO BE SOLVED.   *!CXGCOI38
!*                                                                     *!CXGCOI39
!*          B  AN INPUT/OUTPUT TWO-DIMENSIONAL COMPLEX ARRAY WITH FIRST*!CXGCOI40
!*             DIMENSION EQUAL TO NiMAX AND SECOND DIMENSION AT LEAST   *!CXGCOI41
!*             NRHS. B MUST BE TYPED COMPLEX IN THE CALLING PROGRAM.   *!CXGCOI42
!*             ON ENTRY, THE COLUMNS OF B MUST CONTAIN THE NRHS RIGHT  *!CXGCOI43
!*             HAND SIDE VECTORS OF THE SYSTEM TO BE SOLVED. ON RETURN,*!CXGCOI44
!*             THE COLUMNS OF B WILL CONTAIN THE CORRESPONDING SOLUTION*!CXGCOI45
!*             VECTORS. IF NO RIGHT HAND SIDES ARE TO BE SOLVED FOR,   *!CXGCOI46
!*             (IFAC=-1), B CAN BE A DUMMY PARAMETER AND NEED NOT BE   *!CXGCOI47
!*             DIMENSIONED.                                            *!CXGCOI48
!*                                                                     *!CXGCOI49
!*     CONLIM  AN INPUT REAL PARAMETER SPECIFYING A iMAXIMUM CONDITION  *!CXGCOI50
!*             NUMBER FOR WHICH THE SOLUTION WILL BE ACCEPTED WITHOUT  *!CXGCOI51
!*             IERATIVE IMPROVEMENT. IF NO RIGHT HAND SIDES ARE TO BE  *!CXGCOI52
!*             SOLVED FOR, CONLIM MAY BE A  DUMMY PARAMETER.           *!CXGCOI53
!*                                                                     *!CXGCOI54
!*     IPIVOT  AN INPUT/OUTPUT ONE-DIMENSIONAL INTEGER ARRAY           *!CXGCOI55
!*             DIMENSIONED AT LEAST N.                                 *!CXGCOI56
!*                                                                     *!CXGCOI57
!*                INPUT      IF THE LU DECOMPOSITION OF A IS TO BE     *!CXGCOI58
!*                           INPUT FROM A PREVIOUS CALL TO GLMCOIT,    *!CXGCOI59
!*                           (IFAC=1),THE FIRST N ELEMENTS OF IPIVOT   *!CXGCOI60
!*                           MUST CONTAIN THE PIVOTED STRATEGY USED IN *!CXGCOI61
!*                           PERFORMING THE LU DECOMPOSITION. IF THE   *!CXGCOI62
!*                           LU DECOMPOSITION IS NOT TO BE INPUT,      *!CXGCOI63
!*                           (IFAC NOT EQUAL TO 1), NO INPUT TO        *!CXGCOI64
!*                           IPIVOT IS NEEDED.                         *!CXGCOI65
!*                                                                     *!CXGCOI66
!*                OUTPUT     THE FIRST N ELEMENTS OF IPIVOT WILL       *!CXGCOI67
!*                           CONTAIN THE PIVOTED STRATEGY USED IN PER- *!CXGCOI68
!*                           FORMING THE LU DECOMPOSITION FOR ALL      *!CXGCOI69
!*                           VALUES OF IFAC.                           *!CXGCOI70
!*                                                                     *!CXGCOI71
!*       IFAC  AN INPUT/OUTPUT INTEGER SPECIFYING THE COMPUTATION TO   *!CXGCOI72
!*             BE PERFORMED.                                           *!CXGCOI73
!*                                                                     *!CXGCOI74
!*                INPUT = -1 COMPUTE THE LU DECOMPOSITION OF A, AND THE*!CXGCOI75
!*                           ESTIMATED CONDITION NUMBER ONLY.          *!CXGCOI76
!*                                                                     *!CXGCOI77
!*                      =  0 COMPUTE THE LU DECOMPOSITION OF A, THE    *!CXGCOI78
!*                           ESTIMATED CONDITION NUMBER, AND SOLVE THE *!CXGCOI79
!*                           SYSTEM OF EQUATIONS.                      *!CXGCOI80
!*                                                                     *!CXGCOI81
!*                      =  1 SOLVE  SYSTEM OF EQUATIONS,ASSUMING THE   *!CXGCOI82
!*                           INFORMATION ON THE LU DECOMPOSITION OF A  *!CXGCOI83
!*                           IS STORED IN THE UL AND IPIVOT ARRAYS.    *!CXGCOI84
!*                                                                     *!CXGCOI85
!*               OUTPUT =  1 IN ALL CASES                              *!CXGCOI86
!*                                                                     *!CXGCOI87
!*         UL  AN INPUT/OUTPUT TWO DIMENSIONAL COMPLEX ARRAY WITH FIRST*!CXGCOI88
!*             DIMENSION EQUAL TO NiMAX AND SECOND DIMENSION AT LEAST N.*!CXGCOI89
!*             UL MUST BE TYPED COMPLEX IN THE CALLING PROGRAM.        *!CXGCOI90
!*                                                                     *!CXGCOI91
!*                 INPUT     IF IFAC=1,UL MUST CONTAIN THE LU DECOM-   *!CXGCOI92
!*                           POSITION OF THE A MATRIX. OTHERWISE,UL    *!CXGCOI93
!*                           MAY BE UNDEFINED.                         *!CXGCOI94
!*                                                                     *!CXGCOI95
!*                OUTPUT     UL WILL CONTAIN THE LU DECOMPOSITION OF   *!CXGCOI96
!*                           THE A MATRIX IN ALL CASES.                *!CXGCOI97
!*                                                                     *!CXGCOI98
!*        CON  AN OUTPUT REAL PARAMETER CONTAINING AN ESTIMATE OF THE  *!CXGCOI99
!*             CONDITION NUMBER OF THE A MATRIX. IF THE A MATRIX IS    *!CXGCO100
!*             COMPUTATIONALLY SINGULAR,SO THAT THE CONDITION NUMBER   *!CXGCO101
!*             IS THEORETICALLY INFINITE,CON RETURNS 1.265X10**322     *!CXGCO102
!*             (LARGEST ALLOWABLE REAL NUMBER).                        *!CXGCO103
!*                                                                     *!CXGCO104
!*         WK  A ONE-DIMENSIONAL WORK ARRAY. IF IFAC = -1 OR IF IT IS  *!CXGCO105
!*             KNOWN THAT NO ITERATIVE IMPROVEMENT WILL BE ATTEMPTED,  *!CXGCO106
!*             WK MUST CONTAIN AT LEAST 4N LOCATIONS.                  *!CXGCO107
!*                                                                     *!CXGCO108
!*       IERR  AN OUTPUT INTEGER ERROR PARAMETER                       *!CXGCO109
!*                                                                     *!CXGCO110
!*                       =0  NO ERROR DETECTED                         *!CXGCO111
!*                                                                     *!CXGCO112
!*                       =1  ESTIMATED CONDITION NUMBER IS GREATER THAN*!CXGCO113
!*                           THE RECIPROCAL OF MACHINE PRECISION. NO   *!CXGCO114
!*                           SOLUTION OF EQUATIONS WAS ATTEMPTED AS    *!CXGCO115
!*                           MATRIX IS TOO ILL-CONDITIONED TO OBTAIN  A*!CXGCO116
!*                           MEANINGFUL SOLUTION.                      *!CXGCO117
!*                                                                     *!CXGCO118
!*                       =2  ITERATIVE IMPROVEMENT FAILED TO CONVERGE  *!CXGCO119
!*                           IN TEN ITERATIONS. MATRIX IS COMPUTA-     *!CXGCO120
!*                           TIONALLY SINGULAR.                        *!CXGCO121
!*                                                                     *!CXGCO122
!*    REQUIRED ROUTINES:      QXZ194,QXZ195,QXZ196,QXZ197,QXZ198,QXZ199*!QXZN1946
!*                            QXZ200,QXZ202,QXZ204,QXZ201,QXZ203    *   !QXZN1947
!*                                                                     *!CXGCO125
!*    FORTRAN FUNCTIONS:      ABS,IMAG,MAX,CMPLX,CONJG,REAL         *!CXGCO126
!*                                                                     *!CXGCO127
!*    SOURCE:                 COMPUTER SCIENCES CORPORATION,           *!CXGCO128
!*                            HAMPTON, VA.                             *!CXGCO129
!*                                                                     *!CXGCO130
!*    LANGUAGE:               FORTRAN                                  *!CXGCO131
!*                                                                     *!CXGCO132
!*    DATE RELEASED:          MAY 1,1979                               *!CXGCO133
!*                                                                     *!CXGCO134
!*    LATEST REVISION         MAY 1,1979                               *!CXGCO135
!*                                                                     *!CXGCO136
!***********************************************************************!CXGCO137
      COMPLEX A(NiMAX,*),B(NiMAX,*),UL(NiMAX,*),WK(*)                      
      INTEGER IPIVOT(*)                                                 
      IERR=0                                                            !CXGCO143
      IF(IFAC.EQ.1) GO TO 10                                            !CXGCO144
!              SET UL EQUAL TO A, KEEP A UNMODIFIED                     !CXGCO145
      DO 2 J=1,N                                                        !CXGCO146
      DO 1 I=1,N                                                        !CXGCO147
    1 UL(I,J)=A(I,J)                                                    !CXGCO148
    2 CONTINUE                                                          !CXGCO149
!              CALL SUBROUTINE QXZ196 TO PERFORM THE LU DECOMPOSITION,  !QXZN1948
!              AND PRODUCE THE CONDITION NUMBER                         !CXGCO151
      CALL QXZ196(UL,NiMAX,N,IPIVOT,RCOND,WK)                            
      IF(RCOND.NE.0.0) GO TO 5                                          !CXGCO153
!     CON=O"37767777777777777777"                                       !FTN41839
      con=1.0e30
!     write(6,*)' con = ',con
      GO TO 90                                                          !CXGCO155
    5 CON=1.0/RCOND                                                     !CXGCO156
      T=1.0+RCOND                                                       !CXGCO157
      IF(T.EQ.1.0) GO TO 90                                             !CXGCO158
      IF(IFAC.NE.-1) GO TO 10                                           !CXGCO159
      IFAC=1                                                            !CXGCO160
      RETURN                                                            !CXGCO161
   10 IF(CON.GE.CONLIM) GO TO 20                                        !CXGCO162
!              SOLVE THE SYSTEM OF LINEAR EQUATIONS ONLY, NO ITERATIVE  !CXGCO163
!              IMPROVEMENT OF THE SOLUTIONS REQUIRED                    !CXGCO164
      DO 12 KR=1,NRHS                                                   !CXGCO165
      CALL QXZ198(UL,NiMAX,N,IPIVOT,B(1,KR),0)                           
   12 CONTINUE                                                          !CXGCO167
      IFAC=1                                                            !CXGCO168
      RETURN                                                            !CXGCO169
!              CALL QXZ198 TO SOLVE THE SYSTEM OF LINEAR EQUATIONS      
!              CALL QXZ194 TO PERFORM ITERATIVE IMPROVEMENT OF THE      
!                          SOLUTIONS                                    !CXGCO172
   20 IFAC=1                                                            !CXGCO173
!              SET FIRST HALF OF WK EQUAL TO B, KEEP THE RIGHT HAND     !CXGCO174
!              ARRAY B UNMODIFIED                                       !CXGCO175
      DO 24 KR=1,NRHS                                                   !CXGCO176
      DO 22 I=1,N                                                       !CXGCO177
   22 WK(I)=B(I,KR)                                                     !CXGCO178
      CALL QXZ198(UL,NiMAX,N,IPIVOT,WK,0)                                
      CALL QXZ194(A,UL,NiMAX,N,IPIVOT,WK,B(1,KR),IERR)                   
!              STORE THE IMPROVED SOLUTIONS INTO THE RIGHT HAND ARRAY B !CXGCO181
      DO 23 I=1,N                                                       !CXGCO182
   23 B(I,KR)=WK(I)                                                     !CXGCO183
!              TEST FOR CONVERGANCE                                     !CXGCO184
      IF(IERR.EQ.2) RETURN                                              !CXGCO185
   24 CONTINUE                                                          !CXGCO186
      GO TO 99                                                          !CXGCO187
!              A MATRIX IS SINGULAR                                     !CXGCO188
   90 IERR=1                                                            !CXGCO189
   99 RETURN                                                            !CXGCO190
      END                                                               !CXGCO191

!***********************************************************************!CPROVE 4
      SUBROUTINE QXZ194 (A,UL,NiMAX,N,IPIVOT,WK,B,IERR)                 
!***********************************************************************!CPROVE 4
!*                                                                     *!CPROVE 5
!* PURPOSE:                                                            *!CPROVE 6
!*         TO PERFORM ITERATIVE IMPROVEMENT ON THE SOLUTION TO A       *!CPROVE 7
!*         SYSTEM OF SIMULTANEOUS LINEAR EQUATIONS AS COMPUTED BY      *!CPROVE 8
!*         SUBROUTINE GLMCOIT. ITERATION PROCEEDS UNTIL THE RATIO OF   *!CPROVE 9
!*         THE iMAX-NORM OF THE CORRECTION VECTOR TO THE iMAX-NORM OF    *!CPROVE10
!*         THE SOLUTION VECTOR IS LESS THAN MACHINE PRECISION.         *!CPROVE11
!*         CXPROVE WORKS ON A SINGLE RIGHT HAND SIDE PER CALL.         *!CPROVE12
!*                                                                     *!CPROVE13
!* USE:                                                                *!CPROVE14
!*    CALL QXZ194 (A,UL,NiMAX,N,IPIVOT,WK,B,IERR)                       *!QXZN1958
!*                                                                     *!CPROVE16
!* PARAMETERS:                                                         *!CPROVE17
!*                                                                     *!CPROVE18
!*          A  AN INPUT TWO-DIMENSIONAL COMPLEX ARRAY WITH F&RST       *!CPROVE19
!*             DIMENSION EQUAL TO NiMAX AND SECOND DIMENSION AT LEAST   *!CPROVE20
!*             N. ON ENTRY,A MUST CONTAIN THE COEFFICIENT MATRIX OF    *!CPROVE21
!*             THE SYSTEM OF EQUATIONS.                                *!CPROVE22
!*                                                                     *!CPROVE23
!*         UL  AN INPUT TWO-DIMENSIONAL COMPLEX ARRAY WITH FIRST       *!CPROVE24
!*             DIMENSION EQUAL TO NiMAX AND SECOND DIMENSION AT LEAST   *!CPROVE25
!*             N. ON ENTRY, UL MUST CONTAIN THE LU DECOMPOSITION OF    *!CPROVE26
!*             THE A MATRIX.                                           *!CPROVE27
!*                                                                     *!CPROVE28
!*       NiMAX  AN INPUT INTEGER SPECIFYING THE FIRST DIMENSION OF THE  *!CPROVE29
!*             A AND UL ARRAYS IN THE CALLING PROGRAM.  NiMAX MUST BE   *!CPROVE30
!*             GREATER THAN OR EQUAL TO N.                             *!CPROVE31
!*                                                                     *!CPROVE32
!*          N  AN INPUT INTEGER SPECIFYING THE ORDER OF THE A MATRIX.  *!CPROVE33
!*                                                                     *!CPROVE34
!*     IPIVOT  AN INPUT ONE-DIMENSIONAL INTEGER ARRAY DIMENSIONED AT   *!CPROVE35
!*             LEAST N. ON ENTRY, IPIVOT MUST CONTAIN THE PIVOTING     *!CPROVE36
!*             INFORMATION GENERATED DURING THE LU DECOMPOSITION OF A. *!CPROVE37
!*                                                                     *!CPROVE38
!*         WK  A ONE-DIMENSIONAL COMPLEX WORK ARRAY DIMENSIONED AT     *!CPROVE39
!*             LEAST 2XN. THE FIRST N LOCATIONS OF WK CONTAIN THE LAT- *!CPROVE40
!*             EST APPROXIMATION TO THE SYSTEM OF EQUATIONS, AND THE   *!CPROVE41
!*             NEXT N LOCATIONS CONTAIN THE LATEST RESIDUALS.          *!CPROVE42
!*                                                                     *!CPROVE43
!*          B  AN INPUT ONE-DIMENSIONAL COMPLEX ARRAY DIMENSIONED AT   *!CPROVE44
!*             LEAST N.  B MUST CONTAIN THE RIGHT HAND SIDE OF THE     *!CPROVE45
!*             SYSTEM OF EQUATIONS BEING SOLVED.                       *!CPROVE46
!*                                                                     *!CPROVE47
!*       IERR  AN OUTPUT INTEGER ERROR PARAMETER.                      *!CPROVE48
!*             =0  NO ERROR DETECTED                                   *!CPROVE49
!*                                                                     *!CPROVE50
!*             =2  ITERATIVE IMPROVEMENT FAILED TO CONVERGE IN TEN     *!CPROVE51
!*                 ITERATIONS.  COEFFICIENT MATRIX IS COMPUTATIONALLY  *!CPROVE52
!*                 SINGULAR.                                           *!CPROVE53
!*                                                                     *!CPROVE54
!*    REQUIRED ROUTINES:  QXZ188,QXZ191,QXZ185,QXZ186                  *!QXZN1959
!*                                                                     *!CPROVE56
!*    FORTRAN FUNCTIONS:  ABS,DBLE,MOD,SNGL                            *!CPROVE57
!*                                                                     *!CPROVE58
!*    SOURCE:             COMPUTER SCIENCES CORPORATION,               *!CPROVE59
!*                        HAMPTON, VA.                                 *!CPROVE60
!*                                                                     *!CPROVE61
!*    LANGUAGE:           FORTRAN                                      *!CPROVE62
!*                                                                     *!CPROVE63
!*    DATE RELEASED:      MAY 1,1979                                   *!CPROVE64
!*                                                                     *!CPROVE65
!*    LATEST REVISION:    MAY 1,1979                                   *!CPROVE66
!*                                                                     *!CPROVE67
!***********************************************************************!CPROVE68
!                                                                       !CPROVE69
      COMPLEX A(NiMAX,*),B(NiMAX,*),UL(NiMAX,*),WK(*),SUM,EPS               
      INTEGER IPIVOT(*)                                                  
      DOUBLE PRECISION AR,AI,WKR,WKI,BR,BI,SUMR,SUMI,TEMPR,TEMPI         
      ITiMAX=10                                                           
      INDEX=N+1                                                         !CPROVE74
      EPS=(10.0E-14,10.0E-14)                                           !CPROVE75
!                                                                       !CPROVE76
!              CALL FUNCTION QXZ195 TO COMPUTE THE INFINITY NORM        !QXZN1960
!              OF THE VECTOR WK WHICH CONTAINS THE FIRST SOLUTION X1    !CPROVE78
!                                                                       !CPROVE79
      CXNORM=QXZ195(N,WK(1))                                            
!                                                                       !CPROVE81
!          TEST IF THE SYSTEM OF LINEAR EQUATIONS HAVE ZERO SOLUTIONS   !CPROVE82
!                                                                       !CPROVE83
      IF(CXNORM) 1,90,1                                                 !CPROVE84
!                                                                       !CPROVE85
!              TO PERFORM ITERATIVE IMPROVEMENT OF THE SOLUTIONS        !CPROVE86
!                                                                       !CPROVE87
    1 DO 20 ITER=1,ITiMAX                                                !CPROVE88
      DO 6 I=1,N                                                        !CPROVE89
      SUMR=0.0D0                                                        !CPROVE90
      SUMI=0.0D0                                                        !CPROVE91
!                                                                       !CPROVE92
!              DOUBLE PRECISION COMPUTATION OF A TIMES WK (A*X1)        !CPROVE93
!                                                                       !CPROVE94
      DO 4 K=1,N                                                        !CPROVE95
      AR=DBLE(REAL(A(I,K)))                                             !CPROVE96
      AI=DBLE(IMAG(A(I,K)))                                            !CPROVE97
      WKR=DBLE(REAL(WK(K)))                                             !CPROVE98
      WKI=DBLE(IMAG(WK(K)))                                            !CPROVE99
      TEMPR=AR*WKR-AI*WKI                                               !CPROV100
      TEMPI=AI*WKR+AR*WKI                                               !CPROV101
      SUMR=SUMR+TEMPR                                                   !CPROV102
      SUMI=SUMI+TEMPI                                                   !CPROV103
    4 CONTINUE                                                          !CPROV104
!                                                                       !CPROV105
!              DOUBLE PRECISION COMPUTATION OF THE RESIDUAL. R=B-AX     !CPROV106
!              AND STORE THE SINGLE PRECISION RESIDUAL INTO THE SECOND  !CPROV107
!              HALF OF THE WK-ARRAY.                                    !CPROV108
!                                                                       !CPROV109
      BR=DBLE(REAL(B(I,1)))                                             
      BI=DBLE(IMAG(B(I,1)))                                            
      SUMR=BR-SUMR                                                      !CPROV112
      SUMI=BI-SUMI                                                      !CPROV113
      SUM=CMPLX(SUMR,SUMI)     
      WK(INDEX+I-1)=SUM                                                 !CPROV115
    6 CONTINUE                                                          !CPROV116
!                                                                       !CPROV117
!              SOLVE FOR DX. A*DX=R = UL*DX=R                           !CPROV118
!                                                                       !CPROV119
      CALL QXZ198(UL,NiMAX,N,IPIVOT,WK(INDEX),0)                         
!                                                                       !CPROV121
!              QXZ187 THE ACCURACY OF THE SOLUTIONS                     !QXZN1963
!                                                                       !CPROV123
      DO 8 I=1,N                                                        !CPROV124
    8 WK(I)=WK(I)+WK(INDEX+I-1)                                         !CPROV125
!                                                                       !CPROV126
!              COMPUTE THE INFINITY NORM OF THE DX VECTOR               !CPROV127
!              AND TEST FOR CONVERGANCE                                 !CPROV128
!                                                                       !CPROV129
      CDXNORM=QXZ195(N,WK(INDEX))                                       
      T=CDXNORM/CXNORM                                                  !CPROV131
      IF(T.LE.REAL(EPS) )GO TO 90                                       
   20 CONTINUE                                                          !CPROV133
!                                                                       !CPROV134
!              IF DID NOT CONVERGE AFTER ITiMAX ITERATIONS               !CPROV135
!              SET IERR=2 RETURN                                        !CPROV136
!                                                                       !CPROV137
      IERR=2                                                            !CPROV138
   90 RETURN                                                            !CPROV139
      END                                                               
!**********************************************************************
FUNCTION QXZ195(N,AA)                                             
!***********************************************************************CFNORM 4
!*                                                                     *CFNORM 5
!* PURPOSE:                                                            *CFNORM 6
!*         TO FIND THE iMAX-NORM OF A VECTOR                            *CFNORM 7
!*                                                                     *CFNORM 8
!* USE:                                                                *CFNORM 9
!*         X = QXZ188(N,AA)                                            *!QXZN1968
!*                                                                     *CFNORM11
!* PARAMETERS:                                                         *CFNORM12
!*                                                                     *CFNORM13
!*          N  AN INPUT INTEGER SPECIFYING THE NUMBER OF ELEMENTS IN   *CFNORM14
!*             THE VECTOR.                                             *CFNORM15
!*                                                                     *CFNORM16
!*         AA  AN INPUT ONE-DIMENSIONAL COMPLEX ARRAY DIMENSIONED AT   *CFNORM17
!*             LEAST N, WHICH CONTAINS THE VECTOR WHOSE NORM IS TO     *CFNORM18
!*             BE FOUND.                                               *CFNORM19
!*                                                                     *CFNORM20
!*          X  AN OUTPUT COMPLEX PARAMETER WHICH CONTAINS THE iMAX NORM *CFNORM21
!*             OF THE VECTOR.                                          *CFNORM22
!*                                                                     *CFNORM23
!*    REQUIRED ROUTINES:  NONE                                         *CFNORM24
!*                                                                     *CFNORM25
!*    FORTRAN FUNCTIONS:  ABS                                          *CFNORM26
!*                                                                     *CFNORM27
!*    SOURCE              COMPUTER SCIENCES CORPORATION,               *CFNORM28
!*                        HAMPTON, VA.                                 *CFNORM29
!*                                                                     *CFNORM30
!*    LANGUAGE:           FORTRAN                                      *CFNORM31
!*                                                                     *CFNORM32
!*    DATE RELEASED:      MAY 1,1979                                   *CFNORM33
!*                                                                     *CFNORM34
!*    LATEST REVISION:    MAY 1,1979                                   *CFNORM35
!*                                                                     *CFNORM36
!***********************************************************************CFNORM37
      COMPLEX AA(*),CiMAX                                                
      CiMAX=(0.,0.)                                                      
      DO 20 I=1,N                                                       
      IF(ABS(CiMAX)-ABS(AA(I))) 10,20,20                               
   10 CiMAX=AA(I)                                                        
   20 CONTINUE                                                          
      QXZ195=ABS(CiMAX)                                                 
      RETURN                                                            
      END                                                               

!**********************************************************************
      SUBROUTINE QXZ196(A,LDA,N,IPVT,RCOND,Z)                           
!**********************************************************************
      INTEGER LDA,N,IPVT(*)                                             
      COMPLEX A(LDA,*),Z(*),TEMP                                        
      REAL RCOND                                                        
!                                                                       !CGECO1 6
!     QXZ196 FACTORS A COMPLEX MATRIX BY GAUSSIAN ELIMINATION           !QXZN1973
!     AND ESTIMATES THE CONDITION OF THE MATRIX.                        !CGECO1 8
!                                                                       !CGECO1 9
!     IF  RCOND  IS NOT NEEDED, QXZ197 IS SLIGHTLY FASTER.              !QXZN1974
!     TO SOLVE  A*X = B , FOLLOW QXZ196 BY QXZ198.                      !QXZN1975
!     TO COMPUTE  INVERSE(A)*! , FOLLOW QXZ196 BY QXZ198.               !QXZN1976
!     TO COMPUTE  DETERMINANT(A) , FOLLOW QXZ196 BY CGEDI.              !QXZN1977
!     TO COMPUTE  INVERSE(A) , FOLLOW QXZ196 BY CGEDI.                  !QXZN1978
!                                                                       !CGECO115
!     ON ENTRY                                                          !CGECO116
!                                                                       !CGECO117
!        A       COMPLEX(LDA, N)                                        !CGECO118
!                THE MATRIX TO BE FACTORED.                             !CGECO119
!                                                                       !CGECO120
!        LDA     INTEGER                                                !CGECO121
!                THE LEADING DIMENSION OF THE ARRAY  A .                !CGECO122
!                                                                       !CGECO123
!        N       INTEGER                                                !CGECO124
!                THE ORDER OF THE MATRIX  A .                           !CGECO125
!                                                                       !CGECO126
!     ON RETURN                                                         !CGECO127
!                                                                       !CGECO128
!        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS         !CGECO129
!                WHICH WERE USED TO OBTAIN IT.                          !CGECO130
!                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE       !CGECO131
!                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER          !CGECO132
!                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.       !CGECO133
!                                                                       !CGECO134
!        IPVT    INTEGER(N)                                             !CGECO135
!                AN INTEGER VECTOR OF PIVOT INDICES.                    !CGECO136
!                                                                       !CGECO137
!        RCOND   REAL                                                   !CGECO138
!                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .        !CGECO139
!                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS       !CGECO140
!                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE             !CGECO141
!                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND . !CGECO142
!                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION     !CGECO143
!                           1.0 + RCOND .EQ. 1.0                        !CGECO144
!                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING           !CGECO145
!                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF         !CGECO146
!                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE          !CGECO147
!                UNDERFLOWS.                                            !CGECO148
!                                                                       !CGECO149
!        Z       COMPLEX(N)                                             !CGECO150
!                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.  !CGECO151
!                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS      !CGECO152
!                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT           !CGECO153
!                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .                    !CGECO154
!                                                                       !CGECO155
!     LINPACK. THIS VERSION DATED 03/11/78 .                            !CGECO156
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LABS.     !CGECO157
!                                                                       !CGECO158
!     REVISED BY COMPUTER SCIENCES CORPORATION                          !CGECO159
!     MAY 1,1979                                                        !CGECO160
!                                                                       !CGECO161
!     SUBROUTINES AND FUNCTIONS                                         !CGECO162
!                                                                       !CGECO163
!     LINPACK !CGEFA1                                                    !CGECO164
!     BLAS QXZ199,QXZ200,QXZ202,QXZ204                                  !QXZN1979
!     FORTRAN ABS,IMAG,MAX,CMPLX,CONJG,REAL                          !CGECO166
!                                                                       !CGECO167
!     INTERNAL VARIABLES                                                !CGECO168
!                                                                       !CGECO169
      COMPLEX QXZ200,EK,T,WK,WKM                                        
      REAL ANORM,S,QXZ204,SM,YNORM                                      
      INTEGER INFO,J,K,KB,KP1,L                                         
      COMPLEX ZDUM,ZDUM1,ZDUM2,CSIGN1                                   
      REAL CABS1                                                        
      CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(IMAG(ZDUM))                  
      CSIGN1(ZDUM1,ZDUM2) = CABS1(ZDUM1)*(ZDUM2/CABS1(ZDUM2))           
!                                                                       !CGECO178
!     COMPUTE 1-NORM OF A                                               !CGECO179
!                                                                       !CGECO180
      ANORM = 0.0E0                                                     
      DO 10 J = 1, N                                                    
         ANORM = MAX(ANORM,QXZ204(N,A(1,J)))                          
   10 CONTINUE                                                          
!                                                                       !CGECO185
!     FACTOR                                                            !CGECO186
!                                                                       !CGECO187
      CALL QXZ197(A,LDA,N,IPVT,INFO)                                    
!                                                                       !CGECO189
!     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .              !CGECO190
!     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  CTRANS(A)*Y = E . !CGECO191
!     CTRANS(A)  IS THE CONJUGATE TRANSPOSE OF A .                      !CGECO192
!     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE iMAXIMUM LOCAL           !CGECO193
!     GROWTH IN THE ELEMENTS OF W  WHERE  CTRANS(U)*W = E .             !CGECO194
!     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.            !CGECO195
!                                                                       !CGECO196
!     SOLVE CTRANS(U)*W = E                                             !CGECO197
!                                                                       !CGECO198
      EK = (1.0E0,0.0E0)                                                !CGECO199
      DO 20 J = 1, N                                                    !CGECO100
         Z(J) = (0.0E0,0.0E0)                                           !CGECO101
   20 CONTINUE                                                          !CGECO102
      DO 100 K = 1, N                                                   !CGECO103
         TEMP = -Z(K)                                                   
         IF (CABS1(Z(K)) .NE. 0.0E0) EK = CSIGN1(EK,TEMP)               
         IF (CABS1(EK-Z(K)) .LE. CABS1(A(K,K))) GO TO 30                
            S = CABS1(A(K,K))/CABS1(EK-Z(K))                            !CGECO106
            CALL QXZ202(N,S,Z)                                          
            EK = CMPLX(S,0.0E0)*EK                                      !CGECO108
   30    CONTINUE                                                       !CGECO109
         WK = EK - Z(K)                                                 !CGECO110
         WKM = -EK - Z(K)                                               !CGECO111
         S = CABS1(WK)                                                  !CGECO112
         SM = CABS1(WKM)                                                !CGECO113
         IF (CABS1(A(K,K)) .EQ. 0.0E0) GO TO 40                         !CGECO114
            WK = WK/CONJG(A(K,K))                                       !CGECO115
            WKM = WKM/CONJG(A(K,K))                                     !CGECO116
         GO TO 50                                                       !CGECO117
   40    CONTINUE                                                       !CGECO118
            WK = (1.0E0,0.0E0)                                          !CGECO119
            WKM = (1.0E0,0.0E0)                                         !CGECO120
   50    CONTINUE                                                       !CGECO121
         KP1 = K + 1                                                    !CGECO122
         IF (KP1 .GT. N) GO TO 90                                       !CGECO123
            DO 60 J = KP1, N                                            !CGECO124
               SM = SM + CABS1(Z(J)+WKM*CONJG(A(K,J)))                  !CGECO125
               Z(J) = Z(J) + WK*CONJG(A(K,J))                           !CGECO126
               S = S + CABS1(Z(J))                                      !CGECO127
   60       CONTINUE                                                    !CGECO128
            IF (S .GE. SM) GO TO 80                                     !CGECO129
               T = WKM - WK                                             !CGECO130
               WK = WKM                                                 !CGECO131
               DO 70 J = KP1, N                                         !CGECO132
                  Z(J) = Z(J) + T*CONJG(A(K,J))                         !CGECO133
   70          CONTINUE                                                 !CGECO134
   80       CONTINUE                                                    !CGECO135
   90    CONTINUE                                                       !CGECO136
         Z(K) = WK                                                      !CGECO137
  100 CONTINUE                                                          !CGECO138
      S = 1.0E0/QXZ204(N,Z(1))                                          
      CALL QXZ202(N,S,Z)                                                
!                                                                       !CGECO141
!     SOLVE CTRANS(L)*Y = W                                             !CGECO142
!                                                                       !CGECO143
      DO 120 KB = 1, N                                                  !CGECO144
         K = N + 1 - KB                                                 !CGECO145
         IF (K .LT. N) Z(K) = Z(K) + QXZ200(N-K,A(K+1,K),Z(K+1))        
!        IF (CABS1(REAL(Z(K)) ).LE. 1.0E0) GO TO 110                    
         if (abs(real(z(k))) .le. 1.0e0)goto 110
            S = 1.0E0/CABS1(Z(K))                                       !CGECO148
            CALL QXZ202(N,S,Z)                                          
  110    CONTINUE                                                       !CGECO150
         L = IPVT(K)                                                    !CGECO151
         T = Z(L)                                                       !CGECO152
         Z(L) = Z(K)                                                    !CGECO153
         Z(K) = T                                                       !CGECO154
  120 CONTINUE                                                          !CGECO155
      S = 1.0E0/QXZ204(N,Z(1))                                          
      CALL QXZ202(N,S,Z)                                                
!                                                                       !CGECO158
      YNORM = 1.0E0                                                     !CGECO159
!                                                                       !CGECO160
!     SOLVE L*V = Y                                                     !CGECO161
!                                                                       !CGECO162
      DO 140 K = 1, N                                                   !CGECO163
         L = IPVT(K)                                                    !CGECO164
         T = Z(L)                                                       !CGECO165
         Z(L) = Z(K)                                                    !CGECO166
         Z(K) = T                                                       !CGECO167
         IF (K .LT. N) CALL QXZ199(N-K,T,A(K+1,K),Z(K+1))               !QXZN1991
!        IF (CABS1(REAL(Z(K)) ).LE. 1.0E0) GO TO 130                    F5!CXGCO3
         if(abs(real(z(k))) .le. 1.0e0)goto 130
            S = 1.0E0/CABS1(Z(K))                                       !CGECO170
            CALL QXZ202(N,S,Z)                                          !QXZN1992
            YNORM = S*YNORM                                             !CGECO172
  130    CONTINUE                                                       !CGECO173
  140 CONTINUE                                                          !CGECO174
      S = 1.0E0/QXZ204(N,Z(1))                                          
      CALL QXZ202(N,S,Z)                                                !QXZN1994
      YNORM = S*YNORM                                                   !CGECO177
!                                                                       !CGECO178
!     SOLVE  U*Z = V                                                    !CGECO179
!                                                                       !CGECO180
      DO 160 KB = 1, N                                                  !CGECO181
         K = N + 1 - KB                                                 !CGECO182
         IF (CABS1(Z(K)) .LE. CABS1(A(K,K))) GO TO 150                  
            S = CABS1(A(K,K))/CABS1(Z(K))                               !CGECO184
            CALL QXZ202(N,S,Z)                                          !QXZN1995
            YNORM = S*YNORM                                             !CGECO186
  150    CONTINUE                                                       !CGECO187
         IF (CABS1(A(K,K)) .NE. 0.0E0) Z(K) = Z(K)/A(K,K)               !CGECO188
         IF (CABS1(A(K,K)) .EQ. 0.0E0) Z(K) = (1.0E0,0.0E0)             !CGECO189
         T = -Z(K)                                                      !CGECO190
       CALL QXZ199(K-1,T,A(1,K),Z(1))                                   !QXZN1996
  160 CONTINUE                                                          !CGECO192
!     MAKE ZNORM = 1.0                                                  !CGECO193
      S = 1.0E0/QXZ204(N,Z(1))                                          
      CALL QXZ202(N,S,Z)                                                !QXZN1998
      YNORM = S*YNORM                                                   !CGECO196
!                                                                       !CGECO197
      IF (ANORM .NE. 0.0E0) RCOND = YNORM/ANORM                         !CGECO198
      IF (ANORM .EQ. 0.0E0) RCOND = 0.0E0                               !CGECO199
      RETURN                                                            
      END      
	  
!**********************************************************************
      SUBROUTINE QXZ197(A,LDA,N,IPVT,INFO)                              !!FTN41865
!**********************************************************************
      INTEGER LDA,N,IPVT(*),INFO                                        !!FTN41866
      COMPLEX A(LDA,*)                                                  !!FTN41867
!                                                                       !CGEFA1 5
!     QXZ197 FACTORS A COMPLEX MATRIX BY GAUSSIAN ELIMINATION.          QXZN2002
!                                                                       !CGEFA1 7
!     QXZ197 IS USUALLY CALLED BY QXZ196, BUT IT CAN BE CALLED          QXZN2003
!     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.          !CGEFA1 9
!     (TIME FOR QXZ196) = (1 + 9/N)*(TIME FOR QXZ197) .                 QXZN2004
!                                                                       !CGEFA111
!     ON ENTRY                                                          !CGEFA112
!                                                                       !CGEFA113
!        A       COMPLEX(LDA, N)                                        !CGEFA114
!                THE MATRIX TO BE FACTORED.                             !CGEFA115
!                                                                       !CGEFA116
!        LDA     INTEGER                                                !CGEFA117
!                THE LEADING DIMENSION OF THE ARRAY  A .                !CGEFA118
!                                                                       !CGEFA119
!        N       INTEGER                                                !CGEFA120
!                THE ORDER OF THE MATRIX  A .                           !CGEFA121
!                                                                       !CGEFA122
!     ON RETURN                                                         !CGEFA123
!                                                                       !CGEFA124
!        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS         !CGEFA125
!                WHICH WERE USED TO OBTAIN IT.                          !CGEFA126
!                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE       !CGEFA127
!                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER          !CGEFA128
!                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.       !CGEFA129
!                                                                       !CGEFA130
!        IPVT    INTEGER(N)                                             !CGEFA131
!                AN INTEGER VECTOR OF PIVOT INDICES.                    !CGEFA132
!                                                                       !CGEFA133
!        INFO    INTEGER                                                !CGEFA134
!                = 0  NORMAL VALUE.                                     !CGEFA135
!                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR       !CGEFA136
!                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES        !CGEFA137
!                     INDICATE THAT QXZ198 OR CGEDI WILL DIVIDE BY ZERO QXZN2005
!                     IF CALLED.  USE  RCOND  IN QXZ196 FOR A RELIABLE  QXZN2006
!                     INDICATION OF SINGULARITY.                        !CGEFA140
!                                                                       !CGEFA141
!     LINPACK. THIS VERSION DATED 03/11/78 .                            !CGEFA142
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LABS.     !CGEFA143
!                                                                       !CGEFA144
!     REVISED BY COMPUTER SCIENCES CORPORATION                          !CGEFA145
!     MAY 1,1979                                                        !CGEFA146
!                                                                       !CGEFA147
!     SUBROUTINES AND FUNCTIONS                                         !CGEFA148
!                                                                       !CGEFA149
!     BLAS CAXPY,CSCAL,ICAiMAX                                           !CGEFA150
!     FORTRAN ABS,IMAG,REAL                                            !CGEFA151
!                                                                       !CGEFA152
!     INTERNAL VARIABLES                                                !CGEFA153
!                                                                       !CGEFA154
      COMPLEX T                                                         !CGEFA155
      INTEGER QXZ203,J,K,KP1,L,NM1                                      
!                                                                       !CGEFA157
      COMPLEX ZDUM                                                      !CGEFA158
      REAL CABS1                                                        !CGEFA159
      CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(IMAG(ZDUM))                  !!FTN41868
!                                                                       !CGEFA161
!     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING                        !CGEFA162
!                                                                       !CGEFA163
      INFO = 0                                                          !CGEFA164
      NM1 = N - 1                                                       !CGEFA165
      IF (NM1 .LT. 1) GO TO 70                                          !CGEFA166
      DO 60 K = 1, NM1                                                  !CGEFA167
         KP1 = K + 1                                                    !CGEFA168
!                                                                       !CGEFA169
!        FIND L = PIVOT INDEX                                           !CGEFA170
!                                                                       !CGEFA171
         L = QXZ203(N-K+1,A(K,K)) + K - 1                               
         IPVT(K) = L                                                    !CGEFA173
!                                                                       !CGEFA174
!        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED          !CGEFA175
!                                                                       !CGEFA176
         IF (CABS1(A(L,K)) .EQ. 0.0E0) GO TO 40                         !CGEFA177
!                                                                       !CGEFA178
!           INTERCHANGE IF NECESSARY                                    !CGEFA179
!                                                                       !CGEFA180
            IF (L .EQ. K) GO TO 10                                      !CGEFA181
               T = A(L,K)                                               !CGEFA182
               A(L,K) = A(K,K)                                          !CGEFA183
               A(K,K) = T                                               !CGEFA184
   10       CONTINUE                                                    !CGEFA185
!                                                                       !CGEFA186
!           COMPUTE MULTIPLIERS                                         !CGEFA187
!                                                                       !CGEFA188
            T = -(1.0E0,0.0E0)/A(K,K)                                   !CGEFA189
            CALL QXZ201(N-K,T,A(K+1,K))                                 
!                                                                       !CGEFA191
!           ROW ELIMINATION WITH COLUMN INDEXING                        !CGEFA192
!                                                                       !CGEFA193
            DO 30 J = KP1, N                                            !CGEFA194
               T = A(L,J)                                               !CGEFA195
               IF (L .EQ. K) GO TO 20                                   !CGEFA196
                  A(L,J) = A(K,J)                                       !CGEFA197
                  A(K,J) = T                                            !CGEFA198
   20          CONTINUE                                                 !CGEFA199
               CALL QXZ199(N-K,T,A(K+1,K),A(K+1,J))                     
   30       CONTINUE                                                    !CGEFA101
         GO TO 50                                                       !CGEFA102
   40    CONTINUE                                                       !CGEFA103
            INFO = K                                                    !CGEFA104
   50    CONTINUE                                                       !CGEFA105
   60 CONTINUE                                                          !CGEFA106
   70 CONTINUE                                                          !CGEFA107
      IPVT(N) = N                                                       !CGEFA108
      IF (CABS1(A(N,N)) .EQ. 0.0E0) INFO = N                            !CGEFA109
      RETURN                                                            !CGEFA110
      END                                                               !!FTN41869

!**********************************************************************
      SUBROUTINE QXZ198(A,LDA,N,IPVT,B,JOB)                             !!FTN41870
!**********************************************************************
      INTEGER LDA,N,IPVT(*),JOB                                         !!FTN41871
      COMPLEX A(LDA,*),B(*)                                             !!FTN41872
!                                                                       !CGESL1 5
!     QXZ198 SOLVES THE COMPLEX SYSTEM                                  QXZN2014
!     A * X = B  OR  CTRANS(A) * X = B                                  !CGESL1 7
!     USING THE FACTORS COMPUTED BY QXZ196 OR QXZ197.                   QXZN2015
!                                                                       !CGESL1 9
!     ON ENTRY                                                          !CGESL110
!                                                                       !CGESL111
!        A       COMPLEX(LDA, N)                                        !CGESL112
!                THE OUTPUT FROM QXZ196 OR QXZ197.                      QXZN2016
!                                                                       !CGESL114
!        LDA     INTEGER                                                !CGESL115
!                THE LEADING DIMENSION OF THE ARRAY  A .                !CGESL116
!                                                                       !CGESL117
!        N       INTEGER                                                !CGESL118
!                THE ORDER OF THE MATRIX  A .                           !CGESL119
!                                                                       !CGESL120
!        IPVT    INTEGER(N)                                             !CGESL121
!                THE PIVOT VECTOR FROM QXZ196 OR QXZ197.                QXZN2017
!                                                                       !CGESL123
!        B       COMPLEX(N)                                             !CGESL124
!                THE RIGHT HAND SIDE VECTOR.                            !CGESL125
!                                                                       !CGESL126
!        JOB     INTEGER                                                !CGESL127
!                = 0         TO SOLVE  A*X = B ,                        !CGESL128
!                = NONZERO   TO SOLVE  CTRANS(A)*X = B  WHERE           !CGESL129
!                            CTRANS(A)  IS THE CONJUGATE TRANSPOSE.     !CGESL130
!                                                                       !CGESL131
!     ON RETURN                                                         !CGESL132
!                                                                       !CGESL133
!        B       THE SOLUTION VECTOR  X .                               !CGESL134
!                                                                       !CGESL135
!     ERROR CONDITION                                                   !CGESL136
!                                                                       !CGESL137
!        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A   !CGESL138
!        ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES SINGULARITY  !CGESL139
!        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER       !CGESL140
!        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE     !CGESL141
!        CALLED CORRECTLY AND IF QXZ196 HAS SET RCOND .GT. 0.0          QXZN2018
!        OR QXZ197 HAS SET INFO .EQ. 0 .                                QXZN2019
!                                                                       !CGESL144
!     TO COMPUTE  INVERSE(A) * !  WHERE  C  IS A MATRIX                 !CGESL145
!     WITH  P  COLUMNS                                                  !CGESL146
!           CALL QXZ196(A,LDA,N,IPVT,RCOND,Z)                           QXZN2020
!           IF (RCOND IS TOO SMALL) GO TO ...                           !CGESL148
!           DO 10 J = 1, P                                              !CGESL149
!              CALL QXZ198(A,LDA,N,IPVT,C(1,J),0)                       QXZN2021
!        10 CONTINUE                                                    !CGESL151
!                                                                       !CGESL152
!     LINPACK. THIS VERSION DATED 03/11/78 .                            !CGESL153
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LABS.     !CGESL154
!                                                                       !CGESL155
!     REVISED BY COMPUTER SCIENCES CORPORATION                          !CGESL156
!     MAY 1,1979                                                        !CGESL157
!                                                                       !CGESL158
!     SUBROUTINES AND FUNCTIONS                                         !CGESL159
!                                                                       !CGESL160
!     BLAS CAXPY,CDOTC                                                  !CGESL161
!     FORTRAN CONJG                                                     !CGESL162
!                                                                       !CGESL163
!     INTERNAL VARIABLES                                                !CGESL164
!                                                                       !CGESL165
      COMPLEX QXZ200,T                                                  
      INTEGER K,KB,L,NM1                                                !CGESL167
!                                                                       !CGESL168
      NM1 = N - 1                                                       !CGESL169
      IF (JOB .NE. 0) GO TO 50                                          !CGESL170
!                                                                       !CGESL171
!        JOB = 0 , SOLVE  A * X = B                                     !CGESL172
!        FIRST SOLVE  L*Y = B                                           !CGESL173
!                                                                       !CGESL174
         IF (NM1 .LT. 1) GO TO 30                                       !CGESL175
         DO 20 K = 1, NM1                                               !CGESL176
            L = IPVT(K)                                                 !CGESL177
            T = B(L)                                                    !CGESL178
            IF (L .EQ. K) GO TO 10                                      !CGESL179
               B(L) = B(K)                                              !CGESL180
               B(K) = T                                                 !CGESL181
   10       CONTINUE                                                    !CGESL182
            CALL QXZ199(N-K,T,A(K+1,K),B(K+1))                          
   20    CONTINUE                                                       !CGESL184
   30    CONTINUE                                                       !CGESL185
!                                                                       !CGESL186
!        NOW SOLVE  U*X = Y                                             !CGESL187
!                                                                       !CGESL188
         DO 40 KB = 1, N                                                !CGESL189
            K = N + 1 - KB                                              !CGESL190
            B(K) = B(K)/A(K,K)                                          !CGESL191
            T = -B(K)                                                   !CGESL192
            CALL QXZ199(K-1,T,A(1,K),B(1))                              
   40    CONTINUE                                                       !CGESL194
      GO TO 100                                                         !CGESL195
   50 CONTINUE                                                          !CGESL196
!                                                                       !CGESL197
!        JOB = NONZERO, SOLVE  CTRANS(A) * X = B                        !CGESL198
!        FIRST SOLVE  CTRANS(U)*Y = B                                   !CGESL199
!                                                                       !CGESL100
         DO 60 K = 1, N                                                 !CGESL101
            T = QXZ200(K-1,A(1,K),B(1))                                 
            B(K) = (B(K) - T)/CONJG(A(K,K))                             !CGESL103
   60    CONTINUE                                                       !CGESL104
!                                                                       !CGESL105
!        NOW SOLVE CTRANS(L)*X = Y                                      !CGESL106
!                                                                       !CGESL107
         IF (NM1 .LT. 1) GO TO 90                                       !CGESL108
         DO 80 KB = 1, NM1                                              !CGESL109
            K = N - KB                                                  !CGESL110
            B(K) = B(K) + QXZ200(N-K,A(K+1,K),B(K+1))                   
            L = IPVT(K)                                                 !CGESL112
            IF (L .EQ. K) GO TO 70                                      !CGESL113
               T = B(L)                                                 !CGESL114
               B(L) = B(K)                                              !CGESL115
               B(K) = T                                                 !CGESL116
   70       CONTINUE                                                    !CGESL117
   80    CONTINUE                                                       !CGESL118
   90    CONTINUE                                                       !CGESL119
  100 CONTINUE                                                          !CGESL120
      RETURN                                                            !CGESL121
      END                                                               !!FTN41873

!**********************************************************************
      SUBROUTINE QXZ199(N,CA,CX,CY)                                     !!FTN41874
!**********************************************************************
!                                                                       CAXPY1 3
!     CONSTANT TIMES A VECTOR PLUS A VECTOR.                            CAXPY1 4
!     JACK DONGARRA, LINPACK, 3/11/78.                                  CAXPY1 5
!                                                                       CAXPY1 6
!     REVISED BY COMPUTER SCIENCES CORPORATION                          CAXPY1 7
!     MAY 1,1979                                                        CAXPY1 8
!                                                                       CAXPY1 9
      COMPLEX CX(*),CY(*),CA                                            !!FTN41875
      INTEGER I,N                                                       
!                                                                       CAXPY112
      IF(N.LE.0)RETURN                                                  
      IF (ABS(REAL(CA)) + ABS(IMAG(CA)) .EQ. 0.0 ) RETURN              
!                                                                       CAXPY115
!                                                                       CAXPY116
!        CODE FOR BOTH INCREMENTS EQUAL TO 1                            CAXPY117
!                                                                       CAXPY118
   20 DO 30 I = 1,N                                                     
        CY(I) = CY(I) + CA*CX(I)                                        
   30 CONTINUE                                                          
      RETURN                                                            
      END                                                               !!

!**********************************************************************
      COMPLEX FUNCTION QXZ200(N,CX,CY)                                  !!FTN41877
!**********************************************************************
!                                                                        
!     FORMS THE DOT PRODUCT OF TWO VECTORS, CONJUGATING THE FIRST        
!     VECTOR.                                                            
!     JACK DONGARRA, LINPACK,  3/11/78.                                  
!                                                                        
!     REVISED BY COMPUTER SCIENCES CORPORATION                           
!     MAY 1,1979                                                         
!                                                                       
      COMPLEX CX(*),CY(*),CTEMP                                         
      INTEGER I,N                                                       
!                                                                       
      CTEMP = (0.0,0.0)                                                 
      QXZ200 = (0.0,0.0)                                                
      IF(N.LE.0)RETURN                                                  
!                                                                       
!        CODE FOR BOTH INCREMENTS EQUAL TO 1                            
!                                                                       
   20 DO 30 I = 1,N                                                     
        CTEMP = CTEMP + CONJG(CX(I))*CY(I)                              
   30 CONTINUE                                                          
      QXZ200 = CTEMP                                                    
      RETURN                                                            
      END                                                               

!**********************************************************************
      SUBROUTINE  QXZ201(N,CA,CX)                                       !!FTN41880
!**********************************************************************
!                                                                       CSCAL1 3
!     SCALES A VECTOR BY A CONSTANT.                                    CSCAL1 4
!     JACK DONGARRA, LINPACK,  3/11/78.                                 CSCAL1 5
!                                                                       CSCAL1 6
!     REVISED BY COMPUTER SCIENCES CORPORATION                          CSCAL1 7
!     MAY 1,1979                                                        CSCAL1 8
!                                                                       CSCAL1 9
      COMPLEX CA,CX(*)                                                  !!
      INTEGER I,N                                                       
!                                                                       CSCAL112
      IF(N.LE.0)RETURN                                                  
!                                                                       CSCAL114
!        CODE FOR INCREMENT EQUAL TO 1                                  CSCAL115
!                                                                       CSCAL116
   20 DO 30 I = 1,N                                                     
        CX(I) = CA*CX(I)                                                
   30 CONTINUE                                                          
      RETURN                                                            
      END                                                               

!**********************************************************************
      SUBROUTINE  QXZ202(N,SA,CX)                                       
!**********************************************************************
!                                                                       
!     SCALES A COMPLEX VECTOR BY A REAL CONSTANT.                       
!     JACK DONGARRA, LINPACK, 3/11/78.                                  
!                                                                       
!     REVISED BY COMPUTER SCIENCES CORPORATION                          
!     MAY 1,1979                                                        
!                                                                       
      COMPLEX CX(*)                                                     !!
      REAL SA                                                           
      INTEGER I,N                                                       
!                                                                       
      IF(N.LE.0)RETURN                                                  
!                                                                       
!        CODE FOR INCREMENT EQUAL TO 1                                  
!                                                                       
   20 DO 30 I = 1,N                                                     
        CX(I) = CMPLX(SA*REAL(CX(I)),SA*IMAG(CX(I)))                   
   30 CONTINUE                                                          
      RETURN                                                            
      END                                                               !!

!**********************************************************************
      INTEGER FUNCTION QXZ203(N,CX)                                     !!FTN41886
!**********************************************************************
!                                                                       ICMAX3
!     FINDS THE INDEX OF ELEMENT HAVING iMAX. ABSOLUTE VALUE.            ICMAX4
!     JACK DONGARRA, LINPACK, 3/11/78.                                  ICMAX5
!                                                                       ICMAX6
!     REVISED BY COMPUTER SCIENCES CORPORATION                          ICMAX7
!     MAY 1,1979                                                        ICMAX8
!                                                                       ICMAX9
      COMPLEX CX(*)                                                     !!FTN41887
      REAL SiMAX                                                         
      INTEGER I,N                                                       
      COMPLEX ZDUM                                                      
      REAL CABS1                                                        
      CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(IMAG(ZDUM))                  !!
!                                                                       
      QXZ203 = 0                                                        
      IF( N .LT. 1 ) RETURN                                             
      QXZ203 = 1                                                        
      IF(N.EQ.1)RETURN                                                  
!                                                                       
!        CODE FOR INCREMENT EQUAL TO 1                                  
!                                                                       
   20 SiMAX = CABS1(CX(1))                                               
      DO 30 I = 2,N                                                     
         IF(CABS1(CX(I)) .LE.SiMAX) GO TO 30                             !CXGCO5
         QXZ203 = I                                                     
         SiMAX = CABS1(CX(I))                                            
   30 CONTINUE                                                          
      RETURN                                                            
      END                                                               !!FTN41890

!**********************************************************************
      REAL FUNCTION QXZ204(N,CX)                                        !!FTN41891
!**********************************************************************
!                                                                       
!     TAKES THE SUM OF THE ABSOLUTE VALUES OF A COMPLEX VECTOR AND      
!     RETURNS A SINGLE PRECISION RESULT.                                
!     JACK DONGARRA, LINPACK, 3/11/78.                                  
!                                                                       
!     REVISED BY COMPUTER SCIENCES CORPORATION                          
!     MAY 1,1979                                                        
!                                                                       
COMPLEX CX(*)                                                     
REAL STEMP                                                        
INTEGER I,N                                                       
!                                                                       
      QXZ204 = 0.0E0                                                    
      STEMP = 0.0E0                                                     
      IF(N.LE.0)RETURN                                                  
!                                                                       
!        CODE FOR INCREMENT EQUAL TO 1                                  
!                                                                       
   20 DO 30 I = 1,N                                                     
        STEMP = STEMP + ABS(REAL(CX(I))) + ABS(IMAG(CX(I)))            
   30 CONTINUE                                                          
      QXZ204 = STEMP                                                    
      RETURN                                                            
      END                                                               

!**********************************************************************
      SUBROUTINE IUNI(NiMAX,N,X,NTAB,Y,IORDER,X0,Y0,IPT,IERR)        
!**********************************************************************
      DIMENSION X(*),Y(NiMAX,*),Y0(*)                                    !!FTN41237
      NM1=N-1                                                           !IUNI 105
      IERR=0                                                            !IUNI 106
      J=1                                                               !IUNI 107
      DELX=X(2)-X(1)                                                    !FTN41239
!                                                                       !IUNI 109
      IF (IORDER .EQ. 0) GO TO 10                                       !IUNI 112
      IF (N.LT. 2) GO TO 20                                             !IUNI 113
      GO TO 50                                                          !IUNI 114
  10  IERR=-1                                                           !IUNI 115
      GO TO 30                                                          !IUNI 116
  20  IERR=-2                                                           !IUNI 117
  30  DO 40 NT=1,NTAB                                                   !IUNI 118
         Y0(NT)=Y(1,NT)                                                 !IUNI 119
  40     CONTINUE                                                       !IUNI 120
      RETURN                                                            !IUNI 121
  50  IF (IPT .GT. -1) GO TO 65                                         !IUNI 122
!                                                                       !IUNI 123
      IF (DELX .EQ. 0) GO TO 190                                        !IUNI 128
      IF (N .EQ. 2) GO TO 65                                            !IUNI 129

!                                                                       !IUNI 133
      DO 60 J=2,NM1                                                     !IUNI 134
         IF (DELX * (X(J+1)-X(J))) 190,190,60                           !IUNI 135
  60     CONTINUE                                                       !IUNI 136
!                                                                       !IUNI 137
  65  IF (IPT .LT. 1) IPT=1                                             !IUNI 140
      IF (IPT .GT. NM1) IPT=NM1                                         !IUNI 141
      IN= SIGN (1.0,DELX *( X0-X(IPT)))                                 !IUNI 142
  70  P= X(IPT) - X0                                                    !IUNI 143
      IF (P* (X(IPT +1)- X0)) 90,180,80                                 !IUNI 144
  80  IPT =IPT +IN                                                      !IUNI 145
!                                                                       !IUNI 146
      IF (IPT.GT.0 .AND. IPT .LT. N) GO TO 70                           !IUNI 149
      IERR=-4                                                           !IUNI 150
      IPT=IPT- IN                                                       !IUNI 151
!                                                                       !IUNI 152
!                                                                       !IUNI 155
  90  IF (IORDER .GT. 1) GO TO 120                                      !IUNI 156
!                                                                       !IUNI 157
      DO 100 NT=1,NTAB                                                  !IUNI 160
          Y0(NT)=Y(IPT,NT)+((Y(IPT+1,NT)- Y(IPT,NT))*(X0-X(IPT)))/(X(IPT+1)-X(IPT))                                      !IUNI 162
 100     CONTINUE                                                       !IUNI 163
      IF (IERR .EQ. -4) IPT=IPT+IN                                      !IUNI 164
      RETURN                                                            !IUNI 165
!                                                                       !IUNI 166
 120  IF (N .EQ. 2) GO TO 200                                           !IUNI 169
!                                                                       !IUNI 170
      IF (IPT .EQ. NM1) GO TO  140                                      !IUNI 174
      IF (IPT .EQ. 1) GO TO 130                                         !IUNI 175
      IF (DELX *(X0-X(IPT-1)).LT.DELX* (X(IPT+2)-X0)) GO TO 140         !IUNI 176
 130  L=IPT                                                             !IUNI 177
      GO TO 150                                                         !IUNI 178
 140  L=IPT -1                                                          !IUNI 179
 150  V1=X(L)-X0                                                        !IUNI 180
      V2=X(L+1)-X0                                                      !IUNI 181
      V3=X(L+2)-X0                                                      !IUNI 182
      DO 160 NT=1,NTAB                                                  !IUNI 183
      YY1=(Y(L,NT) * V2 - Y(L+1,NT) * V1)/(X(L+1) - X(L))               !IUNI 184
      YY2=(Y(L+1,NT)*V3-Y(L+2,NT) *V2)/(X(L+2)-X(L+1))                  !IUNI 185
 160  Y0(NT)=(YY1*V3 -YY2*V1)/(X(L+2)-X(L))                             !IUNI 186
      IF (IERR .EQ. -4) IPT=IPT + IN                                    !IUNI 187
      RETURN                                                            !IUNI 188
 180  IF(P .NE. 0) IPT=IPT +1                                           !IUNI 189
      DO 185 NT=1,NTAB                                                  !IUNI 190
         Y0(NT)=Y(IPT,NT)                                               !IUNI 191
 185     CONTINUE                                                       !IUNI 192
      RETURN                                                            !IUNI 193
!                                                                       !IUNI 194
 190  IERR=J +1                                                         !IUNI 198
      RETURN                                                            !IUNI 199
 200  IERR=-3                                                           !IUNI 200
      RETURN                                                            !IUNI 201
      END                                                               !IUNI 202

!**********************************************************************
      SUBROUTINE SPLDER(MNPTS,N,NCVS,MiMAX,M,X,Y,T,F,DER1,DER2,IW,WK,IERR)
!**********************************************************************
      DIMENSION F(MiMAX,NCVS ),DER1(MiMAX,NCVS),DER2(MiMAX,NCVS)           !SPLDER 5
      DIMENSION WK(*)                                                   !FTN41195
      DIMENSION X(N),Y(MNPTS,NCVS),T(M)                                 !SPLDER 7
      IERR=0                                                            !SPLDER11
      JK=N*NCVS                                                         !SPLDER12
      IT=0                                                              !SPLDER13
      IS=0                                                              !SPLDER14
      IP=0                                                              !SPLDER15
      ID=0                                                              !SPLDER16
      IM=1                                                              !SPLDER17
      IN=1                                                              !SPLDER18
!                                                                       !SPLDER19
!     WK(1 TO N)        H                                               !SPLDER20
!                                                                       !SPLDER21
!     WK(N+1 TO JK+N)         IDLY                                      !SPLDER22
!                                                                       !SPLDER23
!     WK(JK+N+1 TO JK+2N)     ICC                                       !SPLDER24
!                                                                       !SPLDER25
!     WK(JK+2N+1 TO JK+3N)    IH2                                       !SPLDER26
!                                                                       !SPLDER27
!     WK(JK+3N+1 TO JK+4N)    IDSQ                                      !SPLDER28
!                                                                       !SPLDER29
!     WK(JK+4N+1 TO JK+5N)    IDD                                       !SPLDER30
!                                                                       !SPLDER31
!     WK(JK+5N+1 TO 2JK+5N)   IS2                                       !SPLDER32
!                                                                       !SPLDER33
!     WK(2JK+5N+1 TO 3JK+5N) IS3                                        !SPLDER34
!     WK(3JK+5N+1 TO 3JK+6N) IWW                                        !SPLDER35
!                                                                       !SPLDER36
!     WK(3JK+6N+1 TO 3JK+7N) IGG                                        !SPLDER37
!                                                                       !SPLDER38
!     WK(3JK+7N+1 TO 3JK+8N) ISV                                        !SPLDER39
!                                                                       !SPLDER40
!     DIMENSION WK (3JK+8N)                                             !SPLDER41
!                                                                       !SPLDER42
      IDLY=N+1                                                          !SPLDER43
      ICC=JK+N+1                                                        !SPLDER44
      IH2=JK+2*N+1                                                      !SPLDER45
      IDSQ=JK+3*N+1                                                     !SPLDER46
      IDD=JK+4*N+1                                                      !SPLDER47
      IS2=JK+5*N+1                                                      !SPLDER48
      IS3=2*JK+5*N+1                                                    !SPLDER49
      IWW=3*JK+5*N+1                                                    !SPLDER50
      IGG=3*JK+6*N+1                                                    !SPLDER51
      ISV=3*JK+7*N+1                                                    !SPLDER52
      IF(IW-(-1))15,2,15                                                !SPLDER53
    2 N1=N-1                                                            !SPLDER54
      IW=1                                                              !SPLDER55
!                                                                       !SPLDER56
!     IS THE INDEPENDENT VARIABLE ARRAY INCREASING                      !SPLDER57
!                                                                       !SPLDER58
      DO 118 I=1,N1                                                     !SPLDER59
      IF(X(I)-X(I+1)) 118,119,119                                       !SPLDER60
  119 PRINT 120,I,X(I),X(I+1)                                           !SPLDER61
  120 FORMAT(1H0,6X,64HINDEPENDENT VARIABLE ARRAY NOT INCREASING IN SPLDER AT POSITION ,I4,2X,2HX=,2F10.4)                                !FTN41198
      IERR=1                                                            !SPLDER64
      RETURN                                                            !SPLDER65
  118 CONTINUE                                                          !SPLDER66
      DO 101 L=1,NCVS                                                   !SPLDER67
    3 DO 51 I=1,N1                                                      !SPLDER68
      WK(I)=X(I+1)-X(I)                                                 !SPLDER69
      II=I+1                                                            !SPLDER70
      IM1=I-1                                                           !SPLDER71
      WK(IDLY+IP)=(Y(II,L)-Y(I,L))/WK(I)                                !SPLDER72
      IP=IP+1                                                           !SPLDER73
      WK(ICC+IM1)=WK(I)                                                 !SPLDER74
   51 CONTINUE                                                          !SPLDER75
      DO 65 I=2,N1                                                      !SPLDER76
      IM1=I-1                                                           !SPLDER77
      WK(IH2+IM1)=(WK(IM1)+WK(I))*2.                                    !SPLDER78
      !WK(IDSQ+IM1)=(WK(ID LY+IM)-WK(ID LY+IM-1))*6                      !SPLDER79
	  WK(IDSQ+IM1)=(WK(IDLY+IM)-WK(IDLY+IM-1))*6                      !SPLDER79
      IM=IM+1                                                           !SPLDER80
   65 CONTINUE                                                          !SPLDER81
  222 CONTINUE                                                          !SPLDER82
      WK(IH2)=1.                                                        !SPLDER83
      WK(IDSQ-1)=1.                                                     !SPLDER84
      WK(ICC)=0.                                                        !SPLDER85
      WK(N1)=0.                                                         !SPLDER86
      WK(IDSQ)=0.                                                       !SPLDER87
      WK(IDD-1)=0.                                                      !SPLDER88
  223 CONTINUE                                                          !SPLDER89
      KKM1=IBC-1                                                        !SPLDER90
!                                                                       !SPLDER91
!C THIS ROUTINE SOLVES THE TRIDIAGONAL (EXCEPT TWO ELEMENTS)   MATRIX    !SPLDER92
!                                                                       !SPLDER93
      IIP=ISV-1                                                         !SPLDER94
      WK(IWW)=WK(IH2)                                                   !SPLDER95
      WK(ISV)=WK(ICC)/WK(IH2)                                           !SPLDER96
      WK(IGG)=WK(IDSQ)/WK(IWW)                                          !SPLDER97
      DO 100 K=2,N                                                      !SPLDER98
      KM2=K-2                                                           !SPLDER99
      KM1 = K-1                                                         !SPLDE100
      WK(IWW+KM1)=WK(IH2+KM1)-WK(KM1)*WK(ISV+KM2)                       !SPLDE101
      IF (K.EQ.N)   GO TO 5                                             !SPLDE102
    4 WK(ISV+KM1)=WK(ICC+KM1)/WK(IWW+KM1)                               !SPLDE103
    5 WK(IGG+KM1)=(WK(IDSQ+KM1)-WK(KM1)*WK(IGG+KM2))/WK(IWW+KM1)        !SPLDE104
  100 CONTINUE                                                          !SPLDE105
      WK(IS2-1)=WK(ISV-1)                                               !SPLDE106
      IBW=IS2-1                                                         !SPLDE107
      DO 200 K=1,N1                                                     !SPLDE108
      IBC=IBW-K                                                         !SPLDE109
      KK= N-K                                                           !SPLDE110
      KKM1=KK-1                                                         !SPLDE111
      WK(IBC)=WK(IGG+KKM1)-WK(ISV+KKM1)*WK(IBC+1)                       !SPLDE112
  200 CONTINUE                                                          !SPLDE113
!                                                                       !SPLDE114
!                                                                       !SPLDE115
      DO 66 I=1,N                                                       !SPLDE116
      IM1=I-1                                                           !SPLDE117
      WK(IS2+ID)=WK(IDD+IM1)                                            !SPLDE118
      ID=ID+1                                                           !SPLDE119
   66 CONTINUE                                                          !SPLDE120
      WK(N1)=WK(IH2-2)                                                  !SPLDE121
   14 DO 53 I=1,N1                                                      !SPLDE122
      II=I+1                                                            !SPLDE123
      WK(IS3+IS)=(WK(IS2+IN)-WK(IS2+IN-1))/WK(I)                        !SPLDE124
      IS=IS+1                                                           !SPLDE125
      IN=IN+1                                                           !SPLDE126
   53 CONTINUE                                                          !SPLDE127
      IM=IM+2                                                           !SPLDE128
      IP=IP+1                                                           !SPLDE129
      IN=IN+1                                                           !SPLDE130
      IS=IS+1                                                           !SPLDE131
  101 CONTINUE                                                          !SPLDE132
   15 CONTINUE                                                          !SPLDE133
  104 J=0                                                               !SPLDE134
  105 J=J+1                                                             !SPLDE135
   16 I=IW                                                              !SPLDE136
   54 IF(T(J)-X(1)) 58,117,55                                           !SPLDE137
  117 I=1                                                               !SPLDE138
      IW=I                                                              !SPLDE139
      GO TO 17                                                          !SPLDE140
   55 IF(T(J)-X(N)) 162,59,58                                           !SPLDE141
  162 IF(T(J)-X(I)) 160,217,56                                          !SPLDE142
  160 I=I-1                                                             !SPLDE143
      GO TO 162                                                         !SPLDE144
   56 IF(T(J)-X(I)) 60,217,57                                           !SPLDE145
   57 I=I+1                                                             !SPLDE146
       GO TO 56                                                         !SPLDE147
   58 CONTINUE                                                          !SPLDE148
      IERR=2                                                            !SPLDE149
      RETURN                                                            !SPLDE150
   59 I=N                                                               !SPLDE151
   60 CONTINUE                                                          !SPLDE152
      I=I-1                                                             !SPLDE153
  217 IW=I                                                              !SPLDE154
   17 CONTINUE                                                          !SPLDE155
      IM1=I-1                                                           !SPLDE156
      ITD=IDLY+IM1                                                      !SPLDE157
      IT2=IS2+IM1                                                       !SPLDE158
      IT3=IS3+IM1                                                       !SPLDE159
      DO 110 K=1,NCVS                                                   !SPLDE160
      HT1=T(J)-X(I)                                                     !SPLDE161
      II=I+1                                                            !SPLDE162
      HT2=T(J)-X(II)                                                    !SPLDE163
      PROD=HT1 * HT2                                                    !SPLDE164
      DER2(J,K)=WK(IT2)+HT1*WK(IT3)                                     !SPLDE165
      DELSQS=(WK(IT2)+WK(IT2+1)+DER2(J,K))/6.                           !SPLDE166
      F(J,K)=Y(I,K)+HT1*WK(ITD)+PROD*DELSQS                             !SPLDE167
      DER1(J,K)=WK(ITD)+(HT1+HT2)*DELSQS+PROD*WK(IT3)/6.                !SPLDE168
      IT3=IT3+N                                                         !SPLDE169
      IT2=IT2+N                                                         !SPLDE170
      ITD=ITD+N                                                         !SPLDE171
  110 CONTINUE                                                          !SPLDE172
   61 CONTINUE                                                          !SPLDE173
      IF(J.LT.M)GO TO 105                                               !SPLDE174
      RETURN                                                            !SPLDE175
      END                                                               !SPLDE176


