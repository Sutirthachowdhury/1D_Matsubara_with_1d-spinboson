SUBROUTINE DIAG(EVALUES,EVECT,CRV)
 
    ! this subroutine only for IBMSP
           
    !c
    !c     CRV: HERMITIAN MATRIX (INPUT)
    !c     EVECT: EIGENVECTORS (OUTPUT)
    !c     EVALUES: EIGENVALUES (OUTPUT)
    !c
    !c     THIS ROUTINE SOLVES THE EIGENVALUE PROBLEM BY CALLING
    !c     THE FOLLOWING IBM ESSL FUNCTION
    !c
    !c     SSPEV, DSPEV, CHPEV, and ZHPEV--Eigenvalues and, Optionally,
    !c     the Eigenvectors of a Real Symmetric Matrix or a Complex
    !c     Hermitian Matrix
    !c
          use global,only : nstate
          IMPLICIT NONE
    
           CHARACTER JOBZ,UPLO
           INTEGER INFO
           INTEGER N,NAP,LDZ
    
           INTEGER I,J,IND
     
           real*8  AP(nstate*(nstate+1)/2),WORK(3*nstate)
    
           !real*8,allocatable:: EVALUES(:)
           !real*8,allocatable:: CRV(:,:),EVECT(:,:)
          
           !allocate(EVALUES(nstate))
           !allocate(CRV(nstate,nstate),EVECT(nstate,nstate))
           
           real*8 EVALUES(nstate)
           real*8 CRV(nstate,nstate),EVECT(nstate,nstate)
    
           N=nstate  ! electronic part of the ham
    
           NAP=N*(N+1)/2
           LDZ=N
          
           EVALUES=0.
           EVECT=0.
    
           JOBZ='V' ! calculate both eigenvalue and eogenvector
           UPLO='L' ! lower diagonal matrix
    
           IND=0
    
           DO J=1,N
              DO I=J,N
                 IND=IND+1
                 AP(IND)=CRV(I,J)
              END DO
           END DO
    
           CALL DSPEV(JOBZ,UPLO,N,AP,EVALUES,EVECT,LDZ,WORK,INFO)
     
           RETURN
         END SUBROUTINE DIAG
    