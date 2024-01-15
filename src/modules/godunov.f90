module godunov
  use modules, only: rk
  implicit none
  
  contains
  
  
  function flux_godunov(AK,CP,P1,P2,R1,R2,U1,U2,UDOT) result(flux)
  real(kind=rk), intent(in):: AK,CP,P1,P2,R1,R2,U1,U2,UDOT
  real(kind=rk)            :: UBIG,PBIG,RBIG1,RBIG2,DL1,DL2,DP1,DP2
  real(kind=rk)            :: Pf, Uf, Rf, CV, RM
  real(kind=rk)            :: flux(3), w(3)
  
  RM=CP*(1.0-1.0/AK)
  CV=CP-RM
    
  call RASPAD(AK,P1,R1,U1,P2,R2,U2,    &
                PBIG,UBIG,RBIG1,RBIG2,DL1,DL2,DP1,DP2)
  call POTOK(AK,P1,R1,U1,P2,R2,U2,     & 
                PBIG,UBIG,RBIG1,RBIG2, &
                DL1,DL2,DP1,DP2,UDOT,  &
                Pf,Uf,Rf)

  w(1)   = Rf
  w(2)   = Rf*Uf
  w(3)   = Rf*(CV*Pf/(Rf*RM) + Uf**2/2)
  
  flux(1) = Rf*Uf
  flux(2) = Rf*Uf*Uf + Pf
  flux(3) = Rf*Uf* (CP*Pf/(Rf*RM)+Uf*Uf/2) 
  
  flux = flux - UDOT*w
  
  end function flux_godunov
  
  
  subroutine RASPAD(AK,P1,R1,U1,P2,R2,U2, &
                  PBIG,UBIG,RBIG1,RBIG2,DL1,DL2,DP1,DP2)
!
! INPUT 
!
  real(kind=rk),intent(in)    :: AK,P1,R1,U1,P2,R2,U2
!      
! OUTPUT 
!
  real(kind=rk),intent(inout) :: PBIG,UBIG       ! parameters on contact surface
  real(kind=rk),intent(inout) :: RBIG1,RBIG2     ! density right and left of the surface
  real(kind=rk),intent(inout) :: DL1,DL2,DP1,DP2 ! speed of characteristics
!
! LOCALS
!
  logical                     :: FLAG
  real(kind=rk)               :: AK1,AK2,C1,C2,ST1,ST2,A1,A2,C2Z,C1Z,PI,EPSP,EPS
  real(kind=rk)               :: PBIGN

!      REAL,EXTERNAL :: FGOD,DIVFGOD

      EPSP=1.0E-8
      EPS=1.0E-8

      AK1=AK  
      AK2=AK  
      C1=SQRT(AK1*P1/R1)
      C2=SQRT(AK2*P2/R2)

      IF ( U1-U2+2.0*(C1/(AK1-1.0)+C2/(AK2-1.0)) < 0.0 ) THEN
        PBIG=1.0E-05
        RBIG1=0.0
        RBIG2=0.0
        RETURN
      ENDIF  

      FLAG=.FALSE.
      PBIG=(P1*R2*C2+P2*R1*C1+(U1-U2)*R1*C1*R2*C2)/(R1*C1+R2*C2)
      IF (PBIG<0.0) PBIG=1.0E-05
1     PBIGN=PBIG-(FGOD(PBIG,P1,R1,AK1)+FGOD(PBIG,P2,R2,AK2)-(U1-U2))/ &
                (DIVFGOD(PBIG,P1,R1,AK1)+DIVFGOD(PBIG,P2,R2,AK2))
            
      IF (PBIGN<0.0) THEN
        WRITE (*,*) PBIG,FLAG
        PAUSE
        PBIG=1.0E-05
        RBIG1=0.0
        RBIG2=0.0
        IF (FLAG) RETURN
        FLAG=.TRUE.
        GO TO 1
      ENDIF
      IF (ABS(PBIGN/PBIG-1.0) > EPS) THEN
        PBIG=PBIGN
        GO TO 1
      ENDIF

      ST1=0.5*(AK1-1.0)/AK1
      ST2=0.5*(AK2-1.0)/AK2
      IF ( PBIG >= P1-EPS ) THEN
        A1=SQRT(R1*(0.5*(AK1+1.0)*PBIG+0.5*(AK1-1.0)*P1))
      ELSE
        PI=PBIG/P1
        A1=0.5*(AK1-1.0)/AK1*R1*C1*(1.0-PI)/(1.0-PI**ST1)
      ENDIF
      IF ( PBIG >= P2-EPS ) THEN
        A2=SQRT(R2*(0.5*(AK2+1.0)*PBIG+0.5*(AK2-1.0)*P2))
      ELSE
        PI=PBIG/P2
        A2=0.5*(AK2-1.0)/AK2*R2*C2*(1.0-PI)/(1.0-PI**ST2)
      ENDIF

      UBIG=(A1*U1+A2*U2+P1-P2)/(A1+A2)
      
      IF ( P1-EPS<PBIG ) THEN
        RBIG1=R1*A1/(A1-R1*(U1-UBIG))
        DL1=U1-A1/R1
        DL2=DL1
      ELSE
        C1Z=C1+0.5*(AK1-1.0)*(U1-UBIG)
        RBIG1=AK1*PBIG/(C1Z*C1Z)
        DL1=U1-C1
        DL2=UBIG-C1Z
      ENDIF
      IF ( P2-EPS<PBIG ) THEN
        RBIG2=R2*A2/(A2+R2*(U2-UBIG))
        DP1=U2+A2/R2
        DP2=DP1
      ELSE
        C2Z=C2-0.5*(AK2-1.0)*(U2-UBIG)
        RBIG2=AK2*PBIG/(C2Z*C2Z)
        DP2=U2+C2
        DP1=UBIG+C2Z
      ENDIF
      
    end subroutine RASPAD


  function FGOD(PBIG,P,R,AK)
  real(kind=rk),intent(in):: PBIG,P,R,AK
  real(kind=rk)           :: PI,ST,FGOD

      PI=PBIG/P
      ST=0.5*(AK-1.0)/AK
      IF (PI >= 1.0) THEN
        FGOD=SQRT(P/R)*(PI-1.0)/SQRT(0.5*(AK+1.0)*PI+0.5*(AK-1.0))
      ELSE
        FGOD=2.0/(AK-1.0)*SQRT(AK*P/R)*(PI**ST-1.0)
      ENDIF
  end function FGOD


  function DIVFGOD(PBIG,P,R,AK)
  real(kind=rk),intent(in):: PBIG,P,R,AK
  real(kind=rk)           :: PI,ST,DIVFGOD  
      
      PI=PBIG/P
      ST=0.5*(AK-1.0)/AK
      IF (PI.GE.1.0) THEN
        DIVFGOD=((AK+1.0)*PI+(3.0*AK-1.0)) &
                /(4.0*AK*R*SQRT(AK*P/R)*   &
         SQRT((0.5*(AK+1.0)/AK*PI+0.5*(AK-1.0)/AK)**3))
      ELSE
         DIVFGOD=1.0/AK/PBIG*SQRT(AK*P/R)*PI**ST
      ENDIF
  end function DIVFGOD


  subroutine POTOK(AK,P1,R1,U1,P2,R2,U2,      &
                      PBIG,UBIG,RBIG1,RBIG2, &
                      DL1,DL2,DP1,DP2,UDOT,  &
                      Pf,Uf,Rf) 
!
! INPUT 
!
  real(kind=rk),intent(in)   :: AK,P1,R1,U1,P2,R2,U2,UDOT
  real(kind=rk),intent(in)   :: PBIG,UBIG                 ! parameters on contact surface
  real(kind=rk),intent(in)   :: RBIG1,RBIG2               ! density right and left of the surface
  real(kind=rk),intent(in)   :: DL1,DL2,DP1,DP2           ! speed of characteristics
!      
! OUTPUT 
!
  real(kind=rk),intent(inout):: Pf,Uf,Rf
!
! LOCALS
!
  real(kind=rk)              :: AK1,AK2
  real(kind=rk)              :: CZJM,CZJP,CZT

      AK1=AK
      AK2=AK
      Uf=UBIG
      Pf=PBIG
      if (UBIG-UDOT > 0. ) then    

        if (DL1-DL2 .eq. 0.) then  
          if (DL1-UDOT >= 0.) then 
            Rf=R1
            Uf=U1
            Pf=P1
          else                     
            Rf=RBIG1
          endif
        else                       
          if (DL2-UDOT <= 0.) then
            Rf=RBIG1             
          else
            if (DL1-UDOT >= 0.) then  
              Rf=R1              
              Uf=U1
              Pf=P1
            else 
              CZJM=SQRT(AK1*P1/R1)
              CZT=((AK1-1.)*(U1-UDOT)+2.*CZJM)/(AK1+1.)
              Uf=CZT+UDOT
              Pf=P1*(CZT/CZJM)**(2.*AK1/(AK1-1.))
              Rf=AK1*Pf/(CZT*CZT)
            endif      
          endif      
        endif      

      else                        

        if (DP1-DP2 .eq. 0.) then  
          if (DP1-UDOT <= 0.) then
            Rf=R2
            Uf=U2
            Pf=P2
          else                  
            Rf=RBIG2
          endif
        else                       
          if (DP1-UDOT >=0.) then
            Rf=RBIG2             
          else
            if (DP2-UDOT <= 0.) then  
              Rf=R2             
              Uf=U2
              Pf=P2
            else                 
              CZJP=SQRT(AK2*P2/R2)
              CZT=(-(AK2-1.)*(U2-UDOT)+2.*CZJP)/(AK2+1.)
              Uf=-CZT+UDOT
              Pf=P2*(CZT/CZJP)**(2.*AK2/(AK2-1.))
              Rf=AK2*Pf/(CZT*CZT)
            endif      
          endif      
        endif      
      
      endif
      
  end subroutine POTOK
  

end module godunov