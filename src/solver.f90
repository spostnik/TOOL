MODULE Solver
 USE Params, ONLY : qp,db,pi,rhos,do, linkpoints, &
	    wflag, xseed,linkdots,yseed,iseed,Rkm, &
	    Love,MomI,nseed,nDen,BE,rss
 USE Chain
 USE Data
IMPLICIT NONE
 SAVE
 REAL(db), PARAMETER :: KM3inFM3=1.0D+54, toMsc2=KM3inFM3/1.1158D+60
CONTAINS
 SUBROUTINE link_TOV(mout,r2out,xx1,d1,x2,d2,stp,yLout,nn,bb)
  INTEGER,INTENT(OUT) :: stp
  REAL(db),INTENT(IN) :: xx1,x2,d1,d2
  REAL(qp),INTENT(INOUT) :: mout,r2out,yLout,nn,bb
  REAL(qp) :: RHS, r,dmdx,dr2dx, dx, alph, bet,y1,mo,&
	    r2o,xo,fo,yLo,pp,rr,El,nup,Q,W,dyLdx,yo
  REAL(qp) :: xx(linkdots),yy(linkdots),ff(linkdots),zz(linkdots+1,3),x1
  INTEGER :: i
  x1=xx1
  bet=do/EXP(x1)
  alph=(d2-d1)/(x2-x1)
  y1=LOG(do+EXP(x1))+d1
  dx=(x2-x1)/linkdots
  DO i=1,linkdots-1
   xx(i)=x1+dx*i
  ENDDO
  xx(linkdots)=x2
  CALL link_eos(x1,y1,xx,yy,ff,alph,bet,stp)
  IF(nDen.AND.nn/=0) CALL link_n(x1,y1,xx,yy,nn,zz)
  IF(stp<0) RETURN
  mo=mout
  r2o=r2out
  fo=1/(1+EXP(y1-x1))
  xo=x1
  zz(1,1)=Rkm*EXP(r2out/2)
  zz(1,2)=mout*Rkm
  IF(Love) THEN
   yo=y1
   yLo=yLout
   DO i=1,linkdots
    r=EXP(r2o/2)

    pp=rss*EXP(xo)
    RHS=(EXP(r2o/2)-2*mo)/(mo+4*pi*EXP(3*r2o/2)*pp)
    dmdx=-4*pi*(1-fo)*RHS*pp*EXP(3*r2o/2)
    dr2dx=-fo*2*RHS
 
    El=1/(1-2*mo/r)
    rr=rss*EXP(yo)
    nup=2*El*(mo+4*pi*r**3*pp)/r**2
    Q=4*pi*El*(5*rr+9*pp+(rr+pp)*rr*(alph+1/(1+EXP(x1-xo)*bet))/pp)-6*El/r**2-nup**2
    W=El*(1+4*pi*r**2*(pp-rr))
    dyLdx=(yLo*yLo+yLo*W+r**2*Q)*(r-2*mo)*fo/(mo+4*pi*r**3*pp)
    
    dx=xx(i)-xo
    mout=mo+dmdx*dx/2
    r2out=r2o+dr2dx*dx/2
    yLout=yLo+dyLdx*dx/2
    r=EXP(r2out/2)

    pp=rss*EXP((xx(i)+xo)/2)
    RHS=(EXP(r2out/2)-2*mout)/(mout+4*pi*EXP(3*r2out/2)*pp)
    dmdx=-4*pi*(1-(ff(i)+fo)/2)*RHS*pp*EXP(3*r2out/2)
    dr2dx=-(ff(i)+fo)*RHS

    El=1/(1-2*mout/r)
    rr=rss*EXP((yy(i)+yo)/2)
    nup=2*El*(mout+4*pi*r**3*pp)/r**2
    Q=4*pi*El*(5*rr+9*pp+(rr+pp)*rr*(alph+1/(1+EXP(x1-(xx(i)+xo)/2)*bet))/pp)-6*El/r**2-nup**2
    W=El*(1+4*pi*r**2*(pp-rr))
    dyLdx=(yLout**2+yLout*W+r**2*Q)*(r-2*mout)*(ff(i)+fo)/(2*(mout+4*pi*r**3*pp))

    mout=mo+dmdx*dx
    r2out=r2o+dr2dx*dx
    yLout=yLo+dyLdx*dx
    IF(mout<mo.OR.r2out<r2o) THEN
!    WRITE(*,*) 'Error stp=0: decreasing mass or radius profile.'
     stp=0 
     RETURN
    ENDIF
    mo=mout
    r2o=r2out
    yLo=yLout
    fo=ff(i)
    xo=xx(i)
    yo=yy(i)
    zz(i+1,1)=Rkm*EXP(r2out/2)
    zz(i+1,2)=mout*Rkm
   ENDDO
  ELSE
   DO i=1,linkdots
    pp=rss*EXP(xo)
    RHS=(EXP(r2o/2)-2*mo)/(mo+4*pi*EXP(3*r2o/2)*pp)
    dmdx=-4*pi*(1-fo)*RHS*pp*EXP(3*r2o/2)
    dr2dx=-fo*2*RHS
 
    dx=xx(i)-xo
    mout=mo+dmdx*dx/2
    r2out=r2o+dr2dx*dx/2
    r=EXP(r2out/2)

    pp=rss*EXP((xx(i)+xo)/2)
    RHS=(EXP(r2out/2)-2*mout)/(mout+4*pi*EXP(3*r2out/2)*pp)
    dmdx=-4*pi*(1-(ff(i)+fo)/2)*RHS*pp*EXP(3*r2out/2)
    dr2dx=-(ff(i)+fo)*RHS

    mout=mo+dmdx*dx
    r2out=r2o+dr2dx*dx
    IF(mout<mo.OR.r2out<r2o) THEN
!    WRITE(*,*) 'Error stp=0: decreasing mass or radius profile.'
     stp=0 
     RETURN
    ENDIF
    mo=mout
    r2o=r2out
    fo=ff(i)
    xo=xx(i)
    zz(i+1,1)=Rkm*EXP(r2out/2)
    zz(i+1,2)=mout*Rkm
   ENDDO
   yLout=-1
  ENDIF
  IF(BE) THEN
   xo=zz(1,1)**2*zz(1,3)/SQRT(1+2*zz(1,2)/zz(1,1))
   DO i=1,linkdots
    fo=xo
    xo=zz(i+1,1)**2*zz(i+1,3)/SQRT(1+2*zz(i+1,2)/zz(i+1,1))
    bb=bb+(xo+fo)*(zz(i+1,1)-zz(i,1))/2
   ENDDO
  ENDIF
 END SUBROUTINE link_TOV

 SUBROUTINE link_eos(x1,y1,xin,yout,fout,alpha,beta,stp)
  INTEGER,INTENT(OUT) :: stp
  REAL(qp),INTENT(IN) :: xin(:),beta,alpha,x1,y1
  REAL(qp),INTENT(OUT) :: yout(:),fout(:)
  REAL(qp) :: ln1b,dx
  INTEGER :: i
  ln1b=LOG(1+beta)
  DO i=1,SIZE(xin)
   dx=xin(i)-x1
   yout(i)=y1+alpha*dx+LOG(EXP(dx)+beta)-ln1b
   fout(i)=1.0D0/(1+EXP(yout(i)-xin(i)))
   ! check if ad. speed of sound exceeds the speed of light
   IF(EXP(yout(i)-x1)*(alpha/EXP(dx)+1/(EXP(dx)+beta))<1) THEN 
    stp=-i
    EXIT
   ENDIF
  END DO
  stp=i
 END SUBROUTINE link_eos
 SUBROUTINE link_n(x1,y1,xxx,yyy,nn,zz)
  REAL(qp),INTENT(IN) :: xxx(:),yyy(:),x1,y1
  REAL(qp),INTENT(INOUT) :: nn,zz(:,:)
  REAL(qp) :: ii, ip, yyi, yyf
  INTEGER :: i
!  ii=1.0D0/(EXP(x1)+EXP(y1))
  ii=1.0D0/(1+EXP(x1-y1))
  yyf=y1
  zz(1,3)=EXP(nn)
  DO i=1,linkdots
   yyi=yyf
   yyf=yyy(i)
   ip=ii
!   ii=1.0D0/(EXP(xxx(i))+EXP(yyf))
!   nn=nn+(ii+ip)*(EXP(yyf)-EXP(yyi))/2
   ii=1.0D0/(1+EXP(xxx(i)-yyf))
   nn=nn+(ii+ip)*(yyf-yyi)/2
   zz(i+1,3)=EXP(nn)
  END DO
 END SUBROUTINE link_n

 SUBROUTINE center_TOV(pc,rhoc,fc,pout,rhoout,fout,mout,r2out,yLout,nn)
  REAL(qp),INTENT(IN) :: pc,rhoc,pout,fc,fout,rhoout
  REAL(qp),INTENT(OUT) :: mout,r2out,yLout,nn
  REAL(qp) :: dr2dx, dx, rout,dmdx
  dr2dx=-fc*3/(2*pi*(rhoc+3*pc))
  dx=LOG(pout/pc)
  r2out=dr2dx*dx/2
  rout=SQRT(r2out)
  mout=4*pi*rhoc*rout**3/3
  dr2dx=-(fout+fc)*r2out*(rout-2*mout)/(mout+2*pi*rout**3*(pout+pc))
  r2out=dr2dx*dx
  mout=4*pi*rhoc*SQRT(r2out)**3/3
  r2out=LOG(r2out/Rkm**2)
  mout=mout/Rkm
!  r2out=(pc-pout)*3/(2*pi*(pc+rhoc)*(3*pc+rhoc))
!  mout=4*pi*r2out**(3.0D0/2)/3
  IF(Love) THEN
   yLout=2-6*(5*rhoc+9*pc+(rhoc+pc)*(rhoout-rhoc)/(pout-pc))*(fc+fout)*LOG(pout/pc)/(14*(3*pc+rhoc))
  ELSE
   yLout=-1
  ENDIF
  IF(nDen) nn=(1.0D0/(rhoout+pout)+1.0D0/(rhoc+pc))*(rhoout-rhoc)/2
 END SUBROUTINE center_TOV

 SUBROUTINE calc_NS(knt,err)
  INTEGER,INTENT(IN) :: knt
  INTEGER,INTENT(OUT) :: err
  INTEGER :: i, stp, io
  REAL(qp) :: Mns
  mMaxAll=0
  err=100000
  io=MAX(2,2+linkpoints*(knt-2))
!$OMP PARALLEL DO PRIVATE(i,Mns,stp) REDUCTION(MIN:err)
  DO i=io,NNS
   CALL integ_TOV(i,Mns,stp)
!   WRITE(*,* ) i,': stp=',stp
   err=MIN(stp,err)
   IF(stp<=0.OR.(Mns.NE.Mns)) THEN
     Mns=-1
     NSchain(i,1)=-1
     NSchain(i,2)=-1
     IF(Love) NSchain(i,5)=-1
     CYCLE
   ELSE
    IF(MomI) THEN
     CALL calc_I(i)
     NSchain(i,6)=profile(1,5,i)/Rkm ! in M_Sun*km^2
    ENDIF
   ENDIF
!$OMP CRITICAL
   IF(Mns>mMaxAll.AND.Mns<=MR) THEN 
    mMaxAll=Mns
    strmax=i
   ENDIF
!$OMP END CRITICAL
  ENDDO
!$OMP END PARALLEL DO
 END SUBROUTINE calc_NS

 SUBROUTINE integ_TOV(str,mstr,stp)
  INTEGER,INTENT(IN) :: str
  INTEGER,INTENT(OUT) :: stp
  REAL(qp), INTENT(INOUT) :: mstr
  REAL(qp) :: pc,rc,fc,po,co,fo,r2str,ro,xc,yc, &
		xx(5*linkdots),yy(5*linkdots),ff(5*linkdots),dx,xo,yo,ddo,ddc, &
		alph, bet, yLstr,lam,nstr,nn,bstr
  INTEGER :: i,j
!  lnk=(str-2)/linkpoints
  xc=NSchain(str,3)
  yc=NSchain(str,4)
  ddc=yc-LOG(do+EXP(xc))
  pc=rhos*EXP(xc)
  rc=rhos*EXP(yc)
  fc=fchain(str)

  xo=NSchain(str-1,3)
  yo=NSchain(str-1,4)
  ddo=yo-LOG(do+EXP(xo))
  
  IF(xo>xc) THEN
   WRITE(*,*) ' str=',str
   STOP 'xo>xc'
  ENDIF
  bet=do/EXP(xc)
  alph=(ddo-ddc)/(xo-xc)
  dx=(xo-xc)/(5*linkdots)
  DO i=1,5*linkdots-1
   xx(i)=xc+dx*i
  ENDDO
  xx(5*linkdots)=xo
  CALL link_eos(xc,yc,xx,yy,ff,alph,bet,stp)
  IF(stp<0) STOP 'error with seed EOS'
  po=rhos*EXP(xx(1))
  ro=rhos*EXP(yy(1))
  fo=1.0D0/(1+ro/po)
  CALL center_TOV(pc,rc,fc,po,ro,fo,mstr,r2str,yLstr,nstr)
  bstr=0
  IF(BE) bstr=(Rkm**3*EXP(3*r2str/2)*EXP(nstr)/SQRT(1+2*mstr/EXP(r2str/2)))/2
  CALL link_TOV(mstr,r2str,REAL(xx(1),db),REAL(yy(1)-LOG(do+EXP(xx(1))),db),REAL(xo,db),REAL(ddo,db),stp,yLstr,nstr,bstr)
  IF(stp<=0) RETURN
  IF(MomI) THEN
   profile(iseed+str,1,str)=xc
   profile(iseed+str,2,str)=yc
   profile(iseed+str,3,str)=0
   profile(iseed+str,4,str)=0
   profile(iseed+str-1,1,str)=xo
   profile(iseed+str-1,2,str)=yo
   profile(iseed+str-1,3,str)=mstr*Rkm
   profile(iseed+str-1,4,str)=Rkm*EXP(r2str/2)
  ENDIF
  i=str-1
  DO WHILE(i>=2)
   CALL link_TOV(mstr,r2str,NSchain(i,3),NSchain(i,4)-LOG(do+EXP(NSchain(i,3))), &
		NSchain(i-1,3),NSchain(i-1,4)-LOG(do+EXP(NSchain(i-1,3))),stp,yLstr,nstr,bstr)
   IF(stp<=0) RETURN
   IF(MomI) THEN
    profile(iseed+i-1,1,str)=NSchain(i-1,3)
    profile(iseed+i-1,2,str)=NSchain(i-1,4)
    profile(iseed+i-1,3,str)=mstr*Rkm
    profile(iseed+i-1,4,str)=Rkm*EXP(r2str/2)
   ENDIF
   i=i-1
  ENDDO
  IF(NSchain(1,1)==0) THEN ! if it bare SQM star, and surface is reached
   NSchain(str,1)=mstr
   NSchain(str,2)=Rkm*EXP(r2str/2)
   IF(nDen) NSchain(str,7)=nseed*EXP(-nstr) ! dimensionless, multiply by 'n' at surface for units [nseed=1.0]
   IF(BE) NSchain(str,8)=mAc2*4*pi*nseed*EXP(-nstr)*bstr*toMsc2 ! baryon mass in M_Sun*fm^3
   IF(Love) THEN
    CALL calc_lambda(mstr,NSchain(str,2),yLstr,lam)
    NSchain(str,5)=lam
   ENDIF
   RETURN
  ENDIF
  ! join point of chain with seed EOS
  CALL link_TOV(mstr,r2str,NSchain(1,3),NSchain(1,4)-LOG(do+EXP(NSchain(1,3))), &
		xEOS(iseed),yEOS(iseed)-LOG(do+EXP(xEOS(iseed))),stp,yLstr,nstr,bstr)
  IF(MomI) THEN
    profile(iseed,3,str)=mstr*Rkm
    profile(iseed,4,str)=Rkm*EXP(r2str/2)
  ENDIF
  IF(stp<=0) RETURN
  ! seed EOS to the surface
  i=iseed
  nn=0
  IF(BE) nn=nstr
  DO WHILE(i>1)
   CALL link_TOV(mstr,r2str,xEOS(i),yEOS(i)-LOG(do+EXP(xEOS(i))),xEOS(i-1),yEOS(i-1)-LOG(do+EXP(xEOS(i-1))),stp,yLstr,nn,bstr)
   IF(stp<=0) RETURN
   IF(MomI) THEN
    profile(i-1,3,str)=mstr*Rkm
    profile(i-1,4,str)=Rkm*EXP(r2str/2)
   ENDIF
   i=i-1
  ENDDO
  NSchain(str,1)=mstr
  NSchain(str,2)=Rkm*EXP(r2str/2)
  IF(nDen) NSchain(str,7)=nEOS(iseed)*EXP(-nstr)
  IF(nDen) NSchain(str,8)=4*pi*mAc2*nEOS(iseed)*EXP(-nstr)*bstr*toMsc2-mstr
  IF(Love) THEN
   CALL calc_lambda(mstr,NSchain(str,2),yLstr,lam)
   IF(lam.NE.lam) THEN
     NSchain(str,5)=-1
   ELSE
     NSchain(str,5)=lam
   ENDIF
  ENDIF
 END SUBROUTINE integ_TOV
 
 SUBROUTINE init_MR(xss,mmstr,rr2str)
  REAL(db),INTENT(IN) :: xss
  REAL(db),INTENT(OUT) :: mmstr,rr2str
  REAL(qp) :: mstr,r2str
  REAL(qp) :: ps,rs,fs,po,ro,fo,ddo,xo,yo,ys,dseed,yL,lam,nstr,nn,bstr,xs
  REAL(qp) :: xx(5*linkdots),yy(5*linkdots),ff(5*linkdots),alph,bet,dx
  INTEGER :: i,stp,s,j
  
  xs=xss
  IF(xs==0.AND.MINknot==1) xs=xEOS(1)
  IF(xs>xEOS(1).AND.xs<=xEOS(dotsEOS)) THEN
   s=1
   DO WHILE(xEOS(s+1)<xs) 
    s=s+1
   ENDDO
   ys=yEOS(s)+(yEOS(s+1)-yEOS(s))*(xs-xEOS(s))/(xEOS(s+1)-xEOS(s))
   nseed=nEOS(s)+(nEOS(s+1)-nEOS(s))*(xs-xEOS(s))/(xEOS(s+1)-xEOS(s))
  ELSEIF(xs<xEOS(1).OR.xs>xEOS(dotsEOS)) THEN
   STOP 'Seed point xseed is out of bounds of EOS'
  ELSE
   xseed=xs
   ys=yEOS(1)
   nseed=nEOS(1)
   s=-1
  ENDIF
  do=EXP(ys)-EXP(xs)
  yseed=ys
  dseed=ys-LOG(do+EXP(xs))
  IF(s>0) THEN
   iseed=s
  ELSE
   iseed=0
  ENDIF
  IF(MomI) CALL init_profile
  ps=rhos*EXP(xs)
  rs=rhos*EXP(ys)
  fs=1.0D0/(1+rs/ps)
  IF(iseed==0) THEN
   WRITE(*,*) 'Seed density at a star surface is allowed to change'
   NSchain(1,1)=0
   NSchain(1,2)=0
   NSchain(1,3)=xs ! rhos*EXP(xseed)/EDtoKM in MeV/fm^3
   NSchain(1,4)=ys ! rhos*EXP(yseed)/EDtoKM in MeV/fm^3
   nseed=1.0
   NSchain(1,7)=nseed ! orbitrary at the surface, set by SQM surface physics
   fchain(1)=1/(1+EXP(ys-xs))
   WRITE(*,*) 'NSchain(1,:)=',NSchain(1,:)
   IF(Love) THEN
    NSchain(1,5)=0
    dyR=1 ! yR=2-3
   ENDIF
   IF(MomI) THEN
    profile(1,1,:)=xEOS(1)
    profile(1,2,:)=yEOS(1)
    profile(1,3,1)=0
    profile(1,4,1)=0
   ENDIF
   mmstr=0
   rr2str=0
   RETURN
  ELSE
   dyR=0
   IF(MomI) THEN
!$OMP PARALLEL DO PRIVATE(i)
    DO i=1,iseed
     profile(i,1,:)=xEOS(i)
     profile(i,2,:)=yEOS(i)
    ENDDO
!$OMP END PARALLEL DO
    profile(iseed+1,1,:)=xseed
    profile(iseed+1,2,:)=yseed
    profile(iseed+1,3,1)=0
    profile(iseed+1,4,1)=0
   ENDIF
  ENDIF
  xo=xEOS(s)
  yo=yEOS(s)
  ddo=yo-LOG(do+EXP(xo))

  bet=do/EXP(xs)
  alph=(ddo-dseed)/(xo-xs)
  dx=(xo-xs)/(5*linkdots)
  DO i=1,5*linkdots-1
   xx(i)=xs+dx*i
  ENDDO
  xx(5*linkdots)=xo
  CALL link_eos(xs,ys,xx,yy,ff,alph,bet,stp)
  IF(stp<0) STOP 'error with seed EOS'
  po=rhos*EXP(xx(1))
  ro=rhos*EXP(yy(1))
  fo=1.0D0/(1+ro/po)
  CALL center_TOV(ps,rs,fs,po,ro,fo,mstr,r2str,yL,nstr)
  bstr=0
  IF(BE) bstr=Rkm**3*EXP(nstr)*EXP(3*r2str/2)/(2*SQRT(1+2*mstr/EXP(r2str/2)))
  CALL link_TOV(mstr,r2str,REAL(xx(1),db),REAL(yy(1)-LOG(do+EXP(xx(1))),db),REAL(xo,db),REAL(ddo,db),stp,yL,nstr,bstr)
  IF(stp<0) STOP 'error with seed EOS:link_TOV'
  IF(MomI) THEN
   profile(s,3,1)=mstr*Rkm
   profile(s,4,1)=Rkm*EXP(r2str/2)
  ENDIF
  nn=0
  IF(BE) nn=nstr
  DO i=s,2,-1
   CALL link_TOV(mstr,r2str,xEOS(i),yEOS(i)-LOG(do+EXP(xEOS(i))),xEOS(i-1),yEOS(i-1)-LOG(do+EXP(xEOS(i-1))),stp,yL,nn,bstr)
   IF(MomI) THEN
    profile(i-1,3,1)=mstr*Rkm
    profile(i-1,4,1)=Rkm*EXP(r2str/2)
   ENDIF
  ENDDO
  NSchain(1,1)=mstr ! in M_Sun
  NSchain(1,2)=Rkm*EXP(r2str/2) ! in km
  NSchain(1,3)=xs ! rhos*EXP(xseed)/EDtoKM ! in MeV/fm^3
  NSchain(1,4)=ys ! rhos*EXP(yseed)/EDtoKM ! in MeV/fm^3
  IF(nDen) NSchain(1,7)=nEOS(s)*EXP(-nstr)
  IF(BE) NSchain(1,8)=nEOS(s)*EXP(-nstr)*4*pi*mAc2*bstr*toMsc2-mstr
  IF(MomI) THEN
   CALL calc_I(1)
   NSchain(1,6)=profile(1,5,1)/Rkm ! in M_Sun*km^2
  ENDIF
  IF(Love) THEN
   IF(dyR/=0) yL=yL-dyR*rss*EXP(xEOS(1))*4*pi*EXP(3*r2str/2.0D0)/mstr
   CALL calc_lambda(mstr*Rkm,NSchain(1,2),yL,lam)
   IF(lam.NE.lam) THEN
     NSchain(1,5)=-1
   ELSE
     NSchain(1,5)=lam
   ENDIF
  ENDIF
  fchain(1)=1/(1+EXP(ys-xs))
  WRITE(*,*) 'NSchain(1,:)=',NSchain(1,:)
  mmstr=mstr
  rr2str=Rkm**2*EXP(r2str)
 END SUBROUTINE init_MR
 SUBROUTINE calc_lambda(mm,rr,yR,ll)
  REAL(db),INTENT(IN) :: rr
  REAL(qp),INTENT(IN) :: mm,yR
  REAL(qp),INTENT(OUT) :: ll
  REAL(qp) :: bt,k2
  bt=mm/rr
!  IF(bt<0.01) THEN
!   k2=(2-yR)/(3+yR)
!   k2=k2+(yR**2-6*yR-6)*bt/(yR+3)**2
!   k2=k2+(yR**3+34*yR**2-8*yR+12)*bt**2/(7*(yR+3)**3)
!   k2=k2+(yR**4+62*yR**3+84*yR**2+48*yR+36)*bt**4/(7*(yR+3)**4)
!   k2=k2+5*(5*yR**5+490*yR**4+1272*yR**3+1884*yR**2+1476*yR+648)*bt**5/(294*(yR+3)**5)
!   k2=k2*(1-2*bt)**2/2
!  ELSE
   k2=2*bt*(6-3*yR+3*bt*(5*yR-8))
   k2=k2+4*bt**3*(13-11*yR+bt*(3*yR-2)+2*bt**2*(1+yR))
   k2=k2+3*(1-2*bt)**2*(2-yR+2*bt*(yR-1))*LOG(1-2*bt)
   k2=8*bt**5*(1-2*bt)**2*(2+2*bt*(yR-1)-yR)/(5*k2)
!  ENDIF
  ll=k2*2*rr**5/3.0D0
 END SUBROUTINE calc_lambda
 SUBROUTINE calc_I(str)
  INTEGER, INTENT(IN) :: str
  INTEGER :: i
  REAL(qp) :: coefI,wpp,wp,wr,fp,fr,ntp,ntr,ir,ip,jp,jr,sq,dr,drp
  coefI=8*rhos*pi/3.0D0
  i=iseed+str
  wp=1
  profile(i,5,str)=0.0D0
  fp=1.0D0/(1+EXP(profile(i,2,str)-profile(i,1,str)))
  jp=1
  ntp=0
  i=i-1
  IF(i==0) RETURN
  sq=SQRT(1-2*profile(i,3,str)/profile(i,4,str))
  fr=1.0D0/(1+EXP(profile(i,2,str)-profile(i,1,str)))
  ntr=ntp+(profile(i,1,str)-profile(i+1,1,str))*(fp+fr)/2
  jr=EXP(-ntr)*sq
  wr=5*jp*wp/(4*jr+jp)
  ip=0
  ir=coefI*profile(i,4,str)**4*EXP(-ntr)*wr*(EXP(profile(i,1,str))+EXP(profile(i,2,str)))/sq
  profile(i,5,str)=profile(i+1,5,str)+(profile(i,4,str)-profile(i+1,4,str))*(ip+ir)/2
  dr=profile(i,4,str)-profile(i+1,4,str)
  i=i-1
  DO WHILE(i>=1)
   wpp=wp
   wp=wr
   sq=SQRT(1-2*profile(i,3,str)/profile(i,4,str))
   fp=fr
   fr=1.0D0/(1+EXP(profile(i,2,str)-profile(i,1,str)))
   ntp=ntr
   ntr=ntp+(profile(i,1,str)-profile(i+1,1,str))*(fp+fr)/2
   jp=jr
   jr=EXP(-ntr)*sq
   drp=dr
   dr=profile(i,4,str)-profile(i+1,4,str)
   wr=((dr+drp)*(dr*(dr+9*drp)*jp-dr*(3*dr+drp)*jr+2*(drp*(-jp+jr)+dr*(jp+jr))*profile(i,4,str))*wp &
   -dr**2*(dr*(jp-3*jr)+2*(jp+jr)*profile(i,4,str))*wpp) &
    /(dr*drp*(2*dr*(jp+2*jr)+drp*(jp+7*jr))+2*drp*(2*dr*jr+drp*(-jp+jr))*profile(i,4,str))
   ip=ir
   ir=coefI*profile(i,4,str)**4*EXP(-ntr)*wr*(EXP(profile(i,1,str))+EXP(profile(i,2,str)))/sq
   profile(i,5,str)=profile(i+1,5,str)+(profile(i,4,str)-profile(i+1,4,str))*(ip+ir)/2
   i=i-1
  ENDDO
  coefI=wr+profile(1,4,str)*(dr**2*(-wp+wpp)+2*dr*drp*(-wp+wr)+drp**2*(-wp+wr))/(dr*drp*(dr+drp)*3)
  coefI=coefI*sq*EXP(-ntr)
  DO i=iseed+str,1
   profile(i,5,str)=profile(i,5,str)/coefI
  ENDDO
 END SUBROUTINE calc_I
END MODULE Solver