! added in ARTOS ver. 0.2 . 2012_07_24 by SCB  
	subroutine onetpenmg_senm(               &
		ng,reigv,sigr,signf,sigs,dc,chi      &
	   ,pflx,cnti,sz,srcbal,srcmomx,srcmomy  &
	   ,cnto,hflx,errh                       &
	   ,sfli)

	implicit double precision (a-h,o-z)

	! HOPEN(Higher Order Polynomial Expansion Nodal method) solver for hexagon 
	!    using HOPEN solution based on Triangle 
	!    input : boundary incoming current, boundary point flux
	!    output : boundary outgoing current  
	!  
	!    pflx : boundary point flux(input) 
	!    cnti : boundary incoming current(input) 
	!    aflx : node average flux at previous step 
	!    xmom : x-moment at previous step 
	!    ymom : y-moment at previous step 
	!    sfli : inner surface flux
	!    cpfl : center point flux
	!    cnto : boundary outgoing current(output) 
	!
	!                        2     
	!                   2---------3 
	!                   / \  2  / \ 
	!                 1/   2   3   \3 
	!                 /  1  \ /   3 \  
	!                1--1--------4---4 
	!                 \   6 / \   4 / 
	!                 6\   6   5   /4  
	!                   \ /  5  \ /   
	!                   6---------5 
	!                        5             
	!
	!
	parameter (myn=6,mxs=6,mxo=6,mxp=1)

	dimension pflx(NG,mxo),cnti(NG,mxo),cnto(NG,mxo),aflx(NG,myn),hflx(NG)
	dimension xmom(NG,myn),ymom(NG,myn),sfli(NG,mxs),cpfl(NG),osfli(NG,mxs) 
	dimension asrc(NG,myn),xsrc(NG,myn),ysrc(NG,myn),ssrc(NG,mxs),ssrcn(NG,mxs)
	dimension ojsrc(NG,mxo),psrc(NG)
	dimension sigr(NG),signf(NG),dc(NG),sigs(NG,NG),chi(NG)
	dimension ca(NG,3),cx(NG,3),cy(NG),cs(NG,5),cj(NG,4),cp(NG,3),caxy(NG,NG)
	dimension csp(NG,NG,4),cjp(NG,NG,3),cpp(NG,NG,3)
	dimension rcaxy(NG,NG),rcjp(NG,NG),rcsp(NG,NG),rcpp(NG,NG),rcsp1(NG,NG)
	dimension oflx(NG,myn),fac(NG,NG),fac1(NG,NG),fac2(NG,NG),fac3(NG,NG),ttt(NG)
	dimension dat1(NG,NG), dat2(NG,NG)
	dimension srcbal(NG,myn),srcmomx(NG,myn),srcmomy(NG,myn)

	dimension fs(NG,NG), fj(NG,NG), fp(NG,NG), fa(NG,NG)
	dimension fs3j1(NG,NG), fs3j3(NG,NG), fp2j1(NG,NG), fp2j3(NG,NG)
	dimension dum2(6*NG,6*NG), dum1(NG)

	dimension mp1(6),mp2(6),mp3(6),mp4(6),mp5(6)
	dimension mg1(4),mg2(4),mg3(4),mg4(4),mg5(4) 
	data mp1/2,3,4,5,6,1/ 
	data mp2/3,4,5,6,1,2/ 
	data mp3/4,5,6,1,2,3/ 
	data mp4/5,6,1,2,3,4/ 
	data mp5/6,1,2,3,4,5/
	data mg1/1,1,2,2/
	data mg2/1,2,1,2/
	data mg3/1,1,3,3/
	data mg4/2,2,4,4/
	data mg5/3,4,3,4/

	real matA(6*NG,6*NG), rmatA(6*NG,6*NG)
	real,save :: total_time=0.


	rt3=sqrt(3.d0)
	!
#ifdef DEBUG
	do it=1,6
	srcbal(:,it)=1.*it
	srcmomx(:,it)=2.*it
	srcmomy(:,it)=3.*it
	enddo
	pflx=1.5
	cnti=0.5
	reigv=1.
#endif

	do it=1,myn 
	do ig=1,ng
	asrc(ig,it)=srcbal(ig,it)
	xsrc(ig,it)=srcmomx(ig,it)
	ysrc(ig,it)=srcmomy(ig,it)
	enddo
	enddo


	caxy=0.
	do ig=1,ng
		do ig2=1,ng
			if(ig2.ne.ig) then
			caxy(ig,ig2)=caxy(ig,ig2)-sigs(ig2,ig)
			endif
			caxy(ig,ig2)=caxy(ig,ig2)-signf(ig2)*reigv*chi(ig)
		enddo

		bt=dc(ig)/sz/sz  
		caxy(ig,ig)=caxy(ig,ig)+80.*bt+sigr(ig)
		ca(ig,1)=-32.*bt
		ca(ig,2)=-64.*bt
		ca(ig,3)=16.*bt/3.
		cx(ig,1)=8.*bt/3.
		cx(ig,2)=-32.*bt/3.
		cx(ig,3)=16.*bt/9.
		cy(ig)=-8.*bt
		cs(ig,1)=20. 
		cs(ig,2)=-60.
		cs(ig,3)=60.
		cs(ig,4)=-48.
		cs(ig,5)=2.
		tt=dc(ig)/rt3/sz
		cj(ig,1)=20.*tt
		cj(ig,2)=120.*tt
		cj(ig,3)=-1.-48.*tt
		cj(ig,4)=2.*tt
		cp(ig,1)=-15.
		cp(ig,2)=2.
		cp(ig,3)=-6.

		psrc(ig)=0.
		do it=1,mxo
			asrc(ig,it)=asrc(ig,it)-ca(ig,2)*cnti(ig,it) &
				-ca(ig,3)*(pflx(ig,it)+pflx(ig,mp1(it)))
			xsrc(ig,it)=xsrc(ig,it)-cx(ig,2)*cnti(ig,it) &
				+0.5*cx(ig,3)*(pflx(ig,it)+pflx(ig,mp1(it)))
			ysrc(ig,it)=ysrc(ig,it) &
				+8.*bt/3.*(pflx(ig,it)-pflx(ig,mp1(it)))
			ssrc(ig,it)=-2.*(pflx(ig,it) &
				+pflx(ig,mp1(it))+pflx(ig,mp5(it)))
			ojsrc(ig,it)=(-1.+48.*tt)*cnti(ig,it) &
				-tt*(pflx(ig,it)+pflx(ig,mp1(it)))
			psrc(ig)=psrc(ig)-2.*cnti(ig,it)
			oflx(ig,it)=aflx(ig,it)
		enddo
	enddo

	call invag(ng,caxy,rcaxy)
	
	r3=1./3.
	r6=1./6.
	do j=1,ng
		caj1=ca(j,1)
		do i=1,ng
			temp=rcaxy(i,j)
! Left hand side
			cjr=cj(i,1)*temp
			csr=cs(i,1)*temp
			cpr=cp(i,1)*temp

			cjra=cjr*caj1
			csra=csr*caj1
			cpra=cpr*caj1

			cjp(i,j,1)=-0.5*cjra
			cjp(i,j,2)=-4.*cjra
			cjp(i,j,3)=0.5*cjra

			csp(i,j,1)=-4.*csra
			csp(i,j,2)=-0.5*csra
			csp(i,j,3)=-csra
			csp(i,j,4)=0.
			
			cpp(i,j,1)=r6*cpra
			cpp(i,j,2)=-r3*cpra
			cpp(i,j,3)=r3*cpra

! Right hand side 
			do it=1,6
				ojsrc(i,it)=ojsrc(i,it)+cjr*(-asrc(j,it)-6.*xsrc(j,it))	

				ssrc(i,it)=ssrc(i,it)+csr*(    &
					-asrc(j,it)-asrc(j,mp5(it)) &
					+3.*(xsrc(j,it)+xsrc(j,mp5(it))-ysrc(j,it)+ysrc(j,mp5(it))))

	
				psrc(i)=psrc(i)-cpr*xsrc(j,it)
			enddo	
		enddo

! Left hand side
		cjp(j,j,2)=cjp(j,j,2)+cj(j,3)
		cjp(j,j,3)=cjp(j,j,3)+0.1*cj(j,1)

		csp(j,j,1)=csp(j,j,1)-2.4*cs(j,1)
		csp(j,j,4)=0.1*cs(j,1)

		cpp(j,j,2)=cpp(j,j,2)+cp(j,2)
		cpp(j,j,3)=cpp(j,j,3)+cp(j,3)
	enddo

! Condensed Matrix: 13 x 13


	call invag(ng,cjp(:,:,2),rcjp)

! Left hand side
	do j=1,ng
	do i=1,ng
		fs3j1(i,j)=0.
		fs3j3(i,j)=0.
		fp2j1(i,j)=0.
		fp2j3(i,j)=0.
		do n=1,ng
		do m=1,ng
			fs3j1(i,j)=fs3j1(i,j)+csp(i,m,3)*rcjp(m,n)*cjp(n,j,1)
			fs3j3(i,j)=fs3j3(i,j)+csp(i,m,3)*rcjp(m,n)*cjp(n,j,3)
			fp2j1(i,j)=fp2j1(i,j)+cpp(i,m,2)*rcjp(m,n)*cjp(n,j,1)
			fp2j3(i,j)=fp2j3(i,j)+cpp(i,m,2)*rcjp(m,n)*cjp(n,j,3)
		enddo
		enddo
		csp(i,j,1)=csp(i,j,1)-2.*fs3j1(i,j)
		csp(i,j,2)=csp(i,j,2)-fs3j1(i,j)

		csp(i,j,4)=csp(i,j,4)-2.*fs3j3(i,j)

		cpp(i,j,1)=cpp(i,j,1)-2.*fp2j1(i,j)
		cpp(i,j,3)=cpp(i,j,3)-6.*fp2j3(i,j)

! Right hand side
		do it=1,6
		do k=1,ng

		ssrc(i,it)=ssrc(i,it)-csp(i,k,3)*rcjp(k,j)*(ojsrc(j,it)+ojsrc(j,mp5(it)))

		psrc(i)=psrc(i)-cpp(i,k,2)*rcjp(k,j)*ojsrc(j,it)
		enddo
		enddo
	enddo
	enddo

! Condensed Matrix: 7 x 7
	do j=1,ng
	do i=1,ng
	fac(i,j)=csp(i,j,1)+2.*csp(i,j,2)
	enddo
	enddo

	call invag(ng,fac,rcsp)

	do j=1,ng
	do i=1,ng
		fac1(i,j)=0.
		do n=1,ng
			do m=1,ng
				fac1(i,j)=fac1(i,j)+6.*cpp(i,m,1)*rcsp(m,n)*csp(n,j,4)
			enddo
		enddo
		cpp(i,j,3)=cpp(i,j,3)-fac1(i,j)
	enddo
	enddo

	ttt=0.
	do j=1,ng
	do it=1,6
	ttt(j)=ttt(j)+ssrc(j,it)
	enddo
	enddo

	do j=1,ng
	do i=1,ng
		do n=1,ng
			psrc(i)=psrc(i)-cpp(i,n,1)*rcsp(n,j)*ttt(j)
		enddo		
	enddo
	enddo

	call invag(ng,cpp(:,:,3),rcpp)
	call a2xb1_c1(ng,rcpp,psrc,cpfl)

	do it=1,6
	do j=1,ng
	do i=1,ng
	ssrc(i,it)=ssrc(i,it)-csp(i,j,4)*cpfl(j)
	ojsrc(i,it)=ojsrc(i,it)-cjp(i,j,3)*cpfl(j)
	enddo
	enddo
	enddo

! Condensed Matrix: 6 x 6
#ifdef DEBUG
	call cpu_time(time0)
	matA=0.
	do j=1,6
	js=(j-1)*ng+1
	je=j*ng
	matA(js:je,js:je)=csp(:,:,1)
		if(j.ne.6) then
		is=(j-1)*ng+1
		ie=j*ng
		jps=j*ng+1
		jpe=(j+1)*ng


		matA(is:ie,jps:jpe)=csp(:,:,2)
		matA(jps:jpe,is:ie)=csp(:,:,2)		
		endif
	enddo

	matA(1:ng,5*ng+1:6*ng)=csp(:,:,2)
	matA(5*ng+1:6*ng,1:ng)=csp(:,:,2)


	call invag(6*ng,matA,rmatA)
	call cpu_time(time1)
	total_time=total_time+time1-time0
!	print *, 'INVERSE',total_time

	sfli=0.
	do j=1,ng
	do jt=1,6
	jm=(jt-1)*ng+j
		do i=1,ng
		do it=1,6
		im=(it-1)*ng+i
		sfli(i,it)=sfli(i,it)+rmatA(im,jm)*ssrc(j,jt)		
		enddo
		enddo
	enddo
	enddo
#else

	call invag(ng,csp(:,:,1),rcsp1)
	osfli=0.

	do iter=1,10
	ssrcn=ssrc
	do it=1,6
		do j=1,ng
			do i=1,ng
				ssrcn(i,it)=ssrcn(i,it)-csp(i,j,2)*(sfli(j,mp1(it))+sfli(j,mp5(it)))
			enddo
			sfli(j,it)=0.
		enddo

		do j=1,ng
			do i=1,ng
				sfli(i,it)=sfli(i,it)+rcsp1(i,j)*ssrcn(j,it)
			enddo
		enddo
	enddo

	err2=0.
	serr2=0.
	do it=1,6
	do m=1,ng
		err=abs(osfli(m,it)-sfli(m,it))
		err2=err2+err*err
		serr2=serr2+sfli(m,it)*sfli(m,it)	
	enddo
	enddo	
	residual=sqrt(err2/serr2)
	if(residual.lt.1.e-4) exit

	osfli=sfli
	enddo
#endif

! Backward Substation
#ifndef DEBUG
	do it=1,6
		do j=1,ng
		do i=1,ng
		ojsrc(i,it)=ojsrc(i,it)-cjp(i,j,1)*(sfli(j,it)+sfli(j,mp1(it)))
		enddo
		enddo
		
		cnto(:,it)=0.
		do j=1,ng
		do i=1,ng
		cnto(i,it)=cnto(i,it)+rcjp(i,j)*ojsrc(j,it)
		enddo
		enddo

		do i=1,ng
		ysrc(i,it)=ysrc(i,it)-cy(i)*(sfli(i,it)-sfli(i,mp1(it)))
		xsrc(i,it)=xsrc(i,it)-cx(i,1)*(sfli(i,it)+sfli(i,mp1(it))) &
							 -cx(i,2)*cnto(i,it)                   &
							 -cx(i,3)*cpfl(i)
		asrc(i,it)=asrc(i,it)-ca(i,1)*(sfli(i,it)+sfli(i,mp1(it))) &
							 -ca(i,2)*cnto(i,it)                   &
							 -ca(i,3)*cpfl(i)
		enddo

		aflx(:,it)=0.
		xmom(:,it)=0.
		ymom(:,it)=0.
		do j=1,ng
		do i=1,ng
		ymom(i,it)=ymom(i,it)+rcaxy(i,j)*ysrc(j,it)
		xmom(i,it)=xmom(i,it)+rcaxy(i,j)*xsrc(j,it)
		aflx(i,it)=aflx(i,it)+rcaxy(i,j)*asrc(j,it)
		enddo
		enddo
	enddo

#else
	do it=1,6
	do i=1,ng
	ysrc(i,it)=ysrc(i,it)-cy(i)*(sfli(i,it)-sfli(i,mp1(it)))
	xsrc(i,it)=xsrc(i,it)-cx(i,1)*(sfli(i,it)+sfli(i,mp1(it))) &
						 -cx(i,2)*cnto(i,it)                   &
						 -cx(i,3)*cpfl(i)
	asrc(i,it)=asrc(i,it)-ca(i,1)*(sfli(i,it)+sfli(i,mp1(it))) &
						 -ca(i,2)*cnto(i,it)                   &
						 -ca(i,3)*cpfl(i)
	enddo
	enddo

	ymom=0.
	xmom=0.
	aflx=0.
	do it=1,6
	do j=1,ng
	do i=1,ng
	ymom(i,it)=ymom(i,it)+rcaxy(i,j)*ysrc(j,it)
	xmom(i,it)=xmom(i,it)+rcaxy(i,j)*xsrc(j,it)
	aflx(i,it)=aflx(i,it)+rcaxy(i,j)*asrc(j,it)
	enddo
	enddo
	enddo
#endif

	errh=0.
	do ig=1,ng
	do it=1,myn
	err=(oflx(ig,it)-aflx(ig,it))/aflx(ig,it)
	if(err.lt.0.) err=-err
	if(err.gt.errh) errh=err 
	enddo
	enddo

	do ig=1,ng
	hflx(ig)=0
	do it=1,6 
	hflx(ig)=hflx(ig)+aflx(ig,it)
	enddo
	hflx(ig)=hflx(ig)/6.
	enddo

	return
	end subroutine


	subroutine a2xb1(n,a,b,c)

	implicit double precision (a-h,o-z)

	dimension a(n,n), b(n), c(n,n)

	do j=1,n
	do i=1,n
	c(i,j)=a(i,j)*b(j)
	enddo
	enddo

	end subroutine


	subroutine a1xb2(n,a,b,c)

	implicit double precision (a-h,o-z)

	dimension a(n), b(n,n), c(n,n)

	do j=1,n
	do i=1,n
	c(i,j)=a(i)*b(i,j)
	enddo
	enddo

	end subroutine


	subroutine a2xb1_c1(n,a,b,c)

	implicit double precision (a-h,o-z)

	dimension a(n,n), b(n), c(n)

	c=0.
	do j=1,n
	do i=1,n
	c(i)=c(i)+a(i,j)*b(j)
	enddo
	enddo

	end subroutine

	subroutine a2xb2(n,a,b,c)

	implicit double precision (a-h,o-z)

	dimension a(n,n), b(n,n), c(n,n)

	c=0.
	do k=1,n
	do j=1,n
	do i=1,n
		c(i,j)=c(i,j)+a(i,k)*b(k,j)
	enddo
	enddo
	enddo

	end subroutine

	subroutine a2pb2(n,a,b,c)

	implicit double precision (a-h,o-z)

	dimension a(n,n), b(n,n), c(n,n)

	do j=1,n
	do i=1,n
		c(i,j)=a(i,j)+b(i,j)
	enddo
	enddo

	end subroutine

	subroutine a1pb2(n,a,b,c)

	implicit double precision (a-h,o-z)

	dimension a(n),b(n,n), c(n,n)

	do j=1,n
		do i=1,n
			c(i,j)=b(i,j)
		enddo
		c(j,j)=c(j,j)+a(j)	
	enddo

	end subroutine

	subroutine inv_mat(n,a,b)

	implicit double precision (a-h,o-z)

	dimension a(n,n), b(n,n)

	b=0.

	do ip=1,n
	do k=ip+1,n
		if(abs(a(k,ip)).gt.abs(a(ip,ip))) then
			a(ip,k)=a(k,ip)
		endif
		a(ip,:)=a(ip,:)/a(ip,ip)
		do i=1,n
			if(i.ne.ip) then
			b(i,:)=b(i,:)-a(i,ip)*a(ip,:)
			endif
		enddo
	enddo
	enddo

	return
	end subroutine
