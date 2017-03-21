! added in ARTOS ver. 0.2 . 2012_07_24 by SCB  
	subroutine onetpenmg_senm_iter(          &
		ng,reigv,sigr,signf,sigs,dc,chi      &
	   ,pflx,cnti,sz,srcbal,srcmomx,srcmomy  &
	   ,cnto,hflx,errh                       &
	   ,aflx,xmom,ymom,cpfl,sfli             &
	   )

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
	dimension xmom(NG,myn),ymom(NG,myn),sfli(NG,mxs),cpfl(NG)  
	dimension asrc(NG,myn),xsrc(NG,myn),ysrc(NG,myn),ssrc(NG,mxs)
	dimension ojsrc(NG,mxo),psrc(NG)
	dimension asrcn(NG,myn),xsrcn(NG,myn),ysrcn(NG,myn),ssrcn(NG,mxs)
	dimension ojsrcn(NG,mxo),psrcn(NG)

	dimension sigr(NG),signf(NG),dc(NG),sigs(NG,NG),chi(NG)
	dimension ca(NG,3),cx(NG,3),cy(NG),cs(NG,5),cj(NG,4),cp(NG,3),caxy(NG,NG)
	dimension csp(NG,NG,4),cjp(NG,NG,3),cpp(NG,NG,3)
	dimension rcaxy(NG,NG),rcjp(NG,NG),rcsp(NG,NG),rcpp(NG,NG)
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
	real, save :: total_time=0.
	logical, save :: first=.TRUE.


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
! Gauss-Seidel Iteration
	do iter=1,100	
	
	do ig=1,ng
	psrcn(ig)=psrc(ig)
	enddo
	do it=1,6
		do ig=1,ng
			asrcn(ig,it)=asrc(ig,it)-ca(ig,1)*(sfli(ig,it)+sfli(ig,mp1(it))) &
									-ca(ig,2)*cnto(ig,it)-ca(ig,3)*cpfl(ig) 
			xsrcn(ig,it)=xsrc(ig,it)-cx(ig,1)*(sfli(ig,it)+sfli(ig,mp1(it))) &
									-cx(ig,2)*cnto(ig,it)-cx(ig,3)*cpfl(ig)
			ysrcn(ig,it)=ysrc(ig,it)+cy(ig)*(sfli(ig,it)-sfli(ig,mp1(it)))
		enddo
#ifdef DEBUG
		asrcn(1,it)=asrcn(1,it)-caxy(1,2)*aflx(2,it)
		asrcn(2,it)=asrcn(2,it)-caxy(2,1)*aflx(1,it)
		
		aflx(1,it)=asrcn(1,it)/caxy(1,1)
		aflx(2,it)=asrcn(2,it)/caxy(2,2)

		xsrcn(1,it)=xsrcn(1,it)-caxy(1,2)*xmom(2,it)
		xsrcn(2,it)=xsrcn(2,it)-caxy(2,1)*xmom(1,it)

		xmom(1,it)=xsrcn(1,it)/caxy(1,1)
		xmom(2,it)=xsrcn(2,it)/caxy(2,2)

		ysrcn(1,it)=ysrcn(1,it)-caxy(1,2)*ymom(2,it)
		ysrcn(2,it)=ysrcn(2,it)-caxy(2,1)*ymom(1,it)

		ymom(1,it)=ysrcn(1,it)/caxy(1,1)
		ymom(2,it)=ysrcn(2,it)/caxy(2,2)

#else
		call a2xb1_c1(ng,rcaxy,asrcn(:,it),aflx(:,it))
		call a2xb1_c1(ng,rcaxy,xsrcn(:,it),xmom(:,it))
		call a2xb1_c1(ng,rcaxy,ysrcn(:,it),ymom(:,it))
#endif
	enddo

	do it=1,6
		do ig=1,ng
			ssrcn(ig,it)=ssrc(ig,it)-cs(ig,1)*(aflx(ig,it)+aflx(ig,mp5(it)))          &
			                        -cs(ig,2)*(xmom(ig,it)+xmom(ig,mp5(it)))          &
			                        -cs(ig,2)*ymom(ig,it)-cs(ig,3)*ymom(ig,mp5(it))   &
									-cs(ig,5)*cpfl(ig)
			sfli(ig,it)=ssrcn(ig,it)/cs(ig,4)

			ojsrcn(ig,it)=ojsrc(ig,it)-cj(ig,1)*aflx(ig,it)-cj(ig,2)*xmom(ig,it)      &
			                          -cj(ig,4)*cpfl(ig)
			cnto(ig,it)=ojsrcn(ig,it)/cj(ig,3)

			psrcn(ig)=psrcn(ig)-cp(ig,1)*xmom(ig,it)-cp(ig,2)*cnto(ig,it)			
		enddo		
	enddo
	
	do ig=1,ng			
		cpfl(ig)=psrcn(ig)/cp(ig,3)
	enddo

	errh=0.
	do ig=1,ng
	do it=1,myn
	err=(oflx(ig,it)-aflx(ig,it))/aflx(ig,it)
	if(err.lt.0.) err=-err
	if(err.gt.errh) errh=err 
	enddo
	enddo

	oflx=aflx
	
	if(err.lt.1.e-6) exit
	enddo

	do ig=1,ng
		hflx(ig)=0.
		do it=1,6 
		hflx(ig)=hflx(ig)+aflx(ig,it)
		enddo
		hflx(ig)=hflx(ig)/6.
	enddo


!	print *, iter

	return
	end subroutine


