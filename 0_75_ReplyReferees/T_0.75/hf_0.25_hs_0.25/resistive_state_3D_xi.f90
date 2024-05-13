program resistive_state_xi

	implicit real(8) (a-h,o-z)
	implicit integer(4) (i-n)

	parameter ( iteMax = 5000000 )

	complex(8), allocatable :: G(:,:,:)             ! order parameter
	complex(8), allocatable :: F(:,:,:)
	complex(8), allocatable :: Ux(:,:,:)			! link variables
	complex(8), allocatable :: Uy(:,:,:)
	complex(8), allocatable :: Uz(:,:,:)
	real(8)   , allocatable :: hx(:,:,:)            ! local magnetic field
	real(8)   , allocatable :: hy(:,:,:)			
	real(8)   , allocatable :: hz(:,:,:)
	real(8)   , allocatable :: fieldOld(:,:,:)
	real(8)   , allocatable :: hty(:,:)
	real(8)   , allocatable :: htz(:,:)
	real(8)   , allocatable :: thetax(:,:,:)        ! vector potential
	real(8)   , allocatable :: phix(:,:,:)
	real(8)   , allocatable :: thetay(:,:,:)
	real(8)   , allocatable :: phiy(:,:,:)
	real(8)   , allocatable :: thetaz(:,:,:)
	real(8)   , allocatable :: phiz(:,:,:)
	real(8)   , allocatable :: Qx(:,:,:)            ! current density
	real(8)   , allocatable :: Qy(:,:,:)
	real(8)   , allocatable :: Qz(:,:,:)
	real(8)   , allocatable :: vtc(:,:,:)           ! vorticity
	real(8)   , allocatable :: phase(:,:,:)         ! phase of the order parameter
	real(8)   , allocatable :: aJx(:,:,:)           ! current density
	real(8)   , allocatable :: aJy(:,:,:)
	real(8)   , allocatable :: aJz(:,:,:)
	real(8)   , allocatable :: aJxn(:,:)
	real(8)   , allocatable :: densityOld(:,:,:)    ! Copper pair density
	real(8)   , allocatable :: densityNew(:,:,:)
	real(8)   , allocatable :: densityNew2(:,:,:)
	real(8)   , allocatable :: Vt(:,:)				! induced voltage
	real(8)   , allocatable :: aVt(:)				! induced voltage
	complex(8), allocatable :: delP(:,:,:)          ! Delta
	real(8),    allocatable :: Pot(:,:,:)			!scalar potential
	real(8),    allocatable :: Qot(:,:,:)           !scalar potential
	real(8),    allocatable :: PotOld(:,:,:)        !scalar potential
	real(8),    allocatable :: rho(:,:,:)			!divergence of the current density

	real(8)   , allocatable :: hyn(:,:,:)			! local magnetic field at the normal state
	real(8)   , allocatable :: hzn(:,:,:)
	
	real(8)   , allocatable :: aa(:,:,:)            ! auxiliaries 
	real(8)   , allocatable :: TT(:,:,:)            ! local temperature
	real(8)   , allocatable :: SS(:,:,:)
	real(8)   , allocatable :: Tyz(:,:,:)
	real(8)   , allocatable :: W(:,:,:)             ! Power
	real(8)   , allocatable :: ttx(:)
	real(8)   , allocatable :: Vortice_Value(:)
	
	complex(8) cp            
    	
	character(len=20) pxyw,pxzw,pyzw,phxyw,phxzw,phyzw,hxw,hyw,hzw,vw,vtcyw,vtczw,jxw,jxbw,jyw,jzw,htyw,htzw,jnxw,jnyw,txyw,tyzw,txw,potw,txysw,txyiw,txymw
	character(len=10) aJastring,CiCountT,CiCountF
	character(len=20) arquiv

	logical arquiv_gl,arquiv_re,arquiv_ET,arquiv_Iter,arquiv_TV,arquiv_IV,arquiv_prmt
	logical arquiv_CVT,arquiv_rg,bVolt,cframe,arquiv_period,arquiv_Q,arquiv_Qm
	logical arquiv_Tg,arquiv_Tgm, arquiv_Vel

	call cpu_time(time_initial)

	pi  = dacos(-1.0d0) ! pi number
	dpi = 2.0d0*pi

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Opening input and output files !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	open(unit=10,file='gl_in_rs.dat',status='unknown')

	arquiv = 'gl_out.dat'
	inquire(file=arquiv,exist=arquiv_gl)
	if (arquiv_gl) then
		continue
	else
		open(unit=20,file='gl_out.dat',status='new',position='append')
		close(20)
	end if

	! file to restart the code in case it is accidentally stoped
	arquiv = 'restart.dat'
	inquire(file=arquiv,exist=arquiv_re)
	open(unit=100,file='restart.dat',status='unknown',form='unformatted')

	arquiv = 'Iterations.dat'
	inquire(file=arquiv,exist=arquiv_Iter)
	if (arquiv_Iter) then
		continue
	else
		open(unit=160,file='Iterations.dat',status='new',position='append')
		close(160)
	end if
	
	arquiv = 'ExecutionTime.dat'
	inquire(file=arquiv,exist=arquiv_ET)
	if (arquiv_ET) then
		continue
	else
		open(unit=110,file='ExecutionTime.dat',status='new',position='append')
		close(110)
	end if
	
	!used to write applied current-difference potential
    arquiv = 'IV.dat'
    inquire(file=arquiv,exist=arquiv_IV)
    if (arquiv_IV) then
        continue
    else
        open(unit=130,file='IV.dat',status='new',position='append')
        close(130)
    end if

    !used to write applied current-difference potential
    arquiv = 'ConvVt.dat'
    inquire(file=arquiv,exist=arquiv_CVT)
    if (arquiv_CVT) then
        continue
    else
        open(unit=150,file='ConvVt.dat',status='new',position='append')
        close(150)
    end if

	!used to write period
    arquiv = 'period.dat'
    inquire(file=arquiv,exist=arquiv_period)
    if (arquiv_period) then
        continue
    else
        open(unit=140,file='period.dat',status='new',position='append')
        close(140)
    end if

	!used to write applied current-difference potential
    arquiv = 'Qm.dat'
    inquire(file=arquiv,exist=arquiv_Qm)
    if (arquiv_Qm) then
        continue
    else
        open(unit=270,file='Qm.dat',status='new',position='append')
        close(270)
    end if

	!used to write applied the time average global temperature
    arquiv = 'Tgm.dat'
    inquire(file=arquiv,exist=arquiv_Tgm)
    if (arquiv_Tgm) then
        continue
    else
        open(unit=280,file='Tgm.dat',status='new',position='append')
        close(280)
    end if    

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! End of opening input and output files !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!!!!!!!!!!!!
	! Read input !
	!!!!!!!!!!!!!!

    read(10,*)alx,aly,alz             ! dimensions of the simulation box
    read(10,*)ak,T0  			      ! kappa, substrate temperature
    read(10,*)anu,zeta,etaT		      ! nu, zeta, etaT
    read(10,*)aJai,aJaf,daJa          ! initial, final1,final2 and step for Ja
    read(10,*)gama				      ! gama
    read(10,*)bdg                     ! to be used in the normal contact
    read(10,*)hf,hs                   ! coefficient of heat transfer for the three surfaces and for the substrate
    read(10,*)dx,dy,dz                ! mesh-grid
    read(10,*)eta,beta,dtFact         ! eta, beta, reduction factor for dtmax
    read(10,*)iVT,iw                  ! every iVT wrtites the V(t)
    read(10,*)iTest,epsP,epsF,epsPh   ! every iTest tests convergence, stops when convergence is reached
    read(10,*)slp				      ! sleep
    read(10,*)istart                  ! istart, 0 -> starts from initial Ja, other -> restarts from last Ja

	if ((istart.eq.0).and.										&
		(arquiv_gl.or.arquiv_re.or.arquiv_ET.or.arquiv_Iter)	&
	   )	then
		open(unit=120,file='Error.dat',status='unknown')
		write(120,*) 'Dat files already exist.'
		write(120,*) 'Set istart = 1 in gl_in.dat file.'
		write(120,*) 'Or delete output dat files if you want to start from the beginning.'
		stop
	end if
	
	iwconv = iw

	!!!!!!!!!!!!!!!!!!!!!!!!
	! End of reading input !
	!!!!!!!!!!!!!!!!!!!!!!!!
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Create directories for working files (linux use mkdir, windowns md), same with -p !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! comment the lines which does not need to save the respective files
	call system('mkdir -p psi')
	call system('mkdir -p psi/xy')
	call system('mkdir -p psi/xz')
	call system('mkdir -p psi/yz')
	call system('mkdir -p vtc')
	call system('mkdir -p J')
	call system('mkdir -p J/xy')
	call system('mkdir -p J/xz')
	call system('mkdir -p J/yz')
	call system('mkdir -p phase')
	call system('mkdir -p field')
	call system('mkdir -p field/yz')
	call system('mkdir -p ht')
	call system('mkdir -p TV')
	call system('mkdir -p restart')
	call system('mkdir -p T')
	call system('mkdir -p Vel')


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! End  creating directories !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Calculation of the size of the mesh !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	Nx = nint(alx/dx) ! sample surface in the x direction
	Ny = nint(aly/dy) ! sample surface in the y direction
	Nz = nint(alz/dz) ! sample surface in the z direction
	
	Ni = nint(1.0d0/dx)+1

	kk = Nz/2+1
	jj = Ny/2+1
	ii = Nx/2+1
	
	Nm = kk
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! End of calculation of the size of the mesh !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	allocate (G(1:Nx+1,1:Ny+1,1:Nz+1))
	allocate (F(1:Nx+1,1:Ny+1,1:Nz+1))
	allocate (Ux(1:Nx,1:Ny+1,1:Nz+1))
	allocate (Uy(1:Nx+1,1:Ny,1:Nz+1))
	allocate (Uz(1:Nx+1,1:Ny+1,1:Nz))
	allocate (hx(1:Nx+1,1:Ny,1:Nz))
	allocate (hy(1:Nx,1:Ny+1,1:Nz))
	allocate (hz(1:Nx,1:Ny,1:Nz+1))
	allocate (fieldOld(1:Nx,1:Ny,1:Nz+1))
	allocate (hty(1:Ny,1:Nz))
	allocate (htz(1:Ny,1:Nz))
	allocate (thetax(1:Nx,1:Ny+1,1:Nz+1))
	allocate (phix(1:Nx,1:Ny+1,1:Nz+1))
	allocate (thetay(1:Nx+1,1:Ny,1:Nz+1))
	allocate (phiy(1:Nx+1,1:Ny,1:Nz+1))
	allocate (thetaz(1:Nx+1,1:Ny+1,1:Nz))
	allocate (phiz(1:Nx+1,1:Ny+1,1:Nz))
	allocate (Qx(1:Nx,1:Ny+1,1:Nz+1))
	allocate (Qy(1:Nx+1,1:Ny,1:Nz+1))
	allocate (Qz(1:Nx+1,1:Ny+1,1:Nz))
	allocate (vtc(1:Nx,1:Ny,1:Nz+1))
	allocate (phase(1:Nx+1,1:Ny+1,1:Nz+1))
	allocate (aJx(1:Nx,1:Ny,1:Nz))
	allocate (aJy(1:Nx,1:Ny,1:Nz))
	allocate (aJz(1:Nx,1:Ny,1:Nz))
	allocate (aJxn(1:Ny,1:Nz))
	allocate (densityOld(1:Nx+1,1:Ny+1,1:Nz+1))
	allocate (densityNew(1:Nx+1,1:Ny+1,1:Nz+1))
	allocate (densityNew2(1:Nx+1,1:Ny+1,1:Nz+1))
	allocate (Vt(1:Nx,1:Ny))
	allocate (aVt(1:Nx-1))
	allocate (delP(1:Nx+1,1:Ny+1,1:Nz+1))
	allocate (Pot(1:Nx+1,1:Ny+1,1:Nz+1))
	allocate (Qot(1:Nx+1,1:Ny+1,1:Nz+1))
	allocate (PotOld(1:Nx+1,1:Ny+1,1:Nz+1))
	allocate (rho(1:Nx+1,1:Ny+1,1:Nz+1))

	allocate (hyn(1:Nx,1:Ny+1,1:Nz))
	allocate (hzn(1:Nx,1:Ny,1:Nz+1))
	
	allocate (aa(1:Nx+1,1:Ny+1,1:Nz+1))
	allocate (TT(1:Nx+1,1:Ny+1,1:Nz+1))
    allocate (SS(1:Nx+1,1:Ny+1,1:Nz+1))
    allocate (Tyz(1:Nx+1,1:Ny+1,1:Nz+1))   
    allocate (W(1:Nx+1,1:Ny+1,1:Nz+1))
    
    allocate (ttx(1:Nx+1))
    allocate (Vortice_Value(1:Ny/2+1))
	
	!!!!!!!!!!!!!!!!!!!!!!!!
	! Format of the output !
	!!!!!!!!!!!!!!!!!!!!!!!!

10	format(2561F12.6)	! format for writing the output variables

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! End of formating of the output !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Calculation of time steps !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	gama2 = gama**2
	gam2i = (0,1)*gama2

	delta = 1.0d0/dx**2 + 1.0d0/dy**2 + 1.0d0/dz**2
	delta = 2.0d0/delta

    dt1 = eta*delta/4.0d0/dsqrt(1.0d0+gama2)
    dt2 = beta*delta/4.0d0/ak**2
    dt3 = anu*delta/(4.0d0*zeta)

	dt = dmin1(dt1,dt2,dt3)	! calculates the dt which garantees convergence
	dt = dt*dtFact
	
	dtP = dtFact*delta/4.0d0 ! time for potential relaxation
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! End of calculation of time steps !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Calculation of the constants to be used !
	! in the recurrence relations, energy     !
	! and currents                            !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	dxdy    = dx*dy     ! area of a cube face
	dydz    = dy*dz
	dzdx    = dz*dx
	vol     = dx*dy*dz/((alx-dx)*(aly-dy)*(alz-dz))	
	NxNyNz	= (Nx-1)*(Ny-1)*(Nz-1)
	
	dxx = 1.0d0/dx**2
	dyy = 1.0d0/dy**2
	dzz = 1.0d0/dz**2

	rx  = dt/(eta*dx**2)		! to be used in the recurrence relations
	ry  = dt/(eta*dy**2)
	rz  = dt/(eta*dz**2)
    c   = dt/eta
    cp  = (0,1)*dt
	d   = dt/beta
	sxy = (ak**2*dt/beta)*(dx/dy)
	sxz = (ak**2*dt/beta)*(dx/dz)
	syx = (ak**2*dt/beta)*(dy/dx)
	syz = (ak**2*dt/beta)*(dy/dz)
	szx = (ak**2*dt/beta)*(dz/dx)
	szy = (ak**2*dt/beta)*(dz/dy)
	dJn = 0.5d0*beta/dx/dt
	dJnb = 0.5d0*beta/dx
	
    sx   = zeta*dt/(anu*dx**2)
    sy   = zeta*dt/(anu*dy**2)
    sz   = zeta*dt/(anu*dx**2)
	wx   = beta**2/(4.0d0*dt*anu*dx**2)
	wy   = beta**2/(4.0d0*dt*anu*dy**2)
	wz   = beta**2/(4.0d0*dt*anu*dy**2)
	wxp  = dt*beta**2/(4.0d0*anu*dx**2)
	wyp  = dt*beta**2/(4.0d0*anu*dy**2)
	wzp  = dt*beta**2/(4.0d0*anu*dz**2)
	su   = etaT*dt/anu
	wt   = eta/(dt*anu)
	wtg  = dt*gama2/(eta*anu)
	
	etaf = 1.0d0-hf
	etas = 1.0d0-hs
    
    rxp   = dtP/dx**2           ! to be used in the potential recurrence relations
    ryp   = dtP/dy**2
    rzp   = dtP/dz**2
    
    if (bdg.eq.0.0d0) then
		bbdg = 0.0d0
	else
		bbdg = (1.0d0-dx/bdg)
	end if
	
	zt_z = zeta*dxdy*(2.0d0/alz)
	
	CcurrentSx  = 1.0d0/dx		! to be used in the current density
	CcurrentSy  = 1.0d0/dy
	CcurrentSz  = 1.0d0/dz

	CaMxy     = 0.5d0*vol*dx/(ak**2*dy)
	CaMyx     = 0.5d0*vol*dy/(ak**2*dx)	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! End the calculation of the constants !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Calculation the local magnetic field at the normal state !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
	! hyn
	do k = 1,Nz
		do j = 1,Ny+1
			do i = 1,Nx

				yy         = (j-0.5d0*Ny-1.0d0)*dy
				zz         = (k-0.5d0*Nz-1.0d0)*dz+0.5d0*dz
				ymnw       = yy-(aly-dy)/2.0d0
				ymaw       = yy+(aly-dy)/2.0d0
				zmnd       = zz-(alz-dz)/2.0d0
				zmad       = zz+(alz-dz)/2.0d0
				c1         = ymnw**2+zmad**2
				c2         = ymnw**2+zmnd**2
				c3         = ymaw**2+zmad**2
				c4         = ymaw**2+zmnd**2
				hyn(i,j,k) = -(-ymnw*dlog(c1/c2)+ymaw*dlog(c3/c4)	&
							   +2.0d0*zmnd*datan(ymnw/zmnd)		&
							   -2.0d0*zmnd*datan(ymaw/zmnd)		&
							   -2.0d0*zmad*datan(ymnw/zmad)		&
							   +2.0d0*zmad*datan(ymaw/zmad)		&
						      )
						     
			end do
		end do
	end do

    ! hzn
    do k = 1,Nz+1
		do j = 1,Ny
			do i = 1,Nx

				yy         = (j-0.5d0*Ny-1)*dy+0.5d0*dy
				zz         = (k-0.5d0*Nz-1.0d0)*dz
				ymnw       = yy-(aly-dy)/2.0d0
				ymaw       = yy+(aly-dy)/2.0d0
				zmnd       = zz-(alz-dz)/2.0d0
				zmad       = zz+(alz-dz)/2.0d0
				c1         = ymaw**2+zmnd**2
				c2         = ymnw**2+zmnd**2
				c3         = ymaw**2+zmad**2
				c4         = ymnw**2+zmad**2
				hzn(i,j,k) = -( zmnd*dlog(c1/c2)-zmad*dlog(c3/c4)	&
							   -2.0d0*ymnw*datan(zmnd/ymnw)		&
							   +2.0d0*ymnw*datan(zmad/ymnw)		&
							   +2.0d0*ymaw*datan(zmnd/ymaw)		&
							   -2.0d0*ymaw*datan(zmad/ymaw)		&
						      )
						     
			end do
		end do
	end do
	
	hyn = hyn/(4.0d0*pi)	
	hzn = hzn/(4.0d0*pi)

	open(unit=170,file='hyn.dat',status='unknown')
	do k = 1,Nz
		write(170,10)(hyn(Nx/2+1,j,k),j=1,Ny+1)
	end do	
	close(170)
	
	open(unit=180,file='hzn.dat',status='unknown')
	do k = 1,Nz+1
		write(180,10)(hzn(Nx/2+1,j,k),j=1,Ny)
	end do	
	close(180)
	
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! End calculation the local magnetic field at the normal state !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	!!!!!!!!!!!!!!!!!!!!!!
	! Initial conditions !
	!!!!!!!!!!!!!!!!!!!!!!
	
	if (istart.eq.0) then  ! (a)
	    
		G = 0.0d0
		do k = 1,Nz+1
		    do j = 1,Ny+1
				do i = 1,Nx+1
					G(i,j,k)  = dsqrt(1-T0)
					TT(i,j,k) = 0
					SS(i,j,k) = 0
				end do
		    end do
		end do
		hx     = 0.0d0
		hy     = 0.0d0
		hz     = 0.0d0
		thetax = 0.0d0
		Ux     = cdexp(-(0,1)*dcmplx(thetax))
		Qx     = 0.0d0
		thetay = 0.0d0
		Uy     = cdexp(-(0,1)*dcmplx(thetay))
		Qy     = 0.0d0
		thetaz = 0.0d0
		Uz     = cdexp(-(0,1)*dcmplx(thetaz))
		Qz     = 0.0d0
		Pot    = 0.0d0
		Qot    = Pot
		PotOld = Pot
		
		time = 0.0d0
		
		avgPotOld  = 0.0d0
		avgPotOld2 = 0.0d0
		voltIVb = 0.0d0

	else ! (a)

		read(100)G,Ux,Uy,Uz,thetax,thetay,thetaz,hx,hy,hz,Pot,aJai,avgPotOld,avgPotOld2,time
		F = G
		PotOld = Pot
		voltIVb = 0.0d0

	end if ! (a)

	densityOld = cdabs(G)
	fieldOld   = hz
	
	if (istart.eq.0) then
		time = 0.0d0	! total time
	end if

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! End of initial conditions !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Rumping up the current density !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	itMax = iteMax/iTest
	
	aJa = aJai

	do while (aJa.le.aJaf)
	
		if (aJa.eq.0.0d0) then
		
			dt1 = eta*delta/4.0d0
			dt2 = beta*delta/4.0d0/ak**2
			dt3 = anu*delta/(4.0d0*zeta)

			dt = dmin1(dt1,dt2,dt3)	! calculates the dt which garantees convergence
			dt = dt*dtFact
			
			rx  = dt/(eta*dx**2)		! to be used in the recurrence relations
			ry  = dt/(eta*dy**2)
			rz  = dt/(eta*dz**2)
			c   = dt/eta
			cp  = (0,1)*dt
			d   = dt/beta
			sxy = (ak**2*dt/beta)*(dx/dy)
			sxz = (ak**2*dt/beta)*(dx/dz)
			syx = (ak**2*dt/beta)*(dy/dx)
			syz = (ak**2*dt/beta)*(dy/dz)
			szx = (ak**2*dt/beta)*(dz/dx)
			szy = (ak**2*dt/beta)*(dz/dy)
			dJn = 0.5d0*beta/dx/dt
			dJnb = 0.5d0*beta/dx
			
			sx   = zeta*dt/(anu*dx**2)
			sy   = zeta*dt/(anu*dy**2)
			sz   = zeta*dt/(anu*dx**2)
			wx   = beta**2/(4.0d0*dt*anu*dx**2)
			wy   = beta**2/(4.0d0*dt*anu*dy**2)
			wz   = beta**2/(4.0d0*dt*anu*dy**2)
			wxp  = dt*beta**2/(4.0d0*anu*dx**2)
			wyp  = dt*beta**2/(4.0d0*anu*dy**2)
			wzp  = dt*beta**2/(4.0d0*anu*dz**2)
			su   = etaT*dt/anu
			wt   = eta/(dt*anu)
			wtg  = dt*gama2/(eta*anu)
			
		else
		
			dt1 = eta*delta/4.0d0/dsqrt(1.0d0+gama2)
			dt2 = beta*delta/4.0d0/ak**2
			dt3 = anu*delta/(4.0d0*zeta)

			dt = dmin1(dt1,dt2,dt3)	! calculates the dt which garantees convergence
			dt = dt*dtFact
			
			rx  = dt/(eta*dx**2)		! to be used in the recurrence relations
			ry  = dt/(eta*dy**2)
			rz  = dt/(eta*dz**2)
			c   = dt/eta
			cp  = (0,1)*dt
			d   = dt/beta
			sxy = (ak**2*dt/beta)*(dx/dy)
			sxz = (ak**2*dt/beta)*(dx/dz)
			syx = (ak**2*dt/beta)*(dy/dx)
			syz = (ak**2*dt/beta)*(dy/dz)
			szx = (ak**2*dt/beta)*(dz/dx)
			szy = (ak**2*dt/beta)*(dz/dy)
			dJn = 0.5d0*beta/dx/dt
			dJnb = 0.5d0*beta/dx
			
			sx   = zeta*dt/(anu*dx**2)
			sy   = zeta*dt/(anu*dy**2)
			sz   = zeta*dt/(anu*dx**2)
			wx   = beta**2/(4.0d0*dt*anu*dx**2)
			wy   = beta**2/(4.0d0*dt*anu*dy**2)
			wz   = beta**2/(4.0d0*dt*anu*dy**2)
			wxp  = dt*beta**2/(4.0d0*anu*dx**2)
			wyp  = dt*beta**2/(4.0d0*anu*dy**2)
			wzp  = dt*beta**2/(4.0d0*anu*dz**2)
			su   = etaT*dt/anu
			wt   = eta/(dt*anu)
			wtg  = dt*gama2/(eta*anu)
			
		end if
	
		iP = 1
	    jP = 1
	    
	    iPA = 1
	    jPA = 1
	    
	    iPosOld = 0
	    
	    VPhi = 0.0d0
	    
	    iOsc  = 0
		iConv = 0

		iw = iwconv+4	! comment this line for the moving average
		
		voltIV    = 0.0d0
		voltIVold = 0.0d0
		iVoltConv = 0
		bVolt     = .false.

		iVmed = 0
		Vmed  = 0.0d0
		
		aJst  = 0.0d0
		aJnt  = 0.0d0
		
		iCountT = 1
		
		heat_m = 0.0d0
		Tgm    = 0.0d0
		
		aJdx = aJa*dx
		
		write(25,'(F8.5)')aJa    ! change the format if it needs
		rewind(25)
		read(25,*) aJastring
		rewind(25)

		! comment the lines which does not need to save the respective files
		call system('mkdir -p psi/xy/'//aJastring)
		call system('mkdir -p psi/xz/'//aJastring)
		call system('mkdir -p psi/yz/'//aJastring)
		call system('mkdir -p vtc/xy/'//aJastring)
		call system('mkdir -p J/xy/'//aJastring)
		call system('mkdir -p J/xz/'//aJastring)
		call system('mkdir -p J/yz/'//aJastring)
		call system('mkdir -p phase/xy/'//aJastring)
		call system('mkdir -p phase/yz/'//aJastring)
		call system('mkdir -p phase/xz/'//aJastring)
		call system('mkdir -p field/perp/'//aJastring)
		call system('mkdir -p field/yz/'//aJastring)
		call system('mkdir -p T/xy/'//aJastring)
		call system('mkdir -p T/xys/'//aJastring)
		call system('mkdir -p T/xym/'//aJastring)
		call system('mkdir -p T/xyi/'//aJastring)
		call system('mkdir -p T/yz/'//aJastring)
		call system('mkdir -p ht/'//aJastring)
		call system('mkdir -p TV/'//aJastring)
		call system('mkdir -p restart/'//aJastring)
		call system('mkdir -p Tx/'//aJastring)
	    call system('mkdir -p restart/'//aJastring)
	    call system('mkdir -p Q/'//aJastring)
	    call system('mkdir -p Tg/'//aJastring)
	    call system('mkdir -p Vel/'//aJastring)
	    call system('mkdir -p pot/xy/'//aJastring)
		
		arquiv = 'restart.dat'
		inquire(file=arquiv,exist=arquiv_re)
		if (arquiv_re) then
			continue
		else
			open(unit=100,file='restart.dat',status='unknown',form='unformatted')
		end if

		write(100)G,Ux,Uy,Uz,thetax,thetay,thetaz,hx,hy,hz,Pot,aJa,avgPotOld,avgPotOld2,time
		rewind(100)
		
		!used to write applied current-time-potential difference
		arquiv = 'TV.dat'
		inquire(file=arquiv,exist=arquiv_TV)
		if (arquiv_TV) then
			call system('rm TV.dat')
			open(unit=120,file='TV.dat',status='new',position='append')
			close(120)			
		else
			open(unit=120,file='TV.dat',status='new',position='append')
			close(120)
		end if

		!used to write the heat transfer rate
		arquiv = 'Q.dat'
		inquire(file=arquiv,exist=arquiv_Q)
		if (arquiv_Q) then
			call system('rm Q.dat')
			open(unit=240,file='Q.dat',status='new',position='append')
			close(240)			
		else
			open(unit=240,file='Q.dat',status='new',position='append')
			close(240)
		end if

		!used to write the global temperatura
		arquiv = 'Tg.dat'
		inquire(file=arquiv,exist=arquiv_Tg)
		if (arquiv_Tg) then
			call system('rm Tg.dat')
			open(unit=290,file='Tg.dat',status='new',position='append')
			close(290)			
		else
			open(unit=290,file='Tg.dat',status='new',position='append')
			close(290)
		end if
		
		!used to write the position of vortex
		arquiv = 'Vel.dat'
		inquire(file=arquiv,exist=arquiv_Vel)
		if (arquiv_Vel) then
			call system('rm Vel.dat')
			open(unit=300,file='Vel.dat',status='new',position='append')
			close(300)			
		else
			open(unit=300,file='Vel.dat',status='new',position='append')
			close(300)
		end if
		
		iFrame  = 0
		cframe  = .true.
		iCountF = 1
		iSaveF  = 0
		
		do k = 1,Nz         ! yz faces
			do j = 1,Ny
				hx(1,j,k)    = 0.0d0
				hx(Nx+1,j,k) = 0.0d0
			end do
		end do

		do k = 1,Nz         ! zx faces
			do i = 2,Nx
				hx(i,1,k)  = 0.0d0
				hx(i,Ny,k) = 0.0d0
			end do
		end do

		do j = 2,Ny-1           ! xy faces
			do i = 2,Nx
				hx(i,j,1)  = 0.0d0
				hx(i,j,Nz) = 0.0d0
			end do
		end do

		do k = 1,Nz         ! zx faces
			do i = 1,Nx
				hy(i,1,k)    = aJa*hyn(i,1,k)
				hy(i,Ny+1,k) = aJa*hyn(i,Ny+1,k)
			end do
		end do

		do j = 2,Ny         ! xy faces
			do i = 1,Nx
				hy(i,j,1)  = aJa*hyn(i,j,1)
				hy(i,j,Nz) = aJa*hyn(i,j,Nz)
			end do
		end do

		do k = 2,Nz-1       ! yz faces
			do j = 2,Ny
				hy(1,j,k)  = aJa*hyn(1,j,k)
				hy(Nx,j,k) = aJa*hyn(Nx,j,k)
			end do
		end do

		do j = 1,Ny         ! xy faces
			do i = 1,Nx
				hz(i,j,1)    = aJa*hzn(i,j,1)
				hz(i,j,Nz+1) = aJa*hzn(i,j,Nz+1)
			end do
		end do

		do k = 2,Nz         ! yz faces
			do j = 1,Ny
				hz(1,j,k)  = aJa*hzn(1,j,k)
				hz(Nx,j,k) = aJa*hzn(Nx,j,k)
			end do
		end do

		do k = 2,Nz         ! zx faces
			do i = 2,Nx-1
				hz(i,1,k)  = aJa*hzn(i,1,k)
				hz(i,Ny,k) = aJa*hzn(i,Ny,k)
			end do
		end do
		
		!!!!!!!!!!!!!!!!!!!!
		! Time iteractions !
		!!!!!!!!!!!!!!!!!!!!
		
		time = 0.0d0
		
		iTmedia = 1
		iStop   = 0

		iiTime = 1
		
		do while (iiTime.le.itMax)

			!$acc data copy(G(:,:,:),F(:,:,:))			            &
			!$acc copy(Qx(:,:,:),Qy(:,:,:),Qz(:,:,:))               &
			!$acc copy(Ux(:,:,:),Uy(:,:,:),Uz(:,:,:))               &
			!$acc copy(thetax(:,:,:),thetay(:,:,:),thetaz(:,:,:))   &
			!$acc copy(phix(:,:,:),phiy(:,:,:),phiz(:,:,:))			&
			!$acc copy(hx(:,:,:),hy(:,:,:),hz(:,:,:))	            &
			!$acc copy(Vt(:,:),aVt(:))      	                    &
			!$acc copy(delP(:,:,:))			                        &
			!$acc copy(TT(:,:,:),SS(:,:,:),W(:,:,:))                &
			!$acc copy(aJxn(:,:),Pot(:,:,:))                        &
			!$acc copy(Qot(:,:,:),PotOld(:,:,:),rho(:,:,:))			&
			!$acc copy(Tyz(:,:,:))		 						
			
			iCountV  = 0

			iCountH  = 0
			
			iCountTg = 0
			
			iCountVe = 0
			
			iTime = 1
			
			iCircle = 0

			do while (iTime.le.iTest)
			
				time = time+dt

				!!!!!!!!!!!!!!!!!!!!!!!!
				! Recurrence relations !
				!!!!!!!!!!!!!!!!!!!!!!!!
				
				! first TDGL equation
				!$acc kernels loop present(G,delP,Ux,Uy,Uz,TT)
				do k = 2,Nz
					do j = 2,Ny
						do i = 2,Nx
							delP(i,j,k) = +rx*(Ux(i,j,k)*G(i+1,j,k)-2.0d0*G(i,j,k)+dconjg(Ux(i-1,j,k))*G(i-1,j,k))  &
										  +ry*(Uy(i,j,k)*G(i,j+1,k)-2.0d0*G(i,j,k)+dconjg(Uy(i,j-1,k))*G(i,j-1,k)) 	&
										  +rz*(Uz(i,j,k)*G(i,j,k+1)-2.0d0*G(i,j,k)+dconjg(Uz(i,j,k-1))*G(i,j,k-1)) 	&
										  +c*G(i,j,k)*(1-T0-TT(i,j,k)-cdabs(G(i,j,k))**2)
                        end do
					end do
				end do
				
				! first TDGL equation
				!$acc kernels loop present(G,F,delP,Pot)
				do k = 2,Nz
					do j = 2,Ny
						do i = 2,Nx
							F(i,j,k) = G(i,j,k)+dsqrt(1.0d0+gama2*cdabs(G(i,j,k))**2)*delP(i,j,k) &
											   -gama2*G(i,j,k)*dreal(dconjg(G(i,j,k))*delP(i,j,k))/dsqrt(1.0d0+gama2*cdabs(G(i,j,k))**2) &
											   -cp*Pot(i,j,k)*G(i,j,k)
						end do
					end do
				end do
				
				! x component of the Ampere's law
				!$acc kernels loop present(G,Ux,thetax,phix,hy,hz,Qx,Pot)
				do k = 2,Nz
					do j = 2,Ny
						do i = 1,Nx
							Qx(i,j,k)   = dimag(dconjg(G(i,j,k))*Ux(i,j,k)*G(i+1,j,k))
							phix(i,j,k) = thetax(i,j,k)+d*Qx(i,j,k)		                &
													   -sxy*(hz(i,j,k)-hz(i,j-1,k))		&
													   +sxz*(hy(i,j,k)-hy(i,j,k-1))	    &
													   -dt*(Pot(i+1,j,k)-Pot(i,j,k))	
							Ux(i,j,k)   = cdexp(-(0,1)*dcmplx(phix(i,j,k)))
						end do
					end do
				end do
				
				! y component of the Ampere's law
				!$acc kernels loop present(G,Uy,thetay,phiy,hx,hz,Qy,Pot)
				do k = 2,Nz
					do j = 1,Ny
						do i = 2,Nx
							Qy(i,j,k)   = dimag(dconjg(G(i,j,k))*Uy(i,j,k)*G(i,j+1,k))
							phiy(i,j,k) = thetay(i,j,k)+d*Qy(i,j,k)			            &
													   -syz*(hx(i,j,k)-hx(i,j,k-1))		&
													   +syx*(hz(i,j,k)-hz(i-1,j,k))     &
													   -dt*(Pot(i,j+1,k)-Pot(i,j,k))
							Uy(i,j,k)   = cdexp(-(0,1)*dcmplx(phiy(i,j,k)))
						end do
					end do
				end do
				
				! z component of the Ampere's law
				!$acc kernels loop present(G,Uz,thetaz,phiz,hx,hy,Qz,Pot)
				do k = 1,Nz
					do j = 2,Ny
						do i = 2,Nx
							Qz(i,j,k)   = dimag(dconjg(G(i,j,k))*Uz(i,j,k)*G(i,j,k+1))
							phiz(i,j,k) = thetaz(i,j,k)+d*Qz(i,j,k)			            &
													   -szx*(hy(i,j,k)-hy(i-1,j,k))		&
													   +szy*(hx(i,j,k)-hx(i,j-1,k))     &
													   -dt*(Pot(i,j,k+1)-Pot(i,j,k))
							Uz(i,j,k)   = cdexp(-(0,1)*dcmplx(phiz(i,j,k)))
						end do
					end do
				end do
				
				! divergent of the supercurrent density
				!$acc kernels loop present(rho,Qx,Qy,Qz)
				do k = 2,Nz
					do j = 2,Ny
						do i = 2,Nx
							rho(i,j,k) = dxx*(Qx(i,j,k)-Qx(i-1,j,k))+dyy*(Qy(i,j,k)-Qy(i,j-1,k))+dzz*(Qz(i,j,k)-Qz(i,j,k-1))
						end do
					end do
				end do
				
				if (iiTime.eq.1) then

					if ((iTime.eq.1).or.(iTime.eq.2)) then
					
						! heat diffusion equation
						!$acc kernels loop present(TT,SS)
						do k = 2,Nz
							do j = 2,Ny
								do i = 2,Nx							
									TT(i,j,k) = SS(i,j,k)	   
								end do
							end do
						end do 

					else

						!heat diffusion equation
						!$acc kernels loop present(SS,TT,W,phix,phiy,phiz,thetax,thetay,thetaz,Pot,F,G,delP,rho)
						do k = 2,Nz
							do j = 2,Ny
								do i = 2,Nx
									W(i,j,k)  =   wx*(phix(i,j,k)-thetax(i,j,k)+phix(i-1,j,k)-thetax(i-1,j,k)              &
												  +dt*(Pot(i+1,j,k)-Pot(i-1,j,k)))**2                                      &
												  +wy*(phiy(i,j,k)-thetay(i,j,k)+phiy(i,j-1,k)-thetay(i,j-1,k)             &
												  +dt*(Pot(i,j+1,k)-Pot(i,j-1,k)))**2                                      &
												  +wz*(phiz(i,j,k)-thetaz(i,j,k)+phiz(i,j,k-1)-thetaz(i,j,k-1)             &
												  +dt*(Pot(i,j,k+1)-Pot(i,j,k-1)))**2                                      &
												  +(1.0d0/dsqrt(1.0d0+gama2*cdabs(G(i,j,k))**2))*(wt*cdabs(delP(i,j,k))**2  &
												                                                +wtg*rho(i,j,k)**2)
									SS(i,j,k) = TT(i,j,k)+sx*(TT(i+1,j,k)-2.0d0*TT(i,j,k)+TT(i-1,j,k)) &
												+sy*(TT(i,j+1,k)-2.0d0*TT(i,j,k)+TT(i,j-1,k)) &
												+sz*(TT(i,j,k+1)-2.0d0*TT(i,j,k)+TT(i,j,k-1)) &
												+W(i,j,k)
								end do
							end do
						end do 
							
					end if

				else

					!heat diffusion equation
					!$acc kernels loop present(SS,TT,W,phix,phiy,phiz,thetax,thetay,thetaz,Pot,F,G,delP,rho)
					do k = 2,Nz
						do j = 2,Ny
							do i = 2,Nx
								W(i,j,k)  =   wx*(phix(i,j,k)-thetax(i,j,k)+phix(i-1,j,k)-thetax(i-1,j,k)              &
											  +dt*(Pot(i+1,j,k)-Pot(i-1,j,k)))**2                                      &
											  +wy*(phiy(i,j,k)-thetay(i,j,k)+phiy(i,j-1,k)-thetay(i,j-1,k)             &
											  +dt*(Pot(i,j+1,k)-Pot(i,j-1,k)))**2                                      &
											  +wz*(phiz(i,j,k)-thetaz(i,j,k)+phiz(i,j,k-1)-thetaz(i,j,k-1)             &
											  +dt*(Pot(i,j,k+1)-Pot(i,j,k-1)))**2                                      &
											  +(1.0d0/dsqrt(1.0d0+gama2*cdabs(G(i,j,k))**2))*(wt*cdabs(delP(i,j,k))**2  &
																						    +wtg*rho(i,j,k)**2)	      										      
								SS(i,j,k) = TT(i,j,k)+sx*(TT(i+1,j,k)-2.0d0*TT(i,j,k)+TT(i-1,j,k)) &
											+sy*(TT(i,j+1,k)-2.0d0*TT(i,j,k)+TT(i,j-1,k)) &
											+sz*(TT(i,j,k+1)-2.0d0*TT(i,j,k)+TT(i,j,k-1)) &
											+W(i,j,k)
							end do
						end do
					end do 
						
				end if
    
				! h = rot A (x component)
				!$acc kernels loop present(phiy,phiz,hx)
				do k = 2,Nz-1
					do j = 2,Ny-1
						do i = 2,Nx
							hx(i,j,k) = (phiz(i,j+1,k)-phiz(i,j,k)-phiy(i,j,k+1)+phiy(i,j,k))/dydz
						end do
					end do
				end do
				
				! h = rot A (y component)
                !$acc kernels loop present(phix,phiz,hy)
                do k = 2,Nz-1
                    do j = 2,Ny
                        do i = 2,Nx-1
                            hy(i,j,k) = (phix(i,j,k+1)-phix(i,j,k)-phiz(i+1,j,k)+phiz(i,j,k))/dzdx
                        end do
                    end do
                end do
                
                ! h = rot A (z component)
                !$acc kernels loop present(phiy,phix,hz)
                do k = 2,Nz
                    do j = 2,Ny-1
                        do i = 2,Nx-1
                            hz(i,j,k) = (phiy(i+1,j,k)-phiy(i,j,k)-phix(i,j+1,k)+phix(i,j,k))/dxdy
                        end do
                    end do
                end do
                
                if (aJa.gt.0.0d0) then
					call poisson(rho,Pot,Qot,PotOld,Nx,Ny,Nz,epsPh,rxP,ryP,rzP,dtP,aJdx)
				end if
				
				!$acc update device(Pot)

				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! End of recurrence relations !
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! Verification of the presence of a ring !
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

				if (verify_ring.eq.0) then

					if (iOsc.gt.iwconv) then

						iCirclea = 0
						!$acc kernels loop present(F)
						do k = 2,Nz/2
							n1 = nint(dimag(cdlog(F(ii,jj,k))))
							n2 = nint(dimag(cdlog(F(ii,jj,k+1))))
							if (abs(n1-n2).eq.3) then
								iCirclea = 1
							end if
						end do
						
						iCircleb = 0
						!$acc kernels loop present(F)
						do j = 2,Ny/2
							n3 = nint(dimag(cdlog(F(ii,j,kk))))
							n4 = nint(dimag(cdlog(F(ii,j+1,kk))))
							if (abs(n3-n4).eq.3) then
								iCircleb = 1
							end if
						end do
						
						if ((iCirclea.eq.1).and.(iCircleb.eq.1)) then
						
							iCircle = 1
							iFrame  = 1
							
						else
		
							if (iFrame.eq.1) then
								cframe = .false.
								iSaveF = 1
							end if
							
						end if
						
					end if
				
				end if
				
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! End verification of the presence of a ring !
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! Calculation of the current !
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				
				if (iOsc.eq.40) then
					iStop = 1
				end if
				
				if ((iOsc.ge.36).and.(iStop.eq.0)) then
				
					aJntb = 0.0d0
					!$acc kernels loop present(Qx,phix,thetax,Pot)
					do k = 2,Nz
						do j = 2,Ny
							aJst  = aJst+0.5d0*CcurrentSx*(Qx(ii,j,k)+Qx(ii-1,j,k))
							aJntb = aJntb-dJn*(phix(ii,j,k)+phix(ii-1,j,k)-thetax(ii,j,k)-thetax(ii-1,j,k)) &
							             -dJnb*(Pot(ii+1,j,k)-Pot(ii-1,j,k))
						end do
					end do
					aJnt = aJnt+aJntb
					
					iTmedia = iTmedia+1
					
				end if
				
				aJntb = 0.0d0
				!$acc kernels loop present(phix,thetax,Pot,aJxn)
				do k = 2,Nz
					do j = 2,Ny
						aJntb = aJntb-dJn*(phix(ii,j,k)+phix(ii-1,j,k)-thetax(ii,j,k)-thetax(ii-1,j,k)) &
							         -dJnb*(Pot(ii+1,j,k)-Pot(ii-1,j,k))
						aJxn(j,k) = -dJn*(phix(ii,j,k)+phix(ii-1,j,k)-thetax(ii,j,k)-thetax(ii-1,j,k)) &
							        -dJnb*(Pot(ii+1,j,k)-Pot(ii-1,j,k))
					end do
				end do
				
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! End of calculation of the current !
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! Calculation of the induced voltage !
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				
				avgPotA = 0.0d0
				!$acc kernels loop present(phix,thetax)
				do j = 2,Ny
					do i = 1,Nx
						avgPotA = avgPotA-(phix(i,j,kk)-thetax(i,j,kk))/dt
					end do
				end do
				Nwp = Ny-1
				avgPotA = avgPotA/Nwp
				
				avgPot = 0.0d0
				!$acc kernels loop present(Pot)
				do j = 2,Ny
						avgPot = avgPot-(Pot(Nx+1,j,kk)-Pot(1,j,kk))
				end do
				avgPot = avgPot/Nwp
				
				if ((avgPot.lt.avgPotOld).and.(avgPotOld2.lt.avgPotOld).and.(iiTIme.gt.1)) then
					iOsc = iOsc+1
					if (iOsc.ge.iwconv) then
						open(unit=140,file='period.dat',status='old',position='append')
						write(140,*)aja,time
						close(140)
					end if
				end if
				avgPotOld2 = avgPotOld
				avgPotOld  = avgPot
				
				! initial potential
				if ((iPA.eq.1).and.(iOsc.eq.(iw-4))) then	! comment this line for the moving average

					Vi = 0.0d0
					!$acc kernels loop present(thetax)
					do j = 2,Ny
						do i = 1,Nx
							Vi = Vi+thetax(i,j,kk)
						end do
					end do
					Vi = Vi/Nwp
					ti = time
					iPA = -1

				end if
				
				! final potential
				if ((jPA.eq.1).and.(iOsc.eq.iw)) then

					Vf = 0.0d0
					!$acc kernels loop present(phix)
					do j = 2,Ny
						do i = 1,Nx
							Vf = Vf+phix(i,j,kk)
						end do
					end do
					Vf = Vf/Nwp
					tf = time
					
					voltIVA = -(Vf-Vi)/(tf-ti)

					iPA = 1	! comment this line for the moving average

				end if
				
				if ((iOsc.eq.(iw-4)).and.(iP.eq.1)) then
					ti = time
					iP = -1
				end if
				
				if ((iOsc.ge.(iw-4)).and.(iOsc.lt.iw)) then	! comment this line for the moving average

					!$acc kernels loop present(Pot)
					do j = 2,Ny
						VPhi = VPhi-(Pot(Nx+1,j,kk)-Pot(1,j,kk))
					end do

				end if
				
				! final potential
				if ((jP.eq.1).and.(iOsc.eq.iw)) then

					tf = time
					voltIV = VPhi*dt/Nwp/(tf-ti)
					erro = dabs(voltIV-voltIVold)
					open(unit=150,file='ConvVt.dat',status='old',position='append')
					write(150,"(F12.8,1x,F12.8,1x,F12.8,1x,F12.8,1x,I8)")aja,voltIV,voltIVA,erro,iOsc
					close(150)
					Vmed = Vmed+voltIV
					iVmed = iVmed+1
					if (erro.lt.1.0d-5) then
						jP = -1
						iVoltConv = 1
						bVolt = .true.
					end if
					voltIVold = voltIV
					iw = iw+4
					
					iP = 1	! comment this line for the moving average
					VPhi = 0.0d0

				end if

				if (aja.eq.0.0d0) then
					iOsc = 0
				end if                
                    
				iCountV = iCountV+1
                if (iCountV.eq.iVT) then				
					if (iOsc.ge.iwconv) then
						open(unit=120,file='TV.dat',status='old',position='append')
						write(120,"(F12.8,1x,F14.4,1x,F12.8,1x,F12.8,1x,I8)")aja,time,avgPot,avgPotA,iOsc
						close(120)
                    end if
                    iCountV = 0
                end if
                
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! End of calculation of the induced voltage !
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! Calculation of the heat transfer rate !
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

				heat = 0.0d0
				!$acc kernels loop present(SS)
				do j = 2,Ny
					do i = 2,Nx
						heat = heat-(SS(i,j,Nz+1)-SS(i,j,kk))
					end do
				end do
				heat_m = heat_m+heat
				heat = zt_z*heat

                iCountH = iCountH+1
                if (iCountH.eq.iVT) then				
					if (iOsc.ge.iwconv) then
						open(unit=240,file='Q.dat',status='old',position='append')
						write(240,"(F12.8,1x,F14.4,1x,F12.8,1x,I8)")aja,time,heat,iOsc
						close(240)
                    end if
                    iCountH = 0
                end if

                !if (iOsc.ge.iwconv) then
				!	!$acc kernels loop present(SS,TT)
				!	do j = 2,Ny
				!		do i = 2,Nx
				!			heat_m = heat_m-((SS(i,j,Nz+1)-SS(i,j,kk))-(TT(i,j,Nz+1)-TT(i,j,kk)))
				!		end do
				!	end do
                !end if

				Tg = 0.0d0
				!$acc kernels loop present(SS)
				do k = 2,Nz
					do j = 2,Ny
						do i = 2,Nx
							Tg = Tg+SS(i,j,k)
						end do
					end do
				end do
				Tg  = Tg/NxNyNz
				Tgm = Tgm+Tg

                iCountTg = iCountTg+1
                if (iCountTg.eq.iVT) then				
					if (iOsc.ge.iwconv) then
						open(unit=290,file='Tg.dat',status='old',position='append')
						write(290,"(F12.8,1x,F14.4,1x,F12.8,1x,I8)")aja,time,Tg,iOsc
						close(290)
                    end if
                    iCountTg = 0
                end if				

                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! End of calculation of the heat transfer rate !
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! Calculation of the vortice velocity!
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				
				iPos = 0
				!do j = 2,Ny/2
				!	Vortice_value(j)= cdabs(F(Nx/2+1,j,Nz/2+1))
				!end do
				
				do j = 3,Ny/2
					!if ((Vortice_value(j).lt.Vortice_value(j+1)).and.(Vortice_value(j).lt.Vortice_value(j-1))) then
					if ((cdabs(F(Nx/2+1,j,Nz/2+1)).lt.cdabs(F(Nx/2+1,j+1,Nz/2+1))).and.(cdabs(F(Nx/2+1,j,Nz/2+1)).lt.cdabs(F(Nx/2+1,j-1,Nz/2+1)))) then
						iPos = j
					end if
					!if (iPos.gt.0) then
					!if (iPos.gt.iPosOld) then
					!	open(unit=300,file='Vel.dat',status='old',position='append')
					!	write(300,"(F16.8,1x,I8)")time,iPos
					!	close(300)
					!end if 
				end do
				if (iPos.gt.iPosOld) then
					open(unit=300,file='Vel.dat',status='old',position='append')
					write(300,"(F16.8,1x,I8)")time,iPos
					close(300)
				end if
				iPosOld = iPos
				
				if (iPos.gt.0) then
				
					!$acc kernels loop present(ttx,SS)
					do i = 1,Nx+1
						ttx(i) = SS(i,Ny/2+1,Nz/2+1)
					end do
					
				end if
					
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! End of calculation of the vortice velocity!
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! Boundary conditions for the order parameter !
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

				!$acc kernels loop present(F)
				do k = 2,Nz             ! xz faces
					do j = 2,Ny
						F(1,j,k) = bbdg*Ux(1,j,k)*F(2,j,k)
						F(Nx+1,j,k) = bbdg*dconjg(Ux(Nx,j,k))*F(Nx,j,k)
					end do
				end do
				
				!$acc kernels loop present(F,Uy)
				do k = 2,Nz             ! xz faces
					do i = 2,Nx
						F(i,1,k) = Uy(i,1,k)*F(i,2,k)
						F(i,Ny+1,k) = dconjg(Uy(i,Ny,k))*F(i,Ny,k)
					end do
				end do
				
				!$acc kernels loop present(F,Uz)
				do j = 2,Ny             ! xy faces
					do i = 2,Nx
						F(i,j,1) = Uz(i,j,1)*F(i,j,2)
						F(i,j,Nz+1) = dconjg(Uz(i,j,Nz))*F(i,j,Nz)
					end do
				end do
				
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! End of boundary conditions for the order parameter !
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				!Boundary conditions for Temperature!
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				
				!$acc kernels loop present(SS)
				do k = 2,Nz             ! xz faces
					do j = 2,Ny
						SS(1,j,k) = 0.0d0
						SS(Nx+1,j,k) = 0.0d0
					end do
				end do
				
				!$acc kernels loop present(SS)
				do k = 2,Nz             ! xz faces
					do i = 2,Nx
						SS(i,1,k) = etaf*SS(i,2,k)
						SS(i,Ny+1,k) = etaf*SS(i,Ny,k)
					end do
				end do
				
				!$acc kernels loop present(SS)
				do j = 2,Ny             ! xy faces
					do i = 2,Nx
						SS(i,j,1) = etas*SS(i,j,2)
						SS(i,j,Nz+1) = etaf*SS(i,j,Nz)
					end do
				end do
								
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				!End  boundary conditions for temperature!
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

				!!!!!!!!!!!!!!!!!!!!!!!!!!
				! Transferring variables !
				!!!!!!!!!!!!!!!!!!!!!!!!!!

				!$acc kernels loop
				do k = 1,Nz+1
					do j = 1,Ny+1
						do i = 1,Nx+1
							G(i,j,k) = F(i,j,k)
						end do
					end do
				end do
				
				!$acc kernels loop
				do k = 1,Nz+1
					do j = 1,Ny+1
						do i = 1,Nx+1
							TT(i,j,k) = SS(i,j,k)
						end do
					end do
				end do
				
				!$acc kernels loop
				do k = 1,Nz+1
					do j = 1,Ny+1
						do i = 1,Nx
							thetax(i,j,k) = phix(i,j,k)
						end do
					end do
				end do
				
				!$acc kernels loop	
				do k = 1,Nz+1
					do j = 1,Ny
						do i = 1,Nx+1
							thetay(i,j,k) = phiy(i,j,k)
						end do
					end do
				end do
				
				!$acc kernels loop	
				do k = 1,Nz
					do j = 1,Ny+1
						do i = 1,Nx+1
							thetaz(i,j,k) = phiz(i,j,k)
						end do
					end do
				end do

				!$acc kernels loop
				do k = 1,Nz+1
					do j = 1,Ny+1
						do i = 1,Nx+1
							Qot(i,j,k) = Pot(i,j,k)
							PotOld(i,j,k) = Pot(i,j,k)
						end do
					end do
				end do
						
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! End of transferring variables !
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				
				if ((iTime.eq.1).and.(iiTime.eq.1)) then 
				    ti = time
				end if
				
				iTime = iTime+1

			end do
			
			!$acc end data

			if (spl.eq.1) then
				call sleep (2)
			end if
			
			do i = 1,Nx+1       ! cosmetic
				F(i,1,1) = 0.5d0*(F(i,1,2)+F(i,2,1))
				F(i,1,Nz+1) = 0.5d0*(F(i,1,Nz)+F(i,2,Nz+1))
				F(i,Ny+1,1) = 0.5d0*(F(i,Ny+1,2)+F(i,Ny,1))
				F(i,Ny+1,Nz+1) = 0.5d0*(F(i,Ny,Nz)+F(i,Ny,Nz))
			end do
			
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! Calculation of the phase !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			do k = 1,Nz+1
				do j = 1,Ny+1
					do i = 1,Nx+1
						phase(i,j,k) = dimag(cdlog(F(i,j,k)))
					end do
				end do
			end do
			
			do k = 2,Nz
				do j = 2,Ny
					do i = 2,Nx
						aJx(i,j,k) = 0.5d0*CcurrentSx*(Qx(i,j,k)+Qx(i-1,j,k))
						aJy(i,j,k) = 0.5d0*CcurrentSy*(Qy(i,j,k)+Qy(i,j-1,k))
						aJz(i,j,k) = 0.5d0*CcurrentSz*(Qz(i,j,k)+Qz(i,j,k-1))
					end do
				end do
			end do

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! End of calculation of the phase !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			
			densityNew  = cdabs(F)
			
			if (iOsc.ge.iwconv) then
			
				write(15,*)iCountT
				rewind(15)
				read(15,*)CiCountT
				rewind(15)
				
				pxyw  = 'psiXY.'//CiCountT
				pxzw  = 'psiXZ.'//CiCountT
				pyzw  = 'psiYZ.'//CiCountT
				potw  = 'phiXY.'//CiCountT
				phyzw = 'phaseYZ.'//CiCountT
				phxzw = 'phaseXZ.'//CiCountT
				jxw   = 'jxsYZ.'//CiCountT
				jnxw  = 'jxnYZ.'//CiCountT
				tyzw = 'tyz.'//CiCountT
				txysw = 'txys.'//CiCountT
				txymw = 'txym.'//CiCountT
				txyiw = 'txyi.'//CiCountT

				open(unit=16,file=pxyw,status='unknown')
				open(unit=17,file=pxzw,status='unknown')
				open(unit=18,file=pyzw,status='unknown')
				open(unit=19,file=phyzw,status='unknown')
				open(unit=20,file=phxzw,status='unknown')
				open(unit=300,file=jxw,status='unknown')
				open(unit=301,file=jnxw,status='unknown')
				open(unit=302,file=potw,status='unknown')
				open(unit=303,file=tyzw,status='unknown')
				open(unit=304,file=txysw,status='unknown')
				open(unit=305,file=txymw,status='unknown')
				open(unit=306,file=txyiw,status='unknown')
				
				do j = 1,Ny+1
					write(16,10)(densityNew(i,j,Nm),i=1,Nx+1)
				end do
				
				do k = 1,Nz+1
					write(17,10)(densityNew(i,jj,k),i=1,Nx+1)
				end do
				
				do k = 1,Nz+1
					write(18,10)(densityNew(ii,j,k),j=1,Ny+1)
				end do
				
				do j = 1,Ny+1
					write(302,10)(Pot(i,j,Nm),i=1,Nx+1)
				end do
				
				do k = 1,Nz+1
					write(19,10)(phase(ii,j,k),j=1,Ny+1)
				end do

				do k = 1,Nz+1
					write(20,10)(phase(i,jj,k),i=1,Nx+1)
				end do
				
				do k = 2,Nz
					write(300,10)(aJx(Nx,j,k),j=2,Ny)
				end do
				
				do k = 2,Nz
					write(301,10)(aJxn(j,k),j=2,Ny)
				end do

				do k = 1,Nz+1
					write(303,10)(SS(ii,j,k),j=1,Ny+1)
				end do

				do j = 1,Ny+1
					write(304,10)(SS(i,j,Nz+1),i=1,Nx+1)
				end do

				do j = 1,Ny+1
					write(305,10)(SS(i,j,Nz/2+1),i=1,Nx+1)
				end do

				do j = 1,Ny+1
					write(306,10)(SS(i,j,1),i=1,Nx+1)
				end do				
								
				close(16)
				close(17)
				close(18)
				close(19)
				close(20)
				close(300)
				close(301)
				close(302)
				close(303)
				close(304)
				close(305)
				close(306)
								
				call system('mv '//pxyw//' psi/xy/'//aJastring)
				call system('mv '//pxzw//' psi/xz/'//aJastring)
				call system('mv '//pyzw//' psi/yz/'//aJastring)
				call system('mv '//phyzw//' phase/yz/'//aJastring)
				call system('mv '//phxzw//' phase/xz/'//aJastring)
				call system('mv '//jxw//' J/yz/'//aJastring)
				call system('mv '//jnxw//' J/yz/'//aJastring)
				call system('mv '//potw//' pot/xy/'//aJastring)
				call system('mv '//tyzw//' T/yz/'//aJastring)
				call system('mv '//txysw//' T/xys/'//aJastring)
				call system('mv '//txymw//' T/xym/'//aJastring)
				call system('mv '//txyiw//' T/xyi/'//aJastring)
				
				
				
				jxw  = 'jxXY.'//CiCountT
				jyw  = 'jyXY.'//CiCountT
				
				open(unit=21,file=jxw,status='unknown')
				open(unit=22,file=jyw,status='unknown')		
				
				do j = 1,Ny+1
					write(21,10)(aJx(i,j,kk),i=1,Nx+1)
				end do
				
				do j = 1,Ny+1
					write(22,10)(aJy(i,j,kk),i=1,Nx+1)
				end do
				
				close(21)
				close(22)
				
				call system('mv '//jxw//' J/xy/'//aJastring)
				call system('mv '//jyw//' J/xy/'//aJastring)
				
				jxw  = 'jxXZ.'//CiCountT
				jzw  = 'jzXZ.'//CiCountT
				
				open(unit=21,file=jxw,status='unknown')
				open(unit=23,file=jzw,status='unknown')
				
				do k = 1,Nz+1
					write(21,10)(aJx(i,jj,k),i=1,Nx+1)
				end do
				
				do k = 1,Nz+1
					write(23,10)(aJz(i,jj,k),i=1,Nx+1)
				end do
				
				close(21)
				close(23)
				
				call system('mv '//jxw//' J/xz/'//aJastring)
				call system('mv '//jzw//' J/xz/'//aJastring)
				
				iCountT = iCountT+1
				
			end if									

			!!!!!!!!!!!!!!!!!!!!
			! Convergence test !
			!!!!!!!!!!!!!!!!!!!!
			
			! this loop determines the largest difference between |psiNew| and |psiOld| 		    
			diffP = 0.0d0
			do k = 2,Nz
				do j = 2,Ny
					do i = 2,Nx
						dif = dabs(densityNew(i,j,k)-densityOld(i,j,k))
						densityOld(i,j,k) = densityNew(i,j,k)
						if (dif.gt.diffP) then
							diffP = dif
						end if
	    			end do
				end do
    		end do
    		
    		! this loop determines the largest difference between |hzNew| and |hzOld| 
			diffH = 0.0d0
			do k  = 1,Nz+1
				do j = 1,Ny
					do i = 1,Nx
						dif = dabs(hz(i,j,k)-fieldOld(i,j,k))
						fieldOld(i,j,k) = hz(i,j,k)
						if (dif.gt.diffH) then
							diffH = dif
						end if
		    		end do
				end do
    		end do
    		
    		! saves the number of iterations and the respective error
			open(unit=20,file='gl_out.dat',status='old',position='append')

			write(20,"(F12.8,1x,I8,1x,F12.8,1x,F12.8,1x,I8)")aJa,iiTime*iTest,diffP,diffH,iOsc
			
			close(20)

            ! If converged, it will stop
            if ((diffP.lt.epsP).and.(diffH.lt.epsF)) then
                open(unit=20,file='gl_out.dat',status='old',position='append')
                 write(20,*)'--------------- Convergence was reached -------------------'
                 write(20,*)' Curr denst   Iterats  Psi error    hz error          iOsc '
                close(20)
                iConv = 1
                exit
			else if (bVolt) then
			!else if (bVolt.and.(iOsc.ge.(iwconv+32))) then
                open(unit=20,file='gl_out.dat',status='old',position='append')
                write(20,*)'------- Desired number of oscilations was reached ---------'
                write(20,*)' Curr denst   Iterats  Psi error    hz error          iOsc '
                close(20)
                exit
            end if
            if (iOsc.gt.220) then
				voltIV = Vmed/iVmed
                open(unit=20,file='gl_out.dat',status='old',position='append')
                write(20,*)'------- Desired number of oscilations was reached ---------'
                write(20,*)' Curr denst   Iterats  Psi error    hz error          iOsc '
                close(20)
                exit
            end if

			!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! End of convergence test !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!
			
			iiTime = iiTime+1

		end do ! Main loop

		!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! End of time iteractions !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		!!!!!!!!!!!!!!!!!
		! Saving output !
		!!!!!!!!!!!!!!!!!
		
		do k = 2,Nz
			do j = 2,Ny
				do i = 2,Nx
					aJx(i,j,k) = 0.5d0*CcurrentSx*(Qx(i,j,k)+Qx(i-1,j,k))
					aJy(i,j,k) = 0.5d0*CcurrentSy*(Qy(i,j,k)+Qy(i,j-1,k))
					aJz(i,j,k) = 0.5d0*CcurrentSz*(Qz(i,j,k)+Qz(i,j,k-1))
				end do
			end do
		end do
		
		! saves normal and superconducting fields

		densityNew2 = densityNew*densityNew

		txw  = 'ttx.dat'
		open(unit=1050,file=txw,status='unknown')
		do i = 1,Nx+1
			xx = (i - 0.5*Nx-1)*dx
			write(1050,"(F12.8,1x,F12.8)") xx,ttx(i)
		end do	
		close(1050)
		call system('mv '//txw//' Tx/'//aJastring)
		
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Calculation of Induced Voltage !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		! average potential
        if (aja.gt.0.0d0) then
            if (iOsc.le.1) then
                VV  = 0.0d0
                tf  = time
                Prd = tf-ti
            else
                VV  = -(Vf-Vi)
                Prd = tf-ti
            end if
        else
            VV  = 0.0d0
            tf  = time
            Prd = tf-ti
        end if
        
        aJt = aJa*(aly-dy)*(alz-dz)
        
        if (voltIV.eq.0.0d0) then
			
			do j = 2,Ny
				voltIV = voltIV-(Pot(Nx+1,j,kk)-Pot(1,j,kk))
			end do
			voltIV = voltIV/Nwp
			
		end if
        
        open(unit=130,file='IV.dat',status='old',position='append')
        write(130,"(F12.8,1x,F12.8,1x,F12.8,1x,F12.8,1x,F16.8,I8)")aja,aJt,voltIV,voltIVA,Prd,iOsc
        close(130)

		open(unit=270,file='Qm.dat',status='old',position='append')
        write(270,"(F12.8,1x,F16.8,1x,I8)")aja,zt_z*heat_m*dt/time,iOsc
        close(270)

		open(unit=280,file='Tgm.dat',status='old',position='append')
        write(280,"(F12.8,1x,F16.8,1x,I8)")aja,Tgm*dt/time,iOsc
        close(280)
 		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! End of calculation of Induced Voltage !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Calculation of the currents !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		do k = 2,Nz
			do j = 2,Ny
				do i = 2,Nx
					aJx(i,j,k) = 0.5d0*CcurrentSx*(Qx(i,j,k)+Qx(i-1,j,k))
					aJy(i,j,k) = 0.5d0*CcurrentSy*(Qy(i,j,k)+Qy(i,j-1,k))
				end do
			end do
		end do
		
		if (iOsc.gt.24) then
			open(unit=500,file='Current.dat',status='unknown',position='append')
			write(500,"(F12.8,1x,F12.8,1x,F12.8)")aJt,aJst*dy*dz/iTmedia,aJnt*dy*dz/iTmedia
			close(500)
		else
			aJst = 0.0d0
			do k = 2,Nz
				do j = 2,Ny
					aJst = aJst+aJx(ii,j,k)
				end do
			end do
			aJst = aJst*dy*dz
			open(unit=500,file='Current.dat',status='unknown',position='append')
			write(500,"(F12.8,1x,F12.8,1x,F12.8)")aJt,aJst,aJntb*dy*dz
			close(500)
		end if
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! End the calculation of the currents !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Calculation of number of vortices !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		do k = 1,Nz+1
			do j = 1,Ny
				do i = 1,Nx
					vtc(i,j,k) = dimag(cdlog(F(i+1,j,k)*dconjg(F(i,j,k))))
					vtc(i,j,k) = vtc(i,j,k)+dimag(cdlog(F(i+1,j+1,k)*dconjg(F(i+1,j,k))))
					vtc(i,j,k) = vtc(i,j,k)+dimag(cdlog(F(i,j+1,k)*dconjg(F(i+1,j+1,k))))
					vtc(i,j,k) = vtc(i,j,k)+dimag(cdlog(F(i,j,k)*dconjg(F(i,j+1,k))))
					vtc(i,j,k) = nint(vtc(i,j,k)/dpi)
				end do
			end do
		end do

		aNV = 0
		do j = 2,Ny
			do i = 1,Nx
				aNV = aNV+vtc(i,j,kk)
			end do
		end do

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! End of calculation of number of vortices !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Calculation of the phase !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		do k = 1,Nz+1
			do j = 1,Ny+1
				do i = 1,Nx+1
					phase(i,j,k) = dimag(cdlog(F(i,j,k)))
				end do
			end do
		end do

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! End of calculation of the phase !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		ii = Nx/2+1
		do k = 2,Nz
			do j = 2,Ny
				hty1     = 0.5d0*(hy(ii-1,j,k)+hy(ii,j,k))
				hty2     = 0.5d0*(hy(ii-1,j,k-1)+hy(ii,j,k-1))
				hty(j,k) = 0.5d0*(hty1+hty2)
				htz1     = 0.5d0*(hz(ii-1,j-1,k)+hz(ii,j-1,k))
				htz2     = 0.5d0*(hz(ii-1,j,k)+hz(ii,j,k))
				htz(j,k) = 0.5d0*(htz1+htz2)
			end do
		end do
		
		if ((iConv.gt.0).and.(iOsc.lt.iw)) then
		
			call system('rm -r psi/xy/'//aJastring)
			call system('rm -r psi/xz/'//aJastring)
			call system('rm -r psi/yz/'//aJastring)
			call system('rm -r pot/xy/'//aJastring)
			call system('rm -r phase/xy/'//aJastring)
			call system('rm -r phase/yz/'//aJastring)
			call system('rm -r J/xy/'//aJastring)
			call system('rm -r J/xz/'//aJastring)
			call system('rm -r J/yz/'//aJastring)

			! comment the lines which does not need to save the respective files
			call system('mkdir -p psi/xy/'//aJastring)
			call system('mkdir -p psi/xz/'//aJastring)
			call system('mkdir -p psi/yz/'//aJastring)
			call system('mkdir -p phase/xy/'//aJastring)
			call system('mkdir -p phase/yz/'//aJastring)
			call system('mkdir -p J/xy/'//aJastring)
			call system('mkdir -p J/xz/'//aJastring)
			call system('mkdir -p J/yz/'//aJastring)
			call system('mkdir -p vtc/xy/'//aJastring)
			call system('mkdir -p T/xy/'//aJastring)
			call system('mkdir -p pot/xy/'//aJastring)
			
			! comment the lines which do not need to save the respective files
			pxyw   = 'psiXY.dat'
			pxzw   = 'psiXZ.dat'
			pyzw   = 'psiYZ.dat'
			hxw    = 'hx.dat'
			hyw    = 'hy.dat'
			hzw    = 'hz.dat'
			vw     = 'vtc.dat'
			phxyw  = 'phaseXY.dat'
			jxw    = 'Jx.dat'
			jyw    = 'Jy.dat'
			htyw   = 'hty.dat'
			htzw   = 'htz.dat'
			txyw   = 'T.dat'
			jxbw   = 'jxsYZ.dat'
			jnxw   = 'jxnYZ.dat'
			potw   = 'phiXY.dat'

			open(unit=16,file=pxyw,status='unknown')
			open(unit=28,file=pxzw,status='unknown')
			open(unit=29,file=pyzw,status='unknown')
			open(unit=27,file=hxw,status='unknown')
			open(unit=26,file=hyw,status='unknown')
			open(unit=17,file=hzw,status='unknown')
			open(unit=18,file=vw,status='unknown')
			open(unit=19,file=phxyw,status='unknown')
			open(unit=21,file=jxw,status='unknown')
			open(unit=22,file=jyw,status='unknown')
			open(unit=23,file=htyw,status='unknown')
			open(unit=24,file=htzw,status='unknown')
			open(unit=31,file=txyw,status='unknown')
			open(unit=300,file=jxbw,status='unknown')
			open(unit=301,file=jnxw,status='unknown')
			open(unit=302,file=potw,status='unknown')
			
			do j = 1,Ny+1
				write(16,10)(densityNew(i,j,Nm),i=1,Nx+1)
			end do
			
			do k = 1,Nz+1
				write(28,10)(densityNew(i,jj,k),i=1,Nx+1)
			end do
			
			do k = 1,Nz+1
				write(29,10)(densityNew(ii,j,k),j=1,Ny+1)
			end do
			
			do j = 1,Ny+1
				write(302,10)(Pot(i,j,Nm),i=1,Nx+1)
			end do
			
			do k = 1,Nz
				write(27,10)(hx(ii,j,k),j=1,Ny)
			end do
			
			do k = 1,Nz
				write(26,10)(hy(i,jj,k),i=1,Nx)
			end do
			
			do j = 1,Ny
				write(17,10)(hz(i,j,Nm),i=1,Nx)
			end do
			
			do j = 1,Ny
				write(18,10)(vtc(i,j,Nm),i=1,Nx)
			end do
			
			do j = 1,Ny+1
				write(19,10)(phase(i,j,Nm),i=1,Nx)
			end do
			
			do j = 2,Ny
				write(21,10)(aJx(i,j,Nm),i=2,Nx)
			end do

			do j = 2,Ny
				write(22,10)(aJy(i,j,Nm),i=2,Nx)
			end do
			
			do k = 2,Nz
				write(23,10)(hty(j,k),j=2,Ny)
			end do

			do k = 2,Nz
				write(24,10)(htz(j,k),j=2,Ny)
			end do

			do j = 1,Ny+1
				write(31,10)(SS(i,j,Nm),i=1,Nx+1)
			end do
			
			do k = 2,Nz
				write(300,10)(aJx(Nx,j,k),j=2,Ny)
			end do
			
			do k = 2,Nz
				write(301,10)(aJxn(j,k),j=2,Ny)
			end do
			
			close(16)
			close(17)
			close(18)
			close(19)
			close(21)
			close(22)
			close(23)
			close(24)
			close(26)
			close(27)
			close(28)
			close(29)
			close(31)
			close(300)
			close(301)
			close(302)
			
			call system('mv '//pxyw//' psi/xy/'//aJastring)
			call system('mv '//pxzw//' psi/xz/'//aJastring)
			call system('mv '//pyzw//' psi/yz/'//aJastring)
			call system('mv '//jxw //' J/xy/'//aJastring)
			call system('mv '//jyw //' J/xy/'//aJastring)
			call system('mv '//jxbw //' J/yz/'//aJastring)
			call system('mv '//jnxw //' J/yz/'//aJastring)
			call system('mv '//phxyw //' phase/xy/'//aJastring)
			call system('mv '//vw //' vtc/xy/'//aJastring)
			call system('mv '//hxw //' field/perp/'//aJastring)
			call system('mv '//hyw //' field/perp/'//aJastring)
			call system('mv '//hzw //' field/perp/'//aJastring)
			call system('mv '//htyw //' ht/'//aJastring)
			call system('mv '//htzw //' ht/'//aJastring)
			call system('mv '//txyw //' T/xy/'//aJastring)
			call system('mv '//potw//' pot/xy/'//aJastring)
			
			iConv = 0
			
		end if
		
		open(unit=180,file='hz.dat',status='unknown')
		do k = 1,Nz+1
			write(180,10)(hz(Nx/2+1,j,k),j=1,Ny)
		end do	
		close(180)

		open(unit=210,file='hy.dat',status='unknown')
		do k = 1,Nz
			write(210,10)(hy(Nx/2+1,j,k),j=1,Ny+1)
		end do	
		close(210)
		
		call system('mv hz.dat  field/yz/'//aJastring)
		call system('mv hy.dat  field/yz/'//aJastring)
		call system('mv TV.dat  TV/'//aJastring)
		call system('mv Q.dat  Q/'//aJastring)
		call system('mv Tg.dat  Tg/'//aJastring)
		close(100)
		call system('mv restart.dat  restart/'//aJastring)
		call system('mv Vel.dat  Vel/'//aJastring)

		!!!!!!!!!!!!!!!!
		
		! saves a message only if the number of iterations 
		! exceeeds the maximum number itMax
		if (iiTime.ge.itMax) then

			open(unit=20,file='gl_out.dat',status='old',position='append')
			
			write(20,*)'-------------- Convergence was not reached ----------------'
			
			close(20)

		end if
		
		!!!!!!!!!!!!!!!!
		
		! saves the the number of iterations (Iterations.dat)
		! and the excecution time (ExecutionTime.dat), for each current

		open(unit=160,file='Iterations.dat',status='old',position='append')

		write(160,"(F12.8,1x,I8)")aJa,iiTime*iTest

		close(160)

		call cpu_time(time_final)
		total_time   = time_final-time_initial
		time_initial = time_final

		open(unit=110,file='ExecutionTime.dat',status='old',position='append')

		write(110,"(F12.8,1x,F12.8)")aJa,total_time/60.0d0

		close(110)

		!!!!!!!!!!!!!!!!!!!!!!!
		! End of saving ouput !
		!!!!!!!!!!!!!!!!!!!!!!!

		psi2Max = maxval(densityNew2)		
		if (psi2Max.lt.epsP) then 
			exit
		end if	
		
		aJa = aJa+daJa	

	end do ! closes current loop

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! End of rumping up the current density !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
end program resistive_state_xi

subroutine poisson(rho,Pot,Qot,PotOld,Nx,Ny,Nz,convP,rxP,ryP,rzP,dtP,aJdx)
!$acc routine seq
	integer i, j, Nx, Ny, Nz, iTime, iiTime
	real (kind=8) rho(Nx+1,Ny+1,Nz+1),Pot(Nx+1,Ny+1,Nz+1),Qot(Nx+1,Ny+1,Nz+1),PotOld(Nx+1,Ny+1,Nz+1)
	real (kind=8) dtP, rxP, ryP, rzP, convP, dif, diff, aJdx
	
	iiTimeP = 1
	
	do while (iiTimeP.le.1000)
	
	iTimeP = 1
	
	do while (iTimeP.le.20)
	
		do k = 2,Nz
			do j = 2,Ny
				do i = 2,Nx
					Pot(i,j,k) = Qot(i,j,k)+rxP*(Qot(i+1,j,k)-2.0d0*Qot(i,j,k)+Qot(i-1,j,k)) &
									       +ryP*(Qot(i,j+1,k)-2.0d0*Qot(i,j,k)+Qot(i,j-1,k)) &
									       +rzP*(Qot(i,j,k+1)-2.0d0*Qot(i,j,k)+Qot(i,j,k-1)) &
									       -dtP*rho(i,j,k)
				end do
			end do
		end do
		
		do k = 2,Nz
			do j = 2,Ny
				Pot(1,j,k)    = Pot(2,j,k)+aJdx
				Pot(Nx+1,j,k) = Pot(Nx,j,k)-aJdx
			end do
		end do
		
		do k = 2,Nz
			do i = 2,Nx
				Pot(i,1,k)    = Pot(i,2,k)
				Pot(i,Ny+1,k) = Pot(i,Ny,k)
			end do
		end do
		
		do j = 2,Ny
			do i = 2,Nx
				Pot(i,j,1)    = Pot(i,j,2)
				Pot(i,j,Nz+1) = Pot(i,j,Nz)
			end do
		end do
		
		do k = 1,Nz+1
			do j = 1,Ny+1
				do i = 1,Nx+1
					Qot(i,j,k) = Pot(i,j,k)
				end do
			end do
		end do
		
		iTimeP = iTimeP+1
		
	end do
	
	dif = 0.0d0
	do k = 2,Nz	
		do j = 2,Ny
			do i = 2,Nx
				diff = dabs(Pot(i,j,k)-PotOld(i,j,k))
				PotOld(i,j,k) = Pot(i,j,k)
				if (diff.gt.dif) then
					dif = diff
				end if
			end do
		end do
	end do
	
	if (dif.lt.convP) then
		exit
	end if
	
	iiTimeP = iiTimeP+1
	
	end do
	
end subroutine
