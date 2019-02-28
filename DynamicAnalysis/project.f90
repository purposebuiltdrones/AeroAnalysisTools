! Sam Weiss
!Flight Mechanics Part 2

!Compile with: '$ gfortran -fdefault-real-8 eigensolver.f90 project.f90'

program project
implicit none


!Eigen Solver Stuff
!Input Stuff
integer :: IV = 1
integer :: N = 12
real, dimension(:,:), allocatable :: Ar
real, dimension(:,:), allocatable :: Br
!Output stuff
complex, dimension(:), allocatable :: EVAL
complex, dimension(:,:), allocatable :: V
integer :: IERR

!Constants
real, parameter :: pi = 4.0*atan(1.0) 
real, parameter :: g = 9.80665!m/s^2

!Input Variables
real :: V_0 !Flight Airspeed (m/s)
real :: MachNum !Flight Mach Number
real :: rho !Air Density (kg/m^3)
real :: W !Flight Gross Weight (N)
real :: th0 !Flight Climb Angle
real :: phi0 !Flight Bank Angle
real :: th0_deg !Flight Climb Angle
real :: phi0_deg !Flight Bank Angle
real :: Sw !Wing Area (m^2)
real :: bw !Wing Span (m)
real :: ThrustOffset !Thrust Offset, positive below CG (m)
real :: ThrustSlope !Thrust Slope, dT/dV (N-sec/m)
real :: l_ref !Reference Length (m)
real :: V_ref !Reference Airspeed (m/s)
real :: MachNum_ref !Reference Mach Number
real :: rho_ref !Reference Air Density (kg/m^3)
real :: Lift_ref !Reference level flight lift (lbf)
real :: Ixx_ref !Reference Ixx (kg-m^2)
real :: Iyy_ref !Reference Iyy (kg-m^2)
real :: Izz_ref !Reference Izz (kg-m^2)
real :: Ixy_ref !Reference Ixy (kg-m^2)
real :: Ixz_ref !Reference Ixz (kg-m^2)
real :: Iyz_ref !Reference Iyz (kg-m^2)
real :: CD_ref !Reference drag coefficient
real :: CL_M_ref !Reference CL,M
real :: CD_M_ref !Reference CD,M
real :: Cm_M_ref !Reference Cm,M
real :: CLa_ref !Reference CL,alpha (lift)
real :: CDa_ref !Reference CD,alpha
real :: Clla_ref !Reference Cl,alpha (rolling moment)
real :: Cma_ref !Reference Cm,alpha
real :: Cna_ref !Reference Cn,alpha
real :: CLahat_ref !Reference CL,alphahat
real :: CDahat_ref !Reference CD,alphahat
real :: Cmahat_ref !Reference Cm,alphahat
real :: CLuhat_ref !Reference CL,muhat
real :: CDuhat_ref !Reference CD,muhat
real :: Cmuhat_ref !Reference Cm,muhat
real :: CYb_ref !Reference CY,beta
real :: Clb_ref !Reference Cl,beta
real :: Cmb_ref !Reference Cm,beta
real :: Cnb_ref !Reference Cn,beta
real :: CLqbar_ref !Reference CL,qbar
real :: CDqbar_ref !Reference CD,qbar
real :: Cmqbar_ref !Reference Cm,qbar
real :: CYpbar_ref !Reference CY,pbar
real :: Clpbar_ref !Reference Cl,pbar
real :: Cnpbar_ref !Reference Cn,pbar
real :: CYrbar_ref !Reference CY,rbar
real :: Clrbar_ref !Reference Cl,rbar
real :: Cnrbar_ref !Reference Cn,rbar
real :: hx_ref !Reference hx (kg-m^2/s)
real :: hy_ref !Reference hy (kg-m^2/s)
real :: hz_ref !Reference hz (kg-m^2/s)


!Transformed Variables
real :: Ixx
real :: Iyy
real :: Izz
real :: Ixy
real :: Ixz
real :: Iyz
real :: CZ0
real :: CX0
real :: Cm0
real :: CZu
real :: CXu
real :: Cmu
real :: CZa
real :: CXa
real :: Clla
real :: Cma
real :: Cna
real :: CZahat
real :: CXahat
real :: Cmahat
real :: CZuhat
real :: CXuhat
real :: Cmuhat
real :: CYb
real :: Clb
real :: Cmb
real :: Cnb
real :: CZqbar
real :: CXqbar
real :: Cmqbar
real :: CYpbar
real :: Clpbar
real :: Cnpbar
real :: CYrbar
real :: Clrbar
real :: Cnrbar
real :: hx
real :: hy
real :: hz

!Matix Variables
real :: txy, txz, tyx, tyz, tzx, tzy
real :: Bxu_, Bzu_, Bmu_, Bxa_, Bza_, Bma_
real :: nxx, nxy, nxz, nyx, nyy, nyz, nzx, nzy, nzz
real :: Ag
real :: Axu, Azu, Amu, Ayb, Alb, Anb
real :: Axa, Aza, Ama, Ala, Amb, Ana
real :: Ayps, Alps, Anps, Axqs, Azqs, Amqs, Ayrs, Alrs, Anrs


!Internal Variables
real :: theta_0 = 0.0
real :: DeltaAlpha
real :: cw
real :: CD
real :: CDa
real :: CL
real :: CL_ref
real :: tp
real :: Cm
real :: M
complex, dimension(:,:), allocatable :: EigenVec
character(100) :: filename = 'PredatorModMetric'!'Ex9_8_1'
character(4) :: ext_txt = '.txt'
character(100) :: temp ,fmt
integer :: ii,jj,i
    
    
allocate(Ar(N,N)) 
allocate(Br(N,N)) 
allocate(Eval(N)) 
allocate(V(N,N)) 
allocate(EigenVec(N,N)) 



call input_read !Reading Ex9_8_1.txt

CL_ref = W/(0.5*rho*V_ref**2*Sw) !(pg 936)
CL = (W*cos(theta_0))/(0.5*rho*V_0**2*Sw*cos(phi0)) !(9.8.23)
CD = CD_ref+((CDa_ref*CL_ref)/(2.0*CLa_ref))*(CL**2/CL_ref**2-1.0) !(9.8.26)
CDa = CDa_ref*(CL/CL_ref) !(9.8.27)
CXa = CL-CDa !(pg 937)
CZa = -CLa_ref-CD !(pg 937)
DeltaAlpha = ((CL-CL_ref)/CLa_ref)
tp = -DeltaAlpha !Transform phi


call transform !Transforming to body axis

call matix_coeff !Defining Constants for Matrices

call matrix_pop !Populating Matrices A and B

!write(*,'(a,f10.6)') 'tp', tp*180/pi
!write(*,'(a,f14.8)') 'txy: ',txz
!write(*,'(a,f14.8)') 'Ixy: ',Ixz
!write(*,'(a,f16.8)') 'Ixx: ',Ixx
!read(*,*) 

write(*,'(2/,a)') 'Matrix [A]'
do i = 1,12
    write(*,'(12(f8.4))') Ar(i,:)
end do

write(*,'(2/,a)') 'Matrix [B]'
do i = 1,12
    write(*,'(12(f8.4))') Br(i,:)
end do

!write(*,*) Ar(5,1)

write(*,*)
write(*,*) '''Enter'' to continue and see Eigen Values and Vectors'

read(*,*)

call eigensolver(IV,N,Ar,Br,EVAL,V,IERR)  

if(IERR /= 0) then
    write(*,*) 'IERR: ', IERR
    stop
    read(*,*)
end if

call order

write(*,'(2/,a)') '=======  Eigen Values  ======'

do ii = 1,12
    write(*,'(f14.8,a,f14.8)') real(EVAL(ii)), ',', imag(EVAL(ii))
end do

call organize


deallocate(Ar,Br,EVAL,V,EigenVec)



write(*,'(2/,a,/)') 'END OF PROGRAM'

read(*,*)

contains

!=======  Subroutines  =========
subroutine input_read
    implicit none


    character(1) :: y_or_n = 'n'

    write(*,'(a,a15,a)') 'File Name? (  ',filename,' )'
    read(*,'(a)') temp; 
    if(temp /= ' ') then; 
        read(temp,*) filename; 
    end if;
    filename = trim(filename) // trim(ext_txt)

    !write(*,*) filename
    !filename = 'PredatorMod.txt'

    OPEN(12, file=filename,status='OLD', action='READ') 

    read(12,*) V_0 !Flight Airspeed (ft/s)
    read(12,*) MachNum !Flight Mach Number
    read(12,*) rho !Air Density (kg/m^3)
    read(12,*) W !Flight Gross Weight (lbf)
    read(12,*) th0_deg !Flight Climb Angle (deg)
    read(12,*) phi0_deg !Flight Bank Angle (deg)
    read(12,*) Sw !Wing Area (ft^2)
    read(12,*) bw !Wing Span (m)
    read(12,*) ThrustOffset !Thrust Offset, positive below CG (m)
    read(12,*) ThrustSlope !Thrust Slope, dT/dV (lbf-sec/ft)
    read(12,*) l_ref !Reference Length (m)
    read(12,*) V_ref !Reference Airspeed (ft/s)
    read(12,*) MachNum_ref !Reference Mach Number
    read(12,*) rho_ref !Reference Air Density (kg/m^3)
    read(12,*) Lift_ref !Reference level flight lift (lbf)
    read(12,*) Ixx_ref !Reference Ixx (kg-m^2)
    read(12,*) Iyy_ref !Reference Iyy (kg-m^2)
    read(12,*) Izz_ref !Reference Izz (kg-m^2)
    read(12,*) Ixy_ref !Reference Ixy (kg-m^2)
    read(12,*) Ixz_ref !Reference Ixz (kg-m^2)
    read(12,*) Iyz_ref !Reference Iyz (kg-m^2)
    read(12,*) CD_ref !Reference drag coefficient
    read(12,*) CL_M_ref !Reference CL,M
    read(12,*) CD_M_ref !Reference CD,M
    read(12,*) Cm_M_ref !Reference Cm,M
    read(12,*) CLa_ref !Reference CL,alpha (lift)
    read(12,*) CDa_ref !Reference CD,alpha
    read(12,*) Clla_ref !Reference Cl,alpha (rolling moment)
    read(12,*) Cma_ref !Reference Cm,alpha
    read(12,*) Cna_ref !Reference Cn,alpha
    read(12,*) CLahat_ref !Reference CL,alphahat
    read(12,*) CDahat_ref !Reference CD,alphahat
    read(12,*) Cmahat_ref !Reference Cm,alphahat
    read(12,*) CLuhat_ref !Reference CL,muhat
    read(12,*) CDuhat_ref !Reference CD,muhat
    read(12,*) Cmuhat_ref !Reference Cm,muhat
    read(12,*) CYb_ref !Reference CY,beta
    read(12,*) Clb_ref !Reference Cl,beta
    read(12,*) Cmb_ref !Reference Cm,beta
    read(12,*) Cnb_ref !Reference Cn,beta
    read(12,*) CLqbar_ref !Reference CL,qbar
    read(12,*) CDqbar_ref !Reference CD,qbar
    read(12,*) Cmqbar_ref !Reference Cm,qbar
    read(12,*) CYpbar_ref !Reference CY,pbar
    read(12,*) Clpbar_ref !Reference Cl,pbar
    read(12,*) Cnpbar_ref !Reference Cn,pbar
    read(12,*) CYrbar_ref !Reference CY,rbar
    read(12,*) Clrbar_ref !Reference Cl,rbar
    read(12,*) Cnrbar_ref !Reference Cn,rbar
    read(12,*) hx_ref !Referenc hx (kg-m^2/s)
    read(12,*) hy_ref !Referenc hx (kg-m^2/s)
    read(12,*) hz_ref !Referenc hx (kg-m^2/s)
    close(12)

    cw = Sw/bw
    !M = MachNum

    !Converting from rad to deg
    th0 = th0_deg*pi/180.0
    phi0 = phi0_deg*pi/180.0

    write(*,'(a,a1,a)') 'View Input From File? ( ' ,y_or_n, ' )'
    read(*,'(a)') temp;
    if(temp /= ' ') then; 
        read(temp,*) y_or_n;  
    end if;

    if (y_or_n == 'y') then

        fmt = '(f14.8,4x,a)'
        write(*,fmt) V_0, 'Flight airspeed (m/s)' !Flight Airspeed (ft/s)
        write(*,fmt) MachNum, 'Flight Mach number' !Flight Mach Number
        write(*,fmt) rho, 'Flight air density (kg/m^3)' !Air Density (kg/m^3)
        write(*,fmt) W, 'Flight Gross weight (N)' !Flight Gross Weight (lbf)
        write(*,fmt) th0_deg, 'Flight climb angle (degrees)' !Flight Climb Angle
        write(*,fmt) phi0_deg, 'Flight bank angle (degrees)' !Flight Bank Angle
        write(*,fmt) Sw, 'Wing area (m^2)' !Wing Area (ft^2)
        write(*,fmt) bw, 'Wing span(m)' !Wing Span (m)
        write(*,fmt) ThrustOffset, 'Thrust offset, positive below CG (m)' !Thrust Offset, positive below CG (m)
        write(*,fmt) ThrustSlope, 'Thrust slope, dT/dV (N-sec/m)' !Thrust Slope, dT/dV (lbf-sec/ft)
        write(*,fmt) l_ref, 'Reference length (m)' !Reference Length (m)
        write(*,fmt) V_ref, 'Reference airspeed (m/s)' !Reference Airspeed (ft/s)
        write(*,fmt) MachNum_ref, 'Reference Mach number' !Reference Mach Number
        write(*,fmt) rho_ref, 'Reference air density (kg/m^3)' !Reference Air Density (kg/m^3)
        write(*,fmt) Lift_ref, 'Reference level flight lift (N)' !Reference level flight lift (lbf)
        write(*,fmt) Ixx_ref, 'Reference Ixx (kg-m^2)' !Reference Ixx (kg-m^2)
        write(*,fmt) Iyy_ref, 'Reference Iyy (kg-m^2)' !Reference Iyy (kg-m^2)
        write(*,fmt) Izz_ref, 'Reference Izz (kg-m^2)' !Reference Izz (kg-m^2)
        write(*,fmt) Ixy_ref, 'Reference Ixy (kg-m^2)' !Reference Ixy (kg-m^2)
        write(*,fmt) Ixz_ref, 'Reference Ixz (kg-m^2)' !Reference Ixz (kg-m^2)
        write(*,fmt) Iyz_ref, 'Reference Iyz (kg-m^2)' !Reference Iyz (kg-m^2)
        write(*,fmt) CD_ref, 'Reference drag coefficient' !Reference drag coefficient
        write(*,fmt) CL_M_ref, 'Reference CL,M' !Reference CL,M
        write(*,fmt) CD_M_ref, 'Reference CD,M' !Reference CD,M
        write(*,fmt) Cm_M_ref, 'Reference Cm,M' !Reference Cm,M
        write(*,fmt) CLa_ref, 'Reference CL,alpha' !Reference CL,alpha (lift)
        write(*,fmt) CDa_ref, 'Reference CD,alpha' !Reference CD,alpha
        write(*,fmt) Clla_ref, 'Reference Clla,alpha' !Reference Cl,alpha (rolling moment)
        write(*,fmt) Cma_ref, 'Reference Cm,alpha' !Reference Cm,alpha
        write(*,fmt) Cna_ref, 'Reference Cn,alpha' !Reference Cn,alpha
        write(*,fmt) CLahat_ref, 'Reference CL,alphahat' !Reference CL,alphahat
        write(*,fmt) CDahat_ref, 'Reference CD,alphahat' !Reference CD,alphahat
        write(*,fmt) Cmahat_ref, 'Reference Cm,alphahat' !Reference Cm,alphahat
        write(*,fmt) CLuhat_ref, 'Reference CL,uhat' !Reference CL,muhat
        write(*,fmt) CDuhat_ref, 'Reference CD,uhat' !Reference CD,muhat
        write(*,fmt) Cmuhat_ref, 'Reference Cm,uhat' !Reference Cm,muhat
        write(*,fmt) CYb_ref, 'Reference CY,beta' !Reference CY,beta
        write(*,fmt) Clb_ref, 'Reference Cl,beta' !Reference Cl,beta
        write(*,fmt) Cmb_ref, 'Reference Cm,beta' !Reference Cm,beta
        write(*,fmt) Cnb_ref, 'Reference Cn,beta' !Reference Cn,beta
        write(*,fmt) CLqbar_ref, 'Reference CL,qbar' !Reference CL,qbar
        write(*,fmt) CDqbar_ref, 'Reference CD,qbar' !Reference CD,qbar
        write(*,fmt) Cmqbar_ref, 'Reference Cm,qbar' !Reference Cm,qbar
        write(*,fmt) CYpbar_ref, 'Reference CY,pbar' !Reference CY,pbar
        write(*,fmt) Clpbar_ref, 'Reference Cl,pbar' !Reference Cl,pbar
        write(*,fmt) Cnpbar_ref, 'Reference Cn,pbar' !Reference Cn,pbar
        write(*,fmt) CYrbar_ref, 'Reference CY,rbar' !Reference CY,rbar
        write(*,fmt) Clrbar_ref, 'Reference Cl,rbar' !Reference Cl,rbar
        write(*,fmt) Cnrbar_ref, 'Reference Cn,rbar' !Reference Cn,rbar
        write(*,fmt) hx_ref, 'Reference hx (kg-m^2/sec)' !Referenc hx (kg-m^2/s)
        write(*,fmt) hy_ref, 'Reference hy (kg-m^2/sec)' !Referenc hx (kg-m^2/s)
        write(*,fmt) hz_ref, 'Reference hz (kg-m^2/sec)' !Referenc hx (kg-m^2/s)


        write(*,*) '''Enter'' to continue'
        read(*,*)
        
    end if

end subroutine input_read

subroutine transform
    implicit none


    CX0 = -CD
    CZ0 = -CL
    Cm0 = 0.0

    CXu = 0.0!-CD_M_ref*cos(tp)**2-(CDa_ref-CL_M_ref)*cos(tp)*sin(tp)-CLa_ref*sin(tp)**2
    CXa = CXa
    CXuhat = (-CDuhat_ref)*cos(tp)**2-(-CDahat_ref-CLuhat_ref)*cos(tp)*sin(tp)+(-CLahat_ref)*sin(tp)**2
    CXahat = (-CDahat_ref)*cos(tp)**2+(-CDuhat_ref+CLahat_ref)*cos(tp)*sin(tp)-(-CLuhat_ref)*sin(tp)**2
    CXqbar = (-CDqbar_ref)*cos(tp)-(-CLqbar_ref)*sin(tp)
    !CXdebar = ()*cos(tp)-()*sin(tp)
    CZu = 0.0! (-CL_M_ref)*cos(tp)**2+(-CD_M_ref+CLa_ref)*cos(tp)*sin(tp)-(-CDa_ref)*sin(tp)**2
    CZa = CZa
    CZuhat = (-CLuhat_ref)*cos(tp)**2+(-CDuhat_ref+CLahat_ref)*cos(tp)*sin(tp)-(-CDahat_ref)*sin(tp)**2
    CZahat = (-CLahat_ref)*cos(tp)**2+(-CDahat_ref-CLuhat_ref)*cos(tp)*sin(tp)+(-CDahat_ref)*sin(tp)**2
    CZqbar = (-CLqbar_ref)*cos(tp)-(-CDqbar_ref)*sin(tp)
    !CZdebar = ()*cos(tp)-()*sin(tp)
    Cmu = 0.0!(Cm_M_ref)*cos(tp)-(Cma_ref)*sin(tp)
    Cma = (Cma_ref)*cos(tp)+(Cm_M_ref)*sin(tp)
    Cmuhat = (Cmuhat_ref)*cos(tp)-(Cmahat_ref)*sin(tp)
    Cmahat = (Cmahat_ref)*cos(tp)+(Cmuhat_ref)*sin(tp)
    Cmqbar = Cmqbar_ref
    !Cmdebar = ()*cos(tp)-()*sin(tp)
    CYb = CYb_ref
    CYpbar = (CYpbar_ref)*cos(tp)-(CYrbar_ref)*sin(tp)
    CYrbar = (CYrbar_ref)*cos(tp)+(CYpbar_ref)*sin(tp)
    !Cyda = 
    !Cydr = 
    Clb = (Clb_ref)*cos(tp)-(Cnb_ref)*sin(tp)
    Clpbar = (Clpbar_ref)*cos(tp)**2-(Clrbar_ref+Cnpbar_ref)*cos(tp)*sin(tp)+(Cnrbar_ref)*sin(tp)**2
    Clrbar = (Clrbar_ref)*cos(tp)**2+(Clpbar_ref-Cnrbar_ref)*cos(tp)*sin(tp)-(Cnpbar_ref)*sin(tp)**2
    !Clda
    !Cldr
    Cnb = (Cnb_ref)*cos(tp)+(Clb_ref)*sin(tp)
    Cnpbar = (Cnpbar_ref)*cos(tp)**2+(Clpbar_ref-Cnrbar_ref)*cos(tp)*sin(tp)-(Clrbar_ref)*sin(tp)**2
    Cnrbar = (Cnrbar_ref)*cos(tp)**2+(Clrbar_ref+Cnpbar_ref)*cos(tp)*sin(tp)+(Clpbar_ref)*sin(tp)**2
    !Cnda
    !Cndr

    Ixx = (Ixx_ref)*cos(tp)**2+(2.0*Ixz_ref)*cos(tp)*sin(tp)+(Izz_ref)*sin(tp)**2
    Iyy = Iyy_ref
    Izz = (Izz_ref)*cos(tp)**2-(2.0*Ixz_ref)*cos(tp)*sin(tp)+(Ixx_ref)*sin(tp)**2
    Ixy = 0.0
    Ixz = (Ixz_ref)*(cos(tp)**2-sin(tp)**2)+(Izz_ref-Ixx_ref)*cos(tp)*sin(tp)
    Iyz = 0.0
    hx = hx_ref
    hy = hy_ref
    hz = hz_ref

    Clla = 0.0
    Cna = 0.0
    Cmb = 0.0


    !write(*,'(/,3(a10,f6.0))') 'Ixx = ', Ixx,'Iyy = ', Iyy,'Izz = ', Izz
    !write(*,'(3(a10,f6.0))') 'Ixy = ', Ixy,'Ixz = ', Ixz,'Iyz = ', Iyz
    !write(*,'(3(a10,f6.0))') 'hx = ', hx,'hy = ', hy,'hz = ', hz
    !write(*,'(3(a10,f6.3))') 'CX0 = ', CX0,'CZ0 = ', CZ0,'Cm0 = ', Cm0
    !write(*,'(3(a10,f6.3))') 'CXuhat = ', CXuhat,'CZuhat = ', CZuhat,'Cmuhat = ', Cmuhat
    !write(*,'(3(a10,f6.3))') 'CXahat = ', CXahat,'CZahat = ', CZahat,'Cmahat = ', Cmahat
    !write(*,'(3(a10,f6.3))') 'CXu = ', CXu,'CZu = ', CZu,'Cmu = ', Cmu
    !write(*,'(5(a10,f6.3))') 'CXa = ', CXa,'CZa = ', CZa,'Cla = ', Clla,'Cma = ', Cma,'Cna = ', Cna
    !write(*,'(4(a10,f7.4))') 'CYb = ', CYb,'Clb = ', Clb,'Cmb = ', Cmb,'Cnb = ', Cnb
    !write(*,'(3(a10,f7.4))') 'CYpbar = ', CYpbar,'Clpbar = ', Clpbar,'Cnpbar = ', Cnpbar
    !write(*,'(3(a10,f8.3))') 'CXqbar = ', CXqbar,'CZqbar = ', CZqbar,'Cmqbar = ', Cmqbar
    !write(*,'(3(a10,f7.3))') 'CYrbar = ', CYrbar,'Clrbar = ', Clrbar,'Cnrbar = ', Cnrbar

    !write(*,'(a,/)') 'Transformed Values'
    !
    !fmt = '(3(a12,f10.4))'
    !write(*,fmt) 'Ixx = ', Ixx,'Iyy = ', Iyy,'Izz = ', Izz
    !write(*,fmt) 'Ixy = ', Ixy,'Ixz = ', Ixz,'Iyz = ', Iyz
    !write(*,fmt) 'hx = ', hx,'hy = ', hy,'hz = ', hz
    !write(*,fmt) 'CX0 = ', CX0,'CZ0 = ', CZ0,'Cm0 = ', Cm0
    !write(*,fmt) 'CXuhat = ', CXuhat,'CZuhat = ', CZuhat,'Cmuhat = ', Cmuhat
    !write(*,fmt) 'CXahat = ', CXahat,'CZahat = ', CZahat,'Cmahat = ', Cmahat
    !write(*,fmt) 'CXu = ', CXu,'CZu = ', CZu,'Cmu = ', Cmu
    !write(*,fmt) 'CXa = ', CXa,'CZa = ', CZa,'Cma = ', Cma
    !write(*,fmt) 'CYb = ', CYb,'Clb = ', Clb,'Cnb = ', Cnb
    !write(*,fmt) 'Cla = ', Clla,'Cmb = ', Cmb,'Cna = ', Cna
    !write(*,fmt) 'CYpbar = ', CYpbar,'Clpbar = ', Clpbar,'Cnpbar = ', Cnpbar
    !write(*,fmt) 'CXqbar = ', CXqbar,'CZqbar = ', CZqbar,'Cmqbar = ', Cmqbar
    !write(*,fmt) 'CYrbar = ', CYrbar,'Clrbar = ', Clrbar,'Cnrbar = ', Cnrbar


    write(*,'(/,3(a10,es18.10))') 'Ixx = ', Ixx,'Iyy = ', Iyy,'Izz = ', Izz
    write(*,'(3(a10,es18.10))') 'Ixy = ', Ixy,'Ixz = ', Ixz,'Iyz = ', Iyz
    write(*,'(3(a10,es18.10))') 'hx = ', hx,'hy = ', hy,'hz = ', hz
    write(*,'(3(a10,es18.10))') 'CX0 = ', CX0,'CZ0 = ', CZ0,'Cm0 = ', Cm0
    write(*,'(3(a10,es18.10))') 'CXuhat = ', CXuhat,'CZuhat = ', CZuhat,'Cmuhat = ', Cmuhat
    write(*,'(3(a10,es18.10))') 'CXahat = ', CXahat,'CZahat = ', CZahat,'Cmahat = ', Cmahat
    write(*,'(3(a10,es18.10))') 'CXu = ', CXu,'CZu = ', CZu,'Cmu = ', Cmu
    write(*,'(3(a10,es18.10))') 'CXa = ', CXa,'CZa = ', CZa,'Cma = ', Cma
    write(*,'(3(a10,es18.10))') 'CYb = ', CYb,'Clb = ', Clb,'Cnb = ', Cnb
    write(*,'(3(a10,es18.10))') 'Cla = ', Clla,'Cmb = ', Cmb,'Cna = ', Cna
    write(*,'(3(a10,es18.10))') 'CYpbar = ', CYpbar,'Clpbar = ', Clpbar,'Cnpbar = ', Cnpbar
    write(*,'(3(a10,es18.10))') 'CXqbar = ', CXqbar,'CZqbar = ', CZqbar,'Cmqbar = ', Cmqbar
    write(*,'(3(a10,es18.10))') 'CYrbar = ', CYrbar,'Clrbar = ', Clrbar,'Cnrbar = ', Cnrbar



    !write(*,'(a15,f13.7)') 'dummy = ', dummy
    write(*,*)
    write(*,*) '''Enter'' to continue and see [A] and [B] matrices'
    read(*,*)

end subroutine transform

subroutine matix_coeff
    implicit none


    txy = Ixy/Ixx
    txz = Ixz/Ixx
    tyx = Ixy/Iyy
    tyz = Iyz/Iyy
    tzx = Ixz/Izz
    tzy = Iyz/Izz
    Bxu_ = ((rho*Sw*cw)/((4.0*W)/g))*CXuhat
    Bzu_ = ((rho*Sw*cw)/((4.0*W)/g))*CZuhat
    Bmu_ = (rho*Sw*cw**2*l_ref)/(4.0*Iyy)*Cmuhat
    Bxa_ = ((rho*Sw*cw)/((4.0*W)/g))*CXahat
    Bza_ = ((rho*Sw*cw)/((4.0*W)/g))*CZahat
    Bma_ = (rho*Sw*cw**2*l_ref)/(4.0*Iyy)*Cmahat
    Ag = (g*l_ref)/V_0**2
    nxx = Ag*(Ixz*tan(phi0)*sin(phi0)*cos(th0)-Ixy*sin(phi0)*cos(th0))/Ixx
    nxy = (hz*l_ref)/(Ixx*V_0)+Ag*((Izz-Iyy)*sin(phi0)*cos(th0)-2.0*Iyz*tan(phi0)*sin(phi0)*cos(th0)+Ixz*tan(phi0)*sin(th0))/Ixx
    nxz = (hy*l_ref)/(Ixx*V_0)+Ag*((Iyy-Izz)*tan(phi0)*sin(phi0)*cos(th0)-2.0*Iyz*sin(phi0)*cos(th0)+Ixz*tan(phi0)*sin(th0))/Ixx
    nyx = (hz*l_ref)/(Iyy*V_0)+Ag*((Izz-Ixx)*sin(phi0)*cos(th0)+2.0*Iyz*tan(phi0)*sin(phi0)-Iyz*tan(phi0)*sin(phi0)*cos(th0))/Iyy
    nyy = Ag*(Ixy*sin(phi0)*cos(th0)+Iyz*tan(phi0)*sin(th0))/Iyy
    nyz = (hx*l_ref)/(Iyy*V_0)+Ag*((Izz-Ixx)*tan(phi0)*sin(th0)-2.0*Ixz*sin(phi0)*cos(th0)-Ixy*tan(phi0)*sin(phi0)*cos(th0))/Iyy
    nzx = (hy*l_ref)/(Izz*V_0)+Ag*((Iyy-Ixx)*tan(phi0)*sin(phi0)*cos(th0)+2.0*Ixy*tan(phi0)*sin(th0)-Iyz*sin(phi0)*cos(th0))/Izz
    nzy = (hx*l_ref)/(Izz*V_0)+Ag*((Iyy-Ixx)*tan(phi0)*sin(th0)-2.0*Ixy*tan(phi0)*sin(phi0)*cos(th0)-Ixz*sin(phi0)*cos(th0))/Izz
    nzz = Ag*(-Iyz*tan(phi0)*sin(phi0)-Ixz*tan(phi0)*sin(phi0)*cos(th0))/Izz
    Axu = ((rho*Sw*l_ref)/((2.0*W)/g))*(-2.0*CD-M*CD_M_ref) !(pg 795)
    Azu = ((rho*Sw*l_ref)/((2.0*W)/g))*(-2.0*CL-M*CD_M_ref) !(pg 795)
    Amu = ((rho*Sw*cw*l_ref**2)/(2.0*Iyy))*(2.0-M**2)/(1.0-M**2)*Cm0 !(pg 795)
    Axa = (rho*Sw*l_ref)/((2.0*W)/g)*CXa
    Aza = (rho*Sw*l_ref)/((2.0*W)/g)*CZa
    Ama = (rho*Sw*cw*l_ref**2)/(2.0*Iyy)*Cma 
    Ayb = (rho*Sw*l_ref)/((2.0*W)/g)*CYb
    Alb = (rho*Sw*bw*l_ref**2)/(2.0*Ixx)*Clb
    Anb = (rho*Sw*bw*l_ref**2)/(2.0*Izz)*Cnb
    Ala = (rho*Sw*bw*l_ref**2)/(2.0*Ixx)*Clla
    Amb = (rho*Sw*cw*l_ref**2)/(2.0*Iyy)*Cmb
    Ana = (rho*Sw*bw*l_ref**2)/(2.0*Izz)*Cna
    Ayps = (rho*Sw*bw)/((4.0*W)/g)*CYpbar
    Alps = (rho*Sw*bw**2*l_ref)/(4.0*Ixx)*Clpbar
    Anps = (rho*Sw*bw**2*l_ref)/(4.0*Izz)*Cnpbar
    Axqs = (rho*Sw*cw)/((4.0*W)/g)*CXqbar
    Azqs = (rho*Sw*cw)/((4.0*W)/g)*CZqbar
    Amqs = (rho*Sw*cw**2*l_ref)/(4.0*Iyy)*Cmqbar
    Ayrs = (rho*Sw*bw)/((4.0*W)/g)*CYrbar
    Alrs = (rho*Sw*bw**2*l_ref)/(4.0*Ixx)*Clrbar
    Anrs = (rho*Sw*bw**2*l_ref)/(4.0*Izz)*Cnrbar





end subroutine matix_coeff

subroutine matrix_pop
    implicit none

    Ar = 0.0
    Ar(1,1) = Axu
    Ar(1,2) = Ag*sin(phi0)*cos(th0)
    Ar(1,3) = Axa-Ag*tan(phi0)*sin(phi0)*cos(th0)
    Ar(1,5) = Axqs
    Ar(1,11) = -Ag*cos(th0)
    Ar(2,1) = -Ag*sin(phi0)*cos(th0)
    Ar(2,2) = Ayb
    Ar(2,3) = -Ag*tan(phi0)*sin(th0)
    Ar(2,4) = Ayps
    Ar(2,6) = Ayrs-1.0
    Ar(2,10) = Ag*cos(phi0)*cos(th0)
    Ar(2,11) = -Ag*sin(phi0)*sin(th0)
    Ar(3,1) = Azu+Ag*tan(phi0)*sin(phi0)*cos(th0)
    Ar(3,2) = Ag*tan(phi0)*sin(th0)
    Ar(3,3) = Aza
    Ar(3,5) = Azqs+1.0
    Ar(3,10) = -Ag*sin(phi0)*cos(th0)
    Ar(3,11) = -Ag*cos(phi0)*sin(th0)
    Ar(4,2) = Alb
    Ar(4,3) = Ala
    Ar(4,4) = Alps+nxx
    Ar(4,5) = -nxy
    Ar(4,6) = Alrs+nxz
    Ar(5,1) = Amu
    Ar(5,2) = Amb
    Ar(5,3) = Ama
    Ar(5,4) = nyx
    Ar(5,5) = Amqs+nyy
    Ar(5,6) = -nyz
    Ar(6,2) = Anb
    Ar(6,3) = Ana
    Ar(6,4) = Anps-nzx
    Ar(6,5) = nzy
    Ar(6,6) = Anrs+nzz
    Ar(7,1) = cos(th0)
    Ar(7,2) = sin(phi0)*sin(th0)
    Ar(7,3) = cos(phi0)*sin(th0)
    Ar(7,11) = -sin(th0)
    Ar(8,2) = cos(phi0)
    Ar(8,3) = -sin(phi0)
    Ar(8,12) = cos(th0)
    Ar(9,1) = -sin(th0)
    Ar(9,2) = sin(phi0)*cos(th0)
    Ar(9,3) = cos(phi0)*cos(th0)
    Ar(9,11) = -cos(th0)
    Ar(10,4) = 1.0
    Ar(10,5) = sin(phi0)*tan(th0)
    Ar(10,6) = cos(phi0)*tan(th0)
    Ar(10,11) = Ag*tan(phi0)/cos(th0)
    Ar(11,5) = cos(phi0)
    Ar(11,6) = -sin(phi0)
    Ar(11,10) = -Ag*tan(phi0)*cos(th0)
    Ar(12,5) = sin(phi0)/cos(th0)
    Ar(12,6) = cos(phi0)/cos(th0)
    Ar(12,11) = Ag*tan(phi0)*tan(th0)

    Br(1,:)  = (/ 1.0-Bxu_, 0.0, -Bxa_,     0.0, 0.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
    Br(2,:)  = (/ 0.0,      1.0, 0.0,       0.0, 0.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
    Br(3,:)  = (/ -Bzu_,    0.0, 1.0-Bza_,  0.0, 0.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
    Br(4,:)  = (/ 0.0,      0.0, 0.0,       1.0, -txy,  -txz, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
    Br(5,:)  = (/ -Bmu_,    0.0, -Bma_,     -tyx, 1.0,  -tyz, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
    Br(6,:)  = (/ 0.0,      0.0, 0.0,       -tzx, -tzy,  1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
    Br(7,:)  = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
    Br(8,:)  = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0/)
    Br(9,:)  = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0/)
    Br(10,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0/)
    Br(11,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0/)
    Br(12,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0/)

end subroutine matrix_pop

subroutine order
    implicit none

    complex, dimension(12) :: EVAL_temp = 0.0
    complex, dimension(12,12) :: EVEC_temp = 0.0
    complex, dimension(12) :: EVEC_switch = 0.0
    complex :: EVAL_switch = 0.0
    integer :: num_real = 0
    integer :: i, j, k, l, m, n
    j = 1
    i = 1
    k = 1
    l = 1
    n = 1

    do m = 1, 12
        if (abs(imag(EVAL(m))) < 0.0000001 .and. &
            & abs(real(EVAL(m))) > 0.0000001) then
            num_real = num_real + 1
        end if
    end do
        
    do j = 1, 12
    !    write(*,*) 'i:', i
    !    write(*,*) 'j:', j
    !    write(*,*) 'k:', k

        if (abs(real(EVAL(j))) < 0.0000001)then
    !        write(*,*) 'A'
            EVAL_temp(l+8) = EVAL(j)
            EVEC_temp(:,l+8) = V(:,j) 
            l = l + 1 
        else if(abs(imag(EVAL(j))) < 0.000001)then
            ! Real
    !        write(*,*) 'B'
            EVAL_temp(i) = EVAL(j)
            EVEC_temp(:,i) = V(:,j) 
            i = i + 1 
        else
            ! Complex
    !        write(*,*) 'C'
            EVAL_temp(k+num_real) = EVAL(j)
            EVEC_temp(:,k+num_real) = V(:,j)  
            k = k + 1
        end if   
            n = 1
        
            
    end do


    !do k = 1,12
    !    write(*,'(f14.8,a,f14.8)') real(EVAL_temp(k)), ',', imag(EVAL_temp(k))
    !end do

    if (real(EVEC_temp(2,7)) < 0.00000001) then
        EVEC_switch = EVEC_temp(:,7)
        EVEC_temp(:,7) = EVEC_temp(:,5)
        EVEC_temp(:,5) = EVEC_switch
        EVEC_switch = EVEC_temp(:,8)
        EVEC_temp(:,8) = EVEC_temp(:,6)
        EVEC_temp(:,6) = EVEC_switch
        
        EVAL_switch = EVAL_temp(7)
        EVAL_temp(7) = EVAL_temp(5)
        EVAL_temp(5) = EVAL_switch
        EVAL_switch = EVAL_temp(8)
        EVAL_temp(8) = EVAL_temp(6)
        EVAL_temp(6) = EVAL_switch
    end if

    !write(*,*) EVEC_temp(2,7)
    !write(*,*) EVEC_temp(2,5)
    !write(*,*) EVEC_switch

    EVAL = EVAL_temp
    V = EVEC_temp


end subroutine

subroutine organize
    implicit none

    real, dimension(12) :: y
    real, dimension(12,12) :: amp, phase
    !complex, dimension(12,12) :: EVAL_Temp
    real, dimension(1) :: x
    integer :: loc, i, j, m, mark, mark2
    real :: converge = 0.000001 !Convergence Criteria
    complex :: conjugate
    real :: mag


    do i=1,N
        !Separating Vectors
        EigenVec(i,:) = V(:,i)
        !Finding Max Magnitude component of vector
        y = cabs(EigenVec(i,:))
        x = maxloc(y)
        loc = int(x(1))
        !Setting the complex conjugate
        conjugate = conjg(EigenVec(i,loc))
        !Normalizing by complex conjugate of Max Magnitude component of vector
        EigenVec(i,:) = EigenVec(i,:) * conjugate
        mag = sqrt(sum(cabs(EigenVec(i,:))**2))
        EigenVec(i,:) = EigenVec(i,:)/mag
        
        
        
        
        do jj=1,N
            amp(i,jj) = cabs(EigenVec(i,jj)); !Magnitude Array
            
            if (amp(i,jj) < 0.0000001) then
                phase(i,jj) = 0.0
            else
                phase(i,jj) = atan2(imag(EigenVec(i,jj)),real(EigenVec(i,jj)))*180.0/pi !Phase Array
            end if
            
            if (amp(i,jj) < 0.0000001) then
                phase(i,jj) = 0.0
            end if
            
            
        end do
        
        
    end do

    !write(*,*)
    !write(*,*) phase
    !write(*,*)

    m = 1
    do j=1,N
        
        if (abs(real(EVAL(m+1))) > abs(real(EVAL(m))) - 0.00000001 &
            & .and. abs(real(EVAL(m+1))) < abs(real(EVAL(m))) + 0.00000001) then
            write(*,'(/,a,f13.7,a,f13.7,a)') '=====  EIGEN VALUE: ', real(EVAL(m)), ' +/-',abs(aimag(EVAL(m))), '  ====='
            m = m + 2
            mark = 1
        else
            write(*,'(/,a,f13.7,a)') '=====  EIGEN VALUE: ', real(EVAL(m)), '  ====='
            m = m + 1
            mark = 2
        end if
            
    !    write(*,*) m
    !    read(*,*)
        

        if (abs(aimag(EVAL(m-1))) < converge) then !checking if complex
            !REAL Eigen Value (Non-oscillatory mode)
            if (abs(real(EVAL(m-1))) < converge) then
                !Eigen Value is Zero (Rigid Body Mode)
                write(*,*) '    (Rigid Body Mode) '
                mark2 = 1
                
            else if (real(EVAL(m-1)) < 0.0) then
                !Eigen Value is Negative (Non-oscillatory Convergent Mode)
                write(*,*) '    (Non-oscillatory Convergent Mode)'
                write(*,'(a30,f13.7)') 'Damping Rate (1/s) = ', DampRate(EVAL(m-1),V_0,l_ref)
                write(*,'(a30,f13.7)') '99% Damping Time (s) = ', DampTime99(EVAL(m-1),V_0,l_ref)

            else
                !Eigen Value is Positive (Non-oscillatory Divergent Mode)
                write(*,*) '    (Non-oscillatory Divergent Mode)'
                write(*,'(a30,f13.7)') 'Damping Rate (1/s) = ', DampRate(EVAL(m-1),V_0,l_ref)
                write(*,'(a30,f13.7)') 'Doubling Time = ', DoubleTime(EVAL(m-1),V_0,cw)
                write(*,'(a30,f13.7)') 'Damping Ratio = ', DampRatio(EVAL(m-2),EVAL(m-1))
                write(*,'(a30,f13.7)') 'Damped Natural Freq. (1/s) = ', DampNrulFreq(EVAL(m-2),V_0,l_ref)
                
            end if
            
        else
            !Complex Eigen Value (Oscillatory mode)
            if (abs(real(EVAL(m-1))) < converge) then
                !Eigen Value is Zero (Undamped Mode)
                write(*,*) '    (Oscillatory Undamped Mode)'
                
            else if (real(EVAL(m-1)) < 0.0) then
                !Eigen Value is Negative (Damped Mode)
                write(*,*) '    (Oscillatory Damped Mode)'
                write(*,'(a30,f13.7)') 'Damping Rate (1/s) = ', DampRate(EVAL(m-1),V_0,l_ref)
                write(*,'(a30,f13.7)') '99% Damping Time (s) = ', DampTime99(EVAL(m-1),V_0,l_ref)
                write(*,'(a30,f13.7)') 'Damping Ratio = ', DampRatio(EVAL(m-2),EVAL(m-1))
                write(*,'(a30,f13.7)') 'Damped Natural Freq. (1/s) = ', DampNrulFreq(EVAL(m-1),V_0,l_ref)
                write(*,'(a30,f13.7)') 'Period (s) = ', Period(EVAL(m-1),V_0,l_ref)

            else
                !Eigen Value is Positive (Divergent Mode)
                write(*,*) '    (Oscillatory Divergent Mode)'
                write(*,'(a24,f13.7)') 'Damping Rate = ', DampRate(EVAL(m-1),V_0,l_ref)
                write(*,'(a24,f13.7)') 'Doubling Time = ', DoubleTime(EVAL(m-1),V_0,l_ref)
                write(*,'(a24,f13.7)') 'Damping Ratio = ', DampRatio(EVAL(m-2),EVAL(m-1))
                write(*,'(a24,f13.7)') 'Damped Natural Freq. = ', DampNrulFreq(EVAL(m-1),V_0,l_ref)
                write(*,'(a24,f13.7)') 'Period = ', Period(EVAL(m-1),V_0,l_ref)
                
                
            end if
            
        end if
        
        
    !    write(*,*) 'm = ',m
        if (abs(EVAL(m-1)) > 0.00000001 .and. mark2 /= 1) then
            if (mark == 2) then
                !Real only
                fmt = "(a12,f10.6)"
                write(*,*)
                write(*,'(a21)')'   Eigen Vector'
                write(*,fmt) 'u/Vo = ',real(EigenVec(m-1,1))
                write(*,fmt) 'Beta = ',real(EigenVec(m-1,2))
                write(*,fmt) 'Alpha = ',real(EigenVec(m-1,3))
                write(*,fmt) 'p*lr/Vo = ',real(EigenVec(m-1,4))
                write(*,fmt) 'q*lr/Vo = ',real(EigenVec(m-1,5))
                write(*,fmt) 'r*lr/Vo = ',real(EigenVec(m-1,6))
                write(*,fmt) 'xc/lr = ',real(EigenVec(m-1,7))
                write(*,fmt) 'yc/lr = ',real(EigenVec(m-1,8))
                write(*,fmt) 'zc/lr = ',real(EigenVec(m-1,9))
                write(*,fmt) 'Phi = ',real(EigenVec(m-1,10))
                write(*,fmt) 'Theta = ',real(EigenVec(m-1,11))
                write(*,fmt) 'Psi = ',real(EigenVec(m-1,12))
                write(*,*)
            else
                !Imaginary
                write(*,*)
                fmt = "(a8,f10.6,f10.6,a2,f10.6,f8.2,a4)"
                write(*,'(a26,a14,a10)')'   Eigen Vector','Amplitude','Phase'
                write(*,fmt) 'u/Vo = ',real(EigenVec(m-1,1)), imag(EigenVec(m-1,1)),'i ', amp(m-1,1), phase(m-1,1),'deg'
                write(*,fmt) 'Beta = ',real(EigenVec(m-1,2)), imag(EigenVec(m-1,2)),'i ', amp(m-1,2), phase(m-1,2),'deg'
                write(*,fmt) 'Alpha = ',real(EigenVec(m-1,3)), imag(EigenVec(m-1,3)),'i ', amp(m-1,3), phase(m-1,3),'deg'
                write(*,fmt) 'p*lr/Vo = ',real(EigenVec(m-1,4)), imag(EigenVec(m-1,4)),'i ', amp(m-1,4), phase(m-1,4),'deg'
                write(*,fmt) 'q*lr/Vo = ',real(EigenVec(m-1,5)), imag(EigenVec(m-1,5)),'i ', amp(m-1,5), phase(m-1,5),'deg'
                write(*,fmt) 'r*lr/Vo = ',real(EigenVec(m-1,6)), imag(EigenVec(m-1,6)),'i ', amp(m-1,6), phase(m-1,6),'deg'
                write(*,fmt) 'xc/lr = ',real(EigenVec(m-1,7)), imag(EigenVec(m-1,7)),'i ', amp(m-1,7), phase(m-1,7),'deg'
                write(*,fmt) 'yc/lr = ',real(EigenVec(m-1,8)), imag(EigenVec(m-1,8)),'i ', amp(m-1,8), phase(m-1,8),'deg'
                write(*,fmt) 'zc/lr = ',real(EigenVec(m-1,9)), imag(EigenVec(m-1,9)),'i ', amp(m-1,9), phase(m-1,9),'deg'
                write(*,fmt) 'Phi = ',real(EigenVec(m-1,10)), imag(EigenVec(m-1,10)),'i ', amp(m-1,10), phase(m-1,10),'deg'
                write(*,fmt) 'Theta = ',real(EigenVec(m-1,11)), imag(EigenVec(m-1,11)),'i ', amp(m-1,11), phase(m-1,11),'deg'
                write(*,fmt) 'Psi = ',real(EigenVec(m-1,12)), imag(EigenVec(m-1,12)),'i ', amp(m-1,12), phase(m-1,12),'deg'
                write(*,*)
                
            end if
        end if



        if (m >= N) then
            exit
        end if
        
    !    read(*,*)
        
    end do


end subroutine organize

!=======  Functions  =========
function DampRate(E,V_0,l_ref)
    real :: DampRate, V_0, l_ref
    complex :: E
    DampRate = -(real(E)*V_0)/l_ref
end function DampRate

function DampTime99(E,V_0,l_ref)
    real :: DampTime99, V_0, l_ref
    complex :: E
    DampTime99 = (log(0.01)/((real(E)*V_0)/l_ref))
end function DampTime99

function DampRatio(E1,E2)
    real :: DampRatio
    complex :: E1, E2
    DampRatio = -((E1)+(E2))/(2.0*sqrt((E1)*(E2)))
end function DampRatio

function DampNrulFreq(E,V_0,l_ref)
    real :: DampNrulFreq, V_0, l_ref
    complex :: E
    DampNrulFreq = abs(real(aimag(E)))*(V_0/l_ref)
end function DampNrulFreq

function NrulFreq(E1,E2,V_0,l_ref)
    real :: NrulFreq, V_0, l_ref
    complex :: E1, E2
    NrulFreq = real(sqrt(E1*E2)*(V_0)/l_ref)
end function NrulFreq

function Period(E,V_0,l_ref)
    real :: Period, V_0, l_ref
    complex :: E 
    Period = (2.0*pi)/DampNrulFreq(E,V_0,l_ref)
end function Period

function DoubleTime(E,V_0,l_ref)
    real :: DoubleTime, V_0, l_ref
    complex :: E
    DoubleTime = (log(2.0)/((real(E)*V_0)/l_ref))
end function DoubleTime



end program project
