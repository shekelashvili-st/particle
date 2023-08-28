Program Main
  implicit none
  character(*), parameter:: InputFile='input.txt',OutputFile='data.plt',sol_file='fields.plt', &
							input_particle='params.nml'							                ! names of input and output files
  character MeshFile*30       ! name of file with computational mesh
  integer::	i,j,ni,nj,m,nm,sch=3,sols = 0, iu
  integer, parameter:: IO = 12 ! input-output unit
  real :: rtmp, pi=4*atan(1.0_8)
  real,allocatable,dimension(:,:):: X,Y,P,CellVolume,DivV,DivV_t,DivV_res,lapP,lapP_t,lapP_res, &
&									DivVP,DivVP_t,DivVP_res,CurlV,CurlV_t,CurlV_res  	 ! scalar arrays
  real,allocatable,dimension(:,:,:):: CellCenter,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector,&
&									  GradP,GradP_t,GradP_res,V  ! vector arrays
!Решение уравнения переноса
  real:: rho_f, rho_p, mu, d_p, r_0(3), v_0(3), r(3), v_p(3), F(3), U(2)
  real:: J_p, m_p, dt, t
  integer:: nt, Ip, Jp, St
  character(:),allocatable:: frmt_screen, frmt_file
  character(*), parameter::  particle_out='traj.dat'

!===  READ INPUT FILE ===
  namelist /params_particles/ rho_f, rho_p, mu, d_p, J_p, r_0, v_0, dt, nt
  open(newunit=iu, file=input_particle)
  read(iu, nml=params_particles)
  close(iu)
  
  WRITE(*,*) 'Read input file: ', InputFile
  OPEN(IO,FILE=InputFile)
  READ(IO,*) MeshFile  ! read name of file with computational mesh
  READ(IO,*) sch 	   ! 1 - linear, 2 - FOU, 3 - SOU
  READ(IO,*) sols 	   ! 0 - initialize fields from func, 1 - initialize from solution file
  CLOSE(IO)

!===   READ NODES NUMBER (NI,NJ) FROM FILE WITH MESH ===
  WRITE(*,*) 'Read nodes number from file: ', MeshFile
  OPEN(IO,FILE = MeshFile)
  READ(IO,*) NI,NJ, rtmp
  WRITE(*,*) 'NI, NJ = ',NI,NJ

!=== ALLOCATE ALL ARRAYS ===
  WRITE(*,*) 'Allocate arrays'       
  allocate(X(NI,NJ)) ! mesh nodes X-coordinates
  allocate(Y(NI,NJ)) ! mesh nodes Y-coordinates
  allocate(P(0:NI,0:NJ))   ! Pressure
  allocate(CellVolume(NI-1,NJ-1))   ! Cell Volumes    
  allocate(CellCenter(0:NI,0:NJ,2)) ! Cell Centers
  allocate(IFaceCenter( NI,NJ-1,2)) ! Face Centers for I-faces
  allocate(IFaceVector( NI,NJ-1,2)) ! Face Vectors for I-faces
  allocate(JFaceCenter( NI-1,NJ,2)) ! Face Centers for J-faces
  allocate(JFaceVector( NI-1,NJ,2)) ! Face Vectors for J-faces
  allocate(GradP(0:NI,0:NJ,2),GradP_t(0:NI,0:NJ,2),GradP_res(0:NI,0:NJ,2))  		! Pressure gradients array
  allocate(V(0:NI,0:NJ,2),DivV(0:NI,0:NJ),DivV_t(0:NI,0:NJ),DivV_res(0:NI,0:NJ), &
&			             DivVP(0:NI,0:NJ),DivVP_t(0:NI,0:NJ),DivVP_res(0:NI,0:NJ))	! Velocity vector and divergence arrays
  allocate(lapP(0:NI,0:NJ),lapP_t(0:NI,0:NJ),lapP_res(0:NI,0:NJ))
  allocate(CurlV(0:NI,0:NJ),CurlV_t(0:NI,0:NJ),CurlV_res(0:NI,0:NJ))
  gradP = 0
  divV = 0
  divVP = 0
  lapP = 0
  CurlV = 0

!===  READ GRID ===
  WRITE(*,*) 'Read mesh from file: ', MeshFile
  READ(IO,*) ((X(I,J),Y(I,J),rtmp,I=1,NI),J=1,NJ)
  CLOSE(IO)

!=== CALCULATE METRIC ===
  WRITE(*,*) 'Calculate metric'       
  Call B_CalcMetric(NI,NJ,X,Y,CellCenter,CellVolume,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector) 

!=== INITIATE FIELDS ===
  if (sols == 1) then
	P = 0; V = 0 
	open(io,file=sol_file)
	read(io,*)
	read(io,*)
	read(io,*) ((rtmp,rtmp,V(i,j,1),V(i,j,2),rtmp,P(i,j),rtmp,rtmp, i=0,NI), J=0,NJ)
	close(io)
  end if

!=== CALCULATE GRADIENT ===
  WRITE(*,*) 'Calculate derivatives'
  do i=1,20
  Call B_CalcGradient(NI,NJ,X,Y,P,GradP,CellVolume,CellCenter,    	  &
&											IFaceVector,JFaceVector,  &
&											IFaceCenter,JFaceCenter)
  enddo
!===CALCULATE DIVERGENCE ===
  call B_CalcDiv(NI,NJ,X,Y,V,DivV,CellVolume,CellCenter,    &
&											IFaceVector,JFaceVector,  &
&											IFaceCenter,JFaceCenter)

  call B_CalcDivphi(NI,NJ,X,Y,V,P,gradP,DivVP,CellVolume,CellCenter,    &
&											IFaceVector,JFaceVector,  &
&											IFaceCenter,JFaceCenter,sch)

!===CALCULATE LAPLACIAN ===
  call B_CalcLap(NI,NJ,X,Y,p,gradP,lapP,CellVolume,CellCenter,    &
&											IFaceVector,JFaceVector,  &
&											IFaceCenter,JFaceCenter)

!===CALCULATE CURL ===
  call B_CalcCurl(NI,NJ,X,Y,V,CurlV,CellVolume,CellCenter,    &
&											IFaceVector,JFaceVector,  &
&											IFaceCenter,JFaceCenter)

!===SOLVE PARTICLE DYNAMICS
!Calculate Stokes number
print*, 'Sk=', rho_p * d_p**2 /(18*mu) !характерная скорость и размер 1

!1st order Euler
r = r_0; v_p=v_0
F = 0
m_p = pi* d_p**3/6 * rho_p
print*, 'm_p=', m_p

frmt_screen = '(a,i6,1x,5(a,ES15.8,1x),2(a,i6,x))'
frmt_file = '(i6,x,5(ES23.16,x))'
open(newunit=iu, file=particle_out)
do i=1, nt
	call C_Location(r(1), r(2), NI, NJ, X, Y, CellVolume, Ip, Jp)
	U(1) = V(IP, JP, 1)
	U(2) = V(IP, JP, 2)
	
	!Расчёт сил
	call C_CalcForce(rho_f,rho_p,mu,d_p,v_p,u,F)
	F = 0
	
	F(1:2) = F(1:2)
	F(3) = F(3)/J_p
	t = t + dt
	r = r_0 + dt*v_p
	v_p = v_0 + dt*F
	!Определение пересечения с границей
	call C_Location(r(1), r(2), NI, NJ, X, Y, CellVolume, Ip, Jp)
	
	if (Ip==-1 .or. Jp==-1) then
		call C_Boundary(x,y,IFaceVector,JFaceVector,NI,NJ, &
		r_0(1),r_0(2),v_0(1),v_0(2),v_0(3),r(1),r(2),v_p(1),v_p(2),v_p(3),Ip,Jp,St)
		write(*,*) 'Status=', St
		
		! St:
		!      2   
		!   3     4
		!      1
		
	end if
	
	r_0 = r
	v_0 = v_p
	
	write(*,frmt_screen) 'i=', i, 't=', t, 'x=', r(1), 'y=', r(2), &
					  'F_x=', F(1), 'F_y=', F(2), 'Ip=', Ip, 'Jp=', Jp
	write(iu,frmt_file) i, t, r(1), r(2), F(1), F(2)		
end do
close(iu)


!=== OUTPUT FIELDS ===
  WRITE(*,*) 'Output fields to file: ', OutputFile       
  Open(IO,FILE=OutputFile)
  Call B_OutputFields(IO,NI,NJ,X,Y,P,V,GradP,divV,lapP,divVP,curlV)
  Close(IO)


	contains
	

END PROGRAM Main  
