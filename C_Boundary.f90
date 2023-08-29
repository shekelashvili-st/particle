Subroutine C_Boundary(x,y,IFaceVector,JFaceVector,NI,NJ, &
		x_m,y_m,u_m,v_m,w_m,x_m1,y_m1,u_m1,v_m1,w_m1,Ip1,Jp1,St)
implicit none
 integer::NI,NJ,I,J
 real,dimension(NI,NJ):: X,Y
 real:: IFaceVector(NI,NJ-1,2), JFaceVector(NI-1,NJ,2)
 real :: x_m1,y_m1,u_m1,v_m1,w_m1,x_m,y_m,u_m,v_m,w_m
 real:: x_new, y_new
 integer:: Ip1, Jp1, St,icr
 real:: xcros,ycros
 real:: TV(2), NV(2), D, Vn, Vtau
 
 St = 0
 
 	! St:
	!      2   
	!   3     4
	!      1
	
 do i=1,ni-1
	
	call Cross_edges(x_m,y_m,x_m1,y_m1,x(i,1),y(i,1),x(i+1,1),y(i+1,1),xcros,ycros,icr)
	if (icr==1) then
		st=1
		call reflect(x_m, y_m, x_m1, y_m1, x(i, 1), y(i, 1), x(i+1, 1), y(i+1, 1), x_new, y_new)
		write(*,*) 'Particle crossed boundary', St, 'at', xcros,ycros
		x_m1 = x_new
		y_m1 = y_new		
		!Случай зеркального отражения
		D = sqrt(JFaceVector(I, 1, 1)**2 + JFaceVector(I, 1, 2)**2)
			
		NV(1) = -JFaceVector(i, 1, 1) / D
		NV(2) = -JFaceVector(i, 1, 2) / D
		TV(1) = -NV(2)
		TV(2) = NV(1)
			
		Vn = u_m * NV(1) + v_m * NV(2)
		Vtau = u_m * TV(1) + v_m * TV(2)
			
		u_m1 = Vtau * TV(1) - Vn * NV(1)
		v_m1 = Vtau * TV(2) - Vn * NV(2)
		w_m1 = w_m		
		
		return
	end if
	
	call Cross_edges(x_m,y_m,x_m1,y_m1,x(i,nj),y(i,nj),x(i+1,nj),y(i+1,nj),xcros,ycros,icr)
	if (icr==1) then
		st=2
!		call reflect(x_m, y_m, x_m1, y_m1, x(i, nj), y(i, nj), x(i+1, nj), y(i+1, nj), x_new, y_new)
		write(*,*) 'Particle crossed boundary', St, 'at', xcros,ycros
!		x_m1 = x_new
!		y_m1 = y_new		
!		!Случай зеркального отражения
!		D = sqrt(JFaceVector(I, nj, 1)**2 + JFaceVector(I, nj, 2)**2)
!			
!		NV(1) = -JFaceVector(i, nj, 1) / D
!		NV(2) = -JFaceVector(i, nj, 2) / D
!		TV(1) = -NV(2)
!		TV(2) = NV(1)
!			
!		Vn = u_m * NV(1) + v_m * NV(2)
!		Vtau = u_m * TV(1) + v_m * TV(2)
!			
!		u_m1 = Vtau * TV(1) - Vn * NV(1)
!		v_m1 = Vtau * TV(2) - Vn * NV(2)
!		w_m1 = w_m		
		
		return
	end if

	call Cross_edges(x_m,y_m,x_m1,y_m1,x(1,j),y(1,j),x(1,j+1),y(1,j+1),xcros,ycros,icr)
	if (icr==1) then
		st=3
		call reflect(x_m, y_m, x_m1, y_m1, x(1, j), y(1, j), x(1, j+1), y(1, j+1), x_new, y_new)
		write(*,*) 'Particle crossed boundary', St, 'at', xcros,ycros
		x_m1 = x_new
		y_m1 = y_new		
		!Случай зеркального отражения, ПЛОСКОСТЬ СИММЕТРИИ
		D = sqrt(IFaceVector(1, j, 1)**2 + IFaceVector(1, j, 2)**2)
			
		NV(1) = -IFaceVector(1, j, 1) / D
		NV(2) = -IFaceVector(1, j, 2) / D
		TV(1) = -NV(2)
		TV(2) = NV(1)
			
		Vn = u_m * NV(1) + v_m * NV(2)
		Vtau = u_m * TV(1) + v_m * TV(2)
			
		u_m1 = Vtau * TV(1) - Vn * NV(1)
		v_m1 = Vtau * TV(2) - Vn * NV(2)
		!ПЛОСКОСТЬ СИММЕТРИИ
		w_m1 = -w_m		
		
		return
	end if

	call Cross_edges(x_m,y_m,x_m1,y_m1,x(ni,j),y(ni,j),x(ni,j+1),y(ni,j+1),xcros,ycros,icr)
	if (icr==1) then
		st=4
		call reflect(x_m, y_m, x_m1, y_m1, x(ni, j), y(ni, j), x(ni, j+1), y(ni, j+1), x_new, y_new)
		write(*,*) 'Particle crossed boundary', St, 'at', xcros,ycros
		x_m1 = x_new
		y_m1 = y_new		
		!Случай зеркального отражения, ПЛОСКОСТЬ СИММЕТРИИ
		D = sqrt(IFaceVector(ni, j, 1)**2 + IFaceVector(ni, j, 2)**2)
			
		NV(1) = -IFaceVector(ni, j, 1) / D
		NV(2) = -IFaceVector(ni, j, 2) / D
		TV(1) = -NV(2)
		TV(2) = NV(1)
			
		Vn = u_m * NV(1) + v_m * NV(2)
		Vtau = u_m * TV(1) + v_m * TV(2)
			
		u_m1 = Vtau * TV(1) - Vn * NV(1)
		v_m1 = Vtau * TV(2) - Vn * NV(2)
		!ПЛОСКОСТЬ СИММЕТРИИ
		w_m1 = -w_m		
		
		return
	end if	
	
 end do

End Subroutine 
	
 subroutine Cross_edges(x_m,y_m,x_m1,y_m1,x_3,y_3,x_4,y_4,xcros,ycros,icr)
 real:: x_m,y_m,x_m1,y_m1
 real:: x_3,y_3,x_4,y_4
 real:: xcros, ycros, eps
 integer:: icr
 real:: r(2),s(2),denom, num_u,num_t
 real:: q(2),p(2)
 eps = 1e-10
 
 q(1) = x_m
 q(2) = y_m
 p(1) = x_3
 p(2) = y_3
 
 r(1) = x_m1-x_m
 r(2) = y_m1-y_m
 s(1) = x_4-x_3
 s(2) = y_4-y_3
 
 
 icr = 0
 denom = r(1)*s(2)-r(2)*s(1)
 num_t = (q(1)-p(1))*s(2) - (q(2)-p(2))*s(1)
 num_u = -((p(1)-q(1))*r(2) - (p(2)-q(2))*r(1))
 
 if (abs(denom)<eps) then
	if(abs(num_u)>eps) return !parallel lines
	Icr = 1
	Xcros = x_m1
	Ycros = y_m1
	return
 end if
 
 t = num_t/denom
 u = -num_u/denom
 
 if ((t>0 .and. t<1) .and. (u>0 .and. u<1)) then
	icr = 1
	xcros = p(1)+t*r(1)
	ycros = p(2)+t*r(2)
 end if
 
 end subroutine
 
subroutine reflect(x1, y1, x2, y2, x3, y3, x4, y4, x2_new, y2_new)
real:: x1, y1, x2, y2, x3, y3, x4, y4
real:: x2_new, y2_new
real:: a, b, c, denom

	a = y4 - y3
	b = x3 - x4
	c = -(a*x4+b*y4)
	denom = a*a + b*b
		
	x2_new = ((b*b-a*a)*x2 - 2.0*a*b*y2 - 2*c*a)/denom
	y2_new = (-2.0*a*b*x2 + (a*a-b*b)*y2  - 2*c*b)/denom
	
end subroutine

