Subroutine C_Boundary(x,y,IFaceVector,JFaceVector,NI,NJ, &
		x_m,y_m,u_m,v_m,w_m,x_m1,y_m1,u_m1,v_m1,w_m1,Ip1,Jp1,St)
implicit none
 integer::NI,NJ,I,J
 real,dimension(NI,NJ):: X,Y
 real:: IFaceVector(NI,NJ-1,2), JFaceVector(NI-1,NJ,2)
 real :: x_m1,y_m1,u_m1,v_m1,w_m1,x_m,y_m,u_m,v_m,w_m
 integer:: Ip1, Jp1, St,icr
 real:: xcros,ycros
 
 St = 0
 
 do i=1,ni-1
	
	call Cross_edges(x_m,y_m,x_m1,y_m1,x(i,1),y(i,1),x(i+1,1),y(i+1,1),xcros,ycros,icr)
	if (icr==1) then
		st=1
		write(*,*) 'Particle crossed boundary', St, 'at', xcros,ycros
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
 eps = 1e-4
 
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

