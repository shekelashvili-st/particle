Subroutine C_Location(x_p, y_p, NI, NJ, X, Y, CellVolume, Ip, Jp)
implicit none
 integer::NI,NJ,I,J
 real,dimension(NI,NJ):: X,Y
 real,dimension(NI-1,NJ-1):: CellVolume ! scalar arrays
 real :: x_p, y_p, S, eps
 integer:: Ip, Jp
 
 Ip = -1; Jp = -1
 eps = 1e-10
 
 outer: do j=1, nj-1
	do i=1, ni-1
		S = tri_S(x(i,j),y(i,j),x(i,j+1),y(i,j+1),x_p,y_p) + &
			tri_S(x(i+1,j),y(i+1,j),x(i+1,j+1),y(i+1,j+1),x_p,y_p) + &
			tri_S(x(i,j),y(i,j),x(i+1,j),y(i+1,j),x_p,y_p) + &
			tri_S(x(i,j+1),y(i,j+1),x(i+1,j+1),y(i+1,j+1),x_p,y_p)
		
		if (abs(S-CellVolume(i,j))<=eps) then
			Ip = I
			Jp = J
			exit outer
		end if
	enddo
 enddo outer
	
 contains
 
 function tri_S(x1,y1,x2,y2,x3,y3)
 real:: tri_S, x1, y1, x2, y2, x3, y3
 real:: a, b, c, p
 a = sqrt((x2-x1)**2+(y2-y1)**2)
 b = sqrt((x3-x1)**2+(y3-y1)**2)
 c = sqrt((x3-x2)**2+(y3-y2)**2)
 p = (a+b+c)/2
 
 tri_S = sqrt(p*(p-a)*(p-b)*(p-c))
 end function
 
End Subroutine 