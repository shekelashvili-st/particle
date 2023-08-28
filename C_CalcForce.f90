Subroutine C_CalcForce(rho_f,rho_p,mu,d_p,v_p,u,F)
implicit none
	real:: rho_f, rho_p, mu, d_p, v_p(3), u(2), F(3)
	real:: pi,m_p,V_r(2),Re_p, f_coef
	
	pi=4*atan(1.0_8)
	
	
	!Относительная скорость, число Рейнольдса
	V_r = u-v_p(1:2)
	Re_p = rho_f * D_p * norm2(V_r)/mu
	
	
	
	!Сила Стокса
	f_coef = (1 + 0.179*Re_p**0.5 + 0.013*Re_p)
	F(1:2) = f_coef*mu*18*V_r / (rho_p*D_p**2)
	
	
End Subroutine 