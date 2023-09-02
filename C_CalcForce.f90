Subroutine C_CalcForce(rho_f,rho_p,mu,d_p,v_p,u,F,&
						v_p0,dt,&
						w_f,&
						fres,fadd,fsaff)
implicit none
	real:: rho_f, rho_p, mu, d_p, v_p(3), u(2), F(3)
	real::fres(2),fadd(2),fsaff(2)
	real:: pi,m_p,V_r(2),Re_p, f_coef
	!added mass
	real::v_p0(3), dt, f_add(2)
	!Saffman
	real::K=6.46,w_f,f_saff(2)
	
	pi=4*atan(1.0_8)
	
	
	!Относительная скорость, число Рейнольдса
	V_r = u-v_p(1:2)
	Re_p = rho_f * D_p * norm2(V_r)/mu
	
	
	!Сила Стокса
	f_coef = (1 + 0.179*Re_p**0.5 + 0.013*Re_p)
	F(1:2) = f_coef*mu*18*V_r / (rho_p*D_p**2)
	fres = F(1:2)
	
	m_p = rho_p * 4./3 * pi * (D_p/2)**3
	!Сила присоединенных масс
	f_add = -2.0_8/3*pi*(D_p/2)**3*rho_f * (v_p(1:2)-v_p0(1:2))/dt / m_p
	fadd = f_add
	
	!Сила Сэффмана
	f_saff(1) = sqrt(2.0_8)*K/4*d_p**2 * sqrt(mu*rho_f) * (u(2)-v_p(2))* w_f/sqrt(abs(w_f)) / m_p
	f_saff(2) = -sqrt(2.0_8)*K/4*d_p**2 * sqrt(mu*rho_f) * (u(1)-v_p(1)) * w_f/sqrt(abs(w_f)) / m_p
	fsaff=f_saff
	
	!Суммарная сила
	F(1:2) = F(1:2) + f_add + f_saff
	
End Subroutine 