

Program P3
	Implicit none
	real(8), parameter :: mass = 1, rho=0.7
	real(8), parameter :: k_b = 1.380649e-23
	integer, parameter :: N=125, Nsteps_ini = 2, Nsteps_prod = 500000
	real(8), dimension(N, 3) :: r, r_ini, vel, vel_ini, r_out
	integer :: step, i, dt_index
	real(8) :: pot, K_energy, L, cutoff, sigma, M, a, Temp, inst_temp, dt, absV, p
	real(8), dimension(7) :: dt_list
	integer, allocatable :: seed(:)
	integer :: nn
	external inst_temp

	dt_list = (/ 10e-1, 10e-2, 10e-3, 10e-4, 10e-5, 10e-6, 10e-7 /)

	call random_seed(size=nn)
	allocate(seed(nn))
	seed = 123456789    ! putting arbitrary seed to all elements
	call random_seed(put=seed)
	deallocate(seed)

	L = (N/rho)**(1./3.)
	Temp = 100
	sigma = Temp**(1.d0/2.d0)

	M = N**(1./3.)
	a = L/M

	print*, L, M, a

	cutoff = L

	! """"
	! ii) Initialize system and run simulation using velocity Verlet
	! 
	! """"
	

	! Initialize bimodal distrubution: v_i = +- sqrt(k_b T / m)
	absV = (k_b * Temp / mass)**(1./2.)

	open(22, file="vel_ini.dat")
	
	call initialize_positions(N, rho, r_ini)
	call initialize_velocities(N, absV, vel_ini)

	do i=1,N
		write(22,*) vel_ini(i,:)
	end do

	close(22)

	open(44, file="energy_verlet.dat")

	! Apply Verlet algorithm
	do dt_index = 1, 7
		dt = dt_list(dt_index)
		write(44,*) ""
		write(44,*) ""
		write(44,*) "dt = ", dt

		r = r_ini
		vel = vel_ini

		do step = 1,Nsteps_ini
			call time_step_vVerlet(r, vel, pot, N, L, cutoff, dt)
		!	print*, r, vel
			call kinetic_energy(vel, K_energy, N)
			! Momentum p = m*v
			p = (2*mass*K_energy)**(1./2.)
			!print*, step, pot, K_energy, pot+K_energy, p
			write(44,*) step, pot, K_energy, pot+K_energy, p
			if (mod(step, 1000).eq.0) then
				print*, real(step)/Nsteps_ini
			end if
		end do

	Temp = inst_temp(N, K_energy)
	write(*,*) "dt = ", dt
	write(*,*) "T = ", Temp

	end do

	close(44)

	open(23, file="vel_fin_Verlet.dat")
	do i = 1, N
		write(23,*) vel(i,:)
	end do
	close(23)

	! """"
	! iii) Initialize system and run simulation using Euler
	! 
	! """"

	! Initialize again, now to apply Euler method
	! Apply Euler algorithm

	open(45, file="energy_euler.dat")

	do dt_index = 1, 7
		dt = dt_list(dt_index)
		write(45,*) ""
		write(45,*) ""
		write(45,*) "dt = ", dt

		r = r_ini
		vel = vel_ini

		! Initialize velocities and positions

		do step = 1,Nsteps_ini
			call time_step_Euler_pbc(r, r_out, vel, N, L, cutoff, dt, pot)
			call kinetic_energy(vel, K_energy, N)
			! Momentum p = m*v
			p = (2*mass*K_energy)**(1./2.)
			write(45,*) step, pot, K_energy, pot+K_energy, p
			if (mod(step, 1000).eq.0) then
				print*, real(step)/Nsteps_ini
			end if
			r = r_out
		end do

	Temp = inst_temp(N, K_energy)
	write(*,*) "dt = ", dt
	write(*,*) "T = ", Temp

	end do

	close(45)
	
	open(24, file="vel_fin_Euler.dat")
	do i = 1, N
		write(24,*) vel(i,:)
	end do
	close(24)
	

End Program


Subroutine therm_Andersen(vel, nu, sigma, N)
	Implicit none
	integer :: i, N
	real(8) :: rand, nu, sigma
	real(8), dimension(N, 3) :: vel
	real(8), dimension(2) :: xnums

 	do i = 1, N
 		call random_number(rand)
 		if (rand.lt.nu) then
 			call BM(2, xnums, sigma)
 			!print*, "xnums: ", xnums
 			vel(i, 1) = xnums(1)
 			vel(i, 2) = xnums(2)
 			call BM(2, xnums, sigma)
 			vel(i, 3) = xnums(1)
		end if
	end do
!	print*, vel
End Subroutine


Subroutine BM(ndat,xnums,sigma)
    Implicit none
    Integer ::  ndat, i
    real(8), dimension(ndat) :: xnums
    real(8) :: r, phi, x1, x2, sigma
    real(8), parameter :: pi = 4.d0*atan(1.d0)
!     ATENCIÓ! Es generen 2ndat numeros
    Do i = 1, ndat, 2
        r = sqrt(-2.d0*log(1.d0-rand()))
        phi = 2.d0*pi*rand()
        x1 = r*cos(phi)
        x2 = r*sin(phi)
        if (i.ne.ndat) then ! Ens assegurem que no haguem acabat la llista
    	    xnums(i) = x1*sigma
            xnums(i+1) = x2*sigma
        endif
    end do
    return
end Subroutine

Subroutine initialize_positions(N, rho, r)
	Implicit none
	integer, intent(in) :: N
	real(8), dimension(N, 3), intent(out) :: r
	real(8) :: L, a, x, y, z
	integer :: M, i, j, k, particle
	real(8), intent(in) :: rho

	L = (N/rho)**(1./3.)

	M = N**(1./3.)
	
	a = L/M

	! Set the position of every particle
	particle = 1
	do i = 1, M
		do j = 1, M
			do k = 1, M
				x = i*a
				y = j*a
				z = k*a
				r(particle, :) = (/x, y, z/)
				particle = particle + 1
			end do
		end do
	end do
End Subroutine

Subroutine time_step_vVerlet(r, vel, pot, N, L, cutoff, dt)
	Implicit none
	integer, intent(in) :: N
	real(8), dimension(N, 3) :: r, F, vel
	real(8) :: dt, L, pot, cutoff
	integer :: i, k

	Call find_force_LJ(r, N, L, cutoff, F, pot)

	do i = 1, N
		do k = 1, 3
 			r(i, k) = r(i, k) + vel(i, k) * dt + 0.5*F(i, k)*dt*dt
			call pbc1(r(i,k), L)
			vel(i, k) = vel(i, k) + F(i, k)* 0.5 * dt
		end do
 	end do

	Call find_force_LJ(r, N, L, cutoff, F, pot)

	print*, F
	do i = 1, N
		do k = 1, 3
			vel(i, k) = vel(i, k) + F(i, k)* 0.5 * dt
		end do
 	end do

End Subroutine


Subroutine pbc1(x, L)
	Implicit none
	real(8) :: x, L

	if (x.ge.L/2.) then
		x = x - L
 	else if (x.le.(-L/2.)) then
 		x = x + L
 	end if

	Return
End Subroutine

Subroutine find_force_LJ(r, N, L, cutoff, F, pot)
	Implicit none
	real(8), dimension(N, 3), intent(in) :: r
	real(8), intent(in) :: L, cutoff
	real(8) :: dx, dy, dz, d
	integer :: i, j, k
	integer, intent(in) :: N
	real(8), dimension(N, 3), intent(out) :: F
	real(8), intent(out) :: pot
	external pbc1

	pot = 0.d0
!	print*, L

	do i = 1, N
		do k = 1, 3
			F(i,k) = 0.d0
		end do
	end do

	do i = 1, N
		do j = i+1, N
			dx = r(i, 1) - r(j, 1)
			dy = r(i, 2) - r(j, 2)
			dz = r(i, 3) - r(j, 3)
			call pbc1(dx, L)
			call pbc1(dy, L)
			call pbc1(dz, L)
	!		print*, i, j, "dx:", dx, dy, dz
			d = (dx**2+dy**2+dz**2)**(1.d0/2.d0)
			if (d.le.cutoff) then
				F(i,1) = F(i,1) + (48.d0 / d**14 - 24.d0 / d**8)*dx
				F(i,2) = F(i,2) + (48.d0 / d**14 - 24.d0 / d**8)*dy
				F(i,3) = F(i,3) + (48.d0 / d**14 - 24.d0 / d**8)*dz
				F(j,1) = F(j,1) - (48.d0 / d**14 - 24.d0 / d**8)*dx
				F(j,2) = F(j,2) - (48.d0 / d**14 - 24.d0 / d**8)*dy
				F(j,3) = F(j,3) - (48.d0 / d**14 - 24.d0 / d**8)*dz
				pot = pot + 4.d0*( 1.d0/ d**12 - 1.d0 /d**6) - 4.d0*( 1/ cutoff**12 - 1.d0 /cutoff**6)
	!			print*, j, F(i,1), d, dx, (48 / d**14 - 24 / d**8)*dx
			end if 
!			print*, i, j, dx, dy, dz, d, F(j, :)

		end do
	!	print*, F(i,1)
	end do

	!do i = 1, N
	!	print*, F(i, :)
	!end do


End Subroutine

Subroutine kinetic_energy(vel, K_energy, N)
	Implicit none
	integer, intent(in) :: N
	real(8), dimension(N, 3) :: vel
	integer :: i, k
	real(8) :: K_energy

	K_energy = 0
	! for each particle
	do i = 1, N
		! loop over coordinates
		do k = 1, 3
			K_energy = K_energy + 0.5*vel(i,k)*vel(i,k) 
		end do
	end do
End Subroutine

Function inst_temp(N, K_energy)
	Implicit none
	integer :: N, N_f
	real(8), parameter :: k_b = 1.380649e-23
	real(8) :: K_energy, inst_temp

	N_f = 3*N - 3
	inst_temp = 2.d0/(N_f * k_b)*K_energy
	Return
End Function


Subroutine time_step_Euler_pbc(r_in, r_out, vel, N, L, cutoff, dt, pot)
 	Implicit none
 	integer, intent(in) :: N
 	real(8), dimension(N, 3), intent(in) :: r_in
 	real(8), dimension(N, 3), intent(out) :: r_out
 	real(8), dimension(N, 3) :: vel, F
 	real(8), intent(in) :: L, dt
 	integer :: i
 	real(8) :: cutoff, pot

 	Call find_force_LJ(r_in, N, L, cutoff, F, pot)


 	do i = 1, N
		r_out(i, 1) = r_in(i, 1) + vel(i, 1) * dt + 0.5*F(i, 1)*dt*dt
		r_out(i, 2) = r_in(i, 2) + vel(i, 2) * dt + 0.5*F(i, 2)*dt*dt
		r_out(i, 3) = r_in(i, 3) + vel(i, 3) * dt + 0.5*F(i, 3)*dt*dt
		vel(i, 1) = vel(i, 1) + F(i, 1) * dt
		vel(i, 2) = vel(i, 2) + F(i, 2) * dt
		vel(i, 3) = vel(i, 3) + F(i, 3) * dt
		call pbc1(r_out(i,1), L)
		call pbc1(r_out(i,2), L)
		call pbc1(r_out(i,3), L)
	end do

End Subroutine

Subroutine initialize_velocities(N, absV, vel)
	Implicit none
	integer, intent(in) :: N
	real(8) :: absV, rnd
	integer :: i, j
	real(8), dimension(N, 3) :: vel

	do i=1,N
		do j = 1,3
			call random_number(rnd)
			if (rnd.ge.0.5) then
				vel(i, j) = +absV
			else if (rnd.lt.0.5) then
				vel(i, j) = -absV
			end if
		end do
	end do


End Subroutine