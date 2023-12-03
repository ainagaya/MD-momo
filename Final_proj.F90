

Program P_final_a
	Implicit none
	real(8), parameter :: mass = 1, rho=0.7, epsilon=1 , sigma=1 
!	real(8), parameter :: k_b = 1.380649e-23
	integer, parameter :: N=125, Nsteps_prod = 500000
	real(8), dimension(N, 3) :: r, r_ini, vel, vel_ini, r_out
	integer :: step, i, dt_index, Nsteps
	real(8) :: pot, K_energy, L, cutoff, M, a, Temp, inst_temp, dt, absV, p, tini, tfin
	real(8), dimension(3) :: dt_list
	integer, allocatable :: seed(:)
	integer :: nn
	external inst_temp

	dt_list = (/ 1e-3, 1e-4, 1e-5/)

	call random_seed(size=nn)
	allocate(seed(nn))
	seed = 123456789    ! putting arbitrary seed to all elements
	call random_seed(put=seed)
	deallocate(seed)

	L = (N/rho)**(1./3.)
	Temp = 100

	M = N**(1./3.)
	a = L/(M)

	print*, L, M, a

	cutoff = 2.5

	! """"
	! ii) Initialize system and run simulation using velocity Verlet
	! 
	! """"

	! Initialize bimodal distrubution: v_i = +- sqrt(T' / m)
	absV = (Temp / mass)**(1./2.)
	print*, absV

	open(22, file="vel_ini.dat")
	
	call initialize_positions(N, rho, r_ini)

	call initialize_velocities(N, absV, vel_ini)

	do i=1,N
		write(22,*) vel_ini(i,:)
	end do

	close(22)

	open(44, file="energy_verlet.dat")
	open(77, file="Temperatures.dat")

	tini = 0
	tfin = 1

	! Apply Verlet algorithm
	do dt_index = 1, 3
		dt = dt_list(dt_index)
		write(44,*) ""
		write(44,*) ""
		write(44,*) "dt = ", dt

		Nsteps = int((tfin - tini)/dt)


		! We roll back to the initial positions and velocities to initialize
		r = r_ini
		vel = vel_ini

		do step = 1,Nsteps

			call time_step_vVerlet(r, vel, pot, N, L, cutoff, dt)
		!	print*, r, vel
			call kinetic_energy(vel, K_energy, N)
			! Momentum p = m*v
			call momentum(vel, p, N)
			!p = (2*mass*K_energy)**(1./2.)
			!print*, step, pot, K_energy, pot+K_energy, p
			write(44,*) step*dt, pot, K_energy, pot+K_energy, p
		!	print*, K_energy
			if (mod(step, 1000).eq.0) then
				print*, real(step)/Nsteps
			end if
		end do


	
	Temp = inst_temp(N, K_energy)
	write(77,*) "Verlet", dt, Temp


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

	dt_list = (/ 1e-3, 1e-4, 1e-5/)


	! Initialize again, now to apply Euler method
	! Apply Euler algorithm

	open(45, file="energy_euler.dat")

	do dt_index = 1, 3
		dt = dt_list(dt_index)
		write(45,*) ""
		write(45,*) ""
		write(45,*) "dt = ", dt

		r = r_ini
		vel = vel_ini

		Nsteps = int((tfin - tini)/dt)

		! Initialize velocities and positions

		do step = 1,Nsteps
			call time_step_Euler_pbc(r, r_out, vel, N, L, cutoff, dt, pot)
			call kinetic_energy(vel, K_energy, N)
			! Momentum p = m*v
			!p = (2*mass*K_energy)**(1./2.)
			call momentum(vel, p, N)

			write(45,*) step*dt, pot, K_energy, pot+K_energy, p
			if (mod(step, 1000).eq.0) then
				print*, real(step)/Nsteps
			end if
			r = r_out
		end do

	Temp = inst_temp(N, K_energy)
	write(77,*) "Euler", dt, Temp

	end do

	close(45)
	close(77)
	
	open(24, file="vel_fin_Euler.dat")
	do i = 1, N
		write(24,*) vel(i,:)
	end do
	close(24)
	

End Program

!########################################################################################################

Subroutine initialize_positions(N, rho, r)
! """"
! Calculates the positions r of N particles in a SC structure
! INPUTS: N, rho
! OUTPUT: r(N, 3)
! """"
	Implicit none
	integer, intent(in) :: N
	real(8), intent(in) :: rho
	real(8), dimension(N, 3), intent(out) :: r
	real(8) :: L, a, x, y, z, ini
	integer :: M, i, j, k, particle

	L = (N/rho)**(1./3.)

	M = N**(1./3.)
	
	a = L/(M)

	! Set the position of every particle
	particle = 1
	ini = -L/2.d0
	! ini = 0.d0
	do i = 0, M-1
		do j = 0, M-1
			do k = 0, M-1
				x = ini + i*a
				y = ini + j*a
				z = ini + k*a
				r(particle, :) = (/x, y, z/)
				particle = particle + 1
			end do
		end do
	end do
End Subroutine

!########################################################################################################

Subroutine time_step_vVerlet(r, vel, pot, N, L, cutoff, dt)
	Implicit none
	integer, intent(in) :: N
	real(8), dimension(N, 3) :: r, F, vel
	real(8) :: dt, L, pot, cutoff
	integer :: i
	external pbc, find_force_LJ

	Call find_force_LJ(r, N, L, cutoff, F, pot)


	do i = 1, N
 		r(i, :) = r(i, :) + vel(i, :) * dt + 0.5*F(i, :)*dt*dt
 			
		do while (any(r(i,:) > L/2.) .or. any(r(i,:) < -L/2.))
    		! Apply periodic boundary conditions using the pbc subroutine
    		call pbc(r(i,:), L, size(r(i,:)))
		end do

			
		vel(i, :) = vel(i, :) + F(i, :)* 0.5 * dt
 	end do

	Call find_force_LJ(r, N, L, cutoff, F, pot)

	do i = 1, N
		vel(i, :) = vel(i, :) + F(i, :)* 0.5 * dt
	end do

End Subroutine

!########################################################################################################

Subroutine pbc1(x, L)
	Implicit none
	real(8) :: x, L

	if (x.gt.L/2.) then
		x = x - L
		if (abs(x).gt.1000) then
			stop 
		end if
 	else if (x.lt.(-L/2.)) then
 		x = x + L
		if (abs(x).gt.1000) then
			stop
		end if
 	end if

	Return
End Subroutine

!########################################################################################################

Subroutine pbc(vector, L, D)
	Implicit none
	integer :: D, i
	real(8), dimension(D), intent(inout) :: vector
	real(8), intent(in) :: L


	do i = 1, D
		if (vector(i).gt.L/2.) then
			vector(i) = vector(i) - L
			if (abs(vector(i)).gt.1000) then
				print*, "YOU HAVE INESTABILITIES!"
				stop 
			end if
 		else if (vector(i).lt.(-L/2.)) then
 			vector(i) = vector(i) + L
			if (abs(vector(i)).gt.1000) then
				print*, "YOU HAVE INESTABILITIES!"
				stop
			end if
 		end if
 	end do

	Return
End Subroutine

!########################################################################################################

Subroutine find_force_LJ(r, N, L, cutoff, F, pot)
	Implicit none
	real(8), dimension(N, 3), intent(in) :: r
	real(8), intent(in) :: L, cutoff
	real(8) :: d, f_ij
	real(8), dimension(3) :: d_r
	integer :: i, j, k
	integer, intent(in) :: N
	real(8), dimension(N, 3), intent(out) :: F
	real(8), intent(out) :: pot
	external pbc

	pot = 0.d0

	F = 0.d0

	do i = 1, N
		do j = i+1, N
			d_r(:) = r(i, :) - r(j, :)
!			print*, r(i, :) - r(j, :)

			do while (any(d_r(:).gt.L/2.).or.(any(d_r(:).lt.(-L/2.))))
				call pbc(d_r, L, size(d_r))
			end do


			d = (d_r(1)**2+d_r(2)**2+d_r(3)**2)**(1.d0/2.d0)
			if (d.le.cutoff) then
				f_ij = 48.d0 / d**14 - 24.d0 / d**8
!				print*, d_r
				F(i,:) = F(i,:) + f_ij*d_r(:)
				F(j,:) = F(j,:) - f_ij*d_r(:)


				if (isnan(F(i,1))) then
					print*, i, j
					stop
				end if

				pot = pot + 4.d0*( 1.d0/ d**12 - 1.d0 /d**6) - 4.d0*( 1/ cutoff**12 - 1.d0 /cutoff**6)

			end if 

		end do
	end do


End Subroutine

!########################################################################################################

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
!	real(8), parameter :: k_b = 1.380649e-23
	real(8) :: K_energy, inst_temp

	N_f = 3*N - 3
	! inst_temp = 2.d0/(N_f * k_b)*K_energy
	inst_temp = 2.d0/(N_f)*K_energy
	Return
End Function

!########################################################################################################

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
 		r_out(i, :) = r_in(i, :) + vel(i, :) * dt + 0.5*F(i, :)*dt*dt
 		vel(i, :) = vel(i, :) + F(i, :) * dt

		do while (any(r_out(i,:).gt.L/2.).or.(any(r_out(i,:).lt.(-L/2.))))
			call pbc(r_out, L, size(r_out))
		end do
		

	end do

End Subroutine

!###########################################################################################################

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


!#############################################################

Subroutine momentum(vel, p, N)
	Implicit none
	real(8), dimension(N, 3) :: vel
	real(8), dimension(3) :: total_p
	integer :: N, i
	real(8), intent(out) :: p 

	total_p(:) = 0

	! Accumulate p
	do i = 1,N
		total_p(:) = total_p(:) + vel(i, :)
	end do

	! Produce the module
	p = ( total_p(1)**2 + total_p(2)**2 + total_p(3)**2 )**(1./2.) 


End Subroutine