! Analysis of properties of a Lennard-Jones liquid. Considering that the interaction between argon atoms
! can be described using a Lennard-Jones potential with epsilon = 0.998kJ/mol and σ = 3.4Å (and the atomic mass of
! argon is m = 40g/mol):
! (a) Calculate the kinetic, potential and total energy per particle of a system of (N > 100) argon atoms in a
! fluid state with densities ρ = 0.05, 0.1, 0.2, 0.4, 0.6, 0.8m/σ 3 and a fixed temperature (using a thermal bath)
! kB T = 1.2.
! • Plot a graph with the kinetic, potential and total energy as a function of density. Show your results in
! units of kJ/mol for the energies and g/cm3 for the density. How do you equilibrate your system? How
! many timesteps do you use to calculate your results?
! • Visualize (part of) the equilibrium trajectory and interpret your observations in terms of the phase
! diagram of the Lennard-Jones fluid shown in Fig. 1.
! (b) Calculate the pressure of a system of (N > 100) argon atoms at kB T = 1.2 for densities ρ =
! 0.05, 0.1, 0.2, 0.4, 0.6, 0.8m/σ 3 . Plot your results in Pascal.
!(c) (extra 1) Calculate the mean square displacement over time of argon atoms in a system of N > 100 particles
! with ρ1 = 0.05m/σ 3 and ρ1 = 0.8m/σ 3 at kB T = 1.2. Extract the corresponding diffusion coefficient
! of argon atoms. Show the plot of the mean square displacement (in units of Å2 ) versus time (in units of
! picoseconds). Interpret your results in terms of the phase diagram of the Lennard-Jones fluid shown in
! Fig. 1.
! (d) (extra 2) Obtain the radial distribution function of argon at kB T = 2.0 and ρ = 0.8m/σ 3 . Show the plot
! of the radial distribution function versus distance in units of Ångstrom.

Program P_final_b
	Implicit none
	real(8), parameter :: nu = 0.1, mass = 40, epsilon=0.998 , sigma=3.4 
!	real(8), parameter :: k_b = 1.380649e-23
	integer, parameter :: N=125, Nsteps_ini = 100, Nsteps_prod = 5000
	real(8), dimension(N, 3) :: r, r_ini, vel, vel_ini, r_out
	integer :: step, i, rho_index, Nsteps
	real(8) :: pot, K_energy, L, cutoff, M, a, Temp, inst_temp, dt, absV, p, tini, tfin, rho, sigma_gaussian, T_inst
	real(8), dimension(6) :: rho_list
	integer, allocatable :: seed(:)
	integer :: nn
	external inst_temp

	rho_list = (/ 0.05, 0.1, 0.2, 0.4, 0.6, 0.8 /)
	Temp = 1.2
	sigma_gaussian = Temp**(1.d0/2.d0)
	dt = 1e-4

	call random_seed(size=nn)
	allocate(seed(nn))
	seed = 123456789    ! putting arbitrary seed to all elements
	call random_seed(put=seed)
	deallocate(seed)

	open(44, file="Inicilaization-Verlet.dat")
	open(55, file="thermodynamics.dat")	

	do rho_index = 1, 6

		! system inicialization
		rho = rho_list(rho_index)

		write(44, *) ""
		write(44, *) ""
		write(44, *) "rho = ", rho

		write(55, *) ""
		write(55, *) ""
		write(55, *) "rho = ", rho
		

		L = (N/rho)**(1./3.)

		M = N**(1./3.)
		a = L/M

		print*, L, M, a

		cutoff = L

		call initialize_positions(N, rho, r)

		! Initialize bimodal distrubution: v_i = +- sqrt(T' / m)
		absV = (Temp / mass)**(1./2.)

		call initialize_velocities(N, absV, vel_ini)

		do step = 1,Nsteps_ini
		!	do i = 1, N
		!		print*, i ,vel(i, :)
		!	end do
			call time_step_vVerlet(r, vel, pot, N, L, cutoff, dt)
			call kinetic_energy(vel, K_energy, N)
			write(44,*) step, pot, K_energy, pot+K_energy
			!print*, real(step)/Nsteps
			if (mod(step, 1000).eq.0) then
				print*, real(step)/Nsteps_ini
			end if
		end do

		print*, "PRODUCTION STARTS"
		print*, "-------------------"
		
		do step = 1,Nsteps_prod
			call time_step_vVerlet(r, vel, pot, N, L, cutoff, dt)
			call therm_Andersen(vel, nu, sigma_gaussian, N)
			call kinetic_energy(vel, K_energy, N)
			T_inst = inst_temp(N, K_energy)
			write(55,*) step, pot, K_energy, pot+K_energy, T_inst
			!print*, real(step)/Nsteps
			if (mod(step, 1000).eq.0) then
				print*, real(step)/Nsteps_prod
			end if
		end do

	end do

	close(44)
	close(55)

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


!#################################################################

Subroutine therm_Andersen(vel, nu, sigma_gaussian, N)
	Implicit none
	integer :: i, N
	real(8) :: rand, nu, sigma_gaussian
	real(8), dimension(N, 3) :: vel
	real(8), dimension(2) :: xnums

 	do i = 1, N
 		call random_number(rand)
 		if (rand.lt.nu) then
 			call BM(2, xnums, sigma_gaussian)
 			!print*, "xnums: ", xnums
 			vel(i, 1) = xnums(1)
 			vel(i, 2) = xnums(2)
 			call BM(2, xnums, sigma_gaussian)
 			vel(i, 3) = xnums(1)
		end if
	end do
!	print*, vel
End Subroutine

!#################################################################

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