
module subroutines_md

implicit none

contains


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
				r(particle, :) = (/ x, y, z /)
				print*, particle, r(particle, :)
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

end module