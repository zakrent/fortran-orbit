module params
        implicit none
        real(kind=8), parameter :: G = 6.674e-11               ! Grav constant   [m^3*kg^(−1)⋅s^(−2)]
        real(kind=8), parameter :: M = 6e24                    ! Planet mass     [kg]
        real(kind=8), parameter :: dt = 0.01                 ! Step size       [s]
        integer, parameter :: n = 5e7                  ! Number of steps

        real(kind=8), parameter, dimension(2) :: pos_0 = (/0.0, 7.0e6/)! Initial position [m]
        real(kind=8), parameter, dimension(2) :: vel_0 = (/10.0e3, 0.0/)! Initial velocity [m/s]

end module params

module model_vars
        use params

        implicit none
        real(kind=8), dimension(2) :: cur_pos 
        real(kind=8), dimension(2) :: cur_vel

        real(kind=8), dimension(n, 2) :: pos
end module model_vars

program orbit
        use params
        use model_vars

        implicit none

        print *, 'Program started!'
        print *, 'Parameters loaded to the memory!'
        call init
        print *, 'Initialized model variables!'
        print *, 'Starting calculations...'
        call solve
        print *, 'Calculations complete!'
        call save_data
        print *, 'Calculation results saved!'        

end program orbit

subroutine init
        use params
        use model_vars

        implicit none

        cur_pos = pos_0
        cur_vel = vel_0
end subroutine init

subroutine solve
        use params
        use model_vars

        implicit none

        integer :: i 
        do i = 1,n
                call step
                pos(i,:) = cur_pos
        end do 
end subroutine solve

subroutine step
!TODO: Implement rk4 and add air drag
        use params
        use model_vars

        implicit none

        real(kind=8), dimension(2) :: k_v1, k_v2, k_v3, k_v4
        real(kind=8), dimension(2) :: k_r1, k_r2, k_r3, k_r4

        call rk4_step(cur_vel          , cur_pos          , k_v1, k_r1)
        call rk4_step(cur_vel+k_v1*dt/2, cur_pos+k_r1*dt/2, k_v2, k_r2)
        call rk4_step(cur_vel+k_v2*dt/2, cur_pos+k_r2*dt/2, k_v3, k_r3)
        call rk4_step(cur_vel+k_v3*dt  , cur_pos+k_r3*dt  , k_v4, k_r4)

        cur_pos = cur_pos + dt/6*(k_r1+2*k_r2+2*k_r3+k_r4)
        cur_vel = cur_vel + dt/6*(k_v1+2*k_v2+2*k_v3+k_v4)

end subroutine step

subroutine rk4_step(vin, pin, vout, pout)
        use params
        use model_vars
        implicit none

        real(kind=8), dimension(2), intent(in) :: vin, pin
        real(kind=8), dimension(2), intent(out) :: vout, pout

        real(kind=8), dimension(2) :: dir
        real(kind=8) :: r

        r = hypot(pin(1), pin(2))
        vout = -G*M*pin/((r)**3)

        pout = vin 
end subroutine rk4_step

subroutine save_data
        use params
        use model_vars

        implicit none
        integer :: i
        open(1, file = 'out.dat')
        do i = 1,n
                write(1,*) pos(i,1), pos(i,2)
        end do
end subroutine save_data
