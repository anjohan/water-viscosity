module mod_autocorr
    implicit none

    contains
        subroutine calc_autocorr(pressure, autocorrs, viscosities, mean_autocorr, mean_viscosity, stddev_viscosity, &
                                 window_length, window_distance, num_samples)
            integer, intent(in) :: window_length, num_samples, window_distance
            double precision, intent(in) :: pressure(num_samples, 3)
            double precision, intent(inout) :: mean_autocorr(window_length), mean_viscosity(window_length), &
                                               stddev_viscosity(window_length)
            double precision, intent(inout), dimension(window_length, 3*(num_samples - window_length + 1)/window_distance) &
                                            :: autocorrs, viscosities

            integer :: axis, start_index, num_windows, i

            num_windows = size(autocorrs, 2)

            !$omp parallel do private(axis, i, start_index) reduction(+:mean_autocorr) collapse(2)
            do axis = 1, 3
                do i = 1, num_windows/3
                    start_index = (i-1)*window_distance + 1
                    autocorrs(:, 3*(i-1)+axis) = pressure(start_index, axis)*pressure(start_index:start_index+window_length-1, axis)
                    viscosities(:, 3*(i-1)+axis) = cumsum(autocorrs(:, 3*(i-1)+axis))
                    mean_autocorr(:) = mean_autocorr(:) + autocorrs(:, 3*(i-1)+axis)
                end do
            end do
            !$omp end parallel do

            mean_autocorr(:) = mean_autocorr(:)/num_windows

            mean_viscosity(:) = cumsum(mean_autocorr)

            !$omp parallel do private(axis, i) reduction(+:stddev_viscosity) collapse(2)
            do axis = 1, 3
                do i = 1, num_windows/3
                    stddev_viscosity(:) = stddev_viscosity(:) + (viscosities(:, 3*(i-1)+axis) - mean_viscosity)**2
                end do
            end do
            !$omp end parallel do

            stddev_viscosity = sqrt(stddev_viscosity/num_windows)

        end subroutine

        function cumsum(input) result (output)
            double precision, intent(in) :: input(:)
            double precision :: output(size(input))

            integer :: i

            output(1) = 0
            do i = 2, size(input)
                output(i) = output(i-1) + 0.5d0*(input(i) + input(i-1))
            end do
        end function
end module
