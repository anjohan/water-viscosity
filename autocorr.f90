module mod_autocorr
    implicit none

    contains
        subroutine calc_autocorr(pressure, mean_autocorr, mean_viscosity, stddev_viscosity, &
                                 window_length, window_distance, num_samples)
            integer, intent(in) :: window_length, window_distance, num_samples
            double precision, intent(in) :: pressure(num_samples, 3)
            double precision, intent(inout) :: mean_autocorr(window_length), &
                                               mean_viscosity(window_length), &
                                               stddev_viscosity(window_length)

            double precision :: tmp_corr(window_length)

            integer :: axis, start_index, num_windows, i

            num_windows = 3*(num_samples - window_length + 1) / window_distance

            !$omp parallel do private(axis, i, start_index, tmp_corr) reduction(+:mean_autocorr) collapse(2)
            do axis = 1, 3
                do i = 1, num_windows/3
                    start_index = (i-1)*window_distance + 1
                    tmp_corr = pressure(start_index, axis)*pressure(start_index:start_index+window_length-1, axis)
                    mean_autocorr(:) = mean_autocorr(:) + tmp_corr(:)
                end do
            end do
            !$omp end parallel do

            mean_autocorr(:) = mean_autocorr(:)/num_windows

            mean_viscosity(:) = cumsum(mean_autocorr)

            !$omp parallel do private(axis, i, tmp_corr, start_index) reduction(+:stddev_viscosity) collapse(2)
            do axis = 1, 3
                do i = 1, num_windows/3
                    start_index = (i-1)*window_distance + 1
                    tmp_corr = pressure(start_index, axis)*pressure(start_index:start_index+window_length-1, axis)
                    stddev_viscosity(:) = stddev_viscosity(:) + (cumsum(tmp_corr) - mean_viscosity)**2
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
