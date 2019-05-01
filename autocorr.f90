module mod_autocorr
    implicit none

    contains
        subroutine calc_autocorr(pressure, autocorr, viscosity, stddev, window_length, num_samples)
            integer, intent(in) :: window_length, num_samples
            double precision, intent(in) :: pressure(num_samples, 3)
            double precision, intent(inout) :: autocorr(window_length), viscosity(window_length), &
                                               stddev(window_length)

            integer :: axis, start_index, num_windows
            double precision :: tmp_corr(window_length)

            num_windows = num_samples - window_length + 1

            !$omp parallel do private(axis, start_index) reduction(+:autocorr) collapse(2)
            do axis = 1, 3
                do start_index = 1, num_windows
                    autocorr(:) = autocorr(:) &
                        + pressure(start_index, axis)*pressure(start_index:start_index+window_length-1, axis)
                end do
            end do
            !$omp end parallel do

            autocorr(:) = autocorr(:)/(3*num_windows)

            viscosity(:) = cumsum(autocorr)

            !$omp parallel do private(axis, start_index, tmp_corr) reduction(+:stddev) collapse(2)
            do axis = 1, 3
                do start_index = 1, num_windows
                    tmp_corr(:) = pressure(start_index, axis)*pressure(start_index:start_index+window_length-1, axis)
                    stddev(:) = stddev(:) + (viscosity(:) - cumsum(tmp_corr))**2
                end do
            end do
            !$omp end parallel do

            stddev = sqrt(stddev/(3*num_windows))

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
