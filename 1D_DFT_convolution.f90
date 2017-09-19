!------------------------------------------------------------------------------------------!
! 1D convolution program using the Convolution Theorem via Discrete Fourier Transforms
!
! The FFTW3 library is required for the DFT and IDFT transforms.
!
! h(x) = (f*g)(x)
! f(x) : Function to convolve
! g(x) : Kernel Gaussian function used for the convolution
! h(x) : Function resulting from the convolution
!
! Compile with:
!   gfortran 1D_DFT_convolution.f90 -fdefault-real-8 -fdefault-integer-8 -Ofast -ffree-line-length-none -I/usr/include -lfftw3 -o 1D_DFT_convolution
!
! Plot function with the python script
!   python3 python/1D_plots.py
!
! Distributed under the GNU GENERAL PUBLIC LICENSE.
!
! Useful link: www.aip.de/groups/soe/local/numres/bookfpdf/f13-1.pdf
!
! Author: B. Font Garcia
! September 2017
!------------------------------------------------------------------------------------------!

module FFTW_helper
	implicit none
	include "fftw3.f"
contains
	function FFT_1D(u) result(uk)
		real, intent(in) :: u(:)
		complex          :: uk(size(u)/2+1)
		integer :: planfw, n

		n = size(u)

		call dfftw_plan_dft_r2c_1d (planfw, n, u, uk, FFTW_ESTIMATE)
		call dfftw_execute_dft_r2c (planfw, u, uk)
		call dfftw_destroy_plan (planfw)
		uk = uk/N    ! Scaling required
	end function FFT_1D

	function IFFT_1D(uk) result(u)
		complex, intent(in) :: uk(:)
		real                :: u((size(uk)-1)*2)
		complex :: uk_copy(size(uk))
		integer :: planbw, n

		n = (size(uk)-1)*2
		uk_copy = uk         ! Create a copy because input is overwritten

		call dfftw_plan_dft_c2r_1d (planbw, n, uk_copy, u, FFTW_ESTIMATE)
		call dfftw_execute_dft_c2r (planbw, uk_copy, u)
		call dfftw_destroy_plan (planbw)
	end function IFFT_1d
end module FFTW_helper

program DFT_convolution
	! External dependencies
	use FFTW_helper
	! Definitions
	implicit none
	complex, parameter :: img = (0.d0,1.d0)    ! sqrt(-1) value
	real, parameter    :: pi = 4.*ATAN(1.)
	! Parameters
	integer, parameter :: n = 256		     ! Function f size
	integer, parameter :: m = 16         ! Kernel g size
	real, parameter    :: sigma = 0.3    ! Kernel standard deviation
	! Useful variables
	integer :: i
	real    :: x, f(n), g(m), h(n)

	write(*,*) '-------------------'
	write(*,*) 'Init 1D convolution'
	write(*,*) '-------------------'

	write(*,*) 'Creating arbitrary function and Gaussian kernel...'
	! Create arbitrary 1D function
	f = 0
	do i = 1,n
		if (i.gt.n/3+1.and.i.le.n/3+1+n/3/2) f(i) = 1.
		if (i.ge.n/3+1+n/3/2.and.i.le.2*n/3-1) f(i) = 0.5
	end do
	! Create Gaussian kernel
	g = 0
	do i = 1,m
		x = -1+2./(m-1)*(i-1)    ! Fit kernel within x = [-1,1] interval
		g(i) = 1./(sigma*sqrt(2.*PI))*exp(-0.5*(x/sigma)**2.)
	end do
	g = g/sum(g)    ! Normalize Gaussian kernel

	! Write kernel function
	write(11,'(a)') 'i x g(x)'
	do i = 1,m
		write(11,'(i4,2f12.4)') i, real(i), g(i)
		flush(11)
	end do

	! Convolute f function with g kernel
	write(*,*) 'Convolving in Fourier space:'
	call convolution(f, g, h)

	! Write convolution result
	write(*,*) 'Writing convolution ouput...'
	write(12,'(a)') 'i x f(x) h(x)'
	do i = 1,n
		write(12,'(i4,3f12.4)') i, real(i), f(i), h(i)
		flush(12)
	end do

	write(*,*) 'Job done!'


contains

	subroutine convolution (f, g, h)
		real, intent(in)  :: f(:), g(:)
		real, intent(out) :: h(size(f))

		integer              :: i, l, m, n
		real, allocatable    :: f_pad(:), g_pad(:), h_padx(:)
		complex, allocatable :: f_padk(:), g_padk(:), h_padk(:)

		n = size(f)
		m = size(g)
		l = n+m-1     ! Size of padded arrays

		write(*,*) ' - Padding function and kernel...'
		allocate(f_pad(-l/2:l/2-1), source=0.)
		allocate(g_pad(-l/2:l/2-1), source=0.)
		allocate(h_padx(-l/2:l/2-1), source=0.)
		allocate(f_padk(l/2+1))
		allocate(g_padk(l/2+1))
		allocate(h_padk(l/2+1))

		f_pad(-n/2:n/2-1) = f       ! Padded and centered function
		g_pad(n/2-1:) = g(1:m/2)    ! Padded and rearrange kernel array to obtain sorted array after convolution
		g_pad(:-n/2) = g(m/2+1:)    ! See this useful link: www.aip.de/groups/soe/local/numres/bookfpdf/f13-1.pdf

		write(*,*) ' - Performing DFTS...'
		f_padk = FFT_1D(f_pad)      ! Transform function to Fourier space
		g_padk = FFT_1D(g_pad)      ! Transform kernel to Fourier space

		write(*,*) ' - Pointwise multiplication...'
		h_padk = f_padk*g_padk*l    ! Scale with l

		write(*,*) ' - Performing IDFT...'
		h_padx = IFFT_1D(h_padk)    ! Transform back to physical space

		h = h_padx(-n/2:n/2-1)      ! Trim padded parts

		! Write convolution result (padded arrays)
		write(13,'(a)') 'i x f_pad(x) g_pad(x) h_pad(x)'
		do i=-l/2,l/2-1
			write(13,'(i4,4f12.4)') i, real(i), f_pad(i), g_pad(i), h_padx(i)
			flush(13)
		end do
	end subroutine convolution
end program DFT_convolution
