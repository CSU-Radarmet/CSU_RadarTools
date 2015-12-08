c_____________________________________________________________
c234567891123456789212345678931234567894123456789512345678961234567897123
      subroutine LSE(a, b, x, y, n)

cccc    This is a Least Square Estimate subroutine to fit a linear
cccc    equation for (xi,yi) (i=1,...,n), so that
cccc                            yi = a * xi + b
cccc    INPUTs: x(i), y(i), n, (i=1,...,n ).
cccc    OUTPUTs: a ( slope ), b ( intercept ).
cccc    Li Liu   Sep. 23, 92

      real x(500), y(500), a, b
      real xsum, ysum, xxsum, xysum, det

      xsum = 0.
      ysum = 0.
      xxsum = 0.
      xysum = 0.
      total = float(n)
      do i = 1,n
        if (x(i).gt.1.e35.or.y(i).gt.1.e35) then
          total = total-1.
        else
          xsum =  xsum + x(i)
          ysum =  ysum + y(i)
          xxsum = xxsum + x(i)*x(i)
          xysum = xysum + x(i)*y(i)
        endif
      enddo
      det = total * xxsum - xsum**2
      a = ( total*xysum - xsum*ysum ) / det
      b = ( ysum*xxsum - xsum*xysum ) / det
      return
      end
c-----------------------------------------------------------------

      subroutine calc_kdp_ray_fir(ngates, dp, dz, rng, thsd, nf,
     +bad, fir_order, fir_gain, fir_coeff, kd_lin, dp_lin, sd_lin,
     +std_gate)
c     """
c     Arguments
c     ---------
c     dp = 1D ray of differential phase
c     dz = 1D ray of reflectivity
c     rng = 1D ray of range
c     thsd = Scalar or 1D ray of diff phase stddev thresholds
c     nf = Number of times to filter the data
c     bad = Bad/missing data value
c     fir = Dictionary containing FIR filter parameters
c     std_gate = Number of gates to use for diff phase stddev calc

c     Returns
c     -------
c     kd_lin = Specific differential phase (deg/km, 1D array)
c     dp_lin = Filtered differential phase (deg, 1D array)
c     sd_lin = Standard deviation of diff. phase (deg, 1D array)
c     """
c     # Define needed variables
      integer*4, intent(in) :: ngates
      real*4, intent(in) :: dp(:)
      real*4, intent(in) :: dz(:)
      real*4, intent(in) :: rng(:)
      real*4, intent(in) :: thsd(:)
      integer*4, intent(in) :: nf
      real*4, intent(in) :: bad
      integer, intent(in) :: fir_order
      real*4, intent(in) :: fir_gain
      real*4, intent(in) :: fir_coeff(:)
      integer*4, intent(in) :: std_gate
      real*4, intent(out) :: kd_lin(ngates)
      real*4, intent(out) :: dp_lin(ngates)
      real*4, intent(out) :: sd_lin(ngates)

c     Internal
      real*4 xx(500), y(ngates), yy(500), z(ngates)
      integer*4 half_std_win, i, mloop, nadp, half_fir_win
      integer*4 index1, index2, j, N, half_nadp
      REAL*4 X, A, V, W

c     # Half window size for calculating stdev phase
      half_std_win = (std_gate - 1) / 2
c     Half window size for FIR filtering
      half_fir_win = fir_order / 2

c     #####################################################################
c     # Calculate standard deviation of phidp
      do i = 1, ngates
        kd_lin(i) = bad
        sd_lin(i) = 100.0
        y(i) = bad
        z(i) = dp(i)
        index1 = i - half_std_win
        index2 = i + half_std_win
        if (index1 .ge. 1 .and. index2 .le. ngates) then
          N = 0
          A = 0.0
          V = 0.0
          do j = index1, index2
            if (dp(j) .ne. bad) then
c             Standard deviation algorithm
              X = dp(j)
              IF (N.LE.0) W = X
              N = N + 1
              D = X - W
              V = (N - 1)*(D - A)**2 /N + V
              A = (D - A)/N + A
            endif
          enddo
          if (N .gt. half_std_win) then
            sd_lin(i) = SQRT(V/N)
          endif
        endif
      enddo

c     # ------------- MAIN LOOP of Phidp Adaptive Filtering ------------------
c     # FIR FILTER SECTION
      do mloop = 1, nf
        do i = (half_fir_win+1), (ngates-half_fir_win)
          if ((sd_lin(i) .le. thsd(i)) .and. (z(i) .ne. bad)) then
            icnt = icnt + 1
            index1 = i - half_fir_win
            index2 = i + half_fir_win
            N = 0
            do j = index1, index2
              if ((sd_lin(j) .le. thsd(j)) .and. (z(j) .ne. bad)) then
                N = N + 1
                yy(N) = z(j)
                xx(N) = rng(j)
              endif
            enddo

c           Now fill in gaps if they aren't too big
            if (REAL(N) .gt. (0.8 * REAL(fir_order))) then
              if (N .lt. (fir_order + 1)) then
                call LSE(aa, bb, xx, yy, N)
                do j = index1, index2
                  if (z(j) .eq. bad) then
                    z(j) = aa * rng(j) + bb
                  endif
                enddo
              endif
c             Now do the FIR filtering
              A = 0.0
              do j = index1, index2
                A = A + fir_coeff(j-index1+1) * z(j)
              enddo
              y(i) = A * fir_gain
            endif
          endif
        enddo
        do i = 1, ngates
          z(i) = y(i) ! Enables re-filtering of processed phase
        enddo
      enddo
      do i = 1, ngates
        dp_lin(i) = z(i)
      enddo

c     # *****************END LOOP for Phidp Adaptive Filtering******************

c     # CALCULATE KDP
c     # Default value for nadp is half_fir_win, but varies based on Zh
      do i = 1, ngates
        if (dz(i) .ne. bad) then
          if (dz(i) .ge. 45.0) nadp = half_fir_win
          if ((dz(i) .ge. 35.0) .and. (dz(i) .lt. 45.0)) then
            nadp = 2 * half_fir_win
          endif
          if (dz(i) .lt. 35.0) nadp = 3 * half_fir_win
          half_nadp = nadp / 2
          index1 = i - half_nadp
          index2 = i + half_nadp
          N = 0
          do j = index1, index2
            if (index1 .ge. 1 .and. index2 .le. ngates) then
              if (dp_lin(j) .ne. bad) then
                N = N + 1
                yy(N) = dp_lin(j)
                xx(N) = rng(j)
              endif
            endif
          enddo
          if (REAL(N) .gt. (0.8 * REAL(nadp))) then
            call LSE(aa, bb, xx, yy, N)
            kd_lin(i) = 0.5 * aa
          endif
        endif
      enddo

c     # *******************END KDP CALCULATION****************************
      return
      end

c-----------------------------------------------------------------
c     Beta function calculator

      subroutine hid_beta_f(ngates, x_arr, a, b, m, beta)

      integer*4, intent(in) :: ngates
      real*4, intent(in) :: x_arr(:)
      real*4, intent(in) :: a
      real*4, intent(in) :: b
      real*4, intent(in) :: m

      real*4, intent(out) :: beta(ngates)
      integer*4 i

c      write(*, *) 'called hid_beta_f'
      do i = 1, ngates
        beta(i) = 1.0/(1.0 + (((x_arr(i) - m)/a)**2.0)**b)
      enddo

      return
      end
