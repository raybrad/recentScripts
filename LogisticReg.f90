MODULE lsq

!  Module for unconstrained linear least-squares calculations.
!  The algorithm is suitable for updating LS calculations as more
!  data are added.   This is sometimes called recursive estimation.
!  Only one dependent variable is allowed.
!  Based upon Applied Statistics algorithm AS 274.
!  Translation from Fortran 77 to Fortran 90 by Alan Miller.
!  A function, VARPRD, has been added for calculating the variances
!  of predicted values, and this uses a subroutine BKSUB2.

!  Version 1.14, 19 August 2002 - ELF90 compatible version
!  Author: Alan Miller
!  e-mail : amiller @ bigpond.net.au
!  WWW-pages: http://www.ozemail.com.au/~milleraj
!             http://users.bigpond.net.au/amiller/

!  Bug fixes:
!  1. In REGCF a call to TOLSET has been added in case the user had
!     not set tolerances.
!  2. In SING, each time a singularity is detected, unless it is in the
!     variables in the last position, INCLUD is called.   INCLUD assumes
!     that a new observation is being added and increments the number of
!     cases, NOBS.   The line:  nobs = nobs - 1 has been added.
!  3. row_ptr was left out of the DEALLOCATE statement in routine startup
!     in version 1.07.
!  4. In COV, now calls SS if rss_set = .FALSE.  29 August 1997
!  5. In TOLSET, correction to accomodate negative values of D.  19 August 2002

!  Other changes:
!  1. Array row_ptr added 18 July 1997.   This points to the first element
!     stored in each row thus saving a small amount of time needed to
!     calculate its position.
!  2. Optional parameter, EPS, added to routine TOLSET, so that the user
!     can specify the accuracy of the input data.
!  3. Cosmetic change of lsq_kind to dp (`Double precision')
!  4. Change to routine SING to use row_ptr rather than calculate the position
!     of first elements in each row.

!  The PUBLIC variables are:
!  dp       = a KIND parameter for the floating-point quantities calculated
!             in this module.   See the more detailed explanation below.
!             This KIND parameter should be used for all floating-point
!             arguments passed to routines in this module.

!  nobs    = the number of observations processed to date.
!  ncol    = the total number of variables, including one for the constant,
!            if a constant is being fitted.
!  r_dim   = the dimension of array r = ncol*(ncol-1)/2
!  vorder  = an integer vector storing the current order of the variables
!            in the QR-factorization.   The initial order is 0, 1, 2, ...
!            if a constant is being fitted, or 1, 2, ... otherwise.
!  initialized = a logical variable which indicates whether space has
!                been allocated for various arrays.
!  tol_set = a logical variable which is set when subroutine TOLSET has
!            been called to calculate tolerances for use in testing for
!            singularities.
!  rss_set = a logical variable indicating whether residual sums of squares
!            are available and usable.
!  d()     = array of row multipliers for the Cholesky factorization.
!            The factorization is X = Q.sqrt(D).R where Q is an ortho-
!            normal matrix which is NOT stored, D is a diagonal matrix
!            whose diagonal elements are stored in array d, and R is an
!            upper-triangular matrix with 1's as its diagonal elements.
!  rhs()   = vector of RHS projections (after scaling by sqrt(D)).
!            Thus Q'y = sqrt(D).rhs
!  r()     = the upper-triangular matrix R.   The upper triangle only,
!            excluding the implicit 1's on the diagonal, are stored by
!            rows.
!  tol()   = array of tolerances used in testing for singularities.
!  rss()   = array of residual sums of squares.   rss(i) is the residual
!            sum of squares with the first i variables in the model.
!            By changing the order of variables, the residual sums of
!            squares can be found for all possible subsets of the variables.
!            The residual sum of squares with NO variables in the model,
!            that is the total sum of squares of the y-values, can be
!            calculated as rss(1) + d(1)*rhs(1)^2.   If the first variable
!            is a constant, then rss(1) is the sum of squares of
!            (y - ybar) where ybar is the average value of y.
!  sserr   = residual sum of squares with all of the variables included.
!  row_ptr() = array of indices of first elements in each row of R.
!
!--------------------------------------------------------------------------

!     General declarations

IMPLICIT NONE

INTEGER, SAVE                :: nobs, ncol, r_dim
INTEGER, ALLOCATABLE, SAVE   :: vorder(:), row_ptr(:)
LOGICAL, SAVE                :: initialized = .false.,                  &
                                tol_set = .false., rss_set = .false.

! Note. dp is being set to give at least 12 decimal digit
!       representation of floating point numbers.   This should be adequate
!       for most problems except the fitting of polynomials.   dp is
!       being set so that the same code can be run on PCs and Unix systems,
!       which will usually represent floating-point numbers in `double
!       precision', and other systems with larger word lengths which will
!       give similar accuracy in `single precision'.

INTEGER, PARAMETER           :: dp = SELECTED_REAL_KIND(12,60)
REAL (dp), ALLOCATABLE, SAVE :: d(:), rhs(:), r(:), tol(:), rss(:)
REAL (dp), SAVE              :: zero = 0.0_dp, one = 1.0_dp, vsmall
REAL (dp), SAVE              :: sserr, toly

PUBLIC                       :: dp, nobs, ncol, r_dim, vorder, row_ptr, &
                                initialized, tol_set, rss_set,          &
                                d, rhs, r, tol, rss, sserr
PRIVATE                      :: zero, one, vsmall


CONTAINS

SUBROUTINE startup(nvar, fit_const)

!     Allocates dimensions for arrays and initializes to zero
!     The calling program must set nvar = the number of variables, and
!     fit_const = .true. if a constant is to be included in the model,
!     otherwise fit_const = .false.
!
!--------------------------------------------------------------------------

IMPLICIT NONE
INTEGER, INTENT(IN)  :: nvar
LOGICAL, INTENT(IN)  :: fit_const

!     Local variable
INTEGER   :: i

vsmall = 10. * TINY(zero)

nobs = 0
IF (fit_const) THEN
  ncol = nvar + 1
ELSE
  ncol = nvar
END IF

IF (initialized) DEALLOCATE(d, rhs, r, tol, rss, vorder, row_ptr)
r_dim = ncol * (ncol - 1)/2
ALLOCATE( d(ncol), rhs(ncol), r(r_dim), tol(ncol), rss(ncol), vorder(ncol),  &
          row_ptr(ncol) )

d = zero
rhs = zero
r = zero
sserr = zero

IF (fit_const) THEN
  DO i = 1, ncol
    vorder(i) = i-1
  END DO
ELSE
  DO i = 1, ncol
    vorder(i) = i
  END DO
END IF ! (fit_const)

! row_ptr(i) is the position of element R(i,i+1) in array r().

row_ptr(1) = 1
DO i = 2, ncol-1
  row_ptr(i) = row_ptr(i-1) + ncol - i + 1
END DO
row_ptr(ncol) = 0

initialized = .true.
tol_set = .false.
rss_set = .false.

RETURN
END SUBROUTINE startup




SUBROUTINE includ(weight, xrow, yelem)

!     ALGORITHM AS75.1  APPL. STATIST. (1974) VOL.23, NO. 3

!     Calling this routine updates D, R, RHS and SSERR by the
!     inclusion of xrow, yelem with the specified weight.

!     *** WARNING  Array XROW is overwritten.

!     N.B. As this routine will be called many times in most applications,
!          checks have been eliminated.
!
!--------------------------------------------------------------------------


IMPLICIT NONE
REAL (dp),INTENT(IN)                    :: weight, yelem
REAL (dp), DIMENSION(:), INTENT(IN OUT) :: xrow

!     Local variables

INTEGER     :: i, k, nextr
REAL (dp)   :: w, y, xi, di, wxi, dpi, cbar, sbar, xk

nobs = nobs + 1
w = weight
y = yelem
rss_set = .false.
nextr = 1
DO i = 1, ncol

!     Skip unnecessary transformations.   Test on exact zeroes must be
!     used or stability can be destroyed.

  IF (ABS(w) < vsmall) RETURN
  xi = xrow(i)
  IF (ABS(xi) < vsmall) THEN
    nextr = nextr + ncol - i
  ELSE
    di = d(i)
    wxi = w * xi
    dpi = di + wxi*xi
    cbar = di / dpi
    sbar = wxi / dpi
    w = cbar * w
    d(i) = dpi
    DO k = i+1, ncol
      xk = xrow(k)
      xrow(k) = xk - xi * r(nextr)
      r(nextr) = cbar * r(nextr) + sbar * xk
      nextr = nextr + 1
    END DO
    xk = y
    y = xk - xi * rhs(i)
    rhs(i) = cbar * rhs(i) + sbar * xk
  END IF
END DO ! i = 1, ncol

!     Y * SQRT(W) is now equal to the Brown, Durbin & Evans recursive
!     residual.

sserr = sserr + w * y * y

RETURN
END SUBROUTINE includ



SUBROUTINE regcf(beta, nreq, ifault)

!     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

!     Modified version of AS75.4 to calculate regression coefficients
!     for the first NREQ variables, given an orthogonal reduction from
!     AS75.1.
!
!--------------------------------------------------------------------------

IMPLICIT NONE
INTEGER, INTENT(IN)                  :: nreq
INTEGER, INTENT(OUT)                 :: ifault
REAL (dp), DIMENSION(:), INTENT(OUT) :: beta

!     Local variables

INTEGER   :: i, j, nextr

!     Some checks.

ifault = 0
IF (nreq < 1 .OR. nreq > ncol) ifault = ifault + 4
IF (ifault /= 0) RETURN

IF (.NOT. tol_set) CALL tolset()

DO i = nreq, 1, -1
  IF (SQRT(d(i)) < tol(i)) THEN
    beta(i) = zero
    d(i) = zero
    ifault = -i
  ELSE
    beta(i) = rhs(i)
    nextr = row_ptr(i)
    DO j = i+1, nreq
      beta(i) = beta(i) - r(nextr) * beta(j)
      nextr = nextr + 1
    END DO ! j = i+1, nreq
  END IF
END DO ! i = nreq, 1, -1

RETURN
END SUBROUTINE regcf



SUBROUTINE tolset(eps)

!     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

!     Sets up array TOL for testing for zeroes in an orthogonal
!     reduction formed using AS75.1.

REAL (dp), INTENT(IN), OPTIONAL :: eps

!     Unless the argument eps is set, it is assumed that the input data are
!     recorded to full machine accuracy.   This is often not the case.
!     If, for instance, the data are recorded to `single precision' of about
!     6-7 significant decimal digits, then singularities will not be detected.
!     It is suggested that in this case eps should be set equal to
!     10.0 * EPSILON(1.0)
!     If the data are recorded to say 4 significant decimals, then eps should
!     be set to 1.0E-03
!     The above comments apply to the predictor variables, not to the
!     dependent variable.

!     Correction - 19 August 2002
!     When negative weights are used, it is possible for an alement of D
!     to be negative.

!     Local variables.
!
!--------------------------------------------------------------------------

!     Local variables

INTEGER    :: col, row, pos
REAL (dp)  :: eps1, ten = 10.0, total, work(ncol)

!     EPS is a machine-dependent constant.

IF (PRESENT(eps)) THEN
  eps1 = MAX(ABS(eps), ten * EPSILON(ten))
ELSE
  eps1 = ten * EPSILON(ten)
END IF

!     Set tol(i) = sum of absolute values in column I of R after
!     scaling each element by the square root of its row multiplier,
!     multiplied by EPS1.

work = SQRT(ABS(d))
DO col = 1, ncol
  pos = col - 1
  total = work(col)
  DO row = 1, col-1
    total = total + ABS(r(pos)) * work(row)
    pos = pos + ncol - row - 1
  END DO
  tol(col) = eps1 * total
END DO

tol_set = .TRUE.
RETURN
END SUBROUTINE tolset




SUBROUTINE sing(lindep, ifault)

!     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

!     Checks for singularities, reports, and adjusts orthogonal
!     reductions produced by AS75.1.

!     Correction - 19 August 2002
!     When negative weights are used, it is possible for an alement of D
!     to be negative.

!     Auxiliary routines called: INCLUD, TOLSET
!
!--------------------------------------------------------------------------

INTEGER, INTENT(OUT)                :: ifault
LOGICAL, DIMENSION(:), INTENT(OUT)  :: lindep

!     Local variables

REAL (dp)  :: temp, x(ncol), work(ncol), y, weight
INTEGER    :: pos, row, pos2

ifault = 0

work = SQRT(ABS(d))
IF (.NOT. tol_set) CALL tolset()

DO row = 1, ncol
  temp = tol(row)
  pos = row_ptr(row)         ! pos = location of first element in row

!     If diagonal element is near zero, set it to zero, set appropriate
!     element of LINDEP, and use INCLUD to augment the projections in
!     the lower rows of the orthogonalization.

  lindep(row) = .FALSE.
  IF (work(row) <= temp) THEN
    lindep(row) = .TRUE.
    ifault = ifault - 1
    IF (row < ncol) THEN
      pos2 = pos + ncol - row - 1
      x = zero
      x(row+1:ncol) = r(pos:pos2)
      y = rhs(row)
      weight = d(row)
      r(pos:pos2) = zero
      d(row) = zero
      rhs(row) = zero
      CALL includ(weight, x, y)
                             ! INCLUD automatically increases the number
                             ! of cases each time it is called.
      nobs = nobs - 1
    ELSE
      sserr = sserr + d(row) * rhs(row)**2
    END IF ! (row < ncol)
  END IF ! (work(row) <= temp)
END DO ! row = 1, ncol

RETURN
END SUBROUTINE sing



SUBROUTINE ss()

!     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

!     Calculates partial residual sums of squares from an orthogonal
!     reduction from AS75.1.
!
!--------------------------------------------------------------------------

!     Local variables

INTEGER    :: i
REAL (dp)  :: total

total = sserr
rss(ncol) = sserr
DO i = ncol, 2, -1
  total = total + d(i) * rhs(i)**2
  rss(i-1) = total
END DO

rss_set = .TRUE.
RETURN
END SUBROUTINE ss



SUBROUTINE cov(nreq, var, covmat, dimcov, sterr, ifault)

!     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

!     Calculate covariance matrix for regression coefficients for the
!     first nreq variables, from an orthogonal reduction produced from
!     AS75.1.

!     Auxiliary routine called: INV
!
!--------------------------------------------------------------------------

INTEGER, INTENT(IN)                   :: nreq, dimcov
INTEGER, INTENT(OUT)                  :: ifault
REAL (dp), INTENT(OUT)                :: var
REAL (dp), DIMENSION(:), INTENT(OUT)  :: covmat, sterr

!     Local variables.

INTEGER                :: dim_rinv, pos, row, start, pos2, col, pos1, k
REAL (dp)              :: total
REAL (dp), ALLOCATABLE :: rinv(:)

!     Check that dimension of array covmat is adequate.

IF (dimcov < nreq*(nreq+1)/2) THEN
  ifault = 1
  RETURN
END IF

!     Check for small or zero multipliers on the diagonal.

ifault = 0
DO row = 1, nreq
  IF (ABS(d(row)) < vsmall) ifault = -row
END DO
IF (ifault /= 0) RETURN

!     Calculate estimate of the residual variance.

IF (nobs > nreq) THEN
  IF (.NOT. rss_set) CALL ss()
  var = rss(nreq) / (nobs - nreq)
ELSE
  ifault = 2
  RETURN
END IF

dim_rinv = nreq*(nreq-1)/2
ALLOCATE ( rinv(dim_rinv) )

CALL INV(nreq, rinv)
pos = 1
start = 1
DO row = 1, nreq
  pos2 = start
  DO col = row, nreq
    pos1 = start + col - row
    IF (row == col) THEN
      total = one / d(col)
    ELSE
      total = rinv(pos1-1) / d(col)
    END IF
    DO K = col+1, nreq
      total = total + rinv(pos1) * rinv(pos2) / d(k)
      pos1 = pos1 + 1
      pos2 = pos2 + 1
    END DO ! K = col+1, nreq
    covmat(pos) = total * var
    IF (row == col) sterr(row) = SQRT(covmat(pos))
    pos = pos + 1
  END DO ! col = row, nreq
  start = start + nreq - row
END DO ! row = 1, nreq

DEALLOCATE(rinv)
RETURN
END SUBROUTINE cov



SUBROUTINE inv(nreq, rinv)

!     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

!     Invert first nreq rows and columns of Cholesky factorization
!     produced by AS 75.1.
!
!--------------------------------------------------------------------------

INTEGER, INTENT(IN)                  :: nreq
REAL (dp), DIMENSION(:), INTENT(OUT) :: rinv

!     Local variables.

INTEGER    :: pos, row, col, start, k, pos1, pos2
REAL (dp)  :: total

!     Invert R ignoring row multipliers, from the bottom up.

pos = nreq * (nreq-1)/2
DO row = nreq-1, 1, -1
  start = row_ptr(row)
  DO col = nreq, row+1, -1
    pos1 = start
    pos2 = pos
    total = zero
    DO k = row+1, col-1
      pos2 = pos2 + nreq - k
      total = total - r(pos1) * rinv(pos2)
      pos1 = pos1 + 1
    END DO ! k = row+1, col-1
    rinv(pos) = total - r(pos1)
    pos = pos - 1
  END DO ! col = nreq, row+1, -1
END DO ! row = nreq-1, 1, -1

RETURN
END SUBROUTINE inv



SUBROUTINE partial_corr(in, cormat, dimc, ycorr, ifault)

!     Replaces subroutines PCORR and COR of:
!     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

!     Calculate partial correlations after the variables in rows
!     1, 2, ..., IN have been forced into the regression.
!     If IN = 1, and the first row of R represents a constant in the
!     model, then the usual simple correlations are returned.

!     If IN = 0, the value returned in array CORMAT for the correlation
!     of variables Xi & Xj is:
!       sum ( Xi.Xj ) / Sqrt ( sum (Xi^2) . sum (Xj^2) )

!     On return, array CORMAT contains the upper triangle of the matrix of
!     partial correlations stored by rows, excluding the 1's on the diagonal.
!     e.g. if IN = 2, the consecutive elements returned are:
!     (3,4) (3,5) ... (3,ncol), (4,5) (4,6) ... (4,ncol), etc.
!     Array YCORR stores the partial correlations with the Y-variable
!     starting with YCORR(IN+1) = partial correlation with the variable in
!     position (IN+1).
!
!--------------------------------------------------------------------------

INTEGER, INTENT(IN)                  :: in, dimc
INTEGER, INTENT(OUT)                 :: ifault
REAL (dp), DIMENSION(:), INTENT(OUT) :: cormat, ycorr

!     Local variables.

INTEGER    :: base_pos, pos, row, col, col1, col2, pos1, pos2
REAL (dp)  :: rms(in+1:ncol), sumxx, sumxy, sumyy, work(in+1:ncol)

!     Some checks.

ifault = 0
IF (in < 0 .OR. in > ncol-1) ifault = ifault + 4
IF (dimc < (ncol-in)*(ncol-in-1)/2) ifault = ifault + 8
IF (ifault /= 0) RETURN

!     Base position for calculating positions of elements in row (IN+1) of R.

base_pos = in*ncol - (in+1)*(in+2)/2

!     Calculate 1/RMS of elements in columns from IN to (ncol-1).

IF (d(in+1) > zero) rms(in+1) = one / SQRT(d(in+1))
DO col = in+2, ncol
  pos = base_pos + col
  sumxx = d(col)
  DO row = in+1, col-1
    sumxx = sumxx + d(row) * r(pos)**2
    pos = pos + ncol - row - 1
  END DO ! row = in+1, col-1
  IF (sumxx > zero) THEN
    rms(col) = one / SQRT(sumxx)
  ELSE
    rms(col) = zero
    ifault = -col
  END IF ! (sumxx > zero)
END DO ! col = in+1, ncol-1

!     Calculate 1/RMS for the Y-variable

sumyy = sserr
DO row = in+1, ncol
  sumyy = sumyy + d(row) * rhs(row)**2
END DO ! row = in+1, ncol
IF (sumyy > zero) sumyy = one / SQRT(sumyy)

!     Calculate sums of cross-products.
!     These are obtained by taking dot products of pairs of columns of R,
!     but with the product for each row multiplied by the row multiplier
!     in array D.

pos = 1
DO col1 = in+1, ncol
  sumxy = zero
  work(col1+1:ncol) = zero
  pos1 = base_pos + col1
  DO row = in+1, col1-1
    pos2 = pos1 + 1
    DO col2 = col1+1, ncol
      work(col2) = work(col2) + d(row) * r(pos1) * r(pos2)
      pos2 = pos2 + 1
    END DO ! col2 = col1+1, ncol
    sumxy = sumxy + d(row) * r(pos1) * rhs(row)
    pos1 = pos1 + ncol - row - 1
  END DO ! row = in+1, col1-1

!     Row COL1 has an implicit 1 as its first element (in column COL1)

  pos2 = pos1 + 1
  DO col2 = col1+1, ncol
    work(col2) = work(col2) + d(col1) * r(pos2)
    pos2 = pos2 + 1
    cormat(pos) = work(col2) * rms(col1) * rms(col2)
    pos = pos + 1
  END DO ! col2 = col1+1, ncol
  sumxy = sumxy + d(col1) * rhs(col1)
  ycorr(col1) = sumxy * rms(col1) * sumyy
END DO ! col1 = in+1, ncol-1

ycorr(1:in) = zero

RETURN
END SUBROUTINE partial_corr




SUBROUTINE vmove(from, to, ifault)

!     ALGORITHM AS274 APPL. STATIST. (1992) VOL.41, NO. 2

!     Move variable from position FROM to position TO in an
!     orthogonal reduction produced by AS75.1.
!
!--------------------------------------------------------------------------

INTEGER, INTENT(IN)    :: from, to
INTEGER, INTENT(OUT)   :: ifault

!     Local variables

REAL (dp)  :: d1, d2, x, d1new, d2new, cbar, sbar, y
INTEGER    :: m, first, last, inc, m1, m2, mp1, col, pos, row

!     Check input parameters

ifault = 0
IF (from < 1 .OR. from > ncol) ifault = ifault + 4
IF (to < 1 .OR. to > ncol) ifault = ifault + 8
IF (ifault /= 0) RETURN

IF (from == to) RETURN

IF (.NOT. rss_set) CALL ss()

IF (from < to) THEN
  first = from
  last = to - 1
  inc = 1
ELSE
  first = from - 1
  last = to
  inc = -1
END IF

DO m = first, last, inc

!     Find addresses of first elements of R in rows M and (M+1).

  m1 = row_ptr(m)
  m2 = row_ptr(m+1)
  mp1 = m + 1
  d1 = d(m)
  d2 = d(mp1)

!     Special cases.

  IF (d1 < vsmall .AND. d2 < vsmall) GO TO 40
  x = r(m1)
  IF (ABS(x) * SQRT(d1) < tol(mp1)) THEN
    x = zero
  END IF
  IF (d1 < vsmall .OR. ABS(x) < vsmall) THEN
    d(m) = d2
    d(mp1) = d1
    r(m1) = zero
    DO col = m+2, ncol
      m1 = m1 + 1
      x = r(m1)
      r(m1) = r(m2)
      r(m2) = x
      m2 = m2 + 1
    END DO ! col = m+2, ncol
    x = rhs(m)
    rhs(m) = rhs(mp1)
    rhs(mp1) = x
    GO TO 40
  ELSE IF (d2 < vsmall) THEN
    d(m) = d1 * x**2
    r(m1) = one / x
    r(m1+1:m1+ncol-m-1) = r(m1+1:m1+ncol-m-1) / x
    rhs(m) = rhs(m) / x
    GO TO 40
  END IF

!     Planar rotation in regular case.

  d1new = d2 + d1*x**2
  cbar = d2 / d1new
  sbar = x * d1 / d1new
  d2new = d1 * cbar
  d(m) = d1new
  d(mp1) = d2new
  r(m1) = sbar
  DO col = m+2, ncol
    m1 = m1 + 1
    y = r(m1)
    r(m1) = cbar*r(m2) + sbar*y
    r(m2) = y - x*r(m2)
    m2 = m2 + 1
  END DO ! col = m+2, ncol
  y = rhs(m)
  rhs(m) = cbar*rhs(mp1) + sbar*y
  rhs(mp1) = y - x*rhs(mp1)

!     Swap columns M and (M+1) down to row (M-1).

  40 pos = m
  DO row = 1, m-1
    x = r(pos)
    r(pos) = r(pos-1)
    r(pos-1) = x
    pos = pos + ncol - row - 1
  END DO ! row = 1, m-1

!     Adjust variable order (VORDER), the tolerances (TOL) and
!     the vector of residual sums of squares (RSS).

  m1 = vorder(m)
  vorder(m) = vorder(mp1)
  vorder(mp1) = m1
  x = tol(m)
  tol(m) = tol(mp1)
  tol(mp1) = x
  rss(m) = rss(mp1) + d(mp1) * rhs(mp1)**2
END DO

RETURN
END SUBROUTINE vmove



SUBROUTINE reordr(list, n, pos1, ifault)

!     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2

!     Re-order the variables in an orthogonal reduction produced by
!     AS75.1 so that the N variables in LIST start at position POS1,
!     though will not necessarily be in the same order as in LIST.
!     Any variables in VORDER before position POS1 are not moved.

!     Auxiliary routine called: VMOVE
!
!--------------------------------------------------------------------------

INTEGER, INTENT(IN)               :: n, pos1
INTEGER, DIMENSION(:), INTENT(IN) :: list
INTEGER, INTENT(OUT)              :: ifault

!     Local variables.

INTEGER    :: next, i, l, j

!     Check N.

ifault = 0
IF (n < 1 .OR. n > ncol+1-pos1) ifault = ifault + 4
IF (ifault /= 0) RETURN

!     Work through VORDER finding variables which are in LIST.

next = pos1
i = pos1
10 l = vorder(i)
DO j = 1, n
  IF (l == list(j)) GO TO 40
END DO
30 i = i + 1
IF (i <= ncol) GO TO 10

!     If this point is reached, one or more variables in LIST has not
!     been found.

ifault = 8
RETURN

!     Variable L is in LIST; move it up to position NEXT if it is not
!     already there.

40 IF (i > next) CALL vmove(i, next, ifault)
next = next + 1
IF (next < n+pos1) GO TO 30

RETURN
END SUBROUTINE reordr



SUBROUTINE hdiag(xrow, nreq, hii, ifault)

!     ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2
!
!                         -1           -1
! The hat matrix H = x(X'X) x' = x(R'DR) x' = z'Dz

!              -1
! where z = x'R

! Here we only calculate the diagonal element hii corresponding to one
! row (xrow).   The variance of the i-th least-squares residual is (1 - hii).
!--------------------------------------------------------------------------

INTEGER, INTENT(IN)                  :: nreq
INTEGER, INTENT(OUT)                 :: ifault
REAL (dp), DIMENSION(:), INTENT(IN)  :: xrow
REAL (dp), INTENT(OUT)               :: hii

!     Local variables

INTEGER    :: col, row, pos
REAL (dp)  :: total, wk(ncol)

!     Some checks

ifault = 0
IF (nreq > ncol) ifault = ifault + 4
IF (ifault /= 0) RETURN

!     The elements of xrow.inv(R).sqrt(D) are calculated and stored in WK.

hii = zero
DO col = 1, nreq
  IF (SQRT(d(col)) <= tol(col)) THEN
    wk(col) = zero
  ELSE
    pos = col - 1
    total = xrow(col)
    DO row = 1, col-1
      total = total - wk(row)*r(pos)
      pos = pos + ncol - row - 1
    END DO ! row = 1, col-1
    wk(col) = total
    hii = hii + total**2 / d(col)
  END IF
END DO ! col = 1, nreq

RETURN
END SUBROUTINE hdiag



FUNCTION varprd(x, nreq) RESULT(fn_val)

!     Calculate the variance of x'b where b consists of the first nreq
!     least-squares regression coefficients.
!
!--------------------------------------------------------------------------

INTEGER, INTENT(IN)                  :: nreq
REAL (dp), DIMENSION(:), INTENT(IN)  :: x
REAL (dp)                            :: fn_val

!     Local variables

INTEGER    :: ifault, row
REAL (dp)  :: var, wk(nreq)

!     Check input parameter values

fn_val = zero
ifault = 0
IF (nreq < 1 .OR. nreq > ncol) ifault = ifault + 4
IF (nobs <= nreq) ifault = ifault + 8
IF (ifault /= 0) THEN
  WRITE(*, '(1x, a, i4)') 'Error in function VARPRD: ifault =', ifault
  RETURN
END IF

!     Calculate the residual variance estimate.

var = sserr / (nobs - nreq)

!     Variance of x'b = var.x'(inv R)(inv D)(inv R')x
!     First call BKSUB2 to calculate (inv R')x by back-substitution.

CALL BKSUB2(x, wk, nreq)
DO row = 1, nreq
  IF(d(row) > tol(row)) fn_val = fn_val + wk(row)**2 / d(row)
END DO

fn_val = fn_val * var

RETURN
END FUNCTION varprd



SUBROUTINE bksub2(x, b, nreq)

!     Solve x = R'b for b given x, using only the first nreq rows and
!     columns of R, and only the first nreq elements of R.
!
!--------------------------------------------------------------------------

INTEGER, INTENT(IN)                  :: nreq
REAL (dp), DIMENSION(:), INTENT(IN)  :: x
REAL (dp), DIMENSION(:), INTENT(OUT) :: b

!     Local variables

INTEGER    :: pos, row, col
REAL (dp)  :: temp

!     Solve by back-substitution, starting from the top.

DO row = 1, nreq
  pos = row - 1
  temp = x(row)
  DO col = 1, row-1
    temp = temp - r(pos)*b(col)
    pos = pos + ncol - col - 1
  END DO
  b(row) = temp
END DO

RETURN
END SUBROUTINE bksub2

END MODULE lsq

MODULE Logistic_Regression
! A package for fitting linear logistic models by iteratively re-weighted
! least squares.

! The model fitted is that the probability of a `success' is:

!      p = 1 - 1/[1 + exp(b0 + b1.X1 + b2.X2 + ... + bk.Xk)]

! where X1, X2, ... , Xk are the predictor variables, and the coefficients
! b0, b1, b2, ..., bk are to be determined.

! N.B. The residual variance used here in estimating standard errors is the
!      larger of 1.0 and that from the weighted least squares calculations;
!      it is not the theoretical residual variance (1.0) assuming a binomial
!      distribution about the logistic curve.   If the fit of the logistic
!      is poor, the standard errors of the coefficients in the logistic will
!      be much larger.

! The calculation of chi-squared was corrected on 9 August 2003.
! My thanks to Zhang Guan

! By Alan Miller
! amiller @ bigpond.net.au
! users.bigpond.net.au/amiller/

! Latest revision - 9 August 2003

USE lsq
IMPLICIT NONE

PUBLIC  :: logistic

CONTAINS


SUBROUTINE logistic(ngroups, x, k, s, n, chisq, devnce, ndf, beta, se_beta,  &
                    ier, cov_beta, fit, stdres)

! Input arguments:

! ngroups     The number of groups of observations.

! x(:,:)      x(i,j) contains the value of the j-th predictor for group i.
!             The X-values are the same for all the n(i) trials in group i.

! k           The number of predictors.   N.B. A constant will be fitted in
!             the model; do not include it in the k predictor variables.

! s(:), n(:)  s(i) is the number of `successes' out of the n(i) `trials' in
!             the i-th group.   In many applications, each n(i) = 1, that is
!             each case will be treated individually, and then s(i) = 0 or 1
!             for each of the two possible outcomes.

INTEGER, INTENT(IN)    :: ngroups, k, s(:), n(:)
REAL (dp), INTENT(IN)  :: x(:,:)

! Output arguments:

! chisq       The value of chi-squared on exit, when a model has been fitted.

! devnce      The deviance on exit, when a model has been fitted.

! ndf         Number of degrees of freedom.

! beta(0:)    The fitted coefficients in the logistic model.

! se_beta(0:) Approximate standard errors of the beta coefficients.

! ier         Error indicator
!             = 0 successful termination
!             = 1 if ngroups < 2 or ndf < 0
!             = 2 if any n(i) < 0
!             = 3 if any r(i) < 0
!             = 4 if any r(i) > n(i)
!             = 5 if any X-variable is constant
!             = 6 if a singularity is detected
!             = 7 if any beta(i) is tending to +/- infinity
!             = 8 failure to converge in 20 iterations

REAL (dp), INTENT(OUT)  :: chisq, devnce, beta(0:), se_beta(0:)
INTEGER, INTENT(OUT)    :: ndf, ier

! Optional output arguments:

! cov_beta(0:,0:)     Approximate covariance matrix of the fitted coefficients.

! fit(:)      The fitted probabilities of a success for each group.

! stdres(:)   Vector of standardized residuals.

REAL (dp), INTENT(OUT), OPTIONAL  :: cov_beta(0:,0:), fit(:), stdres(:)


! Local variables

INTEGER    :: i, iter, j, ncov, pos
REAL (dp)  :: propn(ngroups), p(ngroups), wt(ngroups), xrow(0:k), db(0:k), &
              bnew(0:k), dev_new, xb, pnew(ngroups), wnew(ngroups), a, b,  &
              range(k), var, e(ngroups), hii
LOGICAL    :: lindep(0:k)
REAL (dp), ALLOCATABLE  :: covmat(:)

! Initial checks

ier = 0
ndf = ngroups - k - 1
IF (ngroups < 2 .OR. ndf < 0) THEN
  ier = 1
  RETURN
END IF
IF (ANY(n(1:ngroups) < 0)) THEN
  ier = 2
  RETURN
END IF
IF (ANY(s(1:ngroups) < 0)) THEN
  ier = 3
  RETURN
END IF
IF (ANY(s(1:ngroups) > n(1:ngroups))) THEN
  ier = 4
  RETURN
END IF

! Calculate ranges of the X-variables to use in testing whether any beta
! is tending to +/- infinity.  Also test that no variable is constant.

DO i = 1, k
  a = MAXVAL(x(1:ngroups,i))
  b = MINVAL(x(1:ngroups,i))
  range(i) = a - b
  IF (range(i) < EPSILON(0.0_dp) * (ABS(a) + ABS(b))) THEN
    ier = 5
    RETURN
  END IF
END DO

! Start with all beta's = 0 and weights = 1.

beta(0:k) = 0.0_dp
wt(1:ngroups) = 1.0_dp
p(1:ngroups) = 0.5_dp

! propn stores the sample proportions, i.e. s(i) / n(i)
propn(1:ngroups) = REAL(s(1:ngroups), KIND=dp) / n(1:ngroups)
iter = 1

! Start of iterative cycle

DO
  CALL startup(k, .TRUE.)
  DO i = 1, ngroups
    IF (iter == 1) THEN
      xrow(0) = 0.25_dp
      xrow(1:k) = 0.25_dp*x(i, 1:k)
    ELSE
      xrow(0) = p(i)*(1.0_dp - p(i))
      xrow(1:k) = p(i)*(1.0_dp - p(i))*x(i, 1:k)
    END IF
    CALL includ(wt(i), xrow, propn(i)-p(i))
  END DO

! Test for a singularity

  CALL sing(lindep, ier)
  IF (ier /= 0) THEN
    DO i = 1, k
      IF (lindep(i)) THEN
        WRITE(*, '(a, i6, a)') ' Variable number ', i,  &
                               ' is linearly dependent upon earlier variables'
      END IF
    END DO
    ier = 6
    RETURN
  END IF

  CALL regcf(db, k+1, ier)
  10 bnew = beta(0:k) + db

! Calculate new p(i)'s, weights & deviance

  dev_new = 0.0_dp
  DO i = 1, ngroups
    xb = DOT_PRODUCT( x(i,1:k), bnew(1:k) ) + bnew(0)
    xb = EXP(xb)
    pnew(i) = xb / (1.0_dp + xb)
    wnew(i) = REAL(n(i), KIND=dp) / (pnew(i)*(1.0_dp - pnew(i)))
    IF (iter == 1) wnew(i) = SQRT(wnew(i))
    IF (s(i) > 0) dev_new = dev_new + s(i)*LOG(propn(i)/pnew(i))
    IF (s(i) < n(i)) dev_new = dev_new +   &
                           (n(i)-s(i))*LOG((1.0_dp-propn(i))/(1.0_dp-pnew(i)))
  END DO
  dev_new = 2 * dev_new

! If deviance has increased, reduce the step size.

  IF (iter > 2) THEN
    IF (dev_new > devnce*1.0001_dp) THEN
      db = 0.5_dp * db
      GO TO 10
    END IF
  END IF

! Replace betas, weights & p's with new values

  beta(0:k) = bnew(0:k)
  wt = wnew
  p(1:ngroups) = pnew

! Test for convergence

  IF (iter > 2 .AND. devnce - dev_new < 0.0001_dp) EXIT
  devnce = dev_new
  iter = iter + 1
  IF (iter > 20) THEN
    ier = 8
    RETURN
  END IF

! Test for a very large beta

  DO i = 1, k
    IF (ABS(beta(i))*range(i) > 30.0_dp) THEN
      WRITE(*, '(a, i4, a)') ' Coefficient for variable no.', i,  &
                             ' tending to infinity'
      ier = 7
      RETURN
    END IF
  END DO

END DO

e = n(1:ngroups)*p(1:ngroups)
chisq = SUM( (s(1:ngroups) - e)**2 / (e * (1.0_dp - p(1:ngroups))) )
devnce = dev_new

! Calculate the approximate covariance matrix for the beta's, if ndf > 0.

IF (ndf > 0) THEN
  ncov = (k+1)*(k+2)/2
  ALLOCATE( covmat(ncov) )
  CALL cov(k+1, var, covmat, ncov, se_beta, ier)
  IF (var < 1.0_dp) THEN
    covmat = covmat / var
    se_beta = se_beta / SQRT(var)
  END IF

  IF(PRESENT(cov_beta)) THEN
    pos = 1
    DO i = 0, k
      cov_beta(i,i) = covmat(pos)
      pos = pos + 1
      DO j = i+1, k
        cov_beta(i,j) = covmat(pos)
        cov_beta(j,i) = covmat(pos)
        pos = pos + 1
      END DO
    END DO
  END IF
END IF

IF(PRESENT(fit)) fit(1:ngroups) = p

IF (PRESENT(stdres)) THEN
  DO i = 1, ngroups
    xrow(0) = p(i)*(1.0_dp - p(i))
    xrow(1:k) = p(i)*(1.0_dp - p(i))*x(i, 1:k)
    CALL hdiag(xrow, k+1, hii, ier)
    stdres(i) = (s(i) - n(i)*p(i)) /   &
                SQRT(n(i)*p(i)*(1.0_dp - p(i))*(1.0_dp - hii))
  END DO
END IF

IF (ALLOCATED(covmat)) DEALLOCATE( covmat )

RETURN
END SUBROUTINE logistic



FUNCTION chi_squared(ndf, chi2) RESULT(prob)
! Calculate the chi-squared distribution function
! ndf  = number of degrees of freedom
! chi2 = chi-squared value
! prob = probability of a chi-squared value <= chi2 (i.e. the left-hand
!        tail area)

INTEGER, INTENT(IN)    :: ndf
REAL (dp), INTENT(IN)  :: chi2
REAL (dp)              :: prob

! Local variables
REAL (dp) :: half = 0.5_dp, x, p

x = half * chi2
p = half * REAL(ndf)
prob = gammad(x, p)
RETURN

END FUNCTION chi_squared



FUNCTION gammad(x, p) RESULT(gamma_prob)

!  ALGORITHM AS239  APPL. STATIST. (1988) VOL. 37, NO. 3

!  Computation of the Incomplete Gamma Integral

!  Auxiliary functions required: ALNORM = algorithm AS66 (included) & LNGAMMA

!  Converted to be compatible with ELF90 by Alan Miller
!  N.B. The return parameter IFAULT has been removed as ELF90 allows only
!  one output parameter from functions.   An error message is issued instead.

! This revision - 15 October 1996

REAL (dp), INTENT(IN) :: x, p
REAL (dp)             :: gamma_prob

!     Local variables

REAL (dp) :: pn1, pn2, pn3, pn4, pn5, pn6, tol = 1.d-14, oflo = 1.d+37,  &
             xbig = 1.d+8, arg, c, rn, a, b, one = 1._dp, zero = 0._dp, an, &
             two = 2._dp, elimit = -88._dp, plimit = 1000._dp, three = 3._dp, &
             nine = 9._dp

gamma_prob = zero

!      Check that we have valid values for X and P

IF (p <= zero .OR. x < zero) THEN
  WRITE(*, *)'Error: Function gammad.  1st argument < 0 or 2nd argument <= 0'
  RETURN
END IF
IF (x == zero) RETURN

!      Use a normal approximation if P > PLIMIT

IF (p > plimit) THEN
  pn1 = three * SQRT(p) * ((x / p) ** (one / three) + one / (nine * p) - one)
  gamma_prob = alnorm(pn1, .false.)
  RETURN
END IF

!      If X is extremely large compared to P then set gamma_prob = 1

IF (x > xbig) THEN
  gamma_prob = one
  RETURN
END IF

IF (x <= one .OR. x < p) THEN

!      Use Pearson's series expansion.
!      (Note that P is not large enough to force overflow in LNGAMMA)

  arg = p * LOG(x) - x - lngamma(p + one)
  c = one
  gamma_prob = one
  a = p
  DO
    a = a + one
    c = c * x / a
    gamma_prob = gamma_prob + c
    IF (c < tol) EXIT
  END DO
  arg = arg + LOG(gamma_prob)
  gamma_prob = zero
  IF (arg >= elimit) gamma_prob = EXP(arg)

ELSE

!      Use a continued fraction expansion

  arg = p * LOG(x) - x - lngamma(p)
  a = one - p
  b = a + x + one
  c = zero
  pn1 = one
  pn2 = x
  pn3 = x + one
  pn4 = x * b
  gamma_prob = pn3 / pn4
  DO
    a = a + one
    b = b + two
    c = c + one
    an = a * c
    pn5 = b * pn3 - an * pn1
    pn6 = b * pn4 - an * pn2
    IF (ABS(pn6) > zero) THEN
      rn = pn5 / pn6
      IF (ABS(gamma_prob - rn) <= MIN(tol, tol * rn)) EXIT
      gamma_prob = rn
    END IF

    pn1 = pn3
    pn2 = pn4
    pn3 = pn5
    pn4 = pn6
    IF (ABS(pn5) >= oflo) THEN

  !      Re-scale terms in continued fraction if terms are large

      pn1 = pn1 / oflo
      pn2 = pn2 / oflo
      pn3 = pn3 / oflo
      pn4 = pn4 / oflo
    END IF
  END DO
  arg = arg + LOG(gamma_prob)
  gamma_prob = one
  IF (arg >= elimit) gamma_prob = one - EXP(arg)
END IF

RETURN
END FUNCTION gammad

FUNCTION alnorm(x, upper) RESULT(norm_prob)

!  Algorithm AS66 Applied Statistics (1973) vol.22, no.3

!  Evaluates the tail area of the standardised normal curve
!  from x to infinity if upper is .true. or
!  from minus infinity to x if upper is .false.

REAL (dp), INTENT(IN)  :: x
LOGICAL, INTENT(IN)    :: upper
REAL (dp)              :: norm_prob


! Local variables
REAL (dp) :: zero = 0.0_dp, one = 1.0_dp, half = 0.5_dp
REAL (dp) :: con = 1.28_dp, z, y, ltone = 7.0_dp, utzero = 18.66_dp
REAL (dp) :: p = 0.398942280444_dp, q = 0.39990348504_dp,   &
             r = 0.398942280385_dp, a1 = 5.75885480458_dp,  &
             a2 = 2.62433121679_dp, a3 = 5.92885724438_dp,  &
             b1 = -29.8213557807_dp, b2 = 48.6959930692_dp, &
             c1 = -3.8052D-8, c2 = 3.98064794D-4,         &
             c3 = -0.151679116635_dp, c4 = 4.8385912808_dp, &
             c5 = 0.742380924027_dp, c6 = 3.99019417011_dp, &
             d1 = 1.00000615302_dp, d2 = 1.98615381364_dp,  &
             d3 = 5.29330324926_dp, d4 = -15.1508972451_dp, &
             d5 = 30.789933034_dp
LOGICAL   :: up

up = upper
z = x
IF(z >=  zero) GO TO 10
up = .NOT. up
z = -z
10 IF(z <= ltone .OR. up .AND. z <= utzero) GO TO 20
norm_prob = zero
GO TO 40
20 y = half*z*z
IF(z > con) GO TO 30

norm_prob = half - z*(p-q*y/(y+a1+b1/(y+a2+b2/(y+a3))))
GO TO 40
30 norm_prob = r*EXP(-y)/(z+c1+d1/(z+c2+d2/(z+c3+d3/(z+c4+d4/(z+c5+d5/(z+c6))))))
40 IF(.NOT. up) norm_prob = one - norm_prob
RETURN

END FUNCTION alnorm

FUNCTION lngamma(z) RESULT(lanczos)

!  Uses Lanczos-type approximation to ln(gamma) for z > 0.
!  Reference:
!       Lanczos, C. 'A precision approximation of the gamma
!               function', J. SIAM Numer. Anal., B, 1, 86-96, 1964.
!  Accuracy: About 14 significant digits except for small regions
!            in the vicinity of 1 and 2.

!  Programmer: Alan Miller
!              1 Creswick Street, Brighton, Vic. 3187, Australia
!  Latest revision - 14 October 1996

REAL(dp), INTENT(IN) :: z
REAL(dp)             :: lanczos

! Local variables

REAL(dp)  :: a(9) = (/ 0.9999999999995183_dp, 676.5203681218835_dp, &
                      -1259.139216722289_dp, 771.3234287757674_dp, &
                      -176.6150291498386_dp, 12.50734324009056_dp, &
                      -0.1385710331296526_dp, 0.9934937113930748D-05, &
                       0.1659470187408462D-06 /), zero = 0._dp,   &
                       one = 1._dp, lnsqrt2pi = 0.9189385332046727_dp, &
                       half = 0.5_dp, sixpt5 = 6.5_dp, seven = 7._dp, tmp
INTEGER   :: j

IF (z <= zero) THEN
  WRITE(*, *)'Error: zero or -ve argument for lngamma'
  RETURN
END IF

lanczos = zero
tmp = z + seven
DO j = 9, 2, -1
  lanczos = lanczos + a(j)/tmp
  tmp = tmp - one
END DO
lanczos = lanczos + a(1)
lanczos = LOG(lanczos) + lnsqrt2pi - (z + sixpt5) + (z - half)*LOG(z + sixpt5)
RETURN

END FUNCTION lngamma


END MODULE Logistic_Regression

program fitFermi
use Logistic_Regression
implicit none

integer::i,j
character*20::arg(3),FileName
integer:: n
integer:: Varn(:),Vars(:)
real*8 ::ein,eFermi,occ,b,a,corr,temp
real*8,allocatable :: energy0(:),energy(:),x(:)
real*8, parameter :: boltz = 1.3806506d-23
real*8, parameter :: ev2j  = 1.60217646d-19
real*8::chisq, devnce, beta(0:1), se_beta(0:1)
integer::ndf, ier

eFermi=-4.418
n=1

CALL getarg(1, arg(1))      
CALL getarg(2, arg(2))      
CALL getarg(3, arg(3))      

!input must be *.* real format or there will be error
read(arg(1),'(F12.5)') eFermi
read(arg(2),'(I5)') n
read(arg(3),*) FileName 

write(6,*) 'fermi energy:',eFermi
write(6,*) 'line No. :',n
write(6,*) 'processing ', FileName

allocate(energy0(n),energy(n),x(n,2),Varn(n),Vars(n))
energy0=0.d0
energy =0.d0
Varn=1
Vars=1
open(14,file='diagMO.base')
do i=1,n
read(14,*) energy0(i),temp
enddo

open(16,file=trim(FileName))
do i=1,n
    read(16,*) energy(i), occ
    if(occ>2.0) then
         write(6,*) 'occ approx' , occ
!        cycle
         occ=2.d0-1.d-10
    elseif(occ<0.d0) then
         write(6,*) 'occ approx' , occ
         occ=1.d-10
         !cycle
    endif
x(i,1)=energy0(i)
x(i,2)=occ
!
! fermidirac = 1.d0 / ( 1.d0 + exp(ein / kT) )
!
enddo
close(16)

CALL logistic(n, x, 2, Vars, Varn, chisq, devnce, ndf, beta, se_beta, ier)  
write(6,*) 'chisq',chisq
write(6,*) 'devnce',devnce
write(6,*) 'ndf',ndf
write(6,*) 'beta',beta
write(6,*) 'se_beta'
write(6,*) 'ier',ier

!     Update the average

! open(15, file ='temp.dat')
! write(15,'(2(E20.12,3X))') x , y 
! close(15)

end program fitFermi

