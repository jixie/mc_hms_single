! CONSTANTS.INC

! Note: All angles are in radians
!All distances are in cm
!All energies (momenta,masses) are in MeV
!All deltas are in percent
!All densities (thicknesses) are in g/cm3 (g/cm2)
!All B fields are in kG
!All roads lead to infinity

! ... some physical constants
        real*8 Me, Me2, Mp, Mp2, Mn, Mn2, Mpi, Mpi2, Mk, Mk2, Md, Md2
        real*8 Mlambda, Msigma0, Msigma_minus
        real*8 amu, hbarc, pi, twopi, alpha, alpi, degrad, euler
        parameter (Me=0.51099906)
        parameter (Me2=Me**2)
        parameter (Mp=938.27231)
        parameter (Mp2=Mp**2)
        parameter (Mn=939.56563)
        parameter (Mn2=Mn**2)
        parameter (Mpi=139.56995)
        parameter (Mpi2=Mpi**2)
        parameter (Mk=493.677)
        parameter (Mk2=Mk**2)
        parameter (Md=1875.613)
        parameter (Md2=Md**2)
        parameter (Mlambda=1115.68)
        parameter (Msigma0=1192.64)
        parameter (Msigma_minus=1197.45)
        parameter (amu=931.49432)
        parameter (hbarc=197.327053)
        parameter (pi=3.141592653589793)
        parameter (twopi=2.*pi)
        parameter (alpha=1./137.0359895)
        parameter (alpi=alpha/pi)
        parameter (degrad=180./pi)
        parameter (euler=0.577215665)
