c ======================================================================
c User Subroutine VUMAT for Abaqus visco elastic material
c All rights of reproduction or distribution in any form are reserved.
c By Irfan Habeeb CN (Technion - IIT), cnirfan@gmail.com
c ======================================================================
      subroutine vumat(
C Read only -
     1     nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2     stepTime, totalTime, dt, cmname, coordMp, charLength,
     3     props, density, strainInc, relSpinInc,
     4     tempOld, stretchOld, defgradOld, fieldOld,
     5     stressOld, stateOld, enerInternOld, enerInelasOld,
     6     tempNew, stretchNew, defgradNew, fieldNew,
C Write only -
     7     stressNew, stateNew, enerInternNew, enerInelasNew)
c
      include 'vaba_param.inc'
c
      dimension props(nprops), density(nblock), coordMp(nblock,*),
     1  charLength(nblock), strainInc(nblock,ndir+nshr),
     2  relSpinInc(nblock,nshr), tempOld(nblock),
     3  stretchOld(nblock,ndir+nshr),
     4  defgradOld(nblock,ndir+nshr+nshr),
     5  fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6  stateOld(nblock,nstatev), enerInternOld(nblock),
     7  enerInelasOld(nblock), tempNew(nblock),
     8  stretchNew(nblock,ndir+nshr),
     8  defgradNew(nblock,ndir+nshr+nshr),
     9  fieldNew(nblock,nfieldv),
     1  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     2  enerInternNew(nblock), enerInelasNew(nblock)
c
      character*80 cmname

      integer i, k1, k2

      real*8 E, nu, eta, lambda, mu,
     3  dde1(6,6), dde2(6,6), dde3(6,6)

c material properties
      E = props(1)            ! Young's modulus
      nu = props(2)           ! Poisson's ratio
      eta = props(3)          ! visco-elastic coef.
    
      do 10 i = 1, nblock

c lame's parameters
        lambda = E*nu/((1.0d0+nu)*(1.0d0-2.0d0*nu))
        mu = E/(2.0d0*(1.0d0+nu))

c stiffness matrix with viscous effects
        dde1 = 0.d0
        dde2 = 0.d0
        dde3 = 0.d0

        if (stepTime .gt. 0.d0) then
          do k1 = 1, ndir
            do k2 = 1, ndir
              dde1(k1, k2) = lambda*(1.0d0 + 3.0d0*eta/(E*dt))
              dde2(k1, k2) = lambda
            end do 
            dde1(k1, k1) = lambda*(1.0d0 + 3.0d0*eta/(E*dt)) +
     1        2.0d0*mu*(1.0d0 + eta/(mu*dt))
            dde2(k1, k1) = lambda + 2.0d0*mu
            dde3(k1, k1) = -1.0d0
          end do 
c shear stress
          do k1 = ndir+1, ndir+nshr
            dde1(k1, k1) = mu*(1.0d0 + eta/(mu*dt))
            dde2(k1, k1) = mu
            dde3(k1, k1) = -1.0d0
          end do 
        else 
          do k1 = 1, ndir
            do k2 = 1, ndir
              dde1(k1, k2) = lambda
              dde2(k1, k2) = lambda
            end do 
            dde1(k1, k1) = lambda + 2.0d0*mu
            dde2(k1, k1) = lambda + 2.0d0*mu
            dde3(k1, k1) = -1.0d0
          end do 
        do k1 = ndir+1, ndir+nshr
          dde1(k1, k1) = mu
          dde2(k1, k1) = mu
          dde3(k1, k1) = -1.0d0
        end do 
        end if 

c updating the state variable (old stress)
        stateNew(i, 1:ndir+nshr) = stressOld(i, 1:ndir+nshr)

c stress increment
        do k1 = 1, 6
          do k2 = 1, 6
            stateNew(i,k1) = stateNew(i,k1) + 
     1        dde1(k1,k2)*strainInc(i,k2) +
     2        dde2(k1,k2) * stateOld(i,k2+12) +
     3        dde3(k1,k2) * stressOld(i, k2)
          end do 
        end do

c updating stress from statevariables
        do j = 1, 6
          stressNew(i, j) = stateNew(i, j)
        end do

c Total strain stateNew(i,13-18) = total strain
        stateNew(i, 13:18) = stateOld(i, 13:18) + strainInc(i, 1:6)
        
  10  continue

      return
      end