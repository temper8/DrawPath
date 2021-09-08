      subroutine configlhcd
      implicit double precision (a-h,o-z)
      integer icon,item,itemi,iazef
      double precision con,tem,abtor,rm,x0,z0,rh1,rh,rha,amy
     &,delta,ell,gamma,cdl,cly,cgm,cmy,coeffs,zero,drhodr
      include 'for/parameter.inc'
      include 'for/const.inc'
      include 'for/status.inc'
      external polin,polin1
      dimension con(NRD),tem(NRD),temi(NRD),azef(NRD),afld(NRD)
      dimension rh(NRD),delta(NRD),ell(NRD),gamma(NRD)
      dimension rha(NRD),drhodr(NRD),amy(NRD)
      dimension cdl(10),cly(10),cgm(10),cmy(10),coeffs(10)
      dimension adelta(NRD),aell(NRD),agamma(NRD)
      dimension ddelta(NRD),dell(NRD),dgamma(NRD),fmy(NRD)
      parameter(zero=0.d0,ipsy=5)
cc*********************************************************************
cc   ipsy = number of polinomial decomposition coefficients
cc           used for interpolation of Zakharov's moments.
cc*********************************************************************
cc*********************************************************************
cc    Co-ordinates used in ray-tracing:
cc         (x-x0)/rm=r*cos(teta)-delta-gamma*sin^2(teta)
cc         (z-z0)/rm=ell*r*sin(teta)
cc    Definitions:
cc    (x0,z0) - magnetic axis position, centimeters
cc    rm      - minor radius in mid-plane, cenrimeters
cc    r(rho_ASTRA),delta(r),gamma(r),ell(r) - dimensionless functions
cc    rho_ASTRA=sqrt(Phi_tor/GP/BTOR)
cc    Interval for r:  0.<= r <=1.
cc*********************************************************************

!!!      UPDWN=zero
      ipsy1=ipsy-1
      inpt=NA1          ! ASTRA radial grid number
      do i=1,inpt
       rh(i)=AMETR(i)/ABC
       rha(i)=RHO(i)/ABC  !/ABC instead of /ROC is not a mistake!
       delta(i)=(SHIF(1)-SHIF(i))/ABC  !FRTC Shafr. shift. defin.
       ell(i)=ELON(i)
       gamma(i)=rh(i)*TRIA(i)
!!!       afld(i)=ULON(i)/RTOR/GP2
       afld(i)=UPL(i)/RTOR/GP2
       con(i)=NE(i)
       if(con(i).gt.zero) then
        icon=i
       else
        con(i)=con(icon)
       end if
       tem(i)=dble(TE(i))
       if(tem(i).gt.zero) then
        item=i
       else
        tem(i)=tem(item)
       end if
       temi(i)=dble(TI(i))
       if(temi(i).gt.zero) then
        itemi=i
       else
        temi(i)=temi(itemi)
       end if
       azef(i)=dble(ZEF(i))
       if(azef(i).gt.zero) then
        iazef=i
       else
        azef(i)=azef(iazef)
       end if
      end do

      rh1=rh(1)          !saving the first ASTRA radial grid element
      rh(1)=zero         !shifting the first element to zero
      rha(1)=zero        !shifting the first element to zero
      delta(1)=zero      !putting delta(rh=0.)=0.
      gamma(1)=zero      !putting gamma(rh=0.)=0.

      abtor=1.d4*BTOR*RTOR/(RTOR+SHIF(1)) !B_tor_(magnetic axis), Gauss
      rm=1.d2*ABC                       !minor radius in mid-plane, cm
      x0=1.d2*(RTOR+SHIF(1))     !x-coordinate of the magnetic axis, cm
      z0=1.d2*UPDWN              !z-coordinate of the magnetic axis, cm


cccc   shift as a function of "minor radius":
       call approx(rh,delta,inpt,polin1,ipsy1,coeffs)
       cdl(1)=zero
       do k=2,ipsy
        cdl(k)=coeffs(k-1)
       end do

cccc   triangularity as a function of "minor radius":
       call approx(rh,gamma,inpt,polin1,ipsy1,coeffs)
       cgm(1)=zero
       do k=2,ipsy
        cgm(k)=coeffs(k-1)
       end do

cccc   ellipticity as a function of "minor radius":
       call approx(rh,ell,inpt,polin,ipsy,cly)

cccc   "poloidal magnetic field":
       call diff(rh,rha,inpt,drhodr)
       do i=2,inpt
        amy(i)=1.d4*BTOR*MU(i)*rha(i)*drhodr(i)
       end do
       amy(1)=zero
!! amy=(btor/q)*rho*(drho/dr) is a function of "minor radius" r=rh(i).
!! Poloidal magnetic field: B_pol=amy(r)*sqrt(g22/g), where g is
!! determinant of 3D metric tensor and g22 is the (22) element of
!! the tensor, normalized on ABC^4 and ABC^2, correspondingly.
!!
!!  Polinomial approximation of the amy(r):
       inpt2=inpt-3
       call approx(rh,amy,inpt2,polin1,ipsy1,coeffs)
       cmy(1)=zero
       do k=2,ipsy
        cmy(k)=coeffs(k-1)
       end do
cc************************************************************************
c   Saving data necessary for ray-tracing
c   This file can be used for running ray-tracing without ASTRA
c   for a given plasma slice
cc************************************************************************
       open(97, file='lhcd/config/inputlhcd.dat')
       write(97,5) TIME
       write(97,5) QLH     !input LH power, MW
       write(97,*) inpt    !number of ASTRA grid points
       write(97,*) ipsy    !decomposition coeffs. number
       write(97,5) rm      !minor radius in mid-plane, cm
       write(97,5) x0      !magnetic axis x-coordinate, cm
       write(97,5) z0      !magnetic axis z-coordinate, cm
       write(97,5) rh1     !smallest ASTRA  rh=AMETR(1)/ABC
       write(97,5) abtor   !B_tor at magnetic axis, Gauss
       do i=1,ipsy
        write(97,5) cdl(i),cly(i),cgm(i),cmy(i)
       end do
       do i=1,inpt
        write(97,5) rh(i),con(i),tem(i),temi(i),azef(i),afld(i)
       end do
       close(97)
5     format(7(e12.5,1x))
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine diff(x,y,n,dy)
      implicit real*8 (a-h,o-z)
      dimension y(*),x(*),dy(*)
      dy(1)=(y(2)-y(1))/(x(2)-x(1))
      do k=2,n-1
       dy(k)=(y(k+1)-y(k-1))/(x(k+1)-x(k-1))
      end do
      dy(n)=(y(n)-y(n-1))/(x(n)-x(n-1))
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function polin(k,x)
      implicit real*8 (a-h,o-z)
      polin=1.d0
      if(k.gt.1) then
       polin=x**(k-1)
      end if
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function polin1(k,x)
      implicit real*8 (a-h,o-z)
      polin1=x**k
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine approx(x,y,n,f,m,b)
c
c     y(i)=y(x(i))  the data to be approximated
c     n             number of points in the input data
c     m             number of coefficients of decomposition
c                   over base functions "f(k,x)" :
c                          y(x)=SUM_1^m [b(k)*f(k,x)]
c     b(i)          found decomposition coefficients
c
      implicit real*8 (a-h,o-z)
      parameter(zero=0.d0, np=20)
      dimension y(n),x(n)
      dimension a(np,np),indx(np),b(np)

       if(m.gt.np) then
         write(*,*)'index error subroutine "approx"'
         stop
       end if

         do j=1,m
          do k=1,j
           a(k,j)=zero
            do i=1,n
             a(k,j)=a(k,j)+f(j,x(i))*f(k,x(i))
            end do
          end do
         end do
                 do k=2,m
                  do j=1,k-1
                    a(k,j)=a(j,k)
                  end do
                 end do

          do k=1,m
           b(k)=zero
            do i=1,n
             b(k)=b(k)+y(i)*f(k,x(i))
            end do
          end do

        call ludcmp(a,m,np,indx,d)
        call lubksb(a,m,np,indx,b)

      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine ludcmp(a,n,np,indx,d)
      implicit real*8 (a-h,o-z)
      parameter (nmax=100, tiny=1.d-20, zero=0.d0)
      dimension a(np,np),indx(n),vv(nmax)
      d=1.d0
      do 12 i=1,n
        aamax=zero
        do 11 j=1,n
          if (dabs(a(i,j)).gt.aamax) aamax=dabs(a(i,j))
11      continue
        if (aamax.eq.zero) pause 'singular matrix.'
        vv(i)=1.d0/aamax
12    continue
      do 19 j=1,n
        if (j.gt.1) then
          do 14 i=1,j-1
            sum=a(i,j)
            if (i.gt.1)then
              do 13 k=1,i-1
                sum=sum-a(i,k)*a(k,j)
13            continue
              a(i,j)=sum
            endif
14        continue
        endif
        aamax=zero
        do 16 i=j,n
          sum=a(i,j)
          if (j.gt.1)then
            do 15 k=1,j-1
              sum=sum-a(i,k)*a(k,j)
15          continue
            a(i,j)=sum
          endif
          dum=vv(i)*dabs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(j.ne.n) then
          if(a(j,j).eq.zero) a(j,j)=tiny
          dum=1.d0/a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      if(a(n,n).eq.zero) a(n,n)=tiny
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine lubksb(a,n,np,indx,b)
      implicit real*8 (a-h,o-z)
      parameter(zero=0.d0)
      dimension a(np,np),indx(n),b(n)
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.zero) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        if(i.lt.n)then
          do 13 j=i+1,n
            sum=sum-a(i,j)*b(j)
13        continue
        endif
        b(i)=sum/a(i,i)
14    continue
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine lock(xa,n,x,klo,khi,ierr)
      implicit real*8 (a-h,o-z)
      dimension xa(*)
      parameter(tiny=1.d-14)
      klo=0
      khi=0
      dx1=x-xa(1)
      dx2=x-xa(n)
      if(dx1*dx2/dabs(xa(n)-xa(1)).gt.tiny) then
       ierr=1
       return
      end if
      ierr=0
      klo=1
      khi=n
      do while(khi-klo.gt.1)
       k=(khi+klo)/2
       if(xa(k).gt.x)then
         khi=k
       else
         klo=k
       endif
      end do
      if(khi.eq.klo) ierr=1
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine linf(x,y,t,fout,klo,khi)
      implicit real*8 (a-h,o-z)
      dimension x(*),y(*)
       dout=(y(khi)-y(klo))/(x(khi)-x(klo))
       fout=y(klo)+dout*(t-x(klo))
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
