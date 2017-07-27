program read_band
!*      This program is to get the suitable data for ploting bands of energy
       real kx(1000),ky(1000),kz(1000)
       real y(1000,100)
       character ch1,ch2,ch3,ch4
       open(1,file='EIGENVAL')
       open(2,file='bands.dat')
!*      This file intend to generate the high symmetry point
       open(3,file='spoint.dat')
       open(4,file='kpoints.dat')
       read(1,*)
       read(1,*)
       read(1,*)
       read(1,*)
       read(1,*)
       read(1,*)ch1,nks,nband
       write(*,*) "please input the fermi energy"
       read(*,*)fermi
!c       read(1,*)emin,emax
       do i=1,nks
         read(1,*)
         read(1,*)
          do j=1,nband
           read(1,*)ch3,y(i,j)
          enddo
       enddo
       read(4,*)
       do i=1,nks
          read(4,*)kx(i),ky(i),kz(i)
       enddo
!*      This part intend to chose the energy range
!c       do i=1,nks
!c         do j=1,nband
!c           if (y(i,j).lt.emin) then
!c             y(i,j)=emin
!c           else if (y(i,j).gt.emax) then
!c             y(i,j)=emax
!c           endif
!c         enddo
!c       enddo
       x=0.0
       x_tot=0.0
       do i=1,nks
         if(i.gt.1) then
           x=sqrt((kx(i)-kx(i-1))**2+(ky(i)-ky(i-1))**2+(kz(i)-kz(i-1))**2)
         endif
       x_tot=x_tot+x
         if(mod((i-1),20).eq.0 .or. i.eq.nks) then
           write(3,'(100F10.4)')x_tot
         endif
         write(2,'(100F10.4)')x_tot,(y(i,j)-fermi,j=1,nband)
       enddo
       stop
 end program read_band
