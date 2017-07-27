program TwoDconverter
     implicit none
     integer::nline,i,j,Nenergy,Nwavenum,istat,system
     real*8::x,y,z,denergy,energy,dtmp,wavenum,ramanint,ramanint2,wavemin,wavemax,energymin,energymax
     character(len=100) :: command,charEnergy,tempchar,job_typ,connect_typ,ext_typ,deform_typ,ramanSpecShape
     logical::calDeform,calIR,calRaman

     integer::line_input,nLeadAtoms,nAtoms
     real*8::omega,dep,crit,lambda_raman
     integer,allocatable::atn(:)
     real*8,allocatable ::xyz(:,:)
     logical::doRaman
      
     omega=0.d0
     line_input=10001
     dep=0.1d0
     crit=1.0d-10	  
     !in in.2D copy from in.coord and add some other namelist
     job_typ='mpiso'
     connect_typ='free'
     ext_typ='nogamma' 
     deform_typ='td'
     ramanSpecShape='lorentz'
     calDeform=.false.
     calIR=.false. 
     calRaman=.true.
     namelist /task/ job_typ,connect_typ,ext_typ,deform_typ,ramanSpecShape,calDeform,calIR,calRaman
     rewind 5
     read(5,task,end=90)
     90 continue
     write(6,*) 'connect_typ',connect_typ

     namelist /fft/omega,line_input,dep,crit
     rewind 5
     read(5,fft,end=100)
     100 continue
	
     lambda_raman=514.5d0 !m
     doRaman=.false.
     energymin=0.5
     energymax=3.0
     denergy=0.01d0
     rewind 5
     namelist /ramanpara/lambda_raman,doRaman,energymin,energymax,denergy
     read(5,ramanpara,end=101)
     101 continue
     lambda_raman=1239.84187d0/omega
     nLeadAtoms=0
     nAtoms=0
     rewind 5
     namelist /coordinate/ nLeadAtoms,nAtoms
     read(5,coordinate,end=102)
     102 continue
    
     allocate(atn(nAtoms),xyz(3,nLeadAtoms+nAtoms))
     readcoor: do i=1,nLeadAtoms+nAtoms
    		 read(5,*) j,x,y,z 
    		 atn(i)=j
    		 xyz(1,i)=x
    		 xyz(2,i)=y
    		 xyz(3,i)=z
     enddo readcoor
     
    
!!!!!!!!!
    Nenergy=int((energymax-energymin)/denergy)+1
!!!!!!!!!!
if (doRaman==.true.) then


     do i=1,Nenergy
    	 energy=energymin+dble(i-1)*denergy
         lambda_raman=1239.84187d0/energy
         open(16,file='in.raman')
         if(connect_typ(1:4)=='free') then
          write(16,'(A)') "$task job_typ='mpiso' connect_typ='free' ext_typ='nogamma' deform_typ='td' calDeform=.false. calIR=.false. calRaman=.true. ramanSpecShape='lorentz' $end"
         elseif(connect_typ(1:5)=='bound') then
          write(16,'(A)') "$task job_typ='mpiso' connect_typ='bound' ext_typ='nogamma' deform_typ='td' calDeform=.false. calIR=.false. calRaman=.true. ramanSpecShape='lorentz' $end"
         endif
          write(16,'(A,F9.4,A)')"$ramanPara distortion=0.01d0 atomSta=1 atomEda=11 lambdaMeasure=",lambda_raman," temperature=298.15d0 theta=90.0 $end"
          write(16,'(A,I5,A,F5.2,A,F5.2,A,E12.5,A)')"$fft line_input=",line_input," dep=",dep," omega=",energy," crit=",crit," $end"
          write(16,'(A,I3,A,I3,A)')"$coordinate  nLeadAtoms=",nLeadAtoms," nAtoms=",nAtoms," $end"
	  do j=1,nLeadAtoms+nAtoms
	     write(16,'(I2,3(F13.6))') Atn(j),xyz(1:3,j)
	  enddo
         close(16)
	 istat=system('~/program/dftbraman/raman<in.raman>out.raman')
	 if (istat/=0 ) then
		 write(6,*) 'istat/=0'
		 stop
         else
		 write(6,*) 'ok in energy',energy
	 endif
	 
	 write(tempchar,'(F9.4)') energy
	 charEnergy='2D/raman_lorentz_spectrum.'//trim(adjustl(tempchar))//'eV'
	 command="mv raman_lorentz_spectrum.dat "//trim(charEnergy)
	 write(6,*) trim(command)
	 istat=system(trim(command))
	 if (istat/=0 ) then
		 write(6,*) 'istat/=0'
		 stop
         else
		 write(6,*) 'ok in mv'
	 endif
     enddo
endif
!!!!!!!!!!
     Nwavenum=8000
     wavemin=300
     wavemax=3500
     open(18,file='2Dlorentz_spectrum.dat')
     do i=1,Nenergy
	 energy=energymin+dble(i-1)*denergy
	 write(tempchar,'(F9.4)') energy
	 charEnergy='2D/raman_lorentz_spectrum.'//trim(adjustl(tempchar))//'eV'
	 write(6,*) trim(charEnergy)
	 open(17,file=trim(charEnergy))
	 do j=1,Nwavenum
	 read(17,*) wavenum,ramanint,ramanint2
	 if (wavenum>=wavemin .and. wavenum<=wavemax) then
	 write(18,'(4(e11.4,3x))') energy,wavenum,ramanint,ramanint2
         endif
	 enddo
	 write(18,*)
	 close(17)
     enddo

     close(18)

end program TwoDconverter
