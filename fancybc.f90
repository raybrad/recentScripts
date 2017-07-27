module fancyBC
  implicit none
  private

  !public subroutines
  public :: bndyc, coef, coef_gate
  public :: useFancyBC,initFancyBC,fancyRHS,fancyRHS_3term

  integer,parameter :: dp=8             !< double precision
  logical,save :: useFancyBC=.false.    !< use FancyBC or not
  integer,save :: nConstRegion          !< number of metallic region (const potential)
  real(kind=dp),allocatable,dimension(:,:),save :: constRegion    !< metallic regions: constRegion(1:6,nConstRegion)
  !
  interface fancyRHS
    module procedure fancyRHS
    module procedure fancyRHS_3term
  end interface fancyRHS

  contains
  subroutine initFancyBC(option,poisson)
    use parameters,only: BOHR
    use dataType,only: taskOption,PEData
    type(taskOption) :: option
    type(PEData) :: poisson
    !
    real(kind=dp) :: xa,xb,yc,yd,ze,zf,xcenter,ycenter,zcenter,LeadTop,LeadBot
    real(kind=dp) :: leadArea,gateArea
    !
    integer :: i,j,k
    !
    leadArea = option%LeadLength_y * option%LeadLength_z
    gateArea = option%GateLength_t * option%GateLength_l
    !
    if(leadArea > 0d0 .or. gateArea > 0d0 )then
      useFancyBC=.true.
    else
      useFancyBC=.false.
      return
    endif
    xa=poisson%xa
    xb=poisson%xb
    yc=poisson%yc
    yd=poisson%yd
    ze=poisson%ze
    zf=poisson%zf
    xcenter = (xa + xb) /2
    ycenter = (yc + yd) /2
    zcenter = (ze + zf) /2
    allocate(constRegion(6,6))

    if( leadArea > 0d0 )then
      !-------------------------------------
      ! define leads const potential regions
      !-------------------------------------
      if(option%LeadLength_x <= 0d0)option%LeadLength_x = poisson%dlx
      nConstRegion=2
      constRegion(1,1)  = xa
      constRegion(2,1)  = xa + option%LeadLength_x
      constRegion(3,1)  = ycenter - option%LeadLength_y/2
      constRegion(4,1)  = ycenter + option%LeadLength_y/2
      constRegion(5,1)  = zcenter - option%LeadLength_z/2
      constRegion(6,1)  = zcenter + option%LeadLength_z/2

      constRegion(1,2)  = xb - option%LeadLength_x
      constRegion(2,2)  = xb
      constRegion(3,2)  = ycenter - option%LeadLength_y/2
      constRegion(4,2)  = ycenter + option%LeadLength_y/2
      constRegion(5,2)  = zcenter - option%LeadLength_z/2
      constRegion(6,2)  = zcenter + option%LeadLength_z/2
    else
      nConstRegion=0
    endif

    if( gateArea > 0d0 )then
      !------------------------------------
      ! define gate const potential region
      !------------------------------------
      !no need!if(option%OxThick<0d0 .and. option%GateLength_h<0d0)then
      !no need!  write(6,*)'Both OxThick & GateLength_h < 0'
      !no need!  write(6,*)'GateLength_h is now set to ',max(poisson%dly,poisson%dlz) 
      !no need!endif
      nConstRegion=nConstRegion+1
      i=nConstRegion
      constRegion(1,i)  = xcenter - option%GateLength_l/2
      constRegion(2,i)  = xcenter + option%GateLength_l/2
      select case(option%gateDirect)
      case('Y-',"y-")  ! negative y direction
        constRegion(3,i)  = yc 
        if(option%OxThick >= 0d0)then
          !use OxThick to calculate GateLength_h 
          LeadTop = ycenter - option%LeadLength_y/2
          constRegion(4,i)  = LeadTop - option%OxThick
        else if(option%GateLength_h >= 0d0)then
          !use GateLength_h 
          constRegion(4,i)  = yc+option%GateLength_h
        else
          write(6,*)'ERROR in defining GateLength_h:',option%GateLength_h
        endif
        constRegion(5,i)  = zcenter - option%GateLength_t/2
        constRegion(6,i)  = zcenter + option%GateLength_t/2
      case('Y+',"y+")  ! positive y direction
        if(option%OxThick >= 0d0)then
          LeadBot = ycenter + option%LeadLength_y/2
          constRegion(3,i)  = LeadBot + option%OxThick
        else if(option%GateLength_h >= 0d0)then
          constRegion(3,i)  = yd - option%GateLength_h
        else
          write(6,*)'ERROR in defining GateLength_h:',option%GateLength_h
        endif
        constRegion(4,i)  = yd 
        constRegion(5,i)  = zcenter - option%GateLength_t/2
        constRegion(6,i)  = zcenter + option%GateLength_t/2
      case('Z-',"z-")  ! -ve z direction
        constRegion(5,i)  = ze
        if(option%OxThick >= 0d0)then
          LeadTop = zcenter - option%LeadLength_z/2
          constRegion(6,i)  = LeadTop - option%OxThick
        else if(option%GateLength_h >= 0d0)then
          constRegion(6,i)  = ze+option%GateLength_h
        else
          write(6,*)'ERROR in defining GateLength_h:',option%GateLength_h
        endif
        constRegion(3,i)  = ycenter - option%GateLength_t/2
        constRegion(4,i)  = ycenter + option%GateLength_t/2
      case('Z+',"z+")  ! +ve z direction
        if(option%OxThick >= 0d0)then
          LeadBot = zcenter + option%LeadLength_z/2
          constRegion(5,i)  = LeadBot + option%OxThick
        else if(option%GateLength_h >= 0d0)then
          constRegion(5,i)  = zf - option%GateLength_h
        else
          write(6,*)'ERROR in defining GateLength_h:',option%GateLength_h
        endif
        constRegion(6,i)  = zf 
        constRegion(3,i)  = ycenter - option%GateLength_t/2
        constRegion(4,i)  = ycenter + option%GateLength_t/2
      case default
        write(6,*)'ERROR in option%gateDirect:',option%gateDirect
        stop
      end select

    endif

    write(6,*)'Finish setting up FancyBC:'
    write(6,*)'  - nConstRegion:',nConstRegion
    do i=1,nConstRegion
      write(6,*)'  - Region ',i
      write(6,*)'    x: ',constRegion(1,i)*BOHR,constRegion(2,i)*BOHR 
      write(6,*)'    y: ',constRegion(3,i)*BOHR,constRegion(4,i)*BOHR
      write(6,*)'    z: ',constRegion(5,i)*BOHR,constRegion(6,i)*BOHR
    enddo

  end subroutine initFancyBC

  !> modify poisson%rhs to be the const potential values for <=3 terminal devices
  subroutine fancyRHS_3term(poisson,lead1V,lead2V,gateV)
    use dataType,only: PEData
    type(PEData) :: poisson
    real(kind=dp) :: lead1V,lead2V,gateV
    !
    logical :: debug=.true.
    real(kind=dp),allocatable,dimension(:) ::   constV
    integer :: nConst
    integer :: istat
    !
    nConst=nConstRegion
    allocate(constV(nConst),STAT=istat)
    if(istat/=0)then
      write(6,*)'allocate constV failed in fancyRHS_3term'
      stop
    endif

    if(nConst==1)then
      constV(1)=gateV
    else if(nConst==2)then
      constV(1)=lead1V
      constV(2)=lead2V
    else if(nConst==3)then
      constV(1)=lead1V
      constV(2)=lead2V
      constV(3)=gateV
    endif
    write(6,*)'fancyRHS: const V:',constV
    call fancyRHS(poisson,constV,nConst)
    deallocate(constV,STAT=istat)
  end subroutine

  !> modify poisson%rhs to be the const potential values in constRegions
  subroutine fancyRHS(poisson,constValue,nConst)
    use dataType,only: PEData
    type(PEData) :: poisson
    real(kind=dp) :: constValue(nConst)
    integer :: nConst
    !
    integer :: i,j,k,ix,iy,iz
    real(kind=dp) :: x_x,y_y,z_z
    real(kind=dp) :: xa,xb,yc,yd,ze,zf

    if(nConst /= nConstRegion)then
      write(6,*)'ERROR in fancyRHS:',nConst
      stop
    endif
    xa=poisson%xa
    xb=poisson%xb
    yc=poisson%yc
    yd=poisson%yd
    ze=poisson%ze
    zf=poisson%zf

    do ix=1,poisson%nx
      do iy=1,poisson%ny
        do iz=1,poisson%nz
          x_x=xa + (xb-xa)*(ix-1)/(poisson%nx-1)
          y_y=yc + (yd-yc)*(iy-1)/(poisson%ny-1)
          z_z=ze + (zf-ze)*(iz-1)/(poisson%nz-1)
          do i=1,nConstRegion
            if( x_x >= constRegion(1,i) .and. x_x <= constRegion(2,i) )then
              if( y_y >= constRegion(3,i) .and. y_y <= constRegion(4,i) )then
                if( z_z >= constRegion(5,i) .and. z_z <= constRegion(6,i) )then
                  poisson%rhs(ix,iy,iz)=constValue(i)
                endif
              endif
            endif
          enddo

        enddo
      enddo
    enddo


  end subroutine
  subroutine bndyc(kbdy,xory,yorz,alfa,gbdy)

    integer :: kbdy
    real(kind=dp) :: xory,yorz,alfa,gbdy

    alfa = 0.d0
    gbdy = 0.d0

  end subroutine bndyc

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  Subroutine coef(x,y,z,cxx,cyy,czz,cx,cy,cz,ce)

    real(kind=dp) :: x,y,z,cxx,cyy,czz,cx,cy,cz,ce

    cxx = 1.d0
    cyy = 1.d0
    czz = 1.d0
    cx = 0.d0          
    cy = 0.d0  
    cz = 0.d0
    ce = 0.d0    

    !if(any(localBC.gt.0)) call coef_local(x,y,z,cxx,cyy,czz,cx,cy,cz,ce)

  end subroutine coef
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine coef_gate(x,y,z,cxx,cyy,czz,cx,cy,cz,ce)

    real(kind=dp) :: cxx,cyy,czz,cx,cy,cz,ce
    real(kind=dp) :: x,y,z,x_x,y_y,z_z
    integer :: i_x,i_y,i_z,i

    ! Del^2 V = rhs
    cxx = 1.d0
    cyy = 1.d0
    czz = 1.d0
    cx = 0.d0          
    cy = 0.d0  
    cz = 0.d0
    ce = 0.d0   
    x_x=x 
    y_y=y 
    z_z=z 

    do i=1,nConstRegion
      if( x_x >= constRegion(1,i) .and. x_x <= constRegion(2,i) )then
        if( y_y >= constRegion(3,i) .and. y_y <= constRegion(4,i) )then
          if( z_z >= constRegion(5,i) .and. z_z <= constRegion(6,i) )then
            ! V = rhs
            cxx = 0.d0
            cyy = 0.d0
            czz = 0.d0
            cx = 0.d0          
            cy = 0.d0  
            cz = 0.d0
            ce = 1.d0   
          endif
        endif
      endif
    enddo

  end subroutine coef_gate
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



end module fancyBC
