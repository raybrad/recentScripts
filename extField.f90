!call extField(sys%nAt,sys%nOrbs,sys%XYZ,sys%ind,tdMat%overl,tdMat%inv_overl,tdPara,tt,iExt,Vext)
subroutine extField(nAt,nOrbs,XYZ,ind,overl,inv_overl,tdPara,tt,iExt,Vext)
use dataType
use parameters, only: PI,HBARINV,BOHR
implicit none
integer,intent(in)::nAt,nOrbs,iExt
integer,intent(in)::ind(nAt)
real*8, intent(in)::XYZ(3,nAt),overl(nOrbs,nOrbs),inv_overl(nOrbs,nOrbs),tt
real*8, intent(out)::Vext(nOrbs,nOrbs)
type(tdParameters), intent(in) :: tdPara

integer::ii,ia,ja,li,lj,nDirect,ipot
real*8 ::dipmax(nOrbs,nOrbs,3),field(3),strength(3),freq(3),ref_pos(3)

 do ii=1,3
  ref_pos(ii)=0.d0
      do ia=tdPara%dcensta,tdPara%dceneda
         ref_pos(ii)= ref_pos(ii)+XYZ(ii,ia)
      enddo
         ref_pos(ii)= ref_pos(ii)/(tdPara%dceneda-tdPara%dcensta+1)
	 write(6,'(A20,X,F9.4,X,F9.4 )') 'Dipmax reference Center(A)',ref_pos(ii)*BOHR
 enddo


do ii =1,3
    do ia = 1, nAt
        do ja = 1, nAt
  	       do li = ind(ia)+1, ind(ia+1)
	       do lj = ind(ja)+1, ind(ja+1)
	           dipmax(li,lj,ii) = overl(li,lj)*(0.5d0*(XYZ(ii,ia)+XYZ(ii,ja))-ref_pos(ii))
!	           dipmax(li,lj,ii) = overl(li,lj)*(0.5d0*(XYZ(ii,ia)+XYZ(ii,ja)))
	       enddo
	       enddo
    	 enddo
     enddo
enddo

do ii=1,3
call dtrmm('l', 'u', 't', 'n', nOrbs, nOrbs, 1.d0, inv_overl, nOrbs, dipmax(1,1,ii), nOrbs)
call dtrmm('r', 'u', 'n', 'n', nOrbs, nOrbs, 1.d0, inv_overl, nOrbs, dipmax(1,1,ii), nOrbs)
enddo

if (tdPara%ipot==-1) then
     if (abs(tt-0.d0)<=1.d-9) then
         field(1:3)=tdPara%strength(1:3)
     else
         field(1:3)=0.d0
     endif
elseif (tdPara%ipot==0) then
     field(1:3) = tdPara%strength(1:3)*dsin(tdPara%freq(1:3)*tt*HBARINV) 
elseif (tdPara%ipot == 1) then
     if (tt<=7*tdPara%tlas) then         !here freq is Energy eV (omega=E/hbar=2pi*f) -> omega*tt
	  field(1:3)=tdPara%strength(1:3)*dsin(tdPara%freq(1:3)*tt*HBARINV)*exp(-((tt-3*tdPara%tlas)/tdPara%tlas)**2.0)
     else
 	  field(1:3)=0.d0
     endif
endif
if (iExt) then
open(16,file='extField.dat',access='append')
write(16,'(4(e12.5,2x))') tt,field(1),field(2),field(3)
close(16)
endif
!do ii=1,3
!Vext(1:nOrbs,1:nOrbs)=Vext(1:nOrbs,1:nOrbs)+dipmax(1:nOrbs,1:nOrbs,ii)*field(ii)
!enddo

do ii=1,3
    do ia = tdPara%sta,tdPara%eda
        do ja = tdPara%sta,tdPara%eda
  	       do li = ind(ia)+1, ind(ia+1)
	       do lj = ind(ja)+1, ind(ja+1)
		Vext(li,lj)=Vext(li,lj)+dipmax(li,lj,ii)*field(ii)
	       enddo
	       enddo
	 enddo
     enddo
enddo
!
end subroutine extField
