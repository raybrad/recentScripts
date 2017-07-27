program read_dos
integer:: nlines,natoms,i,j,nepsi,natoms1
character*80 :: f(100),f0
real*4:: emin,emax,nume,Ef,a
real*4:: E1(5000),E2(5000),E4(5000),dos(5000),intdos(5000),dos_s(5000,5000),dos_pn1(5000,5000),dos_p0(5000,5000),dos_pp1(5000,5000),E(5000,5000),E3(5000,5000)
real*4:: dos_dn2(5000,5000),dos_dn1(5000,5000),dos_d0(5000,5000),dos_dp1(5000,5000),dos_dp2(5000,5000)
real*4:: dos_p(5000,5000),dos_d(5000,5000),delta
real*4:: intedos_ds1,intedos_dp1,intedos_dd1,intedos_ds2,intedos_dp2,intedos_dd2,intedos_all
real*4:: dos_s_total_2(5000),dos_p_total_2(5000),dos_d_total_2(5000),dos_s_total_1(5000),dos_p_total_1(5000),dos_d_total_1(5000)
character*6::int_to_char

open(100,file='DOSCAR')
open(200,file='dos_total.dat')
open(300,file='dos_total_1.dat')
open(400,file='dos_total_2.dat')

f0='_dos.dat'
write(*,*) "please input the Ef"
read(*,*)Ef
write(*,*) "please input the number of atoms in the POSCAR file"
read(*,*)natoms
write(*,*) "please input the number of the first atoms in the POSCAR file"
read(*,*)natoms1

do i=10,natoms+9
 f(i)=trim(int_to_char(i))//trim(f0)
 open(i,file=trim(f(i)))
enddo

read(100,*)
read(100,*)
read(100,*)
read(100,*)
read(100,*)  !跳过DOSCAR文件前面5行，从第6行开始起读
read(100,*) emax,emin,neps,Ef,a

delta=0.0
delta=(emax-emin)/1000  
intedos_all=0.0
do i=1,neps
   read(100,*)E1(i),dos(i),intdos(i)    !#第7行数据
   E2(i)=E1(i)-Ef
   intedos_all=intedos_all+dos(i)*delta
   write(200,10)E2(i),dos(i),intdos(i),intedos_all  !#这两个数值有细微差别
enddo

do i=10,natoms+9
   read(100,*)    !#空一行 即从1008行开始读取数据
   do j=1,neps
       read(100,*)E(i,j),dos_s(i,j),dos_p(i,j),dos_d(i,j)  
       E3(i,j)=E(i,j)-Ef  !#减去费米能级的修正    
       write(i,10)E3(i,j),dos_s(i,j),dos_p(i,j),dos_d(i,j)  
   enddo
 enddo

dos_s_total_1(j)=0.0
dos_p_total_1(j)=0.0
dos_d_total_1(j)=0.0
dos_s_total_2(j)=0.0
dos_p_total_2(j)=0.0
dos_d_total_2(j)=0.0

intedos_ds1=0.0
intedos_dp1=0.0
intedos_dd1=0.0
intedos_ds2=0.0
intedos_dp2=0.0
intedos_dd2=0.0


 do j=1,neps
    do i=10,natoms1+9
       dos_s_total_1(j)=dos_s_total_1(j)+dos_s(i,j)
       dos_p_total_1(j)=dos_p_total_1(j)+dos_p(i,j)
       dos_d_total_1(j)=dos_d_total_1(j)+dos_d(i,j)
       E4(j)=E3(i,j)
       intedos_ds1=intedos_ds1+dos_s_total_1(j)*delta
       intedos_dp1=intedos_dp1+dos_p_total_1(j)*delta
       intedos_dd1=intedos_dd1+dos_d_total_1(j)*delta
    enddo
       write(300,20)E4(j),dos_s_total_1(j),dos_p_total_1(j),dos_d_total_1(j),intedos_ds1,intedos_dp1,intedos_dd1 
       !#intedos_dp1 怎么会大于总的态密度积分值呢？ 费解！！！
 enddo

 do j=1,neps
     do i=natoms1+10,natoms+9
       dos_s_total_2(j)=dos_s_total_2(j)+dos_s(i,j)
       dos_p_total_2(j)=dos_p_total_2(j)+dos_p(i,j)
       dos_d_total_2(j)=dos_d_total_2(j)+dos_d(i,j)
       E4(j)=E3(i,j)
       intedos_ds2=intedos_ds2+dos_s_total_2(j)*delta
       intedos_dp2=intedos_dp2+dos_p_total_2(j)*delta
       intedos_dd2=intedos_dd2+dos_d_total_2(j)*delta
     enddo
       write(400,20)E4(j),dos_s_total_2(j),dos_p_total_2(j),dos_d_total_2(j),intedos_ds2,intedos_dp2,intedos_dd2
 enddo
10 format(1x,f8.4,1x,f8.4,1x,f8.4,1x,f8.4)
20 format(1x,f8.4,1x,f8.4,1x,f8.4,1x,f8.4,1x,f8.4,1x,f8.4,1x,f8.4)

end program read_dos
 
!-----------------------------------------------------------------------
FUNCTION int_to_char( int )
!-----------------------------------------------------------------------
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: int
     CHARACTER (LEN=6)   :: int_to_char
     IF ( int < 10 ) THEN
      WRITE( UNIT = int_to_char , FMT = "(I1)" ) int
     ELSE IF ( int < 100 ) THEN
      WRITE( UNIT = int_to_char , FMT = "(I2)" ) int
     ELSE IF ( int < 1000 ) THEN
      WRITE( UNIT = int_to_char , FMT = "(I3)" ) int
     ELSE IF ( int < 10000 ) THEN
      WRITE( UNIT = int_to_char , FMT = "(I4)" ) int
     ELSE
      WRITE( UNIT = int_to_char , FMT = "(I5)" ) int
     END IF
        RETURN
END FUNCTION int_to_char


