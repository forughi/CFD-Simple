!******************************************
!************ CFD Course Project **********
!************** A. F. Forughi *************
!**************** (6/2011) ****************
!***** Sharif University of Technology ****
!***     FVM Method: Staggered Grid     ***
!***          SIMPLE algorithm          ***
!***     Solver: BiCGSTAB (Sparse)      ***
!******************************************


implicit none
allocatable ps(:,:),pp(:,:),us(:,:),vs(:,:),apu(:,:),apv(:,:)
allocatable umatv(:),umati(:),umatj(:),umrhs(:),x(:),vmatv(:),vmati(:),vmatj(:),vmrhs(:)
allocatable pmatv(:),pmati(:),pmatj(:),pmrhs(:),pxx(:),pyy(:)
allocatable uwn(:),uws(:),uww(:),uwe(:),pbcn(:),pbcs(:),pbcw(:),pbce(:)
real(8) lx,ly,dx,dy,ps,pp,us,vs,rho,mu,alphau,alphav,alphap,umatv,umrhs,x,bires,vmatv,vmrhs,pmatv,pmrhs,b
real(8) fe,fw,fn,fs,de,dw,ds,dn,ae,aw,an,as,ap,bound,apu,apv
real(8) uwn,uws,uww,uwe,spw
real(8) idx,jdy,val,mi,ma,contres,creslimit
integer i,j,ti,ni,nj,umati,umatj,unnzero,vmati,vmatj,vnnzero,pmati,pmatj,pnnzero,un,vn,pn,k,bimaxit,pxx,pyy,simpit
integer pbcn,pbcs,pbcw,pbce
character side

open(1,file="result.plt")
open(2,file="residual.plt")
open(3,file="profiles.plt")

print*," CFD Project, A. F. Forughi (6/2011)"
print*," "
print*,"Please enter the dimensions:"
print*,"Lx:"
read*,lx
print*,"Ly:"
read*,ly
print*,"dx=dy:"
read*,dx
!lx=1.0d0
!ly=1.0d0
!dx=0.025d0

dy=dx
rho=1.0d0
mu=1.00d0
print*,"-----------------------------------------------------------"

ni=int(lx/dx) !no. of nodes on x
nj=int(ly/dy) !no. of nodes on y

allocate (ps(1:ni+2,1:nj+2),pp(1:ni+2,1:nj+2),us(1:ni+3,1:nj+3),vs(1:ni+3,1:nj+3))

ps=0.0d0;pp=0.0d0;us=0.0d0;vs=0.0d0
alphau=0.7d0 ; alphav=0.7d0 ; alphap=0.3d0

unnzero=5*(ni-3)*(nj-2)+4*(ni-3)*2+4*(nj-2)*2+3*4!none zero
vnnzero=5*(ni-2)*(nj-3)+4*(ni-2)*2+4*(nj-3)*2+3*4
pnnzero=5*(ni-2)*(nj-2)+4*(ni-2)*2+4*(nj-2)*2+3*4
un=(ni-1)*(nj) ; vn=(ni)*(nj-1) ;pn=ni*nj !n*n
allocate (umatv(1:unnzero),umati(1:unnzero),umatj(1:unnzero),umrhs(1:un),x(1:pn))
allocate (vmatv(1:vnnzero),vmati(1:vnnzero),vmatj(1:vnnzero),vmrhs(1:vn),apu(3:ni+1,2:nj+1),apv(2:ni+1,3:nj+1))
allocate (pmatv(1:pnnzero),pmati(1:pnnzero),pmatj(1:pnnzero),pmrhs(1:pn),pxx(1:pn),pyy(1:pn))
allocate (uwn(1:ni+3),uws(1:ni+3),uww(1:nj+3),uwe(1:nj+3),pbcn(2:ni+1),pbcs(2:ni+1),pbcw(2:nj+1),pbce(2:nj+1))
apu=0.0d0 ; apv=0.0d0 ; pmrhs=0.0d0
uwn=0.0d0 ;uws=0.0d0 ;uww=0.0d0 ;uwe=0.0d0
pbcn=0 ;pbcs=0 ;pbcw=0 ;pbce=0


!Set ICs:
!P IC:
ps=0.0d0
!U IC:
us(3:ni+1,2:nj+1)=-1.0e-20
!V IC:
vs(2:ni+1,3:nj+1)=-1.0e-20


!Set BCs:
!Moved wall:
print*,"                          Default BC is static Wall"
print*," "
print*,"Moving Wall Boundary Condition(s):"
do
    print*,"Enter Side letter: n:North, s:South, w:West, e:East or another to skip this BC"
	read*,side
	if ((side=="n").or.(side=="s").or.(side=="e").or.(side=="w")) then
		print*,"Enter Velocity:"
		read*,val
		if     (side=="n") then
			print*,"Enter X min:"
			read*,mi
			print*,"Enter X max:"
			read*,ma
			uwn(int(mi/dx)+3:int(ma/dx)+1)=val
		elseif (side=="s") then
			print*,"Enter X min:"
			read*,mi
			print*,"Enter X max:"
			read*,ma
			uws(int(mi/dx)+3:int(ma/dx)+1)=val
		elseif (side=="e") then
			print*,"Enter Y min:"
			read*,mi
			print*,"Enter Y max:"
			read*,ma
			uwe(int(mi/dy)+3:int(ma/dy)+1)=val
		elseif (side=="w") then
			print*,"Enter Y min:"
			read*,mi
			print*,"Enter Y max:"
			read*,ma
			uww(int(mi/dy)+3:int(ma/dy)+1)=val
		else
			exit
		endif
	else
		exit
	endif
enddo
print*,"-----------------------------------------------------------"

!uwn(:)=0.0d0 ;uws(:)=0.0d0 ;uww(:)=0.0d0 ;uwe(:)=0.0d0 !Wall velocity on some where

uwn(1)=uwn(3);uwn(2)=uwn(3);uwn(ni+3)=uwn(ni+1);uwn(ni+2)=uwn(ni+1) !for writing data on wall
uws(1)=uws(3);uws(2)=uws(3);uws(ni+3)=uws(ni+1);uws(ni+2)=uws(ni+1)
uww(1)=uww(3);uww(2)=uww(3);uww(nj+3)=uww(nj+1);uww(nj+2)=uww(nj+1)
uwe(1)=uwe(3);uwe(2)=uwe(3);uwe(nj+3)=uwe(nj+1);uwe(nj+2)=uwe(nj+1)

!Velocity Inlet:
print*,"Velocity Inlet Boundary Condition(s):"
do
	print*,"Enter Side letter: n:North, s:South, w:West, e:East or another to skip this BC"
	read*,side
	if ((side=="n").or.(side=="s").or.(side=="e").or.(side=="w")) then
		print*,"Enter Velocity:"
		read*,val
		if     (side=="n") then
			print*,"Enter X min:"
			read*,mi
			print*,"Enter X max:"
			read*,ma
			vs(int(mi/dx)+2:int(ma/dx)+1,nj+2)=-val
		elseif (side=="s") then
			print*,"Enter X min:"
			read*,mi
			print*,"Enter X max:"
			read*,ma
			vs(int(mi/dx)+2:int(ma/dx)+1,2)=val
		elseif (side=="e") then
			print*,"Enter Y min:"
			read*,mi
			print*,"Enter Y max:"
			read*,ma
			us(ni+2,int(mi/dy)+2:int(ma/dy)+1)=-val
		elseif (side=="w") then
			print*,"Enter Y min:"
			read*,mi
			print*,"Enter Y max:"
			read*,ma
			us(2,int(mi/dy)+2:int(ma/dy)+1)=val
		else
			exit
		endif
	else
		exit
	endif
enddo
print*,"-----------------------------------------------------------"

!us(2,:)=+0.0d0    !W
!us(ni+2,:)=-0.0d0 !E
!vs(:,2)=+0.0d0    !S
!vs(:,nj+2)=-0.0d0 !N


!Pressure:
print*,"Constant Peressure Condition(s):"
do
	print*,"Enter Side letter: n:North, s:South, w:West, e:East or another to skip this BC"
	read*,side
	if ((side=="n").or.(side=="s").or.(side=="e").or.(side=="w")) then
		print*,"Enter Pressure:"
		read*,val
		if     (side=="n") then
			print*,"Enter X min:"
			read*,mi
			print*,"Enter X max:"
			read*,ma
			pbcn(int(mi/dx)+1:int(ma/dx)+1)=1
			ps(int(mi/dx)+1:int(ma/dx)+1,nj+1)=val
		elseif (side=="s") then
			print*,"Enter X min:"
			read*,mi
			print*,"Enter X max:"
			read*,ma
			pbcs(int(mi/dx)+1:int(ma/dx)+1)=1
			ps(int(mi/dx)+1:int(ma/dx)+1,2)=val
		elseif (side=="e") then
			print*,"Enter Y min:"
			read*,mi
			print*,"Enter Y max:"
			read*,ma
			pbce(int(mi/dy)+1:int(ma/dy)+1)=1
			ps(ni+1,int(mi/dy)+1:int(ma/dy)+1)=val
		elseif (side=="w") then
			print*,"Enter Y min:"
			read*,mi
			print*,"Enter Y max:"
			read*,ma
			pbcw(int(mi/dy)+1:int(ma/dy)+1)=1
			ps(2,int(mi/dy)+1:int(ma/dy)+1)=val
		else
			exit
		endif
	else
		exit
	endif
enddo
print*,"-----------------------------------------------------------"

!Active Peressure BC:
!pbcn(:)=0 ;pbcs(:)=0 ;pbce(:)=0 ;pbcw(:)=0 ! 1=Active  , 0=Deactive
!Set Value of Pressure BC
!ps(2,:)=+0.0d0    !W
!ps(ni+1,:)=+0.0d0 !E
!ps(:,2)=+0.0d0    !S
!ps(:,nj+1)=+0.0d0 !N

print*," "
print*,"Enter Convergance Conditions:"
print*,"Enter No. of Maximum iteration="
read*,simpit
print*,"Enter minimum Continuity Residual threshold="
read*,creslimit


bires=1.0e-20     !BiCGSTAB Minimum Residual
bimaxit=300       !BiCGSTAB Maximum Iterations
!creslimit=1.0e-10 !Continuity Minimum Residual
!simpit=300        !Continuity Maximum Iterations

do ti=1,simpit !SIMPLE Algorithm Loop

	!U momentum:
	k=0 !counter
	do j=2,nj+1
		do i=3,ni+1
			fe=rho*dy*(us(i,j)+us(i+1,j))/2
			fw=rho*dy*(us(i,j)+us(i-1,j))/2
			fs=rho*dx*(vs(i,j)+vs(i-1,j))/2
			fn=rho*dx*(vs(i,j+1)+vs(i-1,j+1))/2
			de=mu*dy/dx ; dw=mu*dy/dx ; ds=mu*dx/dy ; dn=mu*dx/dy
			!Hybrid:
			ae=max(-fe,(de-fe/2.0),0.0d0)!
			as=max(+fs,(ds+fs/2.0),0.0d0)!
			aw=max(+fw,(dw+fw/2.0),0.0d0)!
			an=max(-fn,(dn-fn/2.0),0.0d0)!
			ap=ae+as+aw+an+(fe-fw+fn-fs)

			!U Matrixer:
			bound=0.0

			if (i/=3) then !W
				k=k+1
				umatv(k)=-aw
				umati(k)=(j-2)*(ni-1)+(i-2-1)
				umatj(k)=(j-2)*(ni-1)+(i-2)
			else
				bound=bound+aw*us(i-1,j) !for the velocity inlet
			endif

			if (i/=ni+1) then !E
				k=k+1
				umatv(k)=-ae
				umati(k)=(j-2)*(ni-1)+(i-2+1)
				umatj(k)=(j-2)*(ni-1)+(i-2)
			else
				bound=bound+ae*us(i+1,j) !for the velocity inlet
			endif

			if (j/=nj+1) then !N
				k=k+1
				umatv(k)=-an
				umati(k)=(j-2+1)*(ni-1)+(i-2)
				umatj(k)=(j-2)*(ni-1)+(i-2)
			else
				ap=ap+mu*dx/(0.5*dy) !moving wall !*((us(i,j)-uw)/us(i,j))
				bound=mu*(uwn(i))*dx/(0.5*dy)
			endif

			if (j/=2) then !S
				k=k+1
				umatv(k)=-as
				umati(k)=(j-2-1)*(ni-1)+(i-2)
				umatj(k)=(j-2)*(ni-1)+(i-2)
			else
				ap=ap+mu*dx/(0.5*dy) !static wall
				bound=mu*(uws(i))*dx/(0.5*dy)
			endif

			k=k+1
			umatv(k)=ap/alphau!-sp
			umati(k)=(j-2)*(ni-1)+(i-2)
			umatj(k)=(j-2)*(ni-1)+(i-2)
			apu(i,j)=ap !for the pressure correction

			umrhs((j-2)*(ni-1)+(i-2))=(ps(i-1,j)-ps(i,j))*dy+((1-alphau)*ap/alphau)*us(i,j)+bound!+sc

		enddo
	enddo ! The End of U momentum



	!V momentum:

	k=0 !counter
	do j=3,nj+1
		do i=2,ni+1
			fe=rho*dy*(us(i+1,j)+us(i+1,j-1))/2
			fw=rho*dy*(us(i,j)+us(i,j-1))/2
			fs=rho*dx*(vs(i,j)+vs(i,j-1))/2
			fn=rho*dx*(vs(i,j)+vs(i,j+1))/2
			de=mu*dy/dx ; dw=mu*dy/dx ; ds=mu*dx/dy ; dn=mu*dx/dy
			!Hybrid:
			ae=max(-fe,(de-fe/2.0),0.0d0)!
			as=max(+fs,(ds+fs/2.0),0.0d0)!
			aw=max(+fw,(dw+fw/2.0),0.0d0)!
			n=max(-fn,(dn-fn/2.0),0.0d0)!
			ap=ae+as+aw+an+(fe-fw+fn-fs)


			!Matrixer:
			bound=0.0
			if (i/=2) then !W
				k=k+1
				vmatv(k)=-aw
				vmati(k)=(j-3)*(ni)+(i-1-1)
				vmatj(k)=(j-3)*(ni)+(i-1)
			else
				ap=ap+mu*dy/(0.5*dx) !stat wall
				bound=mu*(uww(j))*dy/(0.5*dx)
			endif

			if (i/=ni+1) then !E
				k=k+1
				vmatv(k)=-ae
				vmati(k)=(j-3)*(ni)+(i-1+1)
				vmatj(k)=(j-3)*(ni)+(i-1)
			else
				ap=ap+mu*dy/(0.5*dx) !stat wall
				bound=mu*(uwe(j))*dy/(0.5*dx)
				endif

			if (j/=nj+1) then !N
				k=k+1
				vmatv(k)=-an
				vmati(k)=(j-3+1)*(ni)+(i-1)
				vmatj(k)=(j-3)*(ni)+(i-1)
			else
				bound=bound+an*vs(i,j+1)
			endif

			if (j/=3) then !S
				k=k+1
				vmatv(k)=-as
				vmati(k)=(j-3-1)*(ni)+(i-1)
				vmatj(k)=(j-3)*(ni)+(i-1)
			else
				bound=bound+as*vs(i,j-1)
			endif

			k=k+1
			vmatv(k)=ap/alphav!-sp
			vmati(k)=(j-3)*(ni)+(i-1)
			matj(k)=(j-3)*(ni)+(i-1)

			apv(i,j)=ap

			vmrhs((j-3)*(ni)+(i-1))=(ps(i,j-1)-ps(i,j))*dx+((1-alphav)*ap/alphav)*vs(i,j)+bound!+sc

		enddo
	enddo ! End of V momentum


	!Solve U & V, then set results in U* & V*:
2	format(a3,$)
	write(*,2) "U "
	call bicgstab(umatv,umati,umatj,unnzero,umrhs,un,x,bires,bimaxit)
	do k=1,un
		us(mod(k-1,ni-1)+1+2,int((k-1)/(ni-1))+2)=x(k)
	enddo

	write(*,2) "V "
	call bicgstab(vmatv,vmati,vmatj,vnnzero,vmrhs,vn,x,bires,bimaxit)
	do k=1,vn
		vs(mod(k-1,ni)+1+1,int((k-1)/(ni))+3)=x(k)
	enddo

	!Boundary Velocity Correction for Peressure BCs:
	do j=2,nj+1
		if (pbcw(j)==1) then !W
			us(2,j)=us(3,j)+vs(2,j+1)-vs(2,j)
		endif
		if (pbce(j)==1) then !E
			us(ni+2,j)=us(ni+1,j)+vs(ni+1,j)-vs(ni+1,j+1)
		endif
	enddo
	do i=2,ni+1
		if (pbcn(i)==1) then !N
			vs(i,nj+2)=vs(i,nj+1)+us(i,nj+1)-us(i+1,nj+1)
		endif
		if (pbcs(i)==1) then !S
			vs(i,2)=vs(i,3)+us(i,3)-us(i+1,3)
		endif
	enddo


	!Peressure Correction Equations (pp):
	k=0 !counter
	do j=2,nj+1
		do i=2,ni+1
			aw=0.0d0;ae=0.0d0;as=0.0d0;an=0.0d0;ap=0.0d0
			if (i/=2)    aw=rho*dy*dy*alphau/apu(i,j)
			if (i/=ni+1) ae=rho*dy*dy*alphau/apu(i+1,j)
			if (j/=2)    as=rho*dx*dx*alphav/apv(i,j)
			if (j/=nj+1) an=rho*dx*dx*alphav/apv(i,j+1)

			if ((i==2).and.(pbcw(j)==1)) then !P BC on W
				k=k+1
				ap=0.01d0
				pmatv(k)=ap
				pmati(k)=(j-2)*ni+(i-1)
				pmatj(k)=(j-2)*ni+(i-1)
				pmrhs((j-2)*ni+(i-1))=0.0d0
				cycle
			endif

			if ((i==ni+1).and.(pbce(j)==1)) then !P BC on E
				k=k+1
				ap=0.01d0
				pmatv(k)=ap
				pmati(k)=(j-2)*ni+(i-1)
				pmatj(k)=(j-2)*ni+(i-1)
				pmrhs((j-2)*ni+(i-1))=0.0d0
				cycle
			endif

			if ((j==2).and.(pbcs(j)==1)) then !P BC on S
				k=k+1
				ap=0.01d0
				pmatv(k)=ap
				pmati(k)=(j-2)*ni+(i-1)
				pmatj(k)=(j-2)*ni+(i-1)
				pmrhs((j-2)*ni+(i-1))=0.0d0
				cycle
			endif

			if ((j==nj+1).and.(pbcn(j)==1)) then !P BC on N
				k=k+1
				ap=0.01d0
				pmatv(k)=ap
				pmati(k)=(j-2)*ni+(i-1)
				pmatj(k)=(j-2)*ni+(i-1)
				pmrhs((j-2)*ni+(i-1))=0.0d0
				cycle
			endif

			if (i/=2) then
				k=k+1
				pmatv(k)=-aw
				pmati(k)=(j-2)*ni+(i-1-1)
				pmatj(k)=(j-2)*ni+(i-1)
			endif

			if (i/=ni+1) then
				k=k+1
				pmatv(k)=-ae
				pmati(k)=(j-2)*ni+(i-1+1)
				pmatj(k)=(j-2)*ni+(i-1)
			endif

			if (j/=2) then
				k=k+1
				pmatv(k)=-as
				pmati(k)=(j-2-1)*ni+(i-1)
				pmatj(k)=(j-2)*ni+(i-1)
			endif

			if (j/=nj+1) then
				k=k+1
				pmatv(k)=-an
				pmati(k)=(j-2+1)*ni+(i-1)
				pmatj(k)=(j-2)*ni+(i-1)
			endif

			ap=ap+ae+aw+an+as
			k=k+1
			pmatv(k)=ap
			pmati(k)=(j-2)*ni+(i-1)
			pmatj(k)=(j-2)*ni+(i-1)


			pmrhs((j-2)*ni+(i-1))=rho*dy*(us(i,j)-us(i+1,j))+rho*dx*(vs(i,j)-vs(i,j+1))

		enddo
	enddo!The end of pressure Eqns.


	!Solve Pressure Eqns.:
	write(*,2) "P "
	call bicgstab(pmatv,pmati,pmatj,k,pmrhs,pn,x,bires,bimaxit) !pnnzero->k (for press BCs) !! DON'T change k
	do k=1,pn
		pp(mod(k-1,ni)+1+1,int((k-1)/(ni))+1+1)=x(k)
	enddo


	!Continuity Residual Calculation:
	contres=sqrt(dot_product(pmrhs,pmrhs))/pn
	print*,"SIMPLE Iter=",ti,"Continuity Res.=",contres
	write(2,*),ti,contres


	!New Pressure with Corrections (PP):
	do j=2,nj+1
		do i=2,ni+1
			ps(i,j)=ps(i,j)+pp(i,j)*alphap
		enddo
	enddo

	!New U from pressure field:
	do j=2,nj+1
		do i=3,ni+1
			us(i,j)=us(i,j)+(dy*alphau/apu(i,j))*(pp(i-1,j)-pp(i,j))
		enddo
	enddo



	!New V from pressure field:
	do j=3,nj+1
		do i=2,ni+1
			vs(i,j)=vs(i,j)+(dx*alphav/apv(i,j))*(pp(i,j-1)-pp(i,j))
		enddo
	enddo


	if (contres<creslimit) exit !Convergance Checker (by the Continuity Residual)
enddo !!SIMPLE Loop


!Writing Results:
!Extrapolate the P, U & V For Walls:
ps(1,:)=1.5*ps(2,:)-0.5*ps(3,:)
ps(ni+2,:)=1.5*ps(ni+1,:)-0.5*ps(ni,:)
ps(:,1)=1.5*ps(:,2)-0.5*ps(:,3)
ps(:,nj+2)=1.5*ps(:,nj+1)-0.5*ps(:,nj)
us(1,:)=us(2,:);us(ni+3,:)=us(ni+2,:)
vs(:,1)=vs(:,2);vs(:,nj+3)=vs(:,nj+2)
us(:,1)=uws(:);us(:,nj+2)=uwn(:)
do j=1,nj+2
	vs(1,j)=uww(j)
	vs(ni+2,j)=uwe(j)
enddo



!Writing to File:
write(1,*) "VARIABLES= X,Y,U,V,P" !Tecplot header
write(1,*) "ZONE I=",ni+2,",J=",nj+2,",F=POINT"
do j=1,nj+2
	do i=1,ni+2
		idx=0.0d0;jdy=0.0d0!first mesh & last mesh correction
		if (i>1) idx=-dx/2.0d0
		if (i==ni+2) idx=-dx
		if (j>1) jdy=-dy/2.0d0
		if (j==nj+2) jdy=-dy
		write(1,*) (i-1)*dx+idx,(j-1)*dy+jdy,(us(i,j)+us(i+1,j))/2.0,(vs(i,j)+vs(i,j+1))/2.0,ps(i,j)
	enddo
	write(3,*) (j-1)*dy+jdy,"	",(us(ni/2,j)+us(ni/2+1,j))/2.0
enddo
close(1) !"Result.plt"

close(2) !"Residual.plt"

pause "Press any key to exit"
stop
end


![A].b calculater with sparse efficient algorithm (used in BiCGSTAB subroutine):
subroutine matvec(a,ai,aj,nnzero,b,n,c)
	implicit none
	INTEGER, INTENT(inout) :: n !n*n
	INTEGER, INTENT(inout) :: nnzero !none zero
	REAL*8,  INTENT(inout), DIMENSION(1:nnzero) :: a  !A matrix
	INTEGER, INTENT(inout), DIMENSION(1:nnzero) :: ai  !A matrix i
	INTEGER, INTENT(inout), DIMENSION(1:nnzero) :: aj  !A matrix j
	REAL*8,  INTENT(inout), DIMENSION(1:n) :: b !b vector (right hand side)
	REAL*8,  INTENT(inout), DIMENSION(1:n) :: c !Answer vector
	integer i
	c=0.0d0
	do i=1,nnzero
		c(aj(i))=c(aj(i))+a(i)*b(ai(i))
	enddo
endsubroutine


!To solve [A].[x]=[b] system of equations by BiCGSTAB method:
subroutine bicgstab(a,ai,aj,nnzero,b,n,x,bires,bimaxit)
	implicit none
	INTEGER, INTENT(inout) :: n  !Matrix size (n*n)
	INTEGER, INTENT(inout) :: nnzero !No. of none zero elements in [A]
	real*8,  INTENT(inout) :: bires !Max. Residual to converge
	INTEGER, INTENT(inout) :: bimaxit !Max. iterations
	REAL*8,  INTENT(inout), DIMENSION(1:nnzero) :: a  !A matrix
	INTEGER, INTENT(inout), DIMENSION(1:nnzero) :: ai  !A matrix i
	INTEGER, INTENT(inout), DIMENSION(1:nnzero) :: aj  !A matrix j
	REAL*8,  INTENT(inout), DIMENSION(1:n) :: b !Right hand side
	REAL*8,  INTENT(inout), DIMENSION(1:n) :: x !Answer

	allocatable temp(:),r0(:),r(:),p(:),s(:),rr(:)
	integer i
	real(8) temp,r0,r,p,alpha,s,w,rr,beta,res
	allocate (temp(1:n),r0(1:n),r(1:n),p(1:n),s(1:n),rr(1:n))
	i=0

	!BiCGSTAB Algorithm:
	!1:
	x=0.0 !First Guess

	call matvec(a,ai,aj,nnzero,x,n,temp)
	r0=b-temp

	!2:
	p=r0
	r=r0

	!3:
	do
		i=i+1

		!4:
		call matvec(a,ai,aj,nnzero,p,n,temp)
		alpha=(dot_product(r,r0))/(dot_product(temp,r0))

		!5:
		s=r-alpha*temp

		!6:
		call matvec(a,ai,aj,nnzero,s,n,temp)
		w=(dot_product(temp,s))/(dot_product(temp,temp))

		!7:
		x=x+alpha*p+w*s

		!8:
		rr=s-w*temp

		!9:
		beta=((dot_product(rr,r0))/(dot_product(r,r0)))*(alpha/w)
	
		!10:
		call matvec(a,ai,aj,nnzero,p,n,temp)
		p=rr+beta*(p-w*temp)

		!11:
		r=rr

		!12:
		res=dsqrt(dot_product(r,r))/n


		if ((res<bires).or.(i>bimaxit)) then !Convergance or Max iterations Checking
			print*,"BiCGSTAB Res.=",res ,"by",i,"iters"
			exit
		endif
	enddo

endsubroutine
