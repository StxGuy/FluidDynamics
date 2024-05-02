!*************************************************
!*** Lattice Boltzman Fluid Dynamics Simulator ***
!***   By: Prof. Carlo daCunha, Ph.D. (2024)   ***
!*************************************************
!
! Inspired by: https://physics.weber.edu/schroeder/fluids/
! Theory.....: https://www.ndsu.edu/fileadmin/physics.ndsu.edu/Wagner/LBbook.pdf
! 
! To be done: Convert to 2D arrays and use GPU


! Module for fluid dynamics
module FluidModule
    implicit none
    
    ! Calculate some parameters to speed up calculations
    real,private,parameter		:: one9th   = 1.0/9
    real,private,parameter		:: four9ths = 4.0/9
    real,private,parameter		:: one36th  = 1.0/36    
    
    ! Class
    type Fluid
        integer,private     :: Lx,Ly
        double precision,allocatable    :: North(:),South(:)
        double precision,allocatable    :: East(:),West(:)
        double precision,allocatable    :: NorthWest(:),NorthEast(:)
        double precision,allocatable    :: SouthWest(:),SouthEast(:)
        double precision,allocatable    :: n0(:)
        logical,allocatable :: barrier(:)
        real                :: omega

        contains

        procedure           :: Stream
        procedure           :: Collide
        procedure           :: setEq
        procedure           :: initFluid
        procedure           :: setBoundaries
        procedure           :: saveFluid
        procedure           :: mass
        
        final               :: FluidDestructor
    end type Fluid

    interface Fluid
        procedure           :: FluidConstructor
    end interface

    contains

    ! Constructor for the Lattice Boltzmann procedure
    function FluidConstructor(nx,ny,viscosity) result(self)
        integer,intent(in)  :: nx,ny
        real,intent(in)     :: viscosity
        type(Fluid)         :: self
        integer             :: sz,y,idx
        integer,parameter   :: barrierSize = 8

        
        self%Lx = nx
        self%Ly = ny
        sz = nx*ny
        self%omega = 1.0/(3*viscosity + 0.5)

        allocate(self%North(sz))
        allocate(self%South(sz))
        allocate(self%East(sz))
        allocate(self%West(sz))
        allocate(self%NorthWest(sz))
        allocate(self%NorthEast(sz))
        allocate(self%SouthWest(sz))
        allocate(self%SouthEast(sz))
        allocate(self%barrier(sz))
        allocate(self%n0(sz))

        ! Barrier
        self%barrier = .false.
        do y = ishft(ny,-1)-barrierSize,ishft(ny,-1)+barrierSize
            idx = int(real(ny)/real(3.0)) + (y-1)*self%Lx
			self%barrier(idx) = .true.
		end do
		
		call self%initFluid()

		write(*,*) ""
		write(*,*) "*************************************************"
		write(*,*) "*** Lattice Boltzman Fluid Dynamics Simulator ***"
		write(*,*) "***   By: Prof. Carlo daCunha, Ph.D. (2024)   ***"
		write(*,*) "*************************************************"
		write(*,*) ""
		write(*,*) "Simulation:"
    	write(*,*) "  Grid.......: ",nX, nY
    	write(*,*) "  Viscosity..: ",viscosity
    	write(*,*) "  Omega......: ",self%omega
    	!write(*,*) "  Fluid Speed: ",ux, uy
    	write(*,*) ""        		
		
    end function

    ! Destructor for the Lattice Boltzmann procedure
    subroutine FluidDestructor(self)
        type(Fluid) :: self

        deallocate(self%North)
        deallocate(self%South)
        deallocate(self%East)
        deallocate(self%West)
        deallocate(self%NorthWest)
        deallocate(self%NorthEast)
        deallocate(self%SouthWest)
        deallocate(self%SouthEast)
        deallocate(self%barrier)
        deallocate(self%n0)
    end subroutine

    ! Collision procedure that relaxes the distribution towards equilibrium 
    subroutine Collide(self)
        class(Fluid),intent(inout)  :: self
        integer                     :: x,y,idx
        double precision            :: rho,ux,uy
        double precision            :: one9thrho,one36thrho,four9thsrho
        double precision            :: ux2,uy2,u2,uxuy2,ux3,uy3,u215

        do y = 2,(self%Ly-1)
        do x = 2,(self%Lx-1)
            idx = x + (y-1)*self%Lx
            rho = self%n0(idx)+self%North(idx)+self%South(idx)+self%West(idx)+self%East(idx)
            rho = rho + self%NorthEast(idx)+self%NorthWest(idx)
            rho = rho + self%SouthEast(idx)+self%SouthWest(idx)

            ux = self%East(idx)+self%NorthEast(idx)+self%SouthEast(idx)
            ux = ux-(self%West(idx)+self%NorthWest(idx)+self%SouthWest(idx))
            ux = ux/rho

            uy = self%North(idx)+self%NorthEast(idx)+self%NorthWest(idx)
            uy = uy-(self%South(idx)+self%SouthEast(idx)+self%SouthWest(idx))
            uy = uy/rho

            one9thrho = one9th*rho
            one36thrho = one36th*rho
            four9thsrho = four9ths*rho
            ux3 = 3*ux
            uy3 = 3*uy
            ux2 = ux*ux
            uy2 = uy*uy
            uxuy2 = 2*ux*uy
            u2 = ux2 + uy2
            u215 = 1.5*u2

            self%n0(idx)        = (1.0-self%omega)*self%n0(idx)       + self%omega*four9thsrho*(1                       -u215)
            self%East(idx)      = (1.0-self%omega)*self%East(idx)     +   self%omega*one9thrho*(1+ux3    +4.5*ux2       -u215)
            self%West(idx)      = (1.0-self%omega)*self%West(idx)     +   self%omega*one9thrho*(1-ux3    +4.5*ux2       -u215)
            self%North(idx)     = (1.0-self%omega)*self%North(idx)    +   self%omega*one9thrho*(1    +uy3+4.5*uy2       -u215)
            self%South(idx)     = (1.0-self%omega)*self%South(idx)    +   self%omega*one9thrho*(1    -uy3+4.5*uy2       -u215)
            self%NorthEast(idx) = (1.0-self%omega)*self%NorthEast(idx) + self%omega*one36thrho*(1+ux3+uy3+4.5*(u2+uxuy2)-u215)
            self%SouthEast(idx) = (1.0-self%omega)*self%SouthEast(idx) + self%omega*one36thrho*(1+ux3-uy3+4.5*(u2-uxuy2)-u215)
            self%NorthWest(idx) = (1.0-self%omega)*self%NorthWest(idx) + self%omega*one36thrho*(1-ux3+uy3+4.5*(u2-uxuy2)-u215)
            self%SouthWest(idx) = (1.0-self%omega)*self%SouthWest(idx) + self%omega*one36thrho*(1-ux3-uy3+4.5*(u2+uxuy2)-u215)
        end do
        end do

        ! Boundary conditions
        do y = 2,(self%Ly-2)
            idx = self%Lx + (y-1)*self%Lx
            self%West(idx)      = self%West(idx-1)
            self%NorthWest(idx) = self%NorthWest(idx-1)
            self%SouthWest(idx) = self%SouthWest(idx-1)
        end do
    end subroutine    
    
    ! Stream components to neighborhood
    subroutine Stream(self)
        class(Fluid), intent(inout) :: self
        integer                     :: idx, x, y

        ! Northwest Corner
        do y = (self%Ly-1),2,-1
        do x = 2,(self%Lx-1)
            idx = x + (y-1)*self%Lx
            self%North(idx) = self%North(idx - self%Lx)
            self%NorthWest(idx) = self%NorthWest(idx + 1 - self%Lx)
        end do
        end do

        ! Northeast Corner
        do y = (self%Ly-1),2,-1
        do x = (self%Lx-1),2,-1
            idx = x + (y-1)*self%Lx
            self%East(idx) = self%East(idx-1)
            self%NorthEast(idx) = self%NorthEast(idx - 1 - self%Lx)
        end do
        end do

        ! Southeast Corner
        do y = 2, (self%Ly-1)
        do x = (self%Lx-1),2,-1
            idx = x + (y-1)*self%Lx
            self%South(idx) = self%South(idx + self%Lx)
            self%SouthEast(idx) = self%SouthEast(idx - 1 + self%Lx)
        end do
        end do

        ! Southwest Corner
        do y = 1, (self%Ly-1)
        do x = 1, (self%Lx-1)
            idx = x + (y-1)*self%Lx
            self%West(idx) = self%West(idx+1)
            self%SouthWest(idx) = self%SouthWest(idx + 1 + self%Lx)
        end do
        end do

        ! Bounce-back
        do y = 2,(self%Ly-1)
        do x = 2,(self%Lx-1)
            idx = x + (y-1)*self%Lx
            if (self%barrier(idx) .eqv. .true.) then
                self%East(idx + 1)                = self%West(idx)
                self%West(idx - 1)                = self%East(idx)
                self%North(idx + self%Lx)         = self%South(idx)
                self%South(idx - self%Lx)         = self%North(idx)
                self%NorthEast(idx + 1 + self%Lx) = self%SouthWest(idx)
                self%NorthWest(idx - 1 + self%Lx) = self%SouthEast(idx)
                self%SouthEast(idx + 1 - self%Lx) = self%NorthWest(idx)
                self%SouthWest(idx - 1 - self%Lx) = self%NorthEast(idx)
            end if
        end do
        end do
    end subroutine

    ! Auxiliary procedure to monitor the total mass of the system
    function mass(self) result(rho)
        class(Fluid),intent(in)     :: self
        double precision            :: rho
        integer                     :: x,y,idx
        
        rho = 0.0
        do y = 2,(self%Ly-1)
        do x = 2,(self%Lx-1)
            idx = x + (y-1)*self%Lx
            rho = rho + self%n0(idx)+self%North(idx)+self%South(idx)+self%West(idx)+self%East(idx)
            rho = rho + self%NorthEast(idx)+self%NorthWest(idx)
            rho = rho + self%SouthEast(idx)+self%SouthWest(idx)
        end do
        end do
        
        rho = rho/((self%Lx-2)*(self%Ly-2))
    end
        
    ! Subroutine to set a cell to new values of speed and mass
    subroutine setEq(self,x,y,newux,newuy,newrho)
        class(Fluid),intent(inout)  :: self
        integer,intent(in)          :: x,y
        real,intent(in)             :: newux,newuy
        real,intent(in)             :: newrho
        integer                     :: idx
        double precision            :: one9thrho,one36thrho,four9thsrho
        double precision            :: ux3,uy3,ux2,uy2,uxuy2,u2,u215

        idx = x + (y-1)*self%Lx

        one9thrho = one9th*newrho
        one36thrho = one36th*newrho
        four9thsrho = four9ths*newrho
        ux3 = 3*newux
        uy3 = 3*newuy
        ux2 = newux*newux
        uy2 = newuy*newuy
        uxuy2 = 2*newux*newuy
        u2 = ux2 + uy2
        u215 = 1.5*u2

        self%n0(idx)        = four9thsrho*(1                       -u215)
        self%East(idx)      =   one9thrho*(1+ux3    +4.5*ux2       -u215)
        self%West(idx)      =   one9thrho*(1-ux3    +4.5*ux2       -u215)
        self%North(idx)     =   one9thrho*(1    +uy3+4.5*uy2       -u215)
        self%South(idx)     =   one9thrho*(1    -uy3+4.5*uy2       -u215)
        self%NorthEast(idx) =  one36thrho*(1+ux3+uy3+4.5*(u2+uxuy2)-u215)
        self%SouthEast(idx) =  one36thrho*(1+ux3-uy3+4.5*(u2-uxuy2)-u215)
        self%NorthWest(idx) =  one36thrho*(1-ux3+uy3+4.5*(u2-uxuy2)-u215)
        self%SouthWest(idx) =  one36thrho*(1-ux3-uy3+4.5*(u2+uxuy2)-u215)
    end
    
    ! Initialize fluid
    subroutine initFluid(self)
        class(Fluid),intent(inout)  :: self
        integer                     :: x,y
        
        do y = 1,self%Ly
        do x = 1,self%Lx
            call self%setEq(x,y,0.1,0.0,1.0)
        end do
        end do
    end subroutine
    
    ! Set boundaries of the simulation
    subroutine setBoundaries(self)
        class(Fluid),intent(inout)  :: self
        integer                     :: x,y
        
        do x = 1,self%Lx
            call self%setEq(x,1,0.1,0.0,1.0)
            call self%setEq(x,self%Ly-2,0.1,0.0,1.0)
        end do
        do y = 2,self%Ly-1
            call self%setEq(1,y,0.1,0.0,1.0)
            call self%setEq(self%Lx-2,y,0.1,0.0,1.0)
        end do
    end subroutine
    
    ! Save a set of still frames that can be put together into an animation later
    subroutine saveFluid(self,i)
		implicit none
		
		character(len=20)				:: fname
		character(len=10),parameter     :: fname1 = "video/test"
		character(len=4),parameter      :: fname2 = ".txt"
		character(len=4)				:: str
		integer,intent(in)				:: i
		integer							:: x,y,idx
		double precision                :: rho, ux, uy	
		logical							:: existent
		class(Fluid),intent(in)   		:: self
		double precision                :: speed, maxspeed
		
				
		! Filename string handling
		write(str,'(i0)') i
		fname = fname1 // trim(adjustl(str)) // fname2
			
		! File procedure
		inquire(file=fname,exist=existent)
		
		if (existent) then
			open(10,file=fname,status="old")
		else
			open(10,file=fname,status="new")
		end if

		! Save data
		maxspeed = 0.0
		do x = 2,(self%Lx-1)
            do y = 2,(self%Ly-1)
                idx = x + (y-1)*self%Lx
                
                rho = self%n0(idx)+self%North(idx)+self%South(idx)+self%West(idx)+self%East(idx)
                rho = rho + self%NorthEast(idx)+self%NorthWest(idx)
                rho = rho + self%SouthEast(idx)+self%SouthWest(idx)

                ux = self%East(idx)+self%NorthEast(idx)+self%SouthEast(idx)
                ux = ux-(self%West(idx)+self%NorthWest(idx)+self%SouthWest(idx))
                ux = ux/rho

                uy = self%North(idx)+self%NorthEast(idx)+self%NorthWest(idx)
                uy = uy-(self%South(idx)+self%SouthEast(idx)+self%SouthWest(idx))
                uy = uy/rho
                
                speed = sqrt(ux*ux + uy*uy)
                if (speed > maxspeed) then
                    maxspeed = speed
                end if
				write(10,'(E20.4)',advance='no') speed
			end do
			write(10,*)
		end do
			
		close(10)
		!write(*,*) maxspeed
	end subroutine	
end module

!==========================================================!
! MAIN PROGRAM
!==========================================================!
program LatticeBoltzmann
    use FluidModule

    type(Fluid) :: F
    integer     :: it

    F = Fluid(200,80,0.02)
    
    call F%setBoundaries()
    do it = 1,9999
        call F%collide()
        call F%stream()
        
        call F%saveFluid(it)
    end do
    
    
end program
