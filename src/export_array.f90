      module export_array_mod
      implicit none

      logical :: export_mid_plane = .false.
      logical :: export_mid_line  = .false.
      private
      public :: exp_3D_1C_GF ! 3D Fields
      public :: exp_2D_1C_GF ! 2D Fields
      public :: exp_1D_1C_GF ! 1D Fields

      contains

      ! ***********************************************************************
      ! ****************************** 3D FIELDS ******************************
      ! ***********************************************************************

      subroutine exp_3D_1C_GF(g,DL,t,pad,un,u)
        implicit none
        type(grid),intent(in) :: g
        type(data_location),intent(in) :: DL
        integer,intent(in) :: un,pad,t
        type(grid_field),intent(in) :: u
        type(array),dimension(3) :: h
        integer,dimension(3) :: s
        integer :: i,j,k
        s = shape(f)
        do k=1,s(3); do j=1,s(2); do i=1,s(1)
          write(un,*) f(i,j,k)
        enddo; enddo; enddo
        do i=1,3; call delete(h(i)); enddo
      end subroutine


      end module