module map_editor_simple_ops_mod
  ! use quiet_utils
  use map_editor_utils
  use healpix_types
  implicit none


contains

  subroutine operate_on_single_map(map1, operation, value)
    implicit none

    character(len=*),                   intent(in)    :: operation
    real(dp),                           intent(in)    :: value
    real(dp),         dimension(0:,1:), intent(inout) :: map1
    
    real(dp)     :: missval = -1.6375d30
    integer(i4b) :: i

    select case (trim(operation))
       case ('scale')
          where (abs((map1-missval)/missval) > 1d-5) 
             map1 = value * map1
          end where
       case ('a2t')
          where (abs((map1-missval)/missval) > 1d-5) 
             map1 = ant2thermo(real(value,dp)) * map1
          end where
       case ('t2a')
          write(*,*) ant2thermo(real(value,dp)),'= at2 for freq =',value
          where (abs((map1-missval)/missval) > 1d-5) 
             map1 = map1/ant2thermo(real(value,dp))
          end where
       case ('add_offset')
          where (abs((map1-missval)/missval) > 1d-5)
             map1 = map1 + value
          end where
       case ('log')
          where (map1 > 0.)
             map1 = log10(map1)
          elsewhere
             map1 = missval
          end where
       case ('ln')
          where (map1 > 0.) 
             map1 = log(map1)
          end where
       case ('exp')
          where (abs((map1-missval)/missval) > 1d-5) 
             map1 = exp(map1)
          end where
       case ('abs')
          where (abs((map1-missval)/missval) > 1d-5)
             map1 = abs(map1)
          end where
       case ('asinh')
          where (abs((map1-missval)/missval) > 1d-5)
             map1 = asinh(map1)
          end where
       case ('sqrt')
          where (abs((map1-missval)/missval) > 1d-5 .and. map1 > 0.) 
             map1 = sqrt(map1)
          elsewhere
             map1 = 0.d0
          end where
       case ('inv')
          where (abs((map1-missval)/missval) > 1d-5 .and. map1 /= 0.) 
             map1 = 1./map1
          end where
       case ('hitcount2rms')
          where (abs((map1-missval)/missval) > 1d-5 .and. map1 /= 0.) 
             map1 = value/sqrt(map1)
          elsewhere
             map1 = missval
          end where
       case ('max_scalar')
          where (abs((map1-missval)/missval) > 1d-5)
             map1 = max(value, map1)
          elsewhere
             map1 = missval
          end where
       case ('min_scalar')
          where (abs((map1-missval)/missval) > 1d-5)
             map1 = min(value, map1)
          elsewhere
             map1 = missval
          end where
       case ('missing2mask')
          where ((abs(map1-missval))/abs(missval)<1.d-5)
             map1 = 0.d0
          elsewhere
             map1 = 1.d0
          end where
       case ('missing2val')
          where ((abs(map1-missval))/abs(missval)<1.d-5)
             map1 = real(value,dp)
          elsewhere
             map1 = map1
          end where
       case ('zero2mask')
          where (map1 == 0.d0)
             map1 = 0.d0
          elsewhere
             map1 = 1.d0
          end where
       case ('QU2P')
          if (size(map1,2) == 3) then
             do i = 0, size(map1,1)-1
                if (all(abs((map1(i,2:3)-missval)/missval) > 1d-5)) then
                   map1(i,1)   = sqrt(map1(i,2)**2 + map1(i,3)**2)
                else
                   map1(i,1)   = missval
                end if
                !map1(i,2:3) = missval
             end do
          else
             map1 = 0.d0
          end if
       case ('rms2mask')
          where (abs(map1) >= value .or. abs((map1-missval)/missval) < 1d-5)
             map1 = 0.d0
          elsewhere
             map1 = 1.d0
          end where
!!$          if (size(map1(1,:)) == 3) then
!!$             do i = 0, size(map1(:,1))-1
!!$                if (any(map1(i,2:3) == 0.d0)) map1(i,2:3) = 0.d0
!!$             end do
!!$          end if
       case ('amp2mask')
          where (map1 < value .or. abs((map1-missval)/missval) < 1d-5) 
             map1 = 0.d0
          elsewhere
             map1 = 1.d0
          end where
          if (size(map1(1,:)) == 3) then
             do i = 0, size(map1(:,1))-1
                if (any(map1(i,2:3) == 0.d0)) map1(i,2:3) = 0.d0
             end do
          end if
       case ('hitcount2mask')
          write(*,*) 'Number of exposed pixels  = ', count(abs((map1-missval)/missval) > 1d-5)
          where (map1 <= maxval(map1) / value**2 .or. abs((map1-missval)/missval) < 1d-5) 
             map1 = 0.d0
          elsewhere
             map1 = 1.d0
          end where
          write(*,*) 'Number of unmasked pixels = ', sum(map1)
       case ('invert_mask')
          where (map1 < 0.5d0 .and. abs((map1-missval)/missval) > 1d-5)
             map1 = 1.d0
          elsewhere
             map1 = 0.d0
          end where
    end select

  end subroutine operate_on_single_map

  subroutine operate_on_two_maps(map1, map2, operation)
    implicit none

    character(len=*),                   intent(in)    :: operation
    real(dp),         dimension(0:,1:), intent(inout) :: map1
    real(dp),         dimension(0:,1:), intent(in)    :: map2

    integer(i4b) :: i, j, npix, nmaps
    real(dp)     :: missval = -1.6375d30

    npix  = size(map1(:,1))
    nmaps = min(size(map1,2), size(map2,2))

    do i = 1, nmaps
       do j = 0, npix-1
          if (abs((map1(j,i)-missval)/missval) > 1d-5 .and. abs((map2(j,i)-missval)/missval) > 1d-5) then
             select case (trim(operation))
             case ('add')
                map1(j,i) = map1(j,i) + map2(j,i)
             case ('subtract')
                map1(j,i) = map1(j,i) - map2(j,i)
             case ('multiply')
                map1(j,i) = map1(j,i) * map2(j,i)
             case ('divide')
                map1(j,i) = map1(j,i) / map2(j,i)
             case ('half_sum')
                map1(j,i) = 0.5 * (map1(j,i) + map2(j,i))
             case ('half_diff')
                map1(j,i) = 0.5 * (map1(j,i) - map2(j,i))
             case ('max')
                map1(j,i) = max(map1(j,1), map2(j,2))
             case ('min')
                map1(j,i) = min(map1(j,i), map2(j,i))
             case ('splice')
                map1(j,i) = map2(j,i)
             end select
          else
             if (trim(operation) /= 'splice') map1(j,i) = missval
          end if
       end do
    end do

  end subroutine operate_on_two_maps

end module map_editor_simple_ops_mod
