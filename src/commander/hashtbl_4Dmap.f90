  !!!!!!!!!!!!!!!!!!!!
  ! Hashtable module for pixel summation!
  !!!!!!!!!!!!!!!!!!!!

! Based on module by Izaak Beekman, modified by Hans Kristian Eriksen, 2019

! Module implementing an OO hash table (dictionary) in Fortran 2003.
! Compiles and runs with accompanying test program under the Intel 
! Fortran Compiler, version 11.1.046

! Copyright (c) Izaak Beekman 2010

    ! This program is free software: you can redistribute it and/or modify
    ! it under the terms of the GNU Lesser General Public License as published by
    ! the Free Software Foundation, either version 3 of the License, or
    ! (at your option) any later version.

    ! This program is distributed in the hope that it will be useful,
    ! but WITHOUT ANY WARRANTY; without even the implied warranty of
    ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    ! GNU Lesser General Public License for more details.

    ! You should have received a copy of the GNU Lesser General Public License
    ! along with this program.  If not, see <http://www.gnu.org/licenses/>.

MODULE hashtbl_4Dmap
  IMPLICIT NONE ! Use strong typing
  INTEGER, PARAMETER :: tbl_size = 50000

  TYPE sllist_4Dmap
     TYPE(sllist_4Dmap), POINTER :: child => NULL()
     INTEGER :: key(2) = -1
     REAL    :: weight
     REAL, ALLOCATABLE, DIMENSION(:)   :: m_tot
   CONTAINS
     PROCEDURE :: put  => put_sll_4Dmap
     PROCEDURE :: get  => get_sll_4Dmap
     PROCEDURE :: free => free_sll_4Dmap
     PROCEDURE :: gen_len  => get_n_sub
  END TYPE sllist_4Dmap

  TYPE hash_tbl_4Dmap_sll
     TYPE(sllist_4Dmap), DIMENSION(:), ALLOCATABLE :: vec
     INTEGER                                 :: vec_len = 0
     LOGICAL                                 :: is_init = .FALSE.
   CONTAINS
     PROCEDURE :: init => init_hash_tbl_4Dmap_sll
     PROCEDURE :: put  => put_hash_tbl_4Dmap_sll
     PROCEDURE :: get  => get_hash_tbl_4Dmap_sll
     PROCEDURE :: free => free_hash_tbl_4Dmap_sll
     PROCEDURE :: get_n_elements
     PROCEDURE :: linearize 
  END TYPE hash_tbl_4Dmap_sll

  PUBLIC :: hash_tbl_4Dmap_sll, sllist_4Dmap

CONTAINS

  RECURSIVE SUBROUTINE put_sll_4Dmap(list, key, m)
    CLASS(sllist_4Dmap),    INTENT(inout) :: list
    INTEGER              :: key(2)
    REAL, DIMENSION(:)   :: m

    IF (list%key(1) /= -1) THEN
       IF (all(list%key == key)) THEN
          list%weight = list%weight + 1.
          list%m_tot  = list%m_tot  + m
       ELSE
          IF (.NOT. ASSOCIATED(list%child)) ALLOCATE(list%child)
          CALL put_sll_4Dmap(list%child, key, m)
       END IF
    ELSE
       ALLOCATE(list%m_tot(size(m)))
       list%key    = key
       list%weight = 1.
       list%m_tot  = m
       nullify(list%child)
    END IF
  END SUBROUTINE put_sll_4Dmap


  RECURSIVE SUBROUTINE get_sll_4Dmap(list, key, weight, m_tot)
    CLASS(sllist_4Dmap),              INTENT(in)  :: list
    INTEGER,                          INTENT(in)  :: key(2)
    REAL,                             INTENT(out) :: weight
    REAL,             DIMENSION(:),   INTENT(out) :: m_tot

    IF (list%key(1) /= -1 .AND. all(list%key == key)) THEN
       weight = list%weight
       m_tot  = list%m_tot
    ELSE IF (ASSOCIATED(list%child)) THEN ! keep going
       CALL get_sll_4Dmap(list%child, key, weight, m_tot)
    ELSE ! At the end of the list, no key found
       write(*,*) 'HASHTBL_4DMAP error -- entry not found = ', key
       RETURN
    END IF
  END SUBROUTINE get_sll_4Dmap

  RECURSIVE FUNCTION get_n_sub(list)
    CLASS(sllist_4Dmap),              INTENT(in)  :: list
    INTEGER                                       :: get_n_sub

    if (.NOT. ASSOCIATED(list%child)) then
       get_n_sub = 1
    else
       get_n_sub = get_n_sub(list%child) + 1
    end if
  END FUNCTION get_n_sub


  RECURSIVE SUBROUTINE free_sll_4Dmap(list)
    CLASS(sllist_4Dmap), INTENT(inout) :: list
    IF (ASSOCIATED(list%child)) THEN
       CALL free_sll_4Dmap(list%child)
       DEALLOCATE(list%child)
    END IF
    list%child => NULL()
    IF (ALLOCATED(list%m_tot))  DEALLOCATE(list%m_tot)
  END SUBROUTINE free_sll_4Dmap

  SUBROUTINE init_hash_tbl_4Dmap_sll(tbl, tbl_len)
    CLASS(hash_tbl_4Dmap_sll),   INTENT(inout) :: tbl
    INTEGER,     OPTIONAL, INTENT(in)    :: tbl_len

    IF (ALLOCATED(tbl%vec)) DEALLOCATE(tbl%vec)
    IF (PRESENT(tbl_len)) THEN
       ALLOCATE(tbl%vec(0:tbl_len-1))
       tbl%vec_len = tbl_len
    ELSE
       ALLOCATE(tbl%vec(0:tbl_size-1))
       tbl%vec_len = tbl_size
    END IF
    tbl%is_init = .TRUE.
  END SUBROUTINE init_hash_tbl_4Dmap_sll

  SUBROUTINE put_hash_tbl_4Dmap_sll(tbl, key, m)
    CLASS(hash_tbl_4Dmap_sll),          INTENT(inout) :: tbl
    INTEGER,                            INTENT(in)    :: key(2)
    REAL,                DIMENSION(:),  INTENT(in)    :: m
    INTEGER :: hash

    hash = MOD(key(1),tbl%vec_len)
    CALL tbl%vec(hash)%put(key, m)
  END SUBROUTINE put_hash_tbl_4Dmap_sll


  SUBROUTINE get_hash_tbl_4Dmap_sll(tbl, key, weight, m_tot)
    CLASS(hash_tbl_4Dmap_sll),        INTENT(in)  :: tbl
    INTEGER,                          INTENT(in)  :: key(2)
    REAL,                             INTENT(OUT) :: weight
    REAL,             DIMENSION(:),   INTENT(OUT) :: m_tot
    INTEGER :: hash

    hash = MOD(key(1),tbl%vec_len)
    CALL tbl%vec(hash)%get(key, weight, m_tot)
  END SUBROUTINE get_hash_tbl_4Dmap_sll

  SUBROUTINE free_hash_tbl_4Dmap_sll(tbl)
    CLASS(hash_tbl_4Dmap_sll), INTENT(inout) :: tbl    
    INTEGER     :: i, low, high

    low  = LBOUND(tbl%vec,dim=1)
    high = UBOUND(tbl%vec,dim=1) 
    IF (ALLOCATED(tbl%vec)) THEN
       DO i=low,high
          CALL tbl%vec(i)%free()
       END DO
       DEALLOCATE(tbl%vec)
    END IF
    tbl%is_init = .FALSE.
  END SUBROUTINE free_hash_tbl_4Dmap_sll

  function get_n_elements(tbl)
    CLASS(hash_tbl_4Dmap_sll), INTENT(in) :: tbl    
    INTEGER                               :: get_n_elements
    INTEGER     :: i, n, low, high
    CLASS(sllist_4Dmap), pointer :: list

    low  = LBOUND(tbl%vec,dim=1)
    high = UBOUND(tbl%vec,dim=1) 
    n    = 0
    IF (ALLOCATED(tbl%vec)) THEN
       DO i = low, high
          if (tbl%vec(i)%key(1) == -1) cycle
          n = n + get_n_sub(tbl%vec(i))
       END DO
    END IF
    get_n_elements = n
  END function get_n_elements

  subroutine linearize(tbl, pixel, ipsi, weight, m_tot)
    CLASS(hash_tbl_4Dmap_sll), INTENT(in)  :: tbl    
    INTEGER, dimension(:),     INTENT(out) :: pixel, ipsi
    REAL,    dimension(:),     INTENT(out) :: weight
    REAL,    dimension(:,:),   INTENT(OUT) :: m_tot
    INTEGER     :: i, n, low, high
    CLASS(sllist_4Dmap), pointer :: list

    low  = LBOUND(tbl%vec,dim=1)
    high = UBOUND(tbl%vec,dim=1) 
    n    = 0
    IF (ALLOCATED(tbl%vec)) THEN
       DO i = low, high
          if (tbl%vec(i)%key(1) /= -1) then
             n          = n+1
             pixel(n)   = tbl%vec(i)%key(1)
             ipsi(n)    = tbl%vec(i)%key(2)
             weight(n)  = tbl%vec(i)%weight
             m_tot(n,:) = tbl%vec(i)%m_tot
             if (associated(tbl%vec(i)%child)) then
                list => tbl%vec(i)%child
                do while (associated(list))
                   n          = n+1
                   pixel(n)   = list%key(1)
                   ipsi(n)    = list%key(2)
                   weight(n)  = list%weight
                   m_tot(n,:) = list%m_tot
                   list => list%child
                end do
             end if
          end if
       END DO
    END IF

  END subroutine linearize

END MODULE hashtbl_4Dmap
