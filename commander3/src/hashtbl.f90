!================================================================================
!
! Copyright (C) 2020 Institute of Theoretical Astrophysics, University of Oslo.
!
! This file is part of Commander3.
!
! Commander3 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! Commander3 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Commander3. If not, see <https://www.gnu.org/licenses/>.
!
!================================================================================

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

  !!!!!!!!!!!!!!!!!!!!
  ! Hashtable module !
  !!!!!!!!!!!!!!!!!!!!

MODULE hashtbl
  IMPLICIT NONE ! Use strong typing
  INTEGER, PARAMETER :: tbl_size = 50000

  TYPE sllist
     TYPE(sllist), POINTER :: child => NULL()
     CHARACTER(len=:), ALLOCATABLE :: key, val
   CONTAINS
     PROCEDURE :: put  => put_sll
     PROCEDURE :: get  => get_sll
     PROCEDURE :: free => free_sll
  END TYPE sllist

  TYPE hash_tbl_sll
     TYPE(sllist), DIMENSION(:), ALLOCATABLE :: vec
     INTEGER                                 :: vec_len = 0
     LOGICAL                                 :: is_init = .FALSE.
   CONTAINS
     PROCEDURE :: init => init_hash_tbl_sll
     PROCEDURE :: put  => put_hash_tbl_sll
     PROCEDURE :: get  => get_hash_tbl_sll
     PROCEDURE :: free => free_hash_tbl_sll
  END TYPE hash_tbl_sll

  PUBLIC :: hash_tbl_sll

CONTAINS

  RECURSIVE SUBROUTINE put_sll(list,key,val)
    CLASS(sllist),    INTENT(inout) :: list
    CHARACTER(len=*), INTENT(in)    :: key, val
    INTEGER                         :: keylen, vallen

    keylen = LEN(key)
    vallen = LEN(val)
    IF (ALLOCATED(list%key)) THEN
       IF (list%key /= key) THEN
          IF ( .NOT. ASSOCIATED(list%child) ) ALLOCATE(list%child)
          CALL put_sll(list%child,key,val)
       ELSE
          list%val = val
       END IF
    ELSE
       IF (.NOT. ALLOCATED(list%key)) &
            ALLOCATE(CHARACTER(len=keylen) :: list%key)
       list%key = key
       IF (ALLOCATED(list%val)) DEALLOCATE(list%val)
       ALLOCATE(CHARACTER(len=vallen) :: list%val)
       list%val = val
    END IF
  END SUBROUTINE put_sll


  RECURSIVE SUBROUTINE get_sll(list,key,val)
    CLASS(sllist),                 INTENT(in)    :: list
    CHARACTER(len=*),              INTENT(in)    :: key
    CHARACTER(len=:), ALLOCATABLE, INTENT(out)   :: val
    INTEGER                                      :: vallen

    vallen = 0
    if (.not. allocated(list%key)) then
      write(*,*) trim(key)//' not defined'
      stop
    end if
    IF (ALLOCATED(list%key) .AND. (list%key == key)) THEN
       vallen = LEN(list%val)
       IF (ALLOCATED(val)) DEALLOCATE(val)
       ALLOCATE(CHARACTER(len=vallen) :: val)
       val = list%val
    ELSE IF(ASSOCIATED(list%child)) THEN ! keep going
       CALL get_sll(list%child,key,val)
    ELSE ! At the end of the list, no key found
       IF (ALLOCATED(val)) DEALLOCATE(val) ! Exit indication
       RETURN
    END IF
  END SUBROUTINE get_sll


  RECURSIVE SUBROUTINE free_sll(list)
    CLASS(sllist), INTENT(inout) :: list
    IF (ASSOCIATED(list%child)) THEN
       CALL free_sll(list%child)
       DEALLOCATE(list%child)
    END IF
    list%child => NULL()
    IF (ALLOCATED(list%key)) DEALLOCATE(list%key)
    IF (ALLOCATED(list%val)) DEALLOCATE(list%val)
  END SUBROUTINE free_sll

  SUBROUTINE init_hash_tbl_sll(tbl,tbl_len)
    CLASS(hash_tbl_sll),   INTENT(inout) :: tbl
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
  END SUBROUTINE init_hash_tbl_sll

  ! The first part of the hashing procedure using the string
  ! collating sequence
  ELEMENTAL FUNCTION sum_string(str) RESULT(sig)
    CHARACTER(len=*), INTENT(in)   :: str
    INTEGER                        :: sig
    CHARACTER, DIMENSION(LEN(str)) :: tmp
    INTEGER :: i

    FORALL (i=1:LEN(str))
       tmp(i) = str(i:i)
    END FORALL
    sig = SUM(ICHAR(tmp))
  END FUNCTION sum_string


  SUBROUTINE put_hash_tbl_sll(tbl,key,val)
    CLASS(hash_tbl_sll), INTENT(inout) :: tbl
    CHARACTER(len=*),    INTENT(in)    :: key, val
    INTEGER                            :: hash

    hash = MOD(sum_string(key),tbl%vec_len)
    CALL tbl%vec(hash)%put(key=key,val=val)
  END SUBROUTINE put_hash_tbl_sll


  SUBROUTINE get_hash_tbl_sll(tbl,key,val)
    CLASS(hash_tbl_sll),           INTENT(in)    :: tbl
    CHARACTER(len=*),              INTENT(in)    :: key
    CHARACTER(len=:), ALLOCATABLE, INTENT(out)   :: val
    INTEGER                                      :: hash

    hash = MOD(sum_string(key),tbl%vec_len)
    CALL tbl%vec(hash)%get(key=key,val=val)
  END SUBROUTINE get_hash_tbl_sll


  SUBROUTINE free_hash_tbl_sll(tbl)
    CLASS(hash_tbl_sll), INTENT(inout) :: tbl    
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
  END SUBROUTINE free_hash_tbl_sll

END MODULE hashtbl
