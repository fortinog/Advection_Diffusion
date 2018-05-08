module type_defs
!---------------------------------------------------------------------
!  Module which defines precisions
!  Reference: http://fortranwiki.org/fortran/show/Real+precision
!  Ref: https://kiwi.atmos.colostate.edu/fortran/docs/fortran2012.key-6.pdf
!---------------------------------------------------------------------

    integer, parameter :: sp = selected_real_kind(6, 37)
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, parameter :: qp = selected_real_kind(33, 4931)
    
end module type_defs
