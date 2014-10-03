MODULE ITT_FORTRAN

#ifdef VTUNE_PROFILE
    USE, INTRINSIC :: ISO_C_BINDING

    INTERFACE
    
        SUBROUTINE FORTRAN_ITT_RESUME() &
            BIND(C, NAME='fortran_itt_resume')
!            BIND(C, NAME='__itt_resume_ptr__3_0')
        END SUBROUTINE FORTRAN_ITT_RESUME
     
        SUBROUTINE FORTRAN_ITT_PAUSE() &
            BIND(C, NAME='fortran_itt_pause')
!            BIND(C, NAME='__itt_pause_ptr__3_0')
        END SUBROUTINE
    END INTERFACE

#endif

END MODULE

