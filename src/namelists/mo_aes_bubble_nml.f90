! @brief Read configuration parameters as Fortran namelist from an external file.
!
!
! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

MODULE mo_aes_bubble_nml

  USE mo_aes_bubble_config, ONLY: aes_bubble_config, init_aes_bubble_config, &
       &                          eval_aes_bubble_config, print_aes_bubble_config
  USE mo_process_nml,       ONLY: process_nml
  
  IMPLICIT NONE                                                                      
  PRIVATE                                                                            
  PUBLIC :: process_aes_bubble_nml                                                    
                                                                                     
  NAMELIST /aes_bubble_nml/ aes_bubble_config                             
                                                                                     
CONTAINS                                                                             
                                                                                     
  SUBROUTINE process_aes_bubble_nml(filename)                                         
    !                                                                                
    CHARACTER(LEN=*), INTENT(in) :: filename                                         
    !                                                                                
    CALL init_aes_bubble_config                                                       
    !                                                                                
    CALL process_nml(filename, 'aes_bubble_nml', nml_read, nml_write)                 
    !                                                                                
    CALL eval_aes_bubble_config
    CALL print_aes_bubble_config
  CONTAINS                                                                           
    !                                                                                
    SUBROUTINE nml_read(funit)                                                       
      INTEGER, INTENT(in) :: funit                                                   
      READ(funit, NML=aes_bubble_nml)                                                 
    END SUBROUTINE nml_read                                                          
    !                                                                                
    SUBROUTINE nml_write(funit)                                                      
      INTEGER, INTENT(in) :: funit                                                   
      WRITE(funit, NML=aes_bubble_nml)                                                
    END SUBROUTINE nml_write                                                         
    !                                                                                
  END SUBROUTINE process_aes_bubble_nml
                                              
END MODULE mo_aes_bubble_nml
