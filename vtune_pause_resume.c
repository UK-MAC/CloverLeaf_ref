#include "ittnotify.h"
 
void fortran_itt_resume()
{
    __itt_resume();
}
 
void fortran_itt_pause()
{
    __itt_pause();
}

