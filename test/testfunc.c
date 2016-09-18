/* This file exists to provide serial, mock versions of the common utility functions for testing purposes*/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

void endrun(int ierr, const char * fmt, ...)
{
    va_list va;
    va_start(va, fmt);
    vprintf(fmt, va);
    va_end(va);
    printf("This would have ended the run!");
}

void message(int ierr, const char * fmt, ...)
{
    va_list va;
    va_start(va, fmt);
    vprintf(fmt, va);
    va_end(va);
}

void *mymalloc_fullinfo(const char *varname, size_t n, const char *func, const char *file, int line)
{
    return malloc(n);
}
