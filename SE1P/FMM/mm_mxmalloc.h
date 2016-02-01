/* Copyright (C) 2004 Free Software Foundation, Inc.

   This file is part of GCC.

   GCC is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3 or later, or (at your option)
   any later version.

   GCC is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with GCC; see the file COPYING.  If not, write to
   the Free Software Foundation, 59 Temple Place - Suite 330,
   Boston, MA 02111-1307, USA.  */

/* As a special exception, if you include this header file into source
   files compiled by GCC, this header file does not by itself cause
   the resulting executable to be covered by the GNU General Public
   License.  This exception does not however invalidate any other
   reasons why the executable file might be covered by the GNU General
   Public License.  */

#ifndef _MM_MEXMALLOC_H_INCLUDED
#define _MM_MEXMALLOC_H_INCLUDED

#include <stdlib.h>
#include <errno.h>
#include <stdio.h>

/*__inline__ is gcc-specific, turn it off if on Windows*/
#ifndef _WIN32
__inline__ 
#endif
double*
_mm_mxMalloc (size_t size, size_t align)
{
  double * malloc_ptr;
  double * aligned_ptr;

  /* Error if align is not a power of two.  */
  if (align & (align - 1))
    {
      printf("Error in _mm_mxmalloc : alignment not power of two.");
      errno = EINVAL;
      return ((void*) 0);
    }

  if (size == 0) {
    printf("Error in _mm_mxmalloc : cannot allocate zero bytes.");    
    return ((void *) 0);
  }
 /* Assume malloc'd pointer is aligned at least to sizeof (void*).
    If necessary, add another sizeof (void*) to store the value
    returned by malloc. Effectively this enforces a minimum alignment
    of sizeof double. */
    if (align < 2 * sizeof (double *))
      align = 2 * sizeof (double *);

  malloc_ptr = (double*) malloc (size + align);
  if (!malloc_ptr)
    return ((double *) 0);

  /* Align  We have at least sizeof (void *) space below malloc'd ptr. */
  aligned_ptr = (double *) (((size_t) malloc_ptr + align)
                          & ~((size_t) (align) - 1));

  /* Store the original pointer just before p.  */
  ((double **) aligned_ptr) [-1] = malloc_ptr;

  return aligned_ptr;
}


/*__inline__ is gcc-specific, turn it off if on Windows*/
#ifndef _WIN32
__inline__ 
#endif
double*
_mm_mxCalloc (size_t size, size_t typesize ,size_t align)
{
  double * calloc_ptr;
  double * aligned_ptr;

  /* Error if align is not a power of two.  */
  if (align & (align - 1))
    {
      printf("Error in _mm_mxcalloc : alignment not power of two.");
      errno = EINVAL;
      return ((double*) 0);
    }

  if (size == 0) {
      printf("Error in _mm_mxcalloc : cannot allocate zero bytes.");    
    return ((double *) 0);
  }
 /* Assume malloc'd pointer is aligned at least to sizeof (void*).
    If necessary, add another sizeof (void*) to store the value
    returned by malloc. Effectively this enforces a minimum alignment
    of sizeof double. */
    if (align < 2 * sizeof (double *))
      align = 2 * sizeof (double *);

  calloc_ptr = calloc (size + align/typesize,typesize);
  if (!calloc_ptr)
    return ((double *) 0);

  /* Align  We have at least sizeof (void *) space below malloc'd ptr. */
  aligned_ptr = (double *) (((size_t) calloc_ptr + align)
                          & ~((size_t) (align) - 1));

  /* Store the original pointer just before p.  */
  ((double **) aligned_ptr) [-1] = calloc_ptr;

  return aligned_ptr;
}

/*__inline__ is gcc-specific, turn it off if on Windows*/
#ifndef _WIN32
__inline__ 
#endif
void
_mm_mxFree (double * aligned_ptr)
{
  if (aligned_ptr)
    free (((double **) aligned_ptr) [-1]);
}

#endif /* _MM_MEXMALLOC_H_INCLUDED */


