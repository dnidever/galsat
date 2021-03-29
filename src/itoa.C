/*
 * COPYRIGHT:   See COPYING in the top level directory
 * PROJECT:     ReactOS system libraries
 * FILE:        lib/crtdll/stdlib/itoa.c
 * PURPOSE:     converts a integer to ascii
 * PROGRAMER:   
 * UPDATE HISTORY:
 *              1995: Created
 *              1998: Added ltoa Boudewijn Dekker
 */
/* Copyright (C) 1995 DJ Delorie, see COPYING.DJ for details */
#include <errno.h>
#include <stdlib.h>
//#include <file.h>

char *
itoa(int value, char *string)
{
 
  int radix = 10;
  char tmp[33];
  char *tp = tmp;
  int i;
  unsigned v;
  int sign;
  char *sp;

  //if (radix > 36 || radix <= 1)
  //{
  // __set_errno(EDOM);
  //  return 0;
  //}

  sign = (radix == 10 && value < 0);
  if (sign)
    v = -value;
  else
    v = (unsigned)value;
  while (v || tp == tmp)
  {
    i = v % radix;
    v = v / radix;
    if (i < 10)
      *tp++ = i+'0';
    else
      *tp++ = i + 'a' - 10;
  }

  if (string == 0)
    string = (char *)malloc((tp-tmp)+sign+1);
  sp = string;

  if (sign)
    *sp++ = '-';
  while (tp > tmp)
    *sp++ = *--tp;
  *sp = 0;
  return string;
}


#include <windows.h>
#include <stdlib.h>


void* malloc(size_t _size)
{
   return(HeapAlloc(GetProcessHeap(),0,_size));
}

void free(void* _ptr)
{
   HeapFree(GetProcessHeap(),0,_ptr);
}

void* calloc(size_t _nmemb, size_t _size)
{
   return(HeapAlloc(GetProcessHeap(),HEAP_ZERO_MEMORY, _nmemb*_size));
}

void* realloc(void* _ptr, size_t _size)
{
   return(HeapReAlloc(GetProcessHeap(),0,_ptr,_size));
}
#undef alloca
void *alloca(size_t s)
{
	register unsigned int as = s;

// alloca(0) should not return the stack pointer
	if ( s == 0 )
		return NULL;

	
	if ( (s & 0xfffffffc)  != 0 )
		as += 4;
		
	as &= 0xfffffffc;
	
	__asm__ __volatile__(
	"mov %0, %%edx  	\n"
//	"popl %%ebp		\n"
	"leave			\n"
	"popl  %%ecx		\n"
        "subl  %%edx, %%esp	\n"
        "movl  %%esp, %%eax	\n"
        "addl  $20, %%eax        \n"//4 bytes + 16 bytes = arguments
        "push  %%ecx		\n"
        "ret			\n"
        :
        :"g" ( as)
        );
        
       
        return NULL;
}

void *_alloca(size_t s)
{
	register unsigned int as = s;

// alloca(0) should not return the stack pointer
	if ( s == 0 )
		return NULL;

	
	if ( (s & 0xfffffffc)  != 0 )
		as += 4;
		
	as &= 0xfffffffc;
	
	__asm__ __volatile__(
	"mov %0, %%edx  	\n"
//	"popl %%ebp		\n"
	"leave			\n"
	"popl  %%ecx		\n"
        "subl  %%edx, %%esp	\n"
        "movl  %%esp, %%eax	\n"
        "addl  $20, %%eax        \n"//4 bytes + 16 bytes = arguments
        "push  %%ecx		\n"
        "ret			\n"
        :
        :"g" ( as)
        );
        
       
        return NULL;
}
