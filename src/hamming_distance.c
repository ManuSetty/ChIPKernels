
#include "Biostrings_interface.h"

/************************************************************************
 Adding code for y (Hamming distance)

*******************************************************/
static int nmismatch_at_Pshift_fixedPfixedS(const cachedCharSeq *P,
				     const cachedCharSeq *S, int Pshift, int max_nmis)
{
	int nmis, i, j;
	const char *p, *s;

	nmis = 0;
	for (i = 0, j = Pshift, p = P->seq, s = S->seq + Pshift;
	     i < P->length;
	     i++, j++, p++, s++)
	{
		if (j >= 0 && j < S->length && *p == *s)
			continue;
		if (nmis++ >= max_nmis)
			break;
	}
	return nmis;
}


static int nmismatch_at_Pshift_fixedPnonfixedS(const cachedCharSeq *P,
    const cachedCharSeq *S, int Pshift, int max_nmis)
{
  int nmis, i, j;
  const char *p, *s;

  nmis = 0;
  for (i = 0, j = Pshift, p = P->seq, s = S->seq + Pshift;
       i < P->length;
       i++, j++, p++, s++)
  {
    if (j >= 0 && j < S->length && ((*p) & ~(*s)) == 0)
      continue;
    if (nmis++ >= max_nmis)
      break;
  }
  return nmis;
}



SEXP XStringSet_dist_hamming_xy (SEXP x, SEXP y, SEXP max_nmis, SEXP fixed)
{

  static int block_size = 250000;

  cachedCharSeq x_j, y_i;
  cachedXStringSet X, Y;
  int X_length, Y_length, i, j, val, index, prev_index=-1;
  SEXP ans;
  int *ans_elt;
  int fix;

  // Extract fixed
  fix = LOGICAL (fixed)[0];
  
  // Allocate first chunk of memory  
  int max_nmis0 = INTEGER (max_nmis)[0];
  int *accumulator = (int *)Calloc (block_size, int);
  int acc_count = 0;

  X = cache_XStringSet(x);
  Y = cache_XStringSet(y);
  
  X_length = get_cachedXStringSet_length(&X);
  Y_length = get_cachedXStringSet_length(&Y);
  x_j = get_cachedXStringSet_elt(&X, 0);

  for (i = 0; i < Y_length; i++) {
    
    y_i = get_cachedXStringSet_elt(&Y, i);
    for (j = 0; j < X_length; j++) {
      
      // Reallocate memory if necessary
      if (acc_count % block_size >= block_size - 3) {
	index = (int)ceil (acc_count/block_size);
	if (index != prev_index) {
	  accumulator = (int *)Realloc(accumulator, (index + 2) * block_size, int);
	  prev_index = index;
	}

      }

      if (y_i.length != x_j.length) {
	Free (accumulator);
	error("Hamming distance requires equal length strings");
      }

      // Determine score
      x_j = get_cachedXStringSet_elt(&X, j);
      if (fix)
        val = nmismatch_at_Pshift_fixedPfixedS (&y_i, &x_j, 0, max_nmis0);
      else
        val = nmismatch_at_Pshift_fixedPnonfixedS (&x_j, &y_i, 0, max_nmis0);
      
      // Accumulate only if distance > 0
      if (val <= max_nmis0 ) {
	accumulator[acc_count++] = j + 1;
	accumulator[acc_count++] = i + 1;
	accumulator[acc_count++] = val;
      }
    }
    
  }

  // Prepare return value
  PROTECT (ans = allocVector (INTSXP, acc_count));
  ans_elt = INTEGER (ans);
  for (i = 0;i < acc_count; i++)
    ans_elt[i] = accumulator[i];

  UNPROTECT(1);
  Free (accumulator);
  return ans;
}
