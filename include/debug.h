
#ifndef VECORIZATION_PARALLELOGRAMS_DEBUG_H
#define VECORIZATION_PARALLELOGRAMS_DEBUG_H

#define SAYF

#ifdef SAYF

#define SAY(...) fprintf(stdout, __VA_ARGS__)

#else

#define SAY(X...)

#endif


#endif //VECORIZATION_PARALLELOGRAMS_DEBUG_H
