//
// Created by maratdin7 on 26.04.2020.
//

#ifndef VECORIZATION_PARALLELOGRAMS_DEBUG_H
#define VECORIZATION_PARALLELOGRAMS_DEBUG_H

#define DEBUG

#ifdef DEBUG

#define SAY(...) fprintf(stderr, __VA_ARGS__)

#else

#define SAY(X...)

#endif
#endif //VECORIZATION_PARALLELOGRAMS_DEBUG_H
