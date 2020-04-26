//
// Created by maratdin7 on 26.04.2020.
//

#ifndef VECORIZATION_PARALLELOGRAMS_DEBUG_H
#define VECORIZATION_PARALLELOGRAMS_DEBUG_H

#define DEBUG

#ifdef DEBUG

#define SAY(x...) fprintf(stderr, x)

#else

#define SAY(X...)

#endif
#endif //VECORIZATION_PARALLELOGRAMS_DEBUG_H
