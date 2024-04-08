#ifndef CRTS_HPP_INCLUDED
#define CRTS_HPP_INCLUDED

#if defined(__sw_host__) || defined(__sw_slave__)
#ifdef __cplusplus
extern "C" {
#include <crts.h>
}
#else
#include <crts.h>
#endif
#endif

#ifdef __sw_slave__
#include <simd.h>
#endif

#endif