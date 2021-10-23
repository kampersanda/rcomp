/**
 * @file defs.hpp
 */
#pragma once

// #define ENABLE_DEBUG_PRINT
// #define ENABLE_STAT_MONITOR

#ifdef ENABLE_DEBUG_PRINT
#define DEBUG_PRINT(x) x
#else
#define DEBUG_PRINT(x)
#endif

#ifdef ENABLE_STAT_MONITOR
#define STAT_MONITOR(x) x
#else
#define STAT_MONITOR(x)
#endif
