
#include <cstdlib>
#include <cstdarg>
#include <iostream>
#include <getopt.h>

#include "RunningStat.hpp"
#include "discbins.hpp"
#include "kstring.h"
#include "sam.h"
#include "histogram.hpp"
#include "gatlib.hpp"

const char COL_RESET[] = "\x1b[0m";

// Foreground colors are in form of 3x, bacground are 4x
const char RED[]     = "\x1b[31m\x1b[1m";
const char GREEN[]   = "\x1b[32m\x1b[1m";
const char YELLOW[]  = "\x1b[33m\x1b[1m";
const char BLUE[]    = "\x1b[34m\x1b[1m";
const char MAGENTA[] = "\x1b[35m\x1b[1m";
const char CYAN[]    = "\x1b[36m\x1b[1m";
