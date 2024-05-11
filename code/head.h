#pragma once

#ifdef DEBUG_INFO 
this_is_deprecated
#endif
#ifdef DEBUG_TRACE
this_is_deprecated
#endif


#if defined(WIN32)
#elif defined(__CYGWIN__) // cygwin
#include <sys/time.h>
#else //linux
#include <sys/mman.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#endif


#include <iostream>
#include <set>
#include <list>
#include <sstream>
#include <cmath>
#include <queue>
#include <fstream>
#include <string>
#include <cstdio>
#include <functional>
#include <algorithm>
#include <climits>
#include <cstring>
#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <map>
#include <deque>
#include <chrono>
#include <ctime>
#include <ratio>
#include <sys/mman.h>  //included only in Linux, thus put it under the condition of Linux
//#include <sys/stat.h>
#include <fcntl.h>
#include <unordered_set>
#include <set>
//#include <assert.h>

///  Self-defined head files
#include "Algorithm.h"
#include "Argument.h"
#include "CommonFunc.h"
#include "graph.h"
// #include "heap.h"
// #include "Inf_Graph.h"
// #include "Memory.h"
// #include "MemoryUsage.h"
#include "RRcollection.h"
#include "serialize.h"
#include "FileCtrl.h"
// #include "Timer.h"
#include "../dSFMT/dSFMT.h"