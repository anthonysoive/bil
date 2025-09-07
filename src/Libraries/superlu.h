#include "BilConfig.h"

#ifdef HAVE_SUPERLU
#include "superluseq.h.in"
#endif

#ifdef HAVE_SUPERLUMT
#include "superlumt.h.in"
#endif

#ifdef HAVE_SUPERLUDIST
#include "superludist.h.in"
#endif
