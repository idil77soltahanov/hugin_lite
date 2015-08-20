#ifdef HUGIN_LITE
#ifndef _FAST_LIBPANO13_H
#define _FAST_LIBPANO13_H

// "size_t" (used by pano13/panorama.h) is defined in <string>, so we need to include it
#include <string>

extern "C" {
#include <pano13/panorama.h>
#include <pano13/filter.h>
}

namespace libhl {

void hl_SetMakeParams(struct fDesc *stack, struct MakeParams *mp, Image *im , Image *pn, int color);
int hl_sphere_tp_erect(double x_dest, double y_dest, double* x_src, double* y_src, void* params);
int hl_erect_stereographic(double x_dest, double y_dest, double* lon, double* lat, void* params);
int hl_persp_sphere(double x_dest, double y_dest, double* x_src, double* y_src, void* params);

} // namespace

#endif // _H
#endif // HUGIN_LITE
