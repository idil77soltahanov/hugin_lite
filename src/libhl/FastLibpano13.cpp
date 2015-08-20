#ifdef HUGIN_LITE

#include "FastLibpano13.h"
#include "FastMath.h"

#define distanceparam (*((double*)params))

namespace libhl {

/* Replace some transform functions with faster version */
void hl_SetMakeParams(struct fDesc *stack, struct MakeParams *mp, Image *im , Image *pn, int color)
{
    int i;
    double a,b; // field of view in rad
    double tx, ty, tpara; // temporary variables
    /* Joost Nieuwenhuijse, 3 feb 2005: Fix for cropping bug
    If a script containing the 'C' crop parameter was stitched by PTStitcher,
    it would fail if the cropping area is partially outside the source image.

    For 'inside' cropping, PTStitcher apparently pre-crops the images, such that
    *im contains the cropped area of the source image.
    For 'outside' cropping, PTStitcher apparently does nothing. The cropping area
    is stored in im->selection, and im->cp.cutFrame is set, but this information
    was not used at all.

    This is fixed here: All processing is now done based on the width&height of the
    cropped area (instead of the width&height of the image). And an additional horizontal
    and vertical offset are added to compensate for the shift of the center of the
    crop area relative to the center of the image.
    */
    int image_selection_width=im->width;
    int image_selection_height=im->height;
    mp->im = im;
    mp->pn = pn;
    if(im->cP.horizontal) {
        mp->horizontal = im->cP.horizontal_params[color];
    } else {
        mp->horizontal = 0;
    }
    if(im->cP.vertical) {
        mp->vertical = im->cP.vertical_params[color];
    } else {
        mp->vertical = 0;
    }
    if((im->selection.left != 0) || (im->selection.top != 0) || (im->selection.bottom != 0) || (im->selection.right != 0)) {
        if(im->cP.cutFrame) {
          image_selection_width  = im->selection.right  - im->selection.left;
          image_selection_height = im->selection.bottom - im->selection.top;
          mp->horizontal += (im->selection.right  + im->selection.left - (int32_t)im->width)/2.0;
          mp->vertical += (im->selection.bottom + im->selection.top  - (int32_t)im->height)/2.0;
        }
    }

    a = DEG_TO_RAD(im->hfov);    // field of view in rad
    b = DEG_TO_RAD(pn->hfov);

    SetMatrix(- DEG_TO_RAD(im->pitch),
             0.0,
             - DEG_TO_RAD(im->roll),
             mp->mt,
             0);

    /* Pablo d'Angelo, April 2006.
    * Added more output projection types. Broke mp->distance and mp->scale factor calculation
    * into separate parts, making it easier to add new projection types
    */
    // calculate distance
    switch (pn->format)
    {
        case _rectilinear:
            mp->distance        = (double) pn->width / (2.0 * tan(b/2.0));
            break;
        case _equirectangular:
        case _fisheye_ff:
        case _fisheye_circ:
        case _panorama:
        case _lambert:
        case _mercator:
        case _millercylindrical:
        case _sinusoidal:
        case _mirror:
            // horizontal pixels per degree
            mp->distance        = ((double) pn->width) / b;
            break;
        case _panini:
            tpara = 1;
            panini_erect(b/2.0, 0.0, &tx, &ty, & tpara);
            mp->distance = pn->width/(2.0*tx);
            break;
        case _equipanini:
            tpara = 1;
            equipanini_erect(b/2.0, 0.0, &tx, &ty, & tpara);
            mp->distance = pn->width/(2.0*tx);
            break;
        case _panini_general:
            // call setup_panini_general() to set distanceparam
            pn->precomputedCount = 0;	// clear old settings
            setup_panini_general(mp);
            // should abort now if it returns NULL
            break;
        case _architectural:
            tpara = 1;
            arch_erect(b/2.0, 0.0, &tx, &ty, & tpara);
            mp->distance = pn->width/(2.0*tx);
            break;
        case _lambertazimuthal:
            tpara = 1;
            lambertazimuthal_erect(b/2.0, 0.0, &tx, &ty, & tpara);
            mp->distance = pn->width/(2.0*tx);
            break;
        case _hammer:
            tpara = 1;
            hammer_erect(b/2.0, 0.0, &tx, &ty, & tpara);
            mp->distance = pn->width/(2.0*tx);
            break;
        case _stereographic:
            tpara = 1;
            stereographic_erect(b/2.0, 0.0, &tx, &ty, & tpara);
            mp->distance = pn->width/(2.0*tx);
            break;
        case _trans_mercator:
            tpara = 1;
            transmercator_erect(b/2.0, 0.0, &tx, &ty, &tpara);
            mp->distance = pn->width/(2.0*tx);
            break;
        case _albersequalareaconic:
            mp->distance = 1.0;
            //albersequalareaconic_erect(1.924913116, -PI/2.0, &tx, &ty, mp); //b/2.0
            albersequalareaconic_distance(&tx, mp);
            mp->distance = pn->width/(2.0*tx);
            break;
        case _equisolid:
            mp->distance  = (double) pn->width / (4.0 * sin(b/4.0));
            break;
        case _orthographic:
            mp->distance  = (double) pn->width / (2.0 * sin(b/2.0));
            break;
        case _thoby:
            mp->distance  = (double) pn->width / (2.0 * THOBY_K1_PARM * sin(b * THOBY_K2_PARM/2.0));
            break;
        case _biplane:
            biplane_distance(pn->width,b,mp);
            break;
        case _triplane:
            triplane_distance(pn->width,b,mp);
            break;
        default:
            // unknown
            PrintError ("SetMakeParams: Unsupported panorama projection");
            // no way to report an error back to the caller...
            mp->distance = 1;
            break;
    }

    // calculate final scaling factor, that reverses the mp->distance
    // scaling and applies the required output scaling factor
    // printf("im format %d\n", im->format);
    switch (im->format)
    {
        case _rectilinear:
            // calculate distance for this projection
            mp->scale[0] = (double) image_selection_width / (2.0 * tan(a/2.0)) / mp->distance;
            break;
        case _equirectangular:
        case _panorama:
        case _fisheye_ff:
        case _fisheye_circ:
        case _mercator:
        case _sinusoidal:
            mp->scale[0] = ((double) image_selection_width) / a / mp->distance;
            break;
        case _equisolid:
        case _mirror:
            mp->scale[0] = (double) image_selection_width / (4.0 * sin(a/4.0)) / mp->distance;
            break;
        case _orthographic:
            {
                //generate monotonic scale function to help optimizer
                int t=(int)ceil((a-PI)/(2.0*PI));
                mp->scale[0] = (double) image_selection_width / (2.0 * (2 * t + pow(-1.0, t) * sin(a/2.0))) / mp->distance;
            };
            break;
        case _thoby:
            mp->scale[0] = (double) image_selection_width / (2.0 * THOBY_K1_PARM * sin(a * THOBY_K2_PARM /2.0)) / mp->distance;
            break;
        case _stereographic:
            mp->scale[0] = (double) image_selection_width / (4.0 * tan(a/4.0)) / mp->distance;
            break;
        default:
            PrintError ("SetMakeParams: Unsupported input image projection");
            // no way to report an error back to the caller...
            mp->scale[0] = 1;
            break;
    }
    mp->scale[1] = mp->scale[0];

    // printf("new params: mp->distance: %lf, mp->scale: %lf\n\n", mp->distance, mp->scale[0]);

    mp->shear[0] = im->cP.shear_x / image_selection_height;
    mp->shear[1] = im->cP.shear_y / image_selection_width;
    mp->rot[0] = mp->distance * PI; // 180 in screenpoints
    mp->rot[1] = -im->yaw *  mp->distance * PI / 180.0; // rotation angle in screenpoints

    mp->tilt[0] = DEG_TO_RAD(im->cP.tilt_x);
    mp->tilt[1] = DEG_TO_RAD(im->cP.tilt_y);
    mp->tilt[2] = DEG_TO_RAD(im->cP.tilt_z);
    mp->tilt[3] = im->cP.tilt_scale;

    mp->trans[0] = im->cP.trans_x;
    mp->trans[1] = im->cP.trans_y;
    mp->trans[2] = im->cP.trans_z;
    mp->trans[3] = DEG_TO_RAD(im->cP.trans_yaw);
    mp->trans[4] = DEG_TO_RAD(im->cP.trans_pitch);

    mp->test[0] = im->cP.test_p0;
    mp->test[1] = im->cP.test_p1;
    mp->test[2] = im->cP.test_p2;
    mp->test[3] = im->cP.test_p3;

    // panoAdjustPrintMakeParams("SetmakeParms", mp, im);

    mp->perspect[0] = (void*)(mp->mt);
    mp->perspect[1] = (void*)&(mp->distance);

    for(i=0; i<4; i++)
        mp->rad[i]  = im->cP.radial_params[color][i];
    mp->rad[5] = im->cP.radial_params[color][4];

    if((im->cP.correction_mode & 3) == correction_mode_radial)
        mp->rad[4] = ((double)(image_selection_width < image_selection_height ? image_selection_width : image_selection_height)) / 2.0;
    else
        mp->rad[4] = ((double) image_selection_height) / 2.0;

    // Joost: removed, see above
    // mp->horizontal  = im->cP.horizontal_params[color];
    // mp->vertical  = im->cP.vertical_params[color];

    i = 0;


    // Building the stack
    //
    // - Convert from panorama projection to equirectangular
    // - Rotate horizontally
    // - Convert to spherical from equirectangular
    // - Apply perspective correction (pitch and roll) in spherical coordinates
    // - Convert to image format (rectilinear, pano, equirectangular)
    // - Scale output image
    // - Do radial correction
    // - Do tilt
    // - Do vertical shift
    // - Do horizontal shift
    // - Do shear


    //////////////////////////////////////////////////////////////////////
    // Convert from output projection to spherical coordinates
    //
    if(pn->format == _rectilinear) { // rectilinear panorama
        SetDesc(stack[i], erect_rect, &(mp->distance)); i++; // Convert rectilinear to equirect
    } else if(pn->format == _panorama) {
        SetDesc(stack[i], erect_pano, &(mp->distance)); i++; // Convert panoramic to equirect
    } else if(pn->format == _fisheye_circ || pn->format == _fisheye_ff) {
        // the sphere coordinates are actually equivalent to the equidistant fisheye projection
        SetDesc(stack[i], erect_sphere_tp, &(mp->distance)); i++; // Convert fisheye to equirect
    } else if(pn->format == _equisolid) {
        SetDesc(stack[i], sphere_tp_equisolid, &(mp->distance)); i++; // Convert fisheye equisolid to spherical
        SetDesc(stack[i], erect_sphere_tp, &(mp->distance)); i++; // Convert spherical to equirect
    } else if(pn->format == _mirror) {
        SetDesc(stack[i], sphere_cp_mirror, &(mp->distance)); i++; // Convert mirror to spherical
        SetDesc(stack[i], erect_sphere_cp, &(mp->distance)); i++; // Convert spherical to equirect
    } else if(pn->format == _orthographic) {
        SetDesc(stack[i], sphere_tp_orthographic, &(mp->distance)); i++; // Convert fisheye orthographic to spherical
        SetDesc(stack[i], erect_sphere_tp, &(mp->distance)); i++; // Convert spherical to equirect
    } else if(pn->format == _thoby) {
        SetDesc(stack[i], sphere_tp_thoby, &(mp->distance)); i++; // Convert thoby to spherical
        SetDesc(stack[i], erect_sphere_tp, &(mp->distance)); i++; // Convert spherical to equirect
    } else if(pn->format == _mercator) {
        SetDesc(stack[i], erect_mercator, &(mp->distance)); i++; // Convert mercator to equirect
    } else if(pn->format == _millercylindrical) {
        SetDesc(stack[i], erect_millercylindrical, &(mp->distance)); i++; // Convert miller to equirect
    } else if(pn->format == _panini) {
        SetDesc(stack[i], erect_panini, &(mp->distance)); i++; // Convert panini to sphere
    } else if(pn->format == _equipanini) {
        SetDesc(stack[i], erect_equipanini, &(mp->distance)); i++; // Convert equipanini to sphere
    } else if(pn->format == _panini_general) {
        SetDesc(stack[i], erect_panini_general, mp); i++; // Convert general panini to sphere
    } else if(pn->format == _architectural) {
        SetDesc(stack[i], erect_arch, &(mp->distance)); i++; // Convert arch to sphere
    } else if(pn->format == _lambert) {
        SetDesc(stack[i], erect_lambert, &(mp->distance)); i++; // Convert lambert to equirect
    } else if(pn->format == _lambertazimuthal) {
        SetDesc(stack[i], erect_lambertazimuthal, &(mp->distance)); i++; // Convert lambert to equirect
    } else if(pn->format == _hammer) {
        SetDesc(stack[i], erect_hammer, &(mp->distance)); i++; // Convert hammer to equirect
    } else if(pn->format == _trans_mercator) {
        SetDesc(stack[i], erect_transmercator, &(mp->distance)); i++; // Convert transverse mercator to equirect
    } else if(pn->format == _stereographic) {
        SetDesc(stack[i], hl_erect_stereographic, &(mp->distance)); i++;  // Convert stereographic to equirect
    } else if(pn->format == _sinusoidal) {
        SetDesc(stack[i], erect_sinusoidal, &(mp->distance)); i++; // Convert sinusoidal to equirect
    } else if(pn->format == _albersequalareaconic) {
        SetDesc(stack[i], erect_albersequalareaconic, mp); i++; // Convert albersequalareaconic to equirect
    } else if(pn->format == _biplane) {
        SetDesc(stack[i], erect_biplane, mp); i++;  // Convert biplane to equirect
    } else if(pn->format == _triplane) {
        SetDesc(stack[i], erect_triplane, mp); i++;  // Convert triplane to equirect
    } else if(pn->format == _equirectangular) {
        // no conversion needed
    } else {
        PrintError("Projection type %d not supported. Assuming equirectangular", pn->format);
    }

    if (im->cP.trans) {
        SetDesc(stack[i], plane_transfer_to_camera, mp); i++;
    }

    SetDesc(stack[i], rotate_erect, mp->rot); i++; // Rotate equirect. image horizontally
    SetDesc(stack[i], hl_sphere_tp_erect, &(mp->distance)); i++; // Convert spherical image to equirect.
    SetDesc(stack[i], hl_persp_sphere, mp->perspect); i++; // Perspective Control spherical Image

    //////////////////////////////////////////////////////////////////////
    // Convert from spherical coordinates to input projection
    //
    if(im->format == _rectilinear) { // rectilinear image
        SetDesc(stack[i], rect_sphere_tp, &(mp->distance)); i++; // Convert spherical to rectilinear
    } else if(im->format == _panorama) {                                   //  pamoramic image
        SetDesc(stack[i], pano_sphere_tp, &(mp->distance)); i++; // Convert spherical to pano
    } else if(im->format == _equirectangular) { //  equirectangular image
        SetDesc(stack[i], erect_sphere_tp, &(mp->distance)); i++; // Convert spherical to equirect
    } else if (im->format == _fisheye_circ || im->format == _fisheye_ff) {
        // no conversion needed. It is already in spherical coordinates
    } else if (im->format == _mirror) {
        SetDesc(stack[i], mirror_sphere_tp, &(mp->distance)); i++; // Convert spherical to mirror
    } else if (im->format == _stereographic) {
        SetDesc(stack[i], erect_sphere_tp, &(mp->distance)); i++; // Convert spherical to equirectangular
        SetDesc(stack[i], stereographic_erect, &(mp->distance)); i++; // Convert equirectangular to stereographic
    } else if (im->format == _orthographic) {
        SetDesc(stack[i], orthographic_sphere_tp, &(mp->distance)); i++; // Convert spherical to orthographic
    } else if (im->format == _thoby) {
        SetDesc(stack[i], thoby_sphere_tp,           &(mp->distance)); i++; // Convert spherical to thoby
    } else if (im->format == _equisolid) {
        SetDesc(stack[i], erect_sphere_tp, &(mp->distance)); i++; // Convert spherical to equirectangular
        SetDesc(stack[i], lambertazimuthal_erect, &(mp->distance)); i++; // Convert equirectangular to stereographic
    } else {
        PrintError("Invalid input projection %d. Assumed fisheye.", im->format);
    }


    SetDesc(stack[i], resize, mp->scale); i++; // Scale image

    //////////////////////////////////////////////////////////////////////
    // Apply lens corrections
    //

    if(im->cP.radial) {
        switch(im->cP.correction_mode & 3) {
            case correction_mode_radial: SetDesc(stack[i],radial,mp->rad); i++; break;
            case correction_mode_vertical: SetDesc(stack[i],vertical,mp->rad); i++; break;
            case correction_mode_deregister: SetDesc(stack[i],deregister,mp->rad); i++; break;
        }
    }
    if (im->cP.tilt) {
        SetDesc(stack[i], tiltInverse, mp); i++;
    }

    if (mp->vertical != 0.0) {
        SetDesc(stack[i], vert, &(mp->vertical)); i++;
    }
    if (mp->horizontal != 0.0) {
        SetDesc(stack[i], horiz, &(mp->horizontal)); i++;
    }
    if(im->cP.shear) {
        SetDesc(stack[i], shear, mp->shear); i++;
    }

    stack[i].func = (trfn)NULL;
}

/* Replace atan2 with (approximate) faster version */
int hl_sphere_tp_erect(double x_dest, double y_dest, double* x_src, double* y_src, void* params)
{
    // params: double distanceparam

    register double phi, theta, r,s;
    double v[3];

    phi = x_dest / distanceparam;
    theta = - y_dest / distanceparam + PI / 2;
    if(theta < 0) {
        theta = - theta;
        phi += PI;
    }
    if(theta > PI) {
        theta = PI - (theta - PI);
        phi += PI;
    }


    s = sin(theta);
    v[0] = s * sin(phi); // y' -> x
    v[1] = cos(theta); // z' -> y

    r = sqrt(v[1] * v[1] + v[0] * v[0]);

    theta = distanceparam * hl_atan2(r , s * cos(phi));

    *x_src = theta * v[0] / r;
    *y_src = theta * v[1] / r;

    return 1;
}

/* Replace atan2 with (approximate) faster version */
int hl_erect_stereographic(double x_dest, double y_dest, double* lon, double* lat, void* params)
{
    double rh; /* height above sphere */
    double c; /* angle */
    double sinc,cosc; /* sin of c and cos of c */

    /* Inverse equations
     -----------------*/
    double x = x_dest / distanceparam;
    double y = y_dest / distanceparam;
    rh = sqrt(x * x + y * y);
    c = 2.0 * atan(rh / (2.0 * 1));
    sinc = sin(c);
    cosc = cos(c);
    *lon = 0;
    if (fabs(rh) <= EPSLN) {
        *lat = 0;
        return 0;
    } else {
        *lat = asin((y * sinc) / rh) * distanceparam;

        if ((fabs(cosc) < EPSLN) && (fabs(x) < EPSLN)) {
            return 0;
        } else {
            *lon = hl_atan2((x * sinc), (cosc * rh)) * distanceparam;
        }
    }

    return 1;
}

/* Replace atan2 with (approximate) faster version */
int hl_persp_sphere(double x_dest, double y_dest, double* x_src, double* y_src, void* params)
{
    // params : double Matrix[3][3], double distanceparam

    register double theta,s,r;
    double v[3];

    r = sqrt(x_dest * x_dest + y_dest * y_dest);
    theta = r / *((double*) ((void**)params)[1]);
    if(r == 0.0) {
        s = 0.0;
    } else {
        s = sin(theta) / r;
    }

    v[0] = s * x_dest;
    v[1] = s * y_dest;
    v[2] = cos(theta);

    matrix_inv_mult((double(*)[3]) ((void**)params)[0], v);

    r = sqrt(v[0] * v[0] + v[1] * v[1]);
    if(r == 0.0) {
        theta = 0.0;
    } else {
        theta = *((double*) ((void**)params)[1]) * hl_atan2(r, v[2]) / r;
    }
    *x_src = theta * v[0];
    *y_src = theta * v[1];

    return 1;
}

} // namespace

#endif
