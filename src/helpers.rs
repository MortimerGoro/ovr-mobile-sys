use std::mem;
use super::*;
use super::ovrFrameInit::*;
use super::ovrGraphicsAPI::*;
use super::ovrStructureType::*;

//-----------------------------------------------------------------
// Matrix helper functions.
//-----------------------------------------------------------------

// Use left-multiplication to accumulate transformations.
pub fn ovrMatrix4f_Multiply(a: &ovrMatrix4f, b: &ovrMatrix4f) -> ovrMatrix4f {
    let mut out: ovrMatrix4f = unsafe { mem::uninitialized( )};

    out.M[0][0] = a.M[0][0] * b.M[0][0] + a.M[0][1] * b.M[1][0] + a.M[0][2] * b.M[2][0] + a.M[0][3] * b.M[3][0];
    out.M[1][0] = a.M[1][0] * b.M[0][0] + a.M[1][1] * b.M[1][0] + a.M[1][2] * b.M[2][0] + a.M[1][3] * b.M[3][0];
    out.M[2][0] = a.M[2][0] * b.M[0][0] + a.M[2][1] * b.M[1][0] + a.M[2][2] * b.M[2][0] + a.M[2][3] * b.M[3][0];
    out.M[3][0] = a.M[3][0] * b.M[0][0] + a.M[3][1] * b.M[1][0] + a.M[3][2] * b.M[2][0] + a.M[3][3] * b.M[3][0];

    out.M[0][1] = a.M[0][0] * b.M[0][1] + a.M[0][1] * b.M[1][1] + a.M[0][2] * b.M[2][1] + a.M[0][3] * b.M[3][1];
    out.M[1][1] = a.M[1][0] * b.M[0][1] + a.M[1][1] * b.M[1][1] + a.M[1][2] * b.M[2][1] + a.M[1][3] * b.M[3][1];
    out.M[2][1] = a.M[2][0] * b.M[0][1] + a.M[2][1] * b.M[1][1] + a.M[2][2] * b.M[2][1] + a.M[2][3] * b.M[3][1];
    out.M[3][1] = a.M[3][0] * b.M[0][1] + a.M[3][1] * b.M[1][1] + a.M[3][2] * b.M[2][1] + a.M[3][3] * b.M[3][1];

    out.M[0][2] = a.M[0][0] * b.M[0][2] + a.M[0][1] * b.M[1][2] + a.M[0][2] * b.M[2][2] + a.M[0][3] * b.M[3][2];
    out.M[1][2] = a.M[1][0] * b.M[0][2] + a.M[1][1] * b.M[1][2] + a.M[1][2] * b.M[2][2] + a.M[1][3] * b.M[3][2];
    out.M[2][2] = a.M[2][0] * b.M[0][2] + a.M[2][1] * b.M[1][2] + a.M[2][2] * b.M[2][2] + a.M[2][3] * b.M[3][2];
    out.M[3][2] = a.M[3][0] * b.M[0][2] + a.M[3][1] * b.M[1][2] + a.M[3][2] * b.M[2][2] + a.M[3][3] * b.M[3][2];

    out.M[0][3] = a.M[0][0] * b.M[0][3] + a.M[0][1] * b.M[1][3] + a.M[0][2] * b.M[2][3] + a.M[0][3] * b.M[3][3];
    out.M[1][3] = a.M[1][0] * b.M[0][3] + a.M[1][1] * b.M[1][3] + a.M[1][2] * b.M[2][3] + a.M[1][3] * b.M[3][3];
    out.M[2][3] = a.M[2][0] * b.M[0][3] + a.M[2][1] * b.M[1][3] + a.M[2][2] * b.M[2][3] + a.M[2][3] * b.M[3][3];
    out.M[3][3] = a.M[3][0] * b.M[0][3] + a.M[3][1] * b.M[1][3] + a.M[3][2] * b.M[2][3] + a.M[3][3] * b.M[3][3];

    out
}

// Returns the transpose of a 4x4 matrix.
pub fn ovrMatrix4f_Transpose(a: &ovrMatrix4f) -> ovrMatrix4f {
    let mut out: ovrMatrix4f = unsafe { mem::uninitialized( )};

    out.M[0][0] = a.M[0][0]; out.M[0][1] = a.M[1][0]; out.M[0][2] = a.M[2][0]; out.M[0][3] = a.M[3][0];
    out.M[1][0] = a.M[0][1]; out.M[1][1] = a.M[1][1]; out.M[1][2] = a.M[2][1]; out.M[1][3] = a.M[3][1];
    out.M[2][0] = a.M[0][2]; out.M[2][1] = a.M[1][2]; out.M[2][2] = a.M[2][2]; out.M[2][3] = a.M[3][2];
    out.M[3][0] = a.M[0][3]; out.M[3][1] = a.M[1][3]; out.M[3][2] = a.M[2][3]; out.M[3][3] = a.M[3][3];

    out
}

// Returns a 3x3 minor of a 4x4 matrix.
pub fn ovrMatrix4f_Minor(m: &ovrMatrix4f, r0:usize, r1:usize, r2:usize, c0:usize, c1:usize, c2:usize) -> f32 {
    m.M[r0][c0] * ( m.M[r1][c1] * m.M[r2][c2] - m.M[r2][c1] * m.M[r1][c2] ) -
    m.M[r0][c1] * ( m.M[r1][c0] * m.M[r2][c2] - m.M[r2][c0] * m.M[r1][c2] ) +
    m.M[r0][c2] * ( m.M[r1][c0] * m.M[r2][c1] - m.M[r2][c0] * m.M[r1][c1] )
}
 
// Returns the inverse of a 4x4 matrix.
pub fn ovrMatrix4f_Inverse(m: &ovrMatrix4f) -> ovrMatrix4f
{
    let rcp_det = 1.0 / (m.M[0][0] * ovrMatrix4f_Minor( m, 1, 2, 3, 1, 2, 3 ) -
                         m.M[0][1] * ovrMatrix4f_Minor( m, 1, 2, 3, 0, 2, 3 ) +
                         m.M[0][2] * ovrMatrix4f_Minor( m, 1, 2, 3, 0, 1, 3 ) -
                         m.M[0][3] * ovrMatrix4f_Minor( m, 1, 2, 3, 0, 1, 2 ) );
    let mut out: ovrMatrix4f = unsafe { mem::uninitialized( )};
    out.M[0][0] =  ovrMatrix4f_Minor( m, 1, 2, 3, 1, 2, 3 ) * rcp_det;
    out.M[0][1] = -ovrMatrix4f_Minor( m, 0, 2, 3, 1, 2, 3 ) * rcp_det;
    out.M[0][2] =  ovrMatrix4f_Minor( m, 0, 1, 3, 1, 2, 3 ) * rcp_det;
    out.M[0][3] = -ovrMatrix4f_Minor( m, 0, 1, 2, 1, 2, 3 ) * rcp_det;
    out.M[1][0] = -ovrMatrix4f_Minor( m, 1, 2, 3, 0, 2, 3 ) * rcp_det;
    out.M[1][1] =  ovrMatrix4f_Minor( m, 0, 2, 3, 0, 2, 3 ) * rcp_det;
    out.M[1][2] = -ovrMatrix4f_Minor( m, 0, 1, 3, 0, 2, 3 ) * rcp_det;
    out.M[1][3] =  ovrMatrix4f_Minor( m, 0, 1, 2, 0, 2, 3 ) * rcp_det;
    out.M[2][0] =  ovrMatrix4f_Minor( m, 1, 2, 3, 0, 1, 3 ) * rcp_det;
    out.M[2][1] = -ovrMatrix4f_Minor( m, 0, 2, 3, 0, 1, 3 ) * rcp_det;
    out.M[2][2] =  ovrMatrix4f_Minor( m, 0, 1, 3, 0, 1, 3 ) * rcp_det;
    out.M[2][3] = -ovrMatrix4f_Minor( m, 0, 1, 2, 0, 1, 3 ) * rcp_det;
    out.M[3][0] = -ovrMatrix4f_Minor( m, 1, 2, 3, 0, 1, 2 ) * rcp_det;
    out.M[3][1] =  ovrMatrix4f_Minor( m, 0, 2, 3, 0, 1, 2 ) * rcp_det;
    out.M[3][2] = -ovrMatrix4f_Minor( m, 0, 1, 3, 0, 1, 2 ) * rcp_det;
    out.M[3][3] =  ovrMatrix4f_Minor( m, 0, 1, 2, 0, 1, 2 ) * rcp_det;

    out
}

// Returns a 4x4 identity matrix.
pub fn ovrMatrix4f_CreateIdentity() -> ovrMatrix4f {
    let mut out: ovrMatrix4f = unsafe { mem::uninitialized( )};

    out.M[0][0] = 1.0; out.M[0][1] = 0.0; out.M[0][2] = 0.0; out.M[0][3] = 0.0;
    out.M[1][0] = 0.0; out.M[1][1] = 1.0; out.M[1][2] = 0.0; out.M[1][3] = 0.0;
    out.M[2][0] = 0.0; out.M[2][1] = 0.0; out.M[2][2] = 1.0; out.M[2][3] = 0.0;
    out.M[3][0] = 0.0; out.M[3][1] = 0.0; out.M[3][2] = 0.0; out.M[3][3] = 1.0;

    out
}

// Returns a 4x4 homogeneous translation matrix.
pub fn ovrMatrix4f_CreateTranslation(x: f32, y: f32, z: f32) -> ovrMatrix4f {
    let mut out: ovrMatrix4f = unsafe { mem::uninitialized( )};

    out.M[0][0] = 1.0; out.M[0][1] = 0.0; out.M[0][2] = 0.0; out.M[0][3] = x;
    out.M[1][0] = 0.0; out.M[1][1] = 1.0; out.M[1][2] = 0.0; out.M[1][3] = y;
    out.M[2][0] = 0.0; out.M[2][1] = 0.0; out.M[2][2] = 1.0; out.M[2][3] = z;
    out.M[3][0] = 0.0; out.M[3][1] = 0.0; out.M[3][2] = 0.0; out.M[3][3] = 1.0;

    out
}

// Returns a 4x4 homogeneous rotation matrix.
pub fn ovrMatrix4f_CreateRotation(radiansX: f32, radiansY: f32, radiansZ: f32) -> ovrMatrix4f
{
    let sinX = radiansX.sin();
    let cosX = radiansX.cos();
    let rotationX = ovrMatrix4f {
        M: [[ 1.0,  0.0,   0.0, 0.0 ],
            [ 0.0, cosX, -sinX, 0.0 ],
            [ 0.0, sinX,  cosX, 0.0 ],
            [ 0.0,  0.0,   0.0, 1.0 ]]
    };

    let sinY = radiansY.sin();
    let cosY = radiansY.cos();
    let rotationY = ovrMatrix4f {
        M: [[  cosY, 0.0, sinY, 0.0 ],
            [   0.0, 1.0,  0.0, 0.0 ],
            [ -sinY, 0.0, cosY, 0.0 ],
            [   0.0, 0.0,  0.0, 1.0 ]]
    };

    let sinZ = radiansZ.sin();
    let cosZ = radiansZ.cos();
    let rotationZ = ovrMatrix4f {
        M: [[ cosZ, -sinZ, 0.0, 0.0 ],
            [ sinZ,  cosZ, 0.0, 0.0 ],
            [  0.0,   0.0, 1.0, 0.0 ],
            [  0.0,   0.0, 0.0, 1.0 ]]
    };

    let rotationXY = ovrMatrix4f_Multiply(&rotationY, &rotationX);

    ovrMatrix4f_Multiply(&rotationZ, &rotationXY)
}

// Returns a projection matrix based on the specified dimensions.
// The projection matrix transforms -Z=forward, +Y=up, +X=right to the appropriate clip space for the graphics API.
// The far plane is placed at infinity if far_z <= near_z.
// An infinite projection matrix is preferred for rasterization because, except for
// things *right* up against the near plane, it always provides better precision:
//        "Tightening the Precision of Perspective Rendering"
//        Paul Upchurch, Mathieu Desbrun
//        Journal of Graphics Tools, Volume 16, Issue 1, 2012
pub fn ovrMatrix4f_CreateProjection(min_x: f32,
                                    max_x: f32,
                                    min_y: f32,
                                    max_y: f32,
                                    near_z: f32,
                                    far_z: f32)
                                    -> ovrMatrix4f
{
    let width = max_x - min_x;
    let height = max_y - min_y;
    let offsetZ = near_z;    // set to zero for a [0,1] clip space

    let mut out: ovrMatrix4f = unsafe { mem::uninitialized() };
    if far_z <= near_z {
        // place the far plane at infinity
        out.M[0][0] = 2.0 * near_z / width;
        out.M[0][1] = 0.0;
        out.M[0][2] = ( max_x + min_x ) / width;
        out.M[0][3] = 0.0;

        out.M[1][0] = 0.0;
        out.M[1][1] = 2.0 * near_z / height;
        out.M[1][2] = ( max_y + min_y ) / height;
        out.M[1][3] = 0.0;

        out.M[2][0] = 0.0;
        out.M[2][1] = 0.0;
        out.M[2][2] = -1.0;
        out.M[2][3] = -( near_z + offsetZ );

        out.M[3][0] = 0.0;
        out.M[3][1] = 0.0;
        out.M[3][2] = -1.0;
        out.M[3][3] = 0.0;
    } else {
        // normal projection
        out.M[0][0] = 2.0 * near_z / width;
        out.M[0][1] = 0.0;
        out.M[0][2] = ( max_x + min_x ) / width;
        out.M[0][3] = 0.0;

        out.M[1][0] = 0.0;
        out.M[1][1] = 2.0 * near_z / height;
        out.M[1][2] = ( max_y + min_y ) / height;
        out.M[1][3] = 0.0;

        out.M[2][0] = 0.0;
        out.M[2][1] = 0.0;
        out.M[2][2] = -( far_z + offsetZ ) / ( far_z - near_z );
        out.M[2][3] = -( far_z * ( near_z + offsetZ ) ) / ( far_z - near_z );

        out.M[3][0] = 0.0;
        out.M[3][1] = 0.0;
        out.M[3][2] = -1.0;
        out.M[3][3] = 0.0;
    }
    out
}

// Returns a projection matrix based on the given FOV.
pub fn ovrMatrix4f_CreateProjectionFov(fov_degrees_x: f32,
                                       fov_degrees_y: f32,
                                       offset_x: f32,
                                       offset_y: f32,
                                       near_z: f32,
                                       far_z: f32)
                                       -> ovrMatrix4f
{
    let half_width = near_z * (fov_degrees_x * (VRAPI_PI as f32 / 180.032 * 0.532)).tan();
    let half_height = near_z * (fov_degrees_y * (VRAPI_PI as f32 / 180.032 * 0.532)).tan();

    let min_x = offset_x - half_width;
    let max_x = offset_x + half_width;

    let min_y = offset_y - half_height;
    let max_y = offset_y + half_height;

    ovrMatrix4f_CreateProjection( min_x, max_x, min_y, max_y, near_z, far_z )
}


// Returns the 4x4 rotation matrix for the given quaternion.
pub fn ovrMatrix4f_CreateFromQuaternion(q: &ovrQuatf) -> ovrMatrix4f {
    let ww = q.w * q.w;
    let xx = q.x * q.x;
    let yy = q.y * q.y;
    let zz = q.z * q.z;

    let mut out: ovrMatrix4f = unsafe { mem::uninitialized() };
    out.M[0][0] = ww + xx - yy - zz;
    out.M[0][1] = 2.0 * ( q.x * q.y - q.w * q.z );
    out.M[0][2] = 2.0 * ( q.x * q.z + q.w * q.y );
    out.M[0][3] = 0.0;

    out.M[1][0] = 2.0 * ( q.x * q.y + q.w * q.z );
    out.M[1][1] = ww - xx + yy - zz;
    out.M[1][2] = 2.0 * ( q.y * q.z - q.w * q.x );
    out.M[1][3] = 0.0;

    out.M[2][0] = 2.0 * ( q.x * q.z - q.w * q.y );
    out.M[2][1] = 2.0 * ( q.y * q.z + q.w * q.x );
    out.M[2][2] = ww - xx - yy + zz;
    out.M[2][3] = 0.0;

    out.M[3][0] = 0.0;
    out.M[3][1] = 0.0;
    out.M[3][2] = 0.0;
    out.M[3][3] = 1.0;

    out
}

// Convert a standard projection matrix into a TexCoordsFromTanAngles matrix for
// the primary time warp surface.
pub fn ovrMatrix4f_TanAngleMatrixFromProjection(projection: &ovrMatrix4f) -> ovrMatrix4f {
    /*
        A projection matrix goes from a view point to NDC, or -1 to 1 space.
        Scale and bias to convert that to a 0 to 1 space.

        const ovrMatrix3f m =
        { {
            { projection.M[0][0],                0.0, projection.M[0][2] },
            {                0.0, projection.M[1][1], projection.M[1][2] },
            {                0.0,                0.0,               -1.0 }
        } };
        // Note that there is no Y-flip because eye buffers have 0,0 = left-bottom.
        const ovrMatrix3f s = ovrMatrix3f_CreateScaling( 0.5, 0.5 );
        const ovrMatrix3f t = ovrMatrix3f_CreateTranslation( 0.5, 0.5 );
        const ovrMatrix3f r0 = ovrMatrix3f_Multiply( &s, &m );
        const ovrMatrix3f r1 = ovrMatrix3f_Multiply( &t, &r0 );
        return r1;

        clipZ = ( z * projection[2][2] + projection[2][3] ) / ( projection[3][2] * z )
        z = projection[2][3] / ( clipZ * projection[3][2] - projection[2][2] )
        z = ( projection[2][3] / projection[3][2] ) / ( clipZ - projection[2][2] / projection[3][2] )
    */

    let tanAngleMatrix = ovrMatrix4f {
        M: [[ 0.5 * projection.M[0][0], 0.0, 0.5 * projection.M[0][2] - 0.5, 0.0 ],
            [ 0.0, 0.5 * projection.M[1][1], 0.5 * projection.M[1][2] - 0.5, 0.0 ],
            [ 0.0, 0.0, -1.0, 0.0 ],
            // Store the values to convert a clip-Z to a linear depth in the unused matrix elements.
            [ projection.M[2][2], projection.M[2][3], projection.M[3][2], 1.0 ]]
    };

    tanAngleMatrix
}

// If a simple quad defined as a -1 to 1 XY unit square is transformed to
// the camera view with the given modelView matrix, it can alternately be
// drawn as a time warp overlay image to take advantage of the full window
// resolution, which is usually higher than the eye buffer textures, and
// avoids resampling both into the eye buffer, and again to the screen.
// This is used for high quality movie screens and user interface planes.
//
// Note that this is NOT an MVP matrix -- the "projection" is handled
// by the distortion process.
//
// This utility functions converts a model-view matrix that would normally
// draw a -1 to 1 unit square to the view into a TexCoordsFromTanAngles matrix 
// for an overlay surface.
//
// The resulting z value should be straight ahead distance to the plane.
// The x and y values will be pre-multiplied by z for projective texturing.
pub fn ovrMatrix4f_TanAngleMatrixFromUnitSquare(modelView: &ovrMatrix4f) -> ovrMatrix4f{
    /*
        // Take the inverse of the view matrix because the view matrix transforms the unit square
        // from world space into view space, while the matrix needed here is the one that transforms
        // the unit square from view space to world space.
        const ovrMatrix4f inv = ovrMatrix4f_Inverse( modelView );
        // This matrix calculates the projection onto the (-1, 1) X and Y axes of the unit square,
        // of the intersection of the vector (tanX, tanY, -1) with the plane described by the matrix
        // that transforms the unit square into world space.
        const ovrMatrix3f m =
        { {
            {    inv.M[0][0] * inv.M[2][3] - inv.M[0][3] * inv.M[2][0],
                inv.M[0][1] * inv.M[2][3] - inv.M[0][3] * inv.M[2][1],
                inv.M[0][2] * inv.M[2][3] - inv.M[0][3] * inv.M[2][2] },
            {    inv.M[1][0] * inv.M[2][3] - inv.M[1][3] * inv.M[2][0],
                inv.M[1][1] * inv.M[2][3] - inv.M[1][3] * inv.M[2][1],
                inv.M[1][2] * inv.M[2][3] - inv.M[1][3] * inv.M[2][2] },
            {    - inv.M[2][0],
                - inv.M[2][1],
                - inv.M[2][2] }
        } };
        // Flip the Y because textures have 0,0 = left-top as opposed to left-bottom.
        const ovrMatrix3f f = ovrMatrix3f_CreateScaling( 1.0, -1.0 );
        const ovrMatrix3f s = ovrMatrix3f_CreateScaling( 0.5, 0.5 );
        const ovrMatrix3f t = ovrMatrix3f_CreateTranslation( 0.5, 0.5 );
        const ovrMatrix3f r0 = ovrMatrix3f_Multiply( &f, &m );
        const ovrMatrix3f r1 = ovrMatrix3f_Multiply( &s, &r0 );
        const ovrMatrix3f r2 = ovrMatrix3f_Multiply( &t, &r1 );
        return r2;
    */

    let inv = ovrMatrix4f_Inverse( modelView );
    let coef = if inv.M[2][3] > 0.0 {
        1.0
    } else { 
        -1.0
    };

    let mut m: ovrMatrix4f = unsafe { mem::uninitialized() };
    m.M[0][0] = ( 0.5 * ( inv.M[0][0] * inv.M[2][3] - inv.M[0][3] * inv.M[2][0] ) - 0.5 * inv.M[2][0] ) * coef;
    m.M[0][1] = ( 0.5 * ( inv.M[0][1] * inv.M[2][3] - inv.M[0][3] * inv.M[2][1] ) - 0.5 * inv.M[2][1] ) * coef;
    m.M[0][2] = ( 0.5 * ( inv.M[0][2] * inv.M[2][3] - inv.M[0][3] * inv.M[2][2] ) - 0.5 * inv.M[2][2] ) * coef;
    m.M[0][3] = 0.0;

    m.M[1][0] = ( -0.5 * ( inv.M[1][0] * inv.M[2][3] - inv.M[1][3] * inv.M[2][0] ) - 0.5 * inv.M[2][0] ) * coef;
    m.M[1][1] = ( -0.5 * ( inv.M[1][1] * inv.M[2][3] - inv.M[1][3] * inv.M[2][1] ) - 0.5 * inv.M[2][1] ) * coef;
    m.M[1][2] = ( -0.5 * ( inv.M[1][2] * inv.M[2][3] - inv.M[1][3] * inv.M[2][2] ) - 0.5 * inv.M[2][2] ) * coef;
    m.M[1][3] = 0.0;

    m.M[2][0] = ( -inv.M[2][0] ) * coef;
    m.M[2][1] = ( -inv.M[2][1] ) * coef;
    m.M[2][2] = ( -inv.M[2][2] ) * coef;
    m.M[2][3] = 0.0;

    m.M[3][0] = 0.0;
    m.M[3][1] = 0.0;
    m.M[3][2] = 0.0;
    m.M[3][3] = 1.0;

    m
}

// Convert a standard view matrix into a TexCoordsFromTanAngles matrix for
// the looking into a cube map.
pub fn ovrMatrix4f_TanAngleMatrixForCubeMap(viewMatrix: &ovrMatrix4f) -> ovrMatrix4f {
    let mut m = *viewMatrix;
    // clear translation
    for i in 0..3 {
        m.M[ i ][ 3 ] = 0.0;
    }

    ovrMatrix4f_Inverse(&m)
}

// Utility function to calculate external velocity for smooth stick yaw turning.
// To reduce judder in FPS style experiences when the application framerate is
// lower than the vsync rate, the rotation from a joypad can be applied to the
// view space distorted eye vectors before applying the time warp.
pub fn ovrMatrix4f_CalculateExternalVelocity(viewMatrix: &ovrMatrix4f, yawRadiansPerSecond: f32) -> ovrMatrix4f {
    let angle = yawRadiansPerSecond * ( -1.0 / 60.0 );
    let sinHalfAngle = ( angle * 0.5 ).sin();
    let cosHalfAngle = ( angle * 0.5 ).cos();

    // Yaw is always going to be around the world Y axis
    let mut quat: ovrQuatf = unsafe { mem::uninitialized() };
    quat.x = viewMatrix.M[0][1] * sinHalfAngle;
    quat.y = viewMatrix.M[1][1] * sinHalfAngle;
    quat.z = viewMatrix.M[2][1] * sinHalfAngle;
    quat.w = cosHalfAngle;

    ovrMatrix4f_CreateFromQuaternion( &quat )
}

//-----------------------------------------------------------------
// Default initialization helper functions.
//-----------------------------------------------------------------


// Utility function to default initialize the ovrInitParms.
pub fn vrapi_DefaultInitParms(java: *const ovrJava) -> ovrInitParms {
    let mut parms: ovrInitParms = unsafe { mem::zeroed() };

    parms.Type = VRAPI_STRUCTURE_TYPE_INIT_PARMS;
    parms.ProductVersion = VRAPI_PRODUCT_VERSION as i32;
    parms.MajorVersion = VRAPI_MAJOR_VERSION as i32;
    parms.MinorVersion = VRAPI_MINOR_VERSION as i32;
    parms.PatchVersion = VRAPI_PATCH_VERSION as i32;
    parms.GraphicsAPI = VRAPI_GRAPHICS_API_OPENGL_ES_2;
    parms.Java = unsafe { *java };

    return parms;
}

// Utility function to default initialize the ovrModeParms.
pub fn vrapi_DefaultModeParms(java: *const ovrJava) -> ovrModeParms {
    let mut parms: ovrModeParms = unsafe { mem::zeroed() };

    parms.Type = VRAPI_STRUCTURE_TYPE_MODE_PARMS;
    parms.Flags |= ovrModeFlags::VRAPI_MODE_FLAG_ALLOW_POWER_SAVE as u32;
    parms.Flags |= ovrModeFlags::VRAPI_MODE_FLAG_RESET_WINDOW_FULLSCREEN as u32;
    parms.Java = unsafe { *java };

    parms
}

// Utility function to default initialize the ovrPerformanceParms.
pub fn vrapi_DefaultPerformanceParms() -> ovrPerformanceParms {
    let mut parms: ovrPerformanceParms = unsafe { mem::zeroed() };
    parms.CpuLevel = 2;
    parms.GpuLevel = 2;
    parms.MainThreadTid = 0;
    parms.RenderThreadTid = 0;

    parms
}

// Utility function to default initialize the ovrHeadModelParms.
pub fn vrapi_DefaultHeadModelParms() -> ovrHeadModelParms {
    let mut parms: ovrHeadModelParms = unsafe { mem::zeroed() };

    parms.InterpupillaryDistance    = 0.0640;    // average interpupillary distance
    parms.EyeHeight                 = 1.6750;    // average eye height above the ground when standing
    parms.HeadModelDepth            = 0.0805;
    parms.HeadModelHeight           = 0.0750;

    return parms;
}

// Utility function to default initialize the ovrFrameParms.
pub fn vrapi_DefaultFrameParms(java: *const ovrJava,
                               init: ovrFrameInit,
                               currentTime: f64,
                               textureSwapChain: *mut ovrTextureSwapChain)
                               -> ovrFrameParms
{
    let projectionMatrix = ovrMatrix4f_CreateProjectionFov( 90.0, 90.0, 0.0, 0.0, 0.1, 0.0 );
    let texCoordsFromTanAngles = ovrMatrix4f_TanAngleMatrixFromProjection( &projectionMatrix );

    let mut parms: ovrFrameParms = unsafe { mem::zeroed() };

    parms.Type = VRAPI_STRUCTURE_TYPE_FRAME_PARMS;
    for layer in 0..ovrFrameLayerType::VRAPI_FRAME_LAYER_TYPE_MAX as usize{
        parms.Layers[layer].ColorScale = 1.0;
        for eye in 0..ovrFrameLayerEye::VRAPI_FRAME_LAYER_EYE_MAX as usize {
            parms.Layers[layer].Textures[eye].TexCoordsFromTanAngles = texCoordsFromTanAngles;
            parms.Layers[layer].Textures[eye].TextureRect.width = 1.0;
            parms.Layers[layer].Textures[eye].TextureRect.height = 1.0;
            parms.Layers[layer].Textures[eye].HeadPose.Pose.Orientation.w = 1.0;
            parms.Layers[layer].Textures[eye].HeadPose.TimeInSeconds = currentTime;
        }
    }
    parms.LayerCount = 1;
    parms.MinimumVsyncs = 1;
    parms.ExtraLatencyMode = ovrExtraLatencyMode::VRAPI_EXTRA_LATENCY_MODE_OFF;
    parms.ExternalVelocity.M[0][0] = 1.0;
    parms.ExternalVelocity.M[1][1] = 1.0;
    parms.ExternalVelocity.M[2][2] = 1.0;
    parms.ExternalVelocity.M[3][3] = 1.0;
    parms.PerformanceParms = vrapi_DefaultPerformanceParms();
    parms.Java = unsafe { *java };

    parms.Layers[0].SrcBlend = ovrFrameLayerBlend::VRAPI_FRAME_LAYER_BLEND_ONE;
    parms.Layers[0].DstBlend = ovrFrameLayerBlend::VRAPI_FRAME_LAYER_BLEND_ZERO;
    parms.Layers[0].Flags = 0;

    parms.Layers[1].SrcBlend = ovrFrameLayerBlend::VRAPI_FRAME_LAYER_BLEND_SRC_ALPHA;
    parms.Layers[1].DstBlend = ovrFrameLayerBlend::VRAPI_FRAME_LAYER_BLEND_ONE_MINUS_SRC_ALPHA;
    parms.Layers[1].Flags = 0;

    match init
    {
        VRAPI_FRAME_INIT_DEFAULT =>
        {
            //break;
        },
        VRAPI_FRAME_INIT_BLACK |
        VRAPI_FRAME_INIT_BLACK_FLUSH |
        VRAPI_FRAME_INIT_BLACK_FINAL =>
        {
            parms.Flags = ovrFrameFlags::VRAPI_FRAME_FLAG_INHIBIT_SRGB_FRAMEBUFFER as i32;
            for eye in 0..ovrFrameLayerEye::VRAPI_FRAME_LAYER_EYE_MAX as usize
            {
                parms.Layers[0].Textures[eye].ColorTextureSwapChain = unsafe { mem::transmute(ovrDefaultTextureSwapChain::VRAPI_DEFAULT_TEXTURE_SWAPCHAIN_BLACK as usize) };
            }
        },
        VRAPI_FRAME_INIT_LOADING_ICON |
        VRAPI_FRAME_INIT_LOADING_ICON_FLUSH =>
        {
            parms.LayerCount = 2;
            parms.Flags = ovrFrameFlags::VRAPI_FRAME_FLAG_INHIBIT_SRGB_FRAMEBUFFER as i32;
            parms.Layers[1].Flags = ovrFrameLayerFlags::VRAPI_FRAME_LAYER_FLAG_SPIN as i32;
            parms.Layers[1].SpinSpeed = 1.0;        // rotation in radians per second
            parms.Layers[1].SpinScale = 16.0;        // icon size factor smaller than fullscreen
            for eye in 0..ovrFrameLayerEye::VRAPI_FRAME_LAYER_EYE_MAX as usize
            {
                parms.Layers[0].Textures[eye].ColorTextureSwapChain = unsafe { mem::transmute(ovrDefaultTextureSwapChain::VRAPI_DEFAULT_TEXTURE_SWAPCHAIN_BLACK as usize) };
				parms.Layers[1].Textures[eye].ColorTextureSwapChain = if !textureSwapChain.is_null() {
					textureSwapChain
				} else {
					unsafe { mem::transmute(ovrDefaultTextureSwapChain::VRAPI_DEFAULT_TEXTURE_SWAPCHAIN_LOADING_ICON as usize) }
				};
            }
        },
        VRAPI_FRAME_INIT_MESSAGE |
        VRAPI_FRAME_INIT_MESSAGE_FLUSH =>
        {
            parms.LayerCount = 2;
            parms.Flags = ovrFrameFlags::VRAPI_FRAME_FLAG_INHIBIT_SRGB_FRAMEBUFFER as i32;
            parms.Layers[1].SpinSpeed = 0.0;        // rotation in radians per second
            parms.Layers[1].SpinScale = 2.0;        // message size factor smaller than fullscreen
            for eye in 0..ovrFrameLayerEye::VRAPI_FRAME_LAYER_EYE_MAX as usize
            {
                parms.Layers[0].Textures[eye].ColorTextureSwapChain = unsafe { mem::transmute(ovrDefaultTextureSwapChain::VRAPI_DEFAULT_TEXTURE_SWAPCHAIN_BLACK as usize) };
				parms.Layers[1].Textures[eye].ColorTextureSwapChain = if !textureSwapChain.is_null() {
					textureSwapChain
				} else {
					unsafe { mem::transmute(ovrDefaultTextureSwapChain::VRAPI_DEFAULT_TEXTURE_SWAPCHAIN_LOADING_ICON as usize) }
				};
            }
        }
    }

    if init == VRAPI_FRAME_INIT_BLACK_FLUSH || init == VRAPI_FRAME_INIT_LOADING_ICON_FLUSH || init == VRAPI_FRAME_INIT_MESSAGE_FLUSH
    {
        parms.Flags |= ovrFrameFlags::VRAPI_FRAME_FLAG_FLUSH as i32;
    }
    if init == VRAPI_FRAME_INIT_BLACK_FINAL
    {
        parms.Flags |= ovrFrameFlags::VRAPI_FRAME_FLAG_FLUSH as i32 | ovrFrameFlags::VRAPI_FRAME_FLAG_FINAL as i32;
    }

    return parms;
}

//-----------------------------------------------------------------
// Eye view matrix helper functions.
//-----------------------------------------------------------------

// Apply the head-on-a-stick model if head tracking is not available.
pub fn vrapi_ApplyHeadModel(head_model_params: &ovrHeadModelParms, tracking: &ovrTracking) -> ovrTracking {
    if (tracking.Status & ovrTrackingStatus::VRAPI_TRACKING_STATUS_POSITION_TRACKED as u32) == 0 {
        // Calculate the head position based on the head orientation using a head-on-a-stick model.
        let p = head_model_params;
        let m = ovrMatrix4f_CreateFromQuaternion( &tracking.HeadPose.Pose.Orientation );
        let mut new_tracking  = *tracking;
        new_tracking.HeadPose.Pose.Position.x = m.M[0][1] * p.HeadModelHeight - m.M[0][2] * p.HeadModelDepth;
        new_tracking.HeadPose.Pose.Position.y = m.M[1][1] * p.HeadModelHeight - m.M[1][2] * p.HeadModelDepth - p.HeadModelHeight;
        new_tracking.HeadPose.Pose.Position.z = m.M[2][1] * p.HeadModelHeight - m.M[2][2] * p.HeadModelDepth;
        return new_tracking;
    }

    *tracking
}

// Utility function to get the center eye transform.
pub fn vrapi_GetCenterEyeTransform(_head_model_params: &ovrHeadModelParms,
                                   tracking: &ovrTracking,
                                   input: Option<&ovrMatrix4f>)
                                   -> ovrMatrix4f
{
    //VRAPI_UNUSED( headModelParms );

    // Controller input is expected to be applied relative to the head in neutral position, which means
    // ovrTracking::HeadPose.Pose.Translation should be relative to the center of the head in neutral position.
    let center_eye_rotation = ovrMatrix4f_CreateFromQuaternion(&tracking.HeadPose.Pose.Orientation);
    let center_eye_offset = tracking.HeadPose.Pose.Position;
    let center_eye_translation = ovrMatrix4f_CreateTranslation(center_eye_offset.x, center_eye_offset.y, center_eye_offset.z);
    let center_eye_transform = ovrMatrix4f_Multiply(&center_eye_translation, &center_eye_rotation);

    match input {
        Some(input) => ovrMatrix4f_Multiply(input, &center_eye_transform),
        None => center_eye_transform
    }
}

// Utility function to get the center eye view matrix.
// Pass in NULL for 'input' if there is no additional controller input.
pub fn vrapi_GetCenterEyeViewMatrix(head_model_params: &ovrHeadModelParms,
                                    tracking: &ovrTracking,
                                    input: Option<&ovrMatrix4f>)
                                    -> ovrMatrix4f
{
    let center_eye_transform = vrapi_GetCenterEyeTransform(head_model_params, tracking, input);
    ovrMatrix4f_Inverse(&center_eye_transform)
}

// Utility function to get the eye view matrix based on the center eye view matrix and the IPD.
pub fn vrapi_GetEyeViewMatrix(head_model_params: &ovrHeadModelParms,
                              center_eye_view_matrix: &ovrMatrix4f,
                              eye: i32)
                              -> ovrMatrix4f
{
    let eye_offset = ( if eye > 0 { -0.5 } else { 0.5 } ) * head_model_params.InterpupillaryDistance;
    let eye_offset_matrix = ovrMatrix4f_CreateTranslation(eye_offset, 0.0, 0.0);

    ovrMatrix4f_Multiply(&eye_offset_matrix, center_eye_view_matrix)
}
