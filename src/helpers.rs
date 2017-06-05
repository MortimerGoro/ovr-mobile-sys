use std::mem;
use super::*;
use super::ovrGraphicsAPI::*;
use super::ovrStructureType::*;

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

// Returns a projection matrix based on the given FOV.
pub fn ovrMatrix4f_CreateProjectionFov(fov_degrees_x: f32,
									   fov_degrees_y: f32,
									   offset_x: f32,
									   offset_y: f32,
									   near_z: f32,
									   far_z: f32)
									   -> ovrMatrix4f
{
	let half_width = near_z * (fov_degrees_x * (VRAPI_PI as f32 / 180.0f32 * 0.5f32)).tan();
	let half_height = near_z * (fov_degrees_y * (VRAPI_PI as f32 / 180.0f32 * 0.5f32)).tan();

	let min_x = offset_x - half_width;
	let max_x = offset_x + half_width;

	let min_y = offset_y - half_height;
	let max_y = offset_y + half_height;

	ovrMatrix4f_CreateProjection( min_x, max_x, min_y, max_y, near_z, far_z )
}

// Returns a projection matrix based on the specified dimensions.
// The projection matrix transforms -Z=forward, +Y=up, +X=right to the appropriate clip space for the graphics API.
// The far plane is placed at infinity if far_z <= near_z.
// An infinite projection matrix is preferred for rasterization because, except for
// things *right* up against the near plane, it always provides better precision:
//		"Tightening the Precision of Perspective Rendering"
//		Paul Upchurch, Mathieu Desbrun
//		Journal of Graphics Tools, Volume 16, Issue 1, 2012
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
	let offsetZ = near_z;	// set to zero for a [0,1] clip space

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