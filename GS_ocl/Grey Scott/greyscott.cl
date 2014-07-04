#define timek 1
//Gray Scott conrolling parameters
#define diffU 0.16
#define diffV 0.08
#define F 0.035
#define k 0.06


//--------------------------------------------------------------------------------------
// getBufferIndex
// given a thread DTid and group Gid, calculate the corresponding cell index in the buffer
//--------------------------------------------------------------------------------------
uint getBufferIndexUseDispatch( uint2 dispatch_id, uint2 group_id, uint2 win_dim )
{
    uint2 wrapped = ( dispatch_id - 2 * group_id - 1 + win_dim ) % win_dim;
    return wrapped.x + wrapped.y * win_dim.x;
}

//--------------------------------------------------------------------------------------
// threadInWindow
// returns true if the cell corresponding to thread DTid in group Gid is in the application window
//--------------------------------------------------------------------------------------
bool threadInWindow( uint2 dispatch_id, uint2 group_id, uint2 win_dim )
{
    bool inWindowX = dispatch_id.x - 2 * group_id.x - 1 < win_dim.x && dispatch_id.x > 0;
    bool inWindowY = dispatch_id.y - 2 * group_id.y - 1 < win_dim.y && dispatch_id.y > 0;
    return inWindowX && inWindowY;
}

//--------------------------------------------------------------------------------------
// threadOnGroupEdge
// returns true if the group thread GTid is part of the overlapping edge region
//--------------------------------------------------------------------------------------
bool threadOnGroupEdge( uint2 local_id, uint2 local_size )
{
    return local_id.x == 0 || 
           local_id.x == local_size.x - 1 ||
           local_id.y == 0 ||
           local_id.y == local_size.y - 1;
}

uint getOffset( uint GI, int xOffset, int yOffset, uint local_size_x )
{
    return GI + xOffset + yOffset * local_size_x ;
}


//--------------------------------------------------------------------------------------
// threadIDToCell
// given a thread DTid and group Gid, calculate the corresponding cell position
//--------------------------------------------------------------------------------------
uint2 threadIDToCell( uint2 dispatch_id, uint2 group_id )
{
    return dispatch_id - 2 * group_id - 1;
}

__kernel void gs_rd_kernel(__global float2 *rd_in,
							int width,
							int height,
							__local float2 *s_data,
							__global float2 *rd_out)
{
	uint2 l_size = (uint2)(get_local_size(0), get_local_size(1));
	uint2 g_id = (uint2)(get_group_id(0), get_group_id(1));
	uint2 l_id = (uint2)(get_local_id(0), get_local_id(1));
	uint2 w_dim = (uint2)(width, height);
	uint2 dt_id = (uint2)(l_size*g_id + l_id); //Dispatch Thread ID (DC)

	uint idx = getBufferIndexUseDispatch(dt_id, g_id, w_dim);

	// store the state of the current cell in shared memory for others in work group
    uint local_id = l_id.y*l_size.x + l_id.x;
	s_data[local_id] = rd_in[idx];

	barrier(CLK_LOCAL_MEM_FENCE); /* Wait for others in the work-group */

	if ( threadInWindow(dt_id, g_id, w_dim)  && !threadOnGroupEdge(l_id, l_size ) )
    {
		float2 sum = s_data[ getOffset(local_id, 0, -1, l_size.x) ];
			  sum += s_data[ getOffset(local_id, -1, 0, l_size.x)];
			  sum += s_data[ getOffset(local_id, 1, 0, l_size.x) ];
			  sum += s_data[ getOffset(local_id, 0, 1, l_size.x) ];
			
		float2 lap = sum - 4* s_data[local_id];
		
		float uvv = s_data[local_id].x * s_data[local_id].y * s_data[local_id].y;
		
		float delta_U = diffU*lap.x - uvv + F*(1 - s_data[local_id].x);
		float delta_V = diffV*lap.y + uvv - (k+F)*s_data[local_id].y;

		rd_out[idx] = (float2)(s_data[local_id].x + delta_U*timek , s_data[local_id].y + delta_V*timek );
	}
}

__kernel void output_kernel(__write_only image2d_t im, 
							__global float2 *rd_buff)
{
	int x = (int)get_global_id(0);
	int y = (int)get_global_id(1);

	int idx = x + y*get_image_width(im);

	// Write the output value to image
	write_imagef(im, (int2)(x, y), (float4)(0.1f,(float)(1 - rd_buff[idx].x),(float)(1 - rd_buff[idx].y), 1.0f)); //2D
	//write_imagef(im, (int2)(x, y), (float4)(0.1f,(float)(1 - rd_buff[idx].x),0.0f, 1.0f)); // 1D

}

// unroll
void warpReduce(__local volatile float2* s_data, int tid)
{
	s_data[tid] += s_data[tid+32];
	s_data[tid] += s_data[tid+16];
	s_data[tid] += s_data[tid+ 8];
	s_data[tid] += s_data[tid+ 4];
	s_data[tid] += s_data[tid+ 2];
	s_data[tid] += s_data[tid+ 1];	
}

// half gws
// unroll last warp
__kernel void sum_kernel(__global float2 *rd_in,
							__local float2 *s_data,
							__global float2 *sum_partial)
{
	unsigned int lid = get_local_id(0);

	unsigned int gid = get_group_id(0)*(get_local_size(0)*2) + lid;

	// store global value into local shared memory address
	s_data[lid] = rd_in[gid] + rd_in[gid + get_local_size(0)];

	barrier(CLK_LOCAL_MEM_FENCE); /* Wait for others in the work-group */

	// do reduction to find sum of work group
	for(unsigned int s = get_local_size(0)/2; s > 32; s >>= 1)
	{
		if(lid < s)
		{
			s_data[lid] += s_data[lid + s];
		}

		barrier(CLK_LOCAL_MEM_FENCE); /* Wait for others in the work-group */	
	}

	// unroll last warp
	if(lid < 32) warpReduce(s_data, lid);

	// write partial result for this work group
	if(lid == 0)
	{
		sum_partial[get_group_id(0)] = s_data[0];
	}
}