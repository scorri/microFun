
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
int getBufferIndexUseDispatch( int3 dispatch_id, int3 group_id, int3 win_dim )
{
    int3 wrapped = ( dispatch_id - 2 * group_id - 1 + win_dim ) % win_dim;
    return wrapped.x + wrapped.y * win_dim.x + wrapped.z * win_dim.x * win_dim.y;
}

//--------------------------------------------------------------------------------------
// threadInWindow
// returns true if the cell corresponding to thread DTid in group Gid is in the application window
//--------------------------------------------------------------------------------------
bool threadInWindow( int3 dispatch_id, int3 group_id, int3 win_dim )
{
    bool inWindowX = dispatch_id.x - 2 * group_id.x - 1 < win_dim.x && dispatch_id.x > 0;
    bool inWindowY = dispatch_id.y - 2 * group_id.y - 1 < win_dim.y && dispatch_id.y > 0;
    bool inWindowZ = dispatch_id.z - 2 * group_id.z - 1 < win_dim.z && dispatch_id.z > 0;
    return inWindowX && inWindowY && inWindowZ;
}

//--------------------------------------------------------------------------------------
// threadOnGroupEdge
// returns true if the group thread GTid is part of the overlapping edge region
//--------------------------------------------------------------------------------------
bool threadOnGroupEdge( int3 local_id, int3 local_size )
{
    return local_id.x == 0 || 
           local_id.x == local_size.x - 1 ||
           local_id.y == 0 ||
           local_id.y == local_size.y - 1 ||
		   local_id.z == 0 ||
           local_id.z == local_size.z - 1;
}

int getOffset( int GI, int xOffset, int yOffset, int local_size_x, int zOffset, int local_size_y )
{
    return GI + xOffset + yOffset * local_size_x  + zOffset * local_size_x * local_size_y;
}


bool isNear( int3 a, int3 b )
{
    return  (
                a.x == b.x     ||
                a.x == b.x + 1 ||
                a.x == b.x - 1
            ) && (
                a.y == b.y     ||
                a.y == b.y + 1 ||
                a.y == b.y - 1
            ) && (
                a.z == b.z     ||
                a.z == b.z + 1 ||
                a.z == b.z - 1
            );
}

//--------------------------------------------------------------------------------------
// threadIDToCell
// given a thread DTid and group Gid, calculate the corresponding cell position
//--------------------------------------------------------------------------------------
int3 threadIDToCell( int3 dispatch_id, int3 group_id )
{
    return dispatch_id - 2 * group_id - 1;
}

__kernel void output_kernel(__global float4 *colour, 
							__global float2 *rd_buff)
{
	int idx = get_global_id(0) + get_global_size(0)*get_global_id(1) + get_global_size(0)*get_global_size(1)*get_global_id(2);

	float2 rd_colour = rd_buff[idx];

	float4 c = (float4)(0,0,0,0.0f);
	if(rd_colour.x > 0 && rd_colour.y > 0)
	{
		if(rd_colour.x > 0.1f && rd_colour.y > 0.1f)
		{
			c.w = 0.5f;
			if(rd_colour.x > rd_colour.y)	// U > V output is red
				c.x = 1.0f - rd_colour.x;
			else if(rd_colour.y > rd_colour.x) // V > U output is green	
				c.y = 1.0f - rd_colour.y;
			else
				c.z = 1.0f - rd_colour.x; // U == V output is blue 
		}
		//c = (float4)(0.1f, rd_colour.x, rd_colour.y, 0.5f);
	}

	colour[idx] = c;
	//colour[idx] = (float4)(1.0f, rd_colour.x, rd_colour.y, 0.5f);
}

__kernel 
void rd_kernel(__global float2 *rd_in,
			int4 params,
		__local float2 *s_data,
		__global float2 *rd_out)
{
	int3 l_size = (int3)(get_local_size(0), get_local_size(1), get_local_size(2));
	int3 g_id = (int3)(get_group_id(0), get_group_id(1), get_group_id(2));
	int3 l_id = (int3)(get_local_id(0), get_local_id(1), get_local_id(2));
	int3 w_dim = (int3)(params.xyz); 
	int3 dt_id = (int3)(l_size*g_id + l_id); //Dispatch Thread ID (DC)

	int idx = getBufferIndexUseDispatch(dt_id, g_id, w_dim);

	// store the state of the current cell in shared memory for others in work group
    uint local_id = l_id.x + l_id.y * l_size.x + l_id.z * l_size.x * l_size.y;
	s_data[local_id] = rd_in[idx];

	barrier(CLK_LOCAL_MEM_FENCE); /* Wait for others in the work-group */

	bool inWindow = threadInWindow(dt_id, g_id, w_dim);
	bool onEdge = threadOnGroupEdge( l_id, l_size);
	if ( inWindow  && !onEdge )
    {
		float2 sum = s_data[ getOffset(local_id, 0, -1, l_size.x, 0, l_size.y) ];
			sum += s_data[ getOffset(local_id, -1, 0, l_size.x, 0, l_size.y)];
			sum += s_data[ getOffset(local_id, 1, 0, l_size.x, 0, l_size.y) ];
			sum += s_data[ getOffset(local_id, 0, 1, l_size.x, 0, l_size.y) ];

			sum += s_data[ getOffset(local_id, 0, 0, l_size.x, -1, l_size.y) ];
			sum += s_data[ getOffset(local_id, 0, 0, l_size.x,  1, l_size.y) ];

		float2 lap = sum - 6 * s_data[local_id];

		float uvv = s_data[local_id].x * s_data[local_id].y * s_data[local_id].y;
		
		float delta_U = diffU*lap.x - uvv + F*(1 - s_data[local_id].x);
		float delta_V = diffV*lap.y + uvv - (k+F)*s_data[local_id].y;

		rd_out[idx] = (float2)(s_data[local_id].x + delta_U*timek , s_data[local_id].y + delta_V*timek );
	}
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
__kernel 
void sum_kernel(__global float2 *rd_in,
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