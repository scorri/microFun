
//--------------------------------------------------------------------------------------
// getBufferIndex
// given a thread DTid and group Gid, calculate the corresponding cell index in the buffer
//--------------------------------------------------------------------------------------
int getBufferIndexUseDispatch( int2 dispatch_id, int2 group_id, int2 win_dim )
{
    int2 wrapped = ( dispatch_id - 2 * group_id - 1 + win_dim ) % win_dim;
    return wrapped.x + wrapped.y * win_dim.x;
}

//--------------------------------------------------------------------------------------
// threadInWindow
// returns true if the cell corresponding to thread DTid in group Gid is in the application window
//--------------------------------------------------------------------------------------
bool threadInWindow( int2 dispatch_id, int2 group_id, int2 win_dim )
{
    bool inWindowX = dispatch_id.x - 2 * group_id.x - 1 < win_dim.x && dispatch_id.x > 0;
    bool inWindowY = dispatch_id.y - 2 * group_id.y - 1 < win_dim.y && dispatch_id.y > 0;
    return inWindowX && inWindowY;
}

//--------------------------------------------------------------------------------------
// threadOnGroupEdge
// returns true if the group thread GTid is part of the overlapping edge region
//--------------------------------------------------------------------------------------
bool threadOnGroupEdge( int2 local_id, int2 local_size )
{
    return local_id.x == 0 || 
           local_id.x == local_size.x - 1 ||
           local_id.y == 0 ||
           local_id.y == local_size.y - 1;
}

int getOffset( int GI, int xOffset, int yOffset, int local_size_x )
{
    return GI + xOffset + yOffset * local_size_x ;
}


bool isNear( int2 a, int2 b )
{
    return  (
                a.x == b.x     ||
                a.x == b.x + 1 ||
                a.x == b.x - 1
            ) && (
                a.y == b.y     ||
                a.y == b.y + 1 ||
                a.y == b.y - 1
            );
}

//--------------------------------------------------------------------------------------
// threadIDToCell
// given a thread DTid and group Gid, calculate the corresponding cell position
//--------------------------------------------------------------------------------------
int2 threadIDToCell( int2 dispatch_id, int2 group_id )
{
    return dispatch_id - 2 * group_id - 1;
}

__kernel void conway_kernel(__global int *conway_in,
							int width,
							int height,
							__local int *s_data,
							__global int *conway_out)
{
	int2 l_size = (int2)(get_local_size(0), get_local_size(1));
	int2 g_id = (int2)(get_group_id(0), get_group_id(1));
	int2 l_id = (int2)(get_local_id(0), get_local_id(1));
	int2 w_dim = (int2)(width, height); 
	int2 dt_id = (int2)(l_size*g_id + l_id); //Dispatch Thread ID (DC)

	int idx = getBufferIndexUseDispatch(dt_id, g_id, w_dim);

	// store the state of the current cell in shared memory for others in work group
    uint local_id = l_id.y*l_size.x + l_id.x;
	s_data[local_id] = conway_in[idx];

	barrier(CLK_LOCAL_MEM_FENCE); /* Wait for others in the work-group */

	bool inWindow = threadInWindow(dt_id, g_id, w_dim);
	bool onEdge = threadOnGroupEdge( l_id, l_size);
	if ( inWindow  && !onEdge )
    {
		int sum = s_data[ getOffset(local_id, -1, -1, l_size.x) ];
			sum += s_data[ getOffset(local_id, 0, -1, l_size.x) ];
			sum += s_data[ getOffset(local_id, 1, -1, l_size.x) ];
			
			sum += s_data[ getOffset(local_id, -1, 0, l_size.x) ];
			sum += s_data[ getOffset(local_id, +1, 0, l_size.x) ];
			
			sum += s_data[ getOffset(local_id, -1, 1, l_size.x) ];
			sum += s_data[ getOffset(local_id,  0, 1, l_size.x) ];
			sum += s_data[ getOffset(local_id,  1, 1, l_size.x) ];
			
		int nextState;
		if(sum == 3) nextState = 1;
		else if(sum<2 || sum >3) nextState = 0;
		else nextState = s_data[local_id];

		conway_out[idx] = nextState;
	}
}

__kernel void copy_kernel(__global int *conway_in,
							__global int *conway_out)
{
	int idx = get_global_id(0);

	conway_out[idx] = conway_in[idx];
}
