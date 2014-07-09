// REDUCTION SUM KERNELS

// interleaved addressing
__kernel void sum0(__global float2 *rd_in,
							__local float2 *s_data,
							__global float2 *sum_partial)
{
	unsigned int lid = get_local_id(0);

	unsigned int gid = get_global_id(0);

	// store global value into local shared memory address
	s_data[lid] = rd_in[gid];

	barrier(CLK_LOCAL_MEM_FENCE); /* Wait for others in the work-group */

	// do reduction to find sum of work group
	for(unsigned int s = 1; s < get_local_size(0); s*=2)
	{
		if(lid % (2*s) == 0)
		{
			s_data[lid] += s_data[lid + s];
		}

		barrier(CLK_LOCAL_MEM_FENCE); /* Wait for others in the work-group */	
	}

	// write partial result for this work group
	if(lid == 0)
	{
		sum_partial[get_group_id(0)] = s_data[0];
	}
}

// sequential addressing
__kernel void sum1(__global float2 *rd_in,
							__local float2 *s_data,
							__global float2 *sum_partial)
{
	unsigned int lid = get_local_id(0);

	unsigned int gid = get_global_id(0);

	// store global value into local shared memory address
	s_data[lid] = rd_in[gid];

	barrier(CLK_LOCAL_MEM_FENCE); /* Wait for others in the work-group */

	// do reduction to find sum of work group
	for(unsigned int s = 1; s < get_local_size(0); s*=2)
	{
		int index = 2 * s * lid;

		if(index < get_local_size(0))
		{
			s_data[index] += s_data[index + s];
		}

		barrier(CLK_LOCAL_MEM_FENCE); /* Wait for others in the work-group */	
	}

	// write partial result for this work group
	if(lid == 0)
	{
		sum_partial[get_group_id(0)] = s_data[0];
	}
}

// half gws
// first add during load
__kernel void sum2(__global float2 *rd_in,
							__local float2 *s_data,
							__global float2 *sum_partial)
{
	unsigned int lid = get_local_id(0);

	unsigned int gid = get_group_id(0)*(get_local_size(0)*2) + lid;

	// store global value into local shared memory address
	s_data[lid] = rd_in[gid] + rd_in[gid + get_local_size(0)];

	barrier(CLK_LOCAL_MEM_FENCE); /* Wait for others in the work-group */

	// do reduction to find sum of work group
	for(unsigned int s = get_local_size(0)/2; s > 0; s >>= 1)
	{
		if(lid < s)
		{
			s_data[lid] += s_data[lid + s];
		}

		barrier(CLK_LOCAL_MEM_FENCE); /* Wait for others in the work-group */	
	}

	// write partial result for this work group
	if(lid == 0)
	{
		sum_partial[get_group_id(0)] = s_data[0];
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
__kernel void sum3(__global float2 *rd_in,
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