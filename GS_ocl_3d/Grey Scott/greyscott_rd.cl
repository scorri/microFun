
#define timek 1
//Gray Scott conrolling parameters
#define diffU 0.16
#define diffV 0.08
#define F 0.035
#define k 0.06

// microfun parameters
#define alphan 0.95
#define omega 100
#define betan 0.05
#define rho 1.2 //less than 1 for recycling
#define h 0.05

#define vDoc 5
#define insul 0 
#define kDoc 0.0000001
/*
//numerical parameter for diffusion and other processes
#define timek 0.0020833
//fungal genotype parameters
#define replenishment 0
#define uptake_trait 1.0 //check this
*/

// mouse constant buffer initial values
#define Db 0.18
#define Dv 0.13
#define theta 2.0

typedef float scalar;

struct SharedInfo
{
	scalar ni;
	scalar mb;
	int structure;
};

// definition of the Biomass cell buffer element
struct BiomassBufType
{
    scalar ni;
    // Non-insulated immobile biomass
	scalar mb;
    // Mobile biomass
	scalar ins;
    // Insulated biomass
	scalar ex_re;
    // Replenishment rate (?)
};


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


__kernel 
void microFun(__global struct BiomassBufType *rd_in,
			int4 params,
		__local struct SharedInfo *s_data,
		__global struct BiomassBufType *rd_out)
{
	float alphai=alphan/omega;
	float betai=betan*rho*alphai/alphan;


	int3 l_size = (int3)(get_local_size(0), get_local_size(1), get_local_size(2));
	int3 g_id = (int3)(get_group_id(0), get_group_id(1), get_group_id(2));
	int3 l_id = (int3)(get_local_id(0), get_local_id(1), get_local_id(2));
	int3 w_dim = (int3)(params.xyz); 
	int3 dt_id = (int3)(l_size*g_id + l_id); //Dispatch Thread ID (DC)

	int idx = getBufferIndexUseDispatch(dt_id, g_id, w_dim);



	int s= 0.5;//soil_struct.Load(sampleCoords);
	int temp_structure=0;
	if(s>0)
		temp_structure = 0;
	else
		temp_structure = 1;

	//check if thread is on sealing planes and if so set to solid
	if( dt_id.x - 2 * g_id.x - 1 == 0)
		temp_structure=0;
	if( dt_id.y - 2 * g_id.y - 1 == 0)
		temp_structure=0;
	if( dt_id.z - 2 * g_id.z - 1 == 0)
		temp_structure=0;

	if( dt_id.x - 2 * g_id.x - 1 == w_dim.x-1)
		temp_structure=0;
	if( dt_id.y - 2 * g_id.y - 1 == w_dim.y-1)
		temp_structure=0;
	if( dt_id.z - 2 * g_id.z - 1 == w_dim.z-1)
		temp_structure=0;	


    // Store the state of the current cell in shared memory for other threads in this group to access
	// NOTE: We store both ni and mb state in one pass
	uint local_id = l_id.x + l_id.y * l_size.x + l_id.z * l_size.x * l_size.y;
	s_data[local_id].ni = rd_in[idx].ni;
	s_data[local_id].mb = rd_in[idx].ins;
	s_data[local_id].structure = temp_structure;

	barrier(CLK_LOCAL_MEM_FENCE); /* Wait for others in the work-group */


	// updates ni and mb in this case x and y from rd_in to rd_out
	bool inWindow = threadInWindow(dt_id, g_id, w_dim);
	bool onEdge = threadOnGroupEdge( l_id, l_size);
	if ( inWindow  && !onEdge )
    {
		// calculate the sum of the adjacent cell's states
		struct SharedInfo sum;
		sum.ni = s_data[ getOffset(local_id, 0, -1, l_size.x, 0, l_size.y)].ni
				+ s_data[ getOffset(local_id, -1, 0, l_size.x, 0, l_size.y)].ni
				+ s_data[ getOffset(local_id, 1, 0, l_size.x, 0, l_size.y)].ni
				+ s_data[ getOffset(local_id, 0, 1, l_size.x, 0, l_size.y)].ni
				+ s_data[ getOffset(local_id, 0, 0, l_size.x, -1, l_size.y)].ni
				+ s_data[ getOffset(local_id, 0, 0, l_size.x, 1, l_size.y)].ni;

		sum.mb = s_data[ getOffset(local_id, 0, -1, l_size.x, 0, l_size.y)].mb
				+ s_data[ getOffset(local_id, -1, 0, l_size.x, 0, l_size.y)].mb
				+ s_data[ getOffset(local_id, 1, 0, l_size.x, 0, l_size.y)].mb
				+ s_data[ getOffset(local_id, 0, 1, l_size.x, 0, l_size.y)].mb
				+ s_data[ getOffset(local_id, 0, 0, l_size.x, -1, l_size.y)].mb
				+ s_data[ getOffset(local_id, 0, 0, l_size.x, 1, l_size.y)].mb;

		sum.structure = s_data[ getOffset(local_id, 0, -1, l_size.x, 0, l_size.y)].structure
				+ s_data[ getOffset(local_id, -1, 0, l_size.x, 0, l_size.y)].structure
				+ s_data[ getOffset(local_id, 1, 0, l_size.x, 0, l_size.y)].structure
				+ s_data[ getOffset(local_id, 0, 1, l_size.x, 0, l_size.y)].structure
				+ s_data[ getOffset(local_id, 0, 0, l_size.x, -1, l_size.y)].structure
				+ s_data[ getOffset(local_id, 0, 0, l_size.x, 1, l_size.y)].structure;

		// get current state of this threads cell
        //SharedInfo state = getCellByOffset(GI, 0, 0, 0);
		struct SharedInfo state = s_data[ getOffset(local_id, 0, 0, l_size.x, 0, l_size.y)];

        // STRUCTURE CONSTRAINED DIFFUSION _____  diffuse NI -  diffusion algorithm from GPU gems 2 chapter
		// diffuse ni biomass -  diffusion algorithm from GPU gems 2 chapter
		scalar coeff =( 1-(timek*sum.structure*Db)/(h*h) );
        float nextStateNI = state.ni * coeff + (timek/(h*h))*state.structure*Db*sum.ni;

        // finally, update the new biomass buffer resource with the next-state value determined
		rd_out[ idx ].ni = nextStateNI;

		//float coeff_out= 1- (timek*sum.structure*Dv)/h*h;
		// diffuse mobile biomass -  diffusion algorithm from GPU gems 2 chapter
		coeff = 1-((timek*sum.structure*Dv)/(float)(h*h));
        scalar nextStatemb = state.mb * coeff + (timek/(float)(h*h))*state.structure*Dv * sum.mb;
 
		//IS MOBILE BIOMASS DIFFUING
		// finally, update the new biomass buffer resource with the next-state value determined
		rd_out[ idx ].ins = nextStatemb*(state.ni >0);
	}

	barrier(CLK_LOCAL_MEM_FENCE); /* Wait for others in the work-group */

	if( inWindow && !onEdge)
	{
		// create a variables for storing the current cell's next-state value based on recycling
		scalar dtInsul = 0;
		scalar dtNetMobI = 0;
		scalar dtNetMobN = 0;
		scalar dtUptk = 0;
		scalar dtTotMob = 0;
		

		///////MICROBIAL REACTION////////////////////
		struct BiomassBufType copy_output = rd_out[idx];
		// ni == x
		// ins == y
		// ex_re == z
		// mb == w
		//if((BiomassBufferNew[ idx ].ni!=0)||(BiomassBufferNew[ idx ].ins!=0))
		if( (copy_output.ni != 0) || (copy_output.ins != 0) )
		{
			// UPTAKE/////////
			scalar request = 0.0;
			scalar max = 137.5E-06*1000000;
			//scalar tb = BiomassBufferNew[ idx ].ni+BiomassBufferNew[ idx ].ins;
			scalar tb = copy_output.ni + copy_output.ins;			
			scalar ratio = (float)(max-tb)/(float)max;

			// will this not just be += 0.0??
			request += timek * 0.0 * copy_output.ex_re * copy_output.ins * ratio;

			if(copy_output.ex_re > 0)
				request += timek * ( vDoc/(float) (kDoc+ copy_output.ex_re )) * copy_output.ex_re * copy_output.ni * ratio;

			dtUptk = request;

			//check the requested uptake can be satisfied
			if( copy_output.ex_re < dtUptk )
				dtUptk = copy_output.ex_re;

    			////INSULATION///////
			dtInsul = timek * insul * copy_output.ni;

			//RECYCLING//////////
			// Limiting values
			scalar totBioAvail = copy_output.mb + copy_output.ni + copy_output.ins + dtUptk;
			scalar nonIntRes = copy_output.ni + copy_output.ins;	// this is just tb???

			// upLimMobN_coef+upLimMobI_coef+lowLimTotMob_coef <= 1
			scalar upLimMobN_coef = 0.01;
			scalar upLimMobI_coef = 0.01;
			scalar lowLimTotMob_coef = 0.01;
			scalar upLimMobN = 0;
			scalar upLimMobI = 0;
				
			
			if(copy_output.ni > 0) //BiomassBufferNew[ idx ].ni
			{
	    			upLimMobN = copy_output.ni - dtInsul - upLimMobN_coef * totBioAvail;
			}
			else
			{
				upLimMobN = 0;
			}

			if (copy_output.ins > 0) //BiomassBufferNew[ idx ].ins
			{
				upLimMobI = copy_output.ins + dtInsul - upLimMobI_coef * totBioAvail;
			}
			else
			{
				upLimMobI = 0;
			}

			//scalar lowLimTotMob = -BiomassBufferNew[ idx ].mb - dtUptk + lowLimTotMob_coef*totBioAvail;
			scalar lowLimTotMob = -1 * copy_output.mb - dtUptk + lowLimTotMob_coef * totBioAvail;
					
			// pi_ratio is, in fact, the 'internal ressource concentration'
			//scalar pi_ratio =(float) BiomassBufferNew[ idx ].mb / (float)(BiomassBufferNew[ idx ].ni+BiomassBufferNew[ idx ].ins+0.0000000006);		
			scalar pi_ratio =(float) copy_output.mb / (float)(copy_output.ni+copy_output.ins+0.0000000006);	//copy_output x + y is nonIntRes and tb	???
	
			float alphai=alphan/omega;
			float betai=betan*rho*alphai/alphan; // already calculated at beginning of kernel

			// Request
			if (theta == 2)
			{
				dtNetMobI = timek * copy_output.ins * pi_ratio * (betai-pi_ratio*alphai);	//BiomassBufferNew[ idx ].ins
				dtNetMobN = timek * copy_output.ni * pi_ratio * (betan-pi_ratio*alphan); 	//BiomassBufferNew[ idx ].ni
			}

			scalar dtTotMob = dtNetMobN + dtNetMobI;
				
			// sector dependent projection
			if(   ( dtTotMob < lowLimTotMob ) 
				|| ( dtNetMobN > upLimMobN ) 
				|| ( dtNetMobI > upLimMobI )	)
			{
				if(dtNetMobN > upLimMobN)
				{
					if(dtNetMobI > upLimMobI) 
						{dtNetMobI = upLimMobI; }
					else if(dtTotMob < lowLimTotMob) 
						{dtNetMobI = lowLimTotMob-upLimMobN; }
					else dtNetMobI = (float)upLimMobI-(float)totBioAvail/(float)(dtNetMobN+upLimMobI-lowLimTotMob)*(upLimMobI - dtNetMobI);
					dtNetMobN=upLimMobN;
					
				}
				else if(dtNetMobI > upLimMobI)
				{
					if(dtTotMob < lowLimTotMob) 
						{dtNetMobN = lowLimTotMob-upLimMobI; }
					else dtNetMobN = (float)upLimMobN-(float)totBioAvail/(float)(dtNetMobI+upLimMobN-lowLimTotMob)*(upLimMobN - dtNetMobN);
					dtNetMobI = upLimMobI;
				}
				else if(dtTotMob < lowLimTotMob)
				{
					dtNetMobN = (float)upLimMobN-(float)totBioAvail/(float)(nonIntRes-dtTotMob)*(upLimMobN-dtNetMobN);
					dtNetMobI = (float)upLimMobI-(float)totBioAvail/(float)(nonIntRes-dtTotMob)*(upLimMobI-dtNetMobI);
				}			
			}
		}

		// BUFFER UPDATE /////
		// ni == x
		// ins == y
		// ex_re == z
		// mb == w
		//a temp to store biomass data in to transfer to 3D texture for opixel shader so we can do filtering etc...// rgba xyzw
		struct BiomassBufType biomass;
		biomass.ni = copy_output.ni - dtInsul - dtNetMobN; //BiomassBufferNew[ idx ].ni
		biomass.ins = copy_output.mb + dtUptk + dtNetMobN + dtNetMobI;  //BiomassBufferNew[ idx ].mb
		biomass.ex_re = copy_output.ins + dtInsul- dtNetMobI; //BiomassBufferNew[ idx ].ins
		biomass.mb = copy_output.ex_re - dtUptk; //BiomassBufferNew[ idx ].ex_re

		rd_out[ idx ] = biomass;
		
		//pass contents of unstructured buffer into 3D texture for rendering
		//uint3 cellID = threadIDToCell(DTid, Gid );
		//uint3 output;
		//output.x = Gid.x*DTid.x;
		//output.y = DTid.y*Gid.y ;
		//output.z = DTid.z*Gid.z;		//output isnt used

 		//int s= soil_struct.Load(sampleCoords); // s isnt used
		//float4 test;
		//if(s==0)
		//test=float4(1,0,0,1);
		//else if (s==255)
		//{
		//test=float4(1, 1, 1, 0.003);
		//}					//test isnt used

		//gOutputTex[cellID].rgba = biomass.rgba;
		// do we need the outputtex if we just saved biomass to the output buffer rd_out
	}
}