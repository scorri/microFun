
// microfun parameters
#define timek 0.0020833
//fungal genotype parameters
#define insul 0.0000000000001 //0.5 
#define betan 0.2
#define alphan 0.8
#define rho 1.2 //less than 1 for recycling
#define omega 100
#define replenishment 0
#define h 0.05
#define kDoc 0.000001
#define vDoc 5

// mouse constant buffer initial values
#define Db 0.18
#define Dv 0.13
#define theta 2.0

#pragma OPENCL EXTENSION cl_khr_fp64 : enable
typedef double scalar;

typedef struct tag_SharedInfo
{
	scalar ni;
	scalar mb;
	int structure;
}SharedInfo;

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

// output kernel for BIOMASSBUFTYPE rd buffer
__kernel void micro_output(__global float4 *colour, 
							__global struct BiomassBufType *rd_buff)
{

	int idx = get_global_id(0);
	struct BiomassBufType rd_colour = rd_buff[idx];

	float4 c = (float4)(0, 0, 0, 0);
	scalar norm = rd_colour.ni + rd_colour.ins + rd_colour.mb;
		if(rd_colour.ni + rd_colour.ins > 0.05f)
		{   
			c = (float4)((norm -rd_colour.ni)/norm, (norm -rd_colour.ins)/norm, (norm -rd_colour.mb)/norm, 1);
		}

	colour[idx] = c;
}

// output kernel for BIOMASSBUFTYPE rd buffer
__kernel void micro_structure(__global float4 *colour, 
							__global int *structure_buffer)
{
	int idx = get_global_id(0);
	int structure = structure_buffer[idx];

	float4 c = (float4)(0, 0, 0, 0);

	if(structure > 0)
		c = (float4)(1, 0, 0, 1);
	

	colour[idx] = c;
}

// unroll
void microwarpReduce(__local volatile struct BiomassBufType* s_data, int tid)
{
	//s_data[tid] += s_data[tid+32];
	s_data[tid].ni = s_data[tid].ni + s_data[tid+32].ni;
	s_data[tid].mb = s_data[tid].mb + s_data[tid+32].mb;
	s_data[tid].ins = s_data[tid].ins + s_data[tid+32].ins;
	s_data[tid].ex_re = s_data[tid].ex_re + s_data[tid+32].ex_re;

	//s_data[tid] += s_data[tid+16];
	s_data[tid].ni = s_data[tid].ni + s_data[tid+16].ni;
	s_data[tid].mb = s_data[tid].mb + s_data[tid+16].mb;
	s_data[tid].ins = s_data[tid].ins + s_data[tid+16].ins;
	s_data[tid].ex_re = s_data[tid].ex_re + s_data[tid+16].ex_re;

	//s_data[tid] += s_data[tid+ 8];
	s_data[tid].ni = s_data[tid].ni + s_data[tid+8].ni;
	s_data[tid].mb = s_data[tid].mb + s_data[tid+8].mb;
	s_data[tid].ins = s_data[tid].ins + s_data[tid+8].ins;
	s_data[tid].ex_re = s_data[tid].ex_re + s_data[tid+8].ex_re;

	//s_data[tid] += s_data[tid+ 4];
	s_data[tid].ni = s_data[tid].ni + s_data[tid+4].ni;
	s_data[tid].mb = s_data[tid].mb + s_data[tid+4].mb;
	s_data[tid].ins = s_data[tid].ins + s_data[tid+4].ins;
	s_data[tid].ex_re = s_data[tid].ex_re + s_data[tid+4].ex_re;

	//s_data[tid] += s_data[tid+ 2];
	s_data[tid].ni = s_data[tid].ni + s_data[tid+2].ni;
	s_data[tid].mb = s_data[tid].mb + s_data[tid+2].mb;
	s_data[tid].ins = s_data[tid].ins + s_data[tid+2].ins;
	s_data[tid].ex_re = s_data[tid].ex_re + s_data[tid+2].ex_re;

	//s_data[tid] += s_data[tid+ 1];	
	s_data[tid].ni = s_data[tid].ni + s_data[tid+1].ni;
	s_data[tid].mb = s_data[tid].mb + s_data[tid+1].mb;
	s_data[tid].ins = s_data[tid].ins + s_data[tid+1].ins;
	s_data[tid].ex_re = s_data[tid].ex_re + s_data[tid+1].ex_re;
}

// half gws
// unroll last warp
__kernel 
void micro_sum(__global struct BiomassBufType *rd_in,
		__local struct BiomassBufType *s_data,
		__global struct BiomassBufType *sum_partial)
{
	unsigned int lid = get_local_id(0);

	unsigned int gid = get_group_id(0)*(get_local_size(0)*2) + lid;

	// store global value into local shared memory address
	s_data[lid].ni = rd_in[gid].ni + rd_in[gid + get_local_size(0)].ni;
	s_data[lid].mb = rd_in[gid].mb + rd_in[gid + get_local_size(0)].mb;
	s_data[lid].ins = rd_in[gid].ins + rd_in[gid + get_local_size(0)].ins;
	s_data[lid].ex_re = rd_in[gid].ex_re + rd_in[gid + get_local_size(0)].ex_re;

	barrier(CLK_LOCAL_MEM_FENCE); /* Wait for others in the work-group */

	// do reduction to find sum of work group
	for(unsigned int s = get_local_size(0)/2; s > 32; s >>= 1)
	{
		if(lid < s)
		{
			s_data[lid].ni = s_data[lid].ni + s_data[lid + s].ni;
			s_data[lid].mb = s_data[lid].mb + s_data[lid + s].mb;
			s_data[lid].ins = s_data[lid].ins + s_data[lid + s].ins;
			s_data[lid].ex_re = s_data[lid].ex_re + s_data[lid + s].ex_re;
		}

		barrier(CLK_LOCAL_MEM_FENCE); /* Wait for others in the work-group */	
	}

	// unroll last warp
	if(lid < 32) microwarpReduce(s_data, lid);

	// write partial result for this work group
	if(lid == 0)
	{
		sum_partial[get_group_id(0)] = s_data[0];
	}
}



__kernel 
void micro_diffuse(__global struct BiomassBufType *rd_in,
			int4 params,
		__local SharedInfo *s_data,
		__global struct BiomassBufType *rd_out,
		__global int *soil_struct)
{

	float alphai=alphan/omega;
	float betai=betan*rho*alphai/alphan;


	int3 l_size = (int3)(get_local_size(0), get_local_size(1), get_local_size(2));
	int3 g_id = (int3)(get_group_id(0), get_group_id(1), get_group_id(2));
	int3 l_id = (int3)(get_local_id(0), get_local_id(1), get_local_id(2));
	int3 w_dim = (int3)(params.xyz); 
	int3 dt_id = (int3)(l_size*g_id + l_id); //Dispatch Thread ID (DC)

	int idx = getBufferIndexUseDispatch(dt_id, g_id, w_dim);

	int s = soil_struct[idx];
	int temp_structure=0;
	if(s>0)
		temp_structure = 1;
	else
		temp_structure = 0;

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
	s_data[local_id].mb = rd_in[idx].mb;
	s_data[local_id].structure = temp_structure;

	barrier(CLK_LOCAL_MEM_FENCE); /* Wait for others in the work-group */


	// updates ni and mb in this case x and y from rd_in to rd_out
	bool inWindow = threadInWindow(dt_id, g_id, w_dim);
	bool onEdge = threadOnGroupEdge( l_id, l_size);
	if ( inWindow  && !onEdge )
    {
		// calculate the sum of the adjacent cell's states
		SharedInfo sum;
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
		SharedInfo state = s_data[ getOffset(local_id, 0, 0, l_size.x, 0, l_size.y)];

        // STRUCTURE CONSTRAINED DIFFUSION _____  diffuse NI -  diffusion algorithm from GPU gems 2 chapter
		// diffuse ni biomass -  diffusion algorithm from GPU gems 2 chapter
		scalar coeff =( 1-(timek*sum.structure*Db)/(h*h) );
        scalar nextStateNI = state.ni * coeff + (timek/(h*h))*state.structure*Db*sum.ni;

        // finally, update the new biomass buffer resource with the next-state value determined
		//rd_out[ idx ].ni = nextStateNI;

		//float coeff_out= 1- (timek*sum.structure*Dv)/h*h;
		// diffuse mobile biomass -  diffusion algorithm from GPU gems 2 chapter
		coeff = 1-((timek*sum.structure*Dv)/(h*h));
        scalar nextStatemb = state.mb * coeff + (timek/(h*h))*state.structure*Dv * sum.mb;
 
		//IS MOBILE BIOMASS DIFFUING
		// finally, update the new biomass buffer resource with the next-state value determined
		//rd_out[ idx ].mb = nextStatemb*(state.ni >0);

		struct BiomassBufType output;
		output.ni = nextStateNI;
		output.mb = nextStatemb*(state.ni >0);
		output.ins = rd_in[idx].ins;
		output.ex_re = rd_in[idx].ex_re;

		rd_out[idx] = output;
	}
}

__kernel 
void micro_update(__global struct BiomassBufType *rd_in,
		__global struct BiomassBufType *rd_out)
{

	float alphai=alphan/omega;
	float betai=betan*rho*alphai/alphan;

	// could just reduce dimensions of kernel?
	//int idx = get_global_id(0) + get_global_id(1)*get_global_size(0) + get_global_id(2)*get_global_size(0)*get_global_size(1);
	int idx = get_global_id(0); //1D

	// create a variables for storing the current cell's next-state value based on recycling
	scalar dtInsul=0;
    scalar dtNetMobI=0;
    scalar dtNetMobN=0;
    scalar dtUptk=0;
    scalar dtTotMob=0;

	struct BiomassBufType copyBufferNew = rd_in[ idx ];
	if( (copyBufferNew.ni != 0) || (copyBufferNew.ins != 0) )
	{
		/////////UPTAKE//////////////////////
		// Request
		scalar request=0.0;
		scalar max=137.5E-06*1000000;
		scalar tb=copyBufferNew.ni+copyBufferNew.ins;
		scalar ratio =(max-tb)/max;

		if(copyBufferNew.ex_re>0)
			request += timek*( vDoc/(float)(kDoc+copyBufferNew.ex_re))*copyBufferNew.ex_re*copyBufferNew.ni*ratio;

		 dtUptk = request;

		//check the requested uptake can be satisfied
		if(copyBufferNew.ex_re<dtUptk)
		{
			dtUptk=copyBufferNew.ex_re;
		} 

		/////////INSULATION////
		dtInsul = timek*insul*copyBufferNew.ni;

		// Limiting values to ensure all state varaibles remain positive
		float totBioAvail=copyBufferNew.mb+copyBufferNew.ins+copyBufferNew.ni+dtUptk;
    	float nonIntRes=copyBufferNew.ins+copyBufferNew.ni;

		// upLimMobN_coef+upLimMobI_coef+lowLimTotMob_coef <= 1
		scalar upLimMobN_coef=0.01;
		scalar upLimMobI_coef=0.01;
		scalar lowLimTotMob_coef=0.01;
 		scalar upLimMobN=0;
		scalar upLimMobI=0;

        if(copyBufferNew.ni>0)
		{
    		upLimMobN = copyBufferNew.ni - dtInsul - upLimMobN_coef*totBioAvail;
		}
		else
		{
			upLimMobN=0;
		}

		if (copyBufferNew.ins >0)
		{
			upLimMobI = copyBufferNew.ins + dtInsul - upLimMobI_coef*totBioAvail;
		}
		else
		{
			upLimMobI = 0;
		}

		scalar lowLimTotMob = -copyBufferNew.mb - dtUptk + lowLimTotMob_coef*totBioAvail;
 			
		// pi_ratio is, in fact, the 'internal ressource concentration'
		scalar pi_ratio = copyBufferNew.mb / (copyBufferNew.ni + copyBufferNew.ins + 0.0000000006);
  			
		// Request
		if (theta==2)
		{
			dtNetMobI = timek*copyBufferNew.ins*pi_ratio * (betai-pi_ratio*alphai);
			dtNetMobN = timek*copyBufferNew.ni*pi_ratio * (betan-pi_ratio*alphan);
		}     

		scalar dtTotMob=dtNetMobN+dtNetMobI;

		// sector dependent projection
		if(dtTotMob<lowLimTotMob 
				|| dtNetMobN>upLimMobN 
				|| dtNetMobI>upLimMobI)
		{
			if(dtNetMobN>upLimMobN)
			{
				if(dtNetMobI>upLimMobI) {dtNetMobI=upLimMobI; }
				else if(dtTotMob<lowLimTotMob) {dtNetMobI=lowLimTotMob-upLimMobN; }
				else dtNetMobI=upLimMobI-totBioAvail/(dtNetMobN+upLimMobI-lowLimTotMob)*(upLimMobI-dtNetMobI);
					dtNetMobN=upLimMobN;
					
			}
			else if(dtNetMobI>upLimMobI)
			{
				if(dtTotMob<lowLimTotMob) {dtNetMobN=lowLimTotMob-upLimMobI; }
				else dtNetMobN=upLimMobN-totBioAvail/(dtNetMobI+upLimMobN-lowLimTotMob)*(upLimMobN-dtNetMobN);
				dtNetMobI=upLimMobI;
				
			}
			else if(dtTotMob<lowLimTotMob)
			{
				dtNetMobN=upLimMobN-totBioAvail/(nonIntRes-dtTotMob)*(upLimMobN-dtNetMobN);
				dtNetMobI=upLimMobI-totBioAvail/(nonIntRes-dtTotMob)*(upLimMobI-dtNetMobI);
			}
		}
	}
 

	////////////////////////////////////////////////////////////////////BUFFER UPDATE//////////////////////////////////////////////////////////////////////
	//a temp to store biomass data in to transfer to 3D texture for opixel shader so we can do filtering etc...
	struct BiomassBufType biomass;
	biomass.ins=copyBufferNew.ins + dtInsul - dtNetMobI;
	biomass.ni=copyBufferNew.ni - dtInsul - dtNetMobN;
	biomass.mb=copyBufferNew.mb + dtUptk + dtNetMobN+dtNetMobI;
	biomass.ex_re=copyBufferNew.ex_re - dtUptk;

	rd_out[ idx ] = biomass;

}