/*+==============================================================+\
*                                                                 *
*																 *
*                                                                 *
*  File: Modified from HandsOnLab_ConwayCS.hlsl (Copyright (c) Microsoft Corporation.  All rights reserved.)                                *
*                                                                 *
*  Source code for the Conway compute shader                      *
*                                                                 *
\+==============================================================+*/
#include "constants.h"

//numerical parameter for diffusion and other processes
#define timek 0.0020833
//fungal genotype parameters
#define insul 0 
#define betan 0.05
#define alphan 0.95
#define rho 1.2 //less than 1 for recycling
#define omega 100
#define replenishment 0
#define h 0.05
#define kDoc 0.0000001
#define vDoc 5
#define uptake_trait 1.0 //check this


struct DebugProcessBufType
{
	scalar totUptk;
	scalar totIns;
	scalar totMobN;
	scalar totMobI;
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

// Info stored in the groupshared memory, corresponds to the info in BioMassBufType - maximum is 16 or 32 kilobytes - 4 bytes per float
struct SharedInfo
{
	scalar ni;
	scalar mb;
	int structure;
};
// groupshared memory declaration
groupshared SharedInfo sharedState[THREADS_X*THREADS_Y*THREADS_Z];
// definition of the mouse info and paremeters that may modified by keyboard constant buffer
cbuffer CB0
{
    uint mousex;
    uint mousey;
    uint mousebtn;
    float Db;
    //diffussion coefficient of biomass modified via '5' & '6' key presses
	float Dv;
    //diffussion coefficient of internal resource 
	float theta;
    // immobilisation term
};

// a read-only view of the first biomass buffer resource
StructuredBuffer<BiomassBufType> BiomassBufferOld : register (t0);
// a read-only view of the 3D texture buffer resource
Texture3D<uint> soil_struct : register (t1);

//StructuredBuffer<uint> SoilStructure : register (t1);
// a read-write view of the  biomass buffer resource
RWStructuredBuffer<BiomassBufType> BiomassBufferNew: register (u0);
//3D texture - without this specifier a memory barrier or sync will flush the UAV only within the current group
globallycoherent RWTexture3D<float4> gOutputTex: register (u1); 
// a read-write view of the  biomass buffer resource
RWStructuredBuffer<DebugProcessBufType> Processes: register (u2); //u2 specifies the slot it is bound to?
//create a read only vie of the debug processing buffer
//RWStructuredBuffer<DebugProcessBufType> Processes; 

//--------------------------------------------------------------------------------------
// threadOnSealingPlanes
// returns true if the cell corresponding to thread DTid in group Gid is in the application window/volume
// Due to overlapping thread volumes threads on borders are actually outside simulation bounds
//--------------------------------------------------------------------------------------
bool threadInSealingPlanes( uint3 DTid, uint3 Gid )
{

}

//--------------------------------------------------------------------------------------
// threadInWindow
// returns true if the cell corresponding to thread DTid in group Gid is in the application window/volume
// Due to overlapping thread volumes threads on borders are actually outside simulation bounds
//--------------------------------------------------------------------------------------
bool threadInWindow( uint3 DTid, uint3 Gid )
{
    bool inWindowX = DTid.x - 2 * Gid.x - 1 < BUFFER_SIZE_X && DTid.x > 0;
    bool inWindowY = DTid.y - 2 * Gid.y - 1 < BUFFER_SIZE_Y && DTid.y > 0;
    bool inWindowZ = DTid.z - 2 * Gid.z - 1 < BUFFER_SIZE_Z && DTid.z > 0;
    // Check in z?
	return inWindowX && inWindowY && inWindowZ;
}


//--------------------------------------------------------------------------------------
// threadOnGroupEdge
// returns true if the group thread GTid is part of the overlapping edge region
//--------------------------------------------------------------------------------------
bool threadOnGroupEdge( uint3 GTid )
{
    return GTid.x == 0				|| 
		   GTid.x == THREADS_X - 1	||
		   GTid.y == 0				||
		   GTid.y == THREADS_Y - 1 	||
		   GTid.z == 0				||
		   GTid.z == THREADS_Z - 1;
}


//--------------------------------------------------------------------------------------
// threadIDToCell
// given a thread DTid and group Gid, calculate the corresponding cell position
//--------------------------------------------------------------------------------------
uint3 threadIDToCell( uint3 DTid, uint3 Gid )
{
    uint3 output;
    output.x = DTid.x - 2 * Gid.x - 1;
    output.y = DTid.y - 2 * Gid.y - 1;
    output.z = DTid.z - 2 * Gid.z - 1;
    return output;
}


//--------------------------------------------------------------------------------------
// getBufferIndex
// given a thread DTid and group Gid, calculate the corresponding cell index in the buffer
//does wrapping around i.e. id DTID & GID is 0,0,0 then returns 31, 31, 31
//we can create the buffer index associations in such a way as to ensure that no cross-group memory reads are required after group-shared memory has been updated,
// and that every cell has exactly one thread that will determine its next state.
//--------------------------------------------------------------------------------------
uint getBufferIndex( uint3 DTid, uint3 Gid )
{
    uint3 wrapped;
    wrapped.x = ( DTid.x - 2 * Gid.x - 1  + BUFFER_SIZE_X ) % BUFFER_SIZE_X;
    wrapped.y = ( DTid.y - 2 * Gid.y - 1  + BUFFER_SIZE_Y ) % BUFFER_SIZE_Y;
    wrapped.z = ( DTid.z - 2 * Gid.z - 1            + BUFFER_SIZE_Z ) % BUFFER_SIZE_Z;
    return wrapped.x + (wrapped.y * BUFFER_SIZE_X) + (wrapped.z * BUFFER_SIZE_X * BUFFER_SIZE_Y);
}


//--------------------------------------------------------------------------------------
// isNear
// returns true if two points are identical or immediately adjacent
//--------------------------------------------------------------------------------------
bool isNear( uint2 a, uint3 b )
{
    return  (
				a.x == b.x		||
				a.x == b.x + 1	||
				a.x == b.x - 1
			) && (
				a.y == b.y		||
				a.y == b.y + 1	||
				a.y == b.y - 1
			);
	// TODO - for 3D point selection we're gonna need some raycasting or other
	// Z inferring system
     /*&& (
				a.z == b.z		||
				a.z == b.z + 1	||
				a.z == b.z -1 
				);*/
}


//--------------------------------------------------------------------------------------
// getCellByOffset
// given a shared memory group index GI and an X and Y offset, 
// returns the state of the cell at that offset relative to the current cell
// Only non insulating and mobile biomass info is returned as per SharedInfo definition
//--------------------------------------------------------------------------------------
SharedInfo getCellByOffset( uint GI, int xOffset, int yOffset, int zOffset )
{
    return sharedState[ GI + xOffset + (yOffset * THREADS_X) + (zOffset * THREADS_X * THREADS_Y) ];
}


//--------------------------------------------------------------------------------------
// getAdjacentCellSum
// given a shared memory group index GI, calculates the sum of the states of
//     all adjacent cells, ni and mb only, as described for getCellByOffset
//--------------------------------------------------------------------------------------
SharedInfo getAdjacentCellSum( uint GI)
{
    SharedInfo sum;
	sum.ni = 
			getCellByOffset(GI, 0, -1, 0).ni +
			getCellByOffset(GI, -1, 0, 0).ni +
			getCellByOffset(GI, 1, 0, 0).ni +
			getCellByOffset(GI, 0, 1, 0).ni +
			getCellByOffset(GI, 0, 0, -1).ni +
			getCellByOffset(GI, 0, 0, 1).ni;
	sum.mb = 
			getCellByOffset(GI, 0, -1, 0).mb +
			getCellByOffset(GI, -1, 0, 0).mb +
			getCellByOffset(GI, 1, 0, 0).mb +
			getCellByOffset(GI, 0, 1, 0).mb +
			getCellByOffset(GI, 0, 0, -1).mb +
			getCellByOffset(GI, 0, 0, 1).mb;
	sum.structure =
			getCellByOffset(GI, 0, -1, 0).structure +
			getCellByOffset(GI, -1, 0, 0).structure +
			getCellByOffset(GI, 1, 0, 0).structure +
			getCellByOffset(GI, 0, 1, 0).structure +
			getCellByOffset(GI, 0, 0, -1).structure +
			getCellByOffset(GI, 0, 0, 1).structure;

	return sum;
}


//--------------------------------------------------------------------------------------
//
// main entry point
//--------------------------------------------------------------------------------------


// run each thread group with dimensions THREADS_X * THREADS_Y * THREAD_Z
[numthreads(THREADS_X,THREADS_Y,THREADS_Z)]
void ActivityCS( uint3 Gid : SV_GroupID, uint3 DTid : SV_DispatchThreadID, uint3 GTid : SV_GroupThreadID, uint GI : SV_GroupIndex )
{


	float alphai=alphan/omega;
	float betai=betan*rho*alphai/alphan;

	//a temp to store biomass data in to transfer to 3D texture for opixel shader so we can do filtering etc...
	float4 biomass;

    // calculate the buffer index corresponding to the current thread
	uint idx = getBufferIndex( DTid, Gid );
	
	uint3 cellID = threadIDToCell(DTid, Gid );	

	
	uint4 sampleCoords = uint4(cellID,0);

	int s= soil_struct.Load(sampleCoords);
	int temp_structure=0;
	if(s>0)
	temp_structure = 0;
	else
	temp_structure = 1;

	//check if thread is on sealing planes and if so set to solid

	 if( DTid.x - 2 * Gid.x - 1==0)
	 temp_structure=0;
	 if( DTid.x - 2 * Gid.x - 1==BUFFER_SIZE_X-1)
	 temp_structure=0;

	 if( DTid.y - 2 * Gid.y - 1==0)
	 temp_structure=0;
	 if( DTid.y - 2 * Gid.y - 1==BUFFER_SIZE_Y-1)
	 temp_structure=0;

	 if( DTid.z - 2 * Gid.z - 1==0)
	 temp_structure=0;
	 if( DTid.z - 2 * Gid.z - 1==BUFFER_SIZE_Z-1)
	 temp_structure=0;

    // Store the state of the current cell in shared memory for other threads in this group to access
	// NOTE: We store both ni and mb state in one pass
	sharedState[GI].ni = BiomassBufferOld[ idx ].ni;
    sharedState[GI].mb = BiomassBufferOld[ idx ].mb;
	sharedState[GI].structure = temp_structure;

	
    // wait until all threads in this group reach this point to ensure that sharedState contains accurate data
	GroupMemoryBarrierWithGroupSync();

    // only continue if the current cell is within the window and not in an overlapping region
	if( threadInWindow( DTid, Gid ) && !threadOnGroupEdge( GTid ) )
	{
        // create a variable for storing the current cell's next-state value based on diffusion
		float nextStateNI;

        // calculate the sum of the adjacent cell's states
		SharedInfo sum = getAdjacentCellSum( GI );

		// get current state of this threads cell
        SharedInfo state = getCellByOffset(GI, 0, 0, 0);
		SharedInfo state1 = getCellByOffset(GI, 1, 0, 0);
		SharedInfo state2 = getCellByOffset(GI, -1, 0, 0);
		SharedInfo state3 = getCellByOffset(GI, 0, 1, 0);
		SharedInfo state4 = getCellByOffset(GI, 0, -1, 0);
		SharedInfo state5 = getCellByOffset(GI, 0, 0, 1);
		SharedInfo state6 = getCellByOffset(GI, 0, 0, -1);
		


        // STRUCTURE CONSTRAINED DIFFUSION _____  diffuse NI -  diffusion algorithm from GPU gems 2 chapter
		// diffuse ni biomass -  diffusion algorithm from GPU gems 2 chapter
		scalar coeff =( 1-(timek*sum.structure*Db)/(h*h) );
        nextStateNI = state.ni * coeff + (timek/(h*h))*state.structure*Db*sum.ni;

		 
		
		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		uint4 sampleCoords = uint4(threadIDToCell(DTid, Gid ),0);
 		int s= soil_struct.Load(sampleCoords);
		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


        // finally, update the new biomass buffer resource with the next-state value determined
		BiomassBufferNew[ idx ].ni = nextStateNI;

		//float coeff_out= 1- (timek*sum.structure*Dv)/h*h;
		// diffuse mobile biomass -  diffusion algorithm from GPU gems 2 chapter
		coeff = 1-((timek*sum.structure*Dv)/(float)(h*h));
        scalar nextStatemb = state.mb * coeff + (timek/(float)(h*h))*state.structure*Dv * sum.mb;
 
        
		//IS MOBILE BIOMASS DIFFUING
		// finally, update the new biomass buffer resource with the next-state value determined
		BiomassBufferNew[ idx ].mb = nextStatemb*(state.ni >0);
   }

	GroupMemoryBarrierWithGroupSync();

  if( threadInWindow( DTid, Gid ) && !threadOnGroupEdge( GTid ) )
	{  
	    // create a variables for storing the current cell's next-state value based on recycling
		scalar dtInsul=0;
        scalar dtNetMobI=0;
        scalar dtNetMobN=0;
        scalar dtUptk=0;
        scalar dtTotMob=0;

        //assuming a replenished env based on parameter replenishment unit per computational timestep
		scalar current_re=BiomassBufferNew[ idx ].ex_re;
        BiomassBufferNew[ idx ].ex_re = replenishment*(265.2500e-06*1000000-current_re)+current_re;
		
		////////////////////////////////////////////////////////////////////UPTAKE//////////////////////////////////////////////////////////////////////
		// Request
		scalar request=0.0;
		scalar max=137.5E-06*1000000;
		scalar tb=BiomassBufferNew[ idx ].ni+BiomassBufferNew[ idx ].ins;
		scalar ratio =(float)(max-tb)/(float)max;

		request += timek*0.0*BiomassBufferNew[ idx ].ex_re*BiomassBufferNew[ idx ].ins*ratio;

		if(BiomassBufferNew[ idx ].ex_re>0)
		request += timek*( vDoc/(float)(kDoc+BiomassBufferNew[ idx ].ex_re))*BiomassBufferNew[ idx ].ex_re*BiomassBufferNew[ idx ].ni*ratio;

		//request += timek*(ldan+vDOC/(kDOC+DOC_vxl))*DOC_vxl*nonInsBio_vxl*ratio;
		 dtUptk = request;
				 

		//check the requested uptake can be satisfied
        if(BiomassBufferNew[ idx ].ex_re<dtUptk)
		{
            dtUptk=BiomassBufferNew[ idx ].ex_re;
        } 
		Processes[0].totUptk=100;//+=dtUptk;
		////////////////////////////////////////////////////////////////////INSULATION//////////////////////////////////////////////////////////////////////	
        if((BiomassBufferNew[ idx ].ni>0))//&&(insul>0))
		{
            dtInsul = timek*insul*BiomassBufferNew[ idx ].ni;
        }
		
		////////////////////////////////////////////////////////////////////RECYCLING//////////////////////////////////////////////////////////////////////
		 if((BiomassBufferNew[ idx ].ni>0)||(BiomassBufferNew[ idx ].ins>0))
		 {
          // Limiting values to ensure all state varaibles remain positive
			float totBioAvail=BiomassBufferNew[ idx ].mb+BiomassBufferNew[ idx ].ins+BiomassBufferNew[ idx ].ni+dtUptk;
       		float nonIntRes=BiomassBufferNew[ idx ].ins+BiomassBufferNew[ idx ].ni;

			/////////////////// MODIFICATIONS ///////////////////////////
			
			// upLimMobN_coef+upLimMobI_coef+lowLimTotMob_coef <= 1
			scalar upLimMobN_coef=0.05;
			scalar upLimMobI_coef=0.05;
			scalar lowLimTotMob_coef=0.05;
 
            scalar upLimMobN = BiomassBufferNew[ idx ].ni - dtInsul - upLimMobN_coef*totBioAvail;
            scalar upLimMobI = BiomassBufferNew[ idx ].ins + dtInsul - upLimMobI_coef*totBioAvail;
            scalar lowLimTotMob = -BiomassBufferNew[ idx ].mb - dtUptk + lowLimTotMob_coef*totBioAvail;
            // pi_ratio is, in fact, the 'internal ressource concentration'
			float inverse, pi_ratio;
			
			if(BiomassBufferNew[ idx ].ni+BiomassBufferNew[ idx ].ins>0)
			{
			 inverse=1.0/(float)(BiomassBufferNew[ idx ].ni+BiomassBufferNew[ idx ].ins);
			 pi_ratio = BiomassBufferNew[ idx ].mb*inverse;
			}
			else
			{
				pi_ratio=0;
			}
            // Recycling Request
			if (theta==3)
			{
                dtNetMobI = timek*BiomassBufferNew[ idx ].ins*pi_ratio * (betai-pi_ratio*pi_ratio*alphai);
                dtNetMobN = timek*BiomassBufferNew[ idx ].ni*pi_ratio * (betan-pi_ratio*pi_ratio*alphan);
            }
				if (theta==2)
			{
                dtNetMobI = timek*BiomassBufferNew[ idx ].ins*pi_ratio * (betai-pi_ratio*alphai);
                dtNetMobN = timek*BiomassBufferNew[ idx ].ni*pi_ratio * (betan-pi_ratio*alphan);
            }
						
			dtTotMob = dtNetMobI+dtNetMobN;

				 if(dtTotMob<lowLimTotMob || dtNetMobN>upLimMobN || dtNetMobI>upLimMobI)
				 {
					if(dtNetMobN>upLimMobN)
					{
						if(dtNetMobI>upLimMobI) {dtNetMobI=upLimMobI; /*excess_mobI++;*/}
						else if(dtTotMob<lowLimTotMob) {dtNetMobI=lowLimTotMob-upLimMobN; /*excess_immob++;*/}
						else dtNetMobI=upLimMobI-totBioAvail/(float)(dtNetMobN+upLimMobI-lowLimTotMob)*(upLimMobI-dtNetMobI);
						dtNetMobN=upLimMobN;
						//excess_mobN++;
					}
					else if(dtNetMobI>upLimMobI)
					{
						if(dtTotMob<lowLimTotMob) {dtNetMobN=lowLimTotMob-upLimMobI; /*excess_immob++;*/}
						else dtNetMobN=upLimMobN-totBioAvail/(float)(dtNetMobI+upLimMobN-lowLimTotMob)*(upLimMobN-dtNetMobN);
						dtNetMobI=upLimMobI;
						//excess_mobI++;
					}
					else if(dtTotMob<lowLimTotMob)
					{
						dtNetMobN=upLimMobN-totBioAvail/(float)(nonIntRes-dtTotMob)*(upLimMobN-dtNetMobN);
						dtNetMobI=upLimMobI-totBioAvail/(float)(nonIntRes-dtTotMob)*(upLimMobI-dtNetMobI);
						//excess_immob++;
					}
				}
				 Processes[0].totMobN= 10000.0;//+=dtNetMobN;
			     Processes[0].totMobI=666.66;//+=dtNetMobI;
					
			  //   
		}  

	////////////////////////////////////////////////////////////////////BUFFER UPDATE//////////////////////////////////////////////////////////////////////
		BiomassBufferNew[ idx ].ins = BiomassBufferNew[ idx ].ins + dtInsul- dtNetMobI;
		biomass.b=BiomassBufferNew[ idx ].ins ;
        BiomassBufferNew[ idx ].ni =  BiomassBufferNew[ idx ].ni - dtInsul - dtNetMobN;
		biomass.r=BiomassBufferNew[ idx ].ni;
        BiomassBufferNew[ idx ].mb =  BiomassBufferNew[ idx ].mb + dtUptk + dtNetMobN+dtNetMobI;
		biomass.g= BiomassBufferNew[ idx ].mb ;			
		BiomassBufferNew[ idx ].ex_re =BiomassBufferNew[ idx ].ex_re - dtUptk;
		biomass.a=BiomassBufferNew[ idx ].ex_re ;

		

		//pass contents of unstructured buffer into 3D texture for rendering
		uint3 cellID = threadIDToCell(DTid, Gid );
		 uint3 output;
		output.x = Gid.x*DTid.x;
		output.y = DTid.y*Gid.y ;
		output.z = DTid.z*Gid.z;
			
 		int s= soil_struct.Load(sampleCoords);
		float4 test;
		if(s==0)
		test=float4(1,0,0,1);
		else if (s==255)
		{
		test=float4(1, 1, 1, 0.003);
		}
		
		gOutputTex[cellID].rgba = biomass.rgba;			
	
	}//end if threadinwindow
}
