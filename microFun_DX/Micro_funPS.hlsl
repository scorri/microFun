/*+==============================================================+\
*                                                                 *
*  Copyright (c) Microsoft Corporation.  All rights reserved.     *
*                                                                 *
*  File: HandsOnLab_ConwayPS.hlsl                                 *
*                                                                 *
*  Source code for the Conway pixel shader                        *
*                                                                 *
\+==============================================================+*/
 // NOTE: must be identical to "uThreadsY" in "HandsOnLab.h"

#include "constants.h"
#include "common.h"



// definition of the Conway cell buffer element
struct BiomassBufType
{
    scalar ni;
    scalar mb;
    scalar ins;
    scalar ex_re;
};


// a read-only view of a Biomass cell buffer resource
StructuredBuffer<BiomassBufType> BiomassBuffer: register (t0);
Texture3D<float4> gTex3D: register (t1);
Texture3D<uint> gTexStruct: register (t2);
SamplerState Texsampler: register (s0);
//--------------------------------------------------------------------------------------
// ConwayPS
// main entry point
//--------------------------------------------------------------------------------------

struct PS_INPUT
{
    float4 Pos : SV_POSITION;
    float4 Tex : TEXCOORD0;
};


float4 Micro_funPS( PS_INPUT inPix ) : SV_Target
{
	float4 return_col;
	bool g=false;
#if 0
	float4 s= gTex3D.Sample(Texsampler, inPix.Tex);
	if (s.r > 0.0)
	{
		return s;
	}
	return float4(0,1,0,0.001);
#else
	if(geometry)
	{
	// calculate the buffer index corresponding to this screen coordinate
    uint idx = (inPix.Pos.x-.5)+(inPix.Pos.y-.5)*BUFFER_SIZE_X+(inPix.Pos.z-.5)*BUFFER_SIZE_Y*BUFFER_SIZE_Z;  
	//const float4 bias=	float4(-0.5, -0.5, -0.5, 0);
	//const float4 stride=float4(1, BUFFER_SIZE_X, BUFFER_SIZE_X*BUFFER_SIZE_Y, 0);
	//const uint idx= dot(inPix.Tex, stride);
	float4 s= gTex3D.Sample(Texsampler, inPix.Tex);
	
	//uint4 sampleCoords = uint4(inPix.Tex,0);
	int  st = gTexStruct.Load(inPix.Tex);
	float temp_structure=0.0;
	if(st>0)
	temp_structure = 0.0;
	else
	temp_structure = 1.0;
	
    //max-min biomass values
	float normW= 1.0/(s.r+s.g+s.b);

		// [branch]
	if ((s.r+s.g+s.b) > 0.0)
	{
		// show dominant phase
		if ((s.r > s.g)&&(s.r > s.b))
		return_col = float4(1-s.r*1,0,0, s.r*normW);
			
		if ((s.g > s.r)&&(s.g > s.b))
		return_col = float4(0,1-s.g*1,0,s.g*normW);

		if ((s.b > s.r)&&(s.b > s.g))
		return_col = float4(0,0,1-s.b*1,s.b*normW);
	
	     //set the cell color to correspond to the cell's state
		return_col = float4(normW*s.r,normW*s.g,normW*s.b*r_ins,(normW*s.r+normW*s.g+normW*s.b));
		
	}
	else
	{
		return_col = float4(1, 1, 1, 0.0002);
	}
	
	}
	else
	{

		 return_col= float4(0, 0, 0,1);
		 
	}
	
			
	return return_col;// float4((r_ni*ni*weight),(r_mb*mb*weight),(r_ins*ins*weight)*10000,1.0);
#endif
}
