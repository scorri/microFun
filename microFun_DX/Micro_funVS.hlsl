//--------------------------------------------------------------------------------------
// File: WORKING   PDCWorkshopCS11.hlsl
//
// Shaders for rendering
//
// Copyright (c) Microsoft Corporation. All rights reserved.
//--------------------------------------------------------------------------------------

#include "common.h"

struct PS_INPUT
{
    float4 Pos : SV_POSITION;
    float4 Tex : TEXCOORD0;
};
PS_INPUT VS( float3 Pos : POSITION, float3 tex  : TEXCOORD ) 
{
	//float4x4 mat;
	//mat[0] = float4(1,0,0,0); mat[1] = float4(0,1,0,0); mat[2] = float4(0,0,1,0); mat[3] = float4(0,0,0,1);
	//mat[0] = float4(2.4142132,0,0,0); mat[1] = float4(0,2.4142132,0,0); mat[2] = float4(0,0,1.001,1); mat[3] = float4(0,0,9.0090084,10);
    //result.Pos=mul(float4(Pos, 1),mat);
	PS_INPUT result;
	result.Pos= mul(float4(Pos, 1),wvp);
	result.Tex= float4(tex, 0);
    return result;
}
