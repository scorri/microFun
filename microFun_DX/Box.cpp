//=======================================================================================
// Box.cpp by Frank Luna (C) 2008 All Rights Reserved.
//=======================================================================================

#include "Box.h"
#include "Vertex.h"

Box::Box()
: mNumVertices(0), mNumFaces(0), md3dDevice(0), mVB(0), mIB(0)
{
}
 
Box::~Box()
{
	ReleaseCOM(mVB);
	//ReleaseCOM(mIB);
}

void Box::init(ID3D11Device* device, ID3D11DeviceContext* deviceContext, float scale)
{
	md3dDevice = device;
	m_Context = deviceContext;
 
	mNumVertices = 6;
	//mNumFaces    = 12; // 2 per quad

	//// Create vertex buffer
    Vertex v[6];

	// Fill in the front face vertex data.
	//y aligned axes
	v[0] = Vertex(-1.0f, -1.0f, -1.0f, 0.0f, 1.0f, 1.0f);
	v[1] = Vertex(-1.0f,  1.0f, -1.0f, 0.0f, 0.0f, 1.0f);
	//x aligned axes
	v[2] = Vertex(-1.0f, -1.0f, -1.0f, 0.0f, 1.0f, 1.0f);
	v[3] = Vertex( 1.0f, -1.0f, -1.0f, 1.0f, 1.0f, 1.0f);
	//z aligned axes
	v[4] = Vertex(-1.0f, -1.0f, -1.0f, 0.0f, 1.0f, 1.0f);
	v[5] = Vertex( -1.0f, -1.0f, 1.0f, 0.0f, 1.0f, 1.0f);

 //   
	//// Fill in the front face vertex data.
	//v[0] = Vertex(-1.0f, -1.0f, -1.0f, 0.0f, 1.0f, 1.0f);
	//v[1] = Vertex(-1.0f,  1.0f, -1.0f, 0.0f, 0.0f, 1.0f);
	//v[2] = Vertex( 1.0f,  1.0f, -1.0f, 1.0f, 0.0f, 1.0f);
	//v[3] = Vertex( 1.0f, -1.0f, -1.0f, 1.0f, 1.0f, 1.0f);

	//// Fill in the back face vertex data.
	//v[4] = Vertex(-1.0f, -1.0f, 1.0f, 1.0f, 1.0f, 1.0f);
	//v[5] = Vertex( 1.0f, -1.0f, 1.0f, 0.0f, 1.0f, 1.0f);
	//v[6] = Vertex( 1.0f,  1.0f, 1.0f, 0.0f, 0.0f, 1.0f);
	//v[7] = Vertex(-1.0f,  1.0f, 1.0f, 1.0f, 0.0f, 1.0f);

	//// Fill in the top face vertex data.
	//v[8]  = Vertex(-1.0f, 1.0f, -1.0f, 0.0f, 1.0f, 1.0f);
	//v[9]  = Vertex(-1.0f, 1.0f,  1.0f, 0.0f, 0.0f, 1.0f);
	//v[10] = Vertex( 1.0f, 1.0f,  1.0f, 1.0f, 0.0f, 1.0f);
	//v[11] = Vertex( 1.0f, 1.0f, -1.0f, 1.0f, 1.0f, 1.0f);

	//// Fill in the bottom face vertex data.
	//v[12] = Vertex(-1.0f, -1.0f, -1.0f,  1.0f, 1.0f, 1.0f);
	//v[13] = Vertex( 1.0f, -1.0f, -1.0f,  0.0f, 1.0f, 1.0f);
	//v[14] = Vertex( 1.0f, -1.0f,  1.0f,  0.0f, 0.0f, 1.0f);
	//v[15] = Vertex(-1.0f, -1.0f,  1.0f,  1.0f, 0.0f, 1.0f);

	//// Fill in the left face vertex data.
	//v[16] = Vertex(-1.0f, -1.0f,  1.0f,  0.0f, 1.0f, 1.0f);
	//v[17] = Vertex(-1.0f,  1.0f,  1.0f,  0.0f, 0.0f, 1.0f);
	//v[18] = Vertex(-1.0f,  1.0f, -1.0f,  1.0f, 0.0f, 1.0f);
	//v[19] = Vertex(-1.0f, -1.0f, -1.0f,  1.0f, 1.0f, 1.0f);

	//// Fill in the right face vertex data.
	//v[20] = Vertex( 1.0f, -1.0f, -1.0f, 0.0f, 1.0f, 1.0f);
	//v[21] = Vertex( 1.0f,  1.0f, -1.0f, 0.0f, 0.0f, 1.0f);
	//v[22] = Vertex( 1.0f,  1.0f,  1.0f, 1.0f, 0.0f, 1.0f);
	//v[23] = Vertex( 1.0f, -1.0f,  1.0f, 1.0f, 1.0f, 1.0f);

	// Scale the box.
	for(DWORD i = 0; i < mNumVertices; ++i)
		v[i].pos *= scale;

	// create a vertex buffer with the screen quad coordinates
	D3D11_BUFFER_DESC vbDesc =
	{
		sizeof(v),
		D3D11_USAGE_DEFAULT,
		D3D11_BIND_VERTEX_BUFFER,
		0,
		0
	};


	D3D11_SUBRESOURCE_DATA InitData;
	InitData.pSysMem = v;
	InitData.SysMemPitch = 0;
	InitData.SysMemSlicePitch = 0;
	md3dDevice->CreateBuffer( &vbDesc, &InitData, &mVB ) ;


	// Create the index buffer

	DWORD i[6];

	// Fill in the front face index data
	i[0] = 0; i[1] = 1; 
	i[2] = 2; i[3] = 3; 
	i[4] = 4; i[5] = 5;

	//// Fill in the back face index data
	//i[6] = 4; i[7]  = 5; i[8]  = 6;
	//i[9] = 4; i[10] = 6; i[11] = 7;

	//// Fill in the top face index data
	//i[12] = 8; i[13] =  9; i[14] = 10;
	//i[15] = 8; i[16] = 10; i[17] = 11;

	//// Fill in the bottom face index data
	//i[18] = 12; i[19] = 13; i[20] = 14;
	//i[21] = 12; i[22] = 14; i[23] = 15;

	//// Fill in the left face index data
	//i[24] = 16; i[25] = 17; i[26] = 18;
	//i[27] = 16; i[28] = 18; i[29] = 19;

	//// Fill in the right face index data
	//i[30] = 20; i[31] = 21; i[32] = 22;
	//i[33] = 20; i[34] = 22; i[35] = 23;

	D3D11_BUFFER_DESC ibDesc =
	{
		sizeof(i),
		D3D11_USAGE_DEFAULT,
		D3D11_BIND_INDEX_BUFFER,
		0,
		0
	};

	D3D11_SUBRESOURCE_DATA InitDataIB;
	InitDataIB.pSysMem = i;
	InitDataIB.SysMemPitch = 0;
	InitDataIB.SysMemSlicePitch = 0;
    md3dDevice->CreateBuffer(&ibDesc, &InitDataIB, &mIB);


}

void Box::draw()
{

				// Setup the raster description which will determine how and what polygons will be drawn.
	/*rasterDesc.AntialiasedLineEnable = false;
	rasterDesc.CullMode = D3D11_CULL_NONE;
	rasterDesc.DepthBias = 0;
	rasterDesc.DepthBiasClamp = 0.0f;
	rasterDesc.DepthClipEnable = true;
	rasterDesc.FillMode = D3D11_FILL_WIREFRAME;
	rasterDesc.FrontCounterClockwise = false;
	rasterDesc.MultisampleEnable = false;
	rasterDesc.ScissorEnable = false;
	rasterDesc.SlopeScaledDepthBias = 0.0f;
*/
	 //Create the rasterizer state from the description we just filled out.
	md3dDevice->CreateRasterizerState(&rasterDesc, &m_rasterState);
//		// Now set the rasterizer state.
    m_Context->RSSetState(m_rasterState);
	//m_Context->IASetInputLayout( m_pVertexLayout );
	m_Context->IASetPrimitiveTopology( D3D11_PRIMITIVE_TOPOLOGY_LINELIST);
	UINT stride = 24;
	UINT offset = 0;
	//set to draw no fill


    m_Context->IASetVertexBuffers(0, 1, &mVB, &stride, &offset);
	m_Context->IASetIndexBuffer(mIB,DXGI_FORMAT_R32_UINT, 0);
	m_Context->DrawIndexed(6, 0, 0);
}