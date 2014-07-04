//=======================================================================================
// Box.h by Frank Luna (C) 2008 All Rights Reserved.
//=======================================================================================

#ifndef BOX_H
#define BOX_H

//#include "d3dUtil.h"
//////////////
// INCLUDES //
//////////////
#include <dxgi.h>
#include <d3dcommon.h>
#include <d3d11.h>
#include <d3dx10math.h>


//*****************************************************************************
// Convenience macro for releasing COM objects.
//*****************************************************************************

#define ReleaseCOM(x) { if(x){ x->Release();x = 0; } }

class Box
{
public:

	Box();
	~Box();

	void init(ID3D11Device* device, ID3D11DeviceContext*  m_pContext, float scale);
	void draw();

private:
	DWORD mNumVertices;
	DWORD mNumFaces;

	ID3D11Device* md3dDevice;
	ID3D11DeviceContext*  m_Context;
	ID3D11Buffer* mVB;
	ID3D11Buffer* mIB;
	D3D11_RASTERIZER_DESC rasterDesc;
	ID3D11RasterizerState* m_rasterState;
};

#endif // BOX_H
