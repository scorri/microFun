/*+==============================================================+\
*                                                                 *
*  Copyright (c) Microsoft Corporation.  All rights reserved.     *
*                                                                 *
*  File: HandsOnLab.h                                             *
*                                                                 *
*  The header file for the DirectCompute hands-on lab             *
*                                                                 *
\+==============================================================+*/

#pragma once


/******************************************************************
*                                                                 *
*  Includes                                                       *
*                                                                 *
******************************************************************/

#include <windows.h>
#include <d3d11.h>
#include <d3dcompiler.h>
#include <d3dx11.h>
#include "WinUser.h"
#include "stdio.h"
#include "time.h"
#include "constants.h"
#include "GameTimer.h"

/////////////
// LINKING //
/////////////
#pragma comment(lib, "dxgi.lib")
#pragma comment(lib, "d3d11.lib")
#pragma comment(lib, "d3dx11.lib")
#pragma comment(lib, "d3dx10.lib")

//The next thing we do is include the headers for those libraries that we are linking to this object module as well as headers for DirectX type definitions and such. 

//////////////
// INCLUDES //
//////////////
#include <dxgi.h>
#include <d3dcommon.h>
#include <d3d11.h>
#include <d3dx10math.h>

/******************************************************************
*                                                                 *
*  Constant Lab Parameters                                        *
*                                                                 *
******************************************************************/

const UINT uWindowSizeX = WINDOW_SIZE_X;              // width of the window
const UINT uWindowSizeY = WINDOW_SIZE_Y;              // height of the window
const UINT uWindowArea = uWindowSizeX * uWindowSizeY;


// No. of threads in each thread group
// Ideal for 3D.. 16 x 8 x 8
const UINT uThreadsX = THREADS_X;                  // x dimension of compute shader dispatch group size
const UINT uThreadsY = THREADS_Y;                  // y dimension of compute shader dispatch group size
const UINT uThreadsZ = THREADS_Z;					// z...

// No, of thread groups to dispatch
const UINT uGroupsX = GROUPS_X;
const UINT uGroupsY = GROUPS_Y;
const UINT uGroupsZ = GROUPS_Z;

/******************************************************************
*                                                                 *
*  Macros                                                         *
*                                                                 *
******************************************************************/

#ifndef IFR
#define IFR(expr) {hr = (expr); if (FAILED(hr)){OutputDebugStringA("the following operation failed: " #expr "\n"); return(hr);}}
#endif

#ifndef SAFE_RELEASE
#define SAFE_RELEASE(p)      { if (p) { (p)->Release(); (p)=NULL; } }
#endif


/******************************************************************
*                                                                 *
*  Structure definitions                                          *
*                                                                 *
******************************************************************/
typedef struct
{
	unsigned char min, max;
} ByteInterval;

static const ByteInterval ivlTable[]=
{ 
	{0,0},		// empty
	{1, 1},		// first fungus
	{2, 2},		// second fungus
	{3, 3},		// both fungi (overlap)
	{255, 255}	// solid
};

typedef struct 
{
	float	min;
	float	max;
	float	total;
	int		countNZ;
} VolStat;
typedef struct 
{
    unsigned char	b[BUFFER_SIZE_X][BUFFER_SIZE_Y][BUFFER_SIZE_Z];
} ByteArray;
// a buffer element to store cell state information for the Conway exercise
struct BioBufType
{
	scalar state;
	scalar mb;
	scalar ins;
	scalar ex_re;

};
struct int3
{
	int x, y, z;
	int3 (int _x=0, int _y=0, int _z=0)
	{
		x=_x; y=_y; z=_z;
	}
};
struct float3
{
	float x, y, z;
	float3 (float _x=0, float _y=0, float _z=0)
	{
		x=_x; y=_y; z=_z;
	}
};
// a constant buffer to store mouse information
struct CB_MouseInfo
{
	UINT mousex;
	UINT mousey;
	UINT mousebtn;
	float Db;
	float Dv;
	float theta;
};
// a constant buffer to store what state variables to display 
struct CB_DisplayInfo
{
	D3DXMATRIXA16 wvp;
	int ins;
	int ni;
	int mb;
	int all;
	int ex_re;
	UINT geometry;
	
};
// a constant buffer to store per object info
struct CB_ObjectInfo
{
	D3DXMATRIX wvp;
};
struct DebugProcessBufType
{
	scalar uptake;
	scalar insul;
	scalar mob_i;
	scalar mob_n;
};

struct DebugBufType
{
	scalar state;
	scalar mb;
	scalar ins;
	scalar ex_re;
};
struct DebugBufTypeD
{
	double state;
	double mb;
	double ins;
	double ex_re;
};
/******************************************************************
*                                                                 *
*  DemoApp class definition                                       *
*                                                                 *
******************************************************************/

class DemoApp
{
public:
	DemoApp();                                          // constructor
	~DemoApp();                                         // destructor
	HRESULT Run();                                      // run the demo application

private:
	void Initialize();                                  // initialize the demo application
	void RenderText();

	HRESULT CreateComputeDevice();                      // create a DirectCompute device

	HRESULT DemoApp::CreateStagingBufferResourceAndUAVView ///create a staging buffer and UAV another combination of resources....
		(
	UINT uElementSize,
	UINT uNumElements,
	const void* pInitialData,
	ID3D11Buffer** ppBuffer,
	ID3D11Buffer** ppStagingBuffer,
	ID3D11ShaderResourceView** ppBufferSRV,
	ID3D11UnorderedAccessView** ppBufferUAV
	);
	
	HRESULT CreateStagingBuffer
	(
		const UINT		uElementSize,
		const UINT		uNumElements,
		ID3D11Buffer	**ppStagingBuffer
	);

	HRESULT CreateBufferResourceAndViews(               // create a buffer resource and associated views
		UINT uElementSize,
		UINT uNumElements,
		const void* pInitialData,
		ID3D11Buffer** ppBuffer,
		ID3D11ShaderResourceView** ppBufferSRV,
		ID3D11UnorderedAccessView** ppBufferUAV );

	HRESULT Create3DTextureResourceAndViews(
		UINT width,
		UINT height,
		UINT depth,
		DXGI_FORMAT format,
		ID3D11Texture3D** ppTexture3D,
		ID3D11ShaderResourceView** ppTextureSRV,
		ID3D11UnorderedAccessView** ppTextureUAV);			//create texture3D for visualisation

	HRESULT Create3DStructureTextureResourceAndViews(
		UINT width,
		UINT height,
		UINT depth,
		DXGI_FORMAT format,
		ID3D11Texture3D** ppTexture3D,
		ID3D11ShaderResourceView** ppTextureSRV
		//ID3D11UnorderedAccessView** ppTextureUAV
		);			//create texture 3D for structure

	HRESULT CompileComputeShader(                       // compile a compute shader
		LPCWSTR pSrcFile,
		LPCSTR pFunctionName );

	HRESULT CreateConstantBuffer(                       // create a constant buffer
		UINT uSize,
		ID3D11Buffer** ppBuffer );

		HRESULT DemoApp::UpdateCB();							//update cbs 
   
	 HRESULT UpdateDisplayInfoCB();                        // update the display info constant buffer

	 HRESULT UpdateMouseInfoCB();                        // update the mouse info constant buffer

	HRESULT CreateMycoResources();                    // create the myco resources
	HRESULT MycoComputeFunction();                    // run the myco compute shader
	HRESULT MycoDrawFunction(const int tick);                       // draw the output of the myco compute shader
	void updateScene(float dt);

	void Cleanup();                                     // release any live resources and devices


private:
	
	CB_DisplayInfo                m_DisplayInfo;            // a structure for storing display information
	ID3D11Buffer*               m_cbDisplayInfo;          // a constant buffer for sending the display information to the GPU
	ID3D11Device*               m_pDevice;              // the interface to the DirectCompute device
	ID3D11DeviceContext*        m_pContext;             // the interface to send commands to the device
	ID3D11ComputeShader*        m_pComputeShaderState;       // the DirectCompute function
private:
	CB_MouseInfo                m_MouseInfo;            // a structure for storing mouse information
	ID3D11Buffer*               m_cbMouseInfo;          // a constant buffer for sending the mouse information to the GPU

		
	CB_ObjectInfo                m_ObjectInfo;            // a structure for storing object  information
	ID3D11Buffer*               m_cbObjectInfo;          // a constant buffer for sending the object information to the GPU


	float						m_DebugInfo;

	float						conservation;
	


	//CB_DebugInfo

private:
	ID3D11Buffer*               m_pBiomassBuffer0;       // a general buffer for storing cell state information
	ID3D11Buffer*               m_pBiomassBuffer1;       // a general buffer for storing cell state information
	ID3D11Texture3D*			m_p3DTexBioBuff;
	ID3D11Buffer*				m_pDebugProcessBuffer;
	ID3D11ShaderResourceView*   m_pBiomassBuffer0SRV;    // a read-only view of the first cell state buffer
	ID3D11UnorderedAccessView*  m_pBiomassBuffer0UAV;    // a read-write view of the first cell state buffer
	ID3D11ShaderResourceView*   m_pBiomassBuffer1SRV;    // a read-only view of the second cell state buffer
	ID3D11UnorderedAccessView*  m_pBiomassBuffer1UAV;    // a read-write view of the second cell state buffer
	ID3D11ShaderResourceView*   m_p3DTexSRV;			 // a read-only view of the second cell state buffer
	ID3D11UnorderedAccessView*  m_p3DTexUAV;			 //read-write view of the second cell state buffer
	ID3D11Texture3D*			m_p3DTexStructBuff;
	ID3D11ShaderResourceView*   m_p3DTexStructSRV;			 // a read-only view of the 3D texture for pixels shader

	ID3D11Buffer*				m_pDebugStagingBuffer;
	ID3D11Buffer*				m_pBuffer;
	ID3D11Buffer*               m_pDebugProcessesStagingBuffer;
	ID3D11UnorderedAccessView*  m_DebugProcessesStagingUAV; // a unordered access view of the Debug buffwe


 



private:
	static LRESULT CALLBACK WndProc( HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam ); // window procedure
	HRESULT AddDisplayCapabilities();                   // add display capabilities to the DirectCompute device
	HRESULT CompilePixelShader(                         // compile a pixel shader
		LPCWSTR pSrcFile,
		LPCSTR pFunctionName );
	HRESULT MessageLoop();                              // run the message loop
	HRESULT CreateDisplayWindow();                      // attach a display window and swap chain to the device
	HRESULT CreateDisplayResources();                   // create a render target, and compile and set the vertex shader
	HRESULT ComputeFunction();                          // pre-render compute function call to redirect to exercise-specific function
	HRESULT DrawFunction(const int tick);                             // render function call to redirect to exercise-specific function
	template <class T> void swap( T* &x, T* &y )        // swap two pointers
		{ T* t = x; x = y; y = t; }
	void InitialiseTweakBar();							// Initialise and setup tweakbar for adjusting variables
	void FrameBufferImage(int count);							// Write out textures capturing framebuffer


private:
	HWND                        m_hWnd;                 // window

	IDXGIAdapter*               m_pAdapter;             // DXGI adapter
	IDXGIFactory*               m_pFactory;             // DXGI factory
	IDXGISwapChain*             m_pSwapChain;           // DXGI swap chain

	ID3D11RenderTargetView*     m_pRTV;                 // render target view

	ID3D11VertexShader*         m_pVertexShader;        // vertex shader

	ID3D11InputLayout*          m_pVertexLayout;        // vertex layout
	ID3D11Buffer*               m_pVertexBuffer;        // vertex buffer
	ID3D11Buffer*               m_pVertexBufferZPos;        // vertex buffer
	ID3D11Buffer*               m_pVertexBufferZNeg;        // vertex buffer
	ID3D11Buffer*               m_pVertexBufferXPos;        // vertex buffer
	ID3D11Buffer*               m_pVertexBufferXNeg;        // vertex buffer

	ID3D11PixelShader*          m_pPixelShader;         // pixel shader

	D3D11_RASTERIZER_DESC rasterDesc;
	//D3D11_RENDER_TARGET_BLEND_DESC blendDesc;
	ID3D11RasterizerState* m_rasterState;
	ID3D11BlendState* m_BlendState;
	ID3D11SamplerState* m_SamplerState;
		//added
	D3DXMATRIX m_projectionMatrix, mView;
	GameTimer mTimer;
	bool      mAppPaused;

};
