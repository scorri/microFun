/*+==============================================================+\
*                                                                 *
*  Copyright (c) Microsoft Corporation.  All rights reserved.     *
*                                                                 *
*  File: HandsOnLab.cpp                                           *
*                                                                 *
*  The main source file for the DirectCompute hands-on lab        *
*                                                                 *
\+==============================================================+*/

//TO DO
//Add camera class
//Check conservation - where its going wrong - doubles versus floats -- test doubles on als machine...
//Summarize code structure to date
//periodic boundary conditions - change to reflective..
//not a 1:1 relationship between pixels & threads
//SOIL structure
//REMOVE PERIODIC BOUNDARY CONDITIONS
//LOAD in different DOC maps based on MEPSOM simulations, not have DOC-POM 
//tweakbar is not that responsive...test this
//fps

////FOR PAPER///////////////////////
//typedef between float & double
//compare with CPU simulation
//compare float vs double (difference image)
//look at different sized thread groups

#include "Time.h"
#include "Camera.h"
#include "mico_fun_gpu.h"
#include "AntTweakBar.h"	// We use this for GUI tweaking of variables during runtime
#include "Box.h"

#include <fstream>
#include <iostream>
#include <iomanip>
using namespace std;
// MSVC 
#pragma warning(disable:4244) // narrowing context
#pragma warning(disable:4305) // truncation



static float rot=3.14;
static int counter=0;
ByteArray gByteArray={0,};
scalar carbonMap[BUFFER_SIZE_X][BUFFER_SIZE_Y][BUFFER_SIZE_Z]={0.0};
Box mBox;
fstream MatLab;

struct VOLUMEVERTEX
{
	VOLUMEVERTEX (float x=0, float y=0, float z=0, float tu=0, float tv=0, float tw=0)
	{
		v.x= x; v.y= y; v.z= z; 
		t.x= tu; t.y= tv; t.z= tw; 
		//x=i; y=j; z=k;
		//tu=u; tv=v;tw=w;
	}

	D3DXVECTOR3 v, t;
//FLOAT x, y, z;
//FLOAT tu, tv, tw;
}; // VOLUMEVERTEX
  
  struct QUAD
{
	QUAD()
	{
	}
	VOLUMEVERTEX quad[6];
	
};

 #define MAX_PROXY_SLICES	2*BUFFER_SIZE_Z
//#define MAX_VERTICES	(MAX_SLICES * 6)
#define MAX_PROXY_VERTS (MAX_PROXY_SLICES * 6)

QUAD g_quaddata[MAX_PROXY_SLICES];
VOLUMEVERTEX g_vVertices[MAX_PROXY_VERTS];
VOLUMEVERTEX g_vVerticesZPos[MAX_PROXY_VERTS];
VOLUMEVERTEX g_vVerticesZNeg[MAX_PROXY_VERTS];
VOLUMEVERTEX g_vVerticesXPos[MAX_PROXY_VERTS];
VOLUMEVERTEX g_vVerticesXNeg[MAX_PROXY_VERTS];

/******************************************************************
*                                                                 *
*  main                                                           *
*                                                                 *
*  Application entry point                                        *
*                                                                 *
******************************************************************/

int __cdecl main()
{
	// create an instance of the demo app
	DemoApp app;
	//DXUTSetCallbackKeyboard( OnKeyboard );

	// run the app
	HRESULT hr = app.Run();

	// return the HRESULT
	return hr;
}


/******************************************************************
*                                                                 *
*  DemoApp::DemoApp                                               *
*                                                                 *
*  DemoApp Constructor                                            *
*                                                                 *
******************************************************************/

DemoApp::DemoApp()
{
	// perform app initialization
	Initialize();
}


/******************************************************************
*                                                                 *
*  DemoApp::~DemoApp                                              *
*                                                                 *
*  DemoApp Destructor                                             *
*                                                                 *
******************************************************************/

DemoApp::~DemoApp()
{
	// perform app cleanup
	Cleanup();
	TwTerminate();
}


/******************************************************************
*                                                                 *
*  DemoApp::Initialize                                            *
*                                                                 *
*  Initialize all member variables and set the random seed        *
*                                                                 *
******************************************************************/

void DemoApp::Initialize()
{
	m_pDevice = NULL;
	m_pContext = NULL;
	m_pComputeShaderState = NULL;

	m_cbMouseInfo = NULL;
	m_cbDisplayInfo = NULL;

	m_pBiomassBuffer0 = NULL;
	m_pBiomassBuffer1 = NULL;
	m_p3DTexBioBuff = NULL;

	m_pBiomassBuffer0SRV = NULL;
	m_pBiomassBuffer0UAV = NULL;
	m_pBiomassBuffer1SRV = NULL;
	m_pBiomassBuffer1UAV = NULL;
	m_p3DTexSRV = NULL;
	m_p3DTexUAV = NULL;

	m_pAdapter = NULL;
	m_pFactory = NULL;
	m_pSwapChain = NULL;

	m_pRTV = NULL;

	m_pVertexShader = NULL;

	m_pVertexLayout = NULL;
	m_pVertexBuffer = NULL;
	m_pVertexBufferZPos = NULL;
	m_pVertexBufferZNeg = NULL;
	m_pVertexBufferXPos = NULL;
	m_pVertexBufferXNeg = NULL;
	m_pPixelShader = NULL;


	srand( (UINT)time( NULL ) );

	D3DXMatrixIdentity(&m_projectionMatrix);

	D3DXMatrixIdentity(&mView);
	D3DXMATRIX roty;
	D3DXMatrixIdentity(&roty);

	
	GetCamera().position();

	
	GetCamera().rebuildView();

	mAppPaused=false;

		//To profile a portion of your frame, you need a trio of ID3D11Query objects. Two of them need to have the type D3D11_QUERY_TIMESTAMP, and are used to get the GPU timestamp at the start and end of the block you want to profile. The third needs to have the type D3D11_QUERY_TIMESTAMP_DISJOINT, and it tells you whether your timestamps are invalid as well as the frequency used for converting from ticks to seconds.
	ID3D11Query *pComputeEventQuery=NULL;

	ID3D11Query *pDisjointQuery=NULL;

	ID3D11Query *pBeginFrameQuery=NULL;
}

struct PixRGBAU8
{
	PixRGBAU8 (float tr=0.0, float tg=0.0, float tb=0.0, float ta=1.0)
	{
		r= 0xFF*tr; g= 0xFF*tg; b= 0xFF*tb; a= 0xFF*ta;
	}
	UINT8 r, g, b, a;
}; // PixRGBAU8

struct V3U32
{
	V3U32 (UINT tx=0, UINT ty=0, UINT tz=0) {x= tx; y= ty; z= tz;}
	UINT x, y, z;
}; // V3U32

UINT indexV3U32 (const UINT x, const UINT y, const UINT z, const V3U32& dim)
{
	return(x + dim.x * (y + dim.y * z));
} // indexV3U32
UINT loadRaw (const char * const filePath, UINT8 * const pB, const int count)
{
	ifstream inFile(filePath, ios::in); //inFile object name

	if(!inFile)
	{
		cerr << "ERROR opening vol data file" << endl;
		return(0);
	}	
	UINT8 min= 0xFF, max= 0, v;
	float t=0.0, mean;

	inFile.read((char*)pB, count);
	for (int i=0; i < count; i++)
	{	
		//inFile >> v;
		//pB[i++]= v;
		v= pB[i];
		t+=	 v;
		if (v > max) max=	v;		
		if (v < min) min=	v;
	}			
		
	inFile.close();
	mean= t / count;
	return(count);
} // loadRaw

bool setupStructureData (void)
{
	const ByteInterval	*pI=	ivlTable;
	//const RGBA32		*pC=	clrTable;
	//int					histIvlSum[IVL_ENTRIES+1]=	{0};
//	VolStat				vs;
	ByteArray			*pF= NULL;

	if (loadRaw(".//s1centre(200,200,200).u8", (UINT8*)(gByteArray.b), sizeof(gByteArray)))
	{
		//setPalette(gPalette+pI[4].min, pI[4].max - pI[4].min, pC[4]);
		//setHistogramUB(&(gByteArray.b[0][0][0]), FILE_X_DIM*FILE_Y_DIM*FILE_Z_DIM, gHist);
		return(true);
	}
	else
	cerr << "ERROR opening vol data file" << endl;
	//::MessageBox(0, ("Cannot load structure data() - FAILED"), 0, 0);
	return(false);
} // setupVolData
bool setUpCarbonMap (void)
{
		ifstream inFile(".//carbon.txt.txt", ios::in); //inFile object name
		if(!inFile)
		{
			cerr << "ERROR opening C map data file" << endl;
			//return(0);
			return false;
		}	
		double	min=	1E34f;
		double	max=	-1E34f;
		double	v, tot=0.0;
		int		nz= 0;

		for(int k=0; k<BUFFER_SIZE_Z; k++)
		{			
			for(int j=0; j<BUFFER_SIZE_Y; j++)
			{
				for(int i=0; i<BUFFER_SIZE_X; i++)
				{	
					inFile >> v;
					//cout << v << " " ;
					if (0 != v)
					{
						nz++;
						tot+=	 v;
						if (v > max) max=	(float)v;		
						if (v < min) min=	(float)v;	
						carbonMap[i][j][k]= v*2.652500e-04;//(float)v;
						if (v>0)
						{
						//gByteArray.b[i][j][k]='255';
						}

						//cout << v << " " ;
					}
				
									
				}			
			}
		}
		
		inFile.close();
		return(true);


} // setupVolData

bool cheqPixRGBAU8 (PixRGBAU8 * const pPix, const V3U32& dim)
{
	const PixRGBAU8 clr[2]=
	{
		PixRGBAU8(1.0,0.0,0.0,1.0),
		PixRGBAU8(0.0,0.0,1.0,1.0),
		//PixRGBAU8(1.0,1.0,1.0,1.0)
	};
	V3U32 i= dim, m, t;
	UINT j=0, k;

	if (pPix)
	{
		m.x= m.y= m.z= 1<<5; // mask defines chequerboard power-of-two block size
		while (i.z-- > 0)
		{
			t.z= ((i.z & m.z) > 0); // z axis chequerboard test value
			while (i.y-- > 0)
			{
				t.y= ((i.y & m.y) > 0);
				while (i.x-- > 0)
				{
					t.x=  ((i.x & m.x) > 0);
					k= t.x ^ t.y ^ t.z;
					pPix[indexV3U32(i.x, i.y, i.z, dim)]= clr[k];
				}
				i.x= dim.x;
			}
			i.y= dim.y;
		}
		return(true);
	}
	return(false);
} // cheqPixRGBAU8
/******************************************************************
*                                                                 *
*  DemoApp::Run                                                   *
*                                                                 *
*  Run the application                                            *
*                                                                 *
******************************************************************/
int GenerateproxyGeom (void)
{
	const int sliceCount=	2 * BUFFER_SIZE_Z;
	const int sliceMaxIdx=	sliceCount - 1;
	int s;
	//create proxy geometry for volume visualization Z aligned
	D3DXVECTOR3 vG(0.5 *BUFFER_SIZE_X, 0.5 * BUFFER_SIZE_Y, 0.5 * BUFFER_SIZE_Z);
	D3DXVECTOR3 vBmin(0,0,0);
#if 1
	D3DXVECTOR3 vBmax(1,1,1); // texture sampler
#else
	D3DXVECTOR3 vBmax(BUFFER_SIZE_X-1, BUFFER_SIZE_Y-1, BUFFER_SIZE_Z-1); // buffer index
#endif
	float vBz= vBmin.z, dBz= vBmax.z - vBmin.z;
	float dGz= 2 * vG.z; 

	if (sliceCount > 2)
	{
		dBz/= (sliceCount - 1);
		dGz/= sliceCount;
		vG.z= -vG.z;
	}
	else if (sliceCount <= 1)
	{
		vBz= 0.5 * (vBmin.z + vBmax.z);
		vG.z= 0;
	}
	//fill the quads z aligned
	for(s= sliceCount; s > 0; s--)
	{
		g_quaddata[s].quad[0]=VOLUMEVERTEX( vG.x, -vG.y, vG.z, vBmax.x, vBmin.y, vBz);
		g_quaddata[s].quad[1]=VOLUMEVERTEX(-vG.x, -vG.y, vG.z, vBmin.x, vBmin.y, vBz);
		g_quaddata[s].quad[2]=VOLUMEVERTEX( vG.x,  vG.y, vG.z, vBmax.x, vBmax.y, vBz);
		
		g_quaddata[s].quad[3]=VOLUMEVERTEX(-vG.x, -vG.y, vG.z, vBmin.x, vBmin.y, vBz);
		g_quaddata[s].quad[5]=VOLUMEVERTEX(-vG.x,  vG.y, vG.z, vBmin.x, vBmax.y, vBz);
		g_quaddata[s].quad[4]=VOLUMEVERTEX( vG.x,  vG.y, vG.z, vBmax.x, vBmax.y, vBz);

		vG.z+= dGz;
		vBz+= dBz;
	}
	//fill the vertices array
	int vert=0;
	for(s= sliceCount; s > 0; s--)
	{
		for(int v=0; v<6; v++)
		{
			g_vVertices[vert] = g_quaddata[s].quad[v];	
			//g_vVerticesZPos[vert] = g_quaddata[s].quad[v];
		    g_vVerticesZNeg[vert] = g_quaddata[sliceCount-s].quad[v];//5-v
			vert++;
		}
	}
//X aligned
	float vBx= vBmin.x, dBx= vBmax.x - vBmin.x;
	float dGx= 2 * vG.x; 

	if (sliceCount > 2)
	{
		dBx/= (sliceCount - 1);
		dGx/= sliceCount;
		vG.x= -vG.x;
	}
	else if (sliceCount <= 1)
	{
		vBx= 0.5 * (vBmin.x + vBmax.x);
		vG.x= 0;
	}


	for(s= sliceCount; s > 0; s--)
	{
		g_quaddata[s].quad[0]=VOLUMEVERTEX( vG.x, -vG.y, vG.z, vBx, vBmin.y, vBmax.z);
		g_quaddata[s].quad[1]=VOLUMEVERTEX( vG.x, -vG.y, -vG.z, vBx, vBmin.y, vBmin.z);
		g_quaddata[s].quad[2]=VOLUMEVERTEX( vG.x,  vG.y, vG.z, vBx, vBmax.y, vBmax.z);
		
		g_quaddata[s].quad[3]=VOLUMEVERTEX( vG.x, -vG.y,-vG.z, vBx, vBmin.y, vBmin.z);
		g_quaddata[s].quad[5]=VOLUMEVERTEX( vG.x,  vG.y,-vG.z, vBx, vBmax.y, vBmin.z);
		g_quaddata[s].quad[4]=VOLUMEVERTEX( vG.x,  vG.y, vG.z, vBx, vBmax.y, vBmax.z);

		vG.x+= dGx;
		vBx+= dBx;
	}
	//fill the vertices array
	 vert=0;
	for(s= sliceCount; s > 0; s--)
	{
		for(int v=0; v<6; v++)
		{
		
			g_vVerticesXPos[vert] = g_quaddata[s].quad[v];//fills buffer last slices first
		    g_vVerticesXNeg[vert] = g_quaddata[sliceCount-s].quad[v];//5-vfills buffer ist slices first
			vert++;
		}
	}
	//Device->CreateVertexBuffer(MAX_VERTICES*sizeof(VOLUMEVERTEX), 0, D3DFVF_VOLUMEVERTEX, D3DPOOL_MANAGED,&pVB_Zpos, NULL);
	return(vert);
} // GenerateproxyGeom

HRESULT DemoApp::Run()
{
	HRESULT hr = S_OK;

	// create the DirectCompute device
	IFR( CreateComputeDevice() );
	GenerateproxyGeom();
	// add display capabilities to the DirectCompute device
	IFR( AddDisplayCapabilities() );

	// create the Conway resources
	IFR( CreateMycoResources() );

	// compile the compute shader
	IFR( CompileComputeShader( L"fungal_activityCS.hlsl", "ActivityCS" ) );
	// compile the pixel shader
	IFR( CompilePixelShader( L"Micro_funPS.hlsl", "Micro_funPS" ) );

	// Setup anttweakbar for adjusting variables at runtime
	InitialiseTweakBar();
	
	m_DisplayInfo.ins=1;
	m_DisplayInfo.ni=1;
	m_DisplayInfo.mb=1;
	m_DisplayInfo.geometry=1;
	// run the message loop
	IFR( MessageLoop() );

	return hr;
}

void DemoApp::InitialiseTweakBar()
{
	// Initialize AntTweakBar
	if (!TwInit(TW_DIRECT3D11, m_pDevice))
	{
		MessageBoxA(m_hWnd, TwGetLastError(), "AntTweakBar initialization failed", MB_OK|MB_ICONERROR);
	}
	TwWindowSize(WINDOW_SIZE_X, WINDOW_SIZE_Y);
	TwBar *bar = TwNewBar("Myco");
	TwAddVarRO(bar, "Consv", TW_TYPE_DOUBLE, &(conservation), "");
	TwAddVarRW(bar, "Db", TW_TYPE_FLOAT, &m_MouseInfo.Db, "min=0 max=0.18f step=0.01f");
	TwAddVarRW(bar, "Dv", TW_TYPE_FLOAT, &m_MouseInfo.Dv, "min=0 max=0.18f step=0.01f");
	TwAddVarRW(bar, "Theta", TW_TYPE_FLOAT, &m_MouseInfo.theta, "min=0 max=3 step=0.5f");
	TwAddVarRW(bar, "Ni display", TW_TYPE_BOOLCPP, &m_DisplayInfo.ni, "");
	TwAddVarRW(bar, "Ins display", TW_TYPE_BOOLCPP, &m_DisplayInfo.ins, "");
	TwAddVarRW(bar, "Mb display", TW_TYPE_BOOLCPP, &m_DisplayInfo.mb, "");
	TwAddVarRW(bar, "PAUSE ComputeShader", TW_TYPE_BOOLCPP, &mAppPaused, "");
	
}


/******************************************************************
*                                                                 *
*  DemoApp::CreateComputeDevice                                   *
*                                                                 *
*  Create a basic interface to the GPU for use with DirectCompute *
*                                                                 *
******************************************************************/

HRESULT DemoApp::CreateComputeDevice()
{
	HRESULT hr = S_OK;

	// set the device to be created as single-threaded
	UINT uCreationFlags = D3D11_CREATE_DEVICE_SINGLETHREADED;

	// add debug flags for the debug configuration
//#if defined(DEBUG) || defined(_DEBUG)
//	uCreationFlags |= D3D11_CREATE_DEVICE_DEBUG;
//#endif

	// create a variable to store the feature level the device is created with
	D3D_FEATURE_LEVEL flOut = D3D_FEATURE_LEVEL_9_1;

	IFR( D3D11CreateDevice( NULL,                        // Use default graphics card
		D3D_DRIVER_TYPE_HARDWARE,    // Try to create a hardware accelerated device
		NULL,                        // Do not use external software rasterizer modules
		uCreationFlags,              // Device creation flags i.e. debug support
		NULL,                        // Try to get greatest feature level available
		0,                           // # of elements in the previous array
		D3D11_SDK_VERSION,           // SDK version
		&m_pDevice,                  // Device out
		&flOut,                      // Actual feature level created
		&m_pContext ) );             // Context out


	// check to make sure the computer supports DirectCompute
	D3D11_FEATURE_DATA_D3D10_X_HARDWARE_OPTIONS hwOptions;
	IFR( m_pDevice->CheckFeatureSupport( D3D11_FEATURE_D3D10_X_HARDWARE_OPTIONS, &hwOptions, sizeof( hwOptions ) ) );

		D3D11_FEATURE_DATA_DOUBLES dat={0};
		HRESULT r= m_pDevice->CheckFeatureSupport(D3D11_FEATURE_DOUBLES, &dat, sizeof(dat));
		if ((r))
		{
			  printf( "Device supports double precision\n" );
			//return(dat.DoublePrecisionFloatShaderOps);
		}
		else
		{
			  printf( "Device DOES NOT supports double precision\n" );
		}

	// if the computer does not support DirectCompute, display a warning and exit
	if( flOut < D3D_FEATURE_LEVEL_11_0 && !hwOptions.ComputeShaders_Plus_RawAndStructuredBuffers_Via_Shader_4_x )
	{
		OutputDebugStringA( "This lab requires DirectCompute capability, which this computer does not support.\n" );
		return E_FAIL;
	}


	printf( "Successfully created a DirectCompute device\n" );
	return hr;
}


/******************************************************************
*                                                                 *
*  DemoApp::CreateBufferResourceAndViews                          *
*                                                                 *
*  Create a resource and read-only / read-write views of it       *
*                                                                 *
******************************************************************/

HRESULT DemoApp::CreateBufferResourceAndViews
	(
	UINT uElementSize,
	UINT uNumElements,
	const void* pInitialData,
	ID3D11Buffer** ppBuffer,
	ID3D11ShaderResourceView** ppBufferSRV,
	ID3D11UnorderedAccessView** ppBufferUAV
	)
{
	HRESULT hr = S_OK;

	// create a buffer description and initialize all elements to zero
	D3D11_BUFFER_DESC bufferDesc;
	ZeroMemory( &bufferDesc, sizeof( bufferDesc ) );

	// set the size in bytes of a single element in the buffer
	bufferDesc.StructureByteStride = max(4,uElementSize);

	// set the total size of the buffer in bytes
	bufferDesc.ByteWidth = uElementSize * uNumElements;

	// enable the resource as a structured buffer
	bufferDesc.MiscFlags = D3D11_RESOURCE_MISC_BUFFER_STRUCTURED;

	if( ppBufferSRV || ppBufferUAV )
	{
		// if we are creating any views of the buffer, we will be binding the resource to a shader, so enable such a binding
		bufferDesc.BindFlags |= D3D11_BIND_SHADER_RESOURCE;
	}

	if( !ppBufferUAV )
	{
		// if we are not creating an unordered access (read-write) view of the buffer, then we can enable such optimizations by setting the resource as immutable
		bufferDesc.Usage = D3D11_USAGE_IMMUTABLE;
	}
	else
	{
		// if we are creating a unordered access (read-write) view of the buffer, then enable such a binding
		bufferDesc.BindFlags |= D3D11_BIND_UNORDERED_ACCESS;

		// in this case, we cannot make any additional optimizations, so set the resource usage to default
		bufferDesc.Usage = D3D11_USAGE_DEFAULT;
	}

	if( pInitialData )
	{
//InitData.SysMemPitch=	width*sizeof(char);
//InitData.SysMemSlicePitch=	width*height*sizeof(char);
		// if we have data to initialize the resource with, create a subresource to pass the initial data with, and initialize all elements to zero
		D3D11_SUBRESOURCE_DATA InitialSRData;
		ZeroMemory( &InitialSRData, sizeof( InitialSRData ) );

		// set the location of the initial data
		InitialSRData.pSysMem = pInitialData;
		InitialSRData.SysMemPitch= uElementSize * BUFFER_SIZE_X;
		InitialSRData.SysMemSlicePitch= uElementSize * BUFFER_SIZE_X*BUFFER_SIZE_Y;
		*ppBuffer= NULL;
		// create the buffer resource with this description using the initial data provided
		hr= m_pDevice->CreateBuffer( &bufferDesc, &InitialSRData, ppBuffer );
		IFR( hr );
	}
	else
	{
		// create the buffer resource with this description with no initial data
		IFR( m_pDevice->CreateBuffer( &bufferDesc, NULL, ppBuffer ) );
	}

	if( ppBufferSRV )
	{
		// if we wish to create a shader resource view (read-only), create a view description and initialize all elements to zero
		D3D11_SHADER_RESOURCE_VIEW_DESC viewDescSRV;
		ZeroMemory( &viewDescSRV, sizeof( viewDescSRV ) );

		// elements in the buffer may not be a standard structure, so set the format to unknown
		viewDescSRV.Format = DXGI_FORMAT_UNKNOWN;

		// this will be a buffer view of the resource
		viewDescSRV.ViewDimension = D3D11_SRV_DIMENSION_BUFFER;

		// set the view to start with the first element
		viewDescSRV.Buffer.FirstElement = 0;

		// set the view to span the entire resource
		viewDescSRV.Buffer.NumElements = uNumElements;

		// create the shader resource view of the buffer with this description
		IFR( m_pDevice->CreateShaderResourceView( *ppBuffer, &viewDescSRV, ppBufferSRV ) );
	}

	if( ppBufferUAV )
	{
		// if we wish to create a unordered access view (read-write), create a view description and initialize all elements to zero
		D3D11_UNORDERED_ACCESS_VIEW_DESC viewDescUAV;
		ZeroMemory( &viewDescUAV, sizeof( viewDescUAV ) );

		// elements in the buffer may not be a standard structure, so set the format to unknown
		viewDescUAV.Format = DXGI_FORMAT_UNKNOWN;

		// this will be a buffer view of the resource
		viewDescUAV.ViewDimension = D3D11_UAV_DIMENSION_BUFFER;

		// set the view to start with the first element
		viewDescUAV.Buffer.FirstElement = 0;

		// set the view to span the entire resource
		viewDescUAV.Buffer.NumElements = uNumElements;

		// create the unordered access view of the buffer with this description
		IFR( m_pDevice->CreateUnorderedAccessView( *ppBuffer, &viewDescUAV, ppBufferUAV ) );
	}

	return hr;
} // DemoApp::CreateBufferResourceAndViews

	//create Texture Resource and Views
HRESULT DemoApp::Create3DStructureTextureResourceAndViews
	(
	UINT width,
	UINT height,
	UINT depth,
	DXGI_FORMAT format,
	ID3D11Texture3D** ppTexture3D,
	ID3D11ShaderResourceView** ppTextureSRV
	//ID3D11UnorderedAccessView** ppTextureUAV
	)
{
	HRESULT hr = S_OK;

	// create a texture description and initialize all elements to zero
	D3D11_TEXTURE3D_DESC TexDesc={0,};
	// set the size of texture
	TexDesc.Width = width;
	TexDesc.Height = height;
	TexDesc.Depth = depth;
	TexDesc.MipLevels = 1;
	TexDesc.Format =  format;
	TexDesc.Usage = D3D11_USAGE_DEFAULT;
	TexDesc.BindFlags = D3D11_BIND_SHADER_RESOURCE;
	TexDesc.CPUAccessFlags = 0;
	TexDesc.MiscFlags = 0;

	UINT n= width*height*depth;
	
	//fill the 3D texture with some data
	D3D11_SUBRESOURCE_DATA InitData={0,};

	InitData.pSysMem= &(gByteArray);
	InitData.SysMemPitch=	width*sizeof(char);
	InitData.SysMemSlicePitch=	width*height*sizeof(char);


	//IFR(m_pDevice->CreateTexture3D(&TexDesc,&InitData, &m_p3DTexBioBuff ));		
	IFR(m_pDevice->CreateTexture3D(&TexDesc,&InitData, &m_p3DTexStructBuff ));


	
	//create SRV
	D3D11_SHADER_RESOURCE_VIEW_DESC srvDesc;
	srvDesc.Format = format;
	srvDesc.ViewDimension = D3D11_SRV_DIMENSION_TEXTURE3D;
	srvDesc.Texture3D.MostDetailedMip = 0;
	srvDesc.Texture3D.MipLevels = 1;
	IFR(m_pDevice->CreateShaderResourceView(m_p3DTexStructBuff, &srvDesc, &m_p3DTexStructSRV));

	return hr;
} // DemoApp::CreateBufferResource
//create Texture Resource and Views
HRESULT DemoApp::Create3DTextureResourceAndViews
	(
	UINT width,
	UINT height,
	UINT depth,
	DXGI_FORMAT format,
	ID3D11Texture3D** ppTexture3D,
	ID3D11ShaderResourceView** ppTextureSRV,
	ID3D11UnorderedAccessView** ppTextureUAV
	)
{
	HRESULT hr = S_OK;

	// create a texture description and initialize all elements to zero
	D3D11_TEXTURE3D_DESC TexDesc={0,};
	// set the size of texture
	TexDesc.Width = width;
	TexDesc.Height = height;
	TexDesc.Depth = depth;
	TexDesc.MipLevels = 1;
	TexDesc.Format =  DXGI_FORMAT_R8G8B8A8_UNORM;
	TexDesc.Usage = D3D11_USAGE_DEFAULT;
	TexDesc.BindFlags = D3D11_BIND_SHADER_RESOURCE | D3D11_BIND_UNORDERED_ACCESS;
	TexDesc.CPUAccessFlags = 0;
	TexDesc.MiscFlags = 0;

	UINT n= width*height*depth;
	PixRGBAU8 *pPix= new PixRGBAU8[n];
	if (NULL == pPix)
	{
		return E_OUTOFMEMORY;
	}
	cheqPixRGBAU8(pPix, V3U32(width, height, depth));

	//fill the 3D texture with some data
	D3D11_SUBRESOURCE_DATA InitData={0,};

	InitData.pSysMem= pPix;
	InitData.SysMemPitch=	width*sizeof(*pPix);
	InitData.SysMemSlicePitch=	width*height*sizeof(*pPix);


	//IFR(m_pDevice->CreateTexture3D(&TexDesc,&InitData, &m_p3DTexBioBuff ));		
	IFR(m_pDevice->CreateTexture3D(&TexDesc,NULL, &m_p3DTexBioBuff ));


	
	//create SRV
	D3D11_SHADER_RESOURCE_VIEW_DESC srvDesc;
	srvDesc.Format = DXGI_FORMAT_R8G8B8A8_UNORM; ;
	srvDesc.ViewDimension = D3D11_SRV_DIMENSION_TEXTURE3D;
	srvDesc.Texture3D.MostDetailedMip = 0;
	srvDesc.Texture3D.MipLevels = 1;
	IFR(m_pDevice->CreateShaderResourceView(m_p3DTexBioBuff, &srvDesc, ppTextureSRV));

	//crearte UAV
	D3D11_UNORDERED_ACCESS_VIEW_DESC uavDesc;
	uavDesc.Format =DXGI_FORMAT_R8G8B8A8_UNORM; ;
	uavDesc.ViewDimension = D3D11_UAV_DIMENSION_TEXTURE3D;
	uavDesc.Texture3D.MipSlice = 0;
	uavDesc.Texture3D.FirstWSlice = 0;
	uavDesc.Texture3D.WSize=BUFFER_SIZE_Z;
	IFR(m_pDevice->CreateUnorderedAccessView(m_p3DTexBioBuff, &uavDesc, ppTextureUAV));


	return hr;
} // DemoApp::CreateBufferResource

/******************************************************************
*                                                                 *
*  DemoApp::CreateStagingBuffer                          *
*                                                                 *
*  Create a resource and read-only / read-write views of it CPU_GPU reads       *
* For debugging so messages can be printed to console to ensure mass
* conservation and tag the biomass type being displayed
*                                                                 *
******************************************************************/

HRESULT DemoApp::CreateStagingBuffer
	(
	const UINT		uElementSize,
	const UINT		uNumElements,   
	ID3D11Buffer	**ppStagingBuffer
	)
{
	HRESULT hr= S_OK;

	// create a buffer description and initialize all elements to zero
	D3D11_BUFFER_DESC bufferDesc;

	ZeroMemory( &bufferDesc, sizeof( bufferDesc ) );	// clear everything

	bufferDesc.StructureByteStride = uElementSize;	// single element size in bytes
	bufferDesc.ByteWidth = uElementSize * uNumElements; // total size in bytes

	bufferDesc.BindFlags=		0;	// not bound to any shaders
	bufferDesc.Usage=			D3D11_USAGE_STAGING;
	bufferDesc.CPUAccessFlags=	D3D11_CPU_ACCESS_READ;

	// create the buffer resource with this description with no initial data
	hr= m_pDevice->CreateBuffer( &bufferDesc, NULL, ppStagingBuffer );

	IFR(hr);

	return(S_OK);
} // DemoApp::CreateStagingBuffer

int circle(BioBufType *pBiomassBuffer, float3 centre, int radius, float target_value)
{

	// initialize each cell to a random state
	int cx=centre.x, cy=centre.y, r2=radius*radius; 
	for( int x = centre.x-10; x < centre.x+10; x++ )
	{
		for( int y = centre.y-10; y < centre.y+10; y++ )
		{
			const int i= x + y * BUFFER_SIZE_X;
			const int dx= x - cx;
			const int dy= y - cy;
			const bool inCircle= ((dx * dx + dy * dy) <= r2);

			if (inCircle)
			{

				pBiomassBuffer[i].ex_re= target_value;//rand()%3==0?1:0;

			}
		}
	}
	return 1;
}

/******************************************************************
*                                                                 *
*  DemoApp::CompileComputeShader                                  *
*                                                                 *
*  Compile a compute shader                                       *
*                                                                 *
******************************************************************/

HRESULT DemoApp::CompileComputeShader(
	LPCWSTR pSrcFile,
	LPCSTR pFunctionName )
{
	HRESULT hr = S_OK;

	ID3DBlob* pBlob = NULL;         // used to store the compiled compute shader
	ID3DBlob* pErrorBlob = NULL;    // used to store any compilation errors

	// compile the compute shader from the file specified
	hr = D3DX11CompileFromFile(
		pSrcFile,                   // use the code in this file
		NULL,                       // don't use additional defines
		NULL,                       // don't use additional includes
		pFunctionName,              // compile this function
		"cs_5_0",                   // use compute shader 4.0
		NULL,                       // no compile flags
		NULL,                       // no effect flags
		NULL,                       // don't use a thread pump
		&pBlob,                     // store the compiled shader here
		&pErrorBlob,                // store any errors here
		NULL );                     // no thread pump is used, so no asynchronous HRESULT is needed

	// if there were any errors, display them
	if( pErrorBlob )
		OutputDebugStringA( (char*)pErrorBlob->GetBufferPointer() );

	if( FAILED(hr) )
		return hr;

	if(pFunctionName=="ActivityCS")
	{
		// if the compute shader was compiled successfully, create it on the GPU
		IFR( m_pDevice->CreateComputeShader(
			pBlob->GetBufferPointer(),  // use the compute shader that was compiled here
			pBlob->GetBufferSize(),     // with this size
			NULL,                       // don't use any dynamic linkage
			&m_pComputeShaderState ) );      // store the reference to the compute shader here
	}


	printf( "Successfully compiled the \"%s\" compute shader\n",pFunctionName );
	return hr;
}


/******************************************************************
*                                                                 *
*  DemoApp::AddDisplayCapabilities                                *
*                                                                 *
*  Add display capabilities to the DirectCompute device           *
*                                                                 *
******************************************************************/

HRESULT DemoApp::AddDisplayCapabilities()
{
	HRESULT hr = S_OK;

	// create the window and swap chain
	IFR( CreateDisplayWindow() );

	// create a render target, and compile and set the vertex shader
	IFR( CreateDisplayResources() );

	// show the window on the screen
	ShowWindow( m_hWnd, SW_NORMAL );

	printf( "Successfully added display capabilities to the DirectCompute device\n" );
	return hr;
}


/******************************************************************
*                                                                 *
*  DemoApp::CreateDisplayWindow                                   *
*                                                                 *
*  Create a window to draw in                                     *
*                                                                 *
******************************************************************/

HRESULT DemoApp::CreateDisplayWindow()
{
	HRESULT hr = S_OK;

	// register the window class
	WNDCLASSEX wcex;
	wcex.cbSize         = sizeof(WNDCLASSEX);               // set the size of the class
	wcex.style			= CS_HREDRAW | CS_VREDRAW;          // redraw on any movement
	wcex.lpfnWndProc	= DemoApp::WndProc;                 // set the window process function
	wcex.cbClsExtra		= 0;                                // allocate no extra space in the window class
	wcex.cbWndExtra		= sizeof(LONG_PTR);                 // allocate enough extra space in the window instance for a long pointer
	wcex.hInstance		= NULL;                             // don't link to a handle
	wcex.hIcon			= NULL;                             // don't use a window icon
	wcex.hCursor		= LoadCursor( NULL, IDC_ARROW );    // use the default system cursor
	wcex.hbrBackground	= (HBRUSH)(COLOR_WINDOW+1);         // use a basic background color
	wcex.lpszMenuName	= NULL;                             // don't use a menu
	wcex.lpszClassName	= L"GPU_MicroFun_Demo";              // set this as the window class name
	wcex.hIconSm		= NULL;                             // don't use a small window icon

	// register the window class
	IFR( RegisterClassEx(&wcex) ? S_OK : E_FAIL );


	// create a window rectangle with the dimensions specified in the header
	RECT rc;
	SetRect( &rc, 0, 0, uWindowSizeX, uWindowSizeY );

	// adjust the rectangle dimensions to compensate for the window style
	AdjustWindowRect( &rc, WS_OVERLAPPEDWINDOW & ~WS_THICKFRAME & ~WS_MAXIMIZEBOX, false );

	// create the window
	m_hWnd = CreateWindow(
		L"GPU_MicroFun_Demo",                                    // use this window class name
		L"DirectCompute MicroFun",                          // use this as the window title
		WS_OVERLAPPEDWINDOW & ~WS_THICKFRAME & ~WS_MAXIMIZEBOX, // set the window to disable resizing and maximize
		CW_USEDEFAULT,                                          // let the system position the window
		CW_USEDEFAULT,                                          // let the system position the wondow
		rc.right - rc.left,                                     // set the desired window width
		rc.bottom - rc.top,                                     // set the desired window height
		NULL,                                                   // don't set a parent window
		NULL,                                                   // don't use a menu
		NULL,                                                   // don't link to a handle
		this                                                    // pass the demo app class pointer so the window process can access it
		);

	// verify that the window was created properly
	IFR( m_hWnd ? S_OK : E_FAIL );

	// get the DXGI device, adapter, and factory automatically created with the DirectCompute device
	IDXGIDevice* pDXGIDevice = NULL;
	IFR( m_pDevice->QueryInterface( __uuidof( IDXGIDevice ), ( LPVOID* )&pDXGIDevice ) );    
	IFR( pDXGIDevice->GetAdapter( &m_pAdapter ) );
	IFR( m_pAdapter->GetParent( __uuidof( IDXGIFactory ), (LPVOID*) &m_pFactory ) );
	SAFE_RELEASE( pDXGIDevice );


	// link the DXGI factory with the window we created
	IFR( m_pFactory->MakeWindowAssociation( m_hWnd, 0 ) );


	// attach a swap chain to the factory by first creating a swap chain description and initializing its elements to zero
	DXGI_SWAP_CHAIN_DESC pSCDesc;
	ZeroMemory( &pSCDesc, sizeof( pSCDesc ) );

	// set the height to the window height specified
	pSCDesc.BufferDesc.Height = uWindowSizeY;

	// set the width to the window width specified
	pSCDesc.BufferDesc.Width = uWindowSizeX;

	// set the refresh rate to be 60Hz
	pSCDesc.BufferDesc.RefreshRate.Numerator = 60;
	pSCDesc.BufferDesc.RefreshRate.Denominator = 1;

	// set the color format to be a standard format
	pSCDesc.BufferDesc.Format = DXGI_FORMAT_R8G8B8A8_UNORM_SRGB;

	// don't use multisampling
	pSCDesc.SampleDesc.Count = 1;

	// set the swap chain to be used as a render target
	pSCDesc.BufferUsage = DXGI_USAGE_RENDER_TARGET_OUTPUT;

	// use 2 buffers to avoid flicker
	pSCDesc.BufferCount = 2;

	// link the buffer with the window we created
	pSCDesc.OutputWindow = m_hWnd;

	// the app will run windowed, so set this appropriately
	pSCDesc.Windowed = true;


	// create the swap chain with this description
	IFR( m_pFactory->CreateSwapChain( m_pDevice, &pSCDesc, &m_pSwapChain ) );




	return hr;
}


/******************************************************************
*                                                                 *
*  DemoApp::CreateDisplayResources                                *
*                                                                 *
*  Create a render target and set up the vertex shader            *
*                                                                 *
******************************************************************/

HRESULT DemoApp::CreateDisplayResources()
{
	HRESULT hr = S_OK;

	// get the swap chain back buffer resource and create a render target view of it
	ID3D11Texture2D* pBackBuffer;
	IFR( m_pSwapChain->GetBuffer( 0, __uuidof( *pBackBuffer ), ( LPVOID* )&pBackBuffer ) );
	IFR( m_pDevice->CreateRenderTargetView( pBackBuffer, NULL, &m_pRTV ) );
	SAFE_RELEASE( pBackBuffer );



	// set up a viewport
	D3D11_VIEWPORT vp;
	vp.Width   = uWindowSizeX;
	vp.Height  = uWindowSizeY;
	vp.MinDepth = 0;
	vp.MaxDepth = 1;
	vp.TopLeftX = 0;
	vp.TopLeftY = 0;
	m_pContext->RSSetViewports( 1, &vp );

	// set the render target
	m_pContext->OMSetRenderTargets( 1, &m_pRTV, NULL );

#if 1
	// create screen quad
	float svQuad[24] = {
		-0.5*BUFFER_SIZE_X, 0.5*BUFFER_SIZE_Y,1.0f, 0.0f, BUFFER_SIZE_Y,0.0f,
		0.5*BUFFER_SIZE_X, 0.5*BUFFER_SIZE_Y,1.0f, BUFFER_SIZE_X,BUFFER_SIZE_Y,0.0f,
		-0.5*BUFFER_SIZE_X, -0.5*BUFFER_SIZE_X, 1.0f, 0,0,0.0f,
		0.5*BUFFER_SIZE_X, -0.5*BUFFER_SIZE_X,1.0f, BUFFER_SIZE_X,0,0.0f,
	};
#else
		// create screen quad
	float svQuad[24] = {
		-1.0f, 1.0f, 1.0, -0.5*BUFFER_SIZE_X, +0.5*BUFFER_SIZE_Y,-0.5*BUFFER_SIZE_Z,
		1.0f, 1.0f, 1.0, 0.5*BUFFER_SIZE_X,0.5*BUFFER_SIZE_Y,-0.5*BUFFER_SIZE_Z,
		-1.0f, -1.0f, 1.0, -0.5*BUFFER_SIZE_X,-0.5*BUFFER_SIZE_Y,-0.5*BUFFER_SIZE_Z,
		1.0f, -1.0f, 1.0, 0.5*BUFFER_SIZE_X,-0.5*BUFFER_SIZE_Y,-0.5*BUFFER_SIZE_Z,
	};
#endif
	// create a vertex buffer with the screen quad coordinates
	D3D11_BUFFER_DESC vbDesc =
	{
		sizeof(g_vVertices),
		D3D11_USAGE_DEFAULT,
		D3D11_BIND_VERTEX_BUFFER,
		0,
		0
	};

	
	D3D11_SUBRESOURCE_DATA InitData;
	InitData.pSysMem = g_vVertices;
	InitData.SysMemPitch = 0;
	InitData.SysMemSlicePitch = 0;
	IFR( m_pDevice->CreateBuffer( &vbDesc, &InitData, &m_pVertexBuffer ) );

	D3D11_SUBRESOURCE_DATA InitDataZNeg;
	InitDataZNeg.pSysMem =   g_vVerticesZNeg;
	InitDataZNeg.SysMemPitch = 0;
	InitDataZNeg.SysMemSlicePitch = 0;
	IFR( m_pDevice->CreateBuffer( &vbDesc, &InitDataZNeg, &m_pVertexBufferZNeg ) );

	D3D11_SUBRESOURCE_DATA InitDataXNeg;
	InitDataXNeg.pSysMem =   g_vVerticesXNeg;
	InitDataXNeg.SysMemPitch = 0;
	InitDataXNeg.SysMemSlicePitch = 0;
	IFR( m_pDevice->CreateBuffer( &vbDesc, &InitDataXNeg, &m_pVertexBufferXNeg ) );

		D3D11_SUBRESOURCE_DATA InitDataXPos;
	InitDataXPos.pSysMem =   g_vVerticesXPos;
	InitDataXPos.SysMemPitch = 0;
	InitDataXPos.SysMemSlicePitch = 0;
	IFR( m_pDevice->CreateBuffer( &vbDesc, &InitDataXPos, &m_pVertexBufferXPos ) );
		//projection meatrix
	// Setup the projection matrix.
	//float fieldOfView = (float)3.14159265358/4.0f;
	//D3DXMatrixPerspectiveFovLH(&mProj, 0.25f*PI, aspect, 1.0f, 1000.0f);
	float screenAspect = (float)uWindowSizeX /(float)uWindowSizeX;
	// Create the projection matrix for 3D rendering.
	//D3DXMatrixPerspectiveFovLH(&m_projectionMatrix, 0.25f*3.14, screenAspect, 1.0f, 1000.0f);
	GetCamera().setLens(0.25f*3.14, screenAspect, 1.0f, 1000.0f);

	GetCamera().position() = D3DXVECTOR3(0.0f, 0.0f, -400.0f);
	//GetCamera().rotateY(3.14);
	

//	
//		// Setup the raster description which will determine how and what polygons will be drawn.
//	rasterDesc.AntialiasedLineEnable = false;
//	rasterDesc.CullMode = D3D11_CULL_NONE;
//		/*rasterDesc.DepthBias = 0;
//	rasterDesc.DepthBiasClamp = 0.0f;
//	rasterDesc.DepthClipEnable = true;*/
//	rasterDesc.FillMode = D3D11_FILL_SOLID;
//	/*rasterDesc.FrontCounterClockwise = false;
//	rasterDesc.MultisampleEnable = false;
//	rasterDesc.ScissorEnable = false;
//	rasterDesc.SlopeScaledDepthBias = 0.0f;
//*/
//	 //Create the rasterizer state from the description we just filled out.
//	m_pDevice->CreateRasterizerState(&rasterDesc, &m_rasterState);
//		// Now set the rasterizer state.
//     m_pContext->RSSetState(m_rasterState);
//
	D3D11_BLEND_DESC BlendStateDesc;
    ZeroMemory( &BlendStateDesc, sizeof(BlendStateDesc) );
    BlendStateDesc.RenderTarget[0].BlendEnable = TRUE;
    BlendStateDesc.RenderTarget[0].BlendOp = D3D11_BLEND_OP_ADD;
    BlendStateDesc.RenderTarget[0].SrcBlend = D3D11_BLEND_SRC_ALPHA;
    BlendStateDesc.RenderTarget[0].DestBlend = D3D11_BLEND_INV_SRC_ALPHA;
    BlendStateDesc.RenderTarget[0].BlendOpAlpha = D3D11_BLEND_OP_ADD;
    BlendStateDesc.RenderTarget[0].SrcBlendAlpha = D3D11_BLEND_SRC_ALPHA;
    BlendStateDesc.RenderTarget[0].DestBlendAlpha = D3D11_BLEND_INV_SRC_ALPHA;
    BlendStateDesc.RenderTarget[0].RenderTargetWriteMask = D3D11_COLOR_WRITE_ENABLE_ALL;
    m_pDevice->CreateBlendState( &BlendStateDesc, &m_BlendState  );

	//create a sampler object - put all these state creations in there own function

	D3D11_SAMPLER_DESC SamplerStateDesc;
    ZeroMemory( &SamplerStateDesc, sizeof(SamplerStateDesc) );
	SamplerStateDesc.Filter = D3D11_FILTER_MIN_MAG_MIP_POINT;
	SamplerStateDesc.AddressU = D3D11_TEXTURE_ADDRESS_CLAMP;	
	SamplerStateDesc.AddressV = D3D11_TEXTURE_ADDRESS_CLAMP;
	SamplerStateDesc.AddressW = D3D11_TEXTURE_ADDRESS_CLAMP;
	IFR( m_pDevice->CreateSamplerState( &SamplerStateDesc, &m_SamplerState  ));

	
	ID3DBlob* pBlob = NULL;
	ID3DBlob* pErrorBlob = NULL;

	// compile the vertex shader
	hr = D3DX11CompileFromFile(
		L"Micro_funVS.hlsl",       // use the code in this file
		NULL,                       // don't use additional defines
		NULL,                       // don't use additional includes
		"VS",                       // compile this function
		"vs_4_0",                   // use vertex shader 4.0
		NULL,                       // no compile flags
		NULL,                       // no effect flags
		NULL,                       // don't use a thread pump
		&pBlob,                     // compile the shader here
		&pErrorBlob,                // flush any errors here
		NULL );                     // no thread pump used, so no asynchronous HRESULT is needed

	if( pErrorBlob )
		OutputDebugStringA( (char*)pErrorBlob->GetBufferPointer() );

	if( FAILED(hr) )
		return hr;

	// create the vertex shader
	IFR( m_pDevice->CreateVertexShader( pBlob->GetBufferPointer(), pBlob->GetBufferSize(), NULL, &m_pVertexShader ) );

	// create the vertex shader input layout
	D3D11_INPUT_ELEMENT_DESC layout[] =
	{ 
		{ "POSITION", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 0, D3D11_INPUT_PER_VERTEX_DATA, 0 },
		{ "TEXCOORD", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 12, D3D11_INPUT_PER_VERTEX_DATA, 0 }
	};
	IFR( m_pDevice->CreateInputLayout( layout, 2, pBlob->GetBufferPointer(), pBlob->GetBufferSize(), &m_pVertexLayout ) );

	SAFE_RELEASE( pBlob );
	SAFE_RELEASE( pErrorBlob );

	//create bounding box for proxy geometry
	mBox.init(m_pDevice,m_pContext, 100);


	return hr;
}

/******************************************************************
*                                                                 *
*  DemoApp::CompilePixelShader                                    *
*                                                                 *
*  Compile a pixel shader                                         *
*                                                                 *
******************************************************************/

HRESULT DemoApp::CompilePixelShader(
	LPCWSTR pSrcFile,
	LPCSTR pFunctionName )
{
	HRESULT hr = S_OK;

	ID3DBlob* pBlob = NULL;         // used to store the compiled compute shader
	ID3DBlob* pErrorBlob = NULL;    // used to store any compilation errors

	hr = D3DX11CompileFromFile(
		pSrcFile,                   // use the code in this file
		NULL,                       // don't use additional defines
		NULL,                       // don't use additional includes
		pFunctionName,              // compile this function
		"ps_5_0",                   // use pixel shader 4.0
		NULL,                       // no compile flags
		NULL,                       // no effect flags
		NULL,                       // don't use a thread pump
		&pBlob,                     // store the compiled shader here
		&pErrorBlob,                // store any errors here
		NULL );                     // no thread pump is used, so no asynchronous HRESULT is needed

	if( pErrorBlob )
		OutputDebugStringA( (char*)pErrorBlob->GetBufferPointer() );

	if( FAILED(hr) )
		return hr;

	// if the pixel shader was compiled successfully, create it on the GPU
	IFR( m_pDevice->CreatePixelShader(
		pBlob->GetBufferPointer(),  // use the compute shader that was compiled here
		pBlob->GetBufferSize(),     // with this size
		NULL,                       // don't use any dynamic linkage
		&m_pPixelShader ) );        // store the reference to the pixel shader here

	printf( "Successfully compiled the \"%s\" pixel shader\n",pFunctionName );
	return hr;
}


/******************************************************************
*                                                                 *
*  DemoApp::CreateConstantBuffer                                  *
*                                                                 *
*  Create a constant buffer                                       *
*                                                                 *
******************************************************************/

HRESULT DemoApp::CreateConstantBuffer(
	UINT uSize,
	ID3D11Buffer** ppBuffer )
{
	HRESULT hr = S_OK;

	// create a buffer description and initialize its elements to zero
	D3D11_BUFFER_DESC bufferDesc;
	ZeroMemory( &bufferDesc, sizeof( bufferDesc ) );

	// flag the buffer resource to be a constant buffer
	bufferDesc.BindFlags = D3D11_BIND_CONSTANT_BUFFER;

	// constant buffer sizes must be a multiple of 16 bytes, so round up the size
	if( uSize%16 != 0 ) uSize = 16 * ( uSize / 16 ) + 16;

	// set the buffer size appropriately
	bufferDesc.ByteWidth = uSize;

	// allow CPU write access to the buffer
	bufferDesc.CPUAccessFlags = D3D11_CPU_ACCESS_WRITE;

	// the buffer must be accessible by the CPU for write, and the GPU for read, so the resource must be dynamic
	bufferDesc.Usage = D3D11_USAGE_DYNAMIC;

	// create the buffer resource with this description
	IFR( m_pDevice->CreateBuffer( &bufferDesc, NULL, ppBuffer ) );


	return hr;
}


/******************************************************************
*                                                                 *
*  DemoApp::UpdateMouseInfoCB                                     *
*                                                                 *
*  Update the mouse info constant buffer with the most recent data*
*                                                                 *
******************************************************************/

HRESULT DemoApp::UpdateMouseInfoCB()
{
	HRESULT hr = S_OK;

	// create a mapped subresource and map it to the data in the constant buffer
	D3D11_MAPPED_SUBRESOURCE MappedResource;
	m_pContext->Map( m_cbMouseInfo, 0, D3D11_MAP_WRITE_DISCARD, 0, &MappedResource );

	// copy the local mouse info data to the constant buffer
	memcpy( MappedResource.pData, &m_MouseInfo, sizeof( m_MouseInfo ) );

	// unmap the subresource
	m_pContext->Unmap( m_cbMouseInfo, 0 );


	return hr;
}
HRESULT DemoApp::UpdateCB()
{
	HRESULT hr = S_OK;

	// create a mapped subresource and map it to the data in the constant buffer
	D3D11_MAPPED_SUBRESOURCE MappedResource;
	m_pContext->Map( m_cbObjectInfo, 0, D3D11_MAP_WRITE_DISCARD, 0, &MappedResource );

	// copy the local mouse info data to the constant buffer
	memcpy( MappedResource.pData, &m_ObjectInfo, sizeof( m_ObjectInfo ) );

	// unmap the subresource
	m_pContext->Unmap( m_cbObjectInfo, 0 );

	return hr;
}
/******************************************************************
*                                                                 *
*  DemoApp::UpdateDisplayInfoCB                                     *
*                                                                 *
*  Update the display info constant buffer with the most recent data*
*                                                                 *
******************************************************************/

HRESULT DemoApp::UpdateDisplayInfoCB()
{
	HRESULT hr = S_OK;

	// create a mapped subresource and map it to the data in the constant buffer
	D3D11_MAPPED_SUBRESOURCE MappedResource;
	m_pContext->Map( m_cbDisplayInfo, 0, D3D11_MAP_WRITE_DISCARD, 0, &MappedResource );

	// copy the local mouse info data to the constant buffer
	memcpy( MappedResource.pData, &m_DisplayInfo, sizeof( m_DisplayInfo ) );

	// unmap the subresource
	m_pContext->Unmap( m_cbDisplayInfo, 0 );

	return hr;
}

void layeredCarbon(int pos, int thickness)
{
	double sumC=0.0;
	//x plane
	for( int x=pos-thickness; x<pos+thickness; x++)
	{
		for( int y = 0; y < BUFFER_SIZE_Y; y++ )
		{
			for ( int z = 0; z < BUFFER_SIZE_Z; z++)
			{
				
					if(gByteArray.b[x][y][z]==0)//not structure
					{
					carbonMap[x][y][z]= 5*2.652500e-04*1000000;
				
					sumC+=carbonMap[x][y][z];
					}
			//carbonMap[0][y][z]= 2.652500e-04*1000000;
			}	
		}
	}
	printf("Carbon %g ", sumC);
}
/******************************************************************
*                                                                 *
*  DemoApp::CreateMycoResources                                 *
*                                                                 *
*  Create resources and views        *
*                                                                 *
******************************************************************/

HRESULT DemoApp::CreateMycoResources (void)
{
		HRESULT hr = S_OK;

		//for querying

			//DO I NEED TO KEEP CREATING A QUERY>>>
	//Create an event query to sync with GPU
	qd1.Query=D3D11_QUERY_TIMESTAMP;
	qd1.MiscFlags=0;
	m_pDevice->CreateQuery(&qd1, &pComputeEventQuery);


	//Create an event query to sync with GPU
	qd2.Query=D3D11_QUERY_TIMESTAMP_DISJOINT;
	qd2.MiscFlags=0;
	m_pDevice->CreateQuery(&qd2, &pDisjointQuery);

	//Create an event query to sync with GPU
	qd3.Query=D3D11_QUERY_TIMESTAMP;
	qd3.MiscFlags=0;
	m_pDevice->CreateQuery(&qd3, &pBeginFrameQuery);

		// create an array to store the initial cell conditions of the Conway buffer
	BioBufType* pBiomassBuffer = new BioBufType[BUFFER_SIZE_X * BUFFER_SIZE_Y * BUFFER_SIZE_Z];
	IFR( pBiomassBuffer ? S_OK : E_OUTOFMEMORY );

	setupStructureData();

	//layeredCarbon(50,5);
	//layeredCarbon(100,5);
	//layeredCarbon(150,5);
	setUpCarbonMap();
	////// create an struct to store the debug data
	//DebugBufType* pDebugBuffer= new DebugBufType[BUFFER_SIZE_X*BUFFER_SIZE_Y*BUFFER_SIZE_Z];
	//pDebugBuffer.tbio = 0.0;
	//IFR( pDebugBuffer ? S_OK : E_OUTOFMEMORY );
	 double inn=0;
	 double carbon=0.0;
	// initialize each cell to a random state
	for( int x = 0; x < BUFFER_SIZE_X; x++ )
	{
		for( int y = 0; y < BUFFER_SIZE_Y; y++ )
		{
			for ( int z = 0; z < BUFFER_SIZE_Z; z++)
			{
				int i=x + (y * BUFFER_SIZE_X) + (z * BUFFER_SIZE_X*BUFFER_SIZE_Y);
				
				pBiomassBuffer[i].state =0.0;//rand()%3==0?1:0;
				pBiomassBuffer[i].mb = 0.0;//rand()%3==0?1:0;
				pBiomassBuffer[i].ins = 0.0;//rand()%3==0?1:0;
				pBiomassBuffer[i].ex_re =carbonMap[x][y][z]*1000000;// 2.652500e-04*1000000;//carbonMap[x][y][z];//rand()%3==0?1:0;
				//conservation+=pBiomassBuffer[i].ex_re;
				carbon+=pBiomassBuffer[i].ex_re;

			}
		}
	}

	printf(" TOTAL CARBON %f \n", carbon );

	// Initialise some growth by creating a plane of activity in the volume

	//for( int x = 1; x <  BUFFER_SIZE_X-1; x++)
	//{
	//	for( int y = 1; y <  BUFFER_SIZE_Y-1; y++ )
	//	{
	//		for(int z = 1; z < BUFFER_SIZE_Z-1; z++)
	//		{
	//			//const int i= x +(y * BUFFER_SIZE_X )+ (z * BUFFER_SIZE_Y*BUFFER_SIZE_Z);

					//const int i= indexV3U32(0,y,z, V3U32(BUFFER_SIZE_X,BUFFER_SIZE_Y,BUFFER_SIZE_X));
	
					//pBiomassBuffer[i].mb = 2.5507601265e-07*1000000;//rand()%3==0?1:0;theta
					//inn+=pBiomassBuffer[i].mb;
					//pBiomassBuffer[i].state = 2.5507601265e-07*1000000;//rand()%3==0?1:0;
					//inn+=pBiomassBuffer[i].state;
					//pBiomassBuffer[i].ins = 0.0;//rand()%3==0?1:0;
					//inn+=pBiomassBuffer[i].ins;
					
					//randomly allocate innoculation
					int n=200;
					while (n>0)
					{
						int xr = rand() % 199;  
						int yr = rand() % 199;  
						int zr = rand() % 199;  

						//cout << " XR XY XZ  " << xr << " " << yr << " " << zr <<  endl;

						int i= indexV3U32(xr,yr,zr, V3U32(BUFFER_SIZE_X,BUFFER_SIZE_Y,BUFFER_SIZE_X));

						if (gByteArray.b[i]>0)
						{
							
							pBiomassBuffer[i].mb = 2.5507601265e-07*1000000;//rand()%3==0?1:0;theta
							inn+=pBiomassBuffer[i].mb;
							pBiomassBuffer[i].state = 2.5507601265e-07*1000000;//rand()%3==0?1:0;
							inn+=pBiomassBuffer[i].state;
							pBiomassBuffer[i].ins = 0.0;//rand()%3==0?1:0;
							inn+=pBiomassBuffer[i].ins;
							n--;
						}
						
					}

					
		/*			}
			
		}
	}*/
	printf(" INN %f \n", inn/1000000 );

	//mechanism to add hotspots of resources
	//circle(pBiomassBuffer, float3(450,450,0), 10, 10.0);
	//circle(pBiomassBuffer, float3(350,350,0), 10, 10.0);
	//circle(pBiomassBuffer, float3(400,400,0), 10, 10.0);

	//setting initial values of fungal parameteres - Db can be increased/decreased using 5 and 6 keys
	m_MouseInfo.Db=0.18;
	m_MouseInfo.Dv=0.13;
	m_MouseInfo.theta=2.0;


	//create the debug staging buffer
	CreateStagingBuffer(
	sizeof( DebugBufType ),
	BUFFER_SIZE_X * BUFFER_SIZE_Y * BUFFER_SIZE_Z,
	&m_pDebugStagingBuffer );

	// create two buffers for storing the cell grid
	IFR( CreateBufferResourceAndViews(
		sizeof( BioBufType ),
		BUFFER_SIZE_X * BUFFER_SIZE_Y * BUFFER_SIZE_Z,
		pBiomassBuffer,
		&m_pBiomassBuffer0,
		&m_pBiomassBuffer0SRV,
		&m_pBiomassBuffer0UAV ) );



	IFR( CreateBufferResourceAndViews(
		sizeof( BioBufType ),
		BUFFER_SIZE_X * BUFFER_SIZE_Y * BUFFER_SIZE_Z,
		pBiomassBuffer,
		&m_pBiomassBuffer1,
		&m_pBiomassBuffer1SRV,
		&m_pBiomassBuffer1UAV ) );

	//	// create structure buffers for storing the soil struct grid
	//IFR( CreateBufferResourceAndViews(
	//	sizeof( UINT8),
	//	BUFFER_SIZE_X * BUFFER_SIZE_Y * BUFFER_SIZE_Z,
	//	gByteArray.b,
	//	&m_pSoilStructBuffer0,
	//	&m_pSoilStructBuffer0SRV,
	//	NULL ) );

	//create the shader and unordered access view of the texture that is to be used in the visualisation
	IFR( Create3DTextureResourceAndViews(
		BUFFER_SIZE_X,
		BUFFER_SIZE_Y, 
		BUFFER_SIZE_Z,
		DXGI_FORMAT_R8G8B8A8_UNORM,
		&m_p3DTexBioBuff,
		&m_p3DTexSRV,
		&m_p3DTexUAV ) );

		//create the shader and unordered access view of the texture that is to be used in the visualisation
	IFR( Create3DStructureTextureResourceAndViews(
		BUFFER_SIZE_X,
		BUFFER_SIZE_Y, 
		BUFFER_SIZE_Z,
		DXGI_FORMAT_R8_UINT,
		&m_p3DTexStructBuff,
		&m_p3DTexStructSRV
		) );

	// create a constant buffer for storing mouse information
	IFR( CreateConstantBuffer(
		sizeof( CB_MouseInfo ),
		&m_cbMouseInfo ) );

	IFR( CreateConstantBuffer(sizeof(CB_DisplayInfo), &m_cbDisplayInfo) );

			// create a constant buffer for storing mouse information
	IFR( CreateConstantBuffer(
		sizeof( CB_ObjectInfo ),
		&m_cbObjectInfo ) );

	//create proxy geometry
	printf( "Successfully created data resources\n" );
	return hr;
} // DemoApp::CreateConwayResources

void sumDebug (DebugBufType * const pS, const DebugBufType * const pT, int n, int counter )
{
	DebugBufType sum={0};
	
	float min_ni=1000;
	float max_ni=0;
	float min_i=1000;
	float max_i=0;
	float min_mb=1000;
	float max_mb=0;

	if (pT)
	{
		while (n-- > 0)
		{
			sum.state+=	(float)pT[n].state;
			sum.mb+=	(float)pT[n].mb;
			sum.ins+=	(float)pT[n].ins;
			sum.ex_re+=	(float)pT[n].ex_re;

			if (pT[n].state > max_ni) max_ni=pT[n].state;	
			if (pT[n].state < min_ni) min_ni=pT[n].state;

			if (pT[n].mb > max_mb) max_mb=pT[n].mb;	
			if (pT[n].mb < min_mb) min_mb=pT[n].mb;

			if (pT[n].ins > max_i) max_i=pT[n].ins;	
			if (pT[n].ins < min_i) min_i=pT[n].ins;
		}
	}
	if (pS)
	{
		*pS= sum;
	}
	float test=sum.state+sum.mb+sum.ex_re+sum.ins;
	
	//printf(" ni %g i %g mb %g re %g t %g\n", sum.state/100000,  sum.ins/100000,  sum.mb/100000, sum.ex_re/100000, (test)/100000 );
	//printf(" maxNI %g minNI %g maxMB %g minMB %g maxMIN %g minMAX %g\n", max_ni,  min_ni,  max_mb,  min_mb, max_i,  min_i);
	//print out to file 

	if(counter==1)
		MatLab.open(".//GPUevolutiondoubles.txt",ios::out);

	else
		MatLab.open(".//GPUevolutiondoubles.txt",ios::app);	
	
	MatLab << setiosflags( ios::fixed)<< setprecision(6)<< counter << " " << (sum.state+sum.ins+sum.mb)/100000 << "  " << sum.state/100000 <<" "<< sum.ins/100000 <<" "<< sum.mb/100000 <<" "<<  sum.ex_re/100000 << " " << test/100000 <<" "<<endl;	
		//printf("PRINTING output file");
		MatLab.close();

} // sumDebugf
void DemoApp::CollectTimeStamps()
{
		//printf(" IN COLLECT TIMESTAMPS ");
	// Wait for data to be available
    while (m_pContext->GetData(pDisjointQuery, NULL, 0, 0) == S_FALSE)
    {
        Sleep(1);       // Wait a bit, but give other threads a chance to run
    }
 
    //// Check whether timestamps were disjoint during the last frame
    D3D11_QUERY_DATA_TIMESTAMP_DISJOINT tsDisjoint;
    m_pContext->GetData(pDisjointQuery, &tsDisjoint, sizeof(tsDisjoint), 0);
    if (tsDisjoint.Disjoint)
    {
        return;
    }
 
    // Get all the timestamps
    UINT64 tsBeginFrame, tsComputeClear; // ... etc.
    m_pContext->GetData(pBeginFrameQuery, &tsBeginFrame, sizeof(UINT64), 0);
    m_pContext->GetData(pComputeEventQuery, &tsComputeClear, sizeof(UINT64), 0);
    // ... etc.
 
    //// Convert to real time
    float msComputeClear = float(tsComputeClear - tsBeginFrame) /
                           float(tsDisjoint.Frequency) * 1000.0f;

		printf(" TIME COMPUTE KERNEL %d \n", msComputeClear);

	//cout << " TOME COMPUTE " << msComputeClear << endl;
	//return S_OK;//hr;
}
/******************************************************************
*                                                                 *
*  DemoApp::MycoComputeFunction                                 *
*                                                                 *
*  Prepare and run the Myco compute shader                      *
*                                                                 *
******************************************************************/

HRESULT DemoApp::MycoComputeFunction()
{

	HRESULT hr;
	
	 // Begin disjoint query, and timestamp the beginning of the frame
    m_pContext->Begin(pDisjointQuery);
    m_pContext->End(pBeginFrameQuery);

	// Null views to unbind
	ID3D11UnorderedAccessView* no_UAV = 0;
	ID3D11ShaderResourceView*  no_SRV = 0;


	if(GetAsyncKeyState('8') & 0x8000) {mAppPaused=true; } ;
	if(GetAsyncKeyState('9') & 0x8000) {mAppPaused=false; } ;
	if(GetAsyncKeyState('5') & 0x8000) { if(m_MouseInfo.Db<0.18) {m_MouseInfo.Db+= 0.1f;} else{ m_MouseInfo.Db=0.10;} ;//;m_DisplayInfo.ni==true; m_DisplayInfo.mb==true; }
	printf(" Db %f \n", m_MouseInfo.Db );}
	if(GetAsyncKeyState('6') & 0x8000) { m_MouseInfo.Db=0.0;}///*if(m_MouseInfo.Db>0.0) {m_MouseInfo.Db-= 0.1f;} else{ m_MouseInfo.Db=0.0;*/} ;//;m_DisplayInfo.ni==true; m_DisplayInfo.mb==true; }
	//printf(" Db %f \n", m_MouseInfo.Db );}


	// update the mouse info constant buffer with the most recent info
	UpdateMouseInfoCB();

	// set the compute shader to use the mouse info constant buffer
	m_pContext->CSSetConstantBuffers( 0, 1, &m_cbMouseInfo );

	// set the compute shader to use the mouse info constant buffer
	// m_pContext->CSSetConstantBuffers( 1, 1, &m_cbFungalTraitsInfo );

	// set the fungal acitivity compute shader
	m_pContext->CSSetShader( m_pComputeShaderState, NULL, 0 );


	// give it read access to one Conway buffer
	m_pContext->CSSetShaderResources( 0, 1, &m_pBiomassBuffer1SRV );
			//give it read  access to 3D texture structurew
	m_pContext->CSSetShaderResources( 1, 1, &m_p3DTexStructSRV);
	// give it write access to the other Conway buffer
	UINT UAVInitialCounts = 0;
	m_pContext->CSSetUnorderedAccessViews( 0, 1, &m_pBiomassBuffer0UAV, NULL);
	//give it read write access to 3D texture
	m_pContext->CSSetUnorderedAccessViews( 1, 1, &m_p3DTexUAV, NULL);


	// run the compute shader with enough groups so that there is a thread for every pixel on the screen
	m_pContext->Dispatch( uGroupsX, uGroupsY, uGroupsZ); 
	

	// unbind the resources
	m_pContext->CSSetUnorderedAccessViews( 0, 1, &no_UAV, 0 );
	m_pContext->CSSetUnorderedAccessViews( 1, 1, &no_UAV, 0 );
	m_pContext->CSSetShaderResources( 0, 1,&no_SRV );


	
	// swap the Biomass buffers and views
	swap( m_pBiomassBuffer0, m_pBiomassBuffer1 );
	swap( m_pBiomassBuffer0SRV, m_pBiomassBuffer1SRV );
	swap( m_pBiomassBuffer0UAV, m_pBiomassBuffer1UAV );

	 m_pContext->End(pComputeEventQuery);

	 //m_pContext->End(pBeginFrameQuery);
     m_pContext->End(pDisjointQuery);



	counter++;//
	
	
	//copy back to host memory every 64th tick
	if ((counter ==1)||(counter & 0x3f)==0) //0x3f
	{
		D3D11_MAPPED_SUBRESOURCE MappedResource; 
		m_pContext->CopyResource(  m_pDebugStagingBuffer, m_pBiomassBuffer1 );
		hr= m_pContext->Map(m_pDebugStagingBuffer, 0, D3D11_MAP_READ, 0, &MappedResource);    

		if (S_OK == hr)
		{
			DebugBufType sum={0};
			DebugBufType *pRaw = reinterpret_cast< DebugBufType* >( MappedResource.pData );
			//printf(" Debug %f %f %f %f\n", pRaw[505050].state,  pRaw[505050].ins,  pRaw[505050].mb, pRaw[505050].ex_re );//flat[x+ y*WIDTH+ z*WIDTH*DEPTH]=original[x,y,z];
			//printf(" Debug %f %f %f %f\n", pRaw[1005050].state,  pRaw[1005050].ins,  pRaw[1005050].mb, pRaw[1005050].ex_re );//flat[x+ y*WIDTH+ z*WIDTH*DEPTH]=original[x,y,z];

			sumDebug(&sum, pRaw, BUFFER_SIZE_X*BUFFER_SIZE_Y*BUFFER_SIZE_Z, counter);
		
			//printf( "Groups X Y %i %i", uGroupsX , uGroupsY );
		}
		m_pContext->Unmap( m_pDebugStagingBuffer, 0 );//
		printf( "Tick copying GPU data to CPU %i", counter );
	}
	 //CollectTimeStamps();
	return S_OK;//hr;
} 

void DemoApp::FrameBufferImage(const int tick)
{
	ID3D11Texture2D* pSurface;
    HRESULT hr = m_pSwapChain->GetBuffer( 0, __uuidof( ID3D11Texture2D ), reinterpret_cast< void** >( &pSurface ) );

	WCHAR file[256];
	swprintf(file, 255, L".//image_grabs//test%i.png", tick);	

    if( ( pSurface )&&(tick % 100 == 0))
	{
		//frame 
		D3DX11SaveTextureToFile(
		 m_pContext,
        pSurface,
        D3DX11_IFF_JPG,
        (LPCWSTR)file
        );
    }
} // DemoApp::FrameBufferImage

// HACKED in - DISPLACE
int maxIdxAbsFloat (const float f[], const int n)
{
	if (n > 0)
	{
		int i= 1, im= 0;
		float af, afm= fabs(f[0]);

		while (i < n)
		{
			af= fabs(f[i]);
			if (af > afm)
			{
				afm= af;
				im= i;
			}
			++i;
		}
		return(im);
	}
	return(-1);
} // maxIdxAbsFloat

// ?
#define GEOM_AXIS_X (0)
#define GEOM_AXIS_Y (1)
#define GEOM_AXIS_Z (2) // ?
int getPrincipalAxis (const D3DXVECTOR3& v)
{
	int i= maxIdxAbsFloat(&(v.x),3);
	return(i);
} // getPrincipalAxis

/******************************************************************
*                                                                 *
*  PoozleApp::MycoDrawFunction                                    *
*                                                                 *
*  Draw the compute shader results on the screen                  *
*                                                                 *
******************************************************************/
HRESULT DemoApp::MycoDrawFunction(const int tick)
{
	HRESULT hr = S_OK;
	
	// clear the render target
	float ClearColor[4] = { 1.f, 1.f, 1.0f, 1.0f };
	m_pContext->ClearRenderTargetView( m_pRTV, ClearColor );

	//update scene

	// Update info display based on input 
	if(GetAsyncKeyState('3') & 0x8000) { m_DisplayInfo.ins=1;m_DisplayInfo.ni=0; m_DisplayInfo.mb=0; printf( "IB scalar field\n" );}
	if(GetAsyncKeyState('1') & 0x8000) { m_DisplayInfo.ins=0;m_DisplayInfo.ni=1; m_DisplayInfo.mb=0;printf( "NIB scalar field\n" );}
	if(GetAsyncKeyState('2') & 0x8000) { m_DisplayInfo.ins=0;m_DisplayInfo.ni=0; m_DisplayInfo.mb=1;printf( "MB scalar field\n" );}
	if(GetAsyncKeyState('4') & 0x8000) { m_DisplayInfo.ins=1;m_DisplayInfo.ni=1; m_DisplayInfo.mb=1;printf( "MB + NIB + IB scalar fields\n" );}
	if(GetAsyncKeyState('5') & 0x8000) { m_DisplayInfo.ex_re=1;m_DisplayInfo.ni=0; m_DisplayInfo.mb=0;m_DisplayInfo.ins=0;printf( "external resource fields\n" );}
	UpdateDisplayInfoCB();
	

	if(GetAsyncKeyState('0') & 0x8000) {rot=rot+0.01;}
	if(GetAsyncKeyState('9') & 0x8000) {rot=rot-0.01;}

	D3DXMATRIX roty;
	D3DXMatrixIdentity(&roty);
	D3DXMatrixRotationY(&roty,rot);
	
	
	
	GetCamera().rebuildView();
	D3DXMATRIX view = roty*GetCamera().view(); 
	D3DXMATRIX proj = GetCamera().proj();
	m_DisplayInfo.wvp=(view*proj);
	D3DXMATRIX tmp, tmp2;
	D3DXMatrixTranspose(&tmp, &m_DisplayInfo.wvp);
	m_DisplayInfo.wvp=tmp;
	m_DisplayInfo.geometry=1;

	//D3DXMatrixTranspose(&tmp2, &view);
	//m_DisplayInfo.wvp=tmp;
	//not quite right!!!
	D3DXVECTOR3 ov;
	ov.x= view(0,2);
	ov.y= view(1,2);
	ov.z= view(2,2);
	int axis= getPrincipalAxis(ov);
	printf( "AXIS %u", axis );
	static int lastAxis= -1;
	if (lastAxis != axis)
	{
		lastAxis= axis;
		switch (axis)
		{
			case GEOM_AXIS_X : printf("GEOM_AXIS_X"); break;
			case GEOM_AXIS_Y : printf("GEOM_AXIS_Y"); break;
			case GEOM_AXIS_Z : printf("GEOM_AXIS_Z"); break;
		}
	}

	UpdateDisplayInfoCB();
	//UpdateCB();

	// anttweakbar will crash the graphics driver if we do not explicitly set a vertex shader1
	m_pContext->VSSetShader(m_pVertexShader, NULL, 0);
	// set the pixel shader
	m_pContext->PSSetShader( m_pPixelShader, NULL, 0 );
	// set the compute shader to use the mouse info constant buffer
	m_pContext->PSSetConstantBuffers( 0, 1, &m_cbDisplayInfo );
	m_pContext->VSSetConstantBuffers( 0, 1, &m_cbDisplayInfo );

	//m_pContext->PSSetConstantBuffers( 0, 1, &m_cbObjectInfo );
	// give it read access to the Conway buffer
	m_pContext->PSSetShaderResources(0, 1, &m_pBiomassBuffer0SRV );  
	// give it read access to the 3D texture biomass buffer
	m_pContext->PSSetShaderResources( 1, 1, &m_p3DTexSRV );  
	m_pContext->PSSetShaderResources( 2, 1, &m_p3DTexStructSRV );  
	// blendFactor=D3DXVECTOR4( 0.0f, 0.0f, 0.0f, 0.0f );
	m_pContext->PSSetSamplers(0,1,&m_SamplerState);
	m_pContext->OMSetBlendState( m_BlendState, NULL, 0xFFFFFFFF  );

	//draw bounding geometry
	m_DisplayInfo.geometry=0;
	UpdateDisplayInfoCB();
	mBox.draw();

	
	m_DisplayInfo.geometry=1;
	UpdateDisplayInfoCB();

	// set the GPU to use this vertex buffer, layout, shader, and primitive topology
	UINT stride = 24;
	UINT offset = 0;
	m_pContext->IASetInputLayout( m_pVertexLayout );
	m_pContext->IASetPrimitiveTopology( D3D11_PRIMITIVE_TOPOLOGY_TRIANGLESTRIP );
	
	
	//drawgeometry=true;
	//dray proxy gemmetry
	switch(axis){
		case 0: 
			if(ov.x > 0.0f)
			{
				
			
				m_pContext->IASetVertexBuffers(0, 1, &m_pVertexBufferXNeg, &stride, &offset );
				m_pContext->Draw( MAX_PROXY_VERTS, 0 );//DrawSliceStack_NegativeZ();
							printf("NEG X");
			}
			else
			{
				m_pContext->IASetVertexBuffers( 0, 1, &m_pVertexBufferXPos, &stride, &offset );
				m_pContext->Draw( MAX_PROXY_VERTS, 0 );//DrawSliceStack_NegativeZ();
				printf("POS X");
			}
			break;
	
		case 2:
			if(ov.z > 0.0f)
			{
				m_pContext->IASetVertexBuffers( 0, 1, &m_pVertexBufferZNeg, &stride, &offset );
				m_pContext->Draw( MAX_PROXY_VERTS, 0 );//DrawSliceStack_NegativeZ();

					printf("NEG Z");
			}
			else
			{
				m_pContext->IASetVertexBuffers( 0, 1, &m_pVertexBuffer, &stride, &offset );
				m_pContext->Draw( MAX_PROXY_VERTS, 0 );//DrawSliceStack_NegativeZ();
					printf("POS Z");
				
			}
			break;
		}			


	// Draw tweakbar
	TwDraw();

	FrameBufferImage(tick);

	// unbind the pixel shader resources
	ID3D11ShaderResourceView* aSRViewsNULL[ 1 ] = { NULL };
	m_pContext->PSSetShaderResources( 0, 1, aSRViewsNULL );
	m_pContext->PSSetShaderResources( 1, 1, aSRViewsNULL );



	return hr;
}


/******************************************************************
*                                                                 *
*  DemoApp::ComputeFunction                                       *
*                                                                 *
*  Run the desired compute shader function                        *
*                                                                 *
******************************************************************/

HRESULT DemoApp::ComputeFunction()
{
	HRESULT hr = S_OK;

	// run the desired compute shader function
	IFR( MycoComputeFunction());



	return hr;
}


/******************************************************************
*                                                                 *
*  DemoApp::DrawFunction                                          *
*                                                                 *
*  Run the desired draw function                                  *
*                                                                 *
******************************************************************/

HRESULT DemoApp::DrawFunction(const int tick)
{
	HRESULT hr = S_OK;
	
	
	// run the desired draw function
	IFR( MycoDrawFunction(tick) );

	return hr;
}


/******************************************************************
*                                                                 *
*  DemoApp::WndProc                                               *
*                                                                 *
*  Handle messages to the application window                      *
*                                                                 *
******************************************************************/

LRESULT CALLBACK DemoApp::WndProc( HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam )
{
	LRESULT result = 0;

	if(uMsg==WM_CREATE)
	{
		// on creation of the window, store the pointer to the DemoApp class instance
		LPCREATESTRUCT pcs = (LPCREATESTRUCT)lParam;
		DemoApp *pDemoApp = (DemoApp *)pcs->lpCreateParams;

		::SetWindowLongPtrW(
			hWnd,
			GWLP_USERDATA,
			PtrToUlong(pDemoApp)
			);

		result = 1;
	}
	else
	{
		// if the window received a message other than WM_CREATE, first get a pointer to the parent DemoApp class instance
		DemoApp *pDemoApp = reinterpret_cast<DemoApp *>(static_cast<LONG_PTR>(
			::GetWindowLongPtrW(
			hWnd,
			GWLP_USERDATA
			)));

		// create a flag to store whether the current message was handled or not
		bool wasHandled = false;

		// if the parent class instance was successfully retrieved, try to handle the message
		if( TwEventWin(hWnd, uMsg, wParam, lParam) ) // send event message to AntTweakBar
			return 0; // event has been handled by AntTweakBar
		if(pDemoApp)
		{
			switch(uMsg)
			{
			case WM_PAINT:
				{
					// on window paint, simply validate the window
					ValidateRect(hWnd, NULL);
					result = 0;
					wasHandled = true;
					break;
				}
			case WM_MOUSEMOVE:
				{
					// on mouse movement, update the mouse info data
					pDemoApp->m_MouseInfo.mousex = ( UINT )LOWORD( lParam );
					pDemoApp->m_MouseInfo.mousey = ( UINT )HIWORD( lParam );
					break;
				}
			case WM_LBUTTONDOWN:
				{
					// when the mouse button is clicked, set its state appropriately
					pDemoApp->m_MouseInfo.mousebtn = 1;
					break;
				}
			case WM_LBUTTONUP:
				{
					// when the mouse button is released, set its state appropriately
					pDemoApp->m_MouseInfo.mousebtn = 0;
					break;
				}
			case WM_DESTROY:
				{
					// on destruction of the window, post the quit message
					PostQuitMessage(0);
					result = 1;
					wasHandled = true;
					break;
				}
			case WM_QUIT:
				{
					// when the quit message is received, return 0 to terminate the message loop
					result = 0;
					wasHandled = true;
					break;
				}
			}
		}

		if(!wasHandled)
		{
			// if the message was not handled by our own window procedure, handle it with the default one
			result = DefWindowProc(hWnd, uMsg, wParam, lParam);
		}
	}
	return result;
}
void DemoApp::updateScene(float dt)
{
	// Code computes the average frames per second, and also the 
	// average time it takes to render one frame.

	static int frameCnt = 0;
	static float t_base = 0.0f;

	frameCnt++;

	// Compute averages over one second period.
	if( (mTimer.getGameTime() - t_base) >= 1.0f )
	{
		float fps = (float)frameCnt; // fps = frameCnt / 1
		float mspf = 1000.0f / fps;

		WCHAR message[256];
		 swprintf(message,255,L"MycoCompute V1.3    FPS %.2f", (float)fps);
		 SetWindowText(m_hWnd, message);
		
		// Reset for next average.
		frameCnt = 0;
		t_base  += 1.0f;
	}
}

/******************************************************************
*                                                                 *
*  DemoApp::MessageLoop                                           *
*                                                                 *
*  Run the window message loop                                    *
*                                                                 *
******************************************************************/

HRESULT DemoApp::MessageLoop()
{
	HRESULT hr = S_OK;

	bool bGotMsg;
	MSG msg;
	msg.message = WM_NULL;
	
	mTimer.reset();
	// loop forever unless the quit message is posted
	while( msg.message != WM_QUIT )
	{
		// get the next message, if there is one
		bGotMsg = ( PeekMessage( &msg, NULL, 0U, 0U, PM_REMOVE ) != 0);

		if( bGotMsg )
		{
			// if a message was received, translate and handle it
			TranslateMessage(&msg);
			DispatchMessage(&msg);
		}
		else
		{
			if(counter < 100000)
			{
			mTimer.tick();
			updateScene(mTimer.getDeltaTime());
			// if no message was received, execute the compute function
							
			if (mAppPaused)
			{	// ALWAYS draw ONLY when paused
				DrawFunction(counter);
			}
			else
			{
			
				{
					ComputeFunction();//gShadeRC.getTick(CSRC_ID_COMPUTE)
					//gShadeRC.complete(CSRC_ID_COMPUTE);
					//printf( "Tick %u %u", counter , counter );
				}
				if ( counter%50==0 )
				{
					DrawFunction(counter);
				    //printf( "Tick draw %u %u", counter ,counter );
				}
			}
			}
			else
			{
				

				printf( "TIme elapsed %g", mTimer.getGameTime()  ); 
				msg.message = WM_QUIT;
			}
		
			// and present the image to the screen
			m_pSwapChain->Present(0,0);
		}
	}

	return hr;
}


/******************************************************************
*                                                                 *
*  DemoApp::Cleanup                                               *
*                                                                 *
*  Release any live resources and devices                         *
*                                                                 *
******************************************************************/

void DemoApp::Cleanup()
{
	SAFE_RELEASE( m_pDevice );
	SAFE_RELEASE( m_pContext );
	SAFE_RELEASE( m_pComputeShaderState );

	SAFE_RELEASE( m_cbMouseInfo );

	SAFE_RELEASE( m_pBiomassBuffer0 );
	SAFE_RELEASE( m_pBiomassBuffer1 );

	SAFE_RELEASE( m_pBiomassBuffer0SRV );
	SAFE_RELEASE( m_pBiomassBuffer0UAV );
	SAFE_RELEASE( m_pBiomassBuffer1SRV );
	SAFE_RELEASE( m_pBiomassBuffer1UAV );

	SAFE_RELEASE( m_pAdapter );
	SAFE_RELEASE( m_pFactory );
	SAFE_RELEASE( m_pSwapChain );

	SAFE_RELEASE( m_pRTV );

	SAFE_RELEASE( m_pVertexShader );

	SAFE_RELEASE( m_pVertexLayout );
	SAFE_RELEASE( m_pVertexBuffer );

	SAFE_RELEASE( m_pPixelShader );
}
