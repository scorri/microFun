// output kernel queue sizes
#define BUFFER_SIZE_X 128
#define BUFFER_SIZE_Y 128
#define BUFFER_SIZE_Z 64

// for determining results buffer sizes
#define THREAD_X 8
#define THREAD_Y 8
#define THREAD_Z 2

#define RD_X 16
#define RD_Y 16
#define RD_Z 4

// number of points in conway/colour buffers
#define NUM_POINTS BUFFER_SIZE_X*BUFFER_SIZE_Y*BUFFER_SIZE_Z


struct float2
{
	float x; float y;
	float2 (float _x=0, float _y=0)
	{
		x=_x; y=_y;
	}
};


struct float4
{
	float x, y, z, w;
	float4 (float _x=0, float _y=0, float _z=0, float _w=0)
	{
		x=_x; y=_y; z=_z; w=_w;
	}
};

struct int4
{
	int x, y, z, w;
	int4 (int _x=0, int _y=0, int _z=0, int _w=0)
	{
		x=_x; y=_y; z=_z; w=_w;
	}
};



#include <iostream>
#include <fstream>
#include <sstream>
#include "time.h"
#include "math.h"

#ifdef _WIN32
#include <windows.h>
#include <sys/types.h>
#include <sys/timeb.h>
#endif
#include <GL/glew.h>
#include <GL/glut.h>

#ifdef _WIN32
#include <GL/wglew.h>
#endif

#ifdef __GNUC__
#include <GL/glx.h>
#include <sys/time.h>
#endif

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#include <CL/cl_gl.h>
#endif

bool v_sync = false;
bool data_out = false;

//GLuint tex = 0;
GLuint m_vbo = 0;
GLuint m_cbo = 0;

// screen dimensions
int imWidth = 0;
int imHeight = 0;

// view params
int ox, oy;
int buttonState = 0;
float camera_trans[] = {0, 0, -3};
float camera_rot[]   = {0, 0, 0};
float camera_trans_lag[] = {0, 0, -3};
float camera_rot_lag[] = {0, 0, 0};
const float inertia = 0.1f;
float modelView[16];

///
// OpenCL variables, objects
//
cl_kernel kernel;
cl_kernel rd_kernel;
cl_kernel sum_kernel;
cl_mem rd_in;
cl_mem rd_out;
cl_mem rd_results;
cl_mem cl_colour;
cl_context context;
cl_command_queue commandQueue;
cl_program program;

cl_int compute_rd();
cl_int compute_output();
cl_int compute_sum();
void reset_rd(int, int, int);

// Return the current wall-clock time in milliseconds.
// (1000 milliseconds = 1 second.)
long long int get_time()
{
#ifdef _WIN32
	// Code for Windows
	struct _timeb timebuffer;
	_ftime64_s(&timebuffer);
	return (timebuffer.time * 1000LL) + timebuffer.millitm;
#else
	// Code for POSIX operating systems
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (tv.tv_sec * 1000LL) + (tv.tv_usec / 1000LL);
#endif
}

cl_int checkSum(cl_mem rd_gpu, const int w, const int h, const int d)
{
	cl_int errNum;
	float2* rd_cpu = new float2[w*h*d];

	// read gpu buffer into cpu memory
	errNum = clEnqueueReadBuffer(
		commandQueue, 
		rd_gpu, 
		CL_TRUE, 
		0, 
		sizeof(float2)*w*h*d, 
		rd_cpu, 
		0, 
		NULL, 
		NULL);
	if(errNum != CL_SUCCESS)
	{
		std::cerr << errNum << std::endl;
		std::cerr << "Error reading buffer" << std::endl;
	}

	// Calculate sum
	float2 sum;
	sum.x = 0.0f;
	sum.y = 0.0f;
	for(int x = 0; x < w; x++)
	{
		for(int y = 0; y < h; y++)
		{
			for(int z = 0; z < d; z ++)
			{
				sum.x += rd_cpu[x+y*w+z*w*h].x;
				sum.y += rd_cpu[x+y*w+z*w*h].y;
			}
		}
	}

	// output result
	std::cout << "Sum of U component : " << sum.x << std::endl;
	std::cout << "Sum of V component : " << sum.y << std::endl;
	
	delete [] rd_cpu;

	return errNum;
}


void reshape(int w, int h) 
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(60.0, (float) w / (float) h, 0.1, 100.0);

    glMatrixMode(GL_MODELVIEW);
    glViewport(0, 0, w, h);
}

///
// Render the vertex buffer object (VBO) contents
//
void render_vbo( int num_points ) 
{

	glPointSize(3.0);
	glEnable( GL_POINT_SMOOTH );

	/* use data from a buffer by binding it before the gl*Pointer() call */
	glBindBuffer(GL_ARRAY_BUFFER, m_cbo);
	glColorPointer(4, GL_FLOAT, 0, 0);

	glBindBuffer(GL_ARRAY_BUFFER, m_vbo);
	glVertexPointer(4, GL_FLOAT, 0, 0);
	
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	
	/* this turns on the use of glVertexPointer vertex positions for glDraw*() */
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_COLOR_ARRAY);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	
	// draw points
	glDrawArrays(GL_POINTS, 0, num_points);

	glDisable(GL_BLEND);

	/* it is good practice to leave client states disabled */
	glDisableClientState(GL_COLOR_ARRAY);
	glDisableClientState(GL_VERTEX_ARRAY);
}

void updateFPS(int t_ms)
{
	float fps = 1.0f/(t_ms*0.001f);
	char window_name[256];
	sprintf(window_name,"Grey Scott RD - FPS %.2f", fps);
	glutSetWindowTitle(window_name);
}

void renderScene(void)
{
	static int framecounter = 0;
	
	long long int start;
	if(framecounter == 0)	
		start = get_time();

	compute_rd();
	compute_output();

	// output U and V summed values
	if(data_out)
	{
		// rd_out is to be rendered so it is output
		std::cout << "Output UV data" << std::endl;
		std::cout << "CPU" << std::endl;
		checkSum(rd_out, BUFFER_SIZE_X, BUFFER_SIZE_Y, BUFFER_SIZE_Z);
		std::cout << "GPU" << std::endl;
		compute_sum();
		std::cout << std::endl;
		data_out = false;
	}
/*
	try
	{	
		// copy cl_colour to cl_volume 
		std::vector<cl::Memory> i;
		i.push_back(cl_volume);
		i.push_back(cl_colour);

		commandQueue.finish();

		// acquire the shared resources
		commandQueue.enqueueAcquireGLObjects(&i);	

		cl::size_t<3> origin;
		origin.push_back(0);
		origin.push_back(0);
		origin.push_back(0);

		cl::size_t<3> region;
		region.push_back(BUFFER_SIZE_X);
		region.push_back(BUFFER_SIZE_Y);
		region.push_back(BUFFER_SIZE_Z);

		// queue the copy instruction
		commandQueue.enqueueCopyBufferToImage(
			cl_colour,
			cl_volume,
			0,
			origin,
			region);
		
		// now release shared resources
		commandQueue.enqueueReleaseGLObjects(&i);
		commandQueue.finish();

	}
	catch (cl::Error err) 
	{
			 std::cerr
			 << "ERROR: "
			 << err.what()
			 << "("
			 << err.err()
			 << ")"
			 << std::endl;
			
		// return EXIT_FAILURE;
	}
*/


	glClearColor(0.1f, 0.1f, 0.1f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// view transform
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

    for (int c = 0; c < 3; ++c)
    {
        camera_trans_lag[c] += (camera_trans[c] - camera_trans_lag[c]) * inertia;
        camera_rot_lag[c] += (camera_rot[c] - camera_rot_lag[c]) * inertia;
    }

    glTranslatef(camera_trans_lag[0], camera_trans_lag[1], camera_trans_lag[2]);
    glRotatef(camera_rot_lag[0], 1.0, 0.0, 0.0);
    glRotatef(camera_rot_lag[1], 0.0, 1.0, 0.0);

    glGetFloatv(GL_MODELVIEW_MATRIX, modelView);

    // cube
    glColor3f(1.0, 1.0, 1.0);
    glutWireCube(2.0);
	
	// render vbo as points
	render_vbo( NUM_POINTS );	

	glFinish();
	glutSwapBuffers();
	std::swap<cl_mem>(rd_in, rd_out);

	if(framecounter == 0)
	{
		long long int end = get_time();
		//std::cout << "Time taken for frame: " << (end - start) <<" ms"<< std::endl; 
		updateFPS(end - start);
	}
	
	framecounter++;

	if(framecounter > 50)
		framecounter = 0;
}

void KeyboardGL(unsigned char key, int x, int y)
{
    switch(key) 
    {
        case '\033': // escape quits
        case '\015': // Enter quits    
        case 'Q':    // Q quits
        case 'q':    // q (or escape) quits
            // Cleanup up and quit
                exit(0);
            break;
		case 's':
			std::swap<cl_mem>(rd_in, rd_out); //step through
			break;
		case 'r':
			reset_rd(BUFFER_SIZE_X, BUFFER_SIZE_Y, BUFFER_SIZE_Z); //reset buffers
			break;
		case 'o':
			data_out = !data_out;
			break;
    }
}

void MousePressGL(int button, int state, int x, int y)
{
	int mods;

	if (state == GLUT_DOWN)
    {
        buttonState |= 1<<button;
    }
    else if (state == GLUT_UP)
    {
        buttonState = 0;
    }	

    mods = glutGetModifiers();

    if (mods & GLUT_ACTIVE_SHIFT)
    {
        buttonState = 2;
    }
    else if (mods & GLUT_ACTIVE_CTRL)
    {
        buttonState = 3;
    }

    ox = x;
    oy = y;
}

void MouseMotionGL(int x, int y)
{
    float dx, dy;
    dx = (float)(x - ox);
    dy = (float)(y - oy);

	if (buttonState == 3)
	{
	    // left+middle = zoom
	    camera_trans[2] += (dy / 100.0f) * 0.5f * fabs(camera_trans[2]);
	}
	else if (buttonState & 2)
	{
	    // middle = translate
	    camera_trans[0] += dx / 100.0f;
	    camera_trans[1] -= dy / 100.0f;
	}
	else if (buttonState & 1)
	{
	    // left = rotate
	    camera_rot[0] += dy / 5.0f;
	    camera_rot[1] += dx / 5.0f;
	}

    ox = x;
    oy = y;
}

void create_vbo(GLuint *vbo, int num_points, void* data = 0)
{
	/* request unique buffer handles */
	glGenBuffers(1, vbo);

	/* buffer vertices */
	glBindBuffer(GL_ARRAY_BUFFER, *vbo);
	//glBufferData(GL_ARRAY_BUFFER, size, vertices, GL_DYNAMIC_DRAW);
	glBufferData(GL_ARRAY_BUFFER, sizeof(float4)*num_points, data, GL_DYNAMIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

}

void delete_vbo(GLuint *vbo)
{
	glBindBuffer(1, *vbo);
	glDeleteBuffers(1, vbo);

	*vbo = 0;
}

void initGlut(int argc, char *argv[], int wWidth, int wHeight)
{
    glutInit(&argc, argv);
	
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowPosition(100,100);
    glutInitWindowSize(wWidth, wHeight);
    glutCreateWindow("GL interop");

	// set callback functions
    glutDisplayFunc(renderScene);
    glutIdleFunc(renderScene);
    glutReshapeFunc(reshape);
	glutMouseFunc(MousePressGL);
	glutMotionFunc(MouseMotionGL);
    glutKeyboardFunc(KeyboardGL);
        
    glewInit();
    if (glewIsSupported("GL_VERSION_2_1"))
        printf("Ready for OpenGL 2.1\n");
    else {
                printf("Warning: Detected that OpenGL 2.1 not supported\n");
    }

#ifdef _WIN32
	//disable vsync for faster rendering
	wglSwapIntervalEXT(v_sync); 
#elif defined( __GNUC__)	
	//glXSwapIntervalEXT(v_sync);
#elif defined(__APPLE__) 
                //todo
#endif

	glEnable(GL_DEPTH_TEST);
	glutReportErrors();

}


///
//  Round up to the nearest multiple of the group size
//
size_t RoundUp(int groupSize, int globalSize)
{
    int r = globalSize % groupSize;
    if(r == 0)
    {
        return globalSize;
    }
    else
    {
        return globalSize + groupSize - r;
    }
}

cl_int compute_sum()
{
	cl_int errNum;
	
    errNum = clSetKernelArg(sum_kernel, 0, sizeof(cl_mem), &rd_out);
    errNum |= clSetKernelArg(sum_kernel, 1, sizeof(float2) * THREAD_X * THREAD_Y * THREAD_Z, NULL);
	errNum |= clSetKernelArg(sum_kernel, 2, sizeof(cl_mem), &rd_results);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error setting kernel arguments." << std::endl;
    }

	// set work group sizes	
	size_t localWorkSize[1] = { THREAD_X*THREAD_Y*THREAD_Z } ;
    size_t globalWorkSize[1] =  { BUFFER_SIZE_X*(BUFFER_SIZE_Y/2)*BUFFER_SIZE_Z };

	cl_event prof_event;
    errNum = clEnqueueNDRangeKernel(commandQueue, sum_kernel, 1, NULL,
                                    globalWorkSize, localWorkSize,
                                    0, NULL, &prof_event);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error queuing conway kernel for execution." << std::endl;
    }
	//std::cout << "Time to find partial sums : " << completeEvent(prof_event, commandQueue) << " ms\n";


		// calculate remaining results on CPU 
	checkSum(rd_results, BUFFER_SIZE_X/THREAD_X, BUFFER_SIZE_Y/THREAD_Y/2, BUFFER_SIZE_Z/THREAD_Z);

	return 0;
}

cl_int compute_output()
{
	cl_int errNum;

    errNum = clSetKernelArg(kernel, 0, sizeof(cl_mem), &cl_colour);
	errNum |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &rd_in);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error setting kernel arguments." << std::endl;
    }

		// define work groups
    		size_t globalWorkSize[3] = { BUFFER_SIZE_X, BUFFER_SIZE_Y, BUFFER_SIZE_Z };
    		size_t localWorkSize[3] = { THREAD_X, THREAD_Y, THREAD_Z };

    		glFinish();

	errNum = clEnqueueAcquireGLObjects(commandQueue, 1, &cl_colour, 0, NULL, NULL);//&prof_ev );
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error - couldn't acquire GL tex." << std::endl;
    }
	clFinish(commandQueue);

    errNum = clEnqueueNDRangeKernel(commandQueue, kernel, 3, NULL,
                                    globalWorkSize, localWorkSize,
                                    0, NULL, NULL);//&prof_ev);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error queuing texture kernel for execution." << std::endl;
	}
	clFinish(commandQueue);

	errNum = clEnqueueReleaseGLObjects(commandQueue, 1, &cl_colour, 0, NULL, NULL);//&prof_ev );
	clFinish(commandQueue);
	return 0;

}

cl_int compute_rd()
{
	cl_int errNum;

		int4 rd_params = int4(BUFFER_SIZE_X, BUFFER_SIZE_Y, BUFFER_SIZE_Z, 0.0f); // additional parameter for input (required for padding??)

		errNum = clSetKernelArg(rd_kernel, 0, sizeof(cl_mem), &rd_in);
		errNum |= clSetKernelArg(rd_kernel, 1, sizeof(int4), &rd_params); 
		errNum |= clSetKernelArg(rd_kernel, 2, sizeof(float2) * RD_X * RD_Y * RD_Z, NULL);
		errNum |= clSetKernelArg(rd_kernel, 3, sizeof(cl_mem), &rd_out);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error setting kernel arguments." << std::endl;
    }

		int globalX = ((BUFFER_SIZE_X - 1)/(RD_X - 2) + 1) * RD_X;
		int globalY = ((BUFFER_SIZE_Y - 1)/(RD_Y - 2) + 1) * RD_Y;
		int globalZ = ((BUFFER_SIZE_Z - 1)/(RD_Z - 2) + 1) * RD_Z;


		size_t localWorkSize[3] = { RD_X, RD_Y, RD_Z } ;
		size_t globalWorkSize[3] =  { globalX, globalY, globalZ };


		// queue kernel for execution
	cl_event prof_event;
    errNum = clEnqueueNDRangeKernel(commandQueue, rd_kernel, 3, NULL,
                                    globalWorkSize, localWorkSize,
                                    0, NULL, &prof_event);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error queuing conway kernel for execution." << std::endl;
    }

    return 0;
}

void reset_rd(int rdWidth, int rdHeight, int rdDepth)
{
	cl_int errNum;
    float2 *pGSBuffer = new float2[rdWidth * rdHeight * rdDepth];
 
   // initialize each cell to a random state
    for( int x = 0; x < rdWidth; x++ )
    {
        for( int y = 0; y < rdHeight; y++ )
        {
			for( int z = 0; z < rdDepth; z++ )
			{
				pGSBuffer[ x + y * rdWidth + z * rdWidth * rdHeight ].x = 1.0;
				pGSBuffer[ x + y * rdWidth + z * rdWidth * rdHeight ].y = 0.0;
			}
        }
    }

	  // initialize each cells - pattern formation is sensitive to initial state, vlues below are known to work
    for( int x = rdWidth/3; x <  2*rdWidth/3; x++ )
    {
        for( int y = rdHeight/3; y < 2*rdHeight/3; y++ )
        {
			for( int z = rdDepth/3; z < 2*rdDepth/3; z++ )
			{
				pGSBuffer[ x + y * rdWidth + z * rdWidth * rdHeight ].x =0.5;
				pGSBuffer[ x + y * rdWidth + z * rdWidth * rdHeight ].y =0.25;
			}
	    }
    }
	
	// Now perturb the entire grid. Bound the values by [0,1]
    for (int k = 0; k < rdWidth * rdHeight * rdDepth; ++k)
    {
        if ( pGSBuffer[k].x < 1.0 )
        {
            float rRand = .02f*(float)rand() / RAND_MAX - .01f;
            pGSBuffer[k].x += rRand * pGSBuffer[k].x;
        }

        if ( pGSBuffer[k].y < 1.0 )
        {
            float rRand = .02f*(float)rand() / RAND_MAX - .01f;
            pGSBuffer[k].y += rRand * pGSBuffer[k].y;
        }
    }
/*
	// output initial sum to check
	float2 sum;
	sum.x = 0.0;
	sum.y = 0.0;
	for(int x =0; x<rdWidth; x++)
	{
		for(int y = 0; y < rdHeight; y++)
		{
			for(int z = 0; z < rdDepth; z++)
			{
				sum.x += pGSBuffer[x+y*rdWidth+z*rdWidth*rdHeight].x;
				sum.y += pGSBuffer[x+y*rdWidth+z*rdWidth*rdHeight].y;
			}
		}
	}

	// output result
	std::cout << "created initial buffer with sum.." << std::endl;
	std::cout << "Sum of U component : " << sum.x << std::endl;
	std::cout << "Sum of V component : " << sum.y << std::endl;
*/


		// create cl buffers
		rd_in = clCreateBuffer(
			context,
			CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			sizeof(float2)*rdWidth*rdHeight*rdDepth,
			pGSBuffer,&errNum);
		if( errNum != CL_SUCCESS )
	{
		std::cerr<< "Failed creating rd input buffer" << std::endl;
	}
		rd_out = clCreateBuffer(
			context,
			CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
			sizeof(float2)*rdWidth*rdHeight*rdDepth,
			pGSBuffer,
			&errNum);
	if( errNum != CL_SUCCESS )
	{
		std::cerr<< "Failed creating rd output buffer" << std::endl;
	}
		rd_results = clCreateBuffer(
			context,
			CL_MEM_READ_WRITE,
			sizeof(float2)*(rdWidth/THREAD_X)*((rdHeight/THREAD_Y)/2)*(rdDepth/THREAD_Z),
			NULL,
			&errNum);
	if( errNum != CL_SUCCESS )
	{
		std::cerr<< "Failed creating results buffer" << std::endl;
	}

    delete [] pGSBuffer;
}

void Cleanup()
{
    if (commandQueue != 0)
        clReleaseCommandQueue(commandQueue);

    if (program != 0)
        clReleaseProgram(program);

    if (context != 0)
        clReleaseContext(context);

        if( kernel != 0 ) 
                clReleaseKernel(kernel);
        
		if( rd_kernel != 0 ) 
                clReleaseKernel(rd_kernel);

		if( sum_kernel != 0 ) 
                clReleaseKernel(sum_kernel);

        if( rd_in != 0 )
                clReleaseMemObject(rd_in);

		if( rd_out != 0 )
                clReleaseMemObject(rd_out);

		if( rd_results != 0 )
                clReleaseMemObject(rd_results);

		if( cl_colour != 0 )
                clReleaseMemObject(cl_colour);
		
		// after we have released the OpenCL references, we can delete the underlying OpenGL objects
/*		if( tex != 0 ) 
		{
				glBindBuffer(GL_TEXTURE_RECTANGLE_ARB, tex );
				glDeleteBuffers(1, &tex);
		}
*/
	if( m_vbo ) 
	{
		delete_vbo( &m_vbo );
	}
	if( m_cbo ) 
	{
		delete_vbo( &m_cbo );
	}
		system("pause");
	exit(0);
}


bool isSupported(const char* extension, std::string list)
{
	size_t found = list.find(extension);
	if(found == std::string::npos)
		return false;

	return true;
}

cl_context CreateContext()
{
    cl_int errNum;
    cl_uint numPlatforms;
    cl_platform_id firstPlatformId;
    cl_context context = NULL;

    // First, select an OpenCL platform to run on.  For this example, we
    // simply choose the first available platform.  Normally, you would
    // query for all available platforms and select the most appropriate one.
    errNum = clGetPlatformIDs(1, &firstPlatformId, &numPlatforms);
    if (errNum != CL_SUCCESS || numPlatforms <= 0)
    {
        std::cerr << "Failed to find any OpenCL platforms." << std::endl;
        return NULL;
    }

	size_t extensionSize;
	errNum = clGetPlatformInfo( firstPlatformId, CL_PLATFORM_EXTENSIONS, 0, NULL, &extensionSize );
	char* list = (char*)alloca(sizeof(char)*extensionSize);
	errNum = clGetPlatformInfo( firstPlatformId, CL_PLATFORM_EXTENSIONS, extensionSize, list, NULL );
	if(!isSupported("cl_khr_gl_sharing", (std::string)list))
	{
        std::cerr << "sharing extension not available" << std::endl;
        return NULL;
    }

	if(!isSupported("cl_khr_3d_image_writes", (std::string)list))
	{
		std::cerr << "3d image writes not supported" << std::endl;
		//return 0;
	}

	// Pick first platform with gl sharing
    cl_context_properties contextProperties[] =
    {
#ifdef _WIN32
        CL_CONTEXT_PLATFORM,
        (cl_context_properties) firstPlatformId,
                CL_GL_CONTEXT_KHR,
                (cl_context_properties)wglGetCurrentContext(),
                CL_WGL_HDC_KHR,
                (cl_context_properties)wglGetCurrentDC(),
#elif defined( __GNUC__)
                CL_CONTEXT_PLATFORM, (cl_context_properties) firstPlatformId, 
                CL_GL_CONTEXT_KHR, (cl_context_properties)glXGetCurrentContext(), 
                CL_GLX_DISPLAY_KHR, (cl_context_properties)glXGetCurrentDisplay(), 
#elif defined(__APPLE__) 
                //todo
#endif
        0



    };

    context = clCreateContextFromType(contextProperties, CL_DEVICE_TYPE_GPU,
                                      NULL, NULL, &errNum);
    if (errNum != CL_SUCCESS)
    {
        std::cout << "Could not create GPU context, trying CPU..." << std::endl;
        context = clCreateContextFromType(contextProperties, CL_DEVICE_TYPE_CPU,
                                          NULL, NULL, &errNum);
        if (errNum != CL_SUCCESS)
        {
            std::cerr << "Failed to create an OpenCL GPU or CPU context." << std::endl;
            return NULL;
        }
    }

    return context;
}


std::string ReadSource(const char* fileName)
{
    std::ifstream kernelFile(fileName, std::ios::in);
    if (!kernelFile.is_open())
    {
        std::cerr << "Failed to open file for reading: " << fileName << std::endl;
        return NULL;
    }

    std::ostringstream oss;
    oss << kernelFile.rdbuf();

    std::string srcStdStr = oss.str();

	return srcStdStr;
}

///
//  Create a command queue on the first device available on the
//  context
//
cl_command_queue CreateCommandQueue(cl_context context, cl_device_id *device)
{
    cl_int errNum;
    cl_device_id *devices;
    cl_command_queue commandQueue = NULL;
    size_t deviceBufferSize = -1;

    // First get the size of the devices buffer
    errNum = clGetContextInfo(context, CL_CONTEXT_DEVICES, 0, NULL, &deviceBufferSize);
    if (errNum != CL_SUCCESS)
    {
        std::cerr << "Failed call to clGetContextInfo(...,GL_CONTEXT_DEVICES,...)";
        return NULL;
    }

    if (deviceBufferSize <= 0)
    {
        std::cerr << "No devices available.";
        return NULL;
    }

    // Allocate memory for the devices buffer
    devices = new cl_device_id[deviceBufferSize / sizeof(cl_device_id)];
    errNum = clGetContextInfo(context, CL_CONTEXT_DEVICES, deviceBufferSize, devices, NULL);
    if (errNum != CL_SUCCESS)
    {
        std::cerr << "Failed to get device IDs";
        return NULL;
    }

    // In this example, we just choose the first available device.  In a
    // real program, you would likely use all available devices or choose
    // the highest performance device based on OpenCL device queries
    commandQueue = clCreateCommandQueue(context, devices[0], CL_QUEUE_PROFILING_ENABLE, &errNum);
    if (commandQueue == NULL)
    {
        std::cerr << "Failed to create commandQueue for device 0";
        return NULL;
    }
    if (errNum != CL_SUCCESS)
    {
        std::cerr << "Failed to create command queue for profiling";
        return NULL;
    }

    *device = devices[0];
    delete [] devices;
    return commandQueue;
}

///
//  Create an OpenCL program from the kernel source file
//
cl_program CreateProgram(cl_context context, cl_device_id device, const char* fileName)
{
    cl_int errNum;
    cl_program program;

    std::ifstream kernelFile(fileName, std::ios::in);
    if (!kernelFile.is_open())
    {
        std::cerr << "Failed to open file for reading: " << fileName << std::endl;
        return NULL;
    }

    std::ostringstream oss;
    oss << kernelFile.rdbuf();

    std::string srcStdStr = oss.str();
    const char *srcStr = srcStdStr.c_str();
    program = clCreateProgramWithSource(context, 1,
                                        (const char**)&srcStr,
                                        NULL, NULL);
    if (program == NULL)
    {
        std::cerr << "Failed to create CL program from source." << std::endl;
        return NULL;
    }

    errNum = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
    if (errNum != CL_SUCCESS)
    {
        // Determine the reason for the error
        char buildLog[16384];
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG,
                              sizeof(buildLog), buildLog, NULL);

        std::cerr << "Error in kernel: " << std::endl;
        std::cerr << buildLog;
        clReleaseProgram(program);
        return NULL;
    }

    return program;
}
/*
void initTexture( int width, int height, int depth )
{
    // make a texture for output
    glGenTextures(1, &tex);              // texture 
    //glTexEnvi( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE,  GL_REPLACE );
    glBindTexture(GL_TEXTURE_3D, tex);
    glTexImage3D(GL_TEXTURE_3D, 0, GL_RGBA, width,
            height, depth, 0, GL_RGBA, GL_FLOAT, NULL );
}
*/
int main(int argc, char** argv)
{
	srand( time( NULL ) );
	printf("simulating for %d values..\n", NUM_POINTS);
	imWidth = 512; 
	imHeight = 512;
	cl_device_id device = 0;
	initGlut(argc, argv, imWidth, imHeight);
	//initTexture(BUFFER_SIZE_X, BUFFER_SIZE_Y, BUFFER_SIZE_Z);

	float4* colours = new float4[NUM_POINTS];

	float i = 1.0f/NUM_POINTS;
    for( int x = 0; x < NUM_POINTS; x++ )
    {
        colours[ x ].x = i;
		colours[ x ].y = 0;
		colours[ x ].z = 0;
		colours[ x ].w = 1;
	
		i += 1.0f/NUM_POINTS;
    }
	float4 *position_buffer = new float4[NUM_POINTS]; 

	float space_x = 2.0f/(BUFFER_SIZE_X-1);
	float space_y = 2.0f/(BUFFER_SIZE_Y-1);
	float space_z = 2.0f/(BUFFER_SIZE_Z-1);
	for(int x = 0; x < BUFFER_SIZE_X; x++)
	{
		for(int y = 0; y < BUFFER_SIZE_Y; y++)
		{
			for(int z = 0; z < BUFFER_SIZE_Z; z++)
			{
				position_buffer[x + y*BUFFER_SIZE_X + z*BUFFER_SIZE_X*BUFFER_SIZE_Y] = float4(-1 + x* space_x, -1 + y* space_y, -1 + z* space_z, 1);
			}
		}
	}

	create_vbo(&m_vbo, NUM_POINTS, position_buffer);
	create_vbo(&m_cbo, NUM_POINTS, colours);

	delete [] colours;
	delete [] position_buffer;

		std::cout << "OpenGL resources created..\nInitialising OpenCL..\n";
	// Create an OpenCL context with GL sharing 
    context = CreateContext();
    if (context == NULL)
    {
        std::cerr << "Failed to create OpenCL context." << std::endl;
        return 1;
    }

    // Create a command-queue on the first device available
    // on the created context
    commandQueue = CreateCommandQueue(context, &device);
    if (commandQueue == NULL)
    {
        Cleanup();
        return 1;
    }
    // Create OpenCL program from GLinterop.cl kernel source
    program = CreateProgram(context, device, "greyscott_rd.cl");
    if (program == NULL)
    {
        Cleanup();
        return 1;
    }
			
			// create cl buffer with rd data
			reset_rd(BUFFER_SIZE_X, BUFFER_SIZE_Y, BUFFER_SIZE_Z);
			cl_int errNum;

	cl_colour = clCreateFromGLBuffer(
		context,
		CL_MEM_WRITE_ONLY,
		m_cbo,
		&errNum);
	if( errNum != CL_SUCCESS )
	{
		std::cerr<< "Failed creating colour output buffer" << std::endl;
	}

			/*
			cl_volume = cl::Image3DGL(
				context,
				CL_MEM_WRITE_ONLY,
				GL_TEXTURE_3D,
				0,
				tex);
			*/
			// create kernel
    rd_kernel = clCreateKernel(program, "rd_kernel", NULL);
    if (rd_kernel == NULL)
    {
        std::cerr << "Failed to create rd kernel" << std::endl;
        Cleanup();
        return 1;
    }

	// Create OpenCL kernel
    kernel = clCreateKernel(program, "output_kernel", NULL);
    if (kernel == NULL)
    {
        std::cerr << "Failed to create output kernel" << std::endl;
        Cleanup();
        return 1;
    }

	// Create OpenCL kernel
    sum_kernel = clCreateKernel(program, "sum_kernel", NULL);
    if (sum_kernel == NULL)
    {
        std::cerr << "Failed to create sum kernel" << std::endl;
        Cleanup();
        return 1;
    }



	std::cout << "Initial UV data" << std::endl;
		std::cout << "CPU" << std::endl;
		checkSum(rd_out, BUFFER_SIZE_X, BUFFER_SIZE_Y, BUFFER_SIZE_Z);
	//compute_rd();
	//compute_output();

	std::cout << "Beginning Simulation..." << std::endl;
	glutMainLoop();

	Cleanup();

	return 0;
}
