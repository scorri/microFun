// output kernel queue sizes
#define BUFFER_SIZE_X 200
#define BUFFER_SIZE_Y 200
#define BUFFER_SIZE_Z 200

// for determining results buffer sizes
#define THREAD_X 8
#define THREAD_Y 4
#define THREAD_Z 4

#define RD_X 16
#define RD_Y 8
#define RD_Z 8

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

typedef double scalar;
struct BioBufType
{
	scalar state;
	scalar mb;
	scalar ins;
	scalar ex_re;
	void print()
	{
		scalar test = state + mb + ex_re + ins;
		printf("ni: %g\ti: %g\tmb: %g\tre: %g\tt %g\n", state/1000000, ins/1000000, mb/1000000, ex_re/1000000, test/1000000);
	}
};

struct SharedInfo
{
	scalar ni;
	scalar mb;
	int structure;
};

bool v_sync = false;
bool data_out = false;
long long int sim_start;

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

// OpenCL variables, objects
cl_kernel kernel;
cl_kernel diffuse_kernel;
cl_kernel update_kernel;
cl_kernel sum_kernel;
cl_kernel structure_kernel;
cl_mem bio_in;
cl_mem bio_out;
cl_mem structure_tex;
cl_mem bio_results;
cl_mem cl_colour;
cl_context context;
cl_command_queue commandQueue;
cl_program program;

//microFun resources
typedef struct 
{
    unsigned char	b[BUFFER_SIZE_X][BUFFER_SIZE_Y][BUFFER_SIZE_Z];
} ByteArray;

scalar carbonMap[BUFFER_SIZE_X][BUFFER_SIZE_Y][BUFFER_SIZE_Z]={0.0};
ByteArray gByteArray={0,};
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

// forward declarations
cl_int compute_diffuse();
cl_int compute_update();
cl_int compute_microOut();
cl_int compute_microSum();
cl_int compute_structure();

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

// calculates the sum of the contents of rd_gpu buffer and prints
cl_int checkMicroSum(cl_mem rd_gpu, const int w, const int h, const int d)
{
	cl_int errNum;
	BioBufType* rd_cpu = new BioBufType[w*h*d];

	// read gpu buffer into cpu memory
	errNum = clEnqueueReadBuffer(
		commandQueue, 
		rd_gpu, 
		CL_TRUE, 
		0, 
		sizeof(BioBufType)*w*h*d, 
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
	BioBufType sum;
	sum.state = 0.0;
	sum.mb = 0.0;
	sum.ins = 0.0;
	sum.ex_re = 0.0;
	for(int z = 0; z < d; z++)
	{
		for(int y = 0; y < h; y++)
		{
			for(int x = 0; x < w; x ++)
			{
				sum.state += rd_cpu[x+y*w+z*w*h].state;
				sum.mb += rd_cpu[x+y*w+z*w*h].mb;
				sum.ins += rd_cpu[x+y*w+z*w*h].ins;
				sum.ex_re += rd_cpu[x+y*w+z*w*h].ex_re;
			}
		}
	}

	// output result
	sum.print();

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

// Render the vertex buffer object (VBO) contents
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

	// enable blend function so that alpha channel can be used to output
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

void draw()
{
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
}

void renderScene(void)
{
	static int framecounter = 0;
	
	long long int start;
	if(framecounter == 0)	
		start = get_time();
	
	// update colour output buffer first 
	compute_microOut();
			
	compute_diffuse();
	compute_update();	
	//compute_structure();	// disable all other compute_() to output structure buffer

	draw();

	if(framecounter == 0)
	{
		long long int end = get_time();
		updateFPS(end - start);
	}
	
	framecounter++;

	if(framecounter > 50)
	{
		framecounter = 0;
	}

	if(data_out)
	{
		long long int sim_end = get_time();

		std::cout << "Output UV data at " << (sim_end - sim_start)/1000 << " s into simulation" << std::endl;
		std::cout << "CPU" << std::endl;
		checkMicroSum(bio_in, BUFFER_SIZE_X, BUFFER_SIZE_Y, BUFFER_SIZE_Z);
		std::cout << "GPU" << std::endl;
		compute_microSum();
		std::cout << std::endl;	
		data_out = !data_out;
	}
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

//  Round up to the nearest multiple of the group size
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

// calculates the partial sum of components for the input buffer (bio_in)
cl_int compute_microSum()
{
	cl_int errNum;
	
    errNum = clSetKernelArg(sum_kernel, 0, sizeof(cl_mem), &bio_in);
    errNum |= clSetKernelArg(sum_kernel, 1, sizeof(BioBufType) * THREAD_X * THREAD_Y * THREAD_Z, NULL);
	errNum |= clSetKernelArg(sum_kernel, 2, sizeof(cl_mem), &bio_results);
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

	// finish computation on CPU 
	checkMicroSum(bio_results, BUFFER_SIZE_X/THREAD_X, BUFFER_SIZE_Y/THREAD_Y/2, BUFFER_SIZE_Z/THREAD_Z);

	return 0;
}

// updates the colour buffer using bio_in
cl_int compute_microOut()
{
	cl_int errNum;

    errNum = clSetKernelArg(kernel, 0, sizeof(cl_mem), &cl_colour);
	errNum |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &bio_in);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error setting kernel arguments." << std::endl;
    }

	// define work groups
    size_t globalWorkSize[1] = { NUM_POINTS };
    size_t localWorkSize[1] = { 32 };

    glFinish();

	errNum = clEnqueueAcquireGLObjects(commandQueue, 1, &cl_colour, 0, NULL, NULL);//&prof_ev );
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error - couldn't acquire GL tex." << std::endl;
    }
	clFinish(commandQueue);

    errNum = clEnqueueNDRangeKernel(commandQueue, kernel, 1, NULL,
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

// compute update using uptake, insulation, etc 
// second part of original kernel
cl_int compute_update()
{
	cl_int errNum;

	errNum = clSetKernelArg(update_kernel, 0, sizeof(cl_mem), &bio_out);
	errNum |= clSetKernelArg(update_kernel, 1, sizeof(cl_mem), &bio_in);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error setting kernel arguments." << std::endl;
    }

	size_t localWorkSize[1] = { 32 } ;
	size_t globalWorkSize[1] =  { NUM_POINTS };

	// queue kernel for execution
	cl_event prof_event;
    errNum = clEnqueueNDRangeKernel(commandQueue, update_kernel, 1, NULL,
                                    globalWorkSize, localWorkSize,
                                    0, NULL, &prof_event);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error queuing conway kernel for execution." << std::endl;
    }

    return 0;
}

// compute diffusion
// first part of original kernel
cl_int compute_diffuse()
{
	cl_int errNum;

	int4 rd_params = int4(BUFFER_SIZE_X, BUFFER_SIZE_Y, BUFFER_SIZE_Z, 0.0f); // additional parameter for input (required for padding??)

	errNum = clSetKernelArg(diffuse_kernel, 0, sizeof(cl_mem), &bio_in);
	errNum |= clSetKernelArg(diffuse_kernel, 1, sizeof(int4), &rd_params); 
	errNum |= clSetKernelArg(diffuse_kernel, 2, sizeof(SharedInfo) * RD_X * RD_Y * RD_Z, NULL);
	errNum |= clSetKernelArg(diffuse_kernel, 3, sizeof(cl_mem), &bio_out);
	errNum |= clSetKernelArg(diffuse_kernel, 4, sizeof(cl_mem), &structure_tex);
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
    errNum = clEnqueueNDRangeKernel(commandQueue, diffuse_kernel, 3, NULL,
                                    globalWorkSize, localWorkSize,
                                    0, NULL, &prof_event);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error queuing conway kernel for execution." << std::endl;
    }

    return 0;
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

		if( sum_kernel != 0 ) 
                clReleaseKernel(sum_kernel);

        if( diffuse_kernel != 0 ) 
                clReleaseKernel(diffuse_kernel);

		if( update_kernel != 0 ) 
                clReleaseKernel(update_kernel);

		if( structure_kernel != 0 ) 
                clReleaseKernel(structure_kernel);

        if( bio_in != 0 )
                clReleaseMemObject(bio_in);

		if( bio_out != 0 )
                clReleaseMemObject(bio_out);
		
		if( structure_tex != 0 )
                clReleaseMemObject(structure_tex);

		if( bio_results != 0 )
                clReleaseMemObject(bio_results);

		if( cl_colour != 0 )
                clReleaseMemObject(cl_colour);

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

// check extension is supported
bool isSupported(const char* extension, std::string list)
{
	size_t found = list.find(extension);
	if(found == std::string::npos)
		return false;

	return true;
}

// create an opencl context with gl sharing
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

// read .cl file
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

//  Create a command queue on the first device available on the
//  context with profiling enabled
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

//  Create an OpenCL program from the kernel source file
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

bool setUpCarbonMap (void)
{
	std::ifstream inFile("carbon.txt.txt", std::ios::in); //inFile object name
	if(!inFile)
	{
		std::cerr << "ERROR opening C map data file" << std::endl;
		return false;
	}	
	scalar	min=	1E34f;
	scalar	max=	-1E34f;
	scalar	v, tot=0.0;
	int		nz= 0;

	for(int k=0; k<BUFFER_SIZE_X; k++)
	{			
		for(int j=0; j<BUFFER_SIZE_Y; j++)
		{
			for(int i=0; i<BUFFER_SIZE_Z; i++)
			{	
				inFile >> v;
				{
					nz++;
					tot+=	 v;
					if (v > max) max=	(float)v;		
					if (v < min) min=	(float)v;	
					carbonMap[k][j][i]=265.2500e-06;
					// FIXME: order of array indices above?
				}								
			}			
		}
	}
	
	inFile.close();
	return true;
}

UINT loadRaw (const char * const filePath, UINT8 * const pB, const int count)
{
	std::ifstream inFile(filePath, std::ios::in);

	if(!inFile)
	{
		std::cerr << "ERROR opening vol data file" << std::endl;
		return(0);
	}	
	UINT8 min= 0xFF, max= 0, v;
	float t=0.0, mean;

	inFile.read((char*)pB, count);
	for (int i=0; i < count; i++)
	{	
		v= pB[i];
		t+=	 v;
		if (v > max) max=	v;		
		if (v < min) min=	v;
	}			
		
	inFile.close();
	mean= t / count;
	return(count);
} 


bool setupStructureData (void)
{
	const ByteInterval	*pI=	ivlTable;
	ByteArray			*pF= NULL;

	if (loadRaw("s1centre(200,200,200).u8", (UINT8*)(gByteArray.b), sizeof(gByteArray)))
	{
		return true;
	}
	else
	std::cerr << "ERROR opening vol data file" << std::endl;
	return false;
}

struct V3U32
{
	V3U32 (UINT tx=0, UINT ty=0, UINT tz=0) {x= tx; y= ty; z= tz;}
	UINT x, y, z;
};

UINT indexV3U32 (const UINT x, const UINT y, const UINT z, const V3U32& dim)
{
	return(x + dim.x * (y + dim.y * z));

}

void createMicroResources()
{
	cl_int errNum;

	printf("Creating microFun resources..\n");

    BioBufType *pBiomassBuffer = new BioBufType[NUM_POINTS];
	
	//setUpCarbonMap();
	//setupStructureData();

	scalar carbon=0.0;

	// initialize each cell to a random state
	for( int z = 0; z < BUFFER_SIZE_Z; z++ )
	{
		for( int y = 0; y < BUFFER_SIZE_Y; y++ )
		{
			for ( int x = 0; x < BUFFER_SIZE_X; x++)
			{
				int i=x + (y * BUFFER_SIZE_X) + (z * BUFFER_SIZE_X*BUFFER_SIZE_Y);
				
				pBiomassBuffer[i].state = 0.0;//rand()%3==0?1:0;
				pBiomassBuffer[i].mb = 0.0;//rand()%3==0?1:0;
				pBiomassBuffer[i].ins = 0.0;//rand()%3==0?1:0;
				pBiomassBuffer[i].ex_re = 265.25;//carbonMap[x][y][z]*1000000;//carbonMap[x][y][z];//rand()%3==0?1:0;
				//conservation+=pBiomassBuffer[i].ex_re;
				carbon+=pBiomassBuffer[i].ex_re;

			}
		}
	}

	double inn=0;
	printf(" TOTAL CARBON %f \n", carbon/1000000 );
	scalar conservation=0.0;
	scalar inn_vxl=4.755111745e-06;//0.05/(BUFFER_SIZE_Y-2,BUFFER_SIZE_X-2);
	int por=0;
	// Initialise some growth by creating a sphere of activity in the volume
	int radius = 6;
	int r2 = radius * radius; 
	int cx=BUFFER_SIZE_X/2, cy=BUFFER_SIZE_Y/2, cz = BUFFER_SIZE_Z/2;

	printf("range: %d %d\n", cx - 10, cx + 10);
	for( int x = cx - 10; x < cx + 10; x++ )
	{
		for( int y = cx -10; y <  cx + 10; y++ )
		{
			for(int z = cz -10; z < cz + 10; z++)
			{
				const int i= x +(y * BUFFER_SIZE_X )+ (z * BUFFER_SIZE_Y*BUFFER_SIZE_Z);
				const int dx= x - cx;
				const int dy= y - cy;
				const int dz = z - cz;
				const bool inSphere = ((dx * dx + dy * dy + dz * dz) <= r2);

				if (inSphere)
				{

					pBiomassBuffer[i].mb = 0.05*inn_vxl*1000000;//rand()%3==0?1:0;
					conservation+=pBiomassBuffer[i].mb;
					pBiomassBuffer[i].state = 0.90*inn_vxl*1000000;//rand()%3==0?1:0;
					conservation+=pBiomassBuffer[i].state;
					pBiomassBuffer[i].ins =  0.05*inn_vxl*1000000;//rand()%3==0?1:0;
					conservation+=pBiomassBuffer[i].ins;
					pBiomassBuffer[i].ex_re = 265.2500e-06*1000000;
					//conservation+=pBiomassBuffer[i].ex_re;
				}
			}
		}
	}
	printf(" cons %f \n", conservation );
	printf(" INN %f \n", inn/1000000 );
	printf(" PORES %i \n", por);

	// output initial sum to check
	BioBufType sum;
	sum.state = 0.0;
	sum.mb = 0.0;
	sum.ins = 0.0;
	sum.ex_re = 0.0;
	for(int x =0; x<BUFFER_SIZE_X; x++)
	{
		for(int y = 0; y < BUFFER_SIZE_Y; y++)
		{
			for(int z = 0; z < BUFFER_SIZE_Z; z++)
			{
				sum.state += pBiomassBuffer[x+y*BUFFER_SIZE_X+z*BUFFER_SIZE_X*BUFFER_SIZE_Y].state;
				sum.mb += pBiomassBuffer[x+y*BUFFER_SIZE_X+z*BUFFER_SIZE_X*BUFFER_SIZE_Y].mb;
				sum.ins += pBiomassBuffer[x+y*BUFFER_SIZE_X+z*BUFFER_SIZE_X*BUFFER_SIZE_Y].ins;
				sum.ex_re += pBiomassBuffer[x+y*BUFFER_SIZE_X+z*BUFFER_SIZE_X*BUFFER_SIZE_Y].ex_re;
			}
		}
	}

	// output results of sum
	printf("After initialisation - CPU..\n");
	sum.print();

	// create OpenCL buffer
	bio_in = clCreateBuffer(
		context,
		CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
		sizeof(BioBufType)*NUM_POINTS,
		pBiomassBuffer,
		&errNum);
	if( errNum != CL_SUCCESS )
	{
		std::cerr<< "Failed creating biomass input buffer." << std::endl;
	}

	// create OpenCL buffer
	bio_out = clCreateBuffer(
			context,
			CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
			sizeof(BioBufType)*NUM_POINTS,
			pBiomassBuffer,
			&errNum);
	if( errNum != CL_SUCCESS )
	{
		std::cerr<< "Failed creating biomass output buffer" << std::endl;
	}

	// create OpenCL buffer
	bio_results = clCreateBuffer(
			context,
			CL_MEM_READ_WRITE,
			sizeof(BioBufType)*(BUFFER_SIZE_X/THREAD_X)*((BUFFER_SIZE_Y/THREAD_Y)/2)*(BUFFER_SIZE_Z/THREAD_Z),
			NULL,
			&errNum);
	if( errNum != CL_SUCCESS )
	{
		std::cerr<< "Failed creating biomass results buffer" << std::endl;
	}

	// soil structure is just a smaller cube than the simulation area
	// could set as a maze like structure
	int* soil_struct_init = new int [NUM_POINTS];
	for(int x =0; x<BUFFER_SIZE_X; x++)
	{
		for(int y = 0; y < BUFFER_SIZE_Y; y++)
		{
			for(int z = 0; z < BUFFER_SIZE_Z; z++)
			{
				const int i= x +(y * BUFFER_SIZE_X )+ (z * BUFFER_SIZE_Y*BUFFER_SIZE_Z);
			
				soil_struct_init[i] = 0;
			}
		}
	}
	for(int x = 50; x<150; x++)
	{
		for(int y = 50; y < 150; y++)
		{
			for(int z = 50; z < 150; z++)
			{
				const int i= x +(y * BUFFER_SIZE_X )+ (z * BUFFER_SIZE_Y*BUFFER_SIZE_Z);
			
				soil_struct_init[i] = 1;
			}
		}
	}

	// create OpenCL buffer
	structure_tex = clCreateBuffer(
			context,
			CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, 
			sizeof(int)*NUM_POINTS,
			soil_struct_init,
			&errNum);
	if( errNum != CL_SUCCESS )
	{
		std::cerr<< "Failed creating structure buffer" << std::endl;
	}

    delete [] pBiomassBuffer;
}

// output structure to colour buffer to render
cl_int compute_structure()
{
	cl_int errNum;

    errNum = clSetKernelArg(structure_kernel, 0, sizeof(cl_mem), &cl_colour);
	errNum |= clSetKernelArg(structure_kernel, 1, sizeof(cl_mem), &structure_tex);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error setting kernel arguments." << std::endl;
    }

	// define work groups
    size_t globalWorkSize[1] = { NUM_POINTS };
    size_t localWorkSize[1] = { 50 };

    glFinish();

	errNum = clEnqueueAcquireGLObjects(commandQueue, 1, &cl_colour, 0, NULL, NULL);//&prof_ev );
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error - couldn't acquire GL tex." << std::endl;
    }
	clFinish(commandQueue);

    errNum = clEnqueueNDRangeKernel(commandQueue, structure_kernel, 1, NULL,
                                    globalWorkSize, localWorkSize,
                                    0, NULL, NULL);//&prof_ev);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error queuing structure kernel for execution." << std::endl;
	}
	clFinish(commandQueue);

	errNum = clEnqueueReleaseGLObjects(commandQueue, 1, &cl_colour, 0, NULL, NULL);//&prof_ev );
	clFinish(commandQueue);
	return 0;
}

// create colour and position buffer for render output
void initGL(int argc, char** argv)
{
	imWidth = 512; 
	imHeight = 512;

	initGlut(argc, argv, imWidth, imHeight);

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
}

int main(int argc, char** argv)
{
	printf("simulating for %d values..\n", NUM_POINTS);

	cl_device_id device = 0;

	initGL(argc, argv);

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
    program = CreateProgram(context, device, "microFun.cl");
    if (program == NULL)
    {
        Cleanup();
        return 1;
    }
/*		
    // Make sure the device supports images, otherwise exit
    cl_bool imageSupport = CL_FALSE;
    clGetDeviceInfo(device, CL_DEVICE_IMAGE_SUPPORT, sizeof(cl_bool),
                    &imageSupport, NULL);
    if (imageSupport != CL_TRUE)
    {
        std::cerr << "OpenCL device does not support images." << std::endl;
        Cleanup();
        return 1;
    }
	size_t return_val;
	clGetDeviceInfo(device, CL_DEVICE_IMAGE3D_MAX_WIDTH, sizeof(return_val), &return_val, NULL);
	std::cout << "max 3d width " << return_val << std::endl;
	clGetDeviceInfo(device, CL_DEVICE_IMAGE3D_MAX_HEIGHT, sizeof(return_val), &return_val, NULL);
	std::cout << "max 3d height " << return_val << std::endl;
	clGetDeviceInfo(device, CL_DEVICE_IMAGE3D_MAX_DEPTH, sizeof(return_val), &return_val, NULL);
	std::cout << "max 3d depth " << return_val << std::endl;
*/
	createMicroResources();

	cl_int errNum;

	// create clgl sharing resource
	cl_colour = clCreateFromGLBuffer(
		context,
		CL_MEM_WRITE_ONLY,
		m_cbo,
		&errNum);
	if( errNum != CL_SUCCESS )
	{
		std::cerr<< "Failed creating colour output buffer" << std::endl;
	}

	// Create OpenCL kernel
    kernel = clCreateKernel(program, "micro_output", NULL);
    if (kernel == NULL)
    {
        std::cerr << "Failed to create output kernel" << std::endl;
        Cleanup();
        return 1;
    }

	// Create OpenCL kernel
    sum_kernel = clCreateKernel(program, "micro_sum", NULL);
    if (sum_kernel == NULL)
    {
        std::cerr << "Failed to create sum kernel" << std::endl;
        Cleanup();
        return 1;
    }

	// Create OpenCL kernel
    diffuse_kernel = clCreateKernel(program, "micro_diffuse", NULL);
    if (diffuse_kernel == NULL)
    {
        std::cerr << "Failed to create sum kernel" << std::endl;
        Cleanup();
        return 1;
    }

	// Create OpenCL kernel
    update_kernel = clCreateKernel(program, "micro_update", NULL);
    if (update_kernel == NULL)
    {
        std::cerr << "Failed to create sum kernel" << std::endl;
        Cleanup();
        return 1;
    }

	// Create OpenCL kernel
    structure_kernel = clCreateKernel(program, "micro_structure", NULL);
    if (structure_kernel == NULL)
    {
        std::cerr << "Failed to create structure kernel" << std::endl;
        Cleanup();
        return 1;
    }
/*
	printf("OpenCL resources created...\nInitial RD data\nCPU\n");
	checkMicroSum(bio_in, BUFFER_SIZE_X, BUFFER_SIZE_Y, BUFFER_SIZE_Z);
	printf("GPU\n");
	compute_microSum();
*/

	printf("Beginning Simulation...\n");
	sim_start = get_time();

	glutMainLoop();

	Cleanup();

	return 0;
}
