//
// Book:      OpenCL(R) Programming Guide
// Authors:   Aaftab Munshi, Benedict Gaster, Timothy Mattson, James Fung, Dan Ginsburg
// ISBN-10:   0-321-74964-2
// ISBN-13:   978-0-321-74964-2
// Publisher: Addison-Wesley Professional
// URLs:      http://safari.informit.com/9780132488006/
//            http://www.openclprogrammingguide.com
//


// GLinterop.cpp
//
//    This is a simple example that demonstrates basic OpenCL setup and
//    use.

#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <fstream>
#include <sstream>
#include "time.h"

#ifdef _WIN32
#include <windows.h>
#endif
#include <GL/glew.h>
#include <GL/wglew.h>
#include <GL/glut.h>

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#include <CL/cl_gl.h>
#endif


const int width = 512;
const int height = 512;

const int ThreadsX = 16;
const int ThreadsY = 16;

const int iterations = 100;

// view params
int ox, oy;
int buttonState = 0;
float camera_trans[] = {0, 0, -3};
float camera_rot[]   = {0, 0, 0};
float camera_trans_lag[] = {0, 0, -3};
float camera_rot_lag[] = {0, 0, 0};
const float inertia = 0.1f;
float modelView[16];

std::vector<float> results;

struct float2
{
	float x;
	float y;
};

///
// OpenGL variables, objects
//
GLuint tex = 0;

///
// OpenCL variables, objects
//
cl_kernel kernel;
cl_kernel tex_kernel;
cl_kernel sum_kernel;
cl_mem rd_in;
cl_mem rd_out;
cl_mem rd_results;
cl_mem cl_tex_mem;
cl_context context;
cl_command_queue commandQueue;
cl_program program;

bool data_out = false;

///
// Forward declarations
void Cleanup();
cl_int computeVBO();
cl_int computeTexture();


///
// Display the texture in the window
//
void displayTexture(int w, int h)
{
        glEnable(GL_TEXTURE_RECTANGLE_ARB);
        glBindTexture(GL_TEXTURE_RECTANGLE_ARB, tex );
        glBegin(GL_QUADS);
                glTexCoord2f(0, 0);
                glVertex2f(-1, -1);
                glTexCoord2f(0, h);
                glVertex2f(-1, 1);
                glTexCoord2f(w, h);
                glVertex2f(1, 1);
                glTexCoord2f(w, 0);
                glVertex2f(1, -1);
        glEnd();
        glDisable(GL_TEXTURE_RECTANGLE_ARB);
}

void reshape(int w, int h) 
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(60.0, (float) w / (float) h, 0.1, 100.0);

    glMatrixMode(GL_MODELVIEW);
    glViewport(0, 0, w, h);
}

float computeStats(float ms)
{
	float numberoperations = (float) width * (float)height;
	numberoperations *= 1e-6 * 1000 / ms;
	//std::cout << "Number of operations per second : " << numberoperations << " million"<< std::endl;
	return numberoperations;
}

cl_int checkSum(cl_mem rd_gpu, const int w, const int h)
{
	cl_int errNum;
	float2* rd_cpu = new float2[w*h];

	// read gpu buffer into cpu memory
	errNum = clEnqueueReadBuffer(
		commandQueue, 
		rd_gpu, 
		CL_TRUE, 
		0, 
		sizeof(float2)*w*h, 
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
			sum.x += rd_cpu[x+y*w].x;
			sum.y += rd_cpu[x+y*w].y;
		}
	}

	// output result
	std::cout << "Sum of U component : " << sum.x << std::endl;
	std::cout << "Sum of V component : " << sum.y << std::endl;
	
	delete [] rd_cpu;

	return errNum;
}

float completeEvent(cl_event ev, cl_command_queue cq)
{
	cl_int errNum;

	cl_ulong ev_start_time = (cl_ulong)0;
	cl_ulong ev_end_time = (cl_ulong)0;
	
	clFinish(cq);
	errNum = clGetEventProfilingInfo(ev, CL_PROFILING_COMMAND_QUEUED, 				
		sizeof(cl_ulong), &ev_start_time, NULL);
	errNum |= clGetEventProfilingInfo(ev, CL_PROFILING_COMMAND_END, 					
		sizeof(cl_ulong), &ev_end_time, NULL);
	if(errNum != CL_SUCCESS)
		std::cout << "Error profiling event : " << errNum << std::endl;

	float ms = (float)((ev_end_time - ev_start_time)*1.0e-6);

	return ms;
}

cl_int computeSum()
{
	cl_int errNum;
	
    errNum = clSetKernelArg(sum_kernel, 0, sizeof(cl_mem), &rd_out);
    errNum |= clSetKernelArg(sum_kernel, 1, sizeof(float2) * ThreadsX * ThreadsY, NULL);
	errNum |= clSetKernelArg(sum_kernel, 2, sizeof(cl_mem), &rd_results);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error setting kernel arguments." << std::endl;
    }

    size_t localWorkSize[1] = { ThreadsX*ThreadsY } ;
    size_t globalWorkSize[1] =  {  width*height/2 };

	cl_event prof_event;
    errNum = clEnqueueNDRangeKernel(commandQueue, sum_kernel, 1, NULL,
                                    globalWorkSize, localWorkSize,
                                    0, NULL, &prof_event);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error queuing conway kernel for execution." << std::endl;
    }
	std::cout << "Time to find partial sums : " << completeEvent(prof_event, commandQueue) << " ms\n";

	checkSum(rd_results, width/ThreadsX, height/ThreadsY/2);
    return 0;
}

///
// Main rendering call for the scene
//
void renderScene(void)
{
	computeVBO();
	computeTexture();

	// output U and V summed values
	if(data_out)
	{
		// rd_out is to be rendered so it is output
		std::cout << "Output UV data" << std::endl;
		std::cout << "CPU" << std::endl;
		checkSum(rd_out, width, height);
		std::cout << "GPU" << std::endl;
		computeSum();
		std::cout << std::endl;
		data_out = false;
	}

	glClearColor(0.0f, 0.0f, 1.0f, 1.0f );
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

    // render simulation area
    glColor3f(1.0, 1.0, 1.0);
    glutWireCube(2.0);

	// render simulation data
	displayTexture(width, height);
	glutSwapBuffers();	
	std::swap(rd_in, rd_out);
}

///
// Keyboard events handler
//
void KeyboardGL(unsigned char key, int x, int y)
{
    switch(key) 
    {
        case '\033': // escape quits
        case '\015': // Enter quits    
        case 'Q':    // Q quits
        case 'q':    // q (or escape) quits
            // Cleanup up and quit
                Cleanup();
            break;
		case 'o':	// O outputs profile data
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
void initGlut(int argc, char *argv[], int wWidth, int wHeight)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
    //glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowPosition(100,100);
    glutInitWindowSize(wWidth, wHeight);
    glutCreateWindow("Grey Scott");

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

	const GLubyte* strVersion=0;

	if ( (strVersion=glGetString(GL_VERSION)) ) {
	printf("OpenGL version: %s\n", (char*)strVersion);
	}
	else { 
	printf("glGetString returned 0\n");
	}

	strVersion = 0;

	wglSwapIntervalEXT(false); //disable vsync for faster rendering
	glEnable(GL_DEPTH_TEST);
	glutReportErrors();
}

void initTexture( int width, int height )
{
    // make a texture for output
    glGenTextures(1, &tex);              // texture 
    glTexEnvi( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE,  GL_REPLACE );
    glBindTexture(GL_TEXTURE_RECTANGLE_ARB, tex);
    glTexImage2D(GL_TEXTURE_RECTANGLE_ARB, 0, GL_RGBA32F_ARB, width,
            height, 0, GL_LUMINANCE, GL_FLOAT, NULL );
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

cl_int computeVBO()
{
	cl_int errNum;

    errNum = clSetKernelArg(kernel, 0, sizeof(cl_mem), &rd_in);
	errNum |= clSetKernelArg(kernel, 1, sizeof(int), &width);
	errNum |= clSetKernelArg(kernel, 2, sizeof(int), &height);
    errNum |= clSetKernelArg(kernel, 3, sizeof(float2) * ThreadsX * ThreadsY, NULL);
	errNum |= clSetKernelArg(kernel, 4, sizeof(cl_mem), &rd_out);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error setting kernel arguments." << std::endl;
    }

	int globalX = ((width - 1)/(ThreadsX - 2) + 1) * ThreadsX;
	int globalY = ((height - 1)/(ThreadsY - 2) + 1) * ThreadsY;

    size_t localWorkSize[2] = { ThreadsX, ThreadsY } ;
    size_t globalWorkSize[2] =  {  globalX, globalY };

	//std::cout << "compute vbo \n";
	//std::cout << globalX << " " << globalY  << std::endl;
	//std::cout << ThreadsX << " " << ThreadsY << std::endl;

	cl_event prof_event;
    errNum = clEnqueueNDRangeKernel(commandQueue, kernel, 2, NULL,
                                    globalWorkSize, localWorkSize,
                                    0, NULL, &prof_event);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error queuing conway kernel for execution." << std::endl;
    }

	results.push_back( completeEvent(prof_event, commandQueue) );

    return 0;
}


///
// Use OpenCL to compute the colors on the texture background
// Bascially the same functionality as the VBO, execpt operating on a texture object
cl_int computeTexture()
{
	cl_int errNum;

    errNum = clSetKernelArg(tex_kernel, 0, sizeof(cl_mem), &cl_tex_mem);
	errNum = clSetKernelArg(tex_kernel, 1, sizeof(cl_mem), &rd_out);

	// define work groups
	size_t localWorkSize[2] = { ThreadsX, ThreadsY } ;               
    size_t globalWorkSize[2] = {  width, height };

	//std::cout << "compute tex \n";
	///std::cout << width << " " << height  << std::endl;
	//std::cout << ThreadsX << " " << ThreadsY << std::endl;

   glFinish();

	errNum = clEnqueueAcquireGLObjects(commandQueue, 1, &cl_tex_mem, 0, NULL, NULL);//&prof_ev );
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error - couldn't acquire GL tex." << std::endl;
    }
	clFinish(commandQueue);


    errNum = clEnqueueNDRangeKernel(commandQueue, tex_kernel, 2, NULL,
                                    globalWorkSize, localWorkSize,
                                    0, NULL, NULL);//&prof_ev);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error queuing texture kernel for execution." << std::endl;
	}
	clFinish(commandQueue);

	errNum = clEnqueueReleaseGLObjects(commandQueue, 1, &cl_tex_mem, 0, NULL, NULL);//&prof_ev );
	clFinish(commandQueue);

	return 0;
}

///
//  Cleanup any created OpenCL resources
//
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
        
		if( tex_kernel != 0 ) 
                clReleaseKernel(tex_kernel);

		if( sum_kernel != 0 ) 
                clReleaseKernel(sum_kernel);

        if( rd_in != 0 )
                clReleaseMemObject(rd_in);

		if( rd_out != 0 )
                clReleaseMemObject(rd_out);

		if( rd_results != 0 )
                clReleaseMemObject(rd_results);

		if( cl_tex_mem != 0 )
                clReleaseMemObject(cl_tex_mem);
		
		// after we have released the OpenCL references, we can delete the underlying OpenGL objects
		if( tex != 0 ) 
		{
				glBindBuffer(GL_TEXTURE_RECTANGLE_ARB, tex );
				glDeleteBuffers(1, &tex);
		}
		
//		system("pause");
	exit(0);
}


bool isSupported(const char* extension, std::string list)
{
	size_t found = list.find(extension);
	if(found == std::string::npos)
		return false;

	return true;
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

///
//  Create an OpenCL context on the first available platform using
//  either a GPU or CPU depending on what is available.
//
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


    // Next, create an OpenCL context on the platform.  Attempt to
    // create a GPU-based context, and if that fails, try to create
    // a CPU-based context.
    cl_context_properties contextProperties[] =
    {
        CL_CONTEXT_PLATFORM,
        (cl_context_properties)firstPlatformId,
                CL_GL_CONTEXT_KHR,
                (cl_context_properties)wglGetCurrentContext(),
                CL_WGL_HDC_KHR,
                (cl_context_properties)wglGetCurrentDC(),
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

void CreateIOResources()
{
	cl_int errNum;

	printf("Number of RD Points :\t%d\n", width*height);
    float2 *pGSBuffer = new float2[width * height];
 
   // initialize each cell to a random state
    for( int x = 0; x < width; x++ )
    {
        for( int y = 0; y < height; y++ )
        {
            pGSBuffer[ x + y * width ].x = 1.0;
		    pGSBuffer[ x + y * width ].y = 0.0;
        }
    }

	  // initialize each cells - pattern formation is sensitive to initial state, vlues below are known to work
    for( int x = width/3; x <  2*width/3; x++ )
    {
        for( int y = height/3; y < 2*height/3; y++ )
        {
            pGSBuffer[ x + y * width ].x =0.5;
		    pGSBuffer[ x + y * width ].y =0.25;
	     }
    }
	
	// Now perturb the entire grid. Bound the values by [0,1]
    for (int k = 0; k < width * height; ++k)
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
	for(int x =0; x<width; x++)
	{
		for(int y = 0; y < height; y++)
		{
			sum.x += pGSBuffer[x+y*width].x;
			sum.y += pGSBuffer[x+y*width].y;
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
		sizeof(float2)*width*height,
		pGSBuffer,
		&errNum);
	if( errNum != CL_SUCCESS )
	{
		std::cerr<< "Failed creating rd input buffer." << std::endl;
	}

	rd_out = clCreateBuffer(
		context,
		CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
		sizeof(float2)*width*height,
		pGSBuffer,
		&errNum);
	if( errNum != CL_SUCCESS )
	{
		std::cerr<< "Failed creating rd output buffer" << std::endl;
	}

    delete [] pGSBuffer;

	// create buffer to store preliminary sum results
	// each work group will produce one results value
	// (divide by 2 as the kernel performs first add during load)
	rd_results = clCreateBuffer(
		context,
		CL_MEM_READ_WRITE,
		sizeof(float2)*(width/ThreadsX)*(height/ThreadsY)/2,
		NULL,
		&errNum);
	if( errNum != CL_SUCCESS )
	{
		std::cerr<< "Failed creating results buffer" << std::endl;
	}

	// create image
	cl_tex_mem = clCreateFromGLTexture2D(
		context,
		CL_MEM_READ_WRITE,
		GL_TEXTURE_RECTANGLE_ARB,
		0,
		tex,
		&errNum);
	if( errNum != CL_SUCCESS )
	{
		std::cerr<< "Failed creating shared GL image." << std::endl;
	}

}

cl_int checkWGSize(cl_kernel k, cl_device_id d)
{
	cl_int errNum;
	cl_int info;

	// query kernel for max work group size
	errNum = clGetKernelWorkGroupInfo(k, 
		d, 
		CL_KERNEL_WORK_GROUP_SIZE, 
		NULL, 
		&info, 
		NULL);
	if( errNum != CL_SUCCESS )
	{
		std::cerr<< "Problem with get info. " << errNum << std::endl;
		return errNum;
	}

	cl_int threads = ThreadsX * ThreadsY;
	printf("Max work group size for kernel %d\n", info);
	printf("Attempted work group size %d\n\n", threads);

	if(threads > info)
		return -54;			// invalid work group size
	else
		return errNum;
}

///
//      main() for GLinterop example
//
int main(int argc, char** argv)
{
	printf("Grey Scott Benchmark OpenCL\n\n");

	cl_device_id device = 0;

	srand( (UINT)time( NULL ) );

	initGlut(argc, argv, width, height);
	initTexture(width,height);

	printf("GL resources created\n\n");

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
    program = CreateProgram(context, device, "greyscott.cl");
    if (program == NULL)
    {
        Cleanup();
        return 1;
    }

	// create cl input/output buffers
	CreateIOResources();
 
	// Output some useful data
	printf("Image size :\t%d %d\n", width, height);
	printf("Local Work Size :\t%d %d\n", ThreadsX, ThreadsY);
	printf("Global Work Size :\t%d %d\n\n", ((width - 1)/(ThreadsX - 2) + 1) * ThreadsX,
									  ((height - 1)/(ThreadsY - 2) + 1) * ThreadsY );

	// Create OpenCL kernel
    kernel = clCreateKernel(program, "gs_rd_kernel", NULL);
    if (kernel == NULL)
    {
        std::cerr << "Failed to create rd kernel" << std::endl;
        Cleanup();
        return 1;
    }
	cl_int errNum = checkWGSize(kernel, device);
	if( errNum != CL_SUCCESS )
	{
		std::cerr<< "Error " << errNum << std::endl;
		system("pause");
	}

	// Create OpenCL kernel
    tex_kernel = clCreateKernel(program, "output_kernel", NULL);
    if (tex_kernel == NULL)
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

	// check initial sum
	//std::cout << "checking initial sum value against initial sum" << std::endl;
	//checkSum(rd_in, width, height);

		printf("CL resources created. Entering main loop... \n");
        glutMainLoop();

	return 0;
}