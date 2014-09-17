
#include <iostream>
#include <fstream>
#include <sstream>
#include "time.h"

#ifdef _WIN32
#include <windows.h>
#endif

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

bool show_debug = true;
bool show_debug_all = false;

const int imWidth = 16;
const int imHeight = 16;
const int ThreadsX = 4;
const int ThreadsY = 4;

//
// OpenCL variables, objects
//
cl_kernel kernel;
cl_kernel copy_kernel;

// input/output buffers
cl_mem conway_in;
cl_mem conway_out;
cl_mem conway_check;	// output compare buffer

cl_mem sub_buffer_in[3];	// views into conway_in
cl_mem sub_buffer_out;		// view into conway_out for results of middle

cl_mem ghost_zones[2];		// top and bottom ghost zone view into conway_out
cl_mem check_zones[2];		// top and bottom ghost zone view into conway_check

cl_mem temp_zones[2];	// output sections for some to be copied to output and ghost transmit

cl_context context;
cl_command_queue commandQueue;
cl_program program;

void CL_CALLBACK transmit_top_ghost(cl_event my_event, cl_int cmd_exec_status, void *user_data)
{
    // handle the callback here.
	int* output = (int*) user_data;

	/*
	// To Profile Ghost Zone Read
	cl_int errNum;
	cl_ulong ev_start_time = (cl_ulong)0;
	cl_ulong ev_end_time = (cl_ulong)0;
	
	errNum = clGetEventProfilingInfo(my_event, CL_PROFILING_COMMAND_QUEUED, 				
		sizeof(cl_ulong), &ev_start_time, NULL);
	errNum |= clGetEventProfilingInfo(my_event, CL_PROFILING_COMMAND_END, 					
		sizeof(cl_ulong), &ev_end_time, NULL);
	if(errNum != CL_SUCCESS)
		std::cout << "Error profiling event : " << errNum << std::endl;

	double event_time = (ev_end_time - ev_start_time)*1.0e-6;	// ms
	std::cout << event_time << " ms" << std::endl;
	*/

	std::cout << "CALLBACK triggered - send UP" << std::endl; 
	 // initialize each cell to a random state
    for( int y = 0; y < imWidth; y++ )
    {
		std::cout << output[y] << " ";
	}
	std::cout << std::endl; 

	delete [] output;	
	clReleaseEvent(my_event);
}

void CL_CALLBACK receive_top_ghost(cl_event my_event, cl_int cmd_exec_status, void *user_data)
{
    // handle the callback here.

	/*
	// To Profile Ghost Zone Write
	cl_int errNum;
	cl_ulong ev_start_time = (cl_ulong)0;
	cl_ulong ev_end_time = (cl_ulong)0;
	
	errNum = clGetEventProfilingInfo(my_event, CL_PROFILING_COMMAND_QUEUED, 				
		sizeof(cl_ulong), &ev_start_time, NULL);
	errNum |= clGetEventProfilingInfo(my_event, CL_PROFILING_COMMAND_END, 					
		sizeof(cl_ulong), &ev_end_time, NULL);
	if(errNum != CL_SUCCESS)
		std::cout << "Error profiling event : " << errNum << std::endl;

	double event_time = (ev_end_time - ev_start_time)*1.0e-6;	// ms
	std::cout << event_time << " ms" << std::endl;
	*/

	std::cout << "CALLBACK triggered - receive top" << std::endl; 

	clReleaseEvent(my_event);
}

void CL_CALLBACK receive_bottom_ghost(cl_event my_event, cl_int cmd_exec_status, void *user_data)
{
    // handle the callback here.

	/*
	// To Profile Ghost Zone Write
	cl_int errNum;
	cl_ulong ev_start_time = (cl_ulong)0;
	cl_ulong ev_end_time = (cl_ulong)0;
	
	errNum = clGetEventProfilingInfo(my_event, CL_PROFILING_COMMAND_QUEUED, 				
		sizeof(cl_ulong), &ev_start_time, NULL);
	errNum |= clGetEventProfilingInfo(my_event, CL_PROFILING_COMMAND_END, 					
		sizeof(cl_ulong), &ev_end_time, NULL);
	if(errNum != CL_SUCCESS)
		std::cout << "Error profiling event : " << errNum << std::endl;

	double event_time = (ev_end_time - ev_start_time)*1.0e-6;	// ms
	std::cout << event_time << " ms" << std::endl;
	*/

	std::cout << "CALLBACK triggered - receive bottom" << std::endl; 

	clReleaseEvent(my_event);
}

// same as transmit_top_ghost - could include tag in user data to determine where to send
void CL_CALLBACK transmit_bottom_ghost(cl_event my_event, cl_int cmd_exec_status, void *user_data)
{
    // handle the callback here.
	int* output = (int*) user_data;

	/*
	// To Profile Ghost Zone Read
	cl_int errNum;
	cl_ulong ev_start_time = (cl_ulong)0;
	cl_ulong ev_end_time = (cl_ulong)0;
	
	errNum = clGetEventProfilingInfo(my_event, CL_PROFILING_COMMAND_QUEUED, 				
		sizeof(cl_ulong), &ev_start_time, NULL);
	errNum |= clGetEventProfilingInfo(my_event, CL_PROFILING_COMMAND_END, 					
		sizeof(cl_ulong), &ev_end_time, NULL);
	if(errNum != CL_SUCCESS)
		std::cout << "Error profiling event : " << errNum << std::endl;

	double event_time = (ev_end_time - ev_start_time)*1.0e-6;	// ms
	std::cout << event_time << " ms" << std::endl;
	*/

	std::cout << "CALLBACK triggered - send DOWN" << std::endl; 
	// initialize each cell to a random state
    for( int y = 0; y < imWidth; y++ )
	{
		std::cout << output[y] << " ";
	}
	std::cout << std::endl; 

	delete [] output;	
	clReleaseEvent(my_event);
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

	return (float)((ev_end_time - ev_start_time)*1.0e-6);
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

        if( copy_kernel != 0 ) 
                clReleaseKernel(copy_kernel);

        if( conway_in != 0 )
                clReleaseMemObject(conway_in);

		if( conway_out != 0 )
                clReleaseMemObject(conway_out);

		if( conway_check != 0 )
                clReleaseMemObject(conway_check);

		if( temp_zones[0] != 0 )
                clReleaseMemObject(temp_zones[0]);

		if( temp_zones[1] != 0 )
                clReleaseMemObject(temp_zones[1]);
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
    //Source = srcStdStr.c_str();

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

cl_context CreateContext()
{
	cl_int errNum;
	cl_uint numPlatforms;
	cl_platform_id* platformIDs;
	//cl_platform_id firstPlatformId;
	//cl_context context = NULL;
	cl_context* context;

	errNum = clGetPlatformIDs(0, NULL, &numPlatforms);
	if(errNum != CL_SUCCESS || numPlatforms <=0)
	{
		std::cerr << "Failed to find any OpenCl platforms." << std::endl;
		return NULL;
	}
	//std::cout << "Total number of Platforms available : " << numPlatforms << std::endl;

	platformIDs = (cl_platform_id*) alloca(sizeof(cl_platform_id)*numPlatforms);
	
	errNum = clGetPlatformIDs(numPlatforms, platformIDs, NULL);
	if(errNum != CL_SUCCESS)
	{
		std::cerr << errNum << " - Error getting platforms" << std::endl;
		return NULL;
	}

	// Assign memory to context pointer for all platform
	context = (cl_context *)alloca(sizeof(cl_context) * numPlatforms);
	for(cl_uint i = 0; i < numPlatforms; i++)
	{
		//std::cout << "Trying to create context for platform " << i << ".." << std::endl;
		cl_context_properties properties[] =
		{
			CL_CONTEXT_PLATFORM, (cl_context_properties)platformIDs[i], 0
		};

		context[i] = clCreateContextFromType(properties, CL_DEVICE_TYPE_ALL, NULL, NULL, &errNum);
		if(errNum != CL_SUCCESS)
		{
			std::cerr << errNum << " - Error" << std::endl;
			return NULL;
		}
	}

	return context[0];
}

int* get_new_ghost_row()
{
	// create state buffer
	int *pConwayBuffer = new int[imWidth];

	std::cout << "Receiving incoming ghost zone row.." << std::endl; 
	 // initialize each cell to a random state
    for( int x = 0; x < imWidth; x++ )
    {
            pConwayBuffer[x] = rand()%3==0?1:0;
if(show_debug)
{
			std::cout << pConwayBuffer[x] << " ";
}
    }
if(show_debug)
{
	std::cout << std::endl; 
}
	return pConwayBuffer;
}

// simulate receiving top ghost zone from MPI
void receive_top()
{
	cl_int errNum;

	// receive ghost zone
	int* ghost = get_new_ghost_row();
	//print_new_ghost_row(ghost);

	cl_event write_event;
	errNum = clEnqueueWriteBuffer(commandQueue,
		 ghost_zones[0],
		 CL_FALSE,
		 0,
		 sizeof(int)*imWidth,
		 ghost,
		 0,
		 NULL,
		 &write_event);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error writing ghost zone to sub-buffer." << std::endl;
    }

	errNum = clSetEventCallback(write_event, CL_COMPLETE, &receive_top_ghost, NULL);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error setting event callback (write)." << std::endl;
    }


	errNum = clEnqueueWriteBuffer(commandQueue,
		 check_zones[0],
		 CL_TRUE,
		 0,
		 sizeof(int)*imWidth,
		 ghost,
		 0,
		 NULL,
		 NULL);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error writing ghost zone to sub-buffer." << std::endl;
    }

	delete [] ghost;
}

// simulate receiving bottom ghost zone from MPI
// same as receive_top just with different sub-buffer
// still using blocking
void receive_bottom()
{
	cl_int errNum;

	// receive ghost zone
	int* ghost = get_new_ghost_row();
	//print_new_ghost_row(ghost);

	cl_event write_event;
	errNum = clEnqueueWriteBuffer(commandQueue,
		 ghost_zones[1],
		 CL_FALSE,
		 0,
		 sizeof(int)*imWidth,
		 ghost,
		 0,
		 NULL,
		 &write_event);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error writing ghost zone to sub-buffer." << std::endl;
    }

	errNum = clSetEventCallback(write_event, CL_COMPLETE, &receive_bottom_ghost, NULL);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error setting event callback (write)." << std::endl;
    }

	errNum = clEnqueueWriteBuffer(commandQueue,
		 check_zones[1],
		 CL_TRUE,
		 0,
		 sizeof(int)*imWidth,
		 ghost,
		 0,
		 NULL,
		 NULL);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error writing ghost zone to sub-buffer." << std::endl;
    }

	delete [] ghost;
}

void print_new_ghost_row(int* pConwayBuffer)
{
	std::cout << "Printing incoming ghost zone row.." << std::endl; 
	 // initialize each cell to a random state
    for( int x = 0; x < imWidth; x++ )
    {
		std::cout << pConwayBuffer[x ] << " ";
    }
	std::cout << std::endl; 
}

void CreateIOResources()
{
	cl_int errNum;

	// create state buffer
	int *pConwayBuffer = new int[imWidth*imHeight];

	//std::cout << "Creating conway input buffer.." << std::endl; 
	 // initialize each cell to a random state
    for( int x = 0; x < imWidth; x++ )
    {
        for( int y = 0; y < imHeight; y++ )
        {
            pConwayBuffer[y + x * imWidth ] = rand()%3==0?1:0;
			//pConwayBuffer[ y + x * imWidth ] = y + x * imWidth;	// give index to debug
			//std::cout << pConwayBuffer[y + x * imWidth] << " ";
		}
		//std::cout << std::endl; 
    }
	//std::cout << std::endl; 

	// create cl buffers
	conway_in = clCreateBuffer(
		context,
		CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
		sizeof(int)*imWidth*imHeight,
		pConwayBuffer,
		&errNum);
	if( errNum != CL_SUCCESS )
	{
		std::cerr << errNum << std::endl;
		std::cerr<< "Failed creating memory from buffer." << std::endl;
	}
	delete [] pConwayBuffer;

if(show_debug_all)	// check input buffer is set as pConwayBuffer
{
	int* input = new int[imWidth*imHeight];
	errNum = clEnqueueReadBuffer(commandQueue, 
		conway_in,
		CL_TRUE,
		0,
		sizeof(int) * imWidth * imHeight,
		input,
		0,
		NULL,
		NULL);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error reading top." << std::endl;
    }

	std::cout << "Checking input buffer.." << std::endl; 
	 // initialize each cell to a random state
    for( int x = 0; x < imWidth; x++ )
    {
        for( int y = 0; y < imHeight; y++ )
        {
			std::cout << input[y + x * imWidth] << " ";		// instead of to screen we could compare with pConwayBuffer directly and output any errors
		}
		std::cout << std::endl; 
    }
	std::cout << std::endl; 
	delete [] input;
}

	// output buffer
	conway_out = clCreateBuffer(
		context,
		CL_MEM_READ_WRITE,
		sizeof(int)*imWidth*imHeight,
		NULL,
		&errNum);
	if( errNum != CL_SUCCESS )
	{
		std::cerr << errNum << std::endl;
		std::cerr<< "Failed creating memory from buffer." << std::endl;
	}
	// output buffer
	conway_check = clCreateBuffer(
		context,
		CL_MEM_READ_WRITE,
		sizeof(int)*imWidth*imHeight,
		NULL,
		&errNum);
	if( errNum != CL_SUCCESS )
	{
		std::cerr << errNum << std::endl;
		std::cerr<< "Failed creating memory from buffer." << std::endl;
	}	

	// temp zones contain the full amount of memory needed for compute top and bottom ghost 	
	temp_zones[0] = clCreateBuffer(
		context,
		CL_MEM_READ_WRITE,
		sizeof(int)*imWidth*(imHeight/4),
		NULL,
		&errNum);
	if( errNum != CL_SUCCESS )
	{
		std::cerr << errNum << std::endl;
		std::cerr<< "Failed creating memory from buffer." << std::endl;
	}
	temp_zones[1] = clCreateBuffer(
		context,
		CL_MEM_READ_WRITE,
		sizeof(int)*imWidth*(imHeight/4),
		NULL,
		&errNum);
	if( errNum != CL_SUCCESS )
	{
		std::cerr << errNum << std::endl;
		std::cerr<< "Failed creating memory from buffer." << std::endl;
	}

}

cl_int compute_check()
{
	cl_int errNum;
	int width = imWidth;
	int height = imHeight;

if(show_debug_all)
{
	// read input buffer to output to screen
	int* input = new int[width*height];
	errNum = clEnqueueReadBuffer(commandQueue, 
		conway_in,
		CL_TRUE,
		0,
		sizeof(int) * width * height,
		input,
		0,
		NULL,
		NULL);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error reading input buffer." << std::endl;
    }

	clFinish(commandQueue);
	std::cout << "CHECK BEFORE" << std::endl; 
	 // initialize each cell to a random state
    for( int x = 0; x < imHeight; x++ )
    {
        for( int y = 0; y < imWidth; y++ )
        {
			std::cout << input[y + x * imWidth] << " ";
		}
		std::cout << std::endl; 
    }
	std::cout << std::endl; 
	delete [] input;
}

	// set kernel arguments to compute,
	// sub buffers require no need for offset
    errNum = clSetKernelArg(kernel, 0, sizeof(cl_mem), &conway_in);
	errNum |= clSetKernelArg(kernel, 1, sizeof(int), &width);
	errNum |= clSetKernelArg(kernel, 2, sizeof(int), &height);
    errNum |= clSetKernelArg(kernel, 3, sizeof(int) * ThreadsX * ThreadsY, NULL);
	errNum |= clSetKernelArg(kernel, 4, sizeof(cl_mem), &conway_check);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error setting kernel args." << std::endl;
    }

	int globalX = ((imWidth - 1)/(ThreadsX - 2) + 1) * ThreadsX;
	int globalY = ((imHeight - 1)/(ThreadsY - 2) + 1) * ThreadsY;

    size_t localWorkSize[2] = { ThreadsX, ThreadsY } ;
    size_t globalWorkSize[2] =  {  globalX, globalY };

if(show_debug)
{
	std::cout << "Compute Check Work Sizes" << std::endl; 
	std::cout << "Global size " << globalX << ", " << globalY << std::endl;
	std::cout << "Local size " << ThreadsX << ", " << ThreadsY << std::endl;
	std::cout << std::endl;
}

	//cl_event prof_event;
    errNum = clEnqueueNDRangeKernel(commandQueue, kernel, 2, NULL,
                                    globalWorkSize, localWorkSize,
                                    0, NULL, NULL);// &prof_event);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error queuing conway kernel top for execution." << std::endl;
    }

/*
	std::cout << "\tcompute - " << completeEvent(prof_event, commandQueue) << " ms" << std::endl;
	clReleaseEvent(prof_event);
*/
if(show_debug_all)
{
	// read output to screen
	int* output = new int[width*height];
	//cl_event prof_read;
	errNum = clEnqueueReadBuffer(commandQueue, 
		conway_check,
		CL_TRUE,
		0,
		sizeof(int) * width * height,
		output,
		0,
		NULL,
		NULL);
		//&prof_read);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error reading output." << std::endl;
    }
	
	//std::cout << "\tread - " << completeEvent(prof_read, commandQueue) << " ms" << std::endl;
	//clFinish(commandQueue);
	std::cout << "CHECK AFTER" << std::endl; 
	 // initialize each cell to a random state
    for( int x = 0; x < imHeight; x++ )
    {
        for( int y = 0; y < imWidth; y++ )
        {
			std::cout << output[y + x * imWidth] << " ";
		}
		std::cout << std::endl; 
    }
	std::cout << std::endl; 

	//clReleaseEvent(prof_read);
	delete [] output;
}

    return 0;
}

void print_sub_0()
{
if(show_debug_all)
{
	cl_int errNum;
	int width = imWidth;
	int height = (imHeight/4);

	// read input buffer to output to screen
	int* input = new int[width*height];
	errNum = clEnqueueReadBuffer(commandQueue, 
		sub_buffer_in[0],
		CL_TRUE,
		0,
		sizeof(int) * width * height,
		input,
		0,
		NULL,
		NULL);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error reading top input." << std::endl;
    }

	std::cout << "SUB 0" << std::endl; 
	 // initialize each cell to a random state
    for( int x = 0; x < imHeight/4; x++ )
    {
        for( int y = 0; y < imWidth; y++ )
        {
			std::cout << input[y + x * imWidth] << " ";
		}
		std::cout << std::endl; 
    }
	std::cout << std::endl; 

	delete [] input;
}
}

void print_sub_1()
{
if(show_debug_all)
{
	cl_int errNum;
	int width = imWidth;
	int height = ((imHeight/2) + 4);

	int* middle = new int[ width * height];
	errNum = clEnqueueReadBuffer(commandQueue, 
		sub_buffer_in[1],
		CL_TRUE,
		0,
		sizeof(int) * width * height,
		middle,
		0,
		NULL,
		NULL);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error reading middle." << std::endl;
    }

	std::cout << "SUB 1" << std::endl; 
	 // initialize each cell to a random state
    for( int x = 0; x < height; x++ )
    {
        for( int y = 0; y < width; y++ )
        {
			std::cout << middle[y + x * width] << " ";
		}
		std::cout << std::endl; 
    }
	std::cout << std::endl; 

	delete [] middle;
}
}

void print_sub_2()
{
if(show_debug_all)
{
	cl_int errNum;
	int width = imWidth;
	int height = (imHeight/4);

	int* input = new int[width*height];
	errNum = clEnqueueReadBuffer(commandQueue, 
		sub_buffer_in[2],
		CL_TRUE,
		0,
		sizeof(int) * width * height,
		input,
		0,
		NULL,
		NULL);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error reading bottom." << std::endl;
    }

	std::cout << "SUB 2" << std::endl; 
	 // initialize each cell to a random state
    for( int x = 0; x < imHeight/4; x++ )
    {
        for( int y = 0; y < imWidth; y++ )
        {
			std::cout << input[y + x * imWidth] << " ";
		}
		std::cout << std::endl; 
    }
	std::cout << std::endl; 

	delete [] input;
}
}

cl_int compute_top()
{
	cl_int errNum;
	int width = imWidth;
	int height = (imHeight/4);

	print_sub_0();

	// set kernel arguments to compute,
	// sub buffers require no need for offset
    errNum = clSetKernelArg(kernel, 0, sizeof(cl_mem), &sub_buffer_in[0]);
	errNum |= clSetKernelArg(kernel, 1, sizeof(int), &imWidth);
	errNum |= clSetKernelArg(kernel, 2, sizeof(int), &imHeight);
    errNum |= clSetKernelArg(kernel, 3, sizeof(int) * ThreadsX * ThreadsY, NULL);
	errNum |= clSetKernelArg(kernel, 4, sizeof(cl_mem), &temp_zones[0]);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error setting kernel args." << std::endl;
    }

	int globalX = ((width - 1)/(ThreadsX - 2) + 1) * ThreadsX;
	int globalY = ((height - 1)/(ThreadsY - 2) + 1) * ThreadsY;

    size_t localWorkSize[2] = { ThreadsX, ThreadsY } ;
    size_t globalWorkSize[2] =  {  globalX, globalY };

if(show_debug)
{
	std::cout << "Compute Top Work Sizes" << std::endl; 
	std::cout << "Global size " << globalX << ", " << globalY << std::endl;
	std::cout << "Local size " << ThreadsX << ", " << ThreadsY << std::endl;
	std::cout << std::endl;
}

	//cl_event prof_event;
	// compute
    errNum = clEnqueueNDRangeKernel(commandQueue, kernel, 2, NULL,
                                    globalWorkSize, localWorkSize,
                                    0, NULL, NULL);//&prof_event);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error queuing conway kernel top for execution." << std::endl;
    }

	//std::cout << "\tcompute - " << completeEvent(prof_event, commandQueue) << " ms" << std::endl;

	// since there are ghost zones involved, we only require the middle of the buffer,
	// the first and last rows are not needed
	// set offset and size as required!
	int* output = new int[width];
	cl_event read_event;
	errNum = clEnqueueReadBuffer(commandQueue, 
		temp_zones[0],
		CL_FALSE,
		sizeof(int) * width,
		sizeof(int) * width,
		output,
		0,
		NULL,
		&read_event);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error reading top output." << std::endl;
    }
/*
	// this is the same as using blocking
	errNum = clWaitForEvents(1, &read_event);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error waiting on events" << std::endl;
    }
*/
	errNum = clSetEventCallback(read_event, CL_COMPLETE, &transmit_top_ghost, output);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error setting event callback (read)." << std::endl;
    }

	//std::cout << "\tread - " << completeEvent(prof_read, commandQueue) << " ms" << std::endl;

	//clReleaseEvent(read_event);	// cant release if waiting for callback?
	//delete [] output;	// can i delete output if im waiting on an event callback?

	//clReleaseEvent(prof_event);

    return 0;
}


cl_int compute_middle()
{
	cl_int errNum;
	int width = imWidth;
	int height = (imHeight/2) + 4;

if(show_debug_all)
{
	// read input buffer to output to screen
	int* input = new int[width*height];
	errNum = clEnqueueReadBuffer(commandQueue, 
		sub_buffer_in[1],
		CL_TRUE,
		0,
		sizeof(int) * width * height,
		input,
		0,
		NULL,
		NULL);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error reading middle input." << std::endl;
    }

	clFinish(commandQueue);
	std::cout << "MIDDLE BEFORE" << std::endl; 
	 // initialize each cell to a random state
    for( int x = 0; x < height; x++ )
    {
        for( int y = 0; y < imWidth; y++ )
        {
			std::cout << input[y + x * imWidth] << " ";
		}
		std::cout << std::endl; 
    }
	std::cout << std::endl; 

	delete [] input;
}

	// set kernel arguments to compute,
	// sub buffers require no need for offset
    errNum = clSetKernelArg(kernel, 0, sizeof(cl_mem), &sub_buffer_in[1]);
	errNum |= clSetKernelArg(kernel, 1, sizeof(int), &imWidth);
	errNum |= clSetKernelArg(kernel, 2, sizeof(int), &imHeight);
    errNum |= clSetKernelArg(kernel, 3, sizeof(int) * ThreadsX * ThreadsY, NULL);
	errNum |= clSetKernelArg(kernel, 4, sizeof(cl_mem), &sub_buffer_out);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error setting kernel args." << std::endl;
    }

	int globalX = ((width - 1)/(ThreadsX - 2) + 1) * ThreadsX;
	int globalY = ((height - 1)/(ThreadsY - 2) + 1) * ThreadsY;

if(show_debug)
{
	std::cout << "Compute Middle Work Sizes" << std::endl; 
	std::cout << "Global size " << globalX << ", " << globalY << std::endl;
	std::cout << "Local size " << ThreadsX << ", " << ThreadsY << std::endl;
	std::cout << std::endl; 
}

    size_t localWorkSize[2] = { ThreadsX, ThreadsY } ;
    size_t globalWorkSize[2] =  {  globalX, globalY };

	cl_event prof_event;
	// compute
    errNum = clEnqueueNDRangeKernel(commandQueue, kernel, 2, NULL,
                                    globalWorkSize, localWorkSize,
                                    0, NULL, &prof_event);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error queuing conway kernel middle for execution." << std::endl;
    }

	//std::cout << "\tcompute - " << completeEvent(prof_event, commandQueue) << " ms" << std::endl;

if(show_debug_all)
{	// dont need to read middle buffer
	int* output = new int[width*(height-2)];
	cl_event prof_read;
	errNum = clEnqueueReadBuffer(commandQueue, 
		sub_buffer_out,
		CL_TRUE,
		sizeof(int) * width,
		sizeof(int) * width * (height - 2),
		output,
		0,
		NULL,
		&prof_read);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error reading middle output." << std::endl;
    }
	
	//std::cout << "\tread - " << completeEvent(prof_read, commandQueue) << " ms" << std::endl;
	std::cout << "MIDDLE AFTER" << std::endl; 
	 // initialize each cell to a random state
    for( int x = 0; x < (height-2); x++ )
    {
        for( int y = 0; y < width; y++ )
        {
			std::cout << output[y + x * width] << " ";
		}
		std::cout << std::endl; 
    }
	std::cout << std::endl; 

	clReleaseEvent(prof_read);

	delete [] output;
}

	clReleaseEvent(prof_event);

    return 0;
}

cl_int compute_bottom()
{
	cl_int errNum;
	int width = 	imWidth;
	int height =  (imHeight/4);

	print_sub_2();

    errNum = clSetKernelArg(kernel, 0, sizeof(cl_mem), &sub_buffer_in[2]);
	errNum |= clSetKernelArg(kernel, 1, sizeof(int), &imWidth);
	errNum |= clSetKernelArg(kernel, 2, sizeof(int), &imHeight);
    errNum |= clSetKernelArg(kernel, 3, sizeof(int) * ThreadsX * ThreadsY, NULL);
	errNum |= clSetKernelArg(kernel, 4, sizeof(cl_mem), &temp_zones[1]);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error setting kernel args." << std::endl;
    }

	int globalX = ((width - 1)/(ThreadsX - 2) + 1) * ThreadsX;
	int globalY = ((height - 1)/(ThreadsY - 2) + 1) * ThreadsY;

if(show_debug)
{
	std::cout << "Compute Bottom Work Sizes" << std::endl; 
	std::cout << "Global size " << globalX << ", " << globalY << std::endl;
	std::cout << "Local size " << ThreadsX << ", " << ThreadsY << std::endl;
	std::cout << std::endl;
}

    size_t localWorkSize[2] = { ThreadsX, ThreadsY } ;
    size_t globalWorkSize[2] =  {  globalX, globalY };

	cl_event prof_event;
    errNum = clEnqueueNDRangeKernel(commandQueue, kernel, 2, NULL,
                                    globalWorkSize, localWorkSize,
                                    0, NULL, &prof_event);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error queuing conway kernel bottom for execution." << std::endl;
    }

	//std::cout << "\tcompute - " << completeEvent(prof_event, commandQueue) << " ms" << std::endl;

	// read required part of buffer
	int* output = new int[width];
	cl_event read_event;
	errNum = clEnqueueReadBuffer(commandQueue, 
		temp_zones[1],
		CL_FALSE,
		sizeof(int) * width,
		sizeof(int) * width,
		output,
		0,
		NULL,
		&read_event);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error reading bottom output." << std::endl;
    }
	//std::cout << "\tread - " << completeEvent(prof_read, commandQueue) << " ms" << std::endl;

	errNum = clSetEventCallback(read_event, CL_COMPLETE, &transmit_bottom_ghost, output);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error setting event callback (read)." << std::endl;
    }
/*
if(show_debug)
{
	std::cout << "BOTTOM ghost zone to send" << std::endl; 
	// initialize each cell to a random state
    for( int y = 0; y < width; y++ )
	{
		std::cout << output[y] << " ";
	}
	std::cout << std::endl; 
}
*/
	clReleaseEvent(prof_event);
	//clReleaseEvent(prof_read);

	//delete [] output;

    return 0;
}

void update_buffer()
{
	cl_int errNum;
	int width = imWidth;
	int height = imHeight/4;

	// write required temp zone section to output buffer 
	errNum = clEnqueueCopyBuffer(commandQueue,
		 temp_zones[0],	// source buffer
		 conway_out, // destination buffer
		 sizeof(int) * width,	// source offset
		 sizeof(int) * width,		// destination offset
		 sizeof(int) * width * (height - 2),	// size to copy
		 0,
		 NULL,
		 NULL);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error copy temp zone to output buffer." << std::endl;
    }

	// write required temp zone section to output buffer 
	errNum = clEnqueueCopyBuffer(commandQueue,
		 temp_zones[1],	// source buffer
		 conway_out, // destination buffer
		 sizeof(int) * width,	// source offset
		 sizeof(int) * width * (3*height + 1),		// destination offset
		 sizeof(int) * width * (height - 2),	// size to copy
		 0,
		 NULL,
		 NULL);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error copy temp zone to output buffer." << std::endl;
    }
}

void create_subs()
{
	cl_int errNum;

	// Top chunk
	int rows_per_chunk = (imHeight/4);
	int chunk_size = sizeof(int) * imWidth * rows_per_chunk;

	cl_buffer_region region0 =
	{
		0,
		chunk_size
	};
	sub_buffer_in[0] = clCreateSubBuffer(conway_in,
		CL_MEM_READ_ONLY,
		CL_BUFFER_CREATE_TYPE_REGION,
		&region0,
		&errNum);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error creating sub buffer" << std::endl;
    }

	// create ghost zone for incoming ghost zone (1 row)
	cl_buffer_region ghost_top =
	{
		0,
		sizeof(int) * imWidth
	};
	ghost_zones[0] = clCreateSubBuffer(conway_out,
		CL_MEM_READ_WRITE,
		CL_BUFFER_CREATE_TYPE_REGION,
		&ghost_top,
		&errNum);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error creating sub buffer" << std::endl;
    }
	check_zones[0] = clCreateSubBuffer(conway_check,
		CL_MEM_READ_WRITE,
		CL_BUFFER_CREATE_TYPE_REGION,
		&ghost_top,
		&errNum);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error creating sub buffer" << std::endl;
    }

	int width = 	imWidth;
	int height =  (imHeight/4);

	print_sub_0();

		// Bottom chunk
	cl_buffer_region region2 =
	{
		3 * chunk_size,
		chunk_size
	};

	sub_buffer_in[2] = clCreateSubBuffer(conway_in,
		CL_MEM_READ_ONLY,
		CL_BUFFER_CREATE_TYPE_REGION,
		&region2,
		&errNum);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error creating sub buffer" << std::endl;
    }

	// create ghost zone for incoming ghost zone (1 row)
	cl_buffer_region ghost_bottom =
	{
		sizeof(int) * imWidth * (imHeight -1),
		sizeof(int) * imWidth
	};
	ghost_zones[1] = clCreateSubBuffer(conway_out,
		CL_MEM_READ_WRITE,
		CL_BUFFER_CREATE_TYPE_REGION,
		&ghost_bottom,
		&errNum);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error creating sub buffer" << std::endl;
    }
	check_zones[1] = clCreateSubBuffer(conway_check,
		CL_MEM_READ_WRITE,
		CL_BUFFER_CREATE_TYPE_REGION,
		&ghost_bottom,
		&errNum);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error creating sub buffer" << std::endl;
    }

	print_sub_2();

	// Middle chunk (2 center chunks with ghost zones for overlap)
	cl_buffer_region region1 =
	{
		sizeof(int) * imWidth * (rows_per_chunk - 2),
		2 * sizeof(int) * imWidth * (rows_per_chunk + 4)
	};

	sub_buffer_in[1] = clCreateSubBuffer(conway_in,
		CL_MEM_READ_ONLY,
		CL_BUFFER_CREATE_TYPE_REGION,
		&region1,
		&errNum);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error creating sub buffer" << std::endl;
    }
	sub_buffer_out = clCreateSubBuffer(conway_out,
		CL_MEM_READ_WRITE,
		CL_BUFFER_CREATE_TYPE_REGION,
		&region1,
		&errNum);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error creating sub buffer" << std::endl;
    }

	print_sub_1();
}

void print_output()
{
	cl_int errNum;
	int* output = new int[imWidth*imHeight];
	errNum = clEnqueueReadBuffer(commandQueue, 
		conway_out,
		CL_TRUE,
		0,
		sizeof(int) * imWidth * imHeight,
		output,
		0,
		NULL,
		NULL);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error reading output buffer." << std::endl;
    }
	clFinish(commandQueue);
	std::cout << "Displaying Output Buffer.." << std::endl; 
	 // initialize each cell to a random state
    for( int x = 0; x < imWidth; x++ )
    {
        for( int y = 0; y < imHeight; y++ )
        {
			std::cout << output[y + x * imWidth] << " ";		// instead of to screen we could compare with pConwayBuffer directly and output any errors
		}
		std::cout << std::endl; 
    }
	std::cout << std::endl; 
	delete [] output;
}

void print_input()
{
	cl_int errNum;
	int* input = new int[imWidth*imHeight];
	errNum = clEnqueueReadBuffer(commandQueue, 
		conway_in,
		CL_TRUE,
		0,
		sizeof(int) * imWidth * imHeight,
		input,
		0,
		NULL,
		NULL);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error reading input buffer." << std::endl;
    }
	clFinish(commandQueue);
	std::cout << "Displaying Input Buffer.." << std::endl; 
	 // initialize each cell to a random state
    for( int x = 0; x < imWidth; x++ )
    {
        for( int y = 0; y < imHeight; y++ )
        {
			std::cout << input[y + x * imWidth] << " ";		// instead of to screen we could compare with pConwayBuffer directly and output any errors
		}
		std::cout << std::endl; 
    }
	std::cout << std::endl; 
	delete [] input;
}

void check_compute()
{
	cl_int errNum;
	int width = imWidth;
	int height = imHeight;

	// read output buffer
	int* output = new int[width*height];
	//cl_event prof_read;
	errNum = clEnqueueReadBuffer(commandQueue, 
		conway_out,
		CL_TRUE,
		0,
		sizeof(int) * width * height,
		output,
		0,
		NULL,
		NULL);
		//&prof_read);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error reading output." << std::endl;
    }

	// read check buffer
	int* check = new int[width*height];
	//cl_event prof_read;
	errNum = clEnqueueReadBuffer(commandQueue, 
		conway_check,
		CL_TRUE,
		0,
		sizeof(int) * width * height,
		check,
		0,
		NULL,
		NULL);
		//&prof_read);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error reading output." << std::endl;
    }
	
	// check for errors
	int errors = 0;
    for( int x = 0; x < height; x++ )
    {
        for( int y = 0; y < width; y++ )
        {
			if(output[y + x * width] != check[y + x * width])
			{
				std::cout << "error location found - " << x << ", " << y << std::endl;
				errors++;
			}
		}
    }
	// if errors output both to screen
	if(errors > 0)
	{
		std::cout << errors << " ERRORS FOUND!" << std::endl; 
		std::cout << "PRINTING OUTPUT" << std::endl; 
		 // initialize each cell to a random state
		for( int x = 0; x < imHeight; x++ )
		{
			for( int y = 0; y < imWidth; y++ )
			{
				std::cout << output[y + x * imWidth] << " ";
			}
			std::cout << std::endl; 
		}
		std::cout << std::endl; 

		std::cout << "PRINTING CHECK" << std::endl; 
		 // initialize each cell to a random state
		for( int x = 0; x < imHeight; x++ )
		{
			for( int y = 0; y < imWidth; y++ )
			{
				std::cout << check[y + x * imWidth] << " ";
			}
			std::cout << std::endl; 
		}
		std::cout << std::endl; 
		std::cout << "CHECK ERRORS!" << std::endl; 
		system("pause");
	}
	else
		std::cout << "No errors found " << std::endl;

	delete [] output;
	delete [] check;
}

cl_int swap_io()
{
	cl_int errNum;
	int width = imWidth;
	int height = imHeight;

	// set kernel arguments to compute,
	// sub buffers require no need for offset
    errNum = clSetKernelArg(copy_kernel, 0, sizeof(cl_mem), &conway_in);
	errNum |= clSetKernelArg(copy_kernel, 1, sizeof(cl_mem), &conway_out);

    size_t localWorkSize[2] = { ThreadsX, ThreadsY } ;
    size_t globalWorkSize[2] =  {  width, height };

if(show_debug)
{
	std::cout << "Swap IO Work Sizes" << std::endl; 
	std::cout << "Global size " << width << ", " << height << std::endl;
	std::cout << "Local size " << ThreadsX << ", " << ThreadsY << std::endl;
	std::cout << std::endl;
}

	//cl_event prof_event;
	// compute
    errNum = clEnqueueNDRangeKernel(commandQueue, copy_kernel, 2, NULL,
                                    globalWorkSize, localWorkSize,
                                    0, NULL, NULL);//&prof_event);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error queuing conway kernel top for execution." << std::endl;
    }

	//std::cout << "\tcompute - " << completeEvent(prof_event, commandQueue) << " ms" << std::endl;

	//clReleaseEvent(prof_event);

    return 0;
}

///
//      main() 
//
int main(int argc, char** argv)
{
	std::cout<<"starting application .. \n";

	cl_device_id device = 0;

	srand( (UINT)time( NULL ) );


	//imWidth = RoundUp(32, 400); 
	//imHeight = RoundUp(32, 300);

	// Create an OpenCL context
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
    program = CreateProgram(context, device, "life.cl");
    if (program == NULL)
    {
        Cleanup();
        return 1;
    }


	// create cl input/output buffers
	CreateIOResources();
   
	create_subs();

	// Create OpenCL kernel
    kernel = clCreateKernel(program, "conway_kernel", NULL);
    if (kernel == NULL)
    {
        std::cerr << "Failed to create conway kernel" << std::endl;
        Cleanup();
        return 1;
    }

	// Create OpenCL kernel
    copy_kernel = clCreateKernel(program, "copy_kernel", NULL);
    if (copy_kernel == NULL)
    {
        std::cerr << "Failed to create copy kernel" << std::endl;
        Cleanup();
        return 1;
    }

	// insert simulation loop here
	for(int j = 0; j < 100; j++)
	{
		// compute pass to calculate correct output buffer
		// conway_check will not have updated ghost zones
		compute_check();

		compute_top();
		receive_top();

		compute_bottom();
		receive_bottom();

		compute_middle();	

		update_buffer();

		check_compute();

if(show_debug)
{
		print_output();
}

	//	std::cout << "input before" << std::endl; 
	//	print_input();
	//	std::cout << "swapping input and output" << std::endl; 
	//	std::swap<cl_mem>(conway_in, conway_out);
	//	std::swap<cl_mem>(sub_buffer_in[1], sub_buffer_out);
	swap_io();
	//	std::cout << "input after" << std::endl; 
	//	print_input();

	//	std::cout << "sub buffers" << std::endl; 
	//	print_sub_0();
	//	print_sub_1();
	//	print_sub_2();

		system("pause");
	}

    std::cout << std::endl;
    std::cout << "Executed program succesfully." << std::endl;
    Cleanup();

	//system("pause");
    return 0;
}
