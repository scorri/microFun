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

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

const int width = 512;
const int height = 512;
const int ThreadsX = 16;
const int ThreadsY = 16;
const int iterations = 1000;

struct float2
{
	float x;
	float y;
};

struct profile
{
	int threads;
	double median;
	double mean;
	double stddev;
	double band;
	profile(int _t, double _med, double _mean, double sd)
	{
		threads = _t; median = _med, mean = _mean; stddev = sd; band = 0.0;
	}
	void print()
	{
		// calculate GB/sec
		band = width*height*sizeof(float2) * 1e-6 / (median);
		printf("Benchmark Summary\n");
		//printf("Number of elements: %d\n", width*height);
		printf("ThreadsX: %d\n", threads);
		printf("Median: %.2f ms\n", median);	
		printf("Mean: %.2f ms\n", mean);
		printf("Standard Deviation: %.2f\n", stddev);
		printf("Effective Bandwidth: %.2f GB/sec\n", band); 
	}
};

std::vector<float> results;
std::vector<profile> records;

///
// OpenCL variables, objects
//
cl_kernel kernel;
cl_mem rd_in;
cl_mem rd_results;
cl_context context;
cl_command_queue commandQueue;
cl_program program;

bool data_out = false;

///
// Forward declarations
void Cleanup();

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
	std::cout << "Results" << std::endl;
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

cl_int computeSum(cl_kernel k, bool record, int div)
{
	cl_int errNum;
	
    errNum = clSetKernelArg(k, 0, sizeof(cl_mem), &rd_in);
    errNum |= clSetKernelArg(k, 1, sizeof(float2) * ThreadsX * ThreadsY, NULL);
	errNum |= clSetKernelArg(k, 2, sizeof(cl_mem), &rd_results);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error setting kernel arguments." << std::endl;
    }

    size_t localWorkSize[1] = { ThreadsX*ThreadsY } ;
    size_t globalWorkSize[1] =  {  width*height/div };

	cl_event prof_event;
    errNum = clEnqueueNDRangeKernel(commandQueue, k, 1, NULL,
                                    globalWorkSize, localWorkSize,
                                    0, NULL, &prof_event);
    if (errNum != CL_SUCCESS)
    {
		std::cerr << errNum << std::endl;
        std::cerr << "Error queuing conway kernel for execution." << std::endl;
    }

	if(record)
		results.push_back( completeEvent(prof_event, commandQueue) );
	else
		checkSum(rd_results, width/ThreadsX, height/ThreadsY/div);
    
	return 0;
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

        if( rd_in != 0 )
                clReleaseMemObject(rd_in);

		if( rd_results != 0 )
                clReleaseMemObject(rd_results);

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

cl_context CreateContext()
{
	cl_int errNum;
	cl_uint numPlatforms;
	cl_platform_id firstPlatformId;
	cl_context context = NULL;

	errNum = clGetPlatformIDs(1, &firstPlatformId, &numPlatforms);
	if(errNum != CL_SUCCESS || numPlatforms <=0)
	{
		std::cerr << "Failed to find any OpenCl platforms." << std::endl;
		return NULL;
	}

	cl_context_properties contextProperties[] =
	{
		CL_CONTEXT_PLATFORM,
		(cl_context_properties)firstPlatformId,
		0
	};
	context = clCreateContextFromType(contextProperties, 
									CL_DEVICE_TYPE_GPU,
									NULL, NULL, &errNum);

	if(errNum != CL_SUCCESS)
	{
		std::cout << "Could not create GPU context, trying CPU.." << std::endl;
		context = clCreateContextFromType(contextProperties,
									CL_DEVICE_TYPE_CPU,
									NULL, NULL, &errNum);
		if(errNum != CL_SUCCESS)
		{
			std::cerr << "Failed to create an OpenCL GPU or CPU context.";
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

    delete [] pGSBuffer;

	// create buffer to store preliminary sum results
	// each work group will produce one results value
	rd_results = clCreateBuffer(
		context,
		CL_MEM_READ_WRITE,
		sizeof(float2)*(width/ThreadsX)*(height/ThreadsY),
		NULL,
		&errNum);
	if( errNum != CL_SUCCESS )
	{
		std::cerr<< "Failed creating results buffer" << std::endl;
	}
}

void recordStats(std::vector<float>& results)
{
	// Median	
	std::sort( results.begin(), results.end());
	double med = 0.0;
	if(results.size()/2 == 0)
		med = results[ results.size()/2 ];
	else
	{
		med = (results[ results.size()/2 ] + results[ results.size()/2 - 1])/2.0; 
	}

	// Mean
	double sum = std::accumulate(std::begin(results), std::end(results), 0.0);
	double m =  sum / results.size();

	// Standard deviation
	double accum = 0.0;
	std::for_each (std::begin(results), std::end(results), [&](const double d) {
		accum += (d - m) * (d - m);
	});
	double stdev = sqrt(accum / (results.size()-1));

	profile p = profile(ThreadsX*ThreadsY, med, m, stdev);
	p.print();

	// record stats
	records.push_back( p );
	
}

cl_int benchmarkReduce(const char* kernelName, int div)
{
	// Create OpenCL kernel
    kernel = clCreateKernel(program, kernelName, NULL);
    if (kernel == NULL)
    {
        std::cerr << "Failed to create rd kernel" << std::endl;
        Cleanup();
        return 1;
    }

	std::cout << "\nPre benchmarking results" << std::endl;
	computeSum(kernel, false, div);

	for(int i = 0; i < iterations; i++)
	{
		computeSum(kernel, true, div);
	}

	std::cout << "Post benchmarking results" << std::endl;
	computeSum(kernel, false, div);

	// release kernel
	clReleaseKernel(kernel);

	// clear results
	recordStats(results);
	results.clear();
}


///
//      main() for GLinterop example
//
int main(int argc, char** argv)
{
	printf("Reduction Benchmark OpenCL\n\n");

	// Create an OpenCL context with GL sharing 
    context = CreateContext();
    if (context == NULL)
    {
        std::cerr << "Failed to create OpenCL context." << std::endl;
        return 1;
    }

	cl_device_id device = 0;

    // Create a command-queue on the first device available
    // on the created context
    commandQueue = CreateCommandQueue(context, &device);
    if (commandQueue == NULL)
    {
        Cleanup();
        return 1;
    }

    // Create OpenCL program from GLinterop.cl kernel source
    program = CreateProgram(context, device, "reduce.cl");
    if (program == NULL)
    {
        Cleanup();
        return 1;
    }

	// create cl input/output buffers
	CreateIOResources();

	// CPU results
	checkSum(rd_in, width, height);

	benchmarkReduce("sum0", 1);
	benchmarkReduce("sum1", 1);
	benchmarkReduce("sum2", 2);
	benchmarkReduce("sum3", 2);

	std::cout << "Benchmarking complete!" << std::endl;
	Cleanup();

	return 0;
}