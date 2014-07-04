#ifndef __MYCO_CONSTANTS_H__
#define __MYCO_CONSTANTS_H__

// Window dimensions and depth of volume to simulate
#define WINDOW_SIZE_X 500
#define WINDOW_SIZE_Y 500 

//compute buffer dimenions
#define BUFFER_SIZE_X 200
#define BUFFER_SIZE_Y 200
#define BUFFER_SIZE_Z 200
	
// No. of threads to dispatch peer group - NOTE: Max. thread count for a group is 1024
#define THREADS_X 16      
#define THREADS_Y 8       
#define THREADS_Z 8

//// No, of thread groups to dispatch
//#define GROUPS_X  (BUFFER_SIZE_X + THREADS_X-1)/(THREADS_X-1)
//#define GROUPS_Y  (BUFFER_SIZE_Y + THREADS_Y-1)/(THREADS_Y-1)
//#define GROUPS_Z  (BUFFER_SIZE_Z + THREADS_Z-1)/(THREADS_Z-1)
// No, of thread groups to dispatch
#define GROUPS_X  ((BUFFER_SIZE_X-1)/(THREADS_X-2)+1)
#define GROUPS_Y  ((BUFFER_SIZE_Y-1)/(THREADS_Y-2)+1)
#define GROUPS_Z  ((BUFFER_SIZE_Z-1)/(THREADS_Z-2)+1)

typedef double scalar;

#endif