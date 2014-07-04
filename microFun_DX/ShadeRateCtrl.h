#ifndef SHADE_RATE_CTRL_H
#define SHADE_RATE_CTRL_H

#define CSRC_ID_COMPUTE (1)
#define CSRC_ID_DRAW	(0)
#define MAX_SHADERS (2)

#ifndef uint
typedef unsigned int uint;
#endif // uint

struct DelayCounts
{
	int count[MAX_SHADERS];
};

class CShadeRateCtrl
{
public :
	CShadeRateCtrl (void) {;}
	CShadeRateCtrl (uint compute, uint draw)
	{
		refDelay.count[CSRC_ID_COMPUTE]= compute;//max(1,compute);
		refDelay.count[CSRC_ID_DRAW]=	draw;//max(1,draw);

		activeDelay= refDelay;
		exec.count[CSRC_ID_COMPUTE]= exec.count[CSRC_ID_DRAW]= 0;
		nTick= exec;
	} // CShadeRateCtrl

	uint tick (void)
	{
		uint t= 0, i= MAX_SHADERS;
		while (i-- > 0)
		{
			if (--(activeDelay.count[i]) <= 0)
			{
				activeDelay.count[i]= refDelay.count[i];	// reset delay
				exec.count[i]++;	// add number of executions
				t+= exec.count[i];	// total
			}
		}
		return(t);
	} // tick

	// test whether ready to execute
	bool ready (const uint shaderIdx) const
	{
		if (shaderIdx < MAX_SHADERS)
		{
			//printf("exec[%u]=%u\n", shaderIdx, exec.count[shaderIdx]);
			return(exec.count[shaderIdx] <= 0);
		}
		return(false);
	} // ready

	// signal execution complete, testing whether more executions required
	bool complete (const uint shaderIdx)
	{
		if (shaderIdx < MAX_SHADERS)
		{
			++(nTick.count[shaderIdx]);
			return(--(exec.count[shaderIdx]) <= 0);
		}
		return(false);
	} // ready

	uint getTick (const uint shaderIdx)
	{
		if (shaderIdx < MAX_SHADERS)
		{
			return(nTick.count[shaderIdx]);
		}
		return(0);
	} // getTick

protected :
	DelayCounts refDelay;	// reference value used to reinstate active counts
	DelayCounts activeDelay;	// active counts, used to trigger execution
	DelayCounts exec;	// execution count (normally 1 or zero)
	DelayCounts nTick;
}; // CShadeRateCtrl

#endif // SHADE_RATE_CTRL_H