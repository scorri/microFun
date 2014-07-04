// definition of the display info constant buffer
// definition of the displaty info constant buffer
cbuffer updates	: register(b0)
{
	float4x4 wvp	: packoffset(c0);
    int r_ins		: packoffset(c4.x);
    int r_ni		: packoffset(c4.y);
    int r_mb		: packoffset(c4.z) ;
    int r_all_bio	: packoffset(c4.w);
    int r_re_ex		: packoffset(c5.x);
	int geometry	: packoffset(c5.y);

};