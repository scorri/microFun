

#ifndef VERTEX_H
#define VERTEX_H

struct Vertex
{
	Vertex(){}
	Vertex(float x, float y, float z, 
		float u, float v, float w)
		: pos(x,y,z),texC(u,v,w){}

	D3DXVECTOR3 pos;
	D3DXVECTOR3 texC;
};

#endif // VERTEX_H

