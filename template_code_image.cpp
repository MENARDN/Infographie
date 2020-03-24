#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define _USE_MATH_DEFINES
#include "math.h"

#include <iostream>     // std::cout
#include <algorithm>    // std::max

#include <omp.h>

#include <random>

#include <map>
#include <list>

std::default_random_engine engine[8];
std::uniform_real_distribution<double> distrib(0, 1);

class Vector {
public:
	Vector(double x = 0, double y = 0, double z = 0) {
		coord[0] = x;
		coord[1] = y;
		coord[2] = z;
	};
	double coord[3];

	double operator[](int i)  const {
		return coord[i];
	};

	double& operator[](int i) {
		return coord[i];
	};

	double Norm2() const {
		return coord[0] * coord[0] + coord[1] * coord[1] + coord[2] * coord[2];
	}
	
	void normalize() {
		double N = sqrt(Norm2());
		coord[0] = coord[0] / N;
		coord[1] = coord[1] / N;
		coord[2] = coord[2] / N;
	}
};

double dot(const Vector& A, const Vector& B) {
	return A[0] * B[0] + A[1] * B[1] + A[2] * B[2];
}

Vector operator+(const Vector& A, const Vector& B) {
	return Vector(A[0] + B[0], A[1] + B[1], A[2] + B[2]);
}

Vector operator-(const Vector& A, const Vector& B) {
	return Vector(A[0] - B[0], A[1] - B[1], A[2] - B[2]);
}

Vector operator*(const double& a, const Vector& B) {
	return Vector(a*B[0], a*B[1], a*B[2]);
}

Vector operator*(const Vector& A, const double& b) {
	return Vector(A[0]*b, A[1]*b, A[2]*b);
}

Vector operator*(const Vector& A, const Vector& B) {
	return Vector(A[0] * B[0], A[1] * B[1], A[2] * B[2]);
}

Vector operator/(const Vector& A, const double& b) {
	return Vector(A[0] / b, A[1] / b, A[2] / b);
}
Vector cross(const Vector& A, const Vector& B) {
	return Vector(A[1] * B[2] - A[2] * B[1], A[2] * B[0] - A[0] * B[2], A[0] * B[1] - A[1] * B[0]);
}

class Ray {
public:
	Ray(const Vector& C, Vector u) : C(C), u(u) {};
	Vector C, u;
};


class Object {
public:
	Object() {}
	Object(bool diffuse, bool transp, Vector albedo, Vector Emi) : diffuse(diffuse), transp(transp), albedo(albedo), Emi(Emi) {};

	virtual bool intersect(const Ray& r, Vector& P, Vector& N, double& t) = 0;

	void set_emi(Vector& E) {
		Emi = E;
	}

	bool diffuse, transp;
	Vector albedo;
	Vector Emi = Vector(0, 0, 0);
};

class Sphere: public Object{
public:
	Sphere(const Vector& O, double R, const Vector& albedo, bool diffuse, bool transp, Vector Emi) :
		Object(diffuse,transp,albedo,Emi), R(R), O(O) {};


	bool intersect(const Ray& r, Vector &P, Vector &N, double &t) {
		double a = 1;
		double b = 2 * dot(r.u, r.C - O);
		double c = (r.C - O).Norm2() - R * R;

		double delta = b * b - 4 * a * c;
		if (delta < 0) return false;
		double sqrtDelta = sqrt(delta);
		double t1 = (-b + sqrtDelta) / (2 * a);
		double t0 = (-b - sqrtDelta) / (2 * a);
		if (t1 > t0) {
			t = t0;
		}
		else {
			t = t1;
		}
		P = r.C + t * r.u;
		N = P-O;
		N.normalize();
		return true;
	}

	Vector O;
	double R;
};

Vector random_cos(const Vector& N) {
	//orthogonal frame:
	Vector T1;
	if (std::abs(N[0]) <= std::abs(N[1]) && std::abs(N[0]) <= std::abs(N[2])) {
		T1 = Vector(0, -N[2], N[1]);
	}
	else if (std::abs(N[1]) <= std::abs(N[0]) && std::abs(N[1]) <= std::abs(N[2])) {
		T1 = Vector(-N[2],0,N[0]);
	}
	else {
		T1 = Vector(-N[1],N[0],0);
	}
	T1.normalize();
	Vector T2 = cross(N, T1);

	double r1 = distrib(engine[omp_get_thread_num()]);
	double r2 = distrib(engine[omp_get_thread_num()]);
	Vector V;
	V[0] = cos(2 * M_PI * r1) * sqrt(1 - r2);
	V[1] = sin(2 * M_PI * r1) * sqrt(1 - r2);
	V[2] = sqrt(r2);

	return V[0] * T1 + V[1] * T2 + V[2] * N;
}

class Bbox {
public:
	Bbox() {};
	bool intersect(const Ray& r, double& t) {
		bool inside = true;
		Vector pos, plane;
		for (int i = 0; i < 3; i++) {
			if (r.C[i] < m[i]) {
				inside = false;
				pos[i] = -1;
				plane[i] = m[i];
			}
			else if (r.C[i] > M[i])  {
				inside = false;
				pos[i] = 1;
				plane[i] = M[i];
			}
			else {
				pos[i] = 0;
			}
		}
		if (inside) {
			t = 0;
			return true;
		}

		Vector T;
		for (int i = 0; i < 3; i++) {
			if (pos[i] != 0 && r.u[i] != 0) {
				T[i] = (plane[i] - r.C[i]) / r.u[i];
			}
			else {
				T[i] = -1;
			}
		}

		int whichPlane = 0;
		for (int i = 1; i < 3; i++)
			if (T[whichPlane] < T[i])
				whichPlane = i;

		if (T[whichPlane] < 0.) return (false);

		Vector P;

		for (int i = 0; i < 3; i++)
			if (whichPlane != i) {
				P[i] = r.C[i] + T[whichPlane] * r.u[i];
				if (P[i] < m[i] || P[i] > M[i])
					return (false);
			}
			else {
				P[i] = plane[i];
			}
		t = T[whichPlane];
		return true;
	}

	Vector m, M;
};

class BVH {
public:
	BVH* fg, * fd;
	int i0, i1;
	Bbox b;
};

class TriangleIndices {
public:
	TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk) {
	};
	int vtxi, vtxj, vtxk;
	int uvi, uvj, uvk;
	int ni, nj, nk;
	int faceGroup;
};

class Triangle {
public:
	Triangle() {}
	Triangle(Vector a, Vector b, Vector c) {
		A = a;
		B = b;
		C = c;
	}
	bool intersect(const Ray& r, Vector& P, Vector& N, double& t, double& beta, double& gamma) {

		N = -1 * cross(B - A, C - A);
		N.normalize();
		double num = -1 * dot(r.C - A, N);
		double den = dot(r.u, N);
		if (den == 0) {
			return false;
		}
		t = num / den;
		if (t < 0) {
			return false;
		}
		P = r.C + t * r.u;
		Vector v0 = C - A;
		Vector v1 = B - A;
		Vector v2 = P - A;
		double dot00 = dot(v0, v0);
		double dot01 = dot(v0, v1);
		double dot02 = dot(v0, v2);
		double dot11 = dot(v1, v1);
		double dot12 = dot(v1, v2);
		double inv_denom = 1 / (dot00 * dot11 - dot01 * dot01);
		gamma = (dot11 * dot02 - dot01 * dot12) * inv_denom;
		beta = (dot00 * dot12 - dot01 * dot02) * inv_denom;
		return (beta >= 0 && gamma >= 0 && (beta + gamma) < 1);

	}

	Vector A, B, C;
};

class Mesh : public Object {
public:
	Mesh() {}
	Mesh(const char* obj, double scaling, const Vector& offset, bool diffuse, bool transp, Vector albedo, Vector Emi) : Object(diffuse, transp, albedo, Emi) {
		readOBJ(obj);
		for (int i = 0; i < vertices.size(); i++) {
			vertices[i] = vertices[i] * scaling + offset;
		}
		std::cout << indices.size() << std::endl;

		build_bvh(0, indices.size(),&bvh);


	}

	void add_texture(const char* filename) {

		textures.resize(textures.size() + 1);
		w.resize(w.size() + 1);
		h.resize(h.size() + 1);

		FILE* f;
		f = fopen(filename, "rb");
		unsigned char info[54];
		fread(info, sizeof(unsigned char), 54, f); // read the 54-byte header

		w[w.size() - 1] = *(int*)&info[18]; // extract image height and width from header
		h[h.size() - 1] = *(int*)&info[22];

		int size = 3 * w[w.size() - 1] * h[h.size() - 1];
		textures[textures.size() - 1].resize(size); // allocate 3 bytes per pixel
		fread(&textures[textures.size() - 1][0], sizeof(unsigned char), size, f); // read the rest of the data at once
		fclose(f);

		for (int i = 0; i < size; i += 3) {
			std::swap(textures[textures.size() - 1][i], textures[textures.size() - 1][i + 2]);
		}
	}

	void readOBJ(const char* obj) {

		char matfile[255];
		char grp[255];

		FILE* f;
		f = fopen(obj, "r");

		std::map<std::string, int> groupNames;
		int curGroup = -1;
		while (!feof(f)) {
			char line[255];
			if (!fgets(line, 255, f)) break;

			std::string linetrim(line);
			linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
			strcpy(line, linetrim.c_str());

			if (line[0] == 'u' && line[1] == 's') {
				sscanf(line, "usemtl %[^\n]\n", grp);
				if (groupNames.find(std::string(grp)) != groupNames.end()) {
					curGroup = groupNames[std::string(grp)];
				}
				else {
					curGroup = groupNames.size();
					groupNames[std::string(grp)] = curGroup;
				}
			}
			if (line[0] == 'm' && line[1] == 't' && line[2] == 'l') {
				sscanf(line, "mtllib %[^\n]\n", matfile);
			}
			if (line[0] == 'v' && line[1] == ' ') {
				Vector vec;
				Vector col;
				if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[2], &vec[1], &col[0], &col[1], &col[2]) == 6) {
					vertices.push_back(vec);
					vertexcolors.push_back(col);
				}
				else {
					sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[2], &vec[1]);  // helmet
																				 //vec[2] = -vec[2]; //car2
					vertices.push_back(vec);
				}
			}
			if (line[0] == 'v' && line[1] == 'n') {
				Vector vec;
				sscanf_s(line, "vn %lf %lf %lf\n", &vec[0], &vec[2], &vec[1]); //girl
				normals.push_back(vec);
			}
			if (line[0] == 'v' && line[1] == 't') {
				Vector vec;
				sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
				uvs.push_back(vec);
			}
			if (line[0] == 'f') {
				TriangleIndices t;
				int i0, i1, i2, i3;
				int j0, j1, j2, j3;
				int k0, k1, k2, k3;
				int nn;

				char* consumedline = line + 1;
				int offset;
				t.faceGroup = curGroup;
				nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
				if (nn == 9) {
					if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
					if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
					if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
					if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
					if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
					if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
					if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
					if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
					if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;

					indices.push_back(t);
				}
				else {
					nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
					if (nn == 6) {
						if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
						if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
						if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
						if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
						if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
						if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
						indices.push_back(t);
					}
					else {
						nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
						if (nn == 3) {
							if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
							indices.push_back(t);
						}
						else {
							nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
							if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
							if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
							if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
							if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
							indices.push_back(t);
						}
					}
				}


				consumedline = consumedline + offset;

				while (true) {
					if (consumedline[0] == '\n') break;
					if (consumedline[0] == '\0') break;
					nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
					TriangleIndices t2;
					t2.faceGroup = curGroup;
					if (nn == 3) {
						if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
						if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
						if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
						if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
						if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
						if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
						if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
						if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
						if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;
						indices.push_back(t2);
						consumedline = consumedline + offset;
						i2 = i3;
						j2 = j3;
						k2 = k3;
					}
					else {
						nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
						if (nn == 2) {
							if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
							if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
							if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
							if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
							if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
							if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
							consumedline = consumedline + offset;
							i2 = i3;
							j2 = j3;
							indices.push_back(t2);
						}
						else {
							nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
							if (nn == 2) {
								if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
								if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
								if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
								if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
								if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
								if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;
								consumedline = consumedline + offset;
								i2 = i3;
								k2 = k3;
								indices.push_back(t2);
							}
							else {
								nn = sscanf(consumedline, "%u%n", &i3, &offset);
								if (nn == 1) {
									if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
									if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
									if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
									consumedline = consumedline + offset;
									i2 = i3;
									indices.push_back(t2);
								}
								else {
									consumedline = consumedline + 1;
								}
							}
						}
					}
				}

			}


		}
		fclose(f);
	}

	Bbox computeBbox(int i0, int i1) {
		Bbox bbox;
		bbox.m = Vector(1E9, 1E9, 1E9);
		bbox.M = Vector(-1E9, -1E9, -1E9);
		for (int i = i0; i < i1; i++) {
			TriangleIndices tid = indices[i];
			for (int j = 0; j < 3; j++) {
				bbox.m[j] = std::min(bbox.m[j], vertices[tid.vtxi][j]);
				bbox.M[j] = std::max(bbox.M[j], vertices[tid.vtxi][j]);
				bbox.m[j] = std::min(bbox.m[j], vertices[tid.vtxj][j]);
				bbox.M[j] = std::max(bbox.M[j], vertices[tid.vtxj][j]);
				bbox.m[j] = std::min(bbox.m[j], vertices[tid.vtxk][j]);
				bbox.M[j] = std::max(bbox.M[j], vertices[tid.vtxk][j]);
			}
		}
		return bbox;

	}

	void build_bvh(int i0, int i1, BVH* bvh) {
		bvh->i0 = i0;
		bvh->i1 = i1;
		bvh->b = computeBbox(i0, i1);
		bvh->fg = NULL;
		bvh->fd = NULL;


		Vector Diagbox = bvh->b.M - bvh->b.m;

		//longest axis
		int axis;
		if (Diagbox[0] >= Diagbox[1] && Diagbox[0] >= Diagbox[2]) {
			axis = 0;
		} else {
			if (Diagbox[1] >= Diagbox[0] && Diagbox[1] >= Diagbox[2]) {
				axis = 1;
			}
			else {
				axis = 2;
			}
		}

		double middleAxis = bvh->b.m[axis] + Diagbox[axis] * 0.5;
		int pivot = i0;

		for (int i = i0; i < i1; i++) {
			double centerTriangleAxis = (vertices[indices[i].vtxi][axis] + vertices[indices[i].vtxj][axis] + vertices[indices[i].vtxk][axis]) / 3.;
			if (centerTriangleAxis < middleAxis) {
				std::swap(indices[pivot], indices[i]);
				pivot++;
			}
		}

		if (i1 - i0 > 3 && pivot > i0 + 1 && pivot < i1) {
			bvh->fg = new BVH;
			bvh->fd = new BVH;
			build_bvh(i0, pivot, bvh->fg);
			build_bvh(pivot, i1, bvh->fd);
		}
	}

	bool intersect(const Ray& r, Vector& P, Vector& N, double& t) {
		
		double tbox;
		if (!bvh.b.intersect(r, tbox)) {
			return false;
		}

		std::list<BVH*> nodes;
		nodes.push_back(&bvh);

		bool has_inter = false;
		t = 1E9;

		while (!nodes.empty()) {
			BVH* curNode = nodes.front();
			nodes.pop_front();
			if (curNode->fg) {
				double tleft, tright;
				if (curNode->fg->b.intersect(r, tleft)) {
					if (tleft < t) {
						nodes.push_front(curNode->fg);
					}
				}
				if (curNode->fd->b.intersect(r, tright)) {
					if (tright < t) {
						nodes.push_front(curNode->fd);
					}
				}
			}
			else {
				for (int i = curNode->i0; i < curNode->i1; i++) {
					Triangle tri(vertices[indices[i].vtxi], vertices[indices[i].vtxj], vertices[indices[i].vtxk]);
					Vector localP, localN;
					double localt, beta, gamma;
					if (tri.intersect(r, localP, localN, localt, beta, gamma)) {
						has_inter = true;
						if (localt < t) {
							t = localt;
							P = localP;
							double alpha = 1 - beta - gamma;
							N = alpha * normals[indices[i].ni] + beta * normals[indices[i].nj] + gamma * normals[indices[i].nk];
							//N = localN;
						}
					}
				}
			}
		}
		return has_inter;
	}

	BVH bvh;
	std::vector<TriangleIndices> indices;
	std::vector<Vector> vertices;
	std::vector<Vector> normals;
	std::vector<Vector> uvs;
	std::vector<Vector> vertexcolors;

	std::vector<std::vector<unsigned char> > textures;
	std::vector<int> w, h;


};

class Scene {
public:
	Scene() {}
	std::vector<Object*> Objects;
	Vector Light_O, Light_Emi;
	double Light_R;

	void add_objet(Object* S) {
		Objects.push_back(S);
	}

	void set_light(Sphere* S) {
		Objects.push_back(S);
		Light_O = S->O;
		Light_R = S->R;
		Light_Emi = S->Emi;
	}

	bool intersect(Ray r, Vector& P, Vector& N, Vector& albedo, bool& dif, bool& transp) {
		double tmin = -1;
		bool intersected = false;
		bool onscreen = false;
		for (int unsigned i = 0; i < Objects.size(); i++)
		{
			double t;
			Vector localP, localN;
			if (Objects[i]->intersect(r,localP,localN,t)) {
				if ((t < tmin || onscreen == false) && t > 0) {
					tmin = t;
					onscreen = true;
					P = localP;
					N = localN;
					albedo = Objects[i]->albedo;
					dif = Objects[i]->diffuse;
					transp = Objects[i]->transp;
				}

				intersected = true;
			}
		}
		return intersected;
	}

	Vector getColor(Ray r, int num_bounce) {
		Vector P, N, albedo, Emi;
		bool dif = true, transp = true;
		bool has_inter = intersect(r, P, N, albedo, dif, transp);
		if (has_inter) {
			if (dif) {

				Vector PL = P - Light_O;
				PL.normalize();
				Vector w = random_cos(PL);
				w.normalize();
				Vector xi = Light_O + Light_R * w;
				Vector wi = xi - P;
				double d2 = wi.Norm2();
				wi.normalize();
				double costheta = std::max(0., dot(N, wi));
				double costhetaprime = std::max(0., dot(w, -1 * wi));
				double costhetasecond = std::max(0., dot(PL, w));
				Ray r2(P + 0.0001*N, wi);
				Vector P2, N2, albedo2;
				bool has_obs = intersect(r2, P2, N2, albedo2, dif, transp);
				bool lighted = true;
				if (has_obs) {
					double obsdist2 = (P2 - P).Norm2();
					if (obsdist2 < d2-0.001) {
						lighted = false;
					}
				}
				Vector I;
				if (lighted) {
					double p = costhetasecond / (M_PI * 4 * M_PI * Light_R * Light_R);
					I = ((Light_Emi * albedo / M_PI * costheta * costhetaprime) / (d2 * p));
				}
				else {
					I = Vector(0., 0., 0.);
				}
				if (num_bounce != 0) {
					Vector wj = random_cos(N);
					wj.normalize();
					Vector indirect = getColor(Ray(P + 0.0001 * N, wj), num_bounce - 1);
					I = I + indirect * (albedo / M_PI);
				}
				return I;

			}
			else if (transp && num_bounce != 0){
				double n1 = 1;
				double n2 = 1.7;
				Vector N_trans(N);
				if (dot(r.u,N)>0) {
					n1 = 1.7;
					n2 = 1;
					N_trans = -1*N;
				}
				double inside_sqrt = 1 - pow(n1 / n2,2) * (1 - pow(dot(N_trans,r.u),2));
				if (inside_sqrt > 0) {
					double k0 = pow((n1 - n2), 2) / pow((n1 + n2), 2);
					double R;
					Vector newu = (n1 / n2) * (r.u - dot(r.u, N_trans) * N_trans) - N_trans * sqrt(inside_sqrt);
					newu.normalize();
					if (n1 < n2) {
						R = k0 + (1 - k0) * pow(1 - dot(N_trans, -1*r.u), 5);
					}
					else {
						R = k0 + (1 - k0) * pow(1 - dot(N_trans, newu), 5);
					}
					double r1 = distrib(engine[omp_get_thread_num()]);
					if (r1 > R) {
						Ray newr(P - 0.0001 * N_trans, newu);
						return getColor(newr, num_bounce - 1);
					}
					else {
						newu = r.u - 2 * dot(N_trans, r.u) * N_trans;
						newu.normalize();
						Ray newr(P + 0.0001 * N_trans, newu);
						return getColor(newr, num_bounce - 1);
					}

				}
				else {
					Vector newu = r.u - 2 * dot(N_trans, r.u) * N_trans;
					newu.normalize();
					Ray newr(P + 0.0001 * N_trans, newu);
					return getColor(newr, num_bounce - 1);
				}

			}
			else if (num_bounce != 0) {
				Vector newu = r.u - 2 * dot(r.u, N) * N;
				newu.normalize();
				Ray newr(P + 0.0001 * N, newu);
				return getColor(newr, num_bounce - 1);
			}
			else {
				return(Vector(0., 0., 0.));
			}
		}
		else {
			return(Vector(0., 0., 0.));
		}
	}



};

int main() {
	int W = 512;
	int H = 512;
	int N_rays = 500;

	double fov = 60 * M_PI / 180;
	Vector C(0, 0, 55);

	double distance_mise_au_point = 55;

	Vector L(-10, 20, 40);
	double l = 4e7;
	Scene mainscene;

	//Sphere(const Vector& O, double R, const Vector& albedo, bool diffuse, bool transp, Vector Emi)

	//Lightsource
	Sphere Light = Sphere(L, 2, Vector(1., 1., 1.), false, false, (l / (M_PI * 2 * 2)) * Vector(1, 1, 1));
	mainscene.set_light(&Light);

	

	//Spheres
	Sphere Sp1 = Sphere(Vector(10., 0., 10.), 10, Vector(1., 1., 1.), true, false, Vector(0, 0, 0));
	Sphere Sp2 = Sphere(Vector(-10., 0., -10.), 10, Vector(1., 1., 1.), true, false, Vector(0, 0, 0));
	Sphere Sp3 = Sphere(Vector(0., 30., 0.), 12, Vector(1., 1., 1.), false, false, Vector(0, 0, 0));
	Sphere Sp4 = Sphere(Vector(20, 20., 10.), 5, Vector(1., 1., 1.), false, true, Vector(0, 0, 0));
	//mainscene.add_objet(&Sp1);
	//mainscene.add_objet(&Sp2);
	//mainscene.add_objet(&Sp3);
	//mainscene.add_objet(&Sp4);

	

	//Walls
	Sphere W1 = Sphere(Vector(0., 0, -1000), 940, Vector(1., 0., 0.), true, false, Vector(0, 0, 0));
	Sphere W2 = Sphere(Vector(0., 1000, 0), 940, Vector(0., 0., 1.), true, false, Vector(0, 0, 0));
	Sphere W3 = Sphere(Vector(0., -1020, 0), 990, Vector(0., 1., 0.), true, false, Vector(0, 0, 0));
	Sphere W4 = Sphere(Vector(0., 0., 1000), 940, Vector(1., 1., 0.), true, false, Vector(0, 0, 0));
	Sphere W5 = Sphere(Vector(990, 0., 0), 940, Vector(1., 0., 1.), true, false, Vector(0, 0, 0));
	Sphere W6 = Sphere(Vector(-990, 0., 0), 940, Vector(0., 1., 1.), true, false, Vector(0, 0, 0));
	mainscene.add_objet(&W1);
	mainscene.add_objet(&W2);
	mainscene.add_objet(&W3);
	mainscene.add_objet(&W4);
	mainscene.add_objet(&W5);
	mainscene.add_objet(&W6);

	Mesh Girl = Mesh("girl.obj", 25, Vector(0, -15, 10), true, false, Vector(1, 1, 1), Vector(0, 0, 0));
	mainscene.add_objet(&Girl);

	std::vector<unsigned char> image(W*H * 3, 0);

	int num_core = 8;
	#pragma omp parallel for
	for (int i = 0; i < H; i++) {
			std::cout << "Thread " << omp_get_thread_num() << ": " << (100*i*(num_core) / H) - ((omp_get_thread_num())*100) << "%" << std::endl;
		for (int j = 0; j < W; j++) {
			double d = W / (2 * tan(fov / 2));
			
			Vector I(0, 0, 0);
			for (int k = 0; k < N_rays; k++) {
				double r1 = distrib(engine[omp_get_thread_num()]);
				double r2 = distrib(engine[omp_get_thread_num()]);
				double R = sqrt(-2 * log(r1));
				double x1 = R * cos(2 * M_PI * r2) * 0.25;
				double x2 = R * sin(2 * M_PI * r2) * 0.25;
				Vector u((j + x1 - 0.5 - W / 2 ),( -i + x2 - 0.5 + H / 2 ), -d);
				u.normalize();

				double r1b = distrib(engine[omp_get_thread_num()]);
				double r2b = distrib(engine[omp_get_thread_num()]);
				double x1b = R * cos(2 * M_PI * r2) * 0.25;
				double x2b = R * sin(2 * M_PI * r2) * 0.25;
				Vector Cprime = C + Vector(x1b, x2b, 0);
				Vector uprime = C + distance_mise_au_point * u - Cprime;
				uprime.normalize();

				Ray rini(Cprime, uprime);
				I = I + mainscene.getColor(rini, 1);
			}
			I = I / N_rays;
			image[(i * W + j) * 3 + 0] = std::min(255, std::max(0, int(pow(I[0], 0.45))));
			image[(i * W + j) * 3 + 1] = std::min(255, std::max(0, int(pow(I[1], 0.45))));
			image[(i * W + j) * 3 + 2] = std::min(255, std::max(0, int(pow(I[2], 0.45))));
		}
	}
	stbi_write_png("image.png", W, H, 3, &image[0], 0);

	return 0;
}