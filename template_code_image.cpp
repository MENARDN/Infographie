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

class Sphere {
public:
	Sphere(const Vector& O, double R, const Vector& albedo, bool diffuse, bool transp, Vector Emi) : 
		O(O), R(R), albedo(albedo), diffuse(diffuse), transp(transp), Emi(Emi) {};



	bool intersect(const Ray& r, Vector &P, Vector &N) {
		double a = 1;
		double b = 2 * dot(r.u, r.C - O);
		double c = (r.C - 0).Norm2() - R * R;

		double delta = b * b - 4 * a * c;
		if (delta < 0) return false;
		double sqrtDelta = sqrt(delta);
		double t1 = (-b + sqrtDelta) / (2 * a);
		double t0 = (-b - sqrtDelta) / (2 * a);
		double t;
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

	void set_emi(Vector& E) {
		Emi = E;
	}

	Vector O;
	Vector albedo;
	double R;
	bool diffuse, transp;
	Vector Emi = Vector(0,0,0);
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


class Scene {
public:
	Scene() {}
	std::vector<Sphere> Spheres;

	void add_sphere(const Sphere& S) {
		Spheres.push_back(S);
	}

	bool intersect(Ray r, Vector& P, Vector& N, Vector& albedo, bool& dif, bool& transp) {
		double tmin = -1;
		bool intersected = false;
		bool onscreen = false;
		for (int unsigned i = 0; i < Spheres.size(); i++)
		{
			double a = 1;
			double b = 2 * dot(r.u, r.C - Spheres[i].O);
			double c = (r.C - Spheres[i].O).Norm2() - Spheres[i].R * Spheres[i].R;

			double delta = b * b - 4 * a * c;
			if (delta > 0) {
				double sqrtDelta = sqrt(delta);
				double t1 = (-b + sqrtDelta) / (2 * a);
				double t0 = (-b - sqrtDelta) / (2 * a);
				double t;
				if (t1 > t0) {
					t = t0;
				}
				else {
					t = t1;
				}
				if ((t < tmin || onscreen == false) && t > 0) {
					tmin = t;
					onscreen = true;
					P = r.C + t * r.u;
					N = P - Spheres[i].O;
					N.normalize();
					albedo = Spheres[i].albedo;
					dif = Spheres[i].diffuse;
					transp = Spheres[i].transp;
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
				Vector PL = P - Spheres[0].O;
				PL.normalize();
				Vector w = random_cos(PL);
				w.normalize();
				Vector xi = Spheres[0].O + Spheres[0].R * w;
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
					double p = costhetasecond / (M_PI * 4 * M_PI * Spheres[0].R * Spheres[0].R);
					I = ((Spheres[0].Emi * albedo / M_PI * costheta * costhetaprime) / (d2 * p));
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
					Vector newu = (n1 / n2) * (r.u - dot(r.u,N_trans)*N_trans) - N_trans * sqrt(inside_sqrt);
					newu.normalize();
					Ray newr(P - 0.0001 * N_trans, newu);
					return getColor(newr, num_bounce - 1);
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
	int N_rays = 100;

	double fov = 60 * M_PI / 180;
	Vector C(0, 0, 55);

	double distance_mise_au_point = 55;

	Vector L(-10, 20, 40);
	double l = 4e7;
	Scene mainscene;

	//Lightsource
	mainscene.add_sphere(Sphere(L, 2, Vector(1., 1., 1.), false, false,(l/(M_PI*2*2))*Vector(1,1,1)));

	//Spheres
	mainscene.add_sphere(Sphere(Vector(10., 0., 10.), 10, Vector(1., 1., 1.), true, false, Vector(0, 0, 0)));
	mainscene.add_sphere(Sphere(Vector(-10., 0., -10.), 10, Vector(1., 1., 1.), true, false, Vector(0, 0, 0)));
	mainscene.add_sphere(Sphere(Vector(0., 30., 0.), 12, Vector(1., 1., 1.), false, false, Vector(0, 0, 0)));
	mainscene.add_sphere(Sphere(Vector(20, 20., 10.), 5, Vector(1., 1., 1.), false, true, Vector(0, 0, 0)));


	//Walls
	mainscene.add_sphere(Sphere(Vector(0., 0, -1000), 940, Vector(1., 0., 0.), true, false, Vector(0, 0, 0)));
	mainscene.add_sphere(Sphere(Vector(0., 1000, 0), 940, Vector(0., 0., 1.), true, false, Vector(0, 0, 0)));
	mainscene.add_sphere(Sphere(Vector(0., -1000, 0), 990, Vector(0., 1., 0.), true, false, Vector(0, 0, 0)));
	mainscene.add_sphere(Sphere(Vector(0., 0., 1000), 940, Vector(1., 1., 1.), true, false, Vector(0, 0, 0)));
	mainscene.add_sphere(Sphere(Vector(990, 0., 0), 940, Vector(1., 1., 1.), true, false, Vector(0, 0, 0)));
	mainscene.add_sphere(Sphere(Vector(-990, 0., 0), 940, Vector(1., 1., 1.), true, false, Vector(0, 0, 0)));

	std::vector<unsigned char> image(W*H * 3, 0);
	#pragma omp parallel for
	for (int i = 0; i < H; i++) {
		std::cout << i / H << std::endl;
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
				I = I + mainscene.getColor(rini, 2);
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