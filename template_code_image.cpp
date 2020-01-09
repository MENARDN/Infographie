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

Vector operator/(const Vector& A, const double& b) {
	return Vector(A[0] / b, A[1] / b, A[2] / b);
}

class Ray {
public:
	Ray(const Vector& C, Vector u) : C(C), u(u) {};
	Vector C, u;
};

class Sphere {
public:
	Sphere(const Vector& O, double R) : O(O), R(R) {};


	bool intersect(const Ray& r, Vector &P, Vector &N) {
		double a = 1;
		double b = 2 * dot(r.u, r.C - O);
		double c = (r.C - 0).Norm2() - R * R;

		double delta = b * b - 4 * a * c;
		if (delta < 0) return false;
		double sqrtDelta = sqrt(delta);
		double t1 = (-b + sqrtDelta) / (2 * a);
		double t0 = (-b + sqrtDelta) / (2 * a);
		double t;
		if (t1 > t0) {
			t = t0;
		}
		else {
			t = t1;
		}
		P = r.C + t * r.u;
		N = O-P;
		N.normalize();
		return true;
	}

	Vector O;
	double R;
};

class Scene {
public:
	Scene() {}
};

int main() {
	int W = 512;
	int H = 512;
	double fov = 60 * M_PI / 180;
	Sphere s(Vector(0., 0., 0.), 10);
	Vector C(0, 0, 55);
	Vector albedo(1., 0.5, 0);
	Vector L(-10, 20, 40);
	double l = 10000000;

	std::vector<unsigned char> image(W*H * 3, 0);
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			double d = W / (2 * tan(fov/2));
			Vector u(j - W / 2, -i + H / 2, -d);
			u = u + C;
			u.normalize();
			Ray r(C, u);
			Vector P, N;
			bool has_inter = s.intersect(r,P,N);

			if (has_inter) {
				double lightdist = sqrt((L - P).Norm2());

				Vector I = (l * albedo / M_PI * dot(N, (L - P) / lightdist)) / (lightdist * lightdist);
				image[(i * W + j) * 3 + 0] = std::min(255, std::max(0,int(I[0])));
				image[(i * W + j) * 3 + 1] = std::min(255, std::max(0, int(I[1])));
				image[(i * W + j) * 3 + 2] = std::min(255, std::max(0, int(I[2])));
			}
			else {
				image[(i * W + j) * 3 + 0] = 0;
				image[(i * W + j) * 3 + 1] = 0;
				image[(i * W + j) * 3 + 2] = 0;
			}
		}
	}
	stbi_write_png("image.png", W, H, 3, &image[0], 0);

	return 0;
}