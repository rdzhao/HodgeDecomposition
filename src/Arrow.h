#pragma once

/* Create 3d arrow */

#include <iostream>
#include <cmath>
#include <vector>

#include "glm/glm.hpp"

using namespace glm;

const double PI = 3.14159265358979323846264338;

class Arrow {
public:
	Arrow() {
		radius = 1;
		tail = vec3(0, 0, 0);
		head = vec3(1, 0, 0);
		create();
	}

	Arrow(double r, vec3 t, vec3 h) {
		radius = r;
		tail = t;
		head = h;
		create();
	}

	std::vector<float>& getVertices() {
		return vertices;
	}

	std::vector<int>& getFaces() {
		return faces;
	}

private:
	void create(){
		// number of lines to partition a circle
		int num = 16; 

		vec3 d = head - tail;

		if(glm::dot(d,d)<0.00000001)
			return;

		//find a vertical vector
		vec3 f1(0, 0, 0);
		if (d.x != 0) {
			f1.x = d.y;
			f1.y = -d.x;
		}
		else if (d.y != 0) {
			f1.y = d.z;
			f1.z = -d.y;
		}
		else {
			f1.z = d.x;
			f1.x = -d.z;
		}
		//std::cout << f1.x << " " << f1.y << " " << f1.z << std::endl;
		f1 /= sqrt(dot(f1, f1));
		//std::cout << f1.x << " " << f1.y << " " << f1.z << std::endl;
		f1 *= radius;
		//std::cout << f1.x<<" " << f1.y << " " << f1.z << std::endl;

		vec3 f2 = cross(d, f1);
		f2 /= sqrt(dot(f2, f2));
		f2 *= radius;
		//std::cout << f2.x << " " << f2.y << " " << f2.z << std::endl;
		// f1 f2 are frames
		
		// create tail center vertex
		//vertices.push_back(tail);
		vertices.push_back(tail.x);
		vertices.push_back(tail.y);
		vertices.push_back(tail.z);

		// create tail circle vertices
		for (int i = 0; i < num; ++i) {
			float cosine, sine;
			cosine = static_cast<float>(cos((i*1.0 / num) * 2 * PI));
			sine = static_cast<float>(sin((i*1.0 / num) * 2 * PI));
			vec3 pos = cosine * f1 + sine * f2 + tail;
			//vertices.push_back(cosine*f1 + sine*f2 + tail);
			vertices.push_back(pos.x);
			vertices.push_back(pos.y);
			vertices.push_back(pos.z);
		}

		// create head small circle vertices
		for (int i = 0; i < num; ++i) {
			float cosine, sine;
			cosine = static_cast<float>(cos((i*1.0 / num) * 2 * PI));
			sine = static_cast<float>(sin((i*1.0 / num) * 2 * PI));
			vec3 pos = cosine * f1 + sine * f2 + head;
			//vertices.push_back(cosine*f1 + sine*f2 + head);
			vertices.push_back(pos.x);
			vertices.push_back(pos.y);
			vertices.push_back(pos.z);
		}

		// create head big circle vertices
		float ratio = 2;
		for (int i = 0; i < num; ++i) {
			float cosine, sine;
			cosine = static_cast<float>(cos((i*1.0 / num) * 2 * PI));
			sine = static_cast<float>(sin((i*1.0 / num) * 2 * PI));
			vec3 pos = ratio * (cosine * f1 + sine * f2) + head;
			//vertices.push_back(cosine*f1 + sine*f2 + head);
			vertices.push_back(pos.x);
			vertices.push_back(pos.y);
			vertices.push_back(pos.z);
		}

		// create sharp vertex
		vec3 p = head + d * 0.5f;
		//vertices.push_back(head + d * 0.5f);
		vertices.push_back(p.x);
		vertices.push_back(p.y);
		vertices.push_back(p.z);

		// create faces
		// bot faces
		for (int i = 0; i < num; ++i) {
			int zz = i + 2;
			if (zz > num)
				zz %= num;
			//std::cout << zz << std::endl;
			//faces.push_back(vec3(0, zz, i + 1));
			faces.push_back(0);
			faces.push_back(zz);
			faces.push_back(i+1);
		}
		
		// around faces
		for (int i = 0; i < num; ++i) {
			int zz1 = i + 2;
			if (zz1 > num)
				zz1 %= num;

			//faces.push_back(vec3(i + 1, zz1, num + i + 1));
			faces.push_back(i + 1);
			faces.push_back(zz1);
			faces.push_back(num + i + 1);
			//faces.push_back(vec3(zz1, zz1 + num, num + i + 1));
			faces.push_back(zz1);
			faces.push_back(zz1 + num);
			faces.push_back(num + i + 1);
		}

		// ring faces
		for (int i = 0; i < num; ++i) {
			int zz1 = i + 2;
			if (zz1 > num)
				zz1 %= num;

			//faces.push_back(vec3(i + 1, zz1, num + i + 1));
			faces.push_back(i + 1 + num);
			faces.push_back(zz1 + num);
			faces.push_back(2 * num + i + 1);
			//faces.push_back(vec3(zz1, zz1 + num, num + i + 1));
			faces.push_back(zz1 + num);
			faces.push_back(zz1 + 2 * num);
			faces.push_back(2 * num + i + 1);
		}

		// top faces
		for (int i = 0; i < num; ++i) {
			int zz1 = i + 2;
			if (zz1 > num)
				zz1 %= num;

			//faces.push_back(vec3(2 * num + 1, num + i + 1, zz1 + num));
			faces.push_back(3 * num + 1);
			faces.push_back(2*num + i + 1);
			faces.push_back(zz1 + 2 * num);
		}
	}

	double radius;
	vec3 tail;
	vec3 head;

	//model
	std::vector<float> vertices;
	std::vector<int> faces;
};