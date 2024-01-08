#pragma once

#include"vec3.h"

//ray class:origin and direction compose a ray
class ray {
public:
	ray(){}
	ray(const point3& origin, const vec3& direction) : orig(origin), dir(direction) {}

	point3 origin() const { return orig; };
	vec3 direction() const { return dir; }

	//the result of at function is the point the ray get at time t
	point3 at(double t) const {
		return orig + t * dir;
	}
private:
	point3 orig;
	vec3 dir;
};