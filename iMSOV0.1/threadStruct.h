#pragma once

#include "varray.h"
#include "XVec.h"
#include "pointXd.h"

using namespace base;

//���̼߳���������Ҫ�ṹ��

//�� �ṹ��
struct threadParamSpline {
	varray<varray<Vec3>> l2;
	varray<varray<Vec3>> p3d2;
};
//�� �ṹ��
struct threadParamSplineVOL {
	varray<varray<varray<Vec3>>> q3;
	varray<varray<varray<Vec3>>> l3;
};
struct threadParam {
	varray<varray<point3d>> l2;
	varray<varray<point3d>> p3d2;
};
struct threadParamVOL {
	varray<varray<varray<point3d>>> q3;
	varray<varray<varray<point3d>>> l3;
};