#include "stdafx.h"
#include "globalFunc.h"

//template<class _T>
void Global::ReduceVarrayDim(const varray<varray<point4d>>& hightDimVarray, varray<point4d>& lowDimVarray, const bool clear)
{
	if (clear)
		lowDimVarray.clear();
	if (hightDimVarray.size() == 0)return;

	for (int i = 0; i < hightDimVarray.size(); ++i)
		for (int j = 0; j < hightDimVarray[i].size(); ++j)
			lowDimVarray.push_back(hightDimVarray[i][j]);
}

template<class _T>
void Global::Transpose(varray<varray<_T>>& Varr)
{
	varray<varray<_T>> tvarr;
	tvarr.resize(Varr[0].size());
	for (int i = 0; i < Varr.size(); ++i)
		for (int j = 0; j < Varr[i].size(); ++j)
			tvarr[j].push_back(Varr[i][j]);
	Varr = tvarr;
}

double Global::Bin(const int x, const int y)
{
	int x1 = x, y1 = y, d1 = x - y;
	for (int i = x - 1; i > 1; --i)
		x1 *= i;
	for (int i = y - 1; i > 1; --i)
		y1 *= i;
	for (int i = x - y - 1; i > 1; --i)
		d1 *= i;

	if (x1 == 0)x1 = 1;
	if (y1 == 0)y1 = 1;
	if (d1 == 0)d1 = 1;

	return 1.0*x1 / (y1*d1);
}

//模板函数分文件编译不通过
//template<typename _T>
//varray<_T> Global::transVectorToVarray(vector<_T> a)
//{
//	varray<_T> b;
//	for (int i = 0; i != a.size(); ++i)
//	{
//		b.push_back(a[i]);
//	}
//	return b;
//}


varray<double> Global::transVectorToVarray(vector<double> a)
{
	varray<double> b;
	for (int i = 0; i != a.size(); ++i)
	{
		b.push_back(a[i]);
	}
	return b;
}
varray<point4d> Global::transVectorToVarray(vector<point4d> a)
{
	varray<point4d> b;
	for (int i = 0; i != a.size(); ++i)
	{
		b.push_back(a[i]);
	}
	return b;
}
varray<point3d> Global::transVectorToVarray(vector<point3d> a)
{
	varray<point3d> b;
	for (int i = 0; i != a.size(); ++i)
	{
		b.push_back(a[i]);
	}
	return b;
}


vector<double> Global::transVarrayToVector(varray<double> a)
{
	vector<double> b;
	for (int i = 0; i != a.size(); ++i)
	{
		b.push_back(a[i]);
	}
	return b;
}

vector<point4d> Global::transVarrayToVector(varray<point4d> a)
{
	vector<point4d> b;
	for (int i = 0; i != a.size(); ++i)
	{
		b.push_back(a[i]);
	}
	return b;
}

vector<point3d> Global::transVarrayToVector(varray<point3d> a)
{
	vector<point3d> b;
	for (int i = 0; i != a.size(); ++i)
	{
		b.push_back(a[i]);
	}
	return b;
}

void Global::Normalization(varray<point3d>& pts, const double a, const double b)
{
	size_t sz = pts.size();
	double pmax = 0.0;
	for (int i = 0; i < sz; ++i)
	{
		for (int j = 0; j < 3; ++j)
			pmax = abs(pts[i][j]) > pmax ? abs(pts[i][j]) : pmax;
	}
	for (int i = 0; i < sz; ++i)
		pts[i] = (pts[i] * (b - a) / pmax).Trans(a, a, a);
}

void Global::Normalization(varray<varray<point3d>>& pts, const double a, const double b)
{
	double pmax = 0.0;
	for (int i = 0; i < pts.size(); ++i)
		for (int j = 0; j < pts[i].size(); ++j)
			for (int k = 0; k < 3; ++k)
				pmax = abs(pts[i][j][k]) > pmax ? abs(pts[i][j][k]) : pmax;

	for (int i = 0; i < pts.size(); ++i)
		for (int j = 0; j < pts[i].size(); ++j)
			pts[i][j] = (pts[i][j] * (b - a) / pmax).Trans(a, a, a);
}

void Global::GetCenter(varray<point3d>& pts, point3d & center)
{
	point3d sum;
	for (int i = 0; i < pts.size(); ++i)
		sum += pts[i];
	center = sum / pts.size();
}

template<class T>
static void Global::TransPts(varray<T> &pts, const T &vec, bool mode, const T &center /*= T()*/)
{
	if (mode)
	{
		for (int i = 0; i < pts.size(); ++i)
			pts[i] = pts[i] + (vec - center);
	}
	else
	{
		for (int i = 0; i < pts.size(); ++i)
			pts[i] = pts[i] + vec;
	}
}

Vector4d Global::P4dToV4d(const point4d P4d)
{
	Vector4d L2;
	L2(0) = P4d.x;
	L2(1) = P4d.y;
	L2(2) = P4d.z;
	L2(3) = P4d.w;
	return L2;
}

Vec4 Global::Point4dToVec4(const point4d p4d)
{
	Vec4 vec4;
	vec4.x = p4d.x;
	vec4.y = p4d.y;
	vec4.z = p4d.z;
	vec4.w = p4d.w;
	return vec4;
}

varray<Vec4> Global::Point4dToVec4(const varray<point4d> VP4d)
{
	varray<Vec4>VVec4;
	for (int i = 0; i < VP4d.size(); i++)
	{
		VVec4.push_back(Point4dToVec4(VP4d[i]));
	}
	return VVec4;
}



point4d Global::V4dToP4d(const Vector4d V4d)
{
	point4d L2;
	L2.x = V4d(0);
	L2.y = V4d(1);
	L2.z = V4d(2);
	L2.w = V4d(3);
	return L2;
}

void Global::Vector4dToUnit(Vector4d &v)
{
	Vector4d Unit_v;
	Unit_v(0) = v(0) / sqrt(pow(v(0), 2) + pow(v(1), 2) + pow(v(2), 2));
	Unit_v(1) = v(1) / sqrt(pow(v(0), 2) + pow(v(1), 2) + pow(v(2), 2));
	Unit_v(2) = v(2) / sqrt(pow(v(0), 2) + pow(v(1), 2) + pow(v(2), 2));
	v = Unit_v;
}

void Global::AcrossB(const Vector4d & a, const Vector4d & b, Vector4d & c)
{
	c(0) = a(1) * b(2) - a(2) * b(1);
	c(1) = a(2) * b(0) - a(0) * b(2);
	c(2) = a(0) * b(1) - a(1) * b(0);
}

double Global::angleAB(const Vector4d & va, const Vector4d & vb)
{

	point3d a(va[0], va[1], va[2]);
	point3d b(vb[0], vb[1], vb[2]);
	double al = a.Magnitude();
	double bl = b.Magnitude();
	if (al == 0)
		return -1;
	if (bl == 0)
		return -1;

	double cosang = a.Dot(b) / (al*bl);
	double res = acos(cosang) * 180 / PI;
	return res;
}

double Global::CurveLength(const varray<point3d>& L)
{
	double length = 0;
	for (int q = 0; q < L.size() - 1; q++)//计算曲线的长度
	{
		length += sqrt(pow(L[q + 1].x - L[q].x, 2) + pow(L[q + 1].y - L[q].y, 2) + pow(L[q + 1].z - L[q].z, 2));
	}
	return length;
}

//template<class numtype>
bool Global::IsInSet(const double a, const vector<double>& S)
{
	for (int i = 0; i < S.size(); ++i)
	{
		if (a == S[i])
			return true;
	}
	return false;
}

void Global::CalCoordinatesTransMat(const varray<point3d>& objCoordinates, const varray<point3d>& targetCoordinates, Matrix4d& mat)
{
	Matrix4d tran1;//世界坐标到objCoordinates的变换矩阵
	tran1 << objCoordinates[0].x, objCoordinates[0].y, objCoordinates[0].z, 0,
		objCoordinates[1].x, objCoordinates[1].y, objCoordinates[1].z, 0,
		objCoordinates[2].x, objCoordinates[2].y, objCoordinates[2].z, 0,
		objCoordinates[3].x, objCoordinates[3].y, objCoordinates[3].z, 1;


	Matrix4d tran2;//世界坐标到targetCoordinates的变换矩阵
	tran2 << targetCoordinates[0].x, targetCoordinates[0].y, targetCoordinates[0].z, 0,
		targetCoordinates[1].x, targetCoordinates[1].y, targetCoordinates[1].z, 0,
		targetCoordinates[2].x, targetCoordinates[2].y, targetCoordinates[2].z, 0,
		targetCoordinates[3].x, targetCoordinates[3].y, targetCoordinates[3].z, 1;

	mat = tran2.transpose()*(tran1.transpose().inverse());
}

void Global::CalTransMat(const Vector4d & targetPts, const Vector4d & targetVec,	const Vector4d & objPts, const Vector4d & objVec, Matrix4d & mat)
{
	point4d v_tar, v_obj;
	v_tar = V4dToP4d(targetVec);
	v_obj = V4dToP4d(objVec);

	//平移
	Matrix4d T0, T;

	T0 << 1, 0, 0, 0,//交点到坐标原点的平移矩阵
		0, 1, 0, 0,
		0, 0, 1, 0,
		-objPts(0), -objPts(1), -objPts(2), 1;

	T << 1, 0, 0, 0,//坐标原点到交点的平移矩阵
		0, 1, 0, 0,
		0, 0, 1, 0,
		targetPts(0), targetPts(1), targetPts(2), 1;

	//两向量夹角
	double ang = v_obj.Angle(v_tar)*PI / 180;

	//旋转轴
	point4d axis = (v_obj.Cross(v_tar)).Normalize();
	Vector4d axisVec;
	axisVec = P4dToV4d(axis);

	//罗德里格斯变换

	Matrix4d I, A, R;
	I << 1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0,
		0, 0, 0, 1;

	A << 0, -axisVec(2), axisVec(1), 0,
		axisVec(2), 0, -axisVec(0), 0,
		-axisVec(1), axisVec(0), 0, 0,
		0, 0, 0, 1;

	R = cos(ang)*I + (1 - cos(ang))*(axisVec*axisVec.transpose()) + sin(ang)*A.transpose();

	for (int i = 0; i < 4; ++i)
	{
		R(i, 3) = 0;
		R(3, i) = 0;
	}
	R(3, 3) = 1;

	mat = T.transpose()*R.transpose()*T0.transpose();
}

void Global::CalTransMat(const point4d & targetPts, const point4d & targetVec, const point4d & objPts, const point4d & objVec, Matrix4d & mat)
{
	Vector4d vtp, vtv, vop, vov;
	vtp = P4dToV4d(targetPts);
	vtv = P4dToV4d(targetVec);
	vop = P4dToV4d(objPts);
	vov = P4dToV4d(objVec);
	CalTransMat(vtp, vtv, vop, vov, mat);
}

void Global::TransByMat(varray<Vector4d>& varr, const Matrix4d & mat)
{
	for (int i = 0; i < varr.size(); ++i)
	{
		varr[i] = mat * varr[i];
	}
}

void Global::TransByMat(varray<point4d>& varr, const Matrix4d& mat)
{
	for (int i = 0; i < varr.size(); ++i)
	{
		Vector4d v4d = P4dToV4d(varr[i]);
		v4d[3] = 1;
		v4d = mat * v4d;
		v4d[3] = varr[i].w;
		varr[i] = V4dToP4d(v4d);
	}
}

void Global::TransByMat(varray<Vec4>& varr, const Matrix4d & mat)
{

	for (int i = 0; i < varr.size(); ++i)
	{
		point4d p2;
		p2 = { varr[i].x, varr[i].y, varr[i].z, varr[i].w };
		auto p = mat * P4dToV4d(p2);
		varr[i].x = p(0);
		varr[i].y = p(1);
		varr[i].z = p(2);
		varr[i].w = p(3);
	}
}

void Global::TransByMat(vector<point4d>& varr, const Matrix4d& mat)
{
	for (int i = 0; i < varr.size(); ++i)
	{
		Vector4d v4d = P4dToV4d(varr[i]);
		v4d[3] = 1;
		v4d = mat * v4d;
		v4d[3] = varr[i].w;
		varr[i] = V4dToP4d(v4d);
	}
}

//static int Sgn(double a)
//{
//
//	if (a < 0)
//		return -1;
//	if (a == 0)
//		return 0;
//	return 1;
//}
int Global::Sgn(double a)
{

	if (a < 0)
		return -1;
	if (a == 0)
		return 0;
	return 1;
}

template<class _T>
void Global::SetBeginDir(varray<_T>& varr, const size_t beginIdx, const bool dir /*= true*/)
{
	varray<_T> res;
	typename varray<_T>::iterator it0, it;
	it = it0 = varr.begin() + beginIdx;
	do
	{
		res.push_back(*it);
		if (dir)
		{
			++it;
			if (it == varr.end())
				it = varr.begin();
		}
		else
		{
			--it;
			if (it == varr.begin())
			{
				res.push_back(*it);
				it = varr.end() - 1;
			}
		}
	} while (it != it0);
	varr = res;
}

int Global::GCD(int a, int b)
{
	int c;
	if (a < b)
	{
		a = a + b;
		b = a - b;
		a = a - b;
	}
	c = a % b;
	while (a%b != 0)
	{
		a = b;
		b = c;
		c = a % b;
	}
	return b;
}

//template<class double>
void Global::KnotUnify(const varray<double>& KnotsA, const varray<double>& KnotsB, varray<double>& NewKnots)
{
	NewKnots.clear();
	int Anum = 0, Bnum = 0;
	// 如果Anum或Bnum有一个遍历完了就停止
	while (Anum < KnotsA.size() && Bnum < KnotsB.size())
	{
		if (KnotsA[Anum] < KnotsB[Bnum])
		{
			NewKnots.push_back(KnotsA[Anum++]);
			continue;
		}
		else if (KnotsA[Anum] > KnotsB[Bnum])
		{
			NewKnots.push_back(KnotsB[Bnum++]);
			continue;
		}
		else
		{
			NewKnots.push_back(KnotsB[Bnum++]);
			Anum++;
		}
	}
	// 如果 Anum< KnotsA.size()，说明Bnum都添加进数组NewKnots了，只剩下比Bnum大的Anum，逐个添加
	while (Anum < KnotsA.size())
	{
		NewKnots.push_back(KnotsA[Anum++]);
	}
	while (Bnum < KnotsB.size())
	{
		NewKnots.push_back(KnotsB[Bnum++]);
	}
}

//template<class _T>
void Global::KnotsDiff(const varray<double>& KnotsL, const varray<double>& KnotsS, varray<double>& diffKnots)
{
	diffKnots.clear();
	if (KnotsL.size() <= KnotsS.size())
		return;

	int Ln = 0, Sn = 0;
	while (Ln < KnotsL.size() && Sn < KnotsS.size())
	{
		if (KnotsS[Sn] != KnotsL[Ln])
			diffKnots.push_back(KnotsL[Ln]);
		else
			++Sn;
		++Ln;
	}
	while (Ln < KnotsL.size())
	{
		diffKnots.push_back(KnotsL[Ln++]);
	}
}

int Global::Project2Plane(const point3d& v, const point3d& p, const point3d& pts, point3d& res, const point3d& v1/* = point3d()*/)
{
	point3d v2 = v1;
	if (v1 == point3d())
	{
		v2 = -v;
	}
	double a = v.Dot(v2);
	if (a == 0)
		return 0;
	point3d dp = p - pts;
	double t = dp.Dot(v) / a;
	res = pts + v1 * t;
	return 1;
}

int Global::Project2Plane(const point3d & v, const point3d & p, const varray<point4d>& Cpts, varray<point4d>& ProPts, const point3d & v1/* = point3d()*/)
{
	ProPts.clear();
	for (int i = 0; i < Cpts.size(); ++i)
	{
		point3d res;
		int isP = Project2Plane(v, p, Cpts[i], res, v1);
		if (isP)
			ProPts.push_back(res);
		else
			return 0;
	}
	return 1;
}