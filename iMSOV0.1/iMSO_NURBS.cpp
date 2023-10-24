#include "iMSO_NURBS.h"
// NURBS.cpp : 定义 DLL 应用程序的导出函数。
#include "stdafx.h"


NurbsLine::NurbsLine()
{
	_u_Degree = 1;
	_u_Knots = make_shared<vector<double>>();
	_u_Num = 2;
	_ControlPts = make_shared<vector<point4d>>();
}

NurbsLine::NurbsLine(int u_num, int u_degree, vector<double> u_knots, vector<point4d> Cp)
	:_u_Degree(u_degree)
	, _u_Num(u_num)
{
	if (!isParameterCorrect(Cp.size(), u_degree, u_knots.size()))
	{
		throw std::exception();
	}

	_u_Knots = make_shared<vector<double>>();
	*_u_Knots = u_knots;

	_ControlPts = make_shared<vector<point4d>>();
	*_ControlPts = Cp;
}

NurbsLine::NurbsLine(const NurbsLine &nurbs)
{
	_u_Degree = nurbs._u_Degree;
	_u_Knots = make_shared<vector<double>>();
	for (auto iter = nurbs._u_Knots->begin(); iter != nurbs._u_Knots->end(); ++iter)
		_u_Knots->push_back(*iter);
	_u_Num = nurbs.GetUNum();
	_ControlPts = make_shared<vector<point4d>>();
	for (auto iter = nurbs._ControlPts->begin(); iter != nurbs._ControlPts->end(); ++iter)
		_ControlPts->push_back(*iter);
}

NurbsLine& NurbsLine::operator=(const NurbsLine& nurbs)
{
	if (this == &nurbs)
	{
		return *this;
	}
	_u_Degree = nurbs._u_Degree;
	_u_Knots = make_shared<vector<double>>();
	for (auto iter = nurbs._u_Knots->begin(); iter != nurbs._u_Knots->end(); ++iter)
		_u_Knots->push_back(*iter);
	_u_Num = nurbs.GetUNum();
	_ControlPts = make_shared<vector<point4d>>();
	for (auto iter = nurbs._ControlPts->begin(); iter != nurbs._ControlPts->end(); ++iter)
		_ControlPts->push_back(*iter);
	return *this;
}


NurbsLine::~NurbsLine()
{
}

int NurbsLine::GetUDegree()const
{
	return _u_Degree;
}

shared_ptr<vector<double>> NurbsLine::GetUKonts()const
{
	return _u_Knots;
}

int NurbsLine::GetUNum()const
{
	return _u_Num;
}

shared_ptr<vector<point4d>> NurbsLine::GetControlPointer()const
{
	return _ControlPts;
}

int NurbsSurface::GetUDegree()const
{
	return _u_Degree;
}

shared_ptr<vector<double>> NurbsSurface::GetUKonts()const
{
	return _u_Knots;
}

int NurbsSurface::GetUNum()const
{
	return _u_Num;
}

shared_ptr<vector<point4d>> NurbsSurface::GetControlPointer()const
{
	return _ControlPts;
}

int NurbsVol::GetUDegree()const
{
	return _u_Degree;
}

shared_ptr<vector<double>> NurbsVol::GetUKonts()const
{
	return _u_Knots;
}

int NurbsVol::GetUNum()const
{
	return _u_Num;
}

shared_ptr<vector<point4d>> NurbsVol::GetControlPointer()const
{
	return _ControlPts;
}

bool NurbsLine::SetControlPoint(const vector<point4d> cp)
{
	if (cp.empty())
	{
		return false;
	}
	_ControlPts->clear();
	for (auto iter = cp.begin(); iter != cp.end(); ++iter)
	{
		_ControlPts->push_back(*iter);
	}
	return true;
}

bool NurbsLine::SetControlPoint(shared_ptr<vector<point4d>> controlPts)
{
	if (controlPts->empty())
		return false;
	_ControlPts = controlPts;
	return true;
}

bool NurbsLine::SetUDegree(const int degree)
{
	if (degree < 0)
		return false;
	_u_Degree = degree;
	return true;
}

bool NurbsLine::SetUKonts(const vector<double> knots)
{
	if (knots.empty())
	{
		return false;
	}
	_u_Knots->clear();
	for (auto iter = knots.begin(); iter != knots.end(); ++iter)
	{
		_u_Knots->push_back(*iter);
	}
	return true;
}

bool NurbsLine::SetUKonts(shared_ptr<vector<double>> u_knots)
{
	if (u_knots->empty())
		return false;
	_u_Knots = u_knots;
	return true;
}

bool NurbsLine::SetUNum(const int num)
{
	if (num < 0)
		return false;
	_u_Num = num;
	return true;
}

bool NurbsSurface::SetControlPoint(const vector<point4d> cp)
{
	if (cp.empty())
	{
		return false;
	}
	_ControlPts->clear();
	for (auto iter = cp.begin(); iter != cp.end(); ++iter)
	{
		_ControlPts->push_back(*iter);
	}
	return true;
}

bool NurbsSurface::SetControlPoint(shared_ptr<vector<point4d>> controlPts)
{
	if (controlPts->empty())
		return false;
	_ControlPts = controlPts;
	return true;
}

bool NurbsSurface::SetUDegree(const int degree)
{
	if (degree < 0)
		return false;
	_u_Degree = degree;
	return true;
}

bool NurbsSurface::SetUKonts(const vector<double> knots)
{
	if (knots.empty())
	{
		return false;
	}
	_u_Knots->clear();
	for (auto iter = knots.begin(); iter != knots.end(); ++iter)
	{
		_u_Knots->push_back(*iter);
	}
	return true;
}

bool NurbsSurface::SetUKonts(shared_ptr<vector<double>> u_knots)
{
	if (u_knots->empty())
		return false;
	_u_Knots = u_knots;
	return true;
}

bool NurbsSurface::SetUNum(const int num)
{
	if (num < 0)
		return false;
	_u_Num = num;
	return true;
}
bool NurbsVol::SetControlPoint(const vector<point4d> cp)
{
	if (cp.empty())
	{
		return false;
	}
	_ControlPts->clear();
	for (auto iter = cp.begin(); iter != cp.end(); ++iter)
	{
		_ControlPts->push_back(*iter);
	}
	return true;
}

bool NurbsVol::SetControlPoint(shared_ptr<vector<point4d>> controlPts)
{
	if (controlPts->empty())
		return false;
	_ControlPts = controlPts;
	return true;
}

bool NurbsVol::SetUDegree(const int degree)
{
	if (degree < 0)
		return false;
	_u_Degree = degree;
	return true;
}

bool NurbsVol::SetUKonts(const vector<double> knots)
{
	if (knots.empty())
	{
		return false;
	}
	_u_Knots->clear();
	for (auto iter = knots.begin(); iter != knots.end(); ++iter)
	{
		_u_Knots->push_back(*iter);
	}
	return true;
}

bool NurbsVol::SetUKonts(shared_ptr<vector<double>> u_knots)
{
	if (u_knots->empty())
		return false;
	_u_Knots = u_knots;
	return true;
}

bool NurbsVol::SetUNum(const int num)
{
	if (num < 0)
		return false;
	_u_Num = num;
	return true;
}

bool NurbsLine::isParameterCorrect(int num, int degree, int knots)
{

	if ((num + degree + 1) == knots)
	{
		return true;
	}
	else
	{
		return false;
	}
}

void NurbsLine::CreatLineWithTwoPoints(const point4d & p1, const point4d & p2, int degree)
{
	this->_ControlPts->clear();
	this->_u_Knots->clear();

	this->_u_Degree = 1;
	this->_u_Knots->push_back(0);
	this->_u_Knots->push_back(0);
	this->_u_Knots->push_back(1);
	this->_u_Knots->push_back(1);
	this->_ControlPts->push_back(point4d(p1));
	this->_ControlPts->push_back(point4d(p2));
	if (degree > 1) {
		this->DegreeElevate(degree);
	}
}

int NurbsLine::FindSpan(const double x, const int degree, const int CtrlPtsNum, const vector<double>& knots)//varray or vector
{
	int left = degree, right = CtrlPtsNum, mid = 0;
	if (x >= knots[right])
		return right - 1;
	mid = (left + right) / 2;
	while (x < knots[mid] || x >= knots[mid + 1])
	{
		if (x < knots[mid])
			right = mid;
		else
			left = mid;
		mid = (right + left) / 2;
	}
	return mid;
}

int NurbsLine::FindSpan(const double x, const int degree, const int CtrlPtsNum, const varray<double>& knots)//varray or vector
{
	int left = degree, right = CtrlPtsNum, mid = 0;
	if (x >= knots[right])
		return right - 1;
	mid = (left + right) / 2;
	while (x < knots[mid] || x >= knots[mid + 1])
	{
		if (x < knots[mid])
			right = mid;
		else
			left = mid;
		mid = (right + left) / 2;
	}
	return mid;
}

void NurbsLine::BasisFuns(const double u, const int k, const int degree, const vector<double>& knots,vector<double>& N)
{
	int p = degree;
	varray<double>  left, right;
	right.resize(p + 1, 0);
	left.resize(p + 1, 0);
	N.resize(p + 1, 0);

	N[0] = 1.0;
	double saved, temp;
	for (int j = 1; j <= p; j++)
	{
		left[j] = u - knots[k + 1 - j];
		right[j] = knots[k + j] - u;
		saved = 0;
		for (int r = 0; r < j; r++)
		{
			temp = N[r] / (right[r + 1] + left[j - r]);
			N[r] = saved + right[r + 1] * temp;
			saved = left[j - r] * temp;
		}
		N[j] = saved;
	}
	left.clear();
	right.clear();
}
void NurbsLine::BasisFuns(const double u, const int k, const int degree, const varray<double>& knots, varray<double>& N)
{
	int p = degree;
	varray<double>  left, right;
	right.resize(p + 1, 0);
	left.resize(p + 1, 0);
	N.resize(p + 1, 0);

	N[0] = 1.0;
	double saved, temp;
	for (int j = 1; j <= p; j++)
	{
		left[j] = u - knots[k + 1 - j];
		right[j] = knots[k + j] - u;
		saved = 0;
		for (int r = 0; r < j; r++)
		{
			temp = N[r] / (right[r + 1] + left[j - r]);
			N[r] = saved + right[r + 1] * temp;
			saved = left[j - r] * temp;
		}
		N[j] = saved;
	}
	left.clear();
	right.clear();
}


int NurbsSurface::FindSpan(const double x, const int degree, const int CtrlPtsNum, const vector<double>& knots)const//varray or vector
{
	int left = degree, right = CtrlPtsNum, mid = 0;
	if (x >= knots[right])
		return right - 1;
	mid = (left + right) / 2;
	while (x < knots[mid] || x >= knots[mid + 1])
	{
		if (x < knots[mid])
			right = mid;
		else
			left = mid;
		mid = (right + left) / 2;
	}
	return mid;
}

void NurbsSurface::BasisFuns(const double u, const int k, const int degree, const vector<double>& knots,
	vector<double>& N)const
{
	int p = degree;
	varray<double>  left, right;
	right.resize(p + 1, 0);
	left.resize(p + 1, 0);
	N.resize(p + 1, 0);

	N[0] = 1.0;
	double saved, temp;
	for (int j = 1; j <= p; j++)
	{
		left[j] = u - knots[k + 1 - j];
		right[j] = knots[k + j] - u;
		saved = 0;
		for (int r = 0; r < j; r++)
		{
			temp = N[r] / (right[r + 1] + left[j - r]);
			N[r] = saved + right[r + 1] * temp;
			saved = left[j - r] * temp;
		}
		N[j] = saved;
	}
	left.clear();
	right.clear();
}

int NurbsVol::FindSpan(const double x, const int degree, const int CtrlPtsNum, const vector<double>& knots)const//varray or vector
{
	int left = degree, right = CtrlPtsNum, mid = 0;
	if (x >= knots[right])
		return right - 1;
	mid = (left + right) / 2;
	while (x < knots[mid] || x >= knots[mid + 1])
	{
		if (x < knots[mid])
			right = mid;
		else
			left = mid;
		mid = (right + left) / 2;
	}
	return mid;
}

void NurbsVol::BasisFuns(const double u, const int k, const int degree, const vector<double>& knots,
	vector<double>& N)const
{
	int p = degree;
	varray<double>  left, right;
	right.resize(p + 1, 0);
	left.resize(p + 1, 0);
	N.resize(p + 1, 0);

	N[0] = 1.0;
	double saved, temp;
	for (int j = 1; j <= p; j++)
	{
		left[j] = u - knots[k + 1 - j];
		right[j] = knots[k + j] - u;
		saved = 0;
		for (int r = 0; r < j; r++)
		{
			temp = N[r] / (right[r + 1] + left[j - r]);
			N[r] = saved + right[r + 1] * temp;
			saved = left[j - r] * temp;
		}
		N[j] = saved;
	}
	left.clear();
	right.clear();
}

float NurbsLine::OneBasisFun(int p, int m, const varray<double>& U, int i, double u)
{
	double *N = (double*)malloc(sizeof(double) * (p + 1));
	int j;
	int k;
	double saved, Uleft, Uright, temp;
	double Nip;

	if ((i == 0 && u == U[0]) ||
		(i == (m - p - 1) && u == U[m]))
	{
		return(1.0);
	}

	if (u < U[i] || u >= U[i + p + 1])
	{
		return(0.0);
	}
	for (j = 0; j <= p; j++)
	{
		if (u >= U[i + j] && u < U[i + j + 1]) N[j] = 1.0;
		else N[j] = 0.0;
	}
	for (k = 1; k <= p; k++)
	{
		if (N[0] == 0.0) saved = 0.0;
		else saved = ((u - U[i]) * N[0]) / (U[i + k] - U[i]);
		for (j = 0; j < (p - k + 1); j++)
		{
			Uleft = U[i + j + 1];
			Uright = U[i + j + k + 1];
			if (N[j + 1] == 0.0)
			{
				N[j] = saved; saved = 0.0;
			}
			else
			{
				temp = N[j + 1] / (Uright - Uleft);
				N[j] = saved + (Uright - u) * temp;
				saved = (u - Uleft) * temp;
			}
		}
	}
	Nip = N[0];
	free(N);
	return Nip;
}


void NurbsLine::AllBasisFuns(const double u, const int k, const int degree,
	const vector<double>& knots, vector<vector<double>>& ndu)
{
	int order = degree + 1;
	double saved, temp;
	varray<double> left, right;

	left.resize(order);
	right.resize(order);
	ndu.resize(order);

	ndu[0][0] = 1;
	for (int j = 1; j <= degree; j++)
	{
		left[j] = u - knots[k + 1 - j];
		right[j] = knots[k + j] - u;
		saved = 0.0f;
		for (int r = 0; r < j; r++)
		{
			ndu[j][r] = right[r + 1] + left[j - r];
			temp = ndu[r][j - 1] / ndu[j][r];
			ndu[r][j] = saved + right[r + 1] * temp;
			saved = left[j - r] * temp;
		}
		ndu[j][j] = saved;
	}
}

float NurbsSurface::OneBasisFun(int p, int m, const varray<double>& U, int i, double u)
{
	double *N = (double*)malloc(sizeof(double) * (p + 1));
	int j;
	int k;
	double saved, Uleft, Uright, temp;
	double Nip;

	if ((i == 0 && u == U[0]) ||
		(i == (m - p - 1) && u == U[m]))
	{
		return(1.0);
	}

	if (u < U[i] || u >= U[i + p + 1])
	{
		return(0.0);
	}
	for (j = 0; j <= p; j++)
	{
		if (u >= U[i + j] && u < U[i + j + 1]) N[j] = 1.0;
		else N[j] = 0.0;
	}
	for (k = 1; k <= p; k++)
	{
		if (N[0] == 0.0) saved = 0.0;
		else saved = ((u - U[i]) * N[0]) / (U[i + k] - U[i]);
		for (j = 0; j < (p - k + 1); j++)
		{
			Uleft = U[i + j + 1];
			Uright = U[i + j + k + 1];
			if (N[j + 1] == 0.0)
			{
				N[j] = saved; saved = 0.0;
			}
			else
			{
				temp = N[j + 1] / (Uright - Uleft);
				N[j] = saved + (Uright - u) * temp;
				saved = (u - Uleft) * temp;
			}
		}
	}
	Nip = N[0];
	free(N);
	return Nip;
}

void NurbsSurface::AllBasisFuns(const double u, const int k, const int degree,
	const vector<double>& knots, vector<vector<double>>& ndu)
{
	int order = degree + 1;
	double saved, temp;
	varray<double> left, right;

	left.resize(order);
	right.resize(order);
	ndu.resize(order);

	ndu[0][0] = 1;
	for (int j = 1; j <= degree; j++)
	{
		left[j] = u - knots[k + 1 - j];
		right[j] = knots[k + j] - u;
		saved = 0.0f;
		for (int r = 0; r < j; r++)
		{
			ndu[j][r] = right[r + 1] + left[j - r];
			temp = ndu[r][j - 1] / ndu[j][r];
			ndu[r][j] = saved + right[r + 1] * temp;
			saved = left[j - r] * temp;
		}
		ndu[j][j] = saved;
	}
}

float NurbsVol::OneBasisFun(int p, int m, const varray<double>& U, int i, double u)
{
	double *N = (double*)malloc(sizeof(double) * (p + 1));
	int j;
	int k;
	double saved, Uleft, Uright, temp;
	double Nip;

	if ((i == 0 && u == U[0]) ||
		(i == (m - p - 1) && u == U[m]))
	{
		return(1.0);
	}

	if (u < U[i] || u >= U[i + p + 1])
	{
		return(0.0);
	}
	for (j = 0; j <= p; j++)
	{
		if (u >= U[i + j] && u < U[i + j + 1]) N[j] = 1.0;
		else N[j] = 0.0;
	}
	for (k = 1; k <= p; k++)
	{
		if (N[0] == 0.0) saved = 0.0;
		else saved = ((u - U[i]) * N[0]) / (U[i + k] - U[i]);
		for (j = 0; j < (p - k + 1); j++)
		{
			Uleft = U[i + j + 1];
			Uright = U[i + j + k + 1];
			if (N[j + 1] == 0.0)
			{
				N[j] = saved; saved = 0.0;
			}
			else
			{
				temp = N[j + 1] / (Uright - Uleft);
				N[j] = saved + (Uright - u) * temp;
				saved = (u - Uleft) * temp;
			}
		}
	}
	Nip = N[0];
	free(N);
	return Nip;
}

void NurbsVol::AllBasisFuns(const double u, const int k, const int degree,
	const vector<double>& knots, vector<vector<double>>& ndu)
{
	int order = degree + 1;
	double saved, temp;
	varray<double> left, right;

	left.resize(order);
	right.resize(order);
	ndu.resize(order);

	ndu[0][0] = 1;
	for (int j = 1; j <= degree; j++)
	{
		left[j] = u - knots[k + 1 - j];
		right[j] = knots[k + j] - u;
		saved = 0.0f;
		for (int r = 0; r < j; r++)
		{
			ndu[j][r] = right[r + 1] + left[j - r];
			temp = ndu[r][j - 1] / ndu[j][r];
			ndu[r][j] = saved + right[r + 1] * temp;
			saved = left[j - r] * temp;
		}
		ndu[j][j] = saved;
	}
}

void NurbsLine::DerBasisFuns(const double u, const int k, const int degree, const int n, const vector<double>& knots,
	vector<vector<double>>& basisDu)
{
	int p = degree;
	int rn = n > p ? p : n;
	varray<double> left(p + 1, 0), right(p + 1, 0);
	varray<varray<double>>a(2), ndu(p + 1);
	basisDu.resize(n + 1);
	for (int i = 0; i < n + 1; ++i)
		basisDu[i].resize(k + 1, 0);
	for (int i = 0; i < 2; ++i)
		a[i].resize(p + 1, 0);
	for (int i = 0; i < p + 1; ++i)
		ndu[i].resize(p + 1, 0);

	ndu[0][0] = 1.0;
	int r, s1, s2, rk, pk, j1, j2;
	double d;
	double  saved, temp;
	for (int j = 1; j <= p; j++)
	{
		left[j] = u - knots[k + 1 - j];
		right[j] = knots[k + j] - u;
		saved = 0.0;
		for (r = 0; r < j; r++)
		{
			ndu[j][r] = right[r + 1] + left[j - r];
			temp = ndu[r][j - 1] / ndu[j][r];
			ndu[r][j] = saved + right[r + 1] * temp;
			saved = left[j - r] * temp;
		}
		ndu[j][j] = saved;
	}
	for (int j = 0; j <= p; ++j)
		basisDu[0][j] = ndu[j][p];
	for (r = 0; r <= p; r++)
	{
		a[0][0] = 1.0;
		s1 = 0; s2 = 1;
		for (int i = 1; i <= n; ++i)
		{
			d = 0.0;
			rk = r - rn;
			pk = p - rn;
			if (r >= rn)
			{
				a[s2][0] = a[s1][0] / ndu[pk + 1][rk];
				d = a[s2][0] * ndu[rk][pk];
			}
			if (rk >= -1) j1 = 1;
			else j1 = (-rk);
			if (r - 1 <= pk)  j2 = rn - 1;
			else j2 = p - r;
			for (int j = j1; j <= j2; j++)
			{
				a[s2][j] = (a[s1][j] - a[s1][j - 1])
					/ ndu[pk + 1][rk + j];
				d += a[s2][j] * ndu[rk + j][pk];
			}
			if (r <= pk)
			{
				a[s2][rn] = -a[s1][rn - 1] / ndu[pk + 1][r];
				d += a[s2][rn] * ndu[r][pk];

			}
			basisDu[i][r] = d;
			s1 = s1 + s2; s2 = s1 - s2; s1 = s1 - s2;
		}
	}
	r = p;
	for (int i = 1; i <= n; i++)
	{
		for (int j = 0; j <= p; j++)
			basisDu[i][j] *= r;
		r *= (p - k);
	}

}

void NurbsLine::DerBasisFuns(const double u, const int k, const int degree, const int n, const varray<double>& knots, varray<varray<double>>& basisDu)
{
	int p = degree;
	int rn = n > p ? p : n;
	varray<double> left(p + 1, 0), right(p + 1, 0);
	varray<varray<double>>a(2), ndu(p + 1);
	basisDu.resize(n + 1);
	for (int i = 0; i < n + 1; ++i)
		basisDu[i].resize(k + 1, 0);
	for (int i = 0; i < 2; ++i)
		a[i].resize(p + 1, 0);
	for (int i = 0; i < p + 1; ++i)
		ndu[i].resize(p + 1, 0);

	ndu[0][0] = 1.0;
	int r, s1, s2, rk, pk, j1, j2;
	double d;
	double  saved, temp;
	for (int j = 1; j <= p; j++)
	{
		left[j] = u - knots[k + 1 - j];
		right[j] = knots[k + j] - u;
		saved = 0.0;
		for (r = 0; r < j; r++)
		{
			ndu[j][r] = right[r + 1] + left[j - r];
			temp = ndu[r][j - 1] / ndu[j][r];
			ndu[r][j] = saved + right[r + 1] * temp;
			saved = left[j - r] * temp;
		}
		ndu[j][j] = saved;
	}
	for (int j = 0; j <= p; ++j)
		basisDu[0][j] = ndu[j][p];
	for (r = 0; r <= p; r++)
	{
		a[0][0] = 1.0;
		s1 = 0; s2 = 1;
		for (int i = 1; i <= n; ++i)
		{
			d = 0.0;
			rk = r - rn;
			pk = p - rn;
			if (r >= rn)
			{
				a[s2][0] = a[s1][0] / ndu[pk + 1][rk];
				d = a[s2][0] * ndu[rk][pk];
			}
			if (rk >= -1) j1 = 1;
			else j1 = (-rk);
			if (r - 1 <= pk)  j2 = rn - 1;
			else j2 = p - r;
			for (int j = j1; j <= j2; j++)
			{
				a[s2][j] = (a[s1][j] - a[s1][j - 1])
					/ ndu[pk + 1][rk + j];
				d += a[s2][j] * ndu[rk + j][pk];
			}
			if (r <= pk)
			{
				a[s2][rn] = -a[s1][rn - 1] / ndu[pk + 1][r];
				d += a[s2][rn] * ndu[r][pk];

			}
			basisDu[i][r] = d;
			s1 = s1 + s2; s2 = s1 - s2; s1 = s1 - s2;
		}
	}
	r = p;
	for (int i = 1; i <= n; i++)
	{
		for (int j = 0; j <= p; j++)
			basisDu[i][j] *= r;
		r *= (p - k);
	}
}

point3d NurbsLine::GetLinePoint(const double u)
{
	int m_id;
	vector<double> Nu;
	Nu.resize(_u_Degree + 1);
	point3d res(0, 0, 0);
	double SumW = 0;
	m_id = FindSpan(u, _u_Degree, _ControlPts->size(), *_u_Knots);
	BasisFuns(u, m_id, _u_Degree, *_u_Knots, Nu);

	for (int i = 0; i <= _u_Degree; i++)
	{
		res += Nu[i] * (*_ControlPts)[m_id - _u_Degree + i] * (*_ControlPts)[m_id - _u_Degree + i].w;
		SumW += Nu[i] * (*_ControlPts)[m_id - _u_Degree + i].w;
	}
	res /= SumW;
	return res;
}

void NurbsLine::CalLinePoint(const int Unum, vector<point3d>& linePts)
{
	linePts.clear();
	double du = 1.0 / Unum;
	for (int i = 0; i <= Unum; i++)         //采用等分节点矢量插值B样条曲线
	{
		double u = du * i;
		point3d t = GetLinePoint(u);
		linePts.push_back(t);
	}
}

void NurbsLine::PtsDerivsAu(const double u, const int n, vector<point4d>& Der)
{
	Der.clear();

	int span = 0, p = _u_Degree, num = _ControlPts->size(), tmp;
	int du = n < p ? n : p;
	//计算基函数
	vector<vector<double>> ndu;
	span = FindSpan(u, p, num, *_u_Knots);
	AllBasisFuns(u, span, p, *_u_Knots, ndu);

	varray<varray<point4d>> PKK;
	varray<point4d> PK;
	for (int i = 0; i < num; i++)
	{
		point4d pt;
		pt = (*_ControlPts)[i] * (*_ControlPts)[i].w;
		pt.w = (*_ControlPts)[i].w;
		PK.push_back(pt);
	}
	PKK.push_back(PK);

	for (int k = 1; k <= du; k++)
	{
		PK.clear();
		tmp = p - k + 1;
		for (int i = 0; i < num - k; i++)
		{
			point4d pt;
			pt = tmp * (PKK[k - 1][i + 1] - PKK[k - 1][i]) / ((*_u_Knots)[i + p + 1] - (*_u_Knots)[i + k]);
			pt.w = tmp * (PKK[k - 1][i + 1].w - PKK[k - 1][i].w) / ((*_u_Knots)[i + p + 1] - (*_u_Knots)[i + k]);
			PK.push_back(pt);
		}
		PKK.push_back(PK);
	}

	//计算导数
	for (int k = 0; k <= du; k++)
	{
		point4d pt;
		double sw = 0;
		for (int j = 0; j <= p - k; j++)
		{
			pt += ndu[j][p - k] * PKK[k][span - p + j];
			sw += ndu[j][p - k] * PKK[k][span - p + j].w;
		}
		pt.w = sw;
		Der.push_back(pt);
	}
}

void NurbsLine::PtsDerivs(const double u, const int n, vector<point4d>& Der)
{
	Der.clear();
	int du = n < _u_Degree ? n : _u_Degree;
	vector<point4d> DerAu;
	PtsDerivsAu(u, du, DerAu);
	for (int k = 0; k <= du; k++)
	{
		for (int i = 1; i <= k; i++)
			DerAu[k] = DerAu[k] - Global::Bin(k, i) * DerAu[i].w*Der[k - i];

		point4d pt;
		pt = DerAu[k] / DerAu[0].w;
		if (k == 0)
		{
			pt.w = DerAu[0].w;
		}
		Der.push_back(pt);
	}
}

void NurbsLine::CurveDerivs(const int n, const double step, vector<vector<point3d>>& Der)
{
	Der.clear();
	for (double u = 0; u <= 1; u += step)
	{
		vector<point4d> CK;
		PtsDerivs(u, n, CK);
		vector<point3d> CK3;
		for (int i = 0; i < CK.size(); ++i)
			CK3.push_back(CK[i]);
		Der.push_back(CK3);
	}
}

void NurbsLine::KnotInsert(const double u, const int r)
{
	if (r <= 0)return;

	vector<point4d> newbpt;
	vector<double> newKonts;
	int n = _ControlPts->size() - 1;//原始控制点数-1
	newKonts.resize(_u_Knots->size() + r);
	newbpt.resize(newKonts.size() - _u_Degree - 1);
	int a, m;
	m = n + _u_Degree + 1;
	a = FindSpan(u, _u_Degree, _ControlPts->size(), *_u_Knots);
	for (int j = 0; j <= a - _u_Degree; j++)
		newbpt[j] = (*_ControlPts)[j];
	for (int j = a; j <= n; j++)
		newbpt[j + r] = (*_ControlPts)[j];
	for (int j = 0; j <= a; j++)
		newKonts[j] = (*_u_Knots)[j];
	for (int j = a + _u_Degree + 1; j <= m; j++)
		newKonts[j + r] = (*_u_Knots)[j];
	int i = a + _u_Degree;
	int k = a + _u_Degree + r;
	while (u <= (*_u_Knots)[i] && i > a)
	{
		newbpt[k - _u_Degree - 1] = (*_ControlPts)[i - _u_Degree - 1];
		newKonts[k] = (*_u_Knots)[i];
		k--;
		i--;
	}
	newbpt[k - _u_Degree - 1] = newbpt[k - _u_Degree];
	//插入节点r次，逆序
	for (int j = r - 1; j >= 0; j--)
	{
		for (int L = 1; L <= _u_Degree; L++)
		{
			int ind = k - _u_Degree + L;
			double alfa = newKonts[k + L] - u, beta;
			if (alfa <= 0.0)
			{
				newbpt[ind - 1] = newbpt[ind];
			}
			else
			{
				alfa = alfa / (newKonts[k + L] - (*_u_Knots)[i - _u_Degree + L]);
				double nw = alfa * newbpt[ind - 1].w + (1.0f - alfa)*newbpt[ind].w;
				beta = alfa * newbpt[ind - 1].w / nw;
				newbpt[ind - 1] = beta * newbpt[ind - 1] + (1.0f - beta)*newbpt[ind];
				newbpt[ind - 1].w = nw;
			}
		}
		newKonts[k] = u;
		k--;
	}
	*_u_Knots = newKonts;
	*_ControlPts = newbpt;

}

int NurbsLine::KnotRemove(const double u, const int r)
{
	int sp = FindSpan(u, _u_Degree, (*_ControlPts).size(), (*_u_Knots));
	size_t n = (*_ControlPts).size() - 1;
	size_t m = n + _u_Degree + 1;
	size_t ord = _u_Degree + 1;

	int s = 0;
	for (int i = 0; i < (*_u_Knots).size(); ++i)
	{
		if (u == (*_u_Knots)[i])
			s++;
		if (u < (*_u_Knots)[i])
			break;
	}

	double wmin = 1, Pmax = 0;
	for (int i = 0; i < (*_ControlPts).size(); ++i)
	{
		if ((*_ControlPts)[i].w < wmin)
			wmin = (*_ControlPts)[i].w;
		double pd = (*_ControlPts)[i].Magnitude();
		if (pd > Pmax)
			Pmax = pd;
	}
	double TOL = 999;//wmin / (1 + Pmax);//容差

	int fout = (2 * sp - s - _u_Degree) / 2;
	int last = sp - s;
	int first = sp - _u_Degree;
	int t = 0;
	for (t = 0; t < r; ++t)
	{
		varray<point4d> temp;
		temp.resize(2 * _u_Degree + 1);
		int off = first - 1;
		temp[0] = (*_ControlPts)[off];
		temp[last + 1 - off] = (*_ControlPts)[last + 1];
		int i = first, j = last, ii = 1, jj = last - off, remflag = 0;
		while (j - i > t)
		{
			double alfi = (u - (*_u_Knots)[i]) / ((*_u_Knots)[i + ord + t] - (*_u_Knots)[i]);
			double alfj = (u - (*_u_Knots)[j - t]) / ((*_u_Knots)[j + ord] - (*_u_Knots)[j - t]);
			temp[ii] = ((*_ControlPts)[i] - (1.0 - alfi)*temp[ii - 1]) / alfi;
			temp[ii].w = ((*_ControlPts)[i].w - (1.0 - alfi)*temp[ii - 1].w) / alfi;
			temp[jj] = ((*_ControlPts)[j] - alfj * temp[jj + 1]) / (1.0 - alfj);
			temp[jj].w = ((*_ControlPts)[j].w - alfj * temp[jj + 1].w) / (1.0 - alfj);

			i++;
			ii++;
			j--;
			jj--;
		}

		if (j - i < t)
		{
			point4d dp = temp[ii - 1] - temp[jj + 1];
			dp.w = temp[ii - 1].w - temp[jj + 1].w;
			double ds = sqrt(dp.SquareMagnitude() + dp.w*dp.w);
			if (ds <= TOL)
				remflag = 1;
		}
		else
		{
			double alfi = (u - (*_u_Knots)[i]) / ((*_u_Knots)[i + ord + t] - (*_u_Knots)[i]);
			point4d dp = (*_ControlPts)[i] - (alfi*temp[ii + t + 1] + (1.0 - alfi)*temp[ii - 1]);
			dp.w = (*_ControlPts)[i].w - (alfi*temp[ii + t + 1].w + (1.0 - alfi)*temp[ii - 1].w);
			double ds = sqrt(dp.SquareMagnitude() + dp.w*dp.w);
			if (ds <= TOL)
				remflag = 1;
		}
		if (remflag == 0)
			break;
		else
		{
			i = first;
			j = last;
			while (j - i > t)
			{
				(*_ControlPts)[i] = temp[i - off];
				(*_ControlPts)[j] = temp[j - off];
				i++;
				j--;
			}
		}
		first--;
		last++;
	}

	if (t == 0)
		return 0;
	for (int k = sp + 1; k <= m; ++k)
		(*_u_Knots)[k - t] = (*_u_Knots)[k];
	int j = fout, i = j;
	for (int k = 1; k < t; ++k)
		if (k % 2 == 1)
			i++;
		else
			j--;
	for (int k = i + 1; k <= n; ++k)
	{
		(*_ControlPts)[j] = (*_ControlPts)[k];
		j++;
	}
	for (int k = 0; k < t; ++k)
	{
		(*_u_Knots).pop_back();
		(*_ControlPts).pop_back();
	}
	return t;
}

void NurbsLine::KnotsRefine(const vector<double>& u)
{
	if (u.size() <= 0)return;

	vector<point4d> newbpt;
	vector<double> newKonts;

	int n = (*_ControlPts).size() - 1;//原始控制点数-1
	int r = u.size() - 1;//插入的节点数-1
	newKonts.resize((*_u_Knots).size() + u.size());
	newbpt.resize(newKonts.size() - _u_Degree - 1);
	int a, b, m;
	m = n + _u_Degree + 1;
	a = FindSpan(u[0], _u_Degree, (*_ControlPts).size(), (*_u_Knots));
	b = FindSpan(u[r], _u_Degree, (*_ControlPts).size(), (*_u_Knots));
	b = b + 1;
	for (int j = 0; j <= a - _u_Degree; j++)
		newbpt[j] = (*_ControlPts)[j];
	for (int j = b - 1; j <= n; j++)
		newbpt[j + r + 1] = (*_ControlPts)[j];
	for (int j = 0; j <= a; j++)
		newKonts[j] = (*_u_Knots)[j];
	for (int j = b + _u_Degree; j <= m; j++)
		newKonts[j + r + 1] = (*_u_Knots)[j];

	int i = b + _u_Degree - 1;
	int k = b + _u_Degree + r;
	for (int j = r; j >= 0; j--)
	{
		while (u[j] <= (*_u_Knots)[i] && i > a)
		{
			newbpt[k - _u_Degree - 1] = (*_ControlPts)[i - _u_Degree - 1];
			newKonts[k] = (*_u_Knots)[i];
			k--;
			i--;
		}
		newbpt[k - _u_Degree - 1] = newbpt[k - _u_Degree];
		for (int L = 1; L <= _u_Degree; L++)
		{
			int ind = k - _u_Degree + L;
			double alfa = newKonts[k + L] - u[j], beta;
			if (alfa <= 0.0)
			{
				newbpt[ind - 1] = newbpt[ind];
			}
			else
			{
				alfa = alfa / (newKonts[k + L] - (*_u_Knots)[i - _u_Degree + L]);
				double nw = alfa * newbpt[ind - 1].w + (1.0f - alfa)*newbpt[ind].w;
				beta = alfa * newbpt[ind - 1].w / nw;
				newbpt[ind - 1] = beta * newbpt[ind - 1] + (1.0f - beta)*newbpt[ind];
				newbpt[ind - 1].w = nw;
			}
		}
		newKonts[k] = u[j];
		k--;
	}
	(*_u_Knots) = newKonts;
	(*_ControlPts) = newbpt;
}

void NurbsLine::DegreeElevate(const int degree)
{
	int t = degree - _u_Degree;
	if (t < 1)
		return;

	vector<double> Uh;//升阶后的节点矢量
	vector<point4d> Qw;//升阶后的控制点

	{
		int sz = 1;
		for (int i = 0; i < (*_u_Knots).size() - 1; ++i)
		{
			if ((*_u_Knots)[i] != (*_u_Knots)[i + 1])
				++sz;
		}
		Uh.resize((*_u_Knots).size() + 2 * t + (sz - 2)*t);
		Qw.resize(Uh.size() - _u_Degree - 1 - t);
	}

	int n = (*_ControlPts).size() - 1;
	int m = n + _u_Degree + 1;
	int ph = _u_Degree + t, ph2 = ph / 2;

	varray<varray<double>> bezalfs;
	bezalfs.resize(ph + 1);
	for (int i = 0; i < bezalfs.size(); ++i)
		bezalfs[i].resize(_u_Degree + 1);

	varray<double> alphas(_u_Degree - 1, 0);

	varray<point4d> bpts(_u_Degree + 1, point4d()), ebpts(ph + 1, point4d()), Nextbpts(_u_Degree - 1, point4d());

	//计算升阶系数
	bezalfs[0][0] = bezalfs[ph][_u_Degree] = 1.0;
	for (int i = 1; i <= ph2; ++i)
	{
		double inv = 1.0 / Global::Bin(ph, i);
		int mpi = _u_Degree < i ? _u_Degree : i;
		for (int j = 0 > i - t ? 0 : i - t; j <= mpi; ++j)
			bezalfs[i][j] = inv * Global::Bin(_u_Degree, j)*Global::Bin(t, i - j);
	}
	for (int i = ph2 + 1; i < ph; ++i)
	{
		int mpi = _u_Degree < i ? _u_Degree : i;
		for (int j = 0 > i - t ? 0 : i - t; j <= mpi; ++j)
			bezalfs[i][j] = bezalfs[ph - i][_u_Degree - j];
	}
	int mh = ph, kind = ph + 1;
	int r = -1;
	int a = _u_Degree, b = _u_Degree + 1, cind = 1;
	double ua = (*_u_Knots)[0];
	Qw[0] = (*_ControlPts)[0];

	for (int i = 0; i <= ph; ++i)
		Uh[i] = ua;
	for (int i = 0; i <= _u_Degree; ++i)
		bpts[i] = (*_ControlPts)[i];

	while (b < m)
	{
		int i = b;
		while (b < m && (*_u_Knots)[b] == (*_u_Knots)[b + 1])
			++b;
		int mul = b - i + 1;
		mh += mul + t;
		double ub = (*_u_Knots)[b];
		int oldr = r;
		r = _u_Degree - mul;
		int lbz = 0, rbz = 0;
		//插入节点r次
		if (oldr > 0)
			lbz = (oldr + 2) / 2;
		else
			lbz = 1;
		if (r > 0)
			rbz = ph - (r + 1) / 2;
		else
			rbz = ph;
		//插入节点获得贝塞尔曲线段
		if (r > 0)
		{
			double numer = ub - ua;
			varray<double> alfs;
			alfs.resize(_u_Degree - mul);
			for (int k = _u_Degree; k > mul; --k)
				alfs[k - mul - 1] = numer / ((*_u_Knots)[a + k] - ua);
			for (int j = 1; j <= r; ++j)
			{
				int save = r - j, s = mul + j;
				for (int k = _u_Degree; k >= s; --k)
				{
					bpts[k] = alfs[k - s] * bpts[k] + (1 - alfs[k - s])*bpts[k - 1];
					bpts[k].w = alfs[k - s] * bpts[k].w + (1 - alfs[k - s])*bpts[k - 1].w;
				}
				Nextbpts[save] = bpts[_u_Degree];
			}
		}

		//对贝塞尔曲线升阶
		for (i = lbz; i <= ph; ++i)
		{
			ebpts[i] = point4d(0, 0, 0, 0);
			int mpi = _u_Degree < i ? _u_Degree : i;
			for (int j = 0 > i - t ? 0 : i - t; j <= mpi; ++j)
			{
				ebpts[i] += bezalfs[i][j] * bpts[j];
				ebpts[i].w += bezalfs[i][j] * bpts[j].w;
			}
		}

		//消去节点u=(*_u_Knots)[a]
		if (oldr > 1)
		{
			int first = kind - 2, last = kind;
			double den = ub - ua;
			double bet = (ub - Uh[kind - 1]) / den;
			for (int tr = 1; tr < oldr; ++tr)
			{
				i = first;
				int j = last;
				int kj = j - kind + 1;
				while (j - i > tr)//计算消去后的控制点
				{
					if (i < cind)
					{
						double alf = (ub - Uh[i]) / (ua - Uh[i]);
						Qw[i] = alf * Qw[i] + (1 - alf)*Qw[i - 1];
						Qw[i].w = alf * Qw[i].w + (1 - alf)*Qw[i - 1].w;
					}
					if (j >= lbz)
					{
						if (j - tr <= kind - ph + oldr)
						{
							double gam = (ub - Uh[j - tr]) / den;
							ebpts[kj] = gam * ebpts[kj] + (1 - bet)*ebpts[kj + 1];
							ebpts[kj].w = gam * ebpts[kj].w + (1 - bet)*ebpts[kj + 1].w;
						}
						else
						{
							ebpts[kj] = bet * ebpts[kj] + (1 - bet)*ebpts[kj + 1];
							ebpts[kj].w = bet * ebpts[kj].w + (1 - bet)*ebpts[kj + 1].w;
						}
					}
					++i; --j; --kj;
				}
				--first; ++last;
			}
		}//消去节点u=(*_u_Knots)[a]结束
		 //载入节点ua
		if (a != _u_Degree)
			for (i = 0; i < ph - oldr; ++i)
			{
				Uh[kind] = ua;
				++kind;
			}
		//存入控制点
		for (int j = lbz; j <= rbz; ++j)
		{
			Qw[cind] = ebpts[j];
			++cind;
		}
		//为下次循环准备
		if (b < m)
		{
			for (int j = 0; j < r; ++j)
				bpts[j] = Nextbpts[j];
			for (int j = r; j <= _u_Degree; ++j)
				bpts[j] = (*_ControlPts)[b - _u_Degree + j];
			a = b; ++b; ua = ub;
		}
		else//尾端节点
			for (i = 0; i <= ph; ++i)
				Uh[kind + i] = ub;
	}
	//删除容器多余元素及空间
	int nh = mh - ph;//控制点数量
	Qw.erase(Qw.begin() + nh, Qw.end());

	int nk = nh + _u_Degree + t + 1;//节点矢量数量
	Uh.erase(Uh.begin() + nk, Uh.end());

	(*_ControlPts) = Qw;
	(*_u_Knots) = Uh;
	_u_Degree = (*_u_Knots).size() - (*_ControlPts).size() - 1;
}

void NurbsLine::Segmentation(const double u, vector<NurbsLine>& lines)
{
	lines.clear();
	if (u <= 0 || u >= 1)
		return;

	lines.resize(2);

	NurbsLine l0(*this);
	int insN = _u_Degree;
	int span = 0;

	if (Global::IsInSet(u, (*_u_Knots)))
	{
		span = FindSpan(u, _u_Degree, (*_ControlPts).size(), (*_u_Knots));
		for (int i = span; i > 0; --i)
		{
			if (u == (*_u_Knots)[i])
				--insN;
			else
				break;
		}
	}
	if (insN > 0)
		l0.KnotInsert(u, insN);
	lines[0]._u_Degree = l0._u_Degree;
	lines[1]._u_Degree = l0._u_Degree;
	//计算新节点矢量
	for (int i = 0; i < l0._u_Degree + 1; ++i)
	{
		lines[0]._u_Knots->push_back(0);
		lines[1]._u_Knots->push_back(0);
	}
	span = l0.FindSpan(u, l0._u_Degree, (*l0._ControlPts).size(), (*l0._u_Knots));
	for (int i = l0._u_Degree + 1; i <= span - l0._u_Degree; ++i)
	{
		double newu = (*l0._u_Knots)[i] / u;
		(*lines[0]._u_Knots).push_back(newu);
	}
	for (int i = span + 1; i < (*l0._ControlPts).size(); ++i)
	{
		double newu = 1 - (1 - (*l0._u_Knots)[i]) / (1 - u);
		(*lines[1]._u_Knots).push_back(newu);
	}
	for (int i = 0; i < l0._u_Degree + 1; ++i)
	{
		(*lines[0]._u_Knots).push_back(1);
		(*lines[1]._u_Knots).push_back(1);
	}
	//控制点
	int i = 0;
	for (i; i < (*lines[0]._u_Knots).size() - lines[0]._u_Degree - 1; ++i)
		(*lines[0]._ControlPts).push_back((*l0._ControlPts)[i]);
	for (--i; i < (*l0._ControlPts).size(); ++i)
		(*lines[1]._ControlPts).push_back((*l0._ControlPts)[i]);
}

/*NURBS曲线分解为BEZIER
	Qw：BEZIER曲线段控制点*/
void NurbsLine::Decompose(varray<varray<point4d>>& Qw)
{
	int n = _ControlPts->size() - 1;
	int m = n + _u_Degree + 1;
	int a = _u_Degree, b = _u_Degree + 1;
	int nb = 0;
	for (int i = 1; i < _u_Knots->size(); ++i)
		if ((*_u_Knots)[i] != (*_u_Knots)[i - 1])
			nb++;
	Qw.resize(nb);
	for (int i = 0; i < nb; ++i)
		Qw[i].resize(_u_Degree + 1);

	nb = 0;
	for (int i = 0; i <= _u_Degree; ++i)
		Qw[nb][i] = (*_ControlPts)[i];

	while (b < m)
	{
		int i = b;
		while (b < m && (*_u_Knots)[b + 1] == (*_u_Knots)[b])
			b++;
		int mult = b - i + 1;
		if (mult < _u_Degree)
		{
			double numer = (*_u_Knots)[b] - (*_u_Knots)[a];
			varray<double> alphas(_u_Degree - 1);
			for (int j = _u_Degree; j > mult; --j)
				alphas[j - mult - 1] = numer / ((*_u_Knots)[a + j] - (*_u_Knots)[a]);
			int r = _u_Degree - mult;
			for (int j = 1; j <= r; ++j)
			{
				int save = r - j;
				int s = mult + j;
				for (int k = _u_Degree; k >= s; --k)
				{
					double alpha = alphas[k - a];
					Qw[nb][k] = alpha * Qw[nb][k] + (1 - alpha)*Qw[nb][k - 1];
					Qw[nb][k].w = alpha * Qw[nb][k].w + (1 - alpha)*Qw[nb][k - 1].w;
				}
				if (b < m)
					Qw[nb + 1][save] = Qw[nb][_u_Degree];
			}
		}
		nb++;
		if (b < m)
		{
			for (int i = _u_Degree - mult; i <= _u_Degree; ++i)
				Qw[nb][i] = (*_ControlPts)[b - _u_Degree + i];
			a = b;
			b++;
		}
	}
}

void NurbsLine::CruveReverse()
{
	int knotsLth = _u_Knots->size();
	int cptLth = _ControlPts->size();
	shared_ptr<vector<double>> tempKnots;
	shared_ptr<vector<point4d>> tempPts;
	for (int i = 0; i < _u_Degree + 1; i++)
	{
		tempKnots->push_back((*_u_Knots)[i]);
	}
	int imax = knotsLth - 1 - 2 * _u_Degree - 1;  //m-2p-1

	if (imax > 0)
	{
		for (int i = imax; i >= 1; i--)
		{
			tempKnots->push_back((*_u_Knots)[_u_Degree + i] * (-1.0) + (*_u_Knots)[0] + (*_u_Knots)[knotsLth - 1]);
		}
	}
	for (int i = 0; i < _u_Degree + 1; i++)
	{
		tempKnots->push_back((*_u_Knots)[knotsLth - 1]);
	}
	for (int i = 0; i < cptLth; i++)
	{
		tempPts->push_back((*_ControlPts)[cptLth - 1 - i]);
	}

	if (tempKnots->size() != knotsLth)
	{
		std::cout << "\n节点矢量数量不正确" << std::endl;
	}
	if (tempPts->size() != cptLth)
	{
		std::cout << "\n控制点数量不正确" << std::endl;
	}
	_ControlPts = tempPts;
	_u_Knots = tempKnots;
	return;
}

NurbsSurface::NurbsSurface()
{
	_u_Degree = 1;
	_u_Knots = make_shared<vector<double>>();
	_u_Num = 2;
	_ControlPts = make_shared<vector<point4d>>();

	_v_Degree = 1;
	_v_Num = 2;
	_v_Knots = make_shared<vector<double>>();
}

NurbsSurface::NurbsSurface(int u_num, int v_num, int u_degree, int v_degree,
	vector<double> u_knots, vector<double> v_knots, vector<point4d> Cp)
	: _u_Num(u_num)
	, _u_Degree(u_degree) 
	, _v_Degree(v_degree)
	, _v_Num(v_num)
{
	_u_Knots = make_shared<vector<double>>();
	*_u_Knots = u_knots;

	_ControlPts = make_shared<vector<point4d>>();
	*_ControlPts = Cp;

	_v_Knots = make_shared<vector<double>>();
	*_v_Knots = u_knots;
}

NurbsSurface::~NurbsSurface()
{

}

int NurbsSurface::GetVDegree()const
{
	return _v_Degree;
}

shared_ptr<vector<double>> NurbsSurface::GetVKonts()const
{
	return _v_Knots;
}

int NurbsSurface::GetVNum()const
{
	return _v_Num;
}

point4d NurbsSurface::getControlPoint(int u, int v)
{
	return _ControlPts->at(u + v * _u_Num);
}

bool NurbsSurface::SetVDegree(const int v_degree)
{
	if (v_degree < 0)
		return false;
	_v_Degree = v_degree;
	return true;
}

bool NurbsSurface::SetVKonts(const vector<double> knots)
{
	if (knots.empty())
	{
		return true;
	}
	_v_Knots->clear();
	for (auto iter = knots.begin(); iter != knots.end(); ++iter)
	{
		_v_Knots->push_back(*iter);
	}
	return true;
}


bool NurbsSurface::SetVKonts(shared_ptr<vector<double>> vKnots)
{
	if (vKnots->empty())
	{
		return false;
	}
	_v_Knots = vKnots;
	return true;
}

bool NurbsSurface::SetVNum(const int v_num)
{
	if (v_num < 0)
		return false;
	_v_Num = v_num;
	return true;
}

int NurbsVol::GetVDegree()const
{
	return _v_Degree;
}

shared_ptr<vector<double>> NurbsVol::GetVKonts()const
{
	return _v_Knots;
}

int NurbsVol::GetVNum()const
{
	return _v_Num;
}

point4d NurbsVol::getControlPoint(int u, int v)
{
	return _ControlPts->at(u + v * _u_Num);
}

bool NurbsVol::SetVDegree(const int v_degree)
{
	if (v_degree < 0)
		return false;
	_v_Degree = v_degree;
	return true;
}

bool NurbsVol::SetVKonts(const vector<double> knots)
{
	if (knots.empty())
	{
		return true;
	}
	_v_Knots->clear();
	for (auto iter = knots.begin(); iter != knots.end(); ++iter)
	{
		_v_Knots->push_back(*iter);
	}
	return true;
}


bool NurbsVol::SetVKonts(shared_ptr<vector<double>> vKnots)
{
	if (vKnots->empty())
	{
		return false;
	}
	_v_Knots = vKnots;
	return true;
}

bool NurbsVol::SetVNum(const int v_num)
{
	if (v_num < 0)
		return false;
	_v_Num = v_num;
	return true;
}

NurbsSurface::NurbsSurface(const NurbsSurface& ns)
{
	_u_Degree = ns._u_Degree;
	_u_Knots = make_shared<vector<double>>();
	for (auto iter = ns._u_Knots->begin(); iter != ns._u_Knots->end(); ++iter)
		_u_Knots->push_back(*iter);
	_u_Num = ns._u_Num;

	_v_Degree = ns._v_Degree;
	_v_Knots = make_shared<vector<double>>();
	for (auto iter = ns._v_Knots->begin(); iter != ns._v_Knots->end(); ++iter)
		_v_Knots->push_back(*iter);
	_v_Num = ns._v_Num;

	_ControlPts = make_shared<vector<point4d>>();
	for (auto iter = ns._ControlPts->begin(); iter != ns._ControlPts->end(); ++iter)
		_ControlPts->push_back(*iter);

}

NurbsSurface& NurbsSurface::operator=(const NurbsSurface &ns)
{
	if (this == &ns)
	{
		return *this;
	}
	_u_Degree = ns._u_Degree;
	_u_Knots = make_shared<vector<double>>();
	for (auto iter = ns._u_Knots->begin(); iter != ns._u_Knots->end(); ++iter)
		_u_Knots->push_back(*iter);
	_u_Num = ns._u_Num;

	_v_Degree = ns._v_Degree;
	_v_Knots = make_shared<vector<double>>();
	for (auto iter = ns._v_Knots->begin(); iter != ns._v_Knots->end(); ++iter)
		_v_Knots->push_back(*iter);
	_v_Num = ns._v_Num;

	_ControlPts = make_shared<vector<point4d>>();
	for (auto iter = ns._ControlPts->begin(); iter != ns._ControlPts->end(); ++iter)
		_ControlPts->push_back(*iter);
	return *this;
}

point3d NurbsSurface::GetSurFacePoint(const double u, const double v)const
{
	int m_uid = 0, m_vid = 0;
	vector<double> Nu, Nv;
	Nu.resize(this->GetUDegree() + 1);
	Nv.resize(_v_Degree + 1);
	m_uid = FindSpan(u, this->GetUDegree(), this->GetUNum(), (*this->GetUKonts()));
	BasisFuns(u, m_uid, this->GetUDegree(), (*this->GetUKonts()), Nu);
	m_vid = FindSpan(v, _v_Degree, _v_Num, (*_v_Knots));
	BasisFuns(v, m_vid, _v_Degree, (*_v_Knots), Nv);

	point3d res(0, 0, 0);
	double pw = 0.0;
	for (int j = 0; j <= _v_Degree; j++)
	{
		for (int i = 0; i <= this->GetUDegree(); i++)
		{
			//point4d bpti = (*this->GetUKonts())[(m_uid - this->GetUDegree() + i)*_v_Num + m_vid - _v_Degree + j];
			point4d bpti = _ControlPts->at((m_uid - _u_Degree + i)*_v_Num + m_vid - _v_Degree + j);
			res += Nu[i] * Nv[j] * bpti*bpti.w;
			pw += Nu[i] * Nv[j] * bpti.w;
		}
	}
	res /= pw;
	return res;
}

threadParam NurbsSurface::CalQuads(const int Unum, const int Vnum)const
{
	varray<varray<point3d>> quads;
	varray<varray<point3d>> lines;
	quads.clear();
	lines.clear();
	varray<varray<point3d>> L_u;
	double du = 1.0 / Unum;
	double dv = 1.0 / Vnum;
	for (int i = 0; i <= Unum; i++)
	{
		double u = du * i;
		varray<point3d> line;
		for (int j = 0; j <= Vnum; j++)
		{
			double v = dv * j;
			point3d t = GetSurFacePoint(u, v);
			line.push_back(t);
		}
		L_u.push_back(line);
	}

	for (int i = 0; i < L_u.size() - 1; ++i)
	{
		for (int j = 0; j < L_u[i].size() - 1; ++j)
		{
			varray<point3d> quad;
			quad.push_back(L_u[i][j]);
			quad.push_back(L_u[i + 1][j]);
			quad.push_back(L_u[i + 1][j + 1]);
			quad.push_back(L_u[i][j + 1]);
			quads.push_back(quad);
		}
	}
	threadParam param = { quads,lines };
	return param;
}

void NurbsSurface::CoonsInterpolate(const varray<varray<point4d>>& EdgeCtrlPts)
{
	varray<varray<point4d>> coonsPatchCtrlpts;
	int upnum = EdgeCtrlPts[0].size();//u方向数量
	int vpnum = EdgeCtrlPts[1].size();//v方向数量
	point4d P00, P10, P01, P11;//角点
	P00 = EdgeCtrlPts[0][0];
	P10 = EdgeCtrlPts[0][upnum - 1];
	P01 = EdgeCtrlPts[2][0];
	P11 = EdgeCtrlPts[2][upnum - 1];
	for (int i = 0; i < upnum; ++i)
	{
		double tui = 1.0*i / (upnum - 1);
		varray<point4d> Ctrlpts_i;
		for (int j = 0; j < vpnum; ++j)
		{
			double tvj = 1.0*j / (vpnum - 1);
			point4d Pij = (1 - tui)*EdgeCtrlPts[1][j] + tui * EdgeCtrlPts[3][j]
				+ (1 - tvj)*EdgeCtrlPts[0][i] + tvj * EdgeCtrlPts[2][i]
				- (1 - tvj)*((1 - tui)*P00 + tui * P10)
				- tvj * ((1 - tui)*P01 + tui * P11);
			Pij.w = (1 - tui)*EdgeCtrlPts[1][j].w + tui * EdgeCtrlPts[3][j].w
				+ (1 - tvj)*EdgeCtrlPts[0][i].w + tvj * EdgeCtrlPts[2][i].w
				- (1 - tvj)*((1 - tui)*P00.w + tui * P10.w)
				- tvj * ((1 - tui)*P01.w + tui * P11.w);

			Ctrlpts_i.push_back(Pij);
		}
		coonsPatchCtrlpts.push_back(Ctrlpts_i);
	}
	varray<point4d> temp = Global::transVectorToVarray((*this->GetControlPointer()));
	Global::ReduceVarrayDim(coonsPatchCtrlpts, temp, true);
}

void NurbsSurface::SetSurface(const int uDegree, const int vDegree, const int uNum, const int vNum,
	const vector<double>& uKnots, const vector<double>& vKnots)
{
	this->SetUDegree(uDegree);
	SetVDegree(vDegree);
	this->SetUNum(uNum);
	this->SetVNum(vNum);
	this->SetUKonts(uKnots);
	this->SetVKonts(vKnots);
}

bool NurbsSurface::isTwoSurfaceSame(NurbsSurface sf)
{
	float minLentorence = minDisBetweenSurfce;//最小允许误差
	float len = GetDistanceBetweenTwoSurface(sf);
	if (len > minLentorence)
	{
		return false;
	}
	else
	{
		return true;
	}
}

float NurbsSurface::GetDistanceBetweenTwoSurface(NurbsSurface sf)
{
	float meanMinlen = 0;
	for (int i = 0; i < _v_Num; i++)
	{
		for (int j = 0; j < _u_Num; j++)
		{
			point3d vt = getControlPoint(j, i);
			float minlen = 1000000;//magicNumber
			int uminid = -1;
			int vminid = -1;
			for (int k = 0; k < sf._v_Num; k++)
			{
				for (int l = 0; l < sf._u_Num; l++)
				{
					if ((vt - sf.getControlPoint(l, k)).Magnitude() < minlen)
					{
						vminid = k;
						uminid = l;
						minlen = (vt - sf.getControlPoint(l, k)).Magnitude();
					}
				}
			}
			meanMinlen += minlen;
		}
	}
	meanMinlen /= (_v_Num*_u_Num);
	return meanMinlen;
}


void NurbsSurface::CoonsInterpolate(const varray<NurbsLine>& EdgeLines)
{
	if (EdgeLines.size() != 4)
		return;
	SetSurface(EdgeLines[0]._u_Degree, EdgeLines[1]._u_Degree, EdgeLines[0]._u_Num, EdgeLines[1]._u_Num,
		*EdgeLines[0].GetUKonts(), *EdgeLines[1].GetUKonts());
	varray<varray<point4d>> cpt;
	for (int i = 0; i < EdgeLines.size(); ++i)
		cpt.push_back(Global::transVectorToVarray(*(EdgeLines[i].GetControlPointer())));
	CoonsInterpolate(cpt);
}

void NurbsSurface::DegreeElevate(const int Udegree, const int Vdegree)
{
	int tu = Udegree - this->GetUDegree();
	int tv = Vdegree - this->GetVDegree();
	//U方向
	if (tu > 0)
	{
		varray<varray<point4d>> NewConPoint;
		NurbsLine line;

		for (int i = 0; i < this->GetVNum(); ++i)
		{
			line.SetUDegree(this->GetUDegree());
			line.SetUKonts(*this->GetUKonts());
			line.GetControlPointer()->clear();
			for (int j = 0; j < this->GetUNum(); ++j)
				line.GetControlPointer()->push_back((*this->GetControlPointer())
					[this->GetVNum()*j + i]);
			line.DegreeElevate(Udegree);
			NewConPoint.push_back(Global::transVectorToVarray(*line.GetControlPointer()));
		}
		this->SetUDegree(Udegree);
		this->SetUKonts(*line.GetUKonts());
		this->SetUNum(line.GetControlPointer()->size());
		this->GetControlPointer()->resize(this->GetUNum()*this->GetVNum());
		for (int i = 0; i < NewConPoint.size(); ++i)
			for (int j = 0; j < NewConPoint[i].size(); ++j)
				(*this->GetControlPointer())[this->GetVNum()*j + i] = NewConPoint[i][j];
	}
	//V方向
	if (tv > 0)
	{
		varray<varray<point4d>> NewConPoint;
		NurbsLine line;

		for (int i = 0; i < this->GetVNum(); ++i)
		{
			line.SetUDegree(this->GetVDegree());
			line.SetUKonts(*this->GetVKonts());
			line.GetControlPointer()->clear();
			for (int j = 0; j < this->GetVNum(); ++j)
				line.GetControlPointer()->push_back((*this->GetControlPointer())
					[this->GetVNum()*j + i]);
			line.DegreeElevate(Vdegree);
			NewConPoint.push_back(Global::transVectorToVarray(*line.GetControlPointer()));
		}
		this->SetVDegree(Vdegree);
		this->SetVKonts(*line.GetUKonts());
		this->SetVNum(line.GetControlPointer()->size());
		this->GetControlPointer()->clear();
		for (int i = 0; i < NewConPoint.size(); ++i)
			for (int j = 0; j < NewConPoint[i].size(); ++j)
				this->GetControlPointer()->push_back(NewConPoint[i][j]);
	}
}

void NurbsSurface::KnotsRefine(const varray<double>& Uknot, const varray<double>& Vknot)
{
	//U方向
	if (Uknot.size() > 0)
	{
		varray<varray<point4d>> NewConPoint;
		NurbsLine line;

		for (int i = 0; i < this->GetVNum(); ++i)
		{
			line.SetUDegree(this->GetUDegree());
			line.SetUKonts(*this->GetUKonts());
			line.GetControlPointer()->clear();
			for (int j = 0; j < this->GetUNum(); ++j)
				line.GetControlPointer()->push_back(this->GetControlPointer()->at(this->GetVNum()*j + i));
			line.KnotsRefine(*this->GetUKonts());
			NewConPoint.push_back(Global::transVectorToVarray(*line.GetControlPointer()));
		}
		this->SetUKonts(*line.GetUKonts());
		this->SetUNum(line.GetControlPointer()->size());
		line.GetControlPointer()->resize(this->GetVNum()*this->GetUNum());
		for (int i = 0; i < NewConPoint.size(); ++i)
			for (int j = 0; j < NewConPoint[i].size(); ++j)
				(*_ControlPts)[_v_Num*j + i] = NewConPoint[i][j];
	}
	//V方向
	if (Vknot.size() > 0)
	{
		varray<varray<point4d>> NewConPoint;
		NurbsLine line;

		for (int i = 0; i < _u_Num; ++i)
		{
			line._u_Degree = _v_Degree;
			line._u_Knots = _v_Knots;
			line._ControlPts->clear();
			for (int j = 0; j < _v_Num; ++j)
				line._ControlPts->push_back(_ControlPts->at(_v_Num*i + j));
			line.KnotsRefine(*_v_Knots);
			NewConPoint.push_back(Global::transVectorToVarray(*line._ControlPts));
		}
		_v_Knots = line._u_Knots;
		_v_Num = line._ControlPts->size();
		_ControlPts->clear();
		for (int i = 0; i < NewConPoint.size(); ++i)
			for (int j = 0; j < NewConPoint[i].size(); ++j)
				_ControlPts->push_back(NewConPoint[i][j]);
	}
}

inline int NurbsSurface::CtrlPtsIdx(const int uIdx, const int vIdx)
{
	return _v_Num * uIdx + vIdx;
}

void NurbsSurface::OrderCtrlPts()
{
	vector<point4d> ControlPts;
	std::swap(_u_Num, _v_Num);
	std::swap(_u_Degree, _v_Degree);
	std::swap(_u_Knots, _v_Knots);
	for (int i = 0; i < _v_Num; ++i)
	{
		for (int j = 0; j < _u_Num; ++j)
		{
			ControlPts.push_back(_ControlPts->at(CtrlPtsIdx(j, i)));
		}
	}
	_ControlPts->clear();
	for (int i = 0; i != ControlPts.size(); ++i)
	{
		_ControlPts->push_back(ControlPts[i]);
	}
}

void NurbsSurface::OrderCtrlPts(NurbsSurface& sf)
{
	sf._ControlPts->clear();
	for (int i = 0; i < _v_Num; ++i)
	{
		for (int j = 0; j < _u_Num; ++j)
		{
			sf._ControlPts->push_back(_ControlPts->at(CtrlPtsIdx(j, i)));
		}
	}
	std::swap(sf._u_Num, sf._v_Num);
	std::swap(sf._u_Degree, sf._v_Degree);
	std::swap(sf._u_Knots, sf._v_Knots);
}

void NurbsSurface::GetEdgeLines(varray<NurbsLine>& EdgeLines) {
	EdgeLines.clear();
	EdgeLines.resize(4);
	NurbsLine nl;
	int num_pt = this->_ControlPts->size();
	//第一条,U,V=0
	nl._u_Degree = this->_u_Degree;
	nl._u_Knots = this->_u_Knots;
	for (int i = 0; i < this->_u_Num; ++i) {
		shared_ptr<vector<point4d>> pp = this->_ControlPts;
		vector<point4d> pp1 = *pp;
		nl._ControlPts->push_back(pp1[i]);
	}
	EdgeLines[0] = nl;

	//第二条,U,V=1
	nl._ControlPts->clear();
	for (int i = num_pt - this->_u_Num; i < num_pt; ++i) {
		shared_ptr<vector<point4d>> pp = this->_ControlPts;
		vector<point4d> pp1 = *pp;
		nl._ControlPts->push_back(pp1[i]);
	}
	EdgeLines[1] = nl;

	//第三条,V,U=0
	nl._u_Degree = this->_v_Degree;
	nl._u_Knots = this->_v_Knots;
	nl._ControlPts->clear();
	for (int i = 0; i < this->_v_Num; ++i) {
		shared_ptr<vector<point4d>> pp = this->_ControlPts;
		vector<point4d> pp1 = *pp;
		nl._ControlPts->push_back(pp1[i*this->_u_Num]);
	}
	EdgeLines[2] = nl;

	//第四条,V,U=1
	nl._ControlPts->clear();
	for (int i = 0; i < this->_v_Num; ++i) {
		shared_ptr<vector<point4d>> pp = this->_ControlPts;
		vector<point4d> pp1 = *pp;
		nl._ControlPts->push_back(pp1[(i + 1)*this->_u_Num - 1]);
	}
	EdgeLines[3] = nl;
}

NurbsVol::NurbsVol()
{
	_u_Degree = 1;
	_u_Knots = make_shared<vector<double>>();
	_u_Num = 2;
	_ControlPts = make_shared<vector<point4d>>();

	_v_Degree = 1;
	_v_Num = 2;
	_v_Knots = make_shared<vector<double>>();

	_w_Degree = 1;
	_w_Num = 1;

	_w_Knots = make_shared<vector<double>>();
}

NurbsVol::~NurbsVol()
{
}

NurbsVol::NurbsVol(const NurbsVol & nvlos)
{
	_u_Degree = nvlos._u_Degree;
	_u_Knots = make_shared<vector<double>>();
	for (auto iter = nvlos._u_Knots->begin(); iter != nvlos._u_Knots->end(); ++iter)
		_u_Knots->push_back(*iter);
	_u_Num = nvlos._u_Num;

	_v_Degree = nvlos._v_Degree;
	_v_Knots = make_shared<vector<double>>();
	for (auto iter = nvlos._v_Knots->begin(); iter != nvlos._v_Knots->end(); ++iter)
		_v_Knots->push_back(*iter);
	_v_Num = nvlos._v_Num;

	_w_Degree = nvlos._w_Degree;
	_w_Knots = make_shared<vector<double>>();
	for (auto iter = nvlos._w_Knots->begin(); iter != nvlos._w_Knots->end(); ++iter)
		_w_Knots->push_back(*iter);
	_w_Num = nvlos._w_Num;

	_ControlPts = make_shared<vector<point4d>>();
	for (auto iter = nvlos._ControlPts->begin(); iter != nvlos._ControlPts->end(); ++iter)
		_ControlPts->push_back(*iter);
}

NurbsVol& NurbsVol::operator=(const NurbsVol& nvlos)
{
	if (this == &nvlos)
	{
		return *this;
	}
	_u_Degree = nvlos._u_Degree;
	_u_Knots = make_shared<vector<double>>();
	for (auto iter = nvlos._u_Knots->begin(); iter != nvlos._u_Knots->end(); ++iter)
		_u_Knots->push_back(*iter);
	_u_Num = nvlos._u_Num;

	_v_Degree = nvlos._v_Degree;
	_v_Knots = make_shared<vector<double>>();
	for (auto iter = nvlos._v_Knots->begin(); iter != nvlos._v_Knots->end(); ++iter)
		_v_Knots->push_back(*iter);
	_v_Num = nvlos._v_Num;

	_w_Degree = nvlos._w_Degree;
	_w_Knots = make_shared<vector<double>>();
	for (auto iter = nvlos._w_Knots->begin(); iter != nvlos._w_Knots->end(); ++iter)
		_w_Knots->push_back(*iter);
	_w_Num = nvlos._w_Num;

	_ControlPts = make_shared<vector<point4d>>();
	for (auto iter = nvlos._ControlPts->begin(); iter != nvlos._ControlPts->end(); ++iter)
		_ControlPts->push_back(*iter);
	return *this;
}

int NurbsVol::GetWDegree()const
{
	return _w_Degree;
}

shared_ptr<vector<double>> NurbsVol::GetWKonts()const
{
	return _w_Knots;
}

int NurbsVol::GetWNum()const
{
	return _w_Num;
}

point4d NurbsVol::getControlPoint(int u, int v, int w)
{
	return _ControlPts->at(u + v * _u_Num + w * (_u_Num*_v_Num));
}

bool NurbsVol::SetWDegree(const int w_degree)
{
	if (w_degree < 0)
	{
		return false;
	}
	_w_Degree = w_degree;
	return true;
}

bool NurbsVol::SetWKonts(const vector<double> knots)
{
	if (knots.empty())
	{
		return false;
	}
	_w_Knots->clear();
	for (auto iter = knots.begin(); iter != knots.end(); ++iter)
	{
		_w_Knots->push_back(*iter);
	}
	return true;
}

bool NurbsVol::SetWKonts(shared_ptr<vector<double>> w_knots)
{
	if (w_knots.use_count() == 0)
	{
		return false;
	}
	_w_Knots = w_knots;
	return true;
}

bool NurbsVol::SetWNum(const int w_num)
{
	if (w_num < 0)
	{
		return false;
	}
	_w_Num = w_num;
	return true;
}

void NurbsVol::SetVol(const int uDegree, const int vDegree, const int wDegree, const int uNum, const int vNum, const int wNum,
	const vector<double>& uKnots, const vector<double>& vKnots, const vector<double>& wKnots)
{
	_u_Degree = uDegree;
	_v_Degree = vDegree;
	_w_Degree = wDegree;
	_u_Num = uNum;
	_v_Num = vNum;
	_w_Num = wNum;
	for (int i = 0; i != uKnots.size(); ++i)
	{
		_u_Knots->push_back(uKnots[i]);
	}
	for (int i = 0; i != vKnots.size(); ++i)
	{
		_v_Knots->push_back(vKnots[i]);
	}
	for (int i = 0; i != wKnots.size(); ++i)
	{
		_w_Knots->push_back(vKnots[i]);
	}
}

point3d NurbsVol::GetVolPoint(const double u, const double v, const double w)const
{
	int r = FindSpan(u, _u_Degree, _u_Num, *_u_Knots);
	vector<double> uNknot;
	BasisFuns(u, r, _u_Degree, *_u_Knots, uNknot);

	int s = FindSpan(v, _v_Degree, _v_Num, *_v_Knots);
	vector<double> vNknot;
	BasisFuns(v, s, _v_Degree, *_v_Knots, vNknot);

	int t = FindSpan(w, _w_Degree, _w_Num, *_w_Knots);
	vector<double> wNknot;
	BasisFuns(w, t, _w_Degree, *_w_Knots, wNknot);

	point3d val;
	double pw = 0.0;
	int ii, jj, kk;
	ii = jj = kk = 0;
	for (int k = 0; k <= _w_Degree; k++)
	{
		kk = t - _w_Degree + k;
		for (int j = 0; j <= _u_Degree; j++)
		{
			jj = r - _u_Degree + j;
			for (int i = 0; i <= _v_Degree; i++)
			{
				ii = s - _v_Degree + i;
				point4d bpti = _ControlPts->at(kk*_u_Num*_v_Num + jj * _v_Num + ii);
				val += uNknot[j] * vNknot[i] * wNknot[k] * bpti*bpti.w;
				pw += uNknot[j] * vNknot[i] * wNknot[k] * bpti.w;
			}
		}
	}
	val /= pw;
	return val;
}

threadParamVOL NurbsVol::CalQuads(const int Unum, const int Vnum, const int Wnum,
	varray<varray<varray<point3d>>>& quads, varray<varray<varray<point3d>>>& lines) const
{
	quads.clear();
	lines.clear();
	int a = 0;
	int num[3]{ Unum ,Vnum ,Wnum };//曲面细分度三个方向的num
	for (int i = 0; i < 6; ++i)
	{
		varray<varray<point3d>> Q, L;
		CalIsoSurface(a, i % 2, num[(a + 1) % 3], num[(a + 2) % 3], Q, L);
		quads.push_back(Q);
		lines.push_back(L);
		if (i % 2)a += 1;
	}
	threadParamVOL vol = { quads,lines };

	return vol;
}

threadParamVOL NurbsVol::CalQuads2(const int Unum, const int Vnum, const int Wnum)const
{
	varray<varray<varray<point3d>>>quads, lines;
	quads.clear();
	lines.clear();
	int a = 0;
	int num[3]{ Unum ,Vnum ,Wnum };//曲面细分度三个方向的num
	for (int i = 0; i < 6; ++i)
	{
		varray<varray<point3d>> Q, L;
		CalIsoSurface(a, i % 2, num[(a + 1) % 3], num[(a + 2) % 3], Q, L);
		quads.push_back(Q);
		lines.push_back(L);
		if (i % 2)a += 1;
	}
	threadParamVOL vol = { quads,lines };

	return vol;
}

inline	int NurbsVol::CtrlPtsIdx(const int uIdx, const int vIdx, const int wIdx)
{
	return _u_Num * _v_Num*wIdx + _v_Num * uIdx + vIdx;
}

void NurbsVol::OrderCtrlPts()
{
	NurbsVol vol;
	OrderCtrlPts(vol);
	_ControlPts = vol._ControlPts;
}

void NurbsVol::OrderCtrlPts(NurbsVol& vol)
{
	vol._ControlPts->clear();
	for (int k = 0; k < _w_Num; ++k)
	{
		for (int i = 0; i < _v_Num; ++i)
		{
			for (int j = 0; j < _u_Num; ++j)
			{
				vol._ControlPts->push_back(_ControlPts->at(CtrlPtsIdx(j, i, k)));
			}
		}
	}
	std::swap(vol._u_Num, vol._v_Num);
	std::swap(vol._u_Degree, vol._v_Degree);
	std::swap((*vol._u_Knots), (*vol._v_Knots));
}

void NurbsVol::KnotsRefine(varray<double> Uknot, varray<double> Vknot, varray<double> Wknot)
{

	int u = Uknot.size(), v = Vknot.size(), w = Wknot.size();

	if (u > 0 && v > 0)
	{
		vector<vector<point4d>> NewCpts;
		NurbsSurface Suv;
		for (int k = 0; k < _w_Num; k++)
		{
			Suv.SetSurface(_u_Degree, _v_Degree, _u_Num, _v_Num, *_u_Knots, *_v_Knots);
			for (int j = 0; j < _u_Num; j++)
			{
				for (int i = 0; i < _v_Num; i++)
				{
					Suv._ControlPts->push_back(_ControlPts->at(i + j * _v_Num + k * _u_Num*_v_Num));
				}
			}

			Suv.KnotsRefine(Uknot, Vknot);
			NewCpts.push_back(*Suv._ControlPts);
			Suv._ControlPts->clear();
		}
		_ControlPts->clear();
		//_u_Knots = Suv._u_Knots; 共享智能指针，所以不能直接复制，代码结束后会消失
		for (int i = 0; i != Suv._u_Knots->size(); ++i)
		{
			_u_Knots->push_back(Suv._u_Knots->at(i));
		}
		//_v_Knots = Suv._v_Knots;
		for (int i = 0; i != Suv._u_Knots->size(); ++i)
		{
			_u_Knots->push_back(Suv._u_Knots->at(i));
		}
		_v_Num = Suv._v_Num;
		_u_Num = Suv._u_Num;
		for (int i = 0; i < NewCpts.size(); i++)
			for (int j = 0; j < NewCpts[i].size(); j++)
				_ControlPts->push_back(NewCpts[i][j]);
	}

	if (w > 0 && v > 0)
	{

		vector<vector<point4d>> NewCpts;
		NurbsSurface Suv;
		for (int k = 0; k < _u_Num; k++)
		{
			Suv.SetSurface(_w_Degree, _v_Degree, _w_Num, _v_Num, *_w_Knots, *_v_Knots);
			for (int j = 0; j < _w_Num; j++)
			{
				for (int i = 0; i < _v_Num; i++)
				{
					Suv._ControlPts->push_back(_ControlPts->at(i + j * _v_Num * _u_Num + k * _v_Num));
				}
			}
			Vknot.clear();
			Suv.KnotsRefine(Wknot, Vknot);
			NewCpts.push_back(*Suv._ControlPts);
			Suv._ControlPts->clear();
		}
		_ControlPts->clear();
		for (int i = 0; i != Suv._u_Knots->size(); ++i)
		{
			_w_Knots->push_back(Suv._u_Knots->at(i));
		}
		_w_Num = Suv._u_Num;

		for (int k = 0; k < _u_Num; k++) {
			for (int j = 0; j < NewCpts.size(); j++) {
				for (int i = 0; i < _v_Num; i++)
				{
					_ControlPts->push_back(NewCpts[j][i + k * _v_Num]);
				}
			}
		}

	}
}

void NurbsVol::DegreeElevate(const int Udegree, const int Vdegree, const int Wdegree)
{
	int tu = Udegree - _u_Degree;
	int tv = Vdegree - _u_Degree;
	int tw = Wdegree - _w_Degree;

	//U方向
	if (tu > 0)
	{
		vector<vector<vector<point4d>>> NewConPoint;
		NurbsLine line;
		for (int k = 0; k < _w_Num; ++k)
		{
			vector<vector<point4d>> NewConPointSF;
			for (int j = 0; j < _v_Num; ++j)
			{
				line._u_Degree = _u_Degree;
				for (auto iter = _u_Knots->begin(); iter != _u_Knots->end(); ++iter)
				{
					line._ControlPts->push_back(*iter);
				}
				line._ControlPts->clear();
				for (int i = 0; i < _u_Num; ++i)
					line._ControlPts->push_back(_ControlPts->at(CtrlPtsIdx(i, j, k)));
				line.DegreeElevate(Udegree);
				NewConPointSF.push_back(*line._ControlPts);
			}
			NewConPoint.push_back(NewConPointSF);//[w][v][u]
		}
		_u_Degree = Udegree;
		for (auto iter = line._u_Knots->begin(); iter != line._u_Knots->begin(); ++iter)
		{
			_u_Knots->push_back(*iter);
		}
		_u_Num = line._ControlPts->size();
		_ControlPts->resize(_u_Num*_v_Num*_w_Num);
		for (int i = 0; i < NewConPoint.size(); ++i)
		{
			for (int j = 0; j < NewConPoint[i].size(); ++j)
			{
				for (int k = 0; k < NewConPoint[i][j].size(); ++k)
				{
					int idx = CtrlPtsIdx(k, j, i);
					(*_ControlPts)[idx] = NewConPoint[i][j][k];
				}
			}
		}
	}
	//V方向
	if (tv > 0)
	{
		vector<vector<vector<point4d>>> NewConPoint;
		NurbsLine line;
		for (int k = 0; k < _w_Num; ++k)
		{
			vector<vector<point4d>> NewConPointSF;
			for (int i = 0; i < _u_Num; ++i)
			{
				line._u_Degree = _u_Degree;
				for (auto iter = _v_Knots->begin(); iter != _v_Knots->end(); ++iter)
				{
					line._u_Knots->push_back(*iter);
				}
				line._ControlPts->clear();
				for (int j = 0; j < _v_Num; ++j)
					line._ControlPts->push_back(_ControlPts->at(CtrlPtsIdx(i, j, k)));
				line.DegreeElevate(Vdegree);
				NewConPointSF.push_back(*line._ControlPts);
			}
			NewConPoint.push_back(NewConPointSF);//[w][u][v]
		}
		_u_Degree = Vdegree;
		for (auto iter = line._u_Knots->begin(); iter != line._u_Knots->begin(); ++iter)
		{
			_v_Knots->push_back(*iter);
		}
		_v_Num = line._ControlPts->size();
		_ControlPts->resize(_u_Num*_v_Num*_w_Num);
		for (int i = 0; i < NewConPoint.size(); ++i)
		{
			for (int j = 0; j < NewConPoint[i].size(); ++j)
			{
				for (int k = 0; k < NewConPoint[i][j].size(); ++k)
				{
					int idx = CtrlPtsIdx(j, k, i);
					(*_ControlPts)[idx] = NewConPoint[i][j][k];
				}
			}
		}
	}
	//W方向
	if (tw > 0)
	{
		vector<vector<vector<point4d>>> NewConPoint;
		NurbsLine line;
		for (int i = 0; i < _u_Num; ++i)
		{
			vector<vector<point4d>> NewConPointSF;
			for (int j = 0; j < _v_Num; ++j)
			{
				line._u_Degree = _w_Degree;
				for (auto iter = _w_Knots->begin(); iter != _w_Knots->end(); ++iter)
				{
					line._u_Knots->push_back(*iter);
				}
				line._ControlPts->clear();
				for (int k = 0; k < _w_Num; ++k)
					line._ControlPts->push_back(_ControlPts->at(CtrlPtsIdx(i, j, k)));
				line.DegreeElevate(Wdegree);
				NewConPointSF.push_back(*line._ControlPts);
			}
			NewConPoint.push_back(NewConPointSF);//[u][v][w]
		}
		_w_Degree = Wdegree;
		for (auto iter = line._u_Knots->begin(); iter != line._u_Knots->begin(); ++iter)
		{
			_w_Knots->push_back(*iter);
		}
		for (auto iter = line._u_Knots->begin(); iter != line._u_Knots->begin(); ++iter)
		{
			_v_Knots->push_back(*iter);
		}
		_w_Num = line._ControlPts->size();
		_ControlPts->resize(_u_Num*_v_Num*_w_Num);
		for (int i = 0; i < NewConPoint.size(); ++i)
		{
			for (int j = 0; j < NewConPoint[i].size(); ++j)
			{
				for (int k = 0; k < NewConPoint[i][j].size(); ++k)
				{
					int idx = CtrlPtsIdx(i, j, k);
					(*_ControlPts)[idx] = NewConPoint[i][j][k];
				}
			}
		}
	}
}

/*扫描生成Nurbs体模型
	pathT：扫描路径
	nurbsSF：起始截面
	K：截面实例数量,一般取路径控制点数量减1
	*/
void NurbsVol::CreateSweepNurbsVol(const NurbsLine& pathT, const NurbsSurface & nurbsSF, const int K)
{
	//u方向节点矢量，控制点数量
	_u_Degree = nurbsSF._u_Degree;
	_u_Num = nurbsSF._u_Num;
	for (auto iter = nurbsSF._u_Knots->begin(); iter != nurbsSF._u_Knots->end(); ++iter)
		_u_Knots->push_back(*iter);

	//v方向节点矢量，控制点数量
	_v_Degree = nurbsSF._v_Degree;
	_v_Num = nurbsSF._v_Num;
	for (auto iter = nurbsSF._v_Knots->begin(); iter != nurbsSF._v_Knots->end(); ++iter)
		_v_Knots->push_back(*iter);

	_w_Degree = pathT._u_Degree;
	_w_Knots->clear();
	_ControlPts->clear();

	varray<double> pos;
	InsLocation(pathT, K, pos);//同时完成W方向节点矢量
	varray<varray<point4d>> TranMat;
	LocalCoordinates(pathT, pos, TranMat);
	varray<varray<point4d>> nurbsSFctrlPts;
	varray<varray<varray<double>>> SFw;
	varray<varray<double>> sws;
	for (int i = 0; i < nurbsSF._u_Num; ++i)
	{
		varray<point4d> vctrl;
		varray<double> sw;
		for (int j = 0; j < nurbsSF._v_Num; ++j)
		{
			vctrl.push_back(nurbsSF._ControlPts->at(i*nurbsSF._v_Num + j));
			sw.push_back(nurbsSF._ControlPts->at(i*nurbsSF._v_Num + j).w);
		}
		nurbsSFctrlPts.push_back(vctrl);
		sws.push_back(sw);
	}
	varray<varray<varray<point4d>>> allNurbsSF;
	MatrixTran(nurbsSFctrlPts, TranMat, allNurbsSF);
	_w_Num = allNurbsSF.size();
	for (int i = 0; i < _w_Num; ++i)
	{
		SFw.push_back(sws);
	}
	SweepSurface(allNurbsSF, SFw, pos);
}

/*平移扫描生成Nurbs体模型
pathT：扫描路径
nurbsSF：起始截面*/
void NurbsVol::CreateTransSweepNurbsVol(const NurbsLine& pathT, const NurbsSurface & nurbsSF)
{
	//u方向节点矢量，控制点数量
	_u_Degree = nurbsSF._u_Degree;
	_u_Num = nurbsSF._u_Num;
	for (auto iter = nurbsSF._u_Knots->begin(); iter != nurbsSF._u_Knots->end(); ++iter)
		_u_Knots->push_back(*iter);

	//v方向节点矢量，控制点数量
	_v_Degree = nurbsSF._v_Degree;
	_v_Num = nurbsSF._v_Num;
	for (auto iter = nurbsSF._v_Knots->begin(); iter != nurbsSF._v_Knots->end(); ++iter)
		_v_Knots->push_back(*iter);

	_w_Degree = pathT._u_Degree;
	_w_Num = pathT._ControlPts->size();
	for (auto iter = pathT._u_Knots->begin(); iter != pathT._u_Knots->end(); ++iter)
		_w_Knots->push_back(*iter);

	_ControlPts->clear();
	for (int i = 0; i < pathT._ControlPts->size(); i++)
	{
		for (int j = 0; j < nurbsSF._ControlPts->size(); j++)
		{
			point4d ControlPoint;
			ControlPoint = nurbsSF._ControlPts->at(j) + pathT._ControlPts->at(i) - pathT._ControlPts->at(0);
			ControlPoint.w = nurbsSF._ControlPts->at(j).w*pathT._ControlPts->at(i).w;
			_ControlPts->push_back(ControlPoint);
		}
	}
}

//放样
//path:路径
//surfaces:放样截面
void NurbsVol::LoftingNurbsVol(const NurbsLine& path, const varray<NurbsSurface>& surfaces)
{
	_ControlPts->clear();
	_u_Knots->clear();
	_v_Knots->clear();
	_w_Knots->clear();

	_w_Degree = path._u_Degree;
	_w_Num = path._ControlPts->size();

	//最高次
	int maxUDegree, maxVDegree;
	MaxDegree(surfaces, maxUDegree, maxVDegree);
	//升阶
	varray<NurbsSurface> UnitSurfaces = surfaces;
	for (int i = 0; i < UnitSurfaces.size(); ++i)
		UnitSurfaces[i].DegreeElevate(maxUDegree, maxVDegree);
	//统一节点矢量
	varray<double> UKnots, VKnots;
	KnotsUnify(UnitSurfaces, UKnots, VKnots);
	for (int i = 0; i < UnitSurfaces.size(); ++i)
	{
		varray<double> diffUKnots, diffVKnots;

		Global::KnotsDiff(UKnots, Global::transVectorToVarray(*UnitSurfaces[i]._u_Knots), diffUKnots);
		Global::KnotsDiff(VKnots, Global::transVectorToVarray(*UnitSurfaces[i]._v_Knots), diffVKnots);
		//节点插入
		UnitSurfaces[i].KnotsRefine(diffUKnots, diffVKnots);
	}
	_u_Knots = UnitSurfaces[0]._u_Knots;
	_v_Knots = UnitSurfaces[0]._v_Knots;
	_u_Num = UnitSurfaces[0]._u_Num;
	_v_Num = UnitSurfaces[0]._v_Num;
	_u_Degree = UnitSurfaces[0]._u_Degree;
	_v_Degree = UnitSurfaces[0]._v_Degree;

	//varray<double> DIS, UK, U_W; double Dis = 0;
	//for (int i = 0; i<path._ControlPts.size() - 1; i++)
	//{
	//	point3d pt0, pt1; double d;
	//	pt0 = path._ControlPts[i];
	//	pt1 = path._ControlPts[i + 1];
	//	d = sqrt(pow((pt1.x - pt0.x), 2) + pow((pt1.y - pt0.y), 2) + pow((pt1.z - pt0.z), 2));
	//	DIS.push_back(d);
	//	Dis += d;
	//}
	//UK.push_back(0);
	//for (int i = 0; i<DIS.size() - 1; i++)
	//{
	//	UK.push_back(UK[UK.size() - 1] + DIS[i] / Dis);
	//}
	//UK.push_back(1);
	////节点向量
	//for (int i = 0; i < _w_Degree + 1; i++)
	//{
	//	U_W.push_back(0);
	//}
	//for (int j = 1; j<path._ControlPts.size() - _w_Degree; j++)
	//{
	//	double uw = 0;
	//	for (int i = j; i <= j + _w_Degree - 1; i++)
	//	{
	//		uw += UK[i];
	//	}
	//	U_W.push_back(uw / (j + _w_Degree - 1));
	//}
	//for (int i = 0; i<_w_Degree + 1; i++)
	//{
	//	U_W.push_back(1);
	//}
	//for (int i = 0; i<U_W.size(); i++)
	//	_w_Knots.push_back(U_W[i]);
	//for (int i = 0; i<surfaces.size(); i++)
	//	for (int j = 0; j < surfaces[i]._ControlPts.size(); j++)
	//		_ControlPts.push_back(surfaces[i]._ControlPts[j]);

	int K = (surfaces.size() - 1) > (path._ControlPts->size() - 1) ? (surfaces.size() - 1) : (path._ControlPts->size() - 1);
	varray<double> pos;
	InsLocation(path, K, pos);//同时完成W方向节点矢量
	varray<varray<point4d>> TranMat;
	LocalCoordinates(path, pos, TranMat);

	varray<varray<point4d>> nurbsSFctrlPts;
	for (int i = 0; i < UnitSurfaces[0]._u_Num; ++i)
	{
		varray<point4d> vctrl;
		for (int j = 0; j < UnitSurfaces[0]._v_Num; ++j)
			vctrl.push_back(UnitSurfaces[0]._ControlPts->at(i*UnitSurfaces[0]._v_Num + j));
		nurbsSFctrlPts.push_back(vctrl);
	}
	varray<varray<varray<point4d>>> allNurbsSF;
	MatrixTran(nurbsSFctrlPts, TranMat, allNurbsSF);
	_w_Num = allNurbsSF.size();
}

void NurbsVol::LoftingNurbsVol(const NurbsLine & path, const NurbsSurface & surfaces0, const NurbsSurface & surfaces1)
{
	_ControlPts->clear();
	_u_Knots->clear();
	_v_Knots->clear();
	_w_Knots->clear();

	_w_Degree = path._u_Degree;
	_w_Num = path._ControlPts->size();

	varray<NurbsSurface> surfaces;
	surfaces.push_back(surfaces0);
	surfaces.push_back(surfaces1);
	//最高次
	MaxDegree(surfaces, _u_Degree, _v_Degree);
	//升阶
	for (int i = 0; i < surfaces.size(); ++i)
		surfaces[i].DegreeElevate(_u_Degree, _v_Degree);
	//统一节点矢量
	varray<double> a = Global::transVectorToVarray(*_u_Knots);
	varray<double> b = Global::transVectorToVarray(*_v_Knots);
	KnotsUnify(surfaces, a, b);
	for (int i = 0; i < surfaces.size(); ++i)
	{
		varray<double> diffUKnots, diffVKnots;

		Global::KnotsDiff(Global::transVectorToVarray(*_u_Knots), Global::transVectorToVarray(*surfaces[i]._u_Knots), diffUKnots);
		Global::KnotsDiff(Global::transVectorToVarray(*_v_Knots), Global::transVectorToVarray(*surfaces[i]._v_Knots), diffVKnots);
		//节点插入
		surfaces[i].KnotsRefine(diffUKnots, diffVKnots);
	}
	_u_Num = surfaces[0]._u_Num;
	_v_Num = surfaces[0]._v_Num;

	int K = path._ControlPts->size() - 1;
	varray<double> pos;
	InsLocation(path, K, pos);//同时完成W方向节点矢量
	varray<varray<point4d>> TranMat01, TranMat10;
	LocalCoordinates(path, pos, TranMat01);

	//反向
	/*NurbsLine path10 = path;
	path10._ControlPts.reverse();
	for (int i = pos01.size() - 1; i >= 0; --i)
		pos10.push_back(1 - pos01[i]);
	LocalCoordinates(path10, pos10, TranMat10);*/
	TranMat10 = TranMat01;
	TranMat10.reverse();
	for (int i = 0; i < TranMat10.size(); ++i)
	{
		TranMat10[i][2] *= -1;
	}
	//原始权重
	varray<varray<double>> SFw0, SFw1;
	//控制点二维数组
	varray<varray<point4d>> nurbsSFctrlPts0, nurbsSFctrlPts1;

	for (int i = 0; i < surfaces[0]._u_Num; ++i)
	{
		varray<point4d> vctrl;
		varray<double> sw;
		for (int j = 0; j < surfaces[0]._v_Num; ++j)
		{
			vctrl.push_back(surfaces[0]._ControlPts->at(i*surfaces[0]._v_Num + j));
			sw.push_back(surfaces[0]._ControlPts->at(i*surfaces[0]._v_Num + j).w);
		}
		nurbsSFctrlPts0.push_back(vctrl);
		SFw0.push_back(sw);
	}
	for (int i = 0; i < surfaces[1]._u_Num; ++i)
	{
		varray<point4d> vctrl;
		varray<double> sw;
		for (int j = 0; j < surfaces[1]._v_Num; ++j)
		{
			vctrl.push_back(surfaces[1]._ControlPts->at(i*surfaces[1]._v_Num + j));
			sw.push_back(surfaces[0]._ControlPts->at(i*surfaces[0]._v_Num + j).w);
		}
		nurbsSFctrlPts1.push_back(vctrl);
		SFw1.push_back(sw);
	}
	//实例截面计算
	varray<varray<varray<point4d>>> allNurbsSF0, allNurbsSF1, allNurbsSF;
	varray<varray<varray<double>>>  SFw;
	MatrixTran(nurbsSFctrlPts0, TranMat01, allNurbsSF0);
	MatrixTran(nurbsSFctrlPts1, TranMat10, allNurbsSF1);

	//同一位置插值线性插值
	allNurbsSF.resize(pos.size());
	SFw.resize(pos.size());
	for (int i = 0; i < pos.size(); ++i)
	{
		allNurbsSF[i].resize(allNurbsSF0[i].size());
		SFw[i].resize(SFw0.size());
		for (int j = 0; j < allNurbsSF0[i].size(); ++j)
		{
			allNurbsSF[i][j].resize(allNurbsSF0[i][j].size());
			SFw[i][j].resize(SFw0[j].size());
			int sf0ij = allNurbsSF0[i][j].size();
			for (int k = 0; k < sf0ij; ++k)
			{
				point4d dp = allNurbsSF1[pos.size() - 1 - i][j][sf0ij - 1 - k] - allNurbsSF0[i][j][k];
				double dw = abs(allNurbsSF1[pos.size() - 1 - i][j][sf0ij - 1 - k].w - allNurbsSF0[i][j][k].w);
				point4d pts = allNurbsSF0[i][j][k]
					+ pos[i] * dp.Magnitude()*dp.Normalize();
				allNurbsSF[i][j][k] = pts;
				allNurbsSF[i][j][k].w = allNurbsSF0[i][j][k].w + pos[i] * dw;
				SFw[i][j][k] = SFw0[j][k] + pos[i] * abs(SFw1[j][k] - SFw0[j][k]);
			}
		}
	}
	_w_Num = allNurbsSF.size();
	SweepSurface(allNurbsSF, SFw, pos);
}



//计算等参面
//uvw:0=u,1=v,2=w
//t:参数
//num:以u-v-w顺序
//L:等参面四边形面片集
void NurbsVol::CalIsoSurface(const int uvw, const double t, const int num1, const int num2,
	varray<varray<point3d>>& quads, varray<varray<point3d>>& lines)const
{
	quads.clear();
	lines.clear();
	varray<varray<point3d>> L;
	double u[3]{ 0,0,0 };
	u[uvw] = t;
	int a = (uvw + 1) % 3;
	int b = (uvw + 2) % 3;

	double da = 1.0 / num1;
	double db = 1.0 / num2;
	for (int i = 0; i <= num1; i++)
	{
		u[a] = da * i;
		varray<point3d> line;
		for (int j = 0; j <= num2; j++)
		{
			u[b] = db * j;
			point3d t = GetVolPoint(u[0], u[1], u[2]);
			line.push_back(t);
		}
		L.push_back(line);
	}

	for (int i = 0; i < L.size() - 1; ++i)
	{
		for (int j = 0; j < L[i].size() - 1; ++j)
		{
			varray<point3d> quad;
			quad.push_back(L[i][j]);
			quad.push_back(L[i + 1][j]);
			quad.push_back(L[i + 1][j + 1]);
			quad.push_back(L[i][j + 1]);
			quads.push_back(quad);
		}
	}

	lines.push_back(L[0]);
	lines.push_back(L[L.size() - 1]);
	varray<point3d> l0, l1;
	for (int i = 0; i < L.size(); ++i)
	{
		l0.push_back(L[i][0]);
		l1.push_back(L[i][L[i].size() - 1]);
	}
	lines.push_back(l0);
	lines.push_back(l1);
}

/*扫描时的截面实例位置
pathT：扫描路径
K：截面实例数量
pos：截面实例位置
NewKnots：新节点矢量*/
void NurbsVol::InsLocation(const NurbsLine & pathT, int K, varray<double>& pos)
{
	pos.clear();
	int q = pathT._u_Degree;
	int ktv = pathT._u_Knots->size();
	int nsect = K + 1;
	double vsum;

	for (auto iter = pathT._u_Knots->begin(); iter != pathT._u_Knots->end(); ++iter)
		_w_Knots->push_back(*iter);

	if (ktv <= nsect + q)  //细化节点矢量
	{
		int m = nsect + q - ktv + 1;
		for (int i = 0; i < m; ++i)
		{
			//最长节点区间
			auto idx = _w_Knots->begin();
			double maxlen = 0;
			for (auto it = _w_Knots->begin() + pathT._u_Degree;
				it != _w_Knots->begin() + (_w_Knots->size() - pathT._u_Degree - 2); ++it)
			{
				double dk = *(it + 1) - *it;
				if (dk >= maxlen)
				{
					maxlen = dk;
					idx = it;
				}
			}
			//插入中点
			double mid = (*idx + *(idx + 1)) / 2;
			_w_Knots->insert(idx + 1, mid);
		}
	}
	else if (ktv > nsect + q + 1) //增加实例
		nsect = ktv - q - 1;

	//实例位置
	pos.push_back(_w_Knots->at(0));
	for (int k = 1; k < nsect - 1; k++)
	{
		vsum = 0;
		for (int l = k + 1; l < k + q + 1; l++)
		{
			vsum += _w_Knots->at(l);
		}
		pos.push_back(vsum / q);
	}
	pos.push_back(_w_Knots->at(_w_Knots->size() - 1));
}

/*计算局部坐标系
pathT：扫描路径
pos：截面实例位置
TranMat：pos处的局部坐标系*/
void NurbsVol::LocalCoordinates(const NurbsLine & pathT, const varray<double> & pos,
	varray<varray<point4d>> & TranMat)
{
	TranMat.clear();
	int q = pathT._u_Degree;
	//计算弗朗内特标
	//导矢
	vector<vector<point4d>> DersBasis;
	DersBasis.clear();
	for (int i = 0; i < pos.size(); i++)
	{
		vector<point4d> NuDerBasi;
		//转换const->non_const
		NurbsLine & path = const_cast<NurbsLine &>(pathT);
		path.PtsDerivs(pos[i], 2, NuDerBasi);
		if (NuDerBasi.size() < 3 || NuDerBasi[2] == point4d())//不存在2阶导
		{
			point4d vDer = NuDerBasi[1].Normalize().Abs();
			int midx = 0;
			for (int j = 1; j < 3; ++j)
			{
				if (vDer[j] > vDer[midx])
				{
					vDer[midx] = 0;
					midx = j;
				}
				else
					vDer[j] = 0;
			}
			if (vDer.Angle(NuDerBasi[1]) == 0 || vDer.Angle(NuDerBasi[1]) == 180)
			{
				vDer[midx] = 0;
				vDer[(midx + 1) % 3] = 1;
			}
			point4d vt = NuDerBasi[1].Cross(vDer).Normalize();
			if (NuDerBasi.size() < 3)
				NuDerBasi.push_back(vt);
			else
				NuDerBasi[2] = vt;
		}
		DersBasis.push_back(NuDerBasi);
	}
	for (int i = 0; i < pos.size(); i++)
	{
		point4d BV = (DersBasis[i][1].Cross(DersBasis[i][2])).Normalize();
		point4d NV = (BV.Cross(DersBasis[i][1])).Normalize();
		DersBasis[i][1] = DersBasis[i][1].Normalize();
		varray<point4d> Stor;
		Stor.push_back(BV);
		Stor.push_back(NV);
		Stor.push_back(DersBasis[i][1]);
		Stor.push_back(DersBasis[i][0]);
		TranMat.push_back(Stor);
	}
}

/*沿扫描路径对截面进行变换
nurbsSF：截面控制点
TranMat：变换矩阵
allNurbsSF：得到的所有截面实例控制点*/
void NurbsVol::MatrixTran(const varray<varray<point4d>> & nurbsSF, const varray<varray<point4d>>& TranMat,
	varray<varray<varray<point4d>>>& allNurbsSF)
{
	using Eigen::Matrix4d;
	using Eigen::Vector4d;

	allNurbsSF.clear();
	varray<point3d> OriCoordinates;
	for (int i = 0; i < TranMat[0].size(); ++i)
		OriCoordinates.push_back(TranMat[0][i]);

	for (int i = 0; i < TranMat.size(); i++)  //变换矩阵
	{
		Matrix4d mat;
		varray<point3d> newCoordinates;
		newCoordinates.push_back(TranMat[i][0]);
		newCoordinates.push_back(TranMat[i][1]);
		newCoordinates.push_back(TranMat[i][2]);
		newCoordinates.push_back(TranMat[i][3]);

		Global::CalCoordinatesTransMat(OriCoordinates, newCoordinates, mat);

		varray<varray<point4d>> SectPoints;
		for (int j = 0; j < nurbsSF.size(); j++)
		{
			varray<point4d> SectPoint;
			for (int k = 0; k < nurbsSF[j].size(); k++)
			{
				Vector4d bas = Global::P4dToV4d(nurbsSF[j][k]);
				bas[3] = 1;
				Vector4d torm = mat * bas;
				point4d pt = Global::V4dToP4d(torm) * TranMat[i][3].w;
				pt.w = TranMat[i][3].w;
				SectPoint.push_back(pt);
			}
			SectPoints.push_back(SectPoint);
		}
		allNurbsSF.push_back(SectPoints);
	}
}

//扫描生成Nurbs体
void NurbsVol::SweepSurface(const varray<varray<varray<point4d>>>& allNurbsSF,
	const varray<varray<varray<double>>>& SFw, const varray<double> &pos)
{
	using Eigen::MatrixXd;

	//解线性方程组，插值控制点坐标    
	int p = _w_Degree;
	MatrixXd Nbase(allNurbsSF.size(), allNurbsSF.size());
	MatrixXd NpointX(allNurbsSF.size(), _u_Num*_v_Num);
	MatrixXd NpointY(allNurbsSF.size(), _u_Num*_v_Num);
	MatrixXd NpointZ(allNurbsSF.size(), _u_Num*_v_Num);
	MatrixXd NpointW(allNurbsSF.size(), _u_Num*_v_Num);
	MatrixXd ConTropointX(allNurbsSF.size(), _u_Num*_v_Num);
	MatrixXd ConTropointY(allNurbsSF.size(), _u_Num*_v_Num);
	MatrixXd ConTropointZ(allNurbsSF.size(), _u_Num*_v_Num);
	MatrixXd ConTropointW(allNurbsSF.size(), _u_Num*_v_Num);
	for (int i = 0; i < pos.size(); i++)  //截面个数
	{
		//计算基函数
		vector<vector<double>> ndu;
		int span = FindSpan(pos[i], p, allNurbsSF.size(), *_w_Knots);
		AllBasisFuns(pos[i], span, p, *_w_Knots, ndu);
		for (int l = 0; l < span - p; l++)
		{
			Nbase(i, l) = 0;
		}
		for (int m = 0; m <= p; m++)
		{
			double npm = ndu[m][p];
			Nbase(i, span - p + m) = npm;
		}
		for (int n = span + 1; n < allNurbsSF.size(); n++)
		{
			Nbase(i, n) = 0;
		}
	}

	for (int i = 0; i < allNurbsSF.size(); i++) //第i个截面    
	{
		for (int j = 0; j < allNurbsSF[i].size(); j++)         //UV方向
		{
			for (int k = 0; k < allNurbsSF[i][j].size(); k++)
			{
				NpointX(i, j*allNurbsSF[i][j].size() + k) = allNurbsSF[i][j][k].x;          //体上的点
				NpointY(i, j*allNurbsSF[i][j].size() + k) = allNurbsSF[i][j][k].y;          //体上的点
				NpointZ(i, j*allNurbsSF[i][j].size() + k) = allNurbsSF[i][j][k].z;          //体上的点
				NpointW(i, j*allNurbsSF[i][j].size() + k) = allNurbsSF[i][j][k].w;          //体上的点
			}
		}
	}
	MatrixXd nBinv = Nbase.inverse();
	ConTropointX = Nbase.inverse()*NpointX;
	ConTropointY = Nbase.inverse()*NpointY;
	ConTropointZ = Nbase.inverse()*NpointZ;
	ConTropointW = Nbase.inverse()*NpointW;

	for (int i = 0; i < allNurbsSF.size(); i++) //第i个截面,w方向    
	{
		for (int j = 0; j < allNurbsSF[i].size(); j++) //u方向
		{
			for (int k = 0; k < allNurbsSF[i][j].size(); k++) //v方向
			{
				point4d control;
				control.w = ConTropointW(i, j*allNurbsSF[i][j].size() + k);
				control.x = ConTropointX(i, j*allNurbsSF[i][j].size() + k) / control.w;
				control.y = ConTropointY(i, j*allNurbsSF[i][j].size() + k) / control.w;
				control.z = ConTropointZ(i, j*allNurbsSF[i][j].size() + k) / control.w;
				control.w *= SFw[i][j][k];
				_ControlPts->push_back(control);
			}
		}
	}
}

//取最高次
void NurbsVol::MaxDegree(const varray<NurbsSurface>& surfaces, int& uDegree, int& vDegree)
{
	uDegree = 0;
	vDegree = 0;
	for (int i = 0; i < surfaces.size(); ++i)
	{
		uDegree = surfaces[i]._u_Degree > uDegree ? surfaces[i]._u_Degree : uDegree;
		vDegree = surfaces[i]._v_Degree > vDegree ? surfaces[i]._v_Degree : vDegree;
	}
}

//节点矢量并集
void NurbsVol::KnotsUnify(const varray<NurbsSurface>& surfaces, varray<double>& NewUKnots, varray<double>& NewVKnots)
{
	NewUKnots.clear();
	NewVKnots.clear();
	varray<double> temUKnots, temVKnots;
	for (int i = 0; i < surfaces.size(); ++i)
	{
		Global::KnotUnify(Global::transVectorToVarray(*surfaces[i]._u_Knots), temUKnots, NewUKnots);
		temUKnots = NewUKnots;
		Global::KnotUnify(Global::transVectorToVarray(*surfaces[i]._v_Knots), temVKnots, NewVKnots);
		temVKnots = NewVKnots;
	}
}

float NurbsVol::GetJacobianValue(float u, float v, float w)
{
	float nip, njq, nkr;
	varray<double> Nips, Njqs, Nkrs, Nip1s, Njq1s, Nkr1s;
	int i, j, k, i0, j0, k0, i1, j1, k1, i2, j2, k2;
	varray<double> u1Knots = Global::transVectorToVarray(*this->_u_Knots);
	u1Knots.pop_front();
	u1Knots.pop_back();
	varray<double> v1Knots = Global::transVectorToVarray(*this->_v_Knots);
	v1Knots.pop_front();
	v1Knots.pop_back();
	varray<double> w1Knots = Global::transVectorToVarray(*this->_w_Knots);
	w1Knots.pop_front();
	w1Knots.pop_back();
	for (i = 0; i < _u_Num; i++)
	{
		nip = OneBasisFun(_u_Degree, _u_Knots->size() - 1, Global::transVectorToVarray(*this->_u_Knots), i, u);
		Nips.push_back(nip);
	}
	for (i = 0; i < _u_Num - 1; i++)
	{
		nip = OneBasisFun(_u_Degree - 1, u1Knots.size() - 1, u1Knots, i, u);
		Nip1s.push_back(nip);
	}
	for (j = 0; j < _v_Num; j++)
	{
		njq = OneBasisFun(_v_Degree, _v_Knots->size() - 1, Global::transVectorToVarray(*this->_v_Knots), j, v);
		Njqs.push_back(njq);
	}
	for (j = 0; j < _v_Num - 1; j++)
	{
		njq = OneBasisFun(_v_Degree - 1, v1Knots.size() - 1, v1Knots, j, v);
		Njq1s.push_back(njq);
	}
	for (k = 0; k < _w_Num; k++)
	{
		nkr = OneBasisFun(_w_Degree, _w_Knots->size() - 1, Global::transVectorToVarray(*this->_w_Knots), k, w);
		Nkrs.push_back(nkr);
	}
	for (k = 0; k < _w_Num - 1; k++)
	{
		nkr = OneBasisFun(_w_Degree - 1, w1Knots.size() - 1, w1Knots, k, w);
		Nkr1s.push_back(nkr);
	}

	float nip2, nip1, nip0, njq2, njq1, njq0, nkr2, nkr1, nkr0;
	Vec3 vi, vj, vk;
	float Jac = 0;
	float PM;
	float degreeProduct = _u_Degree * _v_Degree*_w_Degree;
	float pm0, pm1, pm2;
	Vec3 crossVjvkproduct;
	for (k2 = 0; k2 < _w_Num - 1; k2++)
	{
		nkr2 = Nkr1s.at(k2);
		if (abs(nkr2) < 1.0e-6)
			break;
		pm2 = nkr2;
		for (j2 = 0; j2 < _v_Num; j2++)
		{
			njq2 = Njqs.at(j2);
			if (abs(njq2) < 1.0e-6)
				break;
			pm2 *= njq2;
			for (i2 = 0; i2 < _u_Num; i2++)
			{
				nip2 = Nips.at(i2);
				if (abs(nip2) < 1.0e-6)
					break;
				point3d tempPoint = getControlPoint(i2, j2, k2 + 1) - getControlPoint(i2, j2, k2);
				vk.x = tempPoint.x;
				vk.y = tempPoint.y;
				vk.z = tempPoint.z;
				pm2 *= nip2;
				for (k1 = 0; k1 < _w_Num; k1++)
				{
					nkr1 = Nkrs.at(k1);
					if (abs(nkr1) < 1.0e-6)
						break;
					pm1 = nkr1;
					for (j1 = 0; j1 < _v_Num - 1; j1++)
					{
						njq1 = Njq1s.at(j1);
						if (abs(njq1) < 1.0e-6)
							break;
						pm1 *= njq1;
						for (i1 = 0; i1 < _u_Num; i1++)
						{
							nip1 = Nips.at(i1);
							if (abs(nip1) < 1.0e-6)
								break;
							pm1 *= nip1;
							//vj = GetControlPoint(i1, j1 + 1, k1) - GetControlPoint(i1, j1, k1);
							point3d tempPoint = getControlPoint(i1, j1 + 1, k1) - getControlPoint(i1, j1, k1);
							vj.x = tempPoint.x;
							vj.y = tempPoint.y;
							vj.z = tempPoint.z;
							crossVjvkproduct = vj.CrossVecX(vk);
							pm1 = nip1 * njq1*nkr1;
							for (k0 = 0; k0 < _w_Num; k0++)
							{
								nkr0 = Nkrs.at(k0);
								if (abs(nkr0) < 1.0e-6)
									break;
								pm0 = nkr0;
								for (j0 = 0; j0 < _v_Num; j0++)
								{
									njq0 = Njqs.at(j0);
									if (abs(njq0) < 1.0e-6)
										break;
									pm0 *= njq0;
									for (i0 = 0; i0 < _u_Num - 1; i0++)
									{
										nip0 = Nip1s.at(i0);
										if (abs(nip0) < 1.0e-6)
											break;
										pm0 *= nip0;
										//vi = GetControlPoint(i0 + 1, j0, k0) - GetControlPoint(i0, j0, k0);
										point3d tempPoint = getControlPoint(i0 + 1, j0, k0) - getControlPoint(i0, j0, k0);
										vi.x = tempPoint.x;
										vi.y = tempPoint.y;
										vi.z = tempPoint.z;
										PM = pm0 * pm1*pm2;
										PM *= degreeProduct;
										PM *= vi.dotmultiple(crossVjvkproduct);
										Jac += PM;
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return Jac;
}

float NurbsVol::GetFirstPartialDirivatives(float u, float v, float w, Vec3& pu, Vec3& pv, Vec3& pw)
{
	int i, j, k;
	Vec3 tup, tvp, twp, deltu, deltv, deltw;
	float nip, njq, nkr;
	float uvwscale;
	varray<double> Nips, Njqs, Nkrs, Nip1s, Njq1s, Nkr1s;
	varray<double> u1Knots = Global::transVectorToVarray(*this->_u_Knots); u1Knots.pop_front(); u1Knots.pop_back();
	varray<double> v1Knots = Global::transVectorToVarray(*this->_v_Knots); v1Knots.pop_front(); v1Knots.pop_back();
	varray<double> w1Knots = Global::transVectorToVarray(*this->_w_Knots); w1Knots.pop_front(); w1Knots.pop_back();
	for (i = 0; i < _u_Num; i++)
	{
		nip = OneBasisFun(_u_Degree, _u_Knots->size() - 1, Global::transVectorToVarray(*this->_u_Knots), i, u);
		Nips.push_back(nip);
	}
	for (i = 0; i < _u_Num - 1; i++)
	{
		nip = OneBasisFun(_u_Degree - 1, u1Knots.size() - 1, u1Knots, i, u);
		Nip1s.push_back(nip);
	}
	for (j = 0; j < _v_Num; j++)
	{
		njq = OneBasisFun(_v_Degree, _v_Knots->size() - 1, Global::transVectorToVarray(*this->_v_Knots), j, v);
		Njqs.push_back(njq);
	}
	for (j = 0; j < _v_Num - 1; j++)
	{
		njq = OneBasisFun(_v_Degree - 1, v1Knots.size() - 1, v1Knots, j, v);
		Njq1s.push_back(njq);
	}
	for (k = 0; k < _w_Num; k++)
	{
		nkr = OneBasisFun(_w_Degree, _w_Knots->size() - 1, Global::transVectorToVarray(*this->_w_Knots), k, w);
		Nkrs.push_back(nkr);
	}
	for (k = 0; k < _w_Num - 1; k++)
	{
		nkr = OneBasisFun(_w_Degree - 1, w1Knots.size() - 1, w1Knots, k, w);
		Nkr1s.push_back(nkr);
	}

	tup = Vec3(0.f, 0.f, 0.f);
	for (k = 0; k < _w_Num; k++)
	{
		nkr = Nkrs.at(k);
		if (abs(nkr) < 1.0e-6)
			break;
		for (j = 0; j < _v_Num; j++)
		{
			njq = Njqs.at(j);
			if (abs(njq) < 1.0e-6)
				break;
			for (i = 0; i < _u_Num - 1; i++)
			{
				nip = Nip1s.at(i);
				if (abs(nip) < 1.0e-6)
					break;
				point3d tempPoint = getControlPoint(i + 1, j, k) - getControlPoint(i, j, k);
				deltu.x = tempPoint.x;
				deltu.y = tempPoint.y;
				deltu.z = tempPoint.z;
				uvwscale = (float)_u_Degree / (_u_Knots->at(i + _u_Degree + 1) - _u_Knots->at(i + 1));
				deltu = deltu * uvwscale;
				tup += nip * njq*nkr*deltu;
			}
		}
	}

	tvp = Vec3(0.f, 0.f, 0.f);
	for (k = 0; k < _w_Num; k++)
	{
		nkr = Nkrs.at(k);
		if (abs(nkr) < 1.0e-6)
			continue;
		for (j = 0; j < _v_Num - 1; j++)
		{
			njq = Njq1s.at(j);
			if (abs(njq) < 1.0e-6)
				continue;
			uvwscale = (float)_v_Degree / (_v_Knots->at(j + _v_Degree + 1) - _v_Knots->at(j + 1));
			for (i = 0; i < _u_Num; i++)
			{
				nip = Nips.at(i);
				if (abs(nip) < 1.0e-6)
					continue;
				point3d tempPoint = getControlPoint(i, j + 1, k) - getControlPoint(i, j, k);
				deltv.x = tempPoint.x;
				deltv.y = tempPoint.y;
				deltv.z = tempPoint.z;
				//deltv = (GetControlPoint(i, j + 1, k) - GetControlPoint(i, j, k));
				deltv = deltv * uvwscale;
				tvp += nip * njq*nkr*deltv;
			}
		}
	}

	twp = Vec3(0.f, 0.f, 0.f);
	for (k = 0; k < _w_Num - 1; k++)
	{
		nkr = Nkr1s.at(k);
		if (abs(nkr) < 1.0e-6)
			continue;
		uvwscale = (float)_w_Degree / (_w_Knots->at(k + _w_Degree + 1) - _w_Knots->at(k + 1));
		for (j = 0; j < _v_Num; j++)
		{
			njq = Njqs.at(j);
			if (abs(njq) < 1.0e-6)
				continue;
			for (i = 0; i < _u_Num; i++)
			{
				nip = Nips.at(i);
				if (abs(nip) < 1.0e-6)
					continue;
				//deltw = (GetControlPoint(i, j, k + 1) - GetControlPoint(i, j, k));
				point3d tempPoint = getControlPoint(i, j, k + 1) - getControlPoint(i, j, k);
				deltw.x = tempPoint.x;
				deltw.y = tempPoint.y;
				deltw.z = tempPoint.z;
				deltw = deltw * uvwscale;
				twp += nip * njq*nkr*deltw;
			}
		}
	}
	pu = tup;
	pv = tvp;
	pw = twp;
	float KK = tup.x*tup.x + tup.y*tup.y + tup.z*tup.z;
	KK += tvp.x*tvp.x + tvp.y*tvp.y + tvp.z*tvp.z;
	KK += twp.x*twp.x + twp.y*twp.y + twp.z*twp.z;
	return KK;
}

void   NurbsVol::TestMeshQuanity(int segmentNum)
{
	m_minorJocbianNum = 0;
	m_isovalAdd = 0;
	m_minJocbian = 1.0e6;
	m_maxJocbian = -1.0e6;

	float u, v, w;
	u = 0.f; v = 0.f; w = 0.f;
	int i, j, k;
	varray<float> Jacs;
	varray<float> Isoparaval;
	Vec3 pu, pv, pw;
	float steplength = 1.f / segmentNum;
	for (k = 0; k <= segmentNum; k++)
	{
		w = steplength * k;
		for (j = 0; j <= segmentNum; j++)
		{
			v = steplength * j;
			for (i = 0; i <= segmentNum; i++)
			{
				u = steplength * i;
				Jacs.push_back(GetJacobianValue(u, v, w));
				Isoparaval.push_back(GetFirstPartialDirivatives(u, v, w, pu, pv, pw));
			}
		}
	}

	for (i = 0; i < Jacs.size(); i++)
	{
		if (Jacs.at(i) < 0)
			m_minorJocbianNum++;
		if (Jacs.at(i) > m_maxJocbian)
			m_maxJocbian = Jacs.at(i);
		else if (Jacs.at(i) < m_minJocbian)
			m_minJocbian = Jacs.at(i);
	}

	for (i = 0; i < Isoparaval.size(); i++)
	{
		m_isovalAdd += Isoparaval.at(i);
	}
}

NurbsSurface NurbsVol::GetSingleSurfaces(int dir) const
{
	NurbsSurface sf;
	assert(dir >= 1 && dir <= 6);
	assert(this->_ControlPts->size());
	//if (dir < 1 || dir>6)
	//	return {};
	if (dir == 1) {
		//UV面,W=0
		sf._u_Degree = this->_u_Degree;
		sf._v_Degree = this->_v_Degree;
		sf._u_Num = this->_u_Num;
		sf._v_Num = this->_v_Num;
		sf._u_Knots = this->_u_Knots;
		sf._v_Knots = this->_v_Knots;
		for (int i = 0; i < this->_v_Num; i++) {
			for (int j = 0; j < this->_u_Num; j++) {
				int ii = i * this->_u_Num + j;
				point4d pp = (*this->_ControlPts)[ii];
				sf._ControlPts->push_back(pp);
				/*sf._ControlPts->push_back(this._ControlPts[i*this->_u_Num + j]);*/
			}
		}
	}

	else if (dir == 2) {
		//UW面，V=0
		sf._u_Degree = this->_u_Degree;
		sf._v_Degree = this->_w_Degree;
		sf._u_Num = this->_u_Num;
		sf._v_Num = this->_w_Num;
		sf._u_Knots = this->_u_Knots;
		sf._v_Knots = this->_w_Knots;
		for (int i = 0; i < this->_w_Num; i++) {
			for (int j = 0; j < this->_u_Num; j++) {
				int ii = i * (this->_u_Num*this->_v_Num) + j;
				point4d pp = (*this->_ControlPts)[ii];
				sf._ControlPts->push_back(pp);
			}
		}
	}

	else if (dir == 5) {
		//VW面，U=0
		sf._u_Degree = this->_v_Degree;
		sf._v_Degree = this->_w_Degree;
		sf._u_Num = this->_v_Num;
		sf._v_Num = this->_w_Num;
		sf._u_Knots = this->_v_Knots;
		sf._v_Knots = this->_w_Knots;
		for (int i = 0; i < this->_w_Num; i++) {
			for (int j = 0; j < this->_v_Num; j++) {
				int ii = i * (this->_u_Num*this->_v_Num) + j * _u_Num;
				point4d pp = (*this->_ControlPts)[ii];
				sf._ControlPts->push_back(pp);
				/*sf._ControlPts->push_back(this->_ControlPts[i*(this->_u_Num*this->_v_Num) + j * _u_Num]);*/
			}
		}
	}

	else if (dir == 6) {
		//UV面,W=1
		sf._u_Degree = this->_u_Degree;
		sf._v_Degree = this->_v_Degree;
		sf._u_Num = this->_u_Num;
		sf._v_Num = this->_v_Num;
		sf._u_Knots = this->_u_Knots;
		sf._v_Knots = this->_v_Knots;
		int pos = this->_u_Num*this->_v_Num*(this->_w_Num - 1);		//第一个点的位置
		for (int i = 0; i < this->_v_Num; i++) {
			for (int j = 0; j < this->_u_Num; j++) {
				int ii = pos + i * this->_u_Num + j;
				point4d pp = (*this->_ControlPts)[ii];
				sf._ControlPts->push_back(pp);
				/*sf._ControlPts->push_back(this->_ControlPts[pos + i * this->_u_Num + j]);*/
			}
		}
	}

	return sf;
}

