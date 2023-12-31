//  Real-time Simulation Library
//	Vec2tor.h: Vec2クラス宣言部
//
//  	Copyright(C) 2001,  T. Ishii
//////////////////////////////////////////////////////////////////////
#pragma once
#include "XIOStream.h"
#include "assert.h"

#define PI 3.1415926535

namespace base {

	//二次元ベクトルクラス
	class  Vec2 {
	public:
		double x, y;
	public:
		//構築
		Vec2();
		//構築
		Vec2(double _x, double _y);

		//代入
		void Set(double _x, double _y);

		//ゼロにする
		void Zero();

		//要素
		const double& e(int i) const { return (&x)[i]; }

		//要素
		double& e(int i) { return (&x)[i]; }

		//要素
		double& operator[] (int i) { return (&x)[i]; }

		//要素
		const double& operator[] (int i) const { return (&x)[i]; }

		//四則演算
		Vec2& operator+=(const Vec2 &a);
		Vec2& operator-=(const Vec2 &a);
		Vec2& operator*=(const double &a);
		Vec2& operator/=(const double &a);

		//大きさの２乗
		double SquareMagnitude() const;

		//大きさ
		double Magnitude() const;

		//正規化
		Vec2 Normalize();

		// 成分の最大を返す
		double Min() const;

		// 成分の最小を返す
		double Max() const;

		//絶対値
		Vec2 Abs();

		double Dot(const Vec2& a)
		{
			return a.x*x + a.y*y;
		}

		double Vec2::Angle(const Vec2 p)
		{
			double a, b, c;
			a = this->Dot(p);
			b = this->Magnitude();
			c = p.Magnitude();
			if (b == 0 || c == 0)
				return -1;

			double d = a / (b*c);
			double res = acos(d) * 180 / PI;
			return res;
		}

		//ベクトル同士の演算
		friend Vec2 operator+(const Vec2& a, const Vec2& b);
		friend Vec2 operator-(const Vec2& a, const Vec2& b);
		friend Vec2 operator/(const Vec2& a, const Vec2& b);
		// 符合反転
		friend Vec2 operator-(const Vec2& a);

		//スカラーとの演算
		friend Vec2 operator*(const Vec2& a, const double& b);
		friend Vec2 operator/(const Vec2& a, const double& b);
		friend Vec2 operator*(const double& b, const Vec2& a);

		//内積(外積)
		friend double    Dot(const Vec2& a, const Vec2& b);
		//デターミナント
		friend double    Det(const Vec2& a, const Vec2& b);
		//外積
		friend Vec2     CrossVecX(const Vec2& a);

		//最大化
		friend Vec2 Minimize(const Vec2& v1, const Vec2& v2);
		//最小化
		friend Vec2 Maximize(const Vec2& v1, const Vec2& v2);

		//比較
		friend int operator == (const Vec2& v1, const Vec2& v2);

		//ストリーム出力
		friend BaseOStream& operator << (BaseOStream& os, const Vec2& a);
		//ストリーム入力
		friend BaseIStream& operator >> (BaseIStream& is, Vec2& a);
	};
	//////////////////////////////////////////////////////////////////////
	//実装部
	//  高速化のためにすべてインライン展開される
	class  Vec3 {
	public:
		double x, y, z;
	public:
		//構築
		Vec3();

		//構築
		Vec3(double _x, double _y, double _z);

		//代入
		void Set(double _x, double _y, double _z);

		//ゼロ
		void Zero();

		//要素
		const double& e(int i) const {
			return (&x)[i];
		}
		//要素
		double& e(int i) {
			return (&x)[i];
		}

		//要素
		double& operator[] (int i) { return (&x)[i]; }

		//要素
		const double& operator[] (int i) const { return (&x)[i]; }
		//四則演算
		Vec3& operator+=(const Vec3 &a);
		Vec3& operator+=(const double &a)
		{
			x += a; y += a; z += a;
			return  (*this);
		};
		Vec3& operator-=(const Vec3 &a);
		Vec3& operator*=(const double &a);
		Vec3& operator/=(const double &a);
		bool operator < (const Vec3& p) const
		{
			double erf = 0.001;					//ｵ羞ﾄﾎ�ｲ�ﾅﾐｶﾏ耙ﾖｵ
			if (x - p.x <0 && abs(x - p.x)>erf)
				return true;
			else if (abs(x - p.x) <= erf && y - p.y < 0 && abs(y - p.y)>erf)
				return true;
			else if (abs(x - p.x) <= erf && abs(y - p.y) <= erf && z - p.z <0 && abs(z - p.z)>erf)
				return true;
			return false;
		}
		//大きさの二乗
		double SquareMagnitude() const;
		//大きさ
		double Magnitude() const;

		//正規化
		Vec3 Normalize() const;

		// 成分の最大を返す
		double Min() const;
		// 成分の最小を返す
		double Max() const;

		bool Equals(const Vec3& v, double off) const;

		//絶対値
		Vec3 Abs();

		//ベクトル同士の演算
		friend Vec3 operator+(const Vec3& a, const Vec3& b);
		friend Vec3 operator-(const Vec3& a, const Vec3& b);
		friend Vec3 operator/(const Vec3& a, const Vec3& b);
		// 符合反転
		friend Vec3 operator-(const Vec3& a);

		//スカラーとの演算
		friend Vec3 operator*(const Vec3& a, const double& b);
		friend Vec3 operator/(const Vec3& a, const double& b);
		friend Vec3 operator*(const double& b, const Vec3& a);

		//内積
		friend double  Dot(const Vec3& a, const Vec3& b);

		double Dot(Vec3& a) {
			return a.x*x + a.y*y + a.z*z;
		}
		//外積
		friend Vec3 CrossVecX(const Vec3& a, const Vec3& b);

		//最大化
		friend Vec3 Minimize(const Vec3& v1, const Vec3& v2);
		//最小化
		friend Vec3 Maximize(const Vec3& v1, const Vec3& v2);

		//比較
		friend int operator == (const Vec3& v1, const Vec3& v2);

		//ストリーム出力
		friend BaseOStream& operator << (BaseOStream& os, const Vec3& a);
		//ストリーム入力
		friend BaseIStream& operator >> (BaseIStream& is, Vec3& a);

		// Added by D. XXX - Apr. 15, 2006
		Vec3	multiple(Vec3 point) const;
		double	dotmultiple(Vec3 point) const;
		double	GetLength();
		void	unit();
		Vec3 Trans(double dx, double dy, double dz);
		Vec3 RotateX(double ang);
		Vec3 RotateY(double ang);
		Vec3 RotateZ(double ang);
		Vec3 RotAxis(Vec3 vect, double ang);
		Vec3 CrossVecX(const Vec3& b);
		Vec3 Cross(const Vec3& b) { return CrossVecX(b); }
		// Compute the global transfomation angle of a vector
		void ComputeGlbAng(double *xrot, double *yrot);

		// scale the vector to 0~1;
		void Scale_to_01()
		{
			assert(!(x == 0 && y == 0 && z == 0));
			double mag = 1 / (double)sqrt(x*x + y * y + z * z);
			x *= mag;
			y *= mag;
			z *= mag;
			*this *= (double)0.5;
			*this += (double)0.5;
		}
		//ﾁｽﾏ�ﾁｿｼﾐｽﾇ｣ｨｶﾈ｣ｩ
		double Vec3::Angle(const Vec3 p)
		{
			double a, b, c;
			a = this->dotmultiple(p);
			b = this->Magnitude();
			c = p.Magnitude();
			if (b == 0 || c == 0)
				return -1;

			double d = a / (b*c);
			if (d > 1)d = 1;
			else if (d < -1)d = -1;
			double res = acos(d) * 180 / PI;
			return res;
		}

	};

	class  IntVec3
	{
	public:
		int x;
		int y;
		int z;
		IntVec3() :x(-1), y(-1), z(-1) {};
		IntVec3(int _x, int _y, int _z) :x(_x), y(_y), z(_z) {};
	};
	typedef Vec3 Vec;

	class  IntVec4
	{
	public:
		int x;
		int y;
		int z;
		int w;
		IntVec4() :x(-1), y(-1), z(-1), w(1) {};
		IntVec4(int _x, int _y, int _z, int _w = 1) :x(_x), y(_y), z(_z), w(_w) {};
	};
	//////////////////////////////////////////////////////////////////////
	//実装部
	//  高速化のためにすべてインライン展開される
	class  Vec4 : public Vec3 {
	public:
		double w;
	public:
		//構築
		Vec4();
		Vec4(double _x, double _y, double _z, double _w = 1.0);
		Vec4(const Vec4& v) {
			x = v.x;
			y = v.y;
			z = v.z;
			w = v.w;
		}
		Vec4(Vec3& v) {
			x = v.x;
			y = v.y;
			z = v.z;
			w = 1;
		}
		Vec4(const Vec3& v) {
			x = v.x;
			y = v.y;
			z = v.z;
			w = 1;
		}

		//代入
		void Set(double _x, double _y, double _z, double _w = 1.0);
		void Zero();

		//要素
		const double& e(int i) const {
			return (&x)[i];
		}
		double& e(int i) {
			return (&x)[i];
		}

		//要素
		double& operator[] (int i) { return (&x)[i]; }

		//要素
		const double& operator[] (int i) const { return (&x)[i]; }

		//四則演算
		Vec4& operator+=(const Vec4 &a);
		Vec4& operator-=(const Vec4 &a);
		Vec4& operator*=(const double &a);
		Vec4& operator/=(const double &a);
		bool operator < (const Vec4& p) const
		{
			double erf = 0.000001;					//ｵ羞ﾄﾎ�ｲ�ﾅﾐｶﾏ耙ﾖｵ
			if (x - p.x <0 && abs(x - p.x)>erf)
				return true;
			else if (abs(x - p.x) <= erf && y - p.y < 0 && abs(y - p.y)>erf)
				return true;
			else if (abs(x - p.x) <= erf && abs(y - p.y) <= erf && z - p.z <0 && abs(z - p.z)>erf)
				return true;
			return false;
		}
		//正規化
		Vec4 Normalize() const;
		//絶対値
		Vec4 Abs()  const;

		//ベクトル同士の演算
		friend Vec4 operator+(const Vec4& a, const Vec4& b);
		friend Vec4 operator-(const Vec4& a, const Vec4& b);
		friend Vec4 operator/(const Vec4& a, const Vec4& b);
		// 符合反転
		friend Vec4 operator-(const Vec4& a);

		//スカラーとの演算
		friend Vec4 operator*(const Vec4& a, const double& b);
		friend Vec4 operator/(const Vec4& a, const double& b);
		friend Vec4 operator*(const double& b, const Vec4& a);


		//最大化
		friend Vec4 Minimize(const Vec4& v1, const Vec4& v2);
		//最小化
		friend Vec4 Maximize(const Vec4& v1, const Vec4& v2);

		//比較
		friend int operator == (const Vec4& v1, const Vec4& v2);
		//ｲ豕ﾋ
		friend Vec4 CrossVecX(const Vec4& a, const Vec4& b);

		Vec4	multiple(Vec4 point) const;
		Vec4 Trans(double dx, double dy, double dz);
		Vec4 RotateX(double ang);
		Vec4 RotateY(double ang);
		Vec4 RotateZ(double ang);
		Vec4 RotAxis(Vec4 vect, double ang);
		Vec4 CrossVecX(const Vec4& b);



	};

	//////////////////////////////////////////////////////////////////////
	//実装部
	//  高速化のためにすべてインライン展開される
	inline
		Vec2::Vec2()
	{
		x = 0; y = 0;
	}

	inline
		Vec2::Vec2(double _x, double _y)
	{
		x = _x; y = _y;
	}

	inline void
		Vec2::Set(double _x, double _y)
	{
		x = _x; y = _y;
	}

	inline void
		Vec2::Zero()
	{
		x = 0.0; y = 0.0;
	}

	//代入演算
	inline Vec2&
		Vec2::operator+=(const Vec2 &a)
	{
		x = x + a.x; y = y + a.y;
		return  (*this);
	}

	inline Vec2&
		Vec2::operator-=(const Vec2 &a)
	{
		x = x - a.x; y = y - a.y;
		return  (*this);
	}

	inline Vec2&
		Vec2::operator*=(const double &a)
	{
		x = a * x; y = a * y;
		return  (*this);
	}

	inline Vec2&
		Vec2::operator/=(const double &a)
	{
		x = x / a; y = y / a;
		return  (*this);
	}

	//大きさ
	inline double
		Vec2::SquareMagnitude() const
	{
		return x * x + y * y;
	}

	inline double
		Vec2::Magnitude()  const
	{
		return sqrt(SquareMagnitude());
	}


	//正規化
	inline Vec2
		Vec2::Normalize()
	{
		double len = Magnitude();
		Vec2 nor;
		nor.x = x / len;
		nor.y = y / len;
		return nor;
	}

	// 成分の最大,最小を返す
	inline double
		Vec2::Min()  const
	{
		double ret = x;
		if (y < ret) ret = y;
		return ret;
	}

	inline double
		Vec2::Max()  const
	{
		double ret = x;
		if (ret < y) ret = y;
		return ret;
	}

	//正規化
	inline Vec2
		Vec2::Abs()
	{
		Vec2 r;
		r.x = fabs(x);
		r.y = fabs(y);
		return r;
	}

	///////////////////////////////////////////////////////////////////
	//フレンド関数

	//ベクトル演算
	inline Vec2
		operator+(const Vec2& a, const Vec2& b)
	{
		Vec2 ret = a;
		ret += b;

		return ret;
	}

	inline Vec2
		operator-(const Vec2& a, const Vec2& b)
	{
		Vec2 ret = a;
		ret -= b;

		return ret;
	}

	inline Vec2
		operator/(const Vec2& a, const Vec2& b)
	{
		Vec2 ret = a;
		ret -= b;

		return Vec2(a.x / b.x, a.y / b.y);
	}

	// 符合反転
	inline Vec2
		operator-(const Vec2& a)
	{
		return Vec2(-a.x, -a.y);
	}

	//スカラーとの演算
	inline Vec2
		operator*(const Vec2& a, const double& b)
	{
		Vec2 ret = a;
		ret *= b;

		return ret;
	}

	inline Vec2
		operator/(const Vec2& a, const double& b)
	{
		Vec2 ret = a;
		ret /= b;

		return ret;
	}

	inline Vec2
		operator*(const double& b, const Vec2& a)
	{
		Vec2 ret = a;
		ret *= b;

		return ret;
	}

	//内積
	inline double
		Dot(const Vec2& a, const Vec2& b)
	{
		return a.x*b.x + a.y*b.y;
	}

	//デターミナント
	inline double
		Det(const Vec2& a, const Vec2& b)
	{
		return a.x*b.y - a.y*b.x;
	}

	//外積
	inline Vec2
		CrossVecX(const Vec2& a)
	{
		return Vec2(a.y, -a.x);
	}


	//最小化
	inline Vec2
		Minimize(const Vec2& v1, const Vec2& v2)
	{
		return Vec2(v1.x < v2.x ? v1.x : v2.x,
			v1.y < v2.y ? v1.y : v2.y);
	}

	//最大化
	inline Vec2
		Maximize(const Vec2& v1, const Vec2& v2)
	{
		return Vec2(v1.x >= v2.x ? v1.x : v2.x,
			v1.y >= v2.y ? v1.y : v2.y);
	}

	//比較
	inline int
		operator == (const Vec2& v1, const Vec2& v2)
	{
		return v1.x == v2.x && v1.y == v2.y;
	}

	//ストリーム
	inline
		BaseOStream& operator << (BaseOStream& os, const Vec2& a)
	{
		if (os.ascii()) {
			os << a.x << _T(" ") << a.y;
		}
		else {
			os.write((char*)&a, sizeof(Vec2));
		}
		return os;
	}

	inline
		BaseIStream& operator >> (BaseIStream& is, Vec2& a)
	{
		if (is.ascii()) {
			_TCHAR str[256];
			is.get_token(str, 256);
			a.x = (double)_tstof(str);
			is.get_token(str, 256);
			a.y = (double)_tstof(str);
		}
		else {
			is.read((char*)&a, sizeof(Vec2));
		}
		return is;
	}
	inline
		Vec3::Vec3()
	{
		x = 0; y = 0; z = 0;
	}

	inline
		Vec3::Vec3(double _x, double _y, double _z)
	{
		x = _x; y = _y; z = _z;
	}

	inline void
		Vec3::Set(double _x, double _y, double _z)
	{
		x = _x; y = _y; z = _z;
	}

	inline void
		Vec3::Zero()
	{
		x = 0.0; y = 0.0; z = 0.0;
	}

	//代入演算
	inline Vec3&
		Vec3::operator+=(const Vec3 &a)
	{
		x = x + a.x; y = y + a.y; z = z + a.z;
		return  (*this);
	}

	inline Vec3&
		Vec3::operator-=(const Vec3 &a)
	{
		x = x - a.x; y = y - a.y; z = z - a.z;
		return  (*this);
	}

	inline Vec3&
		Vec3::operator*=(const double &a)
	{
		x = a * x; y = a * y; z = a * z;
		return  (*this);
	}

	inline Vec3&
		Vec3::operator/=(const double &a)
	{
		x = x / a; y = y / a; z = z / a;
		return  (*this);
	}

	//大きさ
	inline double
		Vec3::SquareMagnitude()  const
	{
		return x * x + y * y + z * z;
	}

	inline double
		Vec3::Magnitude()  const
	{
		return sqrt(SquareMagnitude());
	}


	//正規化
	inline Vec3
		Vec3::Normalize() const
	{
		double len = Magnitude();
		Vec3 nor(0.0, 0.0, 0.0);
		if (len != 0)
		{
			nor.x = x / len;
			nor.y = y / len;
			nor.z = z / len;
		}
		return nor;
	}

	// 成分の最大,最小を返す
	inline double
		Vec3::Min()  const
	{
		double ret = x;
		if (y < ret) ret = y;
		if (z < ret) ret = z;
		return ret;
	}

	inline double
		Vec3::Max()  const
	{
		double ret = x;
		if (ret < y) ret = y;
		if (ret < z) ret = z;
		return ret;
	}

	inline bool
		Vec3::Equals(const Vec3& ver, double off) const
	{
		if ((*this - ver).Magnitude() < off) {
			return true;
		}
		else {
			return false;
		}
	}

	//正規化
	inline Vec3
		Vec3::Abs()
	{
		Vec3 r;
		r.x = fabs(x);
		r.y = fabs(y);
		r.z = fabs(z);
		return r;
	}

	///////////////////////////////////////////////////////////////////
	//フレンド関数

	//ベクトル演算
	inline Vec3
		operator+(const Vec3& a, const Vec3& b)
	{
		Vec3 ret = a;
		ret += b;

		return ret;
	}

	inline Vec3
		operator-(const Vec3& a, const Vec3& b)
	{
		Vec3 ret = a;
		ret -= b;

		return ret;
	}

	inline Vec3
		operator/(const Vec3& a, const Vec3& b)
	{
		Vec3 ret = a;
		ret -= b;

		return Vec3(a.x / b.x, a.y / b.y, a.z / b.z);
	}

	// 符合反転
	inline Vec3
		operator-(const Vec3& a)
	{
		return Vec3(-a.x, -a.y, -a.z);
	}

	//スカラーとの演算
	inline Vec3
		operator*(const Vec3& a, const double& b)
	{
		Vec3 ret = a;
		ret *= b;

		return ret;
	}

	inline Vec3
		operator/(const Vec3& a, const double& b)
	{
		Vec3 ret = a;
		ret /= b;

		return ret;
	}

	inline Vec3
		operator*(const double& b, const Vec3& a)
	{
		Vec3 ret = a;
		ret *= b;

		return ret;
	}

	//内積
	inline double
		Dot(const Vec3& a, const Vec3& b)
	{
		return a.x*b.x + a.y*b.y + a.z*b.z;
	}

	inline double
		Dot(const Vec4& a, const Vec4& b)
	{
		return a.x*b.x + a.y*b.y + a.z*b.z;
	}

	//外積
	inline Vec3
		CrossVecX(const Vec3& a, const Vec3& b)
	{
		Vec3 r;
		r.x = a.y * b.z - a.z * b.y;
		r.y = a.z * b.x - a.x * b.z;
		r.z = a.x * b.y - a.y * b.x;

		return r;
	}

	inline Vec4
		CrossVecX(const Vec4& a, const Vec4& b)
	{
		Vec4 r;
		r.x = a.y * b.z - a.z * b.y;
		r.y = a.z * b.x - a.x * b.z;
		r.z = a.x * b.y - a.y * b.x;

		return r;
	}

	// return the cross product of x
	inline double
		Cross_x(const Vec3& a, const Vec3& b)
	{
		return a.y * b.z - a.z * b.y;
	}

	// return the cross product of y
	inline double
		Cross_y(const Vec3& a, const Vec3& b)
	{
		return a.z * b.x - a.x * b.z;
	}

	// return the cross product of z
	inline double
		Cross_z(const Vec3& a, const Vec3& b)
	{
		return a.x * b.y - a.y * b.x;
	}

	//最小化
	inline Vec3
		Minimize(const Vec3& v1, const Vec3& v2)
	{
		return Vec3(v1.x < v2.x ? v1.x : v2.x,
			v1.y < v2.y ? v1.y : v2.y,
			v1.z < v2.z ? v1.z : v2.z);
	}

	//最大化
	inline Vec3
		Maximize(const Vec3& v1, const Vec3& v2)
	{
		return Vec3(v1.x >= v2.x ? v1.x : v2.x,
			v1.y >= v2.y ? v1.y : v2.y,
			v1.z >= v2.z ? v1.z : v2.z);
	}

	//比較
	inline int
		operator == (const Vec3& v1, const Vec3& v2)
	{
		return v1.x == v2.x && v1.y == v2.y && v1.z == v2.z;
	}

	//ストリーム
	inline
		BaseOStream& operator << (BaseOStream& os, const Vec3& a)
	{
		if (os.ascii()) {
			os << a.x << _T(" ") << a.y << _T(" ") << a.z;
		}
		else {
			os.write((char*)&a, sizeof(Vec3));
		}
		return os;
	}

	inline
		BaseIStream& operator >> (BaseIStream& is, Vec3& a)
	{
		if (is.ascii()) {
			_TCHAR str[256];
			is.get_token(str, 256);
			a.x = (double)_tstof(str);
			is.get_token(str, 256);
			a.y = (double)_tstof(str);
			is.get_token(str, 256);
			a.z = (double)_tstof(str);
		}
		else {
			is.read((char*)&a, sizeof(Vec3));
		}
		return is;
	}


	// Following is added by D. XXX - Apr. 15, 2006
	inline
		Vec3 Vec3::multiple(Vec3 point) const
	{
		double xx, yy, zz;

		xx = y * point.z - z * point.y;
		yy = z * point.x - x * point.z;
		zz = x * point.y - y * point.x;

		return Vec3(xx, yy, zz);
	}

	inline
		double Vec3::dotmultiple(Vec3 point) const
	{
		return x * point.x + y * point.y + z * point.z;
	}

	inline
		Vec3 Vec3::CrossVecX(const Vec3& b)
	{
		Vec3 r;
		r.x = y * b.z - z * b.y;
		r.y = z * b.x - x * b.z;
		r.z = x * b.y - y * b.x;

		return r;
	}

	inline
		double Vec3::GetLength()
	{
		return (double)sqrt(x*x + y * y + z * z);
	}

	inline
		Vec3	Vec3::Trans(double dx, double dy, double dz) {
		Vec3 v;

		v.x = x + dx;
		v.y = y + dy;
		v.z = z + dz;

		return v;
	}

	inline
		void Vec3::unit()
	{
		double len;
		len = (double)sqrt(x*x + y * y + z * z);

		if (len != 0)
		{
			x /= len;
			y /= len;
			z /= len;
		}
	}

	inline
		Vec3	Vec3::RotateX(double ang) {
		Vec3 v;
		double sa, ca;

		sa = (double)sin(ang);
		ca = (double)cos(ang);

		v.x = x;
		v.y = y * ca - z * sa;
		v.z = y * sa + z * ca;
		return v;
	}

	inline
		Vec3	Vec3::RotateY(double ang) {
		Vec3 v;

		double sa, ca;
		sa = (double)sin(ang);
		ca = (double)cos(ang);

		v.x = x * ca + z * sa;
		v.y = y;
		v.z = -x * sa + z * ca;
		return v;
	}

	inline
		Vec3	Vec3::RotateZ(double ang) {
		Vec3 v;

		double sa, ca;
		sa = (double)sin(ang);
		ca = (double)cos(ang);

		v.x = x * ca - y * sa;
		v.y = x * sa + y * ca;
		v.z = z;
		return v;
	}

	inline
		Vec3 Vec3::RotAxis(Vec3 vect, double ang)
	{
		Vec3 p;
		double xa, ya;

		p = Vec3(x, y, z);

		vect.ComputeGlbAng(&xa, &ya);

		p = p.RotateX(xa);
		p = p.RotateY(ya);
		p = p.RotateZ(ang);
		p = p.RotateY(-ya);
		p = p.RotateX(-xa);

		return p;
	}

	inline
		Vec4 Vec4::multiple(Vec4 point) const
	{
		double xx, yy, zz;

		xx = y * point.z - z * point.y;
		yy = z * point.x - x * point.z;
		zz = x * point.y - y * point.x;

		return Vec4(xx, yy, zz, point.w);
	}
	inline
		Vec4 Vec4::CrossVecX(const Vec4& b)
	{
		Vec4 r;
		r.x = y * b.z - z * b.y;
		r.y = z * b.x - x * b.z;
		r.z = x * b.y - y * b.x;
		r.w = b.w;
		return r;
	}


	inline
		Vec4	Vec4::Trans(double dx, double dy, double dz) {
		Vec4 v;

		v.x = x + dx;
		v.y = y + dy;
		v.z = z + dz;
		v.w = w;
		return v;
	}


	inline
		Vec4	Vec4::RotateX(double ang) {
		Vec4 v;
		double sa, ca;

		sa = (double)sin(ang);
		ca = (double)cos(ang);

		v.x = x;
		v.y = y * ca - z * sa;
		v.z = y * sa + z * ca;
		v.w = w;
		return v;
	}

	inline
		Vec4	Vec4::RotateY(double ang) {
		Vec4 v;

		double sa, ca;
		sa = (double)sin(ang);
		ca = (double)cos(ang);

		v.x = x * ca + z * sa;
		v.y = y;
		v.z = -x * sa + z * ca;
		v.w = w;
		return v;
	}

	inline
		Vec4	Vec4::RotateZ(double ang) {
		Vec4 v;

		double sa, ca;
		sa = (double)sin(ang);
		ca = (double)cos(ang);

		v.x = x * ca - y * sa;
		v.y = x * sa + y * ca;
		v.z = z;
		v.w = w;
		return v;
	}

	inline
		Vec4 Vec4::RotAxis(Vec4 vect, double ang)
	{
		Vec4 p;
		double xa, ya;

		p = Vec4(x, y, z, w);

		vect.ComputeGlbAng(&xa, &ya);

		p = p.RotateX(xa);
		p = p.RotateY(ya);
		p = p.RotateZ(ang);
		p = p.RotateY(-ya);
		p = p.RotateX(-xa);

		return p;
	}



	// Compute the global transfomation angle of a vector
	inline
		void Vec3::ComputeGlbAng(double *xrot, double *yrot)
	{
		double tmp;
		double pi = 3.14159265358979323846;
		double err = 0.000001;

		// rotate x-axis
		if (fabs(z) < err&&fabs(y) < err) *xrot = 0;
		else if (fabs(z) < err&&y > 0) *xrot = (double)(pi / 2.0);
		else if (fabs(z) < err&&y < 0) *xrot = (double)(-pi / 2.0);
		else {

			if (fabs(z) < err) {
				return;
			}

			*xrot = (double)atan(y / z);
		}

		if (z < 0) {
			if (y > 0) *xrot = (double)(*xrot + pi);
			else    *xrot = (double)(*xrot - pi);
		}

		// rotate y-axis
		if (fabs(y*sin(*xrot) + z * cos(*xrot)) < err&&fabs(x) < err) *yrot = 0;
		else if (fabs(y*sin(*xrot) + z * cos(*xrot)) < err&&x > 0) *yrot = (double)(-pi / 2.0);
		else if (fabs(y*sin(*xrot) + z * cos(*xrot)) < err&&x < 0) *yrot = (double)(pi / 2.0);
		else {
			tmp = y * (double)sin(*xrot) + z * (double)cos(*xrot);

			if (fabs(tmp) < err) {
				return;
			}

			*yrot = (double)atan(-x / tmp);
		}

		if ((y*sin(*xrot) + z * cos(*xrot)) < 0) {
			if (x < 0) *yrot = (double)(*yrot + pi);
			else    *yrot = (double)(*yrot - pi);
		}

		return;
	}


	inline
		Vec4::Vec4()
	{
		x = 0; y = 0; z = 0; w = 1;
	}

	inline
		Vec4::Vec4(double _x, double _y, double _z, double _w)
	{
		x = _x; y = _y; z = _z; w = _w;
	}

	inline void
		Vec4::Set(double _x, double _y, double _z, double _w)
	{
		x = _x; y = _y; z = _z; w = _w;
	}

	inline void
		Vec4::Zero()
	{
		x = 0.0f; y = 0.0f; z = 0.0f; w = 1.0f;
	}

	//代入演算
	inline Vec4&
		Vec4::operator+=(const Vec4 &a)
	{
		x = x + a.x; y = y + a.y; z = z + a.z; w = w + a.w;
		return  (*this);
	}

	inline Vec4&
		Vec4::operator-=(const Vec4 &a)
	{
		x = x - a.x; y = y - a.y; z = z - a.z; w = w - a.w;
		return  (*this);
	}

	inline Vec4&
		Vec4::operator*=(const double &a)
	{
		x = a * x; y = a * y; z = a * z; w = a * w;
		return  (*this);
	}

	inline Vec4&
		Vec4::operator/=(const double &a)
	{
		x = x / a; y = y / a; z = z / a; w = w / a;
		return  (*this);
	}

	//正規化
	inline Vec4
		Vec4::Normalize() const
	{
		double len = Magnitude();
		Vec4 nor;
		nor.x = x / len;
		nor.y = y / len;
		nor.z = z / len;
		nor.w = w;
		return nor;
	}


	//正規化
	inline Vec4
		Vec4::Abs()  const
	{
		Vec4 r;
		r.x = fabs(x);
		r.y = fabs(y);
		r.z = fabs(z);
		r.w = fabs(w);
		return r;
	}

	///////////////////////////////////////////////////////////////////
	//フレンド関数

	//ベクトル演算
	inline Vec4
		operator+(const Vec4& a, const Vec4& b)
	{
		Vec4 ret = a;
		ret += b;

		return ret;
	}

	inline Vec4
		operator-(const Vec4& a, const Vec4& b)
	{
		Vec4 ret = a;
		ret -= b;

		return ret;
	}

	inline Vec4
		operator/(const Vec4& a, const Vec4& b)
	{
		Vec4 ret = a;
		ret -= b;

		return Vec4(a.x / b.x, a.y / b.y, a.z / b.z, a.w / b.w);
	}

	// 符合反転
	inline Vec4
		operator-(const Vec4& a)
	{
		return Vec4(-a.x, -a.y, -a.z, -a.w);
	}

	//スカラーとの演算
	inline Vec4
		operator*(const Vec4& a, const double& b)
	{
		Vec4 ret = a;
		ret *= b;

		return ret;
	}

	inline Vec4
		operator/(const Vec4& a, const double& b)
	{
		Vec4 ret = a;
		ret /= b;

		return ret;
	}

	inline Vec4
		operator*(const double& b, const Vec4& a)
	{
		Vec4 ret = a;
		ret *= b;

		return ret;
	}


	//最小化
	inline Vec4
		Minimize(const Vec4& v1, const Vec4& v2)
	{
		return Vec4(v1.x < v2.x ? v1.x : v2.x,
			v1.y < v2.y ? v1.y : v2.y,
			v1.z < v2.z ? v1.z : v2.z,
			v1.w < v2.w ? v1.w : v2.w
		);
	}

	//最大化
	inline Vec4
		Maximize(const Vec4& v1, const Vec4& v2)
	{
		return Vec4(v1.x >= v2.x ? v1.x : v2.x,
			v1.y >= v2.y ? v1.y : v2.y,
			v1.z >= v2.z ? v1.z : v2.z,
			v1.w >= v2.w ? v1.w : v2.w
		);
	}

	//比較
	inline int
		operator == (const Vec4& v1, const Vec4& v2)
	{
		return v1.x == v2.x && v1.y == v2.y && v1.z == v2.z && v1.w == v2.w;
	}


}


