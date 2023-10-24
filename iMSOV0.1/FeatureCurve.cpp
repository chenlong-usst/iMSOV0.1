#include "FeatureCurve.h"
#include <cmath>

//void FeatureCurve::Rolate(NurbsLine &NL, double angles, int mode)
//{
//			if (mode == 1)
//		{
//			for (int i = 0; i < NL._ControlPts->size(); i++)
//			{
//				double w = (*NL._ControlPts)[i].w;
//				(*NL._ControlPts)[i] = (*NL._ControlPts)[i].RotateX(angles);
//				(*NL._ControlPts)[i].w = w;
//			}
//		}
//		if (mode == 2)
//		{
//			for (int i = 0; i < NL._ControlPts->size(); i++)
//			{
//				double w = (*NL._ControlPts)[i].w;
//				(*NL._ControlPts)[i] = (*NL._ControlPts)[i].RotateY(angles);
//				(*NL._ControlPts)[i].w = w;
//			}
//		}
//		if (mode == 3)
//		{
//			for (int i = 0; i < NL._ControlPts->size(); i++)
//			{
//				double w = (*NL._ControlPts)[i].w;
//				(*NL._ControlPts)[i] = (*NL._ControlPts)[i].RotateZ(angles);
//				(*NL._ControlPts)[i].w = w;
//			}
//		}
//}

FeatureCurve_Line::FeatureCurve_Line(double x1, double y1, double x2, double y2)
{
	//point4d inputpoint1(1, 1, 1, 1);
	//point4d inputpoint2(2, 2, 1, 1);
	point4d inputpoint1(x1, y1, 1, 1);
	point4d inputpoint2(x2, y2, 1, 1);

	//线段端点赋值
	point4d SP = point4d(inputpoint1.x, inputpoint1.y, 0, 1);
	point4d EP = point4d(inputpoint2.x, inputpoint2.y, 0, 1);

	//绘制nurbs格式线
	NurbsLine nl1;
	nl1.CreatLineWithTwoPoints(SP, EP);
	

	//计算4个occ格式顶点
	gp_Pnt SP0 = { inputpoint1.x, inputpoint1.y, 0 };
	gp_Pnt EP0 = { inputpoint2.x, inputpoint2.y, 0 };

	//绘制occ格式线
	ol1 = BRepBuilderAPI_MakeEdge(SP0, EP0);

	//存储2条线
	FC.nl = nl1;
	//FC.ol = ol1;//这里赋值还有问题，需要重载运算符

}

FeatureCurve_Rectangle::FeatureCurve_Rectangle(double x1, double y1, double x2, double y2, InitType type)
{
	//point4d inputpoint1(-10, -5, 0, 1);
	//point4d inputpoint2(20, 10, 0, 1);
	point4d inputpoint1(x1, y1, 1, 1);
	point4d inputpoint2(x2, y2, 1, 1);

	//计算4个顶点
	point4d SP = point4d(inputpoint1.x, inputpoint1.y, 0, 1);
	point4d EP = point4d(inputpoint2.x, inputpoint2.y, 0, 1);
	point4d SP1 = point4d(inputpoint1.x, inputpoint2.y, 0, 1);
	point4d EP1 = point4d(inputpoint2.x, inputpoint1.y, 0, 1);

	//绘制4条nurbs格式线
	NurbsLine nl1, nl2, nl3, nl4;
	//ctg 这里改了下顺序，逆时针连线
	nl1.CreatLineWithTwoPoints(SP, EP1);
	nl2.CreatLineWithTwoPoints(EP1, EP);
	nl3.CreatLineWithTwoPoints(EP, SP1);
	nl4.CreatLineWithTwoPoints(SP1, SP);
	//分别存储8条线
	FeatureCurve l1, l2, l3, l4;
	l1.nl = nl1;
	l2.nl = nl2;	
	l3.nl = nl3;	
	l4.nl = nl4;
	//l1.ol = ol1;//这里赋值还有问题，需要重载运算符
	//l2.ol = ol2; 
	//l3.ol = ol3; 
	//l4.ol = ol4;

	FC.push_back(l1);
	FC.push_back(l2);
	FC.push_back(l3);
	FC.push_back(l4);

	if (type == Edges)InitEdges(inputpoint1.x, inputpoint1.y, inputpoint2.x, inputpoint2.y);
	else if(type == Wire)InitWire(inputpoint1.x, inputpoint1.y, inputpoint2.x, inputpoint2.y);
	
}





void FeatureCurve_Rectangle::InitEdges(double x1, double y1, double x2, double y2)
{
	//计算4个occ格式顶点
	gp_Pnt SP0 = { x1, y1, 0 };
	gp_Pnt EP0 = { x2, y2, 0 };
	gp_Pnt SP01 = { x1, y2, 0 };
	gp_Pnt EP01 = { x2, y1, 0 };

	//绘制4条occ格式线
	TopoDS_Edge ol1 = BRepBuilderAPI_MakeEdge(SP0, EP01);
	TopoDS_Edge ol2 = BRepBuilderAPI_MakeEdge(EP01, EP0);
	TopoDS_Edge ol3 = BRepBuilderAPI_MakeEdge(EP0, SP01);
	TopoDS_Edge ol4 = BRepBuilderAPI_MakeEdge(SP01, SP0);

	aEdges.push_back(ol1);
	aEdges.push_back(ol2);
	aEdges.push_back(ol3);
	aEdges.push_back(ol4);

	cout << "TopoDS_Edges was created. " << endl;
}



void FeatureCurve_Rectangle::InitWire(double x1, double y1, double x2, double y2)
{
	//计算4个occ格式顶点
	gp_Pnt SP0 = { x1, y1, 0 };
	gp_Pnt EP0 = { x2, y2, 0 };
	gp_Pnt SP01 = { x1, y2, 0 };
	gp_Pnt EP01 = { x2, y1, 0 };

	//绘制4条occ格式线
	TopoDS_Edge ol1 = BRepBuilderAPI_MakeEdge(SP0, EP01);
	TopoDS_Edge ol2 = BRepBuilderAPI_MakeEdge(EP01, EP0);
	TopoDS_Edge ol3 = BRepBuilderAPI_MakeEdge(EP0, SP01);
	TopoDS_Edge ol4 = BRepBuilderAPI_MakeEdge(SP01, SP0);

	BRepBuilderAPI_MakeWire aMakeWire;
	aMakeWire.Add(ol1);
	aMakeWire.Add(ol2);
	aMakeWire.Add(ol3);
	aMakeWire.Add(ol4);
	//线框
	Standard_ASSERT_VOID(aMakeWire.IsDone(), "Added edge isn't connectible!");
	aWire = aMakeWire.Wire();
	cout << "TopoDS_Wire was created. Vertices :" << endl;
	//打印 遍历顶点
	TopTools_IndexedMapOfShape aVertices;
	TopExp::MapShapes(aWire, TopAbs_VERTEX, aVertices);
	for (TopTools_IndexedMapOfShape::Iterator anIt(aVertices); anIt.More(); anIt.Next())
	{
		TopoDS_Vertex aVertex = TopoDS::Vertex(anIt.Value());
		gp_Pnt aPnt = BRep_Tool::Pnt(aVertex);
		cout << "[ " << aPnt.X() << ", " << aPnt.Y() << ", " << aPnt.Z() << " ]" << endl;
	}
}



//FeatureCurve_Circle_arc::FeatureCurve_Circle_arc()
//{
//	point4d inputpoint1(1, 1, 1, 1);
//	point4d inputpoint2(2, 2, 1, 1);
//
//	point4d p(0, r, 0);
//	point4d p0, p1, p2;
//	double ang2 = 0.5*ang;
//	double Op1 = r / cos(ang2);
//	double w1 = cos(ang2);
//	p0 = p.RotateZ(-ang2);
//	p2 = p.RotateZ(ang2);
//	p1 = point4d(0, Op1, 0);
//
//	FC.nl._u_Degree = 2;
//	FC.nl._u_Knots->push_back(0);
//	FC.nl._u_Knots->push_back(0);
//	FC.nl._u_Knots->push_back(0);
//	FC.nl._u_Knots->push_back(1);
//	FC.nl._u_Knots->push_back(1);
//	FC.nl._u_Knots->push_back(1);
//	FC.nl._ControlPts->push_back(point4d(p0.x, p0.y, p0.z, 1));
//	FC.nl._ControlPts->push_back(point4d(p1.x, p1.y, p1.z, w1));
//	FC.nl._ControlPts->push_back(point4d(p2.x, p2.y, p2.z, 1));
//
//}

FeatureCurve_Circle::FeatureCurve_Circle() 
{
	point2d inputpoint1(6, 14);
	point2d inputpoint2(2, 0);

	//圆心、半径
	point2d cp = inputpoint1;
	double cr = sqrt((inputpoint2.x - inputpoint1.x) * (inputpoint2.x - inputpoint1.x) + (inputpoint2.y - inputpoint1.y) * (inputpoint2.y - inputpoint1.y));

	//occ格式圆
	/*gp_Ax2d ocp = { inputpoint1.x, inputpoint1.y, 0 };
	Handle(Geom2d_Circle) OC = new Geom2d_Circle(ocp, cr);*/

	//nurbs格式圆
	FeatureCurve nl1, nl2, nl3, nl4;
	point2d p(cp.x, cp.y+ cr);
	point4d p0, p1, p2;
	double Op1 = cr / cos(0.25*PI);
	double w1 = cos(0.25*PI);
	point2d p00 = p.Rotate(-0.25*PI, cp);
	point2d p02 = p.Rotate(0.25*PI, cp);
	point2d p01 = point2d(cp.x, cp.y + Op1);
	nl1.nl._u_Degree = 2;
	nl1.nl._u_Knots->push_back(0);
	nl1.nl._u_Knots->push_back(0);
	nl1.nl._u_Knots->push_back(0);
	nl1.nl._u_Knots->push_back(1);
	nl1.nl._u_Knots->push_back(1);
	nl1.nl._u_Knots->push_back(1);
	nl1.nl._ControlPts->push_back(point4d(p00.x, p00.y, 0, 1));
	nl1.nl._ControlPts->push_back(point4d(p01.x, p01.y, 0, w1));
	nl1.nl._ControlPts->push_back(point4d(p02.x, p02.y, 0, 1));
	FC.push_back(nl1);
	
	/*point2d p10 = p00.Rotate(-0.5*PI, cp);
	point2d p12 = p02.Rotate(-0.5*PI, cp);
	point2d p11 = p01.Rotate(-0.5*PI, cp);*/
	point2d p10 = p.Rotate(0.75*PI, cp);
	point2d p12 = p.Rotate(0.25*PI, cp);
	point2d p11 = point2d(cp.x - Op1, cp.y);
	nl2.nl._u_Degree = 2;
	nl2.nl._u_Knots->push_back(0);
	nl2.nl._u_Knots->push_back(0);
	nl2.nl._u_Knots->push_back(0);
	nl2.nl._u_Knots->push_back(1);
	nl2.nl._u_Knots->push_back(1);
	nl2.nl._u_Knots->push_back(1);
	nl2.nl._ControlPts->push_back(point4d(p10.x, p10.y, 0, 1));
	nl2.nl._ControlPts->push_back(point4d(p11.x, p11.y, 0, w1));
	nl2.nl._ControlPts->push_back(point4d(p12.x, p12.y, 0, 1));
	FC.push_back(nl2);

	/*point2d p20 = p10.Rotate(-0.5*PI, cp);
	point2d p22 = p12.Rotate(-0.5*PI, cp);
	point2d p21 = p11.Rotate(-0.5*PI, cp);*/
	point2d p20 = p.Rotate(1.25*PI, cp);
	point2d p22 = p.Rotate(0.75*PI, cp);
	point2d p21 = point2d(cp.x, cp.y - Op1);
	nl3.nl._u_Degree = 2;
	nl3.nl._u_Knots->push_back(0);
	nl3.nl._u_Knots->push_back(0);
	nl3.nl._u_Knots->push_back(0);
	nl3.nl._u_Knots->push_back(1);
	nl3.nl._u_Knots->push_back(1);
	nl3.nl._u_Knots->push_back(1);
	nl3.nl._ControlPts->push_back(point4d(p20.x, p20.y, 0, 1));
	nl3.nl._ControlPts->push_back(point4d(p21.x, p21.y, 0, w1));
	nl3.nl._ControlPts->push_back(point4d(p22.x, p22.y, 0, 1));
	FC.push_back(nl3);

	/*point2d p30 = p20.Rotate(-0.5*PI, cp);
	point2d p32 = p22.Rotate(-0.5*PI, cp);
	point2d p31 = p21.Rotate(-0.5*PI, cp);*/
	point2d p30 = p.Rotate(1.75*PI, cp);
	point2d p32 = p.Rotate(1.25*PI, cp);
	point2d p31 = point2d(cp.x + Op1, cp.y);
	nl4.nl._u_Degree = 2;
	nl4.nl._u_Knots->push_back(0);
	nl4.nl._u_Knots->push_back(0);
	nl4.nl._u_Knots->push_back(0);
	nl4.nl._u_Knots->push_back(1);
	nl4.nl._u_Knots->push_back(1);
	nl4.nl._u_Knots->push_back(1);
	nl4.nl._ControlPts->push_back(point4d(p30.x, p30.y, 0, 1));
	nl4.nl._ControlPts->push_back(point4d(p31.x, p31.y, 0, w1));
	nl4.nl._ControlPts->push_back(point4d(p32.x, p32.y, 0, 1));
	FC.push_back(nl4);

}