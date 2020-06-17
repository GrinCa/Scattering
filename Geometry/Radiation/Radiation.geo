// Gmsh project created on Wed Jun 10 09:02:15 2020
//+

sizemesh = 0.07;
coeff = 1/3;


Point(1) = {0, 0.5, 0.5, sizemesh};
//+
Point(2) = {0, 0.5, -0.5, sizemesh};
//+
Point(3) = {0, -0.5, -0.5, sizemesh};
//+
Point(4) = {0, -0.5, 0.5, sizemesh};
//+
Point(5) = {1, 0.5, 0.5, sizemesh};
//+
Point(6) = {1, 0.5, -0.5, sizemesh};
//+
Point(7) = {1, -0.5, -0.5, sizemesh};
//+
Point(8) = {1, -0.5, 0.5, sizemesh};
//+
Point(9) = {0, 0.25, 0.25, sizemesh};
//+
Point(10) = {0, 0.25, -0.25, sizemesh};
//+
Point(11) = {0, -0.25, -0.25, sizemesh};
//+
Point(12) = {0, -0.25, 0.25, sizemesh};
//+
Line(1) = {9, 10};
//+
Line(2) = {10, 11};
//+
Line(3) = {11, 12};
//+
Line(4) = {12, 9};
//+
Line(5) = {1, 4};
//+
Line(6) = {4, 3};
//+
Line(7) = {3, 2};
//+
Line(8) = {2, 1};
//+
Line(9) = {4, 8};
//+
Line(10) = {8, 7};
//+
Line(11) = {7, 3};
//+
Line(12) = {8, 5};
//+
Line(13) = {5, 1};
//+
Line(14) = {5, 6};
//+
Line(15) = {6, 2};
//+
Line(16) = {7, 6};
//+
Line Loop(1) = {13, 5, 9, 12};
//+
Plane Surface(1) = {1};
//+
Line Loop(2) = {13, -8, -15, -14};
//+
Plane Surface(2) = {2};
//+
Line Loop(3) = {7, -15, -16, 11};
//+
Plane Surface(3) = {3};
//+
Line Loop(4) = {6, -11, -10, -9};
//+
Plane Surface(4) = {4};
//+
Line Loop(5) = {12, 14, -16, -10};
//+
Plane Surface(5) = {5};
//+
Line Loop(6) = {5, 6, 7, 8};
//+
Line Loop(7) = {1, 2, 3, 4};
//+
Plane Surface(6) = {6, 7};
//+
Plane Surface(7) = {7};
//+
Point(13) = {0, 1, 1, sizemesh/coeff};
//+
Point(14) = {0, 1, -1, sizemesh/coeff};
//+
Point(15) = {0, -1, -1, sizemesh/coeff};
//+
Point(16) = {0, -1, 1, sizemesh/coeff};
//+
Point(17) = {1.5, 1, 1, sizemesh/coeff};
//+
Point(18) = {1.5, 1, -1, sizemesh/coeff};
//+
Point(19) = {1.5, -1, -1, sizemesh/coeff};
//+
Point(20) = {1.5, -1, 1, sizemesh/coeff};
//+
Line(17) = {14, 13};
//+
Line(18) = {13, 16};
//+
Line(19) = {16, 15};
//+
Line(20) = {15, 14};
//+
Line(21) = {14, 18};
//+
Line(22) = {18, 17};
//+
Line(23) = {17, 13};
//+
Line(24) = {16, 20};
//+
Line(25) = {20, 17};
//+
Line(26) = {20, 19};
//+
Line(27) = {19, 18};
//+
Line(28) = {19, 15};
//+
Line Loop(8) = {19, -28, -26, -24};
//+
Plane Surface(8) = {8};
//+
Line Loop(9) = {18, 24, 25, 23};
//+
Plane Surface(9) = {9};
//+
Line Loop(10) = {17, -23, -22, -21};
//+
Plane Surface(10) = {10};
//+
Line Loop(11) = {21, -27, 28, 20};
//+
Plane Surface(11) = {11};
//+
Line Loop(12) = {25, -22, -27, -26};
//+
Plane Surface(12) = {12};
//+
Line Loop(13) = {20, 17, 18, 19};
//+
Plane Surface(13) = {6, 13};
//+
Surface Loop(1) = {6, 7, 3, 2, 1, 4, 5};
//+
Volume(1) = {1};
//+
Surface Loop(2) = {2, 1, 4, 3, 5, 10, 13, 8, 11, 12, 9};
//+
Volume(2) = {2};
//+
Physical Volume("acoustic",1) = {1};
//+
Physical Volume("PML",2) = {2};
//+
Physical Surface("plate",3) = {7};

