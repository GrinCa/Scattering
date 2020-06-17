// Gmsh project created on Wed May 20 10:50:33 2020

sizemesh = 0.4;

coeff = 2;

PML_thickness = 1;

l = 2;

lh = 0.5;

//+
Point(1) = {l, lh, 0, sizemesh/coeff};
//+
Point(2) = {-l, lh, 0, sizemesh/coeff};
//+
Point(3) = {-l, -lh, 0, sizemesh/coeff};
//+
Point(4) = {l, -lh, 0, sizemesh/coeff};
//+
Point(5) = {l+PML_thickness, -lh-PML_thickness, 0, sizemesh};
//+
Point(6) = {l+PML_thickness, lh+PML_thickness, 0, sizemesh};
//+
Point(7) = {-l-PML_thickness, lh+PML_thickness, 0, sizemesh};
//+
Point(8) = {-l-PML_thickness, -lh-PML_thickness, 0, sizemesh};
//+
Point(9) = {l, lh, l, sizemesh/coeff};
//+
Point(10) = {-l, lh, l, sizemesh/coeff};
//+
Point(11) = {-l, -lh, l, sizemesh/coeff};
//+
Point(12) = {l, -lh, l, sizemesh/coeff};
//+
Point(13) = {l+PML_thickness, -lh-PML_thickness, l+PML_thickness, sizemesh};
//+
Point(14) = {l+PML_thickness, lh+PML_thickness, l+PML_thickness, sizemesh};
//+
Point(15) = {-l-PML_thickness, lh+PML_thickness, l+PML_thickness, sizemesh};
//+
Point(16) = {-l-PML_thickness, -lh-PML_thickness, l+PML_thickness, sizemesh};
//+
Point(17) = {0, lh, 0, sizemesh/coeff};
//+
Point(18) = {0, -lh, 0, sizemesh/coeff};
//+
Point(19) = {0, lh, l, sizemesh/coeff};
//+
Point(20) = {0, -lh, l, sizemesh/coeff};
//+

//+
Line(2) = {4, 1};
//+

//+
Line(4) = {2, 3};
//+
Line(5) = {8, 5};
//+
Line(6) = {5, 6};
//+
Line(7) = {6, 7};
//+
Line(8) = {7, 8};
//+
Line(9) = {8, 16};
//+
Line(10) = {16, 13};
//+
Line(11) = {13, 5};
//+
Line(12) = {7, 15};
//+
Line(13) = {15, 14};
//+
Line(14) = {14, 6};
//+
Line(15) = {14, 13};
//+
Line(16) = {15, 16};
//+
Line(17) = {3, 11};
//+
Line(18) = {11, 10};
//+
Line(19) = {10, 2};
//+
Line(20) = {4, 12};
//+
Line(21) = {12, 9};
//+
Line(22) = {9, 1};
//+

//+
Line(23) = {11, 20};
//+
Line(24) = {20, 12};
//+
Line(25) = {9, 19};
//+
Line(26) = {19, 10};
//+
Line(27) = {2, 17};
//+
Line(28) = {17, 1};
//+
Line(29) = {4, 18};
//+
Line(30) = {18, 3};
//+
Line(31) = {18, 17};
//+
Line(32) = {17, 19};
//+
Line(33) = {19, 20};
//+
Line(34) = {20, 18};
//+
Line Loop(1) = {30, -4, 27, -31};
//+
Plane Surface(1) = {1};
//+
Line Loop(2) = {31, 28, -2, 29};
//+
Plane Surface(2) = {2};
//+
Line Loop(3) = {2, -22, -21, -20};
//+
Plane Surface(3) = {3};
//+
Line Loop(4) = {21, 25, 33, 24};
//+
Plane Surface(4) = {4};
//+
Line Loop(5) = {26, -18, 23, -33};
//+
Plane Surface(5) = {5};
//+
Line Loop(6) = {18, 19, 4, 17};
//+
Plane Surface(6) = {6};
//+
Line Loop(7) = {31, 32, 33, 34};
//+
Plane Surface(7) = {7};
//+
Line Loop(8) = {30, 17, 23, 34};
//+
Plane Surface(8) = {8};
//+
Line Loop(9) = {29, -34, 24, -20};
//+
Plane Surface(9) = {9};
//+
Line Loop(10) = {27, 32, 26, 19};
//+
Plane Surface(10) = {10};
//+
Line Loop(11) = {28, -22, 25, -32};
//+
Plane Surface(11) = {11};
//+
Line Loop(13) = {6, -14, 15, 11};
//+
Plane Surface(13) = {13};
//+
Line Loop(14) = {15, -10, -16, 13};
//+
Plane Surface(14) = {14};
//+
Line Loop(15) = {16, -9, -8, 12};
//+
Plane Surface(15) = {15};
//+
Line Loop(16) = {5, -11, -10, -9};
//+
Plane Surface(16) = {16};
//+
Line Loop(17) = {7, 12, 13, 14};
//+
Plane Surface(17) = {17};
//+
Line Loop(18) = {8, 5, 6, 7};
//+
Line Loop(19) = {27, 28, -2, 29, 30, -4};
//+
Plane Surface(18) = {18, 19};
//+
Surface Loop(1) = {8, 6, 5, 10, 7, 1};
//+
Volume(1) = {1};
//+
Surface Loop(2) = {9, 4, 3, 11, 2, 7};
//+
Volume(2) = {2};
//+
Surface Loop(3) = {16, 18, 15, 14, 13, 17, 9, 4, 3, 11, 8, 6, 5, 10};
//+
Volume(3) = {3};
//+
Physical Volume("BGL") = {1};
//+
Physical Volume("BGR") = {2};
//+
Physical Volume("PML") = {3};
//+
Physical Surface("PMLext") = {16, 13, 14, 17, 15, 18};
//+
Physical Surface("wall") = {2, 1};


//+
Physical Surface("PML_acoustic") = {11, 10, 6, 4, 5, 3, 8, 9};

//+
Physical Surface("BGL_BGR") = {7};
