namespace netgen
{
const char * quadrules[] = {
"rule \"Free Quad (1)\"\n",\
"\n",\
"quality 1\n",\
"\n",\
"mappoints\n",\
"(0, 0);\n",\
"(1, 0);\n",\
"\n",\
"maplines\n",\
"(1, 2) del;\n",\
"\n",\
"newpoints\n",\
"(1, 1) { 1 X2 } { };\n",\
"(0, 1) { } { };\n",\
"\n",\
"newlines\n",\
"(3, 2);\n",\
"(4, 3);\n",\
"(1, 4);\n",\
"\n",\
"freearea\n",\
"(0, 0);\n",\
"(1, 0) { 1 X2 } { };\n",\
"(1.5, 1.5) { 1.5 X2 } { };\n",\
"(-0.5, 1.5) { -0.5 X2 } { };\n",\
"\n",\
"elements\n",\
"(1, 2, 3, 4);\n",\
"\n",\
"endrule\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"rule \"Free Quad (5)\"\n",\
"\n",\
"quality 5\n",\
"\n",\
"mappoints\n",\
"(0, 0);\n",\
"(1, 0);\n",\
"\n",\
"maplines\n",\
"(1, 2) del;\n",\
"\n",\
"newpoints\n",\
"(1, 1) { 1 X2 } { };\n",\
"(0, 1) { } { };\n",\
"\n",\
"newlines\n",\
"(3, 2);\n",\
"(4, 3);\n",\
"(1, 4);\n",\
"\n",\
"freearea\n",\
"(0, 0);\n",\
"(1, 0) { 1 X2 } { };\n",\
"(1.5, 1.5) { 1.5 X2 } { };\n",\
"(-0.5, 1.5) { -0.5 X2 } { };\n",\
"\n",\
"freearea2\n",\
"(0, 0);\n",\
"(1, 0) { 1 X2 } { };\n",\
"(1, 1) { 1 X2 } { };\n",\
"(0, 1) { } { };\n",\
"\n",\
"elements\n",\
"(1, 2, 3, 4);\n",\
"\n",\
"endrule\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"rule \"Quad Right (1)\"\n",\
"\n",\
"quality 1\n",\
"\n",\
"mappoints\n",\
"(0, 0);\n",\
"(1, 0);\n",\
"(1, 1);\n",\
"\n",\
"maplines\n",\
"(1, 2) del;\n",\
"(2, 3) del;\n",\
"\n",\
"newpoints\n",\
"(0, 1) { } { 1 y3 };\n",\
"\n",\
"newlines\n",\
"(1, 4);\n",\
"(4, 3);\n",\
"\n",\
"freearea\n",\
"(0, 0);\n",\
"(1, 0) { 1 X2 } { };\n",\
"(1, 1) { 1 X3 } { 1 Y3 };\n",\
"(-0.5, 1.5) { } { 1.5 Y3 };\n",\
"\n",\
"freearea2\n",\
"(0, 0);\n",\
"(1, 0) { 1 X2 } { };\n",\
"(1, 1) { 1 X3 } { 1 Y3 };\n",\
"(0, 1) { } { 1 Y3 };\n",\
"\n",\
"elements\n",\
"(1, 2, 3, 4);\n",\
"\n",\
"endrule\n",\
"\n",\
"\n",\
"\n",\
"rule \"Quad P Right (2)\"\n",\
"\n",\
"quality 2\n",\
"\n",\
"mappoints\n",\
"(0, 0);\n",\
"(1, 0);\n",\
"(1, 1);\n",\
"\n",\
"maplines\n",\
"(1, 2) del;\n",\
"\n",\
"newpoints\n",\
"(0, 1) { -1 X2, 1 X3 } { 1 Y3 };\n",\
"\n",\
"newlines\n",\
"(1, 4);\n",\
"(4, 3);\n",\
"(3, 2);\n",\
"\n",\
"freearea\n",\
"(0, 0);\n",\
"(1, 0) { 1 X2 } { };\n",\
"(1.2, 0.5) { 0.7 X2, 0.5 X3 } { 0.5 Y3 };\n",\
"(1, 1) { 1 X3 } { 1 Y3 };\n",\
"(-0.5, 1.5) { -2 X2, 1.5 X3 } { 1.5 Y3 };\n",\
"\n",\
"freearea2\n",\
"(0, 0);\n",\
"(1, 0) { 1 X2 } { };\n",\
"(1, 0.5) { 0.5 X2, 0.5 X3 } { 0.5 Y3 };\n",\
"(1, 1) { 1 X3 } { 1 Y3 };\n",\
"(0, 1) { -1 X2, 1 X3 } { 1 Y3 };\n",\
"\n",\
"elements\n",\
"(1, 2, 3, 4);\n",\
"\n",\
"\n",\
"orientations\n",\
"(1, 2, 3);\n",\
"\n",\
"\n",\
"endrule\n",\
"\n",\
"\n",\
"rule \"Quad Right PL (2)\"\n",\
"\n",\
"quality 2\n",\
"\n",\
"mappoints\n",\
"(0, 0);\n",\
"(1, 0);\n",\
"(1, 1);\n",\
"(0, 1);\n",\
"\n",\
"maplines\n",\
"(1, 2) del;\n",\
"(2, 3) del;\n",\
"\n",\
"newpoints\n",\
"\n",\
"newlines\n",\
"(1, 4);\n",\
"(4, 3);\n",\
"\n",\
"freearea\n",\
"(0, 0);\n",\
"(1, 0) { 1 X2 } { };\n",\
"(1, 1) { 1 X3 } { 1 Y3 };\n",\
"(0.5, 1.2) { -0.1 X2, 0.6 X3, 0.6 X4 } { -0.1 Y2, 0.6 Y3, 0.6 Y4 };\n",\
"(0, 1) { 1 X4 } { 1 Y4 };\n",\
"(-0.2, 0.5) { -0.1 X2, -0.1 X3, 0.6 X4 } { -0.1 Y2, -0.1 Y3, 0.6 Y4 };\n",\
"\n",\
"freearea2\n",\
"(0, 0);\n",\
"(1, 0) { 1 X2 } { };\n",\
"(1, 1) { 1 X3 } { 1 Y3 };\n",\
"(0.5, 1) { 0.5 X3, 0.5 X4 } { 0.5 Y3, 0.5 Y4 };\n",\
"(0, 1) { 1 X4 } { 1 Y4 };\n",\
"(0, 0.5) { 0.5 X4 } { 0.5 Y4 };\n",\
"\n",\
"elements\n",\
"(1, 2, 3, 4);\n",\
"\n",\
"orientations\n",\
"(1, 2, 3);\n",\
"(1, 3, 4);\n",\
"(1, 2, 4);\n",\
"(4, 2, 3);\n",\
"\n",\
"\n",\
"endrule\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"rule \"Left Quad (1)\"\n",\
"\n",\
"quality 1\n",\
"\n",\
"mappoints\n",\
"(0, 0);\n",\
"(1, 0);\n",\
"(0, 1);\n",\
"\n",\
"maplines\n",\
"(1, 2) del;\n",\
"(3, 1) del;\n",\
"\n",\
"newpoints\n",\
"(1, 1) { 1 X2, 1 X3 } { 1 Y3 };\n",\
"\n",\
"newlines\n",\
"(3, 4);\n",\
"(4, 2);\n",\
"\n",\
"freearea\n",\
"(0, 0);\n",\
"(1, 0) { 1 X2 } { };\n",\
"(1.5, 1.5) { 1.5 X2, 1.5 X3 } { 1.5 Y3 };\n",\
"(0, 1) { 1 X3 } { 1 Y3 };\n",\
"\n",\
"freearea2\n",\
"(0, 0);\n",\
"(1, 0) { 1 X2 } { };\n",\
"(1, 1) { 1 X2, 1 X3 } { 1 Y3 };\n",\
"(0, 1) { 1 X3 } { 1 Y3 };\n",\
"\n",\
"elements\n",\
"(1, 2, 4, 3);\n",\
"\n",\
"endrule\n",\
"\n",\
"\n",\
"\n",\
"rule \"Left P Quad (2)\"\n",\
"\n",\
"quality 2\n",\
"\n",\
"mappoints\n",\
"(0, 0);\n",\
"(1, 0);\n",\
"(0, 1);\n",\
"\n",\
"maplines\n",\
"(1, 2) del;\n",\
"\n",\
"newpoints\n",\
"(1, 1) { 1 X2, 1 X3 } { 1 Y3 };\n",\
"\n",\
"newlines\n",\
"(1, 3);\n",\
"(3, 4);\n",\
"(4, 2);\n",\
"\n",\
"freearea\n",\
"(0, 0);\n",\
"(1, 0) { 1 X2 } { };\n",\
"(1.5, 1.5) { 1.5 X2, 1.5 X3 } { 1.5 Y3 };\n",\
"(0, 1) { 1 X3 } { 1 Y3 };\n",\
"(-0.2, 0.6) { -0.2 X2, 0.6 X3 } { 0.6 Y3 };\n",\
"\n",\
"freearea2\n",\
"(0, 0);\n",\
"(1, 0) { 1 X2 } { };\n",\
"(1, 1) { 1 X2, 1 X3 } { 1 Y3 };\n",\
"(0, 1) { 1 X3 } { 1 Y3 };\n",\
"(0, 0.5) { 0.5 X3 } { 0.5 Y3 };\n",\
"\n",\
"elements\n",\
"(1, 2, 4, 3);\n",\
"\n",\
"endrule\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"rule \"Left Quad RP (2)\"\n",\
"\n",\
"quality 2\n",\
"\n",\
"mappoints\n",\
"(0, 0);\n",\
"(1, 0);\n",\
"(0, 1);\n",\
"(1, 1);\n",\
"\n",\
"maplines\n",\
"(1, 2) del;\n",\
"(3, 1) del;\n",\
"\n",\
"newpoints\n",\
"\n",\
"newlines\n",\
"(3, 4);\n",\
"(4, 2);\n",\
"\n",\
"freearea\n",\
"(0, 0);\n",\
"(1, 0) { 1 X2 } { };\n",\
"(1.2, 0.5) { 0.6 X2, 0.6 X4, -0.1 X3 } { 0.6 Y2, 0.6 Y4, -0.1 Y3 };\n",\
"(1, 1) { 1 X4 } { 1 Y4 };\n",\
"(0.5, 1.2) { -0.1 X2, 0.6 X3, 0.6 X4 } { -0.1 Y2, 0.6 Y3, 0.6 Y4 };\n",\
"(0, 1) { 1 X3 } { 1 Y3 };\n",\
"\n",\
"freearea2\n",\
"(0, 0);\n",\
"(1, 0) { 1 X2 } { };\n",\
"(1, 0.5) { 0.5 X2, 0.5 X4 } { 0.5 Y2, 0.5 Y4 };\n",\
"(1, 1) { 1 X4 } { 1 Y4 };\n",\
"(0.5, 1) { 0.5 X3, 0.5 X4 } { 0.5 Y3, 0.5 Y4 };\n",\
"(0, 1) { 1 X3 } { 1 Y3 };\n",\
"\n",\
"elements\n",\
"(1, 2, 4, 3);\n",\
"\n",\
"orientations\n",\
"(1, 2, 4);\n",\
"(1, 4, 3);\n",\
"\n",\
"endrule\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"rule \"Two left (1)\"\n",\
"\n",\
"quality 1\n",\
"\n",\
"mappoints\n",\
"(0, 0);\n",\
"(1, 0);\n",\
"(1, 1);\n",\
"(0, 1);\n",\
"\n",\
"maplines\n",\
"(1, 2) del;\n",\
"(3, 4) del;\n",\
"(4, 1) del;\n",\
"\n",\
"newpoints\n",\
"\n",\
"newlines\n",\
"(3, 2);\n",\
"\n",\
"freearea\n",\
"(0, 0);\n",\
"(1, 0) { 1 X2 } { };\n",\
"(1.5, 0.5) { 0.75 X2, 0.75 X3, -0.25 X4 } { 0.75 Y3, -0.25 Y4 };\n",\
"(1, 1) { 1 X3 } { 1 Y3 };\n",\
"(0, 1) { 1 X4 } { 1 Y4 };\n",\
"\n",\
"\n",\
"freearea2\n",\
"(0, 0);\n",\
"(1, 0) { 1 X2 } { };\n",\
"(1, 0.5) { 0.5 X2, 0.5 X3 } { 0.5 Y3 };\n",\
"(1, 1) { 1 X3 } { 1 Y3 };\n",\
"(0, 1) { 1 X4 } { 1 Y4 };\n",\
"\n",\
"elements\n",\
"(1, 2, 3, 4);\n",\
"\n",\
"endrule\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"rule \"Two Right (1)\"\n",\
"\n",\
"quality 1\n",\
"\n",\
"mappoints\n",\
"(0, 0);\n",\
"(1, 0);\n",\
"(1, 1);\n",\
"(0, 1);\n",\
"\n",\
"maplines\n",\
"(1, 2) del;\n",\
"(2, 3) del;\n",\
"(3, 4) del;\n",\
"\n",\
"newpoints\n",\
"\n",\
"newlines\n",\
"(1, 4);\n",\
"\n",\
"freearea\n",\
"(0, 0);\n",\
"(1, 0) { 1 X2 } { };\n",\
"(1, 1) { 1 X3 } { 1 Y3 };\n",\
"(0, 1) { 1 X4 } { 1 Y4 };\n",\
"(-0.5, 0.5) { -0.25 X2, -0.25 X3, 0.75 X4 } { -0.25 Y3, 0.75 Y4 };\n",\
"\n",\
"freearea2\n",\
"(0, 0);\n",\
"(1, 0) { 1 X2 } { };\n",\
"(1, 1) { 1 X3 } { 1 Y3 };\n",\
"(0, 1) { 1 X4 } { 1 Y4 };\n",\
"(0, 0.5) { 0.5 X4 } { 0.5 Y4 };\n",\
"\n",\
"elements\n",\
"(1, 2, 3, 4);\n",\
"\n",\
"endrule\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"rule \"Right 120 (1)\"\n",\
"\n",\
"quality 1000\n",\
"\n",\
"mappoints\n",\
"(0, 0);\n",\
"(1, 0);\n",\
"(1.5, 0.866);\n",\
"\n",\
"maplines\n",\
"(1, 2) del;\n",\
"(2, 3) del;\n",\
"\n",\
"newpoints\n",\
"(0.5, 0.866) { 1 X3, -1 X2 } { 1 Y3 };\n",\
"\n",\
"newlines\n",\
"(1, 4);\n",\
"(4, 3);\n",\
"\n",\
"freearea\n",\
"(0, 0);\n",\
"(1, 0) { 1 X2 } { };\n",\
"(1.5, 0.866) { 1 X3 } { 1 Y3 };\n",\
"(1, 1.732) { -2 X2, 2 X3 } { 2 Y3 };\n",\
"(0, 1.732) { -3 X2, 2 X3 } { 2 Y3 };\n",\
"(-0.5, 0.866) { -2 X2, 1 X3 } {1 Y3 };\n",\
"\n",\
"elements\n",\
"(1, 2, 4);\n",\
"(2, 3, 4);\n",\
"\n",\
"endrule\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"rule \"Left 120 (1)\"\n",\
"\n",\
"quality 1000\n",\
"\n",\
"mappoints\n",\
"(0, 0);\n",\
"(1, 0);\n",\
"(-0.5, 0.866);\n",\
"\n",\
"maplines\n",\
"(1, 2) del;\n",\
"(3, 1) del;\n",\
"\n",\
"newpoints\n",\
"(0.5, 0.866) { 1 X3, 1 X2 } { 1 Y3 };\n",\
"\n",\
"newlines\n",\
"(3, 4);\n",\
"(4, 2);\n",\
"\n",\
"freearea\n",\
"(0, 0);\n",\
"(1, 0) { 1 X2 } { };\n",\
"(1.5, 0.866) { 2 X2, 1 X3 } { 1 Y3 };\n",\
"(1, 1.732) { 2 X2, 2 X3 } { 2 Y3 };\n",\
"(0, 1.732) { -1 X2, 2 X3 } { 2 Y3 };\n",\
"(-0.5, 0.866) { 1 X3 } {1 Y3 };\n",\
"\n",\
"elements\n",\
"(1, 2, 4);\n",\
"(2, 3, 4);\n",\
"\n",\
"endrule\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"rule \"Left Right (1)\"\n",\
"\n",\
"quality 1\n",\
"\n",\
"mappoints\n",\
"(0, 0);\n",\
"(1, 0);\n",\
"(1, 1);\n",\
"(0, 1);\n",\
"\n",\
"maplines\n",\
"(1, 2) del;\n",\
"(2, 3) del;\n",\
"(4, 1) del;\n",\
"\n",\
"\n",\
"newlines\n",\
"(4, 3);\n",\
"\n",\
"freearea\n",\
"(0, 0);\n",\
"(1, 0) { 1 X2 } { };\n",\
"(1, 1) { 1 X3 } { 1 Y3 };\n",\
"(0.5, 1.5) { -0.25 X2, 0.75 X3, 0.75 X4 } { 0.75 Y3, 0.75 Y4 };\n",\
"(0, 1) { 1 X4 } { 1 Y4 };\n",\
"\n",\
"freearea2\n",\
"(0, 0);\n",\
"(1, 0) { 1 X2 } { };\n",\
"(1, 1) { 1 X3 } { 1 Y3 };\n",\
"(0.5, 1) { 0.5 X3, 0.5 X4 } { 0.5 Y3, 0.5 Y4 };\n",\
"(0, 1) { 1 X4 } { 1 Y4 };\n",\
"\n",\
"\n",\
"elements\n",\
"(1, 2, 3, 4);\n",\
"\n",\
"endrule\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"rule \"Fill Quad\"\n",\
"\n",\
"quality 1\n",\
"\n",\
"mappoints\n",\
"(0, 0);\n",\
"(1, 0);\n",\
"(1, 1);\n",\
"(0, 1);\n",\
"\n",\
"maplines\n",\
"(1, 2) del;\n",\
"(2, 3) del;\n",\
"(3, 4) del;\n",\
"(4, 1) del;\n",\
"\n",\
"newpoints\n",\
"\n",\
"newlines\n",\
"\n",\
"freearea\n",\
"(0, 0);\n",\
"(1, 0) { 1 X2 } { 1 Y2 };\n",\
"(1, 1) { 1 X3 } { 1 Y3 };\n",\
"(0, 1) { 1 X4 } { 1 Y4 };\n",\
"\n",\
"elements\n",\
"(1, 2, 3, 4);\n",\
"\n",\
"endrule\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"rule \"Fill Triangle\"\n",\
"\n",\
"quality 10\n",\
"\n",\
"mappoints\n",\
"(0, 0);\n",\
"(1, 0);\n",\
"(0.5, 0.86);\n",\
"\n",\
"maplines\n",\
"(1, 2) del;\n",\
"(2, 3) del;\n",\
"(3, 1) del;\n",\
"\n",\
"newpoints\n",\
"\n",\
"newlines\n",\
"\n",\
"freearea\n",\
"(0, 0);\n",\
"(1, 0) { 1 X2 } { 1 Y2 };\n",\
"(0.5, 0.86) { 1 X3 } { 1 Y3 };\n",\
"\n",\
"elements\n",\
"(1, 2, 3);\n",\
"\n",\
"endrule\n",\
"\n",\
"\n",\
"\n",\
"rule \"Right 60 (1)\"\n",\
"\n",\
"quality 10\n",\
"\n",\
"mappoints\n",\
"(0, 0);\n",\
"(1, 0) { 0.5, 0, 1.0 };\n",\
"(0.5, 0.866) { 0.6, 0, 0.8 };\n",\
"\n",\
"maplines\n",\
"(1, 2) del;\n",\
"(2, 3) del;\n",\
"\n",\
"newpoints\n",\
"\n",\
"newlines\n",\
"(1, 3);\n",\
"\n",\
"freearea\n",\
"(0, 0);\n",\
"(1, 0) { 1 X2 } { };\n",\
"(0.5, 0.866) { 1 X3 } { 1 Y3 };\n",\
"(-0.125, 0.6495) { -0.5 X2, 0.75 X3 } { -0.5 Y2, 0.75 Y3 };\n",\
"\n",\
"freearea2\n",\
"(0, 0);\n",\
"(1, 0) { 1 X2 } { };\n",\
"(0.5, 0.866) { 1 X3 } { 1 Y3 };\n",\
"(0.25, 0.433) { 0.5 X3 } { 0.5 Y3 };\n",\
"\n",\
"elements\n",\
"(1, 2, 3);\n",\
"\n",\
"endrule\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"rule \"Vis A Vis (2)\"\n",\
"\n",\
"quality 2\n",\
"\n",\
"mappoints\n",\
"(0, 0);\n",\
"(1, 0);\n",\
"(1, 1);\n",\
"(0, 1);\n",\
"\n",\
"maplines\n",\
"(1, 2) del;\n",\
"(3, 4) del;\n",\
"\n",\
"newpoints\n",\
"\n",\
"newlines\n",\
"(1, 4);\n",\
"(3, 2);\n",\
"\n",\
"freearea\n",\
"(0, 0);\n",\
"(1, 0) { 1 X2 } { };\n",\
"(1.5, 0.5) { 0.75 X2, 0.75 X3, -0.25 X4 } { 0.75 Y3, -0.25 Y4 };\n",\
"(1, 1) { 1 X3 } { 1 Y3 };\n",\
"(0, 1) { 1 X4 } { 1 Y4 };\n",\
"(-0.5, 0.5) { -0.25 X2, -0.25 X3, 0.75 X4 } { -0.25 Y3, 0.75 Y4 };\n",\
"\n",\
"freearea2\n",\
"(0, 0);\n",\
"(1, 0) { 1 X2 } { };\n",\
"(1, 0.5) { 0.5 X2, 0.5 X3 } { 0.5 Y3 };\n",\
"(1, 1) { 1 X3 } { 1 Y3 };\n",\
"(0, 1) { 1 X4 } { 1 Y4 };\n",\
"(0, 0.5) { 0.5 X4 } { 0.5 Y4 };\n",\
"\n",\
"elements\n",\
"(1, 2, 3, 4);\n",\
"\n",\
"orientations\n",\
"(1, 3, 4);\n",\
"(2, 3, 4);\n",\
"(1, 2, 3);\n",\
"(1, 2, 4);\n",\
"\n",\
"endrule\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"rule \"2 h Vis A Vis (1)\"\n",\
"\n",\
"quality 3000\n",\
"\n",\
"mappoints\n",\
"(0, 0);\n",\
"(1, 0);\n",\
"(1, 1.732);\n",\
"(0, 1.732);\n",\
"\n",\
"maplines\n",\
"(1, 2) del;\n",\
"(3, 4) del;\n",\
"\n",\
"newpoints\n",\
"(0.5, 0.866) { 0.25 X3, 0.25 X4 } { 0.25 Y2, 0.25 Y3, 0.25 Y4 };\n",\
"\n",\
"newlines\n",\
"(1, 5);\n",\
"(5, 4);\n",\
"(3, 5);\n",\
"(5, 2);\n",\
"\n",\
"freearea\n",\
"(0, 0);\n",\
"(1, 0) { 1 X2 } { 1 Y2 };\n",\
"(1.5, 0.866) { 0.75 X2, 0.75 X3, -0.25 X4 } { 0.75 Y2, 0.75 Y3, -0.25 Y4 };\n",\
"(1, 1.732) { 1 X3 } { 1 Y3 };\n",\
"(0, 1.732) { 1 X4 } { 1 Y4 };\n",\
"(-0.5, 0.866) { 0.75 X4, -0.25 X2, -0.25 X3 } { 0.75 Y4, -0.25 Y3 };\n",\
"\n",\
"elements\n",\
"(1, 2, 5);\n",\
"(3, 4, 5);\n",\
"\n",\
"endrule\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
"\n",\
0};
}
