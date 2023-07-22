clear all
clc
close
%Elemento 1-2
%Propiedades geometricas del elemento
b=input('Ingresar valor de la base de la viga [mm]: ');
h=input('Ingresar valor de la altura de la viga [mm]: ')
l=input('Ingresar valor de la longitud de la viga [mm]: ')

%Propiedades del material del elemento
E=input('Ingresar el valor del modulo de elasticidad del material de la viga [Gpa): ');

%Condiciones iniciales del elemento
teta=input('Ingresar el valor del angulo de inclinacion del elemento con respecto al eje global [grados]: ');
w=input('Ingresar el valor de la carga distribuida del elemento [kN/mm]: ');
n=input('Ingresar el valor de la carga distribuida axial del elemento [kN/mm]: ');

%Desarrollo del Metodo Matricial
c=cosd(teta);
s=sind(teta);
a=b*h;
In=b*(h^3)/12

x=E*In/(1^3);
y=a*(1^2)/In;

%Submatrices de la matriz general del elemento 12

K11_12=   x*[(y*c^2+12*s^2)   (y-12)*s*c       -6*1*s;
             (y-12)*s*c       (y*s^2+12*c^2)    6*1*c;
             -6*1*s           6*1*c             4*1^2];

K12_12=   x*[-(y*c^2+12*s^2)   -(y-12)*s*c       -6*1*s;
             -(y-12)*s*c       -(y*s^2+12*c^2)    6*1*c;
              6*1*s            -6*1*c             2*1^2];

K21_12=   x*[-(y*c^2+12*s^2)   -(y-12)*s*c        6*1*s;
             -(y-12)*s*c       -(y*s^2+12*c^2)    -6*1*c;
             -6*1*s             6*1*c             2*1^2];

K22_12=   x*[(y*c^2+12*s^2)   (y-12)*s*c        6*1*s;
             (y-12)*s*c       (y*s^2+12*c^2)    -6*1*c;
              6*1*s           -6*1*c             4*1^2];

   %MATRIZ GENERAL DEL ELEMENTO 12
   K12=[K11_12 K12_12;
       K21_12 K22_12];

   %Fuerzas de empotramiento perfecto
   FEM_12=[-n*l/2;
       w*l/2;
       (w*l^2)/l/2;
       -n*l/2;
       w*l/2;
       -(w*l^2/12)];


   %MATRIZ GENERAL ENSAMBLADA
   K=[K11_12 K12_12
       K21_12 K22_12]

   %ANALISIS DE LA ESTRUCTURA POR METODO MATRICIAL
   %{F}=[K]*{d}+{FEM}

   syms F1X F1Y M1Z F2X F2Y M2Z
   M2Z=0;
   %Vector fuerza
   F=[F1X
       F1Y
       M1Z
       F2X
       F2Y
       M2Z]
   %Vector desplazamientos
   syms u1x v1y r1 u2x v2y r2
   u1x=0;
   v1y=0;
   r1=0;
   u2x=0;
   v2y=0;
   
   d=[u1x
       v1y
       r1
       u2x
       v2y
       r2];

   %Vector FEM
     
     FEM=[-n*l/2;
       w*l/2;
       (w*l^2)/l/2;
       -n*l/2;
       w*l/2;
       -(w*l^2/12)]


%Solucion sistema de ecuaciones

ECU1=F==(K*d)+FEM;

Sol=solve([ECU1],[F1X F1Y M1Z F2X F2Y r2]);

F1X=vpa(Sol. F1X,2)
F1Y=vpa(Sol. F1Y,2)
M1Z=vpa(Sol. M1Z,2)
F2X=vpa(Sol. F2X,2)
F2Y=vpa(Sol. F2Y,2)
r2=vpa(Sol. r2,2)
