
#v11.1UW--UNWEIGHED LPF FOR ELECTRONIC COEFFICIENTS IN NEWTRAJv2
#v11.1--Use generic LSFv2 in NEWTRAJv2 to use different basis for projecting momentum (and possibly fitting other things)
#v11 Use MLSF in NEWTRAJv2 to use different basis for projecting momentum
#v10.1 Added NFrFilt to NEWTRAJv2 LPFFT calls--allows different number of fourier functions for filtering
#v10--Updated error calculation in LPF(computes absolute value instead of just square (since it can be complex) and added it as an argument (to track different errors) 
#v9--New SolRZero Routine and Udegen
#v8--Renamed LPF to LPFFT(LPF Fourier Transform, uses the Fourier proj subroutine) and added LPFG (Uses generic basis functions which are passed to the routine)
#   --Updated NEWTRAJv2 with LPFFT routine
#   --Added LPFG to NEWTRAJv2 for momentum fits (NOTE THIS USES POTENTIAL BASIS AS GLOBAL VARIABLE, WILL NEED TO BE CHANGED IF DIFFERENT BASIS IS USED)
#v7--added LPF (low pass filter) and NEWTRAJv2 (includes LPF routine when computing half-step coefficients)

SMTHFIT:=proc(CF,F,omega,y,dy,yleft,Ny)
description "Internal routine to compute coefficients of basis functions for fitting of g. NB number of basis function, w is vector of basis funcion, omega is DISCRETIZED weighting function, y is variable, dy step-size, yleft is the coordinate of leftmost point, Ny is the number of equidistant points in y. CF is the output vector of coefficients. Numeric integration performed via trapezoidal rule";
global NBF,w,Digits:
local S,b,wy,i,j,k,yright,yn,CL:
b:=vector(NBF,0): wy:=vector(NBF): S:=matrix(NBF,NBF,0): CL:=vector(NBF):
yright:=evalf(yleft+dy*(Ny-1)):

#Compute S and b
for k from 1 to Ny do
 yn:=evalf(yleft+dy*(k-1)):
 for i from 1 to NBF do
  wy[i]:=evalf(subs(y=yn,w[i])):
  for j from 1 to i do
   S[i,j]:=evalf(S[i,j]+dy*omega[k]*wy[i]*wy[j]):
   od:
  b[i]:=evalf(b[i]+dy*omega[k]*F[k]*wy[i]):
  od:
 od:

  #subtracting half of endpoints for trapezoidal
  for i from 1 to NBF do
  wy[i]:=evalf(subs(y=yleft,w[i])):
  for j from 1 to i do
   S[i,j]:=evalf(S[i,j]-dy/2*omega[1]*wy[i]*wy[j]):
   od:
  b[i]:=evalf(b[i]-dy/2*omega[1]*F[1]*wy[i]):
  od:
  for i from 1 to NBF do
  wy[i]:=evalf(subs(y=yright,w[i])):
  for j from 1 to i do
   S[i,j]:=evalf(S[i,j]-dy/2*omega[Ny]*wy[i]*wy[j]):
   od:
  b[i]:=evalf(b[i]-dy/2*omega[Ny]*F[Ny]*wy[i]):
  od:

 for i from 1 to NBF do
  for j from i+1 to NBF do
   S[i,j]:=evalf(S[j,i]):
  od: od:

if abs(LinearAlgebra[Determinant](S))<10^(-Digits/2) then 
 print("WARNING IN SMTHFIT: Small Determinant of S", Determinant(S));
fi:

#Compute coeffs and assign them to LinearAlgebra Vector
 CL:=linalg[multiply](linalg[inverse](S),b):
 for i from 1 to NBF do CF[i]:=CL[i]: od: 
end proc:
SMTHFITTRAJ:=proc(CF,F,omega,yn,y,Ny,NBF,w)
description "Internal routine to compute coefficients of basis functions for fitting of F described in trajectories. NBF number of basis function, w is vector of basis funcion, y is variable, Ny is the number of trajectories. CF is the output vector of coefficients. Numeric integration done by summing over trajectory weights (omega), Trajectory posions given in yn";
local S,b,wy,i,j,k,yright,CL:
global Digits:
b:=vector(NBF,0): wy:=vector(NBF): S:=matrix(NBF,NBF,0): CL:=vector(NBF):
yright:=evalf(yleft+dy*(Ny-1)):

#Compute S and b
for k from 1 to Ny do
 for i from 1 to NBF do
  wy[i]:=evalf(subs(y=yn[k],w[i])):
  for j from 1 to i do
   S[i,j]:=evalf(S[i,j]+omega[k]*wy[i]*wy[j]):
   od:
   b[i]:=evalf(b[i]+omega[k]*F[k]*wy[i]):
  od:
 od:
#print("here1");
 for i from 1 to NBF do
  for j from i+1 to NBF do
   S[i,j]:=evalf(S[j,i]):
  od: od:
#print("here2");

if abs(LinearAlgebra[Determinant](S))<10^(-Digits/2) then 
 print("WARNING IN SMTHFITTRAJ: Small Determinant of S", Determinant(S));
fi:

#Compute coeffs and assign them to LinearAlgebra Vector
 CL:=linalg[multiply](linalg[inverse](S),b):
 for i from 1 to NBF do CF[i]:=CL[i]: od: 
#print("here3");
end proc:
NCMOMFIT:=proc(CF,omega,yn,y,Ny,NBF,w,dw)
#Compute Nonclassical momentum in a basis (dw)
#Does not require Amplitudes on non-classical momentum, uses derivatives of basis functions instead (dw)
#Weight function omega assumed to be Trajectory weights for use with random grid
local S,b,wy,dwy,i,j,k,CL:
global Digits:
b:=vector(NBF,0): wy:=vector(NBF): dwy:=vector(NBF): S:=matrix(NBF,NBF,0): CL:=vector(NBF):

#Compute S and b
for k from 1 to Ny do
 for i from 1 to NBF do
  wy[i]:=evalf(subs(y=yn[k],w[i])):
  dwy[i]:=evalf(subs(y=yn[k],dw[i])):
  for j from 1 to i do
   S[i,j]:=evalf(S[i,j]+omega[k]*wy[i]*wy[j]):
   od:
   b[i]:=evalf(b[i]-omega[k]*dwy[i]/2):
  od:
 od:
#print("here1");
 for i from 1 to NBF do
  for j from i+1 to NBF do
   S[i,j]:=evalf(S[j,i]):
  od: od:
#print("here2");

#print(S):

if abs(linalg[det](S))<10^(-Digits/2) then 
 print("WARNING IN NCMOMFIT: Small Determinant of S", linalg[det](S));
fi:

#print("Done Determinant Check"):

#Compute coeffs and assign them to LinearAlgebra Vector
 CL:=linalg[multiply](linalg[inverse](S),b):
 for i from 1 to NBF do CF[i]:=CL[i]: od: 
#print("here3");
end proc:
SMTHVDYN:=proc(CF,Traj,D2av,S0,R,Nbas,NTraj,dNBF,dw,ider):
#Compute Vdyn in a basis w
#It's gradient is computed according to exact Eq. 49
#Basis derivtives dw are needed for fits of gradient
#Returns expansion coefficients of Vdyn
local i1,i2,j,Ri,omega:
global AvF,GVdyn:

#Compute Average force
AvF:=Vector(NTraj):
 for j from 1 to NTraj do
  Ri:=Traj[1][j]: 
  AvF[j]:=AVA0FORCE(R,Ri,Nbas,Nnuc,Znuc,BAPrim,BCPrim,BNPrim,BasMolR,S0,Traj[5][j],conjugate(Traj[5][j]),ider):
 od:


GVdyn:=Vector(NTraj):
for j from 1 to NTraj do GVdyn[j]:=AvF[j]+4*Traj[3][j]+2*D2av[j]: od:

for j from 1 to NTraj do omega[j]:=Traj[4][j]: od:
SMTHFITTRAJ(CF,GVdyn,omega,Traj[1],R,NTraj,dNBF,dw):

end proc:
#Computation of d(Chi)/Chi, through smoothed Ln, with border region fits, boundary points get special treatment
#Fits in the border region are with the weighting functions exp(-OmegaPar*(y-y0)^2) and |Chi|^2.  y0 is the boundary coordinate
LogDer:=proc(Fin,FOut,N,y,dy,ymin,OmegaPar)
global NBF,w,yj:
local i,j,ij,Nhalf,yn,LF,DLF,DRF,CFL,CFR,fa,fb,tm,ymax:
ymax:=ymin+dy*(N-1):  Nhalf:=floor(N/2)-1:

#Compute Ln of Fin, store in LF, exclude boundary points
LF:=Vector(N-2): DLF:=Vector(N-2): DRF:=Vector(Nhalf): CFL:=Vector(NBF): CFR:=Vector(NBF):
LF[1]:=ln(Fin[2]): 
for i from 2 to N-2 do 
  fa:=ln(Fin[i+1]): 
  tm:=Im(fa)-Im(LF[i-1]): j:=round(-tm/2/Pi): LF[i]:=evalf(fa+I*j*2*Pi):  
od:
###print(Im(LF));
###print(Re(LF));

#Checking LF
###P1:=plot(yj,Re(LF)): P2:=plot(yj,Im(LF)): print(P1,P2):


###print("VR1",LF[1..5]);
#Compute numerical derivatives of LF
FinDer3(N-2,LF,DLF,dy,1):
###print("VR2",DLF[1..5]);

#reuse LF to store total weight function, left and right side in one pass
for i from 1 to N-2 do
 yn:=evalf(ymin+dy*i):
 LF[i]:=evalf(Re(Fin[i+1]*conjugate(Fin[i+1]))*(exp(-OmegaPar*(ymin-yn)^2)+exp(-OmegaPar*(ymax-yn)^2))):
od:
# Fit left side, overwrite the numerical LF with weight function
SMTHFIT(CFL,DLF,LF,y,dy,ymin+dy,Nhalf):
###print("VR-L",CFL);
# Fit right side, shifting content of LF and writing shifted DLF into DRF
for i from 1 to Nhalf do LF[i]:=evalf(LF[N-2-Nhalf+i]): DRF[i]:=DLF[N-2-Nhalf+i]: od:
SMTHFIT(CFR,DRF,LF,y,dy,ymin+dy*(N-1-Nhalf),Nhalf):
###print("VR-R",CFR);

FOut[1]:=add(subs(y=ymin,w[i]*CFL[i]),i=1..NBF):
FOut[N]:=add(subs(y=ymax,w[i]*CFR[i]),i=1..NBF):
for i from 2 to N-1 do 
 yn:=evalf(ymin+dy*(i-1)): ###fa:=evalf(exp(-OmegaPar*((ymin-yn)^2))): fb:=evalf(exp(-OmegaPar*((ymax-yn)^2))):
 FOut[i]:=add(subs(y=yn,w[i]*CFL[i]),i=1..NBF):
od:

end proc:
CFNDIFFFD:=proc(N,NB,CA,D2C,D1C,dy)
#Compute nuclear derivatives of NB electronic basis coefficients on N sized grid w/ 7-point finite difference.
#Outputs coefficient in NB x N array
local DIN,i,j,D2V,D1V:

DIN:=Vector(N):
D1V:=Vector(N): D2V:=Vector(N):
for i from 1 to NB do
    for j from 1 to N do DIN[j]:=CA[i,j]: od: 
    FinDer(N,DIN,D2V,dy,2): FinDer(N,DIN,D1V,dy,1): 
    for j from 1 to N do D2C[i,j]:=D2V[j]: D1C[i,j]:=D1V[j]: od:
   od:

end proc:
CFNDIFFFD3:=proc(N,NB,CA,D2C,D1C,dy)
#Compute nuclear derivatives of NB electronic basis coefficients on N sized grid w/ 7-point finite difference.
#Outputs coefficient in NB x N array
local DIN,i,j,D2V,D1V:

DIN:=Vector(N):
D1V:=Vector(N): D2V:=Vector(N):
for i from 1 to NB do
    for j from 1 to N do DIN[j]:=CA[i,j]: od: 
    FinDer3(N,DIN,D2V,dy,2): FinDer3(N,DIN,D1V,dy,1): 
    for j from 1 to N do D2C[i,j]:=D2V[j]: D1C[i,j]:=D1V[j]: od:
   od:

end proc:
FinDer3:=proc(N,F,Fd,dx,der)
#3 pt centered difference
local i: 

if (der=1) then 
 for i from 2 to N-1 do 
  Fd[i]:=(F[i+1]-F[i-1])/2/dx:
 od:
 Fd[1]:=(F[2]-F[1])/dx:
 Fd[N]:=(F[N]-F[N-1])/dx:
elif (der=2) then 
 for i from 2 to N-1 do
  Fd[i]:=(F[i-1]+F[i+1]-2*F[i])/dx^2:
 od:
 Fd[1]:=(F[3]-2*F[2]+F[1])/dx^2:
 Fd[N]:=(F[N-2]-2*F[N-1]+F[N])/dx^2:
else 
 print("Error, final option (order of derivative) must be 1 or 2"):
fi:
end proc:
#Computation of operator action through projection.  F is Wf in the input, operator-acted in the output.
#operator is in Mat in the BasHerm basis
MatProjSub:=proc(F,N,ymin,NProj,Mat,Bas,dx,y)
##General projection of an operator for any basis, NProj is # of basis functions, Bas is set of basis function, BOG is set of basis 
#funtions defined on an equidistant grid. F is the function, Mat is the matrix of the operator on the Basis, N is the number of #gridpoints, dy is the grid spacing. Numerical integration performed via trapazoidal rule.
#Basis functions are put on the gird internally, to allow varying of grid size during propagation
local i,j,w1,P1,P2,Norm,BOG:
BOG:=Matrix(N,NProj):
for i from 1 to NProj do
 for j from 1 to N do BOG[j,i]:=evalf(subs(y=ymin+dx*(j-1),Bas[i])):
od: od:
P1:=Vector(NProj): P2:=Vector(NProj):
# compute projections of F onto Bas
for i from 1 to NProj do
 w1:=add(F[j]*BOG[j,i],j=2..N-1):
 P1[i]:=evalf((w1+F[N]*BOG[N,i]/2+F[1]*BOG[1,i]/2)*dx):
 od:
####print(here);
###print(P1):
#Norm:=add(P1[i]*conjugate(P1[i]),i=1..NProj);
#print("Matrix projection norm error is",1-Norm);
for i from 1 to NProj do P2[i]:=add(Mat[i,j]*P1[j],j=1..NProj): od:
for i from 1 to N do F[i]:=evalf(add(P2[j]*BOG[i,j],j=1..NProj)): od:
end proc:
MINWF:=proc(NBas,H,S,CF)
#Solve for wf coefficients by solving stationary schrodinger equation in a basis(or grid)
#Renormalizes wf
#Returns energy
local E,w1,gstate,Nrm,i,i1,i2:

#Compute Eigenvectors and find GS
w1:=Eigenvectors(H,S):
gstate:=min[index](Re(w1[1])):
E:=w1[1][gstate]:


#Pull-out new coefficients
for i from 1 to NBas do
 CF[i]:=w1[2][i,gstate]:
od:

#Renormalize wf
Nrm:=evalf(add(add(S[i1,i2]*CF[i1]*CF[i2],i1=1..NBas),i2=1..NBas)):
for i from 1 to NBas do CF[i]:=1/sqrt(Nrm)*CF[i]: od:

##print("EigenValues",w1[1]):
E:
end proc:
EXSWF:=proc(istate,NBas,H,S,CF,Print)
#Solve for wf coefficients by solving stationary schrodinger equation in a basis(or grid)
#Renormalizes wf
#Returns Energy and Coefficients of chosen state (istate=1 is ground, istate=2 is 1st exited, etc.)
#Uses seperate subroutine to find index of chosen excites state
local E,w1,gstate,Nrm,i,i1,i2:

#Compute Eigenvectors and find GS
w1:=Eigenvectors(H,S):
if Print=1 then print(w1); fi:
gstate:=EXSTATE(istate,NBas,Re(w1[1]),Print):
E:=Re(w1[1][gstate]):

#Pull-out new coefficients
for i from 1 to NBas do
 CF[i]:=w1[2][i,gstate]:
od:

#Renormalize wf
Nrm:=evalf(add(add(S[i1,i2]*CF[i1]*CF[i2],i1=1..NBas),i2=1..NBas)):
for i from 1 to NBas do CF[i]:=1/sqrt(Nrm)*CF[i]: od:

E:
end proc:
EXSTATE:=proc(istate,NBas,E0,Print):
#Finds index of istate (ground istate=1, first excited istate=2, etc.)
local state,n,EX0,i:

if istate=1 then 
 state:=min[index](E0):
 if Print=1 then print("Ground State is Index",state): fi:
 state:
 
else 
 EX0:=Vector(NBas):
 for i from 1 to NBas do EX0[i]:=E0[i]: od:
 n:=1:
  while n<istate+1 do 
   state:=min[index](EX0):
   EX0[state]:=abs(EX0[state])*10000000:
   n:=n+1:
  od:
 if Print=1 then print("Excited State",istate-1,"is Index",state): fi:
 state: 
fi:

end proc:
FinDer:=proc(N1,Z,Zd,dx,der)
# Computation of first (der=1) or second (der=2) derivatives of Z via 7-point finite difference 
# number of points is 2^m, dx is spacing between the points, Zd is the derivative
local w,i:
if (der=2) then
 w:=evalf(1./180/dx/dx):
 for i from 4 to N1-3 do
  Zd[i]:=evalf((2*Z[i-3]-27*Z[i-2]+270*Z[i-1]-490*Z[i]+270*Z[i+1]-27*Z[i+2]+2*Z[i+3])*w):
 od:
 Zd[3]:=evalf((-13*Z[1]+228*Z[2]-420*Z[3]+200*Z[4]+15*Z[5]-12*Z[6]+2*Z[7])*w):
 Zd[2]:=evalf((137*Z[1]-147*Z[2]-255*Z[3]+470*Z[4]-285*Z[5]+93*Z[6]-13*Z[7])*w):
 Zd[1]:=evalf((812*Z[1]-3132*Z[2]+5265*Z[3]-5080*Z[4]+2970*Z[5]-972*Z[6]+137*Z[7])*w):
 Zd[N1-2]:=evalf((2*Z[N1-6]-12*Z[N1-5]+15*Z[N1-4]+200*Z[N1-3]-420*Z[N1-2]+228*Z[N1-1]-13*Z[N1])*w):
 Zd[N1-1]:=evalf((-13*Z[N1-6]+93*Z[N1-5]-285*Z[N1-4]+470*Z[N1-3]-255*Z[N1-2]-147*Z[N1-1]+137*Z[N1])*w):
 Zd[N1]:=evalf((137*Z[N1-6]-972*Z[N1-5]+2970*Z[N1-4]-5080*Z[N1-3]+5265*Z[N1-2]-3132*Z[N1-1]+812*Z[N1])*w):
elif (der=1) then
  w:=evalf(1./60/dx):
 for i from 4 to N1-3 do
  Zd[i]:=evalf((-Z[i-3]+9*Z[i-2]-45*Z[i-1]+45*Z[i+1]-9*Z[i+2]+Z[i+3])*w):
 od:
 Zd[3]:=evalf((2*Z[1]-24*Z[2]-35*Z[3]+80*Z[4]-30*Z[5]+8*Z[6]-Z[7])*w):
 Zd[2]:=evalf((-10*Z[1]-77*Z[2]+150*Z[3]-100*Z[4]+50*Z[5]-15*Z[6]+2*Z[7])*w):
 Zd[1]:=evalf((-147*Z[1]+360*Z[2]-450*Z[3]+400*Z[4]-225*Z[5]+72*Z[6]-10*Z[7])*w):
 Zd[N1-2]:=evalf((Z[N1-6]-8*Z[N1-5]+30*Z[N1-4]-80*Z[N1-3]+35*Z[N1-2]+24*Z[N1-1]-2*Z[N1])*w):
 Zd[N1-1]:=evalf((-2*Z[N1-6]+15*Z[N1-5]-50*Z[N1-4]+100*Z[N1-3]-150*Z[N1-2]+77*Z[N1-1]+10*Z[N1])*w):
 Zd[N1]:=evalf((10*Z[N1-6]-72*Z[N1-5]+225*Z[N1-4]-400*Z[N1-3]+450*Z[N1-2]-360*Z[N1-1]+147*Z[N1])*w):
else print("STOP! Unknown derivative in FinDer "); return:

fi:
end proc:  
SIGNMATCH:=proc(N,F):
#Subroutine to match signs of a coefficient when they flip on the grid
#N-number of points, F-vector of coefficient on the grid, FOut-output
local i:

for i from 2 to N do 
 if F[i-1]>0 and F[i]<0 then 
  F[i]:=-1*F[i]:
 elif F[i-1]<0 and F[i]>0 then
  F[i]:=-1*F[i]:
 fi:
od:

end proc:
CFNDIFFCHEB:=proc(N,NB,NCheb,CA,D2C,D1C,Dy)
#Compute nuclear derivatives of NB electronic basis coefficients on N sized grid w/ chebyshev
#Outputs coefficient in NB x N array
local DIN,i,j,D2V,D1V:

DIN:=Vector(N):
D1V:=Vector(N): D2V:=Vector(N):
for i from 1 to NB do
    for j from 1 to N do D1V[j]:=CA[i,j]: D2V[j]:=CA[i,j]: od: 
    DiffCheb(1,N,D1V,NCheb,Dy,1): DiffCheb(1,N,D2V,NCheb,Dy,2):
    for j from 1 to N do D2C[i,j]:=D2V[j]: D1C[i,j]:=D1V[j]: od:
   od:

end proc:
CFNDIFFFT:=proc(NTraj,Nbas,NFr,dw0,Traj,D2C,D1C)
#Compute nuclear derivatives of NB electronic basis coefficients on N sized RANDOM grid w/ fourier series
#Outputs 1st and second coefficient derivatives in Nbas x NTraj array
local CF,RV,PF,i,j:

PF:=Array(0..2*NFr,1..Nbas): RV:=Vector(NTraj): CF:=Matrix(Nbas,NTraj):

for i from 1 to NTraj do RV[i]:=Traj[1][i]: od:
for j from 1 to NTraj do for i from 1 to Nbas do CF[i,j]:=Traj[5][j][i]: od: od:

FTProj(NTraj,Nbas,Traj[4],RV,NFr,dw0,CF,PF):
#for j from 1 to NBas do for i from 1 to 2*NFr+1 do print("CF,BF,i",j,i,PF[i-1,j]); od:od:

FTDer(NTraj,Nbas,1,RV,NFr,dw0,D1C,PF): 
FTDer(NTraj,Nbas,2,RV,NFr,dw0,D2C,PF):


end proc:
CFNDIFFFT2:=proc(NTraj,Nbas,NFr,dw0,dshft,Traj,D2C,D1C)
#Compute nuclear derivatives of NB electronic basis coefficients on N sized RANDOM grid w/ fourier series
#Outputs 1st and second coefficient derivatives in Nbas x NTraj array
local CF,RV,PF,i,j:

PF:=Array(0..2*NFr,1..Nbas): RV:=Vector(NTraj): CF:=Matrix(Nbas,NTraj):

for i from 1 to NTraj do RV[i]:=Traj[1][i]: od:
for j from 1 to NTraj do for i from 1 to Nbas do CF[i,j]:=Traj[5][j][i]: od: od:

FTProj2(NTraj,Nbas,Traj[4],RV,NFr,dw0,dshft,CF,PF):
#for j from 1 to NBas do for i from 1 to 2*NFr+1 do print("CF,BF,i",j,i,PF[i-1,j]); od:od:

FTDer2(NTraj,Nbas,1,RV,NFr,dw0,dshft,D1C,PF): 
FTDer2(NTraj,Nbas,2,RV,NFr,dw0,dshft,D2C,PF):


end proc:
Cheb:=proc(imin,imax,C,T,Nmax)
# expansion of C defined on equidistant grid in terms of Chebyshev Polynomilas mapping [-1,1] onto [imin..imax]
# coeffs are returned in T[i] for polynomials of i-1 order
local i,ii,j,x,x1,w1,N,T0,T1,T2,T3,p1,p2,d,Tx:
N:=imax-imin+1: d:=2./(N-1): #number of points and spacing
w1:=(C[imin+2]-2*C[imin+1]+C[imin])/2./d^2: p1:=(C[imin+1]-C[imin])/d-w1*d: #p1 is derivative of f at imin
w1:=(C[imax-2]-2*C[imax-1]+C[imax])/2./d^2: p2:=(C[imax]-C[imax-1])/d+w1*d: #p2 is derivative of f at imax
T2:=(p2-p1)/4: T3:=(p2+p1-(C[imax]-C[imin]))/4:
T0:=(C[imax]+C[imin])/2-T2: T1:=(C[imax]-C[imin])/2-T3:

Tx:=Vector(Nmax): #Tx[1]:=1: Tx[2]:=x: for i from 3 to Nmax do Tx[i]:=2*x*Tx[i-1]-Tx[i-2]: od:
for i from 1 to Nmax do Tx[i]:=simplify(ChebyshevT(i-1,x)): od:
for j from 1 to Nmax do
 #x:=-1+2*(i-1)/(N-1)
 T[j]:=4./Pi/(N-1)*add(subs(x=2*(ii-1)/(N-1)-1.,(C[imin+ii-1]-T0-T1*x-T2*x^2-T3*x^3)*Tx[j]/sqrt(1-x^2)),ii=2..N-1):
od:
T[4]:=evalf(T[4]+T3/4): T[2]:=evalf(T[2]+(T1+3*T3/4.)): T[3]:=evalf(T[3]+T2/2): T[1]:=evalf(T[1]/2+T0+T2/2):
end proc:

ChebD:=proc(T,Nmax,Ndiff,Dx)
# using Chebyshev coefficients T of f(x), with x_max-x_min=Dx, compute T of Ndiff derivative of f
local i,k,c,c1:
if (Ndiff >= Nmax) then for k from 1 to Nmax do T[k]:=0: od:
elif Ndiff=Nmax-1 then T[1]:=2^(Ndiff-1)*Ndiff!*T[Nmax]: for k from 2 to Nmax do T[k]:=0: od:
else
for k from 1 to Ndiff do
 c:=T[Nmax-k]: T[Nmax-k]:=2*(Nmax-k)*T[Nmax-k+1]: T[Nmax-k+1]:=0:
 for i from Nmax-k to 3 by -1 do 
  c1:=T[i-1]: T[i-1]:=evalf(T[i+1]+2*(i-1)*c): c:=evalf(c1): 
 od:
 T[1]:=evalf(T[3]/2+c):
od:
fi:
c:=(2/Dx)^Ndiff: for k from 1 to Nmax-Ndiff do T[k]:=evalf(T[k]*c): od:
end proc:

ChebN:=proc(imin,imax,C,T,Nmax)
# Placing Chebyshev-defined function (Cheb. coeffs in T) back to grid on C
local i,j,x,x1,Tx:
Tx:=Vector(Nmax): #Tx[1]:=1: Tx[2]:=x: for i from 3 to Nmax do Tx[i]:=simplify(2*x*Tx[i-1]-Tx[i-2]): od:
for i from 1 to Nmax do Tx[i]:=simplify(ChebyshevT(i-1,x)): od:
for i from imin to imax do
 x1:=-1.+2*(i-imin)/(imax-imin):
 C[i]:=add(T[j]*evalf(subs(x=x1,Tx[j])),j=1..Nmax):
od:
end proc:

DiffCheb:=proc(imin,imax,C,Nmax,Dx,Ndiff)
# computing derivative of C via Chebyshev expansion
local T:
T:=Vector(Nmax):
Cheb(imin,imax,C,T,Nmax):
ChebD(T,Nmax,Ndiff,Dx):
ChebN(imin,imax,C,T,Nmax):
end proc:
##############################TRAJECTORIES########################################
;
##################PHI VALUES###############
;
DIFFPHIY:=proc(NTraj,NBas,CF,D1C,D2C,D1M,D2M,S0,Mnuc,D1,D2)
#1NDOF
#Compute D1 and D2 of Phi described in a trajectory centered basis
#Numeric derivs of coefficients provided in D1C and D2C (Matrix(NBas,NTraj))
#Retruns array of derivatives of each Phi coefificent at each grid point in D1 and D2 (w/contrtibution of basis fucntions) 
#Basis function derivs provided as arrays of matrices D1M and D2M defined numerically for each point
#Basis fxn overlaps in S0 defined numerically for each point
local DMO,DBF,DCT,DBOC,i1,i2,i,j;

for j from 1 to NTraj do
 for i1 from 1 to NBas do 
  DMO:=-1/2/Mnuc*add(D2C[i2,j]*S0[j][i1,i2],i2=1..NBas):
  DCT:=add(D1C[i2,j]*D1M[j][i1,i2],i2=1..NBas):
  DBF:=add(D2M[j][i1,i2]*CF[j][i2],i2=1..NBas):
  DBOC:=DMO+DCT+DBF:
  D2[i1,j]:=evalf(DBOC):
od:od:
 
for j from 1 to NTraj do
 for i1 from 1 to NBas do 
  DBF:=add(D1M[j][i1,i2]*CF[j][i2],i2=1..NBas):
  DMO:=-1/Mnuc*add(D1C[i2,j]*S0[j][i1,i2],i2=1..NBas):
  DBOC:=DMO+DBF:
  D1[i1,j]:=evalf(DBOC):
od:od: 

end proc:
PPHIAV:=proc(NTraj,NBas,Traj,D1,M,Pphi)
#Compute PPhi_av from x-averaged y derivative (Eq.45) assuming Normalized wf. Outputs av. electronic momentum in vector Pphi
local w1,i,j,PPHIERROR;
global PPhierrFlag:

PPHIERROR:=Vector(NTraj):
for j from 1 to NTraj do
  w1:=I*M*evalf(add(D1[i,j]*conjugate(Traj[5][j][i]),i=1..NBas)):
  Pphi[j]:=Re(w1):
  #if(abs(Im(w1))>10.^(-Digits/2)) then print("Errors in Pphi at y index",j,w1); fi:
  if PPhierrFlag=1 then PPHIERROR[j]:=Im(w1): fi:
od:

if PPhierrFlag=1 then print("PPhi Im errors",PPHIERROR): fi:

end proc:
##################TRAJECTORY INITIALIZATION###############
TRAJPOS:=proc(QT,Qmin,DQ,N):
#1D-Nuclei ONLY: Subroutine to generate initial equidistant trjectory positions for given starting range (Qmin) trajectory spacing (DQ). #Positions output as vector in array of Trajectories (Traj)
local i:

for i from 1 to N do QT[i]:=Qmin+(i-1)*DQ: od:

end proc:
TRAJWT:=proc(CHI,WT,N,dQ):
#1D-Nuclei ONLY: Compute Trajectory Weights from grid based wf
local i,prob:

for i from 2 to N-1 do
 prob:=CHI[i]*conjugate(CHI[i]):
 WT[i]:=evalf(prob*dQ):
od:

 i:=1:
 prob:=CHI[i]*conjugate(CHI[i]):
 WT[i]:=evalf(prob*dQ/2):
 i:=N:
 prob:=CHI[i]*conjugate(CHI[i]):
 WT[i]:=evalf(prob*dQ/2):

end proc:
TRAJINIT:=proc(N,NBas,CHI,C,Q,Qmin,DQ,TINY,PrintLvl):
#1D-Nuclei ONLY:
#Initialize Trajectories from wf CHI on grid of size N
#NNDOF # of nuclear degrees of freedom
#Trajectories truncated based on cutoff in prob density (TINY)
#Phase and nonclassical momenta computed from CHI w/ LogDer
#Electronic coefficients C in basis size NBas on grid (Stored as Array(1..N) of Arrays(1..NBas))
#Returns Truncated Traj
local BIG,QT,WT,LCHI,Traj,i,j:
global beta,NTraj:

BIG:=1E30:
 
QT:=Vector(N): WT:=Vector(N): LCHI:=Vector(N):
Traj:=Array(1..5):


TRAJPOS(QT,Qmin,DQ,N):
if PrintLvl=1 then print("TRAJPOS",QT): fi:

TRAJWT(CHI,WT,N,DQ):
if PrintLvl=1 then print("TRAJWT",WT): fi:

LogDer(CHI,LCHI,N,Q,DQ,Qmin,beta):
if PrintLvl=1 then print("LogDer",LCHI): fi:

for i from 1 to 5 do Traj[i]:=Array(1..N): od:

for i from 1 to N do
 if WT[i]>TINY then
  Traj[1][i]:=QT[i]:
  Traj[2][i]:=Im(LCHI[i]):
  Traj[3][i]:=Re(LCHI[i]):
  Traj[4][i]:=WT[i]:
  Traj[5][i]:=Vector(NBas):
  for j from 1 to NBas do Traj[5][i][j]:=C[i][j]: od:
 else 
  for j from 1 to 5 do Traj[j][i]:=BIG: od:
 fi:
od:

for i from 1 to 5 do 
 Traj[i]:=remove[flatten](x->x=BIG,Traj[i]);
od:

Traj:
  
end proc:
###########CHI VALUES##############
;
QUANTUMFORCE:=proc(NTraj,D2F,D1F,Traj,QF)
description "Internal routine to compute quantum force directly using non-classical momentum (1 NNDOF),EqIn=3 for QF in Eq. 22, EqIn=2 for similar terms in Eq. 23 using wf momentum";
#MULTIPLY BY 1/M PREFACTOR
local i:

for i from 1 to NTraj do QF[i]:=1/2*D2F[i]+Traj[3][i]*D1F[i]: od:

end proc:
##########TIME-DERIVATIVES############
;
PCHIDT:=proc(PCHI,Vdyn,PPhiav,D1P,QF,NTraj,Traj,M)
description "Compute change in chi momentum for given values (NNDOF=1), stores output in array PCHI";
local GVR,Dz,i:

GVR:=Vector(NTraj):
#Compute Grad Vr using Chebyshev
for i from 1 to NTraj do GVR[i]:=Re(Vdyn[i]): od:
Dz:=Traj[1][NTraj]-Traj[1][1]:
DiffCheb(1,NTraj,GVR,NCheb,Dz,1):

for i from 1 to NTraj do
 PCHI[i]:=-GVR[i]+1/M*QF[i]+PPhiav[i]*D1P[i]/M:
od:

end proc:
RCHIDT:=proc(RCHI,Vdyn,PPhiav,D1R,PQF,NTraj,Traj,M)
description "Compute change in chi non-classical momentum for given values (NNDOF=1), stores output in array RCHI";
local GVI,Dz,i:

GVI:=Vector(NTraj):
#Compute Grad Vi using Chebyshev
for i from 1 to NTraj do GVI[i]:=Im(Vdyn[i]): od:
Dz:=Traj[1][NTraj]-Traj[1][1]:
DiffCheb(1,NTraj,GVI,NCheb,Dz,1):

for i from 1 to NTraj do
 RCHI[i]:=GVI[i]-1/M*PQF[i]+PPhiav[i]*D1R[i]/M:
od:

end proc:
TRAJDT:=proc(Traj,NTraj,M,PPhiav,RKTraj)
description "Compute change in Trajectory positions (NNDOF=1), output in Array RKTraj";
local i:

for i from 1 to NTraj do
 RKTraj[i]:=(Traj[2][i]+PPhiav[i])/M:
od:
end proc:
##################Time-Propogator###############
;
RK1PROP:=proc(R,Nbas,NTraj,Traj,S0,K0,D1M,D2M,M,NCheb,h)
#RK1 propogator--Uses RKBLOCK to compute d/dt of relevant EOM for factorized TDSE
local PCHI,RCHI,RKTraj,i,m:
global CFM:

PCHI:=Vector(NTraj): RCHI:=Vector(NTraj): RKTraj:=Vector(NTraj): CFM:=Matrix(Nbas,NTraj):
RKBLOCK(R,Nbas,NTraj,Traj,S0,K0,D1M,D2M,M,NCheb,PCHI,RCHI,RKTraj,CFM):


for i from 1 to NTraj do
 Traj[1][i]:=Traj[1][i]+RKTraj[i]*h: #Update trajectory Pos.
 Traj[2][i]:=Traj[2][i]+PCHI[i]*h: #Update Chi momentum of trajectory
 Traj[3][i]:=Traj[3][i]+RCHI[i]*h: #Update Chi non-classical momentum of trajectory
 for m from 1 to Nbas do
  Traj[5][i][m]:=Traj[5][i][m]+CFM[m,i]*h: #Update Phi coeffs of trajectory
 od:
od:

end proc:
RK1PROPv2:=proc(R,Nbas,NTraj,Traj,S0,K0,D1M,D2M,M,NCheb,h)
#RK1 propogator--Uses RKBLOCK to compute d/dt of relevant EOM for factorized TDSE
local PCHI,RKTraj,i,m:
global CFM:
#print(here);
PCHI:=Vector(NTraj): 
RKTraj:=Vector(NTraj): 
CFM:=Matrix(Nbas,NTraj):
#print(here,PCHI,RKTraj,CFM);
RKBLOCK(R,Nbas,NTraj,Traj,S0,K0,D1M,D2M,M,NCheb,PCHI,RKTraj,CFM):


for i from 1 to NTraj do
 Traj[1][i]:=Traj[1][i]+RKTraj[i]*h: #Update trajectory Pos.
 Traj[2][i]:=Traj[2][i]+PCHI[i]*h: #Update Chi momentum of trajectory
 for m from 1 to Nbas do
  Traj[5][i][m]:=Traj[5][i][m]+CFM[m,i]*h: #Update Phi coeffs of trajectory
 od:
od:

end proc:
NEWTRAJ:=proc(Nbas,S0,NTraj,Traj,TempTraj,RKTraj,PCHI,CFM,h)
#Internal Routine for constructing new Traj array at different time-points for RK propogation 
local m,i,n,eNrm,w0:
global CNrmErr,TrajN:

eNrm:=Vector(NTraj):

for i from 1 to NTraj do
 TempTraj[1][i]:=Traj[1][i]+RKTraj[i]/2*h:
 TempTraj[2][i]:=Traj[2][i]+PCHI[i]/2*h:
 TempTraj[4][i]:=Traj[4][i]:
 for m from 1 to Nbas do TempTraj[5][i][m]:=Traj[5][i][m]+CFM[m,i]/2*h: od:
od:

#Renormalize electornic wf
for i from 1 to NTraj do 
 eNrm[i]:=Re(add(add(TempTraj[5][i][m]*conjugate(TempTraj[5][i][n])*subs(R=TempTraj[1][i],S0[m,n]),m=1..Nbas),n=1..Nbas)):
 for m from 1 to Nbas do TempTraj[5][i][m]:=TempTraj[5][i][m]/sqrt(eNrm[i]): od:
od:

#print("Renormilization Needed",eNrm);

#Average Norm Error
w0:=evalf(TrajN-add(eNrm[i]*TempTraj[4][i],i=1..NTraj)): #TrajN is discretized wf norm defined globally at beginning of WS
CNrmErr:=CNrmErr+w0^2:
#print("CNrmErr",CNrmErr):

end proc:
NEWTRAJv2:=proc(Nbas,S0,NTraj,Traj,TempTraj,RKTraj,PCHI,CFM,h)
#Internal Routine for constructing new Traj array at different time-points for RK propogation 
local m,i,n,j,eNrm,w0,CF,FLTRBas,CFEC,VTTraj1,VTTraj2:
global CNrmErr,TrajN,NFrFilt,dw0,dwshft,NBFm,Mb,NBFec,ECb,LPFErrFT,LPFErrM:

eNrm:=Vector(NTraj): VTTraj1:=Vector(NTraj): VTTraj2:=Vector(NTraj):

CF:=Vector(NBFm): CFEC:=Vector(NBFec):
 
for i from 1 to NTraj do
 TempTraj[1][i]:=Traj[1][i]+RKTraj[i]/2*h:
 TempTraj[2][i]:=Traj[2][i]+PCHI[i]/2*h:
 TempTraj[4][i]:=Traj[4][i]:
 for m from 1 to Nbas do TempTraj[5][i][m]:=Traj[5][i][m]+CFM[m,i]/2*h: od:
od:

#Low pass filter for nuclear momentum

LSFv2(NBFm,Mb,NTraj,CF,TempTraj[2],TempTraj):
LPFErrM:=LPFG(NTraj,TempTraj[1],TempTraj[4],r,NBFm,CF,Mb,TempTraj[2],LPFErrM):

#Low pass filter for coefficients
for j from 1 to NTraj do VTTraj1[j]:=1: od:
for m from 1 to NBas do
for j from 1 to NTraj do VTTraj2[j]:=TempTraj[5][j][m]: od:
LSFUW(NBFec,ECb,NTraj,CFEC,VTTraj2,TempTraj):
LPFErrFT:=LPFG(NTraj,TempTraj[1],VTTraj1,r,NBFec,CFEC,ECb,VTTraj2,LPFErrFT);
for j from 1 to NTraj do TempTraj[5][j][m]:=VTTraj2[j]: od:
od:

for j from 1 to NTraj do VTTraj1[j]:=TempTraj[5][j][1]: od:
for j from 1 to NTraj do VTTraj2[j]:=TempTraj[5][j][2]: od:
#P1:=display(plot(TempTraj[1],abs(VTTraj1),color=blue,legend="|C1|"),plot(TempTraj[1],abs(VTTraj2),color=green,legend="|C2|")):

#Renormalize electornic wf
for i from 1 to NTraj do 
 eNrm[i]:=Re(add(add(TempTraj[5][i][m]*conjugate(TempTraj[5][i][n])*subs(R=TempTraj[1][i],S0[m,n]),m=1..Nbas),n=1..Nbas)):
 for m from 1 to Nbas do TempTraj[5][i][m]:=TempTraj[5][i][m]/sqrt(eNrm[i]): od:
od:

#P2:=display(plot(TempTraj[1],Re(1-~eNrm),legend="Norm Error")):


#print("Half Step",P1,P2,max(abs(1-~Re(eNrm)))):

#print("Renormilization Needed",eNrm);

#Average Norm Error
w0:=evalf(TrajN-add(eNrm[i]*TempTraj[4][i],i=1..NTraj)): #TrajN is discretized wf norm defined globally at beginning of WS
CNrmErr:=CNrmErr+w0^2:
#print("CNrmErr",CNrmErr):

end proc:
HEL:=proc(R,Ri,Nnuc,Znuc,RA,NBas,BAPrim,BCPrim,BNPrim,BasMolR,S0,K0,Hel):
#1 NUCLEAR DOF
#Compute Electronic Hamiltonian at R=Ri
#Number of nuclei Nnuc, charges in Znuc, positions (in variable R) in RA
#Elec basis function parameters for NBas functions in BAPrim,BCPrim,BNPrim,BasMolR
#Anayltic basis KE and Overlaps S0,K0
#Output is NBasxNBas matrix Hel
local AV,BV,RB,w0,i,i1,i2,j:

AV:=Vector(2): AV[1]:=Vector(3): AV[2]:=Vector(3):
BV:=Vector(2): BV[1]:=Vector(3): BV[2]:=Vector(3): 
RB:=Array(1..Nnuc):
for i from 1 to Nnuc do RB[i]:=subs(R=Ri,RA[i]): od:

for i1 from 1 to NBas do
 for i2 from 1 to i1 do
  for j from 1 to 3 do 
    AV[2][j]:=BL[i1][j]: BV[2][j]:=BL[i2][j]:
    AV[1][j]:=eval(RB[BasMolR[i1]][j]): BV[1][j]:=eval(RB[BasMolR[i2]][j]):
   od:

w0:=add(-Znuc[i]*add(add(I1evalB(BAPrim[i1][k1],AV,BAPrim[i2][k2],BV,RB[i])*BCPrim[i1][k1]*BCPrim[i2][k2],k1=1..BNPrim[i1]),k2=1..BNPrim[i2]),i=1..Nnuc):

Hel[i1,i2]:=evalf(subs(R=Ri,K0[i1,i2])+w0+subs(R=Ri,S0[i1,i2]*1./R)):
od:od:
end proc:
HELv2:=proc(R,Ri,Nnuc,Znuc,RA,NBas,BAPrim,BCPrim,BNPrim,BasMolR,S0,K0,Hel):
#1 NUCLEAR DOF
#Compute Electronic Hamiltonian at R=Ri
#Number of nuclei Nnuc, charges in Znuc, positions (in variable R) in RA
#Elec basis function parameters for NBas functions in BAPrim,BCPrim,BNPrim,BasMolR
#Anayltic basis KE and Overlaps S0,K0
#Output is NBasxNBas matrix Hel
#THIS VERSION DOES NOT INCLUDE NUCLEAR-NUCLEAR INTERACTION
local AV,BV,RB,w0,i,i1,i2,j:

AV:=Vector(2): AV[1]:=Vector(3): AV[2]:=Vector(3):
BV:=Vector(2): BV[1]:=Vector(3): BV[2]:=Vector(3): 
RB:=Array(1..Nnuc):
for i from 1 to Nnuc do RB[i]:=subs(R=Ri,RA[i]): od:

for i1 from 1 to NBas do
 for i2 from 1 to i1 do
  for j from 1 to 3 do 
    AV[2][j]:=BL[i1][j]: BV[2][j]:=BL[i2][j]:
    AV[1][j]:=eval(RB[BasMolR[i1]][j]): BV[1][j]:=eval(RB[BasMolR[i2]][j]):
   od:

w0:=add(-Znuc[i]*add(add(I1evalB(BAPrim[i1][k1],AV,BAPrim[i2][k2],BV,RB[i])*BCPrim[i1][k1]*BCPrim[i2][k2],k1=1..BNPrim[i1]),k2=1..BNPrim[i2]),i=1..Nnuc):

Hel[i1,i2]:=evalf(subs(R=Ri,K0[i1,i2])+w0):
od:od:
end proc:
HSAt:=proc(Nnuc,Znuc,RA,NBas,BAPrim,BCPrim,BNPrim,BasMolR,S0,K0,Hel,Sel):
#A limit of electronic Hamiltonian for R=infinity
#Number of nuclei Nnuc, charges in Znuc, positions (in variable R) in RA
#Elec basis function parameters for NBas functions in BAPrim,BCPrim,BNPrim,BasMolR
#Anayltic basis KE and Overlaps S0,K0
#Output are matrices Hel,Sel
local AV,BV,RB,w0,i,i1,i2,k1,k2,j,nnuc:

AV:=Vector(2): AV[1]:=Vector(3): AV[2]:=Vector(3):
BV:=Vector(2): BV[1]:=Vector(3): BV[2]:=Vector(3): 
RB:=Array(1..Nnuc):
for i from 1 to Nnuc do RB[i]:=subs(R=Ri,RA[i]): od:

for i1 from 1 to NBas do
 for i2 from 1 to i1 do
  if BasMolR[i1]=BasMolR[i2] then
   nnuc:=BasMolR[i1]:
   for j from 1 to 3 do 
    AV[2][j]:=BL[i1][j]: BV[2][j]:=BL[i2][j]:
    AV[1][j]:=eval(RB[BasMolR[i1]][j]): BV[1][j]:=eval(RB[BasMolR[i2]][j]):
   od:
   w0:=-Znuc[nnuc]*add(add(I1evalB(BAPrim[i1][k1],AV,BAPrim[i2][k2],BV,RB[nnuc])*BCPrim[i1][k1]*BCPrim[i2][k2],k1=1..BNPrim[i1]),k2=1..BNPrim[i2]):
  Hel[i1,i2]:=evalf(subs(R=0,K0[i1,i2])+w0): Hel[i2,i1]:=evalf(Hel[i1,i2]):
  Sel[i1,i2]:=evalf(subs(R=0,S0[i1,i2])): Sel[i2,i1]:=evalf(Sel[i1,i2]):
 else
  Hel[i1,i2]:=0: Sel[i1,i2]:=0: Hel[i2,i1]:=0: Sel[i2,i1]:=0:
 fi:
 od:od:
end proc:
#########################################################Analytic Hel Debugging Subroutines#######################################
HELFit:=proc(R,Ri,Nnuc,Znuc,RA,NBas,BAPrim,BCPrim,BNPrim,BasMolR,S0,K0,Hel):
#This is debugging subroutine for 2 basis functions!  Symmetric solution is analytic ground state, antisummetric is Exit higher
#Compute Electronic Hamiltonian at R=Ri
#Number of nuclei Nnuc, charges in Znuc, positions (in variable R) in RA
#Elec basis function parameters for NBas functions in BAPrim,BCPrim,BNPrim,BasMolR
#Anayltic basis KE and Overlaps S0,K0
#Output is NBasxNBas matrix Hel
local Exit,r,V0,s:
if(NBas<>2) then print(" Debugging HEL works only for NBas=2"); return: fi:
V0:=-0.58 + 1./r - 1.66*exp(-(11*r)/20) + 0.4*exp(-(11*r)/10); 
Exit:=0.5:
s:=evalf(subs(R=Ri,S0[1,2])):
Hel[1,1]:=evalf(subs(r=Ri,V0)+Exit/2*(1-s)):
Hel[1,2]:=evalf(subs(r=Ri,V0)*(1+s)-Hel[1,1]):
Hel[2,2]:=evalf(Hel[1,1]):

end proc:

HELv2Fit:=proc(R,Ri,Nnuc,Znuc,RA,NBas,BAPrim,BCPrim,BNPrim,BasMolR,S0,K0,Hel):
#This is debugging subroutine for 2 basis functions!  Symmetric solution is analytic ground state, antisummetric is Exit higher
#Compute Electronic Hamiltonian at R=Ri
#Number of nuclei Nnuc, charges in Znuc, positions (in variable R) in RA
#Elec basis function parameters for NBas functions in BAPrim,BCPrim,BNPrim,BasMolR
#Anayltic basis KE and Overlaps S0,K0
#Output is NBasxNBas matrix Hel
local Exit,r,V0,s:
if(NBas<>2) then print(" Debugging HEL works only for NBas=2"); return: fi:
V0:=-0.58 - 1.66*exp(-(11*r)/20) + 0.4*exp(-(11*r)/10); 
Exit:=0.5:
s:=evalf(subs(R=Ri,S0[1,2])):
Hel[1,1]:=evalf(subs(r=Ri,V0)+Exit/2*(1-s)):
Hel[1,2]:=evalf(subs(r=Ri,V0)*(1+s)-Hel[1,1]):
Hel[2,2]:=evalf(Hel[1,1]):

end proc:

HSAtFit:=proc(Nnuc,Znuc,RA,NBas,BAPrim,BCPrim,BNPrim,BasMolR,S0,K0,Hel,Sel):
#A limit of electronic Hamiltonian for R=infinity
#This is debugging subroutine for 2 basis functions!  Symmetric solution is analytic ground state, antisummetric is Exit higher
#Number of nuclei Nnuc, charges in Znuc, positions (in variable R) in RA
#Elec basis function parameters for NBas functions in BAPrim,BCPrim,BNPrim,BasMolR
#Anayltic basis KE and Overlaps S0,K0
#Output are matrices Hel,Sel
local Exit,r,V0,s:
if(NBas<>2) then print(" Debugging HEL works only for NBas=2"); return: fi:
V0:=-0.58; 
Exit:=0.5:
s:=0:
Hel[1,1]:=evalf(subs(r=Ri,V0)+Exit/2*(1-s)):
Hel[1,2]:=evalf(subs(r=Ri,V0)*(1+s)-Hel[1,1]):
Hel[2,2]:=evalf(Hel[1,1]):
Sel[1,1]:=1: Sel[2,2]:=1: Sel[1,2]:=0:

end proc:
AVA0FORCE:=proc(R,Ri,NBas,Nnuc,Znuc,BAPrim,BCPrim,BNPrim,BasMolR,S0,AC,BC,ider):
#Evaluate average coulomb force in given A0 basis at point R=Ri
#This includes the nuclear-nuclear interaction force and electron nuclear interaction force
#Uses I1evalRder to compute electron nuclear interaction force summing over primitives and there derivatives

local w0,AV,BV,RB,F,i1,i2,k1,k2,j,i:

#Compute electron nuclear interaction force for each bf at each nucleus

#Parameters of nuclear coords and angular momenta of orbitals
AV:=Vector(2): AV[1]:=Vector(3): AV[2]:=Vector(3):
BV:=Vector(2): BV[1]:=Vector(3): BV[2]:=Vector(3):
RB:=Vector(Nnuc):
for i from 1 to Nnuc do RB[i]:=subs(R=Ri,RA[i]): od:

F:=0.:
for i from 1 to Nnuc do
 for i1 from 1 to NBas do
  for i2 from 1 to NBas do
   for j from 1 to 3 do 
    AV[2][j]:=BL[i1][j]: BV[2][j]:=BL[i2][j]:
    AV[1][j]:=eval(RB[BasMolR[i1]][j]): BV[1][j]:=eval(RB[BasMolR[i2]][j]):
   od:

   w0:=-Znuc[i]*add(add(I1evalRder(BAPrim[i1][k1],AV,BAPrim[i2][k2],BV,RB[i],ider)*BCPrim[i1][k1]*BCPrim[i2][k2],k1=1..BNPrim[i1]),k2=1..BNPrim[i2]):
   # print("Nnuc,i1,i2,Force=",i,i1,i2,w0):
   F:=w0+F*AC[i1]*BC[i2]: #Averaging over all A0s
od:od:od:
#print("Electron-Nuclear Force is",F);

#Add x-averaged Nuclear interaction force -grad_R(1/R) #linear molecule

for i1 from 1 to NBas do
 for i2 from 1 to NBas do
  F:=F+Znuc[1]*Znuc[2]*evalf(subs(R=Ri,S0[i1,i2]*1/R^2))*AC[i1]*BC[i2]:
od:od:

end proc:
###################################################INTERPOLATION SUBROUTINES###########################################################################
;
InterpCutoff:=proc(S0,CUT,R):
#Cutoff determined from overlap matrix elements(Analytic in nuclear coordinate R overlap matrix provided in S0), value specified in CUT (10^-CUT)
local i,j,F,Rc:
global Rcut:

Rcut:=0.:
for i from 1 to NBas do 
 for j from 1 to i do 
  ##print(i,j);
  if BasMolR[i]<>BasMolR[j] then
  F:=S0[i,j]=CUT:
  ##print(F);
  Rc:=fsolve(F,R);
  if abs(Rc)>Rcut then Rcut:=abs(Rc): fi:
  fi:
od:od:
print("Dissociation cutoff is",Rcut);

end proc:
InterpBoundFits:=proc(K0,S0,NBas,Nnuc,Znuc,R,dR,DisH):
#Generate analytic functions for long range points outside of Range
#Znuc-Vector of nuclear charges
#Nnuc-# of nuclei
#Output DisH, NBasxNBas matrix containing analytic funtions in R
#Assumes Standard basis used

local Rc,F,AV,BV,Edis,RB,Ealp,RC,Hdiagd,LHdiag,c,alpha,c01,f1,f0,i,j,nnuc,m,n,k1,k2,k:
global Rcut,RA,BCPrim,BNPrim,BAPrim,BL,BasMolR:

AV:=Vector(2): AV[1]:=Vector(3): AV[2]:=Vector(3):
BV:=Vector(2): BV[1]:=Vector(3): BV[2]:=Vector(3): 

Ealp:=Vector(3): RC:=Vector(Nnuc):
Hdiagd:=Vector(3): LHdiag:=Vector(3):

for m from 1 to NBas do
 for n from 1 to m do 
 
if BasMolR[m]=BasMolR[n] then 
#Same Atom Matrix Elements--R independent overlaps
 #Use Edis-Znuc_far/R*SO[m,n]
 
 nnuc:=BasMolR[m]:
## print("on Nucleus",nnuc):

 #Compute Atomic energies for given basis functions
 RB:=subs(R=0,RA[nnuc]):
 for j from 1 to 3 do 
    AV[2][j]:=BL[m][j]: BV[2][j]:=BL[n][j]:
    AV[1][j]:=eval(RB[BasMolR[m]][j]): BV[1][j]:=eval(RB[BasMolR[n]][j]):
 od:
 Edis:=subs(R=0,K0[m,n])-Znuc[nnuc]*add(add(I1evalB(BAPrim[m][k1],AV,BAPrim[n][k2],BV,RB)*BCPrim[m][k1]*BCPrim[n][k2],k1=1..BNPrim[m]),k2=1..BNPrim[n]);
 
 #Add in Coloumb term with overlap and construct overall function
 if nnuc=1 then  
  DisH[m,n]:=Edis-mul(Znuc[i],i=nnuc+1..Nnuc)/R*S0[m,n]:
 elif nnuc=Nnuc then
  DisH[m,n]:=Edis-mul(Znuc[i],i=1..Nnuc-1)/R*S0[m,n]:
 fi:

else
 #Different atom basis functions, now based on overlap matrix 
  #Function exp(-alpha*R), alpha will be computed

  #Compute matrix elements 2 points beyond cutoff
  for i from 1 to 3 do
   RC[1]:=subs(R=Rcut+(i-1)*dR,RA[1]): RC[2]:=subs(R=Rcut+(i-1)*dR,RA[2]):
   for j from 1 to 3 do 
    AV[2][j]:=BL[m][j]: BV[2][j]:=BL[n][j]:
    AV[1][j]:=eval(RC[BasMolR[m]][j]): BV[1][j]:=eval(RC[BasMolR[n]][j]):
   od:
  Ealp[i]:=evalf((subs(R=Rcut+(i-1)*dR,K0[m,n])-Znuc[1]*add(add(I1evalB(BAPrim[m][k1],AV,BAPrim[n][k2],BV,RC[1])*BCPrim[m][k1]*BCPrim[n][k2],k1=1..BNPrim[m]),k2=1..BNPrim[n])-Znuc[2]*add(add(I1evalB(BAPrim[m][k1],AV,BAPrim[n][k2],BV,RC[2])*BCPrim[m][k1]*BCPrim[n][k2],k1=1..BNPrim[m]),k2=1..BNPrim[n])));
 od:

## print(m,n,Ealp);

 for i from 1 to 3 do LHdiag[i]:=ln(abs(Ealp[i])): od:
 FinDer3(3,LHdiag,Hdiagd,dZ,1):
 alpha:=Hdiagd[2];

 f1:=c*exp(alpha*R); f0:=subs(R=Rcut,f1); c01:=solve(f0=Ealp[1]); 
 DisH[m,n]:=subs(c=c01,f1);
fi:
od:od:

end proc:
Interp4ptLagrange:=proc(Z,Zindex,p)
# Interpolate equidistant grid function Z using 4pt lagrange scheme
local w:
  w:=evalf(-p*(p-1)*(p-2)/6*Z[Zindex-1]+(p^2-1)*(p-2)/2*Z[Zindex]-p*(p+1)*(p-2)/2*Z[Zindex+1]+p*(p^2-1)/6*Z[Zindex+2]):
  w:
end proc:  
InterpHam:=proc(NR,NBas,Rmin,R,dR,Ri,Hel):
#Compute interpolated electronic hamiltonian at point Ri
#Uses predefined HGrid (precomputed hamiltonian on the grid of size N)
#DisH predefined analytic hamiltonian for points beyond Rcut
local OnPointFlag,TLFlag,TRFlag,Z,w0,w1,p,i,j,m:
global HGrid,DisH,Rcut:

Z:=Vector(NR):

#Use analytic funtions if Ri is beyond cutoff
if Ri>Rcut then
 for i from 1 to NBas do
  for j from 1 to i do 
   Hel[i,j]:=evalf(subs(R=Ri,DisH[i,j])):
 od:od:

#Interpolation
else
 #Find nearest index for point Ri
  w0:=round((Ri-Rmin)/dR+1):
##  print("Initial value",w0);
 #First check if Ri is equal to a grid point 
  w1:=(Ri-Rmin)/dR+1:
  if (w0-w1=0) then
   OnPointFlag:=1: 
##   print("On Grid Case",w0);
 #Now if it is too close to boundaries for 2 points either side
  elif w0<2 then 
   w0:=w0+1: 
   TLFlag:=1:
##   print("Too far left",w0);
  elif w0>NR-2 then 
   w0:=floor((Ri-Rmin)/dR+1)-1: 
   TRFlag:=1:
##   print("Too far right",w0);
  fi:

 if OnPointFlag=1 then 
  for i from 1 to NBas do 
   for j from 1 to i do 
    Hel[i,j]:=HGrid[w0][i,j]
   od:od:
 else

 if TLFlag=1 then p:=(Rmin+(w0-2)*dR-Ri)/dR:
 elif TRFlag=1 then p:=abs(Rmin+(w0)*dR-Ri)/dR:
 else p:=abs(Rmin+(w0-1)*dR-Ri)/dR:
 fi:
 ##print("p",p);  

 for i from 1 to NBas do 
  for j from 1 to i do
   for m from 1 to NR do Z[m]:=HGrid[m][i,j]: od:
   Hel[i,j]:=Interp4ptLagrange(Z,w0,p):
 od:od:
fi:fi:
 
end proc:
InterpEl:=proc(Z,dint,Zmin,Rcut,NBas,NInt,HInt,DisH,H1)
#4-pt Lagrange interpolation of two matrices, stored on the grid in the triangular form.  Assumes H1, S1 (output matrices) are defined as symmetric
local i0,p,i,j,k,n,PA:
PA:=Vector(4):
if(Z>Rcut) then
 for i from 1 to NBas do
  for j from 1 to i do
   k:=i*(i-1)/2+j: 
   H1[i,j]:=evalf(subs(R=Z,DisH[k])): 
 od: od:
else
 i0:=min(max(2,floor((Z-Zmin)/dint+1)),NInt-2):
 p:=evalf((Z-Zmin)/dint-i0+1):
 PA[1]:=-p*(p-1)*(p-2)/6:
 PA[2]:=(p^2-1)*(p-2)/2:
 PA[3]:=-p*(p+1)*(p-2)/2:
 PA[4]:=p*(p^2-1)/6:
 #print("DB Interp",Z,i0,p,PA);
 for i from 1 to NBas do
  for j from 1 to i do
   k:=i*(i-1)/2+j:
   H1[i,j]:=add(HInt[k,i0-2+n]*PA[n],n=1..4):
 od: od:
fi:
end proc:
InterpElv2:=proc(Z,dint,Zmin,Rcut,NBas,NInt,HInt,DisH,AH,H1)
#4-pt Lagrange interpolation of two matrices, stored on the grid in the triangular form.  Assumes H1, S1 (output matrices) are defined as symmetric
#Includes Analytic long-range (DisH) and short range (AH) bounds
local i0,p,i,j,k,n,PA:
PA:=Vector(4):
if(Z>Rcut) then
 for i from 1 to NBas do
  for j from 1 to i do
   k:=i*(i-1)/2+j: 
   H1[i,j]:=evalf(subs(R=Z,DisH[k])): 
 od: od:
elif(Z<Zmin) then
 for i from 1 to NBas do
  for j from 1 to i do
   k:=i*(i-1)/2+j: 
   H1[i,j]:=evalf(subs(R=Z,AH[k])): 
 od: od:
else
 i0:=min(max(2,floor((Z-Zmin)/dint+1)),NInt-2):
 p:=evalf((Z-Zmin)/dint-i0+1):
 PA[1]:=-p*(p-1)*(p-2)/6:
 PA[2]:=(p^2-1)*(p-2)/2:
 PA[3]:=-p*(p+1)*(p-2)/2:
 PA[4]:=p*(p^2-1)/6:
 #print("DB Interp",Z,i0,p,PA);
 for i from 1 to NBas do
  for j from 1 to i do
   k:=i*(i-1)/2+j:
   H1[i,j]:=add(HInt[k,i0-2+n]*PA[n],n=1..4):
 od: od:
fi:
end proc:
InterpElv3:=proc(Z,dint,Zmin,Rcut,NBas,NInt,HInt,DisH,AH,H1)
#4-pt Lagrange interpolation a full matrix. No assumption of symmetry or asymmetry
#Includes Analytic long-range (DisH) and short range (AH) bounds
local i0,p,i,j,k,n,PA:
PA:=Vector(4):
if(Z>Rcut) then
 for i from 1 to NBas do
  for j from 1 to NBas do
   H1[i,j]:=evalf(subs(R=Z,DisH[i,j])): 
 od: od:
elif(Z<Zmin) then
 for i from 1 to NBas do
  for j from 1 to NBas do
   H1[i,j]:=evalf(subs(R=Z,AH[i,j])): 
 od: od:
else
 i0:=min(max(2,floor((Z-Zmin)/dint+1)),NInt-2):
 p:=evalf((Z-Zmin)/dint-i0+1):
 PA[1]:=-p*(p-1)*(p-2)/6:
 PA[2]:=(p^2-1)*(p-2)/2:
 PA[3]:=-p*(p+1)*(p-2)/2:
 PA[4]:=p*(p^2-1)/6:
 #print("DB Interp",Z,i0,p,PA);
 for i from 1 to NBas do
  for j from 1 to NBas do
   H1[i,j]:=add(HInt[i,j,i0-2+n]*PA[n],n=1..4):
 od: od:
fi:
end proc:
VRLSFv2:=proc(NBF,NTraj,CF,Vr,Traj):
#Least squares fit of TDPES using basis Vb
#Fits done discretley, no integration/weights
local SBF,VE,SM,w0,i,j,k:
global Vb,Acc:

SBF:=Matrix(NBF+1,NBF+1,shape=symmetric):
VE:=Vector(NBF+1): w0:=Vector(NBF+1):

for k from 1 to NTraj do
 for i from 1 to NBF+1 do
  for j from 1 to i do
   SBF[i,j]:=SBF[i,j]+evalf(Traj[4][k]*subs(r=Traj[1][k],Vb[i-1]*Vb[j-1])):
  od:
  VE[i]:=VE[i]+evalf(Traj[4][k]*subs(r=Traj[1][k],Vb[i-1]))*Vr[k]: 
od:od:

if abs(LinearAlgebra[Determinant](SBF))<10^(-Digits/2) then 
 print("WARNING IN VRLSFv2: Small Determinant of S", Determinant(SBF));
fi:

SM:=Matrix(NBF+1,NBF+1):

RefinedMatrixInverse(NBF+1,SBF,SM,Acc):

w0:=SM.VE:

for i from 1 to NBF+1 do CF[i]:=w0[i]: od:

end proc:

VRLSFv2:=proc(NBF,NTraj,CF,Vr,Traj):
#Least squares fit of TDPES using basis Vb
#Fits done discretley, no integration/weights
local SBF,VE,SM,w0,i,j,k:
global Vb,DissocLim,Digits:

SBF:=Matrix(NBF+1,NBF+1,shape=symmetric):
VE:=Vector(NBF+1):

for i from 1 to NBF+1 do
 for j from 1 to i do
  for k from 1 to NTraj do
   if Traj[1][k]<DissocLim then SBF[i,j]:=SBF[i,j]+evalf(subs(r=Traj[1][k],Vb[i-1]*Vb[j-1])): fi:
 od: od:
  for k from 1 to NTraj do
   if Traj[1][k]<DissocLim then VE[i]:=VE[i]+evalf(subs(r=Traj[1][k],Vb[i-1]))*Vr[k]: fi:
od: od:
#print("here1",SBF,VE);
SM:=MatrixInverse(SBF):
w0:=SM.VE:

if abs(LinearAlgebra[Determinant](SBF))<10^(-Digits/2) then 
 print("WARNING IN VRLSFv2: Small Determinant of S", Determinant(SBF));
fi:

#print("here3",w0);

for i from 1 to NBF+1 do CF[i]:=w0[i] od:

end proc:
VRLSF:=proc(NBF,NTraj,CF,Vr,Traj):
#Least squares fit of TDPES using basis Vb(globally declared)
#Fits done discretley, no integration/weights
local SBF,VE,SM,w0,i,j:
global Vb,Digits:

SBF:=Matrix(NBF+1,NBF+1,shape=symmetric):
VE:=Vector(NBF+1):

for i from 1 to NBF+1 do
 for j from 1 to i do
  SBF[i,j]:=add(evalf(subs(r=Traj[1][i1],Vb[i-1]*Vb[j-1])),i1=1..NTraj):
 od: 
 VE[i]:=add(evalf(subs(r=Traj[1][i1],Vb[i-1]))*Vr[i1],i1=1..NTraj): 
od:
#print("here1",SBF,VE);
SM:=MatrixInverse(SBF):
w0:=SM.VE:

#print("here3",w0);

if abs(LinearAlgebra[Determinant](SBF))<10^(-Digits/2) then 
 print("WARNING in VRLSF: Small Determinant of S", Determinant(SBF));
fi:

for i from 1 to NBF+1 do CF[i]:=w0[i] od:

end proc:
SolRZeroOlder:=proc(R,Nnuc,Znuc,NBas,BAPrim,BCPrim,BNPrim,BasMolR,S0,K0,E,C):
#A limit of electronic Hamiltonian for R=0
#Number of nuclei Nnuc, charges in Znuc, positions (in variable R) in RA
#Elec basis function parameters for NBas functions in BAPrim,BCPrim,BNPrim,BasMolR
#Anayltic basis KE and Overlaps S0,K0
local AV,BV,RBi,RB0,w0,w1,i,j,i1,i2,j1,j2,tiny,small,S2L,Hel1,Sel1,Hel0,Sel0,Degen,N0,Ri,C0,E0,iPrint:
tiny:=evalf(10^(-Digits)): if(tiny<1.e-16) then tiny=1.e-16: fi: small:=tiny^(1/4): Ri:=evalf(small*100): Degen:=Vector(NBas): S2L:=Vector(NBas):
iPrint:=2:
#S2L stores index of unpruned b.f. from the pruned one.  

for i from 1 to NBas do Degen[i]:=0: od:
N0:=NBas:
for i from 1 to NBas do
 for j from i+1 to NBas do
  w0:=evalf(subs(R=0,S0[i,j])):
  if(iPrint>2) then print("Ovr",i,j,w0,evalb(w0^2>1.-small)); fi:
  if(w0^2>1.-small) then
   Degen[i]:=j: Degen[j]:=i:
   N0:=eval(N0-1):
  fi:
 od:
od:
if(iPrint>1) then print(NBas-N0," degenerate orbitals",Degen); fi:

Hel1:=Vector(NBas): Sel1:=Vector(NBas): 
Hel0:=Matrix(N0,N0,shape=symmetric): Sel0:=Matrix(N0,N0,shape=symmetric):
C0:=Matrix(N0,N0): E0:=Vector(N0):
j:=1:
for i from 1 to NBas do
 if (Degen[i]<i) then S2L[j]:=i: j:=eval(j+1): fi:
od:
  
AV:=Vector(2): AV[1]:=Vector(3): AV[2]:=Vector(3):
BV:=Vector(2): BV[1]:=Vector(3): BV[2]:=Vector(3): 
RBi:=Array(1..Nnuc): RB0:=Array(1..Nnuc):
for i from 1 to Nnuc do RB0[i]:=subs(R=0,RA[i]): RBi[i]:=subs(R=Ri,RA[i]): od:

for i1 from 1 to N0 do
 j1:=S2L[i1]:
 for i2 from 1 to i1 do
  j2:=S2L[i2]:

  for j from 1 to 3 do 
   AV[2][j]:=BL[j1][j]: BV[2][j]:=BL[j2][j]:
   AV[1][j]:=eval(RB0[BasMolR[j1]][j]): BV[1][j]:=eval(RB0[BasMolR[j2]][j]):
  od:
  w0:=add(-Znuc[i]*add(add(I1evalB(BAPrim[j1][k1],AV,BAPrim[j2][k2],BV,RB0[i])*BCPrim[j1][k1]*BCPrim[j2][k2],k1=1..BNPrim[j1]),k2=1..BNPrim[j2]),i=1..Nnuc):
 Hel0[i1,i2]:=evalf(subs(R=0,K0[j1,j2])+w0): 
 Sel0[i1,i2]:=evalf(subs(R=0,S0[j1,j2])):
 od:

if(Degen[j1]>0) then
 j2:=Degen[j1]:
 for j from 1 to 3 do
  AV[2][j]:=BL[j1][j]: BV[2][j]:=BL[j2][j]: 
  AV[1][j]:=eval(RBi[BasMolR[j1]][j]): BV[1][j]:=eval(RBi[BasMolR[j2]][j]):
 od:
 w0:=add(-Znuc[i]*add(add(I1evalB(BAPrim[j1][k1],AV,BAPrim[j2][k2],BV,RBi[i])*BCPrim[j1][k1]*BCPrim[j2][k2],k1=1..BNPrim[j1]),k2=1..BNPrim[j2]),i=1..Nnuc):
 Hel1[i1]:=evalf(subs(R=Ri,K0[j1,j2])+w0): 
 Sel1[i1]:=evalf(subs(R=Ri,S0[j1,j2])):
 #print("DB6a",AV,BV,RBi);
 #print("DB6",i1,j1,j2,Sel1[i1],Hel1[i1],Ri);
fi:

od:

if(iPrint>1) then
 print("Dimension",N0,S2L):
 print("Overlap"):
 for i from 1 to N0 do 
  for j from 1 to N0 do printf("%10.6f ",Sel0[i,j]): od: printf("\n");
 od:
 print("Hamiltonian"):
 for i from 1 to N0 do 
  for j from 1 to N0 do printf("%10.6f ",Hel0[i,j]): od: printf("\n");
 od:
fi:

EigSort(N0,Hel0,Sel0,E0,C0):
GS(N0,N0,C0,Sel0): 

#Map reduced rank solutions in C0 back into C, splitting degenerate populations indo two, and adding antisymmetric solutions at the end
for i from 1 to N0 do
 for j from 1 to N0 do
  j1:=S2L[j]: if(Degen[j1]>0) then C[i,j1]:=C0[i,j]/2.: C[i,Degen[j1]]:=C0[i,j]/2.: else C[i,j1]:=C0[i,j]: fi:
 od:
E[i]:=E0[i]:
if(iPrint>2) then
 print("Old orb ",i); 
 for j from 1 to N0 do printf("%10.6f ",C0[i,j]): od: printf("\n");
 print("New orb ",i); 
 for j from 1 to NBas do printf("%10.6f ",C[i,j]): od: printf("\n");
fi:
od:
i1:=1:
for i from N0+1 to NBas do
 i2:=0:
 j1:=S2L[i1]: 
 if(Degen[j1]>0 and i2=0) then 
  j2:=Degen[j1]: w0:=evalf((1-Sel1[i1])/Ri^2): w1:=evalf((Hel0[i1,i1]-Hel1[i1])/Ri^2):
  if(iPrint>1) then print("Pair",j1,j2,"dS,dH=",w0,w1); fi:
  E[i]:=evalf(-w1/w0): C[i,j1]:=1./sqrt(2*w0)/R: C[i,j2]:=-C[i,j1]: i2:=1:
 fi:
 i1:=eval(i1+1):
od:
SortV(NBas,E,C):

if(iPrint>0) then 
for i from 1 to NBas do 
 printf("Orb %d E=%10.6f\n",i,E[i]);
 for j from 1 to NBas do if(type(C[i,j],float)) then printf("%10.6f ",C[i,j]): else printf("%10.6f/R ",evalf(simplify(C[i,j]*R))): fi: od: printf("\n");
od:
fi:  

end proc:
SolRZeroOld:=proc(R,Nnuc,Znuc,NBas,BAPrim,BCPrim,BNPrim,BasMolR,S0,K0,E,C):
#A limit of electronic Hamiltonian for R=0
#For the pairs of basis functions that become degenerate at R=0, the sub does a (non-unitary) transform into g and u pair
#In the new basis the H and S are non-degenerate.  Then the solutins are transformed back.
#Number of nuclei Nnuc, charges in Znuc, positions (in variable R) in RA
#Elec basis function parameters for NBas functions in BAPrim,BCPrim,BNPrim,BasMolR
#Anayltic basis KE and Overlaps S0,K0
local AV,BV,RBi,w0,w1,i,j,i1,i2,j1,j2,tiny,small,U,Um,T1,T2,Hel1,Sel1,Hel0,Sel0,Degen,Ri,iPrint:
tiny:=evalf(10^(-Digits)): if(tiny<1.e-16) then tiny=1.e-16: fi: small:=tiny^(1/4): Ri:=evalf(small*100): Degen:=Vector(NBas): 
U:=Matrix(NBas,NBas): Um:=Matrix(NBas,NBas): #Transformation matrix and its inverse
Hel0:=Matrix(NBas,NBas,shape=symmetric): Sel0:=Matrix(NBas,NBas,shape=symmetric):
T1:=Matrix(NBas,NBas): T2:=Matrix(NBas,NBas):
#Hel1:=Matrix(NBas,NBas,shape=symmetric): Sel1:=Matrix(NBas,NBas,shape=symmetric):
iPrint:=0:

for i from 1 to NBas do Degen[i]:=0: od:
for i from 1 to NBas do
 if (Degen[i]<1) then U[i,i]:=1: Um[i,i]:=1: fi:
 for j from i+1 to NBas do
  w0:=evalf(subs(R=0,S0[i,j])):
  w1:=evalf(subs(R=Ri,S0[i,j])):
  Sel0[i,j]:=w1:
  if(iPrint>2) then print("Ovr",i,j,w0,evalb(w0^2>1.-small)); fi:
  if(w0^2>1.-small) then
   Degen[i]:=j: Degen[j]:=i:
   w1:=evalf((1-w1)/Ri^2):
   U[i,i]:=1/2: U[i,j]:=1/sqrt(2*w1)/R: U[j,i]:=1/2: U[j,j]:=-1/sqrt(2*w1)/R: 
   Um[i,i]:=1: Um[i,j]:=1: Um[j,i]:=sqrt(w1/2)*R: Um[j,j]:=-sqrt(w1/2)*R:
  fi:
 od:
 Sel0[i,i]:=evalf(subs(R=Ri,S0[i,i])):
od:
if(iPrint>1) then printf(`Degen: `): for i from 1 to NBas do printf(`%d `,Degen[i]); od: printf(`\n`); fi:

AV:=Vector(2): AV[1]:=Vector(3): AV[2]:=Vector(3):
BV:=Vector(2): BV[1]:=Vector(3): BV[2]:=Vector(3): 
RBi:=Array(1..Nnuc):
for i from 1 to Nnuc do RBi[i]:=subs(R=Ri,RA[i]): od:

for i1 from 1 to NBas do
 for i2 from 1 to i1 do

  for j from 1 to 3 do 
   AV[2][j]:=BL[i1][j]: BV[2][j]:=BL[i2][j]:
   AV[1][j]:=eval(RBi[BasMolR[i1]][j]): BV[1][j]:=eval(RBi[BasMolR[i2]][j]):
  od:
  w0:=add(-Znuc[i]*add(add(I1evalB(BAPrim[i1][k1],AV,BAPrim[i2][k2],BV,RBi[i])*BCPrim[i1][k1]*BCPrim[i2][k2],k1=1..BNPrim[i1]),k2=1..BNPrim[i2]),i=1..Nnuc):
 Hel0[i1,i2]:=evalf(subs(R=Ri,K0[i1,i2])+w0): 
 od:
od:

T1:=evalf(subs(R=Ri,U)):
if(iPrint>2) then
 print("Transformation matrix for R=",Ri):
 for i from 1 to NBas do 
  for j from 1 to NBas do printf("%10.6f ",T1[i,j]): od: printf("\n");
 od:
 print("Original Overlap"):
 for i from 1 to NBas do 
  for j from 1 to NBas do printf("%10.6f ",Sel0[i,j]): od: printf("\n");
 od:
fi:

T2:=Transpose(T1).Sel0.T1: Sel0:=evalf(T2):
T2:=Transpose(T1).Hel0.T1: Hel0:=evalf(T2):

if(iPrint>2) then
 print("Transformed Overlap"):
 for i from 1 to NBas do 
  for j from 1 to NBas do printf("%10.6f ",Sel0[i,j]): od: printf("\n");
 od:
 print("Transformed Hamiltonian"):
 for i from 1 to NBas do 
  for j from 1 to NBas do printf("%10.6f ",Hel0[i,j]): od: printf("\n");
 od:
fi:

EigSort(NBas,Hel0,Sel0,E,C):
GS(NBas,NBas,C,Sel0):
for i from 1 to NBas do for j from 1 to NBas do if (C[i,j]*C[i,j]<tiny) then C[i,j]:=0: fi: od: od:

#print("VR1",C);
T1:=C.Transpose(U):
#print("VR2",T1);
for i from 1 to NBas do for j from 1 to NBas do C[i,j]:=evalf(T1[i,j]): od: od:

if(iPrint>0) then 
for i from 1 to NBas do 
 printf("Orb %d E=%10.6f\n",i,E[i]);
 for j from 1 to NBas do if(type(C[i,j],float)) then printf("%10.6f ",C[i,j]): else printf("%10.6f/R ",evalf(simplify(C[i,j]*R))): fi: od: printf("\n");
od:
fi:  

end proc:
Udegen:=proc(N,Di,Ds,ForwardFlag,T0,TR)
#subroutine to compute a non-unitatry transformation to remove degeneracy from the overlap matrix, when R->0
#returns either direct (ForwardFlag=1) or reciprocal (ForwardFlag=0) trasformation in T0 and TR
#T0 has R-independent terms, TR should be multiplied by R (for reciprocal) or 1/R (for direct)
# Di contains indices of degenerate orbitals, and Ds has quadratic coefficiens sigma for S[i,j]=1-sigma*R^2
local i,j,w0:
for i from 1 to N do for j from 1 to i-1 do T0[i,j]:=0: T0[j,i]:=0: TR[i,j]:=0: TR[j,i]:=0: od: T0[i,i]:=1: TR[i,i]:=0: od:
for i from 1 to N do
 if(Di[i]>i) then
  j:=Di[i]: 
  if(ForwardFlag=1) then
   T0[i,i]:=1/2: TR[i,j]:=evalf(1/sqrt(2*Ds[i])): T0[j,i]:=1/2: T0[j,j]:=0: TR[j,j]:=-TR[i,j]:
  else
   T0[i,i]:=1: T0[i,j]:=1: T0[j,j]:=0: TR[j,i]:=evalf(sqrt(Ds[i]/2)): TR[j,j]:=-TR[j,i]:
  fi:
 fi:
od:
end proc: 

   U[i,i]:=1/2: U[i,j]:=1/sqrt(2*w1)/R: U[j,i]:=1/2: U[j,j]:=-1/sqrt(2*w1)/R: 
   Um[i,i]:=1: Um[i,j]:=1: Um[j,i]:=sqrt(w1/2)*R: Um[j,j]:=-sqrt(w1/2)*R:


SolRZeroT:=proc(R,Nnuc,Znuc,NBas,BAPrim,BCPrim,BNPrim,BasMolR,S0,K0,S,H,H2,do2e):
#A limit of electronic Hamiltonian for R=0
#For the pairs of basis functions that become degenerate at R=0, the sub does a (non-unitary) transform into g(1/2+1/2) and u[1/(2R)-1/(2R)] pair
#In the new basis the H and S are non-degenerate.  Returns a transformation matrix T,transformed S and H, and 2-e integrals H2 if do2e=1
#Number of nuclei Nnuc, charges in Znuc, positions (in variable R) in RA
#Elec basis function parameters for NBas functions in BAPrim,BCPrim,BNPrim,BasMolR
#Anayltic basis KE and Overlaps S0,K0
local AV,BV,RBi,RB0,w0,w1,w2,w3,i,j,i1,i2,j1,j2,tiny,small,U0,UR,Um,T0,T1,T2,T3,Ri,iPrint:
global Degen,DegenS:
tiny:=evalf(10^(-Digits)): if(tiny<1.e-16) then tiny:=1.e-16: fi: small:=tiny^(1/4): Ri:=evalf(small): 
Degen:=Vector(NBas): DegenS:=Vector(NBas): 
U0:=Matrix(NBas,NBas): UR:=Matrix(NBas,NBas): #Um:=Matrix(NBas,NBas): #Transformation matrix and its inverse
# U0 has r-independent value at R->0, UR needs to be multiplied by 1/R
#Hel0:=Matrix(NBas,NBas,shape=symmetric): Sel0:=Matrix(NBas,NBas,shape=symmetric):
T0:=Matrix(NBas,NBas): T1:=Matrix(NBas,NBas): T2:=Matrix(NBas,NBas): T3:=Matrix(NBas,NBas):
#Hel1:=Matrix(NBas,NBas,shape=symmetric): Sel1:=Matrix(NBas,NBas,shape=symmetric):
iPrint:=1:
for i from 1 to NBas do Degen[i]:=0: DegenS[i]:=0: od:
for i from 1 to NBas do for j from 1 to NBas do T0[i,j]:=0: T1[i,j]:=0: T2[i,j]:=0: od: od:
for i from 1 to NBas do
 for j from i+1 to NBas do
  w0:=evalf(subs(R=0,S0[i,j])):
  w1:=evalf(subs(R=Ri,S0[i,j])):
  T0[i,j]:=w0: T0[j,i]:=w0: 
  if(iPrint>2) then print("Ovr",i,j,w0,evalb(w0^2>1.-small)); fi:
  if(w0^2>1.-small) then
   Degen[i]:=j: Degen[j]:=i:
   w2:=evalf((1-w1)/Ri^2):
   DegenS[i]:=w2: DegenS[j]:=w2: 
   T0[i,j]:=1: T0[j,i]:=1: T0[i,i]:=1: T0[j,j]:=1: 
   T2[i,j]:=-w2: T2[j,i]:=-w2:  
  else
   w2:=evalf((w0-w1)/Ri): 
   if(w2^2<small^2) then
    T2[i,j]:=-w2/Ri; T2[j,i]:=-w2/Ri; # at least quadratic dependence
   else 
    T1[i,j]:=-w2: T1[j,i]:=-w2: #linear dependence
   fi:
  fi:
 od:
 T0[i,i]:=1:
od:
if(iPrint>1) then printf(`Degen: `): for i from 1 to NBas do printf(`%d `,Degen[i]); od: printf(`\n`); fi:
Udegen(NBas,Degen,DegenS,1,U0,UR):
#print("DB VR U0,UR",U0,UR);
#print("DB VR T",T0,T1,T2);

#T0 contains R-independent overlap, T1 has O(R), and T2 has O(R^2) terms.  We apply the U transformation and analyze R-dependence
T3:=Transpose(UR).T1.UR: 
w3:=evalf(MatrixNorm(T3,Frobenius)/NBas): 
 if(w3>tiny) then
  print("Potential Error in O(R) terms of the overlap",w3);
 fi:
T3:=Transpose(UR).T0.UR: w3:=evalf(MatrixNorm(T3,Frobenius)/NBas): 
T3:=Transpose(UR).T0.U0: w2:=evalf(MatrixNorm(T3,Frobenius)/NBas):
 if(w3>tiny) then
  print("Potential Error in O(1) terms of the overlap",w3);
 fi:
w0:=Transpose(U0).T0.U0+Transpose(UR).T1.U0+Transpose(U0).T1.UR+Transpose(UR).T2.UR:
for i from 1 to NBas do for j from 1 to i do S[i,j]:=(w0[i,j]+w0[j,i])/2: S[j,i]:=(w0[i,j]+w0[j,i])/2: od: od:



AV:=Vector(2): AV[1]:=Vector(3): AV[2]:=Vector(3):
BV:=Vector(2): BV[1]:=Vector(3): BV[2]:=Vector(3): 
RBi:=Array(1..Nnuc): RB0:=Array(1..Nnuc):
for i from 1 to Nnuc do RBi[i]:=subs(R=Ri,RA[i]): RB0[i]:=subs(R=0,RA[i]): od:

for i1 from 1 to NBas do
 for i2 from 1 to i1 do

  for j from 1 to 3 do 
   AV[2][j]:=BL[i1][j]: BV[2][j]:=BL[i2][j]:
   AV[1][j]:=eval(RBi[BasMolR[i1]][j]): BV[1][j]:=eval(RBi[BasMolR[i2]][j]):
  od:
  w1:=add(-Znuc[i]*add(add(I1evalB(BAPrim[i1][k1],AV,BAPrim[i2][k2],BV,RBi[i])*BCPrim[i1][k1]*BCPrim[i2][k2],k1=1..BNPrim[i1]),k2=1..BNPrim[i2]),i=1..Nnuc):
  for j from 1 to 3 do 
   #AV[2][j]:=BL[i1][j]: BV[2][j]:=BL[i2][j]:
   AV[1][j]:=eval(RB0[BasMolR[i1]][j]): BV[1][j]:=eval(RB0[BasMolR[i2]][j]):
  od:
  w0:=add(-Znuc[i]*add(add(I1evalB(BAPrim[i1][k1],AV,BAPrim[i2][k2],BV,RB0[i])*BCPrim[i1][k1]*BCPrim[i2][k2],k1=1..BNPrim[i1]),k2=1..BNPrim[i2]),i=1..Nnuc):
  w1:=evalf(w1+subs(R=Ri,K0[i1,i2])):
  w0:=evalf(w0+subs(R=0,K0[i1,i2])):
  T0[i1,i2]:=w0: T0[i2,i1]:=w0:
  if(Degen[i1]=i2) then
    w2:=evalf((w0-w1)/Ri^2): T2[i1,i2]:=-w2: T2[i2,i1]:=-w2:  T1[i1,i2]:=0: T1[i2,i1]:=0:
  else
    w2:=evalf((w0-w1)/Ri):
    if(w2^2<small) then
      T2[i1,i2]:=-w2/Ri: T2[i2,i1]:=-w2/Ri:  T1[i1,i2]:=0: T1[i2,i1]:=0:
    else
      T1[i1,i2]:=-w2: T1[i2,i1]:=-w2:  T2[i1,i2]:=0: T2[i2,i1]:=0:
    fi:
  fi:
 od:
od:

#print("DB VR H",T0,T1,T2);

T3:=Transpose(UR).T1.UR: w3:=evalf(MatrixNorm(T3,Frobenius)/NBas): 
 if(w3>tiny) then
  print("Potential Error in O(R) terms of H",w3);
 fi:
T3:=Transpose(UR).T0.UR: w3:=evalf(MatrixNorm(T3,Frobenius)/NBas): 
T3:=Transpose(UR).T0.U0: w2:=evalf(MatrixNorm(T3,Frobenius)/NBas):
 if(w3>tiny) then
  print("Potential Error in O(1) terms of H",w3);
 fi:
w0:=Transpose(U0).T0.U0+Transpose(UR).T1.U0+Transpose(U0).T1.UR+Transpose(UR).T2.UR:
for i from 1 to NBas do for j from 1 to i do H[i,j]:=(w0[i,j]+w0[j,i])/2: H[j,i]:=(w0[i,j]+w0[j,i])/2: od: od:

if(iPrint>2) then
 print("Transformed Overlap"):
 for i from 1 to NBas do 
  for j from 1 to NBas do printf("%10.6f ",S[i,j]): od: printf("\n");
 od:
 print("Transformed Hamiltonian"):
 for i from 1 to NBas do 
  for j from 1 to NBas do printf("%10.6f ",H[i,j]): od: printf("\n");
 od:
fi:

end proc:
LPFG:=proc(NTraj,Y,WT,y,NFLTR,CFLTR,FLTRbas,F,LPFErr)
#Low pass filter for a generic basis
#Y is the set of points (can be random) for which the projection is done
#FLTRbas is a vector of size NFLTR containing the basis functions (in y) which the filter will be performed
#CFLTR is the vecotr of fit coefficients
#Outputs function F projected onto bf and computes cumulative change (F is used to compute error, its values will be replaced but it should be the vector containing the actual fxn for error analysis)
#LPFErr tracks accumulated errror of the filter (weighted with trajectory weights) 
local j,w0,m,w1:

for j from 1 to NTraj do 
  w0:=evalf(add(CFLTR[m]*subs(y=Y[j],FLTRbas[m]),m=1..NFLTR)):
  w1:=LPFErr+sqrt(abs((F[j]-w0))^2)*WT[j]:
  F[j]:=w0:
od:

w1:

end proc:
LPFFT:=proc(NTraj,Nbas,NFr,dw,dshft,Traj,LPFErr)
#Low pass filter for removing high frequency noise from crossing trajectories
#Outputs coefficients projected onto fourier and computes cumulative change
#LPFErr tracks accumulated errror of the filter (weighted with trajectory weights) 
local CF,RV,PF,i,j,w0,m,w1:

PF:=Array(0..2*NFr,1..Nbas): RV:=Vector(NTraj): CF:=Matrix(Nbas,NTraj):

for i from 1 to NTraj do RV[i]:=Traj[1][i]: od:
for j from 1 to NTraj do for i from 1 to Nbas do CF[i,j]:=Traj[5][j][i]: od: od:

FTProj2(NTraj,Nbas,Traj[4],RV,NFr,dw,dshft,CF,PF):
#for j from 1 to Nbas do for i from 1 to 2*NFr+1 do print("CF,BF,i",j,i,PF[i-1,j]); od:od:

##print(PF[2,1],PF[1,2],PF[1,3]);

for j from 1 to NTraj do 
 for i from 1 to Nbas do 
  w0:=evalf(add(PF[m,i]*cos((m)*(dw+dwshft)*RV[j]),m=0..NFr)+add(PF[m+NFr,i]*sin((m)*(dw+dwshft)*RV[j]),m=1..NFr)):
  w1:=LPFErr+sqrt(abs((Traj[5][j][i]-w0))^2)*Traj[4][j]:
  Traj[5][j][i]:=w0:
od:od:

w1:

end proc:
