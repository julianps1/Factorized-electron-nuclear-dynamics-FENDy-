
FFTinit:=proc(m)
description "Initialization of bit-reversal index";
local N,i:
global IndRev:
N:=2^m:
print("Initilaizing of IndRev of size",N,". Counting starts at 1");
IndRev:=Vector(N):
for i from 1 to N do
 IndRev[i]:=Bits:-Join(ListTools:-Reverse(Bits:-Split(i-1,bits=m)))+1:
od:
end proc:

FFTAcc:=proc(m,Z,forw)
description "Numerically accurate FFT. for=1 is FFT, for=0 is inverse FFT";
local sig,i,j,k,s,tz,uz,Z1,N,m2,wm,w:
N:=2^m:
Z1:=Vector(N):
if forw=1 then sig:=1: elif forw=0 then sig:=-1: else print("Error.  The third element should be 1 or zero"); return; fi: 

#bit-reverse copy arrays
for i from 1 to N do Z1[IndRev[i]]:=evalf(Z[i]): od:


#FFT copied from https://en.wikipedia.org/wiki/Cooley-Tukey_FFT_algorithm 
for s from 1 to m do
 m2:=2^(s-1):
 wm:=evalf(exp(-sig*Pi*I/m2)):
 for k from 1 to N by m2*2 do
  w:=1:
  for j from 1 to m2 do
   tz:=w*Z1[k+j+m2-1]:
   uz:=Z1[k+j-1]:
   #print("s,k,j,w,i1,i2,F[i1],F[i2]",s,k,j,w,k+j+m2-1,k+j-1,tz,uz,uz+tz,uz-tz);
   Z1[k+j-1]:=evalf(uz+tz):
   Z1[k+j+m2-1]:=evalf(uz-tz):
   w:=evalf(w*wm):
od: od: od:

w:=sqrt(N):
for i from 1 to N do Z[i]:=evalf(Z1[i]/w) od:
end proc:

FFTDer:=proc(m,Z,dx,der)
# Computation of first (der=1) or second (der=2) derivatives of Z via FFT.  
# number of points is 2^m, dx is spacing between the points, Z is overwritten with the derivative
local w,i,N1:
N1:=2^(m-1):
FFTAcc(m,Z,1):
if (der=1) then
 w:=2*Pi*I/dx/2/N1:
 for i from 0 to N1-1 do Z[i+1]:=evalf(Z[i+1]*w*i): od:
 Z[N1+1]:=0: 
for i from N1+1 to 2*N1-1 do Z[i+1]:=evalf(Z[i+1]*w*(i-2*N1)): od:
elif (der=2) then
 w:=-(2*Pi/dx/2/N1)^2:
 for i from 0 to N1-1 do Z[i+1]:=evalf(Z[i+1]*w*i*i): od:
 for i from N1 to 2*N1-1 do Z[i+1]:=evalf(Z[i+1]*w*(i-2*N1)*(i-2*N1)): od:
fi:
FFTAcc(m,Z,0):
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
FTProj:=proc(N,M,WT,X,nw,dw,F,PF)
# subroutine for projecting functions F onto sin and cos polynomials
# M complex functions are on the 1..N grid, in F[1..M,1..N].  Positions are in X[1..N], weights are in WT[1..N]
# There are nw+1 cos(w_i*x) and nw sin(w_i*x) functions, w_i start from 0 (for cos) or dw (for sin), and increase by dw
# return coefficients are in PF[0..2*nw,1..N], first 0..nw are for cos, then nw+1..2*nw for sin
local i,j,k,n,S,Sm,A,B,wi,wj:
global Digits,Acc:
k:=2*nw+1:
S:=Matrix(k,k,shape=symmetric): A:=Vector(k): B:=Vector(k):
for i from 0 to nw do
 wi:=i*dw:
 for j from 0 to i do
  wj:=j*dw:
   S[i+1,j+1]:=evalf(add(cos(wi*X[n])*cos(wj*X[n])*WT[n],n=1..N)):
 od:
od:
for i from 1 to nw do
 wi:=i*dw:
 for j from 0 to nw do
  wj:=j*dw:
  S[i+nw+1,j+1]:=evalf(add(sin(wi*X[n])*cos(wj*X[n])*WT[n],n=1..N)):
 od:
 for j from 1 to i do
  wj:=j*dw:
  S[i+nw+1,j+nw+1]:=evalf(add(sin(wi*X[n])*sin(wj*X[n])*WT[n],n=1..N)):
 od:
od:

if abs(LinearAlgebra[Determinant](S))<evalf(10^(-Digits/2)) then 
 print("WARNING IN FTProj: Small Determinant of S", Determinant(S));
fi:

Sm:=Matrix(k,k):

#Sm:=LinearAlgebra[MatrixInverse](S):
RefinedMatrixInverse(k,S,Sm,Acc):


for i from 1 to M do
 for j from 0 to nw do
  wj:=j*dw:
  A[j+1]:=evalf(add(cos(wj*X[n])*WT[n]*F[i,n],n=1..N)):
 od:
 for j from 1 to nw do
  wj:=j*dw:
  A[j+nw+1]:=evalf(add(sin(wj*X[n])*WT[n]*F[i,n],n=1..N)):
od:
B:=Sm.A:
for j from 1 to 2*nw+1 do PF[j-1,i]:=evalf(B[j]): od:
od:

end proc:

FTDer:=proc(N,M,Nder,X,nw,dw,F,PF)
# subroutine to compute Nder-th derivative of functions, stored as harmonic expansion in PF, and put it into F
# Needs FTProj to run first to compute PF, for the meanings of variables see FTProj
local i,j,k,wi,r,ff1,ff2,fx1,fx2:
unassign(`r`):
for i from 1 to M do for j from 1 to N do F[i,j]:=0: od: od:
for k from 1 to nw do
 wi:=k*dw:
 if(type(Nder,even)) then 
  ff1:=wi^Nder*(-1)^(Nder/2)*cos(wi*r): 
  ff2:=wi^Nder*(-1)^(Nder/2)*sin(wi*r): 
   else 
  ff1:=wi^Nder*(-1)^((Nder+1)/2)*sin(wi*r):
  ff2:=wi^Nder*(-1)^((Nder+3)/2)*cos(wi*r): 
 fi:
 ##if(k=1) then print("DB",ff1,ff2); fi:
 for j from 1 to N do 
  fx1:=subs(r=X[j],ff1): fx2:=subs(r=X[j],ff2): 
  for i from 1 to M do F[i,j]:=evalf(F[i,j]+PF[k,i]*fx1+PF[k+nw,i]*fx2): od:
 od:
od:
if(Nder=0) then
 for j from 1 to N do
  for i from 1 to M do F[i,j]:=evalf(F[i,j]+PF[0,i]): od:
 od:
fi:

end proc:

FTProj2:=proc(N,M,WT,X,nw,dw,dwshft,F,PF)
# subroutine for projecting functions F onto sin and cos polynomials
# M complex functions are on the 1..N grid, in F[1..M,1..N].  Positions are in X[1..N], weights are in WT[1..N]
# There are nw+1 cos(w_i*x) and nw sin(w_i*x) functions, w_i start from 0 (for cos) or dw (for sin), and increase by dw*dwshft
# dwshft of 0 gives integer spacing (multiples of dw for different functions) or integer+dwshft (negative dshift will decrease spacing)
# return coefficients are in PF[0..2*nw,1..N], first 0..nw are for cos, then nw+1..2*nw for sin
local i,j,k,n,S,Sm,A,B,wi,wj,freq:
global Digits:

k:=2*nw+1:
freq:=Vector(nw):
S:=Matrix(k,k,shape=symmetric): A:=Vector(k): B:=Vector(k):
 
   S[1,1]:=evalf(add(WT[n],n=1..N)): #CONSTANT
for i from 1 to nw do
 wi:=(i)*(dw+dwshft):
 freq[i]:=wi:
   S[i+1,1]:=evalf(add(cos(wi*X[n])*WT[n],n=1..N)): #CONSTANT WITH OTHER COSINES
 for j from 1 to i do
  wj:=(j)*(dw+dwshft):
   S[i+1,j+1]:=evalf(add(cos(wi*X[n])*cos(wj*X[n])*WT[n],n=1..N)):
 od:
od:
##print("Cosine frequencies",freq):

for i from 1 to nw do
 wi:=(i)*(dw+dwshft):
 freq[i]:=wi:
 S[i+nw+1,1]:=evalf(add(sin(wi*X[n])*WT[n],n=1..N)): #CONSTANT COS(0)
 for j from 1 to nw do
  wj:=(j)*(dw+dwshft):
  S[i+nw+1,j+1]:=evalf(add(sin(wi*X[n])*cos(wj*X[n])*WT[n],n=1..N)):
 od:
 for j from 1 to i do
  wj:=(j)*(dw+dwshft):
  S[i+nw+1,j+nw+1]:=evalf(add(sin(wi*X[n])*sin(wj*X[n])*WT[n],n=1..N)):
 od:
od:
##print("Sine frequencies",freq):

if abs(LinearAlgebra[Determinant](S))<evalf(10^(-Digits/2)) then 
 print("WARNING IN FTProj: Small Determinant of S", Determinant(S));
fi:

Sm:=Matrix(k,k):

#Sm:=LinearAlgebra[MatrixInverse](S):
RefinedMatrixInverse(k,S,Sm,Acc):

for i from 1 to M do
  A[1]:=evalf(add(WT[n]*F[i,n],n=1..N)): #CONSTANT
 for j from 1 to nw do
  wj:=(j)*(dw+dwshft):
  freq[j]:=wj:
  A[j+1]:=evalf(add(cos(wj*X[n])*WT[n]*F[i,n],n=1..N)):
 od:
## print("Cosine frequencies for overlap vector",freq):
 for j from 1 to nw do
  wj:=(j)*(dw+dwshft):
  freq[j]:=wj:
  A[j+nw+1]:=evalf(add(sin(wj*X[n])*WT[n]*F[i,n],n=1..N)):
od:
## print("Sine frequencies for overlap vector",freq):
B:=Sm.A:
for j from 1 to 2*nw+1 do PF[j-1,i]:=evalf(B[j]): od:
od:

end proc:

FTDer2:=proc(N,M,Nder,X,nw,dw,dwshft,F,PF)
# subroutine to compute Nder-th derivative of functions, stored as harmonic expansion in PF, and put it into F
# Needs FTProj to run first to compute PF, for the meanings of variables see FTProj
local i,j,k,wi,r,ff1,ff2,fx1,fx2:
unassign(`r`):
for i from 1 to M do for j from 1 to N do F[i,j]:=0: od: od:
for k from 1 to nw do
 wi:=k*(dw+dwshft):
 if(type(Nder,even)) then 
  ff1:=wi^Nder*(-1)^(Nder/2)*cos(wi*r): 
  ff2:=wi^Nder*(-1)^(Nder/2)*sin(wi*r): 
   else 
  ff1:=wi^Nder*(-1)^((Nder+1)/2)*sin(wi*r):
  ff2:=wi^Nder*(-1)^((Nder+3)/2)*cos(wi*r): 
 fi:
 ##if(k=1) then print("DB",ff1,ff2); fi:
 for j from 1 to N do 
  fx1:=subs(r=X[j],ff1): fx2:=subs(r=X[j],ff2): 
  for i from 1 to M do F[i,j]:=evalf(F[i,j]+PF[k,i]*fx1+PF[k+nw,i]*fx2): od:
 od:
od:
if(Nder=0) then
 for j from 1 to N do
  for i from 1 to M do F[i,j]:=evalf(F[i,j]+PF[0,i]): od:
 od:
fi:

end proc:

Legen:=proc(imin,imax,C,T,Nmax)
# expansion of C defined on equidistant grid in terms of Legendre Polynomilas mapping [-1,1] onto [imin..imax]
# coeffs are returned in T[i] for polynomials of i-1 order
local i,ii,j,x,x1,w1,N,d,Tx:
N:=imax-imin+1: d:=2./(N-1): #number of points and spacing
Tx:=Vector(Nmax): for i from 1 to Nmax do Tx[i]:=simplify(LegendreP(i-1,x)): od:
for j from 1 to Nmax do
 #x:=-1+2*(i-1)/(N-1)
 T[j]:=evalf((2*j-1)/(N-1)*(add(subs(x=2*(ii-1)/(N-1)-1.,C[imin+ii-1]*Tx[j]),ii=2..N-1)+
                            1/2.*(C[imin]*subs(x=-1,Tx[j])+C[imax]*subs(x=1,Tx[j])))):
od:
end proc:

LegenD:=proc(T,Nmax,Ndiff,Dx)
# using Legendre coefficients T of f(x), with x_max-x_min=Dx, compute T of Ndiff derivative of f
local i,j,k,c,T1:
T1:=Vector(Nmax+1): for i from 1 to Nmax do T1[i]:=evalf(T[i]): od: T1[Nmax+1]:=0:
for k from 1 to Ndiff do
 for i from 1 to Nmax-k do T[i]:=(2*i-1)*add(T1[i+1+2*j],j=0..floor((Nmax-k-i)/2)): od: T[Nmax-k+1]:=0:
 for i from 1 to Nmax-k+1 do T1[i]:=evalf(T[i]): od:
od:
c:=(2/Dx)^Ndiff: for k from 1 to Nmax-Ndiff do T[k]:=evalf(T[k]*c): od:
end proc:

LegenN:=proc(imin,imax,C,T,Nmax)
# Placing Legendre-defined function (Legen. coeffs in T) back to grid on C
local i,j,x,x1,Tx:
Tx:=Vector(Nmax): for i from 1 to Nmax do Tx[i]:=simplify(LegendreP(i-1,x)): od:
for i from imin to imax do
 x1:=-1.+2*(i-imin)/(imax-imin):
 C[i]:=add(T[j]*evalf(subs(x=x1,Tx[j])),j=1..Nmax):
od:
end proc:

DiffLegen:=proc(imin,imax,C,Nmax,Dx,Ndiff)
# computing derivative of C via Legendre expansion
local T:
T:=Vector(Nmax):
Legen(imin,imax,C,T,Nmax):
LegenD(T,Nmax,Ndiff,Dx):
LegenN(imin,imax,C,T,Nmax):
end proc:
SecDerMillColbInf:=proc(imin,imax,dx,Cin,Cout)
# Miller-Colbert DVR based second derivative for equispaced grid, function vanishes at the end (infinite grid)
local i,j:
for i from imin to imax do
 Cout[i]:=evalf(1/dx^2*(-Pi^2/3*Cin[i]+
 2*add(Cin[j]/(j-i)^2,j=i-1..imin,-2)+
 2*add(Cin[j]/(j-i)^2,j=i+1..imax,2)-
 2*add(Cin[j]/(j-i)^2,j=i-2..imin,-2)-
 2*add(Cin[j]/(j-i)^2,j=i+2..imax,2))):
od:
end proc:

SecDerMillColb:=proc(imin,imax,dx,C,Cout)
# Miller-Colbert DVR based second derivative for equispaced grid, finite grid
local i,j,w1,N:
#to vanish at the ends, we subtract linear function C[imin]+(i-imin)/(imax-imin)*(C[imax]-C[imin])
N:=imax-imin: w1:=evalf(-Pi^2/2/(imax-imin)^2/dx^2):
for i from imin+1 to imax-1 do
 Cout[i]:=w1*evalf((C[i]-C[imin]-(i-imin)/(imax-imin)*(C[imax]-C[imin]))*((2*N^2+1.)/3-1./(sin(Pi*(i-imin)/N))^2)+
  add((C[j]-C[imin]-(j-imin)/(imax-imin)*(C[imax]-C[imin]))*
       (1./(sin(Pi*(i-j)/2/N))^2-1./(sin(Pi*(i+j-2*imin)/2/N))^2),j=i-2..imin+1,-2)+
  add((C[j]-C[imin]-(j-imin)/(imax-imin)*(C[imax]-C[imin]))*
       (1./(sin(Pi*(i-j)/2/N))^2-1./(sin(Pi*(i+j-2*imin)/2/N))^2),j=i+2..imax-1,2)-
  add((C[j]-C[imin]-(j-imin)/(imax-imin)*(C[imax]-C[imin]))*
       (1./(sin(Pi*(i-j)/2/N))^2-1./(sin(Pi*(i+j-2*imin)/2/N))^2),j=i-1..imin+1,-2)-
  add((C[j]-C[imin]-(j-imin)/(imax-imin)*(C[imax]-C[imin]))*
       (1./(sin(Pi*(i-j)/2/N))^2-1./(sin(Pi*(i+j-2*imin)/2/N))^2),j=i+1..imax-1,2)):
 od:

 i:=imin: Cout[i]:=w1*evalf(
   add((C[j]-C[imin]-(j-imin)/(imax-imin)*(C[imax]-C[imin]))*
       (1./(sin(Pi*(i-j)/2/N))^2-1./(sin(Pi*(i+j-2*imin)/2/N))^2),j=i-2..imin+1,-2)+
  add((C[j]-C[imin]-(j-imin)/(imax-imin)*(C[imax]-C[imin]))*
       (1./(sin(Pi*(i-j)/2/N))^2-1./(sin(Pi*(i+j-2*imin)/2/N))^2),j=i+2..imax-1,2)-
  add((C[j]-C[imin]-(j-imin)/(imax-imin)*(C[imax]-C[imin]))*
       (1./(sin(Pi*(i-j)/2/N))^2-1./(sin(Pi*(i+j-2*imin)/2/N))^2),j=i-1..imin+1,-2)-
  add((C[j]-C[imin]-(j-imin)/(imax-imin)*(C[imax]-C[imin]))*
       (1./(sin(Pi*(i-j)/2/N))^2-1./(sin(Pi*(i+j-2*imin)/2/N))^2),j=i+1..imax-1,2)):

 i:=imax: Cout[i]:=w1*evalf(
   add((C[j]-C[imin]-(j-imin)/(imax-imin)*(C[imax]-C[imin]))*
       (1./(sin(Pi*(i-j)/2/N))^2-1./(sin(Pi*(i+j-2*imin)/2/N))^2),j=i-2..imin+1,-2)+
  add((C[j]-C[imin]-(j-imin)/(imax-imin)*(C[imax]-C[imin]))*
       (1./(sin(Pi*(i-j)/2/N))^2-1./(sin(Pi*(i+j-2*imin)/2/N))^2),j=i+2..imax-1,2)-
  add((C[j]-C[imin]-(j-imin)/(imax-imin)*(C[imax]-C[imin]))*
       (1./(sin(Pi*(i-j)/2/N))^2-1./(sin(Pi*(i+j-2*imin)/2/N))^2),j=i-1..imin+1,-2)-
  add((C[j]-C[imin]-(j-imin)/(imax-imin)*(C[imax]-C[imin]))*
       (1./(sin(Pi*(i-j)/2/N))^2-1./(sin(Pi*(i+j-2*imin)/2/N))^2),j=i+1..imax-1,2)):

end proc:


