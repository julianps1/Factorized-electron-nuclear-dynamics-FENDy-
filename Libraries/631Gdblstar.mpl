
# Basis set definition
# This is 6-31G**
# Nbas[Z] is an array for the numbers of basis functions for each Z
# Lbas[Z,k] are the orbital angular momentum values
# Mbas[Z,k,m] are Cartesian components. m=1,2,3 corresponds to (x,y,z)
# Bas[Z,k] are the actual functions
# BasN[Z,k] are the correcting normalization constants for each basis function
# Abas[Z,k,1] is the number of primitive Gaussians in this basis function
# Abas[Z,k,2] is the index of the first exponent in the list Alp (minus 1)
# Abas[Z,k,3] is the index of the first contraction coefficient in the list Gs (minus 1)
Debug:=0:
Zmax:=1: # Worksheet is up to Z=1
Nbas:=Array(1..Zmax): Nbas[1]:=5:
lmt:=-infinity..infinity:
Lbas:=Array(1..Zmax,1..5): Lbas[1,1]:=0: Lbas[1,2]:=0: Lbas[1,3]:=1: Lbas[1,4]:=1: Lbas[1,5]:=1:
Mbas:=Array(1..Zmax,1..5,1..3):
for j from 1 to Zmax do
 k:=1: 
 while (k <= Nbas[j]) do
  if (Lbas[j,k] = 0 ) then
   Mbas[j,k,1]:=0: Mbas[j,k,2]:=0: Mbas[j,k,3]:=0:
   k:=eval(k+1):
  elif (Lbas[j,k] = 1 ) then # p-functions
   Mbas[j,k,1]:=0: Mbas[j,k,2]:=0: Mbas[j,k,3]:=1:
   Mbas[j,k+1,1]:=0: Mbas[j,k+1,2]:=1: Mbas[j,k+1,3]:=0:
   Mbas[j,k+2,1]:=1: Mbas[j,k+2,2]:=0: Mbas[j,k+2,3]:=0:
   k:=eval(k+3):
  fi:
 od:
od:
Abas:=Array(1..Zmax,1..5,1..3):
Abas[1,1,1]:=3: Abas[1,1,2]:=0: Abas[1,1,3]:=0:
Abas[1,2,1]:=1: Abas[1,2,2]:=3: Abas[1,2,3]:=3:
Abas[1,3,1]:=1: Abas[1,3,2]:=4: Abas[1,3,3]:=3:
Abas[1,4,1]:=1: Abas[1,4,2]:=4: Abas[1,4,3]:=3:
Abas[1,5,1]:=1: Abas[1,5,2]:=4: Abas[1,5,3]:=3:
Digits:=30: NL:=Array(0..1):  NL[0]:=(2*alpha/Pi)^(3/4): NL[1]:=2*sqrt(alpha)*NL[0]:
Bas:=Array(1..Zmax,1..5): BasN:=Array(1..Zmax,1..5):
Alp:=Vector(5):
Gs:=Vector(4):
 Alp[1]:=0.1873113696D+02: Alp[2]:=0.2825394365D+01:  Alp[3]:=0.6401216923:
 Alp[4]:=0.1612777588: Alp[5]:=0.1100000000D+01:

 Gs[1]:=0.3349460434D-01: Gs[2]:=0.2347269535: Gs[3]:=0.8137573261:
 Gs[4]:=1.0000000: 
 for j from 1 to Zmax do
  for k from 1 to Nbas[j] do
   Bas[j,k]:=add(subs(alpha=Alp[ii+Abas[j,k,2]],NL[Lbas[j,k]]*exp(-alpha*((x-X0)^2+(y-Y0)^2+(z-Z0)^2))*Gs[ii+Abas[j,k,3]])
  *(x-X0)^Mbas[j,k,1]*(y-Y0)^Mbas[j,k,2]*(z-Z0)^Mbas[j,k,3],ii=1..Abas[j,k,1]):
   w0:=simplify(subs([X0=0,Y0=0,Z0=0],Bas[j,k])):
   w1:=int(int(int(w0^2,x=lmt),y=lmt),z=lmt):
   if (Debug=1) then print(j,k,Abas[j,k,1],Abas[j,k,2],Abas[j,k,3],w0,w1); fi:
   BasN[j,k]:=1/sqrt(w1): 
  od: od:
#if (Debug=1) then print(Bas[1,1]); fi:
#w0:=int(int(int(Bas[1,1]^2,x=lmt),y=lmt),z=lmt): BasN[1,1]:=1/sqrt(w0):
NB:=proc(Z)
 # this procedure returns number of basis functions in the basis
 global Nbas;
 Nbas[Z]:
end:
KB:=proc(Z,k,Xnuc)
 # this procedure returns basis function k for nucleus Z at Xnuc coordinates.  Xnuc must be a Vector(3)
 global Bas,BasN:
 local w0:
 w0:=subs([X0=Xnuc[1],Y0=Xnuc[2],Z0=Xnuc[3]],Bas[Z,k])*BasN[Z,k]:
end:
PB:=proc(Z,k,P)
# This procedure returns parameters of the basis function as elements of the vector P
# P[1]   number of primitives
# P[2]   lx value
# P[3]   ly value
# P[4]   lz value
# P[5]   correcting normalization constant
# P[6..]  values of exponents
# P[6+PB[1]..]  values of contraction coefficients
global Bas,BasN,Abas,Mbas;
local j:
P[1]:=Abas[Z,k,1]:
P[2]:=Mbas[Z,k,1]:
P[3]:=Mbas[Z,k,2]:
P[4]:=Mbas[Z,k,3]:
P[5]:=BasN[Z,k]:
for j from 1 to Abas[Z,k,1] do
 P[5+j]:=Alp[j+Abas[Z,k,2]]:
 P[5+Abas[Z,k,1]+j]:=Gs[j+Abas[Z,k,3]]:
od:
end:

#print("Zmax = ",Zmax);
#NB(1);
#BasMol[1]:=KB(1,1,[1.,0,0]);
#BasMol[2]:=KB(1,2,[1.,0,0]);
#BasMol[3]:=KB(1,3,[1.,0,0]);
#BasMol[4]:=KB(1,4,[1.,0,0]);
#BasMol[5]:=KB(1,5,[1.,0,0]);
#NBas:=NB(1);
#for i from 1 to NBas do 
#w1:=int(int(int(BasMol[i]^2,x=lmt),y=lmt),z=lmt); 
#print(w1):
#od:
