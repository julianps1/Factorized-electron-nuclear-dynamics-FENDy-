
##v2.1 Added MatMult and RefinedMatrixInverse procedures at the bottom--JS:10/31/24
;
GS:=proc(N,n,A,S)
 description "Gram-Schmidt orthonormalization":
 # N orbitals, n basis functions.  A is nxn matrix, top N orbitals are orthonormalized
 # assume overlap matrix S is real, orbitals need not be.
 local i,j1,j2,k,x,y,T:
 T:=Vector(N):

 for i from 1 to N do

  y:=add(add(conjugate(A[i,j1])*S[j1,j2]*A[i,j2],j2=1..n),j1=1..n):
  x:=evalf(1./sqrt(Re(y))): for j1 from 1 to n do A[i,j1]:=evalf(x*A[i,j1]): od:

  for k from 1 to i-1 do 
   y:=add(add(conjugate(A[k,j1])*S[j1,j2]*A[i,j2],j2=1..n),j1=1..n):
   T[k]:=evalf(y):
  od:
  x:=evalf(1-Re(add(conjugate(T[j1])*T[j1],j1=1..i-1))):
  
  for j1 from 1 to n do
   y:=(A[i,j1]-add(T[j2]*A[j2,j1],j2=1..i-1)):
   A[i,j1]:=evalf(y/sqrt(x)):
  od:
 od: 

 end proc:
EigSort1:=proc(N,H,E,C)
# procedure to compute sorted eigenvalues (in E) and eigenvectors in C(state,bf) using linalg
local i,i1,i2,Lst,we,wv,w1,Dig1,Ht:
#print("VR0",N,H);
we:=eigenvalues(H): print("VR",Digits,we);
Ht:=matrix(N,N):
Lst:=[seq(we[i1],i1=1..N)]; w1:=sort(Lst,'output=permutation'):
Dig1:=Digits:  
for i from 1 to N do E[i]:=we[w1[i]]: od:
Digits:=eval(Dig1-2): #need to drop accuracy to make sure kernel exists
for i from 1 to N do 
 for i1 from 1 to N do for i2 from 1 to N do Ht[i1,i2]:=evalf(H[i1,i2]-E[i]*S[i1,i2]): od: od: 
 wv:=kernel(Ht): print(i,E[i],Ht,wv); 
 w1:=add(add(wv[1][i1]*wv[1][i2]*S[i1,i2],i1=1..N),i2=1..N): w1:=evalf(1/sqrt(w1)):
 for i1 from 1 to N do C[i,i1]:=evalf(w1*wv[1][i1]): od:
od:
Digits:=eval(Dig1+2):
end proc:

EigSortG:=proc(N,H,E,C)
# procedure to compute sorted eigenvalues (in E) and ground state eigenvector in C(1,N) using linalg
local i,i1,i2,Lst,we,wv,w1,Dig1,Ht:
#print("VR0",N,H);
we:=eigenvalues(H): #print("VR",Digits,we);
Ht:=matrix(N,N):
Lst:=[seq(we[i1],i1=1..N)]; w1:=sort(Lst,'output=permutation'):
Dig1:=Digits:  
for i from 1 to N do E[i]:=we[w1[i]]: od:
Digits:=eval(Dig1-6): if(Digits>15) then Digits:=15: fi: #need to drop accuracy to make sure kernel exists
i:=1: 
 for i1 from 1 to N do for i2 from 1 to i1-1 do Ht[i1,i2]:=evalf(H[i1,i2]): Ht[i2,i1]:=evalf(H[i2,i1]): od: Ht[i1,i1]:=evalf(H[i1,i1]-E[i]): od: 
 wv:=kernel(Ht): #print(Digits,Ht); print(i,E[i],wv); 
 w1:=add(wv[1][i1]*wv[1][i1],i1=1..N): w1:=evalf(1/sqrt(w1)):
 for i1 from 1 to N do C[i,i1]:=evalf(w1*wv[1][i1]): od:
#od:
Digits:=eval(Dig1):
end proc:

EigSort:=proc(N,H,S,E,C)
# procedure to compute sorted eigenvalues (in E) and eigenvectors in C(state,bf) using LinearAlgebra, later planning with iterative refinement
local i,j,i1,i2,Lst,we,Ct,w1,w2:
we:=Eigenvectors(H,S): #print("VR0",H,S); print("VR1",we[1],we[2]);
Lst:=[seq(Re(we[1][i1]),i1=1..N)]; w1:=sort(Lst,'output=permutation'):
Ct:=Vector(N):
for i from 1 to N do 
 E[i]:=Re(we[1][w1[i]]):
 for j from 1 to N do Ct[j]:=we[2][j,w1[i]]: od: w2:=add(add(Ct[i1]*S[i1,i2]*Ct[i2],i1=1..N),i2=1..N): w2:=evalf(1./sqrt(Re(w2))):
 for j from 1 to N do C[i,j]:=evalf(w2*Re(Ct[j])): od:
od:
#print("VR2",C);
end proc:

EigMatch:=proc(N,H,S,E,C,C0)
# procedure to compute eigenvectors of (H,S) and order them according to the max. overlap with C0
local i,j,i1,i2,Iscr,we,Ct,w1,w2,wn,ix,wx,isign:
Iscr:=Vector(N,0): Ct:=Vector(N):
we:=Eigenvectors(H,S): #print("VR_Match",we,C0);
for i from 1 to N do
 for j from 1 to N do Ct[j]:=Re(add(S[j,i2]*we[2][i2,i],i2=1..N)): od: wn:=add(Ct[i1]*we[2][i1,i],i1=1..N): wn:=evalf(1./sqrt(Re(wn))):
 wx:=0.; ix:=0: isign:=1:
 #print("VR01",wn,Ct);
 for j from 1 to N do
  w1:=add(Ct[i1]*C0[j,i1],i1=1..N): w2:=evalf(Re(w1)^2): #print("Vr ",i,j,wx,w2,ix,Iscr);
  if(w2>wx) then wx:=evalf(w2): ix:=eval(j): if(w1<0) then isign:=-1: else isign:=1: fi: fi:
 od:
 if (Iscr[ix]>0) then 
  print(" Error in assigning orbital ",i,ix," is already assigned to ",Iscr[ix]); 
  print(we);
  print(C0);
  print(S);
 fi: Iscr[ix]:=i:
 #print("VR02",ix,wx);
 E[ix]:=Re(we[1][i]):
 for i1 from 1 to N do C[ix,i1]:=evalf(wn*isign*Re(we[2][i1,i])): od:
od:
#print("VR02",we,Iscr,C0,C);
end proc:

EigMatchDegenU:=proc(N,S,E,C,C0,U)
# procedure to align vectors in C and C0, with degenerate orbitals in blocks.
# Similar to EigMatchDegen, bu doesn't modify C and returns the transformation matrix U (C_new=U*C)
local i,j,j1,i0,i1,i2,we,Ct,Cnrm,w1,iPrint,tiny,iOrb,Ovr,Nd,St,At:
iPrint:=0:
tiny:=evalf(10^(-Digits/2.)):  
if tiny<1.e-8 then tiny:=1.e-8: fi: # criterion for degeneracy:  LinearAlgebra package has accuracy no better than 10^(-16)
Ct:=Vector(N):
iOrb:=Vector(N): for i from 1 to N do iOrb[i]:=i: od: #list of orbital assignments

if(iPrint>2) then
print("Reference orbitals");
for j1 from 1 to N do
 printf (`Orb %d \n`,j1);
 for j from 1 to N do printf(`%8.4f `,C0[j1,j]); od: printf(`\n`);
od:
print("Starting orbitals");
for j1 from 1 to N do
 printf (`Orb %d \n`,j1);
 for j from 1 to N do printf(`%8.4f `,C[j1,j]); od: printf(`\n`);
od:
fi:

#Overlap matrix may be relatively poor for C0 (S should be exact for C).  More accurate to "renormalize" C0 with S, store in Cnrm
Cnrm:=Vector(N):
for i from 1 to N do
 w1:=Re(add(add(S[i1,i2]*C0[i,i2]*conjugate(C0[i,i1]),i2=1..N),i1=1..N)): 
 Cnrm[i]:=evalf(1./sqrt(w1)):
od:
if(iPrint>1) then print("New norms for the ref",Cnrm); fi:

#form of a loop without automatic increments.  "i" is controlled in the procedure
i:=1:
do

#find all states degenerate to "i"
for j from i to N do
 if abs(E[i]-E[j])<tiny then i1:=eval(j): fi: 
od:
Nd:=eval(i1-i+1):

#find (i1-i+1) ket states max overlapping with the same number of bra states, write their indices into iOrb
for j1 from i to i1 do
 for j from 1 to N do Ct[j]:=Re(add(S[j,i2]*C[j1,i2],i2=1..N)): od:
 we:=0: i0:=0:

#printf("VR1a: "); for j from 1 to N do printf(`%d `,iOrb[j]); od: printf(`\n`);

 for i2 from j1 to N do 
  w1:=add(Ct[j]*C0[iOrb[i2],j],j=1..N)*Cnrm[iOrb[i2]]: w1:=evalf(Re(w1)^2):
  if(w1>we) then we:=evalf(w1): i0:=eval(i2): fi:
 od:
 #swap iOrb[j1] and iOrb[i0] (could be the same)
 j:=iOrb[i0]: iOrb[i0]:=iOrb[j1]: iOrb[j1]:=j:
od:

#print("VR2");
St:=Matrix(Nd,Nd,shape=identity): At:=Matrix(Nd,Nd):
if(iPrint>2) then print("State ",i,"block of",i1-i+1,"E=",E[i],E[i1]): for j from i to i1 do printf(`(%d,%d)`,j,iOrb[j]); od: printf(`\n`); fi:
for j1 from i to i1 do
  for j from 1 to N do Ct[j]:=Re(add(S[j,i2]*C0[iOrb[j1],i2],i2=1..N)*Cnrm[iOrb[j1]]): od:
  for i2 from i to i1 do At[j1-i+1,i2-i+1]:=add(Ct[j]*C[i2,j],j=1..N): od:
od:

if(Nd>1) then GS(Nd,Nd,At,St): else At[1,1]:=evalf(At[1,1]/abs(At[1,1])): fi:
for j1 from i to i1 do
 i2:=iOrb[j1]:
 for j from 1 to N do U[i2,j]:=0: od:
 for j from i to i1 do U[i2,j]:=At[j1-i+1,j-i+1]: od:
od:

i:=eval(i+Nd):
until i>N: 

if(iPrint>2) then
print("Transformation Matrix");
for j1 from 1 to N do
 for j from 1 to N do printf(`%8.4f `,U[j1,j]); od: printf(`\n`);
od:
fi:

#Final check
if(iPrint>0) then
 Ovr:=0: w1:=Vector(N):
 for i from 1 to N do
  for j from 1 to N do 
   Ct[j]:=add(C[i1,j]*U[i,i1],i1=1..N):
   w1[j]:=add(C0[i,i1]*S[j,i1],i1=1..N)*Cnrm[i]:
  od:
  we:=add(Ct[j]*w1[j],j=1..N):
  if(iPrint>2) then print("Orb",i,iOrb[i],we,Ovr); fi:
  Ovr:=evalf((1.-we)^2+Ovr):
 od:
 printf(`Orb match is %12.5e \n`,Ovr/N); 
fi:

end proc:
EigMatchDegen:=proc(N,H,S,E,C,C0)
# procedure to compute eigenvectors of (H,S) and order them according to the max. overlap with C0
# If orbitals are degenerate, the linear combinations are chosen
local i,j,j1,i0,i1,i2,we,Ct,Cnrm,w1,iPrint,tiny,iOrb,Ovr,Nd,St,At:
iPrint:=0:
tiny:=evalf(10^(-Digits/2.)):  
if tiny<1.e-8 then tiny:=1.e-8: fi: # criterion for degeneracy:  LinearAlgebra package has accuracy no better than 10^(-16)
Ct:=Vector(N):
iOrb:=Vector(N): for i from 1 to N do iOrb[i]:=i: od: #list of orbital assignments
EigSort(N,H,S,E,C): # Get sorted eigenvectors
GS(N,N,C,S): #Additional G-S orthorormalization due to degeneracy errors

if(iPrint>2) then
print("Reference orbitals");
for j1 from 1 to N do
 printf (`Orb %d \n`,j1);
 for j from 1 to N do printf(`%8.4f `,C0[j1,j]); od: printf(`\n`);
od:
print("Starting orbitals");
for j1 from 1 to N do
 printf (`Orb %d \n`,j1);
 for j from 1 to N do printf(`%8.4f `,C[j1,j]); od: printf(`\n`);
od:
fi:

#Overlap matrix may be relatively poor for C0 (S should be exact for C).  More accurate to "renormalize" C0 with S, store in Cnrm
Cnrm:=Vector(N):
for i from 1 to N do
 w1:=Re(add(add(S[i1,i2]*C0[i,i2]*conjugate(C0[i,i1]),i2=1..N),i1=1..N)): 
 Cnrm[i]:=evalf(1./sqrt(w1)):
od:
if(iPrint>1) then print("New norms for the ref",Cnrm); fi:


#form of a loop without automatic increments.  "i" is controlled in the procedure
i:=1:
do

#find all states degenerate to "i"
for j from i to N do
 if abs(E[i]-E[j])<tiny then i1:=eval(j): fi: 
od:
Nd:=eval(i1-i+1):

#print("VR1",i,i1,Nd);

#find (i1-i+1) ket states max overlapping with the same number of bra states, write their indices into iOrb
for j1 from i to i1 do
 for j from 1 to N do Ct[j]:=Re(add(S[j,i2]*C[j1,i2],i2=1..N)): od:
 we:=0: i0:=0:

#printf("VR1a: "); for j from 1 to N do printf(`%d `,iOrb[j]); od: printf(`\n`);

 for i2 from j1 to N do 
  w1:=add(Ct[j]*C0[iOrb[i2],j],j=1..N)*Cnrm[iOrb[i2]]: w1:=evalf(Re(w1)^2):
  if(w1>we) then we:=evalf(w1): i0:=eval(i2): fi:
 od:
 #swap iOrb[j1] and iOrb[i0] (could be the same)
 j:=iOrb[i0]: iOrb[i0]:=iOrb[j1]: iOrb[j1]:=j:
od:

#print("VR2");
St:=Matrix(Nd,Nd,shape=identity): At:=Matrix(Nd,Nd):
if(iPrint>2) then print("State ",i,"block of",i1-i+1,"E=",E[i],E[i1]): for j from i to i1 do printf(`(%d,%d)`,j,iOrb[j]); od: printf(`\n`); fi:
for j1 from i to i1 do
  for j from 1 to N do Ct[j]:=Re(add(S[j,i2]*C0[iOrb[j1],i2],i2=1..N)*Cnrm[iOrb[j1]]): od:
  for i2 from i to i1 do At[j1-i+1,i2-i+1]:=add(Ct[j]*C[i2,j],j=1..N): od:
od:

if(Nd>1) then GS(Nd,Nd,At,St): else At[1,1]:=evalf(At[1,1]/abs(At[1,1])): fi:

if(Nd>1 and 1=2) then
print("before");
for j1 from i to i1 do
 printf (`Orb %d \n`,j1);
 for j from 1 to N do printf(`%8.4f `,C[j1,j]); od: printf(`\n`);
od:
fi:


for j from 1 to N do
 for i2 from 1 to Nd do Ct[i2]:=C[i2+i-1,j]: od:
 for i2 from 1 to Nd do C[i2+i-1,j]:=add(At[i2,i0]*Ct[i0],i0=1..Nd): od:
od:

if(Nd>1 and 1=2) then
print("after");
for j1 from i to i1 do
 printf (`Orb %d \n`,j1);
 for j from 1 to N do printf(`%8.4f `,C[j1,j]); od: printf(`\n`);
od:
fi:

i:=eval(i+Nd):
until i>N: 

if(iPrint>1) then printf("Orb match: "); for j from 1 to N do printf(`%d `,iOrb[j]); od: printf(`\n`); fi:

#print("VR3a Overlap matrix");
#for i from 1 to N do 
# for j from 1 to N do printf(`%8.4f `,S[i,j]); od: printf(`\n`);
#od:
#print("VR3a Orbs ");
#for i from 1 to N do 
# printf (`Orb %d \n`,i);
# for j from 1 to N do printf(`%8.4f `,C[i,j]); od: printf(`\n`);
# for j from 1 to N do printf(`%8.4f `,C0[iOrb[i],j]); od: printf(`\n`);
#od:

#Orbitals are transformed, now they need to be reordered according to iOrb
# first, flipping iOrb
for i from 1 to N do Ct[i]:=iOrb[i]: od:
for i from 1 to N do for j from 1 to N do if Ct[j]=i then iOrb[i]:=j: fi: od: od:
# second, reorder
for j from 1 to N do
 for i from 1 to N do Ct[i]:=C[i,j]: od:
 for i from 1 to N do C[i,j]:=Ct[iOrb[i]]: od:
od:

#Final check
Ovr:=0:
for i from 1 to N do
 we:=add(add(C[i,i1]*C0[i,i2]*S[i1,i2],i1=1..N),i2=1..N)*Cnrm[i]:
 if(iPrint>2) then print("Orb",i,iOrb[i],we,Ovr); fi:
 Ovr:=evalf((1.-we)^2+Ovr):
od:
if(iPrint>0) then printf(`Orb match is %12.5e \n`,Ovr/N); fi:
end proc:

IterSolv:=proc(N,H,S,E,C)
#iterative solution of generalized eigenfunction equation
# assuming starting vectors in C, H and S are Matrices NxN
# Final Accuracy is Acc, max step size is StepMax
#criterion is sum of suares of all off-diagonal elements
local i,j,i1,i2,j1,j2,It,Acc1,G,P,T,Acc,StepMax,MaxIter,w0,w1,FlagPrint:
#G and P are H and S matrices in current vector representation, respectively
MaxIter:=100: FlagPrint:=0:
Acc:=10^(-Digits+3)*N^2: 
StepMax:=1/10:
if(FlagPrint>0) then printf ("Printing in InterSolv is ON.  Accuracy=%12.4e, MaxStep=%8.4f \n",Acc,StepMax); fi:
G:=Matrix(N,N,shape=hermitian): P:=Matrix(N,N,shape=hermitian): T:=Vector(N):
It:=0: Acc1:=Acc*10:
while (It<MaxIter and Acc1>Acc) do
 It:=eval(It+1):
 for i1 from 1 to N do
  for i2 from 1 to N do T[i2]:=evalf(add(S[i2,j1]*C[i1,j1],j1=1..N)): od:
  w0:=add(T[j1]*C[i1,j1],j1=1..N):
  w0:=evalf(1./sqrt(w0)): for j1 from 1 to N do C[i1,j1]:=evalf(w0*C[i1,j1]): T[j1]:=evalf(w0*T[j1]): od: #make sure the vector norm is 1
  for i2 from 1 to i1 do P[i1,i2]:=evalf(add(T[j1]*C[i2,j1],j1=1..N)): od:
  for i2 from 1 to N do T[i2]:=evalf(add(H[i2,j1]*C[i1,j1],j1=1..N)): od:
  for i2 from 1 to i1 do G[i1,i2]:=evalf(add(T[j1]*C[i2,j1],j1=1..N)): od:
  E[i1]:=evalf(G[i1,i1]):
 od:
 Acc1:=evalf(sqrt(add(add(P[i1,i2]^2+G[i1,i2]^2,i2=1..i1-1),i1=1..N))):
 w0:=0: for i1 from 1 to N do 
  for i2 from 1 to i1-1 do w1:=(G[i1,i2]-E[i1]*P[i1,i2])/(E[i1]-E[i2]): if(w0<w1^2) then w0:=evalf(w1^2):  fi: od:
  for i2 from i1+1 to N do w1:=(G[i1,i2]-E[i1]*P[i1,i2])/(E[i1]-E[i2]): if(w0<w1^2) then w0:=evalf(w1^2):  fi: od:
 od:
 if(sqrt(w0)>StepMax) then w0:=evalf(StepMax/sqrt(w0)): else w0:=1: fi:
 if(FlagPrint>0) then printf("It=%5d Step=%8.4f Acc=%12.4e \n",It,w0,Acc1); fi:
 if(FlagPrint>1 or It>70) then print("H"); printf("%10.6f\n",G); print("S"); printf("%10.6f\n",P); print("C"); printf("%10.6f\n",C); fi:
 for i1 from 1 to N do
  for i2 from 1 to N do T[i2]:=add(w0*(G[j1,i2]-E[i2]*P[j1,i2])/(E[i2]-E[j1])*C[j1,i1],j1=1..i2-1)+
                               add(w0*(G[j1,i2]-E[i2]*P[j1,i2])/(E[i2]-E[j1])*C[j1,i1],j1=i2+1..N)+C[i2,i1]: od:
  for i2 from 1 to N do C[i2,i1]:=evalf(T[i2]): od:
 od:
od:
end proc:

SortV:=proc(N,E,C)
#sorting E and C[vec,bf], in place
local i1,i2,Lst,w0,T:
T:=Vector(N):
Lst:=[seq(E[i1],i1=1..N)]; w0:=sort(Lst,'output=permutation'):
for i1 from 1 to N do T[i1]:=E[w0[i1]]: od: for i1 from 1 to N do E[i1]:=evalf(T[i1]): od:
for i2 from 1 to N do
 for i1 from 1 to N do T[i1]:=C[w0[i1],i2]: od: for i1 from 1 to N do C[i1,i2]:=evalf(T[i1]): od:
od:
end proc:

MatchV:=proc(N,S,E,C,C0)
#matching solutions in C to C0, rearranging E as well
local i,i0,j,i1,T,w0,w1,is:
T:=Vector(N): 
for i from 1 to N do
 for j from 1 to N do T[j]:=Re(add(S[j,i1]*C0[i,i1],i1=1..N)): od: 
 w0:=0: i0:=0:
 for j from i to N do w1:=add(T[i1]*C[j,i1],i1=1..N): if(w0<w1*w1) then w0:=w1*w1: i0:=j: if(w1<0) then is:=-1: else is:=1: fi: fi: od:
 if(i0>i) then 
  w0:=E[i0]: E[i0]:=evalf(E[i]): E[i]:=evalf(w0): 
  for j from 1 to N do w0:=C[i0,j]: C[i0,j]:=evalf(C[i,j]): C[i,j]:=evalf(w0*is): od:
 fi:
od:
end proc:

MatMult:=proc(N1,N2,M,A,B,C)
#Explicitly performs the multiplication A.B=C
#N1,N2 are the dimensions of C(and correspond to dimensions of A and B)
#M is the common dimension between A and B 
local i,j,m:

for i from 1 to N1 do 
 for j from 1 to N2 do 
  C[i,j]:=add(A[i,m]*B[m,j],m=1..M):
od:od:

end proc:

RefinedMatrixInverse:=proc(N,A,AI,Acc)
#A is the NxN matrix  whose inverse you want to take 
#AI is the returned inverse accurate up to 10^(-Acc)
local NMax,w0,wI,wc,we,w1,i,j,n:
global Digits:

NMax:=1000:
n:=0:


w0:=MatrixInverse(A):

wI:=IdentityMatrix(N):
wc:=Matrix(N,N):
we:=Matrix(N,N):
w1:=Matrix(N,N):

do

n:=n+1:

MatMult(N,N,N,w0,A,w1):

#print(w1);

for i from 1 to N do
 for j from 1 to N do
  we[i,j]:=(wI[i,j]-w1[i,j]);
od:od:

#print(we);

MatMult(N,N,N,we,w0,wc):

#print(wc);

for i from 1 to N do
 for j from 1 to N do
  w0[i,j]:=(w0[i,j]+wc[i,j]);
od:od:

until(max(abs(we))<10^(-Acc) or n=NMax):

if n=NMax then print("MAX CYCLES REACHED: FINAL ERROR WAS",we): fi:

for i from 1 to N do 
 for j from 1 to N do 
  AI[i,j]:=w0[i,j]:
od:od:



end proc:

