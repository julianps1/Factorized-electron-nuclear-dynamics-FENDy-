
MaxInd:=proc(N,V)
# returns index of largest element of V, out of N
local x,i,j:
x:=-10^Digits: j:=0: for i from 1 to N do if (x<V[i]) then x:=evalf(V[i]): j:=eval(i): fi: od:
j:
end proc:

PrimNorm:=proc(a,n)
#normalization of a primitive Gaussian
simplify(a^((1+2*n)/4)*(2/Pi)^(1/4)*2^n/sqrt(doublefactorial(2*n-1))):
end proc:
I2init:=proc(lmax)
# Generating equations of integrals over MO primitives

# Using convention Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,Dx,Dy,Dz
#A=1al,B=a1r,C=a2l,D=a2r
global I2Arr,a,b,c,d,Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,Dx,Dy,Dz,LMX,IntFile2:
local i1,s,dS,dT,Px,Py,Pz,Qx,Qy,Qz,w0,wds,wdt,w1,
ax,ay,az,bx,gy,bz,cx,cy,cz,dx,dy,dz,kmx:
LMX:=lmax: 
#IntFile:="/home/rassolov/Maple/BasisSets/Int2e.m":
#IntFile2:="Int2e.m":
I2Arr:=Array(0..lmax,0..lmax,0..lmax,0..lmax,0..lmax,0..lmax,0..lmax,0..lmax,0..lmax,0..lmax,0..lmax,0..lmax):
dS:=Vector(12): dT:=Vector(12): #dS derivative of the overlap, dT is derivative of T that is in Boys
dS[1]:=-2*a*b/(a+b)*(Ax-Bx): dS[2]:=-2*a*b/(a+b)*(Ay-By): dS[3]:=-2*a*b/(a+b)*(Az-Bz):
dS[4]:=-dS[1]: dS[5]:=-dS[2]: dS[6]:=-dS[3]:
dS[7]:=-2*c*d/(c+d)*(Cx-Dx): dS[8]:=-2*c*d/(c+d)*(Cy-Dy): dS[9]:=-2*c*d/(c+d)*(Cz-Dz):
dS[10]:=-dS[7]: dS[11]:=-dS[8]: dS[12]:=-dS[9]:
w1:=2*(a+b)*(c+d)/(a+b+c+d): 
Px:=(a*Ax+b*Bx)/(a+b): Py:=(a*Ay+b*By)/(a+b): Pz:=(a*Az+b*Bz)/(a+b):
Qx:=(c*Cx+d*Dx)/(c+d): Qy:=(c*Cy+d*Dy)/(c+d): Qz:=(c*Cz+d*Dz)/(c+d):
dT[1]:=w1*(Px-Qx)*a/(a+b): dT[2]:=w1*(Py-Qy)*a/(a+b): dT[3]:=w1*(Pz-Qz)*a/(a+b):
dT[4]:=w1*(Px-Qx)*b/(a+b): dT[5]:=w1*(Py-Qy)*b/(a+b): dT[6]:=w1*(Pz-Qz)*b/(a+b):
dT[7]:=w1*(Qx-Px)*c/(c+d): dT[8]:=w1*(Qy-Py)*c/(c+d): dT[9]:=w1*(Qz-Pz)*c/(c+d):
dT[10]:=w1*(Qx-Px)*d/(c+d): dT[11]:=w1*(Qy-Py)*d/(c+d): dT[12]:=w1*(Qz-Pz)*d/(c+d):

I2Arr[0,0,0,0,0,0,0,0,0,0,0,0]:=Array(0..1):
I2Arr[0,0,0,0,0,0,0,0,0,0,0,0][0]:=1:
for ax from 0 to lmax do
 for ay from 0 to lmax-ax do
  for az from 0 to lmax-ax-ay do
for bx from 0 to lmax do
 for gy from 0 to lmax-bx do
  for bz from 0 to lmax-bx-gy do
for cx from 0 to lmax do
 for cy from 0 to lmax-cx do
  for cz from 0 to lmax-cx-cy do
for dx from 0 to lmax do
 for dy from 0 to lmax-dx do
  for dz from 0 to lmax-dx-dy do
  kmx:=ax+ay+az+bx+gy+bz+cx+cy+cz+dx+dy+dz:
  if(dx=lmax and cx=lmax and bx+gy+bz=lmax and ax+ay+az=lmax) then print("VR",kmx,ax,bx,cx,dx," ",ay,gy,cy,dy," ",az,bz,cz,dz); fi:

  if(kmx>0) then
   I2Arr[ax,ay,az,bx,gy,bz,cx,cy,cz,dx,dy,dz]:=Array(0..kmx):

   if(ax>0) then
    for i1 from 0 to kmx do
     if(i1<kmx) then       wds:=I2Arr[ax-1,ay,az,bx,gy,bz,cx,cy,cz,dx,dy,dz][i1]: else wds:=0: fi:
     if(i1>0) then         wdt:=I2Arr[ax-1,ay,az,bx,gy,bz,cx,cy,cz,dx,dy,dz][i1-1]: else wdt:=0: fi:
     if(ax>1 and i1<kmx-1) then w0:=I2Arr[ax-2,ay,az,bx,gy,bz,cx,cy,cz,dx,dy,dz][i1]: else w0:=0: fi:
     I2Arr[ax,ay,az,bx,gy,bz,cx,cy,cz,dx,dy,dz][i1]:=(dS[1]*wds+diff(wds,Ax))/2/a+(ax-1)/2/a*w0-dT[1]/2/a*wdt:
    od:
    
   elif(ay>0) then
    for i1 from 0 to kmx do
     if(i1<kmx) then           wds:=I2Arr[ax,ay-1,az,bx,gy,bz,cx,cy,cz,dx,dy,dz][i1]: else wds:=0: fi:
     if(i1>0) then         wdt:=I2Arr[ax,ay-1,az,bx,gy,bz,cx,cy,cz,dx,dy,dz][i1-1]: else wdt:=0: fi:
     if(ay>1 and i1<kmx-1) then w0:=I2Arr[ax,ay-2,az,bx,gy,bz,cx,cy,cz,dx,dy,dz][i1]: else w0:=0: fi:
     I2Arr[ax,ay,az,bx,gy,bz,cx,cy,cz,dx,dy,dz][i1]:=(dS[2]*wds+diff(wds,Ay))/2/a+(ay-1)/2/a*w0-dT[2]/2/a*wdt:
    od:
    
   elif(az>0) then
     for i1 from 0 to kmx do
     if(i1<kmx) then           wds:=I2Arr[ax,ay,az-1,bx,gy,bz,cx,cy,cz,dx,dy,dz][i1]: else wds:=0: fi:
     if(i1>0) then         wdt:=I2Arr[ax,ay,az-1,bx,gy,bz,cx,cy,cz,dx,dy,dz][i1-1]: else wdt:=0: fi:
     if(az>1 and i1<kmx-1) then w0:=I2Arr[ax,ay,az-2,bx,gy,bz,cx,cy,cz,dx,dy,dz][i1]: else w0:=0: fi:
     I2Arr[ax,ay,az,bx,gy,bz,cx,cy,cz,dx,dy,dz][i1]:=(dS[3]*wds+diff(wds,Az))/2/a+(az-1)/2/a*w0-dT[3]/2/a*wdt:
    od:
    
   elif(bx>0) then
     for i1 from 0 to kmx do
     if(i1<kmx) then           wds:=I2Arr[ax,ay,az,bx-1,gy,bz,cx,cy,cz,dx,dy,dz][i1]: else wds:=0: fi:
     if(i1>0) then         wdt:=I2Arr[ax,ay,az,bx-1,gy,bz,cx,cy,cz,dx,dy,dz][i1-1]: else wdt:=0: fi:
     if(bx>1 and i1<kmx-1) then w0:=I2Arr[ax,ay,az,bx-2,gy,bz,cx,cy,cz,dx,dy,dz][i1]: else w0:=0: fi:
     I2Arr[ax,ay,az,bx,gy,bz,cx,cy,cz,dx,dy,dz][i1]:=(dS[4]*wds+diff(wds,Bx))/2/b+(bx-1)/2/b*w0-dT[4]/2/b*wdt:
    od:
    
   elif(gy>0) then
     for i1 from 0 to kmx do
     if(i1<kmx) then           wds:=I2Arr[ax,ay,az,bx,gy-1,bz,cx,cy,cz,dx,dy,dz][i1]: else wds:=0: fi:
     if(i1>0) then         wdt:=I2Arr[ax,ay,az,bx,gy-1,bz,cx,cy,cz,dx,dy,dz][i1-1]: else wdt:=0: fi:
     if(gy>1 and i1<kmx-1) then w0:=I2Arr[ax,ay,az,bx,gy-2,bz,cx,cy,cz,dx,dy,dz][i1]: else w0:=0: fi:
     I2Arr[ax,ay,az,bx,gy,bz,cx,cy,cz,dx,dy,dz][i1]:=(dS[5]*wds+diff(wds,By))/2/b+(gy-1)/2/b*w0-dT[5]/2/b*wdt:
    od:
    
   elif(bz>0) then
     for i1 from 0 to kmx do
     if(i1<kmx) then           wds:=I2Arr[ax,ay,az,bx,gy,bz-1,cx,cy,cz,dx,dy,dz][i1]: else wds:=0: fi:
     if(i1>0) then         wdt:=I2Arr[ax,ay,az,bx,gy,bz-1,cx,cy,cz,dx,dy,dz][i1-1]: else wdt:=0: fi:
     if(bz>1 and i1<kmx-1) then w0:=I2Arr[ax,ay,az,bx,gy,bz-2,cx,cy,cz,dx,dy,dz][i1]: else w0:=0: fi:
     I2Arr[ax,ay,az,bx,gy,bz,cx,cy,cz,dx,dy,dz][i1]:=(dS[6]*wds+diff(wds,Bz))/2/b+(bz-1)/2/b*w0-dT[6]/2/b*wdt:
    od:
    
   elif(cx>0) then
     for i1 from 0 to kmx do
     if(i1<kmx) then           wds:=I2Arr[ax,ay,az,bx,gy,bz,cx-1,cy,cz,dx,dy,dz][i1]: else wds:=0: fi:
     if(i1>0) then         wdt:=I2Arr[ax,ay,az,bx,gy,bz,cx-1,cy,cz,dx,dy,dz][i1-1]: else wdt:=0: fi:
     if(cx>1 and i1<kmx-1) then w0:=I2Arr[ax,ay,az,bx,gy,bz,cx-2,cy,cz,dx,dy,dz][i1]: else w0:=0: fi:
     I2Arr[ax,ay,az,bx,gy,bz,cx,cy,cz,dx,dy,dz][i1]:=(dS[7]*wds+diff(wds,Cx))/2/c+(cx-1)/2/c*w0-dT[7]/2/c*wdt:
    od:
       
    elif(cy>0) then
     for i1 from 0 to kmx do
     if(i1<kmx) then           wds:=I2Arr[ax,ay,az,bx,gy,bz,cx,cy-1,cz,dx,dy,dz][i1]: else wds:=0: fi:
     if(i1>0) then         wdt:=I2Arr[ax,ay,az,bx,gy,bz,cx,cy-1,cz,dx,dy,dz][i1-1]: else wdt:=0: fi:
     if(cy>1 and i1<kmx-1) then w0:=I2Arr[ax,ay,az,bx,gy,bz,cx,cy-2,cz,dx,dy,dz][i1]: else w0:=0: fi:
     I2Arr[ax,ay,az,bx,gy,bz,cx,cy,cz,dx,dy,dz][i1]:=(dS[8]*wds+diff(wds,Cy))/2/c+(cy-1)/2/c*w0-dT[8]/2/c*wdt:
    od:
       
   elif(cz>0) then
     for i1 from 0 to kmx do
#print("VRa",i1,kmx,dz);
     if(i1<kmx) then           wds:=I2Arr[ax,ay,az,bx,gy,bz,cx,cy,cz-1,dx,dy,dz][i1]: else wds:=0: fi:
     if(i1>0) then         wdt:=I2Arr[ax,ay,az,bx,gy,bz,cx,cy,cz-1,dx,dy,dz][i1-1]: else wdt:=0: fi:
     if(cz>1 and i1<kmx-1) then w0:=I2Arr[ax,ay,az,bx,gy,bz,cx,cy,cz-2,dx,dy,dz][i1]: else w0:=0: fi:
     I2Arr[ax,ay,az,bx,gy,bz,cx,cy,cz,dx,dy,dz][i1]:=(dS[9]*wds+diff(wds,Cz))/2/c+(cz-1)/2/c*w0-dT[9]/2/c*wdt:
    od:
       
   elif(dx>0) then
     for i1 from 0 to kmx do
     if(i1<kmx) then           wds:=I2Arr[ax,ay,az,bx,gy,bz,cx,cy,cz,dx-1,dy,dz][i1]: else wds:=0: fi:
     if(i1>0) then         wdt:=I2Arr[ax,ay,az,bx,gy,bz,cx,cy,cz,dx-1,dy,dz][i1-1]: else wdt:=0: fi:
     if(dx>1 and i1<kmx-1) then w0:=I2Arr[ax,ay,az,bx,gy,bz,cx,cy,cz,dx-2,dy,dz][i1]: else w0:=0: fi:
     I2Arr[ax,ay,az,bx,gy,bz,cx,cy,cz,dx,dy,dz][i1]:=(dS[10]*wds+diff(wds,Dx))/2/d+(dx-1)/2/d*w0-dT[10]/2/d*wdt:
    od:
       
   elif(dy>0) then
     for i1 from 0 to kmx do
     if(i1<kmx) then           wds:=I2Arr[ax,ay,az,bx,gy,bz,cx,cy,cz,dx,dy-1,dz][i1]: else wds:=0: fi:
     if(i1>0) then         wdt:=I2Arr[ax,ay,az,bx,gy,bz,cx,cy,cz,dx,dy-1,dz][i1-1]: else wdt:=0: fi:
     if(dy>1 and i1<kmx-1) then w0:=I2Arr[ax,ay,az,bx,gy,bz,cx,cy,cz,dx,dy-2,dz][i1]: else w0:=0: fi:
     I2Arr[ax,ay,az,bx,gy,bz,cx,cy,cz,dx,dy,dz][i1]:=(dS[11]*wds+diff(wds,Dy))/2/d+(dy-1)/2/d*w0-dT[11]/2/d*wdt:
    od:
       
   elif(dz>0) then
     for i1 from 0 to kmx do
     if(i1<kmx) then           wds:=I2Arr[ax,ay,az,bx,gy,bz,cx,cy,cz,dx,dy,dz-1][i1]: else wds:=0: fi:
     if(i1>0) then         wdt:=I2Arr[ax,ay,az,bx,gy,bz,cx,cy,cz,dx,dy,dz-1][i1-1]: else wdt:=0: fi:
     if(dz>1 and i1<kmx-1) then w0:=I2Arr[ax,ay,az,bx,gy,bz,cx,cy,cz,dx,dy,dz-2][i1]: else w0:=0: fi:
     I2Arr[ax,ay,az,bx,gy,bz,cx,cy,cz,dx,dy,dz][i1]:=(dS[12]*wds+diff(wds,Dz))/2/d+(dz-1)/2/d*w0-dT[12]/2/d*wdt:
    od:
       
   else print("Shouldn't be here");
   fi:

  else I2Arr[ax,ay,az,bx,gy,bz,cx,cy,cz,dx,dy,dz][0]:=1:
 fi:
od: od: od:  od: od: od:  od: od: od:  od: od: od:
save LMX,I2Arr,IntFile2;
end proc:

I1PotInit:=proc()
# creating 1e nuclear potetnial integrals out of 2e integrals over MO primitives
# Using convention Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,Dx,Dy,Dz
# A=1al,B=a1r,C=a2l,D=a2r
global I2Arr,I1Arr,a,b,c,d,Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,Dx,Dy,Dz,Rx,Ry,Rz,LMX,IntFile2,IntFile1:
local i1,w1,ax,ay,az,bx,gy,bz,kmx,lmax:
#IntFile2:="/home/rassolov/Maple/BasisSets/Int2e.m":
#IntFile2:="Int2e.m":
if(FileTools[Exists](IntFile2)) then
 read IntFile2:
 print("Initializing 1e potential integrals with Lmax=",LMX);
else
 print("2e integral file",IntFile2," does not exist.  Run I2Init procedure first");
 return:
fi:
lmax:=LMX:
I1Arr:=Array(0..lmax,0..lmax,0..lmax,0..lmax,0..lmax,0..lmax):
I1Arr[0,0,0,0,0,0]:=Array(0..1):
I1Arr[0,0,0,0,0,0][0]:=1:
for ax from 0 to lmax do
 for ay from 0 to lmax-ax do
  for az from 0 to lmax-ax-ay do
for bx from 0 to lmax do
 for gy from 0 to lmax-bx do
  for bz from 0 to lmax-bx-gy do
  kmx:=ax+ay+az+bx+gy+bz:
  if(bx+gy+bz=lmax and ax+ay+az=lmax) then print("VR",kmx,ax,bx," ",ay,gy," ",az,bz); fi:
  if(kmx>0) then
   I1Arr[ax,ay,az,bx,gy,bz]:=Array(0..kmx):
   for i1 from 0 to kmx do
    w1:=simplify(subs([Dx=Rx,Cx=Rx,Dy=Ry,Cy=Ry,Dz=Rz,Cz=Rz,c=d],I2Arr[ax,ay,az,bx,gy,bz,0,0,0,0,0,0][i1])):
    I1Arr[ax,ay,az,bx,gy,bz][i1]:=limit(w1,d=infinity):
      #if(ax=1) then print("VRa",bx,gy,bz,i1,w1,I1Arr[ax,ay,az,bx,gy,bz][i1]); fi:
   od:
  fi:
 od: od: od: od: od: od:
 #IntFile:="/home/rassolov/Maple/BasisSets/Int1e.m":
 #IntFile:="Int1e.m":
 save LMX,I1Arr,IntFile1;
 end proc:


I2eval:=proc(APr,BPr,CPr,DPr)
# Evaluating integrals over (un)normalized MO primitives
#XPr is a vector of parameters: APr[1]:=a (i.e. alpha), APr[2]:=[Ax,Ay,Az], APr[3]:=[ax,ay,az]
# Using convention Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,Dx,Dy,Dz
#A=1al,B=a1r,C=a2l,D=a2r
global I2Arr,a,b,c,d,Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,Dx,Dy,Dz:
local tiny,i1,Px,Py,Pz,Qx,Qy,Qz,w0,I00,FB,w1,T,
ax,ay,az,bx,gy,bz,cx,cy,cz,dx,dy,dz,kmx,AB,CD,PQ,subst,ast:
tiny:=evalf(sqrt(10^(-Digits))):
subst:=a=APr[1],b=BPr[1],c=CPr[1],d=DPr[1],Ax=APr[2][1],Ay=APr[2][2],Az=APr[2][3],Bx=BPr[2][1],By=BPr[2][2],Bz=BPr[2][3],Cx=CPr[2][1],Cy=CPr[2][2],Cz=CPr[2][3],Dx=DPr[2][1],Dy=DPr[2][2],Dz=DPr[2][3]:
ast:=APr[3][1],APr[3][2],APr[3][3],BPr[3][1],BPr[3][2],BPr[3][3],CPr[3][1],CPr[3][2],CPr[3][3],DPr[3][1],DPr[3][2],DPr[3][3]:
AB:=add((APr[2][i1]-BPr[2][i1])^2,i1=1..3):
CD:=add((CPr[2][i1]-DPr[2][i1])^2,i1=1..3):
w0:=APr[1]+BPr[1]+CPr[1]+DPr[1]:
PQ:=add(((APr[1]*APr[2][i1]+BPr[1]*BPr[2][i1])/(APr[1]+BPr[1])-(CPr[1]*CPr[2][i1]+DPr[1]*DPr[2][i1])/(CPr[1]+DPr[1]))^2,i1=1..3):
I00:=evalf(exp(-APr[1]*BPr[1]/(APr[1]+BPr[1])*AB-CPr[1]*DPr[1]/(CPr[1]+DPr[1])*CD)*2*Pi*w0/(APr[1]+BPr[1])/(CPr[1]+DPr[1])*(Pi/w0)^(3/2)):
T:=evalf((APr[1]+BPr[1])*(CPr[1]+DPr[1])/w0*PQ):
kmx:=add(APr[3][i1]+BPr[3][i1]+CPr[3][i1]+DPr[3][i1],i1=1..3):
FB:=Array(0..kmx):
if(T<tiny) then
 for i1 from 0 to kmx do FB[i1]:=evalf(1/(2*i1+1)-T/(2*i1+3)+T^2/2/(2*i1+5)): od
else
 FB[0]:=evalf(sqrt(Pi/4/T)*erf(sqrt(T))): w0:=exp(-T):
 for i1 from 1 to kmx do FB[i1]:=evalf(((2*i1-1)*FB[i1-1]-w0)/2/T): od:
fi:

#for i1 from 0 to kmx do print("DB",i1,evalf(I00*subs([subst],I2Arr[ast][i1]*FB[i1]))); od:

w0:=evalf(I00*add(subs([subst],I2Arr[ast][i1]*FB[i1]),i1=0..kmx)):
print("VR2e",w0);

w1:=product((4*APr[1])^APr[3][i2]/doublefactorial(2*APr[3][i2]-1)*sqrt(2*APr[1]/Pi),i2=1..3)*
product((4*BPr[1])^BPr[3][i2]/doublefactorial(2*BPr[3][i2]-1)*sqrt(2*BPr[1]/Pi),i2=1..3)*
product((4*CPr[1])^CPr[3][i2]/doublefactorial(2*CPr[3][i2]-1)*sqrt(2*CPr[1]/Pi),i2=1..3)*
product((4*DPr[1])^DPr[3][i2]/doublefactorial(2*DPr[3][i2]-1)*sqrt(2*DPr[1]/Pi),i2=1..3):
w1:=evalf(w0*sqrt(w1)):
end proc:

I1eval:=proc(APr,BPr,RR)
# Evaluating integrals over (un)normalized MO primitives
#XPr is a vector of parameters: APr[1]:=a (i.e. alpha), APr[2]:=[Ax,Ay,Az], APr[3]:=[ax,ay,az]
# Using convention Ax,Ay,Az,Bx,By,Bz
#A=1al,B=a1r
global a,b,Ax,Ay,Az,Bx,By,Bz,Rx,Ry,Rz:
local tiny,i1,w0,I00,FB,w1,ax,ay,az,bx,gy,bz,kmx,AB,subst,ast,PQ,T:
tiny:=evalf(sqrt(10^(-Digits))):
subst:=a=APr[1],b=BPr[1],Ax=APr[2][1],Ay=APr[2][2],Az=APr[2][3],Bx=BPr[2][1],By=BPr[2][2],Bz=BPr[2][3],Rx=RR[1],Ry=RR[2],Rz=RR[3]:
ast:=APr[3][1],APr[3][2],APr[3][3],BPr[3][1],BPr[3][2],BPr[3][3]:
AB:=add((APr[2][i1]-BPr[2][i1])^2,i1=1..3):
PQ:=add(((APr[1]*APr[2][i1]+BPr[1]*BPr[2][i1])/(APr[1]+BPr[1])-RR[i1])^2,i1=1..3):
w0:=APr[1]+BPr[1]:
#I00:=evalf(exp(-APr[1]*BPr[1]/(APr[1]+BPr[1])*AB)*2*Pi*w0/(APr[1]+BPr[1])*(Pi/w0)^(3/2)):
I00:=evalf(exp(-APr[1]*BPr[1]/(APr[1]+BPr[1])*AB)*2*Pi*w0/(APr[1]+BPr[1])/w0):
#T:=evalf((APr[1]+BPr[1])/w0*PQ):
T:=evalf((APr[1]+BPr[1])*PQ):
kmx:=add(APr[3][i1]+BPr[3][i1],i1=1..3):
FB:=Array(0..kmx):
if(T<tiny) then
 for i1 from 0 to kmx do FB[i1]:=evalf(1/(2*i1+1)-T/(2*i1+3)+T^2/2/(2*i1+5)): od
else
 FB[0]:=evalf(sqrt(Pi/4/T)*erf(sqrt(T))): w0:=exp(-T):
 for i1 from 1 to kmx do FB[i1]:=evalf(((2*i1-1)*FB[i1-1]-w0)/2/T): od:
fi:

#for i1 from 0 to kmx do print("DB",i1,evalf(I00*subs([subst],I1Arr[ast][i1]*FB[i1]))); od:

w0:=evalf(I00*add(subs([subst],I1Arr[ast][i1]*FB[i1]),i1=0..kmx)):
print("VR 1e",w0);
w1:=product((4*APr[1])^(APr[3][i2])/doublefactorial(2*APr[3][i2]-1)*sqrt(2*APr[1]/Pi),i2=1..3)*
product((4*BPr[1])^(BPr[3][i2])/doublefactorial(2*BPr[3][i2]-1)*sqrt(2*BPr[1]/Pi),i2=1..3):
w1:=evalf(w0*sqrt(w1)):
end proc:

I1evalB:=proc(A,APr,B,BPr,RR)
# Evaluating integrals over (un)normalized MO primitives
#In contrast o I1eval, alpha and beta values are separate parameters, for convenience of summing over the primitives
#XPr is a vector of parameters: A:=a (i.e. alpha), APr[1]:=[Ax,Ay,Az], APr[2]:=[ax,ay,az]
# Using convention Ax,Ay,Az,Bx,By,Bz
#A=1al,B=a1r
global a,b,Ax,Ay,Az,Bx,By,Bz,Rx,Ry,Rz:
local tiny,i1,i2,w0,I00,FB,w1,ax,ay,az,bx,gy,bz,kmx,AB,subst,ast,PQ,T:
tiny:=evalf(sqrt(10^(-Digits))):
#print("I1ev db1",A,B,APr,BPr,RR);
subst:=a=A,b=B,Ax=APr[1][1],Ay=APr[1][2],Az=APr[1][3],Bx=BPr[1][1],By=BPr[1][2],Bz=BPr[1][3],Rx=RR[1],Ry=RR[2],Rz=RR[3]:
ast:=APr[2][1],APr[2][2],APr[2][3],BPr[2][1],BPr[2][2],BPr[2][3]:
AB:=add((APr[1][i1]-BPr[1][i1])^2,i1=1..3):
PQ:=add(((A*APr[1][i1]+B*BPr[1][i1])/(A+B)-RR[i1])^2,i1=1..3):
w0:=A+B:
#I00:=evalf(exp(-A*B/(A+B)*AB)*2*Pi*w0/(A+B)*(Pi/w0)^(3/2)):
I00:=evalf(exp(-A*B/(A+B)*AB)*2*Pi*w0/(A+B)/w0):
#T:=evalf((A+B)/w0*PQ):
T:=evalf((A+B)*PQ):
kmx:=add(APr[2][i1]+BPr[2][i1],i1=1..3):
#print("I1ev db2",AB,I00,T,kmx);
FB:=Array(0..kmx):
if(T<tiny) then
 for i1 from 0 to kmx do FB[i1]:=evalf(1/(2*i1+1)-T/(2*i1+3)+T^2/2/(2*i1+5)): od
else
 FB[0]:=evalf(sqrt(Pi/4/T)*erf(sqrt(T))): w0:=exp(-T):
 for i1 from 1 to kmx do FB[i1]:=evalf(((2*i1-1)*FB[i1-1]-w0)/2/T): od:
fi:
#print("I1ev db3",FB[0],subst,i1,ast);
#print(I1Arr[ast]);

w0:=evalf(I00*add(subs([subst],I1Arr[ast][i1]*FB[i1]),i1=0..kmx)):
#print("VR 1e",w0,i2,APr[2],BPr[2],APr[2][1],BPr[2][3]);
w1:=(4*A)^(APr[2][1]+APr[2][2]+APr[2][3])/doublefactorial(2*APr[2][1]-1)/doublefactorial(2*APr[2][2]-1)/doublefactorial(2*APr[2][3]-1)
*(4*B)^(BPr[2][1]+BPr[2][2]+BPr[2][3])/doublefactorial(2*BPr[2][1]-1)/doublefactorial(2*BPr[2][2]-1)/doublefactorial(2*BPr[2][3]-1)*
(2*A*2*B)^(3/2)/Pi^3:
#w1:=product((4*A)^(APr[2][i2])/doublefactorial(2*APr[2][i2]-1)*sqrt(2*A/Pi),i2=1..3)*
#product((4*B)^(BPr[2][i2])/doublefactorial(2*BPr[2][i2]-1)*sqrt(2*B/Pi),i2=1..3):
w1:=evalf(w0*sqrt(w1)):
#print("I1ev db4",w1);

end proc:

I1evalRder:=proc(A,APr,B,BPr,RR,ider)
# Evaluating integrals over (un)normalized MO primitives for nuclear  over Coulomb operator in "ider" direction
#In contrast o I1eval, alpha and beta values are separate parameters, for convenience of summing over the primitives
#XPr is a vector of parameters: A:=a (i.e. alpha), APr[1]:=[Ax,Ay,Az], APr[2]:=[ax,ay,az]
# Using convention Ax,Ay,Az,Bx,By,Bz
#A=1al,B=a1r
global a,b,Ax,Ay,Az,Bx,By,Bz,Rx,Ry,Rz:
local tiny,i1,i2,w0,w1,w2,w3,w4,I00,FB,ax,ay,az,bx,gy,bz,kmx,AB,subst,ast,astAp,astAm,astBp,astBm,doA,doB,PQ,T:
tiny:=evalf(sqrt(10^(-Digits))):
#print("I1ev db1",A,B,APr,BPr,RR);
subst:=a=A,b=B,Ax=APr[1][1],Ay=APr[1][2],Az=APr[1][3],Bx=BPr[1][1],By=BPr[1][2],Bz=BPr[1][3],Rx=RR[1],Ry=RR[2],Rz=RR[3]:
ast:=APr[2][1],APr[2][2],APr[2][3],BPr[2][1],BPr[2][2],BPr[2][3]:

#new block to define other primitive Gaussians for nuclear derivatives
APr[2][ider]:=eval(APr[2][ider]+1): astAp:=APr[2][1],APr[2][2],APr[2][3],BPr[2][1],BPr[2][2],BPr[2][3]:
if(APr[2][ider]>1) then 
 doA:=1: APr[2][ider]:=eval(APr[2][ider]-2): astAm:=APr[2][1],APr[2][2],APr[2][3],BPr[2][1],BPr[2][2],BPr[2][3]: 
 APr[2][ider]:=eval(APr[2][ider]+1): 
else
 doA:=0: APr[2][ider]:=eval(APr[2][ider]-1):
fi:
BPr[2][ider]:=eval(BPr[2][ider]+1): astBp:=APr[2][1],APr[2][2],APr[2][3],BPr[2][1],BPr[2][2],BPr[2][3]:
if(BPr[2][ider]>1) then 
 doB:=1: BPr[2][ider]:=eval(BPr[2][ider]-2): astBm:=APr[2][1],APr[2][2],APr[2][3],BPr[2][1],BPr[2][2],BPr[2][3]: 
 BPr[2][ider]:=eval(BPr[2][ider]+1): 
else
 doB:=0: BPr[2][ider]:=eval(BPr[2][ider]-1):
fi:
#print("VR db",ast," ",astAp," ",astBp,astAm,doA,doB);

AB:=add((APr[1][i1]-BPr[1][i1])^2,i1=1..3):
PQ:=add(((A*APr[1][i1]+B*BPr[1][i1])/(A+B)-RR[i1])^2,i1=1..3):
w0:=A+B:
#I00:=evalf(exp(-A*B/(A+B)*AB)*2*Pi*w0/(A+B)*(Pi/w0)^(3/2)):
I00:=evalf(exp(-A*B/(A+B)*AB)*2*Pi*w0/(A+B)/w0):
#T:=evalf((A+B)/w0*PQ):
T:=evalf((A+B)*PQ):
kmx:=add(APr[2][i1]+BPr[2][i1],i1=1..3)+1:
#print("I1ev db2",AB,I00,T,kmx);
FB:=Array(0..kmx):
if(T<tiny) then
 for i1 from 0 to kmx do FB[i1]:=evalf(1/(2*i1+1)-T/(2*i1+3)+T^2/2/(2*i1+5)): od
else
 FB[0]:=evalf(sqrt(Pi/4/T)*erf(sqrt(T))): w0:=exp(-T):
 for i1 from 1 to kmx do FB[i1]:=evalf(((2*i1-1)*FB[i1-1]-w0)/2/T): od:
fi:
#print("I1ev db3",FB[0],subst,i1,ast);
#print(I1Arr[ast]);

#w0:=evalf(I00*add(subs([subst],I1Arr[ast][i1]*FB[i1]),i1=0..kmx-1)):
w1:=evalf(I00*add(subs([subst],I1Arr[astAp][i1]*FB[i1]),i1=0..kmx)):
if(doA=1) then w2:=evalf(I00*add(subs([subst],I1Arr[astAm][i1]*FB[i1]),i1=0..kmx-2)): else w2:=0: fi:
w3:=evalf(I00*add(subs([subst],I1Arr[astBp][i1]*FB[i1]),i1=0..kmx)):
if(doB=1) then w4:=evalf(I00*add(subs([subst],I1Arr[astBm][i1]*FB[i1]),i1=0..kmx-2)): else w4:=0: fi:

w0:=evalf(-(-APr[2][ider]*w2+2*A*w1-BPr[2][ider]*w4+2*B*w3)*PrimNorm(A,APr[2][1])*PrimNorm(A,APr[2][2])*PrimNorm(A,APr[2][3])
*PrimNorm(B,BPr[2][1])*PrimNorm(B,BPr[2][2])*PrimNorm(B,BPr[2][3])):

end proc:

