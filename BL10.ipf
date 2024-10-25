#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include "List_util"
#include "Image_Tool4053"

Menu "BL10" 

	"Convert pixel to angle", p2a()
	"B field correction", Bcor()
End
////////////////////////////////////////////////////////////
Proc p2a(img) : GraphMarquee
string img
prompt img, "Target wave", popup, WaveList("!*ct", ";", "DIMS:3")
variable pa=0.04631
setscale/p x (dimoffset($img,0)-500)*dimdelta($img,0)*0.04631,dimdelta($img,0)*0.04631, $img
end 
///////////////////////////////////////////////////////////////////////
Proc Bcor(img,eng,Arng,outp) : Graph Marqee
string img="f00023"
variable eng=-4.5
string arng=StrVarOrDefault("root","-74.5,-44.2")
string outp="s"
prompt img, "Target wave", popup, WaveList("!*ct", ";", "DIMS:3")
prompt eng, "Cut energy"
prompt arng, "Data range in angle, (start,end)"
prompt outp, "Oupput wave name or extension"

variable x0=ValFromList(arng,0,",")
variable x1=ValFromList(arng,1,",")

duplicate/o $img, $(img+outp)

findmax(eng,x0,x1,$img)
  
string cmd="VolShift(root:"+img+",root:fit_ref"+",\""+"/XZ/D=root:"+img+outp+"/E=-1"+"\")"
print cmd
execute cmd
setscale/p x dimoffset($img,0),dimdelta($img,0), $(img+outp)
end 

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
function findmax(w0,x0,x1,x2)
variable w0
variable x0,x1
wave x2

make/o/n=(dimsize(x2,0),dimsize(x2,2)) D2
setscale/p x dimoffset(x2,0),dimdelta(x2,0),D2
setscale/p y dimoffset(x2,2),dimdelta(x2,2),D2
d2[][]=x2[p](w0)[q]
duplicate/o/r=(x0,x1) d2, d2r
make/o/n=(dimsize(d2r,0)) D1y
setscale/p x dimoffset(d2r,0),dimdelta(d2r,0),D1y
d1y=0
make/o/n=(dimsize(d2r,1)) ref
setscale/p x dimoffset(d2r,1),dimdelta(d2r,1), ref
ref=0

variable dum, dumx
variable x,y

if(x0>0)
for(y=0;y<dimsize(d2r,1);y+=1)
	d1y=d2r[p][y]
	//dum=d1y[dimsize(d2r,0)-1]
	for(x=0;x<5;x+=1)
		dum+=d1y[dimsize(d2r,0)-1-x]
	endfor
	dum/=5
	d1y-=dum
	dum=wavemax(d1y)
	d1y/=dum
	for(x=dimsize(d2r,0)-2;x>1;x-=1)
		if(d1y[x]>0.5)
			dumx=x
			break
		endif
	endfor
	//print dum 
	ref[y]=dimoffset(d1y,0)+dimdelta(d1y,0)*dumx
endfor
else
for(y=0;y<dimsize(d2r,1);y+=1)
	d1y=d2r[p][y]
	//dum=d1y[dimsize(d2r,0)-1]
	for(x=0;x<5;x+=1)
		dum+=d1y[x]
	endfor
	dum/=5
	d1y-=dum
	dum=wavemax(d1y)
	d1y/=dum
	for(x=0;x<dimsize(d1y,0);x+=1)
		if(d1y[x]>0.5)
			dumx=x
			break
		endif
	endfor
	//print dum 
	ref[y]=dimoffset(d1y,0)+dimdelta(d1y,0)*dumx
endfor
endif
//display/k=1 ref
CurveFit/Q/M=2/W=0 poly 6, ref/D
killwaves d2,d2r,d1y
end 