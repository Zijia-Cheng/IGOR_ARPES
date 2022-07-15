#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.





Function ZJC_DrawBz(vector1,vector2)

//This is for creating a surface BZ ploting matrix
//vector1 and vector2 are reciprocal lattice vectors, written in an orthogonal basis
//return a wave named BZ. Then you can plot the BZ using display/appendtograph Bz[][1] vs BZ[][0]
//By CZJ 



wave vector1
wave vector2

Make /N=(8,2) /O nearest_Bz

variable temp = 0
variable i,j

for(i=-1;i<2;i++)
	for(j=-1;j<2;j++)
		if(i==0 && j==0)
		else	
			nearest_BZ[temp][] = i*vector1[q] + j*vector2[q]
			temp++
		endif
	
	endfor
endfor



Make /O /N =2 start_point
Make /O /N =(8,6,3) vertexes
variable temp2
make /O /N = 3 temp3


for(i=0;i<8;i++)
	temp=0
	for(j=0;j<8;j++)
		if((nearest_BZ[i][0] == - nearest_BZ[j][0])&& (nearest_BZ[i][1] == - nearest_BZ[j][1]))
			continue
		elseif(i==j)
			continue
		else
			vertexes[i][temp][0] = (nearest_BZ[j][1]*(nearest_BZ[i][0]^2+nearest_BZ[i][1]^2) - nearest_BZ[i][1]*(nearest_BZ[j][0]^2+nearest_BZ[j][1]^2))/(2*(nearest_BZ[i][0]*nearest_BZ[j][1]-nearest_BZ[j][0]*nearest_BZ[i][1]))
			if(nearest_BZ[i][1]==0)
				vertexes[i][temp][1] = -(nearest_BZ[j][0]/nearest_BZ[j][1])*vertexes[i][temp][0]+(nearest_BZ[j][0]^2/(2*nearest_BZ[j][1])+nearest_BZ[j][1]/2)
			else
				vertexes[i][temp][1] = -(nearest_BZ[i][0]/nearest_BZ[i][1])*vertexes[i][temp][0]+(nearest_BZ[i][0]^2/(2*nearest_BZ[i][1])+nearest_BZ[i][1]/2)
			endif
			
			vertexes[i][temp][2] = sqrt(vertexes[i][temp][0]^2+vertexes[i][temp][1]^2)
			//sort based on the length
			temp2 = temp
			temp3 = 0
			if(temp != 0)
				do
					if(vertexes[i][temp2][2] < vertexes[i][temp2-1][2])
						temp3[] = vertexes[i][temp2-1][p]
						vertexes[i][temp2-1][] = vertexes[i][temp2][r]
						vertexes[i][temp2][] = temp3[r]
						temp2 = temp2 - 1
					else
						break
					endif
				while(temp2>0)
			endif
			temp = temp+1
		endif
	endfor
	
	
	if(i == 0)
		start_point[] = vertexes[0][0][p]
	else
		if(vertexes[i][0][2] < norm(start_point))
			start_point[] = vertexes[i][0][p]
		endif	
	endif
endfor

Make /O /N = (200,2) BZ //This number can becomes larger if there is not enough space
Duplicate /O start_point, vector_start_point
Duplicate /O start_point, vector_end_point

temp = 0
variable temp4=0
variable last_side = -1
variable new_side
do
	
	BZ[temp][] = vector_start_point[q]
	temp2 = 0
	for(i=0;i<8;i++)
		if(i == last_side)
			continue
		else
			for(j =0;j<6;j++)
				if((abs(vertexes[i][j][0]-vector_start_point[0])<0.000001) && (abs(vertexes[i][j][1]-vector_start_point[1])<0.000001))
					temp4=0
					if(temp2 == 0)//initialization
						do
							if((abs(vertexes[i][temp4][0]-vector_start_point[0])>0.000001)||(abs(vertexes[i][temp4][1]-vector_start_point[1])>0.000001))
								vector_end_point[] = vertexes[i][temp4][p]
								new_side = i
								break
							else
								temp4 = temp4+1
							endif
						while(temp4<6)
					else
						do
							if(((abs(vertexes[i][temp4][0]-vector_start_point[0])>0.000001)||(abs(vertexes[i][temp4][1]-vector_start_point[1])>0.000001))&&(vertexes[i][temp4][2]<norm(vector_end_point)))
								vector_end_point[] = vertexes[i][temp4][p]
								new_side = i
								break
							else
								temp4 = temp4+1
							endif
						while(temp4<6)
						
					endif
					temp2 = temp2+1
					break
				endif
			endfor
		endif
	endfor
	last_side = new_side
	BZ[temp+1][] = vector_end_point[q]
	BZ[temp+2][] = Nan
	vector_start_point[] = vector_end_point[p]
	temp = temp+3
	
while((abs(vector_start_point[0]-start_point[0])>0.000001)|| (abs(vector_start_point[1]-start_point[1])>0.000001))




temp4=1
for(i=-1;i<2;i++)
	for(j=-1;j<2;j++)
		if(i==0 && j==0)
		else	
			BZ[temp4*temp,(temp4+1)*temp-1][] = BZ[p-temp4*temp][q]+i*vector1[q] + j*vector2[q]
			temp4++
		endif
	
	endfor
endfor

for(i =(temp4)*temp;i<200;i++)
	BZ[i][] = Nan
endfor


killwaves nearest_Bz,vertexes,start_point, vector_start_point,vector_end_point


End






Function DrawrecBZ(vector1,vector2)
//This is for making rectangle BZs waves for plotting purpose 
//Vector1 and vector2 should be 1d waves with length = 2. Those are reciprocal lattice vector
//Generate BZ and expandBZ


wave vector1
wave vector2

Make /O /N=(12,2) BZ
Make /O /N=(5,2) vertex

vertex[0][] = 0
vertex[1][] = vector1[q]
vertex[2][] = vector1[q]+vector2[q]
vertex[3][] = vector2[q]
vertex[4][] = 0


variable i 
for(i=0;i<4;i++)
	BZ[3*i][] = vertex[i][q]
	BZ[3*i+1][] = vertex[i+1][q]
	BZ[3*i+2][] = Nan
endfor

BZ[][0] = BZ[p][0] - (vector1[0]+vector2[0])/2
BZ[][1] = BZ[p][1] - (vector1[1]+vector2[1])/2

Killwaves vertex


Make /N = (108,2) /O expandBZ

variable j = 0
variable s = 0

Duplicate /O BZ, temp

for(i = -1;i<2;i++)
	for(j= -1;j<2;j++)
	
	temp[][] = temp[p][q]+i*vector1[q]+j*vector2[q]
	
	expandBZ[s*12,s*12+11][] = temp[p-s*12][q]
	s++
	temp = BZ[p][q]
	

	endfor
endfor


killwaves temp
end

