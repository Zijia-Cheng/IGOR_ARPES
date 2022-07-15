#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.






Function Photonenergydep(wavenamelist,photoenergylist,outputwavename)

//The input wavelist should be: the name of 2d cuts, which dim0 is kinetic energy and dim1 is angle(not momentum)
//photoenergylist is the list of photonenergies
//Usually this is used for SSRL 


wave /T wavenamelist
wave photoenergylist
String outputwavename

variable workfunction = 4.365

if(dimsize(wavenamelist,0)!=dimsize(photoenergylist,0))
	print("dimensions of input waves are not right")
	return(0)
endif

variable i=0
for(i=0;i<dimsize(wavenamelist,0);i++)
	cutregulation(wavenamelist[i],photoenergylist[i],workfunction)
endfor

WaveStats photoenergylist
string Max_wavename = wavenamelist[V_maxloc]

Make /O /N=(dimsize($Max_wavename,0),dimsize($Max_wavename,1),dimsize(wavenamelist,0)), $outputwavename
wave temp = $outputwavename

Setscale /P x, dimoffset($Max_wavename,0), dimdelta($Max_wavename,0),$outputwavename
Setscale /P y, dimoffset($Max_wavename,1), dimdelta($Max_wavename,1),$outputwavename
Setscale /I z, Wavemin(photoenergylist), Wavemax(photoenergylist),$outputwavename

Duplicate /O photoenergylist, temp2
sort temp2,temp2 


string wavenamest
variable j,k,l
for(j=0;j<dimsize(wavenamelist,0);j++)
	FindValue /V = (temp2[j]) photoenergylist
	wavenamest = wavenamelist[V_value]
	wave temp3 = $wavenamest
	
	for(k=0;k<dimsize($Max_wavename,0);k++)
		if(indexToscale($Max_wavename,k,0)<min(dimoffset($wavenamest,0),indextoScale($wavenamest,dimsize($wavenamest,0)-1,0))||indexToscale($Max_wavename,k,0)>Max(dimoffset($wavenamest,0),indextoScale($wavenamest,dimsize($wavenamest,0)-1,0)))
			temp[k][][j] = q-q
		else
			for(l=0;l<dimsize($Max_wavename,1);l++)
				if(indexToscale($Max_wavename,l,1)<min(dimoffset($wavenamest,1),indextoScale($wavenamest,dimsize($wavenamest,1)-1,1))||indexToscale($Max_wavename,l,1)>Max(dimoffset($wavenamest,1),indextoScale($wavenamest,dimsize($wavenamest,1)-1,1)))
					temp[k][l][j] = 0
				else
					temp[k][l][j] = temp3(indexToscale($Max_wavename,k,0))(indexToscale($Max_wavename,l,1))
				endif
			endfor
		endif
	endfor
endfor




End


Function cutregulation(wavenamest,photonenergy,workfunction)

string wavenamest
variable photonenergy
variable workfunction

Setscale /P x,dimoffset($wavenamest,0)-(photonenergy - workfunction), dimdelta($wavenamest,0), $wavenamest
Setscale /P y,0.5123*sqrt(photonenergy - workfunction)*Pi/180*(dimoffset($wavenamest,1)),0.5123*sqrt(photonenergy - workfunction)*Pi/180*(dimdelta($wavenamest,1)),$wavenamest

End




Function fixfermilevel_2D(ekcut,startenergy,energyrange,startangle,endangle)
//20220222 This is for aligning the Fermi level of a two dimensional cut
//The cut should have dim0: energy. dim1: angle/momentum
//startenergy should be the guessed Fermi level position for the cut at startangle(scale not index)
//startangle and endangle are angle/momentum range of interest(scale not index)


wave ekcut
variable startenergy,energyrange
variable startangle,endangle

Matrixop /O ekcut_=replaceNaNs(ekcut,0)
wave ekcut_
Setscale /P x,dimoffset(ekcut,0),dimdelta(ekcut,0),ekcut_
Setscale /P y,dimoffset(ekcut,1),dimdelta(ekcut,1),ekcut_


Make /O /N =(2*dimsize(ekcut,0),dimsize(ekcut,1)) $(nameofwave(ekcut)+"_aligned") = 0 //For storing the output result
wave temp1 = $(nameofwave(ekcut)+"_aligned")
Setscale /P x,-dimdelta(ekcut,0)*dimsize(ekcut,0),dimdelta(ekcut,0),temp1
Setscale /P y,dimoffset(ekcut,1),dimdelta(ekcut,1),temp1


make /O /N = (dimsize(ekcut,0)) energyspec = 0 //storing each EDC
Setscale /P x,dimOffset(ekcut,0),dimdelta(ekcut,0),energyspec


variable i,startangle_index = scaleToIndex(ekcut_,startangle,1)
variable endangle_index = scaleToIndex(ekcut_,endangle,1)
variable Fermi

make /O /N = (endangle_index - startangle_index+1) $(nameofwave(ekcut)+"_fl") = 0
wave fl = $(nameofwave(ekcut)+"_fl")
Setscale /I x,startangle,endangle,fl





for(i =scaleToIndex(ekcut_,startangle,1);i<=endangle_index;i++ )
	
	energyspec[] = ekcut_[p][i]
	Duplicate /O /R = (startenergy-energyrange, startenergy+energyrange) energyspec, $"test_12345678"
	
	wave temp = $"test_12345678"
	//method: use differentiate to get the largest change point in the integerated spectrum:
	//seems we need to narrow down the energy range.

	Smooth 3,temp//Should choose whether want to smooth or not
	Differentiate/METH=1 temp /D=energyspec_DIF
	wave energyspec_DIF
	energyspec_DIF[] = abs(energyspec_DIF[p])
	
	Findvalue /V = (Wavemax(energyspec_DIF)) energyspec_DIF
	fermi = Indextoscale(energyspec_DIF,V_value,0)
	
	startenergy = fermi //update the startenergy value for next cut
	fl[i-startangle_index] =fermi 

endfor

CurveFit/M=2/W=0 poly_XOffset 13, fl/D
wave fl_fit = $("fit_"+nameofwave(fl))

for(i =scaleToIndex(ekcut_,startangle,1);i<=endangle_index;i++ )
	
	temp1[][i] = ((x+fl_fit(indexToScale(ekcut_,i,1))-dimoffset(energyspec,0))*(x+fl_fit(indexToScale(ekcut_,i,1))-indextoscale(energyspec,dimsize(energyspec,0)-1,0))<=0)?ekcut_(x+fl_fit(indexToScale(ekcut_,i,1)))[i]:0
endfor

End










Function fixfermilevel(ekcut, startenergy, energyrange)
	//20191210 this is a part of the program used for align the fermi level for photon energy dependence, etc.
	//ekcut should be a 2d wave with dimension 0: angle/momentum; dimension 1: energy.
	//Method: find the largest derivative and that is the fermi energy
	//The fermi level should within startenergy+- energyrange
	
	wave ekcut
	variable startenergy, energyrange
	Make /O /N = (dimsize(ekcut,1)), energyspec = 0
	Variable i,j, fermi
	
	Matrixop /O ekcut_=replaceNaNs(ekcut,0)//Smooth can not handle waves with Nan
	wave ekcut_
	
	// make energy spectrum
	for(i=0;i < dimsize(ekcut,1);i++)
		for(j = 0; j  < dimsize(ekcut,0);j++)
			energyspec[i] =  ekcut_[j][i][0] + energyspec[i]
		endfor		
	endfor
	
	Setscale /P x , dimoffset(ekcut,1), dimdelta(ekcut,1), energyspec 
	Duplicate /O /R = (startenergy-energyrange, startenergy+energyrange) energyspec, $"test_12345678"
	
	wave temp = $"test_12345678"
	//method: use differentiate to get the largest change point in the integerated spectrum:
	//seems we need to narrow down the energy range.
	

	
	Smooth 3,temp//Should choose whether want to smooth or not
	Differentiate/METH=1 temp /D=energyspec_DIF
	wave energyspec_DIF
	
	energyspec_DIF[] = abs(energyspec_DIF[p])
	
	Findvalue /V = (Wavemax(energyspec_DIF)) energyspec_DIF
	fermi = Indextoscale(energyspec_DIF,V_value,0)
	
	//reset scaling of ekcutprin
	Setscale /P y, dimoffset(ekcut,1) - fermi, dimdelta(ekcut,1), ekcut
	killwaves temp
	
	return fermi
	
	
End


Function photonenergydep_ALS(inputwave, outputname,startenergy,energyrange)
		//20191210 by CZJ
		// This is for combiming different photon energy cuts. so we need to rescale kx at each photon energy and also align fermi levels.
		// The input should be a three dimensional wave with dimension : 0 -> angle, 1-> energy, 2-> light energy.
		// output wave's name is the string outputname. Its a three dimensional wave

		
		Wave inputwave
		String outputname
		Variable startenergy, energyrange
		
		
		// prepare different rescaled energy cuts
		
		Variable i
		Make /O /N = (dimsize(inputwave,2)) workfunction
		wave workfunction
		variable HighEn, LowEn, HighMo, LowMo,deltaMo, deltaEn
		
		for(i = 0; i< dimsize(inputwave,2);i++)
			
			String ekcut = "en"+num2str(Indextoscale(inputwave,i,2))
			Duplicate /O /RMD = [][][i] inputwave , $ekcut // remember here
			workfunction[i] = fixfermilevel($ekcut,startenergy,energyrange)
			Setscale /P x, (pi/180)*0.5123*sqrt(Indextoscale(inputwave,i,2) - 4.5)*dimoffset($ekcut,0), (pi/180)*0.5123*sqrt(Indextoscale(inputwave,i,2) - 4.5)*dimdelta($ekcut,0), $ekcut	
			
			if(i == 0)
				HighEn = Max(Dimoffset($ekcut,1), Dimoffset($ekcut,1)+(dimsize($ekcut,1)-1)*dimdelta($ekcut,1))
				LowEn = Min(Dimoffset($ekcut,1), Dimoffset($ekcut,1)+(dimsize($ekcut,1)-1)*dimdelta($ekcut,1))
				HighMo = Max(Dimoffset($ekcut,0), Dimoffset($ekcut,0)+(dimsize($ekcut,0)-1)*dimdelta($ekcut,0))
				LowMo = Min(Dimoffset($ekcut,0), Dimoffset($ekcut,0)+(dimsize($ekcut,0)-1)*dimdelta($ekcut,0))
				deltaMo = abs(dimdelta($ekcut,0))
				deltaEn = abs(dimdelta($ekcut,1))
				
			else
				HighEn = Max(HighEn,Max(Dimoffset($ekcut,1), Dimoffset($ekcut,1)+(dimsize($ekcut,1)-1)*dimdelta($ekcut,1)) )
				LowEn = Min(LowEn,Min(Dimoffset($ekcut,1), Dimoffset($ekcut,1)+(dimsize($ekcut,1)-1)*dimdelta($ekcut,1)) )
				HighMo = Max(HighMo,Max(Dimoffset($ekcut,0), Dimoffset($ekcut,0)+(dimsize($ekcut,0)-1)*dimdelta($ekcut,0)) )
				LowMo = Min(LowMo,Min(Dimoffset($ekcut,0), Dimoffset($ekcut,0)+(dimsize($ekcut,0)-1)*dimdelta($ekcut,0)) )
				deltaMo =Min(deltaMo,abs(dimdelta($ekcut,0)))
				deltaEn = Min(abs(dimdelta($ekcut,1)), deltaEn)
				
				
			endif
			
			
			
		endfor
		
		//Construct a 3D wave for output
		//First need to figure out the size along each direction
		
		Make /O /N = ((ceil((HighMo-LowMo)/deltaMo+1)),(ceil((HighEn-LowEn)/deltaEn+1)), dimsize(inputwave,2)) , $outputname
		wave output = $outputname
		
		//setupscale for outputwave
		
		setscale /I x, LowMo, HighMo, output
		setscale /I y, LowEn,  HighEn, output
		setscale  /P z, dimoffset(inputwave,2), dimdelta(inputwave,2), output
		
		//From output wave find the value in the original wave
		
		Variable j, k 
		
		for(i=0;i<dimsize(output,0);i++)
			for(j = 0; j<dimsize(output,1);j++)
				for(k=0; k<dimsize(output,2);k++)
					String temp = "en"+num2str(Indextoscale(output,k,2))
					wave temp1 = $temp
					if(indextoscale(output,i,0) > Max(Dimoffset($temp,0), Dimoffset($temp,0)+(dimsize($temp,0)-1)*dimdelta($temp,0)) || indextoscale(output,i,0) < Min(Dimoffset($temp,0), Dimoffset($temp,0)+(dimsize($temp,0)-1)*dimdelta($temp,0)) || indextoscale(output,j,1) > Max(Dimoffset($temp,1), Dimoffset($temp,1)+(dimsize($temp,1)-1)*dimdelta($temp,1)) || indextoscale(output,j,1) < Min(Dimoffset($temp,1), Dimoffset($temp,1)+(dimsize($temp,1)-1)*dimdelta($temp,1)))
						output[i][j][k] =0 // out of range
					
					else
					
					output[i][j][k] = temp1(indextoscale(output,i,0))(indextoscale(output,j,1))[0]
					
					endif
					
					
				endfor
			endfor
		endfor
		
		for(i = 0; i< dimsize(inputwave,2);i++)
			ekcut = "en"+num2str(Indextoscale(inputwave,i,2))
			KillWaves $ekcut
			
		
		endfor
		
		
		
		
End

Function photonenergydep_SSRL(inputwave, outputname,startenergy,energyrange,offsetwave)
		//20220528 by CZJ
		// This is for combiming different photon energy cuts. so we need to rescale kx at each photon energy and also align fermi levels.
		// The input should be a three dimensional wave with dimension : 0 -> angle, 1-> energy, 2-> light energy.
		// output wave's name is the string outputname. Its a three dimensional wave
		//Assuming work function to be 4.5 when transfering angle to momentum
		
		Wave inputwave
		String outputname
		Variable startenergy, energyrange
		wave offsetwave //This is the shift of each cut in the angle
		
		
		// prepare different rescaled energy cuts
		
		Variable i
		Make /O /N = (dimsize(inputwave,2)) workfunction
		wave workfunction
		variable HighEn, LowEn, HighMo, LowMo,deltaMo, deltaEn
		
		for(i = 0; i< dimsize(inputwave,2);i++)
			
			String ekcut = "en"+num2str(Indextoscale(inputwave,i,2))
			Duplicate /O /RMD = [][][i] inputwave , $ekcut // remember here
			workfunction[i] = fixfermilevel($ekcut,startenergy,energyrange)
			Setscale /P x, (pi/180)*0.5123*sqrt(Indextoscale(inputwave,i,2) - 4.5)*(dimoffset($ekcut,0)-offsetwave[i]), (pi/180)*0.5123*sqrt(Indextoscale(inputwave,i,2) - 4.5)*dimdelta($ekcut,0), $ekcut	
			
			if(i == 0)
				HighEn = Max(Dimoffset($ekcut,1), Dimoffset($ekcut,1)+(dimsize($ekcut,1)-1)*dimdelta($ekcut,1))
				LowEn = Min(Dimoffset($ekcut,1), Dimoffset($ekcut,1)+(dimsize($ekcut,1)-1)*dimdelta($ekcut,1))
				HighMo = Max(Dimoffset($ekcut,0), Dimoffset($ekcut,0)+(dimsize($ekcut,0)-1)*dimdelta($ekcut,0))
				LowMo = Min(Dimoffset($ekcut,0), Dimoffset($ekcut,0)+(dimsize($ekcut,0)-1)*dimdelta($ekcut,0))
				deltaMo = abs(dimdelta($ekcut,0))
				deltaEn = abs(dimdelta($ekcut,1))
				
			else
				HighEn = Max(HighEn,Max(Dimoffset($ekcut,1), Dimoffset($ekcut,1)+(dimsize($ekcut,1)-1)*dimdelta($ekcut,1)) )
				LowEn = Min(LowEn,Min(Dimoffset($ekcut,1), Dimoffset($ekcut,1)+(dimsize($ekcut,1)-1)*dimdelta($ekcut,1)) )
				HighMo = Max(HighMo,Max(Dimoffset($ekcut,0), Dimoffset($ekcut,0)+(dimsize($ekcut,0)-1)*dimdelta($ekcut,0)) )
				LowMo = Min(LowMo,Min(Dimoffset($ekcut,0), Dimoffset($ekcut,0)+(dimsize($ekcut,0)-1)*dimdelta($ekcut,0)) )
				deltaMo =Min(deltaMo,abs(dimdelta($ekcut,0)))
				deltaEn = Min(abs(dimdelta($ekcut,1)), deltaEn)
				
				
			endif
			
			
			
		endfor
		
		//Construct a 3D wave for output
		//First need to figure out the size along each direction
		
		Make /O /N = ((ceil((HighMo-LowMo)/deltaMo+1)),(ceil((HighEn-LowEn)/deltaEn+1)), dimsize(inputwave,2)) , $outputname
		wave output = $outputname
		
		//setupscale for outputwave
		
		setscale /I x, LowMo, HighMo, output
		setscale /I y, LowEn,  HighEn, output
		setscale  /P z, dimoffset(inputwave,2), dimdelta(inputwave,2), output
		
		//From output wave find the value in the original wave
		
		Variable j, k 
		
		for(i=0;i<dimsize(output,0);i++)
			for(j = 0; j<dimsize(output,1);j++)
				for(k=0; k<dimsize(output,2);k++)
					String temp = "en"+num2str(Indextoscale(output,k,2))
					wave temp1 = $temp
					if(indextoscale(output,i,0) > Max(Dimoffset($temp,0), Dimoffset($temp,0)+(dimsize($temp,0)-1)*dimdelta($temp,0)) || indextoscale(output,i,0) < Min(Dimoffset($temp,0), Dimoffset($temp,0)+(dimsize($temp,0)-1)*dimdelta($temp,0)) || indextoscale(output,j,1) > Max(Dimoffset($temp,1), Dimoffset($temp,1)+(dimsize($temp,1)-1)*dimdelta($temp,1)) || indextoscale(output,j,1) < Min(Dimoffset($temp,1), Dimoffset($temp,1)+(dimsize($temp,1)-1)*dimdelta($temp,1)))
						output[i][j][k] =0 // out of range
					
					else
					
					output[i][j][k] = temp1(indextoscale(output,i,0))(indextoscale(output,j,1))[0]
					
					endif
					
					
				endfor
			endfor
		endfor
		
		for(i = 0; i< dimsize(inputwave,2);i++)
			ekcut = "en"+num2str(Indextoscale(inputwave,i,2))
			KillWaves $ekcut
			
		
		endfor
		
		
		
		
End







Function temperature_dep(inputwave, outputname,photonenergy,startenergy,energyrange)
		//20220424 by CZJ
		// This is for combiming different temperature dependence energy cuts.
		// The input should be a three dimensional wave with dimension : 0 -> angle, 1-> energy, 2-> temperature_index.
		// output wave's name is the string outputname. Its a three dimensional wave

		
		Wave inputwave
		String outputname
		Variable photonenergy,startenergy, energyrange
		wave CL_4f_pos
		
		
		// prepare different rescaled energy cuts
		
		Variable i
		Make /O /N = (dimsize(inputwave,2)) workfunction
		wave workfunction
		variable HighEn, LowEn, deltaEn
		
		for(i = 0; i< dimsize(inputwave,2);i++)
			
			String ekcut = "en"+num2str(Indextoscale(inputwave,i,2))
			Duplicate /O /RMD = [][][i] inputwave , $ekcut // remember here
			Setscale /P y,dimoffset($ekcut,1)-116.494-(CL_4f_pos(i)-CL_4f_pos(101)),dimdelta($ekcut,1),$ekcut
			//workfunction[i] = fixfermilevel($ekcut,startenergy,energyrange)
			Setscale /P x, (pi/180)*0.5123*sqrt(photonenergy - 4.5)*dimoffset($ekcut,0), (pi/180)*0.5123*sqrt(photonenergy - 4.5)*dimdelta($ekcut,0), $ekcut	
			
			if(i == 0)
				HighEn = Max(Dimoffset($ekcut,1), Dimoffset($ekcut,1)+(dimsize($ekcut,1)-1)*dimdelta($ekcut,1))
				LowEn = Min(Dimoffset($ekcut,1), Dimoffset($ekcut,1)+(dimsize($ekcut,1)-1)*dimdelta($ekcut,1))
				deltaEn = abs(dimdelta($ekcut,1))
				
			else
				HighEn = Max(HighEn,Max(Dimoffset($ekcut,1), Dimoffset($ekcut,1)+(dimsize($ekcut,1)-1)*dimdelta($ekcut,1)) )
				LowEn = Min(LowEn,Min(Dimoffset($ekcut,1), Dimoffset($ekcut,1)+(dimsize($ekcut,1)-1)*dimdelta($ekcut,1)) )
				deltaEn = Min(abs(dimdelta($ekcut,1)), deltaEn)
				
				
			endif
			
			
			
		endfor
		
		//Construct a 3D wave for output
		//First need to figure out the size along each direction
		
		Make /O /N = (dimsize(inputwave,0),(ceil((HighEn-LowEn)/deltaEn+1)), dimsize(inputwave,2)) , $outputname
		wave output = $outputname
		
		//setupscale for outputwave
		
		setscale /P x,  (pi/180)*0.5123*sqrt(photonenergy - 4.5)*dimoffset(inputwave,0), (pi/180)*0.5123*sqrt(photonenergy - 4.5)*dimdelta(inputwave,0), output
		setscale /I y, LowEn,  HighEn, output
		setscale  /P z, dimoffset(inputwave,2), dimdelta(inputwave,2), output
		
		//From output wave find the value in the original wave
		
		Variable j, k 
		
		for(i=0;i<dimsize(output,0);i++)
			for(j = 0; j<dimsize(output,1);j++)
				for(k=0; k<dimsize(output,2);k++)
					String temp = "en"+num2str(Indextoscale(output,k,2))
					wave temp1 = $temp
					if(indextoscale(output,i,0) > Max(Dimoffset($temp,0), Dimoffset($temp,0)+(dimsize($temp,0)-1)*dimdelta($temp,0)) || indextoscale(output,i,0) < Min(Dimoffset($temp,0), Dimoffset($temp,0)+(dimsize($temp,0)-1)*dimdelta($temp,0)) || indextoscale(output,j,1) > Max(Dimoffset($temp,1), Dimoffset($temp,1)+(dimsize($temp,1)-1)*dimdelta($temp,1)) || indextoscale(output,j,1) < Min(Dimoffset($temp,1), Dimoffset($temp,1)+(dimsize($temp,1)-1)*dimdelta($temp,1)))
						output[i][j][k] =0 // out of range
					
					else
					
					output[i][j][k] = temp1(indextoscale(output,i,0))(indextoscale(output,j,1))[0]
					
					endif
					
					
				endfor
			endfor
		endfor
		
		for(i = 0; i< dimsize(inputwave,2);i++)
			ekcut = "en"+num2str(Indextoscale(inputwave,i,2))
			KillWaves $ekcut
			
		
		endfor
		
		
		
		
End




Function photon2kz(PhotonEdep,innerpotential, workfunction)


//here photonEdep is a 3D wave:(output from function photonenergydep())
//dim0 is binding energy(in eV) and fermi level should be already aligned
//dim1 is momentum(in inverse A) and Gamma zero needs to be at 0
//dim2 is photon energy(in eV)
//innerpotential is inner potential in eV(should be a postive number)
//workfunction: in eV and is positive
//Output: 3D rescaled wave with name:(nameofwave(PhotonEdep)+"_kzkx")

//CZJ 20210206
//Here we first negelect the effect of binding energy on kz, which is true when binding energy is small comparing to other energy scales.

wave PhotonEdep
variable innerpotential
variable workfunction

//2D normalize inputwave to get rid of possible intensity modulation with photon energy
ZJ_2dnorm(PhotonEdep,2)
wave swave = $(nameofwave(PhotonEdep)+"_norm_2d")

Make /O /N = (dimsize(PhotonEdep,0),dimsize(PhotonEdep,1),5*dimsize(PhotonEdep,2)) $(nameofwave(PhotonEdep)+"_kzkx")
wave temp = $(nameofwave(PhotonEdep)+"_kzkx")

variable kz_min,kz_max,E_max,E_min
//Set the largest Kz and smallest Kz
E_max = max(indexToScale(PhotonEdep,0,2),indextoscale(PhotonEdep,dimsize(PhotonEdep,2)-1,2))
E_min = min(indexToScale(PhotonEdep,0,2),indextoscale(PhotonEdep,dimsize(PhotonEdep,2)-1,2))
kz_max = 0.5123*sqrt(E_max-workfunction+innerpotential)
kz_min = 0.5123*sqrt(E_min-workfunction+innerpotential-max(indexToScale(PhotonEdep,0,1)^2,indexToScale(PhotonEdep,dimsize(PhotonEdep,1)-1,1)^2)*(0.5123)^(-2))


Setscale /P x, dimoffset(PhotonEdep,0),dimdelta(PhotonEdep,0),temp
Setscale /P y, dimoffset(PhotonEdep,1),dimdelta(PhotonEdep,1),temp
Setscale /I z, kz_min,kz_max,temp



//wave assignment, notice that igor pro doesn't use interpolation when using scale for reference in multidimensional waves. So here we need to do it by ourselves. 
multiThread temp[][][] = ((0.5123)^(-2)*(y^2+z^2)+workfunction-innerpotential<=E_max && (0.5123)^(-2)*(y^2+z^2)+workfunction-innerpotential>=E_min)? interpindex3d_z(swave,x,y,(0.5123)^(-2)*(y^2+z^2)+workfunction-innerpotential) : 0



End

ThreadSafe Function interpindex3d_z(swave,x,y,z)

//only interpolate along z
wave swave
variable x,y,z

variable diff,weight

diff = (z - indextoscale(swave,scaletoIndex(swave,z,2),2))
weight = diff/(indextoscale(swave,scaletoIndex(swave,z,2)+sign(diff),2) - indextoScale(swave,scaletoIndex(swave,z,2),2))

return weight*swave(x)(y)[scaletoIndex(swave,z,2)+sign(diff)]+(1-weight)*swave(x)(y)[scaletoIndex(swave,z,2)]


End


Function ZJ_2dnorm(swave,dimension_norm)
//This is a 2D norm on 3 dimensional waves
//Default outputwave name is nameofwave(swave)+"_norm_2d"

wave swave
variable dimension_norm

Duplicate /O swave, $(nameofwave(swave)+"_norm_2d")
wave temp_wave = $(nameofwave(swave)+"_norm_2d")

variable temp_sum = 0
variable i,j,k

if(dimension_norm == 0)
	for(k=0;k<dimsize(swave,0);k++) 
		for(i=0;i<dimSize(swave,1);i++)
			for(j=0;j<dimsize(swave,2);j++)
				temp_sum +=swave[k][i][j]	
			endfor
		endfor
		temp_wave[k][][] = (swave[k][q][r]==0)?0 : swave[k][q][r]/temp_sum
		temp_sum = 0
	endfor	
elseif(dimension_norm == 1)
	for(k=0;k<dimsize(swave,1);k++) 
		for(i=0;i<dimSize(swave,0);i++)
			for(j=0;j<dimsize(swave,2);j++)
				temp_sum +=swave[i][k][j]	
			endfor
		endfor
		temp_wave[][k][] = (swave[p][k][r]==0)?0 : swave[p][k][r]/temp_sum
		temp_sum = 0
	endfor			
	
elseif(dimension_norm == 2)
	for(k=0;k<dimsize(swave,2);k++) 
		for(i=0;i<dimSize(swave,1);i++)
			for(j=0;j<dimsize(swave,0);j++)
				temp_sum +=swave[j][i][k]	
			endfor
		endfor
		temp_wave[][][k] = (swave[p][q][k]==0)?0: swave[p][q][k]/temp_sum
		temp_sum = 0
	endfor		
	
	
else
	print("wrong dimension")

endif



End