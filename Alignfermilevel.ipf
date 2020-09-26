#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.


Function fixfermilevel(ekcut, startenergy, energyrange)
	//20191210 this is a part of the program used for align the fermi level for photon energy dependence, etc.
	//ekcut should be a 2d wave with dimension 0: angle/momentum; dimension 1: energy.
	//Method: find the largest derivative and that is the fermi energy
	//The fermi level should within startenergy+- energyrange
	
	wave ekcut
	variable startenergy, energyrange
	Make /O /N = (dimsize(ekcut,1)), energyspec = 0
	Variable i,j, fermi
	
	
	// make energy spectrum
	for(i=0;i < dimsize(ekcut,1);i++)
		for(j = 0; j  < dimsize(ekcut,0);j++)
			energyspec[i] =  ekcut[j][i] + energyspec[i]
		endfor		
	endfor
	
	Setscale /P x , dimoffset(ekcut,1), dimdelta(ekcut,1), energyspec 
	Duplicate /O /R = (startenergy-energyrange, startenergy+energyrange) energyspec, $"test"
	
	
	//method: use differentiate to get the largest change point in the integerated spectrum:
	//seems we need to narrow down the energy range.
	
	Differentiate/METH=1 $"test"/D=energyspec__DIF;DelayUpdate
	wave energyspec__DIF
	
	energyspec__DIF[] = abs(energyspec__DIF[p])
	
	Findvalue /V = (Wavemax(energyspec__DIF)) energyspec__DIF
	fermi = Indextoscale(energyspec__DIF,V_value,0)
	
	//reset scaling of ekcutprin
	Setscale /P y, dimoffset(ekcut,1) - fermi, dimdelta(ekcut,1), ekcut
	killwaves $"test"
	
	return fermi
	
	
End


Function Photonenergydep1(inputwave, outputname,startenergy,energyrange)
		//20191210 by CZJ
		// This is for combiming different photon energy cuts. so we need to rescale kx at each photon energy and also align fermi levels.
		// The input should be a three dimensional wave with dimension : 0 -> angle, 1-> energy, 2-> light energy.
		// output wave's name is the string outputname. Its a three dimensional wave
		//The fermi level should within startenergy+- energyrange
		
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
			print i
		endfor
		
		for(i = 0; i< dimsize(inputwave,2);i++)
			ekcut = "en"+num2str(Indextoscale(inputwave,i,2))
			KillWaves $ekcut
			
		
		endfor
		
		
		
		
End


Function Photonenergydep2(wavenamelist,photoenergylist,outputwavename)

//The input wavelist should be: the name of 2d cuts, which dim0 is kinetic energy and dim1 is angle(not momentum)
//photoenergylist is the list of photonenergies
//This function is trying to solve the problem if different energy cuts has different energy scaling and angle range.So each cut will be rescaled based on its own range
//NOTE: this function will change the scaling of the waves in the wavenamelist, so strongly suggest to duplicate waves first.

wave /T wavenamelist
wave photoenergylist
String outputwavename

variable workfunction = 4.35

//error handling
if(dimsize(wavenamelist,0)!=dimsize(photoenergylist,0))
	print("dimensions of input waves are not right")
	Abort
endif

//rescale each cut, can be commented if the cuts are already scaled. Note this will scaling on the original wave
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

//perhaps here you need to determine the real photon energy position by fitting
Setscale /P x,dimoffset($wavenamest,0)-(photonenergy - workfunction), dimdelta($wavenamest,0), $wavenamest
Setscale /P y,0.5132*sqrt(photonenergy - workfunction)*Pi/180*(dimoffset($wavenamest,1)),0.5132*sqrt(photonenergy - workfunction)*Pi/180*(dimdelta($wavenamest,1)),$wavenamest

//Needs to more rigoriously fix the fermi energy position

variable E_fermi
Matrixtranspose $wavenamest
E_fermi = fixfermilevel($wavenamest,0,0.5)
Matrixtranspose $wavenamest



End



