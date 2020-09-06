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


Function alignfs(inputwave, outputname,startenergy,energyrange)
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
