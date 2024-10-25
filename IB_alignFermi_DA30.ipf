#pragma rtGlobals=3		// Use modern global access method and strict wave access.

// An attempt at improving Fermi level correction for Fermi surface maps acquired by DA30, suffering from the 'Fermi level distortion'.
//
// 1. IB_bouquetEDC(): Bin bouquets of EDCs (similar to Tyler's DataCondense(), but adapted to the DA30 Fermi surface problem here)
// 2. IB_fit_Fermi_edge_vol(): Fit all EDCs, producing a set of coarse Fermi levels
// 3. IB_unbin_efTrack(): 'Unbin', producing a set of fine Fermi levels that map 1-to-1 to the original Fermi surface
// 4. IB_ef_align(): Correct the Fermi levels for each EDC
//
//
// Example:
//
// binned_data = IB_bouquetEDC(your_raw_data,5,10,0.1)
// efTrack_for_the_binned_data = IB_fit_Fermi_edge_vol(binned_data)
// efTrack_for_your_raw_data = IB_unbin_efTrack(your_raw_data,efTrack_for_the_binned_data)
// aligned_data = IB_ef_align(your_raw_data,efTrack_for_your_raw_data)
//
//
// Dependencies:
// 
// avgDistCurve3D()
// IB_fit_Fermi_edge()
//
//
// Notes:
//
// Assumes that the energy axis is along Dimension 2 and the two angles are along Dimension 0 and Dimension 1.
//
// You will probably need to modify some of the fitting parameters in IB_guess_fermi(), IB_fit_Fermi_edge().
//
// IB_fit_Fermi_edge() could be further improved to recover more data near the boundary of the D30 field of view.
//
// Based on an early version by Tyler, c. 2019.
//
//
// Ilya
//
// 20210124 Created


Function/S IB_bouquetEDC(raw,binning_factorX,binning_factorY,cutoffPercent)

	// raw: three-dimensional wave, presumably a Fermi surface, currently assumed to have the energy axis along Dimension 2
	// binning_factorX: variable, bins along Dimension 0; binning_factorX = 5 means bins {0,1,2,3,4}, {5,6,7,8,9} and so on... The remainder will be grouped into a leftover bin, e.g ... {100,101,102}.
	// binning_factorY: variable, bins along Dimension 1
	// cutoffPercent: variable, threshold for declaring an EDC to be 'dark' (outside of the DA30 field of view)
	
	Wave raw
	Variable binning_factorX,binning_factorY
	Variable cutoffPercent

	Variable verbose = 0

	if (verbose)
		Variable time_id
		Variable time_out
		time_id = StartMSTimer
	endif

	// Collapse fully along the energy dimension
	String collapse_s = GetWavesDataFolder(raw,2) + "_collapse"
	MatrixOp /NTHR=0 /O $collapse_s = sumBeams(raw)
	Wave collapse = $collapse_s

	// Total counts
	String sumWave_s = GetWavesDataFolder(raw,2) + "_sum"
	MatrixOp /NTHR=0 /O $sumWave_s = sum(collapse)
	Wave sumWave = $sumWave_s

	Variable cutoffIntensity = cutoffPercent*sumWave[0]/(DimSize(raw,0)*DimSize(raw,1)) // cutoff number of counts
	
	// Create a copy of the data where we will set EDCs below cutoffIntensity to zero
	String zeroed_s = GetWavesDataFolder(raw,2) + "_zeroed"
	Duplicate /O raw, $zeroed_s
	Wave zeroed = $zeroed_s
	
	// Keep the EDC or not?
	collapse[][] = (collapse[p][q] > cutoffIntensity) ? 1 : 0
	MultiThread zeroed[][][] = (collapse[p][q]) ? zeroed[p][q][r] : 0 

	// Bins in sets of binning_factorX along Dimension 0, keeping the remainder in a last 'incomplete' bin; same for binning_factorY
	String bin_s = GetWavesDataFolder(raw,2) + "_bin"
	Make /O/N=(ceil(DimSize(raw,0)/binning_factorX),ceil(DimSize(raw,1)/binning_factorY),DimSize(raw,2)) $bin_s
	Wave bin = $bin_s

	SetScale /P x, DimOffset(raw,0)+(DimDelta(raw,0)*(binning_factorX-1)/2), DimDelta(raw,0)*binning_factorX, bin
	SetScale /P y, DimOffset(raw,1)+(DimDelta(raw,1)*(binning_factorY-1)/2), DimDelta(raw,1)*binning_factorY, bin
	SetScale /P z, DimOffset(raw,2), DimDelta(raw,2), bin
		
	String bin_EDC_s = bin_s + "_DC"

	Variable i,j

	// Bin it up
	for (i = 0; i < DimSize(bin,0); i += 1)
		for (j = 0; j < DimSize(bin,1); j += 1)

			bin_EDC_s = avgDistCurve3D(GetWavesDataFolder(zeroed,2),2,binning_factorX*i,binning_factorX*i+binning_factorX-1,binning_factorY*j,binning_factorY*j+binning_factorY-1)
			Wave bin_EDC = $bin_EDC_s
			bin[i][j][] = bin_EDC[r]
									
		endfor
	endfor

	if (verbose)
		time_out = StopMSTimer(time_id)
		print "IB_bouquetEDC: ", time_out/10e6, "seconds."
		print "IB_bouquetEDC: ", bin_s
	endif
	
	if (!verbose)
		KillWaves /Z collapse, sumWave, zeroed, bin_EDC
	endif
	
	return bin_s
	
End


// The below snippet does the same thing as
// bin_EDC_s = avgDistCurve3D(GetWavesDataFolder(zeroed,2),2,binning_factorX*i,binning_factorX*i+binning_factorX-1,binning_factorY*j,binning_factorY*j+binning_factorY-1)
// but seems to be ~3 times slower:

//		Variable rs,re,cs,ce
//
//			rs = binning_factorX*i
//			re = binning_factorX*i+binning_factorX-1
//			cs = binning_factorY*j
//			ce = binning_factorY*j+binning_factorY-1
//
//			if (re > DimSize(zeroed,0)-1)
//				re = DimSize(zeroed,0)-1
//			endif
//
//			if (ce > DimSize(zeroed,1)-1)
//				ce = DimSize(zeroed,1)-1
//			endif
//			
//			MatrixOp /NTHR=0 /O $bin_EDC_s = sum(subRange(zeroed,rs,re,cs,ce))


Function/S IB_fit_Fermi_edge_vol(fs)

	// fs: three-dimensional wave, the Fermi surface, presumably the binned output of IB_bouquetEDC() above; assume energy axis along Dimension 2
	//
	// Wrapper script for IB_fit_Fermi_edge()
	
	Wave fs
	
	Variable verbose = 0

	if (verbose)
		Variable time_id
		Variable time_out
		time_id = StartMSTimer
	endif

	// Destination wave for the extracted Fermi levels
	String ef_s = GetWavesDataFolder(fs,2) + "_efTrack"
	Make /O/N=(DimSize(fs,0),DimSize(fs,1)) $ef_s
	Wave ef_w = $ef_s

	SetScale /P x, DimOffset(fs,0), DimDelta(fs,0), ef_w
	SetScale /P y, DimOffset(fs,1), DimDelta(fs,1), ef_w

	// Fit the EDC or not? Avoid attempting to fit border regions which have no overlap with data and are expected to have strictly zero counts.
	String sum_s = GetWavesDataFolder(fs,2) + "_sum"
	MatrixOp /NTHR=0 /O $sum_s = sumBeams(fs)
	Wave sumWave = $sum_s

	// Temporary EDC as a one-dimensional wave
	String EDC_s = GetWavesDataFolder(fs,2) + "_EDC"
	Make /O/N=(DimSize(fs,2)) $EDC_s
	Wave EDC = $EDC_s
	
	SetScale /P x, DimOffset(fs,2), DimDelta(fs,2), EDC

	String fit_message
	Variable i,j

	for (i = 0; i < DimSize(fs,0); i += 1)
//		for (j = 0; j < 1; j += 1)
		for (j = 0; j < DimSize(fs,1); j += 1)

//			if (i == 14 && j == 18) // Troubleshoot one specific EDC

			EDC[] = fs[i][j][p]

			if (sumWave[i][j])

				fit_message = IB_fit_Fermi_edge(EDC)
				
				if (!cmpstr(fit_message,"error"))
					print "IB_fit_Fermi_edge_vol: fitting error for index ", i, j, "of ", GetWavesDataFolder(fs,2), "."
				endif
				
				Wave fitEFcoeff = $(EDC_s + "_coeffF")		
			
				if (fitEFcoeff[0] > DimOffset(EDC,0) && fitEFcoeff[0] < DimOffset(EDC,0) + DimSize(EDC,0)*DimDelta(EDC,0))
								
					eF_w[i][j] = fitEFcoeff[0]
	
				else
					
					eF_w[i][j] = NaN
	
				endif
						
			//	print fitEFcoeff[0]
		
			else
			
				eF_w[i][j] = NaN
			
			endif
			
//			endif // Troubleshoot one specific EDC
			
		endfor
	endfor

	if (verbose)
		time_out = StopMSTimer(time_id)
		print "IB_fit_Fermi_edge_vol:", time_out/10e6, "seconds"
		print "IB_fit_Fermi_edge_vol:", ef_s
	endif
	
	if (!verbose)
		KillWaves /Z sumWave
	endif

	return ef_s

End


Function/S IB_unbin_efTrack(raw,efTrack)

	// raw: three-dimensional wave, presumably a Fermi surface, currently assumed to have the energy axis along Dimension 2
	// efTrack: tracking of the binned Fermi surface, presumably output of IB_fit_Fermi_edge_vol()

	// This is obviously a poorly-defined procedure. Many approaches may be considered: fitting, resampling, interpolation, smoothing...

	Wave raw
	Wave efTrack

	Variable verbose = 0

	if (verbose)
		Variable time_id
		Variable time_out
		time_id = StartMSTimer
	endif

	// _efTrack for the raw DA30 Fermi surface
	String ef_s = GetWavesDataFolder(raw,2) + "_efTrack"
	Make /O/N=(DimSize(raw,0),DimSize(raw,1)) $ef_s
	Wave ef_w = $ef_s

	SetScale /P x, DimOffset(raw,0), DimDelta(raw,0), ef_w
	SetScale /P y, DimOffset(raw,1), DimDelta(raw,1), ef_w


// Using ImageInterpolate:

//	Variable binning_factorX = 10, binning_factorY = 10

//	Variable nx = DimSize(efTrack,0)*binning_factorX
//	Variable ny = DimSize(efTrack,1)*binning_factorY

//	ImageInterpolate /DEST=root:efTrack /D=2 /RESL={nx,ny} Spline efTrack
//	Wave ef_w = root:efTrack

//	SetScale /P x, DimOffset(raw,0), DimDelta(raw,0), ef_w
//	SetScale /P y, DimOffset(raw,1), DimDelta(raw,1), ef_w


// Using Resample:

//	Resample /UP=(binning_factorX) efTrack
//	Resample /UP=(binning_factorY)/DIM=1 efTrack


//	Using Loess?


// Using interp2D:

	ef_w[][] = interp2D(efTrack,DimOffset(ef_w,0)+p*DimDelta(ef_w,0),DimOffset(ef_w,1)+q*DimDelta(ef_w,1)) // these are scaled values


// Using direct unbinning, creating 'plateaus' with the same value for each bin: 

//	ef_w[][] = efTrack[floor(p/binning_factorX)][floor(q/binning_factorY)] // these are scaled values


// Using Smooth:

//	Smooth /B 7, ef_w
//	Smooth /DIM=1/B 7, ef_w


	ef_w[][] = (ef_w[p][q]-DimOffset(raw,2))/DimDelta(raw,2) // convert to index values


// Using polynomial fitting:

//	String ef_fit_s = ef_s + "_fit"
//	Duplicate /O ef_w, $ef_fit_s
//	
//	Make/D/N=15/O W_coef
//	W_coef[0] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}
//	FuncFitMD/NTHR=0 TwoDQaud W_coef ef_w /D=$ef_fit_s
	
	if (verbose)
		time_out = StopMSTimer(time_id)
		print "IB_unbin_efTrack:", time_out/10e6, "seconds"
		print "IB_unbin_efTrack:", ef_s
	endif
	
	return ef_s
	
End


Function/S IB_ef_align(fs,ef)

	// fs: raw three-dimensional wave, presumably a Fermi surface, currently assumed to have the energy axis along Dimension 2
	// ef: tracks the Fermi level of the raw Fermi surface, presumably output of IB_unbin_efTrack()

	Wave fs
	Wave ef

	Variable verbose = 1

	if (verbose)
		Variable time_id
		Variable time_out
		time_id = StartMSTimer
	endif
		
	Variable fermiMax = WaveMax(ef) // indices here, not scaled values
	Variable fermiMin = WaveMin(ef)
	
	//	print fermiMax
	//	print fermiMin

	Variable extra = trunc(fermiMax-fermiMin) + 1 // adding an extra pixel +1 here, should avoid throwing away ANY data (save every pixel!)
	
	//	print extra
	
	// Destination wave for the aligned Fermi surface, enlarged along the energy axis to make room for sliding all EDCs relative to one another.
	String fs_ef_s = GetWavesDataFolder(fs,2) + "_ef"		
	Make /O/N=(DimSize(fs,0),DimSize(fs,1),DimSize(fs,2)+extra) $fs_ef_s
	Wave fs_ef = $fs_ef_s
	
	SetScale /P x, DimOffset(fs,0), DimDelta(fs,0), fs_ef
	SetScale /P y, DimOffset(fs,1), DimDelta(fs,1), fs_ef
	SetScale /P z, -fermiMax*DimDelta(fs,2), DimDelta(fs,2), fs_ef

	// The core alignment step: [ef[p][q] - fermiMax + r] in fs should map to [r] in fs_ef
	MultiThread fs_ef[][][] = ((ef[p][q] - fermiMax + r >= 0) && (ef[p][q] - fermiMax + r <= DimSize(fs,2)-1)) ? interp3DAlongDirection(fs,p,q,ef[p][q]-fermiMax+r) : 0

	if (verbose)
		time_out = StopMSTimer(time_id)
		print "IB_ef_align:", time_out/10e6, "seconds"
		print "IB_ef_align:", fs_ef_s
	endif

	return fs_ef_s

End


ThreadSafe Static Function interp3DAlongDirection(wave_in,i,j,index)

	// Shouldn't Igor have something like this built-in???
	
	// Along dimension 2
	// To do: generalize to any dimension

	// wave_in: three-dimensional wave
	// i: variable, integer index for Dimension 0
	// j: variable, integer index for Dimension 1
	// index: variable, possibly non-integer index for Dimension 2

	// Returns a number

	// Ilya
	
	// 20210124 Created

	Wave wave_in
	Variable i,j,index // points, not scaled values
	
	Variable dim = 2
	Variable wave_size = DimSize(wave_in,dim)

	Variable val
	Variable fracDown, fracUp

	Variable floor_index

	if (numtype(index)) // treat NaNs and Infs separately?
		
		val = NaN
		
	elseif (index < 0)
		
		val = wave_in[i][j][0]
		
	elseif (index > wave_size-1)
	
		val = wave_in[i][j][wave_size-1]
	
	else
	
		// Example: index = 5.6 => fracDown = 5-5.6+1 = -0.6+1 = 0.4 && fracUp = 5.6-5 = 0.6
		// Round number example: index = 4 => fracDown = 4-4+1 = 1 && fracUp = 4-4 = 0

		floor_index = floor(index)
	
		fracDown = floor_index-index+1
		fracUp = index-floor_index
		
		val = fracDown*wave_in[i][j][floor_index] + fracUp*wave_in[i][j][ceil(index)] // linear interpolation
	
	endif
	
	Return val
	
End


// Without multi-threading, actually the above function is still somewhat slower than a loop like this with a temporary "EDC" wave and Igor's internal 1D interpolation
// But the drawback below is cannot multi-thread... This all could probably be further improved.

	for (i = 0; i < DimSize(source,0); i += 1)
		for (j = 0; j < DimSize(source,1); j += 1)

			fermi = ef[i][j]
			
			if (numtype(fermi) == 2) // NaN
		
				dest[i][j][] = 0
				
			else
		
				EDC[] = source[i][j][p]
				dest[i][j][] = ((fermi-fermiMax+r >= 0) && (fermi-fermiMax+r <= DimSize(source,2)-1)) ? EDC[fermi-fermiMax+r] : 0

			endif

		endfor
	endfor


// Random version with ceil() rather than floor()

		// Example: index = 5.6 => fracDown = 6-5.6 = 0.4 && fracUp = 5.6-6+1 = -0.4+1 = 0.6
		// Round number example: index = 4 => fracDown = 4-4 = 0 && fracUp = 4-4+1 = 1
	
//		fracDown = ceil(index)-index
//		fracUp = index-ceil(index)+1
		
//		val = fracDown*wave_in[i][j][floor(index)] + fracUp*wave_in[i][j][ceil(index)] // linear interpolation
	


Function TwoDQaud(w,x,y) : FitFunc
	Wave w
	Variable x
	Variable y

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x,y) = C1+C2*y+C3*y^2+C4*x+C5*x*y+C6*x*y^2+C7*x^2+C8*x^2*y+C9*x^2*y^2
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 2
	//CurveFitDialog/ x
	//CurveFitDialog/ y
	//CurveFitDialog/ Coefficients 9
	//CurveFitDialog/ w[0] = C1
	//CurveFitDialog/ w[1] = C2
	//CurveFitDialog/ w[2] = C3
	//CurveFitDialog/ w[3] = C4
	//CurveFitDialog/ w[4] = C5
	//CurveFitDialog/ w[5] = C6
	//CurveFitDialog/ w[6] = C7
	//CurveFitDialog/ w[7] = C8
	//CurveFitDialog/ w[8] = C9

//	return w[0]+w[1]*y+w[2]*y^2+w[3]*x+w[4]*x*y+w[5]*x*y^2+w[6]*x^2+w[7]*x^2*y+w[8]*x^2*y^2
	return w[0]+w[1]*y+w[2]*y^2+w[3]*x+w[4]*x*y+w[5]*x*y^2+w[6]*x^2+w[7]*x^2*y+w[8]*x^2*y^2+w[9]*x^3+w[10]*y^3+w[11]*x^4+w[12]*x^3*y+w[13]*x*y^3+w[14]*y^4

End



////////////////////////////////////////////////////////////////////////

//Example: smoothderivoutput3d(tmfs1_norm_0,thetest,0.1,3,5,685,830)
//
//Instructions:
//1. Make an outputwave. The dimensions don't have to match, just make a wave with the desired name.
//2. Cutoffpercentage take edcs whos total intensity is less that cutoffpercentage of the average and sets the Fermi level the average Fermi level. Should only apply to the border regions.
//3. Play with SmoothRange and SmoothTimes if problems.
//4. Identify the max and minimum points the Fermi level could be (EfMinPt, EfMaxPt) to eliminate effects from strong deep bands.
//
////////////////////////////////////////////////////////////////////////


//dim0 is scanning angle
//dim1 is analyzer angle
//dim2 is binding energy


// From Tyler, 2019?
// 20210118 Ilya, updated to align by fractional index, so step/plateau artifacts are avoided

Function SmoothDerivOutput3D_ignore(InputWave,OutputWave,CutOffPercentage,SmoothRange, SmoothTimes,EfMinPt,EfMaxPt)
	Wave InputWave
	Wave OutputWave
	Variable CutOffPercentage
	Variable SmoothRange
	Variable SmoothTimes
	Variable EfMinPt
	Variable EfMaxPt
	
	variable v_sum
	wavestats /q InputWave
	Variable CutOffIntensity=v_sum/(dimsize(inputwave,0)*dimsize(inputwave,1))*CutOffPercentage
	
	String EDCName = "SingleEDC"
	Make /O /N=(dimsize(inputwave,2)) $EDCName
	SetScale /P x, Dimoffset(inputwave,2), DimDelta(inputwave,2), $EDCName // (20210118)
	Wave EDC = $EDCName

	 String EfIndexName=NameofWave(InputWave)+"EfIndex"
	 Make /O /N=(dimsize(inputwave,0),dimsize(inputwave,1)) $EfIndexName
	 Wave EfIndex=$EfIndexName
	 
	 setscale /p x, dimoffset(inputwave,0),dimdelta(inputwave,0),efindex
	 setscale /p y, dimoffset(inputwave,1),dimdelta(inputwave,1),efindex
	 
	String TempEfValueWaveName = "TempEf"
	Make /O /N=1 $TempEfValueWaveName
	Wave TempEfValue = $TempEfValueWaveName
	
	variable i,j

	for(i=0;i<dimsize(inputwave,0);i+=1)
	for(j=0;j<dimsize(inputwave,1);j+=1)
		EDC[]=InputWave[i][j][p]
		wavestats /q EDC
		If(v_sum>cutoffintensity)
			//	IB_fit_Fermi_edge(EDC)
			//	Wave fitEFcoeff = $(EDCName + "_coeffF")
			//	EfIndex[i][j] = (fitEFcoeff[0] - DimOffset(EDC,0))/DimDelta(EDC,0)
			//	EfIndex[i][j] = (EfIndex[i][j] >= 0 && EfIndex[i][j] < DimSize(EDC,0)) ? EfIndex[i][j] : NaN
			EfIndex[i][j]=smoothderivoutput_ignore(EDC,SmoothRange,SmoothTimes,EfMinPt,EfMaxPt) // TempefValue
			//EfIndex[i][j]=TempEfValue[0]
		Else
			EfIndex[i][j]=NaN
		EndIf
	endfor
	endfor

	String FitEfIndexName="fit_"+EfIndexName
	duplicate /o /s EFIndex $FitEfIndexName
	Wave FitEfIndex = $FitEfIndexName
	
	Make/D/N=9/O W_coef
	W_coef[0] = {1,1,1,1,1,1,1,1,1}
//	Make/D/N=15/O W_coef
//	W_coef[0] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}
	FuncFitMD/NTHR=0 TwoDQaud W_coef EfIndex /D=$FitEfIndexName
	
	Variable fermiMax = WaveMax(FitEfIndex) // (20210118) no rounding
	Variable fermiMin = WaveMin(FitEfIndex) // (20210118) no rounding
//	Variable fermiMax=round(WaveMax(FitEfIndex))
//	Variable fermiMin=round(WaveMin(FitEfIndex))

	Variable extra = trunc(fermiMax-fermiMin) + 1 // (20210118) adding this extra pixel should avoid throwing away any data
//	Variable extra = fermiMax - fermiMin
	
	Variable new_offset = -fermimax*dimdelta(inputwave,2)
	Variable step_size = dimdelta(inputwave,2)
	
	String photDepShiftStr = NameofWave(OutputWave)
	Make /O/N=(DimSize(Inputwave,0),Dimsize(inputwave,1),DimSize(Inputwave,2)+extra) $photDepShiftStr
       Wave photDepShift = $photDepShiftStr
       setscale /p x, dimoffset(inputwave,0),dimdelta(inputwave,0),photdepshift
       setscale /p y, dimoffset(inputwave,1),dimdelta(inputwave,1),photdepshift
       setscale /p z, -fermimax*dimdelta(inputwave,2),dimdelta(inputwave,2),photdepshift
             
       variable fermi, start
              
       for(i=0;i<dimsize(inputwave,0);i+=1)
       for(j=0;j<dimsize(inputwave,1);j+=1)
  		
       	if(numtype(fitefindex[i][j])==2)
     			wavestats /q fitefindex
     			fermi=v_avg
            	else
       		fermi=fitefindex[i][j]
       	endif
       	start = fermiMax - fermi // (20210118) no rounding, non-integer index
//       	 start=round(fermimax-fermi)
       	if(start<0)
       		print "start is less than zero"
       	endif
       	if(start+dimsize(inputwave,2)-1>dimsize(photdepshift,2)-1)
       		print "troube in middle region"
       		print start
       		print dimsize(inputwave,2)
       		print dimsize(photdepshift,2)
       		print fermi
       		print fermimin
       		return(0)
       	endif

			EDC[]=InputWave[i][j][p]
			photDepShift[i][j][] = ((fermi-fermiMax+r >= 0) && (fermi-fermiMax+r <= DimSize(Inputwave,2)-1)) ? EDC[fermi-fermiMax+r] : 0 // (20210118) interpolation happening here

      	// photDepShift[i][j][,start] = InputWave[i][j][0]
   	   // photDepShift[i][j][start,start+DimSize(Inputwave,2)-1] = InputWave[i][j][r-start]
	  	   // photDepShift[i][j][start+DimSize(inputwave,2)-1,] = InputWave[i][j][DimSize(InputWave,2)-1]

       endfor
       endfor
		
End



//dim0 is analyzer angle
//dim1 is binding energy
Function SmoothDerivOutput2D_ignore(InputCut,OutputCut,CutOffPercentage,SmoothRange,SmoothTimes,EfMinPt,EfMaxPt)
	Wave InputCut
	Wave OutputCut
	variable CutOffPercentage
	Variable SmoothRange
	Variable SmoothTimes
	Variable EfMinPt
	Variable EFMaxPt
	
	variable v_sum
	wavestats /q InputCut
	Variable CutOffIntensity=v_sum/dimsize(inputcut,0)*CutOffPercentage
	
	String EDCName = "SingleEDC"
	Make /O /N=(dimsize(inputcut,1)) $EDCName
	Wave EDC = $EDCName
	
	String EfIndexName = NameofWave(InputCut)+"_EfIndex"
	Make /O /N = (dimsize(inputcut, 0)) $EFIndexname
	Wave EfIndex = $EfIndexName
	
	setscale /p x, dimoffset(inputcut,0),dimdelta(inputcut,0), EfIndex
	
	String TempEfValueWaveName = "TempEf"
	Make /O /N=1 $TempEfValueWaveName
	Wave TempEfValue = $TempEfValueWaveName
	
	variable i
	
	for(i=0;i<dimsize(inputcut,0);i+=1)
		EDC[]=inputcut[i][p]
		wavestats /q EDC
		If(v_sum>CutOffIntensity)
			EfIndex[i] = SmoothDerivOutput_ignore(EDC,SmoothRange,SmoothTimes,EfMinPt,EfMaxPt) // TempEfValue
			//EfIndex[i]=TempEfValue[0]
		Else
			EfIndex[i]=NaN
		EndIf
	endfor
	
	variable fermiMax=WaveMax(EfIndex)
	variable fermiMin=WaveMin(EfIndex)
	Variable extra = fermiMax-fermiMin
	
	String photDepShiftStr = NameofWave(OutputCut)
	Make /O/N=(DimSize(InputCut,0),DimSize(InputCut,1)+extra) $photDepShiftStr
       Wave photDepShift = $photDepShiftStr
       setscale /p x, dimdelta(Inputcut,0),dimdelta(inputcut,0),photDepShift
     	SetScale /P y, -DimDelta(InputCut,1)*fermiMax, DimDelta(InputCut,1), photDepShift
     	
     	variable fermi, start
     	
     	for (i = 0; i < DimSize(InputCut,0); i += 1)
     		If(numtype(EfIndex[i])==2)
     			wavestats /q efindex
     			fermi=v_avg
     		Else
     			fermi = EfIndex[i]
     		EndIf
                start = fermiMax - fermi
                photDepShift[i][,start] = InputCut[i][0]
                photDepShift[i][start,start+DimSize(InputCut,1)-1] = InputCut[i][q-start]
                photDepShift[i][start+DimSize(InputCut,1)-1,] = InputCut[i][DimSize(InputCut,1)-1]
        endfor
	
End


Function SmoothDerivOutput_ignore(InputEDC,SmoothRange,SmoothTimes,EfMinPt,EfMaxPt)
// Function SmoothDerivOutput(InputEDC,Result,SmoothRange,SmoothTimes,EfMinPt,EfMaxPt)

	wave inputEDC
//	wave result //This will be the kinetic energy that corresponds to minimum derivative.
	variable SmoothRange
	Variable SmoothTimes
	Variable EfMinPt
	Variable EfMAxPt
	
	String SmoothInputWaveName = "swave"
	Duplicate /O InputEDC $SmoothInputWaveName
	wave SmoothInputWave = $SmoothInputWaveName
	
	variable i
	for(i=0;i<SmoothTimes;i+=1)
		Smooth SmoothRange, SmoothInputWave
	endfor
	
	String DiffWaveName = "dwave"
	Differentiate SmoothInputWave /D=$DiffWaveName
	Wave DiffWave = $DiffWaveName
	
	String EDCpartname = "windowEDC"
	make /o /n=(EfMaxPt-EfMinPt+1) $EDCpartname
	wave EDCpart = $EDCpartname
	EDCpart[]=DiffWave[EfMinPt+p]
	
	variable v_minloc
	wavestats /q $EDCpartname
	
	//Result[0] = v_minloc+EfMinPt
	
	return v_minloc+EfMinPt
	
End




// From Zijia, version without the quadratic fit 20210119
Function SmoothDerivOutput3D_nonfit_ignore(InputWave,OutputWave,CutOffPercentage,SmoothRange, SmoothTimes,EfMinPt,EfMaxPt)
	Wave InputWave
	Wave OutputWave
	Variable CutOffPercentage
	Variable SmoothRange
	Variable SmoothTimes
	Variable EfMinPt
	Variable EfMaxPt
	
	variable v_sum
	wavestats /q InputWave
	Variable CutOffIntensity=v_sum/(dimsize(inputwave,0)*dimsize(inputwave,1))*CutOffPercentage
	
	String EDCName = "SingleEDC"
	Make /O /N=(dimsize(inputwave,2)) $EDCName
	Wave EDC = $EDCName
	
	 String EfIndexName=NameofWave(InputWave)+"EfIndex"
	 Make /O /N=(dimsize(inputwave,0),dimsize(inputwave,1)) $EfIndexName
	 Wave EfIndex=$EfIndexName
	 
	 setscale /p x, dimoffset(inputwave,0),dimdelta(inputwave,0),efindex
	 setscale /p y, dimoffset(inputwave,1),dimdelta(inputwave,1),efindex
	 
	String TempEfValueWaveName = "TempEf"
	Make /O /N=1 $TempEfValueWaveName
	Wave TempEfValue = $TempEfValueWaveName
	
	variable i,j
	
	for(i=0;i<dimsize(inputwave,0);i+=1)
	for(j=0;j<dimsize(inputwave,1);j+=1)
		EDC[]=InputWave[i][j][p]
		wavestats /q EDC
		If(v_sum>cutoffintensity)
			EfIndex[i][j] = smoothderivoutput_ignore(EDC,SmoothRange,SmoothTimes,EfMinPt,EfMaxPt) // TempEfValue
			// EfIndex[i][j]=TempEfValue[0]
		Else
			EfIndex[i][j]=NaN
		EndIf
	endfor
	endfor

	
	variable fermiMax=round(WaveMax(EfIndex))
	variable fermiMin=round(WaveMin(EfIndex))
	Variable extra = fermiMax-fermiMin
	
	String photDepShiftStr = NameofWave(OutputWave)
	Make /O/N=(DimSize(Inputwave,0),Dimsize(inputwave,1),DimSize(Inputwave,2)+extra) $photDepShiftStr
       Wave photDepShift = $photDepShiftStr
       setscale /p x, dimoffset(inputwave,0),dimdelta(inputwave,0),photdepshift
       setscale /p y, dimoffset(inputwave,1),dimdelta(inputwave,1),photdepshift
       setscale /p z, -fermimax*dimdelta(inputwave,2),dimdelta(inputwave,2),photdepshift
       
       variable fermi, start
       
       for(i=0;i<dimsize(inputwave,0);i+=1)
       for(j=0;j<dimsize(inputwave,1);j+=1)
       	if(numtype(efindex[i][j])==2)
     			wavestats /q efindex
     			fermi=v_avg
            	else
       		fermi=efindex[i][j]
       	endif
       	start=round(fermimax-fermi)
       	if(start<0)
       		print "start is less than zero"
       	endif
       	if(start+dimsize(inputwave,2)-1>dimsize(photdepshift,2)-1)
       		print "troube in middle region"
       		print start
       		print dimsize(inputwave,2)
       		print dimsize(photdepshift,2)
       		print fermi
       		print fermimin
       		return(0)
       	endif
       	photDepShift[i][j][,start] = InputWave[i][j][0]
       	photDepShift[i][j][start,start+DimSize(Inputwave,2)-1] = InputWave[i][j][r-start]
      	photDepShift[i][j][start+DimSize(inputwave,2)-1,] = InputWave[i][j][DimSize(InputWave,2)-1]
       endfor
       endfor
		
End


// Ilya 20210120, naive way to try to get a better result from Tyler's macro, but running it several times over
//Function iterate_SmoothDerivOutput3D(InputWave,CutOffPercentage,SmoothRange, SmoothTimes,EfMinPt,EfMaxPt)
//
//	Wave InputWave
//	Variable CutOffPercentage
//	Variable SmoothRange
//	Variable SmoothTimes
//	Variable EfMinPt
//	Variable EfMaxPt
//
//	String OutputWave_s = GetWavesDataFolder(InputWave,2) + "_ef"
//	print OutputWave_s
//
//End
//
//	Make /O/N=1 $OutputWave_s
//	SmoothDerivOutput3D(InputWave,OutputWave,CutOffPercentage,SmoothRange, SmoothTimes,EfMinPt,EfMaxPt)
//
//End







// Old test version... can be deleted
//
//Function test_bin()
//
//	Wave source = root:NCSS1_5_125eV_t
//
//	MatrixOp /NTHR=0 /O root:NCSS1_5_125eV_t_collapse = sumBeams(source)
//	Wave collapse = root:NCSS1_5_125eV_t_collapse
//
//	MatrixOp /NTHR=0 /O root:NCSS1_5_125eV_t_sum = sum(collapse)
//	Wave sum_v = root:NCSS1_5_125eV_t_sum
//	
//	// Create a copy of the data where we set EDCs below cutoffIntensity to zero
//	String zeroed_s = GetWavesDataFolder(source,2) + "_zeroed"
//	Duplicate /O source, $zeroed_s
//	Wave zeroed = $zeroed_s
//
//	Variable cutoffPercent = 0.1 // cutoff number of counts
//	Variable cutoffIntensity = cutoffPercent*sum_v[0]/(DimSize(source,0)*DimSize(source,1))
//
//	// Variable time_id
//	// Variable time_out
//	// time_id = StartMSTimer
//	
//	// Keep the EDC or not
//	collapse[][] = (collapse[p][q] > cutoffIntensity) ? 1 : 0
//	MultiThread zeroed[][][] = (collapse[p][q]) ? zeroed[p][q][r] : 0 
//
//	// time_out = StopMSTimer(time_id)
//	// print time_out/10e6
//
//	Variable binning_factorX = 5, binning_factorY = 10
//	Make /O/N=(ceil(DimSize(source,0)/binning_factorX),ceil(DimSize(source,1)/binning_factorY),DimSize(source,2)) root:NCSS1_5_125eV_t_bin
//	Wave bin = root:NCSS1_5_125eV_t_bin
//	SetScale /P x, DimOffset(source,0)+(DimDelta(source,0)*(binning_factorX-1)/2), DimDelta(source,0)*binning_factorX, bin
//	SetScale /P y, DimOffset(source,1)+(DimDelta(source,1)*(binning_factorY-1)/2), DimDelta(source,1)*binning_factorY, bin
//	
//
//	String bin_EDC_s
//
//	String ef_s = GetWavesDataFolder(source,2) + "_efTrackBin"
//	Make /O/N=(ceil(DimSize(root:NCSS1_5_125eV_t,0)/binning_factorX),ceil(DimSize(root:NCSS1_5_125eV_t,1)/binning_factorY)) $ef_s
//	Wave ef_w = $ef_s
//	SetScale /P x, DimOffset(source,0)+(DimDelta(source,0)*(binning_factorX-1)/2), DimDelta(source,0)*binning_factorX, ef_w
//	SetScale /P y, DimOffset(source,1)+(DimDelta(source,1)*(binning_factorY-1)/2), DimDelta(source,1)*binning_factorY, ef_w
//
//	MatrixOp /NTHR=0 /O root:NCSS1_5_125eV_t_bin_sum = sumBeams(bin)
//	Wave sum_bin = root:NCSS1_5_125eV_t_bin_sum
//	
//	// Fit the EDC or not? (for border regions which have no overlap with data--strictly zero)
//
//	String fit_message
//	Variable i,j
//
//	for (i = 0; i < DimSize(bin,0); i += 1)
////		for (j = 0; j < 1; j += 1)
//		for (j = 0; j < DimSize(bin,1); j += 1)
//
//			if (i == 21 && j == 4)
//
//			bin_EDC_s = avgDistCurve3D(GetWavesDataFolder(zeroed,2),2,binning_factorX*i,binning_factorX*i+binning_factorX-1,binning_factorY*j,binning_factorY*j+binning_factorY-1)	
//			Wave bin_EDC = $bin_EDC_s
//			bin[i][j][] = bin_EDC[r]
//
//			if (sum_bin[i][j])
//
//				//	eF_w[i][j] = SmoothDerivOutput(bin_EDC,3,5,350,450)
//
//				fit_message = IB_fit_Fermi_edge(bin_EDC)
//				
//				if (!cmpstr(fit_message,"error"))
//					print i, j
//				endif
//				
//				Wave fitEFcoeff = $(bin_EDC_s + "_coeffF")		
//			
//				if (fitEFcoeff[0] > DimOffset(bin_EDC,0) && fitEFcoeff[0] < DimOffset(bin_EDC,0) + DimSize(bin_EDC,0)*DimDelta(bin_EDC,0))
//								
//					eF_w[i][j] = fitEFcoeff[0]
//	
//				else
//					
//					eF_w[i][j] = NaN
//	
//				endif
//						
//			//	print fitEFcoeff[0]
//		
//			else
//			
//				eF_w[i][j] = NaN
//			
//			endif
//			
//			endif
//			
//		endfor
//	endfor
//
//End
