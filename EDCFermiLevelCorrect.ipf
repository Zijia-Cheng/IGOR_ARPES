#pragma rtGlobals=3		// Use modern global access method and strict wave access.

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

Function transposewave(inputwave,outputwave)

//Inputwave is a wave, and outputwave is the wavename
//This is for transposing the SSRL fermi surface map for further FS alignment

wave inputwave
string outputwave

make /O $outputwave
wave temp = $outputwave
matrixop /O temp = transposeVol(inputwave,3)
setscale /p x, dimoffset(inputwave,2),dimdelta(inputwave,2),temp
setscale /p y,dimoffset(inputwave,1),dimdelta(inputwave,1),temp
setscale /P z,dimoffset(inputwave,0),dimdelta(inputwave,0),temp



End


Function SmoothDerivOutput3D(InputWave,OutputWave,CutOffPercentage,SmoothRange, SmoothTimes,EfMinPt,EfMaxPt)
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
			smoothderivoutput(EDC,TempefValue,SmoothRange,SmoothTimes,EfMinPt,EfMaxPt)
			EfIndex[i][j]=TempEfValue[0]
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
	FuncFitMD/NTHR=0 TwoDQaud W_coef EfIndex /D=$FitEfIndexName
	
	variable fermiMax=round(WaveMax(FitEfIndex))
	variable fermiMin=round(WaveMin(FitEfIndex))
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
       	if(numtype(fitefindex[i][j])==2)
     			wavestats /q fitefindex
     			fermi=v_avg
            	else
       		fermi=fitefindex[i][j]
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

Function SmoothDerivOutput3D_nonfit(InputWave,OutputWave,CutOffPercentage,SmoothRange, SmoothTimes,EfMinPt,EfMaxPt)
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
			smoothderivoutput(EDC,TempefValue,SmoothRange,SmoothTimes,EfMinPt,EfMaxPt)
			EfIndex[i][j]=TempEfValue[0]
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








//dim0 is analyzer angle
//dim1 is binding energy
Function SmoothDerivOutput2D(InputCut,OutputCut,CutOffPercentage,SmoothRange,SmoothTimes,EfMinPt,EfMaxPt)
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
			SmoothDerivOutput(EDC, TempEfValue,SmoothRange,SmoothTimes,EfMinPt,EfMaxPt)
			EfIndex[i]=TempEfValue[0]
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

Function SmoothDerivOutput(InputEDC,Result,SmoothRange,SmoothTimes,EfMinPt,EfMaxPt)

	wave inputEDC
	wave result //This will be the kinetic energy that corresponds to minimum derivative.
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
	
	Result[0] = v_minloc+EfMinPt
	
End

//Function TwoDQaud(w,x,y) : FitFunc
//	Wave w
//	Variable x
//	Variable y
//
//	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
//	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
//	//CurveFitDialog/ Equation:
//	//CurveFitDialog/ f(x,y) = C1+C2*y+C3*y^2+C4*x+C5*x*y+C6*x*y^2+C7*x^2+C8*x^2*y+C9*x^2*y^2
//	//CurveFitDialog/ End of Equation
//	//CurveFitDialog/ Independent Variables 2
//	//CurveFitDialog/ x
//	//CurveFitDialog/ y
//	//CurveFitDialog/ Coefficients 9
//	//CurveFitDialog/ w[0] = C1
//	//CurveFitDialog/ w[1] = C2
//	//CurveFitDialog/ w[2] = C3
//	//CurveFitDialog/ w[3] = C4
//	//CurveFitDialog/ w[4] = C5
//	//CurveFitDialog/ w[5] = C6
//	//CurveFitDialog/ w[6] = C7
//	//CurveFitDialog/ w[7] = C8
//	//CurveFitDialog/ w[8] = C9
//
//	return w[0]+w[1]*y+w[2]*y^2+w[3]*x+w[4]*x*y+w[5]*x*y^2+w[6]*x^2+w[7]*x^2*y+w[8]*x^2*y^2
//End
