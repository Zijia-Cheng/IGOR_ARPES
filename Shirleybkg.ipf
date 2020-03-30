#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.



Function shirleybkg(dataset [,A0,threshold,itlimit])
// this is for removing standard shirley background 
//dataset should be a 1D wave which lower binding energy correspond to lower index
//A0:initial value for the scatting amplitude
//threshold: judgement for the itineration.
//itlimit: max itineration numbers

// more information see http://www.qro.cinvestav.mx/~aherrera/reportesInternos/peakShirley.pdf
//By Zijia Cheng, zijiac@princeton.edu
//It is better to use normalized wave. 
//background is in a wave named bg

wave dataset
variable A0,threshold,itlimit

If(ParamIsDefault(A0))
	A0 = 0.001
endif

If(ParamIsDefault(threshold))
	threshold = 0.0005
endif

If(ParamIsDefault(itlimit))
	itlimit = 400
endif


//initialial parameters

Make /O /D /N = (DimSize(dataset,0)) bg = Mean(dataset,pnt2x(dataset,0),pnt2x(dataset,4))
wave bg
Duplicate /O bg, bgold, temp
wave bgold, temp
Setscale /P x, dimoffset(dataset,0),dimdelta(dataset,0),bg,bgold,temp
bgold[] = bg[p] + 5
Variable A = A0
variable ctr=0 // loops number
Variable Aold
Variable indi

Do
	temp[] = bg[p] - bgold[p]
	indi = 2*waveabsmax(temp)/(waveabsmean(bg)+waveabsmean(bgold))
	if((indi)<= threshold || ctr > itlimit )
		break
	endif
	bgold[] = bg[p]
	Aold = A
	Variable i = 0
	temp[] = dataset[p] - bg[p]
	For(i=1;i<DimSize(dataset,0);i++)//fix i =0
	  bg[i] = bg[0]+A*Sum(temp,pnt2x(dataset,0),pnt2x(dataset,i))
	  temp[i] = dataset[i] - bg[i]
	endfor
	A = Aold*(1+(dataset[Dimsize(dataset,0)-1]-bg[Dimsize(dataset,0)-1])/dataset[Dimsize(dataset,0)-1])
	ctr = ctr+1
	
While(1)

If(ctr > itlimit)
	Print("Iteration limit exceed ")
endif
killwaves temp, bgold


end

Function waveabsmean(wave1)
	wave wave1
	variable i,abssum=0
	For(i =0;i<Dimsize(wave1,0);i++)
		abssum = abssum +abs(wave1[i])
	endfor	
	return abssum/(Dimsize(wave1,0))
End

Function waveabsmax(wave1)
	wave wave1
	
	return Max(abs(Wavemax(wave1)),abs(Wavemin(wave1)))
End



Function seriesprocess(dataset)
//goal is to do background substraction on all temperature datas
//

wave dataset
Duplicate /O dataset, $(nameofWave(dataset)+"_nm")
wave temp2 = $(nameofWave(dataset)+"_nm")
wave temp1 = $(normDC2((nameofWave(dataset)),1))
temp2[][] = temp1[dimsize(temp1,0)-1-p][q]

Duplicate /O temp2, $(nameofWave(dataset)+"_bkgrm")
wave dataset_bkgrm = $(nameofWave(dataset)+"_bkgrm")
//Setscale /P x, indextoscale(temp2,0,dimsize(temp2,0)-1),-dimdelta(temp2,0),dataset_bkgrm

Variable i
Make /O /N = (Dimsize(temp1,0)), temp3=0
Setscale /P x, Dimoffset(temp1,0), DimDelta(temp1,0), temp3
For(i=0;i<DimSize(temp1,1);i++)
	//Make /O /N = (Dimsize(temp1,0)), $(nameofwave(dataset)+"_"+num2str(i)+"_")
	//wave temp3 = $(nameofwave(dataset)+"_"+num2str(i)+"_")
	//Setscale /P x, Dimoffset(temp1,0), DimDelta(temp1,0), temp3
	temp3[] = temp2[p][i]
	Shirleybkg(temp3)
	wave bg
	dataset_bkgrm[][i] = temp2[p][i] - bg[p]
	Dowindow /K $(nameofwave(dataset)+"_"+num2str(i)+"_0")
	Display /N = $(nameofwave(dataset)+"_"+num2str(i)+"_") dataset_bkgrm[][i]
	MoveWindow /W = $(nameofwave(dataset)+"_"+num2str(i)+"_0") /I 6.5*(Mod(i,5)),3.5*trunc(i/5),-1,-1

Endfor


End

