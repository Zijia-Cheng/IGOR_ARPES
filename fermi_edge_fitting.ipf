#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.



Function fermi_edge(pw,yw,xw): FitFunc
// For fitting with a fermi-dirac distribution convolued with a gaussian.
//Need to change variable upboard if the fitted wave is too broad.
//The EDC wave should have lower kinatic energy with lower index

wave pw,yw,xw
variable upbord = 100 //upper broading interms of indexes
Make /O /D /FREE /N = (dimsize(yw,0)+upbord) expwave,ywtemp =0 // Create free waves so there is no name confliction
setscale /P x,dimOffset(yw,0)-upbord*deltax(yw),deltax(yw),ywtemp
ywtemp = (pw[0])/(exp((x-pw[1])/pw[2])+1)
expwave[] = exp(-(p*deltax(yw)/pw[3])^2/2)

//normalize
Variable sumexp
sumexp = 2*sum(expwave)
expwave/= sumexp

//We perform convolution between ywtemp and expwave
yw = 0
Variable i,j
For(i=0;i<dimsize(ywtemp,0);i++)
	for(j=0;j<dimsize(ywtemp,0);j++)
		if(i >=upbord)
		 yw[i-upbord] += ywtemp[j]*expwave[abs(i-j)] 
		endif
	endfor
endfor // convolution
yw += pw[4]


End






Function fermi_edge_fitting(spectrum,temperature)
//This is for fitting fermi edge as much automatically as possibly.
//the wave spectrum should be a 1D wave and must be scaled in eV. Lower kinetic energy should have lower index
//variable temperature is in Kelvin
//Zijia


wave spectrum
variable temperature //in K
Make /D /O /N = 4 fitt_pars = 0
Duplicate /O spectrum,$(nameofwave(spectrum)+"_")
wave temp = $(nameofwave(spectrum)+"_")

smooth 5, temp //smooth the data to avoid sudden jump...
Differentiate temp /D = $(nameofwave(spectrum)+"_dif")
wave temp2 = $(nameofwave(spectrum)+"_dif")
wavestats temp2
fitt_pars[0] = V_minloc // x location of the minimum.
fitt_pars[3] = Mean(spectrum,pnt2x(spectrum,dimsize(spectrum,0)-6),pnt2x(spectrum,dimsize(spectrum,0)-1))


//for estimation of half width and height, just assume the fitting function is precise at guessed Ef 

variable yvalue = spectrum(V_minloc) - fitt_pars[3]
fitt_pars[2] = 2*yvalue+fitt_pars[3]

Findvalue /RMD = [x2pnt(spectrum,fitt_pars[0]),*] /T = (abs(temp(V_minloc) - temp(V_minloc+2*deltax(temp)))) /V =(yvalue*2/(1+exp(1))+fitt_pars[3]) spectrum
fitt_pars[1] = abs(pnt2x(spectrum,V_value) - fitt_pars[0])

string name= nameofwave(spectrum)+"dis" //name of the displaying windows
killwindow /Z $name
display /N = $name spectrum
Dowindow /F $name
ModifyGraph width=360,height=144
ModifyGraph /W = $name mode($nameofwave(spectrum)) = 3, marker($nameofwave(spectrum)) = 41

Duplicate /O spectrum, $(nameofwave(spectrum)+"_fit") 
//if you want to choose let constant term be zero, especially when the choosed fitting range doesn't extend to well above fermi level
fitt_pars[3] = 0
Funcfit/TBOX = 768 /H="0001" ZC_Fer,fitt_pars,spectrum /D = $(nameofwave(spectrum)+"_fit")


//Funcfit/TBOX = 768 ZC_Fer,fitt_pars,spectrum /D = $(nameofwave(spectrum)+"_fit")

Appendtograph $(nameofwave(spectrum)+"_fit")

variable kb = 8.617333*10^-5
variable temp3 = kb*temperature/(4*fitt_pars[1])
if(temp3 >=0.25)
	print"Gaussian fitting is wrong,width is even smaller than thermal broadening"
else
	wave /Z temp4 = fermidetable //requires the existence of this wave, which can be generated from function generatedetable()
if(waveexists(temp4) ==0)
	generatedetable()
	wave temp4 = fermidetable
endif
FindValue /T = 0.001 /V=(temp3) temp4 //tolerence should be even smaller of temp3 is appoaching 0.25

//fitting with fermi_edge(pw,yw,xw): FitFunc

Make /O /D /N = 5 fitt_pars2 = 0
fitt_pars2[0] = fitt_pars[2]
fitt_pars2[1] = fitt_pars[0]
fitt_pars2[2] = kb*temperature // fix thermal broadening
fitt_pars2[3] = fitt_pars2[2]/pnt2x(temp4,V_value)
fitt_pars2[4] = fitt_pars[3]

Duplicate /O spectrum,$(nameofwave(spectrum)+"_fit2")
Funcfit/TBOX = 768 /H="00101" fermi_edge,fitt_pars2,spectrum /D = $(nameofwave(spectrum)+"_fit2")
Appendtograph $(nameofwave(spectrum)+"_fit2")

endif




End


Function fermi_edge_fitting_series(spectrum,Twave)
//goal: fitting a series of EDCs' fermi edge in order to get fermi energy and broadening.
//spectrum should be a two dimensional wave with energy be dimension 0.
//Twave is a 1D wave containing temperatures of EDCs.
//Zijia

wave spectrum
wave Twave //1D: dimsize = dimsize(spectrum,1)
Make /D /O /N =(dimsize(spectrum,1),5) $(nameofwave(spectrum)+"_fitpar")
wave temp = $(nameofwave(spectrum)+"_fitpar")

Variable i
For(i = 0;i<dimsize(spectrum,1);i++)
	Make /O /N = (dimsize(spectrum,0)) $(nameofwave(spectrum)+"_"+num2str(i))
	wave temp1 = $(nameofwave(spectrum)+"_"+num2str(i))
	Setscale /P x, dimoffset(spectrum,0),dimdelta(spectrum,0),temp1
	temp1 = spectrum(x)[i]
	fermi_edge_fitting(temp1,Twave[i])
	wave fitt_pars2
	temp[i][] = fitt_pars2[q]
	MoveWindow /W = $(nameofwave(temp1)+"dis") /I 6.5*(Mod(i,5)),3.5*trunc(i/5),-1,-1
	
endfor 

Make /O /N =(dimsize(spectrum,1)) $(nameofwave(spectrum)+"_resolution")
wave temp1 = $(nameofwave(spectrum)+"_resolution")
temp1[] = 2.35*temp[p][3]
display temp[][3] vs Twave[]


killwaves fitt_pars2

end

Function ZC_Fer(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = w_2/(exp((x-w_0)/w_1) + 1)+w_3
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 4
	//CurveFitDialog/ w[0] = w_0
	//CurveFitDialog/ w[1] = w_1
	//CurveFitDialog/ w[2] = w_2
	//curvefitdialog/ w[3] = w_3

	return w[2]/(exp((x-w[0])/w[1]) + 1)+w[3]
End


Function generatedetable()
//generate a dictionary

Make /O /D /N = ((5-0.01)/(0.001)+1) fermidetable = 0
Setscale /P x,0.01,0.001,fermidetable

variable i,j,temp
for(i=0;i<dimsize(fermidetable,0);i++)
	for(j = -200;j<=200;j+=0.01)
		temp = pnt2x(fermidetable,i)
		fermidetable[i]+= temp*exp(-(j*temp)^2/2)/(100*sqrt(2*pi)*(exp(-j)+2+exp(j)))
	endfor

endfor

end