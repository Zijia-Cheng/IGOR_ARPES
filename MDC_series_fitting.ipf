#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include <Multi-peak Fitting 2.0>



Function MDC_series_fitting(inputwave,start_energy,end_energy,bin,start_momen,end_momen,peak_number)



//used for fitting a series of mdcs. Line shape analysis and kink analysis
//inputwave:full data path to 2D wave, with dim0 to be binding energy and dim1 to be momentum
//start_energy,end_energy: energy range for the fitting: in eV. Always start with energy which has the best signal/noise
//bin: index number for the binning, must be integer
//start_momen, end_momen: momentum range for the fitting: in A^-1


wave inputwave
variable start_energy,end_energy,bin
variable start_momen,end_momen
variable peak_number

//Fitting parameters: can adapt
variable minfrac = 0.4//minAutoFindFraction
variable smfact = 10
variable noiseEst = 2


//Creat a new data folder for storing the mdcs
killDataFolder /Z root:$(nameofwave(inputwave)+"_mdc_analysis")
NewDataFolder /O /S root:$(nameofwave(inputwave)+"_mdc_analysis")
string data_folder_name = (nameofwave(inputwave)+"_mdc_analysis")


//Creating mdcs

variable start_energy_index = scaleToIndex(inputwave,start_energy,0)
variable end_energy_index = scaleToIndex(inputwave,end_energy,0)
variable start_momen_index = scaleToIndex(inputwave,start_momen,1)
variable end_momen_index = scaleToIndex(inputwave,end_momen,1)
variable energy_step
variable step_sign

if(start_energy_index <= end_energy_index)
	energy_step = dimdelta(inputwave,0)
	step_sign = 1
else
	energy_step = -dimdelta(inputwave,0)
	step_sign = -1
endif

variable i = 0


//wave fitting_results is the wave for storing the fittting results
variable steps = floor(abs(start_energy_index - end_energy_index)/bin)+1
Make /O /N = (steps,2*peak_number) fitting_results
if(mod(bin,2)==1)
	Setscale /P x,start_energy,energy_step*bin,fitting_results
else
	Setscale /P x,start_energy-energy_step*0.5,energy_step*bin,fitting_results
endif

string mdc_name_list = ""
variable j 
for(i=0;i<steps;i++)
	
	make /O /N = (abs(start_momen_index-end_momen_index)+1) $("mdc"+num2str(i))=0
	wave mdc_temp = $("mdc"+num2str(i)) 
	Setscale /I x, start_momen,end_momen, mdc_temp
	mdc_name_list += "mdc"+num2str(i)+";"
	//Here we ignore the problem of boundary condition, assuming the original data size is much larger than the region of interest
	if(mod(bin,2)==1)
		for(j=-(bin-1)/2;j<=(bin-1)/2;j++)
			mdc_temp += inputwave[start_energy_index+i*bin*step_sign+j](x)  
		endfor
	else
		for(j =-(bin)/2;j<bin/2;j++)
			mdc_temp += inputwave[start_energy_index+i*bin*step_sign+j](x)
		endfor
		
	endif
	
	

endfor

//Fitting
print(MPF2_AutoMPFit("fit_results", "lorentzian", "PeakCoefs%d", "Linear", "baseline", mdc_name_list,"", 4, smFact=smfact,minAutoFindFraction = minfrac,noiseEst = noiseEst))



for(i=0;i<steps;i++)
	
	for(j=0;j<peak_number;j++)
		Wave/SDFR=$("fit_results_"+num2str(i)) temp = $("Peakcoefs"+num2str(j))
		fitting_results[i][2*j] = temp[0]//Peak position
		fitting_results[i][2*j+1] = temp[1]//peak width

	endfor
endfor






End
