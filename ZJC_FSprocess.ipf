//#pragma TextEncoding = "UTF-8"
//#pragma rtGlobals=3		// Use modern global access method and strict wave access.
//
//// 20201023 Fermi surface, preliminary analysis
//Function process_fs0_cut0()
//
//	DFREF saveDFR = GetDataFolderDFR()
//	SetDataFolder root:
//
//	Wave fs = root:fs0_raw_b // source
//
//	String cut_s = "root:fs0_cut0"
//	Make /O/N=(DimSize(fs,1),DimSize(fs,2)) $cut_s
//	Wave cut = $cut_s
//
//	cut[][] = (fs[380][p][q] + fs[381][p][q] + fs[382][p][q] + fs[383][p][q] + fs[384][p][q])/5
//	
//	SetScale /P x, DimOffset(fs,1), DimDelta(fs,1), cut
//	SetScale /P y, DimOffset(fs,2), DimDelta(fs,2), cut
//
//	cut[][] = (numtype(cut[p][q]) == 2) ? 0 : cut[p][q]
//
//	// Horizontal normalization
//	cut_s = IB_norm2D(cut_s,0,startNorm=20,stopNorm=380,fixAvg=1,verbose=0)	
//	Wave cut = $cut_s
//
//	// Vertical normalization
//	cut_s = IB_norm2D(cut_s,1,startNorm=0,stopNorm=60,fixAvg=1,verbose=0)	
//	Wave cut = $cut_s
//
//	// Gaussian blur subtraction (i.e. sharpening)
//	IB_export_wave_hd5("Volumes:Universe:20201021_ALS_BL4_PdTe_Fe3Sn5:GMS_fs0_cut0.h5",cut)
//
//	// ->->- Mathematica magic -<-<-
//
//	// Import from Mathematica
//	cut_s = IB_import_mathematica_wave("Volumes:Universe:20201021_ALS_BL4_PdTe_Fe3Sn5:GMS_fs0_cut0_m.h5",cut)
//	Wave cut = $cut_s
//
//	// Next interpolate, to make image approximately square
//	cut_s = cut_s + "_i"
//	Duplicate /O cut, $cut_s
//	Wave cut = $cut_s
//	Resample /UP=2/DIM=1 cut
//
//	print cut_s
//	SetDataFolder saveDFR
//
//End



Function ZJC_FSprocess(inputwave,polar_dimension, tilt_dimension, polar_center,tilt_center,photonenergy,workfunction)

//This is for angle -> momemtem transfer of fermi surface maps


wave inputwave
variable polar_dimension, tilt_dimension, polar_center,tilt_center,photonenergy
variable workfunction

if(polar_dimension == 0)
	Setscale /P x, 0.5123*Pi/180*sqrt(photonenergy - workfunction)*(dimOffset(inputwave,0)- polar_center),0.5138*Pi/180*sqrt(photonenergy - workfunction)*dimdelta(inputwave,0), inputwave
elseif(polar_dimension == 1)
	Setscale /P y, 0.5123*Pi/180*sqrt(photonenergy - workfunction)*(dimOffset(inputwave,1)- polar_center),0.5138*Pi/180*sqrt(photonenergy - workfunction)*dimdelta(inputwave,1), inputwave
elseif(polar_dimension == 2)
	Setscale /P z, 0.5123*Pi/180*sqrt(photonenergy - workfunction)*(dimOffset(inputwave,2)- polar_center),0.5138*Pi/180*sqrt(photonenergy - workfunction)*dimdelta(inputwave,2), inputwave
else
	print("wrong polar dimension")
	
endif

if(tilt_dimension == 0)
	Setscale /P x, 0.5123*Pi/180*sqrt(photonenergy - workfunction)*(dimOffset(inputwave,0)- tilt_center),0.5138*Pi/180*sqrt(photonenergy - workfunction)*dimdelta(inputwave,0), inputwave
elseif(tilt_dimension == 1)
	Setscale /P y, 0.5123*Pi/180*sqrt(photonenergy - workfunction)*(dimOffset(inputwave,1)- tilt_center),0.5138*Pi/180*sqrt(photonenergy - workfunction)*dimdelta(inputwave,1), inputwave
elseif(tilt_dimension == 2)
	Setscale /P z, 0.5123*Pi/180*sqrt(photonenergy - workfunction)*(dimOffset(inputwave,2)- tilt_center),0.5138*Pi/180*sqrt(photonenergy - workfunction)*dimdelta(inputwave,2), inputwave
else
	print("wrong polar dimension")
	
endif

variable energy_dimension


for(energy_dimension=0;energy_dimension<3;energy_dimension++)
	if((energy_dimension!=tilt_dimension)&&(energy_dimension!=polar_dimension))
		if(energy_dimension == 0)
			Setscale /P x, dimoffset(inputwave,0)-(photonenergy-workfunction),dimdelta(inputwave,0), inputwave
		elseif(energy_dimension == 1)
			Setscale /P y, dimoffset(inputwave,1)-(photonenergy-workfunction),dimdelta(inputwave,1), inputwave
		elseif(energy_dimension == 2)
			Setscale /P z, dimoffset(inputwave,2)-(photonenergy-workfunction),dimdelta(inputwave,2), inputwave
	
endif
		
		
		
		
		
	endif


endfor


end


Function ZJC_cutprocess(inputwave,tilt_dimension, tilt_center, photonenergy,workfunction)

//This is for angle -> momemtem transfer of cuts


wave inputwave
variable tilt_dimension,tilt_center, photonenergy,workfunction

if(tilt_dimension == 0)
	Setscale /P x, 0.5123*Pi/180*sqrt(photonenergy - workfunction)*(dimOffset(inputwave,0)- tilt_center),0.5138*Pi/180*sqrt(photonenergy - workfunction)*dimdelta(inputwave,0), inputwave
elseif(tilt_dimension == 1)
	Setscale /P y, 0.5123*Pi/180*sqrt(photonenergy - workfunction)*(dimOffset(inputwave,1)- tilt_center),0.5138*Pi/180*sqrt(photonenergy - workfunction)*dimdelta(inputwave,1), inputwave
else
	print("wrong polar dimension")
	
endif



if(tilt_dimension == 1)
	Setscale /P x, dimoffset(inputwave,0)-(photonenergy-workfunction),dimdelta(inputwave,0), inputwave
else
	Setscale /P y, dimoffset(inputwave,1)-(photonenergy-workfunction),dimdelta(inputwave,1), inputwave

endif
		


end