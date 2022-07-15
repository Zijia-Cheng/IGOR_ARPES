#pragma TextEncoding = "MacRoman"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.


Function HDF5LoadDataset(datasetName, groupname, filename, desname)

	// by Zijia Cheng

	String datasetName	// Name of dataset to be loaded
	String filename     // Name of the filename where the dataset is
	String groupname    // Name of the groupname where the dataset is
	String desname      // Name of the loaded wave
	
	Variable fileID	// HDF5 file ID will be stored here
   Variable groupID //group ID will be stored here
   
	Variable result = 0	// 0 means no error
	
	// Open the HDF5 file.
	HDF5OpenFile  /R /Z fileID as filename
	if (V_flag != 0)
		Print "HDF5OpenFile failed"
		return -1
	endif
	
	//open the group
	HDF5OpenGroup /Z fileID, groupname, groupID
	if (V_flag != 0)
		Print "HDF5Oloadgroup failed"
		return -1
	endif
	
	
	// Load the HDF5 dataset.
	HDF5LoadData /O /Z /N=$desname groupID, datasetName
	if (V_flag != 0)
		Print "HDF5LoadData failed"
		result = -1
	endif

	// Close the HDF5 group
	HDF5CloseGroup groupID
	
	// Close the HDF5 file.
	HDF5CloseFile fileID

	return result
End



Function LoadHDF5NumericAttribute(pathName, filePath, groupPath, objectName, objectType, attributeName)
	
	//return output
	
	String pathName	// Symbolic path name - can be "" if filePath is a full path
	String filePath	// file name or partial path relative to symbolic path, or full path to file
	String groupPath	// Path to group, such as "/", "/my_group"
	String objectName	// Name of group or dataset
	Variable objectType	// 1=group, 2=dataset
	String attributeName	// Name of attribute
	

	Variable output // output
	Variable result = 0

	// Open the HDF5 file
	Variable fileID	// HDF5 file ID will be stored here
	HDF5OpenFile /P=$pathName /R /Z fileID as filePath
	if (V_flag != 0)
		Print "HDF5OpenFile failed"
		return -1
	endif

	Variable groupID	// HDF5 group ID will be stored here
	HDF5OpenGroup /Z fileID, groupPath, groupID
	if (V_flag != 0)
		Print "HDF5OpenGroup failed"
		HDF5CloseFile fileID
		return -1
	endif
	
	HDF5LoadData /O /A=attributeName /TYPE=(objectType) /N=tempAttributeWave /Z fileID, objectName
	result = V_flag	// 0 if OK or non-zero error code
	
	if (result == 0)
		Wave tempAttributeWave
		if (WaveType(tempAttributeWave) == 0)
			output = NaN	// Attribute is string, not numeric
			result = -1
		else
			output = tempAttributeWave[0]
		endif
		//KillWaves/Z tempAttributeWave
	endif

	// Close the HDF5 group
	HDF5CloseGroup groupID

	// Close the HDF5 file
	HDF5CloseFile fileID
	
	if(result !=0)
		print "error loading attributes" + objectName
	endif

	//print result
	//print "0: no error. Otherwise, error"
	
	return output
	
End


Function SSRL_HDFload_3D(startnumber, endnumber,presur,filepath)
	//For loading both photon energy dependence + maps 
	Variable startnumber, endnumber
	String presur //TMFS1_
	String filepath //"Users:zijiacheng:Desktop:Princeton:ARPES_projects:Co3Sn2S2:20200201:"
	Variable i = 0;
	String name
	
	For( i = Min(startnumber, endnumber); i<= Max(startnumber, endnumber);i++)
		
		
		if(i<10)
			name = presur+"000"+num2str(i)
		elseif(10<=i && i<100)
			name = presur+"00"+num2str(i)
		elseif(i>=100 && i<1000)
			name = presur+"0"+num2str(i)
		else
			print("number out of range")
		endif			
		
		String datapath = filepath +name+".h5"
		HDF5LoadDataset("Count","Data",datapath, name+"_count")
		HDF5LoadDataset("Time","Data",datapath, name+"_time")
		wave count = $(name+"_count")
		wave time0 = $(name+"_time")
		
		Duplicate /O count, $name, temp1
		wave temp = $name
		
		Matrixop /O time0_ = replaceNaNs(time0,0)
		Matrixop /O count_ = replaceNaNs(count,0)
		
		
		
		temp1[][][] = count_[p][q][r]/time0_[p][q][r]
		Matrixop /O temp = replaceNaNs(temp1,0)
		killwaves count, temp1,count_,time0_,time0
		
		Variable offset_energy = LoadHDF5NumericAttribute("",datapath,"/Data/Axes0","/Data/Axes0",1,"Offset")
		Variable delta_energy  = LoadHDF5NumericAttribute("",datapath,"/Data/Axes0","/Data/Axes0",1,"Delta")
		Setscale /P x, offset_energy, delta_energy, temp
		
		Variable offset_mome = LoadHDF5NumericAttribute("",datapath,"/Data/Axes1","/Data/Axes1",1,"Offset")
		Variable delta_mome  = LoadHDF5NumericAttribute("",datapath,"/Data/Axes1","/Data/Axes1",1,"Delta")
		Setscale /P y, offset_mome, delta_mome, temp
		
		Variable offset_polar = LoadHDF5NumericAttribute("",datapath,"/Data/Axes2","/Data/Axes2",1,"Offset")
		Variable delta_polar = LoadHDF5NumericAttribute("",datapath,"/Data/Axes2","/Data/Axes2",1,"Delta")
		Setscale /P z, offset_polar, delta_polar, temp
		
		
		
		

		
	Endfor



End


Function Ultra_HDFload_3D(startnumber, endnumber,presur,filepath)
	//For loading both photon energy dependence + maps for Ultra
	//20220217 by ZJC
	
	Variable startnumber, endnumber
	String presur //TMFS1_
	String filepath //"Users:zijiacheng:Desktop:Princeton:ARPES_projects:Co3Sn2S2:20200201:"
	Variable i = 0;
	String name
	
	For( i = Min(startnumber, endnumber); i<= Max(startnumber, endnumber);i++)
		
		
		if(i<10)
			name = presur+"000"+num2str(i)
		elseif(10<=i && i<100)
			name = presur+"00"+num2str(i)
		elseif(i>=100 && i<1000)
			name = presur+"0"+num2str(i)
		else
			print("number out of range")
		endif			
		
		String datapath = filepath +name+".h5"
		HDF5LoadDataset("Image Data","Electron Analyzer",datapath, name+"_data")
		wave count = $(name+"_data")
		Redimension/S count
		
		//Matrixop /O count_ = replaceNaNs(count,0)
		
		
		
		LoadHDF5NumericAttribute("",datapath,"/Electron Analyzer","/Electron Analyzer/Image Data",2,"Axis0.Scale")
		wave tempAttributeWave
		Setscale /P x, tempAttributeWave[0], tempAttributeWave[1], count
		
		
		LoadHDF5NumericAttribute("",datapath,"/Electron Analyzer","/Electron Analyzer/Image Data",2,"Axis1.Scale")
		Setscale /P y, tempAttributeWave[0], tempAttributeWave[1], count
		
		LoadHDF5NumericAttribute("",datapath,"/Electron Analyzer","/Electron Analyzer/Image Data",2,"Axis2.Scale")
		Setscale /P z, tempAttributeWave[0], tempAttributeWave[1], count
		
		
		
		

		
	Endfor



End


Function Ultra_HDFload_2D(startnumber, endnumber,presur,filepath)
	//For loading cuts for Ultra
	//20220217 by ZJC
	
	Variable startnumber, endnumber
	String presur //TMFS1_
	String filepath //"Users:zijiacheng:Desktop:Princeton:ARPES_projects:Co3Sn2S2:20200201:"
	Variable i = 0;
	String name
	
	For( i = Min(startnumber, endnumber); i<= Max(startnumber, endnumber);i++)
		
		
		if(i<10)
			name = presur+"000"+num2str(i)
		elseif(10<=i && i<100)
			name = presur+"00"+num2str(i)
		elseif(i>=100 && i<1000)
			name = presur+"0"+num2str(i)
		else
			print("number out of range")
		endif			
		
		String datapath = filepath +name+".h5"
		HDF5LoadDataset("Image Data","Electron Analyzer",datapath, name+"_data")
		wave count = $(name+"_data")
		Redimension/S count
		
		//Matrixop /O count_ = replaceNaNs(count,0)
		
		
		
		LoadHDF5NumericAttribute("",datapath,"/Electron Analyzer","/Electron Analyzer/Image Data",2,"Axis0.Scale")
		wave tempAttributeWave
		Setscale /P x, tempAttributeWave[0], tempAttributeWave[1], count
		
		
		LoadHDF5NumericAttribute("",datapath,"/Electron Analyzer","/Electron Analyzer/Image Data",2,"Axis1.Scale")
		Setscale /P y, tempAttributeWave[0], tempAttributeWave[1], count
		


		
	Endfor



End









Function Diamond_HDFload_3D(startnumber, endnumber,filepath)

	//ZJC 20210910
	//Load a series of 3D files with number index(in the name) from start number to endnumber
	Variable startnumber, endnumber // Sixdigits in the file name(usually: like 130119)
	String filepath //"Users:zijiacheng:Desktop:Princeton:ARPES_projects:Co3Sn2S2:20200201:"
	Variable i = 0;
	String name1,name2
	String presur="i05" 
	
	For( i = Min(startnumber, endnumber); i<= Max(startnumber, endnumber);i++)
		
		name1 = presur+"-"+num2istr(i)
		name2 = presur+"_"+num2istr(i)
		String datapath = filepath +name1+".nxs"
		HDF5LoadDataset("data","entry1/analyser",datapath, name2+"_data")
		wave count = $(name2+"_data")
		HDF5LoadDataset("angles","entry1/analyser",datapath, name2+"_angle")
		wave angle = $(name2+"_angle")	
		HDF5LoadDataset("energies","entry1/analyser",datapath, name2+"_energy")
		wave energy = $(name2+"_energy")
		HDF5LoadDataset("sapolar","entry1/analyser",datapath, name2+"_polar")
		wave polar = $(name2+"_polar")					
		
		Setscale /I y, angle[0], angle[dimsize(angle,0)-1], count
		Setscale /I z, energy[0][0], energy[0][dimsize(energy,1)-1], count		
		Setscale /I x, polar[0], polar[dimsize(polar,0)-1], count			
		
		Redimension/S count
		
		killwaves angle,energy,polar
		
		
		
		
				
//		Matrixop /O time0_ = replaceNaNs(time0,0)
//		Matrixop /O count_ = replaceNaNs(count,0)
		
		
//		temp1[][][] = count_[p][q][r]/time0_[p][q][r]
//		Matrixop /O temp = replaceNaNs(temp1,0)
//		killwaves count, temp1,count_,time0_,time0
		
//		Variable offset_energy = LoadHDF5NumericAttribute("",datapath,"/Data/Axes0","/Data/Axes0",1,"Offset")
//		Variable delta_energy  = LoadHDF5NumericAttribute("",datapath,"/Data/Axes0","/Data/Axes0",1,"Delta")
//		Setscale /P x, offset_energy, delta_energy, temp
//		
//		Variable offset_mome = LoadHDF5NumericAttribute("",datapath,"/Data/Axes1","/Data/Axes1",1,"Offset")
//		Variable delta_mome  = LoadHDF5NumericAttribute("",datapath,"/Data/Axes1","/Data/Axes1",1,"Delta")
//		Setscale /P y, offset_mome, delta_mome, temp
//		
//		Variable offset_polar = LoadHDF5NumericAttribute("",datapath,"/Data/Axes2","/Data/Axes2",1,"Offset")
//		Variable delta_polar = LoadHDF5NumericAttribute("",datapath,"/Data/Axes2","/Data/Axes2",1,"Delta")
//		Setscale /P z, offset_polar, delta_polar, temp
		
		
		
		

		
	Endfor



End



Function Diamond_HDFload_2D(startnumber, endnumber,filepath)

	//ZJC 20210910
	Variable startnumber, endnumber
	String filepath //"Users:zijiacheng:Desktop:Princeton:ARPES_projects:Co3Sn2S2:20200201:"
	Variable i = 0;
	String name1,name2
	String presur="i05" 
	
	
	For( i = Min(startnumber, endnumber); i<= Max(startnumber, endnumber);i++)
		
		name1 = presur+"-"+num2istr(i)
		name2 = presur+"_"+num2istr(i)
		String datapath = filepath +name1+".nxs"
		HDF5LoadDataset("data","entry1/analyser",datapath, name2+"_data")
		wave count = $(name2+"_data")
		HDF5LoadDataset("angles","entry1/analyser",datapath, name2+"_angle")
		wave angle = $(name2+"_angle")	
		HDF5LoadDataset("energies","entry1/analyser",datapath, name2+"_energy")
		wave energy = $(name2+"_energy")

		Setscale /I y, angle[0], angle[dimsize(angle,0)-1], count
		Setscale /I z, energy[0][0], energy[0][dimsize(energy,1)-1], count		
		
		
		Redimension/S count
		
		killwaves angle,energy
		
		
		
		
				
//		Matrixop /O time0_ = replaceNaNs(time0,0)
//		Matrixop /O count_ = replaceNaNs(count,0)
		
		
//		temp1[][][] = count_[p][q][r]/time0_[p][q][r]
//		Matrixop /O temp = replaceNaNs(temp1,0)
//		killwaves count, temp1,count_,time0_,time0
		
//		Variable offset_energy = LoadHDF5NumericAttribute("",datapath,"/Data/Axes0","/Data/Axes0",1,"Offset")
//		Variable delta_energy  = LoadHDF5NumericAttribute("",datapath,"/Data/Axes0","/Data/Axes0",1,"Delta")
//		Setscale /P x, offset_energy, delta_energy, temp
//		
//		Variable offset_mome = LoadHDF5NumericAttribute("",datapath,"/Data/Axes1","/Data/Axes1",1,"Offset")
//		Variable delta_mome  = LoadHDF5NumericAttribute("",datapath,"/Data/Axes1","/Data/Axes1",1,"Delta")
//		Setscale /P y, offset_mome, delta_mome, temp
//		
//		Variable offset_polar = LoadHDF5NumericAttribute("",datapath,"/Data/Axes2","/Data/Axes2",1,"Offset")
//		Variable delta_polar = LoadHDF5NumericAttribute("",datapath,"/Data/Axes2","/Data/Axes2",1,"Delta")
//		Setscale /P z, offset_polar, delta_polar, temp
		
		
		
		

		
	Endfor



End





Function SSRL_HDFload_2D(startnumber, endnumber,presur,filepath)
	
	Variable startnumber, endnumber
	string presur//TMFS1_
	string filepath//"Users:zijiacheng:Desktop:Princeton:ARPES_projects:Co3Sn2S2:20200201:"
	string name
	
	
	Variable i = 0;
	
	String /G groupname = ""
	
	
	For( i = Min(startnumber, endnumber); i<= Max(startnumber, endnumber);i++)
		
		if(i<10)
			name = presur+"000"+num2str(i)
		elseif(10<=i && i<100)
			name = presur+"00"+num2str(i)
		elseif(i>=100 && i<1000)
			name = presur+"0"+num2str(i)
		else
			print("number out of range")
		endif		
		
		groupname += name+";"
		String datapath = filepath +name+".h5"
		
		HDF5LoadDataset("Count","Data",datapath, name+"_count")
		HDF5LoadDataset("Time","Data",datapath, name+"_time")
		wave count = $(name+"_count")
		wave time0 = $(name+"_time")
		
		Duplicate /O count, $name, temp1
		wave temp = $name
		
		Matrixop /O time0_ = replaceNaNs(time0,0)
		Matrixop /O count_ = replaceNaNs(count,0)
		
		
		temp1[][] = count_[p][q]/time0_[p][q]
		Matrixop /O temp = replaceNaNs(temp1,0)
		
		Variable offset_energy = LoadHDF5NumericAttribute("",datapath,"/Data/Axes0","/Data/Axes0",1,"Offset")
		Variable delta_energy  = LoadHDF5NumericAttribute("",datapath,"/Data/Axes0","/Data/Axes0",1,"Delta")
		Setscale /P x, offset_energy, delta_energy, temp
		
		Variable offset_mome = LoadHDF5NumericAttribute("",datapath,"/Data/Axes1","/Data/Axes1",1,"Offset")
		Variable delta_mome  = LoadHDF5NumericAttribute("",datapath,"/Data/Axes1","/Data/Axes1",1,"Delta")
		Setscale /P y, offset_mome, delta_mome, temp
		
		killwaves count,  temp1,  count_
		
		
		
		

		
	Endfor



End

Function SSRL_HDFload_photonenergydep(number,presur,filepath,anglecenter)

//anglecenter is the analyser slit angle which gamma zero locates at
//number is the file number
//fermi level correction can be adjusted
//Need to update based on new load HDF5 file format


variable number
string presur
string filepath
variable anglecenter

string name

if(number<10)
	name = presur+"000"+num2str(number)
elseif(10<=number && number<100)
	name = presur+"00"+num2str(number)
elseif(number>=100 && number<1000)
	name = presur+"0"+num2str(number)
else
	print("number out of range")
endif		

String datapath = filepath +name+".h5"


HDF5LoadDataset("Count","Data",datapath, name+"_count")
HDF5LoadDataset("Time","Data",datapath, name+"_time")
wave count = $(name+"_count")
wave time0 = $(name+"_time")
		
Duplicate /O count, $name, temp1
wave temp = $name


Matrixop /O time0_ = replaceNaNs(time0,0)
Matrixop /O count_ = replaceNaNs(count,0)
		
temp1[][][] = count_[p][q][r]
Matrixop /O temp = replaceNaNs(temp1,0)
killwaves count, temp1, count_

variable i = 0

HDF5LoadDataset("Data:Axes0:Offset","MapInfo",datapath, name+"_energyoffset")
wave energyoffset = $(name+"_energyoffset")
//HDF5LoadDataset("Data:Axes1:Offset","Data",datapath, name+"_angleoffset")
//wave angle_offset = $(name+"_angleoffset")
//HDF5LoadDataset("Data:Axes1:Delta","Data",datapath, name+"_angledelta")
//wave angle_delta = $(name+"_angledelta")

Variable delta_energy  = LoadHDF5NumericAttribute("",datapath,"/Data/Axes0","/Data/Axes0",1,"Delta")
Variable startenergy = LoadHDF5NumericAttribute("",datapath,"/Data/Axes2","/Data/Axes2",1,"Offset")
Variable energy_stepsize = LoadHDF5NumericAttribute("",datapath,"/Data/Axes2","/Data/Axes2",1,"Delta")
Variable energy_steps = LoadHDF5NumericAttribute("",datapath,"/Data/Axes2","/Data/Axes2",1,"Count")
Variable angle_offset_c = LoadHDF5NumericAttribute("",datapath,"/Data/Axes1","/Data/Axes1",1,"Offset")
Variable angle_delta_c = LoadHDF5NumericAttribute("",datapath,"/Data/Axes1","/Data/Axes1",1,"Delta")



Make /T /N=(energy_steps) /O $(name+"_namelist")
Make /N =(energy_steps)/O $(name+"_energylist")
wave /T namelist=$(name+"_namelist")
wave energylist = $(name+"_energylist")

string wavenames
variable fermilevel//corection for Ef drift at different photon energy 
for(i=0;i<energy_steps;i++)
	energylist[i] = startenergy+i*energy_stepsize
	wavenames = name+"_"+num2str(energylist[i])+"eV"
	namelist[i] =wavenames
	Duplicate /O /RMD=[][][i] temp,$wavenames
	wave cut2D = $wavenames
	redimension /N=(-1,-1) cut2D
	Setscale /P	x, energyoffset[i],delta_energy,cut2D
	//Setscale /P y, angle_offset[i]-anglecenter,angle_delta[i],cut2D
	Setscale /P y, angle_offset_c-anglecenter,angle_delta_c,cut2D
	fermilevel = fixfermilevel_SSRL(cut2D,energylist[i]-4.365, 0.3)//0.4 can be changed
	Setscale /P	x, energyoffset[i] - fermilevel + energylist[i]-4.365,delta_energy,cut2D
	
	
	

endfor



End


Function fixfermilevel_SSRL(ekcut, startenergy, energyrange)
	//2this is a part of the program used for align the fermi level for photon energy dependence, etc.
	//ekcut should be a 2d wave with dimension 1: angle/momentum; dimension 0: energy.
	//Method: find the largest derivative and that is the fermi energy
	//The fermi level should within startenergy+- energyrange
	
	wave ekcut
	variable startenergy, energyrange
	Make /O /N = (dimsize(ekcut,0)), energyspec = 0
	Variable i,j, fermi
	
	
	// make energy spectrum
	for(i=0;i < dimsize(ekcut,0);i++)
		for(j = 0; j  < dimsize(ekcut,1);j++)
			energyspec[i] =  ekcut[i][j] + energyspec[i]
		endfor		
	endfor
	
	Setscale /P x , dimoffset(ekcut,0), dimdelta(ekcut,0), energyspec 
	if(startenergy+energyrange > dimoffset(energyspec,0)+(dimsize(energyspec,0)-1)*dimDelta(energyspec,0))
		Duplicate /O /R = (startenergy-energyrange,dimoffset(energyspec,0)+(dimsize(energyspec,0)-1)*dimDelta(energyspec,0)) energyspec, $"test_SSRL"
	else
		Duplicate /O /R = (startenergy-energyrange, startenergy+energyrange) energyspec, $"test_SSRL"
	endif
	//method: use differentiate to get the largest change point in the integerated spectrum:
	//seems we need to narrow down the energy range.
	
	Differentiate/METH=1 $"test_SSRL"/D=energyspec__DIF;DelayUpdate
	wave energyspec__DIF
	
	energyspec__DIF[] = abs(energyspec__DIF[p])
	
	Findvalue /V = (Wavemax(energyspec__DIF)) energyspec__DIF
	fermi = Indextoscale(energyspec__DIF,V_value,0)
	
	//reset scaling of ekcutprin
	//Setscale /P y, dimoffset(ekcut,1) - fermi, dimdelta(ekcut,1), ekcut
	killwaves $"test_SSRL"
	
	return fermi
	
	
End













Function combine()
	
	//use for combining cuts
	String /G namelist 
	variable i
	namelist = ""
	
	for(i = 187; i>153;i =i-2)
		String name = "CSS1_"+"00"+num2str(i)
		namelist = namelist + name+";"	
	endfor
	
	concatenate /O namelist, $"tempdep_weyl_cone_2"
	
		
End


Function extractDC3D(inputwave,index,position, bin, outputname)
// for extracting DCs from cuts etc. Inputwave is a three dimenional wave, position is the momentum or energy position(index, not scale), which should be dimension 1(0). 
//The returned wave will be a 2 dimenional wave which has spectrum in 0 dimenion and the other dimenion is the dimention 2 in origianl wave.
//index ==0: EDC; otherwise :MDC

Wave inputwave
Variable position
String outputname
Variable bin,index
variable i = 0

if(index ==0) // EDC

	Make /O /N = (dimsize(inputwave,0), dimsize(inputwave,2)) $outputname
	Setscale /P x, dimoffset(inputwave,0),dimdelta(inputwave,0), $outputname
	Setscale /P y, dimoffset(inputwave,2),dimdelta(inputwave,2), $outputname
	Wave temp = $outputname
	temp[][] = 0
	
	
	if(mod(bin,2) ==0)
		for(i = -bin/2; i<= bin/2; i++)
			temp[][] = temp[p][q]+ inputwave[p][position+i][q]
		
		endfor
		
	else
		for(i = -(bin-1)/2; i<= (bin-1)/2; i++)
			temp[][] = temp[p][q]+ inputwave[p][position+i][q]	
		endfor
	
	endif

else //MDC
	
	Make /O /N = (dimsize(inputwave,1), dimsize(inputwave,2)) $outputname
	Setscale /P x, dimoffset(inputwave,1),dimdelta(inputwave,1), $outputname
	Setscale /P y, dimoffset(inputwave,2),dimdelta(inputwave,2), $outputname
	Wave temp = $outputname
	temp[][] = 0
	
	
	if(mod(bin,2) ==0)
		for(i = -bin/2; i<= bin/2; i++)
			temp[][] = temp[p][q]+ inputwave[position+i][p][q]
		
		endfor
		
	else
		for(i = -(bin-1)/2; i<= (bin-1)/2; i++)
			temp[][] = temp[p][q]+ inputwave[position+i][p][q]
		endfor
	
	endif

endif

End



Function displayEDC(dataset)
// for displaying a list of EDCs. the energy should be dimension 0. dataset is a 2 dimension wave

wave dataset
Variable i
string windname = "tempdep_gap1_bin10"
Dowindow /K $(windname+"0")
For(i = 0;i<Dimsize(dataset,1);i++)
	if(i ==0)
		Display /N =$(windname) dataset[][i]
		Dowindow /F $(windname+num2str(0))
	else
		AppendtoGraph /W = $(windname+num2str(0)) dataset[][i]
		ModifyGraph /W =$(windname+num2str(0))  rgb($(nameofWave(dataset)+"#"+num2str(i))) = (3000*i,65535-3000*i,3000*i)
		Modifygraph /W = $(windname+num2str(0)) offset($(nameofWave(dataset)+"#"+num2str(i))) = {0,i*wavemax(dataset)/1.8}
		
		
	endif

endfor



End





Function dispseries(spectrum)
//display series of EDCs which are offset from each other 
wave spectrum
variable i=0


for(i=0;i<dimsize(spectrum,0);i++)
	if(i==0)
		Dowindow /K $(nameofwave(spectrum)+"_dis")
		display /N = $(nameofwave(spectrum)+"_dis") spectrum[0][]
	else
		appendtograph /W = $(nameofwave(spectrum)+"_dis") spectrum[i][]
		modifygraph offset($(nameofwave(spectrum)+"#"+num2str(i))) = {0,i*WaveMax(spectrum)/4}
	
	endif
	
endfor

end


Function Fermivelocity(dataset)
//this is for studying the fermivelocity and crossing energy dependence with temperature for CSS
//can be used for future general fitting applications
//dependence: IB_iterative_fit



wave dataset
variable index = 17
wave energy_range_weyl
string /G windowlist = ""
Make /O /N = 18 temperature = {20,40,60,80,100,120,130,140,150,160,170,180,190,210,230,250,270,290} // temperature wave
wave /T const1 //fitting constraint


wave fit_pars_1 //this is for storing initial guess for each temperature fitting
//for(index = 0;index < 18;index++)
	variable energystart = energy_range_weyl[0][index]
	variable energyend = energy_range_weyl[1][index] //fitting range for energy 
	Variable steps = 7
	variable stepsize = (energyend - energystart)/(steps-1)
	variable bin = 5
	const1= {"K1 < K4","K4 < 0.5","K1 > 0","K2 > 0.05","K5>0.05","K0>K3","K3 > 0.4*K0","K2 < 0.25","K5 < 0.25"}
	
	String name = "fermi_velocity_2_"+num2str(temperature[index])+"K_"
	String /G namelist = ""
	
	Make /O /N =(7,steps) $(name+"fitpar")
	wave temp3 =  $(name+"fitpar")
	Make /O /N = 7 temp_fitpar
	//wave fitpar_guess = $("fermi_velocity_"+num2str(temperature[index])+"K_fitpar")
	temp_fitpar[] = fit_pars_1[p][index] //this is from first run result
	
	variable i
	
	for(i=0;i<steps;i++)
		string temp = (name+num2str(i))
		namelist =  addlistItem(temp,namelist,";",inf)
		Make /O /N = (dimsize(dataset,1)) $temp = 0
		wave temp1 =  $temp
		Setscale /P x, dimoffset(dataset,1),dimdelta(dataset,1),temp1
		variable bin_temp
		
		if(mod(bin,2) ==0)
		for(bin_temp = -bin/2; bin_temp<= bin/2; bin_temp++)
			temp1[] += dataset(energystart+i*stepsize+bin_temp*DimDelta(dataset,0))[p][index]
		endfor
		else
		for(bin_temp = -(bin-1)/2; bin_temp<= (bin-1)/2; bin_temp++)
			temp1[] += dataset(energystart+i*stepsize+bin_temp*DimDelta(dataset,0))[p][index]
		endfor
		endif
	endfor 
	
	
	if(index < 13 || index == 17)
	windowlist = addlistItem(IB_iterative_fit("namelist",temp_fitpar,"ZC_loroff2sym","lor;lor;shft", fitStart_val = 112, fitStop_val = 225,vanilla = 1),windowlist)
	else
	windowlist = addlistItem(IB_iterative_fit("namelist",temp_fitpar,"ZC_loroff2sym","lor;lor;shft", fitStart_val = 112, fitStop_val = 206,vanilla = 1),windowlist)
	endif
	
	for(i=0;i<steps;i++)
		wave temp2 = $(name+num2str(i)+"_cf_fit")
		temp3[][i] = temp2[p]
	
	endfor
	
	setscale /P y, energystart, stepsize, temp3
	
	//displaying MDC peak positions in one window
	if(index ==0)
	Dowindow /K $("fermi_velocity_renormalize_2")
		display /N = $("fermi_velocity_renormalize_2") temp3[1][]
		Modifygraph /W = $("fermi_velocity_renormalize_2") rgb($(name+"fitpar")) = (65535-2500*index,65535-2500*index,65535-2500*index)
	else
		appendtograph /W = $("fermi_velocity_renormalize_2") temp3[1][]
		Modifygraph /W = $("fermi_velocity_renormalize_2") rgb($(name+"fitpar")) = (65535-2500*index,65535-2500*index,2500*index)
   endif
   
   //displaying fitting results on top of the cut
   string cutname = (nameofwave(dataset)+"_"+num2str(temperature[index])+"K")
   Duplicate /O /R = [][][index,index] dataset, $cutname
   Redimension /N = (dimsize(dataset,0),dimsize(dataset,1)) $cutname 
   //Make /O /N = (steps) $("scale_"+num2str(temperature[index])+"K")
   //wave scale = $("scale_"+num2str(temperature[index])+"K")
   //scale[] = energystart + stepsize*p
   IB_kill_list(winlist(cutname+"*",";","WIN:1"))
   display /N = $(cutname) temp3[1][]
   Modifygraph /W =$(cutname+"0") mode = 3, marker = 41, msize  = 2
   appendimage $cutname
	
	//this is for fitting the lines for each temperature and get velocity and crossing energy
	Make /O /N = (steps) $(nameofwave(temp3)+"_k"), $(nameofwave(temp3)+"_E")
	wave temp5 = $(nameofwave(temp3)+"_k")//k positions
	temp5[] = temp3[1][p] //y wave
	wave temp6 = $(nameofwave(temp3)+"_E")
	temp6[] = energystart+stepsize*p
	
	make /O /N = (18) fervelocity1, EB1
	make /O /N = 18 fervelocity1_error, EB1_error
	CurveFit/TBOX=768 line temp6 /X = temp5 /D 
	wave W_coef,W_sigma
	fervelocity1[index] = W_coef[1]
	fervelocity1_error[index] = W_sigma[1]
	EB1[index] = W_coef[0]
	EB1_error[index] = W_sigma[0]
//endfor

//Dowindow /F $("fermi_velocity_renormalize_2")

//////////////displaying fitting results
//Dowindow  /K $("fermiv1_v_temperature")
//display  /N =$("fermiv1_vs_temperature") fervelocity1[] vs temperature
//ModifyGraph /W =$("fermiv1_vs_temperature")  mode=3,marker=41,msize=2,rgb=(0,0,65535);DelayUpdate
//ErrorBars fervelocity1 Y,wave=(fervelocity1_error,fervelocity1_error)
////
//Dowindow  /K $("EB1_vs_temperature")
//display /N=$("EB1_vs_temperature") EB1[] vs temperature
//ModifyGraph /W = $("EB1_vs_temperature") mode=3,marker=41,msize=2,rgb=(0,20000,65535);DelayUpdate
//ErrorBars EB1 Y,wave=(EB1_error,EB1_error)







end