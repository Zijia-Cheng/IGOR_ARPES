#pragma rtGlobals=3		// Use modern global access method and strict wave access.


Menu "Y_ARPES"
	"Map_InPlane",  Y_MapInPlane_Start()
End




///////////////////////////////////
////  Convert full energy range to k 
///////////////////////////////////


Function	Y_MapIP_A2K_MapTo3D()
	string s_fld = GetDataFolder(1)

	SetDataFolder root:Y_ARPES:MapInPlane
	string/G s_MapInPlane_WorkFld
	string/G s_MapInPlane_ARP3D_name

	SetDataFolder $s_MapInPlane_WorkFld
	wave w3D_ang = $s_MapInPlane_ARP3D_name
	
	SetDataFolder :Misc
	NVAR v_eng, v_eng_width_flg

	variable v_eng_init			= v_eng
	variable v_eng_width_flg_init	= v_eng_width_flg
	v_eng_width_flg = 0	// Set to "No energy integral window" during mapping

	SetDataFolder ::Ang2K
	Duplicate/O w3D_ang, $(s_MapInPlane_ARP3D_name + "_3DEk")
	wave w3D_Ek = $(s_MapInPlane_ARP3D_name + "_3DEk")
	wave w2D_kxky = $(s_MapInPlane_ARP3D_name + "_kxky")
	
	SetDataFolder :Misc
	NVAR v_mesh_kx, v_mesh_ky, v_ky_max, v_ky_min, v_kx_max, v_kx_min, v_k_range_flag
	
	variable v_k_range_flag_init	= v_k_range_flag
	v_k_range_flag = 1	// Manually-fixed k_range
	
	
	Redimension/N=(-1, v_mesh_kx, v_mesh_ky) w3D_Ek
	SetScale/I y v_kx_min, v_kx_max,"", w3D_Ek
	SetScale/I z v_ky_min, v_ky_max,"", w3D_Ek
	
	variable v_eng_start	= DimOffset(w3D_Ek, 0)
	variable v_eng_delta	= DimDelta(w3D_Ek, 0)
	variable v_eng_pnts	= DimSize(w3D_Ek, 0)
	variable i
	for(i=0; i<v_eng_pnts; i+=1)
		v_eng = v_eng_start + i * v_eng_delta
		
		Y_MapIP_2D_eng()
		Y_MapIP_A2K_MapOntoImage()
		w3D_Ek[i][][] = w2D_kxky[q][r]

	endfor
	
	
	//Initialize parameters
	v_eng			= v_eng_init
	v_eng_width_flg 	= v_eng_width_flg_init
	v_k_range_flag		= v_k_range_flag_init
	
	Y_MapIP_2D_eng()
	Y_MapIP_A2K_MapOntoImage()
	
	SetDataFolder $s_fld
End




/////////////////////////////////
////  Convert constant-energy map 
/////////////////////////////////


Function Y_MapIP_A2K()
	string s_fld = GetDataFolder(1)

	SetDataFolder root:Y_ARPES:MapInPlane
	string/G s_MapInPlane_WorkFld
	string/G s_MapInPlane_ARP3D_name

	SetDataFolder $(s_MapInPlane_WorkFld + "Ang2K:Misc")
	
	Y_MapIP_A2K_BoundaryCalc()
	Y_MapIP_A2K_MapOntoImage()
	Y_MapIP_A2K_Graph_SetRange()
	
	SetDataFolder $s_fld
End


Function Y_MapIP_A2K_ChangeType()
	string s_fld = GetDataFolder(1)

	SetDataFolder root:Y_ARPES:MapInPlane
	string/G s_MapInPlane_WorkFld
	string/G s_MapInPlane_ARP3D_name

	SetDataFolder $(s_MapInPlane_WorkFld + "Ang2K:Misc")
	variable/G v_type, v_DA
	
	if(v_type == 0 && v_DA == 0)
		SetVariable setvar_xi 		title=" \\F'Symbol'x\\K(48896,65280,48896)\\Z070\\M\\K(0,0,0)\\F'Times New Roman' : Tilt angle", 		win=P_ARPES_Map_InPlane
		SetVariable setvar_xi0 		title=" \\F'Symbol'x\\F'Times New Roman'\\Z070\\M : Tilt origin", 										win=P_ARPES_Map_InPlane
		SetVariable setvar_beta0 	title=" \\F'Symbol'b\\F'Times New Roman'\\Z070\\K(0,0,0)\\]0 : Cut-angle origin",						win=P_ARPES_Map_InPlane
		SetVariable setvar_chi		title=" \\K(48896,65280,48896)\\F'Symbol'c\\F'Times New Roman'\\Z070\\]0 : Rotary angle",				win=P_ARPES_Map_InPlane
		SetVariable setvar_chi0 	title=" \\K(48896,65280,48896)\\F'Symbol'c\\F'Times New Roman'\\Z070\\]0 : Rotary origin",			win=P_ARPES_Map_InPlane
		
	elseif(v_type == 1 && v_DA == 0)
		SetVariable setvar_xi 		title=" \\F'Symbol'x\\K(48896,65280,48896)\\Z070\\M\\K(0,0,0)\\F'Times New Roman' : Rotary angle", 	win=P_ARPES_Map_InPlane
		SetVariable setvar_xi0 		title=" \\F'Symbol'x\\F'Times New Roman'\\Z070\\M : Rotary origin", 									win=P_ARPES_Map_InPlane
		SetVariable setvar_beta0 	title=" \\F'Symbol'b\\F'Times New Roman'\\Z070\\K(0,0,0)\\]0 : Cut-angle origin",						win=P_ARPES_Map_InPlane
		SetVariable setvar_chi		title=" \\K(48896,65280,48896)\\F'Symbol'c\\F'Times New Roman'\\Z070\\]0 : Rotary angle",				win=P_ARPES_Map_InPlane
		SetVariable setvar_chi0 	title=" \\K(48896,65280,48896)\\F'Symbol'c\\F'Times New Roman'\\Z070\\]0 : Rotary origin",			win=P_ARPES_Map_InPlane
	elseif(v_type == 0 && v_DA == 1)
		SetVariable setvar_xi 		title=" \\F'Symbol'x\\K(48896,65280,48896)\\Z070\\M\\K(0,0,0)\\F'Times New Roman' : Tilt angle", 		win=P_ARPES_Map_InPlane
		SetVariable setvar_xi0 		title=" \\F'Symbol'x\\F'Times New Roman'\\Z070\\M : Tilt origin", 										win=P_ARPES_Map_InPlane
		SetVariable setvar_beta0 	title=" \\K(48896,65280,48896)\\F'Symbol'b\\F'Times New Roman'\\Z070\\]0 : Cut-angle origin",			win=P_ARPES_Map_InPlane
		SetVariable setvar_chi		title=" \\F'Symbol'c\\F'Times New Roman'\\Z07\\K(48896,65280,48896)0\\K(0,0,0)\\]0 : Rotary angle",	win=P_ARPES_Map_InPlane
		SetVariable setvar_chi0 	title=" \\F'Symbol'c\\F'Times New Roman'\\Z070\\K(0,0,0)\\]0 : Rotary origin",						win=P_ARPES_Map_InPlane
	elseif(v_type == 1 && v_DA == 1)
		SetVariable setvar_xi 		title=" \\F'Symbol'x\\K(48896,65280,48896)\\Z070\\M\\K(0,0,0)\\F'Times New Roman' : Gonio angle",		win=P_ARPES_Map_InPlane
		SetVariable setvar_xi0 		title=" \\F'Symbol'x\\F'Times New Roman'\\Z070\\M : Gonio origin", 									win=P_ARPES_Map_InPlane
		SetVariable setvar_beta0 	title=" \\K(48896,65280,48896)\\F'Symbol'b\\F'Times New Roman'\\Z070\\]0 : Cut-angle origin",			win=P_ARPES_Map_InPlane
		SetVariable setvar_chi		title=" \\F'Symbol'c\\F'Times New Roman'\\Z07\\K(48896,65280,48896)0\\K(0,0,0)\\]0 : Rotary angle",	win=P_ARPES_Map_InPlane
		SetVariable setvar_chi0 	title=" \\F'Symbol'c\\F'Times New Roman'\\Z070\\K(0,0,0)\\]0 : Rotary origin",						win=P_ARPES_Map_InPlane
	endif
	
	SetDataFolder $s_fld
End


Function Y_MapIP_A2K_GraphDisp()
	string s_fld = GetDataFolder(1)

	SetDataFolder root:Y_ARPES:MapInPlane
	string/G s_MapInPlane_WorkFld
	string/G s_MapInPlane_ARP3D_name

	SetDataFolder $s_MapInPlane_WorkFld
	wave w2D_ang	= $(s_MapInPlane_ARP3D_name + "_AngCut")
	string s_FldName = GetDataFolder(0)
	SetDataFolder :Ang2K
	wave w2D_k	= $(s_MapInPlane_ARP3D_name + "_kxky")

	string s_WinName = "G_Ang2K_" + s_FldName

	if(strlen(WinList(s_WinName, ";", "WIN:1")))
		DoWindow/F $s_WinName
		if(!(strlen(ListMatch(ImageNameList("", ";"), nameofwave(w2D_ang))) * strlen(ListMatch(ImageNameList("", ";"), nameofwave(w2D_k)))))
			KillWindow $s_WinName
			Y_MapIP_A2K_GraphMake(s_WinName, w2D_ang, w2D_k)

		endif
	else

		Y_MapIP_A2K_GraphMake(s_WinName, w2D_ang, w2D_k)

	endif

	SetDataFolder $s_fld
End


Function Y_MapIP_A2K_Graph_SetRange()
	string s_fld = GetDataFolder(1)

	SetDataFolder root:Y_ARPES:MapInPlane
	string/G s_MapInPlane_WorkFld
	string/G s_MapInPlane_ARP3D_name

	SetDataFolder $s_MapInPlane_WorkFld
	wave w2D_ang	= $(s_MapInPlane_ARP3D_name + "_AngCut")
	string s_FldName = GetDataFolder(0)

	SetDataFolder :Misc
	variable/G v_cut_start, v_cut_end

	SetDataFolder ::Ang2K
	wave w2D_k	= $(s_MapInPlane_ARP3D_name + "_kxky")

	SetDataFolder :Contour
	variable/G v_ky_Boundary_max, v_ky_Boundary_min, v_kx_Boundary_max, v_kx_Boundary_min	
	variable v_alpha_min = DimOffset(w2D_ang, 0)
	variable v_alpha_max = v_alpha_min + DimDelta(w2D_ang, 0) * (DimSize(w2D_ang, 0) - 1)
	variable v_margin_r = 1.06 / 2
	variable v_hwidth_ang	= max(abs(v_alpha_max-v_alpha_min), abs(v_cut_end-v_cut_start)) * v_margin_r
	variable v_hwidth_k		= max(abs(v_kx_Boundary_max-v_kx_Boundary_min), abs(v_ky_Boundary_max-v_ky_Boundary_min)) * v_margin_r
	
	string s_WinName = "G_Ang2K_" + s_FldName

	if(strlen(WinList(s_WinName, ";", "WIN:1")))			
		SetAxis/W=$s_WinName b_ang_alpha 	(v_alpha_max+v_alpha_min)/2 - v_hwidth_ang, 	(v_alpha_max+v_alpha_min)/2 + v_hwidth_ang
		SetAxis/W=$s_WinName l_ang_beta 	(v_cut_start+v_cut_end)/2 - v_hwidth_ang, 		(v_cut_start+v_cut_end)/2 + v_hwidth_ang
		SetAxis/W=$s_WinName b_kx 		(v_kx_Boundary_max+v_kx_Boundary_min)/2 - v_hwidth_k, (v_kx_Boundary_max+v_kx_Boundary_min)/2 + v_hwidth_k
		SetAxis/W=$s_WinName r_ky 		(v_ky_Boundary_max+v_ky_Boundary_min)/2 - v_hwidth_k, (v_ky_Boundary_max+v_ky_Boundary_min)/2 + v_hwidth_k	
	endif
	
	SetDataFolder $s_fld
End


Function Y_MapIP_A2K_GraphMake(s_WinName, w2D_ang, w2D_k)
	string s_WinName
	wave w2D_ang, w2D_k
	
	SetDataFolder :Contour
	
	wave y_Boundary0, x_Boundary0, ky_Boundary0, kx_Boundary0
	wave y_Boundary1, x_Boundary1, ky_Boundary1, kx_Boundary1
	wave y_Boundary2, x_Boundary2, ky_Boundary2, kx_Boundary2
	wave y_Boundary3, x_Boundary3, ky_Boundary3, kx_Boundary3
	
	Display /W=(141.75,303.5,524.25,494.75)/N=$s_WinName/L=l_ang_beta/B=b_ang_alpha y_Boundary0 vs x_Boundary0
	DoWindow/T $s_WinName, s_WinName

	AppendToGraph/L=l_ang_beta/B=b_ang_alpha y_Boundary1 vs x_Boundary1
	AppendToGraph/L=l_ang_beta/B=b_ang_alpha y_Boundary2 vs x_Boundary2
	AppendToGraph/L=l_ang_beta/B=b_ang_alpha y_Boundary3 vs x_Boundary3
	AppendToGraph/R=r_ky/B=b_kx ky_Boundary0 vs kx_Boundary0
	AppendToGraph/R=r_ky/B=b_kx ky_Boundary1 vs kx_Boundary1
	AppendToGraph/R=r_ky/B=b_kx ky_Boundary2 vs kx_Boundary2
	AppendToGraph/R=r_ky/B=b_kx ky_Boundary3 vs kx_Boundary3
	AppendImage/B=b_ang_alpha/L=l_ang_beta w2D_ang
	ModifyImage $nameofwave(w2D_ang) ctab= {0,*,Terrain,1}
	AppendImage/B=b_kx/R=r_ky w2D_k
	ModifyImage $nameofwave(w2D_k) ctab= {0,*,Terrain,1}
	ModifyGraph margin(left)=37,margin(bottom)=34,margin(top)=3,margin(right)=37,gFont="Times New Roman"
	ModifyGraph gfSize=10,width={Aspect,2}
	ModifyGraph lSize=2
	ModifyGraph lStyle=11
	ModifyGraph rgb(y_Boundary1)=(0,52224,0),rgb(y_Boundary2)=(0,0,62976),rgb(y_Boundary3)=(65280,43520,0)
	ModifyGraph rgb(ky_Boundary1)=(0,52224,0),rgb(ky_Boundary2)=(0,0,62976),rgb(ky_Boundary3)=(65280,43520,0)
	ModifyGraph standoff=0
	ModifyGraph axThick=0.8
	ModifyGraph ZisZ=1
	ModifyGraph btLen=4
	ModifyGraph stLen=2
	ModifyGraph freePos(l_ang_beta)=0
	ModifyGraph freePos(b_ang_alpha)=0
	ModifyGraph freePos(r_ky)=0
	ModifyGraph freePos(b_kx)=0
	ModifyGraph axisEnab(l_ang_beta)={0,0.9}
	ModifyGraph axisEnab(b_ang_alpha)={0,0.45}
	ModifyGraph axisEnab(r_ky)={0,0.9}
	ModifyGraph axisEnab(b_kx)={0.55,1}
	TextBox/C/N=text_kx/F=0/H=11/B=1/A=MC/X=27.50/Y=-66.00 "\\f02k\\Bx\\M\\]0 (1/Å)"
	TextBox/C/N=text_ang_alpha/F=0/H=11/B=1/A=MC/X=-27.50/Y=-65.53 "\\F'Symbol'\\f02a \\]0(deg.) // slit"
	TextBox/C/N=text_ang_beta/O=90/F=0/H=11/B=1/A=MC/X=-59.80/Y=-5.00 "\\F'Symbol'\\f02b \\]0(deg.) \\F'Symbol'^\\]0 slit"
	TextBox/C/N=text_ky/O=-90/F=0/H=11/B=1/A=MC/X=59.80/Y=-5.00 "\\f02k\\By\\M\\]0 (1/Å)"
	SetDrawLayer UserFront
	SetDrawEnv linethick= 0.8
	DrawLine 0,0.1,0.45,0.1
	SetDrawEnv linethick= 0.8
	DrawLine 0.45,1,0.45,0.1
	SetDrawEnv linethick= 0.8
	DrawLine 0.55,1,0.55,0.1
	SetDrawEnv linethick= 0.8
	DrawLine 0.55,0.1,1,0.1

End




///////////////////////////////////////////////////////
////  Convert constant-energy map: Angular boundary --> k
///////////////////////////////////////////////////////


Function Y_MapIP_A2K_BoundaryMake()
	string s_fld = GetDataFolder(1)	

	SetDataFolder root:Y_ARPES:MapInPlane
	string/G s_MapInPlane_WorkFld
	string/G s_MapInPlane_ARP3D_name

	SetDataFolder $s_MapInPlane_WorkFld
	wave w2D_AngCut = $(s_MapInplane_ARP3D_name + "_AngCut")

	SetDataFolder :Ang2K
	wave w2D_kxky =  $(s_MapInplane_ARP3D_name + "_kxky")
	
	NewDataFolder/O/S :Contour
	variable/G v_ky_Boundary_min, v_ky_Boundary_max, v_kx_Boundary_min, v_kx_Boundary_max

	Duplicate/O w2D_AngCut, Boundary_alpha, Boundary_beta
	Redimension/N=(-1) Boundary_alpha
	MatrixTranspose Boundary_beta
	Redimension/N=(-1) Boundary_beta
	Duplicate/O Boundary_alpha, 	ky_Boundary2, kx_Boundary2, ky_Boundary0, kx_Boundary0, y_Boundary2, x_Boundary2, y_Boundary0, x_Boundary0
	Duplicate/O Boundary_beta, 		ky_Boundary3, kx_Boundary3, ky_Boundary1, kx_Boundary1, y_Boundary3, x_Boundary3, y_Boundary1, x_Boundary1
	Killwaves Boundary_beta, Boundary_alpha

	x_Boundary0 = x
	y_Boundary1 = x
	x_Boundary2 = x_Boundary0[numpnts(x_Boundary0) - 1- p]
	y_Boundary3 = y_Boundary1[numpnts(y_Boundary1) - 1- p]
	y_Boundary0 = y_Boundary1[0]
	x_Boundary1 = x_Boundary0[numpnts(x_Boundary0) -1]
	y_Boundary2 = y_Boundary1[numpnts(y_Boundary1) -1]
	x_Boundary3 = x_Boundary0[0]

	SetDataFolder $s_fld
End


Function Y_MapIP_A2K_BoundaryCalc()
	string s_fld = GetDataFolder(1)	

	SetDataFolder root:Y_ARPES:MapInPlane
	string/G s_MapInPlane_WorkFld
	string/G s_MapInPlane_ARP3D_name

	SetDataFolder $(s_MapInPlane_WorkFld + "Misc")
	variable/G v_eng
	
	SetDataFolder ::Ang2K
	wave w2D_kxky = $(s_MapInPlane_ARP3D_name + "_kxky")
	
	SetDataFolder :Misc
	variable/G v_type	//0 --> slit parallel to manipulator axis, 1 --> slit perpendicular to manipulator axis
	variable/G v_DA			// --> for future use
	variable/G v_delta, v_xi, v_xi0, v_beta0
	variable/G v_chi, v_chi0	// --> for future use
	variable/G v_hn, v_work	
	NewDataFolder/O/S ::Contour
	variable/G v_ky_Boundary_min, v_ky_Boundary_max, v_kx_Boundary_min, v_kx_Boundary_max
	
	variable i
	for(i=0; i<4; i+=1)
		wave w_y		= $("y_Boundary" + num2str(i))
		wave w_x	 	= $("x_Boundary" + num2str(i))
		wave w_ky	= $("ky_Boundary" + num2str(i))
		wave w_kx	= $("kx_Boundary" + num2str(i))
		
		Y_MapIP_A2K_Contour(w_y, w_x, w_ky, w_kx, v_type, v_DA, v_delta, v_xi-v_xi0, v_beta0, v_chi-v_chi0, v_hn, v_work, v_eng)
		
		if(i==0)
			v_kx_Boundary_min		= wavemin(w_kx)
			v_kx_Boundary_max	= wavemax(w_kx)
			v_ky_Boundary_min		= wavemin(w_ky)
			v_ky_Boundary_max	= wavemax(w_ky)
		else
			v_kx_Boundary_min		= min(v_kx_Boundary_min, wavemin(w_kx))			
			v_kx_Boundary_max	= max(v_kx_Boundary_max, wavemax(w_kx))
			v_ky_Boundary_min		= min(v_ky_Boundary_min, wavemin(w_ky))
			v_ky_Boundary_max	= max(v_ky_Boundary_max, wavemax(w_ky))	
		endif
	endfor
	
	SetDataFolder ::Misc
	variable/G v_k_range_flag
	if(!v_k_range_flag)
		variable/G v_ky_max, v_ky_min, v_kx_max, v_kx_min
		v_kx_max	= v_kx_Boundary_max
		v_kx_min	= v_kx_Boundary_min
		v_ky_max	= v_ky_Boundary_max
		v_ky_min	= v_ky_Boundary_min				
		
		Y_MapIP_Rangekxky()
		
	endif

	SetDataFolder $s_fld
End


Function Y_MapIP_A2K_Contour(w_y, w_x, w_ky, w_kx, v_type, v_DA, v_delta, v_xi, v_beta0, v_chi, v_hn, v_work, v_omega)
	Wave w_y, w_x, w_ky, w_kx
	variable v_type	//0 --> slit parallel to manipulator axis, 1 --> slit perpendicular to manipulator axis
	variable v_DA		//0 --> No deflector, 1 --> With deflector
	variable v_delta 
	variable v_xi	// xi - xi0
	variable v_beta0	// --> used when v_DA = 0
	variable v_chi		// --> used when v_DA = 1
	variable v_hn, v_work, v_omega
	
	
	w_kx[] = 0.513168 * sqrt(v_hn - v_work + v_omega) * Y_MapIP_k_dcos(v_type, v_DA, 0, v_delta, v_xi, v_chi, w_y[p] - v_beta0*(1-v_DA), w_x[p])
	w_ky[] = 0.513168 * sqrt(v_hn - v_work + v_omega) * Y_MapIP_k_dcos(v_type, v_DA, 1, v_delta, v_xi, v_chi, w_y[p] - v_beta0*(1-v_DA), w_x[p])	
End


Function Y_MapIP_k_dcos(v_type, v_DA, v_xy, v_delta, v_xi, v_chi, v_beta, v_alpha)
	variable v_type	//0 --> slit parallel to manipulator axis, 1 --> slit perpendicular to manipulator axis 
	variable v_DA		//0 --> No deflector, 1 --> With deflector
	variable v_xy		//0 --> x, 1 --> y
	variable v_delta, v_xi, v_alpha
	variable v_chi		// --> used when v_DA = 1
	variable v_beta	//beta - beta0
	
	variable v_sd	= sin(	v_delta	/180*pi)
	variable v_cd	= cos(	v_delta	/180*pi)
	variable v_sx	= sin(	v_xi		/180*pi)
	variable v_cx	= cos(	v_xi		/180*pi)
	variable v_sb	= sin(	v_beta	/180*pi)
	variable v_cb	= cos(	v_beta	/180*pi)
	variable v_sa	= sin(	v_alpha	/180*pi)
	variable v_ca	= cos(	v_alpha	/180*pi)
	variable v_sc	= sin(	v_chi	/180*pi)
	variable v_cc	= cos(	v_chi	/180*pi)


//0.513168*sqrt(hn-E_bin-work)*(-cos(phi)*cos(delta)*sin(alpha)+sin(phi)*cos(theta)*cos(delta)*cos(alpha)+sin(theta)*sin(delta)*cos(alpha))  kx
// (v_sx * v_cb * v_cd + v_sb * v_sd)* v_ca - v_cx * v_cd * v_sa		kx

//0.513168*sqrt(hn-E_bin-work)*(-cos(phi)*sin(delta)*sin(alpha)+sin(phi)*cos(theta)*sin(delta)*cos(alpha)-sin(theta)*cos(delta)*cos(alpha))  ky




	
	variable v_return

	if(v_type==0 && v_DA == 0)

		if(v_xy == 0)
			v_return = (v_sd * v_sb + v_cd * v_sx * v_cb) * v_ca - v_cd * v_cx *  v_sa
		elseif(v_xy == 1)
			v_return = (-v_cd * v_sb + v_sd * v_sx * v_cb) * v_ca - v_sd * v_cx *  v_sa
		endif

	elseif(v_type==1 && v_DA == 0)

		if(v_xy == 0)
			v_return = (v_sd * v_sx + v_cd * v_sb * v_cx) * v_ca + (-v_sd * v_cx + v_cd * v_sb * v_sx) * v_sa
		elseif(v_xy == 1)
			v_return = (-v_cd * v_sx + v_sd * v_sb * v_cx) * v_ca + (v_cd * v_cx + v_sd * v_sb * v_sx) * v_sa
		endif

	elseif(v_type==0 && v_DA == 1)

		variable v_polar_01 = sqrt(v_alpha^2 + v_beta^2)*pi/180
		if(v_xy == 0)
//			v_return = sin(v_polar_01) / v_polar_01 * (-v_cd * v_cx * v_alpha + v_sd * v_cc * v_beta - v_cd * v_sx * v_sc * v_beta)*pi/180 + cos(v_polar_01) * (v_sd * v_sc + v_cd * v_sx * v_cc)
			v_return = sinc(v_polar_01) * (-v_cd * v_cx * v_alpha + v_sd * v_cc * v_beta - v_cd * v_sx * v_sc * v_beta)*pi/180 + cos(v_polar_01) * (v_sd * v_sc + v_cd * v_sx * v_cc)
			elseif(v_xy == 1)
//			v_return = sin(v_polar_01/180*pi) / v_polar_01 * (-v_sd * v_cx * v_alpha - v_cd * v_cc * v_beta - v_sd * v_sx * v_sc * v_beta)*pi/180 - cos(v_polar_01) * (v_cd * v_sc - v_sd * v_sx * v_cc)
			v_return = sinc(v_polar_01) * (-v_sd * v_cx * v_alpha - v_cd * v_cc * v_beta - v_sd * v_sx * v_sc * v_beta)*pi/180 - cos(v_polar_01) * (v_cd * v_sc - v_sd * v_sx * v_cc)
		endif

	elseif(v_type==1 && v_DA == 1)
	
		variable v_polar_11 = sqrt(v_alpha^2 + v_beta^2)*pi/180
		if(v_xy == 0)
//			v_return = sin(v_polar_11) / v_polar_11 * (-v_cd * v_cx * v_beta - v_sd * v_cc * v_alpha + v_cd * v_sx * v_sc * v_alpha)*pi/180 + cos(v_polar_11) * (v_sd * v_sc + v_cd * v_sx * v_cc)
			v_return = sinc(v_polar_11) * (-v_cd * v_cx * v_beta - v_sd * v_cc * v_alpha + v_cd * v_sx * v_sc * v_alpha)*pi/180 + cos(v_polar_11) * (v_sd * v_sc + v_cd * v_sx * v_cc)
		elseif(v_xy == 1)
//			v_return = sin(v_polar_11) / v_polar_11 * (-v_sd * v_cx * v_beta + v_cd * v_cc * v_alpha + v_sd * v_sx * v_sc * v_alpha)*pi/180 - cos(v_polar_11) * (v_cd * v_sc - v_sd * v_sx * v_cc)
			v_return = sinc(v_polar_11) * (-v_sd * v_cx * v_beta + v_cd * v_cc * v_alpha + v_sd * v_sx * v_sc * v_alpha)*pi/180 - cos(v_polar_11) * (v_cd * v_sc - v_sd * v_sx * v_cc)
		endif

	endif
	
	return v_return
end




/////////////////////////////////////////////////////////////////
////  Convert constant-energy map: Map angular distribution to k space
/////////////////////////////////////////////////////////////////


Function Y_MapIP_A2K_MapOntoImage()
	string s_fld = GetDataFolder(1)	

	SetDataFolder root:Y_ARPES:MapInPlane
	string/G s_MapInPlane_WorkFld
	string/G s_MapInPlane_ARP3D_name

	SetDataFolder $s_MapInPlane_WorkFld
	wave w2D_ang = $(s_MapInPlane_ARP3D_name + "_AngCut")	

	SetDataFolder :Misc
	variable/G v_eng
	
	SetDataFolder ::Ang2K
	wave w2D_k = $(s_MapInPlane_ARP3D_name + "_kxky")
	
	SetDataFolder :Misc
	variable/G v_type	//0 --> slit parallel to manipulator axis, 1 --> slit perpendicular to manipulator axis
	variable/G v_DA	//0 --> No deflector, 1 --> With deflector
	variable/G v_delta, v_xi, v_xi0
	variable/G v_beta0			// --> used when v_DA = 0
	variable/G v_chi, v_chi0	// --> used when v_DA = 1
	variable/G v_hn, v_work	
	
	Y_MapIP_A2K_ImageCalc(w2D_k, w2D_ang, v_type, v_delta, v_xi-v_xi0, v_beta0, v_chi-v_chi0, v_hn, v_work, v_eng, v_DA)

	SetDataFolder $s_fld
End


Function Y_MapIP_A2K_ImageCalc(w2D_k, w2D_ang, v_type, v_delta, v_xi, v_beta0, v_chi, v_hn, v_work, v_omega, v_DA)
	wave w2D_k, w2D_ang
	variable v_type, v_delta, v_xi
	variable v_beta0	// --> used when v_DA = 0
	variable v_chi		// --> used when v_DA = 1
	variable v_hn, v_work, v_omega, v_DA

	variable v_sd	= sin(	v_delta	/180*pi)
	variable v_cd	= cos(	v_delta	/180*pi)
	variable v_sx	= sin(	v_xi		/180*pi)
	variable v_cx	= cos(	v_xi		/180*pi)
	variable v_sc	= sin(	v_chi	/180*pi)
	variable v_cc	= cos(	v_chi	/180*pi)

	variable v_k = 0.513168 * sqrt(v_hn - v_work + v_omega)

//	kx''' = x / v_k
//	ky''' = y / v_k
//	kz''' = sqrt(1 - kx'''^2 - ky'''^2) = sqrt(v_k^2 -x^2 - y^2) / v_k


//	type = 0		(slit) parallel to (rotary axis)
//	v_alpha	= 180/pi * asin( (-v_cx*v_cd*x-v_cx*v_sd*y+v_sx*sqrt(v_k^2-x^2-y^2))/v_k ) 	
//	v_beta	= v_beta0 + 180/pi * atan( (v_sd*x-v_cd*y)/(v_sx*v_cd*x+v_sx*v_sd*y+v_cx*sqrt(v_k^2-x^2-y^2)) ) 
//	v_alpha	= 180/pi * asin( (v_sx*sqrt(v_k^2-x^2-y^2) - v_cx*(v_cd*x+v_sd*y) ) / v_k )
//	v_beta	= v_beta0 + 180/pi * atan( (v_sd*x-v_cd*y) / (v_sx*(v_cd*x+v_sd*y) + v_cx*sqrt(v_k^2-x^2-y^2)) )

//	type = 1		(slit) perp to (rotary axis)
//	v_alpha	= 180/pi * asin( (v_sx*v_cd*sin( atan( (v_cd*x+v_sd*y)/sqrt(v_k^2-x^2-y^2) ) ) - v_cx*v_sd) * x/v_k   +  (v_sx*v_sd*sin( atan( (v_cd*x+v_sd*y)/sqrt(v_k^2-x^2-y^2) ) ) + v_cx*v_cd) * y/v_k  +  v_sx*cos( atan( (v_cd*x+v_sd*y)/sqrt(v_k^2-x^2-y^2))) * sqrt(v_k^2-x^2 -y^2)/v_k )
//	v_beta	= v_beta0 + 180/pi * atan( (v_cd * x + v_sd * y)/sqrt(v_k^2-x^2-y^2) )
//
//	eqiv
//	v_alpha	= 180/pi * asin( (v_sx*sqrt(v_k^2-(v_sd*x-v_cd*y)^2) - v_cx*(v_sd*x-v_cd*y)) / v_k )
//	v_beta	= v_beta0 + 180/pi * atan( (v_cd * x + v_sd * y)/sqrt(v_k^2-x^2-y^2) )





	if(v_DA == 0 && v_type==0)

		//w2D_k = Interp2D(w2D_ang, 180/pi * asin( (-v_cx*v_cd*x - v_cx*v_sd*y + v_sx*sqrt(v_k^2-x^2-y^2))/v_k ),  v_beta0 + 180/pi * atan( (v_sd*x-v_cd*y)/(v_sx*v_cd*x+v_sx*v_sd*y+v_cx*sqrt(v_k^2-x^2-y^2)) ) )
		w2D_k = Interp2D(w2D_ang, 180/pi * asin( ( v_sx*sqrt(v_k^2-x^2-y^2)-v_cx*(v_cd*x+v_sd*y) ) / v_k ),  v_beta0 + 180/pi * atan(  (v_sd*x-v_cd*y)/(v_sx*v_cd*x+v_sx*v_sd*y+v_cx*sqrt(v_k^2-x^2-y^2)) ) )


	elseif(v_DA == 0 && v_type==1)

		//w2D_k = Interp2D(w2D_ang, 180/pi * asin( (v_sx*v_cd*sin( atan( (v_cd*x+v_sd*y)/sqrt(v_k^2-x^2-y^2) ) ) - v_cx*v_sd) * x/v_k   +  (v_sx*v_sd*sin( atan( (v_cd*x+v_sd*y)/sqrt(v_k^2-x^2-y^2) ) ) + v_cx*v_cd) * y/v_k  +  v_sx*cos( atan( (v_cd*x+v_sd*y)/sqrt(v_k^2-x^2-y^2))) * sqrt(v_k^2-x^2 -y^2)/v_k ),  v_beta0 + 180/pi * atan( (v_cd * x + v_sd * y)/sqrt(v_k^2-x^2-y^2) ) )
		w2D_k = Interp2D(w2D_ang, 180/pi * asin( (v_sx*sqrt(v_k^2-(v_sd*x-v_cd*y)^2) - v_cx*(v_sd*x-v_cd*y)) / v_k ), v_beta0 + 180/pi * atan( (v_cd * x + v_sd * y)/sqrt(v_k^2-x^2-y^2) ) )

	elseif(v_DA == 1)

		variable v_11	= v_cx * v_cd
		variable v_12 	= v_cx * v_sd
		variable v_13	= -v_sx
		variable v_21	= v_sc * v_sx * v_cd - v_cc * v_sd
		variable v_22	= v_sc * v_sx * v_sd + v_cc * v_cd
		variable v_23	= v_sc * v_cx
		variable v_31	= v_cc * v_sx * v_cd + v_sc * v_sd
		variable v_32	= v_cc * v_sx * v_sd - v_sc * v_cd
		variable v_33	= v_cc * v_cx

		// kx	= v_11 * x + v_12 * y + v_13 * sqrt(v_k^2 - x^2 - y^2) 
		// ky	= v_21 * x + v_22 * y + v_23 * sqrt(v_k^2 - x^2 - y^2)
		// kz	= v_31 * x + v_32 * y + v_33 * sqrt(v_k^2 - x^2 - y^2)
		
		if(v_type == 1)
		
			w2D_k = Interp2D(w2D_ang,180/pi*acos((v_31*x+v_32*y+v_33*sqrt(v_k^2-x^2-y^2))/v_k)*(v_21*x+v_22*y+v_23*sqrt(v_k^2-x^2-y^2))/sqrt(v_k^2-(v_31*x+v_32*y+v_33*sqrt(v_k^2-x^2-y^2))^2), -180/pi*acos((v_31*x+v_32*y+v_33*sqrt(v_k^2-x^2-y^2))/v_k)*(v_11*x+v_12*y+v_13*sqrt(v_k^2-x^2-y^2))/sqrt(v_k^2-(v_31*x+v_32*y+v_33*sqrt(v_k^2-x^2-y^2))^2))
		
		elseif(v_type == 0)

			w2D_k = Interp2D(w2D_ang, -180/pi*acos((v_31*x+v_32*y+v_33*sqrt(v_k^2-x^2-y^2))/v_k)*(v_11*x+v_12*y+v_13*sqrt(v_k^2-x^2-y^2))/sqrt(v_k^2-(v_31*x+v_32*y+v_33*sqrt(v_k^2-x^2-y^2))^2), -180/pi*acos((v_31*x+v_32*y+v_33*sqrt(v_k^2-x^2-y^2))/v_k)*(v_21*x+v_22*y+v_23*sqrt(v_k^2-x^2-y^2))/sqrt(v_k^2-(v_31*x+v_32*y+v_33*sqrt(v_k^2-x^2-y^2))^2))

		endif
	endif

End




/////////////////////////////////////////
////  Display slices of ang-ang-eng 3D data
/////////////////////////////////////////


Function Y_ButtonProc_MapIP_2D_GraphDisp(ctrlName) : ButtonControl
	String ctrlName	//button_Disp_Eng, button_Disp_Ang, button_Disp_Cut
	
	string s_slice = ctrlName[12,14]
	Y_MapIP_2D_GraphDisp(s_slice)
	
End


Function Y_MapIP_2D_eng()
	string s_fld = GetDataFolder(1)
	
	SetDataFolder root:Y_ARPES:MapInPlane
	string/G s_MapInPlane_WorkFld
	string/G s_MapInPlane_ARP3D_name

	SetDataFolder $s_MapInPlane_WorkFld
	wave w3D		= $s_MapInPlane_ARP3D_name
	wave w2D		= $(s_MapInPlane_ARP3D_name + "_AngCut")
	wave w1D		= :Misc:w_eng
	NVAR v_val		= :Misc:v_eng
	NVAR v_width		= :Misc:v_eng_width
	NVAR v_width_flg	= :Misc:v_eng_width_flg
	
	w2D = w3D(v_val)[p][q]

	if(v_width_flg)
		variable v_p_s = x2pnt(w1D, v_val - v_width/2)
		variable v_p_e = x2pnt(w1D, v_val + v_width/2)
		
		w2D = 0
		variable i
		for(i=v_p_s; i<(v_p_e+1); i+=1)
			w2D += w3D[i][p][q]
		endfor
		w2D /= (v_p_e - v_p_s + 1)
	endif

	SetDataFolder $s_fld
End


Function Y_MapIP_2D_ang()
	string s_fld = GetDataFolder(1)
	
	SetDataFolder root:Y_ARPES:MapInPlane
	string/G s_MapInPlane_WorkFld
	string/G s_MapInPlane_ARP3D_name

	SetDataFolder $s_MapInPlane_WorkFld
	wave w3D		= $s_MapInPlane_ARP3D_name
	wave w2D		= $(s_MapInPlane_ARP3D_name + "_EngCut")
	wave w1D		= :Misc:w_ang
	NVAR v_val		= :Misc:v_ang
	NVAR v_width		= :Misc:v_ang_width
	NVAR v_width_flg	= :Misc:v_ang_width_flg
	
	w2D = w3D[p](v_val)[q]

	if(v_width_flg)
		variable v_p_s = x2pnt(w1D, v_val - v_width/2)
		variable v_p_e = x2pnt(w1D, v_val + v_width/2)
		
		w2D = 0
		variable i
		for(i=v_p_s; i<(v_p_e+1); i+=1)
			w2D += w3D[p][i][q]
		endfor
		w2D /= (v_p_e - v_p_s + 1)
	endif

	SetDataFolder $s_fld
End


Function Y_MapIP_2D_cut()
	string s_fld = GetDataFolder(1)
	
	SetDataFolder root:Y_ARPES:MapInPlane
	string/G s_MapInPlane_WorkFld
	string/G s_MapInPlane_ARP3D_name

	SetDataFolder $s_MapInPlane_WorkFld
	wave w3D		= $s_MapInPlane_ARP3D_name
	wave w2D		= $(s_MapInPlane_ARP3D_name + "_EngAng")
	wave w1D		= :Misc:w_cut
	NVAR v_val		= :Misc:v_cut
	NVAR v_width		= :Misc:v_cut_width
	NVAR v_width_flg	= :Misc:v_cut_width_flg
	
	w2D = w3D[p][q](v_val)

	if(v_width_flg)
		variable v_p_s = x2pnt(w1D, v_val - v_width/2)
		variable v_p_e = x2pnt(w1D, v_val + v_width/2)
		
		w2D = 0
		variable i
		for(i=v_p_s; i<(v_p_e+1); i+=1)
			w2D += w3D[p][q][i]
		endfor
		w2D /= (v_p_e - v_p_s + 1)
	endif

	SetDataFolder $s_fld
End


Function Y_MapIP_2D_RangeSlider(ctrlName, v_IntWidth)
	string ctrlName
	variable v_IntWidth
	
	if(stringmatch(ctrlName, "*_width_*"))
		string s_var = ctrlName[strlen(ctrlName)-3, strlen(ctrlName)-1]	// = eng or cut or ang
		
		string s_fld = GetDataFolder(1)
		
		SetDataFolder root:Y_ARPES:MapInPlane
		string/G s_MapInPlane_WorkFld
		SetDataFolder $(s_MapInPlane_WorkFld + "Misc:")

		NVAR v_var		= $("v_" + s_var)
		NVAR v_var_start	= $("v_" + s_var + "_start")
		NVAR v_var_end	= $("v_" + s_var + "_end")
		
		// case when invoked from "set cut"
		if(stringmatch(ctrlName, "_width_cut"))
			NVAR v_cut_width 
			v_cut_width = 0
		endif

		variable v_var_min = min(v_var_start, v_var_end) + v_IntWidth/2
		variable v_var_max = max(v_var_start, v_var_end) - v_IntWidth/2
		
		if(v_var < v_var_min)
			v_var = v_var_min
		elseif(v_var > v_var_max)
			v_var = v_var_max
		endif

		string s_command = ""
		if(stringmatch(s_var, "cut"))
			s_command = "Slider slider_cut limits={" + num2str(v_var_min) + ", " + num2str(v_var_max) + ", v_cut_delta}, win = P_ARPES_Map_InPlane"
			Execute s_command
			s_command = "Slider slider_width_cut limits	= {0, abs(v_cut_start - v_cut_end)/2, v_cut_delta * 2}, win = P_ARPES_Map_InPlane"
			Execute s_command
		else
			s_command = "Slider slider_" +  s_var + " limits={" + num2str(v_var_min) + ", " + num2str(v_var_max) + ", 0}, win = P_ARPES_Map_InPlane"
			Execute s_command
		endif
		
		SetDataFolder $s_fld
	endif

End


Function Y_SetVarProc_MapIP_2D_eng(ctrlName,varNum,varStr,varName) : SetVariableControl
	String ctrlName
	Variable varNum
	String varStr
	String varName
	
	Y_MapIP_2D_GraphDisp("Eng")
	Y_MapIP_2D_RangeSlider(ctrlName, varNum)
	Y_MapIP_2D_eng()
	Y_MapIP_A2K()
End


Function Y_SetVarProc_MapIP_2D_ang(ctrlName,varNum,varStr,varName) : SetVariableControl
	String ctrlName
	Variable varNum
	String varStr
	String varName
	
	Y_MapIP_2D_GraphDisp("Ang")
	Y_MapIP_2D_RangeSlider(ctrlName, varNum)
	Y_MapIP_2D_ang()
End


Function Y_SetVarProc_MapIP_2D_cut(ctrlName,varNum,varStr,varName) : SetVariableControl
	String ctrlName
	Variable varNum
	String varStr
	String varName
	
	Y_MapIP_2D_GraphDisp("Cut")
	Y_MapIP_2D_RangeSlider(ctrlName, varNum)
	Y_MapIP_2D_cut()
End


Function Y_SliderProc_MapIP_2D_eng(ctrlName,sliderValue,event) : SliderControl
	String ctrlName
	Variable sliderValue
	Variable event	// bit field: bit 0: value set, 1: mouse down, 2: mouse up, 3: mouse moved

	if(event %& 2)	
		Y_MapIP_2D_GraphDisp("Eng")
	endif

	Y_MapIP_2D_RangeSlider(ctrlName, sliderValue)
	Y_MapIP_2D_eng()
	Y_MapIP_A2K()

	return 0
End


Function Y_SliderProc_MapIP_2D_ang(ctrlName,sliderValue,event) : SliderControl
	String ctrlName
	Variable sliderValue
	Variable event	// bit field: bit 0: value set, 1: mouse down, 2: mouse up, 3: mouse moved

	if(event %& 2)	
		Y_MapIP_2D_GraphDisp("Ang")
	endif	

	Y_MapIP_2D_RangeSlider(ctrlName, sliderValue)
	Y_MapIP_2D_ang()

	return 0
End


Function Y_SliderProc_MapIP_2D_cut(ctrlName,sliderValue,event) : SliderControl
	String ctrlName
	Variable sliderValue
	Variable event	// bit field: bit 0: value set, 1: mouse down, 2: mouse up, 3: mouse moved

	if(event %& 2)	
		Y_MapIP_2D_GraphDisp("Cut")
	endif		

	Y_MapIP_2D_RangeSlider(ctrlName, sliderValue)
	Y_MapIP_2D_cut()

	return 0
End


Function Y_MapIP_2D_GraphDisp(s_slice)
	string s_slice	// = Cut, Eng, Ang
	
	variable v_slice = WhichListItem(s_slice, "Cut;Eng;Ang")		// Cut --> 0; Eng --> 1; Ang --> 2
	string s_slice_name = Stringfromlist(v_slice, "_EngAng;_AngCut;_EngCut") // Cut --> _EngAng; Eng --> _AngCut; Ang --> _EngCut
	
	string s_fld = GetDataFolder(1)

	SetDataFolder root:Y_ARPES:MapInPlane
	string/G s_MapInPlane_WorkFld
	string/G s_MapInPlane_ARP3D_name

	SetDataFolder $s_MapInPlane_WorkFld
	wave w2D	= $(s_MapInPlane_ARP3D_name + s_slice_name)
	string s_FldName = GetDataFolder(0)

	string s_WinName = "G_Slice_"  + s_slice	 + "_" + s_FldName

	if(strlen(WinList(s_WinName, ";", "WIN:1")))
		DoWindow/F $s_WinName
		if(!strlen(ListMatch(ImageNameList("", ";"), nameofwave(w2D))))
			KillWindow $s_WinName
			Y_MapIP_2D_GraphMake(s_WinName, w2D, v_slice)
		endif
	else
			Y_MapIP_2D_GraphMake(s_WinName, w2D, v_slice)
	endif

	SetDataFolder $s_fld
End


Function Y_MapIP_2D_GraphMake(s_WinName, w2D, v_slice)
	string s_WinName
	wave w2D
	variable v_slice
	
	string s_axis_label_bottom	 = StringfromList(v_slice, "Energy (eV);\\F'Symbol'\\f02a \\]0(deg.) // slit;Energy (eV)")
	string s_axis_label_left		 = StringfromList(v_slice, "\\F'Symbol'\\f02a \\]0(deg.) // slit;\F'Symbol'\f02b \]0(deg.) \F'Symbol'^\]0 slit;\F'Symbol'\f02b \]0(deg.) \F'Symbol'^\]0 slit")

	Display /W=(138.75,90.5,333,281.75)/N=$s_WinName as s_WinName
	AppendImage w2D
	ModifyImage $nameofwave(w2D) ctab= {0,*,Terrain,1}
	ModifyGraph margin(left)=37,margin(bottom)=34,margin(top)=9,margin(right)=9,gFont="Times New Roman"
	ModifyGraph gfSize=10,width={Aspect,1}
	ModifyGraph mirror=2
	ModifyGraph standoff=0
	ModifyGraph axThick=0.8
	ModifyGraph ZisZ=1
	ModifyGraph btLen=4
	ModifyGraph stLen=2
	ModifyGraph freePos(left)=0
	ModifyGraph freePos(bottom)=0
	TextBox/C/N=text_bottom/F=0/H=11/B=1/A=MC/X=0.00/Y=-65.53  s_axis_label_bottom
	TextBox/C/N=text_left/O=90/F=0/H=11/B=1/A=MC/X=-70.71/Y=0.00  s_axis_label_left
	Label left "\\u#2"
	Label bottom "\\u#2"
	SetDrawLayer UserFront
End


Function Y_CheckProc_MapIP_Width(ctrlName,checked) : CheckBoxControl
	String ctrlName
	Variable checked
	
	//ctrlName = "setvar_width_XXX"	XXX = eng / ang / cut
	
	string s_setvar	= "setvar_width_"	+ ctrlName[12, 14]
	string s_setslider	= "slider_width_"	+ ctrlName[12, 14]

	SetVariable	$s_setvar disable 		= (2-2*checked), win = P_ARPES_Map_InPlane		//disable == 2  enable == 0
	Slider 		$s_setslider disable	= (2-2*checked), win = P_ARPES_Map_InPlane
End




////////////////////////////////////////
//// Set Range of the kx-ky image  
////////////////////////////////////////


Function Y_MapIP_Rangekxky()
	string s_fld = GetDataFolder(1)
	
	SetDataFolder root:Y_ARPES:MapInPlane
	string/G s_MapInPlane_WorkFld
	string/G s_MapInPlane_ARP3D_name

	SetDataFolder $(s_MapInPlane_WorkFld + "Ang2K")
	wave w2D_kxky = $(s_MapInPlane_ARP3D_name + "_kxky")
	SetDataFolder :Misc
	variable/G v_kx_min, v_kx_max, v_ky_min, v_ky_max
	
	SetScale/I x v_kx_min, v_kx_max,"", w2D_kxky
	SetScale/I y v_ky_min, v_ky_max,"", w2D_kxky
	
	SetDataFolder $s_fld
End




////////////////////////////////////////
//// Set Mesh number of the kx-ky image  
////////////////////////////////////////


Function Y_MapIP_RangeMesh()
	string s_fld = GetDataFolder(1)
	
	SetDataFolder root:Y_ARPES:MapInPlane
	string/G s_MapInPlane_WorkFld
	string/G s_MapInPlane_ARP3D_name

	SetDataFolder $(s_MapInPlane_WorkFld + "Ang2K")
	wave w2D_kxky = $(s_MapInPlane_ARP3D_name + "_kxky")
	NewDataFolder/O/S :Misc
	variable/G v_mesh_kx, v_mesh_ky
	
	Redimension/N=(v_mesh_kx, v_mesh_ky) w2D_kxky
End




///////////////////////////////////////
//// Set range of the A2K-variable slider    
///////////////////////////////////////


Function Y_MapIP_RangeSlider(s_var)
	string s_var	// = beta0, xi, xi0, delta, chi, chi0, work, hn

	string s_fld = GetDataFolder(1)
	
	SetDataFolder root:Y_ARPES:MapInPlane
	string/G s_MapInPlane_WorkFld
	string/G s_MapInPlane_ARP3D_name

	SetDataFolder $(s_MapInPlane_WorkFld + "Ang2K:Misc")
	NVAR v_val_min	= $("v_srange_min_" + s_var)
	NVAR v_val_max	= $("v_srange_max_" + s_var)

	string s_command = ""
	s_command = s_command + "Slider slider_" + s_var + " limits={" + num2str(v_val_min) + ","
	s_command = s_command + num2str(v_val_max) +  ",0}, win=P_ARPES_Map_InPlane" 
	
	Execute s_command

	SetDataFolder $s_fld
End




/////////////////////////////////////
////  Set range for the cut dimension  
/////////////////////////////////////


Function Y_MapIP_RangeCut()
	string s_fld = GetDataFolder(1)
	
	SetDataFolder root:Y_ARPES:MapInPlane
	string/G s_MapInPlane_WorkFld
	string/G s_MapInPlane_ARP3D_name

	SetDataFolder $s_MapInPlane_WorkFld
	wave w3D = $s_MapInPlane_ARP3D_name
	wave w2D_AngCut = $(s_MapInPlane_ARP3D_name + "_AngCut")
	wave w2D_EngCut = $(s_MapInPlane_ARP3D_name + "_EngCut")
	wave w_cut		= :Misc:w_cut
	NVAR v_cut_start	= :Misc:v_cut_start
	NVAR v_cut_delta	= :Misc:v_cut_delta
	NVAR v_cut_end	= :Misc:v_cut_end
	NVAR v_cut_width	= :Misc:v_cut_width

	v_cut_end = v_cut_start + v_cut_delta * (DimSize(w3D, 2) - 1)
	
	SetScale/I y  v_cut_start,  v_cut_end, "", w2D_AngCut
	SetScale/I y  v_cut_start,  v_cut_end, "", w2D_EngCut
	SetScale/I z  v_cut_start,  v_cut_end, "", w3D
	SetScale/I x  v_cut_start,  v_cut_end, "", w_cut
	w_cut = x

	//Update panel
	
	Slider slider_cut limits			= {min(v_cut_start, v_cut_end) + v_cut_width/2,	max(v_cut_start, v_cut_end) - v_cut_width/2, v_cut_delta}, win = P_ARPES_Map_InPlane
	Slider slider_width_cut limits	= {0, abs(v_cut_start - v_cut_end)/2, v_cut_delta * 2}, win = P_ARPES_Map_InPlane
	
	SetDataFolder $s_fld
End




//////////////////////////////
////  Choose/make a 3D wave
//////////////////////////////


Function Y_MapIP_PanelFolderUpdate_A2K(s_WorkFld, s_w3D_name)
	string s_WorkFld		//s_MapInPlane_WorkFld,
	string s_w3D_name		//s_MapInPlane_ARP3D_name

	string s_fld = GetDataFolder(1)
	
	SetDataFolder $s_WorkFld
	wave w3D 		= $s_w3D_name
	wave w2D_AngCut = $(s_w3D_name + "_AngCut")
	
	NewDataFolder/O/S :Ang2K
	Duplicate/O w2D_AngCut, $(s_w3D_name + "_kxky")
	
	NewDataFolder/O/S :Contour
	variable/G v_ky_Boundary_max, v_ky_Boundary_min, v_kx_Boundary_max, v_kx_Boundary_min
	NewDataFolder/O/S ::Misc
	
	Y_MapIP_A2K_InitializeVariables()
	
	variable/G v_mesh_kx = DimSize(w2D_angCut, 0)
	variable/G v_mesh_ky = DimSize(w2D_angCut, 1)
	variable/G v_ky_max, v_ky_min, v_kx_max, v_kx_min
	variable/G v_k_range_flag
	variable/G v_type
	variable/G v_DA	// 0 --> for future use
	variable/G v_hn	
	variable/G v_work	 
	variable/G v_delta, v_xi0, v_xi, v_beta0
	variable/G v_chi, v_chi0	// --> for future use
	variable/G v_srange_min_hn, v_srange_min_work, v_srange_min_delta, v_srange_min_xi0, v_srange_min_xi, v_srange_min_beta0, v_srange_min_chi, v_srange_min_chi0
	variable/G v_srange_max_hn, v_srange_max_work, v_srange_max_delta, v_srange_max_xi0, v_srange_max_xi, v_srange_max_beta0, v_srange_max_chi, v_srange_max_chi0
	
	//Re-assign global values in the panel	
	DoWindow/F P_ARPES_Map_InPlane
	CheckBox	check_DA				variable= v_DA
	Slider		slider_A2K_type			variable= v_type
	Slider		slider_A2K_kRangeFlag		variable= v_k_range_flag
	SetVariable	setvar_kRange_kx_max		value= v_kx_max
	SetVariable	setvar_kRange_kx_min		value= v_kx_min
	SetVariable	setvar_kRange_ky_max		value= v_ky_max
	SetVariable	setvar_kRange_ky_min		value= v_ky_min
	SetVariable	setvar_Mesh_kx 			value= v_mesh_kx
	SetVariable	setvar_Mesh_ky 			value= v_mesh_ky

	string s_var_list = "beta0;xi;xi0;delta;chi;chi0;work;hn"
	variable v_var_list_num = itemsinlist(s_var_list)
	variable i
	for(i=0; i<v_var_list_num; i+=1)
		string s_var = stringfromlist(i, s_var_list)

		//Re-assign global values in the panel
		Execute "Slider slider_" + s_var + " variable= v_" + s_var + ", win = P_ARPES_Map_InPlane"
		Execute "SetVariable setvar_srange_max_" + s_var + " value=  v_srange_max_" + s_var + ", win = P_ARPES_Map_InPlane"
		Execute "SetVariable setvar_srange_min_" + s_var + " value=  v_srange_min_" + s_var + ", win = P_ARPES_Map_InPlane"
		Execute "SetVariable setvar_" + s_var + " value= v_" + s_var + ", win = P_ARPES_Map_InPlane"

		//Updates the slider ranges
		Y_MapIP_RangeSlider(s_var)		
	endfor
			
	SetDataFolder $s_fld
End

	SetVariable	setvar_beta0 				value= v_beta0
	SetVariable	setvar_xi 				value= v_xi
	SetVariable	setvar_xi0 				value= v_xi0
	SetVariable	setvar_delta 				value= v_delta		
	SetVariable	setvar_chi 				value= v_chi		
	SetVariable	setvar_chi0 				value= v_chi0		
	SetVariable	setvar_work 				value= v_work
	SetVariable	setvar_hn	 			value= v_hn

	SetVariable	setvar_srange_min_beta0	value=  v_srange_min_beta0
	SetVariable	setvar_srange_min_xi		value=  v_srange_min_xi
	SetVariable	setvar_srange_min_xi0		value=  v_srange_min_xi0
	SetVariable	setvar_srange_min_delta	value=  v_srange_min_delta
	SetVariable	setvar_srange_min_chi		value=  v_srange_min_chi
	SetVariable	setvar_srange_min_chi0		value=  v_srange_min_chi0
	SetVariable	setvar_srange_min_work	value=  v_srange_min_work
	SetVariable	setvar_srange_min_hn		value=  v_srange_min_hn

	SetVariable	setvar_srange_max_beta0	value=  v_srange_max_beta0
	SetVariable	setvar_srange_max_xi		value=  v_srange_max_xi
	SetVariable	setvar_srange_max_xi0		value=  v_srange_max_xi0
	SetVariable	setvar_srange_max_delta	value=  v_srange_max_delta
	SetVariable	setvar_srange_max_chi		value=  v_srange_max_chi
	SetVariable	setvar_srange_max_chi0	value=  v_srange_max_chi0
	SetVariable	setvar_srange_max_work	value=  v_srange_max_work
	SetVariable	setvar_srange_max_hn		value=  v_srange_max_hn	
	Slider 		slider_beta0 				variable= v_beta0	
	Slider 		slider_xi	 				variable= v_xi	
	Slider 		slider_xi0 				variable= v_xi0	
	Slider 		slider_delta 				variable= v_delta
	Slider 		slider_chi 				variable= v_chi
	Slider 		slider_chi0 				variable= v_chi0		
	Slider 		slider_work				variable= v_work	
	Slider 		slider_hn 				variable= v_hn	

Function Y_MapIP_PanelFolderUpdate(s_WorkFld, s_w3D_name)
	string s_WorkFld		//s_MapInPlane_WorkFld,
	string s_w3D_name		//s_MapInPlane_ARP3D_name

	string s_fld = GetDataFolder(1)
	
	SetDataFolder $s_WorkFld
	wave w3D = $s_w3D_name								//Dim0 --> Eng;   Dim1 --> Ang;  Dim2;  Cut
	
	Duplicate/O w3D, $(s_w3D_name + "_EngAng")
	wave w2D_EngAng	= $(s_w3D_name + "_EngAng")			//Energy dispersion
	Redimension/N=(-1,-1) w2D_EngAng
	
	Duplicate/O w2D_EngAng, $(s_w3D_name + "_EngCut"), $(s_w3D_name + "_AngCut")		
	wave w2D_AngCut	= $(s_w3D_name + "_AngCut")			//Ang-Ang map
	matrixtranspose w2D_AngCut
	Redimension/N=(-1, DimSize(w3D, 2)) w2D_AngCut
	
	wave w2D_EngCut	= $(s_w3D_name + "_EngCut")			//Eng-Cut map
	Redimension/N=(-1, DimSize(w3D, 2)) w2D_EngCut			


	NewDataFolder/O/S :Misc
	variable/G v_cut_width_flg, v_cut_width, v_cut, v_ang_width_flg, v_ang_width, v_ang, v_eng_width_flg, v_eng_width, v_eng
	variable/G v_eng_start		= DimOffset(w3D, 0)
	variable/G v_eng_end		= DimOffset(w3D, 0) + DimDelta(w3D, 0) * (DimSize(w3D, 0) - 1)
	variable/G v_ang_start		= DimOffset(w3D, 1)
	variable/G v_ang_end		= DimOffset(w3D, 1) + DimDelta(w3D, 1) * (DimSize(w3D, 1) - 1)
	variable/G v_cut_start		= DimOffset(w3D, 2)
	variable/G v_cut_end		= DimOffset(w3D, 2) + DimDelta(w3D, 2) * (DimSize(w3D, 2) - 1)
	variable/G v_cut_delta		= DimDelta(w3D, 2)
	
	
	SetScale/I y  v_cut_start,  v_cut_end,"", w2D_AngCut
	SetScale/I y  v_eng_start,  v_eng_end,"", w2D_EngCut
	
	Duplicate/O w2D_EngAng, w_eng
	Redimension/N=(-1) w_eng
	w_eng = x
	
	Duplicate/O w2D_AngCut, w_cut, w_ang
	Redimension/N=(-1) w_ang
	w_ang = x
	
	MatrixTranspose w_cut
	Redimension/N=(-1) w_cut
	w_cut = x
	
	
	////Update panel	

	DoWindow/F P_ARPES_Map_InPlane
	
	Button button_MapIP_WorkingWave title	="\\JL" + s_w3D_name
	Button button_MapIP_WorkingFolder title	="\\JL" + s_WorkFld[5,strlen(s_WorkFld)-2]
	Slider slider_eng limits	=	{v_eng_start + v_eng_width/2,	v_eng_end - v_eng_width/2,	0}
	Slider slider_ang limits	=	{v_ang_start + v_ang_width/2,	v_ang_end - v_ang_width/2,	0}
	Slider slider_cut limits	=	{v_cut_start + v_cut_width/2,	v_cut_end - v_cut_width/2,	v_cut_delta}
	Slider slider_width_eng limits	=	{0,	abs(v_eng_start - v_eng_end)/3,		0},				disable= (2 - 2*v_eng_width_flg)
	Slider slider_width_ang limits	=	{0, 	abs(v_ang_start - v_ang_end)/3,		0},				disable= (2 - 2*v_ang_width_flg)
	Slider slider_width_cut limits	=	{0,	abs(v_cut_start - v_cut_end),		v_cut_delta * 2},	disable= (2 - 2*v_cut_width_flg)


	////Re-assign global values in the panel	

	CheckBox check_width_eng value = v_eng_width_flg,	variable	= v_eng_width_flg
	CheckBox check_width_ang value = v_ang_width_flg, variable	= v_ang_width_flg
	CheckBox check_width_cut value = v_cut_width_flg,	variable	= v_cut_width_flg

	SetVariable setvar_eng value	= v_eng
	SetVariable setvar_ang value	= v_ang
	SetVariable setvar_cut value	= v_cut
	SetVariable setvar_width_eng value	= v_eng_width,		disable=2
	SetVariable setvar_width_ang value	= v_ang_width,		disable=2
	SetVariable setvar_width_cut value	= v_cut_width,		disable=2
	Slider slider_eng		variable	= v_eng
	Slider slider_ang		variable	= v_ang
	Slider slider_cut		variable	= v_cut
	Slider slider_width_eng	variable	= v_eng_width
	Slider slider_width_ang	variable	= v_ang_width
	Slider slider_width_cut	variable	= v_cut_width

	SetVariable setvar_cut_start value	= v_cut_start
	SetVariable setvar_cut_delta value	= v_cut_delta
	SetVariable setvar_cut_end value	= v_cut_end
		
	SetDataFolder $s_fld
End


Function Y_MapIP_A2K_PanelDA()
	string s_fld = GetDataFolder(1)
	
	SetDataFolder root:Y_ARPES:MapInPlane
	string/G s_MapInPlane_WorkFld
	string/G s_MapInPlane_ARP3D_name

	SetDataFolder $(s_MapInPlane_WorkFld) + "Ang2K:Misc:"
	variable/G v_DA
	
	DoWindow/F P_ARPES_Map_InPlane

	Execute "SetVariable setvar_beta0 disable=2 * " 				+ num2str(v_DA) 
	Execute "SetVariable setvar_srange_min_beta0 disable=2 * " 	+ num2str(v_DA)
	Execute "SetVariable setvar_srange_max_beta0 disable=2 * " 	+ num2str(v_DA)
	Execute "Slider slider_beta0 disable=2 * "					+ num2str(v_DA)

	Execute "SetVariable setvar_chi disable=2-2* " 				+ num2str(v_DA) 
	Execute "SetVariable setvar_srange_min_chi disable= 2-2*" 		+ num2str(v_DA)
	Execute "SetVariable setvar_srange_max_chi disable= 2-2*" 	+ num2str(v_DA)
	Execute "Slider slider_chi disable= 2-2*"						+ num2str(v_DA)

	Execute "SetVariable setvar_chi0 disable=2-2*" 				+ num2str(v_DA)
	Execute "SetVariable setvar_srange_min_chi0 disable=2-2*"		+ num2str(v_DA)
	Execute "SetVariable setvar_srange_max_chi0 disable=2-2*"	+ num2str(v_DA)
	Execute "Slider slider_chi0 disable=2-2*"					+ num2str(v_DA)
	
	SetDataFolder $s_fld
End


Function Y_PopMenuProc_MapIP(ctrlName,popNum,popStr) : PopupMenuControl
	String ctrlName
	Variable popNum
	String popStr
	
	string s_fld = GetDataFolder(1)
	
	SetDataFolder root:Y_ARPES:MapInPlane
	variable/G v_MapInPlane_choice 
	
	if(stringmatch(ctrlName, "popup_MapIP_MakeChoose3D"))

		variable/G v_MapInPlane_choice = popNum	//1 --> A 3D wave chosen in browser; 2 --> 2D waves in browser; 3.--> Folders in browser

		DoWindow/F P_ARPES_Map_InPlane
		SetVariable setvar_MapIP_ItemNo disable= (2*(popNum != 3))		//Diables only when "Folders selected in browser" is chosen
		Execute "CreateBrowser"

	endif
		
	SetDataFolder $s_fld
End


Function Y_ButtonProc_MapIP(ctrlName) : ButtonControl
	string ctrlName

	string s_fld = GetDataFolder(1)

	NewDataFolder/O/S root:Y_ARPES
	NewDataFolder/O/S :MapInPlane
	variable/G v_MapInPlane_choice, v_MapInPlane_index
	string/G s_MapInPlane_WorkFld, s_MapInPlane_ARP3D_name
	NewDataFolder/O/S :ARP3D
	NewDataFolder/O :Misc	

	if(stringmatch(ctrlName, "button_MapIP_Make3D"))
	
		string datalist	
		if(v_MapInPlane_Choice == 1)		//--> 1: A 3D wave chosen in browser

			datalist 	= Y_BrowserSelect_Waves_MapIP(";")
			wave w3D = $stringfromlist(0, datalist)

			if(WaveDims(w3D) != 3)
				DoAlert/T="FATAL ERROR IN SELECTION" 0, "\r\r\rThe (1st) selected wave was not a 3D wave"
				DoWindow/F P_ARPES_Map_InPlane
				Button button_MapIP_WorkingWave title	="_none_"
				Button button_MapIP_WorkingFolder title	="_none_"
				Abort
			endif

		else
			if(v_MapInPlane_Choice == 2)	//--> 3: 2D waves chosen in browser
				datalist 	= Y_BrowserSelect_Waves_MapIP(";")
			elseif(v_MapInPlane_Choice == 3)	//--> 3: Folders selected in browser
				datalist 	= Y_BrowserSelect_Folders_MapIP(";", v_MapInPlane_index)
			endif
			
			wave w3D_init = $stringfromlist(0, datalist)
			if(WaveDims(w3D_init) != 2)
				DoAlert/T="FATAL ERROR IN SELECTION" 0, "\r\r\rThe (1st) selected wave was not a 2D wave"
				DoWindow/F P_ARPES_Map_InPlane
				Button button_MapIP_WorkingWave title	="_none_"
				Button button_MapIP_WorkingFolder title	="_none_"
				Abort
			endif
			Duplicate/O w3D_init, $("wARP3D")

			wave w3D = $("wARP3D")
			variable v_NumCuts = itemsinlist(datalist)
			Redimension/N=(-1,-1,v_NumCuts) w3D
			variable i
			for(i=0; i<v_NumCuts; i+=1)
				wave w2D = $stringfromlist(i, datalist)
				w3D[][][i] = w2D[p][q]
			endfor
			
			note w3D, 
			note w3D, ""	
			note w3D, "[Y_MapInPlane 2D waves --> 3D]"
			note w3D, datalist
			
		endif

		s_MapInPlane_WorkFld		= GetWavesDataFolder(w3D, 1)
		s_MapInPlane_ARP3D_name	= nameofwave(w3D)
		
		Y_MapIP_PanelFolderUpdate(s_MapInPlane_WorkFld, s_MapInPlane_ARP3D_name)
		Y_MapIP_PanelFolderUpdate_A2K(s_MapInPlane_WorkFld, s_MapInPlane_ARP3D_name)
		Y_MapIP_A2K_BoundaryMake()
		Y_MapIP_MakePanel_Init(0)
		Y_MapIP_A2K_PanelDA()
		Y_MapIP_A2K_ChangeType()
		
		SetDataFolder $s_fld
	endif
End

Function/S Y_BrowserSelect_Folders_MapIP(SepStr, v_index)
	string SepStr
	variable v_index
	
	string s_fld = GetDataFolder(1)

	string datalist=""
	string s_not_fld = ""
	variable i = 0

	if(stringmatch(GetBrowserSelection(0), "") ==1 )	//When nothing is selected...
		DoAlert/T="FATAL ERROR IN SELECTION" 0, "\r\r\rSelect Folders (not waves) from DataBrowser"
		Abort
	endif	

	do
		if(DataFolderExists(GetBrowserSelection(i)))	//If the i'th object is indeed a folder
			SetDataFolder $GetBrowserSelection(i)
			datalist += GetBrowserSelection(i) + ":" + GetIndexedObjName(":", 1, v_index-1)+ SepStr
		else										//If the i'th object is not a folder
			s_not_fld += num2str(i+1) + ", " 
		endif
		i+=1
	while(stringmatch(GetBrowserSelection(i), "")==0)

	if(stringmatch(datalist, "") == 1)					//When all the selections were not folders...
		DoAlert/T="FATAL ERROR IN SELECTION" 0, "\r\r\rSelect Folders (NOT WAVES!!) from DataBrowser"
		Abort
	elseif(strlen(s_not_fld))		//When items other than folders were selected...
		string s_not_fld_alert = "" 
		s_not_fld_alert += "\rThe following object(s)\r\r    " 
		s_not_fld_alert += s_not_fld[0,strlen(s_not_fld) - 3] + "\r\r"
		s_not_fld_alert += "of the selection in the browser was (were) not folders.\r\r"
		DoAlert/T="Although not fatal..." 0, s_not_fld_alert
	endif
	
	SetDataFolder $s_fld
	return datalist
End


Function/S Y_BrowserSelect_Waves_MapIP(SepStr)
	string SepStr
	
	string s_fld = GetDataFolder(1)

	string datalist=""
	string s_not_wave = ""
	variable i = 0

	if(stringmatch(GetBrowserSelection(0), "") ==1 )	//When nothing is selected...
		DoAlert/T="FATAL ERROR IN SELECTION" 0, "\r\r\rSelect Waves (not folders) from DataBrowser"
		Abort
	endif	

	do
		if(WaveExists($GetBrowserSelection(i)))	//If the i'th object is indeed a wave
			datalist += GetBrowserSelection(i) + SepStr
		else										//If the i'th object is not a wave
			s_not_wave += num2str(i+1) + ", " 
		endif
		
		i+=1
	while(stringmatch(GetBrowserSelection(i), "")==0)

	if(stringmatch(datalist, "") == 1)					//When all the selections were not waves...
		DoAlert/T="FATAL ERROR IN SELECTION" 0, "\r\r\rSelect Waves (NOT FOLDERS!!) from DataBrowser"
		Abort
	elseif(strlen(s_not_wave))		//When items other than folders were selected...
		string s_not_wave_alert = "" 
		s_not_wave_alert += "\rThe following object(s)\r\r    " 
		s_not_wave_alert += s_not_wave[0,strlen(s_not_wave) - 3] + "\r\r"
		s_not_wave_alert += "of the selection in the browser was (were) not waves.\r\r"
		DoAlert/T="Although not fatal..." 0, s_not_wave_alert
	endif
	
	SetDataFolder $s_fld
	return datalist
End




//////////////////////////////////////////////
////  Open only the working folder in the Browser  
//////////////////////////////////////////////


Function Y_ButtonProc_MapIP_OpenFld(ctrlName) : ButtonControl
	String ctrlName

	//If "_none_" is displayed in the button, ABORT.
	ControlInfo/W=P_ARPES_Map_InPlane button_MapIP_WorkingWave
		string s_ButtonInfo1 = S_recreation
	ControlInfo/W=P_ARPES_Map_InPlane button_MapIP_WorkingFolder
		string s_ButtonInfo2 = S_recreation
	if((stringmatch(s_ButtonInfo1, "*_none_*") || stringmatch(s_ButtonInfo1, "*_none_*")))	
		DoAlert/T="FATAL ERROR" 0, "\r\r3D wave is not properly chosen/made."
		Abort
	endif

	string s_fld = GetDataFolder(1)

	NewDataFolder/O/S root:Y_ARPES
	NewDataFolder/O/S :MapInPlane
	string/G s_MapInPlane_WorkFld
	
	Execute "CreateBrowser"						//Brings Browser to top
	Execute "ModifyBrowser collapseAll, expand=0"	//Collapse all folders, then open the root: folder
	
	string n_Fld, s_line
	string n_FldPass = "root:"
	variable v_line
	variable npnts = ItemsInList(s_MapInPlane_WorkFld, ":")
	
	variable i
	for(i=1; i<(npnts); i=i+1)						//Open subfolders
		n_Fld = stringfromlist(i, s_MapInPlane_WorkFld, ":")
		
		n_FldPass = n_FldPass + n_Fld + ":"
		v_line = GetBrowserLine(n_FldPass)
		s_line =  "ModifyBrowser expand=" + num2str(v_line)
		Execute s_line
	endfor

	SetDataFolder $s_fld
End




/////////////////////////////
////  Recreation of the panel
/////////////////////////////


Function Y_MapInPlane_Start()
	if(strlen(WinList("P_ARPES_Map_InPlane", ";", "WIN:64")))
		DoWindow/F P_ARPES_Map_InPlane
	else
		Y_MapIP_MakeAllVariables()
		Y_MapIP_MakePanel()
		Y_MapIP_MakePanel_Init(2)
	endif
End


Function Y_MapIP_MakeAllVariables()
	string s_fld = GetDataFolder(1)
	
	NewDataFolder/O/S root:Y_ARPES
	NewDataFolder/O/S :MapInPlane
	variable/G v_MapInPlane_choice		= 1			//1: FromTopWindowImage; 2: FromSelectedFoldersInBrowser; 3: FromSelectedWavesInBrowser
	variable/G v_MapInPlane_index		= 2			//"v_ImgProc_index"th wave in the folder; Valid only when v_ImageProc_choice == 0
	NewDataFolder/O/S :ARP3D
	NewDataFolder/O/S :Misc
	variable/G v_cut_delta, v_cut_end, v_cut_start, v_ang_end, v_ang_start, v_eng_end, v_eng_start, v_eng, v_eng_width, v_eng_width_flg
	variable/G v_ang, v_ang_width, v_ang_width_flg, v_cut, v_cut_width, v_cut_width_flg
	
	NewDataFolder/O/S ::Ang2K
	NewDataFolder/O/S :Contour
	variable/G v_ky_Boundary_min, v_ky_Boundary_max, v_kx_Boundary_min, v_kx_Boundary_max
	NewDataFolder/O/S ::Misc	
	
	// v_type	// 0 --> slit parallel to manipulator axis; 1 --> slit perpendicular to manipulator axis

	Y_MapIP_A2K_InitializeVariables()

	SetDataFolder $s_fld
End


Function Y_MapIP_A2K_InitializeVariables()
	string s_varlist_name	= ""
	string s_varlist_var	= ""
	s_varlist_name	+= "v_beta0;v_xi;v_xi0;v_delta;v_chi;v_chi0;v_work;v_hn;v_type;"
	s_varlist_name	+= "v_srange_min_beta0;v_srange_max_beta0;v_srange_min_xi;v_srange_max_xi;v_srange_min_xi0;v_srange_max_xi0;v_srange_min_delta;v_srange_max_delta;v_srange_min_chi;v_srange_max_chi;v_srange_min_chi0;v_srange_max_chi0;v_srange_min_work;v_srange_max_work;v_srange_min_hn;v_srange_max_hn;"
	s_varlist_name	+= "v_mesh_kx;v_mesh_ky;v_k_range_flag;v_kx_min;v_kx_max;v_ky_min;v_ky_max;v_DA"
	s_varlist_var		+= "0;0;0;0;0;0;4.18;50;0;"
	s_varlist_var		+= "-80;80;-80;80;-80;80;-200;200;-80;80;-80;80;3.5;5.5;5.5;50;"
	s_varlist_var		+= "10;10;0;-1;1;-1;1;0;"
	
	variable v_varlist_num = itemsinlist(s_varlist_name)
	variable i
	for(i=0; i<v_varlist_num; i+=1)
		NVAR/Z v_var = $stringfromlist(i, s_varlist_name)
		if(!NVAR_Exists(v_var))
			variable/G $stringfromlist(i, s_varlist_name)
			NVAR v_var_exist = $stringfromlist(i, s_varlist_name)
			v_var_exist = str2num(stringfromlist(i, s_varlist_var))
		endif
	endfor
End


Function Y_MapIP_MakePanel_Init(v_Flag_enable)
	variable v_Flag_enable	// 0 --> enable; 2 --> disable
	DoWindow/F P_ARPES_Map_InPlane
	
	string s_list =""
	variable v_list_num
	variable i

	s_list = "eng;ang;cut;A2K_type;A2K_kRangeFlag;"
	s_list += "beta0;xi;xi0;delta;chi;chi0;work;hn"
	v_list_num = itemsinlist(s_list)
	for(i=0; i<v_list_num; i+=1)
		Execute "Slider slider_" + stringfromlist(i, s_list) + " disable = " + num2str(v_Flag_enable)
	endfor
	
	s_list = "Update_A2K_3DMap;"
	v_list_num = itemsinlist(s_list)
	for(i=0; i<v_list_num; i+=1)
		Execute "Button button_" + stringfromlist(i, s_list) + " disable = " + num2str(v_Flag_enable)
	endfor

	s_list = "width_eng;width_ang;width_cut;DA"
	v_list_num = itemsinlist(s_list)
	for(i=0; i<v_list_num; i+=1)
		Execute "CheckBox check_" + stringfromlist(i, s_list) + " disable = " + num2str(v_Flag_enable)
	endfor
	
	s_list = "eng;ang;cut;cut_start;cut_delta;"
	s_list += "kRange_ky_max;kRange_ky_min;kRange_kx_max;kRange_kx_min;"
	s_list += "Mesh_kx;Mesh_ky;"
	s_list += "beta0;srange_min_beta0;srange_max_beta0;"
	s_list += "xi;srange_min_xi;srange_max_xi;"	
	s_list += "xi0;srange_min_xi0;srange_max_xi0;"
	s_list += "delta;srange_min_delta;srange_max_delta;"
	s_list += "chi;srange_min_chi;srange_max_chi;"
	s_list += "chi0;srange_min_chi0;srange_max_chi0;"
	s_list += "work;srange_min_work;srange_max_work;"
	s_list += "hn;srange_min_hn;srange_max_hn;"
	v_list_num = itemsinlist(s_list)

	for(i=0; i<v_list_num; i+=1)
		Execute "SetVariable setvar_" + stringfromlist(i, s_list) + " disable = " + num2str(v_Flag_enable)
	endfor
End


Function Y_MapIP_MakePanel() : Panel
	PauseUpdate; Silent 1		// building window...

	NewPanel /N=P_ARPES_Map_InPlane   /W=(700,71,1011,782)  /K=1  as "P_ARPES_Map_InPlane"		//Note the flag /K=1 (no dialog when killing)

	ModifyPanel cbRGB=(48896,65280,48896)
	SetDrawLayer UserBack
	SetDrawEnv fillfgc= (56576,56576,56576)
	DrawRect 115,449,141,485
	SetDrawEnv xcoord= rel,ycoord= abs,linethick= 0,linefgc= (65280,0,0)
	DrawRect 0,242,1,327
	SetDrawEnv xcoord= rel,ycoord= abs,linethick= 0,linefgc= (65280,0,0),fillfgc= (65280,48896,48896)
	DrawRect 0,0,1,85
	SetDrawEnv linethick= 0
	DrawPoly 186.745928338762,115.75,0.697068,0.75,{360,72,360,144,405,162,477,162,477,90,432,72,360,72}
	SetDrawEnv fillfgc= (48896,65280,65280)
	DrawPoly 255.755700325733,163,0.697068,0.75,{459,135,459,162,414,144,414,117,459,135}
	SetDrawEnv fillfgc= (65280,65280,48896)
	DrawPoly 199.293159609121,155.5,0.522801,0.28125,{405,90,405,162,477,162,477,90,405,90}
	SetDrawEnv fillfgc= (65280,65280,48896)
	DrawPoly 236.934853420196,121.75,0.174267,0.375,{405,90,405,162,477,162,477,90,405,90}
	DrawPoly 93.4625407166122,115.75,0.697068,0.75,{288,72,288,144,333,162,333,90,288,72}
	DrawPoly 68.3680781758957,115.75,0.697068,0.75,{288,72,288,144,333,162,333,90,288,72}
	DrawPoly 55.8208469055375,115.75,0.697068,0.75,{288,72,288,144,333,162,333,90,288,72}
	DrawPoly 43.2736156351792,115.75,0.697068,0.75,{288,72,288,144,333,162,333,90,288,72}
	SetDrawEnv fillpat= 0
	DrawPoly 186.745928338762,115.75,0.697068,0.75,{288,72,288,144,333,162,333,90,288,72}
	SetDrawEnv fillpat= 0
	DrawPoly 218.114006514658,129.25,0.697068,0.75,{405,90,477,90,432,72,360,72,405,90}
	SetDrawEnv fillfgc= (65280,54528,48896)
	DrawPoly 218.114006514658,156.25,0.697068,0.75,{405,90,477,90,432,72,360,72,405,90}
	SetDrawEnv fillfgc= (48896,65280,65280)
	DrawPoly 224.387622149838,115.75,0.697068,0.75,{405,72,405,108,450,126,450,90,405,72}
	SetDrawEnv arrow= 1,arrowlen= 8
	DrawLine 186,109,237,109
	SetDrawEnv arrow= 1,arrowlen= 8
	DrawLine 41,177,72,190
	SetDrawEnv arrow= 1,arrowlen= 8
	DrawLine 37,168,37,120
	SetDrawEnv arrow= 1,arrowlen= 8
	DrawLine 185,177,216,190
	SetDrawEnv arrow= 1,arrowlen= 8
	DrawLine 181,168,181,120
	SetDrawEnv fillfgc= (65280,65280,48896)
	DrawPoly 199.293159609121,121.75,0.522801,0.375,{405,90,405,162,477,162,477,90,405,90}
	DrawLine 218,187,218,197
	DrawLine 268,187,268,197
	SetDrawEnv fillpat= 0
	DrawPoly 218.114006514658,129.25,0.697068,0.75,{405,90,405,162,477,162,477,90,405,90}
	SetDrawEnv fillpat= 0
	DrawPoly 217.416938110749,129.25,0.697068,0.75,{405,90,477,90,432,72,360,72,405,90}
	SetDrawEnv fillfgc= (65280,54528,48896)
	DrawPoly 186.745928338762,142.75,0.697068,0.75,{360,108,360,117,405,135,405,126,360,108}
	SetDrawEnv fillfgc= (65280,54528,48896)
	DrawPoly 218.114006514658,156.25,0.697068,0.75,{405,126,405,135,477,135,477,126,405,126}
	DrawText 181,213,"Start"
	DrawText 288,213,"End"
	SetDrawEnv translate= 45,132,rotate= 25.065,rsabout
	SetDrawEnv translate= 4.39946,127.805,rotate= -25.065,scale= 0.728571,0.860465,rsabout
	SetDrawEnv translate= 4.44916,114.952,rotate= -25.065,scale= 1.07547,0.918919,rsabout
	SetDrawEnv translate= 68.0316,85.2153,rotate= -25.065,scale= 1,1.02941,rsabout
	SetDrawEnv translate= 68.0316,85.2153,rotate= -25.065,scale= 1.07273,1.05714,rsabout
	SetDrawEnv translate= 4.44916,114.952,rotate= -25.065,scale= 0.881356,1,rsabout
	SetDrawEnv translate= 45.1026,119.033,rotate= -4.10184,rsabout
	SetDrawEnv translate= 124.694,156.989,rotate= -20.9632,scale= 1.29091,1.18182,rsabout
	SetDrawEnv translate= 69.0869,219.364,rotate= -20.9632,scale= 1.10145,1.15385,rsabout
	SetDrawEnv translate= 123.628,157.398,rotate= -20.9632,scale= 0.935065,1.06522,rsabout
	SetDrawEnv translate= 97.3147,188.042,rotate= -1.94954,rsabout
	DrawText 60.5251445469633,195.920517342092,"Emission angle"
	DrawText 42,113,"0"
	DrawText 55,113,"1"
	DrawText 67,113,"2"
	DrawText 91,113,"(N-1)"
	DrawText 20,113,"Cut"
	DrawText 76,113,"..."
	SetDrawEnv fname= "Arial",fstyle= 1
	DrawText 6,20,"Choose/make a 3D wave"
	SetDrawEnv xcoord= rel,ycoord= abs,linethick= 3
	DrawLine 0,85,1,85
	SetDrawEnv fname= "Arial",fstyle= 1
	DrawText 73,63,"3D wave:"
	SetDrawEnv fname= "Arial",fstyle= 1
	DrawText 86,80,"Folder:"
	DrawText 197,259,"Integral width"
	DrawText 224,213,"Step size"
	SetDrawEnv fillpat= 0
	DrawBezier 221,189,0.102041,0.148148,{27,405,27,405,36,432,72,432,108,432,117,405,117,405}
	SetDrawEnv fillpat= 0
	DrawBezier 230.183673469388,189,0.102041,0.148148,{27,405,27,405,36,432,72,432,108,432,117,405,117,405}
	SetDrawEnv fillpat= 0
	DrawBezier 239.367346938776,189,0.102041,0.148148,{27,405,27,405,36,432,72,432,108,432,117,405,117,405}
	SetDrawEnv fillpat= 0
	DrawBezier 256.816326530612,189,0.102041,0.148148,{27,405,27,405,36,432,72,432,108,432,117,405,117,405}
	DrawLine 206,202,218,197
	DrawLine 268,197,285,202
	SetDrawEnv arrow= 1,arrowlen= 5,arrowfat= 0.7
	DrawLine 278,145,278,155
	SetDrawEnv arrow= 2,arrowlen= 5,arrowfat= 0.7
	DrawLine 278,165,278,175
	DrawLine 274,156,283,156
	DrawLine 274,163,283,163
	DrawText 274,144,"Width"
	DrawText 48,259,"2D slice"
	SetDrawEnv linethick= 0,fillfgc= (21760,21760,21760)
	DrawPoly 137,155.5,0.425926,0.611111,{72,504,108,504,108,513,126,495,108,477,108,486,72,486,72,504}
	DrawText 250,195,"..."
	SetDrawEnv textrot= 90
	DrawText 163,163,"Energy"
	SetDrawEnv translate= 136,203,rotate= 123.105,rsabout
	SetDrawEnv translate= 113.604,196.115,rotate= -123.105,scale= 1.08108,1.29032,rsabout
	SetDrawEnv translate= 115.119,198.439,rotate= -123.105,scale= 1.05405,0.83871,rsabout
	SetDrawEnv translate= 136.795,206.306,rotate= -2.98428,rsabout
	SetDrawEnv translate= 136.795,206.306,rotate= -2.98428,rsabout
	SetDrawEnv translate= 110.429,187.278,rotate= -9.03297,rsabout
	SetDrawEnv textrot= 90
	DrawText 95.5438118152673,179.284682109076,"Angle"
	DrawText 202,108,"Cut"
	SetDrawEnv textrot= 90
	DrawText 19,163,"Energy"
	SetDrawEnv fname= "Arial",fstyle= 1
	DrawText 7,352,"Angle space  -->  k space"
	SetDrawEnv fsize= 18
	DrawText 163,232,"\\JC\\f01\\F'Symbol'b "
	SetDrawEnv fsize= 18
	DrawText 20,230,"\\JC\\f01\\F'Symbol'a\\]0 // slit"
	SetDrawEnv fsize= 18
	DrawText 5,151,"\\JC\\f01\\F'Symbol'w "
	SetDrawEnv xcoord= rel,ycoord= abs,linethick= 3
	DrawLine 0,327,1,327
	DrawText 157,511,"min"
	DrawText 281,511,"max"
	SetDrawEnv gstart
	SetDrawEnv linethick= 0
	DrawPoly 252.400997946232,419.721873061861,0.515667,0.672043,{287.722685092887,495.522147116049,223.335469539473,495.092310184012,209.50175496149,442.163383033781,269.277732322243,424.950244249075,287.722685092887,495.522147116049}
	SetDrawEnv fillfgc= (48896,65280,48896)
	DrawBezier 252.400997946232,419.721873061861,0.515667,0.672043,{287.722685092887,495.522147116049,265.373847912488,480.276785949687,256.151371527166,444.9908345162,269.277732322243,424.950244249075}
	DrawBezier 219.198657125855,419.433004155922,0.257833,0.504032,{104.670939078946,504.12308024535,80.0164828022167,488.07851422886,66.1827682242335,452.792562795373,77.00350992298,433.551177378375}
	DrawLine 252.195117581284,419.627480827203,218.955108346688,419.308441512967
	DrawLine 211.314110670337,384.788662334792,242.785737712326,371.88375275239
	SetDrawEnv gstop
	DrawLine 42,367,42,429
	DrawLine 103,429,42,429
	DrawRect 51,375,96,420
	DrawLine 200,367,200,429
	DrawLine 263,429,200,429
	SetDrawEnv gstart
	SetDrawEnv linefgc= (0,52224,0)
	DrawLine 216,370,216,424
	SetDrawEnv linefgc= (0,52224,0)
	DrawLine 226,370,226,424
	SetDrawEnv linefgc= (0,52224,0)
	DrawLine 236,370,236,424
	SetDrawEnv linefgc= (0,52224,0)
	DrawLine 247,370,247,424
	SetDrawEnv linefgc= (0,52224,0)
	DrawLine 207,381,256,381
	SetDrawEnv linefgc= (0,52224,0)
	DrawLine 207,392,256,392
	SetDrawEnv linefgc= (0,52224,0)
	DrawLine 207,402,256,402
	SetDrawEnv linefgc= (0,52224,0)
	DrawLine 207,413,256,413
	SetDrawEnv linefgc= (0,52224,0),fillpat= 0
	DrawRect 207,370,256,424
	SetDrawEnv gstop
	DrawLine 263,429,200,429
	SetDrawEnv gstart
	SetDrawEnv fillpat= 0
	DrawBezier 252.400997946232,419.721873061861,0.515667,0.672043,{287.722685092887,495.522147116049,265.373847912488,480.276785949687,256.151371527166,444.9908345162,269.277732322243,424.950244249075}
	SetDrawEnv fillpat= 0
	DrawBezier 219.198657125855,419.433004155922,0.257833,0.504032,{104.670939078946,504.12308024535,80.0164828022167,488.07851422886,66.1827682242335,452.792562795373,77.00350992298,433.551177378375}
	DrawLine 252.195117581284,419.627480827203,218.955108346688,419.308441512967
	DrawLine 211.314110670337,384.788662334792,242.785737712326,371.88375275239
	SetDrawEnv gstop
	DrawLine 207,429,207,436
	DrawLine 256,429,256,436
	DrawLine 200,424,193,424
	DrawLine 200,370,193,370
	DrawText 266,437,"\\JC\\f02kx "
	DrawText 194,364,"\\JC\\f02ky"
	SetDrawEnv linefgc= (30464,30464,30464),fillpat= 0
	DrawBezier 209,367,0.5,0.2,{288,648,288,648,306,621,333,621,360,621,378,648,378,648}
	SetDrawEnv linefgc= (30464,30464,30464),fillpat= 0
	DrawBezier 259.138336553525,372.937523232683,0.4,0.362963,{322.276673107051,566.907049870664,322.276673107051,566.907049870664,340.163214066145,594.07617238441,339.993558528593,634.575372848701,339.82390299104,675.074573312992,321.711154648543,701.904384751633,321.711154648543,701.904384751633}
	SetDrawEnv fname= "Arial",textrgb= (0,52224,0)
	DrawText 257,364,"\\f01Mesh"
	DrawText 29,465,"slit //"
	DrawText 30,484,"slit \\F'Symbol'^"
	DrawText 58,465,"rotary axis"
	DrawText 58,484,"rotary axis"
	DrawText 105,437,"\\JC\\F'Symbol'a\\]0 // slit"
	DrawText 37,367,"\\JC\\F'Symbol'b"
	SetDrawEnv dash= 2
	DrawLine 104,398,35,398
	DrawText 21,406,"\\F'Symbol'b\\]0\\B0 "
	SetDrawEnv arrow= 1,arrowlen= 5,arrowfat= 0.8,fillpat= 0
	DrawBezier 87,390.941176470588,2.18182,1.29412,{45,374,72,364,91,367,111,381}
	SetDrawEnv fname= "Arial",textrgb= (40960,22272,59904)
	DrawText 260,473,"\\f01Auto"
	SetDrawEnv fname= "Arial",textrgb= (40960,22272,59904)
	DrawText 164,473,"\\f01Manual"
	SetDrawEnv fname= "Arial",textrgb= (40960,22272,59904)
	DrawText 164,487,"\\f03k\\f01 range "
	SetDrawEnv fname= "Arial",textrgb= (40960,22272,59904)
	DrawText 252,487,"\\f03k\\f01 range "
	SetDrawEnv fillpat= 0
	DrawRect 6,445,144,489
	DrawText 120,466,"DA"
	SetDrawEnv fname= "Arial"
	DrawText 5,281,"\\f01Energy"
	SetDrawEnv fname= "Arial"
	DrawText 7,300,"\\f01Cut (\\F'Symbol'b\\]0\\f01)"
	SetDrawEnv fname= "Arial"
	DrawText 5,319,"\\f01Ang (\\F'Symbol'a\\]0\\f01)"
	SetVariable setvar_hn,pos={8,513},size={144,18},proc=Y_SetVarProc_MapIP_A2K,title=" Photon energy"
	SetVariable setvar_hn,font="Times New Roman"
	SetVariable setvar_hn,limits={-inf,inf,0},value= root:Y_ARPES:MapInPlane:ARP3D:Ang2K:Misc:v_hn
	SetVariable setvar_srange_min_hn,pos={155,513},size={24,18},proc=Y_SetVarProc_MapIP_SetRange,title=" "
	SetVariable setvar_srange_min_hn,font="Times New Roman"
	SetVariable setvar_srange_min_hn,limits={-inf,inf,0},value= root:Y_ARPES:MapInPlane:ARP3D:Ang2K:Misc:v_srange_min_hn
	Slider slider_hn,pos={181,513},size={96,19},proc=Y_SliderProc_MapIP_A2K
	Slider slider_hn,limits={5.5,50,0},variable= root:Y_ARPES:MapInPlane:ARP3D:Ang2K:Misc:v_hn,vert= 0,ticks= 0
	SetVariable setvar_srange_max_hn,pos={279,513},size={24,18},proc=Y_SetVarProc_MapIP_SetRange,title=" "
	SetVariable setvar_srange_max_hn,font="Times New Roman"
	SetVariable setvar_srange_max_hn,limits={-inf,inf,0},value= root:Y_ARPES:MapInPlane:ARP3D:Ang2K:Misc:v_srange_max_hn
	SetVariable setvar_work,pos={8,533},size={144,18},proc=Y_SetVarProc_MapIP_A2K,title=" Work function"
	SetVariable setvar_work,font="Times New Roman"
	SetVariable setvar_work,limits={-inf,inf,0},value= root:Y_ARPES:MapInPlane:ARP3D:Ang2K:Misc:v_work
	SetVariable setvar_srange_min_work,pos={155,533},size={24,18},proc=Y_SetVarProc_MapIP_SetRange,title=" "
	SetVariable setvar_srange_min_work,font="Times New Roman"
	SetVariable setvar_srange_min_work,limits={-inf,inf,0},value= root:Y_ARPES:MapInPlane:ARP3D:Ang2K:Misc:v_srange_min_work
	Slider slider_work,pos={181,532},size={96,19},proc=Y_SliderProc_MapIP_A2K
	Slider slider_work,limits={3.5,5.5,0},variable= root:Y_ARPES:MapInPlane:ARP3D:Ang2K:Misc:v_work,vert= 0,ticks= 0
	SetVariable setvar_srange_max_work,pos={279,533},size={24,18},proc=Y_SetVarProc_MapIP_SetRange,title=" "
	SetVariable setvar_srange_max_work,font="Times New Roman"
	SetVariable setvar_srange_max_work,limits={-inf,inf,0},value= root:Y_ARPES:MapInPlane:ARP3D:Ang2K:Misc:v_srange_max_work
	SetVariable setvar_xi,pos={8,553},size={144,18},proc=Y_SetVarProc_MapIP_A2K,title=" \\F'Symbol'x\\K(48896,65280,48896)\\Z070\\M\\K(0,0,0)\\F'Times New Roman' : Tilt angle"
	SetVariable setvar_xi,font="Times New Roman"
	SetVariable setvar_xi,limits={-inf,inf,0},value= root:Y_ARPES:MapInPlane:ARP3D:Ang2K:Misc:v_xi
	SetVariable setvar_srange_min_xi,pos={155,553},size={24,18},proc=Y_SetVarProc_MapIP_SetRange,title=" "
	SetVariable setvar_srange_min_xi,font="Times New Roman"
	SetVariable setvar_srange_min_xi,limits={-inf,inf,0},value= root:Y_ARPES:MapInPlane:ARP3D:Ang2K:Misc:v_srange_min_xi
	Slider slider_xi,pos={181,553},size={96,19},proc=Y_SliderProc_MapIP_A2K
	Slider slider_xi,limits={-80,80,0},variable= root:Y_ARPES:MapInPlane:ARP3D:Ang2K:Misc:v_xi,vert= 0,ticks= 0
	SetVariable setvar_srange_max_xi,pos={279,553},size={24,18},proc=Y_SetVarProc_MapIP_SetRange,title=" "
	SetVariable setvar_srange_max_xi,font="Times New Roman"
	SetVariable setvar_srange_max_xi,limits={-inf,inf,0},value= root:Y_ARPES:MapInPlane:ARP3D:Ang2K:Misc:v_srange_max_xi
	SetVariable setvar_xi0,pos={8,573},size={144,18},proc=Y_SetVarProc_MapIP_A2K,title=" \\F'Symbol'x\\F'Times New Roman'\\Z070\\M : Tilt origin"
	SetVariable setvar_xi0,font="Times New Roman"
	SetVariable setvar_xi0,limits={-inf,inf,0},value= root:Y_ARPES:MapInPlane:ARP3D:Ang2K:Misc:v_xi0
	SetVariable setvar_srange_min_xi0,pos={155,573},size={24,18},proc=Y_SetVarProc_MapIP_SetRange,title=" "
	SetVariable setvar_srange_min_xi0,font="Times New Roman"
	SetVariable setvar_srange_min_xi0,limits={-inf,inf,0},value= root:Y_ARPES:MapInPlane:ARP3D:Ang2K:Misc:v_srange_min_xi0
	Slider slider_xi0,pos={181,573},size={96,19},proc=Y_SliderProc_MapIP_A2K
	Slider slider_xi0,limits={-80,80,0},variable= root:Y_ARPES:MapInPlane:ARP3D:Ang2K:Misc:v_xi0,vert= 0,ticks= 0
	SetVariable setvar_srange_max_xi0,pos={279,573},size={24,18},proc=Y_SetVarProc_MapIP_SetRange,title=" "
	SetVariable setvar_srange_max_xi0,font="Times New Roman"
	SetVariable setvar_srange_max_xi0,limits={-inf,inf,0},value= root:Y_ARPES:MapInPlane:ARP3D:Ang2K:Misc:v_srange_max_xi0
	SetVariable setvar_delta,pos={8,593},size={144,18},proc=Y_SetVarProc_MapIP_A2K,title=" \\F'Symbol'd\\F'Times New Roman'\\K(48896,65280,48896)\\Z070\\K(0,0,0)\\]0 : In-plane rotation"
	SetVariable setvar_delta,font="Times New Roman"
	SetVariable setvar_delta,limits={-inf,inf,0},value= root:Y_ARPES:MapInPlane:ARP3D:Ang2K:Misc:v_delta
	SetVariable setvar_srange_min_delta,pos={155,593},size={24,18},proc=Y_SetVarProc_MapIP_SetRange,title=" "
	SetVariable setvar_srange_min_delta,font="Times New Roman"
	SetVariable setvar_srange_min_delta,limits={-inf,inf,0},value= root:Y_ARPES:MapInPlane:ARP3D:Ang2K:Misc:v_srange_min_delta
	Slider slider_delta,pos={181,593},size={96,19},proc=Y_SliderProc_MapIP_A2K
	Slider slider_delta,limits={-200,200,0},variable= root:Y_ARPES:MapInPlane:ARP3D:Ang2K:Misc:v_delta,vert= 0,ticks= 0
	SetVariable setvar_srange_max_delta,pos={279,593},size={24,18},proc=Y_SetVarProc_MapIP_SetRange,title=" "
	SetVariable setvar_srange_max_delta,font="Times New Roman"
	SetVariable setvar_srange_max_delta,limits={-inf,inf,0},value= root:Y_ARPES:MapInPlane:ARP3D:Ang2K:Misc:v_srange_max_delta
	SetVariable setvar_beta0,pos={8,613},size={144,18},proc=Y_SetVarProc_MapIP_A2K,title=" \\F'Symbol'b\\F'Times New Roman'\\Z070\\K(0,0,0)\\]0 : Cut-angle origin"
	SetVariable setvar_beta0,font="Times New Roman"
	SetVariable setvar_beta0,limits={-inf,inf,0},value= root:Y_ARPES:MapInPlane:ARP3D:Ang2K:Misc:v_beta0
	SetVariable setvar_srange_min_beta0,pos={155,613},size={24,18},proc=Y_SetVarProc_MapIP_SetRange,title=" "
	SetVariable setvar_srange_min_beta0,font="Times New Roman"
	SetVariable setvar_srange_min_beta0,limits={-inf,inf,0},value= root:Y_ARPES:MapInPlane:ARP3D:Ang2K:Misc:v_srange_min_beta0
	Slider slider_beta0,pos={181,613},size={96,19},proc=Y_SliderProc_MapIP_A2K
	Slider slider_beta0,limits={-80,80,0},variable= root:Y_ARPES:MapInPlane:ARP3D:Ang2K:Misc:v_beta0,vert= 0,ticks= 0
	SetVariable setvar_srange_max_beta0,pos={279,613},size={24,18},proc=Y_SetVarProc_MapIP_SetRange,title=" "
	SetVariable setvar_srange_max_beta0,font="Times New Roman"
	SetVariable setvar_srange_max_beta0,limits={-inf,inf,0},value= root:Y_ARPES:MapInPlane:ARP3D:Ang2K:Misc:v_srange_max_beta0
	SetVariable setvar_chi,pos={8,633},size={144,18},disable=2,proc=Y_SetVarProc_MapIP_A2K,title=" \\F'Symbol'c\\F'Times New Roman'\\Z07\\K(48896,65280,48896)0\\K(0,0,0)\\]0 : Rotary angle"
	SetVariable setvar_chi,font="Times New Roman"
	SetVariable setvar_chi,limits={-inf,inf,0},value= root:Y_ARPES:MapInPlane:ARP3D:Ang2K:Misc:v_chi
	SetVariable setvar_srange_min_chi,pos={155,633},size={24,18},disable=2,proc=Y_SetVarProc_MapIP_SetRange,title=" "
	SetVariable setvar_srange_min_chi,font="Times New Roman"
	SetVariable setvar_srange_min_chi,limits={-inf,inf,0},value= root:Y_ARPES:MapInPlane:ARP3D:Ang2K:Misc:v_srange_min_chi
	Slider slider_chi,pos={181,633},size={96,19},disable=2,proc=Y_SliderProc_MapIP_A2K
	Slider slider_chi,limits={-80,80,0},variable= root:Y_ARPES:MapInPlane:ARP3D:Ang2K:Misc:v_chi,vert= 0,ticks= 0
	SetVariable setvar_srange_max_chi,pos={279,633},size={24,18},disable=2,proc=Y_SetVarProc_MapIP_SetRange,title=" "
	SetVariable setvar_srange_max_chi,font="Times New Roman"
	SetVariable setvar_srange_max_chi,limits={-inf,inf,0},value= root:Y_ARPES:MapInPlane:ARP3D:Ang2K:Misc:v_srange_max_chi
	SetVariable setvar_chi0,pos={8,653},size={144,18},disable=2,proc=Y_SetVarProc_MapIP_A2K,title=" \\F'Symbol'c\\F'Times New Roman'\\Z070\\K(0,0,0)\\]0 : Rotary origin"
	SetVariable setvar_chi0,font="Times New Roman"
	SetVariable setvar_chi0,limits={-inf,inf,0},value= root:Y_ARPES:MapInPlane:ARP3D:Ang2K:Misc:v_chi0
	SetVariable setvar_srange_min_chi0,pos={155,653},size={24,18},disable=2,proc=Y_SetVarProc_MapIP_SetRange,title=" "
	SetVariable setvar_srange_min_chi0,font="Times New Roman"
	SetVariable setvar_srange_min_chi0,limits={-inf,inf,0},value= root:Y_ARPES:MapInPlane:ARP3D:Ang2K:Misc:v_srange_min_chi0
	Slider slider_chi0,pos={181,653},size={96,19},disable=2,proc=Y_SliderProc_MapIP_A2K
	Slider slider_chi0,limits={-80,80,0},variable= root:Y_ARPES:MapInPlane:ARP3D:Ang2K:Misc:v_chi0,vert= 0,ticks= 0
	SetVariable setvar_srange_max_chi0,pos={279,653},size={24,18},disable=2,proc=Y_SetVarProc_MapIP_SetRange,title=" "
	SetVariable setvar_srange_max_chi0,font="Times New Roman"
	SetVariable setvar_srange_max_chi0,limits={-inf,inf,0},value= root:Y_ARPES:MapInPlane:ARP3D:Ang2K:Misc:v_srange_max_chi0
	PopupMenu popup_MapIP_MakeChoose3D,pos={4,24},size={186,20},bodyWidth=186,proc=Y_PopMenuProc_MapIP
	PopupMenu popup_MapIP_MakeChoose3D,font="Arial"
	PopupMenu popup_MapIP_MakeChoose3D,mode=1,popvalue="Select a 3D wave from Browser",value= #"\"Select a 3D wave from Browser;Select 2D waves in Browser; Select Folders in Browser\""
	Button button_MapIP_Make3D,pos={3,46},size={62,37},proc=Y_ButtonProc_MapIP,title="\\JC\\f01Press to\rstart"
	Button button_MapIP_Make3D,font="Arial"
	Button button_MapIP_WorkingFolder,pos={128,63},size={180,20},proc=Y_ButtonProc_MapIP_OpenFld,title="\\JL_none_"
	Button button_MapIP_WorkingFolder,font="Arial",fStyle=1
	Button button_MapIP_WorkingFolder,fColor=(65280,48896,48896)
	Button button_MapIP_WorkingWave,pos={128,46},size={180,20},proc=Y_ButtonProc_MapIP_OpenFld,title="\\JL_none_"
	Button button_MapIP_WorkingWave,font="Arial",fStyle=1,fColor=(65280,48896,48896)
	Slider slider_eng,pos={90,263},size={103,19},proc=Y_SliderProc_MapIP_2D_eng
	Slider slider_eng,labelBack=(65535,65535,65535)
	Slider slider_eng,limits={-0.775436,0.404816,0},variable= root:Y_ARPES:MapInPlane:ARP3D:Misc:v_eng,vert= 0,ticks= 0
	Slider slider_cut,pos={90,282},size={103,19},proc=Y_SliderProc_MapIP_2D_cut
	Slider slider_cut,labelBack=(65535,65535,65535)
	Slider slider_cut,limits={-10,6,1},variable= root:Y_ARPES:MapInPlane:ARP3D:Misc:v_cut,vert= 0,ticks= 0
	Slider slider_ang,pos={90,302},size={103,19},proc=Y_SliderProc_MapIP_2D_ang
	Slider slider_ang,labelBack=(65535,65535,65535)
	Slider slider_ang,limits={-9.05373,5.83152,0},variable= root:Y_ARPES:MapInPlane:ARP3D:Misc:v_ang,vert= 0,ticks= 0
	CheckBox check_width_eng,pos={200,266},size={16,14},proc=Y_CheckProc_MapIP_Width,title=""
	CheckBox check_width_eng,variable= root:Y_ARPES:MapInPlane:ARP3D:Misc:v_eng_width_flg
	CheckBox check_width_cut,pos={200,285},size={16,14},proc=Y_CheckProc_MapIP_Width,title=""
	CheckBox check_width_cut,variable= root:Y_ARPES:MapInPlane:ARP3D:Misc:v_cut_width_flg
	CheckBox check_width_ang,pos={200,304},size={16,14},proc=Y_CheckProc_MapIP_Width,title=""
	CheckBox check_width_ang,variable= root:Y_ARPES:MapInPlane:ARP3D:Misc:v_ang_width_flg
	Slider slider_width_eng,pos={247,263},size={60,19},disable=2,proc=Y_SliderProc_MapIP_2D_eng
	Slider slider_width_eng,labelBack=(65535,65535,65535)
	Slider slider_width_eng,limits={0,0.393417,0},variable= root:Y_ARPES:MapInPlane:ARP3D:Misc:v_eng_width,vert= 0,ticks= 0
	Slider slider_width_cut,pos={247,282},size={60,19},disable=2,proc=Y_SliderProc_MapIP_2D_cut
	Slider slider_width_cut,labelBack=(65535,65535,65535)
	Slider slider_width_cut,limits={0,16,2},variable= root:Y_ARPES:MapInPlane:ARP3D:Misc:v_cut_width,vert= 0,ticks= 0
	Slider slider_width_ang,pos={247,302},size={60,19},disable=2,proc=Y_SliderProc_MapIP_2D_ang
	Slider slider_width_ang,labelBack=(65535,65535,65535)
	Slider slider_width_ang,limits={0,4.96175,0},variable= root:Y_ARPES:MapInPlane:ARP3D:Misc:v_ang_width,vert= 0,ticks= 0
	Slider slider_A2K_type,pos={8,448},size={19,38},proc=Y_SliderProc_MapIP_ChangeType
	Slider slider_A2K_type,limits={1,0,1},variable= root:Y_ARPES:MapInPlane:ARP3D:Ang2K:Misc:v_type,ticks= 0
	Slider slider_A2K_kRangeFlag,pos={215,465},size={30,19}
	Slider slider_A2K_kRangeFlag,limits={1,0,1},variable= root:Y_ARPES:MapInPlane:ARP3D:Ang2K:Misc:v_k_range_flag,vert= 0,ticks= 0
	SetVariable setvar_MapIP_ItemNo,pos={192,26},size={67,18},disable=2,title=" Item #"
	SetVariable setvar_MapIP_ItemNo,font="Arial"
	SetVariable setvar_MapIP_ItemNo,value= root:Y_ARPES:MapInPlane:v_MapInPlane_index
	SetVariable setvar_cut_start,pos={181,215},size={38,18},proc=Y_SetVarProc_MapIP_SetRange,title=" "
	SetVariable setvar_cut_start,font="Times New Roman"
	SetVariable setvar_cut_start,limits={-inf,inf,0},value= root:Y_ARPES:MapInPlane:ARP3D:Misc:v_cut_start
	SetVariable setvar_cut_delta,pos={225,215},size={38,18},proc=Y_SetVarProc_MapIP_SetRange,title=" "
	SetVariable setvar_cut_delta,font="Times New Roman"
	SetVariable setvar_cut_delta,limits={-inf,inf,0},value= root:Y_ARPES:MapInPlane:ARP3D:Misc:v_cut_delta
	SetVariable setvar_cut_end,pos={270,215},size={38,18},disable=2,title=" "
	SetVariable setvar_cut_end,font="Times New Roman"
	SetVariable setvar_cut_end,limits={-inf,inf,0},value= root:Y_ARPES:MapInPlane:ARP3D:Misc:v_cut_end
	SetVariable setvar_eng,pos={49,264},size={38,18},proc=Y_SetVarProc_MapIP_2D_eng,title=" "
	SetVariable setvar_eng,labelBack=(65280,54528,48896),font="Times New Roman"
	SetVariable setvar_eng,valueBackColor=(65280,54528,48896)
	SetVariable setvar_eng,limits={-inf,inf,0},value= root:Y_ARPES:MapInPlane:ARP3D:Misc:v_eng
	SetVariable setvar_width_eng,pos={216,264},size={28,18},disable=2,proc=Y_SetVarProc_MapIP_2D_eng,title=" "
	SetVariable setvar_width_eng,labelBack=(65280,54528,48896)
	SetVariable setvar_width_eng,font="Times New Roman"
	SetVariable setvar_width_eng,valueBackColor=(65280,54528,48896)
	SetVariable setvar_width_eng,limits={-inf,inf,0},value= root:Y_ARPES:MapInPlane:ARP3D:Misc:v_eng_width
	SetVariable setvar_cut,pos={49,283},size={38,18},proc=Y_SetVarProc_MapIP_2D_cut,title=" "
	SetVariable setvar_cut,labelBack=(48896,65280,65280),font="Times New Roman"
	SetVariable setvar_cut,valueBackColor=(48896,65280,65280)
	SetVariable setvar_cut,limits={-inf,inf,0},value= root:Y_ARPES:MapInPlane:ARP3D:Misc:v_cut
	SetVariable setvar_width_cut,pos={216,283},size={28,18},disable=2,proc=Y_SetVarProc_MapIP_2D_cut,title=" "
	SetVariable setvar_width_cut,labelBack=(48896,65280,65280)
	SetVariable setvar_width_cut,font="Times New Roman"
	SetVariable setvar_width_cut,valueBackColor=(48896,65280,65280)
	SetVariable setvar_width_cut,limits={0,inf,0},value= root:Y_ARPES:MapInPlane:ARP3D:Misc:v_cut_width
	SetVariable setvar_ang,pos={49,302},size={38,18},proc=Y_SetVarProc_MapIP_2D_ang,title=" "
	SetVariable setvar_ang,labelBack=(65280,65280,48896),font="Times New Roman"
	SetVariable setvar_ang,valueBackColor=(65280,65280,48896)
	SetVariable setvar_ang,limits={-inf,inf,0},value= root:Y_ARPES:MapInPlane:ARP3D:Misc:v_ang
	SetVariable setvar_width_ang,pos={216,302},size={28,18},disable=2,proc=Y_SetVarProc_MapIP_2D_ang,title=" "
	SetVariable setvar_width_ang,labelBack=(65280,65280,48896)
	SetVariable setvar_width_ang,font="Times New Roman"
	SetVariable setvar_width_ang,valueBackColor=(65280,65280,48896)
	SetVariable setvar_width_ang,limits={-inf,inf,0},value= root:Y_ARPES:MapInPlane:ARP3D:Misc:v_ang_width
	SetVariable setvar_kRange_ky_max,pos={153,362},size={38,18},proc=Y_SetVarProc_MapIP_SetRange,title=" "
	SetVariable setvar_kRange_ky_max,font="Times New Roman"
	SetVariable setvar_kRange_ky_max,valueBackColor=(51456,44032,58880)
	SetVariable setvar_kRange_ky_max,limits={-inf,inf,0},value= root:Y_ARPES:MapInPlane:ARP3D:Ang2K:Misc:v_ky_max
	SetVariable setvar_kRange_ky_min,pos={153,417},size={38,18},proc=Y_SetVarProc_MapIP_SetRange,title=" "
	SetVariable setvar_kRange_ky_min,font="Times New Roman"
	SetVariable setvar_kRange_ky_min,valueBackColor=(51456,44032,58880)
	SetVariable setvar_kRange_ky_min,limits={-inf,inf,0},value= root:Y_ARPES:MapInPlane:ARP3D:Ang2K:Misc:v_ky_min
	SetVariable setvar_kRange_kx_min,pos={186,438},size={38,18},proc=Y_SetVarProc_MapIP_SetRange,title=" "
	SetVariable setvar_kRange_kx_min,font="Times New Roman"
	SetVariable setvar_kRange_kx_min,valueBackColor=(51456,44032,58880)
	SetVariable setvar_kRange_kx_min,limits={-inf,inf,0},value= root:Y_ARPES:MapInPlane:ARP3D:Ang2K:Misc:v_kx_min
	SetVariable setvar_kRange_kx_max,pos={235,439},size={38,18},proc=Y_SetVarProc_MapIP_SetRange,title=" "
	SetVariable setvar_kRange_kx_max,font="Times New Roman"
	SetVariable setvar_kRange_kx_max,valueBackColor=(51456,44032,58880)
	SetVariable setvar_kRange_kx_max,limits={-inf,inf,0},value= root:Y_ARPES:MapInPlane:ARP3D:Ang2K:Misc:v_kx_max
	SetVariable setvar_Mesh_kx,pos={220,351},size={24,18},proc=Y_SetVarProc_MapIP_SetRange,title=" "
	SetVariable setvar_Mesh_kx,font="Times New Roman"
	SetVariable setvar_Mesh_kx,valueBackColor=(32768,65280,32768)
	SetVariable setvar_Mesh_kx,limits={1,inf,0},value= root:Y_ARPES:MapInPlane:ARP3D:Ang2K:Misc:v_mesh_kx
	SetVariable setvar_Mesh_ky,pos={262,390},size={24,18},proc=Y_SetVarProc_MapIP_SetRange,title=" "
	SetVariable setvar_Mesh_ky,font="Times New Roman"
	SetVariable setvar_Mesh_ky,valueBackColor=(32768,65280,32768)
	SetVariable setvar_Mesh_ky,limits={1,inf,0},value= root:Y_ARPES:MapInPlane:ARP3D:Ang2K:Misc:v_mesh_ky
	Button button_Update_A2K_3DMap,pos={38,681},size={250,23},proc=Y_ButtonProc_MapIP_A2K,title="\\JC\\f01Map 3D Eng-\\F'Symbol'a\\]0\\f01-\\F'Symbol'b\\]0\\f01 --> 3D Eng-kx-ky"
	Button button_Update_A2K_3DMap,fColor=(65535,65535,65535)
	CheckBox check_DA,pos={121,468},size={16,14},proc=Y_CheckProc_MapIP_A2K,title=""
	CheckBox check_DA,variable= root:Y_ARPES:MapInPlane:ARP3D:Ang2K:Misc:v_DA
End




////////////////////////////////////
////  Below are functions for controls
////////////////////////////////////


Function Y_CheckProc_MapIP_A2K(ctrlName,checked) : CheckBoxControl
	String ctrlName
	Variable checked
	
	if(stringmatch(ctrlName, "check_DA"))
	
		Y_MapIP_A2K_PanelDA()
		Y_MapIP_A2K_ChangeType()
		
	endif

End

Function Y_ButtonProc_MapIP_A2K(ctrlName) : ButtonControl
	String ctrlName
	
	if(stringmatch(ctrlName, "button_Update_A2K_3DMap"))
	
		Y_MapIP_A2K_MapTo3D()
	
	elseif(stringmatch(ctrlName, "button_Disp_A2K"))
	
		Y_MapIP_A2K_GraphDisp()
	
	endif

End


Function Y_SetVarProc_MapIP_A2K(ctrlName,varNum,varStr,varName) : SetVariableControl
	String ctrlName
	Variable varNum
	String varStr, varName
	
	Y_MapIP_A2K_GraphDisp()
	Y_MapIP_A2K()

End


Function Y_SliderProc_MapIP_A2K(ctrlName,sliderValue,event) : SliderControl
	String ctrlName
	Variable sliderValue
	Variable event	// bit field: bit 0: value set, 1: mouse down, 2: mouse up, 3: mouse moved

	if(event %& 2)	
		
		Y_MapIP_A2K_GraphDisp()

	endif
	
	Y_MapIP_A2K()

	return 0
End


Function Y_SliderProc_MapIP_ChangeType(ctrlName,sliderValue,event) : SliderControl
	String ctrlName
	Variable sliderValue
	Variable event	// bit field: bit 0: value set, 1: mouse down, 2: mouse up, 3: mouse moved

	if(event %& 2)	
		
		Y_MapIP_A2K_GraphDisp()

	endif
	
	Y_MapIP_A2K()
	Y_MapIP_A2K_ChangeType()

	return 0
End


Function Y_SetVarProc_MapIP_SetRange(ctrlName,varNum,varStr,varName) : SetVariableControl
	String ctrlName
	Variable varNum
	String varStr, varName
	
	if(stringmatch(ctrlName,"setvar_srange_m*"))
	
		string s_selected_var = ctrlName[18,strlen(ctrlName)-1]
		Y_MapIP_RangeSlider(s_selected_var)

	elseif(stringmatch(ctrlName, "setvar_cut_*"))

		Y_MapIP_RangeCut()
		Y_MapIP_2D_RangeSlider("_width_cut", 0)
		Y_MapIP_A2K_BoundaryMake()
		Y_MapIP_A2K()
	
	elseif(stringmatch(ctrlName, "setvar_Mesh_*"))
	
		Y_MapIP_RangeMesh()
		Y_MapIP_Rangekxky()
		Y_MapIP_A2K_MapOntoImage()
		
	elseif(stringmatch(ctrlName, "setvar_kRange_k*")) 

		Y_MapIP_Rangekxky()
		
		string s_fld = GetDataFolder(1)	

		SetDataFolder root:Y_ARPES:MapInPlane
		string/G s_MapInPlane_WorkFld
		SetDataFolder $(s_MapInPlane_WorkFld + "Ang2K:Misc:")
		variable/G v_k_range_flag = 1
		
		SetDataFolder $s_fld
	endif

End