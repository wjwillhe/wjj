load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

undef("ttest2")
function ttest2(x[*][*][*]:float,y[*][*][*]:float,sig:float,n:integer)
local prob_phase,ave_sample,ave_all,var_sample,n
begin
	ave_sample = dim_avg_n_Wrap(x(:,:,:),0)
	ave_all = dim_avg_n_Wrap(y(:,:,:),0)
	var_sample = dim_variance_n_Wrap(x(:,:,:),0)
	;printVarSummary(var_sample)
	
	
	prob_phase = new((/dimsizes(x&lat),dimsizes(x&lon)/),float)
	prob_phase!0 = "lat"
	prob_phase!1 = "lon"
	prob_phase&lat = x&lat
	prob_phase&lon = x&lon
	
	do j = 0,dimsizes(x&lat)-1
		do jj = 0,dimsizes(x&lon)-1		     
				if(ismissing(var_sample(j,jj)))then						
				        var_sample(j,jj)=0
				end if
		end do
	end do
	
	do j = 0,dimsizes(x&lat)-1
		do jj = 0,dimsizes(x&lon)-1		     
				if(var_sample(j,jj).ne.0)then
						prob_phase(j,jj) = abs(ave_sample(j,jj)-ave_all(j,jj))*sqrt(n)/sqrt(var_sample(j,jj))				
				end if
		end do
	end do
    confi=ave_sample
	confi=student_t(prob_phase,n-1)
	;print(ave_sample)
	; print(confi)
	;ave_sample(:,:) = mask(ave_sample(:,:),confi(:,:).gt.sig,False)
	; ave_sample(:,:) = where(confi(:,:).gt.sig,0,ave_sample(:,:))
    ;print(ave_sample)
	return(confi)
	
end

undef("polar_map_circle")
procedure polar_map_circle(wks,plot:graphic,wsize:integer,col:string,offset:numeric)
local degrad,degrees,xcos,xsin,xcenter,ycenter,radius,xc,yc
begin
getvalues plot
"vpXF"       :x
"vpYF"         :y
"vpWidthF"      :w
"vpHeightF"    :h  

end getvalues
degrad = 0.017453292519943
degrees = ispan(0,360,1)
xcos = cos(degrad*degrees)
xsin = sin(degrad*degrees)
xcenter = w/2 +x
ycenter = h/2 +(y-h)
radius = w/2+ offset

xc = xcenter + radius*xcos

yc = ycenter + radius*xcos
lnre = True
lnre@gsLineColor = col
lnre@gsLineThicknessF = wsize
gsn_polyline_ndc(wks,xc,yc,lnre)
end 


begin
data=asciiread("/public/home/wangji19/data/mg1955-2014/three_path/clus1_1532.txt",-1,"integer")
cy=addfile("/public/home/wangji19/data/track_data/TRACK_ff_trs_pos_ERA5_1950_202207_NH_CYC.nc","r")
list_of_files1=systemfunc("ls /public/home/wangji19/data/dhgt/*.nc")
  fz=addfiles(list_of_files1,"r")
list_of_files2=systemfunc("ls /public/home/wangji19/data/du/*.nc")
  fu=addfiles(list_of_files2,"r")
list_of_files4=systemfunc("ls /public/home/wangji19/data/dair/*.nc")
  ft=addfiles(list_of_files4,"r")  
  
npts=dimsizes(data)
print("Number of values: "+npts)
x=data
nnumb=dimsizes(x)
path=ispan(0,nnumb-1,1)
;print(x)
day=ispan(-10,10,1)

point=ispan(-10,10,1)


; letter=(/"0","1","2","3","4","(a)","6",\
		 ; "(b)","8","(c)","10","(d)","12","(e)",\
		 ; "14","15","16","17","18","19","20"/);;
		 
; letter=(/"0","1","2","3","4","(f)","6",\
		 ; "(g)","8","(h)","10","(i)","12","(j)",\
		 ; "14","15","16","17","18","19","20"/);;
		 
letter=(/"0","1","2","3","4","(k)","6",\
		 "(l)","8","(m)","10","(n)","12","(o)",\
		 "14","15","16","17","18","19","20"/);
	

;气旋数据
time1=cy->traj_time 

time1_FillValue=-2147483647
time1=time1*6
time1@units="hours since 1950-01-01 00:00:0.0"
time1@long_name="Time of the trajectory"
Ti1=cd_calendar(time1(x,0),-2)
delete(time1)
delete(cy)
;;;;;;;
 
time=fu[:]->time

hgt1=fz[:]->hgt(:,{875:825},{0:90},:)

u1=fu[:]->uwnd(:,{875:825},{0:90},:)
t1=ft[:]->air(:,{875:825},{0:90},:)
lat1 =hgt1&lat
level=hgt1&level
printVarSummary(u1)
delete(fz)
delete(fu)
delete(ft)


re=6378388
z=new((/nnum,17,361,1441/),float)
z(:,:,:,:)=re*hgt1(:,:,:,:)/(re-hgt1(:,:,:,:))
copy_VarMeta(hgt1(:,:,:,:),z(:,:,:,:))

;;;;求位温
th=new((/nnum,3,361,1441/),float) 
level1=level*100
leveltmp =conform(t1, level1, 1)
th(:,:,:,:) = pot_temp(leveltmp ,t1(:,:,:,:) ,-1, False)

delete(t1)
delete(leveltmp)
delete(hgt1)
lattmp=conform(th, lat1, 2)
egr=new((/nnum,3,361,1441/),float)
eady=new((/nnum,361,1441/),float)
egr(:,:,:,:)=eady_growth_rate(th,u1,z,lattmp,0,1)
eady(:,:,:)=egr(:,{850},:,:)*86400
eady@_FillValue=-2147483647
copy_VarMeta(z(:,:,:,:),egr(:,:,:,:))
copy_VarMeta(egr(:,{850},:,:),eady(:,:,:))

delete(z)
delete(egr)
delete(th)
delete(u1)
                                                       
   TIME2    = cd_calendar(time, 0)          ; type float 
   year    = toint( TIME2(:,0) )
   month   = toint( TIME2(:,1) )
   day     = toint( TIME2(:,2) ) 
   ddd     = day_of_year(year, month, day) 
   yyyyddd = year*1000 + ddd 
	 
    eadyAnom  =  eady
    eadyAnom@_FillValue=-999.0
    eadyClmDay  =  new((/d1,37,144/),float)
    eadyClmDay@_FillValue=-999.0
eadyClmDay    = clmDayTLL(eady, yyyyddd)
   eadyAnom      = calcDayAnomTLL (eady, yyyyddd, eadyClmDay)

 
;;;-----------

ready=new((/nnumb,21,361,1441/),float)
anoeady=new((/nnumb,21,361,1441/),float)
do i=0,nnumb-1
mm=ind(Ti1(i).eq.Time3)
ready(i,:,:,:)=eady(mm-10:mm+10,:,:)
anoeady(i,:,:,:)=eadyAnom(mm-10:mm+10,:,:)
end do
ready!0="path"
ready&path=path
ready!1="day"
ready&day=day
copy_VarMeta(eady(0,:,:),ready(0,0,:,:))
copy_VarMeta(ready(:,:,:,:),anoeady(:,:,:,:))

 eady_n=new((/21,361,1441/),float)
 anoeady_n=new((/21,361,1441/),float)
 eady_n(:,:,:)=dim_avg_n_Wrap(ready(:,:,:,:),0)
 anoeady_n(:,:,:)=dim_avg_n_Wrap(anoeady(:,:,:,:),0)
copy_VarMeta(ready(0,:,:,:),eady_n(:,:,:))
copy_VarMeta(eady_n(:,:,:),anoeady_n(:,:,:))

 eady_std=new((/21,361,1441/),float)
t_eady=new((/21,361,1441/),float)
confi_eady=new((/21,361,1441/),float)
do n=0,20
eady_std(n,:,:)=dim_variance_n_Wrap(anoeady(:,n,:,:),0)
eady_std(n,:,:)=sqrt(eady_std(n,:,:)/nnumb)
eady_std=where(eady_std.eq.0,eady_std@_FillValue,eady_std)
t_eady(n,:,:)=anoeady_n(n,:,:)/eady_std(n,:,:)
confi_eady(n,:,:)=anoeady_n(n,:,:)
confi_eady(n,:,:)=student_t(t_eady(n,:,:),nnumb-1)
end do

 
; ; ; ; eadyave=dim_avg_n_Wrap(dim_avg_n_Wrap(ready(:,9:13,:,:),1),0)
; ; ; ; copy_VarMeta(eady(0,0,:,:),eadyave(:,:))
; ; ; ; outfile=addfile("/public/home/wangji19/data/mg1955-2014/eady_e1532.nc","c") 
  
; ; ; ; ; ; outfile->eady=eady
; ; ; ; ; ; outfile@history="Method 6:eady of east path"
; ; ; ; exit
; Dir = "/public/home/wangji19/nclprogram/"
; wks = gsn_open_wks("eps",Dir+point(n)+"_850hPa_eady_es")   
; wks = gsn_open_wks("eps","/public/home/wangji19/nclprogram/1The_es1532_-1+3_300hPa_uclm&uano")  

wks = gsn_open_wks("eps","/public/home/wangji19/nclprogram/The_clus1_-5+3_850hPa_eady_ano")  

map  = new(21,graphic)   
plot  = new(21,graphic)                ; create graphic array_lowpass
plot1 = new(21,graphic)
plot_p = new(21,graphic)
do n=0,20
gsn_define_colormap(wks,"cmp_b2r" )

res=True
res@gsnAddCyclic = True
; res@gsnMaximize  = False
res@gsnDraw=False
res@gsnFrame=False
res@gsnLeftString=""
res@gsnRightString=""
res@tiMainString="" 

mpres=res
mpres@gsnPolarNH = True
mpres@mpFillOn=False
mpres@mpPerimOn = 	True
mpres@mpPerimDrawOrder = "PostDraw"
mpres@mpPerimLineThicknessF = 0.1
mpres@mpGridLatSpacingF    = 45.                ; change latitude  line spacing
mpres@mpGridLonSpacingF    = 90.
mpres@mpGridLineThicknessF=0.05
mpres@mpGeophysicalLineThicknessF=0.1

mpres@mpCenterLonF=90
mpres@mpMinLatF=20
mpres@mpMaxLatF=90
map(n)=gsn_csm_map_polar(wks,mpres)
 
resc=res
resc@gsnAddCyclic = True
resc@cnFillOn=True
resc@cnLinesOn=False
resc@cnLineLabelsOn=False
resc@cnInfoLabelOn=False
resc@lbLabelBarOn         =  False
resc@lbLabelFontHeightF   =0.01
resc@pmLabelBarOrthogonalPosF  =0.01
resc@pmLabelBarHeightF    =0.1
resc@pmLabelBarWidthF     =0.6


resc@cnLevelSelectionMode = "ManualLevels"       ; set manual contour levels
resc@cnMinLevelValF       = -0.21              ; set min contour level
resc@cnMaxLevelValF       = 0.21               ; set max contour level
resc@cnLevelSpacingF      = 0.03  
plot1(n)=gsn_csm_contour(wks,anoeady_n(n,:,:),resc)
overlay(map(n),plot1(n))

;;;;;;;;;;;;;;;
res2 = res
res2@cnFillOn                = True
res2@cnLinesOn               = False
res2@cnLevelSelectionMode   = "ExplicitLevels"
res2@cnLevels               = (/0.01/) ; significance level
res2@cnFillPalette          ="GMT_gray"
res2@cnFillColors           =(/2,-1/)
res2@cnFillPatterns         = (/17,-1/);;;;
res2@cnFillScales           = (/0.5,1/)
res2@cnFillDotSizeF         = 0.001;0.003
; res2@cnFillScales           = 2
res2@cnMonoFillPattern      = False
res2@cnMonoFillScale        = False
res2@cnInfoLabelOn          = False
res2@cnLineLabelsOn         = False
res2@lbLabelBarOn           = False
plot_p(n)=gsn_csm_contour(wks,confi_eady(n,:,:),res2)
overlay(map(n),plot_p(n))

cnres=res
cnres@gsnAddCyclic = True
cnres@cnLinesOn=True
cnres@cnFillOn=False
cnres@gsnLeftString=letter(n)+" "+"Day"+" "+day(n)
cnres@cnLevelSelectionMode   = "ExplicitLevels"
cnres@cnLevels=(/0.2,0.4,0.6,0.8,1.0,1.2/)
cnres@cnLineColor="black"
cnres@cnLineThicknessF=1;2

cnres@cnLineLabelsOn=True
cnres@cnLineLabelPerimOn =True
cnres@cnLabelDrawOrder     = "PostDraw"
cnres@cnLineLabelFontHeightF=0.008
cnres@cnLineLabelBackgroundColor = "white"
cnres@cnLineLabelInterval=1
cnres@cnLineLabelDensityF=0.5
cnres@cnLabelMasking=True
cnres@cnLineLabelPerimColor = "white"
cnres@cnLineLabelPerimThicknessF = 0.5
cnres@cnLineLabelPlacementMode="Constant";“Computed”
cnres@cnInfoLabelOn=False


plot(n)=gsn_csm_contour(wks,eady_n(n,:,:) ,cnres)
overlay(map(n),plot(n))
  end do
  

  mapdraw = new(5,graphic)
  mapdraw(0:4) = map(5:13:2)
 
pres=True
pres@gsnMaximize                = True;False
pres@txString=""
pres@gsnPanelFigureStringsFontHeightF  =0.015
pres@amJust="TopLeft"
pres@pmLabelBarOrthogonalPosF  = -0.01
pres@gsnPanelLabelBar         =  True;False
pres@pmLabelBarHeightF   =0.03
pres@pmLabelBarWidthF    =0.5
gsn_panel(wks,mapdraw,(/5,1/),pres)

end