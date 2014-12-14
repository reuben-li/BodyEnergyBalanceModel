import math
import textwrap
import matplotlib.pyplot as plt
from numpy import genfromtxt
import csv

def main():
	weather = genfromtxt('weather_h1', delimiter=',') #weather file
	fig = plt.figure(figsize=(6, 9), dpi=100)
	subject = [] #weight,age,gender,MET5,MET10,MET15,colour,fig,weather,height,identifier,clo,tcr_base
	subject.append([67.765,22,'male',4,3.6,3.1,'g',fig,weather,1.781,1,0.57,36.93]) #check 100
	subject.append([57.54,22,'male',6.3,5.6,5.1,'b',fig,weather,1.682,2,0.57,36.91])
	subject.append([70.650,22,'male',5.9,5.4,5.4,'r',fig,weather,1.77,3,0.57,36.73])
	#subject.append([10.0,1,'female','sitting','y',fig,weather,0.8,4])
	#iterate over subjects
	i=0
	while i < len(subject):
		run(*subject[i]) 
		i=i+1
	plt.draw()
	plt.show(block=True)

def input_phy(activity): 	#define input individual variables
	global clo,activity_code
	clo = 1.0 # 0 (naked); summer clothes (0.6); ski (2)
	activity_dict={'rest':1, 'sitting':2, 'walking':3,'light':4}
	activity_code = activity_dict[activity]

def initial():
	global time_unit, simulation_time, time_steps, artificial_weather
	artificial_weather = 0; # 0=false, 1=true
	time_unit=20		#in seconds
	simulation_time=0.5	#in hours
	time_steps= simulation_time * 3600.0 / time_unit
	#~ return {'tsk':tsk,'tcr':tcr}

def pttbl(x): 		#function to convert saturation vapour pressure units
	svp = 6.11*10**((7.5*x)/(237.7+x))	#saturation vapour pressure (hpa)
	svp_m = 0.750061683 * svp 			#svp in  mmhg
	return svp_m

## METABOLISM FUNCTIONS
def bmr(coeff,inter,weight):
	global bmr_W
	# if activity_code == 1:
		# mult = 1.0
	# elif activity_code == 2:
		# mult = 1.5
	# elif activity_code == 3:
		# mult = 3.2
	# elif activity_code == 4:
		# mult = 5
	bmr_MJ = ((coeff*weight)+inter) 
	bmr_W = bmr_MJ * 11.57 		#basal metabolic rate
	mr_W = bmr_W #* mult			#metabolic rate due to activities				
	return mr_W

def metabol(gender,age,weight,met5,met10,met15): 		#function for metabolic rate approximation
	me=0.2 			#mechanical efficiency 20% is practical number
	if gender == 'male':	
		if age < 3:
			rm = bmr(0.249,-0.127,weight)
		elif age < 10:
			rm = bmr(0.095,2.110,weight)
		elif age < 18:
			rm = bmr(0.074,2.754,weight)
		elif age < 30:
			rm = bmr(0.063,2.896,weight)
		elif age < 60:
			rm = bmr(0.049,2.450,weight)
		else: 
			rm = bmr(0.049,2.459,weight)  
	if gender == 'female':
		if age < 3:
			rm = bmr(0.244,-0.130,weight)
		elif age < 10:
			rm = bmr(0.085,2.033,weight)
		elif age < 18:
			rm = bmr(0.056,2.898,weight)
		elif age < 30:
			rm = bmr(0.062,2.036,weight)
		elif age < 60:
			rm = bmr(0.034,3.538,weight)
		else: 
			rm = bmr(0.038,2.755,weight)  
	
	#MET value at each 5 min interval
	
	met=58.2
	
	rm5 = met*met5
	rm10 = met*met10
	rm15 = met*met15
	rmrest = rm
		
	H = (1-me)*rm	# metabolic rate - work done (not converted to internal heat)
	wk = me*rm  # work done
	
	return {'metabolic_rate':rm, 'internal_heat_production':H, "work done":wk,"rm5":rm5,"rm10":rm10,"rm15":rm15,"rmrest":rmrest}
	
def run(weight,age,gender,met5,met10,met15,colour,fig,weather,height,ID,clo,tcr_base):	#main function (physiological inputs)
	
	#activity_code=input_phy(activity)
	#~ input_met()
	#get metabolic rate and activity label
	meta=metabol(gender,age,weight,met5,met10,met15)
	wk=meta['work done']
	rm=meta['metabolic_rate']
	H=meta['internal_heat_production']
	rm5=meta["rm5"]
	rm10=meta["rm10"]
	rm15=meta["rm15"]
	rmrest=meta["rmrest"]
	#intialization
	ta=weather[0][0]
	tmrt=weather[0][2]
	rh_per=weather[0][1]
	v=weather[0][3]
	rh=rh_per/100.0
	
	tsk_set=34.1 #SETPOINTS
	tcr_set=36.6
	#tcr_set=tcr_base
	
	time=0.0
	init=initial()
	#~ tsk=init['tsk']
	#~ tcr=init['tcr']
	tsk=tsk_set
	tcr=tcr_set
	#v=0.5
	ti=1.0
	#convective heat transfer
	if gender == 'male':
		h_b = 0.560 * height/100 #center of mass
	else:
		h_b = 0.543 * height/100
	v = v * math.log(h_b/0.01)/math.log(10/0.01)			#wind velocity at centre of mass	
	htc_c = 4.0*v + 0.35*v * ti - 0.00080*(v*ti)**2 + 3.4
	
	#body metrics
	wt=weight
	ht=height*100.0
	if age < 1:
		adu = wt**0.473 * ht**0.655 * 95.68
	elif age < 6:
		adu = wt**0.423 * ht**0.362 * 381.89
	else:
		adu = wt**0.444 * ht**0.663 * 88.83
	adu_m=adu/10000
	svp = 6.11*10**((7.5*ta)/(237.7+ta))	#saturation vapour pressure (hpa)
	svp_m = 0.750061683 * svp 			#svp in  mmhg
	
	edif=5.0
	# eres=0.0023*rm*(44-rh*pttbl(ta))
	eres=0.0023*rm*(44-rh*pttbl(ta)) #moderated by body area
	ev=eres+edif	
	ersw=0
	edrip=0
	wrsw=0	
	chr=5.23						#radiative heat transfer coefficient Gagge
	chr2=4*0.98*5.67*10**-8*0.7*((273.15+(tsk+tmrt)/2)**3) #parsons 1993
	
	print chr, chr2
	
	chc=htc_c 						#convective heat transfer coefficient
	ctc=chr+chc						#combined heat transfer coefficient
	fcl=1.0/(1.0+0.155*ctc*clo)			#clothing thermal efficiency factor this is actually Fcl differentiated from fcl
	fpcl=1./(1.0+0.143*(chc)*clo) #permeation efficiency factor for clothing
	kmin=5.28						#min. skin heat conductance W/sqm C
	skbfn=6.30
	skbf=skbfn
	mshiv=0
	#arrays to store variables
	atsk=[]
	atcr=[]
	astorage =[]
	xaxis = []
	aev = []
	asignal =[]
	awetted = []
	adry = []
	ata = []
	arh = []
	askbf = []
	aC = []
	aR = []
	atmrt =[]


	
## DUMMY LOOP FOR INITIAL STEADY STATE #################################
	while (time <= 1): 	#start of loop
		C =  fcl * chc * (ta-tsk)  #adu left out cos it will be accounted for
		sbc = 5.670373 * 10**-8
		feff = 0.7
		emiss = 0.97
		R = fcl* feff * emiss * sbc * ((tmrt+273.15)**4 - (tsk+273.15)**4)
		alpha=0.044 + 0.35/ (skbf-0.1386) #mass proportion of skin shell to total body weight
		wt_sk = weight*alpha #skinshell weight 7% (gagge 5%, others 30%)
		wt_cr = weight-wt_sk
		eres=0.0023*rm*(44-rh*pttbl(ta)) 
		hfcr=rm+mshiv-(tcr-tsk)*(kmin+1.163*skbf)-eres-wk #heat flow from core to skin W/m2. 1.163 is the specific heat of blood
	
		hfsk=(tcr-tsk)*(kmin+1.163*skbf)+(R+C)-(ev-eres) #heat flow from skin to air W/m2
		#thermal capacity of skin and core
		tcsk=0.97*wt_sk
		tccr=0.97*wt_cr
		#change in temp/hr
		dtsk=(hfsk*adu_m)/tcsk*simulation_time
		dtcr=(hfcr*adu_m)/tccr*simulation_time
		
		
		#unit of time
		dtim=1.0/60
		u=abs(dtsk)	
		if(u*dtim-0.1>0):
			dtim=0.1/u
		u=abs(dtcr)
		if(u*dtim-0.1>0):
			dtim=0.1/u			
		time=time+dtim
		tsk=tsk+dtsk*dtim    #normalized by time 
		tcr=tcr+dtcr*dtim	
		#control system	
		#shivering
		if age > 1:
			mshiv= 19.4 * max((tsk_set - tsk),0) * max((tcr_set - tcr),0)  #W/m2 Gagge 1971
			
			
		#blood flow
		sksig=(tsk-tsk_set)
		if (sksig>0):
			colds=0.0
			warms=sksig
		else:
			colds=-sksig
			warms=0.0		
		crsig=(tcr-tcr_set)
		if (crsig>0):
			warmc=crsig
			coldc=0	
		else:
			warmc=0.0
			coldc=-crsig
		stric=0.5*colds	#constriction gain heat at skin
		dilat=75*warmc	#dilation loses heat from core
		skbf=(skbfn+dilat)/(1.0+stric)
		
		
			
		if (rm-bmr_W<=0):  #at rest. edited to use bmr_W
			regsw=100*warmc*warms
		else:
			regsw=250*warmc+100*warmc*warms		
		ersw=0.7*regsw*2**((tsk-tsk_set)/3) # this is a power function; 2 != adu_m	
		
		if gender == 'male': 
			swf = 1 #factor for male/female difference (Hoppe 1993)
		else:
			swf = 0.7	
		vhw_w = 675.0 	# Wh/kg latent heat of sweat 
		sw = max(8.47 * 10**-5 * ((alpha* tsk + (1-alpha) * tcr) - tcr_set),0)*swf #Hardy 1978 kg/(s m2)		
		sw_h = sw * 3600  #kg/(h m2) Hoppe
		ersw = sw_h * vhw_w   #W/m2     Hoppe
		if (ersw-500>0):
			exit()		
		
		wrsw=wrsw+ (ersw*adu_m/(0.7*100))*dtim 
		#wrsw = reg. sweat in 100cc units per man for time exp; 0.7 = latent heat of sweat
		emax=2.2*chc*pttbl(tsk)-rh*pttbl(ta)*fpcl
		prsw=ersw/emax
		pwet=0.06+.94*prsw #0.06 = dampness factor of skin; prsw = wettedness due to diff
		#edif=pwet*emax-ersw
		
		svp = 6.11*10**((7.5*ta)/(237.7+ta))	#saturation vapour pressure (hpa)
		vp = svp*rh 							#actual vapour pressure (hpa)
		svp_m = 0.750061683 * svp 				#svp in  mmhg
		vp_m = svp_m*rh	
		edif=-1.694 * 10**-7  * 2430000 * (vp_m - pttbl(tsk))  #terms = permeance coefficient and vhw (J/kg)
		
		
		ev=eres+ersw+edif
			
		if (ersw-emax>0):
			edrip=ersw-emax
			ev=eres+edif
			ersw=emax
			edif=0.0
			prsw=1.0
			pwet=1.0		
		C =  1 * chc * (ta-tsk)
		sbc = 5.670373 * 10**-8
		feff = 0.7
		emiss = 0.97	
		R =  1 * feff * emiss * sbc * ((tmrt+273.15)**4 - (tsk+273.15)**4)		
	
	period=0
	time=0.0
	#tsk=34.1
	
	
#################################################################################################################
## ACTUAL SIMULATION ############################################################################################
#################################################################################################################

	#tcr=tcr_base #use measured oral temperature as base
	#tsk=35
	while (time <= 1): 	#start of loop
		
		me = 0.2
		
		juncture=time*30
		if (juncture <= 5):
			rm = rm5
			
			wk = me*rm
		elif (juncture <=10):
			rm = rm10
			wk = me*rm
		elif (juncture <=15):
			rm = rm15
			wk = me*rm
		else:
			rm = rmrest
			#me = 1.0
			wk = me*rm
		
		print rm
		
		
		
		
		chc=4.0*v + 0.35*v * ti - 0.00080*(v*ti)**2 + 3.4 #Ooka
		#chc=8.6*v**0.5
		
		hfcr=rm+mshiv-(tcr-tsk)*(kmin+1.163*skbf)-eres-wk #heat flow to core W/m2 // +mshiv?
		hfsk=(tcr-tsk)*(kmin+1.163*skbf)+(R+C)-(ev-eres) #heat flow to skin W/m2
		#thermal capacity of skin and core
		tcsk=0.97*wt_sk #added dynamic value
		tccr=0.97*wt_cr
		#change in temp per hr (need to account for simulation time different from an hour)
		offset=simulation_time
		dtsk=(hfsk*adu_m)/tcsk*offset
		dtcr=(hfcr*adu_m)/tccr*offset
		#unit of time
		dtim=1.0/time_steps
	
	
		t = round(time*time_steps*time_unit/60)
		
		u=abs(dtsk)	
		if(u*dtim-0.1>0):
			dtim=0.1/u
		u=abs(dtcr)
		if(u*dtim-0.1>0):
			dtim=0.1/u
			
		time=time+dtim
		tsk=tsk+dtsk*dtim
		tcr=tcr+dtcr*dtim
		#~ print str(dtim) + " " + str(period)
		
## THERMOREGULATION - BLOOD FLOW AND SHIVERING
				
		#shivering
		if age > 1:
			mshiv= 19.4 * max((tsk_set - tsk),0) * max((tcr_set - tcr),0)  #W/m2 Gagge 1971 both terms must be positive else 0.
			#mshiv3=max(69.7*max((tcr-tcr_set),0)*(tsk-tsk_set),0)
			#mshiv= 69.73 * max((tsk_set - tsk),0) * max((tcr_set - tcr),0)
		
		#blood flow
		sksig=(tsk-tsk_set)
		if (sksig>0):
			colds=0.0
			warms=sksig
		else:
			colds=-sksig
			warms=0.0
			
		crsig=(tcr-tcr_set)
		if (crsig>0):
			warmc=crsig
			coldc=0	
		else:
			warmc=0.0
			coldc=-crsig

		stric=0.5*colds	#constriction reduction in heat loss at skin (L/hr/sqm) hands (but not core!)
		dilat=75*warmc	#dilation loses heat from core (L/hr/sqm)
		skbf=(skbfn+dilat)/(1.0+stric)
		#Hardy 1978
		#skbf= max(10.6+36*(tcr-tcr_set)*(tsk-tsk_set) + 0.93*(tsk-tsk_set),0)
					
		
## EVAPORATIVE FLUX - SWEAT, DIFFUSION AND RESPIRATION
		
		vhw_w = 675.0 	# Wh/kg latent heat of sweat 
		
		if (rm-bmr_W<=0):  #account for rest
			regsw=100.0*warmc*warms	#regsw ( g/(hr m2)
		else:
			regsw=250.0*warmc+100.0*warmc*warms 
		ersw=(vhw_w/1000)*regsw*2.0**((tsk-tsk_set)/3.0)
		
		#vhw = 2.23*10**6   # J/kg vaporisation heat of water
		
		if gender == 'male': 
			swf = 1 #factor for male/female difference (Hoppe 1993)
		else:
			swf = 0.7	
		sw = max(8.47 * 10**-5 * ((alpha* tsk + (1-alpha) * tcr) - tcr_set),0)*swf #Hardy 1978 kg/(s m2)		

		sw_h = sw * 3600  #kg/(h m2) Hoppe
		ersw = sw_h * vhw_w   #W/m2     Hoppe
		
		#~ if (ersw-500>0):  # Gagge
			#~ exit()
		
		if (ersw-emax>0): # If sweat produced > maximum cooling potential
			edrip=ersw-emax
			ev=eres+edif	
			ersw=emax
			edif=0.0
			prsw=1.0
			pwet=1.0

		
## EVALUATE NEW DX/DT FOR EACH COMPONENT 
		
		# "Wet" components - sweat, diffusion and respiration
		wrsw=wrsw+ (ersw*2/(0.7*100))*dtim
		emax=2.2*chc*pttbl(tsk)-rh*pttbl(ta)*fpcl
		prsw=ersw/emax
		pwet=0.06+.94*prsw
		#edif=(pwet*emax-ersw)
		
		svp = 6.11*10**((7.5*ta)/(237.7+ta))	#saturation vapour pressure (hpa)
		vp = svp*rh 							#actual vapour pressure (hpa)
		svp_m = 0.750061683 * svp 				#svp in  mmhg
		vp_m = svp_m*rh	
		edif=-1.694 * 10**-7  * 2430000 * (vp_m - pttbl(tsk))  #terms = permeance coefficient and vhw (J/kg)
		
		eres=0.0023*rm*(44-rh*pttbl(ta))
		ev=eres+ersw+edif
		
		# Dynamic changes to body
		alpha=0.044 + 0.35/ (skbf-0.1386) #mass proportion of skin shell to total body weight
		wt_sk = weight*alpha #skinshell weight 7% (gagge 5%, others 30%)
		wt_cr = weight-wt_sk
		
		# Dynamic changes to weather
		
		ta=weather[t][0]
		tmrt=weather[t][2]
		rh_per=weather[t][1]	
		v=weather[t][3]	
		rh=rh_per/100.0
		
		############forced -comment out when using weather file##########################################
		
		if (artificial_weather == 1):
			ta=21
			tmrt=21
			v=0.05
			rh=0.5
		#################################################################################################
		
		
		C =  fcl * chc * (ta-tsk)  #adu left out cos it will be accounted for
		sbc = 5.670373 * 10**-8
		feff = 0.7
		emiss = 0.97
		R =  fcl * feff * emiss * sbc * ((tmrt+273.15)**4 - (tsk+273.15)**4)
		
## RECORD VARIABLES AT END OF EACH TIME STEP
		store=rm-ev+R+C-wk	
		astorage.append(store)		
		atsk.append(tsk)
		atcr.append(tcr)		
		xaxis.append(time*time_steps*time_unit/60)
		aev.append(-ev)
		asignal.append(-skbf+mshiv)
		awetted.append(pwet)
		ata.append(ta)
		arh.append(rh)
		askbf.append(skbf)
		aC.append(C)
		aR.append(R)
		atmrt.append(tmrt)

		chr=5.23
		ctc=chr+chc		
		dry=ctc*(tsk-ta)*fcl
		
		print '[{t:3g}] R:{R:.1f} C:{C:.1f} Ed:{edif:.1f} Er:{eres:.1f} Es:{ersw:.1f} S:{store:.1f} H:{H:.1f} net:{net:.1f} skbf:{skbf:.1f}'\
		.format(t=t,R=R,C=C,edif=-edif,eres=-eres,ersw=-ersw, store=store,H=rm-wk,net=R+C-ev,dry=R+C,th=hfsk+hfcr,skbf=-skbf)
		
		period=period+1
			
	if age > 17:
		bb = "(adult:" + str(age) + "y),"
	elif age > 5:
		bb = "(child:" + str(age) + "y),"
	else:
		bb = "(infant:" + str(age) + "y),"
	
	plt.rcParams["axes.labelsize"] = 10
	plt.rcParams["xtick.labelsize"] = 10
	plt.rcParams["ytick.labelsize"] = 10
	plt.rcParams["legend.fontsize"] = 10
	
	#Plot 1: skin wettedness
	ax1=fig.add_subplot(521)
	ax1.plot(xaxis,awetted,colour)
	plt.xlim([0,time_steps*time_unit/60])
	plt.grid(True)
	#plt.ylim([0,1])
	#~ t=plt.title('Air temp: ' + str(ta)+'$^\circ$C' + ', RH: ' + str(rh*100) +'%')
	#~ box = ax1.get_position()
	#~ t.set_position([box.x0 + box.width*2, box.y0 + box.height*2.2,
    #~ box.width, box.height*0.9])
	plt.ylabel('Skin wettedness')
	
	#Plot 2: C, E
	ax3=fig.add_subplot(522)
	ax3.plot(xaxis,aC,colour+'-',xaxis,aev,colour+'--')
	label="\n".join(textwrap.wrap('C (solid)       E (dashed)', 14))
	plt.ylabel(label)
	plt.grid(True)
	plt.xlim([0,time_steps*time_unit/60])
	ax3.margins(0.04)
	
	#Plot 3: Air temp
	ax3=fig.add_subplot(523)
	ax3.plot(xaxis,ata,'k'+'-',xaxis,atmrt,'k--')
	label="\n".join(textwrap.wrap('Air temp (solid)       Tmrt (dashed)', 16))
	plt.ylabel(label)
	plt.grid(True)
	plt.xlim([0,time_steps*time_unit/60])
	plt.xlabel('Time (mins)')
	ax3.margins(0.04)
	
	#Plot 4: R
	ax2=fig.add_subplot(524)
	ax2.plot(xaxis,aR,colour+'-')
	plt.xlim([0,time_steps*time_unit/60])
	plt.ylabel('R')
	plt.grid(True)
	ax2.margins(0.04)
	
	#Plot 5: RH
	ax2=fig.add_subplot(525)
	ax2.plot(xaxis,arh,'k',label=str(age)+ "y, " + gender+ ", "+ str(met5))
	plt.xlim([0,time_steps*time_unit/60])
	plt.xlabel('Time (mins)')
	plt.ylabel('RH (ratio)')
	plt.grid(True)
	ax2.margins(0.04)
	
	#Plot 6: Thermoregulation
	ax2=fig.add_subplot(526)
	ax2.plot(xaxis,asignal,colour,label=str(age)+ "y, " + gender+ ", "+ str(met5))
	plt.xlim([0,time_steps*time_unit/60])
	plt.xlabel('Time (mins)')
	label="\n".join(textwrap.wrap('Thermoregulation (skbf + mshiv)', 18))
	plt.ylabel(label)
	plt.grid(True)
	ax2.margins(0.04)
	
	#plot#7: skin temp
	axlegend=plt.subplot(514)
	axlegend.plot(xaxis,atsk,colour)
	plt.ylabel('Skin temp. ($^\circ$C)')
	box = axlegend.get_position()
	axlegend.set_position([box.x0, box.y0 +box.height*0.05,box.width, box.height*0.9])
	plt.grid(True)
	plt.xlim([0,time_steps*time_unit/60])
	axlegend.margins(0.04)
	
	#plot#8: core temp
	axlegend=plt.subplot(515)
	axlegend.plot(xaxis,atcr,colour,label=str(age)+ "y, " + gender+ ", "+ str(met5))
	plt.ylabel('Core temp. ($^\circ$C)')
	plt.xlabel('Time (mins)')
	box = axlegend.get_position()
	axlegend.set_position([box.x0, box.y0 + box.height*0.15 , box.width, box.height*0.9])
	axlegend.legend(loc='upper center', bbox_to_anchor=(0.5, -0.35), fancybox=True,  ncol=5)
	plt.grid(True)
	plt.xlim([0,time_steps*time_unit/60])
	plt.ylim(36,38)
	axlegend.margins(0.04)
	
	#print 'Simulating subject' + str(ID) + '...'
	#~ line=0
	#~ f1 = open('output/' + str(ID) + '.csv', 'w')
	#~ writer = csv.writer(f1)
	#~ writer.writerow(['Timestep','Tsk','Tcr','Skin wettedness'])
	#~ while line < len(xaxis):
		#~ writer.writerow([xaxis[line],atsk[line],atcr[line], awetted[line]])
		#~ line+=1
	#~ f1.close()   
	
	#~ for i in [C,R,ersw,edif,eres]:
		#~ 
		#~ print i*adu_m,
	print rm-wk
	
main()






