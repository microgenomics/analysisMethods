export LANG=en_US.UTF-8
set -e

##################PARSE PARAMETERS#############################

workband=0
cfileband=0
localband=0
statusband=0
realdataband=0
simdataband=0

function parsefastaFunction {
	echo 'BEGIN{FS="|"}
{
if($1~">"){
	split($1,array,">")
	$1=array[2]
	for (i=1;i<100;i++){
	band=0;
		if($i != ""){
				if($i==ID){
					ID=$(i+1);
					print ID
					exit 0
				}
		}
				
	}
}
}'
}

for i in "$@"
do
	case $i in
	"--workpath")
		workband=1
	;;
	"--cfile")
		cfileband=1
	;;
	"--simdata")
		simdataband=1
	;;
	"--realdata")
		realdataband=1
	;;
	"--help")
		echo "Usage: bash analysisMethods.bash --workpath [files directory] --cfile [config file] --simdata [simulation table csv format (',' separated)] --realdata [mconf from metasim]"
		echo -e"\nOptions explain:"
		echo "--workpath path where directory tree or your files are (included requirement files)"
		echo "--cfile configuration file"
		echo "--simdata csv with your simulation data obtain from Parse Module provide by SEPA, or your own csv file"
		echo -e "--realdata a file that contain the real values of reads distribution (.mprf of metasim, or some file with format ti-reads [or gi-reads])\n"
		echo "Methods Aviables: R2, RMS_NE, ROC_CURVE. Specify in the config file using the flag 'ANALYSISTYPE'"
		echo -e "For example: ANALYSISTYPE=RMS_NE,R2 \n"
		echo "See the README for more information"
		exit
	;;
	*)
		if [ $((workband)) -eq 1 ];then
			WORKPATH=$i
			EXECUTIONPATH=`pwd`
			statusband=$((statusband+1))
			workband=0
			
		fi
		
		if [ $((cfileband)) -eq 1 ];then
			for parameter in `awk '{print}' $i`
			do
				Pname=`echo "$parameter" |awk 'BEGIN{FS="="}{print $1}'`		
				case $Pname in
					"GENOMESIZEBALANCE")
						GENOMESIZEBALANCE=`echo "$parameter" | awk 'BEGIN{FS="="}{print $2}' | sed "s/,/ /g"`
						#echo "${parameters[$i]}"								
					;;
					"COMMUNITYCOMPLEX")
						COMMUNITYCOMPLEX=`echo "$parameter" | awk 'BEGIN{FS="="}{print $2}' | sed "s/,/ /g"`					
					;;
					"SPECIES")
						SPECIES=`echo "$parameter" | awk 'BEGIN{FS="="}{print $2}' | sed "s/,/ /g"`					
					;;
					"ABUNDANCE")
						ABUNDANCE=`echo "$parameter" | awk 'BEGIN{FS="="}{print $2}' | sed "s/,/ /g"`					
					;;
					"DOMINANCE")
					DOMINANCE=`echo "$parameter" | awk 'BEGIN{FS="="}{print $2}' | sed "s/,/ /g"`					
					;;
					"READSIZE")
						READSIZE=`echo "$parameter" | awk 'BEGIN{FS="="}{print $2}' | sed "s/,/ /g"`					
					;;
					"ABSENT")
						ABSENT=`echo "$parameter" | awk 'BEGIN{FS="="}{print $2}' | sed "s/,/ /g"`					
					;;
					"METHOD")
						METHOD=`echo "$parameter" | awk 'BEGIN{FS="="}{print $2}' | sed "s/,/ /g"`					
					;;
					"TIPERMANENT")
						TIPERMANENT=`echo "$parameter" | awk 'BEGIN{FS="="}{print $2}' | sed "s/,/ /g"`					
					;;
					"TOTALGENOMES")
						TOTALGENOMES=`echo "$parameter" | awk 'BEGIN{FS="="}{print $2}' | sed "s/,/ /g"`		
					;;
					"ANALYSISTYPE")
						ANALYSISTYPE=`echo "$parameter" | awk 'BEGIN{FS="="}{print $2}' | sed "s/,/ /g"`					
					;;
				esac
			done
			statusband=$((statusband+1))
			cfileband=0
		fi
		
		if [ $((realdataband)) -eq 1 ];then
			if [ -f $i ];then
				statusband=$((statusband+1))
				REALDATAFILE=$i
				realdataband=0
				#####################	FETCH ID REAL DATA	########################
				echo "checking your file"
				fields=`awk '{print NF;exit}' $REALDATAFILE`
				case $fields in
					"2")
						echo "file already formated, (cols: ti reads_assigned)"
					;;
					"3")
						echo "ti reads" >> rtmp
						while read line
						do
							ID=`echo "$line" |awk '{print $2}'`
							#parsed file
							case $ID in
								"ti")
									echo "$line" |awk '{print $3, $1}' >> rtmp				
								;;
								"gi")
									abu=`echo "$line" |awk '{print $1}'`
									gi=`echo "$line" |awk '{print $3}'`
									ti=""
									echo "fetching ti by gi: $gi"
									parsefastaFunction > parsefasta.awk
									while [ "$ti" == "" ]
									do
										gi=`curl -s "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$gi&rettype=fasta" |awk -v ID="gi" -f parsefasta.awk`
										ti=`curl -s "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=nuccore&db=taxonomy&id=$gi" |grep "<Id>"|tail -n1 |awk '{print $1}' |cut -d '>' -f 2 |cut -d '<' -f 1`
									done
									rm parsefasta.awk
									echo "$ti $abu" >> rtmp				
								;;
							esac
						done < $REALDATAFILE
						REALDATAFILE="rtmp"
						if [ "$WORKPATH" != "." ];then
							mv $REALDATAFILE ${WORKPATH}
						fi
						echo "Done"
					;;
					*)
						echo "unrecognized file"
						exit
					;;
				esac
				#####################################################################
			else
				echo "ERROR: $i doesn't exist"
				exit
			fi
		fi
		
		if [ $((simdataband)) -eq 1 ];then
			if [ -f $i ];then
				#this is csv file
				statusband=$((statusband+1))
				SIMDATAFILE=$i
				colti=`awk 'BEGIN{FS=","}{if(NR==1){for (i=1;i<=9;i++){if($i ~ "ti"){print i;exit}}}}' $SIMDATAFILE`
				if [ "$colti" == "" ];then
					echo "header 'ti' not found in $SIMDATAFILE, check your csv"
					exit
				else
					awk -v initial=$colti 'BEGIN{FS=","}{for (i=initial;i<=NF;i++){printf "%s ",$i}printf "\n"}' $SIMDATAFILE > stmp
					SIMDATAFILE="stmp"
				
				fi
				simdataband=0
				if [ "$WORKPATH" != "." ];then
					mv $SIMDATAFILE ${WORKPATH}	
				fi

			else
				echo "ERROR: $i doesn't exist"
				exit
			fi
		fi
	;;
	esac
done

if [ $((statusband)) -ge 4 ]; then

function PerdonazoFunction {
##########PERDONAZO METHOD FOR ABSENTS################
if [ "$ABSENT" == "yes" ]; then
	if [ "$TIPERMANENT" != "" ];then
		timayor=`awk 'BEGIN{mayor=-1;ti=1}{if($2>mayor){ti=$1;mayor=$2}}END{print ti}' ti_reads_tmp`
		#make sure you have tifamily.dat
		family=`curl -s "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=$timayor" |awk 'BEGIN{FS="[<|>]";prev=""}{if($2=="ScientificName"){prev=$3}if($3=="family"){printf "%s,",prev}}'`
		FAMILYPERMAMENT=`curl -s "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=$TIPERMANENT" |awk 'BEGIN{FS="[<|>]";prev=""}{if($2=="ScientificName"){prev=$3}if($3=="family"){printf "%s,",prev}}'`			
		if [ "$family" == "$FAMILYPERMANENT" ]; then
			sed -i '' "s/[[:<:]]$timayor[[:>:]]/$TIPERMANENT/g" ti_reads_tmp
			echo "--------------------perdonazo in $filename: yes"
		else
			echo "--------------------perdonazo in $filename: no"
		fi
	else
		echo "tax id permament is needed in configuration file (TIPERMAMENT=xxxxx)"
		rm htmp ti_reads_tmp
		exit
	fi
fi
######################################################
}

function R2function {
	echo "R2function called"
	echo -e "Analysis\nR2" > R2.dat
	totalcol=`awk '{print NF;exit}' $SIMDATAFILE`
	for coli in `seq 2 1 $totalcol`	#col 1 always be tax id, we begin in reads cols >=2
	do
		awk -v coli=$coli '{if(NR>1){print $1, $coli}else{print $coli > "htmp"}}' $SIMDATAFILE > ti_reads_tmp
		cat ti_reads_tmp
		PerdonazoFunction
		firstline=0
		
		while read line
		do
			if [ $((firstline)) -eq 0 ];then
				firstline=1
			else
				tir=`echo "$line" |awk '{print $1}'`
				readr=`echo "$line" | awk '{print $2}'`
				
				tis=`grep "$tir" ti_reads_tmp |awk '{print $1}'`
				reads=`grep "$tis" ti_reads_tmp | awk '{print $2}'`

				if [ "$tis" == "" ]; then
					echo "$readr 0" >> corr
				else
					echo "$reads $readr" >> corr
				fi

			fi
		done < $REALDATAFILE
		getR2Function > getR2.R
		Rscript getR2.R corr |grep "Multiple" |awk '{print $3}' |awk 'BEGIN{FS=","}{print $1}' > R2tmp
		cat htmp R2tmp > Rtmp2
		paste R2.dat Rtmp2 > ftmp
		mv ftmp R2.dat
		rm corr getR2.R
		
	done 
	rm ti_reads_tmp R2tmp htmp Rtmp2
}

function RMSfunction {
	echo "RMSfunction called"
	totalcol=`awk '{print NF;exit}' $SIMDATAFILE`
	echo -e "Analysis\nRRMSE\nAVGRE" > rms.dat
	for coli in `seq 2 1 $totalcol`	#col 1 always be tax id, we begin in reads cols >=2
	do
		awk -v coli=$coli '{if(NR>1){print $1, $coli}else{print $coli > "htmp"}}' $SIMDATAFILE > ti_reads_tmp
		PerdonazoFunction
		#in the next line, we find the ti that match in simulation data file.
		suma=0
		sumprom=0
		firstline=0
		while read line 
		do
			if [ $((firstline)) -eq 0 ];then
				firstline=1
			else
				ti=`echo $line |awk '{print $1}'`
				abu=`echo $line |awk '{print $2}'`
				backupsuma=$suma
				backupprom=$sumprom
				suma=`awk -v realti=$ti -v abu=$abu -v suma=$suma '{if($1==realti){if(abu==0){secure_abu=1}else{secure_abu=abu*2}re=(secure_abu-$2)/secure_abu ; print suma+(re*re);exit}}' ti_reads_tmp`
				sumprom=`awk -v realti=$ti -v abu=$abu -v suma=$sumprom '{if($1==realti){if(abu==0){secure_abu=1}else{secure_abu=abu*2}if((secure_abu-$2)>=0){print suma+((secure_abu-$2)/secure_abu);exit}else{print suma+(((secure_abu-$2)*-1)/secure_abu);exit}}}' ti_reads_tmp`
				if [ "$suma" == "" ];then
					suma=`echo "$backupsuma" |awk '{print $1+1}'`
				fi
				
				if [ "$sumprom" == "" ];then
					sumprom=`echo "$backupprom" |awk '{print $1+1}'`
				fi
			fi
		done < $REALDATAFILE	#mprf parsed file have 'abundance ti' format
		awk -v suma=$suma 'END{print sqrt(suma/NR)}' $REALDATAFILE > rrmsetmp
		awk -v suma=$sumprom 'END{print suma/NR}' $REALDATAFILE > avg

		cat rrmsetmp avg > 2rms
		cat htmp 2rms > rrmse
		paste rms.dat rrmse > rmsvalues
		mv rmsvalues rms.dat
	done
	
	rm ti_reads_tmp htmp rrmsetmp rrmse avg 2rms
}

function ROCfunction {
	
	echo "ROCfunction called"
	echo "fpr tpr" > ROCtmp.dat
	echo "file" > filerocname
	totalcol=`awk '{print NF;exit}' $SIMDATAFILE`
	for coli in `seq 2 1 $totalcol`	#col 1 always be tax id, we begin in reads cols >=2
	do
		awk -v coli=$coli '{if(NR>1){print $1, $coli}else{print $coli > "htmp"}}' $SIMDATAFILE > ti_reads_tmp
		filename=`cat htmp`
		TP=0	#true positive
		TN=0	#true negative
		FP=0	#false positive
		FN=0	#false negative

		PerdonazoFunction
		
		###########TRUE POSITIVE AND FALSE POSITIVE###########
		for tis in `awk '{print $1}' ti_reads_tmp` #$1 is ti
		do
			linetir=`grep -w "$tis" $REALDATAFILE |awk '{print $1, $2}'` #line="" make the script crash, awk print $1, $2; fix it

			if [ "$linetir" != "" ];then
				tir=`echo "$linetir" | awk '{print $1}'`
				readr=`echo "$linetir" | awk '{print $2}'`
				reads=`grep -w "$tir" ti_reads_tmp | awk '{print $2}'`
				resultado=`echo "$reads>=($readr/2)" | bc -l` #bash doesn't work with float
				if [ $((resultado)) -eq 1 ]; then
					TP=`echo  "$TP+1" |bc`
				fi
			else
				FP=`echo "$FP+1" |bc`	
			fi
			
			
		done
		
		###########FALSE NEGATIVE AND TRUE NEGATIVE###########

		for tir in `awk '{print $1}' $REALDATAFILE` #$1 is ti
		do	
			tis=`grep -w "$tir" ti_reads_tmp | awk '{print $1}'`
			if [ "$tis" == "" ]; then
				FN=`echo "$FN+1" |bc`
			fi								
		done
		TN=`wc -l ti_reads_tmp |awk '{print $1}' |awk -v total=$TOTALGENOMES -v tn=$TN '{print tn+(total-$1)}'`

		######################################################
		
		echo -e "\nTP: $TP FP:$FP FN: $FN TN:$TN"
		tpr=`echo "$TP $FN" | awk '{print $1/($1+$2)}'`
		fpr=`echo "$FP $TN" | awk '{print ($1/($1+$2))}'`
		echo "fpr: $fpr tpr: $tpr"
		echo "$filename" >> filerocname
		echo "$fpr $tpr" >> ROCtmp.dat
	done
	
	paste filerocname ROCtmp.dat > ROC.dat

	rm ti_reads_tmp htmp filerocname ROCtmp.dat
	
}
function getR2Function {
	echo 'args <-commandArgs()
myfile<-args[6]

correlaciones<-read.table(myfile,header=FALSE)

summary(lm(V1~V2, data=correlaciones))'
}
	for a in $ANALYSISTYPE
	do
		case $a in
			"R2")
				R2function
			;;
			"RMS_NE")
				RMSfunction
	   		;;
	   		"ROC_CURVE")
	   			ROCfunction
	   		;;
	   		*)
	   			echo "no method aviable for $a"
	   			exit
	   		;;
		esac
	done
			
	rm -f stmp rtmp
else
	echo "Invalid or Missing Parameters, print --help to see the options"
	echo "Usage: bash analysisMethods.bash --workpath [files directory] --cfile [config file] --simdata [table csv] --realdata [mconf from metasim]"
	exit
fi
