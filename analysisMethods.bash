export LANG=en_US.UTF-8
set -ex

##################PARSE PARAMETERS#############################

workband=0
cfileband=0
localband=0
statusband=0
realdataband=0
simdataband=0
metaphlanWorkband=0

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

metawork=`echo "$@" |awk '{if($0 ~ "--metaphlanWork"){print 1}else{print 0}}'`
if [ $((metawork)) -eq 1 ];then
	metaphlanWorkband=1
fi

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
	"--metaphlanWork")
		metaphlanWorkband=1
	;;
	"--help")
		echo "Usage: bash analysisMethods.bash --workpath [files directory] --cfile [config file] --simdata [simulation table csv format (',' separated)] --realdata [mconf from metasim]"
		echo -e"\nOptions explain:"
		echo "--workpath path where directory tree or your files are (included requirement files)"
		echo "--cfile configuration file"
		echo "--simdata csv with your simulation data obtain from Parse Module provide by SEPA, or your own csv file"
		echo -e "--realdata a file that contain the real values of reads distribution (.mprf of metasim, or some file with format ti-reads [or gi-reads])\n"
		echo "Methods Aviables: R2, RMS_NE, ROC. Specify in the config file using the flag 'ANALYSISTYPE'"
		echo -e "For example: ANALYSISTYPE=RMS_NE,R2 \n"
		echo "For ROC_CURVE you must specify a TOTALGENOMES=[integer] flag"
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
						echo "file already formated"
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
						done < <(grep "" $REALDATAFILE)
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

				if [ $((metaphlanWorkband)) -eq 1 ]; then
					echo "metaphlan convertion working"

					#fetch ti for species name
					firstline=0
					while read line
					do
						if [ $((firstline)) -eq 0 ]; then
							firstline=1
							echo "name reads" > newrealtmp
						else
							ti=`echo "$line" |awk '{print $1}'`
							abu=`echo "$line" |awk '{print $2}'`
							echo "fetching lineage from ti: $ti"
							 #fetch the ti by ncbi api
							nofetch=""
							while [ "$nofetch" == "" ]
							do
								curl -s "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=$ti" > tmp.xml
								nofetch=`cat tmp.xml`
							done
							name=`awk 'BEGIN{FS="[<|>]"}{if($2=="ScientificName"){printf "%s\n", $3;exit}}' tmp.xml` #be careful with \n
							lineage=`awk 'BEGIN{FS="[<|>]";prev=""}{if($2=="ScientificName"){prev=$3}if($3=="superkingdom"){printf "%s,",prev}if($3=="phylum"){printf "%s,",prev}if($3=="class"){printf "%s,",prev}if($3=="order"){printf "%s,", prev}if($3=="family"){printf "%s,",prev}if($3=="genus"){printf "%s,",prev}if($3=="species"){printf "%s,",prev}}' tmp.xml`
							cand=`echo "$lineage" |awk '{if($0 ~ "Candidatus"){print "yes"}else{print "no"}}'`
							if [ "$cand" == "yes" ]; then
								name=`echo "$name" |sed "s/Candidatus //g"`
								echo "$name $abu" >> newrealtmp
							else
								name=`echo "$lineage" |awk 'BEGIN{FS=","}{print $7}'`
								echo "$name $abu" >> newrealtmp
							fi
							rm tmp.xml
						fi
					done < <(grep "" $REALDATAFILE)
					rm -f rtmp
					mv newrealtmp rtmp
					REALDATAFILE="rtmp"

					#Grouping same names
					firstline=0
					while read line
					do
						if [ $((firstline)) -eq 0 ]; then
							firstline=1
							echo "name reads" > tmp
						else
							name1=`echo "$line" |awk '{print $1}'`
							name2=`echo "$line" |awk '{print $2}'`
							band=`awk -v name1=$name1 -v name2=$name2 '{if($1==name1 && $2==name2){print "find";exit}}' tmp`
							if [ "$band" != "find" ];then
								suma=`awk -v name1=$name1 -v name2=$name2 'BEGIN{suma=0}{if($1==name1 && $2==name2){suma+=$3}}END{print suma}' $REALDATAFILE`
								echo "$name1 $name2 $suma" >> tmp
							fi
							
						fi
					done < <(grep "" $REALDATAFILE)
					rm rtmp
					mv tmp rtmp
				fi
			
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
				ORIGINALSNAME=`echo "$SIMDATAFILE" |rev |cut -d "/" -f 1 |rev`

				if [ $((metaphlanWorkband)) -eq 1 ]; then
					echo "metaphlan convertion working"
					colti=`awk 'BEGIN{FS=","}{if(NR==1){for (i=1;i<=9;i++){if($i ~ "Species"){print i;exit}}}}' $SIMDATAFILE`			
				else
					colti=`awk 'BEGIN{FS=","}{if(NR==1){for (i=1;i<=9;i++){if($i ~ "ti"){print i;exit}}}}' $SIMDATAFILE`
				fi

				if [ "$colti" == "" ];then
					echo "header 'ti' not found in $SIMDATAFILE, check your csv"
					exit
				else
					if [ $((metaphlanWorkband)) -eq 1 ]; then
						awk -v initial=$colti 'BEGIN{FS=","}{for (i=(initial-1);i<=NF;i++){printf "%s ",$i}printf "\n"}' $SIMDATAFILE > stmp
						SIMDATAFILE="stmp"	
					else
						awk -v initial=$colti 'BEGIN{FS=","}{for (i=initial;i<=NF;i++){printf "%s ",$i}printf "\n"}' $SIMDATAFILE > stmp
						SIMDATAFILE="stmp"							
					fi
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
	if [ $((metaphlanWorkband)) -eq 1 ]; then
		echo -e "Analysis\nR2" > R2.$ORIGINALSNAME
		totalcol=`awk '{print NF;exit}' $SIMDATAFILE`
		for coli in `seq 3 1 $totalcol`	#col 1 and 2 always be metaphlan name, we begin in reads cols >=3
		do
			awk -v coli=$coli '{if(NR>1){print $1, $2, $coli}else{print $coli > "htmp"}}' $SIMDATAFILE > name_reads_tmp
			firstline=0
			PerdonazoFunction

			while read line
			do
				names=""
				if [ $((firstline)) -eq 0 ];then
					firstline=1
				else
					name1=`echo "$line" |awk '{print $1}'`
					name2=`echo "$line" | awk '{print $2}'`
					readr=`echo "$line" | awk '{print $3}'`

					names=`awk -v name1=$name1 -v name2=$name2 '{if($1==name1 && $2==name2){print "find";exit}}' name_reads_tmp`
					reads=`awk -v name1=$name1 -v name2=$name2 '{if($1==name1 && $2==name2){print $3;exit}}' name_reads_tmp`

					if [ "$names" == "" ]; then
						echo "$readr 0" >> corr
					else
						echo "$reads $readr" >> corr
					fi
				fi
			done < <(grep "" $REALDATAFILE)
			getR2Function > getR2.R
			Rscript getR2.R corr |grep "Multiple" |awk '{print $3}' |awk 'BEGIN{FS=","}{print $1}' > R2tmp
			cat htmp R2tmp > Rtmp2
			paste R2.$ORIGINALSNAME Rtmp2 > ftmp
			mv ftmp R2.$ORIGINALSNAME
			rm getR2.R corr
		done 
		rm name_reads_tmp R2tmp htmp Rtmp2
	else
		echo -e "Analysis\nR2" > R2.$ORIGINALSNAME
		totalcol=`awk '{print NF;exit}' $SIMDATAFILE`
		for coli in `seq 2 1 $totalcol`	#col 1 always be tax id, we begin in reads cols >=2
		do
			awk -v coli=$coli '{if(NR>1){print $1, $coli}else{print $coli > "htmp"}}' $SIMDATAFILE > ti_reads_tmp
			firstline=0
			PerdonazoFunction

			while read line
			do
				tis=""
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
			done < <(grep "" $REALDATAFILE)
			getR2Function > getR2.R
			Rscript getR2.R corr |grep "Multiple" |awk '{print $3}' |awk 'BEGIN{FS=","}{print $1}' > R2tmp
			cat htmp R2tmp > Rtmp2
			paste R2.$ORIGINALSNAME Rtmp2 > ftmp
			mv ftmp R2.$ORIGINALSNAME
			rm getR2.R corr
			
		done 
		rm ti_reads_tmp R2tmp htmp Rtmp2
	fi
}

function RMSfunction {
	echo "RMSfunction called"

	if [ $((metaphlanWorkband)) -eq 1 ]; then
		totalcol=`awk '{print NF;exit}' $SIMDATAFILE`
		echo -e "Analysis\nRRMSE\nAVGRE" > RMS.$ORIGINALSNAME
		for coli in `seq 3 1 $totalcol`	#col 1 and 2 always be the name in metaphlan, we begin in reads cols >=3
		do
			awk -v coli=$coli '{if(NR>1){print $1, $2, $coli}else{print $coli > "htmp"}}' $SIMDATAFILE > name_reads_tmp
			PerdonazoFunction
			#in the next line, we find the ti that match in simulation data file.
			suma=0
			sumprom=0
			firstline=0
			PerdonazoFunction

			while read line 
			do
				if [ $((firstline)) -eq 0 ];then
					firstline=1
				else
					name1=`echo $line |awk '{print $1}'`
					name2=`echo $line |awk '{print $2}'`
					abu=`echo $line |awk '{print $3}'`
					backupsuma=$suma
					backupprom=$sumprom
					suma=`awk -v realname1=$name1 -v realname2=$name2 -v abu=$abu -v suma=$suma '{
																			if($1==realname1 && $2==realname2){
																				if(abu==0){
																					secure_abu=1 #to avoid the 0/0
																				}else{
																					secure_abu=abu*2
																				}
																				if((secure_abu-$3)>=0){
																					re=(secure_abu-$3)/secure_abu;
																					if(re>=1){
																						print suma+1
																					}else{
																						print (suma+(re*re));																				
																					}																				
																				}else{
																					re=((secure_abu-$3)*-1)/secure_abu;
																					if(re>=1){
																						print suma+1
																					}else{
																						print (suma+(re*re));																				
																					}
																				}
																				exit
																			}
																		}' name_reads_tmp`
					sumprom=`awk -v realname1=$name1 -v realname2=$name2 -v abu=$abu -v suma=$sumprom '{
																			if($1==realname1 && $2==realname2){
																					if(abu==0){
																						secure_abu=1
																					}else{
																						secure_abu=abu*2
																					}
																					if((secure_abu-$3)>=0){
																						re=(secure_abu-$3)/secure_abu
																						if(re>=1){
																							print suma+1
																						}else{
																							print suma+re;
																						}

																					}else{
																						re=((secure_abu-$3)*-1)/secure_abu
																						if(re>=1){
																							print suma+1
																						}else{
																							print suma+re;
																						}																				
																					}
																					exit
																				}
																			}' name_reads_tmp`
					if [ "$suma" == "" ];then
						suma=`echo "$backupsuma" |awk '{print $1+1}'`
					fi
			
					if [ "$sumprom" == "" ];then
						sumprom=`echo "$backupprom" |awk '{print $1+1}'`
					fi

				fi
			done < <(grep "" $REALDATAFILE)	#mprf parsed file have 'ti abundance' format
			awk -v suma=$suma 'END{print sqrt(suma/(NR-1))}' $REALDATAFILE > rrmsetmp
			awk -v suma=$sumprom 'END{print suma/(NR-1)}' $REALDATAFILE > avg

			cat rrmsetmp avg > 2rms
			cat htmp 2rms > rrmse
			paste RMS.$ORIGINALSNAME rrmse > rmsvalues
			mv rmsvalues RMS.$ORIGINALSNAME
		done
		
		rm name_reads_tmp htmp rrmsetmp rrmse avg 2rms	
	else
		totalcol=`awk '{print NF;exit}' $SIMDATAFILE`
		echo -e "Analysis\nRRMSE\nAVGRE" > RMS.$ORIGINALSNAME
		for coli in `seq 2 1 $totalcol`	#col 1 always be tax id, we begin in reads cols >=2
		do
			awk -v coli=$coli '{if(NR>1){print $1, $coli}else{print $coli > "htmp"}}' $SIMDATAFILE > ti_reads_tmp
			PerdonazoFunction
			#in the next line, we find the ti that match in simulation data file.
			suma=0
			sumprom=0
			firstline=0
			PerdonazoFunction

			while read line 
			do
				if [ $((firstline)) -eq 0 ];then
					firstline=1
				else
					ti=`echo $line |awk '{print $1}'`
					abu=`echo $line |awk '{print $2}'`
					backupsuma=$suma
					backupprom=$sumprom
					suma=`awk -v realti=$ti -v abu=$abu -v suma=$suma '{
																			if($1==realti){
																				if(abu==0){
																					secure_abu=1 #to avoid the 0/0
																				}else{
																					secure_abu=abu*2
																				}
																				if((secure_abu-$2)>=0){
																					re=(secure_abu-$2)/secure_abu;
																					if(re>=1){
																						print suma+1
																					}else{
																						print (suma+(re*re));																				
																					}																				
																				}else{
																					re=((secure_abu-$2)*-1)/secure_abu;
																					if(re>=1){
																						print suma+1
																					}else{
																						print (suma+(re*re));																				
																					}
																				}
																				exit
																			}
																		}' ti_reads_tmp`
					sumprom=`awk -v realti=$ti -v abu=$abu -v suma=$sumprom '{
																			if($1==realti){
																					if(abu==0){
																						secure_abu=1
																					}else{
																						secure_abu=abu*2
																					}
																					if((secure_abu-$2)>=0){
																						re=(secure_abu-$2)/secure_abu
																						if(re>=1){
																							print suma+1
																						}else{
																							print suma+re;
																						}

																					}else{
																						re=((secure_abu-$2)*-1)/secure_abu
																						if(re>=1){
																							print suma+1
																						}else{
																							print suma+re;
																						}																				
																					}
																					exit
																				}
																			}' ti_reads_tmp`
					if [ "$suma" == "" ];then
						suma=`echo "$backupsuma" |awk '{print $1+1}'`
					fi
			
					if [ "$sumprom" == "" ];then
						sumprom=`echo "$backupprom" |awk '{print $1+1}'`
					fi

				fi
			done < <(grep "" $REALDATAFILE)	#mprf parsed file have 'ti abundance' format
			awk -v suma=$suma 'END{print sqrt(suma/(NR-1))}' $REALDATAFILE > rrmsetmp
			awk -v suma=$sumprom 'END{print suma/(NR-1)}' $REALDATAFILE > avg

			cat rrmsetmp avg > 2rms
			cat htmp 2rms > rrmse
			paste RMS.$ORIGINALSNAME rrmse > rmsvalues
			mv rmsvalues RMS.$ORIGINALSNAME
		done
		
		rm ti_reads_tmp htmp rrmsetmp rrmse avg 2rms
	fi
}

function ROCfunction {
	
if [ "$TOTALGENOMES" == "" ]; then
	echo "you need to specify TOTALGENOMES flag, to calculate roc curves"
else
	echo "ROCfunction called"

	if [ $((metaphlanWorkband)) -eq 1 ]; then
		echo "fpr tpr" > ROCtmp.dat
		echo "file" > filerocname
		totalcol=`awk '{print NF;exit}' $SIMDATAFILE`
		for coli in `seq 2 1 $totalcol`	#col 1 always be tax id, we begin in reads cols >=2
		do
			awk -v coli=$coli '{if(NR>1){print $1, $2, $coli}else{print $coli > "htmp"}}' $SIMDATAFILE > name_reads_tmp
			filename=`cat htmp`
			TP=0	#true positive
			TN=0	#true negative
			FP=0	#false positive
			FN=0	#false negative

			PerdonazoFunction
			
			###########TRUE POSITIVE AND FALSE POSITIVE###########
			while read line
			do
				name1=`echo "$line" |awk '{print $1}'`
				name2=`echo "$line" |awk '{print $2}'`
				reads=`echo "$line" |awk '{print $3}'`
				linetir=`awk -v name1=$name1 -v name2=$name2 '{if($1==name1 && $2==name2){print $1, $2, $3;exit}}' $REALDATAFILE` #line="" make the script crash, awk print $1, $2; fix it

				if [ "$linetir" != "" ]; then
					name1=`echo "$linetir" |awk '{print $1}'`
					name2=`echo "$linetir" |awk '{print $2}'`
					readr=`echo "$linetir" |awk '{print $3}'`
					resultado=`echo "$reads $readr" |awk '{if($1 >= $2/2){print "1"}else{print "0"}}'` #bash doesn't work with float
					
					if [ "$resultado" == "1" ]; then
						TP=`echo  "$TP+1" |bc`
					else
						band=`echo "$reads $readr" |awk '{if($1==0 && $2>0){print "fn";exit}}'`
						if [ "$band" == "fn" ]; then
							FN=`echo "$FN+1" |bc`
						fi

					fi

				else
					FP=`echo "$FP+1" |bc`	
				fi
				
			done < <(grep "" name_reads_tmp)	#mprf parsed file have 'ti abundance' format
			
			###########FALSE NEGATIVE AND TRUE NEGATIVE###########

			while read line
			do	
				name1=`echo "$line" |awk '{print $1}'`
				name2=`echo "$line" |awk '{print $2}'`
				names=`awk -v name1=$name1 -v name2=$name2 '{if($1==name1 && $2==name2){print $1, $2;exit}}' name_reads_tmp`
				if [ "$names" == "" ]; then
					FN=`echo "$FN+1" |bc`
				fi								
			done < <(grep "" $REALDATAFILE)	#mprf parsed file have 'ti abundance' format
			TN=`wc -l name_reads_tmp |awk '{print $1}' |awk -v total=$TOTALGENOMES -v tn=$TN '{print tn+(total-$1)}'`

			######################################################
			
			echo -e "\nTP: $TP FP:$FP FN: $FN TN:$TN"
			tpr=`echo "$TP $FN" | awk '{print $1/($1+$2)}'`
			fpr=`echo "$FP $TN" | awk '{print ($1/($1+$2))}'`
			echo "fpr: $fpr tpr: $tpr"
			echo "$filename" >> filerocname
			echo "$fpr $tpr" >> ROCtmp.dat
		done
		
		paste filerocname ROCtmp.dat > ROC.$ORIGINALSNAME

		rm name_reads_tmp htmp filerocname ROCtmp.dat
	else
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

				if [ "$linetir" != "" ]; then
					tir=`echo "$linetir" |awk '{print $1}'`
					readr=`echo "$linetir" |awk '{print $2}'`
					reads=`grep -w "$tir" ti_reads_tmp |awk '{print $2}'`
					resultado=`echo "$reads $readr" |awk '{if($1 >= $2/2){print "1"}else{print "0"}}'` #bash doesn't work with float
					
					if [ "$resultado" == "1" ]; then
						TP=`echo  "$TP+1" |bc`
					else
						band=`echo "$reads $readr" |awk '{if($1==0 && $2>0){print "fn";exit}}'`
						if [ "$band" == "fn" ]; then
							FN=`echo "$FN+1" |bc`
						fi

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
		
		paste filerocname ROCtmp.dat > ROC.$ORIGINALSNAME

		rm ti_reads_tmp htmp filerocname ROCtmp.dat
	fi
fi

}
function getR2Function {
	echo 'args <-commandArgs()
myfile<-args[6]

correlaciones<-read.table(myfile,header=FALSE)

summary(lm(V1~V2, data=correlaciones))'
}

if [ $((statusband)) -ge 4 ]; then

	for a in $ANALYSISTYPE
	do
		case $a in
			"R2")
				R2function
			;;
			"RMS_NE")
				RMSfunction
	   		;;
	   		"ROC")
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
