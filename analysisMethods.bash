if [[ "$@" =~ "--debug" ]]; then
	set -ex
else
	set -e
fi


##################PARSE PARAMETERS#############################

workband=0
cfileband=0
localband=0
statusband=0
realdataband=0
simdataband=0
notiWorkband=0

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

function groupingFunction {
	echo 'args<-commandArgs()
	file<-c(args[6])
	headr<-c(args[7])
	if(headr=="T"){
			df<-read.csv(file, header = T)
			newdf<-aggregate(. ~ Species, df, FUN = sum)

	}else{
			df<-read.csv(file, header = F)
			colnames(df)<-c("COL1","COL2")
			newdf<-aggregate(. ~ COL1, df, FUN = sum)
	}
	write.table(newdf,file,row.names = F,quote = F)' > grp.R
}

metawork=`echo "$@" |awk '{if($0 ~ "--notiWork"){print 1}else{print 0}}'`
if [ $((metawork)) -eq 1 ];then
	notiWorkband=1
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
	"--notiWork")
		notiWorkband=1
	;;
	"--help")
		echo "Usage: bash analysisMethods.bash --workpath [files directory] --cfile [config file] --simdata [simulation table csv format (',' separated)] --realdata [mconf from metasim]"
		echo -e"\nOptions explain:"
		echo "--workpath path where directory tree or your files are (included requirement files)"
		echo "--cfile configuration file"
		echo "--simdata csv with your simulation data obtain from Parse Module provide by SEPA, or your own csv file (e.g. pathoscope_table.csv)"
		echo "--realdata a file that contain the real values of reads distribution (.mprf of metasim, or some file with format ti-reads [or gi-reads])"
		echo "--notiWork work when you haven't a file with ti (metaphlan, constrains and kraken)"
		echo "Methods Aviables: R2, RRMSE (normalized relative RMS error), ROC. Specify in the config file using the flag 'ANALYSISTYPE'"
		echo -e "For example: ANALYSISTYPE=RRMSE,R2\n"
		echo "For ROC_CURVE you must specify a TOTALGENOMES=[integer] flag"
		echo "if ABSENT=YES PERDONAZO method will apply automatically and TIPERMAMENT must be in config file too (TIPERMAMENT=xxxx where xxxx is a taxonomic id of your permament organism)"
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
					"READSIZE")
						READSIZE=`echo "$parameter" |awk 'BEGIN{FS="="}{print $2}' |sed "s/,/ /g"`					
					;;
					"ABSENT")
						ABSENT=`echo "$parameter" |awk 'BEGIN{FS="="}{print $2}' |sed "s/,/ /g"`					
					;;
					"METHOD")
						METHOD=`echo "$parameter" |awk 'BEGIN{FS="="}{print $2}' |sed "s/,/ /g"`					
					;;
					"TIPERMANENT")
						TIPERMANENT=`echo "$parameter" |awk 'BEGIN{FS="="}{print $2}' |sed "s/,/ /g"`					
					;;
					"TOTALGENOMES")
						TOTALGENOMES=`echo "$parameter" |awk 'BEGIN{FS="="}{print $2}' |sed "s/,/ /g"`		
					;;
					"ANALYSISTYPE")
						ANALYSISTYPE=`echo "$parameter" |awk 'BEGIN{FS="="}{print $2}' |sed "s/,/ /g"`					
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
				echo "checking your real data file"
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

				if [ $((notiWorkband)) -eq 1 ]; then
					echo "metaphlan,constrains,kraken convertion working for the real (simulated) data file"

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
							lineage=`echo $lineage |awk '{gsub("\\[|\\]","");print $0}'`
							cand=`echo "$lineage" |awk '{if($0 ~ "Candidatus"){print "YES"}else{print "NO"}}'`
							if [ "$cand" == "YES" ]; then
								name=`echo "$name" |sed "s/Candidatus //g"`
								echo "$name $abu" >> newrealtmp
							else
								name=`echo "$name" |awk '{gsub("\\[|\\]","");print $1, $2}'`
								echo "$name $abu" >> newrealtmp
							fi
							rm tmp.xml
						fi
					done < <(grep "" $REALDATAFILE)
					
					awk '{if(NR>1){print "\""$1" "$2"\""","$3}}' newrealtmp >rtmp
					rm -f newrealtmp
					REALDATAFILE="rtmp"
				
					#Grouping same names
					groupingFunction
					Rscript grp.R rtmp F

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
				BACKUPS=$i
				ORIGINALSNAME=`echo "$SIMDATAFILE" |rev |cut -d "/" -f 1 |rev`

				if [ $((notiWorkband)) -eq 1 ]; then
					echo "metaphlan,constrains,sigma convertion working"
					colti=`awk 'BEGIN{FS=","}{if(NR==1){for (i=1;i<=9;i++){if($i ~ "Name"){print i;exit}}}}' $SIMDATAFILE`			
				else
					colti=`awk 'BEGIN{FS=","}{if(NR==1){for (i=1;i<=9;i++){if($i ~ "ti"){print i;exit}}}}' $SIMDATAFILE`
				fi

				if [ "$colti" == "" ];then
					echo "header 'ti' not found in $SIMDATAFILE, check your csv"
					exit
				else
					#make sim data only with ti numbers
					awk -v initial=$colti -F"," '{for (i=initial;i<=NF;i++){printf "%s ",$i}printf "\n"}' $SIMDATAFILE > stmp
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

function getR2Function {
	echo 'args <-commandArgs()
myfile<-args[6]

correlaciones<-read.table(myfile,header=FALSE)

summary(lm(V1~V2, data=correlaciones))'
}

function PerdonazoFunction {
##########PERDONAZO METHOD FOR ABSENTS################
if [ "$ABSENT" == "YES" ]; then
	if [ "$TIPERMANENT" != "" ];then
		timayor=`awk 'BEGIN{mayor=-1;ti=1}{if($2>mayor){ti=$1;mayor=$2}}END{print ti}' ti_reads_tmp`
		#make sure you have tifamily.dat
		family=`curl -s "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=$timayor" |awk 'BEGIN{FS="[<|>]";prev=""}{if($2=="ScientificName"){prev=$3}if($3=="family"){printf "%s,",prev}}'`
		FAMILYPERMAMENT=`curl -s "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=$TIPERMANENT" |awk 'BEGIN{FS="[<|>]";prev=""}{if($2=="ScientificName"){prev=$3}if($3=="family"){printf "%s,",prev}}'`			
		if [ "$family" == "$FAMILYPERMANENT" ]; then
			sed "s/[[:<:]]$timayor[[:>:]]/$TIPERMANENT/g" ti_reads_tmp >tmp
			rm ti_reads_tmp
			mv tmp ti_reads_tmp
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
	#WARNING WRONG IMPLEMENTED
	folder=`pwd`
	echo "R2function called in $folder"
	if [ $((notiWorkband)) -eq 1 ]; then
		echo -e "Analysis\nR2" > R2.$ORIGINALSNAME
		totalcol=`awk '{print NF;exit}' $SIMDATAFILE`
		for coli in `seq 2 1 $totalcol`	#col 1 and 2 always be metaphlan name, we begin in reads cols >=3
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
					name2=`echo "$line" |awk '{print $2}'`
					readr=`echo "$line" |awk '{print $3}'`

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
					readr=`echo "$line" |awk '{print $2}'`

					tis=`grep "$tir" ti_reads_tmp |awk '{print $1}'`
					reads=`grep "$tis" ti_reads_tmp |awk '{print $2}'`

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
	folder=`pwd`
	echo "RRMSE function called in $folder"

	if [ $((notiWorkband)) -eq 1 ]; then
		totalcol=`awk '{print NF;exit}' $SIMDATAFILE`
		echo -e "Analysis,RRMSE,NAVGRE" > RMS.$ORIGINALSNAME
		for coli in `seq 2 1 $totalcol`	#col 1 and 2 always be the name in metaphlan,contrains and kraken, we begin in reads cols >=3
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

					suma=`awk -v realname1=$name1 -v realname2=$name2 -v abu=$abu -v suma=$suma 'function abs(v) {return v < 0 ? -v : v}{
																			if($1==realname1 && $2==realname2){
																				if(abu==0){
																					exit #to avoid the 0/0
																				}
																				re=abs(($3-abu)/abu);
																				print (suma+(re*re));										
																				exit
																			}
																		}' name_reads_tmp`

					sumprom=`awk -v realname1=$name1 -v realname2=$name2 -v abu=$abu -v suma=$sumprom 'function abs(v) {return v < 0 ? -v : v}{
																			if($1==realname1 && $2==realname2){
																					if(abu==0){
																						exit
																					}
																					re=abs(($3-abu)/abu);
																					print (suma+re);																						
																					exit
																				}
																			}' name_reads_tmp`
					if [ "$suma" == "" ];then
						suma=`echo "$backupsuma"`
					fi
			
					if [ "$sumprom" == "" ];then
						sumprom=`echo "$backupprom"`
					fi

				fi
			done < <(grep "" $REALDATAFILE)	#mprf parsed file have 'ti abundance' format

			awk -v suma=$suma 'END{print sqrt(suma/(NR-1))}' $REALDATAFILE > rrmsetmp
			awk -v suma=$sumprom 'END{print suma/(NR-1)}' $REALDATAFILE > avg

			paste -d "," rrmsetmp avg > 2rms
			paste -d "," htmp 2rms > rrmse
			cat RMS.$ORIGINALSNAME rrmse > rmsvalues
			mv rmsvalues RMS.$ORIGINALSNAME
		done
	else
		totalcol=`awk '{print NF;exit}' $SIMDATAFILE`
		echo -e "Analysis,RRMSE,NAVGRE" > RMS.$ORIGINALSNAME
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


					suma=`awk -v realti=$ti -v abu=$abu -v suma=$suma 'function abs(v) {return v < 0 ? -v : v}{
																			if($1==realti){ 
																				if(abu==0){
																					exit #to avoid the 0/0
																				}
																				re=abs(($2-abu)/abu);
																				print (suma+(re*re));
																				exit
																			}
																		}' ti_reads_tmp`
					sumprom=`awk -v realti=$ti -v abu=$abu -v suma=$sumprom 'function abs(v) {return v < 0 ? -v : v}{
																			if($1==realti){
																					if(abu==0){
																						exit
																					}
																					re=abs(($2-abu)/abu);
																					print (suma+re);
																					exit
																				}
																			}' ti_reads_tmp`
					if [ "$suma" == "" ];then
						suma=`echo "$backupsuma"`
					fi
			
					if [ "$sumprom" == "" ];then
						sumprom=`echo "$backupprom"`
					fi

				fi
			done < <(grep "" $REALDATAFILE)	#mprf parsed file have 'ti abundance' format

			awk -v suma=$suma 'END{print sqrt(suma/(NR-1))}' $REALDATAFILE  > rrmsetmp
			awk -v suma=$sumprom 'END{print suma/(NR-1)}' $REALDATAFILE > avg

			paste -d "," rrmsetmp avg > 2rms
			paste -d "," htmp 2rms > rrmse
			cat RMS.$ORIGINALSNAME rrmse > rmsvalues
			mv rmsvalues RMS.$ORIGINALSNAME
		done
	fi

	rm -f name_reads_tmp htmp rrmsetmp rrmse avg 2rms ti_reads_tmp

}

function ROCfunction {
re='^[0-9]+$'
if ! [[ "$TOTALGENOMES" =~ $re ]] ; then
	echo "you need to specify TOTALGENOMES flag with a valid number to calculate roc curves"
else
	folder=`pwd`
	echo "ROCfunction called in $folder"

	if [ $((notiWorkband)) -eq 1 ]; then
		echo "fpr,tpr" > ROCtmp.dat
		echo "file" > filerocname
		totalcol=`awk '{print NF;exit}' $SIMDATAFILE`
		for coli in `seq 1 1 $totalcol`	#col 1 and 2 always be name in metaphlan,constrains and kraken, we begin in reads cols >=2
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
						band=`echo "$reads $readr" |awk '{if($1==0 && $2>0){print "fp";exit}}'`
						if [ "$band" == "fp" ]; then
							FP=`echo "$FP+1" |bc`
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
			echo "$fpr,$tpr" >> ROCtmp.dat
		done
		
		paste -d "," filerocname ROCtmp.dat > ROC.$ORIGINALSNAME
		sed "2d" ROC.$ORIGINALSNAME > tmp
		rm ROC.$ORIGINALSNAME
		mv tmp ROC.$ORIGINALSNAME

		rm name_reads_tmp htmp filerocname ROCtmp.dat
	else
		echo "fpr,tpr" > ROCtmp.dat
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
						band=`echo "$reads $readr" |awk '{if($1==0 && $2>0){print "fp";exit}}'`
						if [ "$band" == "fp" ]; then
							FP=`echo "$FP+1" |bc`
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
			echo "$fpr,$tpr" >> ROCtmp.dat
		done
		
		paste -d "," filerocname ROCtmp.dat > ROC.$ORIGINALSNAME
		rm ti_reads_tmp htmp filerocname ROCtmp.dat
	fi
fi

}

function SIMOBSfunction {
	folder=`pwd`
	echo "Simulation vs Observed function called in: $folder"

	if [ $((notiWorkband)) -eq 0 ]; then
		echo "convertion ti x lineage working"

		#fetch ti for species name
		###############################FOR REAL DATA
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
				lineage=`echo $lineage |awk '{gsub("\\[|\\]","");print $0}'`
				cand=`echo "$lineage" |awk '{if($0 ~ "Candidatus"){print "YES"}else{print "NO"}}'`
				if [ "$cand" == "YES" ]; then
					name=`echo "$name" |sed "s/Candidatus //g"`
					echo "$name $abu" >> newrealtmp
				else
					name=`echo "$name" |awk '{gsub("\\[|\\]","");print $1, $2}'`
					echo "$name $abu" >> newrealtmp
				fi
				rm tmp.xml
			fi
		done < <(grep "" $REALDATAFILE)
		awk '{if(NR>1){print "\""$1" "$2"\""","$3}}' newrealtmp >rtmp
		rm -f newrealtmp
		REALDATAFILE="rtmp"
		
		#Grouping same names
		echo "Grouping real otus"
		groupingFunction
		Rscript grp.R rtmp F

		####################################FOR OBSERVED DATA (SOFTWARE DETECTION)
		#fetch ti for species name
		awk 'BEGIN{FS=","}{printf "\"%s\",",$7 ;for (i=10;i<NF;i++){printf "%s,",$i}printf "%s\n",$NF}' $BACKUPS > stmp
		SIMDATAFILE="stmp"
		#Grouping same names

		echo "Grouping otus simulated otus"
		groupingFunction
		Rscript grp.R stmp T
		sed "s/\.reads/-reads/g" stmp > tmp
		mv tmp stmp
		rm grp.R

	fi

	totalcol=`awk '{if(NR>1){print NF;exit}}' $SIMDATAFILE`
	awk '{if(NR>1){print $1" "$2","$3}else{print "OTU,REAL"}}' $REALDATAFILE > SIMOBS.$ORIGINALSNAME

	for coli in `seq 3 1 $totalcol`	#col 1 and 2 always be the name we begin in reads cols >=3
	do
		awk -v coli=$coli '{if(NR>1){print $1, $2, $coli}else{print $(coli-1) > "htmp"}}' $SIMDATAFILE > name_reads_tmp

		firstline=0
		PerdonazoFunction
		while read line
		do
			if [ $((firstline)) -eq 0 ];then
				firstline=1
			else
				name1=`echo "$line" |awk '{print $1}'`
				name2=`echo "$line" |awk '{print $2}'`
				readr=`echo "$line" |awk '{print $3}'`
				reads=`awk -v name1=$name1 -v name2=$name2 '{if($1==name1 && $2==name2){print $3;exit}}' name_reads_tmp`
				if [ "$reads" == "" ]; then
					echo "0" >> simobstmp
				else
					echo "$reads" >> simobstmp
				fi
			fi
		done < <(grep "" $REALDATAFILE)	#mprf parsed file have 'ti abundance' format
		cat htmp simobstmp > header
		paste -d "," SIMOBS.$ORIGINALSNAME header > filetmp
		mv filetmp SIMOBS.$ORIGINALSNAME
		rm simobstmp header
	done
	rm htmp name_reads_tmp
}


if [ $((statusband)) -ge 4 ]; then

	for a in $ANALYSISTYPE
	do
		case $a in
			"R2")
				R2function
			;;
			"RRMSE")
				RMSfunction
	   		;;
	   		"SIMOBS")
				SIMOBSfunction
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
