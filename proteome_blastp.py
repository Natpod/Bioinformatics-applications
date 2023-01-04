#!/bin/bash
# Blastp analysis

# This script allow to performs blastp analysis using a urls_file to download the subject proteomes
# this file have to contain two columns, first column contains the species (identifier) for that proteome
# and the url. The script will generate a output project folder which contains two subfolders, data (here
# downloaded proteomes are stored and input fasta files) and results (here blast result, in addition a folder is created
# for every query protein to store blast result, aligment and trees for each protein)



#Function to print help in terminal
ayuda () {
   echo 'Blastp project bash v1.0'
   echo 'Natalia García Sánchez'
   echo -e '\nusage: ejercicio.sh <query_sequences.fa> <ncbi_urls_file> <output_folder> <blast-identity> <blast-coverage>\n'
   echo -e 'query_sequences.fa : a fasta/multifasta file containing protein query sequences for blast'
   echo -e 'ncbi_urls_file     : a text plain file containing species name and url to download fasta protein file'
   echo -e 'output_folder	   : folder in which data and results will be stored'
   echo -e 'blast-identity     : sequence identity cut off value 0-100'
   echo -e 'blast-coverage     : sequence coverage cut off value 0-100\n'
}


#Arguments assignation
query=$1
ncbi_urls_file=$2
project_name=$3
iden=$4 #70
cov=$5 #40


#-------------------------CONTROLES DE ENTRADA EN SCRIPT--------------------------------------------

#Controls help message. Función que necesita usarse con 1 argumento (argumento posicional 1)
saltaayuda(){
	case "$1" in
	#si se evoca la ayuda en el primer argumento, salta la ayuda y sale sin estado de salida de error del programa 
			-h|-help) ayuda; exit 0;;
	#para el resto de primeros argumentos posicionales, se sigue el flujo del programa
			*);;
	esac
}


#Control del número y tipo de argumentos

if [[ "$#" -ne 5 ]]
then
	if [[ "$#" -gt 5 ]]
	then
		# por si el primer argumento de muchos es help
		saltaayuda $1
		# por si no se evocó help
		echo -e "\nerror: too many arguments\n ----------------"
		ayuda
	elif [[ "$#" -lt 5 ]]
	then
		# por si el primer argumento de pocos es help
		saltaayuda $1
		# por si no se evocó help
		echo -e "\nerror: not enough arguments\n --------------"
		ayuda
	fi
#no se presentan el numero de argumentos que se debería y salimos con un estado de error del programa 
exit 1
fi



#-----------Control de la naturaleza de los argumentos de asignación

#----- Control del formato de coverage

# Redondeo de coverage en el caso de ser un decimal (por la propia naturaleza del formato de bash)	
# IMPORTANTE QUE ESTÉ EL FORMATO DE DECIMALES EN . Se asume que el número es positivo, y si no, se someterá al siguiente control igualmente

if [[ "$cov" =~ [0-9]+.[0-9]+? ]]; then
	cov=$(printf %.0f $(echo "$cov"))
fi

#----- Control de cifras de seq identity y coverage

if [[ $iden -gt 100 ]] ||  [[ $iden -lt 0 ]] || [[ $cov -gt 100 ]] || [[ $cov -lt 0 ]]
then
	echo -e "\nerror in usage: coverage or sequence identity cutoff values must be numeric strings ranging 0 to 100\n -------------------"
	#echo "$_" >>$project_name/practica.log
	ayuda
	exit 1
fi


#Control de las características del directorio de salida

if [[ -d $3 ]]
then
	echo -e "WARNING: \t $3 directory already exists. \n If you continue with the operation, the data contained in the directory will be lost...\n Write yes or no [y/n] to continue...\n------------"
	#STDIN para interacción con usuario
	read entrada
	#si no decide abortar la operación, nos movemos a dicho directorio para crear directorios parent(en caso de existir)
	case $entrada in
		Y|y|[Yy]es)rm -r $project_name; mkdir $project_name; cd $project_name; mkdir data;mkdir results;cd ..;;
		N|n|[Nn]o) exit 0;;
		*) echo "no relevant response, aborting operation..."; exit 1;;
	esac

else
	if [[ -f $3 ]]
	then
		# Si es archivo, se sale del script con estado de salida de error.
		echo "error in usage: output_folder is not a directory"
		exit 1
	
	else
		echo -e "\n creating output folders...\n"
		
		#Create project directories and output_files
		mkdir $project_name
		cd $project_name
		mkdir data
		mkdir results
		cd ..
	fi

fi

#Control de archivos con query y subject
if [[ -f $1 ]]
	then
	echo "readable query files"
	cp $1 $project_name/data
else
	# Si no existe el query se sale con un error.
	echo -e "error in usage: query FILE doesn't exist\n ------------"
	echo "$_" >>$project_name/practica.log
	ayuda
	exit 1
fi

#-----------------------DESCARGA Y CREACIÓN DE FICHEROS DICCIONARIOS-------------------------------------------------

#PARSEO DE ARCHIVO CON URLS:
#1. Creación de diccionario de identificadores: "Identifiers.tsv" 
#2. Descarga de PROTEOMAS
#3. Creación de diccionario especies-accesiones de proteinas

# 1. Parseo de proteomas_urls

echo -e "$0 is Running...\n"
echo "Fasta files will be stored in /$project_name/data/ (Downloading and decompressing files...)"

cat $ncbi_urls_file | while read line
do 
	
	linkadescargar=$(echo "$line" | cut -f2)
	archivodescargado=$(echo "$line" | cut -d/ -f11 | sed 's/.gz//')
	especie=$(echo "$line" | cut -f1)

	#Descarga de gz
	wget $linkadescargar -P $project_name/data 2>>$project_name/practica.log

	#Creación de identifiers
	echo -e "$especie\t$project_name/data/$archivodescargado" >>$project_name/data/Identifiers.tsv

done

# PASO INTERMEDIO--> Descompresión de archivos conteniendo proteomas
gzip -d $project_name/data/*.gz 2>>$project_name/practica.log




# 2. Creación del multifasta subject iterando y haciendo cat de archivos de proteomas
# Para ficheros de extensión faa|fasta (no .fa para no confundir con query.fa)

for file in $project_name/data/*.fa?*
do
	
	cat $file >>$project_name/data/proteome.fa

done

# 3. Generate a big file that contains for every species all the fasta headers, necesary for species assgination in blast hits

echo -e "\nGenerating subject file and species_fasta_ID.tsv file.. \nProceeding with subject file testing\n---------------------------------------"



#Parseo de cada proteoma en identifiers para asociarlos con su especie en el mismo archivo
	
cat $project_name/data/Identifiers.tsv | while read especieyarchivo
do
	#Primer campo: Especie
	#Segundo campo: Ruta al proteoma específico para la especie
	
	especie=$(echo "$especieyarchivo" | cut -f1)

	cat $(echo "$especieyarchivo" | cut -f2) | grep "^>" | cut -d " " -f1 | while read IDprot
	do
		
		echo -e "$especie\t$IDprot" >>$project_name/data/species_fasta_ID.tsv
	
	done
done


#Punto de control de archivos que se someterán a blastp. Ya se ha hecho el punto de control del query inicialmente
if [[ -f $project_name/data/proteome.fa ]]
	then
	
	echo -e "\nSubject file exists, \n - Press any key+Enter to continue with blast and muscle alignment \n - or wait to abort operation and keep the proteome \n\n multifasta \t dictionaries for blast/muscle(for use in another script)\n--------------------------------------"
	read -t 15 -p "Prompt in 15 seconds " RESP
	if [[ $? -gt 128 ]] ; then
    	echo -e "\nTimeout, aborting operation"
    	exit 1
	fi

else
	
	echo -e "error in program subject FILES no not exist\n ------------------------------------"
	echo "$_" >>$project_name/practica.log
	ayuda
	exit 1

fi



#------------------------------------BLAST SECTION-----------------------------------------------#

echo -e "\nRunning Blast analysis and result filtering... "

#---------------------#
#PRIMER RESULTADO BLAST. Aseguramos que el formato de salida contiene los ID de las proteínas pato-génicas

blastp -query $project_name/data/query.fa -subject $project_name/data/proteome.fa -outfmt "6 qseqid sseqid pident qcovs sseq" >$project_name/results/blast_result 2>>$project_name/practica.log

#--------------------------------#
#SEGUNDO RESULTADO BLAST FILTRADO. Filtrado para los valores de seqiden y coverage de los argumentos de entrada

cat $project_name/results/blast_result | awk '$3 >= '$iden' && $4 >= '$cov'{print $1"\t"">"$2"\t"$5}' >>$project_name/results/blast_result_filtered



#---------------------------------------------------------------#
#RESULTADOS CON RESPECTO A LOS READS DE CADA PROTEÍNA PATOGÉNICA

# Itero VARIABLE patos -describiendo proteínas pato-génicas- 
# y selecciono dichas proteinas del blast filtrado

# 1. CREACIÓN DIRECTORIO Y ARCHIVO CON proteínas (para iterarlo la próxima vez sin hacer uso de query.fa)
# 2. Almacenamiento de datos de la proteína en cuestión para cada directorio
	

cat $project_name/data/query.fa | grep "^>" | sort | while read patos
do
	
	#Quito el carácter >
	patos=$(echo "$patos" | cut -f1 | sed 's/>//')

	#1.
	echo "$patos" >>$project_name/results/uniq_query_list.txt
	
	cd $project_name/results
	mkdir $patos
	cd ../..

	#2.
	cat $project_name/results/blast_result_filtered | awk -v pat=$patos '$1 == pat{print $2"\t"$3}' >>$project_name/results/$patos/$patos.fa

done

#--------------------------------------------#
#TERCER RESULTADO BLAST FILTRADO CON ESPECIES. 

# PARA ELLO ITERAMOS en el blast filtrado por valor de identificador de proteína, 
# Hacemos grep con el archivo de especiesID y sacamos línea, de la que sacaremos la especie
# El sort permite que la operación sea más limpia y versátil.

cat $project_name/results/blast_result_filtered | sort -k2 | while read lineblast
do
	
	accfastablast=$(echo "$lineblast" | cut -f2)
	cat $project_name/data/species_fasta_ID.tsv | sort -k2 | grep "$accfastablast" | while read linespecies
	do
		
		especieblast=$(echo "$linespecies" | cut -f1)
		echo -e "$especieblast\t$lineblast" >>$project_name/results/blast_result_final
	
	done
done

#------------------------- MEJORA DE RESULTADOS DE READS POR PROTEÍNA PATOGÉNICA -----------------------------#

#Adición de especie en los datos de reads de cada proteína (tras el número de accesión).
# Esto es posible haciendo uso del blast_result_final.

#Itero de nuevo para cada proteína patogénica, los identificadores de proteína 
#y los utilizo con un grep sobre el fichero de blast con las especies 
#para hacer el printing final sobre el mismo archivo

cat $project_name/results/uniq_query_list.txt | while read patos
do
	
	accpat=$(cat $project_name/results/$patos/$patos.fa | sort -k2 | cut -f1)
	cat $project_name/results/blast_result_final | sort -k3 | grep "$accpat" | awk '{print $3"_"$1"\n"$4}' >$project_name/results/$patos/$patos.fa

done

# Contador que me servirá para crear un array de numero de reads para el printing final sobre el terminal.

noreadstotal=0

#--------------------------------------------MUSCLE SECTION----------------------------------------------------#

#Make MUSCLE Phylogenetic trees for each query protein
echo "Making MUSCLE alignment and phylogenetic trees..."


# Este comando nos servirá para almacenar las variables que ocurren en el proceso hijo de while que ocurren al usar un pipeline, para que no se borren las variables como sucede por defecto después de dejar el subshell del pipe

# SIN EMBARGO SOLAMENTE SERÁ FUNCIONAL EN VERSIONES DE BASH SUPERIORES A 4.2
# antes no lo necesitábamos porque almacenábamos los datos directamente en ficheros que no se eliminaban, sin embargo ahora vamos a necesitar salir con una variable de tipo contador

shopt -s lastpipe

# Iteramos por los archivos de cada proteína 
cat $project_name/results/uniq_query_list.txt | while read patos
do
	
	#alineamientos
	muscle -in $project_name/results/$patos/$patos.fa -out $project_name/results/$patos/$patos.aln 2>>$project_name/practica.log
	#arboles filogeneticos
	muscle -maketree -in $project_name/results/$patos/$patos.aln -out $project_name/results/$patos/$patos.aln.nw -cluster neighborjoining 2>>$project_name/practica.log
	
	#aprovecho esta iteración sobre cada dato de proteina patogénica para hacer un recuento de reads, seleccionando las filas con el carácter por el que empiezan los identificadores: ">""
	noreads=$(cat $project_name/results/$patos/$patos.fa | grep "^>" | wc -l)

	#CONTEO TOTAL DE READS
	noreadstotal=$(echo "$noreads" + "$noreadstotal"| bc)

done



###----------------------------------------TERMINAL PRINTING------------------------------------------------####
echo -e "\n*** Done!! :^) ***\n"
echo "Results are available at /$project_name/results/"
echo "Total Blast hits: $noreadstotal"
echo "hits were found for query proteins"

# Final stats to show in terminal 

# Contador para crear el array asociativo de estadśiticas de reads ahora lo uso para expresarlo
cat $project_name/results/uniq_query_list.txt | while read patos
do
	
	echo -n "$patos:"
	noreads=$(cat $project_name/results/$patos/$patos.fa | grep "^>" | wc -l)
	echo -n " $noreads;   "

done

echo -e "\n"
