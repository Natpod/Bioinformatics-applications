
#!/usr/bin/awk -f

BEGIN {


FS=";";
OFS="\t\t"




countpos=0
countneg=0


countcogpos=0
counctcogneg=0
}

#POR CADA FILA DE DATOS EN EL SCRIPT
{

# Conteo de genes inducidos o reprimidos según el valor de expresión en segundo campo
if($2>0){
	countpos+=1
	}else{
	countneg+=1
	}

# Comandos que se cumplen para genes con un sólo COG asociado
if ($6 ~ "^.$"){
	
	# Para genes inducidos con 1 COG se crea un array asociativo (asociado a cog en campo 6). 

	# apos para el conteo de genes inducidos asociados a ese único COG
	# expos para ir sumando de valores de expresión asociados a ese único COG (toma la forma de un contador asociado a cada COG)

	if($2>0){
	apos[$6]+=1
	expos[$6]+=$2
	countcogpos+=1
	
	# Mismo procedimiento de arrays asociativos para genes reprimidos con 1 COG (asociado a cog en campo 6)
	
	}else{
   	if($2<0){
   		aneg[$6]+=1
		exneg[$6]+=$2
		countcogneg+=1
   		}
   		}
   		
		}
}


# En este conjunto de comandos los conteos ya se habrán realizado

END {

#Devolución de genes inducidos con primer contador countpos

print "# Induced genes altogether: "countpos

#Devolución de genes inducidos para el subconjunto con 1 cog con tercer contador countcogpos

print "# Induced genes with a single cog: "countcogpos
print "cogs distribution in induced genes:"
print "COGs","COGs_proportion","Average_express"

# Iteración por cada COG del array asociativo
for (cog in apos){
	# No. total genes inducidos EXCLUSIVAMENTE del cog que se itera
    	cog_value = apos[cog]
    	
	
        print cog, (cog_value / countcogpos), (expos[cog]/ apos[cog])
    	}
    	
 	
print "\n"

#Devolución de genes reprimidos con segundo contador countneg

print "# Repressed genes altogether: "countneg

#Devolución de genes reprimidos para el subconjunto con 1 cog con cuarto contador countcogneg

print "# Repressed genes with a single cog: "countcogneg
print "cogs distribution in induced genes:"
print "COGs","COGs_proportion","Average_express"

# Iteración y cálculo con genes reprimidos asociados a un cog, con sus propios arrays asociativos

for (cog in aneg){
    	cog_value = aneg[cog]

	# Mismo formato de salida que para los genes inducidos
        print cog, (cog_value / countcogneg), (exneg[cog]/ aneg[cog])
    	}
}