#!/bin/env python
#---------------------------------EJERCICIO OPTATIVO PYTHON------------------------------------------------------#


# Programa que detecte sitios de restricción en secuencia importada (sin saltos de línea ni encabezados de FASTA), 
# devuelva el número de fragmentos y su longitud.
# Enzimas de restricción: AbcI y AbcII

#Importo el módulo que me sirve para el reconocimiento de patrones 
import re
# Módulo del sistema para hacer saltar excepciones
import sys

# Imprimo la secuencia que se encuentra en el archivo de texto sobre una variable
file1 = open("gene_sequence.txt","r")

# EL ARCHIVO DEBE DE SER FASTA. Por ello cojo a partir de la segunda línea, es decir, ignoro el encabezado
dna = file1.readlines()
#Cierro el archivo
file1.close()

encabezado=dna[0]
dna = dna[1:]

# Uno las cadenas que he separado para descartar la primera línea
separador=""
dna = separador.join(dna)

# Si no encuentra el encabezado en la secuencia total no será fasta
if not re.search(">", encabezado):
    print('Secuencia dada no es fasta o no es gen, se saldrá del sistema')
    sys.exit(1)

# Quito los carácteres de salto de línea
dna = dna.replace("\n", "")

# Comprobación de que solo tiene bases nitrogenadas convencionales
if re.search("[^ATGC]", dna):
    print('La secuencia no es estrictamente de dna, puede ser una proteína. El programa no funcionará')
    sys.exit(1)

# Contadores de restricciones de enzimas
restr_enz_1, restr_enz_2= 0, 0

# Creo funciones iterables sobre las que me serviré para hacer el diccionario de referencia de posiciones
matches1=re.finditer(r"A[A|T|G|C]TAAT",dna)
matches2=re.finditer(r"GC[A|G]{1}[A|T]{1}TG",dna)

# Diccionario que escriba las posiciones en las que se encuentra el inicio del patrón y a qué enzimas de restricción corresponden los patrones
iteradores={}

for m in matches1:
    restr_enz_1+=1
    iteradores[m.start()] = 'AbcI'

for h in matches2:
    restr_enz_2+=1
    iteradores[h.start()] = 'AbcII'

# Mensaje que devuelve si no encuentra los fragmentos de restricción
if not matches1 and not matches2:
    print('No se encontraron sitios de enzimas de restricción. Secuencia final: ' + str(dna))

# Fragmentos de restricción
numfrag= restr_enz_1+restr_enz_2+1

# Ordeno las posiciones en las que se ha reconocido el patrón y sustituyo el diccionario anterior con un diccionario con las claves ordenadas
posiciones=list(sorted(iteradores.keys()))

referencia=[]

for i in posiciones:
    enz=iteradores[i]
    referencia += [enz]

# DICCIONARIO FINAL ORDENADO: 
iteradores = dict(zip(posiciones, referencia))

# array donde meteré los fragmentos
frags=[]

# Slicing según avance el valor de ilast, que dependerá de la enzima de restricción (cortes en elementos de distinto orden)
ilast=0
for i in list(iteradores.keys()):
    if iteradores[i] == 'AbcI':
        finalfrag=i+3
        frags += [dna[ilast:finalfrag]]
        ilast=i+3
    elif iteradores[i] == 'AbcII':
        finalfrag=i+4
        frags += [dna[ilast:finalfrag]]
        ilast=i+4

# Escritura del último fragmento
frags += [dna[ilast:len(dna)]]

#--------------------ESCRITURA DE RESULTADOS EN EL TERMINAL--------------------------------------------------#
print('El número de fragmentos encontrados es: '+ str(numfrag))
print('\nNúmero de fragmentos sacados: '+ str(len(frags))+'\n')

for i in range(len(frags)):
    long=len(frags[i])
    ordenfrag=i+1
    print('Fragmento no. ' + str(ordenfrag) + ' es: ' + str(frags[i]) + ' y tiene una longitud de: ' + str(long) + '\n')

# Salimos del sistema sin excepción
exit(0)