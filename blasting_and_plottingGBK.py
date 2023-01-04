#!/bin/env python
#---------------------------------EJERCICIO OPTATIVO PYTHON------------------------------------------------------#


# Programa que foormatee un número de genes (proporcionado en argumento) 
# que estén antes y después en el contexto del locus, o utilizando un pblast opcional

import sys
from Bio import SeqIO
from subprocess import Popen, PIPE
import time



# Espicificación de argumentos en línea de comando.
# El primer argumento (después del nombre del script) --> archivo fuente genebank para g_próximos
# El segundo argumento será el locus (query para el programa)
# El tercer argumento será el número de genes vecinos a cada lado del query

input_file = sys.argv[1]
locus = sys.argv[2]
inpnumframe = sys.argv[3]











# Función para print ccon color en los formatos de salida

def esc(code):
    """f\033[31;1;4m se corresponderá a un formato en rojo subrayado
     y f\033[0m a formato normal"""
    
    return f'\033[{code}m'






# Mensaje de ayuda-error que salta en errores

def messagehelp():

    print("~~~~~~~~~~~~~~\nUtilización errónea del archivo (solo ejecutable en linux)\n~~~~~~~~~~~~~~~. Debe de contener\n"+esc('31;1;4')+"- Contestación correcta en la interacción con la máquina\n- BLAST válido \n- En la línea de argumentos:"+ esc(0) +"\n\t nombre archivo \t query (aunque no se utilice) \t número de g_próximos a cada lado \n")
    
    # Salida con mensaje de error
    sys.exit(1)
    return







def make_list_and_dict(input_file):

    """
    Parseamos utilizando biopython un gbk para extraer primero
    una lista (list_of_genes) con todos los locus_tag ordenados
    segun su presencia en el genoma y por otro lado un diccionario
    (dict_of_genes) en el que almacenamos para cada locus_tag su 
    strand y su product (función)
    """

    list_of_genes = []
    dict_of_genes = {}

    # Abrimos el archivo genbank y guardamos la infomración de 
    #       genes, gen:strand y produc
    # En un diccionario y una lista respectivamente


    with open(input_file, "r") as input_handle:
        for record in SeqIO.parse(input_handle, "genbank"):
            for feature in record.features:
                if feature.type == 'CDS':
                    locus = feature.qualifiers['locus_tag'][0] 
                    list_of_genes.append(locus)


                    strand = feature.location.strand
                    try:
                        product = feature.qualifiers['product'][0] 
                    except:
                        product = "NA"
                    dict_of_genes[locus]=[strand,product]


    return  dict_of_genes,list_of_genes     







def get_genomic_neighbourhood_list(locus_tag,list_of_genes,  limite_extremos = inpnumframe):

    """Con esta funcion obtenemos los locus_tags de los genes
    que estan entre los rangos de posición [-n,+n], siendo inpnumframe la variable
    que define dicho n alrededor del gen (locus_tag) que usamos como query 
    (en la línea de argumentos)

    """
    genomic_neighbourhood_list=[]
    locus_tag_index = list_of_genes.index(locus_tag)

    # CONDICIÓN DE CONTROL Si se inserta un número de genes próximos por lado mayor del posible

    if locus_tag_index <= int(limite_extremos):
        
        print("No se podrá ejercer la lista al ser el tercer argumento mayor que el número de elementos de genes vecinos disponibles")
        print("Se efectuará la operación con el número máximo de genes vecinos próximos para ese locus")
        
        limite_extremos=locus_tag_index-1

    # Se itera para un range en número -limite de extremos a límite de extremos para obtener los genes vecinos

    for n in range(-int(limite_extremos), int(limite_extremos)+1):
                    
        genomic_neighbourhood_list.append(list_of_genes[locus_tag_index+n])


    return genomic_neighbourhood_list









def print_result(genomic_neighbourhood_list,dict_of_genes,locus_gene):

    """
    A partir del genomic_neighbourhood_list, que contiene los 
    locus_tag que rodean a nuestro locus_tag query, podemos 
    obtener para todos estos genes su strand, y su funcion 
    usando el diccionario dict_of_genes, y plotemaos 

    1) usando el método format.
    2) usando matplotlib

    En función de lo que nos indique el usuario

    """
    #Para que haya un primer espacio de línea

    print()

    # Variables que se imprimirán en la salida:

    # 1) out_gene : lista que contendrá los formatos en función de 
    #   qué hebra codificante ocupa cada gen

    # 2) products_output: El locus y la función de su proteína

    out_gene=[]
    products_output=["\n"]

    # MÉTODO FORMAT : usando estos strings con las llaves que usa el método
    # format podemos facilmente printear los genes

    strand_negative_format="<[ {} ]"
    strand_positive_format="[ {} ]>"
    

    
    # Diccionario para formatear cada string del gen encontrado con su producto

    # Unpacking en genes próximos de 
    # 1) STRAND para el formateo de la salida (+/-) 2) producto para salida gen:producto
    # Identificación de query para señalarlo en salida
    
    for genelocs in genomic_neighbourhood_list:

        strand_neighbour, producto_neighbour = dict_of_genes[genelocs]

        #Construcción de primera lista con genes en orden
        if strand_neighbour < 0 :

            out_gene.append(strand_negative_format.format(genelocs))

        elif strand_neighbour > 0 :

            out_gene.append(strand_positive_format.format(genelocs))

        #Construcción de segunda lista con los gen:productos
        products_output.append(genelocs+" : "+producto_neighbour)
    




    #ESPECIFICACIÓN DE SI EL PROGRAMA TOMARÁ UN FORMATO SIMPLE DE LA SALIDA O NO

    print('--------------------------------\n¿Formato simple o complejo?\nConteste SIMPLE o complejo\n\t(se recomienda usar el simple si el no. genes vecinos es grande)\n--------------------------------')
    formato = input().upper()
    
    if formato == "SIMPLE":
       
        # # # 1) PRINTING DE LA SALIDA FORMAT. PRINT ROJO SI SE TRATA DEL QUERY

        print() # Primer salto de línea

        # ~~~ Primera parte de la salida ~~~~

        for n_gene_loc in out_gene:

             if locus_gene in n_gene_loc:

                 print(esc('31;1;4') + str(n_gene_loc)+ esc(0) , end="\t")

             else:
                 print(n_gene_loc, end="\t")
        
        # ~~~ Segunda parte de la salida ~~~~

        for output_gen_prod in products_output:

            if locus_gene in output_gen_prod:
                print(esc('31;1;4') + str(output_gen_prod)+ esc(0) + "--> query" , end="\n")
                 
            else:
                print(output_gen_prod, end="\n")
    
    elif formato == "COMPLEJO":

        # # # 2) PRINTING DE LA SALIDA MATPLOTLIB

        # Uso de módulo matplotlib y numpy para la representación. Efectos de texto con path_effects

        import matplotlib.patheffects as path_effects
        import matplotlib.pyplot as plt
        import numpy as np

        # Fijamos un fondo rosa (facecolor) haciendo un call de la figura donde se representará el diagrama

        fig1 = plt.figure(1)
        rect = fig1.patch
        rect.set_facecolor('pink')

        # Para el display de la proteína query 

        font = {'family': 'serif',
        'weight': 'bold',
        'size': 8,
        }

        # Variable que usaré para el display de las proteínas y sus productos en cajas con bbox

        caja = dict(boxstyle="round",
                           ec=(1., 0.5, 0.5),
                           fc=(1., 0.8, 0.8),
                           )

        #1) CADENAS: Dibujo de las cadenas de dna como dos funciones de sen(xmesh) desfasadas

        # Herramienta linspace del módulo numpy para definir rango de representación en x 
        #   (1000 valores)
        # Utilizo numpy para calcular los valores de seno también

        x1=np.linspace(-8,32*np.pi*round(int(inpnumframe)/2),1000)
        y1=np.sin(x1)/20
        x2=-np.pi+x1


        plt.plot(x1,y1, color = "red", label = "Cadena positiva", linewidth= 2.5) # cadena +
        plt.plot(x2,y1, color = "blue", label = "Cadena negativa", linewidth= 2.5) # cadena -
        ax = plt.gca()

        #2) Interespaciado: 

        #   x: Variable que colocará cada texto en orden en el gráfico. Inicialmente xlim izq.
        #   y: A partir de donde aparecen productos. Va variando para que no se superpongan

        # También varían los límites de representación en eje x según el número de genes vecinos.
        # Más adelante se observa la variación del eje y

        ax.set_xlim((ax.get_xlim()[0]/2)*int(inpnumframe),(ax.get_xlim()[1]/2)*int(inpnumframe))

        interespaciadox=ax.get_xlim()[0]-(10*int(inpnumframe))
        interespaciadoy=-(0.1/2)*int(inpnumframe)

        # Flecha que define el orden de los locus en las cajas del producto

        ax.quiver(interespaciadox-8, interespaciadoy, 0, 60, scale=1, minlength=40)

        #3) Printing de cada locus y su función con serie de bucles anidados:

        # Se itera para la lista con los productos,
        # para la lista con la información gráfica de strands, 
        # y para cada uno de los locus disponibles en la lista de genes próximos

        # La condición final es para comprobar si los strings de los locus están en
        # la cadena del producto (str_genprod)/cadena con información sobre strand(n_gene_loc)

        # Se da información de: 1) quién es el query 2) productos y locus 3) cadena

        for str_genprod in products_output:

            for n_gene_loc in out_gene:

                for genelocs in genomic_neighbourhood_list:
                    
                    if genelocs in n_gene_loc and genelocs in str_genprod: 

                        #------------ GENES VECINOS

                        # Strand positiva
                        if "]>" in n_gene_loc and locus_gene not in n_gene_loc:

                            plt.text(interespaciadox, 0.5, str(genelocs), color = "red", fontdict= font) # locus

                            plt.text(0, interespaciadoy, str(str_genprod), color="darkred", bbox=caja) # gen: producto

                            plt.plot(

                            interespaciadox+(4*int(inpnumframe)),0.3, # espaciado (ajustado x inpnumframe)
                            ">", color = "red", mew=2, mec = "pink", ms=10 # formato display

                                    ) # cadena
                        
                        # Strand negativa
                        elif "<[" in n_gene_loc and locus_gene not in n_gene_loc:

                            plt.text(interespaciadox, 0.5, str(genelocs), color = "blue", fontdict = font) # locus

                            plt.text(0, interespaciadoy, str(str_genprod), color="blue", bbox=caja) # gen: producto

                            plt.plot(interespaciadox+(4*int(inpnumframe)),0.3,"<", color ="blue",mew=2, mec = "pink", ms=10) # cadena
                        

                        #------------ SI SE TRATA DEL LOCUS QUERY

                        if locus_gene in n_gene_loc:

                            # Delimitador del locus (limite derecho más grande por el tamaño del texto)

                            plt.axvspan(

                            interespaciadox-(3*int(inpnumframe)),
                            interespaciadox+(7*int(inpnumframe)), # Rango de localización aj
                            alpha=0.1, color="red",  # Transparencia y color
                            label="Proteína query"
                                        ) 
                            # Texto del locus

                            text = plt.text(
                                interespaciadox-3, 0.5, # Localizadores del texto locus

                                str(genelocs), # Texto a imprimir
                                color = "k", 
                                fontdict= font, # Fuente descrita al principio
                                path_effects=[  # Efectos del texto, sombra y glow rosa

                                    path_effects.Stroke(linewidth=3, foreground='pink'),
                                    path_effects.withSimplePatchShadow()
                                    ]
                                            )
                            # Strand positiva
                            if "]>" in n_gene_loc:

                                plt.text(0, interespaciadoy, str(str_genprod), color="darkred", bbox=caja) # gen: producto

                                plt.plot(interespaciadox+(4*int(inpnumframe)),0.3,">", color = "red", mew=2, mec = "pink", ms=10) # info cadena 
                            
                            # Strand negativa
                            if "<[" in n_gene_loc:

                                plt.text(0, interespaciadoy, str(str_genprod), color="darkblue", bbox=caja) # gen: producto

                                plt.plot(interespaciadox+(4*int(inpnumframe)),0.3,"<", color ="blue",mew=2, mec = "pink", ms=10) # info cadena



            interespaciadox = interespaciadox + ((23/2)*int(inpnumframe)) #Distancia de 23/2*numgenesvecinos entre genes
            interespaciadoy = interespaciadoy -(0.1*int(inpnumframe)) #Distancia de 28 entre genes
        
        # Ajuste de valores en el eje y según el número de genes vecinos en el input

        ax.set_ylim(-0.75*int(inpnumframe),0.75*int(inpnumframe))

        # Se quitan los ejes y se deja el superior para una representación más limpia

        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

        # Título y leyenda

        plt.title(
            'Locus vecinos y productos de\n'+str(locus_gene), #título
        
        fontdict=font, fontsize=16, #fuentes y tamaño

        path_effects=[
            path_effects.Stroke(linewidth=3, foreground='pink'),
            path_effects.withSimplePatchShadow()

                    ] #efectos

                )
        plt.legend()

        plt.show()

    else:
        messagehelp()

    return






################################ SI BLAST ############################################################

def locus_retriever_blast():

#-------------------- CREACION DE MULTIFASTA

#Cuando utilicemos blast, para facilidad de uso en el terminal,
#tendremos que utilizar el formato fasta (multifa)en el subject

# Por cada record y feature se imprimirá en formato fasta el locus de la proteína y su secuencia para el alineamiento posterior
# Las proteínas en formato de genebank ya están en su orientación normal de traducción. No hará falta usar el complemento reverse

# CREACIÓN DEL MULTIFASTA PARA EL SUBJECT
    multifa = open('multifasta_subj.fa', 'w')
    with open(input_fichero, "r") as input_handle:
        for record in SeqIO.parse(input_handle, "genbank"):
            for feature in record.features:
                if feature.type == 'CDS':
                    locus = feature.qualifiers['locus_tag'][0]
                    try:
                        sequence =  feature.qualifiers['translation'][0]
                        multifa.write(">"+locus+"\n"+sequence+"\n")
                    except:
                        protein="NA"
    multifa.close()

    # Doy tiempo para que se actualice el archivo en el directorio y se pueda referenciar
    time.sleep(4)

    #------------------------ BLAST

    """En caso de blast dé el resultado del query o ejecute error. 
    Devuelve el query que se ha contrastado para el archivo principal genebank"""

    # Proceso popen con el archivo fasta query que construimos en la rutina princiapl

    proceso = Popen(['blastp', '-query', "prot.fa" , '-subject', "multifasta_subj.fa", '-outfmt','6', '-out', "result.txt"], stdout=PIPE, stderr=PIPE ,shell=True)

    #    #Lectura de error y resultado

    error = proceso.stderr.read()
    resultado_blast = proceso.stdout.read()

    #    #Cierre de procesos en shell

    proceso.stderr.close()
    proceso.stdout.close()

    if not error:

        # Doy tiempo para que se actualice el archivo en el directorio y se pueda referenciar
        time.sleep(4)

        #awk
        cmd = ['awk', '-v', 'OFS="\t"', '$3 > 40 && $4 >= 51{print $2}']

        # Recuperación de locus tras someter al resultado blast al comando cmd 
        # en el terminal de GNU unix

        with open("result.txt", "r+") as f:
            p = Popen(cmd, stdin=f, stdout=PIPE)
            result = p.stdout.read()

        # Quitamos valores b o ' para quedarnos con los locus y \n aliterizados

        result=str(result).strip("'b")

        # Dividimos una vez la cadena con el primer \n que nos encontremos

        # (asumimos que el blast sólo encuentra un solo locus tras el filtrado)
        # Si hay varios locus resultado nos quedamos con el primero de ellos 
        #       (mejores valores de iden y coverage)

        locus_final=str(result).split('\\n',1)[0]

            
    #SI HAY ERROR EN EL BLAST   
    else:
        print(f"Se produjo el siguiente error: {error}")
        messagehelp()


    return locus_final










##########################################################################################

#------------------------ RUTINA PRINCIPAL DEL PROGRAMA ------------------------


# CONTROL DE ENTRADA DE ARGUMENTOS

if len(sys.argv) < 4:
    messagehelp()
else:
    try:

    #------------ Construcción del diccionario de genes y strand,product ---------------

        dict_of_genes,list_of_genes = make_list_and_dict(input_file)
    
    
    #----------------------INTERACCIÓN CON EL USUARIO-----------------------------------

    # ESPECIFICACIÓN DE SI EL PROGRAMA TOMARÁ EL MÉTODO DE BLAST O NO.

        print("\n------------Segundo argumento no es archivo fasta query\n¿Quiere realizar blast para la lista de genes próximos igualmente? \n\nConteste BLAST o NORMAL. No se preocupe de las mayúsculas, sólo se tomará el valor de la secuencia\n------------\n")
        method_used = input().upper()

        if method_used == "NORMAL":

            genomic_neighbourhood_list = get_genomic_neighbourhood_list(locus,list_of_genes)
            print_result(genomic_neighbourhood_list,dict_of_genes, locus)

        elif method_used == "BLAST":

            # SI DIRECTAMENTE SE PROPORCIONA EL FICHERO QUERY

            if locus.endswith(".fa") or locus.endswith(".faa"):

                # Transferencia de contenidos desde 
                # el archivo que se indica en el segundo argument0
                # hasta el archivo query que se someterá a BLAST

                archivoinputquery = open(str(locus), "r+")
                queryfile = open("prot.fa","w")
                queryfile.write(archivoinputquery.read())
                queryfile.close()

            # SI NO SE PROPORCIONA EN ARGUMENTOS DE ENTRADA PERO EL USUARIO LO DECIDE ASÍ

            else:
                # Interacción con el usuario para que escriba la proteína

                print("Indique la secuencia de la proteína\n")
                proteina = input().upper()

                # Creacion de prot.fa
                queryfile = open("prot.fa","w")
                queryfile.write(">\n"+proteina+"\n")
                queryfile.close()

            #-------------- SE DEFINE UN NUEVO LOCUS Y SE EJECUTA EL PROGRAMA --------------------
            locus = locus_retriever_blast()
            genomic_neighbourhood_list = get_genomic_neighbourhood_list(proteina,list_of_genes)
            print_result(genomic_neighbourhood_list,dict_of_genes, locus)


        # SI EL USUARIO NO PRODUCE UNA RESPUESTA RELEVANTE: MENSAJE ERROR Y EXIT 1
        else:
            messagehelp()
    
    # SI EL USUARIO NO INSERTA ARGUMENTOS DE ENTRADA ADECUADOS: MENSAJE ERROR Y EXIT 1
    except:
        messagehelp()
