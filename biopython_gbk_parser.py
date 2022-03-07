#Ejercicio voluntario Biopython
#Paula Robles López, 2021

import sys
from Bio import SeqIO


input_file = sys.argv[1]
locus = sys.argv[2]

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

	#abrimos el archivo y nos quedamos con los features que sean CDS
	#guardamos todos los genes y además en un diccionario aparte
	#guardamos su strand y la función del gen (si existe)
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

	#nos devuleve una lista de genes y un diccionario de genes 
	#siendo key el gen y los  values su funcion y strand)
	return 	dict_of_genes,list_of_genes



def get_genomic_neighbourhood_list(locus_tag,list_of_genes):

	"""
	Con esta funcion obtenemos los locus_tags de los genes
	que estan en la posicion -2,-1,+1 y +2 alrededor del
	gen (locus_tag) que usamos como query

	"""
	#primero identificamos el indice del locus para encontrar sus vecinos
	indice_locus=list_of_genes.index(locus_tag)
	indices=[-2,-1,0,+1,+2] #indices de los vecino y el locus de input
	genomic_neighbourhood_list=[]

	#si el locus es el primer o segundo elemento de la lista 
	#entonces cogemos una sublista de indices para que no nos de
	#genes del final de la lista como proximos, en caso de que interese 
	#(p.ej genoma circular) quitar condicional y descomentar el trozo
	#de código dentro del except
	if indice_locus==0:
		indices=indices[2::]
	elif indice_locus==1:
		indices=indices[1::]
	
	#cogemos los indices vecinos de nuestro locus de nuestra lista de genes
	#obtenemos una lista de genes vecinos que nos interesan
	for i in indices:
		try:
			new_indice=indice_locus+i
			genomic_neighbourhood_list.append(list_of_genes[new_indice])
		except:
			pass
			#en caso de querer que si damos un locus penultimo o ultimo
			#nos de los valores iniciales de la lista como próximos comentar
			#el condicional anterior y descomentar lo siguiente:
			#if list_of_genes[indice_locus]==list_of_genes[-1]:
			#	genomic_neighbourhood_list.append(list_of_genes[0])
			#	genomic_neighbourhood_list.append(list_of_genes[1])
			#	break

			#elif list_of_genes[indice_locus]==list_of_genes[-2]:
			#	genomic_neighbourhood_list.append(list_of_genes[0])

	return genomic_neighbourhood_list



def print_result(genomic_neighbourhood_list,dict_of_genes):

	"""
	A partir del genomic_neighbourhood_list, que contiene los
	locus_tag que rodean a nuestro locus_tag query, podemos
	obtener para todos estos genes su strand, y su funcion
	usando el diccionario dict_of_genes, y plotemaos usando el
	método format

	"""
	print("\n")
	descripciones=[]
	#para cada locus de nuestra lista de genes vecinos
	#si el strand es positivo o negativo hacemos un formateo
	#distinto del print
	#y guardamos el locus y su respectiva funcion tmb para imprimir
	for locusito in genomic_neighbourhood_list:
		if  dict_of_genes[locusito][0]==1:
			positive_strand="[ {} ]>".format(locusito)
			print(positive_strand,end=' ')

		else:
			negative_strand="<[ {} ]".format(locusito)
			print(negative_strand,end=' ')

		descripciones.append("{} : {}".format(locusito,dict_of_genes[locusito][1]))

	print("\n")
	for elemento in descripciones:
		print(elemento)
	print("\n")



#llamamos a las funciones
dict_of_genes,list_of_genes = make_list_and_dict(input_file)
genomic_neighbourhood_list = get_genomic_neighbourhood_list(locus,list_of_genes)
print_result(genomic_neighbourhood_list,dict_of_genes)
