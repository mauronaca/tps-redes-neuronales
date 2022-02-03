# tps-redes-neuronales
Trabajos Prácticos y Ejercicios de Redes Neuronales

Introducción

Las tecnologías actuales de detección de incendios tales como detectores de humo o sensores de calor no suelen ser precisas para espacios grandes o edificaciones complejas. Recientemente se avanzó en nuevas tecnologías más precisas y sensibles en cuanto a la detección basadas en imágenes. 

El proceso de detección basado en imágenes cuenta con 3 pasos importantes, que son el preprocesamiento de las imágenes, la extracción de características o feature extraction y la detección propiamente dicha. El núcleo del algoritmo es el feature extraction, y este paso se puede lograr con alta precisión utilizando redes neuronales convolucionales (CNN) que aprenden automáticamente a extraer características complejas de las imágenes, como pueden ser bordes de objetos o patrones distintivos del fuego. 

Las CNNs cuentan con una capa llamada convolutional layers. Estas capas son las más importantes, son filtros que realizan una transformación a la imagen generando un mapa de features. A los filtros utilizados se los llama kernel, de modo que las capas convolucionales son una serie de distintos kernels con el fin de encontrar un mapa de características de la imagen de entrada. 


Dataset

El dataset utilizado proviene del siguiente repositorio: https://github.com/OlafenwaMoses/FireNET
Link para descargar el dataset: https://drive.google.com/file/d/1dXDxqPJUMEaPeSwb5wl7kO1vpgQRllxW/view?usp=sharing
Otros datasets: https://cvpr.kmu.ac.kr/

Para crear los labels: https://cvat.org 

Contiene imagenes para validacion y entrenamiento de 259 x 194 pixeles. La estructura es la siguiente.
/fire-dataset/validation/images
/fire-dataset/validation/annotations
/fire-dataset/train/images
/fire-dataset/train/annotations

Hay 502 imágenes en total divididas en 412 para training y 90 para validation, esto representa un 82% y 18% respectivamente. 

En annotations se encuentran los archivos .xml asociados a cada imagen y en ellos se detalla la metadata del bounding box. Esto es muy importante ya que YOLO predice el bb y la probabilidad de la clase. En este caso hay una  sola clase que es fire pero puede existir una que sea smoke.

Un ejemplo de un .xml es el siguiente:

<annotation>
	<folder>images</folder>
	<filename>img (103).jpg</filename>
	<path>C:\Users\Moses Olafenwa\Documents\AI\Custom Datasets\FIRE DETECTION\dist\images\img (103).jpg</path>
	<source>
		<database>Unknown</database>
	</source>
	<size>
		<width>188</width>
		<height>268</height>
		<depth>3</depth>
	</size>
	<segmented>0</segmented>
	<object>
		<name>fire</name>
		<pose>Unspecified</pose>
		<truncated>0</truncated>
		<difficult>0</difficult>
		<bndbox>
			<xmin>38</xmin>
			<ymin>199</ymin>
			<xmax>84</xmax>
			<ymax>248</ymax>
		</bndbox>
	</object>
</annotation>

En él se especifica por un lado el tamaño de la imagen y en otra etiqueta el bounding box, en donde min y max corresponden a las coordenadas de la esquina inferior derecha y superior izquierda del bb.

El formato que YOLO interpreta es de la siguiente manera:

<object-class> <x_center> <y_center> <width> <height>

<object-class>: Nombre de la clase
<x_center> <y_center>: Posición del centro del bounding box normalizado al alto y ancho de la imagen
<width> <height>: Alto y ancho del bounding box normalizado.

Se usa el siguiente script para pasar del formato .xml al formato .tx ejecutando pytho xml2yolo.py dentro de la carpeta /fire-dataset/train||validation/annotation:



# -*- coding: utf-8 -*-

from xml.dom import minidom
import os
import glob

lut={}
lut["fire"] =0



def convert_coordinates(size, box):
    dw = 1.0/size[0]
    dh = 1.0/size[1]
    x = (box[0]+box[1])/2.0
    y = (box[2]+box[3])/2.0
    w = box[1]-box[0]
    h = box[3]-box[2]
    x = x*dw
    w = w*dw
    y = y*dh
    h = h*dh
    return (x,y,w,h)


def convert_xml2yolo( lut ):

    for fname in glob.glob("*.xml"):
        
        xmldoc = minidom.parse(fname)
        
        fname_out = (fname[:-4]+'.txt')

        with open(fname_out, "w") as f:

            itemlist = xmldoc.getElementsByTagName('object')
            size = xmldoc.getElementsByTagName('size')[0]
            width = int((size.getElementsByTagName('width')[0]).firstChild.data)
            height = int((size.getElementsByTagName('height')[0]).firstChild.data)

            for item in itemlist:
                # get class label
                classid =  (item.getElementsByTagName('name')[0]).firstChild.data
                if classid in lut:
                    label_str = str(lut[classid])
                else:
                    label_str = "-1"
                    print ("warning: label '%s' not in look-up table" % classid)

                # get bbox coordinates
                xmin = ((item.getElementsByTagName('bndbox')[0]).getElementsByTagName('xmin')[0]).firstChild.data
                ymin = ((item.getElementsByTagName('bndbox')[0]).getElementsByTagName('ymin')[0]).firstChild.data
                xmax = ((item.getElementsByTagName('bndbox')[0]).getElementsByTagName('xmax')[0]).firstChild.data
                ymax = ((item.getElementsByTagName('bndbox')[0]).getElementsByTagName('ymax')[0]).firstChild.data
                b = (float(xmin), float(xmax), float(ymin), float(ymax))
                bb = convert_coordinates((width,height), b)
                #print(bb)

                f.write(label_str + " " + " ".join([("%.6f" % a) for a in bb]) + '\n')

        print ("wrote %s" % fname_out)



def main():
    convert_xml2yolo( lut )


if __name__ == '__main__':
    main()







Entrenamiento de YOLOv4 utilizando la infraestructura de Darknet:

Links varios:

https://stats.stackexchange.com/questions/353607/on-yolo-and-its-loss-function

https://machinelearningspace.com/yolov3-tensorflow-2-part-1/

https://colab.research.google.com/github/luxonis/depthai-ml-training/blob/master/colab-notebooks/YoloV3_V4_tiny_training.ipynb#scrollTo=c7otuBJe4RLS

https://medium.com/@b117020/fire-detection-using-neural-networks-4d52c5cd55c5

https://github.com/AlexeyAB/darknet (v. no original que contiene la v.4 de yolo)
https://www.sciencedirect.com/science/article/pii/S2214157X2030085X (muchos datasets de ejemplo)

La arquitectura de la red tiene la siguiente forma:
yolo-net/
├─ darknet/
│  ├─ darknet
├─ fire-dataset/
│  ├─ test/
│  │  ├─ images/
│  │  ├─ test.txt
│  ├─ train/
│  │  ├─ images/
│  │  ├─ train.txt
│  ├─ validation/
│  │  ├─ images/
│  │  ├─ validation.txt
├─ obj.data
├─ obj.name
├─ yolov4-tiny.cfg

Dentro del dataset, en images/ se encuentran las imágenes correspondientes de entrenamiento o validación e incluye las anotaciones en un .txt en el formato de YOLO.
├─ fire-dataset/
│  ├─ test/
│  │  ├─ images/
│  │  │  ├─ img1.jpg
│  │  │  ├─ img1.txt
|  |  |  ├─ ..
|  |  |  ├─ ..
|  |  |  ├─ ..
│  │  ├─ test.txt



Dentro del directorio darknet/ se encuentra el código fuente para implementar la red YOLO, un archivo Makefile que permite compilar el código, y entre otros archivos más existe el ejecutable que se llama darknet, por el cual se puede entrenar la red y testear con imágenes de entrada entre otras funciones. Además dentro de este directorio se almacenan los resultados obtenidos al ejecutar la red, como por ejemplo los pesos finales y un gráfico que muestra la pérdida en función de la cantidad de iteraciones.


Preparación del Dataset

Al comienzo el dataset solo contiene las imágenes de training, testeo y validación y sus correspondientes anotaciones. Antes de poder comenzar con el entrenamiento se deben crear ciertos archivos que disponen de la información para que YOLO pueda saber donde se ubican las imágenes y sus correspondientes labels, así como también configurar los distintos parámetros de la red.

Agregar más imágenes de fuego en indoor y crear sus labels a partir de la red entrenada con el dataset original.
Agregar imágenes que luzcan similares a un fuego pero que no lo sean (sin label)
Agregar una clase de smoke e imágenes. 

YOLO necesita imágenes que no pertenezcan a la clase que se quiere predecir?

Short answer: No https://stackoverflow.com/questions/55202727/yolo-object-detection-include-images-that-do-not-contain-classes-to-be-predicte

Sin embargo es útil tener preparado un set de imágenes que puedan hacer confundir a la red que contiene un fuego dentro como por ejemplo una foto con el sol brillando, y dejarlo dentro del directorio /fire-dataset/train/images pero sin ningún label, de modo que al momento del entrenamiento la red sepa que esa imagen no contiene la clase de fuego. 


Notebooks

Para extraer los frames de videos: https://colab.research.google.com/drive/1Ja0RjGG_cEtuiGfwJ1cqTZ46HrsA1n8H#scrollTo=xnmVtYdYV4B9

Para dibujar el bounding box en una imagen:
https://colab.research.google.com/drive/1ksGp8PdPhfFbeN-VVooq0HR_perZ-Lxm

La red yolo basada en la arquitectura darknet: https://colab.research.google.com/drive/1cRzasViEzaMzgkQNxuWJMhSYepW14FOY#scrollTo=ODMXegZ2aZQt


