from __future__ import absolute_import, division, print_function
import itertools
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import re
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast import NCBIWWW
from time import time
import os
import tensorflow as tf
from tensorflow import keras
import numpy as np
#import matplotlib.pyplot as plt

print(tf.__version__)
'''
path = '/Volumes/rap/igls/aves/ncbi'

files = []
#r=root, d=directories, f = files
for r, d, f in os.walk(path):
       for file in f:
            if '.fna' in file:
                files.append(os.path.join(r, file))
                print(file)
'''
files = ['/Users/franciscogambondeza/.mounty/rap/igls/reptiles/ncbi/Podarcis_lilfordi/GCA_947686815.1_rPodLil1.2_genomic.fna']
# ii = open('/Volumes/rap/igls/aves/ncbi/V_genes_ncbi_total.fasta','w')
for file in files:
    corte = file.split('/')
    animal = corte[-2]
    especie = file.split(corte[-1])[0] + 'Vs-e'
    if os.path.exists(especie) is False:
        os.mkdir(especie)
    abys = file
    print(animal, file)
    query = '/Users/franciscogambondeza/.mounty/rap/reposteria/aaa/reptiles/query-rep-Vs.txt'
    bas = especie + '/aexones.fasta'
    nnfile = especie + '/nexones.fasta'
    tablafile = especie + '/tabla.txt'

    vvv = "tblastn"
    outfile = tablafile
    blastn_cline = NcbiblastxCommandline(vvv, query=query, subject=abys, evalue=0.1, out=outfile, outfmt=6)
    stdout, stderr = blastn_cline()

    instances = ["YYC", "YFC", "YLC", "YIC", "YHC", "TFC"]
    aa = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

    # tabla obtenida de internet, se obtienen los datos necesarios
    A = open(tablafile, 'r')
    # wgs = ''

    ide = []
    contigs = []
    for line in A:
        par = line.split('\t')
        if '#' not in line:

            ID = par[1]
            a = par[8]
            b = par[9]
            if int(b) - int(a) >= 0:
                frame = 1
                start = a

            else:
                frame = -1
                start = b

            # start = line.split('\t')[10]
            hit = ID, frame, start
            contigs.append(hit)
            if ID not in ide:
                ide.append(ID)

    # Lectura de los contigs. busqueda de exones de tamaño adecuado. Se obtiene archivo de aminoacidos y nucleotidos (bas, nnfile)
    ofile = open(bas, 'w')
    oofile = open(nnfile, 'w')

    for record in SeqIO.parse(file, "fasta"):
        ee = record.id
        f = 0
        r = 0
        if ee in ide:
            seq = record.seq
            seq = seq.upper()
            for m, n, p in contigs:
                if int(n) >= 0 and m == ee:
                    yy = int(p) - 100
                    k = seq[yy:yy + 1000]
                    ini1 = [i.start() for i in re.finditer('AG', str(k))]
                    fin1 = [i.start() for i in re.finditer('CA', str(k))]
                    se = [(i, j) for i, j in itertools.product(ini1, fin1) if j > i + 240 and j < i + 390]
                    DF = []
                    for r, s in se:
                        x = r + yy
                        y = s + yy
                        pis = x, y
                        if pis not in DF:
                            DF.append(pis)
                            exon = seq[x + 2:y]
                            tras = exon[2:]
                            p = tras.translate(to_stop=True)
                            if len(p) >= 80:
                                for xx in instances:
                                    if xx in p[-35:]:
                                        eo = ee.split('-')[0]
                                        nrec = SeqRecord(exon, id=animal + '-' + eo + '-' + str(x + 4) + '-' + str(y),
                                                         description='plus')
                                        arec = SeqRecord(p, id=animal + '-' + eo + '-' + str(x + 4) + '-' + str(y),
                                                         description='plus')
                                        ofile.write(arec.format('fasta'))
                                        oofile.write(nrec.format('fasta'))
                                        f += 1
            print(animal, '----> ', record.id, len(record.seq))
            print('foward squences posibles ', f)
            seq = record.seq.reverse_complement()
            seq = seq.upper()
            ancho = len(seq)

            for m, n, p in contigs:
                if int(n) <= 0 and m == ee:
                    pp = ancho - int(p)
                    yy = int(pp) - 1000
                    k = seq[yy:yy + 2000]
                    ini1 = [i.start() for i in re.finditer('AG', str(k))]
                    fin1 = [i.start() for i in re.finditer('CA', str(k))]
                    se = [(i, j) for i, j in itertools.product(ini1, fin1) if j > i + 240 and j < i + 390]
                    DR = []
                    for r, s in se:
                        x = r + yy
                        y = s + yy
                        pis = x, y
                        if pis not in DR:
                            DR.append(pis)
                            exon = seq[x + 2:y]
                            tras = exon[2:]
                            p = tras.translate(to_stop=True)
                            if len(p) >= 80:
                                for xx in instances:
                                    if xx in p[-35:]:
                                        eo = ee.split('-')[0]
                                        aa = ancho - y
                                        bb = ancho - x
                                        nrec = SeqRecord(exon, id=animal + '-' + eo + '-' + str(aa) + '-' + str(bb),
                                                         description='minus')
                                        arec = SeqRecord(p, id=animal + '-' + eo + '-' + str(aa) + '-' + str(bb),
                                                         description='minus')
                                        ofile.write(arec.format('fasta'))
                                        oofile.write(nrec.format('fasta'))
                                        r += 1

            print('reverse squences posibles ', r)

    ofile.close()
    oofile.close()

    aa = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    Vs = ['back', 'ighv','igkv','iglv','trav','trbv','trgv','trev','igsv']

    x_train = np.load('/Users/franciscogambondeza/.mounty/rap/reposteria/aaa/reptiles/Vs/x_train.npy')
    y_train = np.load('/Users/franciscogambondeza/.mounty/rap/reposteria/aaa/reptiles/Vs/y_train.npy')
    x_test = np.load('/Users/franciscogambondeza/.mounty/rap/reposteria/aaa/reptiles/Vs/x_test.npy')
    y_test = np.load('/Users/franciscogambondeza/.mounty/rap/reposteria/aaa/reptiles/Vs/y_test.npy')

    model = tf.keras.Sequential([
        tf.keras.layers.Dense(128, activation=tf.nn.relu),
        tf.keras.layers.Dropout(0.2),
        tf.keras.layers.Dense(len(Vs), activation=tf.nn.softmax)
    ])

    model.compile(optimizer='rmsprop',
                  loss='sparse_categorical_crossentropy',
                  metrics=['accuracy'])

    model.fit(x_train, y_train, epochs=5)

    test_loss, test_acc = model.evaluate(x_train, y_train)

    print('Test accuracy:', test_acc)

    # os.mkdir('operatoria/ancient-ab/'+especie+'/pos-par')

    ID = []
    X0 = []
    contigs = []
    count = 0
    for record in SeqIO.parse(bas, "fasta"):
        ID.append(record.description)
        par = record.id.split('-')
        if par[1] not in contigs:
            contigs.append(par[1])
        seq = record.seq[:40] + record.seq[-40:]
        pp = record.seq
        fre = []
        for x in range(len(seq)):
            a = seq[x]
            for y in aa:
                if a == y:
                    fre.append(1)
                else:
                    fre.append(0)
        for x in range(len(aa)):
            for y in range(len(aa)):
                i = aa[x] + aa[y]
                q = pp.count(i)
                fre.append(q)
        X0.append(fre)
        count += 1
    print('secuencias  ', count)
    X = np.asarray(X0)
    print(X.shape)

    X0 = np.array(X, dtype=np.int)
    print(X0.shape)
    predictions = model.predict(X0)
    print(predictions[0])
    print('transformacion numerica heha')

    oo = open(especie + '/aexones_asignados.txt', 'w')
    nlp = []
    for a in range(len(predictions)):
        pp = round(max(predictions[a]), 3)
        ale = int(np.argmax(predictions[a]))
        locus = Vs[ale]
        p = ID[a] + '-----------> ' + str(pp) + '-' + Vs[ale] + '\n'
        oo.write(p)
        if locus != 'back':
            name = ID[a]
            data = name, locus, pp
            nlp.append(data)
    oo.close()
    print('asignacion hecha')
    dupli = []
    elegidos = []
    for j in contigs:
        RR = []
        for x, y, z in nlp:
            por = x.split(' ')
            par = por[0].split('-')
            if j == par[1]:
                ran = int(par[2]), int(par[3]), y, z
                RR.append(ran)
        R = sorted(RR, key=lambda x: x[0])
        for a, b, c, d in R:
            if a not in dupli:
                k = range(a - 1, b + 300)
                segmento = []
                for m, bb, cc, dd in R:
                    if m in k and m not in dupli:
                        h = m, bb, cc, dd
                        segmento.append(h)
                        dupli.append(m)
                if len(segmento) >= 1:
                    w = max(segmento, key=lambda x: x[-1])
                    elegidos.append(w)


    print('Secuencias totales ', len(nlp))
    print('Secuencias elegidas ', len(elegidos))

    # obtenemos las secuencias
    ccc = 0
    o = open(especie + '/result-vs.fasta', 'w')

    ll = []
    for rio in SeqIO.parse(bas, 'fasta'):
        par = rio.id.split('-')
        pir = rio.description.split(' ')
        for a, b, c, d in elegidos:
            if int(a) == int(par[2]) and int(b) == int(par[3]) and float(d) >= 0.5:
                rec = SeqRecord(rio.seq, id=pir[0] + '-' + c + '-' + str(d), description=pir[1])
                if rec.id not in ll:
                    ll.append(rec.id)
                    o.write(rec.format('fasta'))
                    ccc += 1
                    # print(rec.id)
                # else:
                # print(rec.id)
    o.close()
    print('Secuencias finales ', ccc)