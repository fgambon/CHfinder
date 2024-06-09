from __future__ import absolute_import, division, print_function
import itertools
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import re
# import cPickle as pickle
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast import NCBIWWW
from time import time
import os
import tensorflow as tf
from tensorflow import keras
import numpy as np


print(tf.__version__)
'''
path = '/Volumes/rap/igls/peces/vgp'
files = []
# r=root, d=directories, f = files
for r, d, f in os.walk(path):
    for file in f:
        if file in est:
                files.append(os.path.join(r, file))
                print(file)
print('archivos ',len(files))
'''
files = ['/Users/franciscogambondeza/.mounty/rap/igls/reptiles/ncbi/Testudines/Dermochelys_coriacea/GCA_009764565.4_rDerCor1.pri.v4_genomic.fna']
for file in files:
    corte = file.split('/')
    animal = corte[-2]
    especie = file.split(corte[-1])[0] + 'chs'
    # especie ='igls_tcrs'
    if os.path.exists(especie) is False:
        os.mkdir(especie)
    print(animal, file)
    qquery = '/Users/franciscogambondeza/.mounty/rap/reposteria/aaa/reptiles/query-rep-chs.txt'

    from Bio.Blast.Applications import NcbiblastxCommandline

    outfile = especie + "/tabla.txt"
    vvv = 'tblastn'
    blastn_cline = NcbiblastxCommandline(vvv, query=qquery, subject=file, evalue=0.01, out=outfile, outfmt=6)
    print(str(blastn_cline))
    stdout, stderr = blastn_cline()

    aa = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

    bas = especie + '/aexones.fasta'
    nnfile = especie + '/nexones.fasta'
    tablafile = especie + '/tabla.txt'
    # tabla obtenida de internet, se obtienen los datos necesarios
    A = open(especie + '/tabla.txt', 'r')
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
    rep = []
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
                    fin1 = [i.start() for i in re.finditer('GT', str(k))]
                    se = [(i, j) for i, j in itertools.product(ini1, fin1) if
                          j > i + 120 and j < i + 390 and ((np.abs(j - i) - 2) % 3 == 0)]
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
                            if len(p) >= 60:
                                nrec = SeqRecord(exon, id=animal + '-' + ee + '-' + str(x + 4) + '-' + str(y),
                                                 description='plus')
                                arec = SeqRecord(p, id=animal + '-' + ee + '-' + str(x + 4) + '-' + str(y),
                                                 description='plus')
                                if arec.id not in rep:
                                    ofile.write(arec.format('fasta'))
                                    oofile.write(nrec.format('fasta'))
                                    rep.append(arec.id)
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
                    fin1 = [i.start() for i in re.finditer('GT', str(k))]
                    se = [(i, j) for i, j in itertools.product(ini1, fin1) if
                          j > i + 120 and j < i + 390 and ((np.abs(j - i) - 2) % 3 == 0)]
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
                            if len(p) >= 60:
                                aa = ancho - y
                                bb = ancho - x
                                nrec = SeqRecord(exon, id=animal + '-' + ee + '-' + str(aa) + '-' + str(bb),
                                                 description='minus')
                                arec = SeqRecord(p, id=animal + '-' + ee + '-' + str(aa) + '-' + str(bb),
                                                 description='minus')
                                if arec.id not in rep:
                                    ofile.write(arec.format('fasta'))
                                    oofile.write(nrec.format('fasta'))
                                    rep.append(arec.id)
                                    r += 1

            print('reverse squences posibles ', r)

    ofile.close()
    oofile.close()
    aa = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

    # chs = ['back','exon2_I','exon3_I','exon4_I','exon2_DA','exon3_DA','exon2_DB','exon3_DB','exon2_DMA','exon3_DMA','exon2_DMB','exon3_DMB','beta2','exon2_CD1','exon3_CD1','exon4_CD1']

    chs = ['back', 'exon1_M', 'exon2_M', 'exon3_M', 'exon4_M',
           'exon1_D', 'exon2_D', 'exon3_D', 'exon4_D', 'exon5_D', 'exon6_D', 'exon7_D', 'exon8_D', 'exon9_D',
           'exon10_D', 'exon11_D', 'exon1_Y', 'exon2_Y', 'exon3_Y', 'exon4_Y', 'exon1_XA', 'exon2_XA', 'exon3_XA',
           'exon4_XA']

    path = '/Users/franciscogambondeza/.mounty/rap/reposteria/aaa/reptiles/chs'
    x_train = np.load(path + '/x_train.npy')
    y_train = np.load(path + '/y_train.npy')
    x_test = np.load(path + '/x_test.npy')
    y_test = np.load(path + '/y_test.npy')

    model = keras.models.Sequential([
        tf.keras.layers.Flatten(input_shape=(2000,)),
        keras.layers.Dense(128.0, activation='relu'),
        tf.keras.layers.Dropout(0.2),
        keras.layers.Dense(len(chs), activation='softmax')
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
        locus = chs[ale]
        p = ID[a] + '-----------> ' + str(pp) + '-' + chs[ale] + '\n'
        oo.write(p)
        if locus != 'back':
            name = ID[a]
            data = name, locus, pp
            nlp.append(data)
    oo.close()
    print('asignacion hecha')

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
        dupli = []
        for a, b, c, d in R:
            if a not in dupli:
                k = range(a - 400, a + 400)
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
    o = open(especie + '/result-chs.fasta', 'w')

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