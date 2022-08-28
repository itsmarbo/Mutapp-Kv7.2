# Sólo código necesario
# Este es el código principal que se anclará a la web. Además de este, existen otros dos scripts auxiliares uno con funciones
# relativas al funcionamiento de este programa y otro que sirve para hacer la traduccion de adn a proteina.


# Se importan las librerías y paquetes necesarios para el algoritmo.
import utils
import pandas as pd

from sklearn.svm import SVC
from sklearn.utils import resample
from sklearn.pipeline import Pipeline
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import VotingClassifier


#Carga de base de datos que sirve para entrenar al algoritmo
wdf = pd.read_csv('corrected_FINAL_DDBB.csv', header = "infer")

#Se realizan una serie de acciones que permiten acondicionar la base de datos para los propósitos de este algoritmo

x1=list(wdf.loc[(wdf["label"] == 1)]["residue_conserv"])
l = []

for i in x1:
    l.append(i)
    
# Se localizan los outliers
new_l = sorted(l)[9:]

# Se eliminan los outliers de los datos originales.
wdf = wdf.drop(wdf.loc[(wdf["label"] == 1) & (wdf["residue_conserv"] <= 0.6197)].index)

# Se buscan valores duplicados
wdf[wdf["mutation"].duplicated()]

# Se buscan datos faltantes
wdf.isnull().sum()


# Los datos vienen etiquetados en dos clases, 0 y 1. Se dividen según a qué clase pertenezcan. La clase 0 pertenece a las
# mutaciones benignas y la clase 1 a las patógenas.
y_be = (wdf.values[:,-1] == 0)
y_pa = (wdf.values[:,-1] == 1)

# En Machine Learning existen dos tipos de datos, los datos de entrenamiento que sirven para entrenar el algoritmo y los datos
# de testeo que sirven para ver si el algoritmo ha aprendido correctamente y cumple con el objetivo esperado.
# Es por ello que se tienen que dividir los datos de entrada en dos conjuntos, 'training set' y 'test set'. En este caso,
# el 75% de los datos corresponderá al 'training set' y el 25% restante al 'test set'.

X = wdf.values[:,2:-1]
y = wdf.values[:,-1].astype('int')

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.25, random_state=1012)

# Se crean dos DataFrame, para el conjunto de datos de entrenamiento. En uno se guardan los datos y variables iniciales
# (X_train_dfos), y en el otro la respuesta correspondiente a esos datos (y_train_dfos). El algoritmo deberá aprender la
# relación entre las variables y la respuesta para poder predecir sobre datos cuya respuesta es desconocida.

X_train_dfos = pd.DataFrame(X_train, columns = ['initial_aa', 
                                                'final_aa', 
                                                'topological_domain', 
                                                'functional_domain', 
                                                'd_size',
                                                'd_hf',
                                                'd_vol',
                                                'd_msa',
                                                'd_charge', 
                                                'd_pol', 
                                                'd_aro', 
                                                'residue_conserv',
                                                'secondary_str',
                                                'pLDDT',
                                                'str_pos',
                                                'MTR'])

y_train_dfos = pd.DataFrame(y_train, columns = ['label'])

# Dentro de los datos de entrenamiento, existen dos categorías, los datos de etiqueta 0 y de etiqueta 1. Los datos están
# descompensados, hay muchos mas datos de la clase 1 que de la 0. Hay que balancearlos, para ello se realiza un oversampling,
# que consiste en producir muestras aleatorias pertenecientes a la clase 0 para que haya más y el algoritmo pueda aprender bien.

# Se unen todos los datos de entrenamiento en un único DataFrame
dfos = pd.concat([X_train_dfos, y_train_dfos], axis=1)

# Se separan la clase mayritaria (la del 1) y la minoritaria (la del 0).
df_majority = dfos[dfos.label==1]
df_minority = dfos[dfos.label==0]

# Se realizan las nuevas muestras de la clase minoritaria
df_minority_upsampled = resample(df_minority, 
                                 replace=True,     # sample with replacement
                                 n_samples=len(df_majority)//2,    
                                 random_state=0)   # reproducible results
 
# Se unen la clase mayoritaria y la nueva clase minoritaria en un mismo DataFrame
df_upsampled = pd.concat([df_majority, df_minority_upsampled])

# Se vuelve a separar este nuevo conjunto de datos, en variables y resultados
X_train_os = df_upsampled.values[:,:-1]
y_train_os = df_upsampled.values[:,-1].astype("int")

# Existen una serie de variables que son de tipo cualitativ, es decir que no son numéricas. Para poder tratarlas computacional-
# mente hace falta escribirlas de forma numérica. Para ello se usan una serie de algoritmos que sirven para codificar estos
# datos.

X_train_enc, X_test_enc, X_train_df, X_test_df = utils.categorical_encoding(X_train_os, X_test)
X_train_df.columns 

# Hasta ahora solo se ha llevado a cabo el tratamiento de los datos. Ahora los datos ya están listos para ser procesados por el
# algoritmo. En este modelo se utiliza un algoritmo 'ensemble', que combina tres algoritmos de clasificación distintos y 
# después predice la respuesta dada por la mayoría de ellos. En este algoritmo se utiliza el algoritmo 'Voting Classifier',
# con el modelo de voto 'Soft Voting'. Dentro de este se encuentran los algoritmos 'Logistic Regression L2', 'Support Vector
# Classifier' y 'Random Forest'.

# Primero se seleccionan las variables mas importantes 
X_train_fs, X_test_fs, featEn, posEn = utils.select_features(X_train_enc,
                                                             X_train_df, 
                                                             y_train_os,
                                                             X_test_enc, 
                                                             n = 45)


# Se implementan los algoritmos mencionados
pipeline_ensemble_soft = Pipeline( [("scaler", StandardScaler()),
                                    ("Ensemble_soft", VotingClassifier( voting = "soft",
                                                                        weights = [1,0.5,1.75],
                                                                        estimators=[("logistic", LogisticRegression(solver = "saga",
                                                                                    penalty = "l2",
                                                                                    max_iter = 10000,
                                                                                    class_weight = {0: 3, 1: 2},
                                                                                    multi_class = "ovr",
                                                                                    C = 2.91,
                                                                                    random_state = 8)),
                                                                                ("SVC", SVC(kernel = "linear", 
                                                                                    class_weight= {0:1, 1:1},
                                                                                    probability=True,
                                                                                    decision_function_shape = "ovr",
                                                                                    degree = 2,
                                                                                    gamma = "auto", 
                                                                                    C = 1,
                                                                                    random_state = 45)),
                                                                                ("RF", RandomForestClassifier(max_depth = 3,
                                                                                    criterion = "log_loss",
                                                                                    max_features = "log2",
                                                                                    oob_score = False,
                                                                                    min_samples_split = 2, # min = 5
                                                                                    class_weight= {0:3, 1:1},
                                                                                    random_state = 45))]))])
pipeline_ensemble_soft.fit(X_train_fs, y_train_os)



# Ahora viene el código correspondiente a la introducción de datos y predicción del algoritmo

def get_prediction(aa, aa_code, aa_mut):
    # Leer secuencia de aminoacidos de archivo
    aa_seq = ""
    with open("protein_aa.dat", "r") as g:
        for line in g.readlines():
            if not (">" in line):
                aa_seq = line.upper()

    # Se crea un DataFrame con la mutación introducida. Esto es así ya que las funciones que se utilizan a continuación están 
    # definidas de tal manera que toman como entrada un DataFrame
    conflictive = pd.DataFrame(columns = ["Mutationppt"])
    c = pd.DataFrame({'Mutationppt':[f"{aa}{aa_code}{aa_mut}"]})
    challenge = pd.concat([conflictive, c])

    # Se utiliza la función siguiente para generar una tabla con todos los descriptores que después se utilizan en el 
    # algoritmo.
    ch_df = utils.KCNQ2_DDBB_generation(challenge)

    # Se crea una lista donde aparece únicamente el nombre de la mutación
    variants_names = list(ch_df["Mutationppt"])

    #Se adapta la tabla a los requerimientos del algoritmo
    ch_df = utils.preprocessing_ch(ch_df)

    #Se convierte la tabla a array numérico para que pueda introducirse en el algoritmo
    X_ch = ch_df.to_numpy()
    X_train_enc, X_ch_enc, X_train_df, X_ch_df = utils.categorical_encoding(X_train_os, X_ch)


    # Se seleccionan los desciptores más importantes como en el entrenamiento
    X_train_fs, X_ch_enc_fs, feat_pred, pos_pred = utils.select_features(X_train_enc,
                                                                     X_train_df, 
                                                                     y_train_os,
                                                                     X_ch_enc,
                                                                     n = 45)

    # Se reliza la predicción de la mutación introducida
    # En la predicción se obtienen dos tipos de resultados. Uno te dice a qué clase pertenece la mutación y el otro la probabilidad 
    # de que la mutación pertenezca a cada clase.
    KCNQ2e_y_ch_p = pipeline_ensemble_soft.predict(X_ch_enc_fs)

    # View probabilities in prediction
    KCNQ2eprob = pipeline_ensemble_soft.predict_proba(X_ch_enc_fs)

    # Se reliza la predicción de la mutación introducida
    # En la predicción se obtienen dos tipos de resultados. Uno te dice a qué clase pertenece la mutación y el otro la probabilidad 
    # de que la mutación pertenezca a cada clase.
    KCNQ2e_y_ch_p = pipeline_ensemble_soft.predict(X_ch_enc_fs)

    # View probabilities in prediction
    KCNQ2eprob = pipeline_ensemble_soft.predict_proba(X_ch_enc_fs)

    # Se saca el mensaje de la predicción
    if KCNQ2e_y_ch_p == 0:
        mut_type = "Benign"
        proba = KCNQ2eprob[0,0]*100
    else:
        mut_type = "Pathogenic"
        proba = KCNQ2eprob[0,1]*100

    if proba <= 60:
        succ_rate = "VERY LOW"
    elif 60 < proba <= 70:
        succ_rate = "LOW"
    elif 70 < proba <= 80:
        succ_rate = "MODERATE"
    elif 80 < proba <= 90:
        succ_rate = "HIGH"
    elif 90 < proba:
        succ_rate = "VERY HIGH"

    message = f"{mut_type} mutation! Success rate: {proba:.4f}%"

    return {"pathogenicity": bool(KCNQ2e_y_ch_p), "percent": proba, "code": f"{aa}{aa_code}{aa_mut}", "message": message}
