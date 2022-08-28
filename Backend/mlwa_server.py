from flask import Flask, request
from flask_cors import CORS
import ML_model
import time
    
#---------------------------For the variables-------------------------------
# To read the DNA sequence
dna_seq = ""
with open("protein_dna.dat", "r") as f:
    for line in f.readlines():
        if not (">" in line):
            dna_seq = line.upper() #In capial letters

# To read the Amino Acid sequence
aa_seq = ""
with open("protein_aa.dat", "r") as g:
    for line in g.readlines():
        if not (">" in line):
            aa_seq = line.upper() #In capital letters
            
# To translate DNA into Amino Acids
translator = {"ATT": "I", "ACT": "T", "AAT": "N", "AGT": "S",
              "ATC": "I", "ACC": "T", "AAC": "N", "AGC": "S",
              "ATA": "I", "ACA": "T", "AAA": "K", "AGA": "R",
              "ATG": "M", "ACG": "T", "AAG": "K", "AGG": "R",
              "CTT": "L", "CCT": "P", "CAT": "H", "CGT": "R",
              "CTC": "L", "CCC": "P", "CAC": "H", "CGC": "R",
              "CTA": "L", "CCA": "P", "CAA": "Q", "CGA": "R",
              "CTG": "L", "CCG": "P", "CAG": "Q", "CGG": "R",
              "GTT": "V", "GCT": "A", "GAT": "D", "GGT": "G",
              "GTC": "V", "GCC": "A", "GAC": "D", "GGC": "G",
              "GTA": "V", "GCA": "A", "GAA": "E", "GGA": "G",
              "GTG": "V", "GCG": "A", "GAG": "E", "GGG": "G",
              "TTT": "F", "TCT": "S", "TAT": "Y", "TGT": "C",
              "TTC": "F", "TCC": "S", "TAC": "Y", "TGC": "C",
              "TTA": "L", "TCA": "S", "TAA": "-", "TGA": "-",
              "TTG": "L", "TCG": "S", "TAG": "-", "TGG": "W"
              }

# List of the three letter code of the amino acids
aa_full_list = ["G","A","L","M","F","W","K","Q","E","S","P","V","I","C","Y","H","R","N","D","T"]

# To create instance of Flask server
app = Flask(__name__)

# Cross Origin Resource Sharing
CORS(app)

#---------------------------Methods-------------------------

# To extract a code from an amino acid sequence
def aa_seq_2_code(seq):
    if (len(seq)) != len(aa_seq):
        raise Exception("There has been an error with the DNA to AA sequence.")
    code = ["M", 1, "M"]
    for i in range(len(seq)):
        if seq[i] != aa_seq[i]:
            code = [seq[i], i + 1, aa_seq[i]]
    return code
    
# DNA utno Amino Acids
def ribosome(dna_s):
    if len(dna_s) % 3 == 0:
        aa_s = []
        for i in range(0, len(dna_s), 3):
            try:
                codon = dna_s[i:i+3]
                aa_s.append(translator[codon])
            except Exception as e:
                print(e)
                print("An error was found in the DNA sequence")
                print(f"Please check codon {dna_s[i:i+3]} at position {i}")
                return [False, ""]
        aa_ss = "".join(aa_s)
        return [True, aa_ss]
    else:
        print("An error was found in the DNA sequence")
        print("The length of the DNA sequence has to be a multiple of 3")
        print(f"Current sequence has a length of {len(dna_s)} nucleotides")
        return [False, ""]

# ------------------------- For DNA --------------------------
# To create a mutated sequence of amino acids based on DNA    
def dna_mutator(code):
    # If the code's length (long or short) is not the one that it should be
    if((len(code) >= 3) and (len(code) <= 5)):
        primera = code[0]
        ultima = code[-1]
        # If the code contains a letter which is not a nucleotide
        #Or if both letters are the same
        if ((primera in "ATGC") and (ultima in "ATGC") and (primera != ultima)):
            numero = int(code[1:-1])
            #If the code's number is out of range
            if ((numero >= 1) and (numero <= len (dna_seq))):
                # If thefirst letter position of the code doesn't match
                if (dna_seq[numero -1 == primera]):
                    new_seq_list = list(dna_seq)
                    new_seq_list[numero - 1] = ultima
                    dna_seq_string = "".join(new_seq_list)
                    translate_ok, aa_seq_string = ribosome(dna_seq_string)
                    if translate_ok:
                        new_aa, new_code, new_mut_aa = aa_seq_2_code(aa_seq_string)
                        predicted_mutation = ML_model.get_prediction(new_aa, new_code + 1, new_mut_aa)
                        return {"pathogenicity": predicted_mutation["pathogenicity"],
                                "percent" : predicted_mutation ["percent"],
                                "code": predicted_mutation ["code"],
                                "model": True,
                                "sequence": aa_seq_string}
                    else:
                        print("An error occured during the translation :(")
                        return {"pathogenicity": False,
                                "percent" : 0,
                                "code": code,
                                "model": False,
                                "sequence": ""}
                else:
                    print("The letter from the code and the sequence at the provided position doesn't match")
                    return {"pathogenicity": False,
                            "percent" : 0,
                            "code": code,
                            "model": False,
                            "sequence": ""}
            else:
                print("The provided number is out of range od the sequence's length")
                return {"pathogenicity": False,
                        "percent" : 0,
                        "code": code,
                        "model": False,
                        "sequence": ""}
        else:
            print("One/Both of the provided letters don't belong to the sequence")
            return {"pathogenicity": False,
                    "percent" : 0,
                    "code": code,
                    "model": False,
                    "sequence": ""}
    else:
        print("The code's length have an error")
        return {"pathogenicity": False,
                "percent" : 0,
                "code": code,
                "model": False,
                "sequence": ""}
#---------------------- For Amino Acids----------------------    
# To create a mutated sequence of amino acids    
def aa_mutator(code):
    # If the code's length (long or short) is not the one that it should be
    if((len(code) >= 3) and (len(code) <= 5)):
        primera = code[0]
        ultima = code[-1]
        # If the code contains a letter which is not an aa
        #Or if both letters are the same
        if ((primera in aa_full_list) and (ultima in aa_full_list) and (primera != ultima)):
            numero = int(code[1:-1])
            #If the code's number is out of range
            if ((numero >= 1) and (numero <= len (aa_seq))):
                # If thefirst letter position of the code doesn't match
                if (aa_seq[numero -1 == primera]):
                    new_seq_list = list(aa_seq)
                    new_seq_list[numero - 1] = ultima
                    predicted_mutation = ML_model.get_prediction(primera, numero + 1, ultima)
                    return {"pathogenicity": predicted_mutation["pathogenicity"],
                                "percent" : predicted_mutation ["percent"],
                                "code": predicted_mutation ["code"],
                                "model": True,
                                "sequence": "".join(new_seq_list)}
                else:
                    print("The letter from the code and the sequence at the provided position doesn't match")
                    return {"pathogenicity": False,
                            "percent" : 0,
                            "code": code,
                            "model": False,
                            "sequence": ""}
            else:
                print("The provided number is out of range od the sequence's length")
                return {"pathogenicity": False,
                        "percent" : 0,
                        "code": code,
                        "model": False,
                        "sequence": ""}
        else:
            print("One/Both of the provided letters don't belong to the sequence")
            return {"pathogenicity": False,
                    "percent" : 0,
                    "code": code,
                    "model": False,
                    "sequence": ""}
    else:
        print("The code's length have an error")
        return {"pathogenicity": False,
                "percent" : 0,
                "code": code,
                "model": False,
                "sequence": ""}        
    
#------------------------Routes------------------------------------

#Ruta para la URL terminada en "/"
@app.route("/")
def hello_world():
    return "<p>Hello, World!</p>"

#Ruta para la URL terminada en "/dnaseq"
@app.route("/dnaseq")
def dna_sequence():
    return dna_seq

#Ruta para la URL terminada en "/aaseq"
@app.route("/aaseq")
def aa_sequence():
    return aa_seq
# Route for the URL /dnamutate
@app.route("/dnamutate", methods=['GET'])
def dna_mutation():
    if request.method == 'GET':
        if "code" in request.args.keys():
            dna_code = request.args["code"]
            mutated_dna_seq = dna_mutator(dna_code)
            time.sleep(3) #To make the frontend wait for 3 seconds
            if len(mutated_dna_seq["sequence"]) > 0:
                return mutated_dna_seq
            else: 
                {"pathogenicity": False, "percent" : 0, "code": dna_code, "model": False, "sequence": mutated_dna_seq}
        else:
            {"pathogenicity": False, "percent" : 0, "code": dna_code, "model": False, "sequence": mutated_dna_seq}
    else:
        {"pathogenicity": False, "percent" : 0, "code": dna_code, "model": False, "sequence": mutated_dna_seq}
 
# Route for the URL /aamutate    
@app.route("/aamutate", methods=['GET'])
def aa_mutation():
    if request.method == 'GET':
        if "code" in request.args.keys():
            aa_code = request.args["code"]
            mutated_aa_seq = aa_mutator(aa_code)
            time.sleep(3) #To make the frontend wait for 3 seconds
            if len(mutated_aa_seq) > 0:
                return mutated_aa_seq
            else:
                return {"pathogenicity": False, "percent" : 0, "code": aa_code, "model": False, "sequence": mutated_aa_seq}        
        else:
            return {"pathogenicity": False, "percent" : 0, "code": aa_code, "model": False, "sequence": mutated_aa_seq}
    else:
        return {"pathogenicity": False, "percent" : 0, "code": aa_code, "model": False, "sequence": mutated_aa_seq}

#------------------------------ For the main program------------------------------
    
if __name__ == "__main__":
    app.run(debug=True)