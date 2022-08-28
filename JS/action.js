const { createApp } = Vue

createApp({
    data() {
        return {
            endpoint: "http://127.0.0.1:5000/",
            dnaSequence: ["A"],
            aaSequence: ["G"],
            dnaPostSequence: ["A"],
            aaPostSequence: ["G"],
            aaMut: {
                "G": ["Glycine", "Glycine.png"],
                "A": ["Alanine", "Alanine.png"],
                "L": ["Leucine", "Leucine.png"],
                "M": ["Methionine", "Methionine.png"],
                "F": ["Phenylalanine", "Phenylalanine.png"],
                "W": ["Tryptophan", "Tryptophan.png"],
                "K": ["Lysine", "Lysine.png"],
                "Q": ["Glutamine", "Glutamine.png"],
                "E": ["Glutamic Acid", "GlutamicAcid.png"],
                "S": ["Serine", "Serine.png"],
                "P": ["Proline", "Proline.png"],
                "V": ["Valine", "Valine.png"],
                "I": ["Isoleucine", "Isoleucine.png"],
                "C": ["Cysteine", "Cysteine.png"],
                "Y": ["Tyrosine", "Tyrosine.png"],
                "H": ["Histidine", "Histidine.png"],
                "R": ["Arginine", "Arginine.png"],
                "N": ["Aspargine", "Aspargine.png"],
                "D": ["Aspartic Acid", "AsparticAcid.png"],
                "T": ["Threonine", "Threonine.png"],
            },
            dnaMut: {
                "A": ["Adenine", "adenine.png"],
                "T": ["Thymine", "thymine.png"],
                "G": ["Guanine", "guanine.png"],
                "C": ["Cytosine", "cytosine.png"]
            },
            dna: 0,
            aa: 0,
            dnaMutCode: "",
            aaMutCode: "",
            patho: "",
            perce: "",
            model: "",
            dnaSlider: false,
            aaSlider: false,
            aaImgs: {},
            dnaImgs: {},
            disableControls: false,
            aboutShow: false,
            howtoShow: true,
            dnaShow: true,
            aaShow: false,
            spinnerShow: false,
            resultsShow: false,
            message: 'Hello World!'
        }
    },
    methods: {
        show (something) {
            this.aboutShow = false;
            this.howtoShow = false;
            this.dnaShow = false;
            this.aaShow = false;
            this. spinnerShow = false;
            this.resultsShow = false;
            if ((something == "dnaShow") || (something == "aaShow")) {
                this.howtoShow = true;
            }
            this[something] = true;
        },
        //Method to actualize the mutation code everytime the slider changes DNA
        updateDNAMut() {
            this.dnaMutCode = this.dnaSequence[this.dna] + (+this.dna + 1) + this.dnaPostSequence[this.dna];
        },
        //Method to actualize the mutation code everytime the mutation changes DNA
        onDNA(nam) {
            this.dnaPostSequence[this.dna] = nam;
            this.dnaMutCode = this.dnaSequence[this.dna] + (+this.dna + 1) + nam;
            var difference = 0;
            for (var i = 0; i < this.dnaSequence.length; i++) {
                if (this.dnaSequence[i] != this.dnaPostSequence[i]) {
                    difference = difference + 1;
                }
            }
            if (difference == 0) {
                this.dnaSlider = false;
            } else {
                this.dnaSlider = true;
            }
        },
        //Method to actualize the mutation code everytime the slider changes AA
        updateAAMut() {
            this.aaMutCode = this.aaSequence[this.aa] + (+this.aa + 1) + this.aaPostSequence[this.aa];
        },
        //Method to actualize the mutation code everytime the mutation changes AA
        onAA(nam) {
            this.aaPostSequence[this.aa] = nam;
            this.aaMutCode = this.aaSequence[this.aa] + (+this.aa + 1) + nam;
            var difference = 0;
            for (var i = 0; i < this.aaSequence.length; i++) {
                if (this.aaSequence[i] != this.aaPostSequence[i]) {
                    difference = difference + 1;
                }
            }
            if (difference == 0) {
                this.aaSlider = false;
            } else {
                this.aaSlider = true;
            }
        },
        //Method to send the code of the mutation to the DNA-backend
        runDNA() {
            //If the code's length is different from the original (shorter or longer)
            if ((this.dnaMutCode.length >=3) && (this.dnaMutCode.length<= 5)) {
                var primera = this.dnaMutCode.slice(0,1); //First letter of the code
                var ultima = this.dnaMutCode.slice(-1); //Last letter of the code
                //If any of the letters in the code does not correspond or if they are the same
                if (Object.keys(this.dnaMut), includes(primera) && Object.keys(this.dnaMut), includes(ultima) && (primera != ultima)) {
                    var numero = this.dnaMutCode.slice(1,-1)
                    //If the code's number is bigger or smaller from the one original in the code
                    if ((+numero >= 1) && (+numero <= this.dnaSequence.length)) {
                        this.disableControls = true;
                        this.spinnerShow = true;
                        axios
                            .get(this.endpoint + 'dnamutate',  //To send the input into a new route
                            {
                                params: {
                                code: this.dnaMutCode  //Include code as an argument
                                }
                            })
                            .then(response => {
                                console.log(response.data);  //If everything works correctly
                                this.spinnerShow = false;
                                if (response.data["model"]) {
                                    this.patho = response.data["pathogenicity"];
                                    this.perce = response.data["percent"];
                                    this.model = response.data["model"];
                                    this.resultsShow = true;
                                }
                            })
                            .catch(error => {
                                this.spinnerShow = false;
                                console.log(error);  //If something does not work correctly
                            })
                        } else {
                            console.log("Error in the number")
                        }
                    } else {
                        console.log("Error in the letters")
                    }
                } else {
                    console.log("Error in the longitude")
                }        
        },
        //Method to send the code of the mutation to the DNA-backend
        runAA() {
            //If the code's length is different from the original (shorter or longer)
            if ((this.aaMutCode.length >=3) && (this.aaMutCode.length<= 5)) {
                var primera = this.aaMutCode.slice(0,1); //First letter of the code
                var ultima = this.aaMutCode.slice(-1); //Last letter of the code
                //If any of the letters in the code does not correspond or if they are the same
                if (Object.keys(this.aaMut), includes(primera) && Object.keys(this.aaMut), includes(ultima) && (primera != ultima)) {
                    var numero = this.aaMutCode.slice(1,-1)
                    //If the code's number is bigger or smaller from the one original in the code
                    if ((+numero >= 1) && (+numero <= this.aaSequence.length)) {
                        this.disableControls = true;
                        this.spinnerShow = true;
                        axios
                            .get(this.endpoint + 'aamutate',
                            {
                                params: {
                                code: this.aaMutCode
                            }
                            })
                            .then(response => {
                                console.log(response.data);
                                this.spinnerShow = false;
                                if (response.data["model"]) {
                                    this.patho = response.data["pathogenicity"];
                                    this.perce = response.data["percent"];
                                    this.model = response.data["model"];
                                    this.resultsShow = true;
                                }
                            })
                            .catch(error => {
                                this.spinnerShow = false;
                                console.log(error);
                            })
                    } else {
                        console.log("Error in the number")
                    }
                } else {
                    console.log("Error in the letters")
                }
            } else {
                console.log("Error in the longitude")
            }
        }
    },
    //This will run before running the complete application
    beforeMount () {
        const getDNASequence = async () => {
            await axios
                .get(this.endpoint + 'dnaseq')
                .then(response => {
                    this.dnaSequence = response.data.split("");
                    this.dnaPostSequence = response.data.split("");
                })
                .catch(error => {
                    this.spinnerShow = false;
                    console.log(error);
                });
        }
        const getAASequence = async () => {
            await axios
                .get(this.endpoint + 'aaseq')
                .then(response => {
                    this.aaSequence = response.data.split("");
                    this.aaPostSequence = response.data.split("");
                })
                .catch(error => {
                    console.log(error);
                });
        }
        getDNASequence();
        getAASequence();
    },
    // To be done after the app has been finished
    mounted() {
        Object.keys(this.aaMut).forEach(element => {
            this.aaImgs[element] = new Image();
            this.aaImgs[element].src = "IMG/" + this.aaMut[element][1];
        });
        Object.keys(this.dnaMut).forEach(element => {
            this.dnaImgs[element] = new Image();
            this.dnaImgs[element].src = "IMG/" + this.dnaMut[element][1];
        });
    }
}).mount('#app')