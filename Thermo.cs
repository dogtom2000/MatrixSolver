using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MatrixSolver
{
    struct DataStructure
    {
        public Dictionary<string, int> components;
        public double molwt;
        public double[] coefHigh;
        public double[] coefLow;
    }

    class Thermo
    {
        string[][] Hdata = new string[][] {
            new string[] {"H                 ",	"L 5/93",	"H ",	" 1",	"  ",	" 0",	"  ",	" 0",	"  ",	" 0",	"G ",	"200.000",	"6000.000",	"      1.00794"},
            new string[] {" 2.50000286E+00",	"-5.65334214E-09",	" 3.63251723E-12",	"-9.19949720E-16",	" 7.95260746E-20",	"",	"",	"",	"",	"",	"",	"",	"",	""},
            new string[] {" 2.54736589E+04",	"-4.46698494E-01",	" 2.50000000E+00",	" 0.00000000E+00",	" 0.00000000E+00",	"",	"",	"",	"",	"",	"",	"",	"",	""},
            new string[] {" 0.00000000E+00",	" 0.00000000E+00",	" 2.54736599E+04",	"-4.46682853E-01",	" 2.62190349E+04",	"",	"",	"",	"",	"",	"",	"",	"",	""},
        };

        string[][] H2data = new string[][] {
            new string[] {"H2                ",	"TPIS78",	"H ",	" 2",	"  ",	" 0",	"  ",	" 0",	"  ",	" 0",	"G ",	"200.000",	"6000.000",	"      2.01588"},
            new string[] {" 2.93286579E+00",	" 8.26607967E-04",	"-1.46402335E-07",	" 1.54100359E-11",	"-6.88804432E-16",	"",	"",	"",	"",	"",	"",	"",	"",	""},
            new string[] {"-8.13065597E+02",	"-1.02432887E+00",	" 2.34433112E+00",	" 7.98052075E-03",	"-1.94781510E-05",	"",	"",	"",	"",	"",	"",	"",	"",	""},
            new string[] {" 2.01572094E-08",	"-7.37611761E-12",	"-9.17935173E+02",	" 6.83010238E-01",	" 0.00000000E+00",	"",	"",	"",	"",	"",	"",	"",	"",	""},
        };

        string[][] H2Odata = new string[][] {
            new string[] {"H2O               ", "L 8/89",   "H ",   " 2",   "O ",   " 1",   "  ",   " 0",   "  ",   " 0",   "G ",   "200.000",  "6000.000", "     18.01528"},
            new string[] {" 2.67703787E+00",    " 2.97318329E-03",  "-7.73769690E-07",  " 9.44336689E-11",  "-4.26900959E-15",  "", "", "", "", "", "", "", "", ""},
            new string[] {"-2.98858938E+04",    " 6.88255571E+00",  " 4.19864056E+00",  "-2.03643410E-03",  " 6.52040211E-06",  "", "", "", "", "", "", "", "", ""},
            new string[] {"-5.48797062E-09",    " 1.77197817E-12",  "-3.02937267E+04",  "-8.49032208E-01",  "-2.90848168E+04",  "", "", "", "", "", "", "", "", ""},
        };

        string[][] Odata = new string[][] {
            new string[] {"O                 ",	"L 1/90",	"O ",	" 1",	"  ",	" 0",	"  ",	" 0",	"  ",	" 0",	"G ",	"200.000",	"6000.000",	"     15.99940"},
            new string[] {" 2.54363697E+00",	"-2.73162486E-05",	"-4.19029520E-09",	" 4.95481845E-12",	"-4.79553694E-16",	"",	"",	"",	"",	"",	"",	"",	"",	""},
            new string[] {" 2.92260120E+04",	" 4.92229457E+00",	" 3.16826710E+00",	"-3.27931884E-03",	" 6.64306396E-06",	"",	"",	"",	"",	"",	"",	"",	"",	""},
            new string[] {"-6.12806624E-09",	" 2.11265971E-12",	" 2.91222592E+04",	" 2.05193346E+00",	" 2.99687009E+04",	"",	"",	"",	"",	"",	"",	"",	"",	""},
        };

        string[][] O2data = new string[][] {
            new string[] {"O2                ",	"TPIS89",	"O ",	" 2",	"  ",	" 0",	"  ",	" 0",	"  ",	" 0",	"G ",	"200.000",	"6000.000",	"     31.99880"},
            new string[] {" 3.66096083E+00",	" 6.56365523E-04",	"-1.41149485E-07",	" 2.05797658E-11",	"-1.29913248E-15",	"",	"",	"",	"",	"",	"",	"",	"",	""},
            new string[] {"-1.21597725E+03",	" 3.41536184E+00",	" 3.78245636E+00",	"-2.99673415E-03",	" 9.84730200E-06",	"",	"",	"",	"",	"",	"",	"",	"",	""},
            new string[] {"-9.68129508E-09",	" 3.24372836E-12",	"-1.06394356E+03",	" 3.65767573E+00",	" 0.00000000E+00",	"",	"",	"",	"",	"",	"",	"",	"",	""},
        };

        public Dictionary<string, DataStructure> BuildSpeciesData()
        {
            var speciesDictionary = new Dictionary<string, DataStructure>();
            speciesDictionary.Add("H", Build(Hdata));
            speciesDictionary.Add("H2", Build(H2data));
            speciesDictionary.Add("H2O", Build(H2Odata));
            speciesDictionary.Add("O", Build(Odata));
            speciesDictionary.Add("O2", Build(O2data));

            return speciesDictionary;
        }
        
        public DataStructure Build(string[][] rawData)
        {
            var newData = new DataStructure();

            newData.components = new Dictionary<string, int>();
            for (var i = 2; i < 10; i += 2)
            {
                if (rawData[0][i].Trim() != "")
                {
                    newData.components.Add(rawData[0][i].Trim(), Convert.ToInt32(rawData[0][i + 1]));
                }
            }
            newData.molwt = Convert.ToDouble(rawData[0][13]);
            newData.coefHigh = new double[] {  Convert.ToDouble(rawData[1][0]),
                                               Convert.ToDouble(rawData[1][1]), 
                                               Convert.ToDouble(rawData[1][2]), 
                                               Convert.ToDouble(rawData[1][3]), 
                                               Convert.ToDouble(rawData[1][4]), 
                                               Convert.ToDouble(rawData[2][0]), 
                                               Convert.ToDouble(rawData[2][1])};
            newData.coefLow = new double[] {   Convert.ToDouble(rawData[2][2]),
                                               Convert.ToDouble(rawData[2][3]),
                                               Convert.ToDouble(rawData[2][4]),
                                               Convert.ToDouble(rawData[3][0]),
                                               Convert.ToDouble(rawData[3][1]),
                                               Convert.ToDouble(rawData[3][2]),
                                               Convert.ToDouble(rawData[3][3])};

            return newData;
        }

    }
}
